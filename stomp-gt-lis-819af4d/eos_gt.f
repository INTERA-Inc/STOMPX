!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE AIRDFL( TX,VISLX,DFLAX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculates binary liquid diffusion coefficient for air in liquid
!     water, using the Wilke and Chang estimation method; pp. 598.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, June, 1993.
!     Last Modified by MD White, Battelle, PNL, June 2, 1993.
!     eos_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/AIRDFL'
      VBA = 2.85D-1*VCRA**1.048D+0
      DFLAX = 7.4D-15*((2.6D+0*WTMW)**5.D-1)*(TX+TABS)/
     &  (VISLX*(VBA**6.D-1))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of AIRDFL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE AIRGSD( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculates component density of air with the ideal gas equation
!     of state.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     eos_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/AIRGSD'
      RHOX = PX/(RCA*(TX+TABS)*0.99966D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of AIRGSD group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE AIRGSH( TX,PX,HX,UX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate the specific enthalpy and internal energy of air as a
!     function of temperature and air partial pressure with a
!     reference state of enthalpy at 0 C.
!
!     Calculate the isobaric specific heat as a function of temperature
!     from an empirical relation. pp. 580. (Sandler)
!
!     Calculate the specific enthalpy by integrating the isobaric
!     specific heat using a two-point Gaussian quadrature method.
!     pp. 101-105. (Carnahan et al.)
!
!     Carnahan, B., H.A. Luther, and J.O. Wilkes. 1969. Applied
!     Numerical Methods. John Wiley and Sons Inc., New York.
!
!     Sandler, S.I. 1989. Chemical and Engineering Thermodynamics.
!     John Wiley and Sons, Inc., New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     eos_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 AX(5),BX(4),CX(4)
!
!----------------------Data Statements---------------------------------!
!
      SAVE AX,BX,CX
      DATA AX /4.184D+3,6.713D+0,4.697D-4,1.147D-6,-4.696D-10/
      DATA BX /11162.5D+0,-12617.8D+0,-3211.64D+0,8374.58D+0/
      DATA CX /7890.68D+0,-8645.83D+0,-2512.88D+0,7377.18D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/AIRGSH'
      PXX = PX
      TK1 = 5.D-1*(-5.7735D-1*TX + TX + 2.D+0*TABS)
      TK2 = 5.D-1*(5.7735D-1*TX + TX + 2.D+0*TABS)
      CP1 = AX(1)*(AX(2)+AX(3)*TK1+AX(4)*(TK1**2)+AX(5)*(TK1**3))/WTMA
      CP2 = AX(1)*(AX(2)+AX(3)*TK2+AX(4)*(TK2**2)+AX(5)*(TK2**3))/WTMA
      HX = TX*(CP1+CP2)/2.D+0
      UX = HX - RCA*(TX+TABS)
!      HX = (BX(1)+BX(2)*EXP(-((TX-BX(3))/BX(4))**2))*1.D+3
!      UX = (CX(1)+CX(2)*EXP(-((TX-CX(3))/CX(4))**2))*1.D+3
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of AIRGSH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE AIRGSC( TX,CPX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate the isobaric specific heat as a function of temperature
!     from an empirical relation. pp. 580. (Sandler)
!
!     Sandler, S.I. 1989. Chemical and Engineering Thermodynamics.
!     John Wiley and Sons, Inc., New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, 4 November 2003.
!     Last Modified by MD White, Battelle, PNNL, 4 November 2003.
!     eos_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 AX(5)
!
!----------------------Data Statements---------------------------------!
!
      SAVE AX
      DATA AX /4.184D+3,6.713D+0,4.697D-4,1.147D-6,-4.696D-10/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/AIRGSC'
      TK = TX+TABS
      CPX = AX(1)*(AX(2)+AX(3)*TK+AX(4)*(TK**2)+AX(5)*(TK**3))/WTMA
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of AIRGSC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE AIRGSK( TX,THKX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate the thermal conductivity as a function of temperature
!     from an empirical relation. pp. 413. (Lide and Kehiaian)
!
!     Lide, D.R., and H.V. Kehiaian.  1994.  CRC Handbook of
!     Thermophysical and Thermochemical Data, CRC Press, Boca Raton.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, 4 November 2003.
!     Last Modified by MD White, Battelle, PNNL, 4 November 2003.
!     eos_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 AX(5)
!
!----------------------Data Statements---------------------------------!
!
      SAVE AX
      DATA AX /0.D+0,0.0965D+0,-9.960D-6,-9.310D-8,8.882D-11/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/AIRGSK'
      TK = TX+TABS
      THKX = 0.D+0
      DO 100 I = 1,5
        THKX = THKX + AX(I)*(TK**(I-1))
  100 CONTINUE
      THKX = 1.D-3*THKX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of AIRGSK group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE AIRGSV( TX,VISX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate the air viscosity from the kinetic theory. pp. 528-530.
!
!     Hirschfelder, J.O., C.F. Curtiss, and R.B. Bird. 1954. Molecular
!     Theory of Gases and Liquids. John Wiley & Sons.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     eos_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/AIRGSV'
      TK = TX+TABS
      VISX = 2.6693D-6*SQRT(WTMA*TK)/(1.5542D+1-(6.88D-3*TK))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of AIRGSD group  ---
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
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     eos_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      TK = TX+TABS
      PB = PX*1.D-5
      TP = TK/SQRT(1.15D+0*TBA*1.15D+0*TBW)
      SIG = 5.D-1*(1.18D+0*(2.85D-1*(VCRW**1.048D+0))**3.333D-1 +
     &  1.18D+0*(2.85D-1*(VCRA**1.048D+0))**3.333D-1 )
      OMG = (A(1)/TP**A(2)) + (A(3)/EXP(A(4)*TP))
     &  + (A(5)/EXP(A(6)*TP))  + (A(7)/EXP(A(8)*TP))
      WTMX = 2.D+0/((1.D+0/WTMA) + (1.D+0/WTMW))
      DIFX = (3.03D+0 - (9.8D-1/SQRT(WTMX)))*1.D-3*(TK**1.5D+0)/
     &  (PB*SQRT(WTMX)*(SIG**2)*(OMG))*1.D-4
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BNDFAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CROUCH
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     The pressurized crack problem, following the displacement
!     discountinuity method of Crouch.
!
!     Crouch, S.L. 1976. Solution of plane elasticity problems by the
!     displacment discountinuity method. International Journal for
!     Numerical Methods in Engineering, 10:301-343.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 March 2015.
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
!----------------------Parameter Statements----------------------------!
!
      PARAMETER(LNDD=40)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AJ(LNDD,LNDD),BJ(LNDD)
      REAL*8 AX(LNDD),XX(LNDD)
      INTEGER IJ(LNDD)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CROUCH'
!
!---  Number of displacment discontinuities  ---
!
      NDD = 40
      PX = 1.D+6
      VX = 0.2D+0
      GX = PX/2.4D-3
      BX = 1.D+0
      DO 10 I = 1,NDD
        AX(I) = BX/REAL(NDD)
        XX(I) = -BX + BX/REAL(NDD) + 2.D+0*REAL(I-1)*BX/REAL(NDD)
        BJ(I) = PX*(REAL(I)/REAL(NDD))*GPI*(1.D+0-VX)/GX
   10 CONTINUE
      DO 30 I = 1,NDD
        DO 20 J = 1,NDD
          AJ(I,J) = AX(J)/(((XX(I)-XX(J))**2)-(AX(J)**2))
   20   CONTINUE
   30 CONTINUE
!
!---  Solve linear system  ---
!
      JP = NDD
      KP = LNDD
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )     
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CROUCH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_B( XLSX,RHOBX,RHOLWX,TX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
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
!     The seawater option uses the formulation of Millero and Huang.
!
!     Millero, F.J., and F. Huang. 2009. The density of seawater as
!     a function of salinity (5 to 70 g/kg) and temperature (273.15
!     to 363.15 K), Ocean Sci. 5:91-100, www.ocean-sci.net/5/91/2009/
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
      REAL*8 CSAX(10)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CCX,CHX,SCX,SAX,VCX,CSAX
      DATA CCX / -3.033405D+0, 10.128163D+0, -8.750567D+0, 2.663107D+0 /
      DATA CHX / -167.219D+0, 448.55D+0, -261.07D+0, -13.644D+0,
     &  13.97D+0, -0.315154D+0, -1.203374D-3, 7.48908D-13,
     &  0.1342489D+0, -3.946963D-3 /
      DATA SCX / -9.9559D+0, 7.0845D+0, 3.9093D+0 /
      DATA SAX / -4.539D-3, -1.638D-4, 2.551D-5 /
      DATA VCX / 3.1975D+0 /
      DATA CSAX / 8.197247E-01,-3.779454E-03,6.821795E-05,-8.009571E-07,
     &  6.158885E-09,-2.001919E-11,-5.808305E-03,5.354872E-05,
     &  -4.714602E-07,5.249266E-04 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_B'
!
!---  Seawater option  ---
!
      IF( ISLC(22).EQ.1 ) THEN
        AX = 0.D+0
        DO M = 1,6
          AX = AX + CSAX(M)*(TX**(M-1))
        ENDDO
        BX = 0.D+0
        DO M = 1,3
          BX = BX + CSAX(M+6)*(TX**(M-1))
        ENDDO
        CX = CSAX(10)
        SRX = (35.16504D+0/35.D+0)*XLSX*1.D+3
        DRHOX = AX*SRX + BX*(SRX**1.5D+0) + CX*(SRX**2)
        RHOBX = RHOLWX + DRHOX
      ELSE
!
!---    Convert mass fraction to molality  ---
!
        GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
!
!---    Convert specific volume units to cm^3/gm  ---
!
        VOX = 1.D+3/RHOLWX
!
!---    Limiting apparent molal volume (cm^3/mol) of NaCl in 
!       solution as the concentration goes to zero  ---
!
        PHIPX = CHX(1) + CHX(2)*VOX + CHX(3)*(VOX**2)
!
!---    Apparent molal volume (cm^3/mol) of NaCl in solution  ---
!
        PHIX = PHIPX + (CHX(4)+CHX(5)*VOX)*((VOX/(VCX-VOX))**2)
     &    *SQRT(GLSX)
!
!---    Brine density (gm/cm^3)  ---
!
        RHOBX = (1.D+3 + GLSX*WTMS)/(1.D+3*VOX + GLSX*PHIX)
!
!---    Convert density units to kg/m^3  ---
!
        RHOBX = 1.D+3*RHOBX
      ENDIF
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
      SUBROUTINE EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,
     &  XMGAX,XMGWX,XMLAX,XMLSX,XMLWX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     This subroutine calculates the equilibrium state between
!     water, air and salt, where brine refers to water + salt without
!     dissolved air.
!
!     Arguments
!
!     XLSX - (in) mass fraction of salt in brine (i.e., salt + water)
!     XLSMX - (in) solubility mass fraction of salt in brine
!     XMLAX - (in) mole fraction of air in aqueous phase
!     XMGAX - (in) mole fraction of air in gas phase
!
!     XGAX - (out) mass fraction of air in gas phase
!     XGWX - (out) mass fraction of water in gas phase
!     XLAX - (out) mass fraction of air in aqueous phase
!     XLSX - (out) mass fraction of salt in aqueous phase
!     XLWX - (out) mass fraction of water in aqueous phase
!     XMGWX - (out) mole fraction of water in gas phase
!     XMLSX - (out) mole fraction of salt in aqueous phase
!     XMLWX - (out) mole fraction of water in aqueous phase
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 June 2010.
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
      SUB_LOG(ISUB_LOG) = '/EQUIL'
      INCG = 1
!
!---  Brine mole fraction of salt   ---
!
      XBSX = XLSX
      XMBSX = (XBSX/WTMS)/(XBSX/WTMS + (1.D+0-XBSX)/WTMW)
!
!---  Aqueous mole fraction of salt  ---
!
      IF( XMLAX/EPSL.GT.EPSL ) THEN
        XMLSX = ((1.D+0/XMLAX)-1.D+0)*XMBSX*XMLAX
      ELSE
        XMLSX = XMBSX
      ENDIF
!
!---  Aqueous mole fraction of water  ---
!
      XMLWX = MAX( 1.D+0-XMLSX-XMLAX,0.D+0 )
!
!---  Aqueous mass fractions  ---
!
      WTMLX = XMLWX*WTMW + XMLAX*WTMA + XMLSX*WTMS
      XLWX = XMLWX*WTMW/WTMLX
      XLAX = XMLAX*WTMA/WTMLX
      XLSX = XMLSX*WTMS/WTMLX
!
!---  Gas mole fraction of water  ---
!
      XMGWX = MAX( 1.D+0-XMGAX,0.D+0 )
!
!---  Gas mass fractions  ---
!
      WTMGX = XMGWX*WTMW + XMGAX*WTMA
      XGWX = XMGWX*WTMW/WTMGX
      XGAX = XMGAX*WTMA/WTMGX
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
      FUNCTION HGBG_GT( BCX,NBS,NBE,M )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Gas hydraulic gradient boundary condition function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
      USE BCVP
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HGBG_GT'
      PGX = PGB(M,NBS) + PATM
      TX = TB(M,NBS)
      TKX = TX + TABS
      IF( TKX.LT.TCRW ) THEN
        INDX = 0
        CALL REGION_4( TX,PSWX,INDX )
        PVWX = PSWX*BCX(5)
      ELSE
        PVWX = PGX*BCX(5)
      ENDIF
      CALL P_IAPWS( TX,PVWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
      PVAX = MAX( PGX-PVWX,0.D+0 )
      CALL AIRGSD( TX,PVAX,RHOGAX )
      XMGAX = PVAX/PGX
      XMGWX = MAX( 1.D+0-XMGAX,0.D+0 )
      WTMGX = XMGWX*WTMW + XMGAX*WTMA
      XGWX = XMGWX*WTMW/WTMGX
      XGAX = XMGAX*WTMA/WTMGX
      RHOGX = XGAX*RHOGAX + XGWX*RHOGWX
      NS = IBCN(NBS)
      NE = IBCN(NBE)
      ISX = ID(NS)
      IEX = ID(NE)
      JSX = JD(NS)
      JEX = JD(NE)
      KSX = KD(NS)
      KEX = KD(NE)
      II = SIGN(1,IEX-ISX)
      JI = SIGN(1,JEX-JSX)
      KI = SIGN(1,KEX-KSX)
!
!---  Convert boundary pressure to node-centroid pressure  ---
!
      IF( IBCD(NBS).EQ.-3 ) THEN
        NPZ = NSZ(NS)
        PGX = PGX + 5.D-1*DZGF(NS)*RHOGX*GRVZ(NPZ)
      ELSEIF( IBCD(NBS).EQ.-2 ) THEN
        NPY = NSY(NS)
        PGX = PGX + 5.D-1*DYGF(NS)*RHOGX*GRVY(NPY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.-1 ) THEN
        NPX = NSX(NS)
        PGX = PGX + 5.D-1*DXGF(NS)*RHOGX*GRVX(NPX)
      ELSEIF( IBCD(NBS).EQ.1 ) THEN
        NQX = NSX(NS)+1
        IF( INBS(4,NS).GT.0 ) NQX = INBS(4,NS)
        PGX = PGX - 5.D-1*DXGF(NS)*RHOGX*GRVX(NQX)
      ELSEIF( IBCD(NBS).EQ.2 ) THEN
        NQY = NSY(NS)+IFLD
        IF( INBS(5,NS).GT.0 ) NQY = INBS(5,NS)
        PGX = PGX - 5.D-1*DYGF(NS)*RHOGX*GRVY(NQY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.3 ) THEN
        NQZ = NSZ(NS)+IJFLD
        IF( INBS(6,NS).GT.0 ) NQZ = INBS(6,NS)
        PGX = PGX - 5.D-1*DZGF(NS)*RHOGX*GRVZ(NQZ)
      ENDIF
!
!---  Loop over nodes in the x direction
!
      I = ISX
      J = JSX
      K = KSX
      NX = ND(I,J,K)
      DO 100 I = ISX,IEX,II
        N = ND(I,J,K)
        NPX = NSX(N)
        GB = (XP(N)-XP(NX))*GRVX(NPX)
        IF( TKX.LT.TCRW ) THEN
          INDX = 0
          CALL REGION_4( TX,PSWX,INDX )
          PVWX = PSWX*BCX(5)
        ELSE
          PVWX = PGX*BCX(5)
        ENDIF
        CALL P_IAPWS( TX,PVWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        PVAX = MAX( PGX-PVWX,0.D+0 )
        CALL AIRGSD( TX,PVAX,RHOGAX )
        XMGAX = PVAX/PGX
        XMGWX = MAX( 1.D+0-XMGAX,0.D+0 )
        WTMGX = XMGWX*WTMW + XMGAX*WTMA
        XGWX = XMGWX*WTMW/WTMGX
        XGAX = XMGAX*WTMA/WTMGX
        RHOGX = XGAX*RHOGAX + XGWX*RHOGWX
        PGX = PGX - RHOGX*GB
        NX = N
  100 CONTINUE
!
!---  Loop over nodes in the y direction
!
      I = IEX
      J = JSX
      K = KSX
      NX = ND(I,J,K)
      DO 110 J = JSX,JEX,JI
        N = ND(I,J,K)
        NPY = NSY(N)
        GB = (YP(N)*RP(I)-YP(NX)*RP(I))*GRVY(NPY)
        IF( TKX.LT.TCRW ) THEN
          INDX = 0
          CALL REGION_4( TX,PSWX,INDX )
          PVWX = PSWX*BCX(5)
        ELSE
          PVWX = PGX*BCX(5)
        ENDIF
        CALL P_IAPWS( TX,PVWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        PVAX = MAX( PGX-PVWX,0.D+0 )
        CALL AIRGSD( TX,PVAX,RHOGAX )
        XMGAX = PVAX/PGX
        XMGWX = MAX( 1.D+0-XMGAX,0.D+0 )
        WTMGX = XMGWX*WTMW + XMGAX*WTMA
        XGWX = XMGWX*WTMW/WTMGX
        XGAX = XMGAX*WTMA/WTMGX
        RHOGX = XGAX*RHOGAX + XGWX*RHOGWX
        PGX = PGX - RHOGX*GB
        NX = N
  110 CONTINUE
!
!---  Loop over nodes in the z direction
!
      I = IEX
      J = JEX
      K = KSX
      NX = ND(I,J,K)
      DO 120 K = KSX,KEX,KI
        N = ND(I,J,K)
        NPZ = NSZ(N)
        GB = (ZP(N)-ZP(NX))*GRVZ(NPZ)
        IF( TKX.LT.TCRW ) THEN
          INDX = 0
          CALL REGION_4( TX,PSWX,INDX )
          PVWX = PSWX*BCX(5)
        ELSE
          PVWX = PGX*BCX(5)
        ENDIF
        CALL P_IAPWS( TX,PVWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        PVAX = MAX( PGX-PVWX,0.D+0 )
        CALL AIRGSD( TX,PVAX,RHOGAX )
        XMGAX = PVAX/PGX
        XMGWX = MAX( 1.D+0-XMGAX,0.D+0 )
        WTMGX = XMGWX*WTMW + XMGAX*WTMA
        XGWX = XMGWX*WTMW/WTMGX
        XGAX = XMGAX*WTMA/WTMGX
        RHOGX = XGAX*RHOGAX + XGWX*RHOGWX
        PGX = PGX - RHOGX*GB
        NX = N
  120 CONTINUE
!
!---  Convert node-centroid pressure to boundary pressure  ---
!
      NE = IBCN(NBE)
      IF( IBCD(NBE).EQ.-3 ) THEN
        NPZ = NSZ(NE)
        PGX = PGX + 5.D-1*DZGF(NE)*RHOGX*GRVZ(NPZ)
      ELSEIF( IBCD(NBE).EQ.-2 ) THEN
        NPY = NSY(NE)
        PGX = PGX + 5.D-1*DYGF(NE)*RHOGX*GRVY(NPY)*RP(ID(NE))
      ELSEIF( IBCD(NBE).EQ.-1 ) THEN
        NPX = NSX(NE)
        PGX = PGX + 5.D-1*DXGF(NE)*RHOGX*GRVX(NPX)
      ELSEIF( IBCD(NBE).EQ.1 ) THEN
        NQX = NSX(NE)+1
        IF( INBS(4,NE).GT.0 ) NQX = INBS(4,NE)
        PGX = PGX - 5.D-1*DXGF(NE)*RHOGX*GRVX(NPX)
      ELSEIF( IBCD(NBE).EQ.2 ) THEN
        NQY = NSY(NE)+IFLD
        IF( INBS(5,NE).GT.0 ) NQY = INBS(5,NE)
        PGX = PGX - 5.D-1*DYGF(NE)*RHOGX*GRVY(NPY)*RP(ID(NE))
      ELSEIF( IBCD(NBE).EQ.3 ) THEN
        NQZ = NSZ(NE)+IJFLD
        IF( INBS(6,NE).GT.0 ) NQZ = INBS(6,NE)
        PGX = PGX - 5.D-1*DZGF(NE)*RHOGX*GRVZ(NPZ)
      ENDIF
      HGBG_GT = PGX - PATM
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HGBG_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE P_IAPWS( TX,PX,RHOGX,RHOLX,HGX,HLX,UGX,ULX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density, enthalpy, and internal energy
!     the Revised Release on the IAPWS Industrial 
!     Formulation 1997 for the Thermodynamic Properties of Water and 
!     Steam, as a function of temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/P_IAPWS'
      TKX = TX + TABS
      INDX = 1
      CALL REGION_4( TSWX,PX,INDX )
      CALL REGION_23( T23X,PX,INDX )
!
!---  Region 1  ---
!
      IF( TX.LE.TSWX .AND. TKX.LE.623.15D+0 ) THEN
        INDX = 0
        CALL REGION_4( TX,PSWX,INDX )
        PMX = MAX( PX,PSWX)
        CALL REGION_1( TX,PMX,RHOLX,HLX,ULX )
        CALL REGION_2( TX,PSWX,RHOGX,HGX,UGX )
!
!---  Region 2 (sub-623.15 K)  ---
!
      ELSEIF( TX.GT.TSWX .AND. TKX.LE.623.15D+0 ) THEN
        CALL REGION_2( TX,PX,RHOGX,HGX,UGX )
        INDX = 0
        CALL REGION_4( TX,PSWX,INDX )
        CALL REGION_1( TX,PSWX,RHOLX,HLX,ULX )
!
!---  Region 2 (sup-623.15 K and super-23 and sub-critical)  ---
!
      ELSEIF( TKX.GT.623.15D+0 .AND. TX.GT.T23X .AND. TKX.LE.TCRW ) THEN     
        CALL REGION_2( TX,PX,RHOGX,HGX,UGX )
        INDX = 0
        CALL REGION_4( TX,PSWX,INDX )
        CALL REGION_3( TX,PSWX,RHOLX,HLX,ULX )
!
!---  Region 2 (sup-623.15 K and supper-23 and super-critical)  ---
!
      ELSEIF( TKX.GT.623.15D+0 .AND. TX.GT.T23X ) THEN     
        CALL REGION_2( TX,PX,RHOGX,HGX,UGX )
        RHOLX = RHOGX
        HLX = HGX
        ULX = UGX
!
!---  Region 3 (super-critical)  ---
!
      ELSEIF( TKX.GT.TCRW .AND. TX.LE.T23X .AND. PX.GT.PCRW ) THEN
        CALL REGION_3( TX,PX,RHOGX,HGX,UGX )
        RHOLX = RHOGX
        HLX = HGX
        ULX = UGX
!
!---  Region 3 (sub-critical gas)  ---
!
      ELSEIF( TX.GT.TSWX .AND. TX.LE.T23X .AND. PX.LE.PCRW
     &  .AND. TKX.LE.TCRW ) THEN
        INDX = 0
        CALL REGION_4( TX,PSWX,INDX )
        PMX = MAX( PX,PSWX)
        CALL REGION_3SUBCRG( TX,PMX,RHOGX,HGX,UGX )
        CALL REGION_3SUBCRL( TX,PSWX,RHOLX,HLX,ULX )
!
!---  Region 3 (sub-critical liquid )  ---
!
      ELSEIF( TKX.GT.623.15D+0 .AND. TKX.LE.TCRW .AND. PX.LE.PCRW ) THEN
        INDX = 0
        CALL REGION_4( TX,PSWX,INDX )
        PMX = MAX( PX,PSWX)
        CALL REGION_3SUBCRL( TX,PMX,RHOLX,HLX,ULX )
        CALL REGION_3SUBCRG( TX,PSWX,RHOGX,HGX,UGX )
!
!---  Region 3 (sub-critical pressure gas )  ---
!
      ELSEIF( TX.GT.TSWX .AND. TX.LE.T23X .AND. PX.LE.PCRW ) THEN
        INDX = 0
        CALL REGION_4( TX,PSWX,INDX )
        PMX = MAX( PX,PSWX)
        CALL REGION_3( TX,PMX,RHOGX,HGX,UGX )
        CALL REGION_3( TX,PSWX,RHOLX,HLX,ULX )
!
!---  Region 5  ---
!
      ELSEIF( TKX.GT.1073.15D+0 ) THEN
        CALL REGION_5( TX,PX,RHOGX,HGX,UGX )
        RHOLX = RHOGX
        HLX = HGX
        ULX = UGX
      ENDIF
!      HLX = 4.184D+3*TX
!      ULX = 4.184D+3*TX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of P_IAPWS group  ---
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
!---  Reduced porosity with salt precipitation  ---
!
      PORDX = MAX( PORDX*(1.D+0-SSX),PORDX*PERM(5,IZN),1.D-12 )
!
!---  Normalized porosity  ---
!
      PORD_NX = MAX( (1.D+0-SSX-PERM(5,IZN))/(1.D+0-PERM(5,IZN)),0.D+0 )
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
      SUBROUTINE REGION_1( TX,PX,RHOX,HX,UX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density, enthalpy, and internal energy for
!     Region 1 per the Revised Release on the IAPWS Industrial 
!     Formulation 1997 for the Thermodynamic Properties of Water and 
!     Steam (The revision only relates to the extension of Region 5 
!     to 50 MPa), as a function of temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(34)
      INTEGER IX(34),JX(34)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA CNX / 0.14632971213167D+0,-0.84548187169114D+0,
     &  -0.37563603672040D+1,0.33855169168385D+1,-0.95791963387872D+0,
     &   0.15772038513228D+0,-0.16616417199501D-1,0.81214629983568D-3,
     &   0.28319080123804D-3,-0.60706301565874D-3,-0.18990068218419D-1,
     &  -0.32529748770505D-1,-0.21841717175414D-1,-0.52838357969930D-4,
     &  -0.47184321073267D-3,-0.30001780793026D-3,0.47661393906987D-4,
     &  -0.44141845330846D-5,-0.72694996297594D-15,-0.31679644845054D-4,
     &  -0.28270797985312D-5,-0.85205128120103D-9,-0.22425281908000D-5,
     &  -0.65171222895601D-6,-0.14341729937924D-12,-0.40516996860117D-6,
     &  -0.12734301741641D-8,-0.17424871230634D-9,-0.68762131295531D-18,
     &   0.14478307828521D-19,0.26335781662795D-22,
     &  -0.11947622640071D-22,0.18228094581404D-23,
     &  -0.93537087292458D-25 /
      DATA IX / 0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,4,5,
     &   8,8,21,23,29,30,31,32 /
      DATA JX / -2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1,3,17,-4,0,6,
     &  -5,-2,10,-8,-11,-6,-29,-31,-38,-39,-40,-41 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_1'
!
!---  Reduced pressure, "pi" and reduced temperature, "tau"  ---
!
      RCX = 0.461526D+3
      PRX = 16.53D+6
      TRX = 1386.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      TAUX = TRX/TKX
!
!---  Partial derivative of the dimensionless Gibbs free energy with
!     respect to the reduced pressure, "pi" and temperature, "tau"  ---
!
      GAMMA_PIX = 0.D+0
      GAMMA_TAUX = 0.D+0
      DO 100 I = 1,34
        GAMMA_PIX = GAMMA_PIX - CNX(I)*REAL(IX(I))*
     &    ((7.1D+0-PIX)**(IX(I)-1))*((TAUX-1.222D+0)**JX(I))
        GAMMA_TAUX = GAMMA_TAUX + CNX(I)*((7.1D+0-PIX)**IX(I))*
     &    REAL(JX(I))*((TAUX-1.222D+0)**(JX(I)-1))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(PIX*GAMMA_PIX*RCX*TKX/PX)
!
!---  Enthalpy, J/kg  ---
!
      HX = RCX*TKX*TAUX*GAMMA_TAUX
!
!---  Internal energy, J/kg  ---
!
      UX = RCX*TKX*(TAUX*GAMMA_TAUX-PIX*GAMMA_PIX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_1 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_2( TX,PX,RHOX,HX,UX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density, enthalpy, and internal energy for
!     Region 2 per the Revised Release on the IAPWS Industrial 
!     Formulation 1997 for the Thermodynamic Properties of Water and 
!     Steam (The revision only relates to the extension of Region 5 
!     to 50 MPa), as a function of temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNOX(9),CNRX(43)
      INTEGER JOX(9),IRX(43),JRX(43)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNOX,CNRX
      SAVE JOX,IRX,JRX
      DATA CNOX / -0.96927686500217D+1,0.10086655968018D+2,
     &  -0.56087911283020D-2,0.71452738081455D-1,-0.40710498223928D+0,
     &   0.14240819171444D+1,-0.43839511319450D+1,-0.28408632460772D+0,
     &   0.21268463753307D-1 /
      DATA CNRX / -0.17731742473213D-2,-0.17834862292358D-1,
     &  -0.45996013696365D-1,-0.57581259083432D-1,-0.50325278727930D-1,
     &  -0.33032641670203D-4,-0.18948987516315D-3,-0.39392777243355D-2,
     &  -0.43797295650573D-1,-0.26674547914087D-4,0.20481737692309D-7,
     &   0.43870667284435D-6,-0.32277677238570D-4,-0.15033924542148D-2,
     &  -0.40668253562649D-1,-0.78847309559367D-9,0.12790717852285D-7,
     &   0.48225372718507D-6,0.22922076337661D-5,-0.16714766451061D-10,
     &  -0.21171472321355D-2,-0.23895741934104D+2,-0.59059564324270D-17,
     &  -0.12621808899101D-5,-0.38946842435739D-1,0.11256211360459D-10,
     &  -0.82311340897998D+1,0.19809712802088D-7,0.10406965210174D-18,
     &  -0.10234747095929D-12,-0.10018179379511D-8,
     &  -0.80882908646985D-10,0.10693031879409D+0,-0.33662250574171D+0,
     &   0.89185845355421D-24,0.30629316876232D-12,
     &  -0.42002467698208D-5,-0.59056029685639D-25,
     &   0.37826947613457D-5,-0.12768608934681D-14,0.73087610595061D-28, 
     &   0.55414715350778D-16,-0.94369707241210D-6 /
      DATA JOX / 0,1,-5,-4,-3,-2,-1,2,3 /
      DATA IRX / 1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,5,6,6,6,7,7,7,
     &  8,8,9,10,10,10,16,16,18,20,20,20,21,22,23,24,24,24 /
      DATA JRX / 0,1,2,3,6,1,2,4,7,36,0,1,3,6,35,1,2,3,7,3,16,35,0,
     &  11,25,8,36,13,4,10,14,29,50,57,20,35,48,21,53,39,26,40,58 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_2'
!
!---  Reduced pressure, "pi" and reduced temperature, "tau"  ---
!
      RCX = 0.461526D+3
      PRX = 1.D+6
      TRX = 540.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      TAUX = TRX/TKX
!
!---  Partial derivative of the ideal-gas part of the dimensionless 
!     Gibbs free energy with respect to the reduced pressure, "pi" 
!     and temperature, "tau"  ---
!
      GAMMAO_PIX = 1.D+0/PIX
      GAMMAO_TAUX = 0.D+0
      DO 100 I = 1,9
        GAMMAO_TAUX = GAMMAO_TAUX + CNOX(I)*REAL(JOX(I))*
     &    (TAUX**(JOX(I)-1))
  100 CONTINUE
!
!---  Partial derivative of the residual part of the dimensionless 
!     Gibbs free energy with respect to the reduced pressure, "pi" 
!     and temperature, "tau"  ---
!
      GAMMAR_PIX = 0.D+0
      GAMMAR_TAUX = 0.D+0
      DO 200 I = 1,43
        GAMMAR_PIX = GAMMAR_PIX + CNRX(I)*REAL(IRX(I))*
     &    (PIX**(IRX(I)-1))*((TAUX-5.D-1)**JRX(I))
        GAMMAR_TAUX = GAMMAR_TAUX + CNRX(I)*REAL(JRX(I))*
     &    (PIX**IRX(I))*((TAUX-5.D-1)**(JRX(I)-1))
  200 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(PIX*(GAMMAO_PIX + GAMMAR_PIX)*RCX*TKX/PX)
!
!---  Enthalpy, J/kg  ---
!
      HX = RCX*TKX*TAUX*(GAMMAO_TAUX + GAMMAR_TAUX)
!
!---  Internal energy, J/kg  ---
!
      UX = RCX*TKX*(TAUX*(GAMMAO_TAUX + GAMMAR_TAUX) - 
     &  PIX*(GAMMAO_PIX + GAMMAR_PIX))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_23( TX,PX,INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Auxiliary equation for the boundary between Regions 2 and 3.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(5)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX
      DATA CNX / 0.34805185628969D+3, -0.11671859879975D+1,
     &  0.10192970039326D-2, 0.57254459862746D+3, 0.13918839778870D+2 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_23'
      PRX = 1.D+6
!
!---  INDX = 0, return pressure  ---
!
      IF( INDX.EQ.0 ) THEN
        TKX = TX + TABS
        IF( TKX.LT.623.15D+0 ) THEN
          PX = 16.5292D+6
        ELSE
          PX = PRX*(CNX(1) + CNX(2)*TKX + CNX(3)*(TKX**2))
        ENDIF
!
!---  INDX != 0, return temperature  ---
!
      ELSE
        IF( PX.LT.16.5292D+6 ) THEN
          TX = 623.15D+0 - TABS
        ELSE
          PIX = PX/PRX
          TKX = CNX(4) + SQRT((PIX-CNX(5))/CNX(3))
          TX = TKX - TABS
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_23 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3BIS( TX,PX,RHOX,HX,UX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density, enthalpy, and internal energy for
!     Region 3 per the Revised Release on the IAPWS Industrial 
!     Formulation 1997 for the Thermodynamic Properties of Water and 
!     Steam (The revision only relates to the extension of Region 5 
!     to 50 MPa), as a function of temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(40),CSLX(6),CSVX(6)
      INTEGER IX(40),JX(40)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,CSLX,CSVX
      SAVE IX,JX
      DATA CNX / 0.10658070028513D+1,-0.15732845290239D+2,
     &  0.20944396974307D+2,-0.76867707878716D+1,0.26185947787954D+1,
     & -0.28080781148620D+1,0.12053369696517D+1,-0.84566812812502D-2,
     & -0.12654315477714D+1,-0.11524407806681D+1,0.88521043984318D+0,
     & -0.64207765181607D+0,0.38493460186671D+0,-0.85214708824206D+0,
     &  0.48972281541877D+1,-0.30502617256965D+1,0.39420536879154D-1,
     &  0.12558408424308D+0,-0.27999329698710D+0,0.13899799569460D+1,
     & -0.20189915023570D+1,-0.82147637173963D-2,-0.47596035734923D+0,
     &  0.43984074473500D-1,-0.44476435428739D+0,0.90572070719733D+0,
     &  0.70522450087967D+0,0.10770512626332D+0,-0.32913623258954D+0,
     & -0.50871062041158D+0,-0.22175400873096D-1,0.94260751665092D-1,
     &  0.16436278447961D+0,-0.13503372241348D-1,-0.14834345352472D-1,
     &  0.57922953628084D-3,0.32308904703711D-2,0.80964802996215D-4,
     & -0.16557679795037D-3,-0.44923899061815D-4 /
      DATA CSLX / 1.99274064D+0,1.09965342D+0,-0.510839303D+0,
     & -1.75493479D+0,-45.5170352D+0,-6.74694450D+5 /
      DATA CSVX / -2.03150240D+0,-2.68302940D+0,-5.38626492D+0,
     & -17.2991605D+0,-44.7586581D+0,-63.9201063D+0 /
      DATA IX / 0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,
     &  5,5,5,6,6,6,7,8,9,9,10,10,11 /
      DATA JX / 0,0,1,2,7,10,12,23,2,6,15,17,0,2,6,7,22,26,0,2,4,16,26,
     &  0,2,4,26,1,3,26,0,2,26,2,26,2,26,0,1,26 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3BIS'
!
!---  Critical temperature and density and specific gas constant  ---
!
      TCRX = 647.096D+0
      PCRX = 22.064D+6
      RHOCRX = 322.D+0
      RCX = 0.461526D+3
!
!---  Reduced temperature, "tau"  ---
!
      TKX = TX + TABS
      TAUX = TCRX/TKX
!
!---  Initial guess of density, using the density computed from Region 1
!     at pressure and temperature = 623.15 K, and density computed
!     from Region 3 at pressure and temperature = the boundary between
!     Regions 2 and 3  ---
!
      IF( PX.GT.PCRX ) THEN
        T1X = 623.15 - TABS
        CALL REGION_1( T1X,PX,RHO1X,HX,UX )
        INDX = 1
        CALL REGION_23( T2X,PX,INDX )
        CALL REGION_2( T2X,PX,RHO2X,HX,UX )
        RHOX = (RHO2X-RHO1X)*(TX-T1X)/(T2X-T1X) + RHO1X
        DELTA1X = RHO1X/RHOCRX
        DELTA2X = RHO2X/RHOCRX
      ELSE
        INDX = 1
        CALL REGION_4( T4X,PX,INDX )
!
!---    Initial guess of density, using the density computed from 
!       Region 1 at pressure and temperature = 623.15 K, and density 
!       computed from density of saturated liquid at temperature  ---
!
        IF( TX.LE.T4X ) THEN
          T1X = 623.15 - TABS
          CALL REGION_1( T1X,PX,RHO1X,HX,UX )
          DELTA1X = RHO1X/RHOCRX
          TAU2X = 1.D+0 - TKX/TCRX
          DELTA2X = 1.D+0 + CSLX(1)*(TAU2X**(1.D+0/3.D+0)) +
     &      CSLX(2)*(TAU2X**(2.D+0/3.D+0)) + 
     &      CSLX(3)*(TAU2X**(5.D+0/3.D+0)) + 
     &      CSLX(4)*(TAU2X**(16.D+0/3.D+0)) + 
     &      CSLX(5)*(TAU2X**(43.D+0/3.D+0)) + 
     &      CSLX(6)*(TAU2X**(110.D+0/3.D+0))
          RHO2X = DELTA2X*RHOCRX
          T2X = T4X
          RHOX = (RHO2X-RHO1X)*(TX-T1X)/(T2X-T1X) + RHO1X
!
!---    Initial guess of density, using the density computed from
!       saturated vapor at temperature, and density computed
!       from Region 3 at pressure and temperature = the boundary between
!       Regions 2 and 3  ---
!
        ELSE
          TAU1X = 1.D+0 - TKX/TCRX
          DELTA1X = EXP(CSVX(1)*(TAU1X**(2.D+0/6.D+0)) +
     &      CSVX(2)*(TAU1X**(4.D+0/6.D+0)) + 
     &      CSVX(3)*(TAU1X**(8.D+0/6.D+0)) + 
     &      CSVX(4)*(TAU1X**(18.D+0/6.D+0)) + 
     &      CSVX(5)*(TAU1X**(37.D+0/6.D+0)) + 
     &      CSVX(6)*(TAU1X**(71.D+0/6.D+0)))
          RHO1X = DELTA1X*RHOCRX
          T1X = T4X
          INDX = 1
          CALL REGION_23( T2X,PX,INDX )
          CALL REGION_2( T2X,PX,RHO2X,HX,UX )
          RHOX = (RHO2X-RHO1X)*(TX-T1X)/(T2X-T1X) + RHO1X
          DELTA2X = RHO2X/RHOCRX
        ENDIF
      ENDIF
!
!---  Reduced density, "delta"  ---
!
      DELTAX = RHOX/RHOCRX
!
!---  Newton-Raphson loop to determine reduced density, "delta"  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 12
        IMSGX = N_DB
        NMSGX = 0
        SUBLOGX = 'REGION_3BIS'
        RLMSGX = 0.D+0
        CHMSG = 'Region 3BIS Convergence Failure: ' // 
     &    'Newton-Raphson @ Node'
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  First and second partial derivatives of the Helmholtz free energy
!     with respect to the reduced density, "delta"  ---
!
      PHI_DELTAX = CNX(1)/DELTAX
      PHI_2DELTAX = -CNX(1)/(DELTAX**2)
      DO 200 I = 2,40
        PHI_DELTAX = PHI_DELTAX + CNX(I)*REAL(IX(I))*
     &    (DELTAX**(IX(I)-1))*(TAUX**JX(I))
        PHI_2DELTAX = PHI_2DELTAX + CNX(I)*REAL(IX(I))*(REAL(IX(I)-1))*
     &    (DELTAX**(IX(I)-2))*(TAUX**JX(I))
  200 CONTINUE
      FX = PX - (DELTAX**2)*PHI_DELTAX*RCX*TKX*RHOCRX
      DFX = -RCX*TKX*RHOCRX*(2.D+0*DELTAX*PHI_DELTAX + 
     &  (DELTAX**2)*PHI_2DELTAX)
      DX = -FX/DFX
      DX = SIGN( MIN(0.125D+0,ABS(DX)),DX )
      DELTAX = DELTAX + DX
      DELTAX = MAX( MIN( DELTAX,DELTA1X),DELTA2X )
      IF( ABS(DX).GT.1.D-6 ) GOTO 100
!
!---  Density  ---
!
      RHOX = DELTAX*RHOCRX
!
!---  First partial derivatives of the Helmholtz free energy
!     with respect to the reduced density, "delta" and reduced
!     temperature, "tau" ---
!
      PHI_DELTAX = CNX(1)/DELTAX
      PHI_TAUX = 0.D+0
      DO 300 I = 2,40
        PHI_DELTAX = PHI_DELTAX + CNX(I)*REAL(IX(I))*
     &    (DELTAX**(IX(I)-1))*(TAUX**JX(I))
        PHI_TAUX = PHI_TAUX + CNX(I)*(DELTAX**IX(I))*REAL(JX(I))*
     &    (TAUX**(JX(I)-1))
  300 CONTINUE
!
!---  Enthalpy, J/kg  ---
!
      HX = RCX*TKX*(TAUX*PHI_TAUX + DELTAX*PHI_DELTAX)
!
!---  Internal energy, J/kg  ---
!
      UX = RCX*TKX*TAUX*PHI_TAUX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3BIS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3( TX,PX,RHOX,HX,UX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density, enthalpy, and internal energy for
!     Region 3 per the Revised Release on the IAPWS Industrial 
!     Formulation 1997 for the Thermodynamic Properties of Water and 
!     Steam (The revision only relates to the extension of Region 5 
!     to 50 MPa), as a function of temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(40)
      REAL*8 C3ABX(5),C3CDX(4),C3GHX(5),C3IJX(5),C3JKX(5),C3MNX(4)
      REAL*8 C3OPX(5),C3QUX(4),C3RXX(4),C3UVX(4),C3WXX(5)
      REAL*8 P3X(14)
      INTEGER IX(40),JX(40)
      INTEGER I3ABX(5),I3CDX(4),I3GHX(5),I3IJX(5),I3JKX(5),I3MNX(4)
      INTEGER I3OPX(5),I3QUX(4),I3RXX(4),I3UVX(4),I3WXX(5)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,C3ABX,C3CDX,C3GHX,C3IJX,C3JKX,C3MNX,C3OPX,C3QUX,C3RXX
      SAVE C3UVX,C3WXX
      SAVE P3X
      SAVE IX,JX,I3ABX,I3CDX,I3GHX,I3JKX,I3MNX,I3OPX,I3QUX,I3RXX
      SAVE I3UVX,I3WXX
      SAVE TCRX,PCRX,RHOCRX,RCX
      DATA TCRX / 647.096D+0 /
      DATA PCRX / 22.064D+6 /
      DATA RHOCRX / 322.D+0 /
      DATA RCX / 0.461526D+3 /
      DATA CNX / 0.10658070028513D+1,-0.15732845290239D+2,
     &  0.20944396974307D+2,-0.76867707878716D+1,0.26185947787954D+1,
     & -0.28080781148620D+1,0.12053369696517D+1,-0.84566812812502D-2,
     & -0.12654315477714D+1,-0.11524407806681D+1,0.88521043984318D+0,
     & -0.64207765181607D+0,0.38493460186671D+0,-0.85214708824206D+0,
     &  0.48972281541877D+1,-0.30502617256965D+1,0.39420536879154D-1,
     &  0.12558408424308D+0,-0.27999329698710D+0,0.13899799569460D+1,
     & -0.20189915023570D+1,-0.82147637173963D-2,-0.47596035734923D+0,
     &  0.43984074473500D-1,-0.44476435428739D+0,0.90572070719733D+0,
     &  0.70522450087967D+0,0.10770512626332D+0,-0.32913623258954D+0,
     & -0.50871062041158D+0,-0.22175400873096D-1,0.94260751665092D-1,
     &  0.16436278447961D+0,-0.13503372241348D-1,-0.14834345352472D-1,
     &  0.57922953628084D-3,0.32308904703711D-2,0.80964802996215D-4,
     & -0.16557679795037D-3,-0.44923899061815D-4 /
      DATA C3ABX / 0.154793642129415D+4,-0.187661219490113D+3,
     &  0.213144632222113D+2,-0.191887498864292D+4,
     &  0.918419702359447D+3 /
      DATA C3CDX / 0.585276966696349D+3,0.278233532206915D+1,
     &  -0.127283549295878D-1,0.159090746562729D-3 /
      DATA C3GHX / -0.249284240900418D+5,0.428143584791546D+4,
     &  -0.269029173140130D+3,0.751608051114157D+1,
     &  -0.787105249910383D-1 /
      DATA C3IJX / 0.584814781649163D+3,-0.616179320924617D+0,
     &  0.260763050899562D+0,-0.587071076864459D-2,
     &  0.515308185433082D-4 /
      DATA C3JKX / 0.617229772068439D+3,-0.770600270141675D+1,
     &  0.697072596851896D+0,-0.157391839848015D-1,
     &  0.137897492684194D-3 /
      DATA C3MNX / 0.535339483742384D+3,0.761978122720128D+1,
     &  -0.158365725441648,0.192871054508108D-2 /
      DATA C3OPX / 0.969461372400213D+3,-0.332500170441278D+3,
     &  0.642859598466067D+2,0.773845935768222D+3,
     &  -0.152313732937084D+4 /
      DATA C3QUX / 0.565603648239126D+3,0.529062258221222D+1,
     &  -0.102020639611016D+0,0.122240301070145D-2 /
      DATA C3RXX / 0.584561202520006D+3,-0.102961025163669D+1,
     &  0.243293362700452D+0,-0.294905044740799D-2 /
      DATA C3UVX / 0.528199646263062D+3,0.890579602135307D+1,
     &  -0.222814134903755,0.286791682263697D-2 /
      DATA C3WXX / 0.728052609145380D+1,0.973505869861952D+2,
     &  0.147370491183191D+2,0.329196213998375D+3,
     &  0.873371668682417D+3 /
      DATA IX / 0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,
     &  5,5,5,6,6,6,7,8,9,9,10,10,11 /
      DATA JX / 0,0,1,2,7,10,12,23,2,6,15,17,0,2,6,7,22,26,0,2,4,16,26,
     &  0,2,4,26,1,3,26,0,2,26,2,26,2,26,0,1,26 /
      DATA I3ABX / 0,1,2,-1,-2 /
      DATA I3CDX / 0,1,2,3 /
      DATA I3GHX / 0,1,2,3,4 /
      DATA I3IJX / 0,1,2,3,4 /
      DATA I3JKX / 0,1,2,3,4 /
      DATA I3MNX / 0,1,2,3 /
      DATA I3OPX / 0,1,2,-1,-2 /
      DATA I3QUX / 0,1,2,3 /
      DATA I3RXX / 0,1,2,3 /
      DATA I3UVX / 0,1,2,3 /
      DATA I3WXX / 0,1,2,-1,-2 /
      DATA P3X / 16.5292D+6,19.00881189173929D+6,20.5D+6,
     &  21.04336732D+6,21.90096265D+6,21.93161551D+6,22.064D+6,
     &  22.11D+6,22.5D+6,23.D+6,23.5D+6,25.D+6,40.D+6,100.D+6 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3'
      TKX = TX + TABS
      PRX = 1.D+6
      PIX = PX/PRX
      IF( TKX.LE.TCRX ) THEN
        INDX = 1
        CALL REGION_4( TSAT97X,PX,INDX )
        TSAT97X = TSAT97X + TABS
      ELSE
        TSAT97X = TCRX
      ENDIF
      INDX = 1
      CALL REGION_23( TB23X,PX,INDX )
      TB23X = TB23X + TABS
!
!---    ab subregion boundary  ---
!
        T3ABX = 0.D+0
        DO 100 I = 1,5
          T3ABX = T3ABX + C3ABX(I)*(LOG(PIX)**I3ABX(I))
  100   CONTINUE
!
!---    cd subregion boundary  ---
!
        T3CDX = 0.D+0
        DO 110 I = 1,4
          T3CDX = T3CDX + C3CDX(I)*(PIX**I3CDX(I))
  110   CONTINUE
!
!---    ef subregion boundary  ---
!
        T3EFX = 3.727888004D+0*(PIX-22.064) + 647.096
!
!---    gh subregion boundary  ---
!
        T3GHX = 0.D+0
        DO 120 I = 1,5
          T3GHX = T3GHX + C3GHX(I)*(PIX**I3GHX(I))
  120   CONTINUE
!
!---    ij subregion boundary  ---
!
        T3IJX = 0.D+0
        DO 130 I = 1,5
          T3IJX = T3IJX + C3IJX(I)*(PIX**I3IJX(I))
  130   CONTINUE
!
!---    jk subregion boundary  ---
!
        T3JKX = 0.D+0
        DO 140 I = 1,5
          T3JKX = T3JKX + C3JKX(I)*(PIX**I3JKX(I))
  140   CONTINUE
!
!---    mn subregion boundary  ---
!
        T3MNX = 0.D+0
        DO 150 I = 1,4
          T3MNX = T3MNX + C3MNX(I)*(PIX**I3MNX(I))
  150   CONTINUE
!
!---    op subregion boundary  ---
!
        T3OPX = 0.D+0
        DO 160 I = 1,5
          T3OPX = T3OPX + C3OPX(I)*(LOG(PIX)**I3OPX(I))
  160   CONTINUE
!
!---    qu subregion boundary  ---
!
        T3QUX = 0.D+0
        DO 170 I = 1,4
          T3QUX = T3QUX + C3QUX(I)*(PIX**I3QUX(I))
  170   CONTINUE
!
!---    rx subregion boundary  ---
!
        T3RXX = 0.D+0
        DO 180 I = 1,4
          T3RXX = T3RXX + C3RXX(I)*(PIX**I3RXX(I))
  180   CONTINUE
!
!---    uv subregion boundary  ---
!
        T3UVX = 0.D+0
        DO 190 I = 1,4
          T3UVX = T3UVX + C3UVX(I)*(PIX**I3UVX(I))
  190   CONTINUE
!
!---    wx subregion boundary  ---
!
        T3WXX = 0.D+0
        DO 200 I = 1,5
          T3WXX = T3WXX + C3WXX(I)*(LOG(PIX)**I3WXX(I))
  200   CONTINUE
!
!---    Region 3a  ---
!
        IF( PX.GT.P3X(13) .AND. PX.LE.P3X(14) .AND. 
     &    TKX.GT.623.15D+0 .AND. TKX.LE.T3ABX ) THEN
          CALL REGION_3A( TX,PX,RHOX )
!
!---    Region 3b  ---
!
        ELSEIF( PX.GT.P3X(13) .AND. PX.LE.P3X(14) .AND. 
     &    TKX.GT.T3ABX ) THEN
          CALL REGION_3B( TX,PX,RHOX )
!
!---    Region 3c  ---
!
        ELSEIF( ( PX.GT.P3X(2) .AND. PX.LE.P3X(13) .AND. 
     &    TKX.LE.T3CDX ) .OR. ( PX.GT.P3X(1) .AND. PX.LE.P3X(2) .AND.
     &    TKX.LE.TSAT97X ) ) THEN
          CALL REGION_3C( TX,PX,RHOX )
!
!---    Region 3d  ---
!
        ELSEIF( PX.GT.P3X(12) .AND. PX.LE.P3X(13) .AND. 
     &    TKX.GT.T3CDX .AND. TKX.LE.T3ABX ) THEN
          CALL REGION_3D( TX,PX,RHOX )
!
!---    Region 3e  ---
!
        ELSEIF( PX.GT.P3X(12) .AND. PX.LE.P3X(13) .AND. 
     &    TKX.GT.T3ABX .AND. TKX.LE.T3EFX ) THEN
          CALL REGION_3E( TX,PX,RHOX )
!
!---    Region 3f  ---
!
        ELSEIF( PX.GT.P3X(12) .AND. PX.LE.P3X(13) .AND. 
     &    TKX.GT.T3EFX .AND. TKX.LE.TB23X ) THEN
          CALL REGION_3F( TX,PX,RHOX )
!
!---    Region 3g  ---
!
        ELSEIF( PX.GT.P3X(11) .AND. PX.LE.P3X(12) .AND. 
     &    TKX.GT.T3CDX .AND. TKX.LE.T3GHX ) THEN
          CALL REGION_3G( TX,PX,RHOX )
!
!---    Region 3h  ---
!
        ELSEIF( PX.GT.P3X(10) .AND. PX.LE.P3X(12) .AND. 
     &    TKX.GT.T3GHX .AND. TKX.LE.T3EFX ) THEN
          CALL REGION_3H( TX,PX,RHOX )
!
!---    Region 3i  ---
!
        ELSEIF( PX.GT.P3X(10) .AND. PX.LE.P3X(12) .AND. 
     &    TKX.GT.T3EFX .AND. TKX.LE.T3IJX ) THEN
          CALL REGION_3I( TX,PX,RHOX )
!
!---    Region 3j  ---
!
        ELSEIF( PX.GT.P3X(9) .AND. PX.LE.P3X(12) .AND. 
     &    TKX.GT.T3IJX .AND. TKX.LE.T3JKX ) THEN
          CALL REGION_3J( TX,PX,RHOX )
!
!---    Region 3k  ---
!
        ELSEIF( PX.GT.P3X(3) .AND. PX.LE.P3X(12) .AND. 
     &    TKX.GT.T3JKX .AND. TKX.LE.TB23X ) THEN
          CALL REGION_3K( TX,PX,RHOX )
!
!---    Region 3l  ---
!
        ELSEIF( PX.GT.P3X(9) .AND. PX.LE.P3X(11) .AND. 
     &    TKX.GT.T3CDX .AND. TKX.LE.T3GHX ) THEN
          CALL REGION_3L( TX,PX,RHOX )
!
!---    Region 3m  ---
!
        ELSEIF( PX.GT.P3X(9) .AND. PX.LE.P3X(10) .AND. 
     &    TKX.GT.T3GHX .AND. TKX.LE.T3MNX ) THEN
          CALL REGION_3M( TX,PX,RHOX )
!
!---    Region 3n  ---
!
        ELSEIF( PX.GT.P3X(9) .AND. PX.LE.P3X(10) .AND. 
     &    TKX.GT.T3MNX .AND. TKX.LE.T3EFX ) THEN
          CALL REGION_3N( TX,PX,RHOX )
!
!---    Region 3o  ---
!
        ELSEIF( PX.GT.P3X(9) .AND. PX.LE.P3X(10) .AND. 
     &    TKX.GT.T3EFX .AND. TKX.LE.T3OPX ) THEN
          CALL REGION_3O( TX,PX,RHOX )
!
!---    Region 3p  ---
!
        ELSEIF( PX.GT.P3X(9) .AND. PX.LE.P3X(10) .AND. 
     &    TKX.GT.T3OPX .AND. TKX.LE.T3IJX ) THEN
          CALL REGION_3P( TX,PX,RHOX )
!
!---    Region 3q  ---
!
        ELSEIF( PX.GT.P3X(4) .AND. PX.LE.P3X(9) .AND. 
     &    TKX.GT.T3CDX .AND. TKX.LE.T3QUX ) THEN
          CALL REGION_3Q( TX,PX,RHOX )
!
!---    Region 3r  ---
!
        ELSEIF( PX.GT.P3X(3) .AND. PX.LE.P3X(9) .AND. 
     &    TKX.GT.T3RXX .AND. TKX.LE.T3JKX ) THEN
          CALL REGION_3R( TX,PX,RHOX )
!
!---    Region 3s  ---
!
        ELSEIF( PX.GT.P3X(2) .AND. PX.LE.P3X(4) .AND. 
     &    TKX.GT.T3CDX .AND. TKX.LE.TSAT97X ) THEN
          CALL REGION_3S( TX,PX,RHOX )
!
!---    Region 3t  ---
!
        ELSEIF( PX.GT.P3X(1) .AND. PX.LE.P3X(3) .AND. 
     &    TKX.GT.TSAT97X .AND. TKX.LE.TB23X ) THEN
          CALL REGION_3T( TX,PX,RHOX )
!
!---    Region 3u  ---
!
        ELSEIF( ( PX.GT.P3X(4) .AND. PX.LE.P3X(6) .AND. 
     &    TKX.GT.T3QUX .AND. TKX.LE.TSAT97X ) .OR. 
     &    ( PX.GT.P3X(6) .AND. PX.LE.P3X(9) .AND. 
     &    TKX.GT.T3QUX .AND. TKX.LE.T3UVX ) ) THEN
          CALL REGION_3U( TX,PX,RHOX )
!
!---    Region 3v  ---
!
        ELSEIF( PX.GT.P3X(8) .AND. PX.LE.P3X(9) .AND. 
     &    TKX.GT.T3UVX .AND. TKX.LE.T3EFX ) THEN
          CALL REGION_3V( TX,PX,RHOX )
!
!---    Region 3w  ---
!
        ELSEIF( PX.GT.P3X(8) .AND. PX.LE.P3X(9) .AND. 
     &    TKX.GT.T3EFX .AND. TKX.LE.T3WXX ) THEN
          CALL REGION_3W( TX,PX,RHOX )
!
!---    Region 3x  ---
!
        ELSEIF( ( PX.GT.P3X(4) .AND. PX.LE.P3X(5) .AND. 
     &    TKX.GT.TSAT97X .AND. TKX.LE.T3RXX ) .OR. 
     &    ( PX.GT.P3X(5) .AND. PX.LE.P3X(9) .AND. 
     &    TKX.GT.T3WXX .AND. TKX.LE.T3RXX ) ) THEN
          CALL REGION_3X( TX,PX,RHOX )
!
!---    Region 3y  ---
!
        ELSEIF( ( PX.GT.P3X(6) .AND. PX.LE.P3X(7) .AND. 
     &    TKX.GT.T3UVX .AND. TKX.LE.TSAT97X ) .OR. 
     &    ( PX.GT.P3X(7) .AND. PX.LE.P3X(8) .AND. 
     &    TKX.GT.T3UVX .AND. TKX.LE.T3EFX ) ) THEN
          CALL REGION_3Y( TX,PX,RHOX )
!
!---    Region 3z  ---
!
        ELSEIF( ( PX.GT.P3X(5) .AND. PX.LE.P3X(7) .AND. 
     &    TKX.GT.TSAT97X .AND. TKX.LE.T3WXX ) .OR. 
     &    ( PX.GT.P3X(7) .AND. PX.LE.P3X(8) .AND. 
     &    TKX.GT.T3EFX .AND. TKX.LE.T3WXX ) ) THEN
          CALL REGION_3Z( TX,PX,RHOX )
!
!---    Region 3 subregion not identified  ---
!
        ELSE
          INDX = 12
          IMSGX = N_DB
          NMSGX = 0
          SUBLOGX = 'REGION_3'
          RLMSGX = 0.D+0
          CHMSGX = 'Region 3 Subregion Not Identified: ' // 
     &      ' @ Node'
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
!
!---  First partial derivatives of the Helmholtz free energy
!     with respect to the reduced density, "delta" and reduced
!     temperature, "tau" ---
!
      DELTAX = RHOX/RHOCRX
      TAUX = TCRX/TKX
      PHI_DELTAX = CNX(1)/DELTAX
      PHI_TAUX = 0.D+0
      DO 300 I = 2,40
        PHI_DELTAX = PHI_DELTAX + CNX(I)*REAL(IX(I))*
     &    (DELTAX**(IX(I)-1))*(TAUX**JX(I))
        PHI_TAUX = PHI_TAUX + CNX(I)*(DELTAX**IX(I))*REAL(JX(I))*
     &    (TAUX**(JX(I)-1))
  300 CONTINUE
!
!---  Enthalpy, J/kg  ---
!
      HX = RCX*TKX*(TAUX*PHI_TAUX + DELTAX*PHI_DELTAX)
!
!---  Internal energy, J/kg  ---
!
      UX = RCX*TKX*TAUX*PHI_TAUX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3SUBCRG( TX,PX,RHOX,HX,UX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density, enthalpy, and internal energy for
!     Region 3 (subcritical gas) per the Revised Release on the IAPWS  
!     Industrial Formulation 1997 for the Thermodynamic Properties of  
!     Water and Steam (The revision only relates to the extension of  
!     Region 5 to 50 MPa), as a function of temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(40)
      REAL*8 C3ABX(5),C3CDX(4),C3GHX(5),C3IJX(5),C3JKX(5),C3MNX(4)
      REAL*8 C3OPX(5),C3QUX(4),C3RXX(4),C3UVX(4),C3WXX(5)
      REAL*8 P3X(14)
      INTEGER IX(40),JX(40)
      INTEGER I3ABX(5),I3CDX(4),I3GHX(5),I3IJX(5),I3JKX(5),I3MNX(4)
      INTEGER I3OPX(5),I3QUX(4),I3RXX(4),I3UVX(4),I3WXX(5)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,C3ABX,C3CDX,C3GHX,C3IJX,C3JKX,C3MNX,C3OPX,C3QUX,C3RXX
      SAVE C3UVX,C3WXX
      SAVE P3X
      SAVE IX,JX,I3ABX,I3CDX,I3GHX,I3JKX,I3MNX,I3OPX,I3QUX,I3RXX
      SAVE I3UVX,I3WXX
      SAVE TCRX,PCRX,RHOCRX,RCX
      DATA TCRX / 647.096D+0 /
      DATA PCRX / 22.064D+6 /
      DATA RHOCRX / 322.D+0 /
      DATA RCX / 0.461526D+3 /
      DATA CNX / 0.10658070028513D+1,-0.15732845290239D+2,
     &  0.20944396974307D+2,-0.76867707878716D+1,0.26185947787954D+1,
     & -0.28080781148620D+1,0.12053369696517D+1,-0.84566812812502D-2,
     & -0.12654315477714D+1,-0.11524407806681D+1,0.88521043984318D+0,
     & -0.64207765181607D+0,0.38493460186671D+0,-0.85214708824206D+0,
     &  0.48972281541877D+1,-0.30502617256965D+1,0.39420536879154D-1,
     &  0.12558408424308D+0,-0.27999329698710D+0,0.13899799569460D+1,
     & -0.20189915023570D+1,-0.82147637173963D-2,-0.47596035734923D+0,
     &  0.43984074473500D-1,-0.44476435428739D+0,0.90572070719733D+0,
     &  0.70522450087967D+0,0.10770512626332D+0,-0.32913623258954D+0,
     & -0.50871062041158D+0,-0.22175400873096D-1,0.94260751665092D-1,
     &  0.16436278447961D+0,-0.13503372241348D-1,-0.14834345352472D-1,
     &  0.57922953628084D-3,0.32308904703711D-2,0.80964802996215D-4,
     & -0.16557679795037D-3,-0.44923899061815D-4 /
      DATA C3ABX / 0.154793642129415D+4,-0.187661219490113D+3,
     &  0.213144632222113D+2,-0.191887498864292D+4,
     &  0.918419702359447D+3 /
      DATA C3CDX / 0.585276966696349D+3,0.278233532206915D+1,
     &  -0.127283549295878D-1,0.159090746562729D-3 /
      DATA C3GHX / -0.249284240900418D+5,0.428143584791546D+4,
     &  -0.269029173140130D+3,0.751608051114157D+1,
     &  -0.787105249910383D-1 /
      DATA C3IJX / 0.584814781649163D+3,-0.616179320924617D+0,
     &  0.260763050899562D+0,-0.587071076864459D-2,
     &  0.515308185433082D-4 /
      DATA C3JKX / 0.617229772068439D+3,-0.770600270141675D+1,
     &  0.697072596851896D+0,-0.157391839848015D-1,
     &  0.137897492684194D-3 /
      DATA C3MNX / 0.535339483742384D+3,0.761978122720128D+1,
     &  -0.158365725441648,0.192871054508108D-2 /
      DATA C3OPX / 0.969461372400213D+3,-0.332500170441278D+3,
     &  0.642859598466067D+2,0.773845935768222D+3,
     &  -0.152313732937084D+4 /
      DATA C3QUX / 0.565603648239126D+3,0.529062258221222D+1,
     &  -0.102020639611016D+0,0.122240301070145D-2 /
      DATA C3RXX / 0.584561202520006D+3,-0.102961025163669D+1,
     &  0.243293362700452D+0,-0.294905044740799D-2 /
      DATA C3UVX / 0.528199646263062D+3,0.890579602135307D+1,
     &  -0.222814134903755,0.286791682263697D-2 /
      DATA C3WXX / 0.728052609145380D+1,0.973505869861952D+2,
     &  0.147370491183191D+2,0.329196213998375D+3,
     &  0.873371668682417D+3 /
      DATA IX / 0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,
     &  5,5,5,6,6,6,7,8,9,9,10,10,11 /
      DATA JX / 0,0,1,2,7,10,12,23,2,6,15,17,0,2,6,7,22,26,0,2,4,16,26,
     &  0,2,4,26,1,3,26,0,2,26,2,26,2,26,0,1,26 /
      DATA I3ABX / 0,1,2,-1,-2 /
      DATA I3CDX / 0,1,2,3 /
      DATA I3GHX / 0,1,2,3,4 /
      DATA I3IJX / 0,1,2,3,4 /
      DATA I3JKX / 0,1,2,3,4 /
      DATA I3MNX / 0,1,2,3 /
      DATA I3OPX / 0,1,2,-1,-2 /
      DATA I3QUX / 0,1,2,3 /
      DATA I3RXX / 0,1,2,3 /
      DATA I3UVX / 0,1,2,3 /
      DATA I3WXX / 0,1,2,-1,-2 /
      DATA P3X / 16.5292D+6,19.00881189173929D+6,20.5D+6,
     &  21.04336732D+6,21.90096265D+6,21.93161551D+6,22.064D+6,
     &  22.11D+6,22.5D+6,23.D+6,23.5D+6,25.D+6,40.D+6,100.D+6 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3SUBCRG'
      TKX = TX + TABS
      PRX = 1.D+6
      PIX = PX/PRX
      INDX = 1
      CALL REGION_23( TB23X,PX,INDX )
      TB23X = TB23X + TABS
!
!---    ef subregion boundary  ---
!
        T3EFX = 3.727888004D+0*(PIX-22.064) + 647.096
!
!---    jk subregion boundary  ---
!
        T3JKX = 0.D+0
        DO 140 I = 1,5
          T3JKX = T3JKX + C3JKX(I)*(PIX**I3JKX(I))
  140   CONTINUE
!
!---    rx subregion boundary  ---
!
        T3RXX = 0.D+0
        DO 180 I = 1,4
          T3RXX = T3RXX + C3RXX(I)*(PIX**I3RXX(I))
  180   CONTINUE
!
!---    wx subregion boundary  ---
!
        T3WXX = 0.D+0
        DO 200 I = 1,5
          T3WXX = T3WXX + C3WXX(I)*(LOG(PIX)**I3WXX(I))
  200   CONTINUE
!
!---    Region 3k  ---
!
        IF( PX.GT.P3X(3) .AND. PX.LE.P3X(12) .AND. 
     &    TKX.GT.T3JKX .AND. TKX.LE.TB23X ) THEN
          CALL REGION_3K( TX,PX,RHOX )
!
!---    Region 3r  ---
!
        ELSEIF( PX.GT.P3X(3) .AND. PX.LE.P3X(9) .AND. 
     &    TKX.GT.T3RXX .AND. TKX.LE.T3JKX ) THEN
          CALL REGION_3R( TX,PX,RHOX )
!
!---    Region 3t  ---
!
        ELSEIF( PX.GT.P3X(1) .AND. PX.LE.P3X(3) .AND. 
     &    TKX.LE.TB23X ) THEN
          CALL REGION_3T( TX,PX,RHOX )
!
!---    Region 3x  ---
!
        ELSEIF( ( PX.GT.P3X(4) .AND. PX.LE.P3X(5) .AND. 
     &    TKX.LE.T3RXX ) .OR. 
     &    ( PX.GT.P3X(5) .AND. PX.LE.P3X(9) .AND. 
     &    TKX.GT.T3WXX .AND. TKX.LE.T3RXX ) ) THEN
          CALL REGION_3X( TX,PX,RHOX )
!
!---    Region 3z  ---
!
        ELSEIF( ( PX.GT.P3X(5) .AND. PX.LE.P3X(7) .AND. 
     &    TKX.LE.T3WXX ) .OR. 
     &    ( PX.GT.P3X(7) .AND. PX.LE.P3X(8) .AND. 
     &    TKX.GT.T3EFX .AND. TKX.LE.T3WXX ) ) THEN
          CALL REGION_3Z( TX,PX,RHOX )
!
!---    Region 3 subregion not identified  ---
!
        ELSE
          INDX = 12
          IMSGX = N_DB
          NMSGX = 0
          SUBLOGX = 'REGION_3SUBCRG'
          RLMSGX = 0.D+0
          CHMSGX = 'Region 3 Subregion Not Identified: ' // 
     &      ' @ Node'
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
!
!---  First partial derivatives of the Helmholtz free energy
!     with respect to the reduced density, "delta" and reduced
!     temperature, "tau" ---
!
      DELTAX = RHOX/RHOCRX
      TAUX = TCRX/TKX
      PHI_DELTAX = CNX(1)/DELTAX
      PHI_TAUX = 0.D+0
      DO 300 I = 2,40
        PHI_DELTAX = PHI_DELTAX + CNX(I)*REAL(IX(I))*
     &    (DELTAX**(IX(I)-1))*(TAUX**JX(I))
        PHI_TAUX = PHI_TAUX + CNX(I)*(DELTAX**IX(I))*REAL(JX(I))*
     &    (TAUX**(JX(I)-1))
  300 CONTINUE
!
!---  Enthalpy, J/kg  ---
!
      HX = RCX*TKX*(TAUX*PHI_TAUX + DELTAX*PHI_DELTAX)
!
!---  Internal energy, J/kg  ---
!
      UX = RCX*TKX*TAUX*PHI_TAUX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3SUBCRG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3SUBCRL( TX,PX,RHOX,HX,UX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density, enthalpy, and internal energy for
!     Region 3 (subcritical liquid) per the Revised Release on the  
!     IAPWS Industrial Formulation 1997 for the Thermodynamic Properties  
!     of Water andSteam (The revision only relates to the extension 
!     of Region 5 to 50 MPa), as a function of temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 8 April 2015.
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
      REAL*8 CNX(40)
      REAL*8 C3ABX(5),C3CDX(4),C3GHX(5),C3IJX(5),C3JKX(5),C3MNX(4)
      REAL*8 C3OPX(5),C3QUX(4),C3RXX(4),C3UVX(4),C3WXX(5)
      REAL*8 P3X(14)
      INTEGER IX(40),JX(40)
      INTEGER I3ABX(5),I3CDX(4),I3GHX(5),I3IJX(5),I3JKX(5),I3MNX(4)
      INTEGER I3OPX(5),I3QUX(4),I3RXX(4),I3UVX(4),I3WXX(5)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,C3ABX,C3CDX,C3GHX,C3IJX,C3JKX,C3MNX,C3OPX,C3QUX,C3RXX
      SAVE C3UVX,C3WXX
      SAVE P3X
      SAVE IX,JX,I3ABX,I3CDX,I3GHX,I3JKX,I3MNX,I3OPX,I3QUX,I3RXX
      SAVE I3UVX,I3WXX
      SAVE TCRX,PCRX,RHOCRX,RCX
      DATA TCRX / 647.096D+0 /
      DATA PCRX / 22.064D+6 /
      DATA RHOCRX / 322.D+0 /
      DATA RCX / 0.461526D+3 /
      DATA CNX / 0.10658070028513D+1,-0.15732845290239D+2,
     &  0.20944396974307D+2,-0.76867707878716D+1,0.26185947787954D+1,
     & -0.28080781148620D+1,0.12053369696517D+1,-0.84566812812502D-2,
     & -0.12654315477714D+1,-0.11524407806681D+1,0.88521043984318D+0,
     & -0.64207765181607D+0,0.38493460186671D+0,-0.85214708824206D+0,
     &  0.48972281541877D+1,-0.30502617256965D+1,0.39420536879154D-1,
     &  0.12558408424308D+0,-0.27999329698710D+0,0.13899799569460D+1,
     & -0.20189915023570D+1,-0.82147637173963D-2,-0.47596035734923D+0,
     &  0.43984074473500D-1,-0.44476435428739D+0,0.90572070719733D+0,
     &  0.70522450087967D+0,0.10770512626332D+0,-0.32913623258954D+0,
     & -0.50871062041158D+0,-0.22175400873096D-1,0.94260751665092D-1,
     &  0.16436278447961D+0,-0.13503372241348D-1,-0.14834345352472D-1,
     &  0.57922953628084D-3,0.32308904703711D-2,0.80964802996215D-4,
     & -0.16557679795037D-3,-0.44923899061815D-4 /
      DATA C3ABX / 0.154793642129415D+4,-0.187661219490113D+3,
     &  0.213144632222113D+2,-0.191887498864292D+4,
     &  0.918419702359447D+3 /
      DATA C3CDX / 0.585276966696349D+3,0.278233532206915D+1,
     &  -0.127283549295878D-1,0.159090746562729D-3 /
      DATA C3GHX / -0.249284240900418D+5,0.428143584791546D+4,
     &  -0.269029173140130D+3,0.751608051114157D+1,
     &  -0.787105249910383D-1 /
      DATA C3IJX / 0.584814781649163D+3,-0.616179320924617D+0,
     &  0.260763050899562D+0,-0.587071076864459D-2,
     &  0.515308185433082D-4 /
      DATA C3JKX / 0.617229772068439D+3,-0.770600270141675D+1,
     &  0.697072596851896D+0,-0.157391839848015D-1,
     &  0.137897492684194D-3 /
      DATA C3MNX / 0.535339483742384D+3,0.761978122720128D+1,
     &  -0.158365725441648,0.192871054508108D-2 /
      DATA C3OPX / 0.969461372400213D+3,-0.332500170441278D+3,
     &  0.642859598466067D+2,0.773845935768222D+3,
     &  -0.152313732937084D+4 /
      DATA C3QUX / 0.565603648239126D+3,0.529062258221222D+1,
     &  -0.102020639611016D+0,0.122240301070145D-2 /
      DATA C3RXX / 0.584561202520006D+3,-0.102961025163669D+1,
     &  0.243293362700452D+0,-0.294905044740799D-2 /
      DATA C3UVX / 0.528199646263062D+3,0.890579602135307D+1,
     &  -0.222814134903755,0.286791682263697D-2 /
      DATA C3WXX / 0.728052609145380D+1,0.973505869861952D+2,
     &  0.147370491183191D+2,0.329196213998375D+3,
     &  0.873371668682417D+3 /
      DATA IX / 0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,
     &  5,5,5,6,6,6,7,8,9,9,10,10,11 /
      DATA JX / 0,0,1,2,7,10,12,23,2,6,15,17,0,2,6,7,22,26,0,2,4,16,26,
     &  0,2,4,26,1,3,26,0,2,26,2,26,2,26,0,1,26 /
      DATA I3ABX / 0,1,2,-1,-2 /
      DATA I3CDX / 0,1,2,3 /
      DATA I3GHX / 0,1,2,3,4 /
      DATA I3IJX / 0,1,2,3,4 /
      DATA I3JKX / 0,1,2,3,4 /
      DATA I3MNX / 0,1,2,3 /
      DATA I3OPX / 0,1,2,-1,-2 /
      DATA I3QUX / 0,1,2,3 /
      DATA I3RXX / 0,1,2,3 /
      DATA I3UVX / 0,1,2,3 /
      DATA I3WXX / 0,1,2,-1,-2 /
      DATA P3X / 16.5292D+6,19.00881189173929D+6,20.5D+6,
     &  21.04336732D+6,21.90096265D+6,21.93161551D+6,22.064D+6,
     &  22.11D+6,22.5D+6,23.D+6,23.5D+6,25.D+6,40.D+6,100.D+6 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3SUBCRL'
      TKX = TX + TABS
      PRX = 1.D+6
      PIX = PX/PRX
!
!---    cd subregion boundary  ---
!
        T3CDX = 0.D+0
        DO 110 I = 1,4
          T3CDX = T3CDX + C3CDX(I)*(PIX**I3CDX(I))
  110   CONTINUE
!
!---    ef subregion boundary  ---
!
        T3EFX = 3.727888004D+0*(PIX-22.064) + 647.096
!
!---    qu subregion boundary  ---
!
        T3QUX = 0.D+0
        DO 170 I = 1,4
          T3QUX = T3QUX + C3QUX(I)*(PIX**I3QUX(I))
  170   CONTINUE
!
!---    uv subregion boundary  ---
!
        T3UVX = 0.D+0
        DO 190 I = 1,4
          T3UVX = T3UVX + C3UVX(I)*(PIX**I3UVX(I))
  190   CONTINUE
!
!---    Region 3c  ---
!
        IF( ( PX.GT.P3X(2) .AND. PX.LE.P3X(13) .AND. 
     &    TKX.LE.T3CDX ) .OR. ( PX.GT.P3X(1) .AND. PX.LE.P3X(2) ) ) THEN
          CALL REGION_3C( TX,PX,RHOX )
!
!---    Region 3q  ---
!
        ELSEIF( PX.GT.P3X(4) .AND. PX.LE.P3X(9) .AND. 
     &    TKX.GT.T3CDX .AND. TKX.LE.T3QUX ) THEN
          CALL REGION_3Q( TX,PX,RHOX )
!
!---    Region 3s  ---
!
        ELSEIF( PX.GT.P3X(2) .AND. PX.LE.P3X(4) .AND. 
     &    TKX.GT.T3CDX ) THEN
          CALL REGION_3S( TX,PX,RHOX )
!
!---    Region 3u  ---
!
        ELSEIF( ( PX.GT.P3X(4) .AND. PX.LE.P3X(6) .AND. 
     &    TKX.GT.T3QUX ) .OR. 
     &    ( PX.GT.P3X(6) .AND. PX.LE.P3X(9) .AND. 
     &    TKX.GT.T3QUX .AND. TKX.LE.T3UVX ) ) THEN
          CALL REGION_3U( TX,PX,RHOX )
!
!---    Region 3y  ---
!
        ELSEIF( ( PX.GT.P3X(6) .AND. PX.LE.P3X(7) .AND. 
     &    TKX.GT.T3UVX ) .OR. 
     &    ( PX.GT.P3X(7) .AND. PX.LE.P3X(8) .AND. 
     &    TKX.GT.T3UVX .AND. TKX.LE.T3EFX ) ) THEN
          CALL REGION_3Y( TX,PX,RHOX )
!
!---    Region 3 subregion not identified  ---
!
        ELSE
          INDX = 12
          IMSGX = N_DB
          NMSGX = 0
          SUBLOGX = 'REGION_3SUBCRL'
          RLMSGX = 0.D+0
          CHMSGX = 'Region 3 Subregion Not Identified: ' // 
     &      ' @ Node'
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
!
!---  First partial derivatives of the Helmholtz free energy
!     with respect to the reduced density, "delta" and reduced
!     temperature, "tau" ---
!
      DELTAX = RHOX/RHOCRX
      TAUX = TCRX/TKX
      PHI_DELTAX = CNX(1)/DELTAX
      PHI_TAUX = 0.D+0
      DO 300 I = 2,40
        PHI_DELTAX = PHI_DELTAX + CNX(I)*REAL(IX(I))*
     &    (DELTAX**(IX(I)-1))*(TAUX**JX(I))
        PHI_TAUX = PHI_TAUX + CNX(I)*(DELTAX**IX(I))*REAL(JX(I))*
     &    (TAUX**(JX(I)-1))
  300 CONTINUE
!
!---  Enthalpy, J/kg  ---
!
      HX = RCX*TKX*(TAUX*PHI_TAUX + DELTAX*PHI_DELTAX)
!
!---  Internal energy, J/kg  ---
!
      UX = RCX*TKX*TAUX*PHI_TAUX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3SUBCRL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3A( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(30)
      INTEGER IX(30),JX(30)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-12,-12,-10,-10,-10,-8,-8,-8,-6,-5,-5,-5,-4,-3,-3,
     &  -3,-3,-2,-2,-2,-1,-1,-1,0,0,1,1,2,2 /
      DATA JX / 5,10,12,5,10,12,5,8,10,1,1,5,10,8,0,1,3,6,0,2,3,0,1,2,
     &  0,1,0,2,0,2 /
      DATA CNX / 0.110879558823853D-2,0.572616740810616D+3,
     &  -0.767051948380852D+5,-0.253321069529674D-1,
     &  0.628008049345689D+4,0.234105654131876D+6,
     &  0.216867826045856D+0,-0.156237904341963D+3,
     &  -0.269893956176613D+5,-0.180407100085505D-3,
     &  0.116732227668261D-2,0.266987040856040D+2,
     &  0.282776617243286D+5,-0.242431520029523D+4,
     &  0.435217323022733D-3,-0.122494831387441D-1,
     &  0.179357604019989D+1,0.442729521058314D+2,
     &  -0.593223489018342D-2,0.453186261685774D+0,
     &  0.135825703129140D+1,0.408748415856745D-1,
     &  0.474686397863312D+0,0.118646814997915D+1,
     &  0.546987265727549D+0,0.195266770452643D+0,
     &  -0.502268790869663D-1,-0.369645308193377D+0,
     &  0.633828037528420D-2,0.797441793901017D-1 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3A'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 100.D+6
      TRX = 760.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.085D+0
      BX = 0.817D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,30
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0024D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3A group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3B( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(32)
      INTEGER IX(32),JX(32)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-12,-10,-10,-8,-6,-6,-6,-5,-5,-5,-4,-4,-4,-3,-3,
     &  -3,-3,-3,-2,-2,-2,-1,-1,0,0,1,1,2,3,4,4 /
      DATA JX / 10,12,8,14,8,5,6,8,5,8,10,2,4,5,0,1,2,3,5,0,2,5,0,2,
     &  0,1,0,2,0,2,0,1 /
      DATA CNX / -0.827670470003621D-1,0.416887126010565D+2,
     &  0.483651982197059D-1,-0.291032084950276D+5,
     &  -0.111422582236948D+3,-0.202300083904014D-1,
     &  0.294002509338515D+3,0.140244997609658D+3,
     &  -0.344384158811459D+3,0.361182452612149D+3,
     &  -0.140699677420738D+4,-0.202023902676481D-2,
     &  0.171346792457471D+3,-0.425597804058632D+1,
     &  0.691346085000334D-5,0.151140509678925D-2,
     &  -0.416375290166236D-1,-0.413754957011042D+2,
     &  -0.506673295721637D+2,-0.572212965569023D-3,
     &  0.608817368401785D+1,0.239600660256161D+2,
     &  0.122261479925384D-1,0.216356057692938D+1,
     &  0.398198903368642D+0,-0.116892827834085D+0,
     &  -0.102845919373532D+0,-0.492676637589284D+0,
     &  0.655540456406790D-1,-0.240462535078530D+0,
     &  -0.269798180310075D-1,0.128369435967012D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3B'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 100.D+6
      TRX = 860.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.280D+0
      BX = 0.779D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,32
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0041D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3C( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(35)
      INTEGER IX(35),JX(35)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-12,-12,-10,-10,-10,-8,-8,-8,-6,-5,-5,-5,-4,-4,-3,
     &  -3,-2,-2,-2,-1,-1,-1,0,0,0,1,1,2,2,2,2,3,3,8 /
      DATA JX / 6,8,10,6,8,10,5,6,7,8,1,4,7,2,8,0,3,0,4,5,0,1,2,0,1,2,
     &  0,2,0,1,3,7,0,7,1 /
      DATA CNX / 3.11967788763030D+0,2.76713458847564D+04,
     &  3.22583103403269D+07,-3.42416065095363D+02,
     &  -8.99732529907377D+05,-7.93892049821251D+07,
     &  9.53193003217388D+01,2.29784742345072D+03,
     &  1.75336675322499D+05,7.91214365222792D+06,
     &  3.19933345844209D-05,-6.59508863555767D+01,
     &  -8.33426563212851D+05,6.45734680583292D-02,
     &  -3.82031020570813D+06,4.06398848470079D-05,
     &  3.10327498492008D+01,-8.92996718483724D-04,
     &  2.34604891591616D+02,3.77515668966951D+03,
     &  1.58646812591361D-02,7.07906336241843D-01,
     &  1.26016225146570D+01,7.36143655772152D-01,
     &  6.76544268999101D-01,-1.78100588189137D+01,
     &  -1.56531975531713D-01,1.17707430048158D+01,
     &  8.40143653860447D-02,-1.86442467471949D-01,
     &  -4.40170203949645D+01,1.23290423502494D+06,
     &  -2.40650039730845D-02,-1.07077716660869D+06,
     &  4.38319858566475D-02 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3C'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 40.D+6
      TRX = 690.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.259D+0
      BX = 0.903D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,35
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0022D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3C group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3D( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(38)
      INTEGER IX(38),JX(38)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-12,-12,-12,-12,-12,-10,-10,-10,-10,-10,-10,-10,
     &  -8,-8,-8,-8,-6,-6,-5,-5,-5,-5,-4,-4,-4,-3,-3,-2,-2,-1,-1,-1,
     &  0,0,1,1,3 /
      DATA JX / 4,6,7,10,12,16,0,2,4,6,8,10,14,3,7,8,10,6,8,
     &  1,2,5,7,0,1,7,2,4,0,1,0,1,5,0,2,0,6,0 /
      DATA CNX / -4.52484847171645D-10,3.15210389538801D-05,
     &  -2.14991352047545D-03,5.08058874808345D+02,
     &  -1.27123036845932D+07,1.15371133120497D+12,
     &  -1.97805728776273D-16,2.41554806033972D-11,
     &  -1.56481703640525D-06,2.77211346836625D-03,
     &  -2.03578994462286D+01,1.44369489909053D+06,
     &  -4.11254217946539D+10,6.23449786243773D-06,
     &  -2.21774281146038D+01,-6.89315087933158D+04,
     &  -1.95419525060713D+07,3.16373510564015D+03,
     &  2.24040754426988D+06,-4.36701347922356D-06,
     &  -4.04213852833996D-04,-3.48153203414663D+02,
     &  -3.85294213555289D+05,1.35203700099403D-07,
     &  1.34648383271089D-04,1.25031835351736D+05,
     &  9.68123678455841D-02,2.25660517512438D+02,
     &  -1.90102435341872D-04,-2.99628410819229D-02,
     &  5.00833915372121D-03,3.87842482998411D-01,
     &  -1.38535367777182D+03,8.70745245971773D-01,
     &  1.71946252068742D+00,-3.26650121426383D-02,
     &  4.98044171727877D+03,5.51478022765087D-03 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3D'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 40.D+6
      TRX = 690.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.559D+0
      BX = 0.939D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,38
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
      SPC_VOLX = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0029D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3D group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3E( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(29)
      INTEGER IX(29),JX(29)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-12,-10,-10,-10,-10,-10,-8,-8,-8,-6,-5,-4,-4,-3,
     &  -3,-3,-2,-2,-2,-2,-1,0,0,1,1,1,2,2 /
      DATA JX / 14,16,3,6,10,14,16,7,8,10,6,6,2,4,2,6,7,0,1,3,4,0,0,
     &  1,0,4,6,0,2 /
      DATA CNX / 7.15815808404721D+08,-1.14328360753449D+11,
     &  3.76531002015720D-12,-9.03983668691157D-05,
     &  6.65695908836252D+05,5.35364174960127D+09,
     &  7.94977402335603D+10,9.22230563421437D+01,
     &  -1.42586073991215D+05,-1.11796381424162D+06,
     &  8.96121629640760D+03,-6.69989239070491D+03,
     &  4.51242538486834D-03,-3.39731325977713D+01,
     &  -1.20523111552278D+00,4.75992667717124D+04,
     &  -2.66627750390341D+05,-1.53314954386524D-04,
     &  3.05638404828265D-01,1.23654999499486D+02,
     &  -1.04390794213011D+03,-1.57496516174308D-02,
     &  6.85331118940253D-01,1.78373462873903D+00,
     &  -5.44674124878910D-01,2.04529931318843D+03,
     &  -2.28342359328752D+04,4.13197481515899D-01,
     &  -3.41931835910405D+01 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3E'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 40.D+6
      TRX = 710.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.587D+0
      BX = 0.918D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,29
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0032D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3E group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3F( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(42)
      INTEGER IX(42),JX(42)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / 0,0,0,0,0,0,1,1,1,1,2,2,3,3,3,4,5,5,6,7,7,10,12,12,
     &  12,14,14,14,14,14,16,16,18,18,20,20,20,22,24,24,28,32 /
      DATA JX / -3,-2,-1,0,1,2,-1,1,2,3,0,1,-5,-2,0,-3,-8,1,-6,
     &  -4,1,-6,-10,-8,-4,-12,-10,-8,-6,-4,-10,-8,-12,-10,-12,
     &  -10,-6,-12,-12,-4,-12,-12 /
      DATA CNX / -2.51756547792325D-08,6.01307193668763D-06,
     &  -1.00615977450049D-03,9.99969140252192D-01,2.14107759236486D+00,
     &  -1.65175571959086D+01,-1.41987303638727D-03,
     &  2.69251915156554D+00,3.49741815858722D+01,
     &  -3.00208695771783D+01,-1.31546288252539D+00,
     &  -8.39091277286169D+00,1.81545608337015D-10,
     &  -5.91099206478909D-04,1.52115067087106D+00,
     &  2.52956470663225D-05,1.00726265203786D-15,
     &  -1.49774533860650D+00,-7.93940970562969D-10,
     &  -1.50290891264717D-04,1.51205531275133D+00,
     &  4.70942606221652D-06,1.95049710391712D-13,
     &  -9.11627886266077D-09,6.04374640201265D-04,
     &  -2.25132933900136D-16,6.10916973582981D-12,
     &  -3.03063908043404D-07,-1.37796070798409D-05,
     &  -9.19296736666106D-04,6.39288223132545D-10,
     &  7.53259479898699D-07,-4.00321478682929D-13,
     &  7.56140294351614D-09,-9.12082054034891D-12,
     &  -2.37612381140539D-08,2.69586010591874D-05,
     &  -7.32828135157839D-11,2.41995578306660D-10,
     &  -4.05735532730322D-04,1.89424143498011D-10,
     &  -4.86632965074563D-10 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3F'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 40.D+6
      TRX = 730.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.587D+0
      BX = 0.891D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,42
        SPC_VOLX  = SPC_VOLX + CNX(I)*((SQRT(PIX-AX))**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
      SPC_VOLX = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0064D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3F group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3G( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(38)
      INTEGER IX(38),JX(38)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-12,-12,-12,-12,-12,-10,-10,-10,-8,-8,-8,-8,-6,
     &  -6,-5,-5,-4,-3,-2,-2,-2,-2,-1,-1,-1,0,0,0,1,1,1,3,5,6,8,10,10 /
      DATA JX / 7,12,14,18,22,24,14,20,24,7,8,10,12,8,22,7,20,22,7,3,
     &  5,14,24,2,8,18,0,1,2,0,1,3,24,22,12,3,0,6 /
      DATA CNX / 4.12209020652996D-05,-1.14987238280587D+06,
     &  9.48180885032080D+09,-1.95788865718971D+17,4.96250704871300D+24,
     &  -1.05549884548496D+28,-7.58642165988278D+11,
     &  -9.22172769596101D+22,7.25379072059348D+29,
     &  -6.17718249205859D+01,1.07555033344858D+04,
     &  -3.79545802336487D+07,2.28646846221831D+11,
     &  -4.99741093010619D+06,-2.80214310054101D+30,
     &  1.04915406769586D+06,6.13754229168619D+27,
     &  8.02056715528378D+31,-2.98617819828065D+07,
     &  -9.10782540134681D+01,1.35033227281565D+05,
     &  -7.12949383408211D+18,-1.04578785289542D+36,
     &  3.04331584444093D+01,5.93250797959445D+09,
     &  -3.64174062110798D+27,9.21791403532461D-01,
     &  -3.37693609657471D-01,-7.24644143758508D+01,
     &  -1.10480239272601D-01,5.36516031875059D+00,
     &  -2.91441872156205D+03,6.16338176535305D+39,
     &  -1.20889175861180D+38,8.18396024524612D+22,
     &  9.40781944835829D+08,-3.67279669545448D+04,
     &  -8.37513931798655D+15 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3G'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 25.D+6
      TRX = 660.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.872D+0
      BX = 0.971D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,38
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
      SPC_VOLX = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0027D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3G group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3H( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(29)
      INTEGER IX(29),JX(29)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-12,-10,-10,-10,-10,-10,-10,-8,-8,-8,-8,-8,-6,
     &  -6,-6,-5,-5,-5,-4,-4,-3,-3,-2,-1,-1,0,1,1 /
      DATA JX / 8,12,4,6,8,10,14,16,0,1,6,7,8,4,6,8,2,3,4,2,4,1,2,
     &  0,0,2,0,0,2 /
      DATA CNX / 5.61379678887577D-02,7.74135421587083D+09,
     &  1.11482975877938D-09,-1.43987128208183D-03,1.93696558764920D+03,
     &  -6.05971823585005D+08,1.71951568124337D+13,
     &  -1.85461154985145D+16,3.87851168078010D-17,
     &  -3.95464327846105D-14,-1.70875935679023D+02,
     &  -2.12010620701220D+03,1.77683337348191D+07,
     &  1.10177443629575D+01,-2.34396091693313D+05,
     &  -6.56174421999594D+06,1.56362212977396D-05,
     &  -2.12946257021400D+00,1.35249306374858D+01,
     &  1.77189164145813D-01,1.39499167345464D+03,
     &  -7.03670932036388D-03,-1.52011044389648D-01,
     &  9.81916922991113D-05,1.47199658618076D-03,
     &  2.02618487025578D+01,8.99345518944240D-01,
     &  -2.11346402240858D-01,2.49971752957491D+01 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3H'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 25.D+6
      TRX = 660.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.898D+0
      BX = 0.983D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,29
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
      SPC_VOLX = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0032D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3H group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3I( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(42)
      INTEGER IX(42),JX(42)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / 0,0,0,1,1,1,1,2,3,3,4,4,4,5,5,5,7,7,8,8,10,12,12,12,
     &  14,14,14,14,18,18,18,18,18,20,20,22,24,24,32,32,36,36 /
      DATA JX / 0,1,10,-4,-2,-1,0,0,-5,0,-3,-2,-1,-6,-1,12,-4,-3,-6,
     &  10,-8,-12,-6,-4,-10,-8,-4,5,-12,-10,-8,-6,2,-12,-10,-12,-12,
     &  -8,-10,-5,-10,-8 /
      DATA CNX / 1.06905684359136D+00,-1.48620857922333D+00,
     &  2.59862256980408D+14,-4.46352055678749D-12,
     &  -5.66620757170032D-07,-2.35302885736849D-03,
     &  -2.69226321968839D-01,9.22024992944392D+00,
     &  3.57633505503772D-12,-1.73942565562222D+01,
     &  7.00681785556229D-06,-2.67050351075768D-04,
     &  -2.31779669675624D+00,-7.53533046979752D-13,
     &  4.81337131452891D+00,-2.23286270422356D+21,
     &  -1.18746004987383D-05,6.46412934136496D-03,
     &  -4.10588536330937D-10,4.22739537057241D+19,
     &  3.13698180473812D-13,1.64395334345040D-24,
     &  -3.39823323754373D-06,-1.35268639905021D-02,
     &  -7.23252514211625D-15,1.84386437538366D-09,
     &  -4.63959533752385D-02,-9.92263100376750D+13,
     &  6.88169154439335D-17,-2.22620998452197D-11,
     &  -5.40843018624083D-08,3.45570606200257D-03,
     &  4.22275800304086D+10,-1.26974478770487D-15,
     &  9.27237985153679D-10,6.12670812016489D-14,
     &  -7.22693924063497D-12,-3.83669502636822D-04,
     &  3.74684572410204D-04,-9.31976897511086D+04,
     &  -2.47690616026922D-02,6.58110546759474D+01 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3I'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 25.D+6
      TRX = 660.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.910D+0
      BX = 0.984D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,42
        SPC_VOLX  = SPC_VOLX + CNX(I)*((SQRT(PIX-AX))**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
      SPC_VOLX = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0041D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3I group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3J( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(29)
      INTEGER IX(29),JX(29)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / 0,0,0,1,1,1,2,2,3,4,4,5,5,5,6,10,12,12,14,14,14,16,
     &  18,20,20,24,24,28,28 /
      DATA JX / -1,0,1,-2,-1,1,-1,1,-2,-2,2,-3,-2,0,3,-6,-8,-3,-10,
     &  -8,-5,-10,-12,-12,-10,-12,-6,-12,-5 /
      DATA CNX / -1.11371317395540D-04,1.00342892423685D+00,
     &  5.30615581928979D+00,1.79058760078792D-06,-7.28541958464774D-04,
     &  -1.87576133371704D+01,1.99060874071849D-03,2.43574755377290D+01,
     &  -1.77040785499444D-04,-2.59680385227130D-03,
     &  -1.98704578406823D+02,7.38627790224287D-05,
     &  -2.36264692844138D-03,-1.61023121314333D+00,
     &  6.22322971786473D+03,-9.60754116701669D-09,
     &  -5.10572269720488D-11,7.67373781404211D-03,
     &  6.63855469485254D-15,-7.17590735526745D-10,
     &  1.46564542926508D-05,3.09029474277013D-12,
     &  -4.64216300971708D-16,-3.90499637961161D-14,
     &  -2.36716126781431D-10,4.54652854268717D-12,
     &  -4.22271787482497D-03,2.83911742354706D-11,
     &  2.70929002720228D+00 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3J'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 25.D+6
      TRX = 670.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.875D+0
      BX = 0.964D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,29
        SPC_VOLX  = SPC_VOLX + CNX(I)*((SQRT(PIX-AX))**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
      SPC_VOLX = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0054D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3J group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3K( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(34)
      INTEGER IX(34),JX(34)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -2,-2,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,2,
     &  5,5,5,6,6,6,6,8,10,12 /
      DATA JX / 10,12,-5,6,-12,-6,-2,-1,0,1,2,3,14,-3,-2,0,1,2,-8,
     &  -6,-3,-2,0,4,-12,-6,-3,-12,-10,-8,-5,-12,-12,-10 /
      DATA CNX / -4.01215699576099D+08,4.84501478318406D+10,
     &  3.94721471363678D-15,3.72629967374147D+04,-3.69794374168666D-30,
     &  -3.80436407012452D-15,4.75361629970233D-07,
     &  -8.79148916140706D-04,8.44317863844331D-01,
     &  1.22433162656600D+01,-1.04529634830279D+02,
     &  5.89702771277429D+02,-2.91026851164444D+13,
     &  1.70343072841850D-06,-2.77617606975748D-04,
     &  -3.44709605486686D+00,2.21333862447095D+01,
     &  -1.94646110037079D+02,8.08354639772825D-16,
     &  -1.80845209145470D-11,-6.96664158132412D-06,
     &  -1.81057560300994D-03,2.55830298579027D+00,
     &  3.28913873658481D+03,-1.73270241249904D-19,
     &  -6.61876792558034D-07,-3.95688923421250D-03,
     &  6.04203299819132D-18,-4.00879935920517D-14,
     &  1.60751107464958D-09,3.83719409025556D-05,
     &  -6.49565446702457D-15,-1.49095328506000D-12,
     &  5.41449377329581D-09 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3K'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 25.D+6
      TRX = 680.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.802D+0
      BX = 0.935D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,34
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0077D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3K group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3L( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(43)
      INTEGER IX(43),JX(43)
!
!----------------------Data Statements---------------------------------!
!
      DATA IX / -12,-12,-12,-12,-12,-10,-10,-8,-8,-8,-8,-8,-8,-8,
     &  -6,-5,-5,-4,-4,-3,-3,-3,-3,-2,-2,-2,-1,-1,-1,0,0,0,0,1,1,
     &  2,4,5,5,6,10,10,14 /
      DATA JX / 14,16,18,20,22,14,24,6,10,12,14,18,24,36,8,4,5,7,
     &  16,1,3,18,20,2,3,10,0,1,3,0,1,2,12,0,16,1,0,0,1,14,4,12,10 /
      DATA CNX / 2.60702058647537D+09,-1.88277213604704D+14,
     &  5.54923870289667D+18,-7.58966946387758D+22,
     &  4.13865186848908D+26,-8.15038000738060D+11,
     &  -3.81458260489955D+32,-1.23239564600519D-02,
     &  2.26095631437174D+07,-4.95017809506720D+11,
     &  5.29482996422863D+15,-4.44359478746295D+22,
     &  5.21635864527315D+34,-4.87095672740742D+54,
     &  -7.14430209937547D+05,1.27868634615495D-01,
     &  -1.00752127917598D+01,7.77451437960990D+06,
     &  -1.08105480796471D+24,-3.57578581169659D-06,
     &  -2.12857169423484D+00,2.70706111085238D+29,
     &  -6.95953622348829D+32,1.10609027472280D-01,
     &  7.21559163361354D+01,-3.06367307532219D+14,
     &  2.65839618885530D-05,2.53392392889754D-02,
     &  -2.14443041836579D+02,9.37846601489667D-01,
     &  2.23184043101700D+00,3.38401222509191D+01,
     &  4.94237237179718D+20,-1.98068404154428D-01,
     &  -1.41415349881140D+30,-9.93862421613651D+01,
     &  1.25070534142731D+02,-9.96473529004439D+02,
     &  4.73137909872765D+04,1.16662121219322D+32,
     &  -3.15874976271533D+15,-4.45703369196945D+32,
     &  6.42794932373694D+32 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3L'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 24.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.908D+0
      BX = 0.989D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,43
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
      SPC_VOLX = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0026D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3L group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3M( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(40)
      INTEGER IX(40),JX(40)
!
!----------------------Data Statements---------------------------------!
!
      DATA IX / 0,3,8,20,1,3,4,5,1,6,2,4,14,2,5,3,0,1,1,1,28,2,16,0,
     &  5,0,3,4,12,16,1,8,14,0,2,3,4,8,14,24 /
      DATA JX / 0,0,0,2,5,5,5,5,6,6,7,8,8,10,10,12,14,14,18,20,20,22,
     &  22,24,24,28,28,28,28,28,32,32,32,36,36,36,36,36,36,36 /
      DATA CNX / 8.11384363481847D-01,-5.68199310990094D+03,
     &  -1.78657198172556D+10,7.95537657613427D+31,
     &  -8.14568209346872D+04,-6.59774567602874D+07,
     &  -1.52861148659302D+10,-5.60165667510446D+11,
     &  4.58384828593949D+05,-3.85754000383848D+13,
     &  4.53735800004273D+07,9.39454935735563D+11,
     &  2.66572856432938D+27,-5.47578313899097D+09,
     &  2.00725701112386D+14,1.85007245563239D+12,
     &  1.85135446828337D+08,-1.70451090076385D+11,
     &  1.57890366037614D+14,-2.02530509748774D+15,
     &  3.68193926183570D+59,1.70215539458936D+17,
     &  6.39234909918741D+41,-8.21698160721956D+14,
     &  -7.95260241872306D+23,2.33415869478510D+17,
     &  -6.00079934586803D+22,5.94584382273384D+24,
     &  1.89461279349492D+39,-8.10093428842645D+45,
     &  1.88813911076809D+21,1.11052244098768D+35,
     &  2.91133958602503D+45,-3.29421923951460D+21,
     &  -1.37570282536696D+25,1.81508996303902D+27,
     &  -3.46865122768353D+29,-2.11961148774260D+37,
     &  -1.28617899887675D+48,4.79817895699239D+64 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3M'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 23.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 1.D+0
      BX = 0.997D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,40
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    (((THETAX-BX)**2.5D-1)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0028D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3M group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3N( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(39)
      INTEGER IX(39),JX(39)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / 0,3,4,6,7,10,12,14,18,0,3,5,6,8,12,0,3,7,12,2,3,4,2,
     &  4,7,4,3,5,6,0,0,3,1,0,1,0,1,0,1 /
      DATA JX / -12,-12,-12,-12,-12,-12,-12,-12,-12,-10,-10,-10,-10,
     &  -10,-10,-8,-8,-8,-8,-6,-6,-6,-5,-5,-5,-4,-3,-3,-3,-2,-1,-1,
     &  0,1,1,2,4,5,6 /
      DATA CNX / 2.80967799943151D-39,6.14869006573609D-31,
     &  5.82238667048942D-28,3.90628369238462D-23,
     &  8.21445758255119D-21,4.02137961842776D-15,
     &  6.51718171878301D-13,-2.11773355803058D-08,
     &  2.64953354380072D-03,-1.35031446451331D-32,
     &  -6.07246643970893D-24,-4.02352115234494D-19,
     &  -7.44938506925544D-17,1.89917206526237D-13,
     &  3.64975183508473D-06,1.77274872361946D-26,
     &  -3.34952758812999D-19,-4.21537726098389D-09,
     &  -3.91048167929649D-02,5.41276911564176D-14,
     &  7.05412100773699D-12,2.58585887897486D-09,
     &  -4.93111362030162D-11,-1.58649699894543D-06,
     &  -5.25037427886100D-01,2.20019901729615D-03,
     &  -6.43064132636925D-03,6.29154149015048D+01,
     &  1.35147318617061D+02,2.40560808321713D-07,
     &  -8.90763306701305D-04,-4.40209599407714D+03,
     &  -3.02807107747776D+02,1.59158748314599D+03,
     &  2.32534272709876D+05,-7.92681207132600D+05,
     &  -8.69871364662769D+10,3.54542769185671D+11,
     &  4.00849240129329D+14 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3N'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 23.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.976D+0
      BX = 0.997D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,39
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
      SPC_VOLX = EXP( SPC_VOLX )
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0031D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3N group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3O( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(24)
      INTEGER IX(24),JX(24)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / 0,0,0,2,3,4,4,4,4,4,5,5,6,7,8,8,8,10,10,14,14,20,20,24 /
      DATA JX / -12,-4,-1,-1,-10,-12,-8,-5,-4,-1,-4,-3,-8,-12,-10,-8,
     &  -4,-12,-8,-12,-8,-12,-10,-12 /
      DATA CNX / 1.28746023979718D-35,-7.35234770382342D-12,
     &  2.89078692149150D-03,2.44482731907223D-01,
     &  1.41733492030985D-24,-3.54533853059476D-29,
     &  -5.94539202901431D-18,-5.85188401782779D-09,
     &  2.01377325411803D-06,1.38647388209306D+00,
     &  -1.73959365084772D-05,1.37680878349369D-03,
     &  8.14897605805513D-15,4.25596631351839D-26,
     &  -3.87449113787755D-18,1.39814747930240D-13,
     &  -1.71849638951521D-03,6.41890529513296D-22,
     &  1.18960578072018D-11,-1.55282762571611D-18,
     &  2.33907907347507D-08,-1.74093247766213D-13,
     &  3.77682649089149D-09,-5.16720236575302D-11 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3O'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 23.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.974D+0
      BX = 0.996D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,24
        SPC_VOLX  = SPC_VOLX + CNX(I)*((SQRT(PIX-AX))**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0034D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3O group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3P( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(27)
      INTEGER IX(27),JX(27)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / 0,0,0,0,1,2,3,3,4,6,7,7,8,10,12,12,12,14,14,14,16,
     &  18,20,22,24,24,36 /
      DATA JX / -1,0,1,2,1,-1,-3,0,-2,-2,-5,-4,-2,-3,-12,-6,-5,-10,
     &  -8,-3,-8,-8,-10,-10,-12,-8,-12 /
      DATA CNX / -9.82825342010366D-05,1.05145700850612D+00,
     &  1.16033094095084D+02,3.24664750281543D+03,
     &  -1.23592348610137D+03,-5.61403450013495D-02,
     &  8.56677401640869D-08,2.36313425393924D+02,
     &  9.72503292350109D-03,-1.03001994531927D+00,
     &  -1.49653706199162D-09,-2.15743778861592D-05,
     &  -8.34452198291445D+00,5.86602660564988D-01,
     &  3.43480022104968D-26,8.16256095947021D-06,
     &  2.94985697916798D-03,7.11730466276584D-17,
     &  4.00954763806941D-10,1.07766027032853D+01,
     &  -4.09449599138182D-07,-7.29121307758902D-06,
     &  6.77107970938909D-09,6.02745973022975D-08,
     &  -3.82323011855257D-11,1.79946628317437D-03,
     &  -3.45042834640005D-04 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3P'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 23.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.972D+0
      BX = 0.997D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,27
        SPC_VOLX  = SPC_VOLX + CNX(I)*((SQRT(PIX-AX))**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0041D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3P group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3Q( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(24)
      INTEGER IX(24),JX(24)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-12,-10,-10,-10,-10,-8,-6,-5,-5,-4,-4,-3,-2,
     &  -2,-2,-2,-1,-1,-1,0,1,1,1 /
      DATA JX / 10,12,6,7,8,10,8,6,2,5,3,4,3,0,1,2,4,0,1,2,0,0,1,3 /
      DATA CNX / -8.20433843259950D+04,4.73271518461586D+10,
     &  -8.05950021005413D-02,3.28600025435980D+01,
     &  -3.56617029982490D+03,-1.72985781433335D+09,
     &  3.51769232729192D+07,-7.75489259985144D+05,
     &  7.10346691966018D-05,9.93499883820274D+04,
     &  -6.42094171904570D-01,-6.12842816820083D+03,
     &  2.32808472983776D+02,-1.42808220416837D-05,
     &  -6.43596060678456D-03,-4.28577227475614D+00,
     &  2.25689939161918D+03,1.00355651721510D-03,
     &  3.33491455143516D-01,1.09697576888873D+00,
     &  9.61917379376452D-01,-8.38165632204598D-02,
     &  2.47795908411492D+00,-3.19114969006533D+03 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3Q'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 23.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.848D+0
      BX = 0.983D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,24
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
      SPC_VOLX = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0022D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3Q group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3R( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(27)
      INTEGER IX(27),JX(27)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -8,-8,-3,-3,-3,-3,-3,0,0,0,0,3,3,8,8,8,8,10,10,10,10,
     &  10,10,10,10,12,14 /
      DATA JX / 6,14,-3,3,4,5,8,-1,0,1,5,-6,-2,-12,-10,-8,-5,-12,-10,
     &  -8,-6,-5,-4,-3,-2,-12,-12 /
      DATA CNX / 1.44165955660863D-03,-7.01438599628258D+12,
     &  -8.30946716459219D-17,2.61975135368109D-01,
     &  3.93097214706245D+02,-1.04334030654021D+04,
     &  4.90112654154211D+08,-1.47104222772069D-04,
     &  1.03602748043408D+00,3.05308890065089D+00,
     &  -3.99745276971264D+06,5.69233719593750D-12,
     &  -4.64923504407778D-02,-5.35400396512906D-18,
     &  3.99988795693162D-13,-5.36479560201811D-07,
     &  1.59536722411202D-02,2.70303248860217D-15,
     &  2.44247453858506D-08,-9.83430636716454D-06,
     &  6.63513144224454D-02,-9.93456957845006D+00,
     &  5.46491323528491D+02,-1.43365406393758D+04,
     &  1.50764974125511D+05,-3.37209709340105D-10,
     &  3.77501980025469D-09 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3R'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 23.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.874D+0
      BX = 0.982D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,27
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0054D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3R group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3S( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(29)
      INTEGER IX(29),JX(29)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-12,-10,-8,-6,-5,-5,-4,-4,-3,-3,-2,-1,-1,-1,
     &  0,0,0,0,1,1,3,3,3,4,4,4,5,14 /
      DATA JX / 20,24,22,14,36,8,16,6,32,3,8,4,1,2,3,0,1,4,28,0,
     &  32,0,1,2,3,18,24,4,24 /
      DATA CNX / -5.32466612140254D+22,1.00415480000824D+31,
     &  -1.91540001821367D+29,1.05618377808847D+16,
     &  2.02281884477061D+58,8.84585472596134D+07,
     &  1.66540181638363D+22,-3.13563197669111D+05,
     &  -1.85662327545324D+53,-6.24942093918942D-02,
     &  -5.04160724132590D+09,1.87514491833092D+04,
     &  1.21399979993217D-03,1.88317043049455D+00,
     &  -1.67073503962060D+03,9.65961650599775D-01,
     &  2.94885696802488D+00,-6.53915627346115D+04,
     &  6.04012200163444D+49,-1.98339358557937D-01,
     &  -1.75984090163501D+57,3.56314881403987D+00,
     &  -5.75991255144384D+02,4.56213415338071D+04,
     &  -1.09174044987829D+07,4.37796099975134D+33,
     &  -6.16552611135792D+45,1.93568768917797D+09,
     &  9.50898170425042D+53 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3S'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 21.D+6
      TRX = 640.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.886D+0
      BX = 0.990D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,29
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
      SPC_VOLX = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0022D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3S group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3T( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(33)
      INTEGER IX(33),JX(33)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / 0,0,0,0,1,1,2,2,2,3,3,4,4,7,7,7,7,7,10,10,10,10,10,
     &  18,20,22,22,24,28,32,32,32,36 /
      DATA JX / 0,1,4,12,0,10,0,6,14,3,8,0,10,3,4,7,20,36,10,12,14,
     &  16,22,18,32,22,36,24,28,22,32,36,36 /
      DATA CNX / 1.55287249586268D+00,6.64235115009031D+00,
     &  -2.89366236727210D+03,-3.85923202309848D+12,
     &  -2.91002915783761D+00,-8.29088246858083D+11,
     &  1.76814899675218D+00,-5.34686695713469D+08,
     &  1.60464608687834D+17,1.96435366560186D+05,
     &  1.56637427541729D+12,-1.78154560260006D+00,
     &  -2.29746237623692D+15,3.85659001648006D+07,
     &  1.10554446790543D+09,-6.77073830687349D+13,
     &  -3.27910592086523D+30,-3.41552040860644D+50,
     &  -5.27251339709047D+20,2.45375640937055D+23,
     &  -1.68776617209269D+26,3.58958955867578D+28,
     &  -6.56475280339411D+35,3.55286045512301D+38,
     &  5.69021454413270D+57,-7.00584546433113D+47,
     &  -7.05772623326374D+64,1.66861176200148D+52,
     &  -3.00475129680486D+60,-6.68481295196808D+50,
     &  4.28432338620678D+68,-4.44227367758304D+71,
     &  -2.81396013562745D+76 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3T'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 20.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      AX = 0.803D+0
      BX = 1.02D+0
      SPC_VOLX = 0.D+0
      DO 100 I = 1,33
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-AX)**IX(I))*
     &    ((THETAX-BX)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0088D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3T group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3U( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(38)
      INTEGER IX(38),JX(38)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-10,-10,-10,-8,-8,-8,-6,-6,-5,-5,-5,-3,-1,-1,
     &  -1,-1,0,0,1,2,2,3,5,5,5,6,6,8,8,10,12,12,12,14,14,14,14 /
      DATA JX / 14,10,12,14,10,12,14,8,12,4,8,12,2,-1,1,12,14,-3,
     &  1,-2,5,10,-5,-4,2,3,-5,2,-8,8,-4,-12,-4,4,-12,-10,-6,6 /
      DATA CNX / 1.22088349258355D+17,1.04216468608488D+09,
     &  -8.82666931564652D+15,2.59929510849499D+19,
     &  2.22612779142211D+14,-8.78473585050085D+17,
     &  -3.14432577551552D+21,-2.16934916996285D+12,
     &  1.59079648196849D+20,-3.39567617303423D+02,
     &  8.84387651337836D+12,-8.43405926846418D+20,
     &  1.14178193518022D+01,-1.22708229235641D-04,
     &  -1.06201671767107D+02,9.03443213959313D+24,
     &  -6.93996270370852D+27,6.48916718965575D-09,
     &  7.18957567127851D+03,1.05581745346187D-03,
     &  -6.51903203602581D+14,-1.60116813274676D+24,
     &  -5.10254294237837D-09,-1.52355388953402D-01,
     &  6.77143292290144D+11,2.76378438378930D+14,
     &  1.16862983141686D-02,-3.01426947980171D+13,
     &  1.69719813884840D-08,1.04674840020929D+26,
     &  -1.08016904560140D+04,-9.90623601934295D-13,
     &  5.36116483602738D+06,2.26145963747881D+21,
     &  -4.88731565776210D-10,1.51001548880670D-05,
     &  -2.27700464643920D+04,-7.81754507698846D+27 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3U'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 23.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      SPC_VOLX = 0.D+0
      DO 100 I = 1,38
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-0.902D+0)**IX(I))*
     &    ((THETAX-0.988D+0)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0026D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3U group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3V( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(39)
      INTEGER IX(39),JX(39)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -10,-8,-6,-6,-6,-6,-6,-6,-5,-5,-5,-5,-5,-5,-4,-4,
     &  -4,-4,-3,-3,-3,-2,-2,-1,-1,0,0,0,1,1,3,4,4,4,5,8,10,12,14 /
      DATA JX / -8,-12,-12,-3,5,6,8,10,1,2,6,8,10,14,-12,-10,-6,10,
     &  -3,10,12,2,4,-2,0,-2,6,10,-12,-10,3,-6,3,10,2,-12,-2,-3,1 /
      DATA CNX / -4.15652812061591D-55,1.77441742924043D-61,
     &  -3.57078668203377D-55,3.59252213604114D-26,
     &  -2.59123736380269D+01,5.94619766193460D+04,
     &  -6.24184007103158D+10,3.13080299915944D+16,
     &  1.05006446192036D-09,-1.92824336984852D-06,
     &  6.54144373749937D+05,5.13117462865044D+12,
     &  -6.97595750347391D+18,-1.03977184454767D+28,
     &  1.19563135540666D-48,-4.36677034051655D-42,
     &  9.26990036530639D-30,5.87793105620748D+20,
     &  2.80375725094731D-18,-1.92359972440634D+22,
     &  7.42705723302738D+26,-5.17429682450605D+01,
     &  8.20612048645469D+06,-1.88214882341448D-09,
     &  1.84587261114837D-02,-1.35830407782663D-06,
     &  -7.23681885626348D+16,-2.23449194054124D+26,
     &  -1.11526741826431D-35,2.76032601145151D-29,
     &  1.34856491567853D+14,6.52440293345860D-10,
     &  5.10655119774360D+16,-4.68138358908732D+31,
     &  -7.60667491183279D+15,-4.17247986986821D-19,
     &  3.12545677756104D+13,-1.00375333864186D+14,
     &  2.47761392329058D+26 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3V'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 23.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      SPC_VOLX = 0.D+0
      DO 100 I = 1,39
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-0.960D+0)**IX(I))*
     &    ((THETAX-0.995D+0)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0031D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3V group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3W( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(35)
      INTEGER IX(35),JX(35)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -12,-12,-10,-10,-8,-8,-8,-6,-6,-6,-6,-5,-4,-4,-3,
     &  -3,-2,-2,-1,-1,-1,0,0,1,2,2,3,3,5,5,5,8,8,10,10 /
      DATA JX / 8,14,-1,8,6,8,14,-4,-3,2,8,-10,-1,3,-10,3,1,2,-8,
     &  -4,1,-12,1,-1,-1,2,-12,-5,-10,-8,-6,-12,-10,-12,-8 /
      DATA CNX / -5.86219133817016D-08,-8.94460355005526D+10,
     &  5.31168037519774D-31,1.09892402329239D-01,
     &  -5.75368389425212D-02,2.28276853990249D+04,
     &  -1.58548609655002D+18,3.29865748576503D-28,
     &  -6.34987981190669D-25,6.15762068640611D-09,
     &  -9.61109240985747D+07,-4.06274286652625D-45,
     &  -4.71103725498077D-13,7.25937724828145D-01,
     &  1.87768525763682D-39,-1.03308436323771D+03,
     &  -6.62552816342168D-02,5.79514041765710D+02,
     &  2.37416732616644D-27,2.71700235739893D-15,
     &  -9.07886213483600D+01,-1.71242509570207D-37,
     &  1.56792067854621D+02,9.23261357901470D-01,
     &  -5.97865988422577D+00,3.21988767636389D+06,
     &  -3.99441390042203D-30,4.93429086046981D-08,
     &  8.12036983370565D-20,-2.07610284654137D-12,
     &  -3.40821291419719D-07,5.42000573372233D-18,
     &  -8.56711586510214D-13,2.66170454405981D-14,
     &  8.58133791857099D-06 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3W'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 23.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      SPC_VOLX = 0.D+0
      DO 100 I = 1,35
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-0.959D+0)**IX(I))*
     &    ((THETAX-0.995D+0)**JX(I))
  100 CONTINUE
      SPC_VOLX = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0039D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3X( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(36)
      INTEGER IX(36),JX(36)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -8,-6,-5,-4,-4,-4,-3,-3,-1,0,0,0,1,1,2,3,3,3,4,5,5,
     &  5,6,8,8,8,8,10,12,12,12,12,14,14,14,14 /
      DATA JX / 14,10,10,1,2,14,-2,12,5,0,4,10,-10,-1,6,-12,0,8,3,
     &  -6,-2,1,1,-6,-3,1,8,-8,-10,-8,-5,-4,-12,-10,-8,-6 /
      DATA CNX / 3.77373741298151D+18,-5.07100883722913D+12,
     &  -1.03363225598860D+15,1.84790814320773D-06,
     &  -9.24729378390945D-04,-4.25999562292738D+23,
     &  -4.62307771873973D-13,1.07319065855767D+21,
     &  6.48662492280682D+10,2.44200600688281D+00,
     &  -8.51535733484258D+09,1.69894481433592D+21,
     &  2.15780222509020D-27,-3.20850551367334D-01,
     &  -3.82642448458610D+16,-2.75386077674421D-29,
     &  -5.63199253391666D+05,-3.26068646279314D+20,
     &  3.97949001553184D+13,1.00824008584757D-07,
     &  1.62234569738433D+04,-4.32355225319745D+10,
     &  -5.92874245598610D+11,1.33061647281106D+00,
     &  1.57338197797544D+06,2.58189614270853D+13,
     &  2.62413209706358D+24,-9.20011937431142D-02,
     &  2.20213765905426D-03,-1.10433759109547D+01,
     &  8.47004870612087D+06,-5.92910695762536D+08,
     &  -1.83027173269660D-05,1.81339603516302D-01,
     &  -1.19228759669889D+03,4.30867658061468D+06 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3X'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 23.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      SPC_VOLX = 0.D+0
      DO 100 I = 1,36
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-0.910D+0)**IX(I))*
     &    ((THETAX-0.988D+0)**JX(I))
  100 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0049D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3X group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3Y( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(20)
      INTEGER IX(20),JX(20)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / 0,0,0,0,1,2,2,2,2,3,3,3,4,4,5,5,8,8,10,12 /
      DATA JX / -3,1,5,8,8,-4,-1,4,5,-8,4,8,-6,6,-2,1,-8,-2,-5,-8 /
      DATA CNX / -5.25597995024633D-10,5.83441305228407D+03,
     &  -1.34778968457925D+16,1.18973500934212D+25,
     &  -1.59096490904708D+26,-3.15839902302021D-07,
     &  4.96212197158239D+02,3.27777227273171D+18,
     &  -5.27114657850696D+21,2.10017506281863D-17,
     &  7.05106224399834D+20,-2.66713136106469D+30,
     &  -1.45370512554562D-08,1.49333917053130D+27,
     &  -1.49795620287641D+07,-3.81881906271100D+15,
     &  7.24660165585797D-05,-9.37808169550193D+13,
     &  5.14411468376383D+09,-8.28198594040141D+04 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3Y'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 22.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      SPC_VOLX = 0.D+0
      DO 100 I = 1,20
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-0.996D+0)**IX(I))*
     &    ((THETAX-0.994D+0)**JX(I))
  100 CONTINUE
      SPC_VOLX  = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0031D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3Y group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_3Z( TX,PX,RHOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density from the Revised Supplementary Release 
!     on Backward Equations for Specific Volume as a Function of 
!     Pressure and Temperature v(p,T) for Region 3a of the IAPWS 
!     Industrial Formulation 1997 for the Thermodynamic 
!     Properties of Water and Steam, as a function of 
!     temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNX(23)
      INTEGER IX(23),JX(23)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNX,IX,JX
      DATA IX / -8,-6,-5,-5,-4,-4,-4,-3,-3,-3,-2,-1,0,1,2,3,3,
     &  6,6,6,6,8,8 /
      DATA JX / 3,6,6,8,5,6,8,-2,5,6,2,-6,3,1,6,-6,-2,-6,-5,-4,
     &  -1,-8,-4 /
      DATA CNX / 2.4400789229065D-11,-4.6305743033124D+06,
     &  7.2880327477771D+09,3.2777630285886D+15,
     &  -1.1059817011841D+09,-3.2389991572996D+12,
     &  9.2381400702325D+15,8.4225008041371D-13,
     &  6.6322143624551D+11,-1.6717018667214D+14,
     &  2.5374935870139D+03,-8.1973155961052D-21,
     &  3.2838058789066D+11,-6.2500479117154D+07,
     &  8.0319795746202D+20,-2.0439701133835D-11,
     &  -3.7839104705594D+03,9.7287654593862D-03,
     &  1.5435572168146D+01,-3.7396286292864D+03,
     &  -6.8285901137457D+10,-2.4848801561454D-04,
     &  3.9453604949707D+06 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_3Z'
!
!---  Reduced pressure, "pi" and reduced temperature, "theta"  ---
!
      PRX = 22.D+6
      TRX = 650.D+0
      TKX = TX + TABS
      PIX = PX/PRX
      THETAX = TKX/TRX
      SPC_VOLX = 0.D+0
      DO 100 I = 1,23
        SPC_VOLX  = SPC_VOLX + CNX(I)*((PIX-0.993D+0)**IX(I))*
     &    ((THETAX-0.994D+0)**JX(I))
  100 CONTINUE
      SPC_VOLX  = SPC_VOLX**4
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(SPC_VOLX*0.0038D+0)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_3Z group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_4( TX,PX,INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density, enthalpy, and internal energy for
!     Region 5 per the Revised Release on the IAPWS Industrial 
!     Formulation 1997 for the Thermodynamic Properties of Water and 
!     Steam (The revision only relates to the extension of Region 4 
!     to 50 MPa), as a function of temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!     Vapor pressure above ice determined from Murphy, D.M, and T. Koop
!     2005. Review of the vapour pressures of ice and supercooled water
!     for atmospheric applications. Q.J.R. Meteorol. Soc. 131:1539-1565
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 February 2015.
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
      REAL*8 CNX(10)
!      REAL*8 SNX(3)
!
!----------------------Data Statements---------------------------------!
!
      DATA CNX / 0.11670521452767D+4, -0.72421316703206D+6, 
     &  -0.17073846940092D+2, 0.12020824702470D+5,
     &  -0.32325550322333D+7, 0.14915108613530D+2,
     &  -0.48232657361591D+4, 0.40511340542057D+6,
     &  -0.23855557567849D+0, 0.65017534844798D+3 /
!      DATA SNX / 0.53822143709235D+2, -0.356303356187D+0,
!     &  0.474775625854D-3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_4'
!
!---  INDX = 0, return pressure  ---
!
      IF( INDX.EQ.0 ) THEN
        TKX = TX+TABS
!
!---    Vapor pressure over ice  ---
!
        IF( TKX.LT.TABS ) THEN
          PX = EXP(9.550426D+0 - 5723.265D+0/TKX + 3.53068D+0*LOG(TKX)
     &      - 0.00728332*TKX)
!
!---    Vapor pressure over water  ---
!
        ELSEIF( TKX.LE.TCRW ) THEN
          THETAX = TKX + CNX(9)/(TKX-CNX(10))
          AX = THETAX**2 + CNX(1)*THETAX + CNX(2)
          BX = CNX(3)*THETAX**2 + CNX(4)*THETAX + CNX(5)
          CX = CNX(6)*THETAX**2 + CNX(7)*THETAX + CNX(8)
          PRX = 1.D+6
          PX = PRX*(2.D+0*CX/(-BX + SQRT(BX**2 - 4.D+0*AX*CX)))**4
        ELSE
!          PX = SNX(1) + SNX(2)*TKX + SNX(3)*(TKX**2)
          PX = PCRW
        ENDIF
!
!---  INDX != 0, return temperature  ---
!
      ELSE
        PRX = 1.D+6
!
!---    Vapor pressure over ice  ---
!
        IF( PX.LT.611.657D+0 ) THEN
          TKX = 114.539D+0 + 98.114D+0*(PX**0.074871D+0)
          TX = TKX - TABS
!
!---    Vapor pressure over water  ---
!
        ELSEIF( PX.LE.PCRW ) THEN
          BETAX = (PX/PRX)**2.5D-1
          EX = BETAX**2 + CNX(3)*BETAX + CNX(6)
          FX = CNX(1)*(BETAX**2) + CNX(4)*BETAX + CNX(7)
          GX = CNX(2)*(BETAX**2) + CNX(5)*BETAX + CNX(8)
          DX = 2.D+0*GX/(-FX - SQRT((FX**2)-(4.D+0*EX*GX)))
          TKX = (CNX(10) + DX - SQRT((CNX(10)+DX)**2 - 
     &      4.D+0*(CNX(9)+CNX(10)*DX)))/2.D+0
          TX = TKX - TABS
        ELSE
          TX = TCRW - TABS
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REGION_5( TX,PX,RHOX,HX,UX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density, enthalpy, and internal energy for
!     Region 5 per the Revised Release on the IAPWS Industrial 
!     Formulation 1997 for the Thermodynamic Properties of Water and 
!     Steam (The revision only relates to the extension of Region 5 
!     to 50 MPa), as a function of temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 February 2015.
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
      REAL*8 CNOX(6),CNRX(6)
      INTEGER JOX(6),IRX(6),JRX(6)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CNOX,CNRX
      SAVE JOX,IRX,JRX
      DATA CNOX / -0.13179983674201D+2,0.68540841634434D+1,
     &  -0.24805148933466D-1,0.36901534980333D+0,-0.31161318213925D+1,
     &  -0.32961626538917D+0/
      DATA CNRX / 0.15736404855259D-2,0.90153761673944D-3,
     &  -0.50270077677648D-2,0.22440037409485D-5,-0.41163275453471D-5,
     &  0.37919454822955D-7/
      DATA JOX / 0,1,-3,-2,-1,2 /
      DATA IRX / 1,1,1,2,2,3 /
      DATA JRX / 1,2,3,3,9,7 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REGION_5'
!
!---  Reduced pressure, "pi" and reduced temperature, "tau"  ---
!
      RCX = 0.461526D+3
      PRX = 1.D+6
      TRX = 1.D+3
      TKX = TX + TABS
      PIX = PX/PRX
      TAUX = TRX/TKX
!
!---  Partial derivative of the ideal-gas part of the dimensionless 
!     Gibbs free energy with respect to the reduced pressure, "pi" 
!     and temperature, "tau"  ---
!
      GAMMAO_PIX = 1.D+0/PIX
      GAMMAO_TAUX = 0.D+0
      DO 100 I = 1,6
        GAMMAO_TAUX = GAMMAO_TAUX + CNOX(I)*REAL(JOX(I))*
     &    (TAUX**(JOX(I)-1))
  100 CONTINUE
!
!---  Partial derivative of the residual part of the dimensionless 
!     Gibbs free energy with respect to the reduced pressure, "pi" 
!     and temperature, "tau"  ---
!
      GAMMAR_PIX = 0.D+0
      GAMMAR_TAUX = 0.D+0
      DO 200 I = 1,6
        GAMMAR_PIX = GAMMAR_PIX + CNRX(I)*REAL(IRX(I))*
     &    (PIX**(IRX(I)-1))*(TAUX**JRX(I))
        GAMMAR_TAUX = GAMMAR_TAUX + CNRX(I)*REAL(JRX(I))*
     &    (PIX**IRX(I))*(TAUX**(JRX(I)-1))
  200 CONTINUE
!
!---  Density, kg/m^3  ---
!
      RHOX = 1.D+0/(PIX*(GAMMAO_PIX + GAMMAR_PIX)*RCX*TKX/PX)
!
!---  Enthalpy, J/kg  ---
!
      HX = RCX*TKX*TAUX*(GAMMAO_TAUX + GAMMAR_TAUX)
!
!---  Internal energy, J/kg  ---
!
      UX = RCX*TKX*(TAUX*(GAMMAO_TAUX + GAMMAR_TAUX) - 
     &  PIX*(GAMMAO_PIX + GAMMAR_PIX))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REGION_5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SCF_GL( BTGLX,PWX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     This subroutine calculates the capillary head scaling factor
!     based on the surface tension of water.
!
!     Vargaftik, N.B., B.N. Volkov, and L.D. Voljak. 1983 International
!     tables of surface tension of water. J. Phys. Chem. Ref. Data.
!     12(3):817-820.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 April 2015.
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
      SUB_LOG(ISUB_LOG) = '/SCF_GL'
      IF( ISLC(7).EQ.0 ) THEN
        BTGLX = 1.D+0
      ELSE
        IF( PWX.GE.PCRW ) THEN
          BTGLX = 0.D+0
        ELSE
          INDX = 1
          CALL REGION_4( TSWX,PWX,INDX)
          TKSWX = MIN( TSWX+TABS,TCRW )
          TRX = (TCRW-TKSWX)/TCRW
          SIGMAX = 235.8D-3*(TRX**1.256D+0)*(1.D+0 - 0.625D+0*TRX)
          BTGLX = SIGMAX/72.75D-3
          BTGLX = MAX( BTGLX,0.D+0 )
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SCF_GL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SOL_BRNS( TX,PX,XLSMX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Solubility of NaCl in liquid water and 
!     supercritical water, for the following phase conditions:
!
!     Supercritical region:
!     Leusbrock, I., S.J. Metz, G. Rexwinkel, and G.F. Versteeg. 2008.
!     Quantitative approaches for the description of solubilities of
!     inorganic compounds in near-critical and supercritical water.
!     The Journal of Supercritical Fluids, 47:117-127.
!
!      Liquid water region:
!      Sawamura, S., N. Egoshi, Y. Setoguchi, and H. Matsuo. 2007. 
!      Solubility of sodium chloride in water under high pressure.
!      Fluid Phase Equilibria, 254:158-162.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 February 2015.
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
      REAL*8 CSPX(4),CSBX(6)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CSPX,CSBX
      SAVE RHOMOX,RHOMCX
      DATA CSPX / -59.93D+0,-2513.7D+0,6.23D+0,4.01D+0 /
      DATA CSBX / 2.86623D+0,-571.361D+0,7.7128D+4,4.0479D-4,
     &  -7.445D-7,6.209D-10 /
!
!---  RHOMOX - ln of the reference molar density
!     RHOMCX - ln of the critical molar density  --- 
!
      DATA RHOMOX / 4.01655053478208D+0 /
      DATA RHOMCX / 2.88334680134435D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SOL_BRNS'
!
!---  Absolute and reduced temperature and pressure  ---
!
      TKX = TX + TABS
      PMX = 1.D-6*PX
      CALL P_IAPWS( TX,PX,RHOGX,RHOLX,HGX,HLX,UGX,ULX )
      RHOMX = LOG(RHOLX/WTMW)
      RHOMRX = MAX(MIN(((RHOMX-RHOMOX)/(RHOMCX-RHOMOX)),1.D+0),0.D+0)
!
!---  High-temperature/low-pressure molality, mol/kg water  ---
!
      BWVX = EXP( CSBX(1) + CSBX(2)/TKX + CSBX(3)/(TKX**2) + 
     &  CSBX(4)*PMX + CSBX(5)*(PMX**2) + CSBX(6)*(PMX**3) )
!
!---  Supercritical molality, mol/kg water  ---
!
      IF( TKX.GT.TCRW ) THEN
        BSCX = EXP( CSPX(1) - CSPX(2)/TKX + CSPX(3)*LOG(TKX) +
     &    CSPX(4)*LOG(RHOMX) )
      ELSE
        BSCX = EXP( CSPX(1) - CSPX(2)/TCRW + CSPX(3)*LOG(TCRW) +
     &    CSPX(4)*LOG(RHOMCX) )
      ENDIF
!
!---  Combination molality, mol/kg water  ---
!
      BX = (1.D+0-RHOMRX)*BWVX + RHOMRX*BSCX
!
!---  Mass fraction of dissolved salt in brine  ---
!
      XLSMX = 1.D+0/(1.D+0 + 1.D+0/(BX*WTMS*1.D-3))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SOL_BRNS group  ---
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
      INDX = 0
      CALL REGION_4( TWX,PSBX,INDX )
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
      IF( RHOBX.LT.EPSL ) THEN
        TKX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
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
      SUB_LOG(ISUB_LOG) = '/VISC_G'
      IF( ISLC(9).EQ.1 ) THEN
        VISGX = VISGI
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
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
      IF( ISLC(9).EQ.1 ) THEN
        VISLX = VISLI
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
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
     & RHOBX.GE.0.755 .AND. RHOBX.LE.1.290 ) THEN
       DPX = 1.D-1
       DBETAX = DPX/PREF
       PIX = PX+DPX
       CALL P_IAPWS( TX,PIX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
       IF( (1.D+0-ABS(RHOLWX/RHOWX)).LT.(1.D+0-ABS(RHOGWX/RHOWX)) ) THEN
         RHOBIX = RHOLWX/RHOREF
       ELSE
         RHOBIX = RHOGWX/RHOREF
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
      SUBROUTINE V_IAPWS
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate water density, enthalpy, and internal energy for
!     Region 3 per the Revised Release on the IAPWS Industrial 
!     Formulation 1997 for the Thermodynamic Properties of Water and 
!     Steam, as a function of temperature and pressure.
!
!     Revised Release on the IAPWS Industrial Formulation 1997 for 
!     the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/V_IAPWS'
      IERR = 0
!
!---  Region 1  ---
!
      INDX = 1
      PX = 100.D+6
      CALL REGION_4( TX,PX,INDX )
!
!---  Region 1  ---
!
      TX = 300.D+0-TABS
      PX = 3.D+6
      CALL REGION_1( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.00215168D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 300.D+0-TABS
      PX = 80.D+6
      CALL REGION_1( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 0.971180894D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 500.D+0-TABS
      PX = 3.D+6
      CALL REGION_1( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.20241800D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 2  ---
!
      TX = 300.D+0-TABS
      PX = 0.0035D+6
      CALL REGION_2( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 0.394913866D+2
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 700.D+0-TABS
      PX = 0.0035D+6
      CALL REGION_2( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 0.923015898D+2
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 700.D+0-TABS
      PX = 30.D+6
      CALL REGION_2( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 0.542946619D-2
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3  ---
!
      TX = 650.D+0-TABS
      PX = 25.5837018D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.D+0/500.D+0
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 650.D+0-TABS
      PX = 22.2930643D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.D+0/200.D+0
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 750.D+0-TABS
      PX = 78.3095639D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.D+0/500.D+0
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3a  ---
!
      TX = 630.D+0-TABS
      PX = 50.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.470853100D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 670.D+0-TABS
      PX = 80.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.503831359D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3b  ---
!
      TX = 710.D+0-TABS
      PX = 50.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.204728587D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 750.D+0-TABS
      PX = 80.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.973692940D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3c  ---
!
      TX = 630.D+0-TABS
      PX = 20.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.761696406D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 650.D+0-TABS
      PX = 30.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.819560617D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3d  ---
!
      TX = 656.D+0-TABS
      PX = 26.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.245587720D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 670.D+0-TABS
      PX = 30.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.506897702D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3e  ---
!
      TX = 661.D+0-TABS
      PX = 26.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.970225962D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 675.D+0-TABS
      PX = 30.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 3.004627086D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3f  ---
!
      TX = 671.D+0-TABS
      PX = 26.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 5.019029401D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 690.D+0-TABS
      PX = 30.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 4.656470142D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3g  ---
!
      TX = 649.D+0-TABS
      PX = 23.6D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.163198378D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 650.D+0-TABS
      PX = 24.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.166044161D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3h  ---
!
      TX = 652.D+0-TABS
      PX = 23.6D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.651081407D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 654.D+0-TABS
      PX = 24.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.967802335D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3i  ---
!
      TX = 653.D+0-TABS
      PX = 23.6D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      TX = 655.D+0-TABS
      PX = 24.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
!
!---  Region 3j  ---
!
      TX = 655.D+0-TABS
      PX = 23.5D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 4.545001142D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 660.D+0-TABS
      PX = 24.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 5.100267704D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3k  ---
!
      TX = 660.D+0-TABS
      PX = 23.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 6.109525997D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 670.D+0-TABS
      PX = 24.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 6.427325645D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3l  ---
!
      TX = 646.D+0-TABS
      PX = 22.6D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.117860851D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 646.D+0-TABS
      PX = 23.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.062374674D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3m  ---
!
      TX = 648.6D+0-TABS
      PX = 22.6D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.533063780D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 649.3D+0-TABS
      PX = 22.8D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.572971781D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3n  ---
!
      TX = 649.D+0-TABS
      PX = 22.6D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      TX = 649.7D+0-TABS
      PX = 22.8D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
!
!---  Region 3o  ---
!
      TX = 649.1D+0-TABS
      PX = 22.6D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 3.131208996D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 649.9D+0-TABS
      PX = 22.8D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 3.221160278D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3p  ---
!
      TX = 649.4D+0-TABS
      PX = 22.6D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 3.715596186D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 650.2D+0-TABS
      PX = 22.8D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 3.221160278D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3q  ---
!
      TX = 640.D+0-TABS
      PX = 21.1D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.970999272D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 643.D+0-TABS
      PX = 21.8D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.043919161D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3r  ---
!
      TX = 644.D+0-TABS
      PX = 21.1D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 5.251009921D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 648.D+0-TABS
      PX = 21.8D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 5.256844741D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3s  ---
!
      TX = 635.D+0-TABS
      PX = 19.1D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.932829079D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 638.D+0-TABS
      PX = 20.0D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 1.985387227D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3t  ---
!
      TX = 626.D+0-TABS
      PX = 17.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 8.483262001D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 640.D+0-TABS
      PX = 20.0D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 6.227528101D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3u  ---
!
      TX = 644.6D+0-TABS
      PX = 21.5D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.268366647D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 646.1D+0-TABS
      PX = 22.D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.2963505530D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3v  ---
!
      TX = 648.6D+0-TABS
      PX = 22.5D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.832373260D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 647.9D+0-TABS
      PX = 22.3D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.811424405D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3w  ---
!
      TX = 647.5D+0-TABS
      PX = 22.15D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 3.694032281D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 648.1D+0-TABS
      PX = 22.3D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 3.622226305D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3x  ---
!
      TX = 648.0D+0-TABS
      PX = 22.11D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 4.528072649D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 649.0D+0-TABS
      PX = 22.3D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 4.556905799D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3y  ---
!
      TX = 646.84D+0-TABS
      PX = 22.0D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.698354719D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 647.05D+0-TABS
      PX = 22.064D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 2.717655648D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 3z  ---
!
      TX = 646.89D+0-TABS
      PX = 22.0D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 3.798732962D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 647.15D+0-TABS
      PX = 22.064D+6
      CALL REGION_3( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 3.701940010D-3
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Region 4  ---
!
      INDX = 0
      TX = 300.D+0-TABS
      CALL REGION_4( TX,PX,INDX )
      PRX = 0.353658941D+4
      IF( ABS(PX-PRX)/PRX.GT.1.D-6 ) IERR = 1
      TX = 500.D+0-TABS
      CALL REGION_4( TX,PX,INDX )
      PRX = 0.263889776D+7
      IF( ABS(PX-PRX)/PRX.GT.1.D-6 ) IERR = 1
      TX = 600.D+0-TABS
      CALL REGION_4( TX,PX,INDX )
      PRX = 0.123443146D+8
      IF( ABS(PX-PRX)/PRX.GT.1.D-6 ) IERR = 1
!
!---  Region 5  ---
!
      TX = 1500.D+0-TABS
      PX = 0.5D+6
      CALL REGION_5( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 0.138455090D+1
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 1500.D+0-TABS
      PX = 30.D+6
      CALL REGION_5( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 0.230761299D-1
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
      TX = 2000.D+0-TABS
      PX = 30.D+6
      CALL REGION_5( TX,PX,RHOX,HX,UX )
      VX = 1.D+0/RHOX
      VRX = 0.311385219D-1
      IF( ABS(VX-VRX)/VRX.GT.1.D-6 ) IERR = 1
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of V_IAPWS group  ---
!
      RETURN
      END

!!----------------------Subroutine--------------------------------------!
!!
!      SUBROUTINE VPL_B( TX,PSBX,PGX,PLX,RHOBX,PVBX )
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
!!     Vapor pressure lowering of brine.
!!
!!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  The simulator
!!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!!     noncondensible gas.  Geothermics, 26(4): 437-464.
!!
!!
!!----------------------Authors-----------------------------------------!
!!
!!     Written by MD White, PNNL, 25 April 2002.
!!
!!----------------------Fortran 90 Modules------------------------------!
!!
!      USE GLB_PAR
!      USE SOLTN
!      USE CONST
!!
!!----------------------Implicit Double Precision-----------------------!
!!
!      IMPLICIT REAL*8 (A-H,O-Z)
!      IMPLICIT INTEGER (I-N)
!!
!!----------------------Executable Lines--------------------------------!
!!
!      ISUB_LOG = ISUB_LOG+1
!      SUB_LOG(ISUB_LOG) = '/VPL_B'
!!
!!---  Kelvin's equation  ---
!!
!      TKX = TX + TABS
!      PCX = PGX - PLX
!!
!!---  Positive capillary pressure  ---
!!
!      IF( PCX/EPSL.GT.EPSL ) THEN
!!
!!---    Minimum reduced vapor pressure  ---
!!
!        IF( PGX.LE.PSBX ) THEN
!          VARX = PSBX*WTMW*EXP( WTMW*PLX/(RHOBX*RCU*TKX) )/
!     &      (RHOBX*RCU*TKX)
!          PVMNX = 0.D+0
!          FCX = 1.D+0
!          DO 10 N = 1,7
!            RNX = REAL(N)
!            FCX = FCX*RNX
!            PVMNX = PVMNX + ((-RNX)**(N-1))/FCX*(VARX**N)
!   10     CONTINUE
!          PVMNX = PVMNX*RHOBX*RCU*TKX/WTMW
!          PVBX = PSBX*EXP( -WTMW*PCX/(RHOBX*RCU*TKX) )
!          PVBX = MAX( PVMNX,PVBX )
!        ELSE
!          PVBX = PSBX*EXP( -WTMW*PCX/(RHOBX*RCU*TKX) )
!        ENDIF
!!
!!---  Zero capillary pressure  ---
!!
!      ELSE
!        PVBX = PSBX
!      ENDIF
!!
!!---  Reset subroutine string sequence  ---
!!
!      ISUB_LOG = ISUB_LOG-1
!!
!!---  End of VPL_B group  ---
!!
!      RETURN
!      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VPL_B( TX,PSBX,PCX,RHOBX,PVBX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
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
!
!---  Positive capillary pressure  ---
!
      IF( PCX/EPSL.GT.EPSL ) THEN
        PVBX = PSBX*EXP( -WTMW*PCX/(RHOBX*RCU*TKX) )
!
!---  Zero capillary pressure  ---
!
      ELSE
        PVBX = PSBX
      ENDIF
      PVBX = MAX( 1.D+0,PVBX )
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
!     Written by MD White, PNNL, 16 July 2010
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
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WEBB_BC'
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
        IMSGX = N_DB
        NMSGX = 0
        SUBLOGX = 'WEBB_BC'
        RLMSGX = 0.D+0
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
!     Written by MD White, PNNL, 16 July 2010
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
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WEBB_VG'
!!
!!---  Initial guess of the matching-point head  ---
!!
!      SMPX = 1.D-2*(1.D+0-SRX) + SRX
!      ESMPX = (SMPX-SRX)/(1.D+0-SRX)
!      HMPX = (1.D+0/ALPHAX)*
!     &  (((1.D+0/(ESMPX**(1.D+0/CMX)))-1.D+0)**(1.D+0/CNX))
!      HMPX = 1.D+1**(4.D-1*LOG10(HDOD))
!!
!!---  Newton-Raphson iteration for the matching-point head  ---
!!
!      NC = 0
!  100 CONTINUE
!      NC = NC + 1
!      DHX = MAX( 1.D-4,1.D-6*HMPX )
!      DHX = SIGN( DHX,5.D-1*HDOD-HMPX )
!      AX = ((1.D+0 + (ALPHAX*HMPX)**CNX)**(-CMX))*(1.D+0-SRX) + SRX
!      BX = HMPX*(LOG(HDOD)-LOG(HMPX))
!      CX = (1.D+0 + (ALPHAX*HMPX)**CNX)**(-CMX-1.D+0)
!      DX = CNX*CMX*((ALPHAX*HMPX)**CNX)*(1.D+0-SRX)/HMPX
!      F1X = AX/BX - CX*DX
!      HMPY = HMPX + DHX
!      AX = ((1.D+0 + (ALPHAX*HMPY)**CNX)**(-CMX))*(1.D+0-SRX) + SRX
!      BX = HMPY*(LOG(HDOD)-LOG(HMPY))
!      CX = (1.D+0 + (ALPHAX*HMPY)**CNX)**(-CMX-1.D+0)
!      DX = CNX*CMX*((ALPHAX*HMPY)**CNX)*(1.D+0-SRX)/HMPY
!      F2X = AX/BX - CX*DX
!      DFX = (F2X-F1X)/DHX
!      DHMPX = -F1X/DFX
!      HMPX = MIN( MAX( HMPX+DHMPX,1.D-6 ),HDOD-1.D+0 )
!!
!!---  No convergence on saturation matching point  ---
!!
!      IF( NC.GT.64 ) THEN
!        INDX = 7
!        IMSG = N_DB
!        CHMSG = 'No Convergence on van Genuchten'
!     &    // ' Matching Point Head @ Node: '
!        CALL WRMSGX( INDX )
!      ENDIF
!      IF( ABS(DHMPX).GT.1.D-4 ) GOTO 100
!!
!!---  Find the capillary matching point saturation  ---
!!
!      ESMPX = (1.D+0 + (ALPHAX*HMPX)**CNX)**(-CMX)
!      SMPX = ESMPX*(1.D+0-SRX) + SRX
!
!---  Use the matrix saturation at 0.4 LOG10(HDOD) as
!     the initial guess  ---
!
      HDX = 1.D+1**(4.D-1*LOG10(HDOD))
      SMPX = (1.D+0/(1.D+0 + (ALPHAX*HDX)**CNX))**CMX
      SMPX = SMPX*(1.D+0-SRX) + SRX
      SMPX = 1.D-1*(1.D+0-SRX) + SRX
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
        INDX = 7
        IMSGX = N_DB
        NMSGX = 0
        SUBLOGX = 'WEBB_VG'
        RLMSGX = 0.D+0
        CHMSGX = 'No Convergence on Saturation '
     &    // 'Matching Point @ Node: '
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
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


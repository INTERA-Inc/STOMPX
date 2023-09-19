!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VOCLQH( TX,HX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate the liquid VOC enthalpy as a function of temperature
!     with a reference state of liquid VOC at 0 C.
!
!     Calculate the ideal gas isobaric specific heat from empirical
!     coefficients. pp. 657 (Reid et al.)
!
!     Calculate the specific heat for liquid VOC from the ideal gas
!     molar specfic heat using the corresponding states method.
!     pp. 140 (Reid et al.)
!
!     Calculate the liquid VOC enthalpy by integrating the specific heat
!     for liquid VOC from a reference temperature of 0 C, with a
!     two-point Gauss-Legendre quadrature method.
!     pp. 101-105 (Carnahan et al.)
!
!     Carnahan, B., H.A. Luther, and J.O. Wilkes. 1969. Applied
!     Numerical Methods. John Wiley and Sons Inc., New York.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York.
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < TX < VOC critical temperature
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     voclqh.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/VOCLQH'
      TK1 = (-0.57735D+0*TX + TX + 2.D+0*TABS)/2.D+0
      TK2 = (0.57735D+0*TX + TX + 2.D+0*TABS)/2.D+0
      TR1 = TK1/TCRO
      TR2 = TK2/TCRO
      CPO1 = CPOC(1) + CPOC(2)*TK1 + CPOC(3)*TK1**2 + CPOC(4)*TK1**3
      CPO2 = CPOC(1) + CPOC(2)*TK2 + CPOC(3)*TK2**2 + CPOC(4)*TK2**3
      DCPL1 = 1.45D+0 + 4.5D-1/(1.D+0-TR1) +
     &  2.5D-1*PAFO*(1.711D+1 +
     &  2.52D+1*((1.D+0-TR1)**(3.333D-1))/TR1 + 1.742D+0/(1.D+0-TR1))
      DCPL2 = 1.45D+0 + 4.5D-1/(1.D+0-TR2) +
     &  2.5D-1*PAFO*(1.711D+1 +
     &  2.52D+1*((1.D+0-TR2)**(3.333D-1))/TR2 + 1.742D+0/(1.D+0-TR2))
      DCSL1 = EXP(2.01D+1*TR1-1.79D+1)
      DCSL2 = EXP(2.01D+1*TR2-1.79D+1)
      CPNR1 = RCO*(DCPL1-DCSL1) + 1.0D+3*CPO1/WTMO
      CPNR2 = RCO*(DCPL2-DCSL2) + 1.0D+3*CPO2/WTMO
      HX = TX*(CPNR1+CPNR2)/2.D+0
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VOCLQH group  ---
!
      RETURN
      END

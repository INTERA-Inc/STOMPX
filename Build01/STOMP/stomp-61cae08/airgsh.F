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
!     airgsh.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!     airgsh.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!     airgsh.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of AIRGSK group  ---
!
      RETURN
      END


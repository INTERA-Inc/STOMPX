
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
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     bndfow.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      TK = TX+TABS
      PB = PX*1.D-5
      TP = TK/SQRT(1.15D+0*TBO*1.15D+0*TBW)
      SIG = 5.D-1*(1.18D+0*(2.85D-1*(VCRO**1.048D+0))**3.333D-1
     &  + 1.18D+0*(2.85D-1*(VCRW**1.048D+0))**3.333D-1)
      OMG = (A(1)/TP**A(2)) + (A(3)/EXP(A(4)*TP))
     &  + (A(5)/EXP(A(6)*TP))  + (A(7)/EXP(A(8)*TP))
      WTMX = 2.D+0/((1.D+0/WTMO) + (1.D+0/WTMW))
      DIFX = (3.03D+0 - (9.8D-1/SQRT(WTMX)))*1.D-3*(TK**1.5D+0)/
     &  (PB*SQRT(WTMX)*(SIG**2)*(OMG))*1.D-4
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BNDFOW group  ---
!
      RETURN
      END

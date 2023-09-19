!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VOCGSH( TX,PX,RHOX,HLQX,HVPX,UVPX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate the enthalpy and internal energy of VOC vapor as a
!     function of temperature and pressure with a reference state of
!     liquid VOC at 0 C.
!
!     Calculate the enthalpy as the sum of the liquid enthalpy and heat
!     of vaporization.
!
!     Calculate the heat of vaporization at the normal boiling point
!     with a method by Chen. pp. 226
!
!     Calculate the heat of vaporization as a function of temperature
!     with the Watson relation and the Viswanath and Kuloor extention.
!     pp. 228
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     vocgsh.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 N
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VOCGSH'
      IF( PX .LT. SMALL ) THEN
        HVPX = 0.D+0
        UVPX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      TR = (TX+TABS)/TCRO
      TBR = TBO/TCRO
      HVPB = RCO*TBO*((3.978D+0*TBR) - 3.958D+0 +
     &  (1.555D+0*LOG(PCRO*1.D-5)))/(1.07D+0-TBR)
      N = (2.64D-3*HVPB/(RCO*TBO) + 8.794D-1)**10
      HVPX = HVPB*(((1.D+0-TR)/(1.D+0-TBR))**N)
      HVPX = HLQX + HVPX
      UVPX = HVPX - PX/RHOX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VOCGSH group  ---
!
      RETURN
      END

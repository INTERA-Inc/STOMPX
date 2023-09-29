!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATGSV( TX,VISX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculates component viscosity of water vapor
!     with the corresponding states method.  pp 397.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     watgsv.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/WATGSV'
      TR = (TX+TABS)/TCRW
      DPMR = 5.246D+1*(DPMW**2)*(PCRW/1.D+5)/TCRW
      IF( DPMR .LT. 2.2D-1 ) THEN
        FPO = 1.D+0
      ELSEIF( DPMR .LT. 7.5D-2 ) THEN
        FPO = 1.D+0 + 3.055D+1*((2.92D-1-ZCRW)**1.72D+0)
      ELSE
        FPO = 1.D+0 + 3.055D+1*((2.92D-1-ZCRW)**1.72D+0)*
     &    ABS(9.6D-1 + 1.D-1*(TR-7.D-1))
      ENDIF
      EPSILON = (1.76D-1*((TCRW/((WTMW**3)*((PCRW/1.D+5)**4)))
     &  **(1.6667D-1)))
      VISX = (0.807D+0*(TR**0.618D+0) - 0.357D+0*EXP(-0.449D+0*TR)
     &  + 0.340D+0*EXP(-4.058D+0*TR) + 0.018D+0)*FPO/EPSILON
      VISX = VISX*1.D-7
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATGSV group  ---
!
      RETURN
      END

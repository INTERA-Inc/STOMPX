!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VOCLQV( TX,VISX,VISCX,IVISX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculate the saturated volatile organic compound liquid viscosity
!     as a function of temperature, with the Lewis and Squires liquid
!     viscosity-temperature correlation or with empirical
!     equations.  pp 441- 445.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < TX < VOC critical temperature
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     voclqv.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 VISCX(4)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VOCLQV'
      TK = TX + TABS
      IF( IVISX .EQ. 1 ) THEN
        VISX = ( VISCX(1)**(-2.661D-1) +
     &    (TX-VISCX(2))/2.33D+2 )**(-3.758D+0)
      ELSEIF( IVISX .EQ. 2 ) THEN
        VISX = VISCX(1)*(TK**VISCX(2))*1.D-3
      ELSEIF( IVISX .EQ. 3 ) THEN
        VISX = EXP( VISCX(1) + VISCX(2)/TK )*1.D-3
      ELSEIF( IVISX .EQ. 4 ) THEN
        VISX = EXP( VISCX(1) + VISCX(2)/TK + VISCX(3)*TK +
     &    VISCX(4)*(TK**2) )*1.D-3
      ELSEIF( IVISX .EQ. 5 ) THEN
        VISX = VISCX(1)
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VOCLQV group  ---
!
      RETURN
      END

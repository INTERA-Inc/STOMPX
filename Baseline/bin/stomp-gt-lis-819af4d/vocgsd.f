!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VOCGSD( TX,PX,RHOX,INDX )
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
!     INDX = 0
!
!     Calculates vapor component density of an volatile organic compound
!     with the ideal gas law equation of state.
!
!     INDX = 1
!
!     Calculates component density of an volatile organic compound
!     with van der Waals equation of state using an iterative
!     solution. pp. 43.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     vocgsd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/VOCGSD'
      IF( INDX .EQ. 0 ) THEN
        RHOX = PX/((TX+TABS)*RCO)
      ELSEIF( INDX .EQ. 1 ) THEN
        IF( PX .LT. SMALL ) THEN
          RHOX = 0.D+0
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
        TK = TX + TABS
        VO = RCO*TK/PX
  100   CONTINUE
        A = 2.7D+1*((RCO*TCRO)**2)/(6.4D+1*PCRO)
        B = RCO*TCRO/(8.D+0*PCRO)
        F = PX - RCO*TK/(VO-B) + A/(VO**2)
        DF = RCO*TK/((VO-B)**2) - 2.D+0*A/(VO**3)
        DVO = -F/DF
        VO = VO + DVO
        IF( ABS(DVO/(VO+SMALL)) .GE. 1.0D-8 ) GOTO 100
        RHOX = 1.D+0/VO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VOCGSD group  ---
!
      RETURN
      END


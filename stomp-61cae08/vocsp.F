!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VOCSP( TX,PX )
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
!     Calculates the saturation pressure of an organic vapor
!     as a function of temperature by an empirical equation. pp. 656.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids. Appendix A.
!     McGraw-Hill, New York, New York
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < TX < VOC critical temperature
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     vocsp.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/VOCSP'
      IF( IVAPO .EQ. 1 ) THEN
        TR1 = 1.D+0-(TX+TABS)/TCRO
        PX = PCRO*EXP((SATOC(1)*TR1 + SATOC(2)*TR1**1.5D+0
     &    + SATOC(3)*TR1**3 + SATOC(4)*TR1**6)/(1.D+0-TR1))
      ELSEIF( IVAPO .EQ. 2 ) THEN
        TK = TX + TABS
        PX = EXP( SATOC(1)-SATOC(2)/TK+SATOC(3)*LOG(TK) )
  100   CONTINUE
        F = LOG(PX) - SATOC(1) + SATOC(2)/TK - SATOC(3)*LOG(TK)
     &    - PX*SATOC(4)/(TK**2)
        DF = 1.D+0/PX - SATOC(4)/(TK**2)
        DPX = -F/DF
        PX = PX + DPX
        IF( ABS(DPX/(PX+SMALL)) .GE. 1.0D-8 ) GOTO 100
        PX = PX*1.D+5
      ELSEIF( IVAPO .EQ. 3 ) THEN
        TK = TX + TABS
        PX = EXP( SATOC(1) - SATOC(2)/(TK + SATOC(3)))*1.D+5
      ELSEIF( IVAPO .EQ. 4 ) THEN
        PX = SATOC(1)
      ELSEIF( IVAPO .EQ. 5 ) THEN
        TK = TX + TABS
        PX = EXP( SATOC(1) - SATOC(2)/(TK + SATOC(3)))
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VOCSP group  ---
!
      RETURN
      END

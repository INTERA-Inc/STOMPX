!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VOCLQD( TX,PX,PSOX,RHOX )
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
!     Calculate the saturated volatile organic compound liquid density
!     as a function of temperature with a modified Rackett technique.
!     Calculate the compressed volatile organic compound liquid density
!     as a function of temperature and pressure following the
!     Hankinson-Brobst-Thomson technique. pp 66-67.
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
!     voclqd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 A1(8),A2(8)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A1,A2
      DATA A1 /-9.070217D+0,6.245326D+1,-1.351102D+2,4.79594D+0,
     &        2.50047D-1,1.14188D+0,8.61488D-2,3.44483D-2/
      DATA A2 /-1.52816D+0,1.43907D+0,-8.1446D-1,1.90454D-1,
     &         -2.96123D-1,3.86914D-1,-4.27258D-2,-4.80645D-2/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VOCLQD'
      TR = (TX+TABS)/TCRO
      TR1 = 1.D+0 - TR
      TRR = (TDRO+TABS)/TCRO
      TRR1 = 1.D+0 - TRR
!
!---  Hankinson-Brobst-Thomson Technique  ---
!
      IF( IRHOO .EQ. 1 ) THEN
        VRO = 1.D+0 + A2(1)*TR1**(3.333D-1) +
     &    A2(2)*TR1**(6.666D-1) + A2(3)*TR1 +
     &    A2(4)*TR1**(1.333D+0)
        VRD = (A2(5) +A2(6)*TR +A2(7)*TR**2 +A2(8)*TR**3)/
     &    (TR-1.00001D+0)
        RHOSX = WTMO/(PCCVO*VRO*(1.D+0-WHBTO*VRD))
      ELSEIF( IRHOO .EQ. 2 ) THEN
        RHOSX = PCRO/(RCO*TCRO*ZRAO**(1.D+0 + TR1**(2.8571D-1)))
      ELSEIF( IRHOO .EQ. 3 ) THEN
        PHI = TR1**(2.8571D-1) - TRR1**(2.8571D-1)
        RHOSX = 1.D+0/((1.D+0/RHORO)*(ZRAO+SMALL)**PHI)
      ELSEIF( IRHOO .EQ. 4 ) THEN
        RHOX = RHORO
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      E = EXP( A1(4) + A1(5)*WHBTO + A1(6)*WHBTO**2 )
      CX = A1(7) + A1(8)*WHBTO
      B = PCRO*(-1.D+0 + A1(1)*TR1**(3.333D-1) +
     &  A1(2)*TR1**(6.666D-1) + A1(3)*TR1 +
     &  E*TR1**(1.333D+0))
      PX = MAX( PX,PSOX )
      RHOX = RHOSX*(1.D+0/(1.D+0 - CX*LOG((B+PX)/(B+PSOX))))
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VOCLQD group  ---
!
      RETURN
      END

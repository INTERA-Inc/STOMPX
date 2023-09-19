!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATGSH( TX,P,RHO,H,U )
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
!     Calculate the water vapor enthalpy and internal energy, as a
!     function of temperature and pressure per the steam table equations
!     as given by the 1967 International Formulation Committee:
!     Formulation for Industrial Use.
!
!     Thermodynamic and Transport Properties of Steam.
!     1967. ASME Steam Tables.
!     The American Society of Mechanical Engineers.
!     United Engineering Center, 345 East 47th Street, New York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < T > 364.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     watgsh.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 B(31),BB(5),L(3)
      INTEGER JN(8),JL(8),JX(8,2),JZ(8,3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE B,BB,L,SL1,JN,JL,JX,JZ
      DATA B /1.683599274D+1,2.856067796D+1,-5.438923329D+1,
     &4.330662834D-1,
     &-6.547711697D-1,8.565182058D-2,6.670375918D-2,1.388983801D+0,
     &8.390104328D-2,2.614670893D-2,-3.373439453D-2,4.520918904D-1,
     &1.069036614D-1,-5.975336707D-1,-8.847535804D-2,5.958051609D-1,
     &-5.159303373D-1,2.075021122D-1,1.190610271D-1,-9.867174132D-2,
     &1.683998803D-1,-5.809438001D-2,6.552390126D-3,5.710218649D-4,
     &1.936587558D+2,-1.388522425D+3,4.126607219D+3,-6.508211677D+3,
     &5.745984054D+3,-2.693088365D+3,5.235718623D+2/
      DATA BB /7.633333333D-1,4.006073948D-1,8.636081627D-2,
     &-8.532322921D-1,3.460208861D-1/
      DATA L/1.574373327D+1,-3.417061978D+1,1.931380707D+1/
      DATA SL1/4.260321148D+0/
      DATA JN/2,3,2,2,3,2,2,2/
      DATA JL/0,0,0,0,0,1,1,2/
      DATA JX/0,0,0,0,0,14,19,54,0,0,0,0,0,0,0,27/
      DATA JZ/13,18,18,25,32,12,24,24,3,2,10,14,28,11,18,
     &14,0,1,0,0,24,0,0,0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATGSH'
      TR = (TX+TABS)/TCRW
      PR = MAX( P/PCRW,1.D-12 )
      CX = EXP(BB(1)*(1.D+0-TR))
      BL = L(1) + L(2)*TR + L(3)*TR*TR
      DBL = L(2) + 2.D+0*L(3)*TR
      DDBL = 2.D+0*L(3)
      HX = B(1)*TR
      DO 110 I = 1,5
      HX = HX - B(I+1)*(I-2)*(TR**(I-1))
  110 CONTINUE
      IB = 6
      DO 130 I = 1,5
        DO 120 J = 1,JN(I)
          IB = IB + 1
          HX = HX-(B(IB)*(1.D+0+JZ(I,J)*BB(1)*TR)*(CX**JZ(I,J)))*(PR**I)
  120   CONTINUE
  130 CONTINUE
      IB = 18
      DO 180 I = 6,8
        ISB = 1
        T1 = 0.D+0
        DO 140 K = 1,JL(I)
          ISB = ISB+1
          T1 = T1 + BB(ISB)*(CX**JX(I,K))
  140   CONTINUE
        T5 = 0.D+0
        DO 170 J = 1,JN(I)
          IB = IB + 1
          T2 = 0.D+0
          ISB = 1
          DO 150 K = 1,JL(I)
            ISB = ISB + 1
            T2 = T2 + JX(I,K)*BB(ISB)*(CX**JX(I,K))
  150     CONTINUE
          T2 = T2*BB(1)*TR
          T3 = PR**(2-I)
          ISB = 1
          DO 160 K = 1,JL(I)
            ISB = ISB+1
            T3 = T3 + BB(ISB)*(CX**JX(I,K))
  160     CONTINUE
          T4 = (1.D+0 + JZ(I,J)*BB(1)*TR) - T2/T3
          T4 = T4*B(IB)*(CX**JZ(I,J))
          T5 = T5 + T4
  170   CONTINUE
      HX = HX - T5/(PR**(2-I))
  180 CONTINUE
      IB = 24
      T1 = PR*(PR/BL)**10
      DO 190 I = 1,7
        IB = IB + 1
        HX = HX +
     &    (1.D+0+TR*(1.D+1*DBL/BL+(I-1)*BB(1)))*B(IB)*(CX**(I-1))*T1
  190 CONTINUE
      H = HX*PCRW*VCRW*1.D-3/WTMW
      IF( P.LT.EPSL.OR.RHO.EQ.0.D+0 )THEN
        U = H
      ELSE
        U = H - P/RHO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATGSH group  ---
!
      RETURN
      END

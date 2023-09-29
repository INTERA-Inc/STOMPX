!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BNDFGC2( TX,PX,DIFX,IGC,JGC )
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
!     Calculates gas binary diffusion coefficient from
!     the Wilke and Lee theory. pp. 587.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
      USE CCP
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
      SAVE A,TP1X,TP2X,SIG1X,SIG2X,WTMX
      DATA A /1.06036D+0,1.5610D-1,1.9300D-1,4.7635D-1,1.03587D+0,
     &        1.52996D+0,1.76474D+0,3.89411D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BNDFGC2'
      TP1X = 1.15D+0*GCPP(2,IGC)
      TP2X = 1.15D+0*GCPP(2,JGC)
      SIG1X = 1.18D+0*(2.85D-1*(GCPP(5,IGC)**1.048D+0))**(1.0D+0/3.0D+0)
      SIG2X = 1.18D+0*(2.85D-1*(GCPP(5,JGC)**1.048D+0))**(1.0D+0/3.0D+0)
      WTMX = 2.D+0/((1.D+0/GCPP(1,IGC)) + (1.D+0/GCPP(1,JGC)))
      TKX = TX+TABS
      PB = PX*1.D-5
      TP = TKX/SQRT(TP1X*TP2X)
      SIG = 5.D-1*(SIG1X+SIG2X)
      OMG = (A(1)/TP**A(2)) + (A(3)/EXP(A(4)*TP))
     &  + (A(5)/EXP(A(6)*TP))  + (A(7)/EXP(A(8)*TP))
      DIFX = (3.03D+0 - (9.8D-1/SQRT(WTMX)))*1.D-3*(TKX**1.5D+0)/
     &  (PB*SQRT(WTMX)*(SIG**2)*(OMG))*1.D-4
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BNDFGC2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BNDFGCW( TX,PX,DIFX,IGC )
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
!     Calculates gas binary diffusion coefficient from
!     the Wilke and Lee theory. pp. 587.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
      USE CCP
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
      SAVE A,TPCX,TPWX,SIGCX,SIGWX,WTMX
      DATA A /1.06036D+0,1.5610D-1,1.9300D-1,4.7635D-1,1.03587D+0,
     &        1.52996D+0,1.76474D+0,3.89411D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BNDFGCW'
      TPCX = 1.15D+0*GCPP(2,IGC)
      TPWX = 1.15D+0*TBW
      SIGCX = 1.18D+0*(2.85D-1*(GCPP(5,IGC)**1.048D+0))**(1.0D+0/3.0D+0)
      SIGWX = 1.18D+0*(2.85D-1*(VCRW**1.048D+0))**(1.0D+0/3.0D+0)
      WTMX = 2.D+0/((1.D+0/GCPP(1,IGC)) + (1.D+0/WTMW))
      TKX = TX+TABS
      PB = PX*1.D-5
      TP = TKX/SQRT(TPCX*TPWX)
      SIG = 5.D-1*(SIGCX+SIGWX)
      OMG = (A(1)/TP**A(2)) + (A(3)/EXP(A(4)*TP))
     &  + (A(5)/EXP(A(6)*TP))  + (A(7)/EXP(A(8)*TP))
      DIFX = (3.03D+0 - (9.8D-1/SQRT(WTMX)))*1.D-3*(TKX**1.5D+0)/
     &  (PB*SQRT(WTMX)*(SIG**2)*(OMG))*1.D-4
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BNDFGCW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BNDFLC( TX,VISLX,VISGCX,DFLCX,IGC )
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
!     Calculates binary liquid diffusion coefficient for gas component
!     in liquid water, using the Wilke and Chang estimation
!      method; pp. 598.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BNDFLC'
      VBC = 2.85D-1*GCPP(5,IGC)**1.048D+0
      DFLCX = 7.4D-15*((2.6D+0*WTMW)**5.D-1)*(TX+TABS)/
     &  (VISLX*(VBC**6.D-1))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BNDFLC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_GC( TX,PX,XMGCX,RHOGX,ISRX,ZGX )
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
!     Density of liquid or vapor NGC components, excluding water vapor
!     with Peneloux volume correction.
! 
!     
!     ISRX liquid or vapor index: 1 - liquid 2 - vapor or supercritical
!       
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE NAPL
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER I_PX(2),I_SX,I_TX(2)
      REAL*8 TRX(30),SAX(30),BIJ(30,30),XMGCX(LNGC)
      REAL*8 YX(LNGC),XX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_GC'
      INCG = 1
      TCRA = 304.1D+0
      PCRA = 73.8D+5
      TKX = TX + TABS
      RCUBX = RCU*1.D-2
      PBX = PX*1.D-5
      
      IF(PX.GT.SMALL) THEN
        IF( (XMGCX(1).LE.0.D+0) .AND. (NGC.EQ.1) ) THEN
          RHOGX = 0.D+0
          RETURN
        ENDIF
        DO IGC = 1, NGC
         YX(IGC) = XMGCX(IGC)/SUM(XMGCX)
        ENDDO
        DO IGC = 1, NGC
          TRX(IGC) = TKX/GCPP(3,IGC)
          GCPP(15,IGC) = 0.37464D+0 + 1.54226D+0*GCPP(7,IGC) -
     &      0.26992D+0*(GCPP(7,IGC)**2.D+0)
          SAX(IGC) =  (1.D+0 + GCPP(15,IGC)*
     &      (1.D+0-SQRT(TRX(IGC))))**2.D+0
        ENDDO
!
!--- Peneloux correction factor ---
!
      CIX = 0.D+0
      CX = 0.D+0
      DO IGC = 1, NGC
! Correct Rackett compressibility factor for temperature and pressure
        ZRAC = GCPP(16,IGC) + GCPP(67,IGC)*(PX/1.d+8-GCPP(68,IGC))
     &          + GCPP(69,IGC)*(TX/8.d+2-GCPP(70,IGC))
        CIX = GCPP(65,IGC)*(GCPP(66,IGC)-ZRAC)*RCUBX*GCPP(3,IGC)/
     &   (GCPP(4,IGC)*1.D-5)
       CX = CX + YX(IGC)*CIX
      ENDDO
      SCMX = PBX*CX/(RCUBX*TKX)
!
!--- Binary interation coefficient ---
!
       DO IGC = 1, NGC
        DO JGC = 1, NGC
         BIJ(IGC,JGC) = BIJNA(IGC,JGC)+BIJTD(IGC,JGC)*TRX(IGC)
        ENDDO
       ENDDO
!
!--- Peng-Robinson temperature dependent parameters
!
          DO IGC = 1, NGC
           GCPP(13,IGC) = SAX(IGC)*0.45724D+0*(RCUBX**2.D+0)*
     &          ((GCPP(3,IGC)**2.D+0)/(GCPP(4,IGC)*1.D-5))*
     &          PBX/(RCUBX*TKX)**2.D+0
           GCPP(14,IGC) = (0.0778D+0*RCUBX*GCPP(3,IGC)/
     &          (GCPP(4,IGC)*1.D-5))*PBX/(RCUBX*TKX)
          ENDDO
!
!--- Peng-Robinson concentration dependent parameters
!
         IF( NGC.EQ.1 ) YX(NGC)=1.D+0
          SAMX = 0.D+0
          SBMX = 0.D+0
          DO IGC = 1, NGC
           DO JGC = 1, NGC
            SAMX = SAMX + YX(IGC)*YX(JGC)*SQRT(GCPP(13,IGC)*
     &          GCPP(13,JGC))*(1.D+0 - BIJ(IGC,JGC))
           ENDDO
          ENDDO
   
          DO IGC = 1, NGC
           SBMX = SBMX + YX(IGC)*GCPP(14,IGC)
          ENDDO
!
!---    Polynomial coefficients for Peng-Robinson equation of state  ---
!
         SCMX = 0.D+0
         CAX = 1.D+0
         CBX = -(1.0D+0 - SBMX - 3.D+0*SCMX)
         CCX = SAMX - 2.D+0*SBMX - 3.D+0*(SBMX**2.D+0) + 2.D+0*SBMX*SCMX
     &         - 2.D+0*SCMX + 3.D+0*(SCMX**2.D+0)
         CDX = -(SAMX*SBMX - (SBMX**2.D+0) - (SBMX**3.D+0) - SAMX*SCMX +
     &         2.D+0*SBMX*SCMX + 3.D+0*(SBMX**2.D+0)*SCMX 
     &         - (SCMX**2.D+0)*SBMX + (SCMX**2.D+0) - (SCMX**3.D+0))
         CALL NICKALLS( CAX,CBX,CCX,CDX,Z1X,Z2X,Z3X )
!
!---    Compressibility factor Z---
!
        ZGX = MAX( Z1X,Z2X,Z3X )
        ZNX = MIN( Z1X,Z2X,Z3X )
!
!---    Molar volume, cm^3/mol  ---
!
        VGX = RCUBX*TKX*ZGX/PBX
        VNX = RCUBX*TKX*ZNX/PBX

        DELTA1 = 1.D+0+SQRT(2.D+0)
        DELTA2 = 1.D+0-SQRT(2.D+0)

        GSTAR = (ZNX-ZGX)+LOG(ABS(((ZGX-SBMX)/(ZNX-SBMX))))+
     +   (1.D+0/(DELTA2-DELTA1))*(SAMX/SBMX)
     +   *LOG(ABS(((ZGX+DELTA2*SBMX)*(ZNX+DELTA1*SBMX)/
     +    ((ZNX+DELTA2*SBMX)*(ZGX+DELTA1*SBMX)))))
!
!---    Gas conditions  ---
!
        IF( NGC.EQ.1 )THEN
          TCIJ = GCPP(3,NGC)
        ELSE
        TCIJ = 0.D+0
        DO IGC = 1, NGC
         TCIJ = TCIJ + XMGCX(IGC)*GCPP(3,IGC)
        ENDDO
        ENDIF
       
        IF( TKX.GE.TCIJ ) THEN
          VX = VGX
          ZX = ZGX
        ELSE
         IF(GSTAR.GE.0.D+0) THEN
           VX = VGX
           ZX = ZGX
           GOTO 901
         ELSE
          VX = VNX
          ZX = ZNX
         GOTO 901
         ENDIF
          W1X = PBX*(VGX-VNX)
          W2X = RCUBX*TKX*LOG((VGX-SBMX)/(VNX-SBMX)) + 
     &      (SAMX/(SQRT(TKX)*SBMX))*LOG(((VGX+SBMX)*VNX)/(VNX+SBMX)*VGX)
!
!---      Gas phase condition  ---
!
          IF( (W1X-W2X)/EPSL.GE.EPSL ) THEN
            VX = VGX
            ZX = ZGX
!
!---      Liquid phase condition  ---
!
          ELSEIF( (W1X-W2X)/EPSL.LT.EPSL ) THEN
            VX = VNX
            ZX = ZNX
!
!---      Gas-liquid phase condition, use gas molar volume  ---
!
          ELSE
            VX = VGX
            ZX = ZGX
          ENDIF
        ENDIF
  901   CONTINUE
!
!--- Molecular weight of mixture without water 
!
        WTMX = 0.D+0
        DO IGC = 1, NGC
         WTMX = WTMX + YX(IGC)*GCPP(1,IGC)
        ENDDO
!
!--- Density of mixture without water ---
!
!
!---  Liquid density  ---
!
       IF( ISRX.EQ.1 ) THEN
          TKX = TX+TABS
          VOLX = ZX*RCU*TKX/PX
          RHOGX = WTMX/(VOLX)
!
!---  Vapor or supercritical density  ---
!
       ELSE
          TKX = TX+TABS
          VOGX = (ZX*RCUBX*TKX/PBX)+CX
          VOGX = VX+CX
          RHOGX = WTMX/(VOGX)
          RHOGX = RHOGX*1.D+3
       ENDIF
      ELSE 
        RHOGX = 0.D+0
      ENDIF
!       
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_GC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_B( TX,PX,XLSX,RHOBX )
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
!     This subroutine calculates the density of NaCl brine as a
!     function of temperature (C), pressure (Pa), and NaCl
!     concentration (molality, mol NaCl/kg H2O).
!
!     Haas Jr., J.L.  1976.  Physical Properties of the Coexisting
!     Phases and Thermochemical Properties of the H2O Component
!     in Boiling NaCl Solutions.  Preliminary Steam Tables for
!     NaCl Solutions.  Geological Survey Bulletin 1421-A.
!
!     Phillips, S.L., H. Ozbek, and L.F. Silvester.  1983.
!     Density of Sodium Chloride Solutions at High Temperatures and
!     Pressures, LBL-16275, Lawrence Berkeley Laboratory, University
!     of California, Berkeley, California.
!
!     Temperature Range: 0 - 350 C
!     Pressure Range:  0.1 - 100 MPa
!     NaCl Concentration Range:  0 - 5 Molal (mol NaCl/kg H2O)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 CCX(4),CHX(10),SCX(3),SAX(3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CCX,CHX,SCX,SAX,VCX
      DATA CCX / -3.033405D+0, 10.128163D+0, -8.750567D+0, 2.663107D+0 /
      DATA CHX / -167.219D+0, 448.55D+0, -261.07D+0, -13.644D+0,
     &  13.97D+0, -0.315154D+0, -1.203374D-3, 7.48908D-13,
     &  0.1342489D+0, -3.946963D-3 /
      DATA SCX / -9.9559D+0, 7.0845D+0, 3.9093D+0 /
      DATA SAX / -4.539D-3, -1.638D-4, 2.551D-5 /
      DATA VCX / 3.1975D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_B'
!
!---  Convert pressure to bar and mass fraction to molality  ---
!
      PBX = 1.D-5*PX
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
!
!---  Compressed or vapor-saturated density of pure water using the
!     ASME formulations  ---
!
      CALL SP_W( TX,PSWX )
      ISRX = 1
      PWX = MAX( PX,PSWX )
      CALL DENS_W( TX,PWX,RHOLWX,RHOVWX,ISRX )
!
!---  Convert density units to gm/cm^3  ---
!
      RHOLWX = 1.D-3*RHOLWX
      VOX = 1.D+0/RHOLWX
!
!---  Limiting apparent molal volume (cm^3/mol) of NaCl in solution as
!     the concentration goes to zero  ---
!
      PHIPX = CHX(1) + CHX(2)*VOX + CHX(3)*(VOX**2)
!
!---  Apparent molal volume (cm^3/mol) of NaCl in solution  ---
!
      PHIX = PHIPX + (CHX(4)+CHX(5)*VOX)*((VOX/(VCX-VOX))**2)*SQRT(GLSX)
!
!---  Brine density (gm/cm^3)  ---
!
      RHOBX = (1.D+3 + GLSX*WTMS)/(1.D+3*VOX + GLSX*PHIX)
!
!---  Convert density units to kg/m^3  ---
!
      RHOBX = 1.D+3*RHOBX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_L( TX,RHOBX,XLCX,RHOLX )
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
!     This subroutine calculates the density of CO2 gas dissolved
!     in NaCl brine as a function of the brine density and dissolved
!     gas mass fraction.
!
!     Alendal, G., and H. Drange.  2001.  Two-phase, near-field,
!     modeling of purposefully released CO2 in the ocean.  Journal
!     of Geophysical Research, 106(C1):1085-1096.
!
!     Anderson, G.M., and D.A. Crerar.  1992.  Thermodynamics in
!     Geochemistry: The Equilibrium Model, Oxford University Press.
!
!     Variable definitions
!
!     PMV_A - molar volume of CO2 (m^3/kmol)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 PMV_GC(LNGC)
      REAL*8 XLCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_L'

      XLAX = 0.D+0
      CXA = 0.D+0
      DO IGC = 1,NGC
        XLAX = XLAX + XLCX(IGC)
!
!---  Partial molar volume of CO2 (m^3/kmol), formulation of
!     Anderson et al. (1992)  ---
!
        PMV_GC(IGC) = 0.D+0
        DO 10 M = 1,5
          PMV_GC(IGC) = PMV_GC(IGC) + GCPP(M+22,IGC)*(TX**(M-1))
   10   CONTINUE
        CXA = CXA + PMV_GC(IGC)*RHOBX*XLCX(IGC)/GCPP(1,IGC)
      ENDDO
!
!---  Poynting correction  ---
!
      RHOLX = RHOBX/(1.D+0 + CXA - XLAX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_L group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_S( TX,PX,RHOSX )
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
!     Density of precipitated NaCl.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas.  Geothermics, 26(4): 437-464.
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_S'
!
!---  Formulation of Battistelli et al.  ---
!
      RHOSX = 2.165D+3*EXP(-1.2D-4*TX + 4.D-11*PX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_S group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_W( TX,PX,RHOLX,RHOVX,ISRX )
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
!     Density (kg/m^3) of pure water as a function of temperature and
!     pressure.
!
!     Pressure Range: 0 - 100 MPa (1000 Bar)
!     Temperature Range:  273.16 K (0.01 C) to 1073.15K (800 C)
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 CAX(23),SAX(12),CBX(31),SBX(5),LFCX(3),SLX
      REAL*8 TX,PX,RHOLX,RHOVX
      INTEGER ISRX,INX(8),IZX(8,3),ITX(8),IXX(8,2)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CAX,SAX,CBX,SBX
      SAVE LFCX,SLX,INX,IZX,ITX,IXX
      DATA CAX / 6.824687741D+3, -5.422063673D+2, -2.096666205D+4,
     &  3.941286787D+4, -6.733277739D+4, 9.902381028D+4,
     &  -1.093911774D+5, 8.590841667D+4, -4.511168742D+4,
     &  1.418138926D+4, -2.017271113D+3, 7.982692717D+0,
     &  -2.616571843D-2, 1.522411790D-3, 2.284279054D-2,
     &  2.421647003D+2, 1.269716088D-10, 2.074838328D-7,
     &  2.174020350D-8, 1.105710498D-9, 1.293441934D+1,
     &  1.308119072D-5, 6.047626338D-14 /
      DATA SAX / 8.438375405D-1, 5.362162162D-4, 1.720000000D+0,
     &  7.342278489D-2, 4.975858870D-2, 6.537154300D-1,
     &  1.150000000D-6, 1.150800000D-5, 1.418800000D-1,
     &  7.002753165D+0, 2.995284926D-4, 2.040000000D-1 /
      DATA CBX / 1.683599274D+1, 2.856067796D+1, -5.438923329D+1,
     &  4.330662834D-1, -6.547711697D-1, 8.565182058D-2,
     &  6.670375918D-2, 1.388983801D+0, 8.390104328D-2,
     &  2.614670893D-2, -3.373439453D-2, 4.520918904D-1,
     &  1.069036614D-1, -5.975336707D-1, -8.847535804D-2,
     &  5.958051609D-1, -5.159303373D-1, 2.075021122D-1,
     &  1.190610271D-1, -9.867174132D-2, 1.683998803D-1,
     &  -5.809438001D-2, 6.552390126D-3, 5.710218649D-4,
     &  1.936587558D+2, -1.388522425D+3, 4.126607219D+3,
     &  -6.508211677D+3, 5.745984054D+3, -2.693088365D+3,
     &  5.235718623D+2/
      DATA SBX / 7.633333333D-1, 4.006073948D-1, 8.636081627D-2,
     &  -8.532322921D-1, 3.460208861D-1 /
      DATA LFCX / 1.574373327D+1, -3.417061978D+1, 1.931380707D+1 /
      DATA SLX / 4.260321148D+0 /
      DATA INX / 2, 3, 2, 2, 3, 2, 2, 2 /
      DATA IZX / 13, 18, 18, 25, 32, 12, 24, 24,
     &  3, 2, 10, 14, 28, 11, 18, 14,
     &  0, 1, 0, 0, 24, 0, 0, 0 /
      DATA ITX / 0, 0, 0, 0, 0, 1, 1, 2 /
      DATA IXX / 0, 0, 0, 0, 0, 14, 19, 54,
     &  0, 0, 0, 0, 0, 0, 0, 27 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_W'
!
!---  Check temperature and pressure  ---
!
      IF( TX.LT.1.D-2 .OR. TX.GT.8.D+2 ) THEN
!        INDX = 13
!        CHMSG = 'Temperature Out of Range, C: '
!        RLMSG = TX
!        CALL WRMSGS( INDX )
        TX = MAX( TX,1.D-2 )
        TX = MIN( TX,8.D+2 )
      ENDIF
      IF( PX.LT.0.D+0 .OR. PX.GT.1.D+8 ) THEN
!        INDX = 13
!        CHMSG = 'Pressure Out of Range, Pa: '
!        RLMSG = PX
!        CALL WRMSGS( INDX )
        PX = MAX( PX,0.D+0 )
        PX = MIN( PX,1.D+8 )
      ENDIF
      IF( PX/EPSL.LT.EPSL ) THEN
        RHOLX = 0.D+0
        RHOVX = 0.D+0
        GOTO 100
      ENDIF
!
!---  Reduced temperature and pressure  ---
!
      RHOLX = 0.D+0
      RHOVX = 0.D+0
      THETAX = (TX+TABS)/TCRW
      BETAX = PX/PCRW
!
!---  L-function  ---
!
      BETALX = LFCX(1) + LFCX(2)*THETAX + LFCX(3)*(THETAX**2)
!
!---  Subregions 1 or 6  ---
!
      IF( ISRX.EQ.1 .OR. ISRX.EQ.6 ) THEN
        CYX = 1.D+0 - SAX(1)*(THETAX**2) - SAX(2)/(THETAX**6)
        CZX = CYX + SQRT((SAX(3)*(CYX**2)) - (2.D+0*SAX(4)*THETAX) +
     &    (2.D+0*SAX(5)*BETAX))
        RVX = CAX(12)*SAX(5)*(CZX**(-5.D+0/17.D+0))
        RVX = RVX +  (CAX(13) + CAX(14)*THETAX + CAX(15)*(THETAX**2) +
     &    CAX(16)*((SAX(6)-THETAX)**10) + CAX(17)/(SAX(7)+(THETAX**19)))
        RVX = RVX - (CAX(18) + 2.D+0*CAX(19)*BETAX +
     &    3.D+0*CAX(20)*(BETAX**2))/(SAX(8)+(THETAX**11))
        RVX = RVX - CAX(21)*(THETAX**18)*(SAX(9)+(THETAX**2))*
     &    (-3.D+0/((SAX(10)+BETAX)**4))
        RVX = RVX + 3.D+0*CAX(22)*(SAX(12)-THETAX)*(BETAX**2)
        RVX = RVX + 4.D+0*CAX(23)*(BETAX**3)/(THETAX**20)
        RHOLX = (1.D+3*WTMW)/(RVX*VCRW)
      ENDIF
!
!---  Subregions 2 or 6  ---
!
      IF( ISRX.EQ.2 .OR. ISRX.EQ.6 ) THEN
        CXX = EXP(SBX(1)*(1.D+0-THETAX))
        RVX = SLX*THETAX/BETAX
        ICBX = 6
        DO 30 I = 1,5
          RVAX = 0.D+0
          DO 20 J = 1,INX(I)
            ICBX = ICBX + 1
            RVAX = RVAX + CBX(ICBX)*(CXX**IZX(I,J))
   20     CONTINUE
          REALX = REAL(I)
          RVX = RVX - REALX*(BETAX**(I-1))*RVAX
   30   CONTINUE
        ICBX = 18
        ISBX = 1
        DO 60 I = 6,8
          RVAX = 0.D+0
          DO 40 J = 1,INX(I)
            ICBX = ICBX + 1
            RVAX = RVAX + CBX(ICBX)*(CXX**IZX(I,J))
   40     CONTINUE
          RVBX = 0.D+0
          DO 50 J = 1,ITX(I)
            ISBX = ISBX + 1
            RVBX = RVBX + SBX(ISBX)*(CXX**IXX(I,J))
   50     CONTINUE
          REALX = REAL(I)
          RVX = RVX - ((REALX-2.D+0)*(BETAX**(1-I))*RVAX)/
     &      (((BETAX**(2-I))+RVBX)**2)
   60   CONTINUE
        ICBX = 24
        RVAX = 0.D+0
        DO 70 I = 0,6
          ICBX = ICBX + 1
          RVAX = RVAX + CBX(ICBX)*(CXX**I)
   70   CONTINUE
        RVX = RVX + 1.1D+1*((BETAX/BETALX)**10)*RVAX
        RHOVX = (1.D+3*WTMW)/(RVX*VCRW)
      ENDIF
!
!---  Subregions 3 or 5  ---
!
      IF( ISRX.EQ.3 .OR. ISRX.EQ.5 ) THEN
        INDX = 14
        CHMSG = 'Steam Table Subregion 3/5: Exceeds Temperature Range: '
        RLMSG = TX
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Subregions 4 or 5  ---
!
      IF( ISRX.EQ.4 .OR. ISRX.EQ.5 ) THEN
        INDX = 14
        CHMSG = 'Steam Table Subregion 4/5: Exceeds Temperature Range: '
        RLMSG = TX
        CALL WRMSGS( INDX )
      ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DIFC_GW( TX,PX,DFGWX )
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
!
!     Calculates the CO2-H2O binary diffusion coefficient from
!     the Wilke and Lee theory. pp. 587
!
!     Liquid molar volume is computed from critical molar volume with
!     the Tyn and Calus method. pp. 53-54
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids. pp. 587.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 CAX(8),CLJ_W(2),CLJ_A(2)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CAX,CLJ_W,CLJ_A
      DATA CAX /1.06036D+0,1.5610D-1,1.9300D-1,4.7635D-1,1.03587D+0,
     &  1.52996D+0,1.76474D+0,3.89411D+0/
      DATA CLJ_W / 3.190008977D+0, 429.18D+0 /
      DATA CLJ_A / 3.795165630D+0, 95.85245D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DIFC_GW'
!
!---  Absolute temperature and pressure in bars  ---
!
      TKX = TX+TABS
      PBX = PX*1.D-5
!
!---  Mean characteristic Lennard-Jones energy and length (angstrom) ---
!
      WTMA = GCPP(1,1)
      EPSX = SQRT(CLJ_W(2)*CLJ_A(2))
      SIGX = 5.D-1*(CLJ_W(1)+CLJ_A(1))
      TPX = TKX/EPSX
      OMGX = (CAX(1)/TPX**CAX(2)) + (CAX(3)/EXP(CAX(4)*TPX))
     &  + (CAX(5)/EXP(CAX(6)*TPX))  + (CAX(7)/EXP(CAX(8)*TPX))
      WTMX = 2.D+0/((1.D+0/WTMA) + (1.D+0/WTMW))
      DFGWX = (3.03D+0 - (9.8D-1/SQRT(WTMX)))*1.D-3*(TKX**1.5D+0)/
     &  (PBX*SQRT(WTMX)*(SIGX**2)*(OMGX))*1.D-4
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIFC_GW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DIFC_LS( TX,XLSX,VISLX,DFLSX )
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
!     Calculates the diffusion coefficient (m^2/s) for NaCl in
!     aqueous solutions, following the method of Nernst-Haskell for
!     dilute solutions and the method of Gordon for concentrated
!     solutions; where the mean ionic activity is computed according
!     to the method of Bromley.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids. pp. 620-621.
!     McGraw-Hill, New York, New York
!
!     Bromley, L.A.  1973.  Thermodynamic properties of strong
!     electrolytes in aqueous solutions.  AIChE Journal, 19(2):313-320
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
!----------------------Data Statements---------------------------------!
!
      SAVE CBX,TRX,VISWRX
      DATA CBX / 0.0547D+0 /
      DATA TRX / 25.D+0 /
      DATA VISWRX / 0.8904339807D-3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DIFC_LS'
!
!---  Convert halite mass fraction to molaity   ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
      TKX = TX+TABS
!
!---  Partial derivative of the natural logarithm  of the mean ionic
!     activity with respect to the molaity at 298 K   ---
!
      IF( GLSX.GT.EPSL ) THEN
        DLNGX = (-0.2555D+0/(SQRT(GLSX)*(1.D+0+SQRT(GLSX))) +
     &    0.2555D+0/((1.D+0+SQRT(GLSX))**2) +
     &    (6.D-2 + 6.D-1*CBX)/((1.D+0+1.5D+0*GLSX)**2) -
     &    3.D+0*(6.D-2 + 6.D-1*CBX)*GLSX/((1.D+0+1.5D+0*GLSX)**3) +
     &    CBX)*2.302585D+0
      ELSE
        DLNGX = -5.883D-6
      ENDIF
!
!---  Diffusion coefficient (m^2/s) for dilute NaCl aqueous solutions
!     at 298 K   ---
!
      DFLSX = 2.254D-9
!
!---  Viscosity of brine as a function of NaCl mass fraction
!     at 298.15 K   ---
!
      CALL VISC_B( TRX,XLSX,VISWRX,VISBRX )
!
!---  Diffusion coefficient for concentrated NaCl aqueous solutions
!     at 298.15 K   ---
!
      DFLSX = DFLSX*(VISWRX/VISBRX)*(1.D+0 + GLSX*(DLNGX))
!
!---  Viscosity of brine as a function of NaCl mass fraction  ---
!
      CALL VISC_B( TKX,XLSX,VISLX,VISBX )
!
!---  Correct diffusion coefficient for temperature  ---
!
      DFLSX = DFLSX*(TKX/2.9815D+2)*(VISBRX/VISBX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIFC_LS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_GC( TX,PX,XMGCX,XMGWX,WTMX,HGAX,UGAX,ISRX )
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
!     Enthalpy and internal energy of CO2.
!
!     Span, R., and W. Wagner.  1996.  A New Equation of State for
!     Carbon Dioxide Covering the Fluid Region from the Triple-Point
!     to 1100 K at Pressures up to 800 MPa.  J. Phys. Chem. Ref. Data
!     25(6):1509-1588.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER I_PX(2),I_SX,I_TX(2)
!--RR
      CHARACTER*132 CHMSGX(1)
      INTEGER NC,JNCG
      REAL*8 XX(2),PC_I,TC_I,PC_J,TC_J,VC_I,VC_J,VCIJ
      REAL*8 TCIJ,PCM,TCM,VCM,OMEGAXM,OMEGAX,OMEGAXR,OMEGAXW
      REAL*8 PRX,B1O,B2O,B3O,B4O,C1O,C2O,C3O,C4O,D1O,D2O
      REAL*8 B1R,B2R,B3R,B4R,C1R,C2R,C3R,C4R,D1R,D2R
      REAL*8 BETAXO,GAMMAXO,BETAXR,GAMMAXR,BLK,CLK,DLK
      REAL*8 V0,V0N,FX,FXD,ZX,HDEP0,HDEPR,HDEP,E1O,E1R,HAX,UAX
      REAL*8 Z0,ZR,HCPAX,HCPWX,HCPMX,DELHFX,DELHFXW
      REAL*8 HCPAXA,HCPAXB,HCPAXC,HCPAXD,HCPWXA,HCPWXB,HCPWXC,HCPWXD
      REAL*8 HDEPMX,HDEPRMX,PRF,TRF,PMAX
      REAL*8 XMGCX(LNGC),BIJ(LNGC,LNGC),VCX(LNGC,LNGC),SAX(30)
      REAL*8 TRX(LNGC+1),YX(LNGC+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_GC'
      INCG = 1
      RHOCRA = 467.6D+0
      TKX = TX+TABS
      PBX = PX*1.D-5
      PMX = PX*1.D-6
      RCUMPX = RCU*1.D-2*1D-1
      RCUBX = RCU*1.D-2
      OMEGAXR = 0.3978
      OMEGAX = 0.239
      OMEGAW = 0.344
      JNCG = 2
      P_TXM = 1.D-4
      XVC = 0.D+0
      VCGMX = 0.D+0
      OMEGAGMX = 0.D+0
      TCGMX = 0.D+0
!
!--- Critical properties. Pressure in MPA and temperature in Kelvin
!     
      PC_I = 73.8*1.D-1
      TC_I = 304.1
      VC_I = 93.9
      PC_J = 221.2*1.D-1
      TC_J = 647.3
      VC_J = 57.1
!      KIJ = 0.920
!
!--- Constants ---
!  
      RTX = RCUBX*TKX
      SQR8X = SQRT(8.D+0)
      SQR2X = SQRT(2.D+0)
        DO IGC = 1, NGC
         YX(IGC) = XMGCX(IGC)/SUM(XMGCX)
        ENDDO

        DO IGC = 1, NGC
          TRX(IGC) = TKX/GCPP(3,IGC)
          GCPP(15,IGC) = 0.37464D+0 + 1.54226D+0*GCPP(7,IGC) -
     &      0.26992D+0*(GCPP(7,IGC)**2.D+0)
          SAX(IGC) =  (1.D+0 + GCPP(15,IGC)*
     &      (1.D+0-SQRT(TRX(IGC))))**2.D+0
        ENDDO
!
!--- Binary interation coefficient ---
!
       DO IGC = 1, NGC
        DO JGC = 1, NGC
         BIJ(IGC,JGC) = BIJNA(IGC,JGC)+BIJTD(IGC,JGC)*TRX(IGC)
        ENDDO
       ENDDO
!
!--- Peng-Robinson temperature dependent parameters
!
          DO IGC = 1, NGC
           GCPP(13,IGC) = SAX(IGC)*0.45724D+0*(RCUBX**2.D+0)*
     &          ((GCPP(3,IGC)**2.D+0)/(GCPP(4,IGC)*1.D-5))*
     &          PBX/(RCUBX*TKX)**2.D+0
           GCPP(14,IGC) = (0.0778D+0*RCUBX*GCPP(3,IGC)/
     &          (GCPP(4,IGC)*1.D-5))*PBX/(RCUBX*TKX)
          ENDDO
!
!--- Peng-Robinson concentration dependent parameters
!
         IF( NGC.EQ.1 ) YX(NGC)=1.D+0
          SAMX = 0.D+0
          SBMX = 0.D+0
          DO IGC = 1, NGC
           DO JGC = 1, NGC
            SAMX = SAMX + YX(IGC)*YX(JGC)*SQRT(GCPP(13,IGC)*
     &          GCPP(13,JGC))*(1.D+0 - BIJ(IGC,JGC))
           ENDDO
          ENDDO
          SUMXA = SAMX*((RCUBX*TKX)**2)/PBX
   
          DO IGC = 1, NGC
           SBMX = SBMX + YX(IGC)*GCPP(14,IGC)
          ENDDO
          SUMXB = SBMX*(RCUBX*TKX)/PBX
       DADT = 0.D+0
       DO IGC = 1, NGC
        DO JGC = 1, NGC
         AXI = GCPP(13,IGC)*((RCUBX*TKX)**2/PBX)
         AXJ = GCPP(13,JGC)*((RCUBX*TKX)**2/PBX)
         DADT = DADT+0.5*(-RCUBX)*SQRT(0.45724D+0/TKX)*(YX(IGC)*
     &          YX(JGC)*(1-BIJNA(IGC,JGC)))*((GCPP(15,IGC)*
     &          SQRT(AXI*GCPP(3,IGC)/(GCPP(4,IGC)*1.D-5)))+
     &          GCPP(15,JGC)*(SQRT(AXJ*GCPP(3,JGC)
     &         /(GCPP(4,JGC)*1.D-5))))
          ENDDO
         ENDDO
! 
!---    Polynomial coefficients for Peng-Robinson equation of state  ---
!
         CAX = 1.D+0
         CBX = -(1.0D+0 - SBMX)
         CCX = SAMX - 2.D+0*SBMX - 3.D+0*(SBMX**2.D+0)
         CDX = -(SAMX*SBMX - (SBMX**2.D+0) - (SBMX**3.D+0))

         CALL NICKALLS( CAX,CBX,CCX,CDX,Z1X,Z2X,Z3X )
!
!---    Compressibility factor Z---
!
        ZGX = MAX( Z1X,Z2X,Z3X )
        ZNX = MIN( Z1X,Z2X,Z3X )
!
!---    Molar volume, cm^3/mol  ---
!
        VGX = RCUBX*TKX*ZGX/PBX
        VNX = RCUBX*TKX*ZNX/PBX
        DELTA1 = 1.D+0+SQRT(2.D+0)
        DELTA2 = 1.D+0-SQRT(2.D+0)
        GSTAR = (ZNX-ZGX)+LOG(ABS(((ZGX-SBMX)/(ZNX-SBMX))))+
     +   (1.D+0/(DELTA2-DELTA1))*(SAMX/SBMX)
     +   *LOG(ABS(((ZGX+DELTA2*SBMX)*(ZNX+DELTA1*SBMX)/
     +    ((ZNX+DELTA2*SBMX)*(ZGX+DELTA1*SBMX)))))
!
!---    Gas conditions  ---
!
      SFX = 2.D+0/3.D+0
      TCIJ = 0.D+0
      SUMYX = 0.D+0
      DO IGC = 1,NGC
        SUMYX = SUMYX + YX(IGC)
      ENDDO
      IF( NGC.EQ.1 .OR. SUMYX.LT.EPSL ) THEN
        TCIJ = GCPP(3,NGC)
      ELSE
      CR1MX = 0.D+0
      DO IGC = 1, NGC
        CR1MX = CR1MX + YX(IGC)*GCPP(5,IGC)
        TCIJ = TCIJ + YX(IGC)*GCPP(5,IGC)*GCPP(3,IGC)
      ENDDO
      TCIJ = TCIJ/CR1MX
      ENDIF
        IF( TKX.GE.TCIJ ) THEN
          VX = VGX
          ZX = ZGX
        ELSE
           IF(GSTAR.GE.0.D+0) THEN
           VX = VGX
           ZX = ZGX
           GOTO 951
         ELSE
          VX = VNX
          ZX = ZNX
         GOTO 951
         ENDIF

          W1X = PBX*(VGX-VNX)
          W2X = RCUBX*TKX*LOG((VGX-SBMX)/(VNX-SBMX)) +
     &      (SAMX/(SQRT(TKX)*SBMX))*LOG(((VGX+SBMX)*VNX)/(VNX+SBMX)*VGX)
!
!---      Gas phase condition  ---
!
          IF( (W2X-W1X)/EPSL.LE.EPSL ) THEN
            VX = VGX
            ZX = ZGX
!
!---      Liquid phase condition  ---
!
          ELSEIF( (W1X-W2X)/EPSL.LT.EPSL ) THEN
            VX = VNX
            ZX = ZNX
!
!---      Gas-liquid phase condition, use gas molar volume  ---
!
          ELSE
            VX = VGX
            ZX = ZGX
          ENDIF
        ENDIF
  951  CONTINUE
!  
!---    Enthalpy departure function at temperature
!       and pressure conditions  ---
!
        HDEP = ((-SUMXA/(2.D+0*SUMXB*SQRT(2.D+0)))+
     &      (TKX*DADT/(2.D+0*SUMXB*
     &       sqrt(2.D+0))))*LOG((ZX-0.414D+0*SBMX)/(ZX+2.414D+0*SBMX))-
     &       RCUBX*TKX*(ZX-1.D+0)
        VXX = ZX*TKX/PBX
!
!--- Ideal gas enthalpy. Using integral of specific heat capacity ---
!
      TSX = 2.9815D+2
      HCPAX = 0.D+0
         DO IGC = 1, NGC
          HCPAX = HCPAX + YX(IGC)*(GCPP(9,IGC)*(TKX-TSX) 
     &                  + GCPP(10,IGC)*(TKX**2.D+0-TSX**2.D+0)/2.D+0 
     &                  + GCPP(11,IGC)*(TKX**3.D+0-TSX**3.D+0)/3.D+0 
     &                  + GCPP(12,IGC)*(TKX**4.D+0-TSX**4.D+0)/4.D+0)
         ENDDO
!
!---Enthalpy calculations at reference temperature and pressure
!
        TSAVE = TX
        PSAVE = PX
        TX = 2.5D+1
        TKX = TX + TABS
        PBX = 101325.D+0*1.D-5
        SQR8X = SQRT(8.D+0)
        SQR2X = SQRT(2.D+0)

        DO IGC = 1, NGC
          TRX(IGC) = TKX/GCPP(3,IGC)
          GCPP(15,IGC) = 0.37464D+0 + 1.54226D+0*GCPP(7,IGC) -
     &      0.26992D+0*(GCPP(7,IGC)**2.D+0)
          SAX(IGC) =  (1.D+0 + GCPP(15,IGC)*
     &      (1.D+0-SQRT(TRX(IGC))))**2.D+0
        ENDDO
!
!
!--- Peng-Robinson temperature dependent parameters
!
          DO IGC = 1, NGC
           GCPP(13,IGC) = SAX(IGC)*0.45724D+0*(RCUBX**2.D+0)*
     &          ((GCPP(3,IGC)**2.D+0)/(GCPP(4,IGC)*1.D-5))*
     &          PBX/(RCUBX*TKX)**2.D+0
           GCPP(14,IGC) = (0.0778D+0*RCUBX*GCPP(3,IGC)/
     &          (GCPP(4,IGC)*1.D-5))*PBX/(RCUBX*TKX)
          ENDDO
!
!--- Peng-Robinson concentration dependent parameters
!
         IF( NGC.EQ.1 ) YX(NGC)=1.D+0
          SAMX = 0.D+0
          SBMX = 0.D+0
          DO IGC = 1, NGC
           DO JGC = 1, NGC
            SAMX = SAMX + YX(IGC)*YX(JGC)*SQRT(GCPP(13,IGC)*
     &          GCPP(13,JGC))*(1.D+0 - BIJ(IGC,JGC))
           ENDDO
          ENDDO
          SUMXA = SAMX*((RCUBX*TKX)**2)/PBX
   
          DO IGC = 1, NGC
           SBMX = SBMX + YX(IGC)*GCPP(14,IGC)
          ENDDO
          SUMXB = SBMX*(RCUBX*TKX)/PBX
       DADT = 0.D+0
       DO IGC = 1, NGC
        DO JGC = 1, NGC
         AXI = GCPP(13,IGC)*((RCUBX*TKX)**2/PBX)
         AXJ = GCPP(13,JGC)*((RCUBX*TKX)**2/PBX)
         DADT = DADT+(-RCUBX)*0.5*SQRT(0.45724D+0/TKX)*
     &           (YX(IGC)*YX(JGC)*
     &         (1-BIJNA(IGC,JGC)))*(GCPP(15,IGC)*(SQRT(AXI*GCPP(3,IGC)/
     &         (GCPP(4,IGC)*1.D-5)))+GCPP(15,JGC)*(SQRT(AXJ*GCPP(3,JGC)
     &         /(GCPP(4,JGC)*1.D-5))))
          ENDDO
         ENDDO
! 
!---    Polynomial coefficients for Peng-Robinson equation of state  ---
!
         CAX = 1.D+0
         CBX = -(1.D+0 - SBMX)
         CCX = SAMX - 2.D+0*SBMX - 3.D+0*(SBMX**2.D+0)
         CDX = -(SAMX*SBMX - (SBMX**2.D+0) - (SBMX**3.D+0))

         CALL NICKALLS( CAX,CBX,CCX,CDX,Z1X,Z2X,Z3X )
!
!---    Compressibility factor Z---
!
        ZGX = MAX( Z1X,Z2X,Z3X )
        ZNX = MIN( Z1X,Z2X,Z3X )
!
!---    Molar volume, cm^3/mol  ---
!
        VGX = RCUBX*TKX*ZGX/PBX
        VNX = RCUBX*TKX*ZNX/PBX
!
!---    Gas conditions  ---
!
      SFX = 2.D+0/3.D+0
      SUMYX = 0.D+0
      DO IGC = 1,NGC
        SUMYX = SUMYX + YX(IGC)
      ENDDO
      IF( NGC.EQ.1 .OR. SUMYX.LT.EPSL ) THEN
        TCIJ = GCPP(3,NGC)
      ELSE
      TCIJ = 0.D+0
      CR1MX = 0.D+0
      DO IGC = 1, NGC
        CR1MX = CR1MX + YX(IGC)*GCPP(5,IGC)
        TCIJ = TCIJ + YX(IGC)*GCPP(5,IGC)*GCPP(3,IGC)
      ENDDO
      TCIJ = TCIJ/CR1MX
      ENDIF
        IF( TKX.GE.TCIJ ) THEN
          VX = VGX
          ZX = ZGX
        ELSE
          W1X = PBX*(VGX-VNX)
          W2X = RCUBX*TKX*LOG((VGX-SBMX)/(VNX-SBMX)) +
     &      (SAMX/(SQRT(TKX)*SBMX))*LOG(((VGX+SBMX)*VNX)/(VNX+SBMX)*VGX)
!
!---      Gas phase condition  ---
!
          IF( (W2X-W1X)/EPSL.GE.EPSL ) THEN
            VX = VGX
            ZX = ZGX
!
!---      Liquid phase condition  ---
!
          ELSEIF( (W1X-W2X)/EPSL.GE.EPSL ) THEN
            VX = VNX
            ZX = ZNX
!
!---      Gas-liquid phase condition, use gas molar volume  ---
!
          ELSE
            VX = VGX
            ZX = ZGX
          ENDIF
        ENDIF
!  
!---    Enthalpy departure function at temperature
!       and pressure conditions  ---
!
        HDEPR = ((-SUMXA/(2.D+0*SUMXB*SQRT(2.D+0)))+
     &      (TKX*DADT/(2.D+0*SUMXB*
     &       sqrt(2.D+0))))*LOG((ZX-0.414D+0*SBMX)/(ZX+2.414D+0*SBMX))-
     &       RCUBX*TKX*(ZX-1.D+0)
  900 CONTINUE
        WTMX = 0.D+0
        DO IGC = 1, NGC
         WTMX = WTMX + YX(IGC)*GCPP(1,IGC)
        ENDDO
         HGAX =   HCPAX - HDEP*(RCUBX/1.D+1)/RCUBX + 
     &            HDEPR*(RCUBX/1.D+1)/RCUBX
         HGAX = HGAX/WTMX
         TX = TSAVE
         PX = PSAVE
         UGAX = HGAX - PBX*VXX/WTMX
!
!---  Convert units to J/kg  ---
!
!  900 CONTINUE
      HGAX = HGAX*1.D+3
      UGAX = UGAX*1.D+3
!
!--- Liquid Enthalpy
!
      IF(ISRX==1) THEN
       DO IGC = 1, NGC
       IF(XMGCX(IGC).GT.EPSL) THEN
        IF(TKX.LT.GCPP(3,IGC)) THEN
            CALL ENTH_VAP( DHVX,TX,IGC )
            HGAX = HGAX - DHVX
          ELSE
            HGAX = HGAX 
          ENDIF
        ENDIF
        ENDDO
      ENDIF
!     WRITE(134,*)'Enthalpy and internal energy', HGCX,UGCX,IGC
      
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_GC group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_VAP( DHVX,TX,IGC )
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
!     Enthalpy of vaporization, using the two-fluid reference
!     model (pp. 221), J/kg
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 CR1(6),CR2(6),TAU(6)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CR1,CR2,OMGR1,OMGR2
      DATA CR1 /6.537D+0,-2.467D+0,-77.521D+0,
     &  59.634D+0,36.009D+0,-14.606D+0/
      DATA CR2 /-0.133D+0,-28.215D+0,-82.958D+0,
     &  99.000D+0,19.105D+0,-2.796D+0/
      DATA OMGR1 /0.21D+0/
      DATA OMGR2 /0.46D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_VAP'
!
!---  Sequence of tau powers  ---
!
      TKX = TX + TABS
      TAU(4) = 1.D+0 - TKX/GCPP(3,IGC)
      TAU(1) = TAU(4)**(THIRD)
      TAU(2) = TAU(4)**(5.D+0/6.D+0)
      TAU(3) = TAU(4)**1.208D+0
      TAU(5) = TAU(4)**2
      TAU(6) = TAU(4)**3
      DHVR1 = 0.D+0
      DHVR2 = 0.D+0
      DO 100 I = 1,6
        DHVR1 = DHVR1 + CR1(I)*TAU(I)
        DHVR2 = DHVR2 + CR2(I)*TAU(I)
        DHVR2 = DHVR2 + CR2(I)*TAU(I)
  100 CONTINUE
      DHVX = DHVR1 + ((GCPP(7,IGC)-OMGR1)/(OMGR2-OMGR1))*DHVR2
!
!---  Enthalpy of vaporization, J/kg  ---
!
      DHVX = DHVX*RCU*GCPP(3,IGC)/GCPP(1,IGC)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_VAP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_B( TX,XLSX,HLWX,HBX )
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
!     This subroutine calculates the enthalpy of NaCl solutions
!     as a function of temperature and NaCl concentration.
!
!     tx - temperature, C
!     xlsx - mass fraction of aqueous NaCl
!     tkbx - thermal conductivity of NaCl brine, W/m K
!     HLWX - Enthalpy of pure water at vapor-saturated conditions
!
!     Michaelides, E.E.  1981. Thermodynamic properties of geothermal
!     fluids.  Geothermal Resources Council, Transactions 5:361-364.
!
!     Gudmundsson, J.S., and H. Thrainsson.  1989.  Power potential of
!     two-phase geothermal wells.  Geothermics 18(3):357-366.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 SAX(12),SBX(3),SCX(4),CAX(6),CBX(15)
!
!----------------------Data Statements---------------------------------!
!
      SAVE SAX,SBX,SCX,CAX,CBX
      DATA SAX / 9633.6D+0, -4080.0D+0, 286.49D+0, 166.58D+0,
     &  68.577D+0, -4.6856D+0, -0.90963D+0, -0.36524D+0,
     &  0.249667D-1, 0.17965D-2, 0.71924D-3, -0.4900D-4 /
      DATA SBX / -0.83624D-3, 0.16792D+0, -25.9293D+0 /
      DATA SCX / 0.12453D-4, -0.45137D-2, 4.81155D+0, -29.578D+0 /
      DATA CAX / 25.19D+0, 0.1973D+0, -6.0114D-04, 8.81505D-7,
     &  -4.76500D-10, -1.923188214D+5 /
      DATA CBX / -104.51D+0, 81.086D+0, -308.22D+0, -1.6952D+0,
     &  -16.65D+0, -8.6385D+0, 0.010618D+0, 0.029634D+0, 0.61366D+0,
     &  -2.4977D+0, -1.9876D-05, -0.00032121D+0, 0.0022773D+0,
     &  -0.015262D+0, 0.081865D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_B'
!
!---  Convert mass fraction to molality and weight percent  ---
!
      TKX = TX + TABS
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
      YLSX = XLSX*1.D+2
!
!---  Enthalpy of mixing  ---
!
      DHMX = 0.D+0
      NC = 0
      DO 20 I = 0,3
        DO 10 J = 0,2
          NC = NC + 1
          DHMX = DHMX + SAX(NC)*(TX**I)*(GLSX**J)
   10   CONTINUE
   20 CONTINUE
      DHMX = (4.184D+3/(1.D+3+WTMS*GLSX))*DHMX
!
!---  Enthalpy of pure sodium chloride (halite)  ---
!
      HSX = 4.184D+3*(SBX(1)*(TX**3) + SBX(2)*(TX**2) +
     &  SBX(3)*TX)/WTMS
!
!---  Enthalpy of brine  ---
!
      HBX = (1.D+0-XLSX)*HLWX + XLSX*HSX + GLSX*DHMX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_L( TX,XLSX,XLCX,DHCX,HBX,HGAX,HLX )
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
!     This subroutine calculates the enthalpy of H2O-NaCl-CO2 solutions
!     as a function of brine enthalpy, gaseous C02 enthalpy,
!     temperature, mass fraction of dissolved NaCl, mass fraction
!     of dissolved CO2 and mass fraction of dissolved CH4.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas.  Geothermics, 26(4): 437-464.
!
!     Himmelblau, D. M.  1959.  Partial molal heats and entropies of
!     solution for gases dissolved in water from the freezing point
!     to near the critical point.  Journal of Physical Chemsitry,
!     63:1803-1808.
!
!     TX - temperature, C
!     XLSX - dissolved NaCl mass fraction
!     XLCX - dissolved gas component mass fraction
!     HBX - enthalpy of brine (NaCl-H2O), J/kg
!     HGAX - enthalpy of gas phase, J/kg
!     HLX - enthalpy of aqueous phase, J/kg
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DHCX(LNGC),HSAX(LNGC),XLCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_L'

      XLAX = 0.D+0
      DO IGC = 1, NGC
       XLAX = XLAX + XLCX(IGC)
      ENDDO
!
!---  Aqueous enthalpy  ---
!
        HLX = MAX(1.D+0-XLAX,0.D+0)*HBX 
        HLX = HLX + XLAX*HGAX
!
!---  Heat of solution for CO2 from Himmelblau (1959)  ---
!
      DO IGC = 1,NGC
        TKX = TX + TABS
        HSAX(IGC) = -RCU*(TKX**2)*DHCX(IGC)/GCPP(1,IGC)
        HLX = HLX + XLCX(IGC)*HSAX(IGC)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_L group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_S( TX,HSX )
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
!     Enthalpy of precipitated NaCl.
!
!     Lide, D.R. and H.V. Kehiaian.  1994.  CRC Handbook of
!     Thermophysical and Thermochemical Data, CRC Press, Inc.,
!     Boca Raton, Florida, pp. 97-98.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 CAX(5)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CAX
      DATA CAX / 25.19D+0, 0.1973D+0, -6.0114D-4, 8.81505D-7,
     &  -4.765D-10 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_S'
!
!---  0 C Reference  ---
!
      TKX = TX+TABS
      HSX = -1.24858D-4
      DO 10 I = 1,5
        REALX = REAL(I)
        HSX = HSX + CAX(I)*(TKX**I)/REALX
   10 CONTINUE
!
!---  Convert from J/mol to J/kg  ---
!
      HSX = 1.D+3*HSX/WTMS
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_S group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_W( TX,PX,HLX,HVX,ISRX )
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
!     Enthalpy and internal energy (J/kg) of pure water as a function
!     of temperature and pressure.
!
!     Pressure Range: 0 - 100 MPa (1000 Bar)
!     Temperature Range:  273.16 K (0.01 C) to 1073.15K (800 C)
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 CAX(23),SAX(12),CBX(31),SBX(5),LFCX(3),SLX
      REAL*8 TX,PX
      INTEGER ISRX,INX(8),IZX(8,3),ITX(8),IXX(8,2)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CAX,SAX,CBX,SBX
      SAVE LFCX,SLX,INX,IZX,ITX,IXX
      DATA CAX / 6.824687741D+3, -5.422063673D+2, -2.096666205D+4,
     &  3.941286787D+4, -6.733277739D+4, 9.902381028D+4,
     &  -1.093911774D+5, 8.590841667D+4, -4.511168742D+4,
     &  1.418138926D+4, -2.017271113D+3, 7.982692717D+0,
     &  -2.616571843D-2, 1.522411790D-3, 2.284279054D-2,
     &  2.421647003D+2, 1.269716088D-10, 2.074838328D-7,
     &  2.174020350D-8, 1.105710498D-9, 1.293441934D+1,
     &  1.308119072D-5, 6.047626338D-14 /
      DATA SAX / 8.438375405D-1, 5.362162162D-4, 1.720000000D+0,
     &  7.342278489D-2, 4.975858870D-2, 6.537154300D-1,
     &  1.150000000D-6, 1.150800000D-5, 1.418800000D-1,
     &  7.002753165D+0, 2.995284926D-4, 2.040000000D-1 /
      DATA CBX / 1.683599274D+1, 2.856067796D+1, -5.438923329D+1,
     &  4.330662834D-1, -6.547711697D-1, 8.565182058D-2,
     &  6.670375918D-2, 1.388983801D+0, 8.390104328D-2,
     &  2.614670893D-2, -3.373439453D-2, 4.520918904D-1,
     &  1.069036614D-1, -5.975336707D-1, -8.847535804D-2,
     &  5.958051609D-1, -5.159303373D-1, 2.075021122D-1,
     &  1.190610271D-1, -9.867174132D-2, 1.683998803D-1,
     &  -5.809438001D-2, 6.552390126D-3, 5.710218649D-4,
     &  1.936587558D+2, -1.388522425D+3, 4.126607219D+3,
     &  -6.508211677D+3, 5.745984054D+3, -2.693088365D+3,
     &  5.235718623D+2/
      DATA SBX / 7.633333333D-1, 4.006073948D-1, 8.636081627D-2,
     &  -8.532322921D-1, 3.460208861D-1 /
      DATA LFCX / 1.574373327D+1, -3.417061978D+1, 1.931380707D+1 /
      DATA SLX / 4.260321148D+0 /
      DATA INX / 2, 3, 2, 2, 3, 2, 2, 2 /
      DATA IZX / 13, 18, 18, 25, 32, 12, 24, 24,
     &  3, 2, 10, 14, 28, 11, 18, 14,
     &  0, 1, 0, 0, 24, 0, 0, 0 /
      DATA ITX / 0, 0, 0, 0, 0, 1, 1, 2 /
      DATA IXX / 0, 0, 0, 0, 0, 14, 19, 54,
     &  0, 0, 0, 0, 0, 0, 0, 27 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_W'
!
!---  Restrict temperature to above 0.01 C and below 800 C  ---
!
      TY = MIN( MAX( TX,1.D-2 ),8.D+2 )
!
!---  Restrict pressure to above 0.0 MPa and below 100 MPa  ---
!
      PY = MIN( MAX( PX,0.D+0 ),1.D+8 )
!
!---  Reduced temperature and pressure  ---
!
      THETAX = (TY+TABS)/TCRW
      BETAX = PY/PCRW
!
!---  L-function  ---
!
      BETALX = LFCX(1) + LFCX(2)*THETAX + LFCX(3)*(THETAX**2)
      BETALPX = LFCX(2) + 2.D+0*LFCX(3)*THETAX
!
!---  Subregions 1 or 6  ---
!
      IF( ISRX.EQ.1 .OR. ISRX.EQ.6 ) THEN
        CYX = 1.D+0 - SAX(1)*(THETAX**2) - SAX(2)/(THETAX**6)
        CYPX = -2.D+0*SAX(1)*THETAX + 6.D+0*SAX(2)/(THETAX**7)
        CZX = SAX(3)*(CYX**2)
        CZX = CZX - 2.D+0*SAX(4)*THETAX
        CZX = CZX + 2.D+0*SAX(5)*BETAX
        IF( (CZX/EPSL).GT.EPSL ) THEN
          CZX = CYX + SQRT(CZX)
        ELSE
          CZX = CYX
        ENDIF
        RHX = CAX(1)*THETAX
        DO 10 I = 1,10
          REALX = REAL(I)
          RHX = RHX - (REALX-2.D+0)*CAX(I+1)*(THETAX**(I-1))
   10   CONTINUE
        RHX = RHX + CAX(12)*(CZX*(1.7D+1*((CZX/2.9D+1)-(CYX/1.2D+1)) +
     &    5.D+0*THETAX*(CYPX/1.2D+1)) + SAX(4)*THETAX -
     &    (SAX(3)-1.D+0)*THETAX*CYX*CYPX)/(CZX**(5.D+0/17.D+0))
        RHX = RHX + (CAX(13) - CAX(15)*(THETAX**2) +
     &    CAX(16)*(9.D+0*THETAX + SAX(6))*((SAX(6)-THETAX)**9) +
     &    CAX(17)*(2.D+1*(THETAX**19) +
     &    SAX(7))/((SAX(7) + (THETAX**19))**2))*BETAX
        RHX = RHX - (1.2D+1*(THETAX**11) + SAX(8))*
     &    (CAX(18)*BETAX + CAX(19)*(BETAX**2) + CAX(20)*(BETAX**3))
     &    /((SAX(8) + (THETAX**11))**2)
        RHX = RHX + CAX(20)*(THETAX**18)*
     &    (1.7D+1*SAX(9) + 1.9D+1*(THETAX**2))*
     &    (((SAX(10) + BETAX)**3) + SAX(11)*BETAX)
        RHX = RHX + CAX(22)*SAX(12)*(BETAX**3)
        RHX = RHX + 2.1D+1*CAX(23)*(BETAX**4)/(THETAX**20)
        HLX = 1.D-3*RHX*PCRW*VCRW/WTMW
      ENDIF
!
!---  Below 0.01 C  ---
!
      IF( TX.LT.1.D-2 ) HLX = HLX - 4.202405*(1.D-2-TX)
!
!---  Subregions 2 or 6  ---
!
      IF( ISRX.EQ.2 .OR. ISRX.EQ.6 ) THEN
        CXX = EXP(SBX(1)*(1.D+0-THETAX))
        RHX = CBX(1)*THETAX
        ICBX = 1
        DO 20 I = 1,5
          ICBX = ICBX + 1
          REALX = REAL(I)
          RHX = RHX - CBX(ICBX)*(REALX-2.D+0)*(THETAX**(I-1))
   20   CONTINUE
        ICBX = 6
        DO 40 I = 1,5
          RHAX = 0.D+0
          DO 30 J = 1,INX(I)
            ICBX = ICBX + 1
            RHAX = RHAX + CBX(ICBX)*(1.D+0 + IZX(I,J)*SBX(1)*THETAX)*
     &        (CXX**IZX(I,J))
   30     CONTINUE
          RHX = RHX - (BETAX**I)*RHAX
   40   CONTINUE
        ICBX = 18
        DO 70 I = 6,8
          RHBX = 0.D+0
          RHCX = 0.D+0
          ISBX = 1
          DO 50 K = 1,ITX(I)
            ISBX = ISBX + 1
            RHBX = RHBX + IXX(I,K)*SBX(ISBX)*(CXX**IXX(I,K))
            RHCX = RHCX + SBX(ISBX)*(CXX**IXX(I,K))
   50     CONTINUE
          RHBX = RHBX*SBX(1)*THETAX
          RHCX = RHCX + (BETAX**(2-I))
          DO 60 J = 1,INX(I)
            ICBX = ICBX + 1
            RHAX = CBX(ICBX)*(CXX**IZX(I,J))
            RHAX = RHAX*((1.D+0 + IZX(I,J)*SBX(1)*THETAX) - RHBX/RHCX )
   60     CONTINUE
          RHX = RHX - SBX(1)*RHAX/RHCX
   70   CONTINUE
        RHAX = BETAX*((BETAX/BETALX)**10)
        RHBX = 1.D+1*BETALPX/BETALX
        ICBX = 24
        DO 80 I = 0,6
          ICBX = ICBX + 1
          REALX = REAL(I)
          RHX = RHX + RHAX*(1.D+0 + THETAX*(RHBX + REALX*SBX(1))*
     &      CBX(ICBX)*(CXX**I))
   80   CONTINUE
        HVX = 1.D-3*RHX*PCRW*VCRW/WTMW
      ENDIF
!
!---  Subregions 3 or 5  ---
!
      IF( ISRX.EQ.3 .OR. ISRX.EQ.5 ) THEN
        INDX = 14
        CHMSG = 'Steam Table Subregion 3/5: Exceeds Temperature Range: '
        RLMSG = TX
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Subregions 4 or 5  ---
!
      IF( ISRX.EQ.4 .OR. ISRX.EQ.5 ) THEN
        INDX = 14
        CHMSG = 'Steam Table Subregion 4/5: Exceeds Temperature Range: '
        RLMSG = TX
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTHG_I( HGX,HRGOX,PVAX,PVWX,QMX,QTHX,RHOGX,
     &  THKGX,TPX,TX,VOLW,XGCX,XGWX,XMGCX )
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
!     H2O-NaCl-CO2 Equation of state.  This subroutine finds the
!     temperature at a given gas enthalpy, CO2 partial pressure,
!     and water partial pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
          REAL*8 XMGCX(LNGC),XGCX(LNGC),HGCX(LNGC),UEGCX(LNGC)
          REAL*8 HGCZ(LNGC),UEGCZ(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTHG_I'
      ISRX = 2
!
!---  Newton scheme  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
      CALL DENS_GC( TX,PVAX,XMGCX,RHOGCX,ISRX,ZGX )
      RHOGX = RHOGWX+RHOGCX
      XGAX = RHOGCX/RHOGX
      XGWX = RHOGWX/RHOGX
      WTMX = 0.D+0
      DO IGC = 1, NGC
         WTMX = XMGCX(IGC)*GCPP(1,IGC)
      ENDDO
      PVX = PVAX
        CALL ENTH_GC( TX,PVX,XMGCX,XMGWX,WTMX,HCGX,UCGX,
     &      ISRX )
      CALL ENTH_W( TX,PVWX,HLWX,HGWX,ISRX )
      HGY = XGAX*HCGX + XGWX*HGWX
      RY = QMX*(HGY-HCGX) + (HGY*RHOGX-HRGOX)*VOLW*DTI + (TX-TPX)*THKGX
      TZ = TX + 1.D-7
      CALL DENS_W( TZ,PVWX,RHOLWZ,RHOGWZ,ISRX )
      CALL DENS_GC( TZ,PVAX,XMGCX,RHOGCZ,ISRX,ZGX )
      RHOGZ = RHOGWZ+RHOGCZ
      XGCZ = RHOGCZ/RHOGZ
      XGWZ = RHOGWZ/RHOGZ
      WTMX = 0.D+0
      DO IGC = 1, NGC
       WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
      ENDDO
      PVX = PVAX
        CALL ENTH_GC( TZ,PVX,XMGCX,XMGWX,WTMX,HCGZ,UCGZ,ISRX)
      CALL ENTH_W( TZ,PVWX,HLWZ,HGWZ,ISRX )
      HGZ = XGAZ*HCGZ + XGWZ*HGWZ
      RZ = QMX*(HGZ-HCGZ) + (HGZ*RHOGZ-HRGOX)*VOLW*DTI + (TZ-TPX)*THKGX
      FX = RY
      DFX = 1.D+7*(RZ-RY)
      DTX = -FX/DFX
      SFX = 2.D+0/3.D+0
      SUMXMGCX = 0.D+0
      DO IGC = 1,NGC
        SUMXMGCX = SUMXMGCX + XMGCX(IGC)
      ENDDO
      IF( NGC.EQ.1 .OR. SUMXMGCX.EQ.0.D+0 ) THEN
        TCIJ = GCPP(3,NGC)
      ELSE
      TCIJ = 0.D+0
      CR1MX = 0.D+0
      DO IGC = 1, NGC
       CR1MX = CR1MX + XMGCX(IGC)*GCPP(5,IGC)
       TCIJ = TCIJ + XMGCX(IGC)*GCPP(5,IGC)*GCPP(3,IGC)
      ENDDO
      TCIJ = TCIJ/CR1MX
      ENDIF
      IF( ABS(TX+TABS-TCIJ).LT.2.D+0 ) THEN
        DTX = SIGN( MIN(ABS(DTX),(ABS(TX+TABS-TCIJ)/8.D-1)),DTX )
      ELSE
        DTX = 8.D-1*DTX
      ENDIF
      TX = TX + DTX
      IF( NC.GT.32 ) THEN
        INDX = 14
        RLMSG = TX
        CHMSG = 'Unconverged Well Temperature'
        CALL WRMSGS( INDX )
      ENDIF
      IF( ABS(DTX).GT.1.D-6 ) GOTO 100
!
!---  Use the inlet enthalpy and outlet temperature to calculate 
!     energy transport into the formation  ---
!
      QTHX = QMX*HGX + (TX-TPX)*THKGX
      CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
      CALL DENS_GC( TX,PVAX,XMGCX,RHOGAX,ISRX,ZGX )
      RHOGX = RHOGWX+RHOGAX
      WTMX = 0.D+0
      DO IGC = 1, NGC
       WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
      ENDDO
      PVX = PVAX
        CALL ENTH_GC( TX,PVX,XMGCX,XMGWX,WTMX,HCGX,UECGX,
     &        ISRX )
      CALL ENTH_W( TX,PVWX,HLWX,HGWX,ISRX )
      HGX =  XGWX*HGWX + XGAX*HCGX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTHG_I group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XMLWX,ZGX,
     &  WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
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
!     This subroutine calculates the equilibrium state between
!     a multiple component gas mixture and water
!
!----------------------Authors-----------------------------------------!
!
!     Written by R Ramanathan, PNNL, 11 June 2012.
!     Modified by DH Bacon, PNNL, 31 October 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE NAPL
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (MAXSSIT=200)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 CLMBX(3),CXIX(3)
      REAL*8 FUGX(LNGC+1),FULX(LNGC+1)
      REAL*8 BIJ(LNGC+1,LNGC+1),YT(LNGC+1),XT(LNGC+1),AXI(LNGC+1,LNGC+1)
      REAL*8 OMEGAX, OMEGAXW, KS(LNGC+1),TRX(LNGC+1),SAX(LNGC+1)
      REAL*8 XTEMP(LNGC+1),BIJQ(LNGC+1,LNGC+1)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),PGCX(LNGC)
      REAL*8 YX(LNGC+1),XX(LNGC+1)
      REAL*8 XGCX(LNGC),XLCX(LNGC),AIX(LNGC+1,LNGC+1)
      REAL*8 XMDGCX(LNGC),XMDLCX(LNGC)
      REAL*8 ERR1X(MAXSSIT),ERR2X(MAXSSIT)
      INTEGER JNCG,ICTX
!
!----------------------Data Statements---------------------------------!
!
      SAVE CLMBX,CXIX
      DATA CLMBX / 2.217D-4,1.074D+0,2.648D+3 /
      DATA CXIX / 1.3D-5,-2.012D+1,5.259D+3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/EQUIL'

      INCG = 1
      ERR1X = 0.D+0
      ERR2X = 0.D+0
!
!---  Temperature in Kelvin, pressure in bar, 
!     gas constant in cm^3 bar/K mol   ---
!
      TKX = TX+TABS
      PBX = PX*1.D-5
      RCUBX = RCU*1.D-2
      TLX = 99.D+0
      THX = 109.D+0
      TKLX = TLX+TABS
      TKHX = THX+TABS
      TCRA = 304.1D+0
      PCRA = 73.8D+0
      TCRAW = 647.14D+0
      PCRAW = 220.64D+0
      OMEGAX = 0.239D+0
      OMEGAXW = 0.344D+0
!
!---  Salt effects   ---
!
      XIX = CXIX(1)*TKX + CXIX(2)/TKX + CXIX(3)/(TKX**2.D+0)
      CLX = CLMBX(1)*TKX + CLMBX(2)/TKX + CLMBX(3)/(TKX**2.D+0)
!
!---  Convert NaCl concentration to molality   ---
!
      XMOLSX = 1.D+3*(XLSX/(1.D+0-XLSX))/WTMS
      XMOLWX = 1.D+3/WTMW
      XMOLNAX = XMOLSX
      XMOLCLX = XMOLSX
      APCX = (1.D+0+(XMOLNAX+XMOLCLX)/XMOLWX)*
     &  EXP(2.D+0*CLX*XMOLNAX + XIX*XMOLCLX*XMOLNAX)
      CWX = XMOLSX
!
!---  Peng-Robinson pure component parameters    ---
!
      DO IGC = 1, NGC
        TRX(IGC) = TKX/GCPP(3,IGC)
        GCPP(15,IGC) = 0.37464D+0 + 1.54226D+0*GCPP(7,IGC) -
     &    0.26992D+0*(GCPP(7,IGC)**2.D+0)
        SAX(IGC) =  (1.D+0 + GCPP(15,IGC)*
     &    (1.D+0-SQRT(TRX(IGC))))**2.D+0
      ENDDO
!
!---  Water
!
      TRWX = TKX/TCRW
      FWW =  0.37464D+0 + 1.54226D+0*OMEGAXW -
     &  0.26992D+0*(OMEGAXW**2.D+0)
      SAWX = ((1.D+0 + FWW*(1.D+0-SQRT(TRWX)))**2.D+0)
      CALL SP_W( TX,PSWX )
      YX(NGC+1) = PSWX/PX
      YXT = MAX( MIN( 1.D+0-YX(NGC+1),1.D+0 ),0.D+0 )
      DO IGC = 1, NGC
        YX(IGC) = XMDGCX(IGC)*(1.D+0-YX(NGC+1))
      ENDDO
      DO IGC = 1, NGC
        XX(IGC) = 1.D-2
      ENDDO
      XX(NGC+1) = 0.D+0
      XX(NGC+1) = 1.D+0 - SUM(XX)
!
!---  Succesive Substitution  ---
!
      ICTX = 0
  111 CONTINUE
      ICTX = ICTX + 1
      IF(ICTX.GT.MAXSSIT) THEN
        INDX = 0
        CHMSG = 'Mutual Solubility Convergence Failure: ' //
     &    'Successive-Substitution '
        CALL WRMSGS( INDX )
        GO TO 800
      ENDIF
!
!---  Non-aqueous Binary Interaction Coefficients
!
      DO IGC = 1, NGC+1
        BIJ(IGC,IGC) = 0.D+0
      ENDDO
      DO IGC = 1, NGC
        DO JGC = IGC, NGC
          BIJ(IGC,JGC) = BIJNA(IGC,JGC)+BIJTD(IGC,JGC)*TRX(IGC)
          BIJ(JGC,IGC)=BIJ(IGC,JGC)
        ENDDO
      ENDDO
      DO IGC = 1, NGC
        BIJ(IGC,NGC+1) = BIJNA(IGC,NGC+1)+BIJTD(IGC,NGC+1)*TRX(IGC)
        BIJ(NGC+1,IGC) = BIJ(IGC,NGC+1)
      ENDDO
!
!---  Peng-Robinson temperature dependent parameters
!
      DO IGC = 1, NGC
        GCPP(13,IGC) = SAX(IGC)*0.45724D+0*(RCUBX**2.D+0)
     &    *(GCPP(3,IGC)**2.D+0/(GCPP(4,IGC)*1.D-5))*PBX/
     &    ((RCUBX*TKX)**2.D+0)
        GCPP(14,IGC) = (0.07780D+0*RCUBX*GCPP(3,IGC)/
     &    (GCPP(4,IGC)*1.D-5))*PBX/(RCUBX*TKX)
      ENDDO
      AWX = SAWX*0.45724D+0*(RCUBX**2.D+0)*(TCRAW**2.D+0/PCRAW)*PBX/
     &  ((RCUBX*TKX)**2.D+0)
      BWX = (0.07780D+0*RCUBX*TCRAW/PCRAW)*PBX/(RCUBX*TKX)
!
!---  Peng-Robinson concentration dependent parameters
!
      SAMX = 0.D+0
      SBMX = 0.D+0
      DO IGC = 1, NGC
        DO JGC = 1, NGC
          SAMX = SAMX + YX(IGC)*YX(JGC)*SQRT(GCPP(13,IGC)*
     &      GCPP(13,JGC))*(1.D+0 - BIJ(IGC,JGC))
          AIX(IGC,JGC) = SQRT(GCPP(13,IGC)*
     &      GCPP(13,JGC))*(1.D+0 - BIJ(IGC,JGC))
        ENDDO
      ENDDO
      DO IGC = 1, NGC
        SBMX = SBMX + YX(IGC)*GCPP(14,IGC)
      ENDDO
!
!---  Add contribution of water vapor ---
!
      DO IGC = 1, NGC
        SAMX = SAMX + 2.D+0*(YX(IGC)*YX(NGC+1)*SQRT(GCPP(13,IGC)*
     &    AWX)*(1.D+0-BIJ(IGC,NGC+1)))
        AIX(IGC,NGC+1) = SQRT(GCPP(13,IGC)*AWX)*(1.D+0-BIJ(IGC,NGC+1))
        AIX(NGC+1,IGC) = AIX(IGC,NGC+1)
      ENDDO
      AIX(NGC+1,NGC+1) = SQRT(AWX*AWX)*(1.D+0-BIJ(NGC+1,NGC+1))
      SAMX = SAMX + YX(NGC+1)*YX(NGC+1)*SQRT(AWX*AWX)*
     &  (1.D+0-BIJ(NGC+1,NGC+1))
      SBMX = SBMX + YX(NGC+1)*BWX
! 
!---  Polynomial coefficients for Peng-Robinson equation of state  ---
!
      CAX = 1.D+0
      CBX = -(1.D+0 - SBMX)
      CCX = SAMX - 2.D+0*SBMX - 3.D+0*(SBMX**2.D+0)
      CDX = -(SAMX*SBMX - (SBMX**2.D+0) - (SBMX**3.D+0))
      CALL NICKALLS( CAX,CBX,CCX,CDX,Z1X,Z2X,Z3X )
!
!---  Compressibility factor Z---
!
      ZGX = MAX( Z1X,Z2X,Z3X )
      ZNX = MIN( Z1X,Z2X,Z3X )
!
!---  Molar volume, cm^3/mol  ---
!
      VGX = RCUBX*TKX*ZGX/PBX
      VNX = RCUBX*TKX*ZNX/PBX
!
!---  Gas conditions  ---
!
      SFX = 2.D+0/3.D+0
      SUMYX = 0.D+0
      DO IGC = 1,NGC
        SUMYX = SUMYX + YX(IGC)
      ENDDO
      IF( NGC.EQ.1 .OR. SUMYX.EQ.0.D+0 ) THEN
        TCIJ = GCPP(3,NGC)
      ELSE
      TCIJ = 0.D+0
      CR1MX = 0.D+0
      DO IGC = 1, NGC
        CR1MX = CR1MX + YX(IGC)*GCPP(5,IGC)
        TCIJ = TCIJ + YX(IGC)*GCPP(5,IGC)*GCPP(3,IGC)
      ENDDO
      TCIJ = TCIJ/CR1MX
      ENDIF
      DELTA1 = 1.D+0+SQRT(2.D+0)
      DELTA2 = 1.D+0-SQRT(2.D+0)
      GSTAR = (ZNX-ZGX)+LOG(ABS(((ZGX-SBMX)/(ZNX-SBMX))))+
     +  (1.D+0/(DELTA2-DELTA1))*(SAMX/SBMX)
     +  *LOG(ABS(((ZGX+DELTA2*SBMX)*(ZNX+DELTA1*SBMX)/
     +  ((ZNX+DELTA2*SBMX)*(ZGX+DELTA1*SBMX)))))

      IF( TKX.GE.TCIJ ) THEN
        VX = VGX
        ZX = ZGX
      ELSE
        IF(GSTAR.GE.0.D+0) THEN
          VX = VGX
          ZX = ZGX
          GOTO 902
        ELSE
          VX = VNX
          ZX = ZNX
          GOTO 902
        ENDIF
!
!---    Gas phase condition  ---
!
        IF( (W2X-W1X)/EPSL.LE.EPSL ) THEN
          VX = VGX
          ZX = ZGX
!
!---      Liquid phase condition  ---
!
        ELSEIF( (W1X-W2X)/EPSL.GT.EPSL ) THEN
          VX = VNX
          ZX = ZNX
!
!---    Gas-liquid phase condition, use gas molar volume  ---
!
        ELSE
          VX = VGX
          ZX = ZGX
        ENDIF
      ENDIF
  902 CONTINUE
      ZVX = ZX
!
!---  Fugacity coefficients  ---
!
      D1X = SQRT(2.D+0) + 1.D+0
      D2X = 1.D+0 - SQRT(2.D+0)
      DO 120 KGC = 1,NGC+1
        SUMX = 0.D+0
        DO 110 IGC = 1,NGC+1
          SUMX = SUMX + YX(IGC)*AIX(IGC,KGC)
  110   CONTINUE
        IF (KGC.LE.NGC)BXX=GCPP(14,KGC)
        IF (KGC.GT.NGC)BXX=BWX
        DLTA3 = SAMX/SBMX/SQRT(8.D+0)
        DLTA4 = LOG(ABS((ZVX+D1X*SBMX)/
     &    (ZVX+D2X*SBMX)))
        FUGX(KGC) = (BXX*(ZVX-1.D+0)/SBMX) - LOG(ABS(ZVX-SBMX))
        FUGX(KGC) = FUGX(KGC) - DLTA3*DLTA4*((2.D+0*SUMX/SAMX) -
     &    (BXX/SBMX))
        FUGX(KGC) = EXP(FUGX(KGC))
  120 CONTINUE
      JNCG = NGC+1
      YT(1:JNCG) = YX(1:JNCG)
!
!---  Aqueous Binary Interaction Coefficients  ---
!
      JNCG = NGC+1
      DO IGC = 1,NGC
        TRX(IGC) = TKX/GCPP(3,IGC)
      ENDDO
      BIJQ(1:JNCG,1:JNCG) = 0.0D+0
      DO IGC = 1, NGC
        DO JGC = 1, NGC
          BIJQ(IGC,JGC) = 0.D+0
        ENDDO
      ENDDO
      DO IGC = 1, NGC
        CWX = XMOLSX
        BIJQ(IGC,NGC+1) = 
     &  GCPP(37,IGC)*(1.d+0+GCPP(41,IGC)*CWX**GCPP(45,IGC))
     & +GCPP(38,IGC)*(1.d+0+GCPP(42,IGC)*CWX**GCPP(46,IGC))*TRX(IGC)
     & +GCPP(39,IGC)*(1.d+0+GCPP(43,IGC)*CWX**GCPP(47,IGC))
     &  *TRX(IGC)**2.d+0
     & +GCPP(40,IGC)*exp(GCPP(44,IGC)*CWX+GCPP(48,IGC)*TRX(IGC))
        BIJQ(NGC+1,IGC) = BIJQ(IGC,NGC+1)
        BIJQ(NGC+1,IGC) = BIJQ(IGC,NGC+1)
      ENDDO
!
!---  Water
!
      BIJQ(NGC+1,NGC+1) = 0.0D+0
      TRX(NGC+1) = TKX/TCRAW
      FWW  = (1.5201D+0 + TRX(NGC+1)*(-6.0350D-1 + 5.0273D-3*
     &  (XMOLSX**1.1D+0) + 8.4155D-2*TRX(NGC+1) + 1.2560D-3*
     &  (TRX(NGC+1)**(-4.D+0))) - 9.8324D-6*XMOLSX*(6.D+0-XMOLSX) - 
     &  1.2584D-4*XMOLSX*(1.D+0 + 1.1D+0*(TRX(NGC+1)**(-1.256D-3)))*
     &  (1.D+0 + (3.3057D+0*(TRX(NGC+1)**2.D+0) - 8.1723D+0*
     &  TRX(NGC+1) + 4.1454D+0)*(-0.107D+0*(XMOLSX**2.D+0) 
     &  + 0.8027D+0*XMOLSX - 1.0D+0)))**2.D+0

      AWX = FWW*0.45724D+0*(RCUBX**2.D+0)*(TCRAW**2.D+0/PCRAW)*PBX/
     &  ((RCUBX*TKX)**2.D+0)
      DO IGC = 1, NGC
        GCPP(13,IGC) = SAX(IGC)*0.45724D+0*(RCUBX**2.D+0)*
     &    (GCPP(3,IGC)**2.D+0/(GCPP(4,IGC)*1.D-5))*PBX
     &    /((RCUBX*TKX)**2.D+0)
      ENDDO
          
      SAMXQ = 0.D+0
      SBMXQ = 0.D+0
      DO IGC = 1, NGC
        DO JGC = 1, NGC
          SAMXQ = SAMXQ + XX(IGC)*XX(JGC)*SQRT(GCPP(13,IGC)*
     &      GCPP(13,JGC))*(1.D+0 - BIJQ(IGC,JGC))
          AIX(IGC,JGC) = SQRT(GCPP(13,IGC)*
     &      GCPP(13,JGC))*(1.D+0 - BIJQ(IGC,JGC))
        ENDDO
      ENDDO
!
!--- Add water contributions
!          
      DO IGC = 1, NGC
        SAMXQ = SAMXQ + 2.D+0*XX(IGC)*XX(NGC+1)*
     &    SQRT(AWX*GCPP(13,IGC))*(1.D+0-BIJQ(IGC,NGC+1))
        AIX(IGC,NGC+1) = SQRT(GCPP(13,IGC)*
     &    AWX)*(1.D+0 - BIJQ(IGC,NGC+1))
        AIX(NGC+1,IGC) = AIX(IGC,NGC+1)
      ENDDO
      AIX(NGC+1,NGC+1) = SQRT(AWX*AWX)*(1.D+0-BIJQ(NGC+1,NGC+1))
      SAMXQ = SAMXQ + XX(NGC+1)*XX(NGC+1)*SQRT(AWX*AWX)*
     &  (1.D+0-BIJQ(NGC+1,NGC+1))
      DO IGC = 1, NGC
        SBMXQ = SBMXQ + XX(IGC)*GCPP(14,IGC)
      ENDDO
      SBMXQ = SBMXQ + XX(NGC+1)*BWX
!
!---  Polynomial coefficients for Peng-Robinson equation of state  ---
!
      CAX = 1.D+0
      CBX = -(1.D+0 - SBMXQ)
      CCX = SAMXQ - 2.D+0*SBMXQ - 3.D+0*(SBMXQ**2.D+0)
      CDX = -(SAMXQ*SBMXQ - (SBMXQ**2.D+0) - (SBMXQ**3.D+0))
      CALL NICKALLS( CAX,CBX,CCX,CDX,Z1X,Z2X,Z3X )
!
!---  Compressibility factor Z---
!
      ZGX = MAX( Z1X,Z2X,Z3X )
      ZNX = MIN( Z1X,Z2X,Z3X )
!
!---  Molar volume, cm^3/mol  ---
!
      VGX = RCUBX*TKX*ZGX/PBX
      VNX = RCUBX*TKX*ZNX/PBX
!
!---  Gas conditions  ---
!
      SFX = 2.D+0/3.D+0
      SUMYX = 0.D+0
      DO IGC = 1,NGC
        SUMYX = SUMYX + YX(IGC)
      ENDDO
      IF( NGC.EQ.1 .OR. SUMYX.LT.EPSL ) THEN
        TCIJ = GCPP(3,NGC)
      ELSE
      TCIJ = 0.D+0
      CR1MX = 0.D+0
      DO IGC = 1, NGC
        CR1MX = CR1MX + YX(IGC)*GCPP(5,IGC)
        TCIJ = TCIJ + YX(IGC)*GCPP(5,IGC)*GCPP(3,IGC)
      ENDDO
        TCIJ = TCIJ/CR1MX
      ENDIF
!diana
      DELTA1 = 1.D+0+SQRT(2.D+0)
      DELTA2 = 1.D+0-SQRT(2.D+0)
      GSTAR = (ZNX-ZGX)+LOG(ABS(((ZGX-SBMX)/(ZNX-SBMX))))+
     +  (1.D+0/(DELTA2-DELTA1))*(SAMX/SBMX)
     +  *LOG(ABS(((ZGX+DELTA2*SBMX)*(ZNX+DELTA1*SBMX)/
     +  ((ZNX+DELTA2*SBMX)*(ZGX+DELTA1*SBMX)))))

      IF( TKX.GE.TCIJ ) THEN
        VX = VGX
        ZX = ZGX
      ELSE
        W1X = PBX*(VGX-VNX)
        W2X = RCUBX*TKX*LOG((VGX-SBMXQ)/(VNX-SBMXQ)) + (SAMXQ/
     &    (SQRT(TKX)*SBMXQ))*LOG(((VGX+SBMXQ)*VNX)/(VNX+SBMXQ)*VGX)
        IF(GSTAR.GE.0.D+0) THEN
          VX = VGX
          ZX = ZGX
          GOTO 903
        ELSE
          VX = VNX
          ZX = ZNX
          GOTO 903
        ENDIF
!
!---    Gas phase condition  ---
!
        IF( (W2X-W1X)/EPSL.LE.EPSL ) THEN
          VX = VGX
          ZX = ZGX
!---    Liquid phase condition  ---
!
        ELSEIF( (W1X-W2X)/EPSL.LT.EPSL ) THEN
          VX = VNX
          ZX = ZNX
!
!---    Gas-liquid phase condition, use gas molar volume  ---
!
        ELSE
          VX = VGX
          ZX = ZGX
        ENDIF
      ENDIF
  903 CONTINUE
      ZX = ZNX
      ZLX = ZX
!
!---  Fugacity coefficients  ---
!
      D1X = SQRT(2.D+0) + 1.D+0
      D2X = 1.D+0 - SQRT(2.D+0)
      DO 140 KGC = 1,NGC+1
        SUMX = 0.D+0
        IF (KGC.LE.NGC)AXX=GCPP(13,KGC)
        IF (KGC.GT.NGC)AXX=AWX
        DO 130 IGC = 1,NGC+1
          SUMX = SUMX + XX(IGC)*AIX(IGC,KGC)
  130   CONTINUE
        IF (KGC.LE.NGC)BXX=GCPP(14,KGC)
        IF (KGC.GT.NGC)BXX=BWX
        DLTA3 = SAMXQ/SBMXQ/SQRT(8.D+0)
        DLTA4 = LOG((ZLX+D1X*SBMXQ)/
     &    (ZLX+D2X*SBMXQ))
        FULX(KGC) = (BXX*(ZLX-1.D+0)/SBMXQ) - LOG(ABS(ZLX-SBMXQ))
        FULX(KGC) = FULX(KGC) - DLTA3*DLTA4*((2.D+0*SUMX/SAMXQ) -
     &    (BXX/SBMXQ))
        FULX(KGC) = EXP(FULX(KGC))
  140 CONTINUE
!
!---  Update mole fractions  ---
!
      JNCG = NGC+1
      XT(1:JNCG) = XX(1:JNCG)
      YT(1:JNCG) = YX(1:JNCG)
      DO IGC = 1, JNCG-1
        XX(IGC) = YX(IGC)*FUGX(IGC)/FULX(IGC)
      ENDDO
      XX(JNCG) = 0.D+0
      XX(JNCG) = 1.D+0 - SUM(XX)
      DO IGC = JNCG,JNCG
        YX(IGC) = XX(IGC)*FULX(IGC)/FUGX(IGC)
      ENDDO
      YX(NGC) = 0.D+0
      YX(NGC) = 1.D+0 - SUM(YX)
      DO IGC = 1, NGC
        IF( YX(IGC).LT.0.D+0 ) YX(IGC)=0.D+0
      ENDDO          
!
!---  Convergence check
!
      DO I = 1, JNCG
        ERR1X(ICTX) = ABS(XT(I) - XX(I))
        ERR2X(ICTX) = ABS(YT(I) - YX(I))
        IF( (ERR1X(ICTX).GT.1.0D-6) .OR. (ERR2X(ICTX).GT.1.0D-6) ) THEN
!
!---  Divergence check
!
          IF( (ICTX.GT.3) ) THEN
            IF( ERR1X(ICTX).GT.ERR1X(ICTX-1) .AND.
     &          ERR1X(ICTX).GT.ERR1X(ICTX-2) .AND.
     &          ERR1X(ICTX).GT.ERR1X(ICTX-3) ) GO TO 800 
          ENDIF
          GOTO 111
        ENDIF
      ENDDO
  800 CONTINUE
!
!---  Equilibrium mole fractions  ---
!
      XMGWX = YX(NGC+1)
      DO IGC = 1,NGC
        XMGCX(IGC) = YX(IGC)
      ENDDO
      XMGCXM = SUM(YX) - YX(NGC+1)
      XMLWX = XX(NGC+1)
      DO IGC = 1, NGC
        XMLCX(IGC) = XX(IGC)
      ENDDO
      XMLCXM = SUM(XX) - XX(NGC+1)
!
!---Molecular weight of mixture without water
!
      WTMX = 0.D+0
      WTMGCX = 0.D+0
      DO IGC = 1, NGC
        IF(NGC.NE.1) THEN
          WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
          WTMGCX = WTMGCX + XMGCX(IGC)*GCPP(1,IGC)
        ELSE
          WTMX = WTMX + GCPP(1,IGC)
          WTMGCX = WTMGCX + GCPP(1,IGC)
        ENDIF
      ENDDO

      IF( XLSX.GT.EPSL ) THEN
        XMLWX = (WTMS*(1.D+0-XLSX) + XLSX*XMLCXM*(WTMS-WTMX)
     &    - XMLCXM*WTMS)/(WTMS*(1.D+0-XLSX) + XLSX*WTMW)
        XMLSX = 1.D+0 - XMLCXM - XMLWX
      ELSE
        XMLWX = 1.D+0 - XMLCXM
        XMLSX = 0.D+0
      ENDIF
!
!---  Water vapor pressure lowering, reduce water mole fraction in 
!     the CO2-rich phase by the ratio of the reduced to 
!     saturated water vapor pressure  ---
!
      XMGWX = XMGWX*PVBX/PSBX
      IF( XMGCXM.NE.0.D+0 ) THEN
      DO IGC = 1, NGC
        XMGCX(IGC) = (1.D+0-XMGWX)*XMGCX(IGC)/XMGCXM
      ENDDO
      ENDIF
!
!---  Limit mole fractions  ---
!
      IF( XMLWX.LT.EPSL ) XMLWX = 0.D+0
      IF( XMGWX.LT.EPSL ) XMGWX = 0.D+0
      DO IGC = 1, NGC
        IF( XMGCX(IGC).LT.EPSL ) XMGCX(IGC) = 0.D+0
        IF( XMLCX(IGC).LT.EPSL ) XMLCX(IGC) = 0.D+0
      ENDDO
!
!---  Gas partial pressures, Pa  ---
!
      PGAX = 0.D+0
      DO IGC = 1,NGC
        PGCX(IGC) = XMGCX(IGC)*PX
        PGAX = PGAX + PGCX(IGC)
      ENDDO
      PGWX = XMGWX*PX
!
!---  Equilibrium mass fractions  ---
!
      WTMX = 0.D+0
      DO IGC = 1, NGC
        WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
      ENDDO
      WTMX = WTMX + XMGWX*WTMW 
      XGWX = XMGWX*WTMW/WTMX
      DO IGC = 1,NGC
        XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/WTMX
      ENDDO

      WTMCL = 0.D+0
      DO IGC = 1,NGC
       WTMCL = WTMCL + XMLCX(IGC)*GCPP(1,IGC)
      ENDDO
      
      WTMLX = XMLWX*WTMW + WTMCL + XMLSX*WTMS
      XLWX = XMLWX*WTMW/WTMLX
      DO IGC = 1, NGC
        XLCX(IGC) = XMLCX(IGC)*GCPP(1,IGC)/WTMLX
      ENDDO

      ZGX = ZVX
      ZNX = ZLX
!
!--- Update dry gas mole fractions ---
!
      DO IGC = 1, NGC
       XMDGCX(IGC) = XMGCX(IGC)/SUM(XMGCX)
       XMDLCX(IGC) = XMLCX(IGC)
      ENDDO
      XMDGWX = XMGWX
      XMDLWX = XMLWX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of EQUIL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_11( TX,PX,PCX,PGX,RHOLSX,RHOLCX,YLSX,XLSX,XLAX,
     &  CHMSGX,XMDGCX )
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
!     Compute dissolved-salt and -CO2 aqueous mass fractions from
!     dissolved-salt and -CO2 aqueous concentrations
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 YLSY(4),XLAY(4),XLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2),DXLCX(LNGC),XLCX(LNGC),XLCYX(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC)
      REAL*8 XLCMX(LNGC),RHOLCX(LNGC),XLCY(LNGC,LSV),XMDGCX(LNGC),
     &XMDLCX(LNGC),XGCX(LNGC)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_11'
!
!---  Pure water density  ---
!
      ISRX = 1
      CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
!
!---  Guess dissolved-salt and -co2 aqueous mass fractions  ---
!
      YLSX = RHOLSX/(RHOLWX+RHOLSX+RHOLAX)
      XLAX = RHOLAX/(RHOLWX+RHOLSX+RHOLAX)
      XLSX = 0.D+0
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,
     &  XLWX,ZGX,WTGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
      CALL SOL_LS( TX,XLSMX )
      DYLSX = 1.D-4*XLSMX
      DXLAX = 0.D+0
      DO IGC = 1, NGC
       DXLCX(IGC) = 1.D-4*XLCMX(IGC)
      ENDDO
      NC = 0
!
!---  Two-variable Newton-Raphson loop (YLS, XLA)  ---
!
  10  CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = CHMSGX(1)
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 2,4
        YLSY(M) = YLSX
         IF( M.EQ.3 ) YLSY(M) = YLSX + DYLSX
        DO IGC = 1, NGC
         XLCY(IGC,M) = XLCX(IGC)
         IF( M.EQ.4 ) XLCY(IGC,M) = XLCX(IGC) + DXLCX(IGC)
        ENDDO
        XLSY(M) = MIN( YLSY(M),XLSMX )
        CALL SP_B( TX,XLSY(M),PSBX )
        PVBX = PSBX
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,
     &  XLWX,ZGX,WTGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
        XLAX = 0.D+0
        DO IGC = 1, NGC
          XLAX = XLAX + XLCX(IGC)
        ENDDO
        XLAY(M) = 0.D+0
        DO IGC = 1, NGC
         XLAY(M) = XLAY(M) + XLCX(IGC)
        ENDDO
        PVAX = MIN( XLAY(M)/XLAX,1.D+0 )*PGAX
        IF( PSBX+PVAX .GT. PGX .AND. M.EQ.2 ) THEN
          INDX = 9
          RLMSG = PSBX+PVAX
          CHMSG = CHMSGX(2)
          CALL WRMSGS( INDX )
        ENDIF
        CALL DENS_B( TX,PX,XLSY(M),RHOBX )
        DO IGC = 1, NGC
         XLCYX(IGC) = XLCY(IGC,M)
        ENDDO
          CALL DENS_L( TX,RHOBX,XLCYX,RHOLX )
        GX(M,1) = XLSY(M) - RHOLSX/RHOLX
        GX(M,2) = XLAY(M) - RHOLAX/RHOLX
   20 CONTINUE
      RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
      RX(1,2) = (GX(4,1)-GX(2,1))/DXLAX
      RX(2,1) = (GX(3,2)-GX(2,2))/DYLSX
      RX(2,2) = (GX(4,2)-GX(2,2))/DXLAX
      RPX(1) = -GX(2,1)
      RPX(2) = -GX(2,2)
      CYLSX = (RPX(2)-RPX(1)*RX(2,2)/(RX(1,2)+SMALL))/
     &  (RX(2,1)-RX(1,1)*RX(2,2)/(RX(1,2)+SMALL))
        CXLAX = (RPX(2)-RPX(1)*RX(2,1)/(RX(1,1)+SMALL))/
     &  (RX(2,2)-RX(1,2)*RX(2,1)/(RX(1,1)+SMALL))
      YLSX = YLSX + CYLSX
      XLAX = XLAX + CXLAX
      IF( ABS(CYLSX).GT.(1.D-6*XLSMX) .OR.
     &  ABS(CXLAX).GT.(1.D-6*XLAMX) ) GOTO 10
      XLSX = MIN( YLSX,XLSMX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_11 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_12( TX,PX,PCX,PGX,RHOLSX,PHILAX,YLSX,XLSX,XLAX,
     &  CHMSGX,XMDGCX )
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
!     Compute dissolved-salt and -CO2 aqueous mass fractions from
!     dissolved-salt aqueous concentration and dissolved-CO2 aqueous
!     relative saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 YLSY(4),XLAY(4),XLSY(4),XMDGCX(LNGC),XMDLCX(LNGC),
     &       XLCYX(LNGC)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      REAL*8 XLCMX(LNGC),XLCY(LNGC,4),XLCZ(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XGCX(LNGC),XLCX(LNGC),DXLCX(LNGC)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_12'
!
!---  Pure-water density  ---
!
      ISRX = 1
      CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
!
!---  Guess dissolved salt aqueous mass fraction  ---
!
      YLSX = RHOLSX/(RHOLWX+RHOLSX)
      XLSX = 0.D+0
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &  WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
      CALL SOL_LS( TX,XLSMX )
      DYLSX = 1.D-4*XLSMX
      DXLAX = 0.D+0
      DO IGC = 1, NGC
        DXLCX(IGC) = 1.D-4*XLCMX(IGC)
        DXLAX = DXLAX + DXLCX(IGC)
      ENDDO
      XLSX = MIN( YLSX,XLSMX )
      CALL SP_B( TX,XLSX,PSBX )
      PVAX = MAX( PX-PSBX,0.D+0 )
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &  WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
      XLAX = 0.D+0
      DO IGC = 1, NGC
       XLAX = XLAX + XLCX(IGC)
      ENDDO
!
!---  Guess dissolved CO2 aqueous mass fraction  ---
!
      XLAX = PHILAX*XLAX
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCZ,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &  WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
      XLAZ = 0.D+0
      DO IGC = 1, NGC
        XLAZ = XLAZ + XLCZ(IGC)
      ENDDO
!      XLAX = 0.D+0
!      DO IGC = 1, NGC
!        XLAX = XLAX + XLCX(IGC)
!      ENDDO
      PVAX = MIN( XLAX/XLAZ,1.D+0 )*PGAX
      IF( PSBX+PVAX.GT.PGX ) THEN
        INDX = 9
        RLMSG = PSBX+PVAX
        CHMSG = CHMSGX(1)
        CALL WRMSGS( INDX )
      ENDIF
      NC = 0
!
!---  Two-variable Newton-Raphson loop (YLS, XLA)  ---
!
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = CHMSGX(2)
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 2,4
        YLSY(M) = YLSX
        DO IGC = 1, NGC
          XLCY(IGC,M) = XLCX(IGC)
        ENDDO
        IF( M.EQ.3 ) YLSY(M) = YLSX + DYLSX
        DO IGC = 1, NGC
          IF( M.EQ.4 ) XLCY(IGC,M) = XLCX(IGC) + DXLCX(IGC)
        ENDDO
        XLSY(M) = MIN( YLSY(M),XLSMX )
        CALL SP_B( TX,XLSY(M),PSBX )
        CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &    XGCX,XGWX,XLCX,XLSY(M),XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,
     &    ZGX,WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
        XLAX = 0.D+0
        DO IGC = 1, NGC
         XLAX = XLAX + XLCX(IGC)
        ENDDO
        PVAX = MIN( XLAY(M)/XLAX,1.D+0 )*PGAX
        IF( PSBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
          INDX = 9
          RLMSG = PSBX+PVAX
          CHMSG = CHMSGX(1)
          CALL WRMSGS( INDX )
        ENDIF
        CALL DENS_B( TX,PX,XLSY(M),RHOBX )
        DO IGC = 1, NGC
          XLCYX(IGC) = XLCY(IGC,M)
        ENDDO
        CALL DENS_L( TX,RHOBX,XLCYX,RHOLX )
        GX(M,1) = XLSY(M) - RHOLSX/RHOLX
        GX(M,2) = XLAY(M) - PHILAX*XLAX
   20 CONTINUE
      DXLAX = 0.D+0
      DO IGC = 1, NGC
        DXLAX = DXLAX + DXLCX(IGC)
      ENDDO
      RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
      RX(1,2) = (GX(4,1)-GX(2,1))/DXLAX
      RX(2,1) = (GX(3,2)-GX(2,2))/DYLSX
      RX(2,2) = (GX(4,2)-GX(2,2))/DXLAX
      RPX(1) = -GX(2,1)
      RPX(2) = -GX(2,2)
      CYLSX = (RPX(2)-RPX(1)*RX(2,2)/(RX(1,2)+SMALL))/
     &  (RX(2,1)-RX(1,1)*RX(2,2)/(RX(1,2)+SMALL))
      CXLAX = (RPX(2)-RPX(1)*RX(2,1)/(RX(1,1)+SMALL))/
     &  (RX(2,2)-RX(1,2)*RX(2,1)/(RX(1,1)+SMALL))
      YLSX = YLSX + CYLSX
      XLAX = XLAX + CXLAX
      IF( ABS(CYLSX).GT.(1.D-6*XLSMX) .OR.
     &  ABS(CXLAX).GT.(1.D-6*XLAMX) ) GOTO 10
      XLSX = MIN( YLSX,XLSMX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_12 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_13( TX,PX,PCX,PGX,RHOLSX,YLSX,XLSX,XLAX,CHMSGX,
     &          XMDGCX )
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
!     Compute dissolved-salt and -CO2 aqueous mass fractions from
!     dissolved-salt aqueous concentration and dissolved-CO2 aqueous
!     mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 YLSY(4),XLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      REAL*8 XLCZ(LNGC),XMDGCX(LNGC),XMDLCX(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XGCX(LNGC),XLCX(LNGC),DXLCX(LNGC)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_13'
!
!---  Pure-water density  ---
!
      ISRX = 1
      CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
!
!---  Guess dissolved-salt aqueous mass fractions  ---
!
      YLSX = RHOLSX/(RHOLWX+RHOLSX)
      CALL SOL_LS( TX,XLSMX )
      DYLSX = 1.D-4*XLSMX
      NC = 0
!
!---  Single-variable Newton-Raphson loop (YLS)  ---
!
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = CHMSGX(1)
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 2,3
        YLSY(M) = YLSX
        IF( M.EQ.3 ) YLSY(M) = YLSX + DYLSX
        XLSY(M) = MIN( YLSY(M),XLSMX )
        CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &    XGCX,XGWX,XLCZ,XLSY(M),XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX, 
     &    WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
        XLAZ = 0.D+0
        DO IGC = 1, NGC
          XLAZ = XLAZ + XLCZ(IGC)
        ENDDO
        PVAX = MIN( XLAX/XLAZ,1.D+0 )*PGAX
        CALL SP_B( TX,XLSY(M),PSBX )
        IF( PSBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
          INDX = 9
          RLMSG = PSBX+PVAX
          CHMSG = CHMSGX(2)
          CALL WRMSGS( INDX )
        ENDIF
        CALL DENS_B( TX,PX,XLSY(M),RHOBX )
        CALL DENS_L( TX,RHOBX,XLCX,RHOLX )
        GX(M,1) = XLSY(M) - RHOLSX/RHOLX
   20 CONTINUE
      RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
      RPX(1) = -GX(2,1)
      CYLSX = RPX(1)/RX(1,1)
      YLSX = YLSX + CYLSX
      IF( ABS(CYLSX).GT.(1.D-6*XLSMX) ) GOTO 10
      XLSX = MIN( YLSX,XLSMX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_13 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_21( TX,PX,PCX,PGX,PHILSX,RHOLAX,YLSX,XLSX,XLAX,
     &  CHMSGX,XMDGCX )
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
!     Compute dissolved-salt and -CO2 aqueous mass fractions from
!     dissolved-salt aqueous relative saturation and
!     dissolved-CO2 aqueous concentration
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XMDGCX(LNGC),XMDLCX(LNGC),XLCYX(LNGC)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      REAL*8 XLCY(LNGC,4),XLMCX(LNGC),XLCMX(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XGCX(LNGC),XLCX(LNGC),DXLCX(LNGC)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_21'
!
!---  Convert dissolved-salt relative saturation to mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      YLSX = PHILSX*XLSMX
      XLSX = MIN( YLSX,XLSMX )
      CALL DENS_B( TX,PX,XLSX,RHOBX )
      CALL SP_B( TX,XLSX,PSBX )
!
!---  Guess dissolved co2 aqueous mass fraction  ---
!
      XLAX = RHOLAX/(RHOBX+RHOLAX)
      XLSX = 0.D+0
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &  WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
      XLAMX = 0.D+0
      DO IGC = 1, NGC
       XLAMX = XLAMX + XLCMX (IGC)
      ENDDO
      DXLAX = 0.D+0
      DO IGC = 1, NGC
        DXLCX(IGC) = 1.D-4*XLCMX(IGC)
        DXLAX = DXLAX + DXLCX(IGC)
      ENDDO
      NC = 0
!
!---  Single-variable Newton-Raphson loop (XLA)  ---
!
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = CHMSGX(1)
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 2,3
        XLAY(M) = XLAX
        IF( M.EQ.3 ) XLAY(M) = XLAX + DXLAX
        CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &    XGCX,XGWX,XLCX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &    WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
        XLAX = 0.D+0
        DO IGC = 1, NGC
         XLAX = XLAX + XLCX(IGC)
        ENDDO
        PVAX = MIN( XLAY(M)/XLAX,1.D+0 )*PGAX
        IF( PSBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
          INDX = 9
          RLMSG = PSBX+PVAX
          CHMSG = CHMSGX(2)
          CALL WRMSGS( INDX )
        ENDIF
        DO IGC = 1, NGC
          XLCYX(IGC) = XLCY(IGC,M)
        ENDDO
        CALL DENS_L( TX,RHOBX,XLCYX,RHOLX )
        GX(M,1) = XLAY(M) - RHOLAX/RHOLX
   20 CONTINUE
      RX(1,1) = (GX(3,1)-GX(2,1))/DXLAX
      RPX(1) = -GX(2,1)
      CXLAX = RPX(1)/RX(1,1)
      XLAX = XLAX + CXLAX
      IF( ABS(CXLAX).GT.(1.D-6*XLAMX) ) GOTO 10
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_21 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_22( TX,PX,PCX,PGX,PHILSX,PHILAX,YLSX,XLSX,XLAX,
     &  CHMSGX,XMDGCX )
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
!     Compute dissolved-salt and -CO2 aqueous mass fractions from
!     dissolved-salt aqueous relative saturation and
!     dissolved-CO2 aqueous relative saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLCZ(LNGC),XMDGCX(LNGC),XMDLCX(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XGCX(LNGC),XLCX(LNGC),DXLCX(LNGC)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_22'
!
!---  Convert dissolved-salt relative saturation to mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      YLSX = PHILSX*XLSMX
      XLSX = MIN( YLSX,XLSMX )
      CALL SP_B( TX,XLSX,PSBX )
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCZ,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &  WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
      XLAZ = 0.D+0
      DO IGC = 1, NGC
       XLAZ = XLAZ + XLCZ(IGC)
      ENDDO
      XLAX = PHILAX*XLAZ
      PVAX = MIN( XLAX/XLAZ,1.D+0 )*PGAX
      IF( PSBX+PVAX.GT.PGX ) THEN
        INDX = 9
        RLMSG = PSBX+PVAX
        CHMSG = CHMSGX(1)
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_22 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_23( TX,PX,PCX,PGX,PHILSX,YLSX,XLSX,XLAX,CHMSGX,
     &            XMDGCX )
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
!     Compute dissolved-salt mass fraction from
!     dissolved-salt aqueous relative saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLCZ(LNGC),XMDGCX(LNGC),XMDLCX(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XGCX(LNGC)
      REAL*8 XLCX(LNGC),DXLCX(LNGC)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_23'
!
!---  Convert dissolved-salt relative saturation to mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      YLSX = PHILSX*XLSMX
      XLSX = MIN( YLSX,XLSMX )
      CALL SP_B( TX,XLSX,PSBX )
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &   XGCX,XGWX,XLCZ,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &   WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
      XLAZ = 0.D+0
      DO IGC = 1, NGC
       XLAZ = XLAZ + XLCZ(IGC)
      ENDDO
      PVAX = MIN( XLAX/XLAZ,1.D+0 )*PGAX
      IF( PSBX+PVAX.GT.PGX ) THEN
        INDX = 9
        RLMSG = PSBX+PVAX
        CHMSG = CHMSGX(1)
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_23 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_31( TX,PX,PCX,PGX,RHOLAX,YLSX,XLSX,XLAX,CHMSGX,
     &                  XMDGCX )
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
!     Compute dissolved-CO2 aqueous mass fractions from
!     dissolved-CO2 aqueous concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XMDGCX(LNGC),XMDLCX(LNGC)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      REAL*8 XLCMX(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XGCX(LNGC),XLCX(LNGC),DXLCX(LNGC)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_31'
!
!---  Limit dissolved-salt mass fraction to solubility limit  ---
!
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
      CALL DENS_B( TX,PX,XLSX,RHOBX )
      CALL SP_B( TX,XLSX,PSBX )
!
!---  Guess dissolved co2 aqueous mass fraction  ---
!
      XLAX = RHOLAX/(RHOBX+RHOLAX)
      XLSX = 0.D+0
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCMX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &  WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
      XLAMX = 0.D+0
      DO IGC = 1, NGC
        XLAMX = XLAMX + XLCMX(IGC)
      ENDDO
      DXLAX = 0.D+0
      DO IGC = 1, NGC
         DXLAX = DXLAX + 1.D-4*XLCMX(IGC)
      ENDDO
      NC = 0
!
!---  Single-variable Newton-Raphson loop (XLA)  ---
!
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = CHMSGX(1)
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 2,3
        XLAY(M) = XLAX
        IF( M.EQ.3 ) XLAY(M) = XLAX + DXLAX
        CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &    XGCX,XGWX,XLCX,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &    WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
        XLAX = 0.D+0
        DO IGC = 1, NGC
         XLAX = XLAX + XLCX(IGC)
        ENDDO
        PVAX = MIN( XLAY(M)/XLAX,1.D+0 )*PGAX
        IF( PSBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
          INDX = 9
          RLMSG = PSBX+PVAX
          CHMSG = CHMSGX(2)
          CALL WRMSGS( INDX )
        ENDIF
        CALL DENS_L( TX,RHOBX,XLCX,RHOLX )
        GX(M,1) = XLAY(M) - RHOLAX/RHOLX
   20 CONTINUE
      RX(1,1) = (GX(3,1)-GX(2,1))/DXLAX
      RPX(1) = -GX(2,1)
      CXLAX = RPX(1)/RX(1,1)
      XLAX = XLAX + CXLAX
      IF( ABS(CXLAX).GT.(1.D-6*XLAMX) ) GOTO 10
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_31 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_32( TX,PX,PCX,PGX,PHILAX,YLSX,XLSX,XLAX,CHMSGX,
     &                  XMDGCX )
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
!     Compute dissolved-CO2 aqueous mass fractions from
!     dissolved-CO2 aqueous relative saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLCZ(LNGC),XMDGCX(LNGC),XMDLCX(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XGCX(LNGC),XLCX(LNGC),DXLCX(LNGC)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_32'
!
!---  Limit dissolved-salt mass fraction to solubility limit  ---
!
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN(YLSX,XLSMX)
      CALL SP_B( TX,XLSX,PSBX )
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCZ,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &  WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
      DO IGC = 1, NGC
        XLCX(IGC) = PHILAX*XLCZ(IGC)
      ENDDO
      XLAZ = 0.D+0
      DO IGC = 1, NGC
       XLAZ = XLAZ + XLCZ(IGC)
      ENDDO
      PVAX = MIN( XLAZ/XLAZ,1.D+0 )*PGAX
      IF( PSBX+PVAX.GT.PGX ) THEN
        INDX = 9
        RLMSG = PSBX+PVAX
        CHMSG = CHMSGX(1)
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_32 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_33( TX,PX,PCX,PGX,YLSX,XLSX,XLAX,CHMSGX,XMDGCX )
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
!     Compute dissolved-salt mass fraction from
!     dissolved-salt aqueous relative saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLCZ(LNGC),XMDGCX(LNGC),XMDLCX(LNGC)
      REAL*8 XMGCX(LNGC),XMLCX(LNGC),XGCX(LNGC),XLCX(LNGC),DXLCX(LNGC)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_33'
!
!---  Limit dissolved-salt mass fraction to solubility limit  ---
!
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
      CALL SP_B( TX,XLSX,PSBX )
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGCX,XGWX,XLCZ,XLSX,XLWX,XMGCX,XMGWX,XMLCX,XMLSX,XLWX,ZGX,
     &  WTMGCX,XMDGCX,XMDGWX,XMDLCX,XMDLWX )
      XLAZ = 0.D+0
      DO IGC = 1, NGC
       XLAZ = XLAZ + XLCZ(IGC)
      ENDDO
      PVAX = MIN( XLAX/XLAZ,1.D+0 )*PGAX
      IF( PSBX+PVAX.GT.PGX ) THEN
        INDX = 9
        RLMSG = PSBX+PVAX
        CHMSG = CHMSGX(1)
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE NICKALLS( CAX,CBX,CCX,CDX,R1X,R2X,R3X )
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
!     This subroutine calculates the roots of a cubic equation using
!     the Nickalls formulation.
!
!     Nickalls, R.W.D. 1993. A new approach to solving the cubic:
!     Cardan's solution revealed.  Math. Gaz. 77:354-359.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
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
      SUB_LOG(ISUB_LOG) = '/NICKALLS'
      INCG = 1
!
!---  Nickalls cubic polynomial solver  ---
!
      XNX = -CBX/(3.D+0*CAX)
      YNX = CAX*(XNX**3) + CBX*(XNX**2) + CCX*XNX + CDX
      YN2X = YNX**2
      D2X = ((CBX**2)-(3.D+0*CAX*CCX))/((3.D+0*CAX)**2)
      IF( D2X.LE.0.D+0 ) THEN
        HX = 0.D+0
      ELSE
        DX = SQRT(D2X)
        HX = -2.D+0*(DX**3)
      ENDIF
      H2X = 4.D+0*(CAX**2)*(D2X**3)
!
!---  YN2X > H2X, 1 real root  ---
!
      IF( YN2X.GT.H2X ) THEN
        R1X = (5.D-1/CAX)*(-YNX+SQRT(YN2X-H2X))
        R2X = (5.D-1/CAX)*(-YNX-SQRT(YN2X-H2X))
        IF( R1X.LT.0.D+0 .AND. R2X.LT.0.D+0 ) THEN
          VX = XNX + (ABS(R1X)**(1.D+0/3.D+0))
     &      + (ABS(R2X)**(1.D+0/3.D+0))
          CHKX = CAX*(VX**3) + CBX*(VX**2) + CCX*VX + CDX
        ENDIF
        R3X = XNX + SIGN((ABS(R1X)**(1.D+0/3.D+0)),R1X)
     &    + SIGN((ABS(R2X)**(1.D+0/3.D+0)),R2X)
        R1X = R3X
        R2X = R3X
!
!---  YN2X < H2X, 3 distinct real roots  ---
!
      ELSEIF( YN2X.LT.H2X ) THEN
        THETAX = (ACOS(MIN(YNX/HX,1.D+0)))/3.D+0
        R1X = XNX + 2.D+0*DX*COS(THETAX)
        R2X = XNX + 2.D+0*DX*COS(((2.D+0*GPI)/3.D+0)+THETAX)
        R3X = XNX + 2.D+0*DX*COS(((4.D+0*GPI)/3.D+0)+THETAX)
!
!---  YN2X = H2X, 3 real roots (two or three equal roots)  ---
!
      ELSE
!
!---    HX /= 0 (two equal roots)  ---
!
        IF( ABS(HX)/EPSL.GT.EPSL ) THEN
          DX = YNX/(2.D+0*CAX)
          DX = SIGN((ABS(DX)**(1.D+0/3.D+0)),DX)
          R1X = XNX + DX
          R2X = XNX + DX
          R3X = XNX - 2.D+0*DX
!
!---    HX = 0 (three equal roots)  ---
!
        ELSE
          R1X = XNX
          R2X = XNX
          R3X = XNX
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of NICKALLS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PERM_R( SSX,PERMRFX,PORDX,IZN )
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
!     Calculation of permeability reduction factor.
!
!     Verma, A., and K. Pruess.  1988.  Thermohydrological Conditions
!     and Silica Redistribution Near High-Level Nuclear Wastes
!     Emplaced in Saturated Geological Formations.  Journal of
!     Geophysical Research, 93(B2):1159-1173.
!
!     Pruess, K., and J. Garcia.  2002.  Multiphase flow dynamics
!     during CO2 disposal into saline aquifers.  Environmental Geology
!     http://link.springer.de/link/service/journals/00254/contents
!     /01/00498/paper/s00254-001-0498-3ch110.html
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PERM_R'
!
!---  Reduced porosity with salt precipitation  ---
!
      PORDX = MAX( PORDX*(1.D+0-SSX),PORDX*PERM(5,IZN),1.D-12 )
!
!---  Normalized porosity  ---
!
      PORD_NX = MAX( (1.D+0-SSX-PERM(5,IZN))/(1.D+0-PERM(5,IZN)),0.D+0 )
!
!---  Tube area ratio  ---
!
      OMEGAX = 1.D+0 + (1.D+0/PERM(4,IZN))/((1.D+0/PERM(5,IZN))-1.D+0)
!
!---  Permeability reduction factor  ---
!
      PERMRFX = (PORD_NX**2)*(1.D+0-PERM(4,IZN)+PERM(4,IZN)/(OMEGAX**2))
     &  /(1.D+0-PERM(4,IZN)+PERM(4,IZN)*
     &  ((PORD_NX/(PORD_NX+OMEGAX-1.D+0))**2))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PERM_R group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PERM_V( SSX,PERMRFX,PORDX,PORM1,PORM0,IZN )
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
!     Calculation of permeability reduction factor.
!
!     Verma, A., and K. Pruess.  1988.  Thermohydrological Conditions
!     and Silica Redistribution Near High-Level Nuclear Wastes
!     Emplaced in Saturated Geological Formations.  Journal of
!     Geophysical Research, 93(B2):1159-1173.
!
!     Pruess, K., and J. Garcia.  2002.  Multiphase flow dynamics
!     during CO2 disposal into saline aquifers.  Environmental Geology
!     http://link.springer.de/link/service/journals/00254/contents
!     /01/00498/paper/s00254-001-0498-3ch110.html
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PERM_V'
!
!---  Reduced porosity with salt precipitation  ---
!
      PORDX = MAX( PORDX*(1.D+0-SSX),PORDX*PERM(5,IZN),1.D-12 )
!
!--- Precipitated mineral saturation
!
      IF( ISLC(43).EQ.1 ) THEN
        SPM = 1.D+0 - PORM1/PORM0
      ENDIF
!
!---  Normalized porosity  ---
!
      PORD_NX = MAX( (1.D+0-SSX-SPM-PERM(5,IZN))/(1.D+0-PERM(5,IZN)),
     &  0.D+0 )
!
!---  Tube area ratio  ---
!
      OMEGAX = 1.D+0 + (1.D+0/PERM(4,IZN))/((1.D+0/PERM(5,IZN))-1.D+0)
!
!---  Permeability reduction factor  ---
!
      PERMRFX = (PORD_NX**2)*(1.D+0-PERM(4,IZN)+PERM(4,IZN)/(OMEGAX**2))
     &  /(1.D+0-PERM(4,IZN)+PERM(4,IZN)*
     &  ((PORD_NX/(PORD_NX+OMEGAX-1.D+0))**2))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PERM_V group  ---
!
      RETURN
      END

      SUBROUTINE RDPF_GC
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
!     Water-NComponent-NaCl-Energy Mode (STOMP-WNSE)
!     Read gas component properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FILES
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,FDUM,RDUM,UNTS
      CHARACTER*512 CHDUM
      CHARACTER*64 GCNMX,GC1NMX,GC2NMX
      LOGICAL FCHK
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDPF_GC'
!
!---  Check for property file  ---
!
      INQUIRE( FILE='comp_prop.dat', FORM=FDUM, EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        INDX = 4
        CHMSG = 'Missing Component Property File: comp_prop.dat'
        CALL WRMSGS( INDX )
      ENDIF
      OPEN(UNIT=26,FILE='comp_prop.dat',STATUS='OLD',FORM='FORMATTED')
      WRITE (IWR,'(/,A)') '~ Component Properties'
!
!---  Loop over number of gas components  ---
!
      DO 500 IGC = 1,NGC
        REWIND(26)
        WRITE (IWR,'(/,2A)') 'Component Name: ',GCNM(IGC)
!
!---    Critical Properties  ---
!
        VARB = 'Critical Properties'
        WRITE (IWR,'(2X,4A,1PE11.4)') TRIM(VARB)
!
!---    Read new input line ---
!---    Skip input lines that begin with '#' or '!'  ---
!
   10   READ(26,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' .OR. CHDUM(1:1).EQ.'!' ) GOTO 10
        CALL LCASE( CHDUM )
        IF( CHDUM(1:9).NE.'~critical' ) GOTO 10
!
!---    Read new input line ---
!
   15   READ(26,'(A)') CHDUM
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---    Read species name
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GCNMX)
        CALL LCASE( GCNMX )
        IF( GCNMX /= GCNM(IGC) ) GOTO 15
!
!---    Read molecular weight, kg/kmol ---
!
        VARB = 'Molecular Weight'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(1,IGC))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE (IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',GCPP(1,IGC)
        INDX = 0
        IUNKG = 1
        IUNMOL = -1
        CALL RDUNIT(UNTS,GCPP(1,IGC),INDX)
!
!---    Read normal boiling point temperature, C ---
!
        VARB = 'Normal Boiling Point, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(2,IGC))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE (IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',GCPP(2,IGC)
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,GCPP(2,IGC),INDX)
        GCPP(2,IGC) = GCPP(2,IGC) + TABS
!
!---    Read critical temperature, C ---
!
        VARB = 'Critical Temperature, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(3,IGC))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE (IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',GCPP(3,IGC)
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,GCPP(3,IGC),INDX)
        GCPP(3,IGC) = GCPP(3,IGC) + TABS
!
!---    Read critical pressure, Pa ---
!
        VARB = 'Critical Pressure, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(4,IGC))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE (IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',GCPP(4,IGC)
        INDX = 0
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        CALL RDUNIT(UNTS,GCPP(4,IGC),INDX)
!
!---    Read critical molar volume ---
!
        VARB = 'Critical Molar Volume, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(5,IGC))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE (IWR,'(4X,3A,1PE11.4)') VARB(1:IVR),
     &    UNTS(1:NCH),': ',GCPP(5,IGC)
        INDX = 0
        IUNM = 3
        IUNMOL = -1
        CALL RDUNIT(UNTS,GCPP(5,IGC),INDX)
!
!---    Convert to units of cm^3/mole  ---
!
        GCPP(5,IGC) = 1.D+3*GCPP(5,IGC)
!
!---    Read critical compressiblity  ---
!
        VARB = 'Critical Compressibility: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(6,IGC))
        WRITE (IWR,'(4X,A,1PE11.4)') VARB(1:IVR),GCPP(6,IGC)
!
!---    Read Pitzer acentric factor ---
!
        VARB = 'Pitzer acentric factor: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(7,IGC))
        WRITE (IWR,'(4X,A,1PE11.4)') VARB(1:IVR),GCPP(7,IGC)
!
!---    Read dipole moment, debyes ---
!
        VARB = 'Dipole moment, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(8,IGC))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE (IWR,'(4X,3A,1PE11.4)') VARB(1:IVR),
     &    UNTS(1:NCH),': ',GCPP(8,IGC)
        INDX = 0
        CALL RDUNIT(UNTS,GCPP(8,IGC),INDX)
!
!---    Isobaric molar specific heat  ---
!
        VARB = 'Isobaric Molar Specific Heat Coefficients: '
!
!---    Read new input line ---
!---    Skip input lines that begin with '#' or '!'  ---
!
   20   READ(26,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' .OR. CHDUM(1:1).EQ.'!' ) GOTO 20
        CALL LCASE( CHDUM )
        IF( CHDUM(1:5).NE.'~heat' ) GOTO 20
!
!---    Read new input line ---
!
   25   READ(26,'(A)') CHDUM
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---    Read species name
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GCNMX)
        CALL LCASE( GCNMX )
        IF( GCNMX /= GCNM(IGC) ) GOTO 25
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(9,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(10,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(11,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(12,IGC))
        WRITE (IWR,'(2X,A,$)') VARB(1:IVR)
        WRITE (IWR,'(A)') ' (Cp = a + bT + cT^2 + dT^3)'
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant a: ',GCPP(9,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant b: ',GCPP(10,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant c: ',GCPP(11,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant d: ',GCPP(12,IGC)
!
!---    Peng-Robinson pure component parameters (a, b, fw)  ---
!
        GCPP(13,IGC) = 0.45724D+0*(RCU**2)*(GCPP(3,IGC)**2)/
     &    GCPP(4,IGC)
        GCPP(14,IGC) = 0.07780D+0*RCU*GCPP(3,IGC)/GCPP(4,IGC)
        GCPP(15,IGC) = 0.37464D+0 + 1.54226D+0*GCPP(7,IGC) -
     &    0.26992D+0*(GCPP(7,IGC)**2)
!
!---    Partial molar volume ---
!
        VARB = 'Partial Molar Volume Coefficients: '
!
!---    Read new input line ---
!---    Skip input lines that begin with '#' or '!'  ---
!
   30   READ(26,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' .OR. CHDUM(1:1).EQ.'!' ) GOTO 30
        CALL LCASE( CHDUM )
        IF( CHDUM(1:8).NE.'~partial' ) GOTO 30
!
!---    Read new input line ---
!
   35   READ(26,'(A)') CHDUM
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---    Read species name
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GCNMX)
        CALL LCASE( GCNMX )
        IF( GCNMX /= GCNM(IGC) ) GOTO 35
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(23,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(24,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(25,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(26,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(27,IGC))
        WRITE (IWR,'(2X,A,$)') VARB(1:IVR)
        WRITE (IWR,'(A)') ' (PMV = a + bT + cT^2 + dT^3 + eT^4)'
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant a: ',GCPP(23,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant b: ',GCPP(24,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant c: ',GCPP(25,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant d: ',GCPP(26,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant e: ',GCPP(27,IGC)

   40   READ(26,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' .OR. CHDUM(1:1).EQ.'!' ) GOTO 40
        CALL LCASE( CHDUM )
        IF( CHDUM(1:18).NE.'~nonaqueous binary' ) GOTO 40
!
!---    Non-aqueous Binary interaction coefficients ---
!
        VARB = 'Non-aqueous Binary Interaction Coefficients: '
        WRITE(IWR,'(2X,A)') VARB

        BIJNA(IGC,IGC) = 0.D+0 
        WRITE (IWR,'(4X,A,"-",A,":",ES10.3)') 
     &    TRIM(GCNM(IGC)),TRIM(GCNM(IGC)),BIJNA(IGC,IGC)
!
!---    Read new input line ---
!
   45   READ(26,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'~' ) GOTO 51
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---    Read species names
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GC1NMX)
        CALL LCASE( GC1NMX )
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GC2NMX)
        CALL LCASE( GC2NMX )
        IF( GC1NMX /= GCNM(IGC) .AND. GC2NMX /= GCNM(IGC) ) GOTO 45

        IF( TRIM(GC1NMX) == 'h2o' .OR. TRIM(GC2NMX) == 'h2o' ) THEN
          CALL RDDPR(ISTART,ICOMMA,CHDUM,BIJNA(IGC,NGC+1))
          WRITE (IWR,'(4X,A,"-",A,":",ES10.3)') 
     &      TRIM(GC1NMX),TRIM(GC2NMX),BIJNA(IGC,NGC+1)
          GOTO 45
        ELSE
          DO JGC = 1,NGC
            IF( (GC1NMX == GCNM(IGC) .AND. GC2NMX == GCNM(JGC)) 
     &        .OR. (GC1NMX == GCNM(JGC) .AND. GC2NMX == GCNM(IGC))
     &          ) THEN
              CALL RDDPR(ISTART,ICOMMA,CHDUM,BIJNA(IGC,JGC))
              WRITE (IWR,'(4X,A,"-",A,":",ES10.3)') 
     &          TRIM(GC1NMX),TRIM(GC2NMX),BIJNA(IGC,JGC)
            ENDIF
          END DO
          GOTO 45
        ENDIF
!
!---    Read new input line ---
!---    Skip input lines that begin with '#' or '!'  ---
!
   50   READ(26,'(A)') CHDUM
   51   IF( CHDUM(1:1).EQ.'#' .OR. CHDUM(1:1).EQ.'!' ) GOTO 50
        CALL LCASE( CHDUM )
        IF( CHDUM(1:16).NE.'~nonaqueous temp' ) GOTO 50
!
!---    Non-aqueous Temperature coefficients ---
!
        VARB = 'Non-aqueous Temperature Dependence: '
        WRITE(IWR,'(2X,A)') VARB
        BIJTD(IGC,IGC) = 0.D+0
        WRITE (IWR,'(4X,A,"-",A,":",ES10.3)')  
     &    TRIM(GCNM(IGC)),TRIM(GCNM(IGC)),BIJTD(IGC,IGC)
!
!---    Read new input line ---
!
   55   READ(26,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'~' ) GOTO 61
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---    Read species names
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GC1NMX)
        CALL LCASE( GC1NMX )
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GC2NMX)
        CALL LCASE( GC2NMX )
        IF( GC1NMX /= GCNM(IGC) .AND. GC2NMX /= GCNM(IGC) ) GOTO 55

        IF( TRIM(GC1NMX) == 'h2o' .OR. TRIM(GC2NMX) == 'h2o' ) THEN
          CALL RDDPR(ISTART,ICOMMA,CHDUM,BIJTD(IGC,NGC+1))
          WRITE (IWR,'(4X,A,"-",A,":",ES10.3)')  
     &      TRIM(GC1NMX),TRIM(GC2NMX),BIJTD(IGC,NGC+1)
          GOTO 55
        ELSE
          DO JGC = 1,NGC
            IF( (GC1NMX == GCNM(IGC) .AND. GC2NMX == GCNM(JGC)) 
     &        .OR. (GC1NMX == GCNM(JGC) .AND. GC2NMX == GCNM(IGC))
     &          ) THEN
              CALL RDDPR(ISTART,ICOMMA,CHDUM,BIJTD(IGC,JGC))
            WRITE (IWR,'(4X,A,"-",A,":",ES10.3)') 
     &        TRIM(GC1NMX),TRIM(GC2NMX),BIJTD(IGC,JGC)
            ENDIF
          END DO
          GOTO 55
        ENDIF

   59   CONTINUE 
!
!---    Read new input line ---
!---    Skip input lines that begin with '#' or '!'  ---
!
   60   READ(26,'(A)') CHDUM
   61   IF( CHDUM(1:1).EQ.'#' .OR. CHDUM(1:1).EQ.'!' ) GOTO 60
        CALL LCASE( CHDUM )
        IF( CHDUM(1:8).NE.'~aqueous' ) GOTO 60
!
!---    Aqueous Binary interaction coefficients ---
!
        VARB = 'Aqueous Binary Interaction Coefficients: '
        WRITE(IWR,'(2X,A)') VARB
!
!---    Read new input line ---
!
   65   READ(26,'(A)') CHDUM
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---    Read species name 
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GCNMX)
        CALL LCASE( GCNMX )
        IF( GCNMX /= GCNM(IGC) ) GOTO 65

        DO INOC = 37,48
          CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(INOC,IGC))
        END DO
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant a1: ',GCPP(37,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant a2: ',GCPP(38,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant a3: ',GCPP(39,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant a4: ',GCPP(40,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant b1: ',GCPP(41,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant b2: ',GCPP(42,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant b3: ',GCPP(43,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant b4: ',GCPP(44,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant c1: ',GCPP(45,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant c2: ',GCPP(46,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant c3: ',GCPP(47,IGC)
        WRITE (IWR,'(4X,A,1PE11.4)') 'Constant c4: ',GCPP(48,IGC)
!
!---    Viscosity coefficients ---
!
        VARB = 'Viscosity Coefficients: '
        WRITE(IWR,'(2X,A)') VARB
!
!---    Read new input line ---
!---    Skip input lines that begin with '#' or '!'  ---
!
   70   READ(26,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' .OR. CHDUM(1:1).EQ.'!' ) GOTO 70
        CALL LCASE( CHDUM )
        IF( CHDUM(1:10).NE.'~viscosity' ) GOTO 70
!
!---    Read new input line ---
!
   75   READ(26,'(A)') CHDUM
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---    Read species name 
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GCNMX)
        CALL LCASE( GCNMX )
        IF( GCNMX /= GCNM(IGC) ) GOTO 75

        VARB = 'Friction Theory Coefficients: '
        WRITE(IWR,'(4X,A)') VARB
        DO INOC = 49,64
          CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(INOC,IGC))
          WRITE(IWR,'(4X,ES15.6)')GCPP(INOC,IGC)
        END DO
!
!---    Density coefficients ---
!
        VARB = 'Density Coefficients: '
        WRITE(IWR,'(2X,A)') VARB
!
!---    Read new input line ---
!---    Skip input lines that begin with '#' or '!'  ---
!
   80   READ(26,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' .OR. CHDUM(1:1).EQ.'!' ) GOTO 80
        CALL LCASE( CHDUM )
        IF( CHDUM(1:8).NE.'~density' ) GOTO 80
!
!---    Read new input line ---
!
   85   READ(26,'(A)') CHDUM
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---    Read species name 
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GCNMX)
        CALL LCASE( GCNMX )
        IF( GCNMX /= GCNM(IGC) ) GOTO 85

        VARB = 'Peneloux Correction Factor Coefficients: '
        GCPP(65,IGC)=0.40768
        GCPP(66,IGC)=0.294410

        VARB = 'Rackett Compressibility Factor: '
        WRITE(IWR,'(4X,A,$)') TRIM(VARB)
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(16,IGC))
        WRITE(IWR,'(ES15.6)') GCPP(16,IGC)

        VARB = 'Volume Shift Coefficients'
        WRITE(IWR,'(4X,A)') TRIM(VARB)
        DO INOC = 67,70
          CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(INOC,IGC))
          WRITE(IWR,'(4X,ES15.6)')GCPP(INOC,IGC)
        END DO
!
!---    Read new input line ---
!---    Skip input lines that begin with '#' or '!'  ---
!
   90   READ(26,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' .OR. CHDUM(1:1).EQ.'!' ) GOTO 90
        CALL LCASE( CHDUM )
!
!---    Thermal conductivity coefficients ---
!
        IF( CHDUM(1:8).NE.'~thermal' ) GOTO 90
        VARB = 'Thermal conductivity Coefficients: '
        WRITE(IWR,'(2X,A)') VARB
!
!---    Read new input line ---
!
   95   READ(26,'(A)') CHDUM
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---    Read species name 
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,GCNMX)
        CALL LCASE( GCNMX )
        IF( GCNMX /= GCNM(IGC) ) GOTO 95

        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(19,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(20,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(21,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(22,IGC))
        WRITE(IWR,*)GCPP(19,IGC),GCPP(20,IGC),GCPP(21,IGC),GCPP(22,IGC)

        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(28,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(29,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(30,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(31,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(32,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(33,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(34,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(35,IGC))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GCPP(36,IGC))
        WRITE(IWR,*) GCPP(28,IGC),GCPP(29,IGC),GCPP(30,IGC),
     &    GCPP(31,IGC),GCPP(32,IGC),GCPP(33,IGC),GCPP(34,IGC),
     &    GCPP(35,IGC),GCPP(36,IGC)
!
!---  Go to next gas component  ---
!
  500 CONTINUE
!
!---  End of RDPF_GC group ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RHOG_I( TX,PX,XGWX,XMGCX,RHOGX )
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
!     H2O-NaCl-CO2 Equation of state.  This subroutine finds the gas
!     density at fixed gas water/CO2 mass fraction, temperature,
!     and pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
                REAL*8 XMGCX(LNGC)        
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RHOG_I'
      ISRX = 2
!
!---  Estimate partial pressures  ---
!
      XGAX = MAX( 1.D+0-XGWX,0.D+0 )
      XMGWX = (XGWX/WTMW)/(XGWX/WTMW + XGAX/WTMA)
      PVWX = XMGWX*PX
!
!---  Newton scheme  ---
!
  100 CONTINUE
      PVAX = MAX( PX-PVWX,0.D+0 )
      CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
      CALL DENS_GC( TX,PVAX,XMGCX,RHOGCX,ISRX,ZGX )
      RHOGX = RHOGWX+RHOGCX
      DPVWY = SIGN( 1.D-1,(5.D-1*PX-PVWX) )
      PVWY = PVWX + DPVWY
      PVAY = MAX( PX-PVWY,0.D+0 )
      CALL DENS_W( TX,PVWY,RHOLWY,RHOGWY,ISRX )
      CALL DENS_GC( TX,PVAY,XMGCX,RHOGCY,ISRX,ZGX )
      RHOGY = RHOGWY+RHOGCY
      FX = XGWX - RHOGWX/RHOGX
      DFX = (RHOGWX/RHOGX - RHOGWY/RHOGY)/DPVWY
      DPVWX = -FX/DFX
      PVWX = PVWX + DPVWX
      IF( ABS(DPVWX).GT.1.D-4 ) GOTO 100
      PVAX = MAX( PX-PVWX,0.D+0 )
      CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
      CALL DENS_GC( TX,PVAX,XMGCX,RHOGCX,ISRX,ZGX )
      RHOGX = RHOGWX+RHOGCX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RHOG_I group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFT_L( TX,XLSX,SFTLX )
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
!     This subroutine calculates the surface tension of aqueous
!     solutions of sodium chloride as a function of temperature
!     and salt concentration.
!
!     Abramzon, A.A., and R.D. Gaukhberg.  1993.  Surface tension of
!     salt solutions.  Russian Journal of Applied Chemistry,
!     66(6):1139-1146.
!
!     Lide, D.R., and H.V. Kehiaian.  1994.  CRC Handbook of
!     Thermophysical and Thermochemical Data.  CRC Press, Boca Raton,
!     Florida.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/SFT_L'
!
!---  Convert mass fraction to molality  ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
!
!---  Pure water vapor surface tension as a function of temperature
!     by D. Lide and H. Kehiaian
!
      SFTWX = 1.D-3*(75.6592D+0 - 1.40959D-1*TX - 2.66317D-4*(TX**2))
!
!---  Function by A. Abramzon and R. Gaukhberg  ---
!
      SFTLX = SFTWX + 1.57D-3*GLSX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFT_L group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SOL_LS( TX,XLSMX )
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
!     This subroutine calculates the mass fraction of NaCl salt in
!     saturated aqueous solutions.
!
!     McKibbin, R., and A. McNabb.  1993.  Modeling the phase
!     boundaries and fluid properties of the system H2O-NaCl at high
!     temperatures and pressures.  Proceedings 15th NZ Geothermal
!     Workshop, University of Auckland, New Zealand.
!
!     Temperature range 0-800 C.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 SFX(3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE SFX
      DATA SFX / 2.6218D-1, 7.2D-5, 1.06D-6 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SOL_LS'
!
!---  Maximum NaCl saturation  ---
!
      XLSMX = SFX(1) + SFX(2)*TX + SFX(3)*(TX**2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SOL_LS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_B( TX,XLSX,PSBX )
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
!     This subroutine calculates the saturation pressure of NaCl brine
!     as a function of temperature and salt concentration.
!
!     tx - temperature, C
!     xlsx - mass fraction of NaCl salt
!     psbx - saturation pressure of brine, Pa
!
!     Haas, J.L., Jr.  1976.  Physical Properties of the Coexisting
!     Phases and Thermochemical Properties of the H2O Component in
!     Boiling NaCl Solutions, U.S. Geological Survey Bulletin, 1421-A,
!     United States Government Printing Office, Washington.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 SAX(3),SBX(5)
!
!----------------------Data Statements---------------------------------!
!
      SAVE SAX,SBX
      DATA SAX / 5.93582D-6, -5.19386D-5, 1.23156D-5 /
      DATA SBX / 1.15420D-6, 1.41254D-7, -1.92476D-8, -1.70717D-9,
     &  1.05390D-10 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SP_B'
!
!---  Convert temperature to Kelvin and mass fraction to molality  ---
!
      TKX = TX + TABS
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
!
!---  Concentration dependent coefficients  ---
!
      AX = 1.D+0
      DO 10 I = 1,3
        AX = AX + SAX(I)*(GLSX**I)
   10 CONTINUE
      BX = 0.D+0
      DO 20 I = 1,5
        BX = BX + SBX(I)*(GLSX**I)
   20 CONTINUE
!
!---  Temperature dependent coefficient  ---
!
      CX = 1.D+0/(AX + BX*TKX)
!
!---  Equivalent pure water temperature  ---
!
      TWX = EXP( CX*LOG(TKX) ) - TABS
      CALL SP_W( TWX,PSBX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_W( TX,PSWX )
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
!     Saturation pressure (Pa) of pure water as a function of
!     temperature.
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 K(9),PSWX,THETAX,THETAY,TX
      INTEGER I
!
!----------------------Data Statements---------------------------------!
!
      SAVE K
      DATA K / -7.691234564D+0, -2.608023696D+1, -1.681706546D+2,
     &  6.423285504D+1, -1.189646225D+2, 4.167117320D+0,
     &  2.097506760D+1, 1.D+9, 6.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SP_W'
!
!---  The K-function (saturation line)  ---
!
      THETAX = (TX+TABS)/TCRW
      THETAY = (1.D+0-THETAX)
      PSWX = 0.D+0
      DO 10 I = 1,5
        PSWX = PSWX + K(I)*(THETAY**I)
   10 CONTINUE
      PSWX = PSWX/((1.D+0 + K(6)*THETAY + K(7)*(THETAY**2))*THETAX)
      PSWX = PSWX - THETAY/(K(8)*(THETAY**2)+K(9))
      PSWX = EXP(PSWX)*PCRW
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SR_W( TX,PX,ISRX )
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
!     Formulation subregion as a function of temperature and pressure.
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 LFCX(3),PSWX,PSRBX,THETAX,TOLX
      INTEGER ISRX
!
!----------------------Data Statements---------------------------------!
!
      SAVE LFCX,TOLX
      DATA LFCX / 1.574373327D+1, -3.417061978D+1, 1.931380707D+1 /
      DATA TOLX / 1.D-2 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SR_W'
!
!---  Subcritical-critical temperature  ---
!
      TKX = TX+TABS
      IF( TKX.LE.TCRW ) THEN
!
!---    The K-function (saturation line)  ---
!
        CALL SP_W( TX,PSWX )
!
!---    Subregions 5/6 and 1/4 boundaries  ---
!
        IF( TX.LE.350.D+0 ) THEN
          IF( (PX-PSWX).GE.TOLX ) THEN
            ISRX = 1
          ELSEIF( (PX-PSWX).LE.-TOLX ) THEN
            ISRX = 2
          ELSE
            ISRX = 6
          ENDIF
        ELSE
          IF( (PX-PSWX).GE.TOLX ) THEN
            ISRX = 4
          ELSEIF( (PX-PSWX).LE.-TOLX ) THEN
            ISRX = 2
          ELSE
            ISRX = 5
          ENDIF
        ENDIF
!
!---  Supercritical temperature  ---
!
      ELSE
!
!---    The L-function (subregions 2/3 boundary)  ---
!
        THETAX = TKX/TCRW
        PSRBX = PCRW*(LFCX(1) + LFCX(2)*THETAX + LFCX(3)*(THETAX**2))
        ISRX = 2
        IF( PX.GT.PSRBX ) ISRX = 3
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SR_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_B( TX,XLSX,TKLWX,TKBX )
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
!     This subroutine calculates the thermal conductivity of pure
!     water as a function of temperature.
!
!     tx - temperature, C
!     xlsx - mass fraction of aqueous NaCl
!     tkbx - thermal conductivity of NaCl brine, W/m K
!
!     Ozbek, H. and S.L. Phillips.  1980.  Thermal conductivity of
!     aqueous sodium chloride solutions from 20 to 330 C.
!     J. Chem. Engr. Data, 25:263-267.
!
!     Temperature Range:  20 - 330 C
!     NaCl Concentration Range:  0 - saturation
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 SCX(6)
!
!----------------------Data Statements---------------------------------!
!
      SAVE SCX
      DATA SCX / 2.3434D-3, -7.924D-6, 3.924D-8, 1.06D-5, -2.D-8,
     &  1.2D-10 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_B'
!
!---  Yusufova correlation  ---
!
      YLSX = 1.D+2*XLSX
      TKBX = 1.D+0 - (SCX(1) + SCX(2)*TX + SCX(3)*(TX**2))*YLSX +
     &  (SCX(4) + SCX(5)*TX + SCX(6)*(TX**2))*(YLSX**2)
      TKBX = TKBX*TKLWX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_G( TX,TKAX,TKWX,XMAX,XMWX,TKGX )
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
!     Calculate the gas thermal conductivity using the modification of
!     Mason and Saxena to the Waasiljewa equation.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York. pp: 530-531.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 XMAX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_G'

      CHIWX = 0.D+0 
      CHIAX = 0.D+0 
      TKGX =  0.D+0 
!
!---  Reduced, inverse thermal conductivity  ---
!
      TKX = TX+TABS
      TKGX = XMWX*TKWX + SUM(XMAX)*TKAX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_G group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_GC( PGX,RHOMGX,THKGX,TX,XMGCX,ZX )
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
!     Low-pressure thermal conductivity of pure gas components
!     computed from a polynomial function (pp. 516)
!     Pressure correction from the excess thermal conductivity
!     correlation (pp. 520) following Stiel and Thodos.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 THKM_T(LNGC),GAMMAX(LNGC),THKGCX(LNGC),XMGCX(LNGC),MMX
      REAL*8 TRX(LNGC) 
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_GC'

      TKX = TX+273.15D+0
      RCUBX = RCU*1.D-2
      PBX = PGX*1.D-5
      THKM_T(1:NGC)=0.D+0
!
!--- Mixture molecular weight---
!
      DO IGC = 1,NGC
        TRX(IGC) = TKX/GCPP(3,IGC)
      ENDDO
      DO IGC = 1, NGC
        IF(XMGCX(IGC)<=1.D-8) XMGCX(IGC)=0.D+0
        IF( SUM(XMGCX).NE.0.D+0 ) THEN
        XMGCX(IGC) = XMGCX(IGC)/SUM(XMGCX)
        ELSE
        XMGCX(IGC) = 0.D+0
        ENDIF
      ENDDO
      WTMX = 0.D+0
      DO IGC = 1, NGC
        WTMX = WTMX + XMGCX(IGC)*GCPP(1,IGC)
      ENDDO
      IF(NGC==1)WTMX = GCPP(1, NGC)
      RMM = 0.D+0
      DO IGC = 1, NGC
        RMM = RMM + XMGCX(IGC)/(GCPP(1,IGC)**0.30D+0)
      ENDDO
!
!---  Pure component thermal conductivities  ---
!
      DO 100 IGC = 1,NGC
!
!---    Low-pressure polynomial correlation  ---
!
        THKGCX(IGC) = GCPP(19,IGC) + GCPP(20,IGC)*TKX
     &    + GCPP(21,IGC)*(TKX**2) + GCPP(22,IGC)*(TKX**3)

 100  CONTINUE
!
!---   Excess thermal conductivity correlations  ---
!
       THKM=0.D+0
       TCMX = 0.D+0
       VCMX = 0.D+0
       WMX = 0.D+0
       MMX = 0.D+0
       DO IGC = 1, NGC
       IF(XMGCX(IGC)>0.D+0) THEN
        IF(NGC==1) THEN
         ZS = 1.D+0
        ELSE
         ZS = XMGCX(IGC)/((GCPP(1,IGC)**0.3D+0)*RMM)
        ENDIF
        GCPP(14,IGC) = (0.0778D+0*RCUBX*GCPP(3,IGC)/
     &                      (GCPP(4,IGC)*1.D-5))
        GCPP(15,IGC) = 0.37464D+0 + 1.54226D+0*GCPP(7,IGC) -
     &          0.26992D+0*(GCPP(7,IGC)**2.D+0)
        GCPP(17,IGC) = (1.D+0 + GCPP(15,IGC)
     &     *(1.D+0-SQRT(TRX(IGC))))**2.D+0
        GCPP(13,IGC) = (0.45724D+0*GCPP(17,IGC)*(RCUBX**2.D+0)*
     &           (GCPP(3,IGC)**2.D+0/(GCPP(4,IGC)*1.D-5)))
        SAMX = GCPP(13,IGC)*PBX/((RCUBX*TKX)**2.D+0)
        SBMX = GCPP(14,IGC)*PBX/(RCUBX*TKX)
! 
!---    Polynomial coefficients for Peng-Robinson equation of state  ---
!
        CAX = 1.D+0
        CBX = -(1.0D+0 - SBMX)
        CCX = SAMX - 2.D+0*SBMX - 3.D+0*(SBMX**2.D+0)
        CDX = -(SAMX*SBMX - (SBMX**2.D+0) - (SBMX**3.D+0))
        CALL NICKALLS( CAX,CBX,CCX,CDX,Z1X,Z2X,Z3X )
!
!---    Compressibility factor Z---
!
        ZGX = MAX( Z1X,Z2X,Z3X )
        ZNX = MIN( Z1X,Z2X,Z3X )
!
!---    Molar volume, cm^3/mol  ---
!
        VGX = RCUBX*TKX*ZGX/PBX
        VNX = RCUBX*TKX*ZNX/PBX
        DELTA1 = 1.D+0+SQRT(2.D+0)
        DELTA2 = 1.D+0-SQRT(2.D+0)
        GSTAR = (ZNX-ZGX)+LOG(ABS(((ZGX-SBMX)/(ZNX-SBMX))))+
     +   (1.D+0/(DELTA2-DELTA1))*(SAMX/SBMX)
     +   *LOG(ABS(((ZGX+DELTA2*SBMX)*(ZNX+DELTA1*SBMX)/
     +    ((ZNX+DELTA2*SBMX)*(ZGX+DELTA1*SBMX)))))
!
!---    Gas conditions  ---
!
        SFX = 2.D+0/3.D+0
        TCIJ = GCPP(3,IGC)
        IF( TKX.GE.TCIJ ) THEN
          VX = VGX
          ZX = ZGX
        ELSE
         IF(GSTAR.GE.0.D+0) THEN
           VX = VGX
           ZX = ZGX
           GOTO 651
         ELSE
          VX = VNX
          ZX = ZNX
         GOTO 651
         ENDIF
        ENDIF
  651  CONTINUE
       IF(XMGCX(IGC)==1.D+0) ZS = 1.D+0
       TCIJX = GCPP(3,IGC)
       VCIJX = GCPP(5,IGC)
       TCMX = TCIJX
       VCMX = VCIJX
       WMX = GCPP(7,IGC)
       MMX = GCPP(1,IGC)
       ZCMX = GCPP(6,IGC)
       PCMX = GCPP(4,IGC)*1.D-5
       PCMX = PCMX*1.D+0
       TOUX = 210.D+0*(TCMX*(MMX**3)/(PCMX**4))**(1.D+0/6.D+0)
 499   CONTINUE
       VX=(ZX*RCUBX*TKX/PBX)
       RHORX = VCMX/VX
       IF( RHORX<0.5) THEN
         THKM = GCPP(28,IGC)*(exp(GCPP(29,IGC)*RHORX)+
     &      GCPP(30,IGC))/(TOUX*(ZCMX**5))
       ELSEIF(RHORX>=0.5 .and. RHORX<2.0) THEN
         THKM = GCPP(31,IGC)*(exp(GCPP(32,IGC)*RHORX)+
     &           GCPP(33,IGC))/(TOUX*(ZCMX**5))
       ELSEIF(RHORX>=2.0) THEN
         THKM = GCPP(34,IGC)*(exp(GCPP(35,IGC)*RHORX)+
     &           GCPP(36,IGC))/(TOUX*(ZCMX**5))
       endif
        THKM_T(IGC) = ZS/(THKM+THKGCX(IGC))
       ENDIF
       enddo
      XLPTHKS = 0.D+0
      THKM_SUM = 0.D+0
      DO IGC = 1, NGC
       IF(XMGCX(IGC)==0.D+0) THKGCX(IGC)=0.D+0
       XLPTHKS = XLPTHKS + XMGCX(IGC)*THKGCX(IGC)
       THKM_SUM = THKM_SUM+THKM_T(IGC)
      ENDDO
      THKGX =  1.D+0/THKM_SUM
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_GC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_W( TX,PX,RHOX,TKX )
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
!     Thermal conductivity (W/m K) of pure water as a function of
!     temperature and density.
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 SAX(4),SBX(3),CBX(2),SDX(4),CCX(6)
!
!----------------------Data Statements---------------------------------!
!
      SAVE SAX,SBX,CBX,SDX,CCX
      SAVE TREF,RHOREF,PREF,TKREF
      DATA TREF / 6.4727D+2 /
      DATA RHOREF / 3.17763D+2 /
      DATA PREF / 2.2115D+7 /
      DATA TKREF / 1.D+0 /
      DATA SAX / 0.0102811D+0, 0.0299621D+0, 0.0156146D+0,
     &  -0.00422464 /
      DATA SBX / -0.397070D+0, 0.400302D+0, 1.060000D+0 /
      DATA CBX / -0.171587D+0, 2.392190D+0 /
      DATA SDX / 0.0701309D+0, 0.0118520D+0, 0.00169937D+0, -1.0200D+0 /
      DATA CCX / 0.642857D+0, -4.11717D+0, -6.17937D+0, 0.00308976D+0,
     &  0.0822994D+0, 10.0932D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_W'
!
!---  Reduced temperature, density, and pressure  ---
!
      THETAX = (TX+TABS)/TREF
      RHOBX = RHOX/RHOREF
      BETAX = PX/PREF
!
!---  Zero term  ---
!
      TKX = 0.D+0
      ISAX = 0
      DO 10 I = 0,3
        ISAX = ISAX + 1
        TKX = TKX + SAX(ISAX)*(THETAX**I)
   10 CONTINUE
      TKX = SQRT(THETAX)*TKX
!
!---  First term  ---
!
      TKX = TKX + SBX(1) + SBX(2)*RHOBX +
     &  SBX(3)*EXP(CBX(1)*((RHOBX+CBX(2))**2))
!
!---  Second term  ---
!
      DTHETAX = ABS(THETAX-1.D+0) + CCX(4)
      CQX = 2.D+0 + CCX(5)/(DTHETAX**(3.D+0/5.D+0))
      IF( DTHETAX.GE.1.D+0 ) THEN
        CSX = 1.D+0/DTHETAX
      ELSE
        CSX = CCX(6)/(DTHETAX**(3.D+0/5.D+0))
      ENDIF
      TKX = TKX + ((SDX(1)/(THETAX**10))+SDX(2))*
     &  (RHOBX**(9.D+0/5.D+0))*
     &  EXP(CCX(1)*(1.D+0-(RHOBX**(14.D+0/5.D+0))))
      TKX = TKX + SDX(3)*CSX*(RHOBX**CQX)*
     &  EXP((CQX/(1.D+0+CQX))*(1.D+0-(RHOBX**(1.D+0+CQX))))
      VARX = CCX(2)*(THETAX**(3.D+0/2.D+0)) + CCX(3)/(RHOBX**5)
      VARX = MAX( VARX,-50.D+0 )
      TKX = TKX + SDX(4)*EXP(VARX)
!
!---  Dimensionalize thermal conductivity  ---
!
      TKX = TKX*TKREF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THKS_LG( THKSX,IZN )
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
!     Calculation of average intrinsic permeability around well,
!     using logarithmic averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/THKS_LG'
!
!---  Average intrinsic permeability  ---
!
      THKSX = 0.D+0
      NTHKSX = 0
      IF( (THKS(1,IZN)/EPSL).GT.EPSL ) THEN
        THKSX = THKSX + LOG(THKS(1,IZN))
        NTHKSX = NTHKSX+1
      ENDIF
      IF( (THKS(2,IZN)/EPSL).GT.EPSL ) THEN
        THKSX = THKSX + LOG(THKS(2,IZN))
        NTHKSX = NTHKSX+1
      ENDIF
      IF( (THKS(3,IZN)/EPSL).GT.EPSL ) THEN
        THKSX = THKSX + LOG(THKS(3,IZN))
        NTHKSX = NTHKSX+1
      ENDIF
      IF( NTHKSX.GT.0 ) THEN
        REALX = REAL(NTHKSX)
        THKSX = EXP( THKSX/REALX )
      ELSE
        THKSX = 0.D+0
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THKS_LG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_GC( TX,PX,ZGX,XMAX,RHOGX,VISCGX,VISGCX )
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
!     Calculation of CO2 viscosity using the formulation of
!     Fenghour et al. in the temperature range 200K < T < 1500K
!     and densities up to 1400 kg/m^3.
!
!     Fenghour, A., W. A. Wakeham, V. Vesovic.  1998.  The viscosity
!     of carbon dioxide.  J. Phys. Chem. Ref. Data, 27(1):31-41.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 SAX(5),SBX(5),BIJ(3,3),XMAX(LNGC),VISCGX(LNGC)
      REAL*8 RNC(LNGC),RKR(LNGC),RKA(LNGC),RKRR(LNGC) 
      REAL*8 TRX(LNGC)
!
!----------------------Data Statements---------------------------------!
!
      SAVE SAX,SBX,ESPX
      DATA SAX / 0.235156D+0, -0.491266D+0,
     &  5.211155D-2, 5.347906D-2, -1.537102D-2 /
      DATA SBX / 0.4071119D-2, 0.7198037D-4,
     &  0.2411697D-16, 0.2971072D-22, -0.1627888D-22 /
      DATA ESPX / 251.196D+0 /
!     
      SAVE RMM,RN0MX,RKRMX,RKAMX,RKRRMX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_GC'
!
!---  Zero-density-limit component  ---
!
      TKX = TX+TABS
!
!---  Critical viscosity of CO2 in cm^3/mol 
!
      VCRAX = 93.9D0
      DO IGC = 1, NGC
       IF (XMAX(IGC).EQ.0) XMAX(IGC) = SMALL
       XMAX(IGC)=XMAX(IGC)/SUM(XMAX)
      ENDDO
!
!---Initialize variables
!
       RMM = 0.D+0
       RN0MX = 0.D+0
       RKRMX = 0.D+0
       RKAMX = 0.D+0
       RKRRMX = 0.D+0
       VISGCX = 0.D+0
       VISCGX(1:NGC) = 0.D+0
!
!--- Critical and other properties
!
      PBX = PX*1.D-5
      RCUBX = RCU*1.D-2
      OMEGAXW = 0.344D+0
!
!--- Convert critical pressure to bars
!
      DO IGC = 1, NGC
        GCPP(4,IGC) = GCPP(4,IGC)*1.D-5
      ENDDO
       BIJ(1:NGC,1:NGC)=0.D+0
!      ENDIF
      DO IGC = 1,NGC
        TRX(IGC) = TKX/GCPP(3,IGC)
      ENDDO
       DO IGC = 1, NGC
        DO JGC = 1, NGC
          BIJ(IGC,JGC) = BIJNA(IGC,JGC)+BIJTD(IGC,JGC)*TRX(IGC)
        ENDDO
       ENDDO
!
!--- Check for mole fraction ---
!
       DO IGC = 1, NGC
        IF(XMAX(IGC).LE.1.0D-7)XMAX(IGC)=0.D+0
       ENDDO

!
!--- Convert critical pressure to bars
!
      DO IGC = 1, NGC
            GCPP(14,IGC) = (0.0778D+0*RCUBX*GCPP(3,IGC)/
     &                      GCPP(4,IGC))
            GCPP(15,IGC) = 0.37464D+0 + 1.54226D+0*GCPP(7,IGC) -
     &          0.26992D+0*(GCPP(7,IGC)**2.D+0)
         GCPP(17,IGC) = (1.D+0 + GCPP(15,IGC)
     &     *(1.D+0-SQRT(TRX(IGC))))**2.D+0
         GCPP(13,IGC) = (0.45724D+0*GCPP(17,IGC)*(RCUBX**2.D+0)*
     &           (GCPP(3,IGC)**2.D+0/GCPP(4,IGC)))
      ENDDO 
!
!--- Peng-Robinson concentration dependent parameters
!
          SAMX = 0.D+0
          SBMX = 0.D+0
          DO IGC = 1, NGC
           DO JGC = 1, NGC
            SAMX = SAMX + (XMAX(IGC)*XMAX(JGC)*SQRT(GCPP(13,IGC)*
     &          GCPP(13,JGC))*(1.D+0 - BIJ(IGC,JGC)))
           ENDDO
          ENDDO
          DO IGC = 1, NGC
           SBMX = SBMX + (XMAX(IGC)*GCPP(14,IGC))
          ENDDO
          AX = SAMX
          BX = SBMX

          SAMX = SAMX*PBX/((RCUBX*TKX)**2.D+0)
          SBMX = SBMX*PBX/(RCUBX*TKX)
! 
!---    Polynomial coefficients for Peng-Robinson equation of state  ---
!
         CAX = 1.D+0
         CBX = -(1.0D+0 - SBMX)
         CCX = SAMX - 2.D+0*SBMX - 3.D+0*(SBMX**2.D+0)
         CDX = -(SAMX*SBMX - (SBMX**2.D+0) - (SBMX**3.D+0))

         CALL NICKALLS( CAX,CBX,CCX,CDX,Z1X,Z2X,Z3X )
!
!---    Compressibility factor Z---
!
        ZGX = MAX( Z1X,Z2X,Z3X )
        ZNX = MIN( Z1X,Z2X,Z3X )
!
!---    Molar volume, cm^3/mol  ---
!
        VGX = RCUBX*TKX*ZGX/PBX
        VNX = RCUBX*TKX*ZNX/PBX
        DELTA1 = 1.D+0+SQRT(2.D+0)
        DELTA2 = 1.D+0-SQRT(2.D+0)
        GSTAR = (ZNX-ZGX)+LOG(ABS(((ZGX-SBMX)/(ZNX-SBMX))))+
     +   (1.D+0/(DELTA2-DELTA1))*(SAMX/SBMX)
     +   *LOG(ABS(((ZGX+DELTA2*SBMX)*(ZNX+DELTA1*SBMX)/
     +    ((ZNX+DELTA2*SBMX)*(ZGX+DELTA1*SBMX)))))
!
!---    Gas conditions  ---
!
        SFX = 2.D+0/3.D+0
        SUMXMAX = 0.D+0
        DO IGC = 1,NGC
          SUMXMAX = SUMXMAX + XMAX(IGC)
        ENDDO
        IF( NGC.EQ.1 .OR. SUMXMAX.LT.EPSL ) THEN
          TCIJ = GCPP(3,NGC)
        ELSE
        TCIJ = 0.D+0
        CR1MX = 0.D+0
        DO IGC = 1, NGC
         CR1MX = CR1MX + XMAX(IGC)*GCPP(5,IGC)
         TCIJ = TCIJ + XMAX(IGC)*GCPP(5,IGC)*GCPP(3,IGC)
        ENDDO
        TCIJ = TCIJ/CR1MX
        ENDIF
        IF( TKX.GE.TCIJ ) THEN
          VX = VGX
          ZX = ZGX
        ELSE
         IF(GSTAR.GE.0.D+0) THEN
           VX = VGX
           ZX = ZGX
           GOTO 651
         ELSE
          VX = VNX
          ZX = ZNX
         GOTO 651
         ENDIF

          W1X = PBX*(VGX-VNX)
          W2X = RCUBX*TKX*LOG((VGX-SBMX)/(VNX-SBMX)) +
     &      (SAMX/(SQRT(TKX)*SBMX))*LOG(((VGX+SBMX)*VNX)/(VNX+SBMX)*VGX)
!
!---      Gas phase condition  ---
!
          IF( (W2X-W1X)/EPSL.LE.EPSL ) THEN
            VX = VGX
            ZX = ZGX
!
!---      Liquid phase condition  ---
!
          ELSEIF( (W1X-W2X)/EPSL.LT.EPSL ) THEN
            VX = VNX
            ZX = ZNX
!
!---      Gas-liquid phase condition, use gas molar volume  ---
!
          ELSE
            VX = VGX
            ZX = ZGX
          ENDIF
        ENDIF
  651  CONTINUE
         IF(NGC/=1) THEN
          VGX = RCUBX*TKX*ZGX/PBX
          RMM = 0.D+0
          DO IGC = 1, NGC
           RMM = RMM + XMAX(IGC)/(GCPP(1,IGC)**0.3D+0)
          ENDDO
          RPRMX = RCUBX*TKX/(VGX - BX)
          RPAMX = -AX/((VGX**2.D+0)+2.D+0*BX*VGX-(BX**2.D+0))
         ENDIF
      DO IGC = 1, NGC
         IF( NGC.EQ.1 ) THEN
          VGX = RCUBX*TKX*ZGX/PBX
         ENDIF
          RPR = RCUBX*TKX/(VGX - GCPP(14,IGC))
          RPA = -GCPP(13,IGC)/(VGX**2.D+0 + 
     &      2.D+0*GCPP(14,IGC)*VGX - GCPP(14,IGC)**2.D+0)
         IF( NGC.EQ.1 ) THEN
          RPRMX = RPR
          RPAMX = RPA
         ENDIF
         TRRX = 1.D+0/TRX(IGC) - 1.D+0
         SIX = RCUBX*GCPP(3,IGC)/GCPP(4,IGC)
!
!    Friction coefficients
!
         XKA00=GCPP(49,IGC)
         XKA10=GCPP(50,IGC)
         XKA11=GCPP(51,IGC)
         XKA20=GCPP(52,IGC)
         XKA21=GCPP(53,IGC)
         XKA22=GCPP(54,IGC)
         XKR00=GCPP(55,IGC)
         XKR10=GCPP(56,IGC)
         XKR11=GCPP(57,IGC)
         XKR20=GCPP(58,IGC)
         XKR21=GCPP(59,IGC)
         XKR22=GCPP(60,IGC)
         XKRR21=GCPP(61,IGC)
         XKAC=GCPP(62,IGC)
         XKRC=GCPP(63,IGC)
         XKRRC=GCPP(64,IGC)

      DELTA_KA = XKA00*TRRX + (XKA10 +
     &  (XKA11*SIX))*(EXP(TRRX)-1.D+0) +
     &  (XKA20 + XKA21*SIX +XKA22*(SIX**2.D+0))*
     &  (EXP(2.D+0*TRRX) - 1.D+0)
      RKA(IGC) = XKAC + DELTA_KA
      RNFA = RKA(IGC)*RPA/GCPP(4,IGC)
      DELTA_KR = XKR00*TRRX + (XKR10 + (XKR11*SIX))*
     &  (EXP(TRRX) - 1.D+0) +
     &  (XKR20 + XKR21*SIX + XKR22*(SIX**2.D+0))*
     &  (EXP(2.D+0*TRRX) - 1.D+0)
      DELTA_KRR = XKRR21*SIX*(EXP(2.D+0*(TRRX + 1.D+0)) 
     &  - 1.0D+0)*(TRRX**2.D+0)
      RKR(IGC) = XKRC + DELTA_KR
      RKRR(IGC) = XKRRC + DELTA_KRR

      RNFR = RKR(IGC)*(RPR/GCPP(4,IGC)) 
     &  + RKRR(IGC)*(RPR/GCPP(4,IGC))**2.D+0
      RNFC = RNFR + RNFA
      RNC(IGC) = 7.9483D+0*(SQRT(GCPP(1,IGC))
     &  *(GCPP(4,IGC)**0.66667D+0)/(GCPP(3,IGC)**0.166667D+0))
      RNF = RNFC*RNC(IGC)

      RMULT = 18.0323D+0*(1.2593D+0*TKX/GCPP(3,IGC))**(-0.76830D+0)
     &  - 7.27371D+0
      OMGAX = (1.16145D+0/(1.2593D+0*TKX/GCPP(3,IGC))**0.14874D+0) +
     &  (0.52487D+0/EXP(0.77320D+0*(1.2593D+0*TKX/GCPP(3,IGC)))) +
     &  (2.16178D+0/EXP(2.43787D+0*1.2593D+0*TKX/GCPP(3,IGC))) -
     &  6.435D-4*((1.2593D+0*TKX/GCPP(3,IGC))**0.14874D+0)*
     &  SIN(3.14D+0*RMULT/180.D+0)
      FC = 1.0D+0-0.2756D+0*GCPP(7,IGC)
      RN0 = 40.785D+0*SQRT(GCPP(1,IGC)*TKX)*FC/(OMGAX*
     &           (GCPP(5,IGC)**0.666667D+0))
      IF(XMAX(IGC).GT.0.0D+0) THEN
        VISCGX(IGC) = (RN0 + RNF)*1D-7
        VISGCX = VISCGX(IGC)
      ELSE
        VISCGX(IGC) = 0.0D+0
        VISGCX = VISCGX(IGC)
      ENDIF
      ENDDO
!
!--- Mixture viscosity ---
!
      IF(NGC>1) THEN
       DO IGC = 1, NGC
        IF(XMAX(IGC).LE.0.D+0) CYCLE
        ZS = XMAX(IGC)/((GCPP(1,IGC)**0.3D+0)*RMM)
        RN0MX = RN0MX + (XMAX(IGC)*LOG(RN0))
        RKRMX = RKRMX + (ZS*RNC(IGC)*RKR(IGC)/GCPP(4,IGC))
        RKAMX = RKAMX + (ZS*RNC(IGC)*RKA(IGC)/GCPP(4,IGC))
        RKRRMX = RKRRMX + (ZS*RNC(IGC)*RKRR(IGC)/(GCPP(4,IGC)**2.D+0))
       ENDDO
        RN0MX=EXP(RN0MX)
        RNFMX = RKRMX*RPRMX + RKAMX*RPAMX + RKRRMX*(RPRMX**2.D+0)
        VISGCX = (RN0MX + RNFMX)*1D-7
        IF(XMAX(NGC).LE.0.0D+0) VISGCX=VISCGX(1) 
       ENDIF
!
!--- Reset critical pressure to pascals
!
      Do IGC = 1, NGC
       GCPP(4,IGC) = GCPP(4,IGC)/1.D-5
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_GC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_B( TX,XLSX,VISWX,VISBX )
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
!     Viscosity (Pa s) of brine as a function of temperature (C),
!     NaCl brine mass fraction and pure water viscosity (Pa s)
!
!     Phillips, S.L., A. Igbene, J.A. Fair, H. Ozbek, and M. Tavana.
!     1981.  A Technical Databook for Geothermal Energy Utilization
!     LBL-12810, UC-66a, Lawrence Berkeley Laboratory, University of
!     California, Berkeley, California.
!
!     Temperature Range 10-350 C
!     Pressure Range 0.1-50 MPa
!     NaCl Concentration Range 0-5 Molal (mol NaCl/kg H2O)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 SAX(5)
!
!----------------------Data Statements---------------------------------!
!
      SAVE SAX
      DATA SAX / 0.0816D+0, 0.0122D+0, 0.000128D+0, 0.000629D+0,
     &  -0.7D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_B'
!
!---  Convert mass fraction to molality  ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
!
!---  Formulation of Phillips et al.  ---
!
      VISBX = VISWX*(1.D+0 + SAX(1)*GLSX + SAX(2)*(GLSX**2) +
     &  SAX(3)*(GLSX**3) + SAX(4)*TX*(1.D+0-EXP(SAX(5)*GLSX)))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_G( VISGCX,VISWX,XMGCX,XMGWX,VISGX )
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
!     Calculate the gas viscosity using the method of Wilke, p.407.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York. pp: 332-337.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VISGCX(LNGC),XMGCX(LNGC),CHIGC(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_G'
!
!---  Water contribution  ---
!
      CHIW = 0.D+0
      DO 10 IGC = 1,NGC
        PHIX = ((1.D+0+SQRT(VISWX/(VISGCX(IGC)+SMALL))*
     &    SQRT(SQRT(GCPP(1,IGC)/WTMW)))**2)/
     &    SQRT(8.D+0*(1.D+0 + WTMW/GCPP(1,IGC)))
        CHIW = CHIW + XMGCX(IGC)*PHIX
   10 CONTINUE
!
!---  Component contributions  ---
!
      DO 30 JGC = 1,NGC
        CHIGC(JGC) = 0.D+0
        PHIX = ((1.D+0+SQRT(VISGCX(JGC)/(VISWX+SMALL))*
     &    SQRT(SQRT(WTMW/GCPP(1,JGC))))**2)/
     &    SQRT(8.D+0*(1.D+0 + GCPP(1,JGC)/WTMW))
        CHIGC(JGC) = CHIGC(JGC) + XMGWX*PHIX
        DO 20 IGC = 1,NGC
          PHIX = ((1.D+0+SQRT(VISGCX(JGC)/(VISGCX(IGC)+SMALL))*
     &      SQRT(SQRT(GCPP(1,IGC)/GCPP(1,JGC))))**2)/
     &      SQRT(8.D+0*(1.D+0 + GCPP(1,JGC)/GCPP(1,IGC)))
          CHIGC(JGC) = CHIGC(JGC) + XMGCX(IGC)*PHIX
   20   CONTINUE
   30 CONTINUE
      VISGX = XMGWX*VISWX/CHIW
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_G group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_L( XMLCX,VISBX,VISCX,VISLX )
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
!     Viscosity (Pa s) of aqueous NaCl solutions with dissolved CO2,
!     following the method of Grunberg and Nissan, calibrated against
!     the data of Kumagai and Yokoyama.
!
!     Kumagai, A., and C. Yokoyama.  1999.  Viscosities of aqueous
!     NaCl solutions containing CO2 at high pressures.  J. Chem. Eng.
!     Data, 44:227-229.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids. pp. 474-475.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
            REAL*8 XMLCX(LNGC), VISCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_L'
!
!---  Formulation of Grunberg and Nissan  ---
!
      XMLAX = 0.D+0
      DO IGC = 1, NGC
        XMLAX = XMLAX + XMLCX(IGC)
      ENDDO
      VISLX = (1.D+0-XMLAX)*LOG(VISBX) 
      DO IGC = 1, NGC
        IF( VISCX(IGC).LE.0.D+0 ) CYCLE
        VISLX = VISLX + XMLCX(IGC)*LOG(VISCX(IGC))
      ENDDO
      VISLX = EXP(VISLX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_L group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_W( TX,PX,RHOWX,VISWX )
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
!     Viscosity (Pa s) of pure water as a function of temperature and
!     density.
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 CHX(46)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CHX,TREF,RHOREF,PREF,VISREF
      DATA TREF / 6.4727D+2 /
      DATA RHOREF / 3.17763D+2 /
      DATA PREF / 2.2115D+7 /
      DATA VISREF / 5.5071D+1 /
      DATA CHX / 1.D+0, 9.78197D-1, 5.79829D-1, -2.02354D-1,
     &  5.132047D-1, 3.205656D-1, 0.D+0, 0.D+0, -7.782567D-1,
     &  1.885447D-1, 2.151778D-1, 7.317883D-1, 1.241044D+0,
     &  1.476783D+0, 0.D+0, 0.D+0, -2.818107D-1, -1.070786D+0,
     & -1.263184D+0, 0.D+0, 0.D+0, 0.D+0, 1.778064D-1,
     &  4.605040D-1, 2.340379D-1, -4.924179D-1, 0.D+0, 0.D+0,
     &  -4.176610D-2, 0.D+0, 0.D+0, 1.600435D-1, 0.D+0, 0.D+0,
     &  0.D+0, -1.578386D-2, 0.D+0, 0.D+0, 0.D+0, 0.D+0, 0.D+0,
     &  0.D+0, 0.D+0, -3.629481D-3,  0.D+0, 0.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_W'
!
!---  Reduced temperature, density, and pressure  ---
!
      THETAX = (TX+TABS)/TREF
      RHOBX = RHOWX/RHOREF
      BETAX = PX/PREF
!
!---  Zero term  ---
!
      VISWX = 0.D+0
      ICH = 0
      DO 10 I = 0,3
        ICH = ICH + 1
        VISWX = VISWX + CHX(ICH)/(THETAX**I)
   10 CONTINUE
      VISWX = SQRT(THETAX)/VISWX
!
!---  First term  ---
!
      VISAX = 0.D+0
      DO 30 I = 0,5
        DO 20 J = 0,6
          ICH = (J*6) + I + 5
          VISAX = VISAX + CHX(ICH)*(((1.D+0/THETAX)-1.D+0)**I)*
     &      ((RHOBX-1.D+0)**J)
   20   CONTINUE
   30 CONTINUE
      VISWX = VISWX*EXP(RHOBX*VISAX)
!
!---  Second term  ---
!
      IF( THETAX.GE.0.9970 .AND. THETAX.LE.1.0082 .AND.
     &  RHOBX.GE.0.755 .AND. RHOBX.LE.1.290 ) THEN
        DPX = 1.D-1
        DBETAX = DPX/PREF
        PIX = PX+DPX
        CALL SR_W( TX,PIX,ISRX )
        CALL DENS_W( TX,PIX,RHOLX,RHOVX,ISRX )
        IF( (1.D+0-ABS(RHOLX/RHOWX)).LT.(1.D+0-ABS(RHOVX/RHOWX)) ) THEN
          RHOBIX = RHOLX/RHOREF
        ELSE
          RHOBIX = RHOVX/RHOREF
        ENDIF
        CHIX = RHOBX*(RHOBIX-RHOBX)/DBETAX
        IF( CHIX.GE.21.93 ) VISWX = VISWX*0.922D+0*(CHIX**0.0263D+0)
      ENDIF
!
!---  Dimensionalize viscosity  ---
!
      VISWX = 1.D-6*VISWX*VISREF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VPL_B( TX,PSBX,PCX,RHOBX,PVBX,XLSX )
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
!     Vapor pressure lowering of brine.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas.  Geothermics, 26(4): 437-464.
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      SUB_LOG(ISUB_LOG) = '/VPL_B'
!
!---  Kelvin's equation  ---
!
      TKX = TX + TABS
      WTMX = XLSX*WTMS + (1.D+0-XLSX)*WTMW
      PVBX = PSBX*EXP( -WTMW*(PCX**1.25D+0)/(RHOBX*RCU*TKX) )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VPL_B group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WEBB_BC( IZN )
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
!     Webb saturation and capillary pressure matching points for
!     the Brooks-Corey capillary pressure-saturation
!     function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/WEBB_BC'
!
!---  Find the saturation matching point  ---
!
      SRX = SCHR(4,IZN)
      PSIX = SCHR(1,IZN)
      CLX = MAX( SCHR(3,IZN),SMALL )
!
!---  Use the matrix saturation at 0.4 LOG10(HDOD) as
!     the initial guess  ---
!
      HDX = 1.D+1**(4.D-1*LOG10(HDOD))
      SMPX = (PSIX/HDX)**CLX
      SMPX = SMPX*(1.D+0-SRX) + SRX
!
!---  Newton-Raphson iteration for the matrix saturation
!     matching point  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      SEMPX = (SMPX-SRX)/(1.D+0-SRX)
      FX = LOG10(HDOD) - LOG10((PSIX/(SEMPX**(1.D+0/CLX))))
     &  - 1.D+0/(LOG(1.D+1)*CLX*(SMPX-SRX))
      DFX = 1.D+0/(LOG(1.D+1)*CLX*(SMPX-SRX))
     &  + 1.D+0/(LOG(1.D+1)*CLX*((SMPX-SRX)**2))
      DSMPX = -FX/DFX
      SMPX = MAX( SMPX+DSMPX,SRX+1.D-12 )
!
!---  No convergence on saturation matching point  ---
!
      IF( NC.GT.32 ) THEN
        INDX = 7
        IMSG = N
        CHMSG = 'No Convergence on Saturation '
     &    // 'Matching Point @ Node: '
        CALL WRMSGS( INDX )
      ENDIF
      IF( ABS(DSMPX).GT.1.D-9 ) GOTO 100
      SCHR(8,IZN) = SMPX
!
!---  Find the capillary head matching point  ---
!
      SEMPX = (SMPX-SRX)/(1.D+0-SRX)
      SCHR(9,IZN) = PSIX/(SEMPX**(1.D+0/CLX))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WEBB_BC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WEBB_VG( IZN )
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
!     Webb saturation and capillary pressure matching points for
!     the van Genuchten capillary pressure-saturation
!     function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/WEBB_VG'
!
!---  Find the matrix saturation matching point  ---
!
      SRX = SCHR(4,IZN)
      ALPHAX = SCHR(1,IZN)
      CNX = SCHR(3,IZN)
      CMX = SCHR(14,IZN)
!
!---  Use the matrix saturation at 0.4 LOG10(HDOD) as
!     the initial guess  ---
!
      HDX = 1.D+1**(4.D-1*LOG10(HDOD))
      SMPX = (1.D+0/(1.D+0 + (ALPHAX*HDX)**CNX))**CMX
      SMPX = SMPX*(1.D+0-SRX) + SRX
!
!---  Newton-Raphson iteration for the matrix saturation
!     matching point  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      SEMPX = (SMPX-SRX)/(1.D+0-SRX)
      ESEMPX = (1.D+0/SEMPX)**(1.D+0/CMX)
      FX1 = LOG10(HDOD)
      FX1 = FX1 - LOG10(((ESEMPX-1.D+0)**(1.D+0/CNX))/ALPHAX)
      FX1 = FX1 - (SMPX/(SMPX-SRX))/(LOG(1.D+1)*CNX*CMX*
     &  (1.D+0-(SEMPX**(1.D+0/CMX))))
      SMPY = SMPX + 1.D-8
      SEMPX = (SMPY-SRX)/(1.D+0-SRX)
      ESEMPX = (1.D+0/SEMPX)**(1.D+0/CMX)
      FX2 = LOG10(HDOD)
      FX2 = FX2 - LOG10(((ESEMPX-1.D+0)**(1.D+0/CNX))/ALPHAX)
      FX2 = FX2 - (SMPY/(SMPY-SRX))/(LOG(1.D+1)*CNX*CMX*
     &  (1.D+0-(SEMPX**(1.D+0/CMX))))
      DFX = (FX2-FX1)/1.D-8
      DSMPX = -FX1/DFX
      SMPX = MAX( SMPX+DSMPX,SRX+1.D-12 )
!
!---  No convergence on matrix saturation matching point  ---
!
      IF( NC.GT.32 ) THEN
        INDX = 7
        IMSG = N
        CHMSG = 'No Convergence on Saturation '
     &    // 'Matching Point @ Node: '
        CALL WRMSGS( INDX )
      ENDIF
      IF( ABS(DSMPX).GT.1.D-9 ) GOTO 100
      SCHR(8,IZN) = SMPX
!
!---  Find the matrix capillary head matching point  ---
!
      SEMPX = (SMPX-SRX)/(1.D+0-SRX)
      SCHR(9,IZN) = (1.D+0/ALPHAX)*
     &  (((1.D+0/(SEMPX**(1.D+0/CMX)))-1.D+0)**(1.D+0/CNX))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WEBB_VG group  ---
!
      RETURN
      END


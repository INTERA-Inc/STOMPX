!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNCN( TX,SLX,PORDX,SX,RHOLX,XLSX,INDX )
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
!     Compute salt (sodium chloride) concentration
!     and aqueous-brine solution density.
!
!     INDX = 0: salt volumetric concentration loaded in SX
!     INDX = 1: salt aqueous concentration loaded in SX
!     INDX = 2: relative salt aqueous concentration loaded in SX
!     INDX = 3: salt mass fraction loaded in SX
!     INDX = 4: salt mass fraction loaded in XLSX
!
!     This subroutine allows for salt concentrations and salt
!     aqueous mass fractions greater than
!     the solubility limit.  The returned brine density, however,
!     will not exceed the density at maximum solubility.  Returned
!     values are volumetric salt concentration in SX, aqueous salt
!     mass fraction in XLSX, and brine density in RHOLX.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July 1995.
!     Last Modified by MD White, Battelle, January 26, 2000.
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
      SUB_LOG(ISUB_LOG) = '/BRNCN'
!
!---  Compute maximum dissolved salt mass fraction ---
!
      CALL BRNS( TX,XLSM )
      CALL BRND( TX,XLSM,RHOLX,RHLSM,DRHLSM )
!
!---  Function for sodium chloride Elder's problem  ---
!
      IF( ISLC(22).EQ.4 ) RHOLX = 1.D+3
!
!---  Compute aqueous density and salt volumetric concentration from
!     the salt aqueous mass fraction  ---
!
      IF( INDX.EQ.4 ) THEN
        CALL BRND( TX,XLSX,RHOLX,RHLSX,DRHLSX )
        RHLSX = MIN( RHLSX,RHLSM )
        SX = XLSX*RHLSX*PORDX*SLX
      ELSEIF( INDX.EQ.3 ) THEN
        XLSX = SX
        CALL BRND( TX,XLSX,RHOLX,RHLSX,DRHLSX )
        RHLSX = MIN( RHLSX,RHLSM )
        IF( XLSX.GE.XLSM ) THEN
          SX = XLSX*RHLSX*PORDX*SLX
        ELSE
          SX = RHOLX*XLSX*PORDX*SLX/(1.D+0-XLSX)
        ENDIF
      ELSEIF( INDX.EQ.2 ) THEN
        XLSX = SX*XLSM
        CALL BRND( TX,XLSX,RHOLX,RHLSX,DRHLSX )
        RHLSX = MIN( RHLSX,RHLSM )
        IF( XLSX.GE.XLSM ) THEN
          SX = XLSX*RHLSX*PORDX*SLX
        ELSE
          SX = RHOLX*XLSX*PORDX*SLX/(1.D+0-XLSX)
        ENDIF
      ELSEIF( INDX.EQ.1 ) THEN
        XLSX = SX/RHOLX
!
!---  Newton iteration to determine the brine aqueous mass
!     fraction and aqueous density  ---
!
  100   CONTINUE
        CALL BRND( TX,XLSX,RHOLX,RHLSX,DRHLSX )
        IF( RHLSX.GE.RHLSM ) THEN
          RHLSX = RHLSM
          DRHLSX = 0.D+0
        ENDIF
        F = XLSX*RHLSX - SX
        DF = (RHLSX + XLSX*DRHLSX)
        DXLSX = -F/DF
        XLSX = XLSX + DXLSX
        IF( ABS(DXLSX/XLSX).GE.1.D-8 ) GOTO 100
        CALL BRND( TX,XLSX,RHOLX,RHLSX,DRHLSX )
        RHLSX = MIN( RHLSX,RHLSM )
        SX = XLSX*RHLSX*PORDX*SLX
      ELSE
        XLSX = SX/(PORDX*SLX*RHOLX+SMALL)
!
!---  Newton iteration to determine the brine aqueous mass
!     fraction and aqueous density  ---
!
  110   CONTINUE
        CALL BRND( TX,XLSX,RHOLX,RHLSX,DRHLSX )
        IF( RHLSX.GE.RHLSM ) THEN
          RHLSX = RHLSM
          DRHLSX = 0.D+0
        ENDIF
        F = XLSX*PORDX*SLX*RHLSX - SX
        DF = PORDX*SLX*(RHLSX + XLSX*DRHLSX)
        DXLSX = -F/DF
        XLSX = XLSX + DXLSX
        IF( ABS(DXLSX/XLSX).GE.1.D-8 ) GOTO 110
        CALL BRND( TX,XLSX,RHOLX,RHLSX,DRHLSX )
        RHLSX = MIN( RHLSX,RHLSM )
        SX = XLSX*RHLSX*PORDX*SLX
      ENDIF
!
!---  Function for sodium chloride Henry's problem  ---
!
      IF( ISLC(22).EQ.3 ) THEN
        IF( INDX.EQ.4 ) THEN
          IF( XLSX.GT.XLSM ) THEN
            RHLSX = RHOLX/(1.D+0-XLSM)
          ELSE
            RHLSX = RHOLX/(1.D+0-XLSX)
          ENDIF
          SX = RHLSX*XLSX*PORDX*SLX
        ELSEIF( INDX.EQ.3 ) THEN
          XLSX = SX
          IF( XLSX.GT.XLSM ) THEN
            RHLSX = RHOLX/(1.D+0-XLSM)
          ELSE
            RHLSX = RHOLX/(1.D+0-XLSX)
          ENDIF
          SX = RHLSX*XLSX*PORDX*SLX
        ELSEIF( INDX.EQ.2 ) THEN
          XLSX = SX*XLSM
          IF( XLSX.GT.XLSM ) THEN
            RHLSX = RHOLX/(1.D+0-XLSM)
          ELSE
            RHLSX = RHOLX/(1.D+0-XLSX)
          ENDIF
          SX = RHOLX*XLSX*PORDX*SLX
        ELSEIF( INDX.EQ.1 ) THEN
          RHLSM = RHOLX/(1.D+0-XLSM)
          RHLSX = MIN( RHLSM,(RHOLX + 0.6829D+0*SX) )
          XLSX = SX/(SX+RHOLX)
          SX = SX*PORDX*SLX
        ELSE
          RHLSM = RHOLX/(1.D+0-XLSM)
          RHLSX = MIN( RHLSM,(RHOLX + 0.6829D+0*SX/(PORDX*SLX+SMALL)) )
          XLSX = SX/(PORDX*SLX*RHLSX+SMALL)
        ENDIF
      ENDIF
!
!---  Reset aqueous density to brine density  ---
!
      RHOLX = RHLSX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNCN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRND( TX,XLSX,RHOLX,RHLSX,DRHLSX )
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
!     Compute brine density from brine solution mass fraction
!     liquid water density.
!
!     Leijnse, A.  1992.  Three-Dimensional Modeling of Coupled Flow
!     and Transport in Porous Media.  Ph.D. Dissertation, Department
!     of Civil Engineering and Geological Sciences, University of
!     Notre Dame, Notre Dame, Indiana.
!
!     Mahiuddin, S., and K. Ismail.  1996.  Temperature and
!     concentration dependence of the viscosity of aqueous sodium
!     nitrate and sodium thiosulfate electrolytic solutions.
!     Fluid Phase Equilibria, 123:231-243.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, August 1995.
!     Last Modified by MD White, Battelle, PNL, April 22, 1995.
!     brncn.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
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
      REAL*8 A(3,2),B(6,2)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A,B
      DATA A /1012.6D+0,628.46D+0,442.01D+0,
     &  1014.5D+0,775.87D+0,449.55D+0 /
      DATA B /0.48785D+0,-1.512D+0,33.036D+0,
     &  -198.96D+0,518.24D+0,-479.19D+0,
     &  0.43609D+0,0.80069D+0,-3.7273D+0,
     &  20.445D+0,-42.401D+0,28.528D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRND'
!
!---  Function for sodium chloride by A. Leijnse  ---
!
      IF( ISLC(22).EQ.0 ) THEN
        RHLSX = RHOLX*EXP(7.D-1*XLSX)
        DRHLSX = 7.D-1*RHOLX*EXP(7.D-1*XLSX)
!
!---  Function for sodium nitrate by S. Mahiuddin and K. Ismail  ---
!
      ELSEIF( ISLC(22).EQ.1 ) THEN
        AX = A(1,1)
        DAX = 0.D+0
        DO 100 N = 1,2
          AX = AX + A(N,1)*(XLSX**N)
          DAX = DAX + A(N,1)*(XLSX**(N-1))
  100   CONTINUE
        BX = B(1,1)
        DBX = 0.D+0
        DO 110 N = 1,5
          BX = BX + B(N,1)*(XLSX**N)
          DBX = DBX + B(N,1)*(XLSX**(N-1))
  110   CONTINUE
        RHLSX = AX + TX*BX
        DRHLSX = DAX + TX*DBX
!
!---  Function for sodium thiosulfate by S. Mahiuddin and K. Ismail  ---
!
      ELSEIF( ISLC(22).EQ.2 ) THEN
        AX = A(1,2)
        DAX = 0.D+0
        DO 200 N = 1,2
          AX = AX + A(N,2)*(XLSX**N)
          DAX = DAX + A(N,2)*(XLSX**(N-1))
  200   CONTINUE
        BX = B(1,2)
        DBX = 0.D+0
        DO 210 N = 1,5
          BX = BX + B(N,2)*(XLSX**N)
          DBX = DBX + B(N,2)*(XLSX**(N-1))
  210   CONTINUE
        RHLSX = AX + TX*BX
        DRHLSX = DAX + TX*DBX
!
!---  Function for sodium chloride Elder's problem  ---
!
      ELSEIF( ISLC(22).EQ.4 ) THEN
        RHLSX = RHOLX + 200.D+0*(XLSX/0.2605D+0)
        DRHLSX = 200.D+0/0.2605D+0
      ENDIF
!
!---  Function for Hanford tank SX109  ---
!
!      RHLSX = RHOLX*EXP(2.0513D+0*XLSX)
!      DRHLSX = 2.0513D+0*RHOLX*EXP(2.0513D+0*XLSX)
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRND group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNH( TX,XLWX,XLAX,XLSX,HGAX,HLX )
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
!     Compute brine enthalpy (NaCl solutions) by the method of Pitzer.
!
!     Lewis, G.B. and M. Randall. 1961.  Thermodynamics.  McGraw-Hill
!     Book Company, New York. pp.381-392.
!
!     Pitzer, K.S., J.C. Peiper, and R.H. Busey.  1984.  Thermodynamic
!     Properties of Aqueous Sodium Chloride Solutions.
!     J. Phys. Chem. Ref. Data, Vol. 13, No. 1.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, November 1995.
!     Last Modified by MD White, Battelle, November 16, 1995.
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
      REAL*8 VAR(4)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNH'
!
!---  Compute the salt mole fraction and molality  ---
!
      TK = TX+TABS
      XMLSX = (XLSX/WTMS)/(XLSX/WTMS + XLAX/WTMA + XLWX/WTMW )
      XMLWX = (XLWX/WTMW)/(XLSX/WTMS + XLAX/WTMA + XLWX/WTMW )
      XMLAX = 1.D+0 - XMLWX - XMLSX
      XGLSX = 1.D+3*XMLSX/(WTMW - XMLSX*WTMW)
!
!---  Compute the infinite dilution enthalpy  ---
!
      HLSO = (77.734D+0-0.60371D+0*TK+0.0015662D+0*(TK**2) -
     &  1.3913D-06*(TK**3))*RCU*TK
!
!---  Compute the excess enthalpy  ---
!
      VAR(1) = -9.3340D+1 + 6.8661D+1*EXP(-3.4452D-1*XGLSX)
     &  + 1.6031D+1*EXP(-3.5850D+0*XGLSX)
      VAR(2) = 6.77174D-1 - 4.7246D-1*EXP(-3.6703D-1*XGLSX)
     &  - 1.2881D-1*EXP(-3.5899D+0*XGLSX)
      VAR(3) = -1.6263D-3 + 1.0999D-3*EXP(-3.9371D-1*XGLSX)
     &  - 3.3829D-4*EXP(-3.6424D+0*XGLSX)
      VAR(4) = 1.3749D-6 - 9.0136D-7*EXP(-4.1024D-1*XGLSX)
     &  - 3.0285D-7*EXP(-3.6789D+0*XGLSX)
      HLSEX = (VAR(1)+VAR(2)*TK+VAR(3)*(TK**2)+VAR(4)*(TK**3))*RCU*TK
!
!---  Compute the brine enthalpy  ---
!
      HLX = (XMLWX*WTMW*HLX + XMLSX*(HLSO+HLSEX+3.883D+3)
     &  + XMLAX*WTMA*HGAX)/(XMLWX*WTMW + XMLAX*WTMA + XMLSX*WTMS)
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNK( XLWX,XLAX,XLSX,THKLX )
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
!     Compute brine thermal conductivity (NaCl solutions) by the
!     method of Riedel.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids. pp 566-567.
!     McGraw-Hill, New York, New York

!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, November 1995.
!     Last Modified by MD White, Battelle, November 16, 1995.
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
      SUB_LOG(ISUB_LOG) = '/BRNK'
!
!---  Compute the salt mole fraction and molality  ---
!
      XMLSX = (XLSX/WTMS)/(XLSX/WTMS + XLAX/WTMA + XLWX/WTMW )
      XGLSX = 1.D+3*XMLSX/(WTMW - XMLSX*WTMW)
!
!---  Compute the brine (NaCl solution) thermal conductivity at
!     reference temperature (20C) conditions and
!     correct to other temperatures  ---
!
      THKLR = 6.034709854D-01 + XGLSX*(-5.466D-5)
      THKLX = THKLR*THKLX/6.034709854D-01
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNK group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNOEC( N,SLX,RKLX,RKLAX,PORDX,TORLX,OECX )
!
!-------------------------Notice---------------------------------------!
!
!     This software is being made available to you for testing and
!     evaluation purposes only, solely for use in performance of work
!     under contracts with the U.S. Department of Energy or other
!     Federal Government agencies.  This software is a pre-release test
!     version which has not yet been authenticated and cleared for
!     publication.  Adherence to this notice may be necessary for the
!     author, Battelle Memorial Institute, to successfully assert
!     copyright in and commercialize this software.  This software
!     is not intended for duplication or distribution to third parties
!     without the permission of the Manager of Software Products at
!     Pacific Northwest Laboratory, Richland, Washington, 99352.
!
!----------------------Description-------------------------------------!
!
!     Compute the osmotic efficiency coefficient.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January 15, 1995.
!     Last Modified by MD White, Battelle, PNL, January 17, 1995.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNOEC'
!
!---  Compute the osmotic efficiency coefficient  ---
!
      IF( PERM(1,IZ(N)).GT.0.D+0 ) THEN
        IF( PERM(2,IZ(N)).GT.0.D+0 ) THEN
          IF( PERM(3,IZ(N)).GT.0.D+0 ) THEN
            PERMX = (PERM(1,IZ(N))*PERM(2,IZ(N))*PERM(3,IZ(N)))
     &        **(THIRD)
          ELSE
            PERMX = SQRT(PERM(1,IZ(N))*PERM(2,IZ(N)))
          ENDIF
        ELSEIF( PERM(3,IZ(N)).GT.0.D+0 ) THEN
          PERMX = SQRT(PERM(1,IZ(N))*PERM(3,IZ(N)))
        ELSE
          PERMX = PERM(1,IZ(N))
        ENDIF
      ELSE
        IF( PERM(2,IZ(N)).GT.0.D+0 ) THEN
          IF( PERM(3,IZ(N)).GT.0.D+0 ) THEN
            PERMX = SQRT(PERM(2,IZ(N))*PERM(3,IZ(N)))
          ELSE
            PERMX = PERM(2,IZ(N))
          ENDIF
        ELSEIF( PERM(3,IZ(N)).GT.0.D+0 ) THEN
          PERMX = PERM(3,IZ(N))
        ENDIF
      ENDIF
      RKLMX = (RKLX(1)*RKLX(2)*RKLX(3))**(THIRD)
      AWFT = SQRT(3.D+0*RKLMX*PERMX/(TORLX*SLX*PORDX))
      OECX = MAX( (AWFT-1.5D-10),SMALL )
      OECX = MIN( 1.32D-10/OECX,1.D+0 )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNOEC group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNS( TX,XLSX )
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
!     Compute the brine solubility in terms of salt mass fraction as
!     a function of temperature.
!
!     Leijnse, A.  1992.  Three-Dimensional Modeling of Coupled Flow
!     and Transport in Porous Media.  Ph.D. Dissertation, Department
!     of Civil Engineering and Geological Sciences, University of
!     Notre Dame, Notre Dame, Indiana.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, August 1995.
!     Last Modified by MD White, Battelle, PNL, April 22, 1995.
!     brncn.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/BRNS'
!
!---  Function for sodium nitrate by S. Mahiuddin and K. Ismail  ---
!
      IF( ISLC(22).EQ.1 ) THEN
        TK = TX + TABS
        XLSX = -6.91775D-9*(TK**2) + (2.1929D-3*(TK)) - 1.75953D-1
!
!---  Function for sodium thiosulfate by S. Mahiuddin and K. Ismail  ---
!
      ELSEIF( ISLC(22).EQ.2 ) THEN
        TK = TX + TABS
        XLSX = 4.5608D-5*(TK**2) + (-2.19878D-2*(TK)) + 2.93319D+0
!
!---  Function for sodium chloride by A. Leijnse  ---
!
      ELSE
        TR = (TX+2.7315D+2)/1.D+2
        XLSX = EXP( -1.30546D+0 + (-1.08431D-1)/TR 
     &    + (-9.25719D-1)*LOG(TR) + 3.43444D-1*TR )
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNS group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNCN( TX,SLX,PORDX,SX,RHOLX,XLSX,YLSX,INDX )
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
!     INDX = 1: salt aqueous concentration loaded in SX
!     INDX = 2: salt relative saturation loaded in SX
!     INDX = 3: salt mass fraction loaded in SX
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
      USE FDVS
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
!---  Set pure water density for Elder's problem  ---
!
      IF( ISALT.EQ.4 ) RHOLX = 1.D+3
!
!---  Compute maximum dissolved salt mass fraction ---
!
      CALL BRNS( TX,XLSM )
!
!---  Compute maximum dissolved salt mass fraction ---
!
      RHLSM = RHOLX
      RHLSX = RHOLX
      CALL BRND( TX,XLSM,RHLSM,DRHLSM )
!
!---  Compute aqueous density and salt volumetric concentration from
!     the salt aqueous mass fraction  ---
!
      IF( INDX.EQ.3 ) THEN
        YLSX = SX
        XLSX = MIN( YLSX,XLSM )
        RHLSX = RHOLX
        CALL BRND( TX,XLSX,RHLSX,DRHLSX )
        SX = XLSX*RHLSX*PORDX*SLX
!
!---  Compute aqueous density and salt volumetric concentration from
!     the salt relative saturation  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
        YLSX = SX*XLSM
        XLSX = MIN( YLSX,XLSM )
        RHLSX = RHOLX
        CALL BRND( TX,XLSX,RHLSX,DRHLSX )
        SX = XLSX*RHLSX*PORDX*SLX
!
!---  Compute aqueous density and salt volumetric concentration from
!     the salt aqueous density  ---
!
      ELSEIF( INDX.EQ.1 ) THEN
!
!---    Function for sodium chloride Henry's problem  ---
!
        IF( ISALT.EQ.3 ) THEN
!          RHLSX = RHOLX + 0.707314D+0*SX
!          SX = (RHLSX-RHOLX)/0.707314D+0
          RHLSX = RHOLX + 0.6829D+0*SX
          SX = (RHLSX-RHOLX)/0.6829D+0
          YLSX = SX/(RHLSX-SX)
          RHLSX = MIN( RHLSX,RHLSM )
!          SX = (RHLSX-RHOLX)/0.707314D+0
          SX = (RHLSX-RHOLX)/0.6829D+0
          XLSX = SX/(RHLSX-SX)
        ELSE
          YLSX = SX/RHOLX
  100     CONTINUE
          RHLSX = RHOLX
          CALL BRND( TX,YLSX,RHLSX,DRHLSX )
          IF( RHLSX.GE.RHLSM ) THEN
            RHLSX = RHLSM
            DRHLSX = 0.D+0
          ENDIF
          F = YLSX*RHLSX - SX
          DF = (RHLSX + YLSX*DRHLSX)
          DYLSX = -F/DF
          YLSX = YLSX + DYLSX
          IF( ABS(DYLSX/XLSM).GE.1.D-8 ) GOTO 100
          RHLSX = RHOLX
          CALL BRND( TX,YLSX,RHLSX,DRHLSX )
          XLSX = MIN( YLSX,XLSM )
          RHLSX = MIN( RHLSX,RHLSM )
          SX = XLSX*RHLSX*PORDX*SLX
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
      SUBROUTINE BRND( TX,XLSX,RHLSX,DRHLSX )
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
!     brn_prp.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVS
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
      IF( ISALT.EQ.0 ) THEN
        RHLSX = RHLSX*EXP(7.D-1*XLSX)
        DRHLSX = 7.D-1*RHLSX*EXP(7.D-1*XLSX)
!
!---  Function for sodium nitrate by S. Mahiuddin and K. Ismail  ---
!
      ELSEIF( ISALT.EQ.1 ) THEN
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
!
!---    ASME Pure-Water Density Correction  ---
!
        RHWX = A(1,1) + TX*B(1,1)
        RHLSX = (AX + TX*BX)*RHLSX/RHWX
        DRHLSX = (DAX + TX*DBX)*RHLSX/RHWX
!
!---  Function for sodium thiosulfate by S. Mahiuddin and K. Ismail  ---
!
      ELSEIF( ISALT.EQ.2 ) THEN
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
!
!---    ASME Pure-Water Density Correction  ---
!
        RHWX = A(1,2) + TX*B(1,2)
        RHLSX = (AX + TX*BX)*RHLSX/RHWX
        DRHLSX = (DAX + TX*DBX)*RHLSX/RHWX
!
!---  Function for sodium chloride Henry's problem  ---
!
      ELSEIF( ISALT.EQ.3 ) THEN
!        RHLSX = RHLSX/(1.D+0 - 0.707314D+0*XLSX/(1.D+0+XLSX))
!        DRHLSX = -RHLSX*0.707314D+0/
!     &    ((1.D+0 + XLSX - 0.707314D+0*XLSX)**2)
         RHLSX = RHLSX/(1.D+0 - 0.6829D+0*XLSX)
         DRHLSX = 0.6829D+0*RHLSX/(1.D+0 - 0.6829D+0*XLSX)
!
!---  Function for sodium chloride Elder's problem  ---
!
      ELSEIF( ISALT.EQ.4 ) THEN
        RHLSX = RHLSX + 200.D+0*(XLSX/0.2605D+0)
        DRHLSX = 200.D+0/0.2605D+0
      ENDIF
!
!---  Function for Hanford tank SX109  ---
!
!      RHLSX = RHLSX*EXP(2.0513D+0*XLSX)
!      DRHLSX = 2.0513D+0*RHLSX*EXP(2.0513D+0*XLSX)
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
      SUBROUTINE BRNOEC( N,SLX,RKLX,PORDX,TORLX,OECX )
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
!     brn_prp.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVS
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
      IF( ISALT.EQ.1 ) THEN
        TK = TX + TABS
        XLSX = -6.91775D-9*(TK**2) + (2.1929D-3*(TK)) - 1.75953D-1
!
!---  Function for sodium thiosulfate by S. Mahiuddin and K. Ismail  ---
!
      ELSEIF( ISALT.EQ.2 ) THEN
        TK = TX + TABS
        XLSX = 4.5608D-5*(TK**2) + (-2.19878D-2*(TK)) + 2.93319D+0
!
!---  Function for sodium chloride by A. Leijnse  ---
!
      ELSE
        TR = (TX+2.7315D+2)/1.D+2
        XLSX = EXP( -1.30546D+0 + (-1.08431D-1)/TR +
     &    (-9.25719D-1)*LOG(TR) + 3.43444D-1*TR )
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNS group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNLQV( TX,XLSX,VISLX )
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
!     Compute aqueous-brine (sodium chloride) liquid viscosity.
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
!     Written by MD White, Battelle, PNL, July 1995.
!     Last Modified by MD White, Battelle, 13 Sept 2001.
!     brn_prp.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVS
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 A_C(3,2),B_C(6,2),T_C(3,2)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A_C,B_C,T_C
      DATA A_C /9.91828D-4,2.98772D-3,4.22179D-3,
     &  1.00716D-3,3.65659D-3,1.76215D-3 /
      DATA B_C /642.237D+0,-689.63D+0,8273.85D+0,
     &  -51214.8D+0,125889.D+0,-104430.D+0,
     &  643.915D+0,121.295D+0,-3997.75D+0,
     &  19616.9D+0,-34793.9D+0,21344.6D+0 /
      DATA T_C /134.711D+0,15.6163D+0,29.475D+0,
     &  134.045D+0,16.4398D+0,136.606D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNLQV'
!
!---  Brine viscosity computed from an empirical relation by Leijnse ---
!
      IF( ISALT.EQ.0 ) THEN
        VISLX = VISLX*(1.D+0 + 1.85D+0*XLSX + 4.1D+0*(XLSX**2)
     &    + 4.45D+1*(XLSX**3))
!
!---  Function for sodium nitrate by S. Mahiuddin and K. Ismail  ---
!
      ELSEIF( ISALT.EQ.1 ) THEN
        A_CX = A_C(1,1)
        DO 100 N = 1,2
          A_CX = A_CX + A_C(N,1)*(XLSX**N)
  100   CONTINUE
        B_CX = B_C(1,1)
        DO 110 N = 1,5
          B_CX = B_CX + B_C(N,1)*(XLSX**N)
  110   CONTINUE
        T_CX = T_C(1,1)
        DO 120 N = 1,2
          T_CX = T_CX + T_C(N,1)*(XLSX**N)
  120   CONTINUE
        VISLX = 1.D-3*A_CX*SQRT(TX+TABS)*EXP(B_CX/(TX+TABS-T_CX))
!
!---  Function for sodium thiosulfate by S. Mahiuddin and K. Ismail  ---
!
      ELSEIF( ISALT.EQ.2 ) THEN
        A_CX = A_C(1,2)
        DO 200 N = 1,2
          A_CX = A_CX + A_C(N,2)*(XLSX**N)
  200   CONTINUE
        B_CX = B_C(1,2)
        DO 210 N = 1,5
          B_CX = B_CX + B_C(N,2)*(XLSX**N)
  210   CONTINUE
        T_CX = T_C(1,2)
        DO 220 N = 1,2
          T_CX = T_CX + T_C(N,2)*(XLSX**N)
  220   CONTINUE
        VISLX = 1.D-3*A_CX*SQRT(TX+TABS)*EXP(B_CX/(TX+TABS-T_CX))
!
!---  Function for sodium chloride Henry's problem  ---
!
      ELSEIF( ISALT.EQ.3 ) THEN
        VISLX = 1.0176D-3
!
!---  Function for sodium chloride Elder's problem  ---
!
      ELSEIF( ISALT.EQ.4 ) THEN
        VISLX = 1.0D-3
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNLQV group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNDBH( TX,PX,APHI )
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
!     Compute the Debye-Huckel parameter as a function of temperature
!     and pressure. Below 373.15 K (100 oC), the parameters
!     apply to a standard pressure of 0.1 MPa (1 bar). Above 373.15 K,
!     they apply to water under its own vapor presssure.
!
!     Bradley, D.J. and K.S.Pitzer. 1983. Thermodynamics of electrolytes
!     XII. Dielectric proprties of water and Debye-Huckel parameters to
!     350 oC and 1kbar. J. Phys. Chem 83:1599,1979;87:3798.
!
!     ELC = electronic charge,
!     ANO = Avagadro's number, molecules/mole
!     BZM = Boltzmann's constant,
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, November 1995.
!     Last Modified by MD White, Battelle, November 16, 1995.
!     Last Modified by AL Ward, Battelle, July 23, 1999.
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
!----------------------Data Statements---------------------------------!
!
      SAVE ELC,ANO,BZM,EPO
      DATA ELC,ANO,BZM / 1.60210D-19,6.02252D+23,1.38054D-23 /
      DATA EPO / 8.854D-12 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNDBH'
!
!---  Convert input temperature to absolute temperature ---
!
        TK = TX+TABS
!
!---  Compute the liquid water density as a function of P and T ---
!
        CALL WATLQD( TX,PX,RHOLX )
!
!---  Compute relative permittivity of pore water at temperature T---
!
        ER = EPR(TX)
!
!---  Compute Debye-Huckel parameter for the osmotic coefficient ---
!
        VAR1 = SQRT(2.D+0*GPI*ANO*RHOLX)
        VAR2 = ((ELC**2)/(4.D+0*GPI*EPO*ER*BZM*TK))**1.5D+0
        APHI = VAR1*VAR2/3.D+0
!
!---  End of BRNDBH group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!
!----------------------Subroutine--------------------------------------!
!
      FUNCTION EPR( TX )
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
!     Compute the the relative dielectric constant as a function of
!     temperature assuming a pressure of (</=100 MPa, or 1 kbar). At
!     T< 373.15 K (100 oC), the parameters apply to a standard pressure
!     of 0.1 MPa (1 bar). At T >373.15 K, they apply to water under its
!     own vapor presssure.
!
!     Bradley, D.J. and K.S.Pitzer. 1983. Thermodynamics of electrolytes
!     XII. Dielectric proprties of water and Debye-Huckel parameters to
!     350 oC and 1kbar. J. Phys. Chem 83:1599,1979;87:3798.
!
!----------------------Authors-----------------------------------------!
!
!     Written by AL Ward, Battelle, PNNL, July 23, 1999.
!     Last Modified by AL Ward, Battelle, July 24, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 U(9)
!
!----------------------Data Statements---------------------------------!
!
      SAVE U
      DATA U / 3.4279D+2, -5.0866D-3, 9.4690D-7, -2.0525D+0, 3.1159D+3,
     &  -1.8289D+2, -8.0325D+3, 4.2142D+6, 2.1417D+0 /
!
!----------------------Executable Lines--------------------------------!
!
!
      TK = TX+TABS
      EPK = U(1)* EXP(U(2)*TK+U(3)*TK*TK)
      CX = U(4)+U(5)/(U(6)+TK)
      BX = U(7)+U(8)/TK+U(9)*TK
      VPX = PATM
      IF( TK.GT.3.73149D+2 ) CALL WATSP( TX,VPX )
      EPR = EPK+CX*LOG((BX+VPX*1.D-5)/(BX+1.D+3))
!
!---  End of EPR group  ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNOSM ( TX,XLSX,APHI,BOC,SOC,PVWX,POSX )
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
!     Compute the the osmotic coefficient as a function of salt mass
!     fraction and temperature using the Pitzer equation
!
!     Pitzer, K J.C. 1991. Ion interaction approach: Theory and data
!     correlation. In Activity coefficeints in electrolyte solutions,
!     2nd Ed.., CRC Press: Boston.
!
!     Li, Z-B,Y-G Li, and J-F Lu. 1999. Surface Tension Model for
!     Concentrated Electrolyte Aqueous Solutions by the
!     Pitzer Equation. Ind. Eng. Chem. Res. 38, 1133-1139
!
!     PTPS(1) = B0
!     PTPS(2) = B1
!     PTPS(3) = B2
!     PTPS(4) = CPHI
!     PTPS(5) = GOC
!     IPTPS(1) = IVPOS
!     IPTPS(2) = IVNEG
!     IPTPS(3) = IZPOS
!     IPTPS(4) = IZNEG
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, November 1995.
!     Last Modified by MD White, Battelle, November 16, 1995.
!     Last Modified by AL Ward, Battelle, July 23, 1999.
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
!----------------------Data Statements---------------------------------!
!
      SAVE SOMV
      DATA SOMV / 1.805D-5 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNOSM'
      TK = TX+TABS
!
!---  Compute stoichiometric and charge balance terms ---
!
      ZPOS = REAL(IPTPS(3))
      ZNEG = REAL(IPTPS(4))
      VPOS = REAL(IPTPS(1))
      VNEG = REAL(IPTPS(2))
      ZMX = ABS(ZPOS*ZNEG)
      VPN = VPOS+VNEG
      VMX = VPOS*VNEG
      OCB = 1.2D+0
!
!---  Compute mole fraction and molality in bulk and surface phases ---
!
      XMLSX = XLSX/WTMS/(XLSX/WTMS + (1.D+0-XLSX)/WTMW)
      XGLSX = 1.D+3*XMLSX/(WTMW - XMLSX*WTMW)
      SXGLSX = PTPS(5)*XGLSX
!
!---  Compute the bulk phase ionic strength and its square root ---
!
      BSI = 5.D-1*(XGLSX*ZPOS*ZPOS + XGLSX*ZNEG*ZNEG)
      DBSI = SQRT(BSI)
      SSI = 0.5D+0*(SXGLSX*ZPOS*ZPOS + SXGLSX*ZNEG*ZNEG)
      DSSI = SQRT(SSI)
!
!---  Compute bulk phase osmotic coefficent ---
!
      IF( XGLSX.GT.EPSL ) THEN
        IF( IPTPS(1).EQ.1 .AND. IPTPS(2).GE.1 .AND. IPTPS(2).LE.3 ) THEN
          A1 = 2.D+0
          BOC = 1.D+0 + ZMX*APHI*DBSI/(1.D+0+OCB*DBSI) +
     &      2.D+0*XGLSX*VMX*(PTPS(1)+PTPS(2)*EXP(-A1*DBSI))/VPN +
     &      2.D+0*(XGLSX**2)*(VMX**1.5D+0)*PTPS(4)/VPN
        ELSEIF( IPTPS(1).GE.2 .AND. IPTPS(2).GE.2 ) THEN
          A1 = 1.4D+0
          A2 = 1.2D+1
          BOC = 1.D+0 + ZMX*APHI*DBSI/(1.D+0+OCB*DBSI) +
     &      2.D+0*XGLSX*VMX*(PTPS(1)+PTPS(2)*EXP(-A1*DBSI))/VPN +
     &      2.D+0*(XGLSX**2)*(VMX**1.5D+0)*PTPS(4)/VPN +
     &      PTPS(3)*EXP(-A2*DBSI)
        END IF
      ELSE
        BOC = 0.D+0
      END IF
!
!---  Compute surface phase osmotic coefficent ---
!
      IF( XGLSX.GT.EPSL ) THEN
        IF( IPTPS(1).EQ.1 .AND. IPTPS(2).GE.1 .AND. IPTPS(2).LE.3 ) THEN
          A1 = 2.D+0
          SOC = 1.D+0 + ZMX*APHI*DSSI/(1.D+0+OCB*DSSI) +
     &      2.D+0*SXGLSX*VMX*(PTPS(1)+PTPS(2)*EXP(-A1*DSSI))/VPN +
     &      2.D+0*(SXGLSX**2)*(VMX**1.5D+0)*PTPS(4)/VPN
        ELSEIF( IPTPS(1).GE.2 .AND. IPTPS(2).GE.2 )THEN
          A1 = 1.4D+0
          A2 = 1.2D+1
          SOC = 1.D+0 + ZMX*APHI*DSSI/(1.D+0+OCB*DSSI) +
     &      2.D+0*SXGLSX*VMX*(PTPS(1)+PTPS(2)*EXP(-A1*DSSI))/VPN +
     &      2.D+0*(SXGLSX**2)*(VMX**1.5D+0)*PTPS(4)/VPN +
     &      PTPS(3)*EXP(-A2*DSSI)
        END IF
      ELSE
        SOC = 0.D+0
      END IF
!
!---  Compute overall osmotic coefficient
!
      OC = BOC - PTPS(5)*SOC
!
!---  Compute the osmotic pressure of the pore water ---
!
      POSX = 0.D+0
      IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 ) THEN
        POSX = -OC*(RCU/1.D+3)*TK*VPN*XGLSX*WTMW/(SOMV*1.D+3)
      END IF
!
!---  Compute the vapor pressure of the pore water ---
!
      AAC = EXP(-3.D+0*OC*WTMW*XGLSX/1.D+3 )
      PVWX =PVWX*AAC
!
!---  End of BRNOSM group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNSTN( TX,XLSX,SOC,BOC,ST )
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
!     Compute the pore water surface tension as a function of salt mass
!     fraction and temperature using the Pitzer equation
!
!     Adamson, A.W. 1976. The physical Chemistry of Surfaces. 3rd. Ed
!     Wiley--Interscience, NY.
!
!     SOMV = water molar volume, m^3/mol
!     STWX = water-air surface tension, N/m
!     PTPS(5) = GOC
!
!----------------------Authors-----------------------------------------!
!
!     Written by AL Ward, Battelle, PNL, July 17, 1999.
!     Last Modified by AL Ward, Battelle, PNL, July 19, 1999.
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
!----------------------Data Statements---------------------------------!
!
      SAVE SOMO,SOMV,STWX,ANO
      DATA SOMO,SOMV,STWX,ANO / 55.50837D+1,1.805D-5,7.2D-2,6.0225D+23 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNSTN'
!
!---  Compute the molar surface area of pure water ---
!
      SAMW = (ANO**(THIRD))*(SOMV**(2.D+0/3.D+0))
!
!---  Compute mole fraction and molality from wt % in the bulk phase ---
!
      XMLSX = (XLSX/WTMS)/(XLSX/WTMS + (1.D+0-XLSX)/WTMW)
      XGLSX = 1.D+3*XMLSX/(WTMW - XMLSX*WTMW)
!
!---  Compute the van't Hoff coefficient ---
!
      TK = TX + TABS
      VPOS = REAL(IPTPS(1))
      VNEG = REAL(IPTPS(2))
      VHOF =((VPOS+VNEG)*(RCU/1.D+3)*TK*XGLSX)/(SOMO*SAMW)
!
!---  Compute the contribution of the bulk liquid phase ---
!
      STB = VHOF*BOC
!
!---  Compute the contribution of the surface phase ---
!
      STS = VHOF*PTPS(5)*SOC
!
!---  Compute the surface tension of solution ---
!
      ST = STWX + (STB - STS)
!
!---  End of BRNSTN group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRN_VPR( PVWX,TX,XLSX )
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
!     Vapor pressure of aqueous solutions of sodium chloride,
!     sodium nitrate, or sodium thiosulfate.
!
!     Apelblat, A., and E. Korin.  1998.  The vapour pressure of
!     saturated aqueous solutions of sodium chloride, sodium bromide,
!     sodium nitrate, sodium nitrite, potassium iodate, and
!     rubidium chloride at temperatures from 227 K to 323 K.
!     Journal of Chemical Thermodynamics, 30:59-71.
!
!     Apelblat, A., and E. Korin.  1998.  The vapour pressure of
!     saturated aqueous solutions of ammonium iodide, potassium iodide,
!     potassium nitrate, strontium chloride, lithium sulphate,
!     sodium thiosulfate, magnesium nitrate, and
!     uranyl nitrate at temperatures from 227 K to 323 K.
!     Journal of Chemical Thermodynamics, 30:459-471.
!
!     ISALT:   0    Sodium Chloride (NaCl)
!              1    Sodium Nitrate (NaNO3)
!              2    Sodium Thiosulfate (Na2S2O3)
!              3    Sodium Chloride (NaCl) (Henry's Problem Functions)
!              4    Sodium Chloride (NaCl) (Elder's Problem Functions)
!              5    Sodium Bromide (NaBr)
!              6    Sodium Nitrite (NaNO2)
!              7    Potassium Iodate (KIO3)
!              8    Rubidium Chloride (RbCl)
!              9    Ammonium Iodide (NH4I)
!             10    Potassium Iodide (KI)
!             11    Potassum Nitrate (KNO2)
!             12    Strontium Chloride (SrCl2)
!             13    Lithium Sulphate (LiSO4)
!             14    Magnesium Nitrate (Mg(NO3)2)
!             15    Uranyl Nitrate UO2(NO3)2
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 September 2001.
!     Last Modified by MD White, PNNL, 18 September 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVS
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 A_C(9),B_C(9),C_C(9)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A_C,B_C,C_C
      DATA A_C / 82.503D+0, 161.958D+0, 139.720D+0, 82.503D+0,
     &  82.503D+0, 105.921D+0, 147.920D+0, 222.615D+0, 137.885D+0 /
      DATA B_C / -8160.1D+0, -11680.3D+0, -10189.0D+0, -8160.1D+0,
     &  -8160.1D+0, -9125.2D+0, -11161.2D+0, -14420.0D+0, -10680.9D+0 /
      DATA C_C / -9.523D+0, -21.401D+0, -18.360D+0, -9.523D+0,
     &  -9.523D+0, -13.110D+0, -19.280D+0, -30.392D+0, -17.767D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRN_VPR'
!
!---  Pure water vapor pressure, saturated aqueous solution
!     mole fractions, and molar saturation index
!
      CALL WATSP( TX,PSWX )
      CALL BRNS( TX,XLSMX )
      WTML = 1.D+0/(XLSMX/WTMS + (1.D+0-XLSMX)/WTMW)
      XMLSMX = XLSMX*WTML/WTMS
      WTML = 1.D+0/(XLSX/WTMS + (1.D+0-XLSX)/WTMW)
      XMLSX = XLSX*WTML/WTMS
      SMLSX = XMLSX/XMLSMX
      TK = TX + TABS
!
!---  Function by A. Apelblat and E. Korin  ---
!
      PSLX = 1.D+3*EXP( A_C(ISALT+1) + B_C(ISALT+1)/TK +
     &    C_C(ISALT+1)*LOG(TK) )
!
!---  Linear interpolation between vapor pressure of pure water and
!     vapor pressure of saturated salt solution using linear
!     interpolation with salt mole fraction  ---
!
      PVWX = SMLSX*PSLX + (1.D+0-SMLSX)*PSWX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRN_VPR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRN_SFT( TX,XLSX,SFTX )
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
!     Surface tension of aqueous solutions of sodium chloride,
!     sodium nitrate, or sodium thiosulfate.
!
!     Abramzon, A.A., and R.D. Gaukhberg.  1993.  Surface tension of
!     salt solutions.  Russian Journal of Applied Chemistry,
!     66(6):1139-1146.
!
!     Lide, D.R., and H.V. Kehiaian.  1994.  CRC Handbook of
!     Thermophysical and Thermochemical Data.  CRC Press, Boca Raton,
!     Florida.
!
!     ISALT:   0    Sodium Chloride (NaCl)
!              1    Sodium Nitrate (NaNO3)
!              2    Sodium Thiosulfate (Na2S2O3)
!              3    Sodium Chloride (NaCl) (Henry's Problem Functions)
!              4    Sodium Chloride (NaCl) (Elder's Problem Functions)
!              5    Sodium Bromide (NaBr)
!              6    Sodium Nitrite (NaNO2)
!              7    Potassium Iodate (KIO3)
!              8    Rubidium Chloride (RbCl)
!              9    Ammonium Iodide (NH4I)
!             10    Potassium Iodide (KI)
!             11    Potassum Nitrate (KNO2)
!             12    Strontium Chloride (SrCl2)
!             13    Lithium Sulphate (LiSO4)
!             14    Magnesium Nitrate (Mg(NO3)2)
!             15    Uranyl Nitrate UO2(NO3)2
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 September 2001.
!     Last Modified by MD White, PNNL, 19 September 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVS
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 A_C(9)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A_C
      DATA A_C / 1.57D+0, 1.23D+0, 15.0D+0, 1.57D+0,
     &  1.57D+0, 1.26D+0, 0.00D+0, 0.00D+0, 1.57D+0 /
!      DATA A_C / 1.57D+0, 1.23D+0, 2.90D+0, 1.57D+0,
!     &  1.57D+0, 1.26D+0, 0.00D+0, 0.00D+0, 1.57D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRN_SFT'
!
!---  Pure water vapor surface tension as a function of temperature
!     by D. Lide and H. Kehiaian
!
      SFTWX = 1.D-3*(75.6592D+0 - 1.40959D-1*TX - 2.66317D-4*(TX**2))
!
!---  Convert mass fraction to molality (mole salt/kg water)
!
      CMOLX = XLSX/(WTMS*(1.D+0-XLSX))
!
!---  Function by A. Abramzon and R. Gaukhberg  ---
!
      SFTX = SFTWX + A_C(ISALT+1)*CMOLX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRN_SFT group  ---
!
      RETURN
      END


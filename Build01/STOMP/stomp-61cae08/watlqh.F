!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATLQH( TX,PX,H )
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
!     Calculate the subcooled or saturated enthalpy, as a function of
!     temperature and pressure per the steam table equations
!     as given by the 1967 International Formulation Committee:
!     Formulation for Industrial Use.
!
!     Thermodynamic and Transport Properties of Steam.
!     1967. ASME Steam Tables.
!     The American Society of Mechanical Engineers.
!     United Engineering Center, 345 East 47th Street, New York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < T > 364.0 !
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 6, 1992.
!     watlqh.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 A(23),SA(12),B(4)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A,SA,B
      DATA A /6.824687741D+3,-5.422063673D+2,-2.096666205D+4,
     &3.941286787D+4,-6.733277739D+4,9.902381028D+4,-1.093911774D+5,
     &8.590841667D+4,-4.511168742D+4,1.418138926D+4,-2.017271113D+3,
     &7.982692717D+0,-2.616571843D-2,1.522411790D-3,2.284279054D-2,
     &2.421647003D+2,1.269716088D-10,2.074838328D-7,2.174020350D-8,
     &1.105710498D-9,1.293441934D+1,1.308119072D-5,6.047626338D-14/
      DATA SA /8.438375405D-1,5.362162162D-4,1.720,7.342278489D-2,
     &4.975858870D-2,6.537154300D-1,1.15D-6,1.1508D-5,1.4188D-1,
     &7.002753165D+0,2.995284926D-4,2.040D-1/
      DATA B / -6.36443D+5,1.14468D+2,3.64104D+0,-3.3561464D+5 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATLQH'
!
!---  Liquid water enthalpy extended below freezing point  ---
!
      TR = MIN( (MAX(1.D-2,TX)+TABS)/TCRW,1.D+0 )
      PR = PX/PCRW
      CY = 1.0D+0 - SA(1)*TR*TR-SA(2)/(TR**6)
      CDY = -2.D+0*SA(1)*TR + 6.D+0*SA(2)/(TR**7)
      CZ = CY +
     &  SQRT(MAX( ZERO,(SA(3)*CY*CY-2.D+0*SA(4)*TR+2.D+0*SA(5)*PR)))
      HX = TR*A(1)
      DO 100 I = 1,10
        HX = HX - (I-2)*A(I+1)*(TR**(I-1))
  100 CONTINUE
      HX = HX + A(12)*(CZ*(1.7D+1*(CZ/2.9D+1 - CY/1.2D+1) +
     &  4.1667D-1*TR*CDY) + SA(4)*TR - (SA(3) - 1.D+0)*TR*CY*CDY)/
     &  (CZ**2.9412D-1)
      HX = HX + (A(13) - A(15)*TR*TR
     &  + A(16)*(9.D+0*TR + SA(6))*((SA(6) - TR)**9)
     &  + A(17)*(2.D+1*(TR**19) + SA(7))/((SA(7) + (TR**19))**2))*PR
      HX = HX - ((1.2D+1*(TR**11) + SA(8))/((SA(8) + (TR**11))**2))*
     &  (A(18)*PR + A(19)*PR*PR + A(20)*PR*PR*PR)
      HX = HX + A(21)*(TR**18)*(1.7D+1*SA(9) + 1.9D+1*TR*TR)*
     &  ((1.D+0/((SA(10) + PR)**3) + SA(11)*PR))
      HX = HX + A(22)*SA(12)*PR*PR*PR + 2.1D+1*A(23)*(PR**4)/(TR**20)
      H = HX*PCRW*VCRW*1.D-3/WTMW
!
!---  Extend the enthalpy function for freezing conditions
!     by incorporating the latent heat of fusion over the
!     temperature range -0.01 C to -1.0 C and using a polynomial
!     fit for the enthalpy of ice, referenced to -1.0 C  ---
!
      IF( TX.LT.1.D-2 ) THEN
        H = H + (MAX(-1.D+0,TX)-1.D-2)*333.7D+3/1.01D+0
      ENDIF
      IF( TX.LT.-1.D+0 ) THEN
        TK = TX + TABS
        H = H + B(1) + B(2)*TK + B(3)*TK**2 - B(4)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATLQH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_I( TX,RHOIX )
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
!     Calculate the density of ice as a function of temperature using a
!     using polynomial fit of data from ASHRAE (1977).
!
!     ASHRAE, 1977; Psychrometric Tables given by ASHRAE Handbook 1977
!     Fundamentals, American Society of Heating, Refrigerating and
!     Air-Conditioning Engineering, Inc., 345 East 47th Street, New
!     York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  -100.0 C < T <= 0.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, PNNL, September 1994.
!     Last Modified by MD White, PNNL, 29 November 2004.
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
      REAL*8 A(3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A
      DATA A / 1.06472D-3,1.23173D-8,3.0203D-10 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_I'
      TKX = TX+TABS
!
!---  Polynomial fit for specific volume  ---
!
      SPVI = A(1) + A(2)*TKX + A(3)*(TKX**2)
      RHOIX = 1.D+0/SPVI
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_I group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_I( TX,HIX )
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
!     Calculate the enthalpy of ice as a function of temperature using a
!     using polynomial fit of data from ASHRAE (1977) with a 
!     reference point of liquid water at 273.15 K.
!
!     ASHRAE, 1977; Psychrometric Tables given by ASHRAE Handbook 1977
!     Fundamentals, American Society of Heating, Refrigerating and
!     Air-Conditioning Engineering, Inc., 345 East 47th Street, New
!     York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  -100.0 C < T <= 0.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, PNNL, September 1994.
!     Last Modified by MD White, PNNL, 29 November 2004.
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
      REAL*8 A(3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A
      DATA A / -6.36443D+5,1.14468D+2,3.64104D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_I'
      TKX = TX+TABS
!
!---  Polynomial fit for reference point of liquid water
!     at 273.15 K (includes the heat of dissociation 333.7 kJ/kg)  ---
!
      HIX = A(1) + A(2)*TKX + A(3)*(TKX**2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_I group  ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ICE( TX,PLX,PGX,TFPX,PIX )
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
!     Compute ice freezing point and ice pressure as a function of
!     temperature, pressure, and dissolved aqueous components.
!
!     Pure water freezing point as a function of pressure from
!     Wagner W., A. Saul, and A. Pruss.  1994.  International equations
!     for the pressure along the melting and along the sublimation
!     curve of ordinary water substance.  Journal of Physical and 
!     Chemical Reference Data, 23:515-527.
!
!     Ice-aqueous radius of curvature from Brun, M., A. Lallemand,
!     J.F. Quinson, and C. Eyraud.  A new method for the simultaneous
!     determination of the size and shape of pores:
!     The thermoporometry.  Thermochimica Acta, 21:59-88.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 2 August 2005.
!     Last Modified by MD White, PNNL, 2 August 2005.
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
      SUB_LOG(ISUB_LOG) = '/ICE'
!
!---  Ice I freezing point (273.16 - 251.165 K) from Wagner et al., 
!     1994, using Newton iteration and a starting guess of 0 C ---
!
      TNX = 273.16D+0
      PNX = 6.11657D+2
      PX = MAX(PLX,PGX)+PATM
      TFPX = TABS
      THETAX = TFPX/TNX
   10 CONTINUE
      PIX = 1.D+0 - 0.626D+6*(1.D+0-(1.D+0/(THETAX**3))) + 
     &  0.197135D+6*(1.D+0-(THETAX**(21.2D+0)))
      DPIX = -1.878D+6/(THETAX**4) - 4.179262D+6*(THETAX**20.2D+0)
      FX = PIX*PNX-PX
      DFX = DPIX*PNX
      DTHETAX = -FX/DFX
      THETAX = THETAX + DTHETAX
      IF( ABS(DTHETAX).GT.1.D-12 ) GOTO 10
      TFPX = THETAX*TNX - TABS
!
!---  Ice-aqueous interfacial radius from Brun et al., 1977  ---
!
      IF( TX.LT.TFPX ) THEN
!        RADX = 1.D-9*((64.67D+0/(TFPX-TX)) + 0.57D+0)
        RADX = 1.D-10*((646.7D+0/(TFPX-TX)) - 2.3D+0)
        PIX = PLX + 2.D+0*0.0267D+0/RADX
      ELSE
        PIX = PLX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ICE group  ---
!
      RETURN
      END


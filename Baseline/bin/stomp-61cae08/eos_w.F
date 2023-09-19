!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE AIRGSD( TX,PX,RHOX )
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
!     Calculates component density of air with the ideal gas equation
!     of state.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     airgsd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/AIRGSD'
      RHOX = PX/(RCA*(TX+TABS))
      ISUB_LOG = ISUB_LOG-1
!
!---  End of AIRGSD group  ---
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION FNHGBL( NBS,NBE,M )
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
!     Aqueous hydraulic gradient boundary condition function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, June 1994.
!     Last Modified by MD White, Battelle, PNL, September 10, 1997.
!     Last Modified by MD White, Battelle, PNNL, July 12, 2000.
!     fnhgbl.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FNHGBL'
!
!---  Set starting-point pressures and temperatures  ---
!
      PLX = PLB(M,NBS)
      PGX = PGB(M,NBS)
      PNX = PLB(M,NBS)
      IF( IEQO.GT.0 ) PNX = PNB(M,NBS)
      TX = TB(M,NBS)
      PX = MAX( PLX,PGX,PNX ) + PATM
      CALL WATLQD( TX,PX,RHOLX )
      NS = IBCN(NBS)
      NE = IBCN(NBE)
      ISX = ID(NS)
      IEX = ID(NE)
      JSX = JD(NS)
      JEX = JD(NE)
      KSX = KD(NS)
      KEX = KD(NE)
      II = SIGN(1,IEX-ISX)
      JI = SIGN(1,JEX-JSX)
      KI = SIGN(1,KEX-KSX)
!
!---  Convert boundary pressure to node-centroid pressure  ---
!
      IF( IBCD(NBS).EQ.-3 ) THEN
        NPZ = NSZ(NS)
        PLX = PLX - 5.D-1*DZGF(NS)*RHOLX*GRVZ(NPZ)
      ELSEIF( IBCD(NBS).EQ.-2 ) THEN
        NPY = NSY(NS)
        PLX = PLX - 5.D-1*DYGF(NS)*RHOLX*GRVY(NPY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.-1 ) THEN
        NPX = NSX(NS)
        PLX = PLX - 5.D-1*DXGF(NS)*RHOLX*GRVX(NPX)
      ELSEIF( IBCD(NBS).EQ.1 ) THEN
        NQX = NSX(NS)+1
        PLX = PLX + 5.D-1*DXGF(NS)*RHOLX*GRVX(NQX)
      ELSEIF( IBCD(NBS).EQ.2 ) THEN
        NQY = NSY(NS)+IFLD
        PLX = PLX + 5.D-1*DYGF(NS)*RHOLX*GRVY(NQY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.3 ) THEN
        NQZ = NSZ(NS)+IJFLD
        PLX = PLX + 5.D-1*DZGF(NS)*RHOLX*GRVZ(NQZ)
      ENDIF
!
!---  Loop over nodes in the x direction
!
      I = ISX
      J = JSX
      K = KSX
      NX = ND(I,J,K)
      DO 100 I = ISX+II,IEX,II
        N = ND(I,J,K)
        NPX = NSX(N)
        GB = (XP(N)-XP(NX))*GRVX(NPX)
        PX = MAX( PLX,PGX,PNX ) + PATM
        CALL WATLQD( TX,PX,RHOLX )
        PLX = PLX - RHOLX*GB
        NX = N
  100 CONTINUE
!
!---  Loop over nodes in the y direction
!
      I = IEX
      J = JSX
      K = KSX
      NX = ND(I,J,K)
      DO 110 J = JSX+JI,JEX,JI
        N = ND(I,J,K)
        NPY = NSY(N)
        GB = (YP(N)*RP(I)-YP(NX)*RP(I))*GRVY(NPY)
        PX = MAX( PLX,PGX,PNX ) + PATM
        CALL WATLQD( TX,PX,RHOLX )
        PLX = PLX - RHOLX*GB
        NX = N
  110 CONTINUE
!
!---  Loop over nodes in the z direction
!
      I = IEX
      J = JEX
      K = KSX
      NX = ND(I,J,K)
      DO 120 K = KSX+KI,KEX,KI
        N = ND(I,J,K)
        NPZ = NSZ(N)
        GB = (ZP(N)-ZP(NX))*GRVZ(NPZ)
        PX = MAX( PLX,PGX,PNX ) + PATM
        CALL WATLQD( TX,PX,RHOLX )
        PLX = PLX - RHOLX*GB
        NX = N
  120 CONTINUE
!
!---  Convert node-centroid pressure to boundary pressure  ---
!
      NE = IBCN(NBE)
      IF( IBCD(NBE).EQ.-3 ) THEN
        NPZ = NSZ(NE)
        PLX = PLX + 5.D-1*DZGF(NE)*RHOLX*GRVZ(NPZ)
      ELSEIF( IBCD(NBE).EQ.-2 ) THEN
        NPY = NSY(NE)
        PLX = PLX + 5.D-1*DYGF(NE)*RHOLX*GRVY(NPY)*RP(ID(NE))
      ELSEIF( IBCD(NBE).EQ.-1 ) THEN
        NPX = NSX(NE)
        PLX = PLX + 5.D-1*DXGF(NE)*RHOLX*GRVX(NPX)
      ELSEIF( IBCD(NBE).EQ.1 ) THEN
        NQX = NSX(NE)+1
        PLX = PLX - 5.D-1*DXGF(NE)*RHOLX*GRVX(NQX)
      ELSEIF( IBCD(NBE).EQ.2 ) THEN
        NQY = NSY(NE)+IFLD
        PLX = PLX - 5.D-1*DYGF(NE)*RHOLX*GRVY(NQY)*RP(ID(NE))
      ELSEIF( IBCD(NBE).EQ.3 ) THEN
        NQZ = NSZ(NE)+IJFLD
        PLX = PLX - 5.D-1*DZGF(NE)*RHOLX*GRVZ(NQZ)
      ENDIF
      FNHGBL = PLX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FNHGBL group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FNHGBLX( XSBCX,YSBCX,ZSBCX,XPBCX,YPBCX,ZPBCX,TSX,
     &  TX,PLX )
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
!     Compute the aqueous pressure at a boundary surface, assuming
!     hydrostatic conditions, from a given pressure at a boundary
!     surface. Use twenty density calculations being points. 
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FNHGBLX'
!
!---  x,y,z deltas  ---
!
      XDX = 5.D-2*(XPBCX - XSBCX)
      YDX = 5.D-2*(YPBCX - YSBCX)
      ZDX = 5.D-2*(ZPBCX - ZSBCX)
      DO N = 1,20
        TDX = (TX-TSX)*5.D-2 + TSX
        CALL WATLQD( TDX,PLX,RHOLX )
        PLX = PLX - RHOLX*(XDX*GRAVX + YDX*GRAVY + ZDX*GRAVZ)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FNHGBLX group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATGSD( TX,PX,RHOX,INDX )
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
!     Calculate the water vapor component density
!     with the ideal gas law equation of state.
!
!     INDX = 1
!
!     Calculate the water vapor component density, as a function of
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
!     values:  0.01 C < T > 364.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     watgsd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 BX(31),SBX(5),L(3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE BX,SBX,L,SL1
      DATA BX /1.683599274D+1,2.856067796D+1,-5.438923329D+1,
     &4.330662834D-1,
     &-6.547711697D-1,8.565182058D-2,6.670375918D-2,1.388983801D+0,
     &8.390104328D-2,2.614670893D-2,-3.373439453D-2,4.520918904D-1,
     &1.069036614D-1,-5.975336707D-1,-8.847535804D-2,5.958051609D-1,
     &-5.159303373D-1,2.075021122D-1,1.190610271D-1,-9.867174132D-2,
     &1.683998803D-1,-5.809438001D-2,6.552390126D-3,5.710218649D-4,
     &1.936587558D+2,-1.388522425D+3,4.126607219D+3,-6.508211677D+3,
     &5.745984054D+3,-2.693088365D+3,5.235718623D+2/
      DATA SBX /7.633333333D-1,4.006073948D-1,8.636081627D-2,
     &-8.532322921D-1,3.460208861D-1/
      DATA L/1.574373327D+1,-3.417061978D+1,1.931380707D+1/
      DATA SL1/4.260321148D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATGSD'
      IF( INDX .EQ. 0 ) THEN
        RHOX = PX/((TX+TABS)*RCW)
      ELSEIF( INDX .EQ. 1 ) THEN
        IF( PX .LT. SMALL ) THEN
          RHOX = 0.D+0
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
        TR = (TX+TABS)/TCRW
        PR = PX/PCRW
        CX = EXP(SBX(1)*(1.D+0-TR))
        BL = L(1) + L(2)*TR + L(3)*TR*TR
!----------------------------------------------------------------------!
!     Compute the vapor density.
!----------------------------------------------------------------------!
        CHX = SL1*TR/PR
        CHX = CHX -(BX(7)*CX**13 + BX(8)*CX**3)
        CHX = CHX -2.D+0*PR*(BX(9)*CX**18 + BX(10)*CX**2 + BX(11)*CX)
        CHX = CHX -3.D+0*(PR**2)*(BX(12)*CX**18 + BX(13)*CX**10)
        CHX = CHX -4.D+0*(PR**3)*(BX(14)*CX**25 + BX(15)*CX**14)
        CHX = CHX -5.D+0*(PR**4)*
     &    (BX(16)*CX**32 + BX(17)*CX**28 + BX(18)*CX**24)
        CHX = CHX -4.D+0*(1.D+0/(PR**5))*(BX(19)*CX**12 +
     &    BX(20)*CX**11)/
     &    (((1.D+0/(PR**4))+(SBX(2)*CX**14))**2)
        CHX = CHX -5.D+0*(1.D+0/(PR**6))*(BX(21)*CX**24 +
     &    BX(22)*CX**18)/
     &    (((1.D+0/(PR**5))+(SBX(3)*CX**19))**2)
        CHX = CHX -6.D+0*(1.D+0/(PR**7))*(BX(23)*CX**24 +
     &    BX(24)*CX**14)/
     &    (((1.D+0/(PR**6))+(SBX(4)*CX**54 + SBX(5)*CX**27))**2)
        T1 = 1.1D+1*((PR/BL)**10)
        CHX = CHX + T1*BX(25)
        CHX = CHX + T1*BX(26)*CX
        CHX = CHX + T1*BX(27)*CX**2
        CHX = CHX + T1*BX(28)*CX**3
        CHX = CHX + T1*BX(29)*CX**4
        CHX = CHX + T1*BX(30)*CX**5
        CHX = CHX + T1*BX(31)*CX**6
        RHOX = 1.D+0/(CHX*VCRW*1.D-3/WTMW)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATGSD group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATLQD( TX,PX,RHO )
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
!     Calculate the subcooled or saturated density, as a function of
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
!     values:  0.01 C < T > 364.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     watlqd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 A(23),SA(12),B(10)
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
      DATA SA /8.438375405D-1,5.362162162D-4,1.720D+0,7.342278489D-2,
     &4.975858870D-2,6.537154300D-1,1.15D-6,1.5108D-5,1.4188D-1,
     &7.002753165D+0,2.995284926D-4,2.040D-1/
      DATA B/9.99667D+2,6.85021D-2,-7.0966D-3,2.76483D-5,-5.4108D-8,
     &5.20175D-7,-7.41396D-9,1.41879D-10,-8.82877D-13,1.92152D-15/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATLQD'
!
!---  Constant aqueous density  ---
!
      IF( ISLC(9).EQ.1 .OR. ISLC(53).EQ.1 ) THEN
        RHO = RHOLI
!
!---  Polynomial formulation  ---
!
      ELSEIF( IOM.NE.3 .AND. IOM.NE.13 ) THEN
        RHO = B(1) + PX*B(6)
        DO 10 M = 1,4
          RHO = RHO + B(M+1)*(TX**M) + PX*B(M+6)*(TX**M)
   10   CONTINUE
!
!---  ASME formulation  ---
!
      ELSE
        TR = MIN( (TX+TABS)/TCRW,1.D+0 )
        PR = PX/PCRW
        YC = 1.0D+0 - SA(1)*TR*TR-SA(2)/(TR**6)
        ZC = YC +
     &    SQRT(MAX( ZERO,(SA(3)*YC*YC-2.D+0*SA(4)*TR+2.D+0*SA(5)*PR)))
        CHI = A(12)*SA(5)/ZC**(2.9412D-1) +A(13) +A(14)*TR +A(15)*TR*TR
     &    +A(16)*(SA(6)-TR)**10 +A(17)/(SA(7)+TR**19)
     &    -(A(18)+2.D+0*A(19)*PR+3.D+0*A(20)*PR*PR)/(SA(8)+TR**11)
     &    -(A(21)*TR**18*(SA(9)+TR*TR)*(-3.D+0/((SA(10)+PR)**4)+SA(11)))
     &    +3.D+0*A(22)*(SA(12)-TR)*PR*PR +4.D+0*A(23)*PR*PR*PR/(TR**20)
        RHO = 1.D+0/(CHI*VCRW*1.D-3/WTMW)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATLQD group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ELC_DEN( RHOLX,CLX,DCFX )
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
!     Compute brine density from brine solution mass fraction and
!     liquid water density.
!
!     Leijnse, A.  1992.  Three-Dimensional Modeling of Coupled Flow
!     and Transport in Porous Media.  Ph.D. Dissertation, Department
!     of Civil Engineering and Geological Sciences, University of
!     Notre Dame, Notre Dame, Indiana.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, October 2000.
!     Last Modified by MD White, Battelle, PNNL, October 9, 2000.
!     watlqd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DCFX(4)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ELC_DEN'
!
!---  Electrolyte aqeuous density function by A. Leijnse  ---
!
      IF( IDF_ELC.EQ.1 ) THEN
        NC = 0
        RHOLSX = RHOLX
   10   CONTINUE
        NC = NC+1
        F = RHOLSX - RHOLX*EXP(DCFX(1)*CLX/RHOLSX)
        DF = 1.D+0 + RHOLX*DCFX(1)*CLX*EXP(DCFX(1)*CLX/RHOLSX)/
     &    (RHOLSX**2)
        DRHOX = -F/DF
        RHOLSX = RHOLSX + DRHOX
        IF( NC.GT.32 ) THEN
          INDX = 3
          CHMSG = 'Convergence Failure on Electrolyte Density'
          CALL WRMSGS( INDX )
        ENDIF
        IF( ABS(DRHOX/RHOLSX).GT.1.D-9 ) GOTO 10
        RHOLX = RHOLSX
!
!---  Fourth-order electrolyte aqeuous density function  ---
!
      ELSEIF( IDF_ELC.EQ.2 ) THEN
        CLXX = CLX/ELC_DUN
        RHOLX = RHOLX*(DCFX(1) + DCFX(2)*CLXX + DCFX(3)*(CLXX**2) +
     &    DCFX(4)*(CLXX**3))
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ELC_DEN group  ---
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION WLQDF( TX )
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
!     Water liquid density function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, 16 May 2001.
!     Last Modified by MD White, Battelle, PNNL, 16 May 2001.
!     watlqd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 B(10)
!
!----------------------Data Statements---------------------------------!
!
      SAVE B
      DATA B/9.99667D+2,6.85021D-2,-7.0966D-3,2.76483D-5,-5.4108D-8,
     &5.20175D-7,-7.41396D-9,1.41879D-10,-8.82877D-13,1.92152D-15/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WLQDF'
!
!---  Water liquid density function  ---
!
      WLQDF = B(1)+B(2)*TX+B(3)*(TX**2)+B(4)*(TX**3)+B(5)*(TX**4)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WLQDF group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION WLQDG( TX )
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
!     Water liquid density function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, 16 May 2001.
!     Last Modified by MD White, Battelle, PNNL, 16 May 2001.
!     watlqd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 B(10)
!
!----------------------Data Statements---------------------------------!
!
      SAVE B
      DATA B/9.99667D+2,6.85021D-2,-7.0966D-3,2.76483D-5,-5.4108D-8,
     &5.20175D-7,-7.41396D-9,1.41879D-10,-8.82877D-13,1.92152D-15/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WLQDG'
!
!---  Water liquid density function  ---
!
      WLQDG = B(6)+B(7)*TX+B(8)*(TX**2)+B(9)*(TX**3)+B(10)*(TX**4)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WLQDG group
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATLQV( TX,PX,PSWX,VIS )
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
!     Calculate the subcooled or saturated viscosity, as a function of
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
!     values:  0.01 C < T > 364.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     watlqv.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 K(4)
!
!----------------------Data Statements---------------------------------!
!
      SAVE K
      DATA K/-2.471D+1, 4.209D+3, 4.527D-02, -3.376D-5/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATLQV'
      TK = MIN( TX+TABS,TCRW )
      DPX = MAX( PX,PSWX ) - PSWX
      PHI = 1.D+0 + 1.0467D+0*(TK-3.05D+2)*DPX/1.D+11
      VIS = PHI*EXP( K(1) + K(2)/TK + K(3)*TK + K(4)*TK*TK )*1.D-3
      IF( ISLC(9).EQ.1 ) VIS = VISLI
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATLQV group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ELC_VIS( VISLX,CLX,VCFX )
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
!     Compute brine viscosity from brine solution mass fraction and
!     liquid water viscosity.
!
!     Leijnse, A.  1992.  Three-Dimensional Modeling of Coupled Flow
!     and Transport in Porous Media.  Ph.D. Dissertation, Department
!     of Civil Engineering and Geological Sciences, University of
!     Notre Dame, Notre Dame, Indiana.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, October 2000.
!     Last Modified by MD White, Battelle, PNNL, October 9, 2000.
!     watlqv.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VCFX(4)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ELC_VIS'
!
!---  Electrolyte aqeuous viscosity function by A. Leijnse  ---
!
      IF( IVF_ELC.EQ.1 ) THEN
        VISLX = VISLX*(VCFX(1) + VCFX(2)*CLX + VCFX(3)*(CLX**2)
     &    + VCFX(4)*(CLX**3))
!
!---  Fourth-order electrolyte aqeuous density function  ---
!
      ELSEIF( IDF_ELC.EQ.2 ) THEN
        CLXX = CLX/ELC_VUN
        VISLX = VISLX*(VCFX(1) + VCFX(2)*CLXX + VCFX(3)*(CLXX**2)
     &    + VCFX(4)*(CLXX**3))
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ELC_VIS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATSP( TX,PX )
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
!     Calculate the saturation pressure of water as a function of
!     temperature per the 1967 International Formulation Committee
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
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     watsp.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 K(9)
!
!----------------------Data Statements---------------------------------!
!
      SAVE K
      DATA K/-7.691234564D+0,-2.608023696D+1,-1.681706546D+2,
     &  6.423285504D+1,-1.189646225D+2,4.167117320D+0,2.097506760D+1,
     &  1.D+9,6.D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATSP'
      TR = MIN( (TX+TABS)/TCRW,1.D+0 )
      TR1 = 1.D+0-TR
      SVM = K(1)*(TR1**1) + K(2)*(TR1**2) + K(3)*(TR1**3)
     &  + K(4)*(TR1**4) + K(5)*(TR1**5)
      PR = EXP((SVM/(TR*(1.D+0+K(6)*TR1+K(7)*TR1*TR1)))
     &  + (TR1/(K(8)*TR1*TR1+K(9))))
      PX = PR*PCRW
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATSP group  ---
!
      RETURN
      END

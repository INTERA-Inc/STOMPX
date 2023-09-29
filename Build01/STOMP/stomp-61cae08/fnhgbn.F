!----------------------Function----------------------------------------!
!
      FUNCTION FNHGBN( NBS,NBE,M )
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
!     NAPL hydraulic gradient boundary condition function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, June 1994.
!     Last Modified by MD White, Battelle, PNL, September 10, 1997.
!     Last Modified by MD White, Battelle, PNNL, July 12, 2000.
!     fnhgbn.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/FNHGBN'
      PLX = PLB(M,NBS)
      PGX = PGB(M,NBS)
      PNX = PLB(M,NBS)
      TX = TB(M,NBS)
      PX = MAX( PLX+PATM,PGX+PATM,PNX+PATM )
      CALL VOCSP( TX,PSOX )
      CALL VOCLQD( TX,PX,PSOX,RHONX )
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
        PNX = PNX - 5.D-1*DZGF(NS)*RHONX*GRVZ(NPZ)
      ELSEIF( IBCD(NBS).EQ.-2 ) THEN
        NPY = NSY(NS)
        PNX = PNX - 5.D-1*DYGF(NS)*RHONX*GRVY(NPY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.-1 ) THEN
        NPX = NSX(NS)
        PNX = PNX - 5.D-1*DXGF(NS)*RHONX*GRVX(NPX)
      ELSEIF( IBCD(NBS).EQ.1 ) THEN
        NQX = NSX(NS)+1
        PNX = PNX + 5.D-1*DXGF(NS)*RHONX*GRVX(NQX)
      ELSEIF( IBCD(NBS).EQ.2 ) THEN
        NQY = NSY(NS)+IFLD
        PNX = PNX + 5.D-1*DYGF(NS)*RHONX*GRVY(NQY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.3 ) THEN
        NQZ = NSZ(NS)+IJFLD
        PNX = PNX + 5.D-1*DZGF(NS)*RHONX*GRVZ(NQZ)
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
        PX = MAX( PNX+PATM,PGX+PATM,PNX+PATM )
        CALL VOCLQD( TX,PX,PSOX,RHONX )
        PNX = PNX - RHONX*GB
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
        PX = MAX( PNX+PATM,PGX+PATM,PNX+PATM )
        CALL VOCLQD( TX,PX,PSOX,RHONX )
        PNX = PNX - RHONX*GB
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
        PX = MAX( PNX+PATM,PGX+PATM,PNX+PATM )
        CALL VOCLQD( TX,PX,PSOX,RHONX )
        PNX = PNX - RHONX*GB
        NX = N
  120 CONTINUE
!
!---  Convert node-centroid pressure to boundary pressure  ---
!
      NE = IBCN(NBE)
      IF( IBCD(NBE).EQ.-3 ) THEN
        NPZ = NSZ(NE)
        PNX = PNX + 5.D-1*DZGF(NE)*RHONX*GRVZ(NPZ)
      ELSEIF( IBCD(NBE).EQ.-2 ) THEN
        NPY = NSY(NE)
        PNX = PNX + 5.D-1*DYGF(NE)*RHONX*GRVY(NPY)*RP(ID(NE))
      ELSEIF( IBCD(NBE).EQ.-1 ) THEN
        NPX = NSX(NE)
        PNX = PNX + 5.D-1*DXGF(NE)*RHONX*GRVX(NPX)
      ELSEIF( IBCD(NBE).EQ.1 ) THEN
        NQX = NSX(NE)+1
        PNX = PNX - 5.D-1*DXGF(NE)*RHONX*GRVX(NQX)
      ELSEIF( IBCD(NBE).EQ.2 ) THEN
        NQY = NSY(NE)+IFLD
        PNX = PNX - 5.D-1*DYGF(NE)*RHONX*GRVY(NQY)*RP(ID(NE))
      ELSEIF( IBCD(NBE).EQ.3 ) THEN
        NQZ = NSZ(NE)+IJFLD
        PNX = PNX - 5.D-1*DZGF(NE)*RHONX*GRVZ(NQZ)
      ENDIF
      FNHGBN = PNX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FNHGBN group
!
      RETURN
      END


!----------------------Function----------------------------------------!
!
      FUNCTION HGBG_EOR( NBS,NBE,M )
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
!     Gas hydraulic gradient boundary condition function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 May 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
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
      SUB_LOG(ISUB_LOG) = '/HGBG_EOR'
      PGX = PGB(M,NBS)
      PVBX = PVWB(M,NBS)
      PVAX = PVAB(M,NBS)
      PVOX = PVOB(M,NBS)
      TX = TB(M,NBS)
      ISRX = 2
      CALL DENS_W( TX,PVBX,RHOLWX,RHOGWX,ISRX )
!      CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
!      CALL DENS_O( TX,PVOX,RHOGOX )
      RHOGX = RHOGWX+RHOGAX+RHOGOX
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
        PGX = PGX + 5.D-1*DZGF(NS)*RHOGX*GRVZ(NPZ)
      ELSEIF( IBCD(NBS).EQ.-2 ) THEN
        NPY = NSY(NS)
        PGX = PGX + 5.D-1*DYGF(NS)*RHOGX*GRVY(NPY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.-1 ) THEN
        NPX = NSX(NS)
        PGX = PGX + 5.D-1*DXGF(NS)*RHOGX*GRVX(NPX)
      ELSEIF( IBCD(NBS).EQ.1 ) THEN
        NQX = NSX(NS)+1
        PGX = PGX - 5.D-1*DXGF(NS)*RHOGX*GRVX(NQX)
      ELSEIF( IBCD(NBS).EQ.2 ) THEN
        NQY = NSY(NS)+IFLD
        PGX = PGX - 5.D-1*DYGF(NS)*RHOGX*GRVY(NQY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.3 ) THEN
        NQZ = NSZ(NS)+IJFLD
        PGX = PGX - 5.D-1*DZGF(NS)*RHOGX*GRVZ(NQZ)
      ENDIF
!
!---  Loop over nodes in the x direction
!
      I = ISX
      J = JSX
      K = KSX
      NX = ND(I,J,K)
      DO 100 I = ISX,IEX,II
        N = ND(I,J,K)
        NPX = NSX(N)
        GB = (XP(N)-XP(NX))*GRVX(NPX)
        ISRX = 2
        CALL DENS_W( TX,PVBX,RHOLWX,RHOGWX,ISRX )
!        CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
!        CALL DENS_O( TX,PVOX,RHOGOX )
        RHOGX = RHOGWX+RHOGAX+RHOGOX
        PGX = PGX - RHOGX*GB
        NX = N
  100 CONTINUE
!
!---  Loop over nodes in the y direction
!
      I = IEX
      J = JSX
      K = KSX
      NX = ND(I,J,K)
      DO 110 J = JSX,JEX,JI
        N = ND(I,J,K)
        NPY = NSY(N)
        GB = (YP(N)*RP(I)-YP(NX)*RP(I))*GRVY(NPY)
        ISRX = 2
        CALL DENS_W( TX,PVBX,RHOLWX,RHOGWX,ISRX )
!        CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
!        CALL DENS_O( TX,PVOX,RHOGOX )
        RHOGX = RHOGWX+RHOGAX+RHOGOX
        PGX = PGX - RHOGX*GB
        NX = N
  110 CONTINUE
!
!---  Loop over nodes in the z direction
!
      I = IEX
      J = JEX
      K = KSX
      NX = ND(I,J,K)
      DO 120 K = KSX,KEX,KI
        N = ND(I,J,K)
        NPZ = NSZ(N)
        GB = (ZP(N)-ZP(NX))*GRVZ(NPZ)
        ISRX = 2
        CALL DENS_W( TX,PVBX,RHOLWX,RHOGWX,ISRX )
!        CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
!        CALL DENS_O( TX,PVOX,RHOGOX )
        RHOGX = RHOGWX+RHOGAX+RHOGOX
        PGX = PGX - RHOGX*GB
        NX = N
  120 CONTINUE
!
!---  Convert node-centroid pressure to boundary pressure  ---
!
      NE = IBCN(NBE)
      IF( IBCD(NBE).EQ.-3 ) THEN
        NPZ = NSZ(NE)
        PGX = PGX + 5.D-1*DZGF(NE)*RHOGX*GRVZ(NPZ)
      ELSEIF( IBCD(NBE).EQ.-2 ) THEN
        NPY = NSY(NE)
        PGX = PGX + 5.D-1*DYGF(NE)*RHOGX*GRVY(NPY)*RP(ID(NE))
      ELSEIF( IBCD(NBE).EQ.-1 ) THEN
        NPX = NSX(NE)
        PGX = PGX + 5.D-1*DXGF(NE)*RHOGX*GRVX(NPX)
      ELSEIF( IBCD(NBE).EQ.1 ) THEN
        NQX = NSX(NE)+1
        PGX = PGX - 5.D-1*DXGF(NE)*RHOGX*GRVX(NQX)
      ELSEIF( IBCD(NBE).EQ.2 ) THEN
        NQY = NSY(NE)+IFLD
        PGX = PGX - 5.D-1*DYGF(NE)*RHOGX*GRVY(NQY)*RP(ID(NE))
      ELSEIF( IBCD(NBE).EQ.3 ) THEN
        NQZ = NSZ(NE)+IJFLD
        PGX = PGX - 5.D-1*DZGF(NE)*RHOGX*GRVZ(NQZ)
      ENDIF
      HGBG_EOR = PGX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HGBG_EOR group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION HGBL_EOR( NBS,NBE,M )
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
!     Written by MD White, PNNL, 13 May 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
      USE BCVS
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/HGBL_EOR'
      PLX = PLB(M,NBS)
      PGX = PGB(M,NBS)
      PGX = PGB(M,NBS)
      PVBX = PVWB(M,NBS)
      PVAX = PVAB(M,NBS)
      PVOX = PVOB(M,NBS)
      XLSX = XLSB(M,NBS)
      XLAX = XLAB(M,NBS)
      XLOX = XLOB(M,NBS)
      PX = MAX( PLX+PATM,PGX+PATM )
      CALL DENS_B( TX,PX,XLSX,RHOBX )
      CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
!
!---  Micro-emulsion  ---
!
      IF( ZLAB(M,NBS).GT.EPSL ) THEN
        CALL SP_A( TX,PSAX )
!
!---    Liquid micro-emulsion  ---
!
        IF( PX.GT.PSAX ) THEN
          ISRX = 1
!          CALL DENS_A( TX,PX,RHONX,ISRX )
          RHOMEX = RHOLX*(1.D+0-ZLAB(M,NBS)) + RHONX*ZLAB(M,NBS)
!
!---    Vapor micro-emulsion  ---
!
        ELSE
          ISRX = 2
          CALL DENS_W( TX,PVBX,RHOLWX,RHOGWX,ISRX )
!          CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
!          CALL DENS_O( TX,PVOX,RHOGOX )
          RHOGX = RHOGWX+RHOGAX+RHOGOX
          RHOMEX = RHOLX*(1.D+0-ZLAB(M,NBS)) + RHOGX*ZLAB(M,NBS)
        ENDIF
      ENDIF
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
        PLX = PLX + 5.D-1*DZGF(NS)*RHOMEX*GRVZ(NPZ)
      ELSEIF( IBCD(NBS).EQ.-2 ) THEN
        NPY = NSY(NS)
        PLX = PLX + 5.D-1*DYGF(NS)*RHOMEX*GRVY(NPY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.-1 ) THEN
        NPX = NSX(NS)
        PLX = PLX + 5.D-1*DXGF(NS)*RHOMEX*GRVX(NPX)
      ELSEIF( IBCD(NBS).EQ.1 ) THEN
        NQX = NSX(NS)+1
        PLX = PLX - 5.D-1*DXGF(NS)*RHOMEX*GRVX(NQX)
      ELSEIF( IBCD(NBS).EQ.2 ) THEN
        NQY = NSY(NS)+IFLD
        PLX = PLX - 5.D-1*DYGF(NS)*RHOMEX*GRVY(NQY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.3 ) THEN
        NQZ = NSZ(NS)+IJFLD
        PLX = PLX - 5.D-1*DZGF(NS)*RHOMEX*GRVZ(NQZ)
      ENDIF
!
!---  Loop over nodes in the x direction
!
      I = ISX
      J = JSX
      K = KSX
      NX = ND(I,J,K)
      DO 100 I = ISX,IEX,II
        N = ND(I,J,K)
        NPX = NSX(N)
        GB = (XP(N)-XP(NX))*GRVX(NPX)
        PX = MAX( PLX+PATM,PGX+PATM )
        CALL DENS_B( TX,PX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
!
!---    Micro-emulsion  ---
!
        IF( ZLAB(M,NBS).GT.EPSL ) THEN
          CALL SP_A( TX,PSAX )
!
!---      Liquid micro-emulsion  ---
!
          IF( PX.GT.PSAX ) THEN
            ISRX = 1
!            CALL DENS_A( TX,PX,RHONX,ISRX )
            RHOMEX = RHOLX*(1.D+0-ZLAB(M,NBS)) + RHONX*ZLAB(M,NBS)
!
!---      Vapor micro-emulsion  ---
!
          ELSE
            ISRX = 2
            CALL DENS_W( TX,PVBX,RHOLWX,RHOGWX,ISRX )
!            CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
!            CALL DENS_O( TX,PVOX,RHOGOX )
            RHOGX = RHOGWX+RHOGAX+RHOGOX
            RHOMEX = RHOLX*(1.D+0-ZLAB(M,NBS)) + RHOGX*ZLAB(M,NBS)
          ENDIF
        ENDIF
        PLX = PLX - RHOMEX*GB
        NX = N
  100 CONTINUE
!
!---  Loop over nodes in the y direction
!
      I = IEX
      J = JSX
      K = KSX
      NX = ND(I,J,K)
      DO 110 J = JSX,JEX,JI
        N = ND(I,J,K)
        NPY = NSY(N)
        GB = (YP(N)*RP(I)-YP(NX)*RP(I))*GRVY(NPY)
        PX = MAX( PLX+PATM,PGX+PATM )
        CALL DENS_B( TX,PX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
!
!---    Micro-emulsion  ---
!
        IF( ZLAB(M,NBS).GT.EPSL ) THEN
          CALL SP_A( TX,PSAX )
!
!---      Liquid micro-emulsion  ---
!
          IF( PX.GT.PSAX ) THEN
            ISRX = 1
!            CALL DENS_A( TX,PX,RHONX,ISRX )
            RHOMEX = RHOLX*(1.D+0-ZLAB(M,NBS)) + RHONX*ZLAB(M,NBS)
!
!---      Vapor micro-emulsion  ---
!
          ELSE
            ISRX = 2
            CALL DENS_W( TX,PVBX,RHOLWX,RHOGWX,ISRX )
!            CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
!            CALL DENS_O( TX,PVOX,RHOGOX )
            RHOGX = RHOGWX+RHOGAX+RHOGOX
            RHOMEX = RHOLX*(1.D+0-ZLAB(M,NBS)) + RHOGX*ZLAB(M,NBS)
          ENDIF
        ENDIF
        PLX = PLX - RHOMEX*GB
        NX = N
  110 CONTINUE
!
!---  Loop over nodes in the z direction
!
      I = IEX
      J = JEX
      K = KSX
      NX = ND(I,J,K)
      DO 120 K = KSX,KEX,KI
        N = ND(I,J,K)
        NPZ = NSZ(N)
        GB = (ZP(N)-ZP(NX))*GRVZ(NPZ)
        PX = MAX( PLX+PATM,PGX+PATM )
        CALL DENS_B( TX,PX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
!
!---    Micro-emulsion  ---
!
        IF( ZLAB(M,NBS).GT.EPSL ) THEN
          CALL SP_A( TX,PSAX )
!
!---      Liquid micro-emulsion  ---
!
          IF( PX.GT.PSAX ) THEN
            ISRX = 1
!            CALL DENS_A( TX,PX,RHONX,ISRX )
            RHOMEX = RHOLX*(1.D+0-ZLAB(M,NBS)) + RHONX*ZLAB(M,NBS)
!
!---      Vapor micro-emulsion  ---
!
          ELSE
            ISRX = 2
            CALL DENS_W( TX,PVBX,RHOLWX,RHOGWX,ISRX )
!            CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
!            CALL DENS_O( TX,PVOX,RHOGOX )
            RHOGX = RHOGWX+RHOGAX+RHOGOX
            RHOMEX = RHOLX*(1.D+0-ZLAB(M,NBS)) + RHOGX*ZLAB(M,NBS)
          ENDIF
        ENDIF
        PLX = PLX - RHOMEX*GB
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
      HGBL_EOR = PLX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HGBL_EOR group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION HGBN_EOR( NBS,NBE,M )
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
!     Written by MD White, PNNL, 13 May 2005.
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
      SUB_LOG(ISUB_LOG) = '/HGBN_EOR'
      PLX = PLB(M,NBS)
      PGX = PGB(M,NBS)
      PNX = PLB(M,NBS)
      TX = TB(M,NBS)
      PX = MAX( PLX+PATM,PGX+PATM,PNX+PATM )
      ISRX = 1
!      CALL DENS_A( TX,PX,RHONX,ISRX )
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
        PNX = PNX + 5.D-1*DZGF(NS)*RHONX*GRVZ(NPZ)
      ELSEIF( IBCD(NBS).EQ.-2 ) THEN
        NPY = NSY(NS)
        PNX = PNX + 5.D-1*DYGF(NS)*RHONX*GRVY(NPY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.-1 ) THEN
        NPX = NSX(NS)
        PNX = PNX + 5.D-1*DXGF(NS)*RHONX*GRVX(NPX)
      ELSEIF( IBCD(NBS).EQ.1 ) THEN
        NQX = NSX(NS)+1
        PNX = PNX - 5.D-1*DXGF(NS)*RHONX*GRVX(NQX)
      ELSEIF( IBCD(NBS).EQ.2 ) THEN
        NQY = NSY(NS)+IFLD
        PNX = PNX - 5.D-1*DYGF(NS)*RHONX*GRVY(NQY)*RP(ID(NS))
      ELSEIF( IBCD(NBS).EQ.3 ) THEN
        NQZ = NSZ(NS)+IJFLD
        PNX = PNX - 5.D-1*DZGF(NS)*RHONX*GRVZ(NQZ)
      ENDIF
!
!---  Loop over nodes in the x direction
!
      I = ISX
      J = JSX
      K = KSX
      NX = ND(I,J,K)
      DO 100 I = ISX,IEX,II
        N = ND(I,J,K)
        NPX = NSX(N)
        GB = (XP(N)-XP(NX))*GRVX(NPX)
        PX = MAX( PNX+PATM,PGX+PATM,PNX+PATM )
        ISRX = 1
!        CALL DENS_A( TX,PX,RHONX,ISRX )
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
      DO 110 J = JSX,JEX,JI
        N = ND(I,J,K)
        NPY = NSY(N)
        GB = (YP(N)*RP(I)-YP(NX)*RP(I))*GRVY(NPY)
        PX = MAX( PNX+PATM,PGX+PATM,PNX+PATM )
        ISRX = 1
!        CALL DENS_A( TX,PX,RHONX,ISRX )
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
      DO 120 K = KSX,KEX,KI
        N = ND(I,J,K)
        NPZ = NSZ(N)
        GB = (ZP(N)-ZP(NX))*GRVZ(NPZ)
        PX = MAX( PNX+PATM,PGX+PATM,PNX+PATM )
        ISRX = 1
 !       CALL DENS_A( TX,PX,RHONX,ISRX )
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
      HGBN_EOR = PNX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HGBN_EOR group
!
      RETURN
      END

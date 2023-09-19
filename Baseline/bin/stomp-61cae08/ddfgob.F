!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFGOB( N,NB )
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
!     Compute advective, diffusive, dispersive gas oil
!     fluxes on bottom boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddfgob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
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
      SUB_LOG(ISUB_LOG) = '/DDFGOB'
      IEQDX = IEQD
      IF( IOM.EQ.4 .OR. IOM.EQ.5 .OR. IOM.EQ.24 ) IEQDX = IEQDO
      K = KD(N)
      IZN = IZ(N)
      NPZ = NSZ(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Diffusion coefficients  ---
!
        DFSP = TORG(MP,N)*(SG(MP,N)-SGT(MP,N))*PORD(MP,N)*DFGO(MP,N)
     &    *RHOMG(MP,N)
        DFSB = TORGB(MP,NB)*SGB(MP,NB)*PORDB(MP,NB)*DFGOB(MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFSB = DIFMN( DFSB,DFSP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVBB( PORD(MP,N),PORDB(MP,NB),SG(MP,N),SGB(MP,NB),
     &      UG,VG,WG,UGX,VGX,WGX,N,M )
          UGX = UGX*UGX
          VGX = VGX*VGX
          WGX = WGX*WGX
          ZGBX = SQRT(UGX+VGX+WGX)
          DPSB = (DPLD(IZN)*WGX + DPTD(IZN)*(UGX+VGX))/(ZGBX+SMALL)
        ELSE
          DPSB = 0.D+0
        ENDIF
!
!---   Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDGB = (DFSB+DPSB)/(5.D-1*DZGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
            VMCB = SGB(MP,NB)*PORDB(MP,NB)
            INDX = 13
            VMCB = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WG(1,NPZ),INDX)
            CRGB = ABS(WG(M,NPZ))*DT/(DZGF(N)*VMCB+SMALL)
!
!---        On the boundary surface  --
!
            WDGO(M,NPZ) = DDGB*(XMGOB(MP,NB)-XMGO(MP,N))
            IF( WG(1,NPZ).GE.0.D+0 ) THEN
              WGO(M,NPZ) = XGOB(MP,NB)*RHOGB(MP,NB)*WG(M,NPZ)
            ELSEIF( WG(1,NPZ).LT.0.D+0 .AND. K.LT.KFLD ) THEN
              NTN = N+IJFLD
              R = ((XMGO(1,N)-XMGO(1,NTN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NTN)))
              THETA = FLIMIT( R,CRGB,ISLC(8) )
              DZF = 5.D-1
              WGO(M,NPZ) = XGOB(MP,NB)*RHOGB(MP,NB)*WG(M,NPZ)*THETA*DZF
     &          + XGO(MP,N)*RHOG(MP,N)*WG(M,NPZ)*(1.D+0-THETA*DZF)
            ELSEIF( WG(1,NPZ).LT.0.D+0 .AND. K.EQ.KFLD ) THEN
              WGO(M,NPZ) = XGO(MP,N)*RHOG(MP,N)*WG(M,NPZ)
            ENDIF
            WGO(M,NPZ) = WGO(M,NPZ) + WDGO(M,NPZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQZ = NPZ+IJFLD
            NTN = N+IJFLD
            MN = MNEG(M)
            IF( WG(1,NQZ).GE.0.D+0 ) THEN
              VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
              VMCT = (SG(MP,NTN)-SGT(MP,NTN))*PORD(MP,NTN)
              INDX = 13
              VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NTN),WG(1,NQZ),INDX)
              CRGT = ABS( WG(M,NQZ) )*DT/(DZGF(NTN)*VMCT+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NTN)-XMGO(1,N)+SMALL))
     &          *((DZGF(NTN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRGT,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NTN)+DZGF(N))
              DZF = 5.D-1
              WGO(M,NQZ) = XGO(MN,NTN)*RHOG(MN,NTN)*WG(M,NQZ)*THETA*DZF
     &          + XGO(MP,N)*RHOG(MP,N)*WG(M,NQZ)*(1.D+0-THETA*DZF)
            ENDIF
            WGO(M,NQZ) = WGO(M,NQZ) + WDGO(M,NQZ)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AG = MAX( WG(M,NPZ),ZERO ) +
     &       DDGB*MAX((ONE-(TENTH*ABS(WG(M,NPZ))/(DDGB+SMALL)))**5,ZERO)
            AGP = MAX( -WG(M,NPZ),ZERO ) +
     &       DDGB*MAX((ONE-(TENTH*ABS(WG(M,NPZ))/(DDGB+SMALL)))**5,ZERO)
            WGO(M,NPZ) = WTMO*(XMGOB(MP,NB)*AG-XMGO(MP,N)*AGP)
            WDGO(M,NPZ) = DDGB*(XMGOB(MP,NB)-XMGO(MP,N))
!
!---      Upwind scheme  --
!
          ELSE
            WDGO(M,NPZ) = DDGB*(XMGOB(MP,NB)-XMGO(MP,N))
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
            INDX = -4
            FGO = DIFMN( FGOB,FGOP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
            WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
          ENDIF
!
!---    Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
            VMCB = SGB(MP,NB)*PORDB(MP,NB)
            INDX = 13
            VMCB = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WG(1,NPZ),INDX)
            CRGB = ABS(WG(M,NPZ))*DT/(DZGF(N)*VMCB+SMALL)
!
!---        On the boundary surface  --
!
            WDGO(M,NPZ) = 0.D+0
            IF( WG(1,NPZ).LT.0.D+0 .AND. K.LT.KFLD ) THEN
              NTN = N+IJFLD
              R = ((XMGO(1,N)-XMGO(1,NTN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NTN)))
              THETA = FLIMIT( R,CRGB,ISLC(8) )
              DZF = 5.D-1
              WGO(M,NPZ) = XGOB(MP,NB)*RHOGB(MP,NB)*WG(M,NPZ)*THETA*DZF
     &          + XGO(MP,N)*RHOG(MP,N)*WG(M,NPZ)*(1.D+0-THETA*DZF)
            ELSEIF( WG(1,NPZ).LT.0.D+0 .AND. K.EQ.KFLD ) THEN
              WGO(M,NPZ) = XGO(MP,N)*RHOG(MP,N)*WG(M,NPZ)
            ENDIF
            WGO(M,NPZ) = WGO(M,NPZ) + WDGO(M,NPZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQZ = NPZ+IJFLD
            NTN = N+IJFLD
            MN = MNEG(M)
            IF( WG(1,NQZ).GE.0.D+0 ) THEN
              VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
              VMCT = (SG(MP,NTN)-SGT(MP,NTN))*PORD(MP,NTN)
              INDX = 13
              VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NTN),WG(1,NQZ),INDX)
              CRGT = ABS( WG(M,NQZ) )*DT/(DZGF(NTN)*VMCT+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NTN)-XMGO(1,N)+SMALL))
     &          *((DZGF(NTN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRGT,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NTN)+DZGF(N))
              DZF = 5.D-1
              WGO(M,NQZ) = XGO(MN,NTN)*RHOG(MN,NTN)*WG(M,NQZ)*THETA*DZF
     &          + XGO(MP,N)*RHOG(MP,N)*WG(M,NQZ)*(1.D+0-THETA*DZF)
            ENDIF
            WGO(M,NQZ) = WGO(M,NQZ) + WDGO(M,NQZ)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            WDGO(M,NPZ) = 0.D+0
            AGP = MAX( -WG(M,NPZ),ZERO )
            WGO(M,NPZ) = -WTMO*XMGO(MP,N)*AGP
!
!---      Upwind scheme  --
!
          ELSE
            WDGO(M,NPZ) = 0.D+0
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOB = 0.D+0
            INDX = -4
            FGO = DIFMN( FGOB,FGOP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
            WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
          ENDIF
!
!---    Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
            VMCB = SGB(MP,NB)*PORDB(MP,NB)
            INDX = 13
            VMCB = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WG(1,NPZ),INDX)
            CRGB = ABS(WG(M,NPZ))*DT/(DZGF(N)*VMCB+SMALL)
!
!---        On the boundary surface  --
!
            WDGO(M,NPZ) = 0.D+0
            IF( WG(1,NPZ).GE.0.D+0 ) THEN
              WGO(M,NPZ) = XGOB(MP,NB)*RHOGB(MP,NB)*WG(M,NPZ)
            ENDIF
            WGO(M,NPZ) = WGO(M,NPZ) + WDGO(M,NPZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQZ = NPZ+IJFLD
            NTN = N+IJFLD
            MN = MNEG(M)
            IF( WG(1,NQZ).GE.0.D+0 ) THEN
              VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
              VMCT = (SG(MP,NTN)-SGT(MP,NTN))*PORD(MP,NTN)
              INDX = 13
              VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NTN),WG(1,NQZ),INDX)
              CRGT = ABS( WG(M,NQZ) )*DT/(DZGF(NTN)*VMCT+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NTN)-XMGO(1,N)+SMALL))
     &          *((DZGF(NTN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRGT,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NTN)+DZGF(N))
              DZF = 5.D-1
              WGO(M,NQZ) = XGO(MN,NTN)*RHOG(MN,NTN)*WG(M,NQZ)*THETA*DZF
     &          + XGO(MP,N)*RHOG(MP,N)*WG(M,NQZ)*(1.D+0-THETA*DZF)
            ENDIF
            WGO(M,NQZ) = WGO(M,NQZ) + WDGO(M,NQZ)*WTMO
!
!---      Patankar scheme  ---
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            WDGO(M,NPZ) = 0.D+0
            AG = MAX( WG(M,NPZ),ZERO )
            WGO(M,NPZ) = WTMO*XMGOB(MP,NB)*AG
!
!---      Upwind scheme  ---
!
          ELSE
            WDGO(M,NPZ) = 0.D+0
            FGOP = 0.D+0
            FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
            INDX = -4
            FGO = DIFMN( FGOB,FGOP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
            WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFGOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFGOS( N,NB )
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
!     Compute advective, diffusive, dispersive gas oil
!     fluxes on south boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddfgob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
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
      SUB_LOG(ISUB_LOG) = '/DDFGOS'
      IEQDX = IEQD
      IF( IOM.EQ.4 .OR. IOM.EQ.5 .OR. IOM.EQ.24 ) IEQDX = IEQDO
      J = JD(N)
      I = ID(N)
      IZN = IZ(N)
      NPY = NSY(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        DFSP = TORG(MP,N)*(SG(MP,N)-SGT(MP,N))*PORD(MP,N)*DFGO(MP,N)
     &    *RHOMG(MP,N)
        DFSS = TORGB(MP,NB)*SGB(MP,NB)*PORDB(MP,NB)*DFGOB(MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFSS = DIFMN(DFSS,DFSP,DYGF(N),DYGF(N),VG(1,NPY),INDX)
!
!---  Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVSB( PORD(MP,N),PORDB(MP,NB),SG(MP,N),SGB(MP,NB),
     &      UG,VG,WG,UGX,VGX,WGX,N,M )
          UGX = UGX*UGX
          VGX = VGX*VGX
          WGX = WGX*WGX
          ZGS = SQRT(UGX+VGX+WGX)
          DPSS = (DPLD(IZN)*VGX + DPTD(IZN)*(UGX+WGX))/(ZGS+SMALL)
        ELSE
          DPSS = 0.D+0
        ENDIF
!
!---    Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDGS = (DFSS+DPSS)/RP(I)/(5.D-1*DYGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
            VMCS = SGB(MP,NB)*PORDB(MP,NB)
            INDX = 13
            VMCS = DIFMN(VMCS,VMCP,DYGF(N),DYGF(N),VG(1,NPY),INDX)
            CRGS = ABS(VG(M,NPY))*DT/(DYGP(NPY)*VMCS+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDGO(M,NPY) = DDGS*(XMGOB(MP,NB)-XMGO(MP,N))
            IF( VG(1,NPY).GE.0.D+0 ) THEN
              VGO(M,NPY) = XGOB(MP,NB)*RHOGB(MP,NB)*VG(M,NPY)
            ELSEIF( VG(1,NPY).LT.0.D+0 .AND. J.LT.JFLD ) THEN
              NNN = N+IFLD
              R = ((XMGO(1,N)-XMGO(1,NNN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NNN)))
              THETA = FLIMIT( R,CRGS,ISLC(8) )
              DYF = 5.D-1
              VGO(M,NPY) = XGOB(MP,NB)*RHOGB(MP,NB)*VG(M,NPY)*THETA*DYF
     &          + XGO(MP,N)*RHOG(MP,N)*VG(M,NPY)*(1.D+0-THETA*DYF)
            ELSEIF( VG(1,NPY).LT.0.D+0 .AND. J.EQ.JFLD ) THEN
              VGO(M,NPY) = XGO(MP,N)*RHOG(MP,N)*VG(M,NPY)
            ENDIF
            VGO(M,NPY) = VGO(M,NPY) + VDGO(M,NPY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQY = NPY+IFLD
            NNN = N+IFLD
            MN = MNEG(M)
            IF( VG(1,NQY).GE.0.D+0 ) THEN
              VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
              VMCN = (SG(MP,NNN)-SGT(MP,NNN))*PORD(MP,NNN)
              INDX = 13
              VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(NNN),VG(1,NQY),INDX)
              CRGN = ABS( VG(M,NQY) )*DT/(DYGF(NNN)*VMCN+SMALL)/RP(I)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NNN)-XMGO(1,N)+SMALL))
     &          *((DYGF(NNN)+DYGF(N))/(DYGF(N)))
              THETA = FLIMIT( R,CRGN,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NNN)+DYGF(N))
              DYF = 5.D-1
              VGO(M,NQY) = XGO(MN,NNN)*RHOG(MN,NNN)*VG(M,NQY)*THETA*DYF
     &          + XGO(MP,N)*RHOG(MP,N)*VG(M,NQY)*(1.D+0-THETA*DYF)
            ENDIF
            VGO(M,NQY) = VGO(M,NQY) + VDGO(M,NQY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AG = MAX( VG(M,NPY),ZERO ) +
     &       DDGS*MAX((ONE-(TENTH*ABS(VG(M,NPY))/(DDGS+SMALL)))**5,ZERO)
            AGP = MAX( -VG(M,NPY),ZERO ) +
     &       DDGS*MAX((ONE-(TENTH*ABS(VG(M,NPY))/(DDGS+SMALL)))**5,ZERO)
            VGO(M,NPY) = WTMO*(XMGOB(MP,NB)*AG-XMGO(MP,N)*AGP)
            VDGO(M,NPY) = DDGS*(XMGOB(MP,NB)*-XMGO(MP,N))
!
!---      Upwind scheme  --
!
          ELSE
            VDGO(M,NPY) = DDGS*(XMGOB(MP,NB)-XMGO(MP,N))
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOS = XGOB(MP,NB)*RHOGB(MP,NB)
            INDX = -4
            FGO = DIFMN( FGOS,FGOP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
            VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
          ENDIF
!
!---   Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
            VMCS = SGB(MP,NB)*PORDB(MP,NB)
            INDX = 13
            VMCS = DIFMN(VMCS,VMCP,DYGF(N),DYGF(N),VG(1,NPY),INDX)
            CRGS = ABS(VG(M,NPY))*DT/(DYGP(NPY)*VMCS+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDGO(M,NPY) = 0.D+0
            IF( VG(1,NPY).LT.0.D+0 .AND. J.LT.JFLD ) THEN
              NNN = N+IFLD
              R = ((XMGO(1,N)-XMGO(1,NNN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NNN)))
              THETA = FLIMIT( R,CRGS,ISLC(8) )
              DYF = 5.D-1
              VGO(M,NPY) = XGOB(MP,NB)*RHOGB(MP,NB)*VG(M,NPY)*THETA*DYF
     &          + XGO(MP,N)*RHOG(MP,N)*VG(M,NPY)*(1.D+0-THETA*DYF)
            ELSEIF( VG(1,NPY).LT.0.D+0 .AND. J.EQ.JFLD ) THEN
              VGO(M,NPY) = XGO(MP,N)*RHOG(MP,N)*VG(M,NPY)
            ENDIF
            VGO(M,NPY) = VGO(M,NPY) + VDGO(M,NPY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQY = NPY+IFLD
            NNN = N+IFLD
            MN = MNEG(M)
            IF( VG(1,NQY).GE.0.D+0 ) THEN
              VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
              VMCN = (SG(MP,NNN)-SGT(MP,NNN))*PORD(MP,NNN)
              INDX = 13
              VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(NNN),VG(1,NQY),INDX)
              CRGN = ABS( VG(M,NQY) )*DT/(DYGF(NNN)*VMCN+SMALL)/RP(I)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NNN)-XMGO(1,N)+SMALL))
     &          *((DYGF(NNN)+DYGF(N))/(DYGF(N)))
              THETA = FLIMIT( R,CRGN,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NNN)+DYGF(N))
              DYF = 5.D-1
              VGO(M,NQY) = XGO(MN,NNN)*RHOG(MN,NNN)*VG(M,NQY)*THETA*DYF
     &          + XGO(MP,N)*RHOG(MP,N)*VG(M,NQY)*(1.D+0-THETA*DYF)
            ENDIF
            VGO(M,NQY) = VGO(M,NQY) + VDGO(M,NQY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            VDGO(M,NPY) = 0.D+0
            AGP = MAX( -VG(M,NPY),ZERO )
            VGO(M,NPY) = -WTMO*XMGO(MP,N)*AGP
!
!---      Upwind scheme  --
!
          ELSE
            VDGO(M,NPY) = 0.D+0
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOS = 0.D+0
            INDX = -4
            FGO = DIFMN( FGOS,FGOP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
            VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
          ENDIF
!
!---   Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
            VMCS = SGB(MP,NB)*PORDB(MP,NB)
            INDX = 13
            VMCS = DIFMN(VMCS,VMCP,DYGF(N),DYGF(N),VG(1,NPY),INDX)
            CRGS = ABS(VG(M,NPY))*DT/(DYGP(NPY)*VMCS+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDGO(M,NPY) = 0.D+0
            IF( VG(1,NPY).GE.0.D+0 ) THEN
              VGO(M,NPY) = XGOB(MP,NB)*RHOGB(MP,NB)*VG(M,NPY)
            ENDIF
            VGO(M,NPY) = VGO(M,NPY) + VDGO(M,NPY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQY = NPY+IFLD
            NNN = N+IFLD
            MN = MNEG(M)
            IF( VG(1,NQY).GE.0.D+0 ) THEN
              VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
              VMCN = (SG(MP,NNN)-SGT(MP,NNN))*PORD(MP,NNN)
              INDX = 13
              VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(NNN),VG(1,NQY),INDX)
              CRGN = ABS( VG(M,NQY) )*DT/(DYGF(NNN)*VMCN+SMALL)/RP(I)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NNN)-XMGO(1,N)+SMALL))
     &          *((DYGF(NNN)+DYGF(N))/(DYGF(N)))
              THETA = FLIMIT( R,CRGN,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NNN)+DYGF(N))
              DYF = 5.D-1
              VGO(M,NQY) = XGO(MN,NNN)*RHOG(MN,NNN)*VG(M,NQY)*THETA*DYF
     &          + XGO(MP,N)*RHOG(MP,N)*VG(M,NQY)*(1.D+0-THETA*DYF)
            ENDIF
            VGO(M,NQY) = VGO(M,NQY) + VDGO(M,NQY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            VDGO(M,NPY) = 0.D+0
            AG = MAX( VG(M,NPY),ZERO )
            VGO(M,NPY) = WTMO*XMGOB(MP,NB)*AG
!
!---      Upwind scheme  --
!
          ELSE
            VDGO(M,NPY) = 0.D+0
            FGOP = 0.D+0
            FGOS = XGOB(MP,NB)*RHOGB(MP,NB)
            INDX = -4
            FGO = DIFMN( FGOS,FGOP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
            VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
          ENDIF
        ENDIF

  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFGOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFGOW( N,NB )
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
!     Compute advective, diffusive, dispersive gas oil
!     fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddfgob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
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
      SUB_LOG(ISUB_LOG) = '/DDFGOW'
      IEQDX = IEQD
      IF( IOM.EQ.4 .OR. IOM.EQ.5 .OR. IOM.EQ.24 ) IEQDX = IEQDO
      I = ID(N)
      IZN = IZ(N)
      NPX = NSX(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        DFSP = TORG(MP,N)*(SG(MP,N)-SGT(MP,N))*PORD(MP,N)*DFGO(MP,N)
      &   *RHOMG(MP,N)
        DFSW = TORGB(MP,NB)*SGB(MP,NB)*PORDB(MP,NB)*DFGOB(MP,NB)
      &   *RHOMGB(MP,NB)
        INDX = 13
        DFSW = DIFMN(DFSW,DFSP,DXGF(N),DXGF(N),UG(1,NPX),INDX)
!
!---  Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVWB( PORD(MP,N),PORDB(MP,NB),SG(MP,N),SGB(MP,NB),
     &      UG,VG,WG,UGX,VGX,WGX,N,M )
          UGX = UGX*UGX
          VGX = VGX*VGX
          WGX = WGX*WGX
          ZGW = SQRT(UGX+VGX+WGX)
          DPSW = (DPLD(IZN)*UGX + DPTD(IZN)*(WGX+VGX))/(ZGW+SMALL)
        ELSE
          DPSW = 0.D+0
        ENDIF
!
!---   Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDGW = (DFSW+DPSW)/(5.D-1*DXGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
            VMCW = SGB(MP,NB)*PORDB(MP,NB)
            INDX = 13
            VMCW = DIFMN(VMCW,VMCP,DXGF(N),DXGF(N),UG(1,NPX),INDX)
            CRGW = ABS(UG(M,NPX))*DT/(DXGF(N)*VMCW+SMALL)
!
!---        On the boundary surface  --
!
            UDGO(M,NPX) = DDGW*(XMGOB(MP,NB)-XMGO(MP,N))
            IF( UG(1,NPX).GE.0.D+0 ) THEN
              UGO(M,NPX) = XGOB(MP,NB)*RHOGB(MP,NB)*UG(M,NPX)
            ELSEIF( UG(1,NPX).LT.0.D+0 .AND. I.LT.IFLD ) THEN
              NEN = N+1
              R = ((XMGO(1,N)-XMGO(1,NEN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NEN)))
              THETA = FLIMIT( R,CRGW,ISLC(8) )
              DXF = 5.D-1
              UGO(M,NPX) = XGOB(MP,NB)*RHOGB(MP,NB)*UG(M,NPX)*THETA*DXF
     &          + XGO(MP,N)*RHOG(MP,N)*UG(M,NPX)*(1.D+0-THETA*DXF)
            ELSEIF( UG(1,NPX).LT.0.D+0 .AND. I.EQ.IFLD ) THEN
              UGO(M,NPX) = XGO(MP,N)*RHOG(MP,N)*UG(M,NPX)
            ENDIF
            UGO(M,NPX) = UGO(M,NPX) + UDGO(M,NPX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQX = NPX+1
            NEN = N+1
            MN = MNEG(M)
            IF( UG(1,NQX).GE.0.D+0 ) THEN
              VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
              VMCE = (SG(MP,NEN)-SGT(MP,NEN))*PORD(MP,NEN)
              INDX = 13
              VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NEN),UG(1,NQX),INDX)
              CRGE = ABS( UG(M,NQX) )*DT/(DXGF(NEN)*VMCE+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NEN)-XMGO(1,N)+SMALL))
     &          *((DXGF(NEN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRGE,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NEN)+DXGF(N))
              DXF = 5.D-1
              UGO(M,NQX) = XGO(MN,NEN)*RHOG(MN,NEN)*UG(M,NQX)*THETA*DXF
     &          + XGO(MP,N)*RHOG(MP,N)*UG(M,NQX)*(1.D+0-THETA*DXF)
            ENDIF
            UGO(M,NQX) = UGO(M,NQX) + UDGO(M,NQX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AG = MAX( UG(M,NPX),ZERO ) +
     &       DDGW*MAX((ONE-(TENTH*ABS(UG(M,NPX))/(DDGW+SMALL)))**5,ZERO)
            AGP = MAX( -UG(M,NPX),ZERO ) +
     &       DDGW*MAX((ONE-(TENTH*ABS(UG(M,NPX))/(DDGW+SMALL)))**5,ZERO)
            UGO(M,NPX) = WTMO*(XMGOB(MP,NB)*AG-XMGO(MP,N)*AGP)
            UDGO(M,NPX) = DDGW*(XMGOB(MP,NB)-XMGO(MP,N))
!
!---      Upwind scheme  --
!
          ELSE
            UDGO(M,NPX) = DDGW*(XMGOB(MP,NB)-XMGO(MP,N))
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOW = XGOB(MP,NB)*RHOGB(MP,NB)
            INDX = -4
            FGO = DIFMN( FGOW,FGOP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
            UGO(M,NPX) = UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
          ENDIF
!
!---   Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
            VMCW = SGB(MP,NB)*PORDB(MP,NB)
            INDX = 13
            VMCW = DIFMN(VMCW,VMCP,DXGF(N),DXGF(N),UG(1,NPX),INDX)
            CRGW = ABS(UG(M,NPX))*DT/(DXGF(N)*VMCW+SMALL)
!
!---        On the boundary surface  --
!
            UDGO(M,NPX) = 0.D+0
            IF( UG(1,NPX).LT.0.D+0 .AND. I.LT.IFLD ) THEN
              NEN = N+1
              R = ((XMGO(1,N)-XMGO(1,NEN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NEN)))
              THETA = FLIMIT( R,CRGW,ISLC(8) )
              DXF = 5.D-1
              UGO(M,NPX) = XGOB(MP,NB)*RHOGB(MP,NB)*UG(M,NPX)*THETA*DXF
     &          + XGO(MP,N)*RHOG(MP,N)*UG(M,NPX)*(1.D+0-THETA*DXF)
            ELSEIF( UG(1,NPX).LT.0.D+0 .AND. I.EQ.IFLD ) THEN
              UGO(M,NPX) = XGO(MP,N)*RHOG(MP,N)*UG(M,NPX)
            ENDIF
            UGO(M,NPX) = UGO(M,NPX) + UDGO(M,NPX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQX = NPX+1
            NEN = N+1
            MN = MNEG(M)
            IF( UG(1,NQX).GE.0.D+0 ) THEN
              VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
              VMCE = (SG(MP,NEN)-SGT(MP,NEN))*PORD(MP,NEN)
              INDX = 13
              VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NEN),UG(1,NQX),INDX)
              CRGE = ABS( UG(M,NQX) )*DT/(DXGF(NEN)*VMCE+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NEN)-XMGO(1,N)+SMALL))
     &          *((DXGF(NEN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRGE,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NEN)+DXGF(N))
              DXF = 5.D-1
              UGO(M,NQX) = XGO(MN,NEN)*RHOG(MN,NEN)*UG(M,NQX)*THETA*DXF
     &          + XGO(MP,N)*RHOG(MP,N)*UG(M,NQX)*(1.D+0-THETA*DXF)
            ENDIF
            UGO(M,NQX) = UGO(M,NQX) + UDGO(M,NQX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            UDGO(M,NPX) = 0.D+0
            AGP = MAX( -UG(M,NPX),ZERO )
            UGO(M,NPX) = -WTMO*XMGO(MP,N)*AGP
!
!---      Upwind scheme  --
!
          ELSE
            UDGO(M,NPX) = 0.D+0
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOW = 0.D+0
            INDX = -4
            FGO = DIFMN( FGOW,FGOP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
            UGO(M,NPX) = UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
          ENDIF
!
!---   Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
            VMCW = SGB(MP,NB)*PORDB(MP,NB)
            INDX = 13
            VMCW = DIFMN(VMCW,VMCP,DXGF(N),DXGF(N),UG(1,NPX),INDX)
            CRGW = ABS(UG(M,NPX))*DT/(DXGF(N)*VMCW+SMALL)
!
!---        On the boundary surface  --
!
            UDGO(M,NPX) = 0.D+0
            IF( UG(1,NPX).GE.0.D+0 ) THEN
              UGO(M,NPX) = XGOB(MP,NB)*RHOGB(MP,NB)*UG(M,NPX)
            ENDIF
            UGO(M,NPX) = UGO(M,NPX) + UDGO(M,NPX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQX = NPX+1
            NEN = N+1
            MN = MNEG(M)
            IF( UG(1,NQX).GE.0.D+0 ) THEN
              VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
              VMCE = (SG(MP,NEN)-SGT(MP,NEN))*PORD(MP,NEN)
              INDX = 13
              VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NEN),UG(1,NQX),INDX)
              CRGE = ABS( UG(M,NQX) )*DT/(DXGF(NEN)*VMCE+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NEN)-XMGO(1,N)+SMALL))
     &          *((DXGF(NEN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRGE,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NEN)+DXGF(N))
              DXF = 5.D-1
              UGO(M,NQX) = XGO(MN,NEN)*RHOG(MN,NEN)*UG(M,NQX)*THETA*DXF
     &          + XGO(MP,N)*RHOG(MP,N)*UG(M,NQX)*(1.D+0-THETA*DXF)
            ENDIF
            UGO(M,NQX) = UGO(M,NQX) + UDGO(M,NQX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            UDGO(M,NPX) = 0.D+0
            AG = MAX( UG(M,NPX),ZERO )
            UGO(M,NPX) = WTMO*XMGOB(MP,NB)*AG
!
!---      Upwind scheme  --
!
          ELSE
            UDGO(M,NPX) = 0.D+0
            FGOP = 0.D+0
            FGOW = XGOB(MP,NB)*RHOGB(MP,NB)
            INDX = -4
            FGO = DIFMN( FGOW,FGOP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
            UGO(M,NPX) = UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFGOW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFGOE( N,NB )
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
!     Compute advective, diffusive, dispersive gas oil
!     fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddfgob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
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
      SUB_LOG(ISUB_LOG) = '/DDFGOE'
      IEQDX = IEQD
      IF( IOM.EQ.4 .OR. IOM.EQ.5 .OR. IOM.EQ.24 ) IEQDX = IEQDO
      I = ID(N)
      IZN = IZ(N)
      NQX = NSX(N)+1
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Diffusion coefficients  ---
!
        DFSP = TORG(MN,N)*(SG(MN,N)-SGT(MN,N))*PORD(MN,N)*DFGO(MN,N)
     &    *RHOMG(MN,N)
        DFSE = TORGB(MN,NB)*SGB(MN,NB)*PORDB(MN,NB)*DFGOB(MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFSE = DIFMN(DFSP,DFSE,DXGF(N),DXGF(N),UG(1,NQX),INDX)
!
!---    Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVEB( PORD(MN,N),PORDB(MN,NB),SG(MN,N),SGB(MN,NB),
     &      UG,VG,WG,UGX,VGX,WGX,N,M )
          UGX = UGX*UGX
          VGX = VGX*VGX
          WGX = WGX*WGX
          ZGE = SQRT(UGX+VGX+WGX)
          DPSE = (DPLD(IZN)*UGX + DPTD(IZN)*(WGX+VGX))/(ZGE+SMALL)
        ELSE
          DPSE = 0.D+0
        ENDIF
!
!---   Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDGE = (DFSE+DPSE)/(5.D-1*DXGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
            VMCE = SGB(MN,NB)*PORDB(MN,NB)
            INDX = 13
            VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(N),UG(1,NQX),INDX)
            CRGE = ABS(UG(M,NQX))*DT/(DXGF(N)*VMCE+SMALL)
!
!---        On the boundary surface  --
!
            UDGO(M,NQX) = DDGE*(XMGO(MN,N)-XMGOB(MN,NB))
            IF( UG(1,NQX).LT.0.D+0 ) THEN
              UGO(M,NQX) = XGOB(MN,NB)*RHOGB(MN,NB)*UG(M,NQX)
            ELSEIF( UG(1,NQX).GE.0.D+0 .AND. I.GT.1 ) THEN
              NWN = N-1
              R = ((XMGO(1,N)-XMGO(1,NWN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NWN)))
              THETA = FLIMIT( R,CRGE,ISLC(8) )
              DXF = 5.D-1
              UGO(M,NQX) = XGOB(MN,NB)*RHOGB(MN,NB)*UG(M,NQX)*THETA*DXF
     &          + XGO(MN,N)*RHOG(MN,N)*UG(M,NQX)*(1.D+0-THETA*DXF)
            ELSEIF( UG(1,NQX).GE.0.D+0 .AND. I.EQ.1 ) THEN
              UGO(M,NQX) = XGO(MN,N)*RHOG(MN,N)*UG(M,NQX)
            ENDIF
            UGO(M,NQX) = UGO(M,NQX) + UDGO(M,NQX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPX = NQX-1
            NWN = N-1
            MP = MPOS(M)
            IF( UG(1,NPX).LT.0.D+0 ) THEN
              VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
              VMCW = (SG(MN,NWN)-SGT(MN,NWN))*PORD(MN,NWN)
              INDX = 13
              VMCW = DIFMN(VMCW,VMCP,DXGF(NWN),DXGF(N),UG(1,NPX),INDX)
              CRGW = ABS( UG(M,NPX) )*DT/(DXGF(NWN)*VMCW+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NWN)-XMGO(1,N)+SMALL))
     &          *((DXGF(NWN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRGW,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NWN)+DXGF(N))
              DXF = 5.D-1
              UGO(M,NPX) = XGO(MN,NWN)*RHOG(MN,NWN)*UG(M,NPX)*THETA*DXF
     &          + XGO(MP,N)*RHOG(MP,N)*UG(M,NPX)*(1.D+0-THETA*DXF)
            ENDIF
            UGO(M,NPX) = UGO(M,NPX) + UDGO(M,NPX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AG = MAX( -UG(M,NQX),ZERO ) +
     &       DDGE*MAX((ONE-(TENTH*ABS(UG(M,NQX))/(DDGE+SMALL)))**5,ZERO)
            AGP = MAX( UG(M,NQX),ZERO ) +
     &       DDGE*MAX((ONE-(TENTH*ABS(UG(M,NQX))/(DDGE+SMALL)))**5,ZERO)
            UGO(M,NQX) = WTMO*(XMGO(MN,N)*AGP-XMGOB(MN,NB)*AG)
            UDGO(M,NQX) = DDGE*(XMGO(MN,N)-XMGOB(MN,NB))
!
!---      Upwind scheme  --
!
          ELSE
            UDGO(M,NQX) = DDGE*(XMGO(MN,N)-XMGOB(MN,NB))
            FGOP = XGO(MN,N)*RHOG(MN,N)
            FGOE = XGOB(MN,NB)*RHOGB(MN,NB)
            INDX = -4
            FGO = DIFMN( FGOP,FGOE,DXGF(N),DXGF(N),UG(1,NQX),INDX )
            UGO(M,NQX) = UG(M,NQX)*FGO + WTMO*UDGO(M,NQX)
          ENDIF
!
!---   Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
            VMCE = SGB(MN,NB)*PORDB(MN,NB)
            INDX = 13
            VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(N),UG(1,NQX),INDX)
            CRGE = ABS(UG(M,NQX))*DT/(DXGF(N)*VMCE+SMALL)
!
!---        On the boundary surface  --
!
            UDGO(M,NQX) = 0.D+0
            IF( UG(1,NQX).GE.0.D+0 .AND. I.GT.1 ) THEN
              NWN = N-1
              R = ((XMGO(1,N)-XMGO(1,NWN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NWN)))
              THETA = FLIMIT( R,CRGE,ISLC(8) )
              DXF = 5.D-1
              UGO(M,NQX) = XGOB(MN,NB)*RHOGB(MN,NB)*UG(M,NQX)*THETA*DXF
     &          + XGO(MN,N)*RHOG(MN,N)*UG(M,NQX)*(1.D+0-THETA*DXF)
            ELSEIF( UG(1,NQX).GE.0.D+0 .AND. I.EQ.1 ) THEN
              UGO(M,NQX) = XGO(MN,N)*RHOG(MN,N)*UG(M,NQX)
            ENDIF
            UGO(M,NQX) = UGO(M,NQX) + UDGO(M,NQX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPX = NQX-1
            NWN = N-1
            MP = MPOS(M)
            IF( UG(1,NPX).LT.0.D+0 ) THEN
              VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
              VMCW = (SG(MN,NWN)-SGT(MN,NWN))*PORD(MN,NWN)
              INDX = 13
              VMCW = DIFMN(VMCW,VMCP,DXGF(NWN),DXGF(N),UG(1,NPX),INDX)
              CRGW = ABS( UG(M,NPX) )*DT/(DXGF(NWN)*VMCW+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NWN)-XMGO(1,N)+SMALL))
     &          *((DXGF(NWN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRGW,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NWN)+DXGF(N))
              DXF = 5.D-1
              UGO(M,NPX) = XGO(MN,NB)*RHOG(MN,NB)*UG(M,NPX)*THETA*DXF
     &          + XGO(MP,N)*RHOG(MP,N)*UG(M,NPX)*(1.D+0-THETA*DXF)
            ENDIF
            UGO(M,NPX) = UGO(M,NPX) + UDGO(M,NPX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            UDGO(M,NQX) = 0.D+0
            AGP = MAX( UG(M,NQX),ZERO )
            UGO(M,NQX) = WTMO*XMGO(MN,N)*AGP
!
!---      Upwind scheme  --
!
          ELSE
            UDGO(M,NQX) = 0.D+0
            FGOP = XGO(MN,N)*RHOG(MN,N)
            FGOE = 0.D+0
            INDX = -4
            FGO = DIFMN( FGOP,FGOE,DXGF(N),DXGF(N),UG(1,NQX),INDX )
            UGO(M,NQX) = UG(M,NQX)*FGO + WTMO*UDGO(M,NQX)
          ENDIF
!
!---   Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
            VMCE = SGB(MN,NB)*PORDB(MN,NB)
            INDX = 13
            VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(N),UG(1,NQX),INDX)
            CRGE = ABS(UG(M,NQX))*DT/(DXGF(N)*VMCE+SMALL)
!
!---        On the boundary surface  --
!
            UDGO(M,NQX) = 0.D+0
            IF( UG(1,NQX).LT.0.D+0 ) THEN
              UGO(M,NQX) = XGOB(MN,NB)*RHOGB(MN,NB)*UG(M,NQX)
            ENDIF
            UGO(M,NQX) = UGO(M,NQX) + UDGO(M,NQX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPX = NQX-1
            NWN = N-1
            MP = MPOS(M)
            IF( UG(1,NPX).LT.0.D+0 ) THEN
              VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
              VMCW = (SG(MN,NWN)-SGT(MN,NWN))*PORD(MN,NWN)
              INDX = 13
              VMCW = DIFMN(VMCW,VMCP,DXGF(NWN),DXGF(N),UG(1,NPX),INDX)
              CRGW = ABS( UG(M,NPX) )*DT/(DXGF(NWN)*VMCW+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NWN)-XMGO(1,N)+SMALL))
     &          *((DXGF(NWN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRGW,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NWN)+DXGF(N))
              DXF = 5.D-1
              UGO(M,NPX) = XGO(MP,N)*RHOG(MP,N)*UG(M,NPX)*THETA*DXF
     &          + XGO(MN,NWN)*RHOG(MN,NWN)*UG(M,NPX)*(1.D+0-THETA*DXF)
            ENDIF
            UGO(M,NPX) = UGO(M,NPX) + UDGO(M,NPX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            UDGO(M,NQX) = 0.D+0
            AG = MAX( -UG(M,NQX),ZERO )
            UGO(M,NQX) = -WTMO*XMGOB(MN,NB)*AG
!
!---      Upwind scheme  --
!
          ELSE
            UDGO(M,NQX) = 0.D+0
            FGOP = 0.D+0
            FGOE = XGOB(MN,NB)*RHOGB(MN,NB)
            INDX = -4
            FGO = DIFMN( FGOP,FGOE,DXGF(N),DXGF(N),UG(1,NQX),INDX )
            UGO(M,NQX) = UG(M,NQX)*FGO + WTMO*UDGO(M,NQX)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFGOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFGON( N,NB )
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
!     Compute advective, diffusive, dispersive gas oil
!     fluxes on north boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddfgob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
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
      SUB_LOG(ISUB_LOG) = '/DDFGON'
      IEQDX = IEQD
      IF( IOM.EQ.4 .OR. IOM.EQ.5 .OR. IOM.EQ.24 ) IEQDX = IEQDO
      J = JD(N)
      I = ID(N)
      IZN = IZ(N)
      NQY = NSY(N)+IFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        DFSP = TORG(MN,N)*(SG(MN,N)-SGT(MN,N))*PORD(MN,N)*DFGO(MN,N)
     &    *RHOMG(MN,N)
        DFSN = TORGB(MN,NB)*SGB(MN,NB)*PORDB(MN,NB)*DFGOB(MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFSN = DIFMN(DFSP,DFSN,DYGF(N),DYGF(N),VG(1,NQY),INDX)
!
!---  Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVNB( PORD(MN,N),PORDB(MN,NB),SG(MN,N),SGB(MN,NB),
     &      UG,VG,WG,UGX,VGX,WGX,N,M )
          UGX = UGX*UGX
          VGX = VGX*VGX
          WGX = WGX*WGX
          ZGN = SQRT(UGX+VGX+WGX)
          DPSN = (DPLD(IZN)*VGX + DPTD(IZN)*(UGX+WGX))/(ZGN+SMALL)
        ELSE
          DPSN = 0.D+0
        ENDIF
!
!---   Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDGN = (DFSN+DPSN)/RP(I)/(5.D-1*DYGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
            VMCN = SGB(MN,NB)*PORDB(MN,NB)
            INDX = 13
            VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(N),VG(1,NQY),INDX)
            CRGN = ABS(VG(M,NQY))*DT/(DYGF(N)*VMCN+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDGO(M,NQY) = DDGN*(XMGO(MN,N)-XMGOB(MN,NB))
            IF( VG(1,NQY).LT.0.D+0 ) THEN
              VGO(M,NQY) = XGOB(MN,NB)*RHOGB(MN,NB)*VG(M,NQY)
            ELSEIF( VG(1,NQY).GE.0.D+0 .AND. J.GT.1 ) THEN
              NSN = N-IFLD
              R = ((XMGO(1,N)-XMGO(1,NSN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NSN)))
              THETA = FLIMIT( R,CRGN,ISLC(8) )
              DYF = 5.D-1
              VGO(M,NQY) = XGOB(MN,NB)*RHOGB(MN,NB)*VG(M,NQY)*THETA*DYF
     &          + XGO(MN,N)*RHOG(MN,N)*VG(M,NQY)*(1.D+0-THETA*DYF)
            ELSEIF( VG(1,NQY).GE.0.D+0 .AND. J.EQ.1 ) THEN
              VGO(M,NQY) = XGO(MN,N)*RHOG(MN,N)*VG(M,NQY)
            ENDIF
            VGO(M,NQY) = VGO(M,NQY) + VDGO(M,NQY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPY = NQY-IFLD
            NSN = N-IFLD
            MP = MPOS(M)
            IF( VG(1,NPY).LT.0.D+0 ) THEN
              VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
              VMCS = (SG(MN,NSN)-SGT(MN,NSN))*PORD(MN,NSN)
              INDX = 13
              VMCS = DIFMN(VMCS,VMCP,DYGF(NSN),DYGF(N),VG(1,NPY),INDX)
              CRGS = ABS( VG(M,NPY) )*DT/(DYGF(NSN)*VMCS+SMALL)/RP(I)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NSN)-XMGO(1,N)+SMALL))
     &          *((DXGF(NSN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRGS,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NSN)+DYGF(N))
              DYF = 5.D-1
              VGO(M,NPY) = XGO(MN,NSN)*RHOG(MN,NSN)*VG(M,NPY)*THETA*DYF
     &          + XGO(MP,N)*RHOG(MP,N)*VG(M,NPY)*(1.D+0-THETA*DYF)
            ENDIF
            VGO(M,NPY) = VGO(M,NPY) + VDGO(M,NPY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AG = MAX( -VG(M,NQY),ZERO ) +
     &       DDGN*MAX((ONE-(TENTH*ABS(VG(M,NQY))/(DDGN+SMALL)))**5,ZERO)
            AGP = MAX( VG(M,NQY),ZERO ) +
     &       DDGN*MAX((ONE-(TENTH*ABS(VG(M,NQY))/(DDGN+SMALL)))**5,ZERO)
            VGO(M,NQY) = WTMO*(XMGO(MN,N)*AGP-XMGOB(MN,NB)*AG)
            VDGO(M,NQY) = DDGN*(XMGO(MN,N)-XMGOB(MN,NB))
!
!---      Upwind scheme  --
!
          ELSE
            VDGO(M,NQY) = DDGN*(XMGO(MN,N)-XMGOB(MN,NB))
            FGOP = XGO(MN,N)*RHOG(MN,N)
            FGON = XGOB(MN,NB)*RHOGB(MN,NB)
            INDX = -4
            FGO = DIFMN( FGOP,FGON,DYGF(N),DYGF(N),VG(1,NQY),INDX )
            VGO(M,NQY) = VG(M,NQY)*FGO + WTMO*VDGO(M,NQY)
          ENDIF
!
!---   Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
            VMCN = SGB(MN,NB)*PORDB(MN,NB)
            INDX = 13
            VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(N),VG(1,NQY),INDX)
            CRGN = ABS(VG(M,NQY))*DT/(DYGF(N)*VMCN+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDGO(M,NQY) = 0.D+0
            IF( VG(1,NQY).GE.0.D+0 .AND. J.GT.1 ) THEN
              NSN = N-IFLD
              R = ((XMGO(1,N)-XMGO(1,NSN))
     &          /(XMGOB(1,NB)-XMGO(1,N))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NSN)))
              THETA = FLIMIT( R,CRGN,ISLC(8) )
              DYF = 5.D-1
              VGO(M,NQY) = XGOB(MN,NB)*RHOGB(MN,NB)*VG(M,NQY)*THETA*DYF
     &          + XGO(MN,N)*RHOG(MN,N)*VG(M,NQY)*(1.D+0-THETA*DYF)
            ELSEIF( VG(1,NQY).GE.0.D+0 .AND. J.EQ.1 ) THEN
              VGO(M,NQY) = XGO(MN,N)*RHOG(MN,N)*VG(M,NQY)
            ENDIF
            VGO(M,NQY) = VGO(M,NQY) + VDGO(M,NQY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPY = NQY-IFLD
            NSN = N-IFLD
            MP = MPOS(M)
            IF( VG(1,NPY).LT.0.D+0 ) THEN
              VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
              VMCS = (SG(MN,NSN)-SGT(MN,NSN))*PORD(MN,NSN)
              INDX = 13
              VMCS = DIFMN(VMCS,VMCP,DYGF(NSN),DYGF(N),VG(1,NPY),INDX)
              CRGS = ABS( VG(M,NPY) )*DT/(DYGF(NSN)*VMCS+SMALL)/RP(I)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NSN)-XMGO(1,N)+SMALL))
     &          *((DXGF(NSN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRGS,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NSN)+DYGF(N))
              DYF = 5.D-1
              VGO(M,NPY) = XGO(MN,NSN)*RHOG(MN,NSN)*VG(M,NPY)*THETA*DYF
     &          + XGO(MP,N)*RHOG(MP,N)*VG(M,NPY)*(1.D+0-THETA*DYF)
            ENDIF
            VGO(M,NPY) = VGO(M,NPY) + VDGO(M,NPY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            VDGO(M,NQY) = 0.D+0
            AGP = MAX( VG(M,NQY),ZERO )
            VGO(M,NQY) = WTMO*XMGO(MN,N)*AGP
!
!---      Upwind scheme  --
!
          ELSE
            VDGO(M,NQY) = 0.D+0
            FGOP = 0.D+0
            FGON = XGOB(MN,NB)*RHOGB(MN,NB)
            INDX = -4
            FGO = DIFMN( FGOP,FGON,DYGF(N),DYGF(N),VG(1,NQY),INDX )
            VGO(M,NQY) = VG(M,NQY)*FGO + WTMO*VDGO(M,NQY)
          ENDIF
!
!---   Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
            VMCN = SGB(MN,NB)*PORDB(MN,NB)
            INDX = 13
            VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(N),VG(1,NQY),INDX)
            CRGN = ABS(VG(M,NQY))*DT/(DYGF(N)*VMCN+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDGO(M,NQY) = 0.D+0
            IF( VG(1,NQY).LT.0.D+0 ) THEN
              VGO(M,NQY) = XGOB(MN,NB)*RHOGB(MN,NB)*VG(M,NQY)
            ENDIF
            VGO(M,NQY) = VGO(M,NQY) + VDGO(M,NQY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPY = NQY-IFLD
            NSN = N-IFLD
            MP = MPOS(M)
            IF( VG(1,NPY).LT.0.D+0 ) THEN
              VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
              VMCS = (SG(MN,NSN)-SGT(MN,NSN))*PORD(MN,NSN)
              INDX = 13
              VMCS = DIFMN(VMCS,VMCP,DYGF(NSN),DYGF(N),VG(1,NPY),INDX)
              CRGS = ABS( VG(M,NPY) )*DT/(DYGF(NSN)*VMCS+SMALL)/RP(I)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NSN)-XMGO(1,N)+SMALL))
     &          *((DXGF(NSN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRGS,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NSN)+DYGF(N))
              DYF = 5.D-1
              VGO(M,NPY) = XGO(MP,N)*RHOG(MP,N)*VG(M,NPY)*THETA*DYF
     &          + XGO(MN,NSN)*RHOG(MN,NSN)*VG(M,NPY)*(1.D+0-THETA*DYF)
            ENDIF
            VGO(M,NPY) = VGO(M,NPY) + VDGO(M,NPY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            VDGO(M,NQY) = 0.D+0
            AG = MAX( -VG(M,NQY),ZERO )
            VGO(M,NQY) = -WTMO*XMGOB(MN,NB)*AG
!
!---      Upwind scheme  --
!
          ELSE
            VDGO(M,NQY) = 0.D+0
            FGOP = 0.D+0
            FGON = XGOB(MN,NB)*RHOGB(MN,NB)
            INDX = -4
            FGO = DIFMN( FGOP,FGON,DYGF(N),DYGF(N),VG(1,NQY),INDX )
            VGO(M,NQY) = VG(M,NQY)*FGO + WTMO*VDGO(M,NQY)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFGON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFGOT( N,NB )
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
!     Compute advective, diffusive, dispersive gas oil
!     fluxes on top boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddfgob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
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
      SUB_LOG(ISUB_LOG) = '/DDFGOT'
      IEQDX = IEQD
      IF( IOM.EQ.4 .OR. IOM.EQ.5 .OR. IOM.EQ.24 ) IEQDX = IEQDO
      K = KD(N)
      IZN = IZ(N)
      NQZ = NSZ(N)+IJFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        DFSP = TORG(MN,N)*(SG(MN,N)-SGT(MN,N))*PORD(MN,N)*DFGO(MN,N)
     &    *RHOMG(MN,N)
        DFST = TORGB(MN,NB)*SGB(MN,NB)*PORDB(MN,NB)*DFGOB(MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFST = DIFMN(DFSP,DFST,DZGF(N),DZGF(N),WG(1,NQZ),INDX)
!
!---  Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVTB( PORD(MN,N),PORDB(MN,NB),SG(MN,N),SGB(MN,NB),
     &      UG,VG,WG,UGX,VGX,WGX,N,M )
          UGX = UGX*UGX
          VGX = VGX*VGX
          WGX = WGX*WGX
          ZGT = SQRT(UGX+VGX+WGX)
          DPST = (DPLD(IZN)*WGX + DPTD(IZN)*(UGX+VGX))/(ZGT+SMALL)
        ELSE
          DPST = 0.D+0
        ENDIF
!
!---   Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDGT = (DFST+DPST)/(5.D-1*DZGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
            VMCT = SGB(MN,NB)*PORDB(MN,NB)
            INDX = 13
            VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(N),WG(1,NQZ),INDX)
            CRGT = ABS(WG(M,NQZ))*DT/(DZGF(N)*VMCT+SMALL)
!
!---        On the boundary surface  --
!
            WDGO(M,NQZ) = DDGT*(XMGO(MN,N)-XMGOB(MN,NB))
            IF( WG(1,NQZ).LT.0.D+0 ) THEN
              WGO(M,NQZ) = XGOB(MN,NB)*RHOGB(MN,NB)*WG(M,NQZ)
            ELSEIF( WG(1,NQZ).GE.0.D+0 .AND. K.GT.1 ) THEN
              NBN = N-IJFLD
              R = ((XMGO(1,N)-XMGO(1,NBN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NBN)))
              THETA = FLIMIT( R,CRGT,ISLC(8) )
              DZF = 5.D-1
              WGO(M,NQZ) = XGOB(MN,NB)*RHOGB(MN,NB)*WG(M,NQZ)*THETA*DZF
     &          + XGO(MN,N)*RHOG(MN,N)*WG(M,NQZ)*(1.D+0-THETA*DZF)
            ELSEIF( WG(1,NQZ).GE.0.D+0 .AND. K.EQ.1 ) THEN
              WGO(M,NQZ) = XGO(MN,N)*RHOG(MN,N)*WG(M,NQZ)
            ENDIF
            WGO(M,NQZ) = WGO(M,NQZ) + WDGO(M,NQZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPZ = NQZ-IJFLD
            NBN = N-IJFLD
            MP = MPOS(M)
            IF( WG(1,NPZ).LT.0.D+0 ) THEN
              VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
              VMCB = (SG(MN,NBN)-SGT(MN,NBN))*PORD(MN,NBN)
              INDX = 13
              VMCB = DIFMN(VMCB,VMCP,DZGF(NBN),DZGF(N),WG(1,NPZ),INDX)
              CRGB = ABS( WG(M,NPZ) )*DT/(DZGF(NBN)*VMCB+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NBN)-XMGO(1,N)+SMALL))
     &          *((DZGF(NBN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRGB,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NBN)+DZGF(N))
              DZF = 5.D-1
              WGO(M,NPZ) = XGO(MN,NBN)*RHOG(MN,NBN)*WG(M,NPZ)*THETA*DZF
     &          + XGO(MP,N)*RHOG(MP,N)*WG(M,NPZ)*(1.D+0-THETA*DZF)
            ENDIF
            WGO(M,NPZ) = WGO(M,NPZ) + WDGO(M,NPZ)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AG = MAX( -WG(M,NQZ),ZERO ) +
     &       DDGT*MAX((ONE-(TENTH*ABS(WG(M,NQZ))/(DDGT+SMALL)))**5,ZERO)
            AGP = MAX( WG(M,NQZ),ZERO ) +
     &       DDGT*MAX((ONE-(TENTH*ABS(WG(M,NQZ))/(DDGT+SMALL)))**5,ZERO)
            WGO(M,NQZ) = WTMO*(XMGO(MN,N)*AGP-XMGOB(MN,NB)*AG)
            WDGO(M,NQZ) = DDGT*(XMGO(MN,N)-XMGOB(MN,NB))
!
!---      Upwind scheme  --
!
          ELSE
            WDGO(M,NQZ) = DDGT*(XMGO(MN,N)-XMGOB(MN,NB))
            FGOP = XGO(MN,N)*RHOG(MN,N)
            FGOT = XGOB(MN,NB)*RHOGB(MN,NB)
            INDX = -4
            FGO = DIFMN( FGOP,FGOT,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
            WGO(M,NQZ) = WG(M,NQZ)*FGO + WTMO*WDGO(M,NQZ)
          ENDIF
!
!---   Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
            VMCT = SGB(MN,NB)*PORDB(MN,NB)
            INDX = 13
            VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(N),WG(1,NQZ),INDX)
            CRGT = ABS(WG(M,NQZ))*DT/(DZGF(N)*VMCT+SMALL)
!
!---        On the boundary surface  --
!
            WDGO(M,NQZ) = 0.D+0
            IF( WG(1,NQZ).GE.0.D+0 .AND. K.GT.1 ) THEN
              NBN = N-IJFLD
              R = ((XMGO(1,N)-XMGO(1,NBN))
     &          /(XMGOB(1,NB)-XMGO(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NBN)))
              THETA = FLIMIT( R,CRGT,ISLC(8) )
              DZF = 5.D-1
              WGO(M,NQZ) = XGOB(MN,NB)*RHOGB(MN,NB)*WG(M,NQZ)*THETA*DZF
     &          + XGO(MN,N)*RHOG(MN,N)*WG(M,NQZ)*(1.D+0-THETA*DZF)
            ELSEIF( WG(1,NQZ).GE.0.D+0 .AND. K.EQ.1 ) THEN
              WGO(M,NQZ) = XGO(MN,N)*RHOG(MN,N)*WG(M,NQZ)
            ENDIF
            WGO(M,NQZ) = WGO(M,NQZ) + WDGO(M,NQZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPZ = NQZ-IJFLD
            NBN = N-IJFLD
            MP = MPOS(M)
            IF( WG(1,NPZ).LT.0.D+0 ) THEN
              VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
              VMCB = (SG(MN,NBN)-SGT(MN,NBN))*PORD(MN,NBN)
              INDX = 13
              VMCB = DIFMN(VMCB,VMCP,DZGF(NBN),DZGF(N),WG(1,NPZ),INDX)
              CRGB = ABS( WG(M,NPZ) )*DT/(DZGF(NBN)*VMCB+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NBN)-XMGO(1,N)+SMALL))
     &          *((DZGF(NBN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRGB,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NBN)+DZGF(N))
              DZF = 5.D-1
              WGO(M,NPZ) = XGO(MN,NBN)*RHOG(MN,NBN)*WG(M,NPZ)*THETA*DZF
     &          + XGO(MP,N)*RHOG(MP,N)*WG(M,NPZ)*(1.D+0-THETA*DZF)
            ENDIF
            WGO(M,NPZ) = WGO(M,NPZ) + WDGO(M,NPZ)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            WDGO(M,NQZ) = 0.D+0
            AGP = MAX( WG(M,NQZ),ZERO )
            WGO(M,NQZ) = WTMO*XMGO(MN,N)*AGP
!
!---      Upwind scheme  --
!
          ELSE
            WDGO(M,NQZ) = 0.D+0
            FGOP = XGO(MN,N)*RHOG(MN,N)
            FGOT = 0.D+0
            INDX = -4
            FGO = DIFMN( FGOP,FGOT,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
            WGO(M,NQZ) = WG(M,NQZ)*FGO + WTMO*WDGO(M,NQZ)
          ENDIF
!
!---   Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = (SG(MN,N)-SGT(MN,N))*PORD(MN,N)
            VMCT = SGB(MN,NB)*PORDB(MN,NB)
            INDX = 13
            VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(N),WG(1,NQZ),INDX)
            CRGT = ABS(WG(M,NQZ))*DT/(DZGF(N)*VMCT+SMALL)
!
!---        On the boundary surface  --
!
            WDGO(M,NQZ) = 0.D+0
            IF( WG(1,NQZ).LT.0.D+0 ) THEN
              WGO(M,NQZ) = XGOB(MN,NB)*RHOGB(MN,NB)*WG(M,NQZ)
            ENDIF
            WGO(M,NQZ) = WGO(M,NQZ) + WDGO(M,NQZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPZ = NQZ-IJFLD
            NBN = N-IJFLD
            MP = MPOS(M)
            IF( WG(1,NPZ).LT.0.D+0 ) THEN
              VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
              VMCB = (SG(MN,NBN)-SGT(MN,NBN))*PORD(MN,NBN)
              INDX = 13
              VMCB = DIFMN(VMCB,VMCP,DZGF(NBN),DZGF(N),WG(1,NPZ),INDX)
              CRGB = ABS( WG(M,NPZ) )*DT/(DZGF(NBN)*VMCB+SMALL)
              R = ((XMGO(1,N)-XMGOB(1,NB))
     &          /(XMGO(1,NBN)-XMGO(1,N)+SMALL))
     &          *((DZGF(NBN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRGB,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NBN)+DZGF(N))
              DZF = 5.D-1
              WGO(M,NPZ) = XGO(MN,NBN)*RHOG(MN,NBN)*WG(M,NPZ)*THETA*DZF
     &          + XGO(MP,N)*RHOG(MP,N)*WG(M,NPZ)*(1.D+0-THETA*DZF)
            ENDIF
            WGO(M,NPZ) = WGO(M,NPZ) + WDGO(M,NPZ)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            WDGO(M,NQZ) = 0.D+0
            AG = MAX( -WG(M,NQZ),ZERO )
            WGO(M,NQZ) = -WTMO*XMGOB(MN,NB)*AG
!
!---      Upwind scheme  --
!
          ELSE
            WDGO(M,NQZ) = 0.D+0
            FGOP = 0.D+0
            FGOT = XGOB(MN,NB)*RHOGB(MN,NB)
            INDX = -4
            FGO = DIFMN( FGOP,FGOT,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
            WGO(M,NQZ) = WG(M,NQZ)*FGO + WTMO*WDGO(M,NQZ)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFGOT group  ---
!
      RETURN
      END

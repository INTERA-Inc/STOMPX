!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFLOB( N,NB )
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
!     Compute advective, diffusive, dispersive aqueous oil
!     fluxes on bottom boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddflob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DDFLOB'
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
        DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &    *RHOML(MP,N)
        DFSB = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLOB(MP,NB)
     &    *RHOMLB(MP,NB)
        INDX = 15
        DFSB = DIFMN( DFSB,DFSP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVBB( PORD(MP,N),PORDB(MP,NB),SL(MP,N),SLB(MP,NB),
     &      UL,VL,WL,ULX,VLX,WLX,N,M )
          ULX = ULX*ULX
          VLX = VLX*VLX
          WLX = WLX*WLX
          ZLB = SQRT(ULX+VLX+WLX)
          DPSB = (DPLD(IZN)*WLX + DPTD(IZN)*(ULX+VLX))/(ZLB+SMALL)
        ELSE
          DPSB = 0.D+0
        ENDIF
!
!---   Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDLB = (DFSB+DPSB)/(5.D-1*DZGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MP,N)*PORD(MP,N)
            VMCB = SLB(MP,NB)*PORDB(MP,NB)
            INDX = 15
            VMCB = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WL(1,NPZ),INDX)
            CRLB = ABS(WL(M,NPZ))*DT/(DZGF(N)*VMCB+SMALL)
!
!---        On the boundary surface  --
!
            WDLO(M,NPZ) = DDLB*(XMLOB(MP,NB)-XMLO(MP,N))
            IF( WL(1,NPZ).GE.0.D+0 ) THEN
              WLO(M,NPZ) = XLOB(MP,NB)*RHOLB(MP,NB)*WL(M,NPZ)
            ELSEIF( WL(1,NPZ).LT.0.D+0 .AND. K.LT.KFLD ) THEN
              NTN = N+IJFLD
              R = ((XMLO(1,N)-XMLO(1,NTN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NTN)))
              THETA = FLIMIT( R,CRLB,ISLC(8) )
              DZF = 5.D-1
              WLO(M,NPZ) = XLOB(MP,NB)*RHOLB(MP,NB)*WL(M,NPZ)*THETA*DZF
     &          + XLO(MP,N)*RHOL(MP,N)*WL(M,NPZ)*(1.D+0-THETA*DZF)
            ELSEIF( WL(1,NPZ).LT.0.D+0 .AND. K.EQ.KFLD ) THEN
              WLO(M,NPZ) = XLO(MP,N)*RHOL(MP,N)*WL(M,NPZ)
            ENDIF
            WLO(M,NPZ) = WLO(M,NPZ) + WDLO(M,NPZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQZ = NPZ+IJFLD
            NTN = N+IJFLD
            MN = MNEG(M)
            IF( WL(1,NQZ).GE.0.D+0 ) THEN
              VMCP = SL(MN,N)*PORD(MN,N)
              VMCT = SL(MP,NTN)*PORD(MP,NTN)
              INDX = 15
              VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NTN),WL(1,NQZ),INDX)
              CRLT = ABS( WL(M,NQZ) )*DT/(DZGF(NTN)*VMCT+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NTN)-XMLO(1,N)+SMALL))
     &          *((DZGF(NTN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRLT,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NTN)+DZGF(N))
              DZF = 5.D-1
              WLO(M,NQZ) = XLO(MN,NTN)*RHOL(MN,NTN)*WL(M,NQZ)*THETA*DZF
     &          + XLO(MP,N)*RHOL(MP,N)*WL(M,NQZ)*(1.D+0-THETA*DZF)
            ENDIF
            WLO(M,NQZ) = WLO(M,NQZ) + WDLO(M,NQZ)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AL = MAX( WL(M,NPZ),ZERO ) +
     &       DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/(DDLB+SMALL)))**5,ZERO)
            ALP = MAX( -WL(M,NPZ),ZERO ) +
     &       DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/(DDLB+SMALL)))**5,ZERO)
            WLO(M,NPZ) = WTMO*(XMLOB(MP,NB)*AL -
     &        XMLO(MP,N)*ALP)
            WDLO(M,NPZ) = DDLB*(XMLOB(MP,NB) -
     &        XMLO(MP,N))
!
!---      Upwind scheme  --
!
          ELSE
            WDLO(M,NPZ) = DDLB*(XMLOB(MP,NB) -
     &        XMLO(MP,N))
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
            INDX = -4
            FLO = DIFMN( FLOB,FLOP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
            WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
          ENDIF
!
!---    Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MP,N)*PORD(MP,N)
            VMCB = SLB(MP,NB)*PORDB(MP,NB)
            INDX = 15
            VMCB = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WL(1,NPZ),INDX)
            CRLB = ABS(WL(M,NPZ))*DT/(DZGF(N)*VMCB+SMALL)
!
!---        On the boundary surface  --
!
            WDLO(M,NPZ) = 0.D+0
            IF( WL(1,NPZ).LT.0.D+0 .AND. K.LT.KFLD ) THEN
              NTN = N+IJFLD
              R = ((XMLO(1,N)-XMLO(1,NTN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NTN)))
              THETA = FLIMIT( R,CRLB,ISLC(8) )
              DZF = 5.D-1
              WLO(M,NPZ) = XLOB(MP,NB)*RHOLB(MP,NB)*WL(M,NPZ)*THETA*DZF
     &          + XLO(MP,N)*RHOL(MP,N)*WL(M,NPZ)*(1.D+0-THETA*DZF)
            ELSEIF( WL(1,NPZ).LT.0.D+0 .AND. K.EQ.KFLD ) THEN
              WLO(M,NPZ) = XLO(MP,N)*RHOL(MP,N)*WL(M,NPZ)
            ENDIF
            WLO(M,NPZ) = WLO(M,NPZ) + WDLO(M,NPZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQZ = NPZ+IJFLD
            NTN = N+IJFLD
            MN = MNEG(M)
            IF( WL(1,NQZ).GE.0.D+0 ) THEN
              VMCP = SL(MN,N)*PORD(MN,N)
              VMCT = SL(MP,NTN)*PORD(MP,NTN)
              INDX = 15
              VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NTN),WL(1,NQZ),INDX)
              CRLT = ABS( WL(M,NQZ) )*DT/(DZGF(NTN)*VMCT+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NTN)-XMLO(1,N)+SMALL))
     &          *((DZGF(NTN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRLT,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NTN)+DZGF(N))
              DZF = 5.D-1
              WLO(M,NQZ) = XLO(MN,NTN)*RHOL(MN,NTN)*WL(M,NQZ)*THETA*DZF
     &          + XLO(MP,N)*RHOL(MP,N)*WL(M,NQZ)*(1.D+0-THETA*DZF)
            ENDIF
            WLO(M,NQZ) = WLO(M,NQZ) + WDLO(M,NQZ)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            WDLO(M,NPZ) = 0.D+0
            ALP = MAX( -WL(M,NPZ),ZERO )
            WLO(M,NPZ) = -WTMO*XMLO(MP,N)*ALP
!
!---      Upwind scheme  --
!
          ELSE
            WDLO(M,NPZ) = 0.D+0
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOB = 0.D+0
            INDX = -4
            FLO = DIFMN( FLOB,FLOP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
            WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
          ENDIF
!
!---    Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MP,N)*PORD(MP,N)
            VMCB = SLB(MP,NB)*PORDB(MP,NB)
            INDX = 15
            VMCB = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WL(1,NPZ),INDX)
            CRLB = ABS(WL(M,NPZ))*DT/(DZGF(N)*VMCB+SMALL)
!
!---        On the boundary surface  --
!
            WDLO(M,NPZ) = 0.D+0
            IF( WL(1,NPZ).GE.0.D+0 ) THEN
              WLO(M,NPZ) = XLOB(MP,NB)*RHOLB(MP,NB)*WL(M,NPZ)
            ENDIF
            WLO(M,NPZ) = WLO(M,NPZ) + WDLO(M,NPZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQZ = NPZ+IJFLD
            NTN = N+IJFLD
            MN = MNEG(M)
            IF( WL(1,NQZ).GE.0.D+0 ) THEN
              VMCP = SL(MN,N)*PORD(MN,N)
              VMCT = SL(MP,NTN)*PORD(MP,NTN)
              INDX = 15
              VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NTN),WL(1,NQZ),INDX)
              CRLT = ABS( WL(M,NQZ) )*DT/(DZGF(NTN)*VMCT+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NTN)-XMLO(1,N)+SMALL))
     &          *((DZGF(NTN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRLT,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NTN)+DZGF(N))
              DZF = 5.D-1
              WLO(M,NQZ) = XLO(MN,NTN)*RHOL(MN,NTN)*WL(M,NQZ)*THETA*DZF
     &          + XLO(MP,N)*RHOL(MP,N)*WL(M,NQZ)*(1.D+0-THETA*DZF)
            ENDIF
            WLO(M,NQZ) = WLO(M,NQZ) + WDLO(M,NQZ)*WTMO
!
!---      Patankar scheme  ---
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            WDLO(M,NPZ) = 0.D+0
            AL = MAX( WL(M,NPZ),ZERO )
            WLO(M,NPZ) = WTMO*XMLOB(MP,NB)*AL
!
!---      Upwind scheme  ---
!
          ELSE
            WDLO(M,NPZ) = 0.D+0
            FLOP = 0.D+0
            FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
            INDX = -4
            FLO = DIFMN( FLOB,FLOP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
            WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFLOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFLOS( N,NB )
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
!     Compute advective, diffusive, dispersive aqueous oil
!     fluxes on south boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddflob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DDFLOS'
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
        DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &    *RHOML(MP,N)
        DFSS = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLOB(MP,NB)
     &    *RHOMLB(MP,NB)
        INDX = 15
        DFSS = DIFMN(DFSS,DFSP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
!
!---  Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVSB( PORD(MP,N),PORDB(MP,NB),SL(MP,N),SLB(MP,NB),
     &      UL,VL,WL,ULX,VLX,WLX,N,M )
          ULX = ULX*ULX
          VLX = VLX*VLX
          WLX = WLX*WLX
          ZLS = SQRT(ULX+VLX+WLX)
          DPSS = (DPLD(IZN)*VLX + DPTD(IZN)*(ULX+WLX))/(ZLS+SMALL)
        ELSE
          DPSS = 0.D+0
        ENDIF
!
!---    Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDLS = (DFSS+DPSS)/RP(I)/(5.D-1*DYGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MP,N)*PORD(MP,N)
            VMCS = SLB(MP,NB)*PORDB(MP,NB)
            INDX = 15
            VMCS = DIFMN(VMCS,VMCP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
            CRLS = ABS(VL(M,NPY))*DT/(DYGP(NPY)*VMCS+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDLO(M,NPY) = DDLS*(XMLOB(MP,NB) -
     &        XMLO(MP,N))
            IF( VL(1,NPY).GE.0.D+0 ) THEN
              VLO(M,NPY) = XLOB(MP,NB)*RHOLB(MP,NB)*VL(M,NPY)
            ELSEIF( VL(1,NPY).LT.0.D+0 .AND. J.LT.JFLD ) THEN
              NNN = N+IFLD
              R = ((XMLO(1,N)-XMLO(1,NNN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NNN)))
              THETA = FLIMIT( R,CRLS,ISLC(8) )
              DYF = 5.D-1
              VLO(M,NPY) = XLOB(MP,NB)*RHOLB(MP,NB)*VL(M,NPY)*THETA*DYF
     &          + XLO(MP,N)*RHOL(MP,N)*VL(M,NPY)*(1.D+0-THETA*DYF)
            ELSEIF( VL(1,NPY).LT.0.D+0 .AND. J.EQ.JFLD ) THEN
              VLO(M,NPY) = XLO(MP,N)*RHOL(MP,N)*VL(M,NPY)
            ENDIF
            VLO(M,NPY) = VLO(M,NPY) + VDLO(M,NPY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQY = NPY+IFLD
            NNN = N+IFLD
            MN = MNEG(M)
            IF( VL(1,NQY).GE.0.D+0 ) THEN
              VMCP = SL(MN,N)*PORD(MN,N)
              VMCN = SL(MP,NNN)*PORD(MP,NNN)
              INDX = 15
              VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(NNN),VL(1,NQY),INDX)
              CRLN = ABS( VL(M,NQY) )*DT/(DYGF(NNN)*VMCN+SMALL)/RP(I)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NNN)-XMLO(1,N)+SMALL))
     &          *((DYGF(NNN)+DYGF(N))/(DYGF(N)))
              THETA = FLIMIT( R,CRLN,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NNN)+DYGF(N))
              DYF = 5.D-1
              VLO(M,NQY) = XLO(MN,NNN)*RHOL(MN,NNN)*VL(M,NQY)*THETA*DYF
     &          + XLO(MP,N)*RHOL(MP,N)*VL(M,NQY)*(1.D+0-THETA*DYF)
            ENDIF
            VLO(M,NQY) = VLO(M,NQY) + VDLO(M,NQY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AL = MAX( VL(M,NPY),ZERO ) +
     &       DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/(DDLS+SMALL)))**5,ZERO)
            ALP = MAX( -VL(M,NPY),ZERO ) +
     &       DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/(DDLS+SMALL)))**5,ZERO)
            VLO(M,NPY) = WTMO*(XMLOB(MP,NB)*AL -
     &        XMLO(MP,N)*ALP)
            VDLO(M,NPY) = DDLS*(XMLOB(MP,NB) -
     &        XMLO(MP,N))
!
!---      Upwind scheme  --
!
          ELSE
            VDLO(M,NPY) = DDLS*(XMLOB(MP,NB) -
     &        XMLO(MP,N))
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOS = XLOB(MP,NB)*RHOLB(MP,NB)
            INDX = -4
            FLO = DIFMN( FLOS,FLOP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
            VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
          ENDIF
!
!---   Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MP,N)*PORD(MP,N)
            VMCS = SLB(MP,NB)*PORDB(MP,NB)
            INDX = 15
            VMCS = DIFMN(VMCS,VMCP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
            CRLS = ABS(VL(M,NPY))*DT/(DYGP(NPY)*VMCS+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDLO(M,NPY) = 0.D+0
            IF( VL(1,NPY).LT.0.D+0 .AND. J.LT.JFLD ) THEN
              NNN = N+IFLD
              R = ((XMLO(1,N)-XMLO(1,NNN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NNN)))
              THETA = FLIMIT( R,CRLS,ISLC(8) )
              DYF = 5.D-1
              VLO(M,NPY) = XLOB(MP,NB)*RHOLB(MP,NB)*VL(M,NPY)*THETA*DYF
     &          + XLO(MP,N)*RHOL(MP,N)*VL(M,NPY)*(1.D+0-THETA*DYF)
            ELSEIF( VL(1,NPY).LT.0.D+0 .AND. J.EQ.JFLD ) THEN
              VLO(M,NPY) = XLO(MP,N)*RHOL(MP,N)*VL(M,NPY)
            ENDIF
            VLO(M,NPY) = VLO(M,NPY) + VDLO(M,NPY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQY = NPY+IFLD
            NNN = N+IFLD
            MN = MNEG(M)
            IF( VL(1,NQY).GE.0.D+0 ) THEN
              VMCP = SL(MN,N)*PORD(MN,N)
              VMCN = SL(MP,NNN)*PORD(MP,NNN)
              INDX = 15
              VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(NNN),VL(1,NQY),INDX)
              CRLN = ABS( VL(M,NQY) )*DT/(DYGF(NNN)*VMCN+SMALL)/RP(I)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NNN)-XMLO(1,N)+SMALL))
     &          *((DYGF(NNN)+DYGF(N))/(DYGF(N)))
              THETA = FLIMIT( R,CRLN,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NNN)+DYGF(N))
              DYF = 5.D-1
              VLO(M,NQY) = XLO(MN,NNN)*RHOL(MN,NNN)*VL(M,NQY)*THETA*DYF
     &          + XLO(MP,N)*RHOL(MP,N)*VL(M,NQY)*(1.D+0-THETA*DYF)
            ENDIF
            VLO(M,NQY) = VLO(M,NQY) + VDLO(M,NQY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            VDLO(M,NPY) = 0.D+0
            ALP = MAX( -VL(M,NPY),ZERO )
            VLO(M,NPY) = -WTMO*XMLO(MP,N)*ALP
!
!---      Upwind scheme  --
!
          ELSE
            VDLO(M,NPY) = 0.D+0
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOS = 0.D+0
            INDX = -4
            FLO = DIFMN( FLOS,FLOP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
            VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
          ENDIF
!
!---   Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MP,N)*PORD(MP,N)
            VMCS = SLB(MP,NB)*PORDB(MP,NB)
            INDX = 15
            VMCS = DIFMN(VMCS,VMCP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
            CRLS = ABS(VL(M,NPY))*DT/(DYGP(NPY)*VMCS+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDLO(M,NPY) = 0.D+0
            IF( VL(1,NPY).GE.0.D+0 ) THEN
              VLO(M,NPY) = XLOB(MP,NB)*RHOLB(MP,NB)*VL(M,NPY)
            ENDIF
            VLO(M,NPY) = VLO(M,NPY) + VDLO(M,NPY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQY = NPY+IFLD
            NNN = N+IFLD
            MN = MNEG(M)
            IF( VL(1,NQY).GE.0.D+0 ) THEN
              VMCP = SL(MN,N)*PORD(MN,N)
              VMCN = SL(MP,NNN)*PORD(MP,NNN)
              INDX = 15
              VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(NNN),VL(1,NQY),INDX)
              CRLN = ABS( VL(M,NQY) )*DT/(DYGF(NNN)*VMCN+SMALL)/RP(I)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NNN)-XMLO(1,N)+SMALL))
     &          *((DYGF(NNN)+DYGF(N))/(DYGF(N)))
              THETA = FLIMIT( R,CRLN,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NNN)+DYGF(N))
              DYF = 5.D-1
              VLO(M,NQY) = XLO(MN,NNN)*RHOL(MN,NNN)*VL(M,NQY)*THETA*DYF
     &          + XLO(MP,N)*RHOL(MP,N)*VL(M,NQY)*(1.D+0-THETA*DYF)
            ENDIF
            VLO(M,NQY) = VLO(M,NQY) + VDLO(M,NQY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            VDLO(M,NPY) = 0.D+0
            AL = MAX( VL(M,NPY),ZERO )
            VLO(M,NPY) = WTMO*XMLOB(MP,NB)*AL
!
!---      Upwind scheme  --
!
          ELSE
            VDLO(M,NPY) = 0.D+0
            FLOP = 0.D+0
            FLOS = XLOB(MP,NB)*RHOLB(MP,NB)
            INDX = -4
            FLO = DIFMN( FLOS,FLOP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
            VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
          ENDIF
        ENDIF

  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFLOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFLOW( N,NB )
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
!     Compute advective, diffusive, dispersive aqueous oil
!     fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddflob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DDFLOW'
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
        DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &    *RHOML(MP,N)
        DFSW = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLOB(MP,NB)
     &    *RHOMLB(MP,NB)
        INDX = 15
        DFSW = DIFMN(DFSW,DFSP,DXGF(N),DXGF(N),UL(1,NPX),INDX)
!
!---  Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVWB( PORD(MP,N),PORDB(MP,NB),SL(MP,N),SLB(MP,NB),
     &      UL,VL,WL,ULX,VLX,WLX,N,M )
          ULX = ULX*ULX
          VLX = VLX*VLX
          WLX = WLX*WLX
          ZLW = SQRT(ULX+VLX+WLX)
          DPSW = (DPLD(IZN)*ULX + DPTD(IZN)*(WLX+VLX))/(ZLW+SMALL)
        ELSE
          DPSW = 0.D+0
        ENDIF
!
!---   Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDLW = (DFSW+DPSW)/(5.D-1*DXGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MP,N)*PORD(MP,N)
            VMCW = SLB(MP,NB)*PORDB(MP,NB)
            INDX = 15
            VMCW = DIFMN(VMCW,VMCP,DXGF(N),DXGF(N),UL(1,NPX),INDX)
            CRLW = ABS(UL(M,NPX))*DT/(DXGF(N)*VMCW+SMALL)
!
!---        On the boundary surface  --
!
            UDLO(M,NPX) = DDLW*(XMLOB(MP,NB) -
     &        XMLO(MP,N))
            IF( UL(1,NPX).GE.0.D+0 ) THEN
              ULO(M,NPX) = XLOB(MP,NB)*RHOLB(MP,NB)*UL(M,NPX)
            ELSEIF( UL(1,NPX).LT.0.D+0 .AND. I.LT.IFLD ) THEN
              NEN = N+1
              R = ((XMLO(1,N)-XMLO(1,NEN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NEN)))
              THETA = FLIMIT( R,CRLW,ISLC(8) )
              DXF = 5.D-1
              ULO(M,NPX) = XLOB(MP,NB)*RHOLB(MP,NB)*UL(M,NPX)*THETA*DXF
     &          + XLO(MP,N)*RHOL(MP,N)*UL(M,NPX)*(1.D+0-THETA*DXF)
            ELSEIF( UL(1,NPX).LT.0.D+0 .AND. I.EQ.IFLD ) THEN
              ULO(M,NPX) = XLO(MP,N)*RHOL(MP,N)*UL(M,NPX)
            ENDIF
            ULO(M,NPX) = ULO(M,NPX) + UDLO(M,NPX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQX = NPX+1
            NEN = N+1
            MN = MNEG(M)
            IF( UL(1,NQX).GE.0.D+0 ) THEN
              VMCP = SL(MN,N)*PORD(MN,N)
              VMCE = SL(MP,NEN)*PORD(MP,NEN)
              INDX = 15
              VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NEN),UL(1,NQX),INDX)
              CRLE = ABS( UL(M,NQX) )*DT/(DXGF(NEN)*VMCE+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NEN)-XMLO(1,N)+SMALL))
     &          *((DXGF(NEN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRLE,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NEN)+DXGF(N))
              DXF = 5.D-1
              ULO(M,NQX) = XLO(MN,NEN)*RHOL(MN,NEN)*UL(M,NQX)*THETA*DXF
     &          + XLO(MP,N)*RHOL(MP,N)*UL(M,NQX)*(1.D+0-THETA*DXF)
            ENDIF
            ULO(M,NQX) = ULO(M,NQX) + UDLO(M,NQX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AL = MAX( UL(M,NPX),ZERO ) +
     &       DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/(DDLW+SMALL)))**5,ZERO)
            ALP = MAX( -UL(M,NPX),ZERO ) +
     &       DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/(DDLW+SMALL)))**5,ZERO)
            ULO(M,NPX) = WTMO*(XMLOB(MP,NB)*AL -
     &        XMLO(MP,N)*ALP)
            UDLO(M,NPX) = DDLW*(XMLOB(MP,NB) -
     &        XMLO(MP,N))
!
!---      Upwind scheme  --
!
          ELSE
            UDLO(M,NPX) = DDLW*(XMLOB(MP,NB) -
     &        XMLO(MP,N))
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOW = XLOB(MP,NB)*RHOLB(MP,NB)
            INDX = -4
            FLO = DIFMN( FLOW,FLOP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
            ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
          ENDIF
!
!---   Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MP,N)*PORD(MP,N)
            VMCW = SLB(MP,NB)*PORDB(MP,NB)
            INDX = 15
            VMCW = DIFMN(VMCW,VMCP,DXGF(N),DXGF(N),UL(1,NPX),INDX)
            CRLW = ABS(UL(M,NPX))*DT/(DXGF(N)*VMCW+SMALL)
!
!---        On the boundary surface  --
!
            UDLO(M,NPX) = 0.D+0
            IF( UL(1,NPX).LT.0.D+0 .AND. I.LT.IFLD ) THEN
              NEN = N+1
              R = ((XMLO(1,N)-XMLO(1,NEN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NEN)))
              THETA = FLIMIT( R,CRLW,ISLC(8) )
              DXF = 5.D-1
              ULO(M,NPX) = XLOB(MP,NB)*RHOLB(MP,NB)*UL(M,NPX)*THETA*DXF
     &          + XLO(MP,N)*RHOL(MP,N)*UL(M,NPX)*(1.D+0-THETA*DXF)
            ELSEIF( UL(1,NPX).LT.0.D+0 .AND. I.EQ.IFLD ) THEN
              ULO(M,NPX) = XLO(MP,N)*RHOL(MP,N)*UL(M,NPX)
            ENDIF
            ULO(M,NPX) = ULO(M,NPX) + UDLO(M,NPX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQX = NPX+1
            NEN = N+1
            MN = MNEG(M)
            IF( UL(1,NQX).GE.0.D+0 ) THEN
              VMCP = SL(MN,N)*PORD(MN,N)
              VMCE = SL(MP,NEN)*PORD(MP,NEN)
              INDX = 15
              VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NEN),UL(1,NQX),INDX)
              CRLE = ABS( UL(M,NQX) )*DT/(DXGF(NEN)*VMCE+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NEN)-XMLO(1,N)+SMALL))
     &          *((DXGF(NEN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRLE,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NEN)+DXGF(N))
              DXF = 5.D-1
              ULO(M,NQX) = XLO(MN,NEN)*RHOL(MN,NEN)*UL(M,NQX)*THETA*DXF
     &          + XLO(MP,N)*RHOL(MP,N)*UL(M,NQX)*(1.D+0-THETA*DXF)
            ENDIF
            ULO(M,NQX) = ULO(M,NQX) + UDLO(M,NQX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            UDLO(M,NPX) = 0.D+0
            ALP = MAX( -UL(M,NPX),ZERO )
            ULO(M,NPX) = -WTMO*XMLO(MP,N)*ALP
!
!---      Upwind scheme  --
!
          ELSE
            UDLO(M,NPX) = 0.D+0
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOW = 0.D+0
            INDX = -4
            FLO = DIFMN( FLOW,FLOP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
            ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
          ENDIF
!
!---   Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MP,N)*PORD(MP,N)
            VMCW = SLB(MP,NB)*PORDB(MP,NB)
            INDX = 15
            VMCW = DIFMN(VMCW,VMCP,DXGF(N),DXGF(N),UL(1,NPX),INDX)
            CRLW = ABS(UL(M,NPX))*DT/(DXGF(N)*VMCW+SMALL)
!
!---        On the boundary surface  --
!
            UDLO(M,NPX) = 0.D+0
            IF( UL(1,NPX).GE.0.D+0 ) THEN
              ULO(M,NPX) = XLOB(MP,NB)*RHOLB(MP,NB)*UL(M,NPX)
            ENDIF
            ULO(M,NPX) = ULO(M,NPX) + UDLO(M,NPX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NQX = NPX+1
            NEN = N+1
            MN = MNEG(M)
            IF( UL(1,NQX).GE.0.D+0 ) THEN
              VMCP = SL(MN,N)*PORD(MN,N)
              VMCE = SL(MP,NEN)*PORD(MP,NEN)
              INDX = 15
              VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NEN),UL(1,NQX),INDX)
              CRLE = ABS( UL(M,NQX) )*DT/(DXGF(NEN)*VMCE+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NEN)-XMLO(1,N)+SMALL))
     &          *((DXGF(NEN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRLE,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NEN)+DXGF(N))
              DXF = 5.D-1
              ULO(M,NQX) = XLO(MN,NEN)*RHOL(MN,NEN)*UL(M,NQX)*THETA*DXF
     &          + XLO(MP,N)*RHOL(MP,N)*UL(M,NQX)*(1.D+0-THETA*DXF)
            ENDIF
            ULO(M,NQX) = ULO(M,NQX) + UDLO(M,NQX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            UDLO(M,NPX) = 0.D+0
            AL = MAX( UL(M,NPX),ZERO )
            ULO(M,NPX) = WTMO*XMLOB(MP,NB)*AL
!
!---      Upwind scheme  --
!
          ELSE
            UDLO(M,NPX) = 0.D+0
            FLOP = 0.D+0
            FLOW = XLOB(MP,NB)*RHOLB(MP,NB)
            INDX = -4
            FLO = DIFMN( FLOW,FLOP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
            ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFLOW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFLOE( N,NB )
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
!     Compute advective, diffusive, dispersive aqueous oil
!     fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddflob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DDFLOE'
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
        DFSP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLO(MN,N)
     &    *RHOML(MN,N)
        DFSE = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLOB(MN,NB)
     &    *RHOMLB(MN,NB)
        INDX = 15
        DFSE = DIFMN(DFSP,DFSE,DXGF(N),DXGF(N),UL(1,NQX),INDX)
!
!---    Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVEB( PORD(MN,N),PORDB(MN,NB),SL(MN,N),SLB(MN,NB),
     &      UL,VL,WL,ULX,VLX,WLX,N,M )
          ULX = ULX*ULX
          VLX = VLX*VLX
          WLX = WLX*WLX
          ZLE = SQRT(ULX+VLX+WLX)
          DPSE = (DPLD(IZN)*ULX + DPTD(IZN)*(WLX+VLX))/(ZLE+SMALL)
        ELSE
          DPSE = 0.D+0
        ENDIF
!
!---   Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDLE = (DFSE+DPSE)/(5.D-1*DXGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MN,N)*PORD(MN,N)
            VMCE = SLB(MN,NB)*PORDB(MN,NB)
            INDX = 15
            VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(N),UL(1,NQX),INDX)
            CRLE = ABS(UL(M,NQX))*DT/(DXGF(N)*VMCE+SMALL)
!
!---        On the boundary surface  --
!
            UDLO(M,NQX) = DDLE*(XMLO(MN,N)-XMLOB(MN,NB))
            IF( UL(1,NQX).LT.0.D+0 ) THEN
              ULO(M,NQX) = XLOB(MN,NB)*RHOLB(MN,NB)*UL(M,NQX)
            ELSEIF( UL(1,NQX).GE.0.D+0 .AND. I.GT.1 ) THEN
              NWN = N-1
              R = ((XMLO(1,N)-XMLO(1,NWN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NWN)))
              THETA = FLIMIT( R,CRLE,ISLC(8) )
              DXF = 5.D-1
              ULO(M,NQX) = XLOB(MN,NB)*RHOLB(MN,NB)*UL(M,NQX)*THETA*DXF
     &          + XLO(MN,N)*RHOL(MN,N)*UL(M,NQX)*(1.D+0-THETA*DXF)
            ELSEIF( UL(1,NQX).GE.0.D+0 .AND. I.EQ.1 ) THEN
              ULO(M,NQX) = XLO(MN,N)*RHOL(MN,N)*UL(M,NQX)
            ENDIF
            ULO(M,NQX) = ULO(M,NQX) + UDLO(M,NQX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPX = NQX-1
            NWN = N-1
            MP = MPOS(M)
            IF( UL(1,NPX).LT.0.D+0 ) THEN
              VMCP = SL(MP,N)*PORD(MP,N)
              VMCW = SL(MN,NWN)*PORD(MN,NWN)
              INDX = 15
              VMCW = DIFMN(VMCW,VMCP,DXGF(NWN),DXGF(N),UL(1,NPX),INDX)
              CRLW = ABS( UL(M,NPX) )*DT/(DXGF(NWN)*VMCW+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NWN)-XMLO(1,N)+SMALL))
     &          *((DXGF(NWN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRLW,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NWN)+DXGF(N))
              DXF = 5.D-1
              ULO(M,NPX) = XLO(MN,NWN)*RHOL(MN,NWN)*UL(M,NPX)*THETA*DXF
     &          + XLO(MP,N)*RHOL(MP,N)*UL(M,NPX)*(1.D+0-THETA*DXF)
            ENDIF
            ULO(M,NPX) = ULO(M,NPX) + UDLO(M,NPX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AL = MAX( -UL(M,NQX),ZERO ) +
     &       DDLE*MAX((ONE-(TENTH*ABS(UL(M,NQX))/(DDLE+SMALL)))**5,ZERO)
            ALP = MAX( UL(M,NQX),ZERO ) +
     &       DDLE*MAX((ONE-(TENTH*ABS(UL(M,NQX))/(DDLE+SMALL)))**5,ZERO)
            ULO(M,NQX) = WTMO*(XMLO(MN,N)*ALP -
     &       XMLOB(MN,NB)*AL)
            UDLO(M,NQX) = DDLE*(XMLO(MN,N) -
     &       XMLOB(MN,NB))
!
!---      Upwind scheme  --
!
          ELSE
            UDLO(M,NQX) = DDLE*(XMLO(MN,N) -
     &       XMLOB(MN,NB))
            FLOP = XLO(MN,N)*RHOL(MN,N)
            FLOE = XLOB(MN,NB)*RHOLB(MN,NB)
            INDX = -4
            FLO = DIFMN( FLOP,FLOE,DXGF(N),DXGF(N),UL(1,NQX),INDX )
            ULO(M,NQX) = UL(M,NQX)*FLO + WTMO*UDLO(M,NQX)
          ENDIF
!
!---   Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MN,N)*PORD(MN,N)
            VMCE = SLB(MN,NB)*PORDB(MN,NB)
            INDX = 15
            VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(N),UL(1,NQX),INDX)
            CRLE = ABS(UL(M,NQX))*DT/(DXGF(N)*VMCE+SMALL)
!
!---        On the boundary surface  --
!
            UDLO(M,NQX) = 0.D+0
            IF( UL(1,NQX).GE.0.D+0 .AND. I.GT.1 ) THEN
              NWN = N-1
              R = ((XMLO(1,N)-XMLO(1,NWN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NWN)))
              THETA = FLIMIT( R,CRLE,ISLC(8) )
              DXF = 5.D-1
              ULO(M,NQX) = XLOB(MN,NB)*RHOLB(MN,NB)*UL(M,NQX)*THETA*DXF
     &          + XLO(MN,N)*RHOL(MN,N)*UL(M,NQX)*(1.D+0-THETA*DXF)
            ELSEIF( UL(1,NQX).GE.0.D+0 .AND. I.EQ.1 ) THEN
              ULO(M,NQX) = XLO(MN,N)*RHOL(MN,N)*UL(M,NQX)
            ENDIF
            ULO(M,NQX) = ULO(M,NQX) + UDLO(M,NQX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPX = NQX-1
            NWN = N-1
            MP = MPOS(M)
            IF( UL(1,NPX).LT.0.D+0 ) THEN
              VMCP = SL(MP,N)*PORD(MP,N)
              VMCW = SL(MN,NWN)*PORD(MN,NWN)
              INDX = 15
              VMCW = DIFMN(VMCW,VMCP,DXGF(NWN),DXGF(N),UL(1,NPX),INDX)
              CRLW = ABS( UL(M,NPX) )*DT/(DXGF(NWN)*VMCW+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NWN)-XMLO(1,N)+SMALL))
     &          *((DXGF(NWN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRLW,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NWN)+DXGF(N))
              DXF = 5.D-1
              ULO(M,NPX) = XLO(MN,NB)*RHOL(MN,NB)*UL(M,NPX)*THETA*DXF
     &          + XLO(MP,N)*RHOL(MP,N)*UL(M,NPX)*(1.D+0-THETA*DXF)
            ENDIF
            ULO(M,NPX) = ULO(M,NPX) + UDLO(M,NPX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            UDLO(M,NQX) = 0.D+0
            ALP = MAX( UL(M,NQX),ZERO )
            ULO(M,NQX) = WTMO*XMLO(MN,N)*ALP
!
!---      Upwind scheme  --
!
          ELSE
            UDLO(M,NQX) = 0.D+0
            FLOP = XLO(MN,N)*RHOL(MN,N)
            FLOE = 0.D+0
            INDX = -4
            FLO = DIFMN( FLOP,FLOE,DXGF(N),DXGF(N),UL(1,NQX),INDX )
            ULO(M,NQX) = UL(M,NQX)*FLO + WTMO*UDLO(M,NQX)
          ENDIF
!
!---   Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MN,N)*PORD(MN,N)
            VMCE = SLB(MN,NB)*PORDB(MN,NB)
            INDX = 15
            VMCE = DIFMN(VMCP,VMCE,DXGF(N),DXGF(N),UL(1,NQX),INDX)
            CRLE = ABS(UL(M,NQX))*DT/(DXGF(N)*VMCE+SMALL)
!
!---        On the boundary surface  --
!
            UDLO(M,NQX) = 0.D+0
            IF( UL(1,NQX).LT.0.D+0 ) THEN
              ULO(M,NQX) = XLOB(MN,NB)*RHOLB(MN,NB)*UL(M,NQX)
            ENDIF
            ULO(M,NQX) = ULO(M,NQX) + UDLO(M,NQX)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPX = NQX-1
            NWN = N-1
            MP = MPOS(M)
            IF( UL(1,NPX).LT.0.D+0 ) THEN
              VMCP = SL(MP,N)*PORD(MP,N)
              VMCW = SL(MN,NWN)*PORD(MN,NWN)
              INDX = 15
              VMCW = DIFMN(VMCW,VMCP,DXGF(NWN),DXGF(N),UL(1,NPX),INDX)
              CRLW = ABS( UL(M,NPX) )*DT/(DXGF(NWN)*VMCW+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NWN)-XMLO(1,N)+SMALL))
     &          *((DXGF(NWN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRLW,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(NWN)+DXGF(N))
              DXF = 5.D-1
              ULO(M,NPX) = XLO(MP,N)*RHOL(MP,N)*UL(M,NPX)*THETA*DXF
     &          + XLO(MN,NWN)*RHOL(MN,NWN)*UL(M,NPX)*(1.D+0-THETA*DXF)
            ENDIF
            ULO(M,NPX) = ULO(M,NPX) + UDLO(M,NPX)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            UDLO(M,NQX) = 0.D+0
            AL = MAX( -UL(M,NQX),ZERO )
            ULO(M,NQX) = -WTMO*XMLOB(MN,NB)*AL
!
!---      Upwind scheme  --
!
          ELSE
            UDLO(M,NQX) = 0.D+0
            FLOP = 0.D+0
            FLOE = XLOB(MN,NB)*RHOLB(MN,NB)
            INDX = -4
            FLO = DIFMN( FLOP,FLOE,DXGF(N),DXGF(N),UL(1,NQX),INDX )
            ULO(M,NQX) = UL(M,NQX)*FLO + WTMO*UDLO(M,NQX)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFLOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFLON( N,NB )
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
!     Compute advective, diffusive, dispersive aqueous oil
!     fluxes on north boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddflob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DDFLON'
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
        DFSP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLO(MN,N)
     &    *RHOML(MN,N)
        DFSN = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLOB(MN,NB)
     &    *RHOMLB(MN,NB)
        INDX = 15
        DFSN = DIFMN(DFSP,DFSN,DYGF(N),DYGF(N),VL(1,NQY),INDX)
!
!---  Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVNB( PORD(MN,N),PORDB(MN,NB),SL(MN,N),SLB(MN,NB),
     &      UL,VL,WL,ULX,VLX,WLX,N,M )
          ULX = ULX*ULX
          VLX = VLX*VLX
          WLX = WLX*WLX
          ZLN = SQRT(ULX+VLX+WLX)
          DPSN = (DPLD(IZN)*VLX + DPTD(IZN)*(ULX+WLX))/(ZLN+SMALL)
        ELSE
          DPSN = 0.D+0
        ENDIF
!
!---   Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDLN = (DFSN+DPSN)/RP(I)/(5.D-1*DYGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MN,N)*PORD(MN,N)
            VMCN = SLB(MN,NB)*PORDB(MN,NB)
            INDX = 15
            VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(N),VL(1,NQY),INDX)
            CRLN = ABS(VL(M,NQY))*DT/(DYGF(N)*VMCN+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDLO(M,NQY) = DDLN*(XMLO(MN,N)-XMLOB(MN,NB))
            IF( VL(1,NQY).LT.0.D+0 ) THEN
              VLO(M,NQY) = XLOB(MN,NB)*RHOLB(MN,NB)*VL(M,NQY)
            ELSEIF( VL(1,NQY).GE.0.D+0 .AND. J.GT.1 ) THEN
              NSN = N-IFLD
              R = ((XMLO(1,N)-XMLO(1,NSN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NSN)))
              THETA = FLIMIT( R,CRLN,ISLC(8) )
              DYF = 5.D-1
              VLO(M,NQY) = XLOB(MN,NB)*RHOLB(MN,NB)*VL(M,NQY)*THETA*DYF
     &          + XLO(MN,N)*RHOL(MN,N)*VL(M,NQY)*(1.D+0-THETA*DYF)
            ELSEIF( VL(1,NQY).GE.0.D+0 .AND. J.EQ.1 ) THEN
              VLO(M,NQY) = XLO(MN,N)*RHOL(MN,N)*VL(M,NQY)
            ENDIF
            VLO(M,NQY) = VLO(M,NQY) + VDLO(M,NQY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPY = NQY-IFLD
            NSN = N-IFLD
            MP = MPOS(M)
            IF( VL(1,NPY).LT.0.D+0 ) THEN
              VMCP = SL(MP,N)*PORD(MP,N)
              VMCS = SL(MN,NSN)*PORD(MN,NSN)
              INDX = 15
              VMCS = DIFMN(VMCS,VMCP,DYGF(NSN),DYGF(N),VL(1,NPY),INDX)
              CRLS = ABS( VL(M,NPY) )*DT/(DYGF(NSN)*VMCS+SMALL)/RP(I)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NSN)-XMLO(1,N)+SMALL))
     &          *((DXGF(NSN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRLS,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NSN)+DYGF(N))
              DYF = 5.D-1
              VLO(M,NPY) = XLO(MN,NSN)*RHOL(MN,NSN)*VL(M,NPY)*THETA*DYF
     &          + XLO(MP,N)*RHOL(MP,N)*VL(M,NPY)*(1.D+0-THETA*DYF)
            ENDIF
            VLO(M,NPY) = VLO(M,NPY) + VDLO(M,NPY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AL = MAX( -VL(M,NQY),ZERO ) +
     &       DDLN*MAX((ONE-(TENTH*ABS(VL(M,NQY))/(DDLN+SMALL)))**5,ZERO)
            ALP = MAX( VL(M,NQY),ZERO ) +
     &       DDLN*MAX((ONE-(TENTH*ABS(VL(M,NQY))/(DDLN+SMALL)))**5,ZERO)
            VLO(M,NQY) = WTMO*(XMLO(MN,N)*ALP -
     &        XMLOB(MN,NB)*AL)
            VDLO(M,NQY) = DDLN*(XMLO(MN,N) -
     &        XMLOB(MN,NB))
!
!---      Upwind scheme  --
!
          ELSE
            VDLO(M,NQY) = DDLN*(XMLO(MN,N) -
     &        XMLOB(MN,NB))
            FLOP = XLO(MN,N)*RHOL(MN,N)
            FLON = XLOB(MN,NB)*RHOLB(MN,NB)
            INDX = -4
            FLO = DIFMN( FLOP,FLON,DYGF(N),DYGF(N),VL(1,NQY),INDX )
            VLO(M,NQY) = VL(M,NQY)*FLO + WTMO*VDLO(M,NQY)
          ENDIF
!
!---   Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MN,N)*PORD(MN,N)
            VMCN = SLB(MN,NB)*PORDB(MN,NB)
            INDX = 15
            VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(N),VL(1,NQY),INDX)
            CRLN = ABS(VL(M,NQY))*DT/(DYGF(N)*VMCN+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDLO(M,NQY) = 0.D+0
            IF( VL(1,NQY).GE.0.D+0 .AND. J.GT.1 ) THEN
              NSN = N-IFLD
              R = ((XMLO(1,N)-XMLO(1,NSN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NSN)))
              THETA = FLIMIT( R,CRLN,ISLC(8) )
              DYF = 5.D-1
              VLO(M,NQY) = XLOB(MN,NB)*RHOLB(MN,NB)*VL(M,NQY)*THETA*DYF
     &          + XLO(MN,N)*RHOL(MN,N)*VL(M,NQY)*(1.D+0-THETA*DYF)
            ELSEIF( VL(1,NQY).GE.0.D+0 .AND. J.EQ.1 ) THEN
              VLO(M,NQY) = XLO(MN,N)*RHOL(MN,N)*VL(M,NQY)
            ENDIF
            VLO(M,NQY) = VLO(M,NQY) + VDLO(M,NQY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPY = NQY-IFLD
            NSN = N-IFLD
            MP = MPOS(M)
            IF( VL(1,NPY).LT.0.D+0 ) THEN
              VMCP = SL(MP,N)*PORD(MP,N)
              VMCS = SL(MN,NSN)*PORD(MN,NSN)
              INDX = 15
              VMCS = DIFMN(VMCS,VMCP,DYGF(NSN),DYGF(N),VL(1,NPY),INDX)
              CRLS = ABS( VL(M,NPY) )*DT/(DYGF(NSN)*VMCS+SMALL)/RP(I)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NSN)-XMLO(1,N)+SMALL))
     &          *((DXGF(NSN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRLS,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NSN)+DYGF(N))
              DYF = 5.D-1
              VLO(M,NPY) = XLO(MN,NSN)*RHOL(MN,NSN)*VL(M,NPY)*THETA*DYF
     &          + XLO(MP,N)*RHOL(MP,N)*VL(M,NPY)*(1.D+0-THETA*DYF)
            ENDIF
            VLO(M,NPY) = VLO(M,NPY) + VDLO(M,NPY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            VDLO(M,NQY) = 0.D+0
            ALP = MAX( VL(M,NQY),ZERO )
            VLO(M,NQY) = WTMO*XMLO(MN,N)*ALP
!
!---      Upwind scheme  --
!
          ELSE
            VDLO(M,NQY) = 0.D+0
            FLOP = 0.D+0
            FLON = XLOB(MN,NB)*RHOLB(MN,NB)
            INDX = -4
            FLO = DIFMN( FLOP,FLON,DYGF(N),DYGF(N),VL(1,NQY),INDX )
            VLO(M,NQY) = VL(M,NQY)*FLO + WTMO*VDLO(M,NQY)
          ENDIF
!
!---   Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MN,N)*PORD(MN,N)
            VMCN = SLB(MN,NB)*PORDB(MN,NB)
            INDX = 15
            VMCN = DIFMN(VMCP,VMCN,DYGF(N),DYGF(N),VL(1,NQY),INDX)
            CRLN = ABS(VL(M,NQY))*DT/(DYGF(N)*VMCN+SMALL)/RP(I)
!
!---        On the boundary surface  --
!
            VDLO(M,NQY) = 0.D+0
            IF( VL(1,NQY).LT.0.D+0 ) THEN
              VLO(M,NQY) = XLOB(MN,NB)*RHOLB(MN,NB)*VL(M,NQY)
            ENDIF
            VLO(M,NQY) = VLO(M,NQY) + VDLO(M,NQY)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPY = NQY-IFLD
            NSN = N-IFLD
            MP = MPOS(M)
            IF( VL(1,NPY).LT.0.D+0 ) THEN
              VMCP = SL(MP,N)*PORD(MP,N)
              VMCS = SL(MN,NSN)*PORD(MN,NSN)
              INDX = 15
              VMCS = DIFMN(VMCS,VMCP,DYGF(NSN),DYGF(N),VL(1,NPY),INDX)
              CRLS = ABS( VL(M,NPY) )*DT/(DYGF(NSN)*VMCS+SMALL)/RP(I)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NSN)-XMLO(1,N)+SMALL))
     &          *((DXGF(NSN)+DXGF(N))/(DXGF(N)))
              THETA = FLIMIT( R,CRLS,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(NSN)+DYGF(N))
              DYF = 5.D-1
              VLO(M,NPY) = XLO(MP,N)*RHOL(MP,N)*VL(M,NPY)*THETA*DYF
     &          + XLO(MN,NSN)*RHOL(MN,NSN)*VL(M,NPY)*(1.D+0-THETA*DYF)
            ENDIF
            VLO(M,NPY) = VLO(M,NPY) + VDLO(M,NPY)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            VDLO(M,NQY) = 0.D+0
            AL = MAX( -VL(M,NQY),ZERO )
            VLO(M,NQY) = -WTMO*XMLOB(MN,NB)*AL
!
!---      Upwind scheme  --
!
          ELSE
            VDLO(M,NQY) = 0.D+0
            FLOP = 0.D+0
            FLON = XLOB(MN,NB)*RHOLB(MN,NB)
            INDX = -4
            FLO = DIFMN( FLOP,FLON,DYGF(N),DYGF(N),VL(1,NQY),INDX )
            VLO(M,NQY) = VL(M,NQY)*FLO + WTMO*VDLO(M,NQY)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFLON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFLOT( N,NB )
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
!     Compute advective, diffusive, dispersive aqueous oil
!     fluxes on top boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddflob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DDFLOT'
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
        DFSP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLO(MN,N)
     &    *RHOML(MN,N)
        DFST = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLOB(MN,NB)
     &    *RHOMLB(MN,NB)
        INDX = 15
        DFST = DIFMN(DFSP,DFST,DZGF(N),DZGF(N),WL(1,NQZ),INDX)
!
!---  Hydraulic dispersion  ---
!
        IF( IDSPD.EQ.1 ) THEN
          CALL ADVTB( PORD(MN,N),PORDB(MN,NB),SL(MN,N),SLB(MN,NB),
     &      UL,VL,WL,ULX,VLX,WLX,N,M )
          ULX = ULX*ULX
          VLX = VLX*VLX
          WLX = WLX*WLX
          ZLT = SQRT(ULX+VLX+WLX)
          DPST = (DPLD(IZN)*WLX + DPTD(IZN)*(ULX+VLX))/(ZLT+SMALL)
        ELSE
          DPST = 0.D+0
        ENDIF
!
!---   Dirichlet boundary  ---
!
        IF( IBCT(IEQDX,NB).EQ.1 .OR. IBCT(IEQDX,NB).EQ.8 .OR.
     &    IBCT(IEQDX,NB).EQ.12 ) THEN
          DDLT = (DFST+DPST)/(5.D-1*DZGF(N))
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MN,N)*PORD(MN,N)
            VMCT = SLB(MN,NB)*PORDB(MN,NB)
            INDX = 15
            VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(N),WL(1,NQZ),INDX)
            CRLT = ABS(WL(M,NQZ))*DT/(DZGF(N)*VMCT+SMALL)
!
!---        On the boundary surface  --
!
            WDLO(M,NQZ) = DDLT*(XMLO(MN,N)-XMLOB(MN,NB))
            IF( WL(1,NQZ).LT.0.D+0 ) THEN
              WLO(M,NQZ) = XLOB(MN,NB)*RHOLB(MN,NB)*WL(M,NQZ)
            ELSEIF( WL(1,NQZ).GE.0.D+0 .AND. K.GT.1 ) THEN
              NBN = N-IJFLD
              R = ((XMLO(1,N)-XMLO(1,NBN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NBN)))
              THETA = FLIMIT( R,CRLT,ISLC(8) )
              DZF = 5.D-1
              WLO(M,NQZ) = XLOB(MN,NB)*RHOLB(MN,NB)*WL(M,NQZ)*THETA*DZF
     &          + XLO(MN,N)*RHOL(MN,N)*WL(M,NQZ)*(1.D+0-THETA*DZF)
            ELSEIF( WL(1,NQZ).GE.0.D+0 .AND. K.EQ.1 ) THEN
              WLO(M,NQZ) = XLO(MN,N)*RHOL(MN,N)*WL(M,NQZ)
            ENDIF
            WLO(M,NQZ) = WLO(M,NQZ) + WDLO(M,NQZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPZ = NQZ-IJFLD
            NBN = N-IJFLD
            MP = MPOS(M)
            IF( WL(1,NPZ).LT.0.D+0 ) THEN
              VMCP = SL(MP,N)*PORD(MP,N)
              VMCB = SL(MN,NBN)*PORD(MN,NBN)
              INDX = 15
              VMCB = DIFMN(VMCB,VMCP,DZGF(NBN),DZGF(N),WL(1,NPZ),INDX)
              CRLB = ABS( WL(M,NPZ) )*DT/(DZGF(NBN)*VMCB+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NBN)-XMLO(1,N)+SMALL))
     &          *((DZGF(NBN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRLB,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NBN)+DZGF(N))
              DZF = 5.D-1
              WLO(M,NPZ) = XLO(MN,NBN)*RHOL(MN,NBN)*WL(M,NPZ)*THETA*DZF
     &          + XLO(MP,N)*RHOL(MP,N)*WL(M,NPZ)*(1.D+0-THETA*DZF)
            ENDIF
            WLO(M,NPZ) = WLO(M,NPZ) + WDLO(M,NPZ)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AL = MAX( -WL(M,NQZ),ZERO ) +
     &       DDLT*MAX((ONE-(TENTH*ABS(WL(M,NQZ))/(DDLT+SMALL)))**5,ZERO)
            ALP = MAX( WL(M,NQZ),ZERO ) +
     &       DDLT*MAX((ONE-(TENTH*ABS(WL(M,NQZ))/(DDLT+SMALL)))**5,ZERO)
            WLO(M,NQZ) = WTMO*(XMLO(MN,N)*ALP -
     &        XMLOB(MN,NB)*AL)
            WDLO(M,NQZ) = DDLT*(XMLO(MN,N) -
     &        XMLOB(MN,NB))
!
!---      Upwind scheme  --
!
          ELSE
            WDLO(M,NQZ) = DDLT*(XMLO(MN,N) -
     &        XMLOB(MN,NB))
            FLOP = XLO(MN,N)*RHOL(MN,N)
            FLOT = XLOB(MN,NB)*RHOLB(MN,NB)
            INDX = -4
            FLO = DIFMN( FLOP,FLOT,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
            WLO(M,NQZ) = WL(M,NQZ)*FLO + WTMO*WDLO(M,NQZ)
          ENDIF
!
!---   Outflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.7 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MN,N)*PORD(MN,N)
            VMCT = SLB(MN,NB)*PORDB(MN,NB)
            INDX = 15
            VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(N),WL(1,NQZ),INDX)
            CRLT = ABS(WL(M,NQZ))*DT/(DZGF(N)*VMCT+SMALL)
!
!---        On the boundary surface  --
!
            WDLO(M,NQZ) = 0.D+0
            IF( WL(1,NQZ).GE.0.D+0 .AND. K.GT.1 ) THEN
              NBN = N-IJFLD
              R = ((XMLO(1,N)-XMLO(1,NBN))
     &          /(XMLOB(1,NB)-XMLO(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NBN)))
              THETA = FLIMIT( R,CRLT,ISLC(8) )
              DZF = 5.D-1
              WLO(M,NQZ) = XLOB(MN,NB)*RHOLB(MN,NB)*WL(M,NQZ)*THETA*DZF
     &          + XLO(MN,N)*RHOL(MN,N)*WL(M,NQZ)*(1.D+0-THETA*DZF)
            ELSEIF( WL(1,NQZ).GE.0.D+0 .AND. K.EQ.1 ) THEN
              WLO(M,NQZ) = XLO(MN,N)*RHOL(MN,N)*WL(M,NQZ)
            ENDIF
            WLO(M,NQZ) = WLO(M,NQZ) + WDLO(M,NQZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPZ = NQZ-IJFLD
            NBN = N-IJFLD
            MP = MPOS(M)
            IF( WL(1,NPZ).LT.0.D+0 ) THEN
              VMCP = SL(MP,N)*PORD(MP,N)
              VMCB = SL(MN,NBN)*PORD(MN,NBN)
              INDX = 15
              VMCB = DIFMN(VMCB,VMCP,DZGF(NBN),DZGF(N),WL(1,NPZ),INDX)
              CRLB = ABS( WL(M,NPZ) )*DT/(DZGF(NBN)*VMCB+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NBN)-XMLO(1,N)+SMALL))
     &          *((DZGF(NBN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRLB,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NBN)+DZGF(N))
              DZF = 5.D-1
              WLO(M,NPZ) = XLO(MN,NBN)*RHOL(MN,NBN)*WL(M,NPZ)*THETA*DZF
     &          + XLO(MP,N)*RHOL(MP,N)*WL(M,NPZ)*(1.D+0-THETA*DZF)
            ENDIF
            WLO(M,NPZ) = WLO(M,NPZ) + WDLO(M,NPZ)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            WDLO(M,NQZ) = 0.D+0
            ALP = MAX( WL(M,NQZ),ZERO )
            WLO(M,NQZ) = WTMO*XMLO(MN,N)*ALP
!
!---      Upwind scheme  --
!
          ELSE
            WDLO(M,NQZ) = 0.D+0
            FLOP = XLO(MN,N)*RHOL(MN,N)
            FLOT = 0.D+0
            INDX = -4
            FLO = DIFMN( FLOP,FLOT,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
            WLO(M,NQZ) = WL(M,NQZ)*FLO + WTMO*WDLO(M,NQZ)
          ENDIF
!
!---   Inflow boundary  ---
!
        ELSEIF( IBCT(IEQDX,NB).EQ.13 .OR. IBCT(IEQDX,NB).EQ.14 ) THEN
!
!---      Leonard-TVD and Roe-Superbee schemes  --
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN          
            VMCP = SL(MN,N)*PORD(MN,N)
            VMCT = SLB(MN,NB)*PORDB(MN,NB)
            INDX = 15
            VMCT = DIFMN(VMCP,VMCT,DZGF(N),DZGF(N),WL(1,NQZ),INDX)
            CRLT = ABS(WL(M,NQZ))*DT/(DZGF(N)*VMCT+SMALL)
!
!---        On the boundary surface  --
!
            WDLO(M,NQZ) = 0.D+0
            IF( WL(1,NQZ).LT.0.D+0 ) THEN
              WLO(M,NQZ) = XLOB(MN,NB)*RHOLB(MN,NB)*WL(M,NQZ)
            ENDIF
            WLO(M,NQZ) = WLO(M,NQZ) + WDLO(M,NQZ)*WTMO
!
!---        On the interior surface adjacent to the boundary
!           if the flow direction is from the boundary  --
!
            NPZ = NQZ-IJFLD
            NBN = N-IJFLD
            MP = MPOS(M)
            IF( WL(1,NPZ).LT.0.D+0 ) THEN
              VMCP = SL(MP,N)*PORD(MP,N)
              VMCB = SL(MN,NBN)*PORD(MN,NBN)
              INDX = 15
              VMCB = DIFMN(VMCB,VMCP,DZGF(NBN),DZGF(N),WL(1,NPZ),INDX)
              CRLB = ABS( WL(M,NPZ) )*DT/(DZGF(NBN)*VMCB+SMALL)
              R = ((XMLO(1,N)-XMLOB(1,NB))
     &          /(XMLO(1,NBN)-XMLO(1,N)+SMALL))
     &          *((DZGF(NBN)+DZGF(N))/(DZGF(N)))
              THETA = FLIMIT( R,CRLB,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(NBN)+DZGF(N))
              DZF = 5.D-1
              WLO(M,NPZ) = XLO(MN,NBN)*RHOL(MN,NBN)*WL(M,NPZ)*THETA*DZF
     &          + XLO(MP,N)*RHOL(MP,N)*WL(M,NPZ)*(1.D+0-THETA*DZF)
            ENDIF
            WLO(M,NPZ) = WLO(M,NPZ) + WDLO(M,NPZ)*WTMO
!
!---      Patankar scheme  --
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            WDLO(M,NQZ) = 0.D+0
            AL = MAX( -WL(M,NQZ),ZERO )
            WLO(M,NQZ) = -WTMO*XMLOB(MN,NB)*AL
!
!---      Upwind scheme  --
!
          ELSE
            WDLO(M,NQZ) = 0.D+0
            FLOP = 0.D+0
            FLOT = XLOB(MN,NB)*RHOLB(MN,NB)
            INDX = -4
            FLO = DIFMN( FLOP,FLOT,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
            WLO(M,NQZ) = WL(M,NQZ)*FLO + WTMO*WDLO(M,NQZ)
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFLOT group  ---
!
      RETURN
      END

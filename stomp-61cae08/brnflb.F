!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNFLB( N,NB )
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
!     Compute salt aqueous-phase fluxes on bottom boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1995.
!     Last Modified by MD White, Battelle, PNL, October 17, 1995.
!     Last Modified by MD White, PNNL, 23 May 2002.
!     brnflb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      SUB_LOG(ISUB_LOG) = '/BRNFLB'
      K = KD(N)
      IZN = IZ(N)
      NPZ = NSZ(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MP,N)+TABS)/TSPRF
          SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
          TCOR = (TB(MP,NB)+TABS)/TSPRF
          SMDLB = DFLSB(MP,NB)*TCOR*(VISRL/VISLB(MP,NB))
          DFFLB = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*SMDLB
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZN))
          DFFLB = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SLB(MP,NB)*PORDB(MP,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
          DFFLB = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLSB(MP,NB)
        ENDIF
        INDX = 18
        DFFLB = DIFMN(DFFLB,DFFLP,DZGF(N),DZGF(N),WL(1,NPZ),INDX)
!
!---  Hydraulic dispersion  ---
!
        IF( IDSPS.EQ.1 ) THEN
          CALL ADVBB( PORD(MP,N),PORDB(MP,NB),SL(MP,N),SLB(MP,NB),
     &      UL,VL,WL,UBX,VBX,WBX,N,M )
          ULBX = UBX*UBX
          VLBX = VBX*VBX
          WLBX = WBX*WBX
          ZLB = SQRT(ULBX+VLBX+WLBX)
          DPLB = (DPLGS(IZN)*WLBX + DPTRS(IZN)*(ULBX+VLBX))/(ZLB+SMALL)
        ELSE
          DPLB = 0.D+0
        ENDIF
!
!---   Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 ) THEN
          DDLB = (DFFLB+DPLB)/(5.D-1*DZGF(N))
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            WDS(M,NPZ) = DDLB*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &        XLS(MP,N)*RHOL(MP,N))
            WS(M,NPZ) = XLSB(1,NB)*RHOLB(1,NB)*WL(1,NPZ)
            IF( WL(1,NPZ).LT.ZERO ) THEN
              NBT = N+IJFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBT)*RHOL(1,NBT))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              WS(M,NPZ) = XLSB(1,NB)*RHOLB(1,NB)*WL(1,NPZ)*THETA
     &          + XLS(1,N)*RHOL(1,N)*WL(1,NPZ)*(1.D+0-THETA)
            ENDIF
            NQZ = NPZ+IJFLD
            IF( WL(1,NQZ).GE.ZERO ) THEN
              NBT = N+IJFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLSB(1,NB)*RHOLB(1,NB))
     &          /(XLS(1,NBT)*RHOL(1,NBT)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DZGF(NBT)+DZGF(N))/DZGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
              WS(M,NQZ) = WS(M,NQZ)
     &          + XLS(1,N)*RHOL(1,N)*WL(1,NQZ)*(-THETA*DZF)
     &          + XLS(1,NBT)*RHOL(1,NBT)*WL(1,NQZ)*THETA*DZF
            ENDIF
            WS(M,NPZ) = WS(M,NPZ) + WDS(M,NPZ)
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( WL(M,NPZ),ZERO ) +
     &       DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/(DDLB+SMALL)))**5,ZERO)
            ALP = MAX( -WL(M,NPZ),ZERO ) +
     &       DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/(DDLB+SMALL)))**5,ZERO)
            WS(M,NPZ) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP)
            WDS(M,NPZ) = DDLB*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &        XLS(MP,N)*RHOL(MP,N))
          ENDIF
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            WS(M,NPZ) = 0.D+0
            IF( WL(1,NPZ).LT.ZERO ) THEN
              NBT = N+IJFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBT)*RHOL(1,NBT))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              WS(M,NPZ) = XLSB(1,NB)*RHOLB(1,NB)*WL(1,NPZ)*THETA
     &          + XLS(1,N)*RHOL(1,N)*WL(1,NPZ)*(1.D+0-THETA)
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            ALP = MAX( -WL(M,NPZ),ZERO )
            WS(M,NPZ) = -XLS(MP,N)*RHOL(MP,N)*ALP
          ENDIF
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            WS(M,NPZ) = XLSB(1,NB)*RHOLB(1,NB)*MAX( WL(1,NPZ),ZERO )
            NQZ = NPZ+IJFLD
            IF( WL(1,NQZ).GE.ZERO ) THEN
              NBT = N+IJFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLSB(1,NB)*RHOLB(1,NB))
     &          /(XLS(1,NBT)*RHOL(1,NBT)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DZGF(NBT)+DZGF(N))/DZGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
              WS(M,NQZ) = WS(M,NQZ)
     &          + XLS(1,N)*RHOL(1,N)*WL(1,NQZ)*(-THETA*DZF)
     &          + XLS(1,NBT)*RHOL(1,NBT)*WL(1,NQZ)*THETA*DZF
            ENDIF
!
!---  Patankar salt transport  ---
!
          ELSE
            AL = MAX( WL(M,NPZ),ZERO )
            WS(M,NPZ) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNFLB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNFLS( N,NB )
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
!     Compute salt aqueous-phase fluxes on south boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1995.
!     Last Modified by MD White, Battelle, PNL, October 17, 1995.
!     Last Modified by MD White, PNNL, 23 May 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      SUB_LOG(ISUB_LOG) = '/BRNFLS'
      J = JD(N)
      I = ID(N)
      IZN = IZ(N)
      NPY = NSY(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MP,N)+TABS)/TSPRF
          SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
          TCOR = (TB(MP,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MP,NB)*TCOR*(VISRL/VISLB(MP,NB))
          DFFLS = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*SMDLP
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZN))
          DFFLS = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SLB(MP,NB)*PORDB(MP,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
          DFFLS = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLSB(MP,NB)
        ENDIF
        INDX = 18
        DFFLS = DIFMN(DFFLS,DFFLP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
!
!---  Hydraulic dispersion
!
        IF( IDSPS.EQ.1 ) THEN
          CALL ADVSB( PORD(MP,N),PORDB(MP,NB),SL(MP,N),SLB(MP,NB),
     &      UL,VL,WL,USX,VSX,WSX,N,M )
          ULSX = USX*USX
          VLSX = VSX*VSX
          WLSX = WSX*WSX
          ZLS = SQRT(ULSX+VLSX+WLSX)
          DPLS = (DPLGS(IZN)*VLSX + DPTRS(IZN)*(ULSX+WLSX))/(ZLS+SMALL)
        ELSE
          DPLS = 0.D+0
        ENDIF
!
!---   Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 ) THEN
          DDLS = (DFFLS+DPLS)/RP(I)/(5.D-1*DYGF(N))
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            VDS(M,NPY) = DDLS*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &        XLS(MP,N)*RHOL(MP,N))
            VS(M,NPY) = XLSB(1,NB)*RHOLB(1,NB)*VL(1,NPY)
            IF( VL(1,NPY).LT.ZERO ) THEN
              NBN = N+IFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBN)*RHOL(1,NBN))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              VS(M,NPY) = XLSB(1,NB)*RHOLB(1,NB)*VL(1,NPY)*THETA
     &          + XLS(1,N)*RHOL(1,N)*VL(1,NPY)*(1.D+0-THETA)
            ENDIF
            NQY = NPY+IFLD
            IF( VL(1,NQY).GE.ZERO ) THEN
              NBN = N+IFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLSB(1,NB)*RHOLB(1,NB))
     &          /(XLS(1,NBN)*RHOL(1,NBN)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DYGF(NBN)+DYGF(N))/DYGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
              VS(M,NQY) = VS(M,NQY)
     &          + XLS(1,N)*RHOL(1,N)*VL(1,NQY)*(-THETA*DYF)
     &          + XLS(1,NBN)*RHOL(1,NBN)*VL(1,NQY)*THETA*DYF
            ENDIF
            VS(M,NPY) = VS(M,NPY) + VDS(M,NPY)
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( VL(M,NPY),ZERO ) +
     &       DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/(DDLS+SMALL)))**5,ZERO)
            ALP = MAX( -VL(M,NPY),ZERO ) +
     &       DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/(DDLS+SMALL)))**5,ZERO)
            VS(M,NPY) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP)
            VDS(M,NPY) = DDLS*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &        XLS(MP,N)*RHOL(MP,N))
          ENDIF
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            VS(M,NPY) = 0.D+0
            IF( VL(1,NPY).LT.ZERO ) THEN
              NBN = N+IFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBN)*RHOL(1,NBN))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              VS(M,NPY) = XLSB(1,NB)*RHOLB(1,NB)*VL(1,NPY)*THETA
     &          + XLS(1,N)*RHOL(1,N)*VL(1,NPY)*(1.D+0-THETA)
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            ALP = MAX( -VL(M,NPY),ZERO )
            VS(M,NPY) = -XLS(MP,N)*RHOL(MP,N)*ALP
          ENDIF
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            VS(M,NPY) = XLSB(1,NB)*RHOLB(1,NB)*MAX( VL(1,NPY),ZERO )
            NQY = NPY+IFLD
            IF( VL(1,NQY).GE.ZERO ) THEN
              NBN = N+IFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLSB(1,NB)*RHOLB(1,NB))
     &          /(XLS(1,NBN)*RHOL(1,NBN)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DYGF(NBN)+DYGF(N))/DYGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
              VS(M,NQY) = VS(M,NQY)
     &          + XLS(1,N)*RHOL(1,N)*VL(1,NQY)*(-THETA*DYF)
     &          + XLS(1,NBN)*RHOL(1,NBN)*VL(1,NQY)*THETA*DYF
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( VL(M,NPY),ZERO )
            VS(M,NPY) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
          ENDIF
        ENDIF

  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNFLS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNFLW( N,NB )
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
!     Compute salt aqueous-phase fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1995.
!     Last Modified by MD White, Battelle, PNL, October 17, 1995.
!     Last Modified by MD White, PNNL, 23 May 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      SUB_LOG(ISUB_LOG) = '/BRNFLW'
      I = ID(N)
      IZN = IZ(N)
      NPX = NSX(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MP,N)+TABS)/TSPRF
          SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
          TCOR = (TB(MP,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MP,NB)*TCOR*(VISRL/VISLB(MP,NB))
          DFFLW = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*SMDLP
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZN))
          DFFLW = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SLB(MP,NB)*PORDB(MP,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
          DFFLW = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLSB(MP,NB)
        ENDIF
        INDX = 18
        DFFLW = DIFMN(DFFLW,DFFLP,DXGF(N),DXGF(N),UL(1,NPX),INDX)
!
!---  Hydraulic dispersion
!
        IF( IDSPS.EQ.1 ) THEN
          CALL ADVWB( PORD(MP,N),PORDB(MP,NB),SL(MP,N),SLB(MP,NB),
     &      UL,VL,WL,UWX,VWX,WWX,N,M )
          ULX = UWX*UWX
          VLX = VWX*VWX
          WLX = WWX*WWX
          ZLW = SQRT(ULX+VLX+WLX)
          DPLW = (DPLGS(IZN)*ULX + DPTRS(IZN)*(WLX+VLX))/(ZLW+SMALL)
        ELSE
          DPLW = 0.D+0
        ENDIF
!
!---   Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 ) THEN
          DDLW = (DFFLW+DPLW)/(5.D-1*DXGF(N))
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            UDS(M,NPX) = DDLW*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &        XLS(MP,N)*RHOL(MP,N))
            US(M,NPX) = XLSB(1,NB)*RHOLB(1,NB)*UL(1,NPX)
            IF( UL(1,NPX).LT.ZERO ) THEN
              NBE = N+1
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBE)*RHOL(1,NBE))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              US(M,NPX) = XLSB(1,NB)*RHOLB(1,NB)*UL(1,NPX)*THETA
     &          + XLS(1,N)*RHOL(1,N)*UL(1,NPX)*(1.D+0-THETA)
            ENDIF
            NQX = NPX+1
            IF( UL(1,NQX).GE.ZERO ) THEN
              NBE = N+1
              R = ((XLS(1,N)*RHOL(1,N)-XLSB(1,NB)*RHOLB(1,NB))
     &          /(XLS(1,NBE)*RHOL(1,NBE)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DXGF(NBE)+DXGF(N))/DXGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
              US(M,NQX) = US(M,NQX)
     &          + XLS(1,N)*RHOL(1,N)*UL(1,NQX)*(-THETA*DXF)
     &          + XLS(1,NBE)*RHOL(1,NBE)*UL(1,NQX)*THETA*DXF
            ENDIF
            US(M,NPX) = US(M,NPX) + UDS(M,NPX)
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( UL(M,NPX),ZERO ) +
     &       DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/(DDLW+SMALL)))**5,ZERO)
            ALP = MAX( -UL(M,NPX),ZERO ) +
     &       DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/(DDLW+SMALL)))**5,ZERO)
            US(M,NPX) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP)
            UDS(M,NPX) = DDLW*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &        XLS(MP,N)*RHOL(MP,N))
          ENDIF
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            US(M,NPX) = 0.D+0
            IF( UL(1,NPX).LT.ZERO ) THEN
              NBE = N+1
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBE)*RHOL(1,NBE))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              US(M,NPX) = XLSB(1,NB)*RHOLB(1,NB)*UL(1,NPX)*THETA
     &          + XLS(1,N)*RHOL(1,N)*UL(1,NPX)*(1.D+0-THETA)
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            ALP = MAX( -UL(M,NPX),ZERO )
            US(M,NPX) = -XLS(MP,N)*RHOL(MP,N)*ALP
          ENDIF
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            US(M,NPX) = XLSB(1,NB)*RHOLB(1,NB)*MAX( UL(1,NPX),ZERO )
            NQX = NPX+1
            IF( UL(1,NQX).GE.ZERO ) THEN
              NBE = N+1
              R = ((XLS(1,N)*RHOL(1,N)-XLSB(1,NB)*RHOLB(1,NB))
     &          /(XLS(1,NBE)*RHOL(1,NBE)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DXGF(NBE)+DXGF(N))/DXGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
              US(M,NQX) = US(M,NQX)
     &          + XLS(1,N)*RHOL(1,N)*UL(1,NQX)*(-THETA*DXF)
     &          + XLS(1,NBE)*RHOL(1,NBE)*UL(1,NQX)*THETA*DXF
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( UL(M,NPX),ZERO )
            US(M,NPX) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNFLW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNFLE( N,NB )
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
!     Compute salt aqueous-phase fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1995.
!     Last Modified by MD White, Battelle, PNL, October 17, 1995.
!     Last Modified by MD White, PNNL, 23 May 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      SUB_LOG(ISUB_LOG) = '/BRNFLE'
      I = ID(N)
      IZN = IZ(N)
      NQX = NSX(N)+1
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MN,N)+TABS)/TSPRF
          SMDLP = DFLS(MN,N)*TCOR*(VISRL/VISL(MN,N))
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*SMDLP
          TCOR = (TB(MN,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MN,NB)*TCOR*(VISRL/VISLB(MN,NB))
          DFFLE = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*SMDLP
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SL(MN,N)*PORD(MN,N)*SDCLS(3,IZN))
          DFFLE = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SLB(MN,NB)*PORDB(MN,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLS(MN,N)
          DFFLE = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLSB(MN,NB)
        ENDIF
        INDX = 18
        DFFLE = DIFMN(DFFLP,DFFLE,DXGF(N),DXGF(N),UL(1,NQX),INDX)
!
!---  Hydraulic dispersion
!
        IF( IDSPS.EQ.1 ) THEN
          CALL ADVEB( PORD(MN,N),PORDB(MN,NB),SL(MN,N),SLB(MN,NB),
     &      UL,VL,WL,UEX,VEX,WEX,N,M )
          ULEX = UEX*UEX
          VLEX = VEX*VEX
          WLEX = WEX*WEX
          ZLE = SQRT(ULEX+VLEX+WLEX)
          DPLE = (DPLGS(IZN)*ULEX + DPTRS(IZN)*(WLEX+VLEX))/(ZLE+SMALL)
        ELSE
          DPLE = 0.D+0
        ENDIF
!
!---   Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 ) THEN
          DDLE = (DFFLE+DPLE)/(5.D-1*DXGF(N))
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            UDS(M,NQX) = DDLE*(XLSB(MN,NB)*RHOLB(MN,NB) -
     &        XLS(MN,N)*RHOL(MN,N))
            US(M,NQX) = XLSB(1,NB)*RHOLB(1,NB)*UL(1,NQX)
            IF( UL(1,NQX).GE.ZERO ) THEN
              NBW = N-1
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBW)*RHOL(1,NBW))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              US(1,NQX) = XLS(1,N)*RHOL(1,N)*UL(1,NQX)*(1.D+0-THETA)
     &          + XLSB(1,NB)*RHOLB(1,NB)*UL(1,NQX)*THETA
            ENDIF
            NPX = NSX(N)
            IF( UL(1,NPX).LT.ZERO ) THEN
              NBW = N-1
              R = ((XLS(1,N)*RHOL(1,N)-XLSB(1,NB)*RHOLB(1,NB))
     &          /(XLS(1,NBW)*RHOL(1,NBW)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DXGF(NBW)+DXGF(N))/DXGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
              US(M,NPX) = US(M,NPX)
     &          + XLS(1,N)*RHOL(1,N)*UL(1,NPX)*(-THETA*DXF)
     &          + XLS(1,NBW)*RHOL(1,NBW)*UL(1,NPX)*THETA*DXF
            ENDIF
            US(M,NQX) = US(M,NQX) + UDS(M,NQX)
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( -UL(M,NQX),ZERO ) +
     &       DDLE*MAX((ONE-(TENTH*ABS(UL(M,NQX))/(DDLE+SMALL)))**5,ZERO)
            ALP = MAX( UL(M,NQX),ZERO ) +
     &       DDLE*MAX((ONE-(TENTH*ABS(UL(M,NQX))/(DDLE+SMALL)))**5,ZERO)
            US(M,NQX) = (XLS(MN,N)*RHOL(MN,N)*ALP -
     &       XLSB(MN,NB)*RHOLB(MN,NB)*AL)
            UDS(M,NQX) = DDLE*(XLS(MN,N)*RHOL(MN,N) -
     &       XLSB(MN,NB)*RHOLB(MN,NB))
          ENDIF
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            US(M,NQX) = 0.D+0
            IF( UL(1,NQX).GE.ZERO ) THEN
              NBW = N-1
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBW)*RHOL(1,NBW))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              US(M,NQX) = XLS(1,N)*RHOL(1,N)*UL(1,NQX)*(1.D+0-THETA)
     &          + XLSB(1,NB)*RHOLB(1,NB)*UL(1,NQX)*THETA
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            ALP = MAX( UL(M,NQX),ZERO )
            US(M,NQX) = XLS(MN,N)*RHOL(MN,N)*ALP
          ENDIF
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            US(M,NQX) = XLSB(1,NB)*RHOLB(1,NB)*MIN( UL(1,NQX),ZERO )
            NPX = NSX(N)
            IF( UL(1,NPX).LT.ZERO ) THEN
              NBW = N-1
              R = ((XLS(1,N)*RHOL(1,N)-XLSB(1,NB)*RHOLB(1,NB))
     &          /(XLS(1,NBW)*RHOL(1,NBW)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DXGF(NBW)+DXGF(N))/DXGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
              US(M,NPX) = US(M,NPX)
     &          + XLS(1,N)*RHOL(1,N)*UL(1,NPX)*(-THETA*DXF)
     &          + XLS(1,NBW)*RHOL(1,NBW)*UL(1,NPX)*THETA*DXF
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( -UL(M,NQX),ZERO )
            US(M,NQX) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNFLE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNFLN( N,NB )
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
!     Compute salt aqueous-phase fluxes on north boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1995.
!     Last Modified by MD White, Battelle, PNL, October 17, 1995.
!     Last Modified by MD White, PNNL, 23 May 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      SUB_LOG(ISUB_LOG) = '/BRNFLN'
      J = JD(N)
      I = ID(N)
      IZN = IZ(N)
      NQY = NSY(N)+IFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MN,N)+TABS)/TSPRF
          SMDLP = DFLS(MN,N)*TCOR*(VISRL/VISL(MN,N))
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*SMDLP
          TCOR = (TB(MN,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MN,NB)*TCOR*(VISRL/VISLB(MN,NB))
          DFFLN = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*SMDLP
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SL(MN,N)*PORD(MN,N)*SDCLS(3,IZN))
          DFFLN = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SLB(MN,NB)*PORDB(MN,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLS(MN,N)
          DFFLN = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLSB(MN,NB)
        ENDIF
        INDX = 18
        DFFLN = DIFMN(DFFLP,DFFLN,DYGF(N),DYGF(N),VL(1,NQY),INDX)
!
!---  Hydraulic dispersion
!
        IF( IDSPS.EQ.1 ) THEN
          CALL ADVNB( PORD(MN,N),PORDB(MN,NB),SL(MN,N),SLB(MN,NB),
     &      UL,VL,WL,UNX,VNX,WNX,N,M )
          ULNX = UNX*UNX
          VLNX = VNX*VNX
          WLNX = WNX*WNX
          ZLN = SQRT(ULNX+VLNX+WLNX)
          DPLN = (DPLGS(IZN)*VLNX + DPTRS(IZN)*(ULNX+WLNX))/(ZLN+SMALL)
        ELSE
          DPLN = 0.D+0
        ENDIF
!
!---   Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 ) THEN
          DDLN = (DFFLN+DPLN)/RP(I)/(5.D-1*DYGF(N))
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            VDS(M,NQY) = DDLN*(XLSB(MN,NB)*RHOLB(MN,NB) -
     &        XLS(MN,N)*RHOL(MN,N))
            VS(M,NQY) = XLSB(1,NB)*RHOLB(1,NB)*VL(1,NQY)
            IF( VL(1,NQY).GE.ZERO ) THEN
              NBS = N-IFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBS)*RHOL(1,NBS))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              VS(M,NQY) = XLS(1,N)*RHOL(1,N)*VL(1,NQY)*(1.D+0-THETA)
     &          + XLSB(1,NB)*RHOLB(1,NB)*VL(1,NQY)*THETA
            ENDIF
            NPY = NSY(N)
            IF( VL(1,NPY).LT.ZERO ) THEN
              NBS = N-IFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLSB(1,NB)*RHOLB(1,NB))
     &          /(XLS(1,NBS)*RHOL(1,NBS)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DYGF(NBS)+DYGF(N))/DYGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
              VS(M,NPY) = VS(M,NPY)
     &          + XLS(1,N)*RHOL(1,N)*VL(1,NPY)*(-THETA*DYF)
     &          + XLS(1,NBS)*RHOL(1,NBS)*VL(1,NPY)*THETA*DYF
            ENDIF
            VS(M,NQY) = VS(M,NQY) + VDS(M,NQY)
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( -VL(M,NQY),ZERO ) +
     &       DDLN*MAX((ONE-(TENTH*ABS(VL(M,NQY))/(DDLN+SMALL)))**5,ZERO)
            ALP = MAX( VL(M,NQY),ZERO ) +
     &       DDLN*MAX((ONE-(TENTH*ABS(VL(M,NQY))/(DDLN+SMALL)))**5,ZERO)
            VS(M,NQY) = (XLS(MN,N)*RHOL(MN,N)*ALP -
     &        XLSB(MN,NB)*RHOLB(MN,NB)*AL)
            VDS(M,NQY) = DDLN*(XLS(MN,N)*RHOL(MN,N) -
     &        XLSB(MN,NB)*RHOLB(MN,NB))
          ENDIF
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            VS(M,NQY) = 0.D+0
            IF( VL(1,NQY).GE.ZERO ) THEN
              NBS = N-IFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBS)*RHOL(1,NBS))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              VS(M,NQY) = XLS(1,N)*RHOL(1,N)*VL(1,NQY)*(1.D+0-THETA)
     &          + XLSB(1,NB)*RHOLB(1,NB)*VL(1,NQY)*THETA
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            ALP = MAX( VL(M,NQY),ZERO )
            VS(M,NQY) = XLS(MN,N)*RHOL(MN,N)*ALP
          ENDIF
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            VS(M,NQY) = XLSB(1,NB)*RHOLB(1,NB)*MIN( VL(1,NQY),ZERO )
            NPY = NSY(N)
            IF( VL(1,NPY).LT.ZERO ) THEN
              NBS = N-IFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLSB(1,NB)*RHOLB(1,NB))
     &          /(XLS(1,NBS)*RHOL(1,NBS)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DYGF(NBS)+DYGF(N))/DYGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
              VS(M,NPY) = VS(M,NPY)
     &          + XLS(1,N)*RHOL(1,N)*VL(1,NPY)*(-THETA*DYF)
     &          + XLS(1,NBS)*RHOL(1,NBS)*VL(1,NPY)*THETA*DYF
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( -VL(M,NQY),ZERO )
            VS(M,NQY) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNFLN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNFLT( N,NB )
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
!     Compute salt aqueous-phase fluxes on top boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1995.
!     Last Modified by MD White, Battelle, PNL, October 17, 1995.
!     Last Modified by MD White, PNNL, 23 May 2002.
!     brnflb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      SUB_LOG(ISUB_LOG) = '/BRNFLT'
      K = KD(N)
      IZN = IZ(N)
      NQZ = NSZ(N)+IJFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MN,N)+TABS)/TSPRF
          SMDLP = DFLS(MN,N)*TCOR*(VISRL/VISL(MN,N))
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*SMDLP
          TCOR = (TB(MN,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MN,NB)*TCOR*(VISRL/VISLB(MN,NB))
          DFFLT = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*SMDLP
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SL(MN,N)*PORD(MN,N)*SDCLS(3,IZN))
          DFFLT = SDCLS(1,IZN)*SDCLS(2,IZN)*
     &      EXP(SLB(MN,NB)*PORDB(MN,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLS(MN,N)
          DFFLT = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLSB(MN,NB)
        ENDIF
        INDX = 18
        DFFLT = DIFMN(DFFLP,DFFLT,DZGF(N),DZGF(N),WL(1,NQZ),INDX)
!
!---  Hydraulic dispersion
!
        IF( IDSPS.EQ.1 ) THEN
          CALL ADVTB( PORD(MN,N),PORDB(MN,NB),SL(MN,N),SLB(MN,NB),
     &      UL,VL,WL,UTX,VTX,WTX,N,M )
          ULTX = UTX*UTX
          VLTX = VTX*VTX
          WLTX = WTX*WTX
          ZLT = SQRT(ULTX+VLTX+WLTX)
          DPLT = (DPLGS(IZN)*WLTX + DPTRS(IZN)*(ULTX+VLTX))/(ZLT+SMALL)
        ELSE
          DPLT = 0.D+0
        ENDIF
!
!---   Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 ) THEN
          DDLT = (DFFLT+DPLT)/(5.D-1*DZGF(N))
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            WDS(M,NQZ) = DDLT*(XLSB(MN,NB)*RHOLB(MN,NB) -
     &        XLS(MN,N)*RHOL(MN,N))
            WS(M,NQZ) = XLSB(1,NB)*RHOLB(1,NB)*WL(1,NQZ)
            IF( WL(1,NQZ).GE.ZERO ) THEN
              NBB = N-IJFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBB)*RHOL(1,NBB))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              WS(M,NQZ) = XLS(1,N)*RHOL(1,N)*WL(1,NQZ)*(1.D+0-THETA)
     &          + XLSB(1,NB)*RHOLB(1,NB)*WL(1,NQZ)*THETA
            ENDIF
            NPZ = NSZ(N)
            IF( WL(1,NPZ).LT.ZERO ) THEN
              NBB = N-IJFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBB)*RHOL(1,NBB))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DZGF(NBB)+DZGF(N))/DZGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
              WS(M,NPZ) = WS(M,NPZ)
     &          + XLS(1,N)*RHOL(1,N)*WL(1,NPZ)*(-THETA*DZF)
     &          + XLS(1,NBB)*RHOL(1,NBB)*WL(1,NPZ)*THETA*DZF
            ENDIF
            WS(M,NQZ) = WS(M,NQZ) + WDS(M,NQZ)
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( -WL(M,NQZ),ZERO ) +
     &       DDLT*MAX((ONE-(TENTH*ABS(WL(M,NQZ))/(DDLT+SMALL)))**5,ZERO)
            ALP = MAX( WL(M,NQZ),ZERO ) +
     &       DDLT*MAX((ONE-(TENTH*ABS(WL(M,NQZ))/(DDLT+SMALL)))**5,ZERO)
            WS(M,NQZ) = WS(M,NQZ) + (XLS(MN,N)*RHOL(MN,N)*ALP -
     &        XLSB(MN,NB)*RHOLB(MN,NB)*AL)
            WDS(M,NQZ) = DDLT*(XLS(MN,N)*RHOL(MN,N) -
     &        XLSB(MN,NB)*RHOLB(MN,NB))
          ENDIF
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            WS(M,NQZ) = 0.D+0
            IF( WL(1,NQZ).GE.ZERO ) THEN
              NBB = N-IJFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBB)*RHOL(1,NBB))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              WS(M,NQZ) = XLS(1,N)*RHOL(1,N)*WL(1,NQZ)*(1.D+0-THETA)
     &          + XLSB(1,NB)*RHOLB(1,NB)*WL(1,NQZ)*THETA
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            ALP = MAX( WL(M,NQZ),ZERO )
            WS(M,NQZ) = XLS(MN,N)*RHOL(MN,N)*ALP
          ENDIF
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            WS(M,NQZ) = XLSB(1,NB)*RHOLB(1,NB)*MIN( WL(1,NQZ),ZERO )
            NPZ = NSZ(N)
            IF( WL(1,NPZ).LT.ZERO ) THEN
              NBB = N-IJFLD
              R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NBB)*RHOL(1,NBB))
     &          /(XLSB(1,NB)*RHOLB(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &          *((DZGF(NBB)+DZGF(N))/DZGF(N))
              THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
              DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
              WS(M,NPZ) = WS(M,NPZ)
     &          + XLS(1,N)*RHOL(1,N)*WL(1,NPZ)*(-THETA*DZF)
     &          + XLS(1,NBB)*RHOL(1,NBB)*WL(1,NPZ)*THETA*DZF
            ENDIF
!
!---  Patankar salt transport  --
!
          ELSE
            AL = MAX( -WL(M,NQZ),ZERO )
            WS(M,NQZ) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNFLT group  ---
!
      RETURN
      END

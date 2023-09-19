!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNFL
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
!     Compute salt aqueous-phase fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1995.
!     Last Modified by MD White, Battelle, PNL, October 17, 1995.
!     Last Modified by MD White, PNNL, 20 May 2002.
!     brnfl.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNFL'
!
!---  X-direction Darcy velocities, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
        DO 200 J = 1,JFLD
          DO 200 I = 2,IFLD
            N = ND(I,J,K)
            NW = N-1
            IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &        INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
            NPX = NSX(N)
            DO 100 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
              IF( IEDLS.EQ.1 ) THEN
                TCOR = (T(MP,N)+TABS)/TSPRF
                SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
                DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
                TCOR = (T(MN,NW)+TABS)/TSPRF
                SMDLP = DFLS(MN,NW)*TCOR*(VISRL/VISL(MN,NW))
                DFFLW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*SMDLP
              ELSEIF( IEDLS.EQ.2 ) THEN
                DFFLP = SDCLS(1,IZ(N))*SDCLS(2,IZ(N))*
     &            EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
                DFFLW = SDCLS(1,IZ(NW))*SDCLS(2,IZ(NW))*
     &            EXP(SL(MN,NW)*PORD(MN,NW)*SDCLS(3,IZ(NW)))
              ELSEIF( IEDLS.EQ.3 ) THEN
                DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
                DFFLW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*DFLS(MN,NW)
              ENDIF
              INDX = 18
              DFFLW = DIFMN(DFFLW,DFFLP,DXGF(NW),DXGF(N),UL(1,NPX),INDX)
!
!---  Hydraulic dispersion  ---
!
              IF( IDSPS.EQ.1 ) THEN
                CALL ADVW( PORD,SL,UL,VL,WL,ULX,VLX,WLX,N,M )
                ULX = ULX*ULX
                VLX = VLX*VLX
                WLX = WLX*WLX
                ZVW = SQRT(ULX+VLX+WLX)
                INDX = 17
                DPLW = DIFMN(DPLGS(IZ(NW)),DPLGS(IZ(N)),
     &            DXGF(NW),DXGF(N),UL(1,NPX),INDX)
                DPTW = DIFMN(DPTRS(IZ(NW)),DPTRS(IZ(N)),
     &            DXGF(NW),DXGF(N),UL(1,NPX),INDX)
                DPLW = (DPLW*ULX + DPTW*(VLX+WLX))/(ZVW+SMALL)
              ELSE
                DPLW = 0.D+0
              ENDIF
!
!---  Salt aqueous flux by advection, diffusion, and dispersion  ---
!
              DDLW = (DFFLW+DPLW)/DXGP(NPX)
              IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  ---
!
                UDS(M,NPX) = DDLW*(XLS(MN,NW)*RHOL(MN,NW) -
     &            XLS(MP,N)*RHOL(MP,N))
                IF( UL(1,NPX).GE.ZERO ) THEN
                  IF( I.GT.2 ) THEN
                    NWW = NW-1
                    R = ((XLS(1,NW)*RHOL(1,NW)-XLS(1,NWW)*RHOL(1,NWW))
     &             /(XLS(1,N)*RHOL(1,N)-XLS(1,NW)*RHOL(1,NW)+SMALL))
     &                *((DXGF(N)+DXGF(NW))/(DXGF(NW)+DXGF(NWW)))
                  ELSE
                    R = 0.D+0
                  ENDIF
                  THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
                  DXF = DXGF(NW)/(DXGF(N)+DXGF(NW))
                  US(M,NPX) = XLS(1,N)*RHOL(1,N)*UL(1,NPX)*THETA*DXF
     &              + XLS(1,NW)*RHOL(1,NW)*UL(1,NPX)*(1.D+0-THETA*DXF)
                ELSEIF( UL(1,NPX).LT.ZERO ) THEN
                  IF( I.LT.IFLD ) THEN
                    NE = N+1
                    R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NE)*RHOL(1,NE))
     &             /(XLS(1,NW)*RHOL(1,NW)-XLS(1,N)*RHOL(1,N)+SMALL))
     &                *((DXGF(NW)+DXGF(N))/(DXGF(N)+DXGF(NE)))
                  ELSE
                    R = 0.D+0
                  ENDIF
                  THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
                  DXF = DXGF(N)/(DXGF(N)+DXGF(NW))
                  US(M,NPX) = XLS(1,NW)*RHOL(1,NW)*UL(1,NPX)*THETA*DXF
     &              + XLS(1,N)*RHOL(1,N)*UL(1,NPX)*(1.D+0-THETA*DXF)
                ENDIF
                US(M,NPX) = US(M,NPX) + UDS(M,NPX)
              ELSE
!
!---  Patankar salt transport  ---
!
                AL = MAX( UL(M,NPX),ZERO ) +
     &            DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/
     &            (DDLW+SMALL)))**5,ZERO)
                ALP = MAX( -UL(M,NPX),ZERO ) +
     &            DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/
     &            (DDLW+SMALL)))**5,ZERO)
                US(M,NPX) = XLS(MN,NW)*RHOL(MN,NW)*AL -
     &            XLS(MP,N)*RHOL(MP,N)*ALP
                UDS(M,NPX) = DDLW*(XLS(MN,NW)*RHOL(MN,NW) -
     &            XLS(MP,N)*RHOL(MP,N))
              ENDIF
  100       CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
      DO 400 I = 1,IFLD
        DO 400 K = 1,KFLD
          DO 400 J = 2,JFLD
            N = ND(I,J,K)
            NS = N-IFLD
            IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &        INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
            NPY = NSY(N)
            DO 300 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
              IF( IEDLS.EQ.1 ) THEN
                TCOR = (T(MP,N)+TABS)/TSPRF
                SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
                DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
                TCOR = (T(MN,NS)+TABS)/TSPRF
                SMDLP = DFLS(MN,NS)*TCOR*(VISRL/VISL(MN,NS))
                DFFLS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*SMDLP
              ELSEIF( IEDLS.EQ.2 ) THEN
                DFFLP = SDCLS(1,IZ(N))*SDCLS(2,IZ(N))*
     &            EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
                DFFLS = SDCLS(1,IZ(NS))*SDCLS(2,IZ(NS))*
     &            EXP(SL(MN,NS)*PORD(MN,NS)*SDCLS(3,IZ(NS)))
              ELSEIF( IEDLS.EQ.3 ) THEN
                DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
                DFFLS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*DFLS(MN,NS)
              ENDIF
              INDX = 18
              DFFLS = DIFMN(DFFLS,DFFLP,DYGF(NS),DYGF(N),VL(1,NPY),INDX)
!
!---  Hydraulic dispersion  ---
!
              IF( IDSPS.EQ.1 ) THEN
                CALL ADVS( PORD,SL,UL,VL,WL,ULSX,VLSX,WLSX,N,M )
                ULSX = ULSX*ULSX
                VLSX = VLSX*VLSX
                WLSX = WLSX*WLSX
                ZVS = SQRT(ULSX+VLSX+WLSX)
                INDX = 17
                DPLS = DIFMN(DPLGS(IZ(NS)),DPLGS(IZ(N)),
     &            DYGF(NS),DYGF(N),VL(1,NPY),INDX)
                DPTS = DIFMN(DPTRS(IZ(NS)),DPTRS(IZ(N)),
     &            DYGF(NS),DYGF(N),VL(1,NPY),INDX)
                DPLS = (DPLS*VLSX + DPTS*(WLSX+ULSX))/(ZVS+SMALL)
              ELSE
                DPLS = 0.D+0
              ENDIF
!
!---  Salt aqueous flux by advection, diffusion, and dispersion  ---
!
              DDLS = (DFFLS+DPLS)/DYGP(NPY)/RP(I)
              IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  ---
!
                VDS(M,NPY) = DDLS*(XLS(MN,NS)*RHOL(MN,NS) -
     &            XLS(MP,N)*RHOL(MP,N))
                IF( VL(1,NPY).GE.ZERO ) THEN
                  IF( J.GT.2 ) THEN
                    NSS = NS-IFLD
                    R = ((XLS(1,NS)*RHOL(1,NS)-XLS(1,NSS)*RHOL(1,NSS))
     &             /(XLS(1,N)*RHOL(1,N)-XLS(1,NS)*RHOL(1,NS)+SMALL))
     &                *((DYGF(N)+DYGF(NS))/(DYGF(NS)+DYGF(NSS)))
                  ELSE
                    R = 0.D+0
                  ENDIF
                  THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
                  DYF = DYGF(NS)/(DYGF(N)+DYGF(NS))
                  VS(M,NPY) = XLS(1,N)*RHOL(1,N)*VL(1,NPY)*THETA*DYF
     &              + XLS(1,NS)*RHOL(1,NS)*VL(1,NPY)*(1.D+0-THETA*DYF)
                ELSEIF( VL(1,NPY).LT.ZERO ) THEN
                  IF( J.LT.JFLD ) THEN
                    NN = N+IFLD
                    R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NN)*RHOL(1,NN))
     &             /(XLS(1,NS)*RHOL(1,NS)-XLS(1,N)*RHOL(1,N)+SMALL))
     &                *((DYGF(NS)+DYGF(N))/(DYGF(N)+DYGF(NN)))
                  ELSE
                    R = 0.D+0
                  ENDIF
                  THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
                  DYF = DYGF(N)/(DYGF(N)+DYGF(NS))
                  VS(M,NPY) = XLS(1,NS)*RHOL(1,NS)*VL(1,NPY)*THETA*DYF
     &              + XLS(1,N)*RHOL(1,N)*VL(1,NPY)*(1.D+0-THETA*DYF)
                ENDIF
                VS(M,NPY) = VS(M,NPY) + VDS(M,NPY)
              ELSE
!
!---  Patankar salt transport  ---
!
                AL = MAX( VL(M,NPY),ZERO ) +
     &            DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/
     &            (DDLS+SMALL)))**5,ZERO)
                ALP = MAX( -VL(M,NPY),ZERO ) +
     &            DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/
     &            (DDLS+SMALL)))**5,ZERO)
                VS(M,NPY) = (XLS(MN,NS)*RHOL(MN,NS)*AL -
     &            XLS(MP,N)*RHOL(MP,N)*ALP)
                VDS(M,NPY) = DDLS*(XLS(MN,NS)*RHOL(MN,NS) -
     &            XLS(MP,N)*RHOL(MP,N))
              ENDIF
  300       CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction Darcy velocities, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
      DO 600 J = 1,JFLD
        DO 600 I = 1,IFLD
          DO 600 K = 2,KFLD
            N = ND(I,J,K)
            NB = N-IJFLD
            IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &        INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
            NPZ = NSZ(N)
            DO 500 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
              IF( IEDLS.EQ.1 ) THEN
                TCOR = (T(MP,N)+TABS)/TSPRF
                SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
                DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
                TCOR = (T(MN,NB)+TABS)/TSPRF
                SMDLP = DFLS(MN,NB)*TCOR*(VISRL/VISL(MN,NB))
                DFFLB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*SMDLP
              ELSEIF( IEDLS.EQ.2 ) THEN
                DFFLP = SDCLS(1,IZ(N))*SDCLS(2,IZ(N))*
     &            EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
                DFFLB = SDCLS(1,IZ(NB))*SDCLS(2,IZ(NB))*
     &            EXP(SL(MN,NB)*PORD(MN,NB)*SDCLS(3,IZ(NB)))
              ELSEIF( IEDLS.EQ.3 ) THEN
                DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
                DFFLB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*DFLS(MN,NB)
              ENDIF
              INDX = 18
              DFFLB = DIFMN(DFFLB,DFFLP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
!
!---  Hydraulic dispersion  ---
!
              IF( IDSPS.EQ.1 ) THEN
                CALL ADVB( PORD,SL,UL,VL,WL,ULBX,VLBX,WLBX,N,M )
                ULBX = ULBX*ULBX
                VLBX = VLBX*VLBX
                WLBX = WLBX*WLBX
                ZVB = SQRT(ULBX+VLBX+WLBX)
                INDX = 17
                DPLB = DIFMN(DPLGS(IZ(NB)),DPLGS(IZ(N)),
     &            DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
                DPTB = DIFMN(DPTRS(IZ(NB)),DPTRS(IZ(N)),
     &            DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
                DPLB = (DPLB*WLBX + DPTB*(ULBX+VLBX))/(ZVB+SMALL)
              ELSE
                DPLB = 0.D+0
              ENDIF
!
!---  Salt aqueous flux by advection, diffusion, and dispersion  ---
!
              DDLB = (DFFLB+DPLB)/DZGP(NPZ)
              IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  ---
!
                WDS(M,NPZ) = DDLB*(XLS(MN,NB)*RHOL(MN,NB) -
     &            XLS(MP,N)*RHOL(MP,N))
                IF( WL(1,NPZ).GE.ZERO ) THEN
                  IF( K.GT. 2 ) THEN
                    NBB = NB-IJFLD
                    R = ((XLS(1,NB)*RHOL(1,NB)-XLS(1,NBB)*RHOL(1,NBB))
     &             /(XLS(1,N)*RHOL(1,N)-XLS(1,NB)*RHOL(1,NB)+SMALL))
     &                *((DZGF(N)+DZGF(NB))/(DZGF(NB)+DZGF(NBB)))
                  ELSE
                    R = 0.D+0
                  ENDIF
                  THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
                  DZF = DZGF(NB)/(DZGF(N)+DZGF(NB))
                  WS(M,NPZ) = XLS(1,N)*RHOL(1,N)*WL(1,NPZ)*THETA*DZF
     &              + XLS(1,NB)*RHOL(1,NB)*WL(1,NPZ)*(1.D+0-THETA*DZF)
                ELSEIF( WL(1,NPZ).LT.ZERO ) THEN
                  IF( K.LT.KFLD ) THEN
                    NT = N+IJFLD
                    R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NT)*RHOL(1,NT))
     &             /(XLS(1,NB)*RHOL(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &                *((DZGF(NB)+DZGF(N))/(DZGF(N)+DZGF(NT)))
                  ELSE
                    R = 0.D+0
                  ENDIF
                  THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &              (1.D+0+2.D+0*R)/3.D+0))
                  DZF = DZGF(N)/(DZGF(N)+DZGF(NB))
                  WS(M,NPZ) = XLS(1,NB)*RHOL(1,NB)*WL(1,NPZ)*THETA*DZF
     &              + XLS(1,N)*RHOL(1,N)*WL(1,NPZ)*(1.D+0-THETA*DZF)
                 ENDIF
                 WS(M,NPZ) = WS(M,NPZ) + WDS(M,NPZ)
              ELSE
!
!---  Patankar salt transport  ---
!
                AL = MAX( WL(M,NPZ),ZERO ) +
     &            DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/
     &            (DDLB+SMALL)))**5,ZERO)
                ALP = MAX( -WL(M,NPZ),ZERO ) +
     &            DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/
     &            (DDLB+SMALL)))**5,ZERO)
                WS(M,NPZ) = (XLS(MN,NB)*RHOL(MN,NB)*AL -
     &            XLS(MP,N)*RHOL(MP,N)*ALP)
                WDS(M,NPZ) = DDLB*(XLS(MN,NB)*RHOL(MN,NB) -
     &            XLS(MP,N)*RHOL(MP,N))
              ENDIF
  500       CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNFL group  ---
!
      RETURN
      END


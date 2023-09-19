!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGA
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
!     Compute diffusive CO2 gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGA'
!
!---  X-direction vapor mole diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( IXP(N).LE.0 .OR. IXP(N-1).LE.0 ) GOTO 200
          NPX = NSX(N)
          DXGA = XMGA(2,NW) - XMGA(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGA(MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*(SG(MN,NW)-SGT(MN,NW))*
     &        DFGA(MN,NW)*RHOMG(MN,NW)
            INDX = 20
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXGA,INDX )
            UDGA(M,NPX) = DFM*(XMGA(MN,NW)
     &        - XMGA(MP,N))/DXGP(NPX)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAW = XGA(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGA = DIFMN( FGAW,FGAP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGA(M,NPX) = UG(M,NPX)*FGA + WTMA*UDGA(M,NPX)
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( IXP(N).LE.0 .OR. IXP(N-IFLD).LE.0 ) GOTO 400
          NPY = NSY(N)
          DXGA = XMGA(2,NS) - XMGA(2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGA(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*(SG(MN,NS)-SGT(MN,NS))*
     &        DFGA(MN,NS)*RHOMG(MN,NS)
            INDX = 20
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXGA,INDX )
            VDGA(M,NPY) = DFM*(XMGA(MN,NS)
     &        - XMGA(MP,N))/(DYGP(NPY)*RP(I))
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAS = XGA(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGA = DIFMN( FGAS,FGAP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGA(M,NPY) = VG(M,NPY)*FGA + WTMA*VDGA(M,NPY)
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 600
          NB = ICM(1,1,N)
          IF( IXP(N).LE.0 .OR. IXP(N-IJFLD).LE.0 ) GOTO 600
          NPZ = NSZ(N)
          DXGA = XMGA(2,NB) - XMGA(2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGA(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*(SG(MN,NB)-SGT(MN,NB))*
     &        DFGA(MN,NB)*RHOMG(MN,NB)
            INDX = 20
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXGA,INDX )
            WDGA(M,NPZ) = DFM*(XMGA(MN,NB)
     &        - XMGA(MP,N))/DZGP(NPZ)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAB = XGA(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGA = DIFMN( FGAB,FGAP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGA(M,NPZ) = WG(M,NPZ)*FGA + WTMA*WDGA(M,NPZ)
  500     CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGA group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLA
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
!     Compute dissolved air molar diffusion rates through the
!     aqueous phase.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLA'
!
!---  X-direction molar diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( IXP(N).LE.0 .OR. IXP(N-1).LE.0 ) GOTO 200
          NPX = NSX(N)
          DXLA = XMLA(2,NW) - XMLA(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*
     &        DFLA(MN,NW)*RHOML(MN,NW)
            INDX = 14
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLA,INDX)
            UDLA(M,NPX) = DFM*(XMLA(MN,NW)-XMLA(MP,N))/DXGP(NPX)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAW = XLA(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLA = DIFMN( FLAW,FLAP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULA(M,NPX) = UL(M,NPX)*FLA + WTMA*UDLA(M,NPX)
  100     CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( IXP(N).LE.0 .OR. IXP(N-IFLD).LE.0 ) GOTO 400
          NPY = NSY(N)
          DXLA = XMLA(2,NS) - XMLA(2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*
     &        DFLA(MN,NS)*RHOML(MN,NS)
            INDX = 14
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLA,INDX )
            VDLA(M,NPY) = DFM*(XMLA(MN,NS)-XMLA(MP,N))/(DYGP(NPY)*RP(I))
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAS = XLA(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLA = DIFMN( FLAS,FLAP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLA(M,NPY) = VL(M,NPY)*FLA + WTMA*VDLA(M,NPY)
  300   CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 600
          NB = ICM(1,1,N)
          IF( IXP(N).LE.0 .OR. IXP(N-IJFLD).LE.0 ) GOTO 600
          NPZ = NSZ(N)
          DXLA = XMLA(2,NB) - XMLA(2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*
     &        DFLA(MN,NB)*RHOML(MN,NB)
            INDX = 14
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLA,INDX)
            WDLA(M,NPZ) = DFM*(XMLA(MN,NB)-XMLA(MP,N))/DZGP(NPZ)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAB = XLA(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLA = DIFMN( FLAB,FLAP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLA(M,NPZ) = WL(M,NPZ)*FLA + WTMA*WDLA(M,NPZ)
  500   CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLA group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLS
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
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
      SUB_LOG(ISUB_LOG) = '/DFFLS'
!
!---  X-direction Darcy velocities, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
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
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
              TCOR = (T(MN,NW)+TABS)/TSPRF
              SMDLP = DFLS(MN,NW)*TCOR*(VISRL/VISL(MN,NW))
              DFCLW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*SMDLP
            ELSEIF( IEDLS.EQ.2 ) THEN
              DFCLP = SDCLS(1,IZ(N))*SDCLS(2,IZ(N))*
     &          EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
              DFCLW = SDCLS(1,IZ(NW))*SDCLS(2,IZ(NW))*
     &          EXP(SL(MN,NW)*PORD(MN,NW)*SDCLS(3,IZ(NW)))
            ELSEIF( IEDLS.EQ.3 ) THEN
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
              DFCLW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*DFLS(MN,NW)
            ENDIF
            INDX = 18
            DFCLW = DIFMN(DFCLW,DFCLP,DXGF(NW),DXGF(N),UL(1,NPX),INDX)
!
!---  Hydraulic dispersion  ---
!
            IF( IDSPS.EQ.1 ) THEN
              CALL ADVW( PORD,SL,UL,VL,WL,ULWX,VLWX,WLWX,N,M )
              ULWX = ULWX*ULWX
              VLWX = VLWX*VLWX
              WLWX = WLWX*WLWX
              ZVW = SQRT(ULWX+VLWX+WLWX)
              INDX = 17
              DPLW = DIFMN(DPLGS(IZ(NW)),DPLGS(IZ(N)),DXGF(NW),DXGF(N),
     &          UL(1,NPX),INDX)
              DPTW = DIFMN(DPTRS(IZ(NW)),DPTRS(IZ(N)),DXGF(NW),DXGF(N),
     &          UL(1,NPX),INDX)
              DPLW = (DPLW*ULWX + DPTW*(VLWX+WLWX))/(ZVW+SMALL)
            ELSE
              DPLW = 0.D+0
            ENDIF
!
!---  Salt aqueous flux by advection, diffusion, and dispersion  ---
!
            DDLW = (DFCLW+DPLW)/DXGP(NPX)
            AL = MAX( UL(M,NPX),ZERO ) +
     &        DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/
     &        (DDLW+SMALL)))**5,ZERO)
            ALP = MAX( -UL(M,NPX),ZERO ) +
     &        DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/
     &        (DDLW+SMALL)))**5,ZERO)
            US(M,NPX) = XLS(MN,NW)*RHOL(MN,NW)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP
            UDS(M,NPX) = DDLW*(XLS(MN,NW)*RHOL(MN,NW) -
     &        XLS(MP,N)*RHOL(MP,N))
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
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
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
              TCOR = (T(MN,NS)+TABS)/TSPRF
              SMDLP = DFLS(MN,NS)*TCOR*(VISRL/VISL(MN,NS))
              DFCLS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*SMDLP
            ELSEIF( IEDLS.EQ.2 ) THEN
              DFCLP = SDCLS(1,IZ(N))*SDCLS(2,IZ(N))*
     &          EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
              DFCLS = SDCLS(1,IZ(NS))*SDCLS(2,IZ(NS))*
     &          EXP(SL(MN,NS)*PORD(MN,NS)*SDCLS(3,IZ(NS)))
            ELSEIF( IEDLS.EQ.3 ) THEN
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
              DFCLS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*DFLS(MN,NS)
            ENDIF
            INDX = 18
            DFCLS = DIFMN(DFCLS,DFCLP,DYGF(NS),DYGF(N),VL(1,NPY),INDX)
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
              DPLS = DIFMN(DPLGS(IZ(NS)),DPLGS(IZ(N)),DYGF(NS),DYGF(N),
     &          VL(1,NPY),INDX)
              DPTS = DIFMN(DPTRS(IZ(NS)),DPTRS(IZ(N)),DYGF(NS),DYGF(N),
     &          VL(1,NPY),INDX)
              DPLS = (DPLS*VLSX + DPTS*(WLSX+ULSX))/(ZVS+SMALL)
            ELSE
              DPLS = 0.D+0
            ENDIF
!
!---  Salt aqueous flux by advection, diffusion, and dispersion  ---
!
            DDLS = (DFCLS+DPLS)/(DYGP(NPY)*RP(I))
            AL = MAX( VL(M,NPY),ZERO ) +
     &        DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/
     &        (DDLS+SMALL)))**5,ZERO)
            ALP = MAX( -VL(M,NPY),ZERO ) +
     &        DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/
     &        (DDLS+SMALL)))**5,ZERO)
            VS(M,NPY) = (XLS(MN,NS)*RHOL(MN,NS)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP)
            VDS(M,NPY) = DDLS*(XLS(MN,NS)*RHOL(MN,NS) -
     &        XLS(MP,N)*RHOL(MP,N))
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction Darcy velocities, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 600
          NB = ICM(1,1,N)
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
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
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
              TCOR = (T(MN,NB)+TABS)/TSPRF
              SMDLP = DFLS(MN,NB)*TCOR*(VISRL/VISL(MN,NB))
              DFCLB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*SMDLP
            ELSEIF( IEDLS.EQ.2 ) THEN
              DFCLP = SDCLS(1,IZ(N))*SDCLS(2,IZ(N))*
     &          EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
              DFCLB = SDCLS(1,IZ(NB))*SDCLS(2,IZ(NB))*
     &          EXP(SL(MN,NB)*PORD(MN,NB)*SDCLS(3,IZ(NB)))
            ELSEIF( IEDLS.EQ.3 ) THEN
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
              DFCLB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*DFLS(MN,NB)
            ENDIF
            INDX = 18
            DFCLB = DIFMN(DFCLB,DFCLP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
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
              DPLB = DIFMN(DPLGS(IZ(NB)),DPLGS(IZ(N)),DZGF(NB),DZGF(N),
     &          WL(1,NPZ),INDX)
              DPTB = DIFMN(DPTRS(IZ(NB)),DPTRS(IZ(N)),DZGF(NB),DZGF(N),
     &          WL(1,NPZ),INDX)
              DPLB = (DPLB*WLBX + DPTB*(ULBX+VLBX))/(ZVB+SMALL)
            ELSE
              DPLB = 0.D+0
            ENDIF
!
!---  Salt aqueous flux by advection, diffusion, and dispersion  ---
!
            DDLB = (DFCLB+DPLB)/DZGP(NPZ)
            AL = MAX( WL(M,NPZ),ZERO ) +
     &        DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/
     &        (DDLB+SMALL)))**5,ZERO)
            ALP = MAX( -WL(M,NPZ),ZERO ) +
     &        DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/
     &        (DDLB+SMALL)))**5,ZERO)
            WS(M,NPZ) = (XLS(MN,NB)*RHOL(MN,NB)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP)
            WDS(M,NPZ) = DDLB*(XLS(MN,NB)*RHOL(MN,NB) -
     &        XLS(MP,N)*RHOL(MP,N))
  500     CONTINUE
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAB( N,NB )
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
!     Compute CO2 mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMGA = XMGAB(2,NB) - XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.26 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGA(M,NPZ) = WG(M,NPZ)*FGA
          ELSE
            WGA(M,NPZ) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQA,NB).EQ.42 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGA(M,NPZ) = WG(M,NPZ)*FGA
          ELSE
            WGA(M,NPZ) = 0.D+0
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGA(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGA,INDX )
          WDGA(M,NPZ) = DFM*(XMGAB(MP,NB)
     &      - XMGA(MP,N))/(5.D-1*DZGF(N))
          WGA(M,NPZ) = WG(M,NPZ)*FGA + WTMA*WDGA(M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAE( N,NB )
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
!     Compute VOC vapor mole diffusion rates on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMGA = XMGA(2,N) - XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.26 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGA(M,NQX) = UG(M,NQX)*FGA
          ELSE
            UGA(M,NQX) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQA,NB).EQ.42 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGA(M,NQX) = UG(M,NQX)*FGA
          ELSE
            UGA(M,NQX) = 0.D+0
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGA(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGA,INDX )
          UDGA(M,NQX) = DFM*(XMGA(MN,N)
     &      - XMGAB(MN,NB))/(5.D-1*DXGF(N))
          UGA(M,NQX) = UG(M,NQX)*FGA + WTMA*UDGA(M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAN( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGA = XMGA(2,N) - XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.26 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGA(M,NQY) = VG(M,NQY)*FGA
          ELSE
            VGA(M,NQY) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQA,NB).EQ.42 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGA(M,NQY) = VG(M,NQY)*FGA
          ELSE
            VGA(M,NQY) = 0.D+0
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGA(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGA,INDX )
          VDGA(M,NQY) = DFM*(XMGA(MN,N)
     &      - XMGAB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
          VGA(M,NQY) = VG(M,NQY)*FGA + WTMA*VDGA(M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAS( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGA = XMGAB(2,NB) - XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.26 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGA(M,NPY) = VG(M,NPY)*FGA
          ELSE
            VGA(M,NPY) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQA,NB).EQ.42 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGA(M,NPY) = VG(M,NPY)*FGA
          ELSE
            VGA(M,NPY) = 0.D+0
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGA(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGA,INDX )
          VDGA(M,NPY) = DFM*(XMGAB(MP,NB)
     &      - XMGA(MP,N))/((5.D-1*DYGF(N))*RP(I))
          VGA(M,NPY) = VG(M,NPY)*FGA + WTMA*VDGA(M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAT( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGA = XMGA(2,N)-XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.26 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGA(M,NQZ) = WG(M,NQZ)*FGA
          ELSE
            WGA(M,NQZ) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQA,NB).EQ.42 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGA(M,NQZ) = WG(M,NQZ)*FGA
          ELSE
            WGA(M,NQZ) = 0.D+0
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGA(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGA,INDX )
          WDGA(M,NQZ) = DFM*(XMGA(MN,N)
     &      - XMGAB(MN,NB))/(5.D-1*DZGF(N))
          WGA(M,NQZ) = WG(M,NQZ)*FGA + WTMA*WDGA(M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAW( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAW'
      I = ID(N)
      NPX = NSX(N)
      DXMGA = XMGAB(2,NB)-XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.26 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGA(M,NPX) = UG(M,NPX)*FGA
          ELSE
            UGA(M,NPX) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQA,NB).EQ.42 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGA(M,NPX) = UG(M,NPX)*FGA
          ELSE
            UGA(M,NPX) = 0.D+0
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGA(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGA,INDX )
          UDGA(M,NPX) = DFM*(XMGAB(MP,NB)
     &      - XMGA(MP,N))/(5.D-1*DXGF(N))
          UGA(M,NPX) = UG(M,NPX)*FGA + WTMA*UDGA(M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAB( N,NB )
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
!     Compute CO2 mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVA
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMLA = XMLAB(2,NB) - XMLA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*(1.D+0-ZLAB(MP,NB)) +
     &      XNAB(MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*(1.D+0+ZLAB(MP,NB)) +
     &      XGAB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLA(M,NPZ) = WL(M,NPZ)*FLA
          ELSE
            WLA(M,NPZ) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLA(M,NPZ) = WL(M,NPZ)*FLA
          ELSE
            WLA(M,NPZ) = 0.D+0
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)*
     &      (1.D+0-ABS(ZLAB(MP,NB)))*RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,NPZ) = DFM*(XMLAB(MP,NB)
     &      - XMLA(MP,N))/(5.D-1*DZGF(N))
          WLA(M,NPZ) = WL(M,NPZ)*FLA + WTMA*WDLA(M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAE( N,NB )
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
!     Compute VOC vapor mole diffusion rates on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVA
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMLA = XMLA(2,N) - XMLAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*(1.D+0-ZLAB(MN,NB)) +
     &      XNAB(MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*(1.D+0+ZLAB(MN,NB)) +
     &      XGAB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULA(M,NQX) = UL(M,NQX)*FLA
          ELSE
            ULA(M,NQX) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULA(M,NQX) = UL(M,NQX)*FLA
          ELSE
            ULA(M,NQX) = 0.D+0
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)*
     &      (1.D+0-ABS(ZLAB(MN,NB)))*RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,NQX) = DFM*(XMLA(MN,N)
     &      - XMLAB(MN,NB))/(5.D-1*DXGF(N))
          ULA(M,NQX) = UL(M,NQX)*FLA + WTMA*UDLA(M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAN( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVA
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMLA = XMLA(2,N) - XMLAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*(1.D+0-ZLAB(MN,NB)) +
     &      XNAB(MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*(1.D+0+ZLAB(MN,NB)) +
     &      XGAB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLA(M,NQY) = VL(M,NQY)*FLA
          ELSE
            VLA(M,NQY) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLA(M,NQY) = VL(M,NQY)*FLA
          ELSE
            VLA(M,NQY) = 0.D+0
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)*
     &      (1.D+0-ABS(ZLAB(MN,NB)))*RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,NQY) = DFM*(XMLA(MN,N)
     &      - XMLAB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
          VLA(M,NQY) = VL(M,NQY)*FLA + WTMA*VDLA(M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAS( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVA
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMLA = XMLAB(2,NB) - XMLA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*(1.D+0-ZLAB(MP,NB)) +
     &      XNAB(MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*(1.D+0+ZLAB(MP,NB)) +
     &      XGAB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLA(M,NPY) = VL(M,NPY)*FLA
          ELSE
            VLA(M,NPY) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLA(M,NPY) = VL(M,NPY)*FLA
          ELSE
            VLA(M,NPY) = 0.D+0
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)*
     &      (1.D+0-ABS(ZLAB(MP,NB)))*RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,NPY) = DFM*(XMLAB(MP,NB)
     &      - XMLA(MP,N))/((5.D-1*DYGF(N))*RP(I))
          VLA(M,NPY) = VL(M,NPY)*FLA + WTMA*VDLA(M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAT( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVA
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMLA = XMLA(2,N)-XMLAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*(1.D+0-ZLAB(MN,NB)) +
     &      XNAB(MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*(1.D+0+ZLAB(MN,NB)) +
     &      XGAB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLA(M,NQZ) = WL(M,NQZ)*FLA
          ELSE
            WLA(M,NQZ) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLA(M,NQZ) = WL(M,NQZ)*FLA
          ELSE
            WLA(M,NQZ) = 0.D+0
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)*
     &      (1.D+0-ABS(ZLAB(MN,NB)))*RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,NQZ) = DFM*(XMLA(MN,N)
     &      - XMLAB(MN,NB))/(5.D-1*DZGF(N))
          WLA(M,NQZ) = WL(M,NQZ)*FLA + WTMA*WDLA(M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAW( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVA
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAW'
      I = ID(N)
      NPX = NSX(N)
      DXMLA = XMLAB(2,NB)-XMLA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*(1.D+0-ZLAB(MP,NB)) +
     &      XNAB(MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*(1.D+0+ZLAB(MP,NB)) +
     &      XGAB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULA(M,NPX) = UL(M,NPX)*FLA
          ELSE
            ULA(M,NPX) = 0.D+0
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULA(M,NPX) = UL(M,NPX)*FLA
          ELSE
            ULA(M,NPX) = 0.D+0
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)*
     &      (1.D+0-ABS(ZLAB(MP,NB)))*RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,NPX) = DFM*(XMLAB(MP,NB)
     &      - XMLA(MP,N))/(5.D-1*DXGF(N))
          ULA(M,NPX) = UL(M,NPX)*FLA + WTMA*UDLA(M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSB( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFLSB'
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
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MP,N)*
     &      PORD(MP,N)*SDCLS(3,IZN))
          DFFLB = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MP,NB)*
     &      PORDB(MP,NB)*SDCLS(3,IZN))
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          ALP = MAX( -WL(M,NPZ),ZERO )
          WS(M,NPZ) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          AL = MAX( WL(M,NPZ),ZERO )
          WS(M,NPZ) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLB = (DFFLB+DPLB)/(5.D-1*DZGF(N))
          AL = MAX( WL(M,NPZ),ZERO ) +
     &     DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/(DDLB+SMALL)))**5,ZERO)
          ALP = MAX( -WL(M,NPZ),ZERO ) +
     &     DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/(DDLB+SMALL)))**5,ZERO)
          WS(M,NPZ) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL
     &     - XLS(MP,N)*RHOL(MP,N)*ALP)
          WDS(M,NPZ) = DDLB*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &      XLS(MP,N)*RHOL(MP,N))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSS( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFLSS'
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
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MP,N)*
     &      PORD(MP,N)*SDCLS(3,IZN))
          DFFLS = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MP,NB)*
     &      PORDB(MP,NB)*SDCLS(3,IZN))
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          ALP = MAX( -VL(M,NPY),ZERO )
          VS(M,NPY) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          AL = MAX( VL(M,NPY),ZERO )
          VS(M,NPY) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLS = (DFFLS+DPLS)/RP(I)/(5.D-1*DYGF(N))
          AL = MAX( VL(M,NPY),ZERO ) +
     &     DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/(DDLS+SMALL)))**5,ZERO)
          ALP = MAX( -VL(M,NPY),ZERO ) +
     &     DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/(DDLS+SMALL)))**5,ZERO)
          VS(M,NPY) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL
     &      - XLS(MP,N)*RHOL(MP,N)*ALP)
          VDS(M,NPY) = DDLS*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &      XLS(MP,N)*RHOL(MP,N))
        ENDIF

  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSW( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFLSW'
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
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MP,N)*
     &      PORD(MP,N)*SDCLS(3,IZN))
          DFFLW = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MP,NB)*
     &      PORDB(MP,NB)*SDCLS(3,IZN))
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
          ULWX = UWX*UWX
          VLWX = VWX*VWX
          WLWX = WWX*WWX
          ZLW = SQRT(ULWX+VLWX+WLWX)
          DPLW = (DPLGS(IZN)*ULWX + DPTRS(IZN)*(WLWX+VLWX))/(ZLW+SMALL)
        ELSE
          DPLW = 0.D+0
        ENDIF
!
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          ALP = MAX( -UL(M,NPX),ZERO )
          US(M,NPX) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          AL = MAX( UL(M,NPX),ZERO )
          US(M,NPX) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLW = (DFFLW+DPLW)/(5.D-1*DXGF(N))
          AL = MAX( UL(M,NPX),ZERO ) +
     &     DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/(DDLW+SMALL)))**5,ZERO)
          ALP = MAX( -UL(M,NPX),ZERO ) +
     &     DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/(DDLW+SMALL)))**5,ZERO)
          US(M,NPX) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL
     &      - XLS(MP,N)*RHOL(MP,N)*ALP)
          UDS(M,NPX) = DDLW*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &      XLS(MP,N)*RHOL(MP,N))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSE( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFLSE'
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
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MN,N)*
     &      PORD(MN,N)*SDCLS(3,IZN))
          DFFLE = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MN,NB)*
     &      PORDB(MN,NB)*SDCLS(3,IZN))
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          ALP = MAX( UL(M,NQX),ZERO )
          US(M,NQX) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          AL = MAX( -UL(M,NQX),ZERO )
          US(M,NQX) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLE = (DFFLE+DPLE)/(5.D-1*DXGF(N))
          AL = MAX( -UL(M,NQX),ZERO ) +
     &     DDLE*MAX((ONE-(TENTH*ABS(UL(M,NQX))/(DDLE+SMALL)))**5,ZERO)
          ALP = MAX( UL(M,NQX),ZERO ) +
     &     DDLE*MAX((ONE-(TENTH*ABS(UL(M,NQX))/(DDLE+SMALL)))**5,ZERO)
          US(M,NQX) = (XLS(MN,N)*RHOL(MN,N)*ALP
     &     - XLSB(MN,NB)*RHOLB(MN,NB)*AL)
          UDS(M,NQX) = DDLE*(XLS(MN,N)*RHOL(MN,N) -
     &     XLSB(MN,NB)*RHOLB(MN,NB))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSN( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFLSN'
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
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MN,N)*
     &      PORD(MN,N)*SDCLS(3,IZN))
          DFFLN = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MN,NB)*
     &      PORDB(MN,NB)*SDCLS(3,IZN))
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          ALP = MAX( VL(M,NQY),ZERO )
          VS(M,NQY) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          AL = MAX( -VL(M,NQY),ZERO )
          VS(M,NQY) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLN = (DFFLN+DPLN)/RP(I)/(5.D-1*DYGF(N))
          AL = MAX( -VL(M,NQY),ZERO ) +
     &     DDLN*MAX((ONE-(TENTH*ABS(VL(M,NQY))/(DDLN+SMALL)))**5,ZERO)
          ALP = MAX( VL(M,NQY),ZERO ) +
     &     DDLN*MAX((ONE-(TENTH*ABS(VL(M,NQY))/(DDLN+SMALL)))**5,ZERO)
          VS(M,NQY) = (XLS(MN,N)*RHOL(MN,N)*ALP -
     &      XLSB(MN,NB)*RHOLB(MN,NB)*AL)
          VDS(M,NQY) = DDLN*(XLS(MN,N)*RHOL(MN,N) -
     &      XLSB(MN,NB)*RHOLB(MN,NB))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLST( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFLST'
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
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MN,N)*
     &      PORD(MN,N)*SDCLS(3,IZN))
          DFFLT = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MN,NB)*
     &      PORDB(MN,NB)*SDCLS(3,IZN))
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQW,NB).EQ.26 ) THEN
          ALP = MAX( WL(M,NQZ),ZERO )
          WS(M,NQZ) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.42 ) THEN
          AL = MAX( -WL(M,NQZ),ZERO )
          WS(M,NQZ) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLT = (DFFLT+DPLT)/(5.D-1*DZGF(N))
          AL = MAX( -WL(M,NQZ),ZERO ) +
     &     DDLT*MAX((ONE-(TENTH*ABS(WL(M,NQZ))/(DDLT+SMALL)))**5,ZERO)
          ALP = MAX( WL(M,NQZ),ZERO ) +
     &     DDLT*MAX((ONE-(TENTH*ABS(WL(M,NQZ))/(DDLT+SMALL)))**5,ZERO)
          WS(M,NQZ) = (XLS(MN,N)*RHOL(MN,N)*ALP
     &      - XLSB(MN,NB)*RHOLB(MN,NB)*AL)
          WDS(M,NQZ) = DDLT*(XLS(MN,N)*RHOL(MN,N) -
     &      XLSB(MN,NB)*RHOLB(MN,NB))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLST group  ---
!
      RETURN
      END

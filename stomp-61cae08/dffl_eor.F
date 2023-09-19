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
!     Diffusive CO2 aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
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
      DO 200 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO 150 NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
          DXMLA = XMLA(2,NW) - XMLA(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*
     &        DFLA(MN,NW)*RHOML(MN,NW)
            INDX = 14
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMLA,INDX)
            UDLA(M,NPX) = DFM*(XMLA(MN,NW)-XMLA(MP,N))/DXGP(NPX)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAW = XLA(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLA = DIFMN( FLAW,FLAP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULA(M,NPX) = UL(M,NPX)*FLA + WTMA*UDLA(M,NPX)
  100     CONTINUE
  150   CONTINUE
  200 CONTINUE
!
!---  Y-direction molar diffusion, excluding boundaries
!
      DO 400 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 400
        DO 350 NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
          DXMLA = XMLA(2,NS) - XMLA(2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*
     &        DFLA(MN,NS)*RHOML(MN,NS)
            INDX = 14
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMLA,INDX )
            VDLA(M,NPY) = DFM*(XMLA(MN,NS)-XMLA(MP,N))/(DYGP(NPY)*RP(I))
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAS = XLA(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLA = DIFMN( FLAS,FLAP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLA(M,NPY) = VL(M,NPY)*FLA + WTMA*VDLA(M,NPY)
  300     CONTINUE
  350   CONTINUE
  400 CONTINUE
!
!---  Z-direction molar diffusion, excluding boundaries
!
      DO 600 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 600
        DO 550 NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
          DXMLA = XMLA(2,NB) - XMLA(2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*
     &        DFLA(MN,NB)*RHOML(MN,NB)
            INDX = 14
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMLA,INDX)
            WDLA(M,NPZ) = DFM*(XMLA(MN,NB)-XMLA(MP,N))/DZGP(NPZ)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAB = XLA(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLA = DIFMN( FLAB,FLAP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLA(M,NPZ) = WL(M,NPZ)*FLA + WTMA*WDLA(M,NPZ)
  500     CONTINUE
  550   CONTINUE
  600 CONTINUE
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
!     Diffusive salt aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
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
      DO 200 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO 150 NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Diffusion coefficients  ---
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
!---        Salt aqueous flux by advection and diffusion  ---
!
            DDLW = DFCLW/DXGP(NPX)
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
  150   CONTINUE
  200 CONTINUE
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      DO 400 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 400
        DO 350 NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
          DO 300 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Diffusion coefficients  ---
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
!---        Salt aqueous flux by advection and diffusion  ---
!
            DDLS = DFCLS/(DYGP(NPY)*RP(I))
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
  350   CONTINUE
  400 CONTINUE
!
!---  Z-direction Darcy velocities, excluding boundaries
!
      DO 600 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 600
        DO 550 NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
          DO 500 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Diffusion coefficients  ---
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
!---        Salt aqueous flux by advection and diffusion  ---
!
            DDLB = DFCLB/DZGP(NPZ)
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
  550   CONTINUE
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLW
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
!     Diffusive water aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
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
      SUB_LOG(ISUB_LOG) = '/DFFLW'
!
!---  X-direction molar diffusion, excluding boundaries
!
      DO 200 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO 150 NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            UDLW(M,NPX) = -(UDLA(M,NPX) + UDS(M,NPX)/WTMS)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWW = XLW(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULW(M,NPX) = UL(M,NPX)*FLW + WTMW*UDLW(M,NPX)
  100     CONTINUE
  150   CONTINUE
  200 CONTINUE
!
!---  Y-direction molar diffusion, excluding boundaries
!
      DO 400 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 400
        DO 350 NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            VDLW(M,NPY) = -(VDLA(M,NPY) + VDS(M,NPY)/WTMS)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWS = XLW(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLW = DIFMN( FLWS,FLWP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLW(M,NPY) = VL(M,NPY)*FLW + WTMW*VDLW(M,NPY)
  300     CONTINUE
  350   CONTINUE
  400 CONTINUE
!
!---  Z-direction molar diffusion, excluding boundaries
!
      DO 600 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 600
        DO 550 NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            WDLW(M,NPZ) = -(WDLA(M,NPZ) + WDS(M,NPZ)/WTMS)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWB = XLW(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLW(M,NPZ) = WL(M,NPZ)*FLW + WTMW*WDLW(M,NPZ)
  500     CONTINUE
  550   CONTINUE
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLW group  ---
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
!     Diffusive CO2 aqeuous fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAB'
      ITFX = MOD(IBCT(2,NB),10)
      K = KD(N)
      NPZ = NSZ(N)
      DXMLA = XMLAB(2,NB) - XMLA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLAP = XLA(MP,N)*RHOL(MP,N)
        FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLA(M,NPZ) = WL(M,NPZ)*FLA
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLA(M,NPZ) = WL(M,NPZ)*FLA
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,NPZ) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/DZGF(N)
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
!     Diffusive CO2 aqueous flux on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAE'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      NQX = NSX(N)+1
      DXMLA = XMLA(2,N) - XMLAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLAP = XLA(MN,N)*RHOL(MN,N)
        FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULA(M,NQX) = UL(M,NQX)*FLA
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULA(M,NQX) = UL(M,NQX)*FLA
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,NQX) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))/DXGF(N)
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
!     Diffusive CO2 aqeuous fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAN'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMLA = XMLA(2,N) - XMLAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLAP = XLA(MN,N)*RHOL(MN,N)
        FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLA(M,NQY) = VL(M,NQY)*FLA
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLA(M,NQY) = VL(M,NQY)*FLA
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,NQY) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))
     &      /(DYGF(N)*RP(I))
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
!     Diffusive CO2 aqeuous fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAS'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMLA = XMLAB(2,NB) - XMLA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLAP = XLA(MP,N)*RHOL(MP,N)
        FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLA(M,NPY) = VL(M,NPY)*FLA
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLA(M,NPY) = VL(M,NPY)*FLA
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,NPY) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))
     &      /(DYGF(N)*RP(I))
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
!     Diffusive CO2 aqeuous fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAT'
      ITFX = MOD(IBCT(2,NB),10)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMLA = XMLA(2,N)-XMLAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLAP = XLA(MN,N)*RHOL(MN,N)
        FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLA(M,NQZ) = WL(M,NQZ)*FLA
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLA(M,NQZ) = WL(M,NQZ)*FLA
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,NQZ) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))/DZGF(N)
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
!     Diffusive CO2 aqeuous fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAW'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      NPX = NSX(N)
      DXMLA = XMLAB(2,NB)-XMLA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLAP = XLA(MP,N)*RHOL(MP,N)
        FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULA(M,NPX) = UL(M,NPX)*FLA
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULA(M,NPX) = UL(M,NPX)*FLA
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,NPX) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/DXGF(N)
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
!     Diffusive salt aqueous fluxes on bottom boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( -WL(M,NPZ),ZERO )
          WS(M,NPZ) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( WL(M,NPZ),ZERO )
          WS(M,NPZ) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLB = DFFLB/(5.D-1*DZGF(N))
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
!     Diffusive salt aqueous fluxes on south boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( -VL(M,NPY),ZERO )
          VS(M,NPY) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( VL(M,NPY),ZERO )
          VS(M,NPY) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLS = DFFLS/RP(I)/(5.D-1*DYGF(N))
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
!     Diffusive salt aqueous fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( -UL(M,NPX),ZERO )
          US(M,NPX) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( UL(M,NPX),ZERO )
          US(M,NPX) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLW = DFFLW/(5.D-1*DXGF(N))
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
!     Diffusive salt aqueous fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( UL(M,NQX),ZERO )
          US(M,NQX) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( -UL(M,NQX),ZERO )
          US(M,NQX) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLE = DFFLE/(5.D-1*DXGF(N))
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
!     Diffusive salt aqueous fluxes on north boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( VL(M,NQY),ZERO )
          VS(M,NQY) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( -VL(M,NQY),ZERO )
          VS(M,NQY) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLN = DFFLN/RP(I)/(5.D-1*DYGF(N))
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
!     Diffusive salt aqueous fluxes on top boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
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
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( WL(M,NQZ),ZERO )
          WS(M,NQZ) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( -WL(M,NQZ),ZERO )
          WS(M,NQZ) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLT = DFFLT/(5.D-1*DZGF(N))
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

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWB( N,NB )
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
!     Diffusive water aqueous flux on bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWB'
      ITFX = MOD(IBCT(2,NB),10)
      K = KD(N)
      NPZ = NSZ(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLW(M,NPZ) = WL(M,NPZ)*FLW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLW(M,NPZ) = WL(M,NPZ)*FLW
          ENDIF
        ELSE
          WDLW(M,NPZ) = -(WDLA(M,NPZ) + WDS(M,NPZ)/WTMS)
          WLW(M,NPZ) = WL(M,NPZ)*FLW + WTMW*WDLW(M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWE( N,NB )
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
!     Diffusive water aqueous flux on east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWE'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      NQX = NSX(N)+1
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULW(M,NQX) = UL(M,NQX)*FLW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULW(M,NQX) = UL(M,NQX)*FLW
          ENDIF
        ELSE
          UDLW(M,NQX) = -(UDLA(M,NQX) + UDS(M,NQX)/WTMS)
          ULW(M,NQX) = UL(M,NQX)*FLW + WTMW*UDLW(M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWN( N,NB )
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
!     Diffusive water aqueous flux on north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWN'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLW(M,NQY) = VL(M,NQY)*FLW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLW(M,NQY) = VL(M,NQY)*FLW
          ENDIF
        ELSE
          VDLW(M,NQY) = -(VDLA(M,NQY) + VDS(M,NQY)/WTMS)
          VLW(M,NQY) = VL(M,NQY)*FLW + WTMW*VDLW(M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWS( N,NB )
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
!     Diffusive water aqueous flux on south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWS'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLW(M,NPY) = VL(M,NPY)*FLW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLW(M,NPY) = VL(M,NPY)*FLW
          ENDIF
        ELSE
          VDLW(M,NPY) = -(VDLA(M,NPY) + VDS(M,NPY)/WTMS)
          VLW(M,NPY) = VL(M,NPY)*FLW + WTMW*VDLW(M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWT( N,NB )
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
!     Diffusive water aqueous flux on top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWT'
      ITFX = MOD(IBCT(2,NB),10)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLW(M,NQZ) = WL(M,NQZ)*FLW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLW(M,NQZ) = WL(M,NQZ)*FLW
          ENDIF
        ELSE
          WDLW(M,NQZ) = -(WDLA(M,NQZ) + WDS(M,NQZ)/WTMS)
          WLW(M,NQZ) = WL(M,NQZ)*FLW + WTMW*WDLW(M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWW( N,NB )
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
!     Diffusive water aqueous flux on west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWW'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      NPX = NSX(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULW(M,NPX) = UL(M,NPX)*FLW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULW(M,NPX) = UL(M,NPX)*FLW
          ENDIF
        ELSE
          UDLW(M,NPX) = -(UDLA(M,NPX) + UDS(M,NPX)/WTMS)
          ULW(M,NPX) = UL(M,NPX)*FLW + WTMW*UDLW(M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWW group  ---
!
      RETURN
      END


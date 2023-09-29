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
!     Written by MD White, PNNL, 13 June 2012
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
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DXMLA = XMLA(2,NW) - XMLA(2,N)
          DO M = 1,ISVF
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
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DXMLA = XMLA(2,NS) - XMLA(2,N)
          DO M = 1,ISVF
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
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DXMLA = XMLA(2,NB) - XMLA(2,N)
          DO M = 1,ISVF
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
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLA group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLN
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
!     Diffusive N2 aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
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
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DFFLN'
!
!---  X-direction molar diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DXMLN = (XMLN(2,NW)-XMLN(2,N))
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLN(MP,N)
     &        *RHOML(MP,N)
            DFSW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*DFLN(MN,NW)
     &        *RHOML(MN,NW)
            INDX = 14
            DFM = DIFMN( DFSW,DFSP,DXGF(NW),DXGF(N),DXMLN,INDX )
            UDLN(M,NPX) = DFM*(XMLN(MN,NW)-XMLN(MP,N))/DXGP(NPX)
            FLNP = XLN(MP,N)*RHOL(MP,N)
            FLNW = XLN(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLN = DIFMN( FLNW,FLNP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULN(M,NPX) = UL(M,NPX)*FLN + WTMN*UDLN(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DXMLN = (XMLN(2,NS)-XMLN(2,N))
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLN(MP,N)
     &        *RHOML(MP,N)
            DFSS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*DFLN(MN,NS)
     &        *RHOML(MN,NS)
            INDX = 14
            DFM = DIFMN( DFSS,DFSP,DYGF(NS),DYGF(N),DXMLN,INDX )
            VDLN(M,NPY) = DFM*(XMLN(MN,NS)-XMLN(MP,N))/(DYGP(NPY)*RP(I))
            FLNP = XLN(MP,N)*RHOL(MP,N)
            FLNS = XLN(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLN = DIFMN( FLNS,FLNP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLN(M,NPY) = VL(M,NPY)*FLN + WTMN*VDLN(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DXMLN = (XMLN(2,NB)-XMLN(2,N))
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLN(MP,N)
     &        *RHOML(MP,N)
            DFSB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*DFLN(MN,NB)
     &        *RHOML(MN,NB)
            INDX = 14
            DFM = DIFMN( DFSB,DFSP,DZGF(NB),DZGF(N),DXMLN,INDX )
            WDLN(M,NPZ) = DFM*(XMLN(MN,NB)-XMLN(MP,N))/DZGP(NPZ)
            FLNP = XLN(MP,N)*RHOL(MP,N)
            FLNB = XLN(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLN = DIFMN( FLNB,FLNP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLN(M,NPZ) = WL(M,NPZ)*FLN + WTMN*WDLN(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLO
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
!     Diffusive CH4 aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLO'
!
!---  X-direction molar diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DXMLO = (XMLO(2,NW)-XMLO(2,N))
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &        *RHOML(MP,N)
            DFSW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*DFLO(MN,NW)
     &        *RHOML(MN,NW)
            INDX = 14
            DFM = DIFMN( DFSW,DFSP,DXGF(NW),DXGF(N),DXMLO,INDX )
            UDLO(M,NPX) = DFM*(XMLO(MN,NW)-XMLO(MP,N))/DXGP(NPX)
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOW = XLO(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLO = DIFMN( FLOW,FLOP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DXMLO = (XMLO(2,NS)-XMLO(2,N))
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &        *RHOML(MP,N)
            DFSS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*DFLO(MN,NS)
     &        *RHOML(MN,NS)
            INDX = 14
            DFM = DIFMN( DFSS,DFSP,DYGF(NS),DYGF(N),DXMLO,INDX )
            VDLO(M,NPY) = DFM*(XMLO(MN,NS)-XMLO(MP,N))/(DYGP(NPY)*RP(I))
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOS = XLO(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLO = DIFMN( FLOS,FLOP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DXMLO = (XMLO(2,NB)-XMLO(2,N))
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &        *RHOML(MP,N)
            DFSB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*DFLO(MN,NB)
     &        *RHOML(MN,NB)
            INDX = 14
            DFM = DIFMN( DFSB,DFSP,DZGF(NB),DZGF(N),DXMLO,INDX )
            WDLO(M,NPZ) = DFM*(XMLO(MN,NB)-XMLO(MP,N))/DZGP(NPZ)
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOB = XLO(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLO = DIFMN( FLOB,FLOP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLO group  ---
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
!     Written by MD White, PNNL, 13 June 2012
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
      SUB_LOG(ISUB_LOG) = '/DFFLS'
!
!---  X-direction Darcy velocities, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DO M = 1,ISVF
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
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DO M = 1,ISVF
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
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction Darcy velocities, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DO M = 1,ISVF
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
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by MD White, PNNL, 13 June 2012
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
      IF( IFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            UDLW(M,NPX) = -(UDLA(M,NPX) + UDLO(M,NPX) + UDLN(M,NPX)
     &        + UDS(M,NPX)/WTMS)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWW = XLW(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULW(M,NPX) = UL(M,NPX)*FLW + WTMW*UDLW(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            VDLW(M,NPY) = -(VDLA(M,NPY) + VDLO(M,NPY) + VDLN(M,NPY)
     &        + VDS(M,NPY)/WTMS)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWS = XLW(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLW = DIFMN( FLWS,FLWP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLW(M,NPY) = VL(M,NPY)*FLW + WTMW*VDLW(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            WDLW(M,NPZ) = -(WDLA(M,NPZ) + WDLO(M,NPZ) + WDLN(M,NPZ)
     &        + WDS(M,NPZ)/WTMS)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWB = XLW(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLW(M,NPZ) = WL(M,NPZ)*FLW + WTMW*WDLW(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by MD White, PNNL, 13 June 2012
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
      K = KD(N)
      NPZ = NSZ(N)
      DXMLA = XMLAB(2,NB) - XMLA(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLAP = XLA(MP,N)*RHOL(MP,N)
        FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLA(M,NPZ) = WL(M,NPZ)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
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
      ENDDO
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
!     Written by MD White, PNNL, 13 June 2012
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
      I = ID(N)
      NQX = NSX(N)+1
      DXMLA = XMLA(2,N) - XMLAB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FLAP = XLA(MN,N)*RHOL(MN,N)
        FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULA(M,NQX) = UL(M,NQX)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
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
      ENDDO
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
!     Written by MD White, PNNL, 13 June 2012
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
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMLA = XMLA(2,N) - XMLAB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FLAP = XLA(MN,N)*RHOL(MN,N)
        FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLA(M,NQY) = VL(M,NQY)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
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
      ENDDO
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
!     Written by MD White, PNNL, 13 June 2012
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
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMLA = XMLAB(2,NB) - XMLA(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLAP = XLA(MP,N)*RHOL(MP,N)
        FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLA(M,NPY) = VL(M,NPY)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
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
      ENDDO
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
!     Written by MD White, PNNL, 13 June 2012
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
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMLA = XMLA(2,N)-XMLAB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FLAP = XLA(MN,N)*RHOL(MN,N)
        FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLA(M,NQZ) = WL(M,NQZ)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
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
      ENDDO
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
!     Written by MD White, PNNL, 13 June 2012
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
      I = ID(N)
      NPX = NSX(N)
      DXMLA = XMLAB(2,NB)-XMLA(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLAP = XLA(MP,N)*RHOL(MP,N)
        FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULA(M,NPX) = UL(M,NPX)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
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
      ENDDO
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
      SUBROUTINE DFFLNB( N,NB )
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
!     Diffusive N2 aqeuous fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
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
      USE FDVH
      USE CONST
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLNB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMLN = XMLNB(2,NB) - XMLN(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLNP = XLN(MP,N)*RHOL(MP,N)
        FLNB = XLNB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLN = DIFMN( FLNB,FLNP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLN(M,NPZ) = WL(M,NPZ)*FLN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLN(M,NPZ) = WL(M,NPZ)*FLN
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLN(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLNB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLN,INDX )
          WDLN(M,NPZ) = 2.D+0*DFM*(XMLNB(MP,NB)-XMLN(MP,N))/DZGF(N)
          WLN(M,NPZ) = WL(M,NPZ)*FLN + WTMN*WDLN(M,NPZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLNB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLNE( N,NB )
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
!     Diffusive N2 aqeuous fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
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
      USE FDVH
      USE CONST
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLNE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMLN = XMLN(2,N) - XMLNB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FLNP = XLN(MN,N)*RHOL(MN,N)
        FLNB = XLNB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLN = DIFMN( FLNP,FLNB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULN(M,NQX) = UL(M,NQX)*FLN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULN(M,NQX) = UL(M,NQX)*FLN
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLN(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLNB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLN,INDX )
          UDLN(M,NQX) = 2.D+0*DFM*(XMLN(MN,N)-XMLNB(MN,NB))/DXGF(N)
          ULN(M,NQX) = UL(M,NQX)*FLN + WTMN*UDLN(M,NQX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLNE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLNN( N,NB )
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
!     Diffusive N2 aqeuous fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
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
      USE FDVH
      USE CONST
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLNN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMLN = XMLN(2,N) - XMLNB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FLNP = XLN(MN,N)*RHOL(MN,N)
        FLNB = XLNB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLN = DIFMN( FLNP,FLNB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLN(M,NQY) = VL(M,NQY)*FLN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLN(M,NQY) = VL(M,NQY)*FLN
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLN(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLNB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLN,INDX )
          VDLN(M,NQY) = 2.D+0*DFM*(XMLN(MN,N)-XMLNB(MN,NB))
     &      /(DYGF(N)*RP(I))
          VLN(M,NQY) = VL(M,NQY)*FLN + WTMN*VDLN(M,NQY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLNN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLNS( N,NB )
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
!     Diffusive CH4 aqueous fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
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
      USE FDVH
      USE CONST
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLNS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMLN = XMLNB(2,NB) - XMLN(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLNP = XLN(MP,N)*RHOL(MP,N)
        FLNB = XLNB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLN = DIFMN( FLNB,FLNP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLN(M,NPY) = VL(M,NPY)*FLN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLN(M,NPY) = VL(M,NPY)*FLN
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLN(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLNB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLN,INDX )
          VDLN(M,NPY) = 2.D+0*DFM*(XMLNB(MP,NB)-XMLN(MP,N))
     &      /(DYGF(N)*RP(I))
          VLN(M,NPY) = VL(M,NPY)*FLN + WTMN*VDLN(M,NPY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLNS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLNT( N,NB )
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
!     Diffusive N2 aqueous fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
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
      USE FDVH
      USE CONST
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLNT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMLN = XMLN(2,N)-XMLNB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FLNP = XLN(MN,N)*RHOL(MN,N)
        FLNB = XLNB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLN = DIFMN( FLNP,FLNB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLN(M,NQZ) = WL(M,NQZ)*FLN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLN(M,NQZ) = WL(M,NQZ)*FLN
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLN(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLNB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLN,INDX )
          WDLN(M,NQZ) = 2.D+0*DFM*(XMLN(MN,N)-XMLNB(MN,NB))/DZGF(N)
          WLN(M,NQZ) = WL(M,NQZ)*FLN + WTMN*WDLN(M,NQZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLNT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLNW( N,NB )
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
!     Diffusive N2 aqueous fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
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
      USE FDVH
      USE CONST
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLNW'
      I = ID(N)
      NPX = NSX(N)
      DXMLN = XMLNB(2,NB)-XMLN(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLNP = XLN(MP,N)*RHOL(MP,N)
        FLNB = XLNB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLN = DIFMN( FLNB,FLNP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULN(M,NPX) = UL(M,NPX)*FLN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULN(M,NPX) = UL(M,NPX)*FLN
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLN(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLNB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLN,INDX )
          UDLN(M,NPX) = 2.D+0*DFM*(XMLNB(MP,NB)-XMLN(MP,N))/DXGF(N)
          ULN(M,NPX) = UL(M,NPX)*FLN + WTMN*UDLN(M,NPX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLNW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOB( N,NB )
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
!     Diffusive CH4 aqeuous fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLOB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMLO = XMLOB(2,NB) - XMLO(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLOP = XLO(MP,N)*RHOL(MP,N)
        FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLO = DIFMN( FLOB,FLOP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLO(M,NPZ) = WL(M,NPZ)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLO(M,NPZ) = WL(M,NPZ)*FLO
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLO,INDX )
          WDLO(M,NPZ) = 2.D+0*DFM*(XMLOB(MP,NB)-XMLO(MP,N))/DZGF(N)
          WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOE( N,NB )
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
!     Diffusive CH4 aqeuous fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLOE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMLO = XMLO(2,N) - XMLOB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FLOP = XLO(MN,N)*RHOL(MN,N)
        FLOB = XLOB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLO = DIFMN( FLOP,FLOB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULO(M,NQX) = UL(M,NQX)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULO(M,NQX) = UL(M,NQX)*FLO
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLO,INDX )
          UDLO(M,NQX) = 2.D+0*DFM*(XMLO(MN,N)-XMLOB(MN,NB))/DXGF(N)
          ULO(M,NQX) = UL(M,NQX)*FLO + WTMO*UDLO(M,NQX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLON( N,NB )
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
!     Diffusive CH4 aqeuous fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLON'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMLO = XMLO(2,N) - XMLOB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FLOP = XLO(MN,N)*RHOL(MN,N)
        FLOB = XLOB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLO = DIFMN( FLOP,FLOB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLO(M,NQY) = VL(M,NQY)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLO(M,NQY) = VL(M,NQY)*FLO
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLO,INDX )
          VDLO(M,NQY) = 2.D+0*DFM*(XMLO(MN,N)-XMLOB(MN,NB))
     &      /(DYGF(N)*RP(I))
          VLO(M,NQY) = VL(M,NQY)*FLO + WTMO*VDLO(M,NQY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOS( N,NB )
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
!     Diffusive CH4 aqueous fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLOS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMLO = XMLOB(2,NB) - XMLO(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLOP = XLO(MP,N)*RHOL(MP,N)
        FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLO = DIFMN( FLOB,FLOP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLO(M,NPY) = VL(M,NPY)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLO(M,NPY) = VL(M,NPY)*FLO
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLO,INDX )
          VDLO(M,NPY) = 2.D+0*DFM*(XMLOB(MP,NB)-XMLO(MP,N))
     &      /(DYGF(N)*RP(I))
          VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOT( N,NB )
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
!     Diffusive CH4 aqueous fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLOT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMLO = XMLO(2,N)-XMLOB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FLOP = XLO(MN,N)*RHOL(MN,N)
        FLOB = XLOB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLO = DIFMN( FLOP,FLOB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLO(M,NQZ) = WL(M,NQZ)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLO(M,NQZ) = WL(M,NQZ)*FLO
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLO,INDX )
          WDLO(M,NQZ) = 2.D+0*DFM*(XMLO(MN,N)-XMLOB(MN,NB))/DZGF(N)
          WLO(M,NQZ) = WL(M,NQZ)*FLO + WTMO*WDLO(M,NQZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOW( N,NB )
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
!     Diffusive CH4 aqueous fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLOW'
      I = ID(N)
      NPX = NSX(N)
      DXMLO = XMLOB(2,NB)-XMLO(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLOP = XLO(MP,N)*RHOL(MP,N)
        FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLO = DIFMN( FLOB,FLOP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULO(M,NPX) = UL(M,NPX)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULO(M,NPX) = UL(M,NPX)*FLO
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLO,INDX )
          UDLO(M,NPX) = 2.D+0*DFM*(XMLOB(MP,NB)-XMLO(MP,N))/DXGF(N)
          ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOW group  ---
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
!     Written by MD White, PNNL, 13 June 2012
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
      DO M = 1,ISVF
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by MD White, PNNL, 13 June 2012
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
      DO M = 1,ISVF
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by MD White, PNNL, 13 June 2012
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
      DO M = 1,ISVF
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by MD White, PNNL, 13 June 2012
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
      DO M = 1,ISVF
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by MD White, PNNL, 13 June 2012
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
      DO M = 1,ISVF
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by MD White, PNNL, 13 June 2012
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
      DO M = 1,ISVF
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by MD White, PNNL, 13 June 2012
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
      K = KD(N)
      NPZ = NSZ(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLW(M,NPZ) = WL(M,NPZ)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLW(M,NPZ) = WL(M,NPZ)*FLW
          ENDIF
        ELSE
          WDLW(M,NPZ) = -(WDLA(M,NPZ) + WDLO(M,NPZ) + WDLN(M,NPZ)
     &      + WDS(M,NPZ)/WTMS)
          WLW(M,NPZ) = WL(M,NPZ)*FLW + WTMW*WDLW(M,NPZ)
        ENDIF
      ENDDO
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
!     Written by MD White, PNNL, 13 June 2012
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
      I = ID(N)
      NQX = NSX(N)+1
      DO M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULW(M,NQX) = UL(M,NQX)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULW(M,NQX) = UL(M,NQX)*FLW
          ENDIF
        ELSE
          UDLW(M,NQX) = -(UDLA(M,NQX) + UDLO(M,NQX) + UDLN(M,NQX)
     &      + UDS(M,NQX)/WTMS)
          ULW(M,NQX) = UL(M,NQX)*FLW + WTMW*UDLW(M,NQX)
        ENDIF
      ENDDO
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
!     Written by MD White, PNNL, 13 June 2012
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
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DO M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLW(M,NQY) = VL(M,NQY)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLW(M,NQY) = VL(M,NQY)*FLW
          ENDIF
        ELSE
          VDLW(M,NQY) = -(VDLA(M,NQY) + VDLO(M,NQY) + VDLN(M,NQY)
     &      + VDS(M,NQY)/WTMS)
          VLW(M,NQY) = VL(M,NQY)*FLW + WTMW*VDLW(M,NQY)
        ENDIF
      ENDDO
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
!     Written by MD White, PNNL, 13 June 2012
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
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLW(M,NPY) = VL(M,NPY)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLW(M,NPY) = VL(M,NPY)*FLW
          ENDIF
        ELSE
          VDLW(M,NPY) = -(VDLA(M,NPY) + VDLO(M,NPY) + VDLN(M,NPY)
     &      + VDS(M,NPY)/WTMS)
          VLW(M,NPY) = VL(M,NPY)*FLW + WTMW*VDLW(M,NPY)
        ENDIF
      ENDDO
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
!     Written by MD White, PNNL, 13 June 2012
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
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DO M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLW(M,NQZ) = WL(M,NQZ)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLW(M,NQZ) = WL(M,NQZ)*FLW
          ENDIF
        ELSE
          WDLW(M,NQZ) = -(WDLA(M,NQZ) + WDLO(M,NQZ) + WDLN(M,NQZ)
     &      + WDS(M,NQZ)/WTMS)
          WLW(M,NQZ) = WL(M,NQZ)*FLW + WTMW*WDLW(M,NQZ)
        ENDIF
      ENDDO
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
!     Written by MD White, PNNL, 13 June 2012
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
      I = ID(N)
      NPX = NSX(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULW(M,NPX) = UL(M,NPX)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULW(M,NPX) = UL(M,NPX)*FLW
          ENDIF
        ELSE
          UDLW(M,NPX) = -(UDLA(M,NPX) + UDLO(M,NPX) + UDLN(M,NPX)
     &      + UDS(M,NPX)/WTMS)
          ULW(M,NPX) = UL(M,NPX)*FLW + WTMW*UDLW(M,NPX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWW group  ---
!
      RETURN
      END


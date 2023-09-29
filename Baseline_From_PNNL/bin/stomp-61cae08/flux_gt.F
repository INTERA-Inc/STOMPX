!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGW_GT
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
!     Compute water vapor mole diffusion rates.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFGW_GT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
!
!---  X-direction vapor mole diffusion, excluding boundaries
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
          DXMGW = XMGW(2,NW)-XMGW(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*(SG(MN,NW)-SGT(MN,NW))*
     &        DFGW(MN,NW)*RHOMG(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMGW,INDX )
            UDGW(M,NPX) = DFM*(XMGW(MN,NW)-XMGW(MP,N))/DXGP(NPX)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWW = XGW(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
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
          DXMGW = XMGW(2,NS)-XMGW(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*(SG(MN,NS)-SGT(MN,NS))*
     &        DFGW(MN,NS)*RHOMG(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGW,INDX )
            VDGW(M,NPY) = DFM*(XMGW(MN,NS)-XMGW(MP,N))/DYGP(NPY)/RP(I)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWS = XGW(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGW = DIFMN( FGWS,FGWP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
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
          DXMGW = XMGW(2,NB)-XMGW(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*(SG(MN,NB)-SGT(MN,NB))*
     &        DFGW(MN,NB)*RHOMG(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGW,INDX )
            WDGW(M,NPZ) = DFM*(XMGW(MN,NB)-XMGW(MP,N))/DZGP(NPZ)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWB = XGW(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGW = DIFMN( FGWB,FGWP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
!
!---  X-direction vapor mole diffusion, excluding boundaries
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
          DXMGW = XMGW(2,NW)*RHOMG(2,NW)-XMGW(2,N)*RHOMG(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*(SG(MN,NW)-SGT(MN,NW))*
     &        DFGW(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMGW,INDX )
            UDGW(M,NPX) = DFM*(XMGW(MN,NW)*RHOMG(MN,NW)
     &        - XMGW(MP,N)*RHOMG(MP,N))/DXGP(NPX)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWW = XGW(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
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
          DXMGW = XMGW(2,NS)*RHOMG(2,NS)-XMGW(2,N)*RHOMG(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*(SG(MN,NS)-SGT(MN,NS))*
     &        DFGW(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGW,INDX )
            VDGW(M,NPY) = DFM*(XMGW(MN,NS)*RHOMG(MN,NS)
     &        - XMGW(MP,N)*RHOMG(MP,N))/DYGP(NPY)/RP(I)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWS = XGW(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGW = DIFMN( FGWS,FGWP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
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
          DXMGW = XMGW(2,NB)*RHOMG(2,NB)-XMGW(2,N)*RHOMG(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*(SG(MN,NB)-SGT(MN,NB))*
     &        DFGW(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGW,INDX )
            WDGW(M,NPZ) = DFM*(XMGW(MN,NB)*RHOMG(MN,NB)
     &        - XMGW(MP,N)*RHOMG(MP,N))/DZGP(NPZ)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWB = XGW(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGW = DIFMN( FGWB,FGWP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLA_GT
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
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFLA_GT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
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
          DXLA = (XMLA(2,NW)-XMLA(2,N))
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &        *RHOML(MP,N)
            DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*DFLA(MN,NW)
     &        *RHOML(MN,NW)
            INDX = 14
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLA,INDX)
            UDLA(M,NPX) = DFM*(XMLA(MN,NW)-XMLA(MP,N))/DXGP(NPX)
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
          DXLA = (XMLA(2,NS)-XMLA(2,N))
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &        *RHOML(MP,N)
            DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*DFLA(MN,NS)
     &        *RHOML(MN,NS)
            INDX = 14
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLA,INDX )
            VDLA(M,NPY) = DFM*(XMLA(MN,NS)-XMLA(MP,N))/DYGP(NPY)/RP(I)
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
          DXLA = (XMLA(2,NB)-XMLA(2,N))
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &        *RHOML(MP,N)
            DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*DFLA(MN,NB)
     &        *RHOML(MN,NB)
            INDX = 14
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLA,INDX)
            WDLA(M,NPZ) = DFM*(XMLA(MN,NB)-XMLA(MP,N))/DZGP(NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
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
          DXLA = (XMLA(2,NW)*RHOML(2,NW)-XMLA(2,N)*RHOML(2,N))
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
            DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*DFLA(MN,NW)
            INDX = 14
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLA,INDX)
            UDLA(M,NPX) = DFM*(XMLA(MN,NW)*RHOML(MN,NW)
     &       - XMLA(MP,N)*RHOML(MP,N))/DXGP(NPX)
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
          DXLA = (XMLA(2,NS)*RHOML(2,NS)-XMLA(2,N)*RHOML(2,N))
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
            DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*DFLA(MN,NS)
            INDX = 14
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLA,INDX )
            VDLA(M,NPY) = DFM*(XMLA(MN,NS)*RHOML(MN,NS)
     &        - XMLA(MP,N)*RHOML(MP,N))/DYGP(NPY)/RP(I)
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
          DXLA = (XMLA(2,NB)*RHOML(2,NB)-XMLA(2,N)*RHOML(2,N))
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
            DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*DFLA(MN,NB)
            INDX = 14
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLA,INDX)
            WDLA(M,NPZ) = DFM*(XMLA(MN,NB)*RHOML(MN,NB)
     &        - XMLA(MP,N)*RHOML(MP,N))/DZGP(NPZ)
          ENDDO
        ENDDO
      ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLA_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLS_GT
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
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFLS_GT'
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
!---  Diffusion coefficients  ---
!
            IF( IEDLS.EQ.1 ) THEN
              TCOR = (T(MP,N)+TABS)/TSPRF
              SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
              TCOR = (T(MN,NW)+TABS)/TSPRF
              SMDLP = DFLS(MN,NW)*TCOR*(VISRL/VISL(MN,NW))
              DFCLW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*SMDLP
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
              CALL ADVW( PORD,SL,UL,VL,WL,ULX,VLX,WLX,N,M )
              ULX = ULX*ULX
              VLX = VLX*VLX
              WLX = WLX*WLX
              ZVW = SQRT(ULX+VLX+WLX)
              INDX = 17
              DPLW = DIFMN(DPLGS(IZ(NW)),DPLGS(IZ(N)),
     &          DXGF(NW),DXGF(N),UL(1,NPX),INDX)
              DPTW = DIFMN(DPTRS(IZ(NW)),DPTRS(IZ(N)),
     &          DXGF(NW),DXGF(N),UL(1,NPX),INDX)
              DPLW = (DPLW*ULX + DPTW*(VLX+WLX))/(ZVW+SMALL)
            ELSE
              DPLW = 0.D+0
            ENDIF
!
!---  Salt aqueous flux by advection, diffusion, and dispersion  ---
!
            DDLW = (DFCLW+DPLW)/DXGP(NPX)
            IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  ---
!
              UDS(M,NPX) = DDLW*(XLS(MN,NW)*RHOL(MN,NW) -
     &          XLS(MP,N)*RHOL(MP,N))
              IF( UL(1,NPX).GE.ZERO ) THEN
                IF( I.GT.2 ) THEN
                  NWW = NW-1
                  R = ((XLS(1,NW)*RHOL(1,NW)-XLS(1,NWW)*RHOL(1,NWW))
     &           /(XLS(1,N)*RHOL(1,N)-XLS(1,NW)*RHOL(1,NW)+SMALL))
     &              *((DXGF(N)+DXGF(NW))/(DXGF(NW)+DXGF(NWW)))
                ELSE
                  R = 0.D+0
                ENDIF
                THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &            (1.D+0+2.D+0*R)/3.D+0))
                DXF = DXGF(NW)/(DXGF(N)+DXGF(NW))
                US(M,NPX) = XLS(1,N)*RHOL(1,N)*UL(1,NPX)*THETA*DXF
     &            + XLS(1,NW)*RHOL(1,NW)*UL(1,NPX)*(1.D+0-THETA*DXF)
              ELSEIF( UL(1,NPX).LT.ZERO ) THEN
                IF( I.LT.IFLD ) THEN
                  NE = N+1
                  R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NE)*RHOL(1,NE))
     &           /(XLS(1,NW)*RHOL(1,NW)-XLS(1,N)*RHOL(1,N)+SMALL))
     &              *((DXGF(NW)+DXGF(N))/(DXGF(N)+DXGF(NE)))
                ELSE
                  R = 0.D+0
                ENDIF
                THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &            (1.D+0+2.D+0*R)/3.D+0))
                DXF = DXGF(N)/(DXGF(N)+DXGF(NW))
                US(M,NPX) = XLS(1,NW)*RHOL(1,NW)*UL(1,NPX)*THETA*DXF
     &            + XLS(1,N)*RHOL(1,N)*UL(1,NPX)*(1.D+0-THETA*DXF)
              ENDIF
!              US(M,NPX) = US(M,NPX) + UDS(M,NPX)
            ELSE
!
!---  Patankar salt transport  ---
!
              AL = MAX( UL(M,NPX),ZERO ) +
     &          DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/
     &          (DDLW+SMALL)))**5,ZERO)
              ALP = MAX( -UL(M,NPX),ZERO ) +
     &          DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/
     &          (DDLW+SMALL)))**5,ZERO)
              US(M,NPX) = XLS(MN,NW)*RHOL(MN,NW)*AL -
     &          XLS(MP,N)*RHOL(MP,N)*ALP
              UDS(M,NPX) = DDLW*(XLS(MN,NW)*RHOL(MN,NW) -
     &          XLS(MP,N)*RHOL(MP,N))
            ENDIF
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
!---  Diffusion coefficients  ---
!
            IF( IEDLS.EQ.1 ) THEN
              TCOR = (T(MP,N)+TABS)/TSPRF
              SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
              TCOR = (T(MN,NS)+TABS)/TSPRF
              SMDLP = DFLS(MN,NS)*TCOR*(VISRL/VISL(MN,NS))
              DFCLS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*SMDLP
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
              DPLS = DIFMN(DPLGS(IZ(NS)),DPLGS(IZ(N)),
     &          DYGF(NS),DYGF(N),VL(1,NPY),INDX)
              DPTS = DIFMN(DPTRS(IZ(NS)),DPTRS(IZ(N)),
     &          DYGF(NS),DYGF(N),VL(1,NPY),INDX)
              DPLS = (DPLS*VLSX + DPTS*(WLSX+ULSX))/(ZVS+SMALL)
            ELSE
              DPLS = 0.D+0
            ENDIF
!
!---  Salt aqueous flux by advection, diffusion, and dispersion  ---
!
            DDLS = (DFCLS+DPLS)/DYGP(NPY)/RP(I)
            IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  ---
!
              VDS(M,NPY) = DDLS*(XLS(MN,NS)*RHOL(MN,NS) -
     &          XLS(MP,N)*RHOL(MP,N))
              IF( VL(1,NPY).GE.ZERO ) THEN
                IF( J.GT.2 ) THEN
                  NSS = NS-IFLD
                  R = ((XLS(1,NS)*RHOL(1,NS)-XLS(1,NSS)*RHOL(1,NSS))
     &           /(XLS(1,N)*RHOL(1,N)-XLS(1,NS)*RHOL(1,NS)+SMALL))
     &              *((DYGF(N)+DYGF(NS))/(DYGF(NS)+DYGF(NSS)))
                ELSE
                  R = 0.D+0
                ENDIF
                THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &            (1.D+0+2.D+0*R)/3.D+0))
                DYF = DYGF(NS)/(DYGF(N)+DYGF(NS))
                VS(M,NPY) = XLS(1,N)*RHOL(1,N)*VL(1,NPY)*THETA*DYF
     &            + XLS(1,NS)*RHOL(1,NS)*VL(1,NPY)*(1.D+0-THETA*DYF)
              ELSEIF( VL(1,NPY).LT.ZERO ) THEN
                IF( J.LT.JFLD ) THEN
                  NN = N+IFLD
                  R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NN)*RHOL(1,NN))
     &           /(XLS(1,NS)*RHOL(1,NS)-XLS(1,N)*RHOL(1,N)+SMALL))
     &              *((DYGF(NS)+DYGF(N))/(DYGF(N)+DYGF(NN)))
                ELSE
                  R = 0.D+0
                ENDIF
                THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &            (1.D+0+2.D+0*R)/3.D+0))
                DYF = DYGF(N)/(DYGF(N)+DYGF(NS))
                VS(M,NPY) = XLS(1,NS)*RHOL(1,NS)*VL(1,NPY)*THETA*DYF
     &            + XLS(1,N)*RHOL(1,N)*VL(1,NPY)*(1.D+0-THETA*DYF)
              ENDIF
!              VS(M,NPY) = VS(M,NPY) + VDS(M,NPY)
            ELSE
!
!---  Patankar salt transport  ---
!
              AL = MAX( VL(M,NPY),ZERO ) +
     &          DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/
     &          (DDLS+SMALL)))**5,ZERO)
              ALP = MAX( -VL(M,NPY),ZERO ) +
     &          DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/
     &          (DDLS+SMALL)))**5,ZERO)
              VS(M,NPY) = (XLS(MN,NS)*RHOL(MN,NS)*AL -
     &          XLS(MP,N)*RHOL(MP,N)*ALP)
              VDS(M,NPY) = DDLS*(XLS(MN,NS)*RHOL(MN,NS) -
     &          XLS(MP,N)*RHOL(MP,N))
            ENDIF
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
!---  Diffusion coefficients  ---
!
            IF( IEDLS.EQ.1 ) THEN
              TCOR = (T(MP,N)+TABS)/TSPRF
              SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
              TCOR = (T(MN,NB)+TABS)/TSPRF
              SMDLP = DFLS(MN,NB)*TCOR*(VISRL/VISL(MN,NB))
              DFCLB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*SMDLP
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
              DPLB = DIFMN(DPLGS(IZ(NB)),DPLGS(IZ(N)),
     &          DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
              DPTB = DIFMN(DPTRS(IZ(NB)),DPTRS(IZ(N)),
     &          DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
              DPLB = (DPLB*WLBX + DPTB*(ULBX+VLBX))/(ZVB+SMALL)
            ELSE
              DPLB = 0.D+0
            ENDIF
!
!---  Salt aqueous flux by advection, diffusion, and dispersion  ---
!
            DDLB = (DFCLB+DPLB)/DZGP(NPZ)
            IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  ---
!
              WDS(M,NPZ) = DDLB*(XLS(MN,NB)*RHOL(MN,NB) -
     &          XLS(MP,N)*RHOL(MP,N))
              IF( WL(1,NPZ).GE.ZERO ) THEN
                IF( K.GT. 2 ) THEN
                  NBB = NB-IJFLD
                  R = ((XLS(1,NB)*RHOL(1,NB)-XLS(1,NBB)*RHOL(1,NBB))
     &           /(XLS(1,N)*RHOL(1,N)-XLS(1,NB)*RHOL(1,NB)+SMALL))
     &              *((DZGF(N)+DZGF(NB))/(DZGF(NB)+DZGF(NBB)))
                ELSE
                  R = 0.D+0
                ENDIF
                THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &            (1.D+0+2.D+0*R)/3.D+0))
                DZF = DZGF(NB)/(DZGF(N)+DZGF(NB))
                WS(M,NPZ) = XLS(1,N)*RHOL(1,N)*WL(1,NPZ)*THETA*DZF
     &            + XLS(1,NB)*RHOL(1,NB)*WL(1,NPZ)*(1.D+0-THETA*DZF)
              ELSEIF( WL(1,NPZ).LT.ZERO ) THEN
                IF( K.LT.KFLD ) THEN
                  NT = N+IJFLD
                  R = ((XLS(1,N)*RHOL(1,N)-XLS(1,NT)*RHOL(1,NT))
     &           /(XLS(1,NB)*RHOL(1,NB)-XLS(1,N)*RHOL(1,N)+SMALL))
     &              *((DZGF(NB)+DZGF(N))/(DZGF(N)+DZGF(NT)))
                ELSE
                  R = 0.D+0
                ENDIF
                THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,
     &            (1.D+0+2.D+0*R)/3.D+0))
                DZF = DZGF(N)/(DZGF(N)+DZGF(NB))
                WS(M,NPZ) = XLS(1,NB)*RHOL(1,NB)*WL(1,NPZ)*THETA*DZF
     &            + XLS(1,N)*RHOL(1,N)*WL(1,NPZ)*(1.D+0-THETA*DZF)
               ENDIF
!               WS(M,NPZ) = WS(M,NPZ) + WDS(M,NPZ)
            ELSE
!
!---  Patankar salt transport  ---
!
              AL = MAX( WL(M,NPZ),ZERO ) +
     &          DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/
     &          (DDLB+SMALL)))**5,ZERO)
              ALP = MAX( -WL(M,NPZ),ZERO ) +
     &          DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/
     &          (DDLB+SMALL)))**5,ZERO)
              WS(M,NPZ) = (XLS(MN,NB)*RHOL(MN,NB)*AL -
     &          XLS(MP,N)*RHOL(MP,N)*ALP)
              WDS(M,NPZ) = DDLB*(XLS(MN,NB)*RHOL(MN,NB) -
     &          XLS(MP,N)*RHOL(MP,N))
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWB_GT( N,NB )
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
!     Compute water vapor mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWB_GT'
      NPZ = NSZ(N)
      IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        K = KD(N)
        DXMGW = XMGWB(2,NB)-XMGW(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &      *DFGW(MP,N)*RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)
     &      *DFGWB(MP,NB)*RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGW,INDX )
          WDGW(M,NPZ) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))
     &      /DZGF(N)
          FGWP = XGW(MP,N)*RHOG(MP,N)
          FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
          INDX = 3
          FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
            IF( WG(1,NPZ).LT.-EPSL ) THEN
              WGW(M,NPZ) = WG(M,NPZ)*FGW
            ELSE
              WGW(M,NPZ) = 0.D+0
            ENDIF
          ELSE
            WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
          ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        K = KD(N)
        DXMGW = XMGWB(2,NB)*RHOMGB(2,NB) - XMGW(2,N)*RHOMG(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGW(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGW,INDX )
          WDGW(M,NPZ) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB)
     &      - XMGW(MP,N)*RHOMG(MP,N))/DZGF(N)
          FGWP = XGW(MP,N)*RHOG(MP,N)
          FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
          INDX = 3
          FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
            IF( WG(1,NPZ).LT.-EPSL ) THEN
              WGW(M,NPZ) = WG(M,NPZ)*FGW
            ELSE
              WGW(M,NPZ) = 0.D+0
            ENDIF
          ELSE
            WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWE_GT( N,NB )
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
!     Compute water vapor mole diffusion rates on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWE_GT'
      NQX = NSX(N)+1
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        I = ID(N)
        DXMGW = XMGW(2,N)-XMGWB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &      *DFGW(MN,N)*RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*
     &      DFGWB(MN,NB)*RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGW,INDX )
          UDGW(M,NQX) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))
     &      /DXGF(N)
          FGWP = XGW(MN,N)*RHOG(MN,N)
          FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
          INDX = 3
          FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
            IF( UG(1,NQX).GT.EPSL ) THEN
              UGW(M,NQX) = UG(M,NQX)*FGW
            ELSE
              UGW(M,NQX) = 0.D+0
            ENDIF
          ELSE
            UGW(M,NQX) = UG(M,NQX)*FGW + WTMW*UDGW(M,NQX)
          ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        I = ID(N)
        DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGW(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGW,INDX )
          UDGW(M,NQX) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N)
     &      - XMGWB(MN,NB)*RHOMGB(MN,NB))/DXGF(N)
          FGWP = XGW(MN,N)*RHOG(MN,N)
          FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
          INDX = 3
          FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
            IF( UG(1,NQX).GT.EPSL ) THEN
              UGW(M,NQX) = UG(M,NQX)*FGW
            ELSE
              UGW(M,NQX) = 0.D+0
            ENDIF
          ELSE
            UGW(M,NQX) = UG(M,NQX)*FGW + WTMW*UDGW(M,NQX)
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWN_GT( N,NB )
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
!     Compute water vapor mole diffusion rates on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWN_GT'
      NQY = NSY(N)+IFLD
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        I = ID(N)
        J = JD(N)
        DXMGW = XMGW(2,N)-XMGWB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &      *DFGW(MN,N)*RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*
     &      DFGWB(MN,NB)*RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGW,INDX )
          VDGW(M,NQY) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))
     &      /(DYGF(N)*RP(I))
          FGWP = XGW(MN,N)*RHOG(MN,N)
          FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
          INDX = 3
          FGW = DIFMN( FGWP,FGWB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
            IF( VG(1,NQY).GT.EPSL ) THEN
              VGW(M,NQY) = VG(M,NQY)*FGW
            ELSE
              VGW(M,NQY) = 0.D+0
            ENDIF
          ELSE
            VGW(M,NQY) = VG(M,NQY)*FGW
          ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        I = ID(N)
        J = JD(N)
        DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGW(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGW,INDX )
          VDGW(M,NQY) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N)
     &      - XMGWB(MN,NB)*RHOMGB(MN,NB))/(DYGF(N)*RP(I))
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWS_GT( N,NB )
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
!     Compute water vapor mole diffusion rates on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWS_GT'
      NPY = NSY(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        I = ID(N)
        J = JD(N)
        DXMGW = XMGWB(2,NB)-XMGW(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &      *DFGW(MP,N)*RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*
     &      DFGWB(MP,NB)*RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGW,INDX )
          VDGW(M,NPY) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))/
     &      (DYGF(N)*RP(I))
          FGWP = XGW(MP,N)*RHOG(MP,N)
          FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
          INDX = 3
          FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
            IF( VG(1,NPY).LT.-EPSL ) THEN
              VGW(M,NPY) = VG(M,NPY)*FGW
            ELSE
              VGW(M,NPY) = 0.D+0
            ENDIF
          ELSE
            VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
          ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        I = ID(N)
        J = JD(N)
        DXMGW = XMGWB(2,NB)*RHOMGB(2,NB) - XMGW(2,N)*RHOMG(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGW(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGW,INDX )
          VDGW(M,NPY) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB)
     &      - XMGW(MP,N)*RHOMG(MP,N))/(DYGF(N)*RP(I))
          FGWP = XGW(MP,N)*RHOG(MP,N)
          FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
          INDX = 3
          FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
            IF( VG(1,NPY).LT.-EPSL ) THEN
              VGW(M,NPY) = VG(M,NPY)*FGW
            ELSE
              VGW(M,NPY) = 0.D+0
            ENDIF
          ELSE
            VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWT_GT( N,NB )
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
!     Compute water vapor mole diffusion rates on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWT_GT'
      NQZ = NSZ(N)+IJFLD
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        K = KD(N)
        DXMGW = XMGW(2,N)-XMGWB(2,NB)
        DO M = 1,ISVF
           MN = MNEG(M)
           DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &       *DFGW(MN,N)*RHOMG(MN,N)
           DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*
     &       DFGWB(MN,NB)*RHOMGB(MN,NB)
           INDX = 12
           DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGW,INDX )
           WDGW(M,NQZ) = 2.D+0*DFM*(XMGW(MN,N)-
     &       XMGWB(MN,NB))/DZGF(N)
           FGWP = XGW(MN,N)*RHOG(MN,N)
           FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
           INDX = 3
           FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---       Dirichlet-Outflow boundary condition  ---
!
           IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
             IF( WG(1,NQZ).GT.EPSL ) THEN
               WGW(M,NQZ) = WG(M,NQZ)*FGW
             ELSE
               WGW(M,NQZ) = 0.D+0
             ENDIF
           ELSE
             WGW(M,NQZ) = WG(M,NQZ)*FGW + WTMW*WDGW(M,NQZ)
           ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        K = KD(N)
        DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGW(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGW,INDX )
          WDGW(M,NQZ) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N)
     &      - XMGWB(MN,NB)*RHOMGB(MN,NB))/DZGF(N)
          FGWP = XGW(MN,N)*RHOG(MN,N)
          FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
          INDX = 3
          FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
            IF( WG(1,NQZ).GT.EPSL ) THEN
              WGW(M,NQZ) = WG(M,NQZ)*FGW
            ELSE
              WGW(M,NQZ) = 0.D+0
            ENDIF
          ELSE
            WGW(M,NQZ) = WG(M,NQZ)*FGW + WTMW*WDGW(M,NQZ)
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWW_GT( N,NB )
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
!     Compute water vapor mole diffusion rates on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWW_GT'
      NPX = NSX(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        I = ID(N)
        DXMGW = XMGWB(2,NB)-XMGW(2,N)
        DO M = 1,ISVF
           MP = MPOS(M)
           DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &       *DFGW(MP,N)*RHOMG(MP,N)
           DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*
     &       DFGWB(MP,NB)*RHOMGB(MP,NB)
           INDX = 12
           DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGW,INDX )
           UDGW(M,NPX) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))/DXGF(N)
           FGWP = XGW(MP,N)*RHOG(MP,N)
           FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
           INDX = 3
           FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---       Dirichlet-Outflow boundary condition  ---
!
           IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
             IF( UG(1,NPX).LT.-EPSL ) THEN
               UGW(M,NPX) = UG(M,NPX)*FGW
             ELSE
               UGW(M,NPX) = 0.D+0
             ENDIF
           ELSE
             UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
           ENDIF
         ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
       I = ID(N)
       DXMGW = XMGWB(2,NB)*RHOMGB(2,NB)-XMGW(2,N)*RHOMG(2,N)
       DO M = 1,ISVF
         MP = MPOS(M)
         DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGW(MP,N)
         DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
         INDX = 12
         DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGW,INDX )
         UDGW(M,NPX) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB)
     &     - XMGW(MP,N)*RHOMG(MP,N))/DXGF(N)
         FGWP = XGW(MP,N)*RHOG(MP,N)
         FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
         INDX = 3
         FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---     Dirichlet-Outflow boundary condition  ---
!
         IF( MOD(IBCT(2,NB),100).EQ.19 ) THEN
           IF( UG(1,NPX).LT.-EPSL ) THEN
             UGW(M,NPX) = UG(M,NPX)*FGW
           ELSE
             UGW(M,NPX) = 0.D+0
           ENDIF
         ELSE
           UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
         ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAB_GT( N,NB )
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
!     Compute dissolved air mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLAB_GT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        K = KD(N)
        DXMLA = XMLAB(2,NB)-XMLA(2,N)
        DO M = 1,ISVF
           MP = MPOS(M)
           DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &       *RHOML(MP,N)
           DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &       *RHOMLB(MP,NB)
           INDX = 14
           DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
           WDLA(M,NSZ(N)) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))
     &      /DZGF(N)
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        K = KD(N)
        DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,NSZ(N)) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB)
     &      - XMLA(MP,N)*RHOML(MP,N))/DZGF(N)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAE_GT( N,NB )
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
!     Compute dissolved air mole diffusion rates on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLAE_GT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        I = ID(N)
        DXMLA = XMLA(2,N)-XMLAB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &     *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,NSX(N)+1) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))
     &     /DXGF(N)
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        I = ID(N)
        DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,NSX(N)+1) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N)
     &      - XMLAB(MN,NB)*RHOMLB(MN,NB))/DXGF(N)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAN_GT( N,NB )
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
!     Compute dissolved air mole diffusion rates on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLAN_GT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        I = ID(N)
        J = JD(N)
        DXMLA = XMLA(2,N)-XMLAB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,NSY(N)+IFLD) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))
     &      /(DYGF(N)*RP(I))
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        I = ID(N)
        J = JD(N)
        DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,NSY(N)+IFLD) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N)
     &      - XMLAB(MN,NB)*RHOMLB(MN,NB))/(DYGF(N)*RP(I))
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAS_GT( N,NB )
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
!     Compute dissolved air mole diffusion rates on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLAS_GT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        I = ID(N)
        J = JD(N)
        DXMLA = XMLAB(2,NB)-XMLA(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,NSY(N)) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/
     &      (DYGF(N)*RP(I))
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        I = ID(N)
        J = JD(N)
        DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,NSY(N)) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB)
     &      -XMLA(MP,N)*RHOML(MP,N))/(DYGF(N)*RP(I))
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAT_GT( N,NB )
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
!     Compute dissolved air mole diffusion rates on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLAT_GT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        K = KD(N)
        DXMLA = XMLA(2,N)-XMLAB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,NSZ(N)+IJFLD) = 2.D+0*DFM*(XMLA(MN,N)-
     &      XMLAB(MN,NB))/DZGF(N)
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        K = KD(N)
        DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,NSZ(N)+IJFLD) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N)
     &      - XMLAB(MN,NB)*RHOMLB(MN,NB))/DZGF(N)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAW_GT( N,NB )
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
!     Compute dissolved air mole diffusion rates on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLAW_GT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        I = ID(N)
        DXMLA = XMLAB(2,NB)-XMLA(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,NSX(N)) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/DXGF(N)
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        I = ID(N)
        DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,NSX(N)) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB)
     &      - XMLA(MP,N)*RHOML(MP,N))/DXGF(N)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSB_GT( N,NB )
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
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSB_GT'
      K = KD(N)
      IZN = IZ(N)
      NPZ = NSZ(N)
      IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
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
!---   Inflow/outflow boundary types  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.1 .OR. MOD(IBCT(2,NB),100).EQ.2 .OR.
     &    MOD(IBCT(2,NB),100).EQ.4 .OR. MOD(IBCT(2,NB),100).EQ.12 .OR.
     &    MOD(IBCT(2,NB),100).EQ.13 ) THEN
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
            IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.6 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.9 ) THEN
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.5 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.7 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            WS(M,NPZ) = XLSB(1,NB)*RHOLB(1,NB)*MAX( WL(1,NPZ),ZERO )
            NQZ = NPZ+IJFLD
            IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSS_GT( N,NB )
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
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSS_GT'
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
!---   Inflow/outflow boundary types  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.1 .OR. MOD(IBCT(2,NB),100).EQ.2 .OR.
     &    MOD(IBCT(2,NB),100).EQ.4 .OR. MOD(IBCT(2,NB),100).EQ.12 .OR.
     &    MOD(IBCT(2,NB),100).EQ.13 ) THEN
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
            IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.6 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.9 ) THEN
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.5 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.7 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            VS(M,NPY) = XLSB(1,NB)*RHOLB(1,NB)*MAX( VL(1,NPY),ZERO )
            NQY = NPY+IFLD
            IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
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

      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSW_GT( N,NB )
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
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSW_GT'
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
!---   Inflow/outflow boundary types  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.1 .OR. MOD(IBCT(2,NB),100).EQ.2 .OR.
     &    MOD(IBCT(2,NB),100).EQ.4 .OR. MOD(IBCT(2,NB),100).EQ.12 .OR.
     &    MOD(IBCT(2,NB),100).EQ.13 ) THEN
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
            IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.6 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.9 ) THEN
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.5 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.7 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            US(M,NPX) = XLSB(1,NB)*RHOLB(1,NB)*MAX( UL(1,NPX),ZERO )
            NQX = NPX+1
            IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSE_GT( N,NB )
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
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSE_GT'
      I = ID(N)
      IZN = IZ(N)
      NQX = NSX(N)+1
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
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
!---   Inflow/outflow boundary types  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.1 .OR. MOD(IBCT(2,NB),100).EQ.2 .OR.
     &    MOD(IBCT(2,NB),100).EQ.4 .OR. MOD(IBCT(2,NB),100).EQ.12 .OR.
     &    MOD(IBCT(2,NB),100).EQ.13 ) THEN
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.6 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.9 ) THEN
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.5 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.7 ) THEN
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSN_GT( N,NB )
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
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSN_GT'
      J = JD(N)
      I = ID(N)
      IZN = IZ(N)
      NQY = NSY(N)+IFLD
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
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
!---   Inflow/outflow boundary types  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.1 .OR. MOD(IBCT(2,NB),100).EQ.2 .OR.
     &    MOD(IBCT(2,NB),100).EQ.4 .OR. MOD(IBCT(2,NB),100).EQ.12 .OR.
     &    MOD(IBCT(2,NB),100).EQ.13 ) THEN
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.6 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.9 ) THEN
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.5 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.7 ) THEN
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLST_GT( N,NB )
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
!     Written by MD White, PNNL, 26 February 2015.
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
      SUB_LOG(ISUB_LOG) = '/DFFLST_GT'
      K = KD(N)
      IZN = IZ(N)
      NQZ = NSZ(N)+IJFLD
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
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
!---   Inflow/outflow boundary types  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.1 .OR. MOD(IBCT(2,NB),100).EQ.2 .OR.
     &    MOD(IBCT(2,NB),100).EQ.4 .OR. MOD(IBCT(2,NB),100).EQ.12 .OR.
     &    MOD(IBCT(2,NB),100).EQ.13 ) THEN
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
            IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.6 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.9 ) THEN
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
        ELSEIF( MOD(IBCT(2,NB),100).EQ.5 .OR. 
     &    MOD(IBCT(2,NB),100).EQ.7 ) THEN
          IF( ISLC(6).EQ.1 ) THEN
!
!---  TVD salt transport  --
!
            WS(M,NQZ) = XLSB(1,NB)*RHOLB(1,NB)*MIN( WL(1,NQZ),ZERO )
            NPZ = NSZ(N)
            IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLST_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVG_GT
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 KGM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVG_GT'
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
            HDGX = PG(MN,NW) - PG(MP,N) - 0.5D+0*GRVX(NPX)*
     &          (RHOG(MN,NW)*DXGF(N)+RHOG(MP,N)*DXGF(NW))
            IF( M.EQ.1 ) HDG = HDGX
!
!---        Permeability reduction factor
!
            PERM_WX = PERMRF(MN,NW)*PERM(1,IZ(NW))
            PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
            INDX = 11
            KGM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDG,INDX)
            IF( PERM_WX/EPSL.LT.EPSL ) KGM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
            INDX = 9
            RKGM = DIFMN(RKG(MN,NW),RKG(MP,N),DXGF(NW),DXGF(N),
     &        HDG,INDX)
            INDX = 6
            VGM = DIFMN(VISG(MN,NW),VISG(MP,N),DXGF(NW),DXGF(N),
     &        HDG,INDX)
            UG(M,NPX) = KGM*RKGM*HDGX/DXGP(NPX)/VGM
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
            HDGY = PG(MN,NS) - PG(MP,N) - 0.5D+0*GRVY(NPY)*
     &          (RHOG(MN,NS)*DYGF(N)+RHOG(MP,N)*DYGF(NS))
            IF( M.EQ.1 ) HDG = HDGY
!
!---        Permeability reduction factor
!
            PERM_SX = PERMRF(MN,NS)*PERM(2,IZ(NS))
            PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
            INDX = 11
            KGM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDG,INDX)
            IF( PERM_SX/EPSL.LT.EPSL ) KGM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
            INDX = 9
            RKGM = DIFMN(RKG(MN,NS),RKG(MP,N),DYGF(NS),DYGF(N),
     &        HDG,INDX)
            INDX = 6
            VGM = DIFMN(VISG(MN,NS),VISG(MP,N),DYGF(NS),DYGF(N),
     &        HDG,INDX)
            VG(M,NPY) = KGM*RKGM*HDGY/DYGP(NPY)/VGM/RP(I)
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
            HDGZ = PG(MN,NB) - PG(MP,N) - 0.5D+0*GRVZ(NPZ)*
     &          (RHOG(MN,NB)*DZGF(N)+RHOG(MP,N)*DZGF(NB))
            IF( M.EQ.1 ) HDG = HDGZ
!
!---        Permeability reduction factor
!
            PERM_BX = PERMRF(MN,NB)*PERM(3,IZ(NB))
            PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
            INDX = 11
            KGM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDG,INDX)
            IF( PERM_BX/EPSL.LT.EPSL ) KGM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
            INDX = 9
            RKGM = DIFMN(RKG(MN,NB),RKG(MP,N),DZGF(NB),DZGF(N),
     &        HDG,INDX)
            INDX = 6
            VGM = DIFMN(VISG(MN,NB),VISG(MP,N),DZGF(NB),DZGF(N),
     &        HDG,INDX)
            WG(M,NPZ) = KGM*RKGM*HDGZ/DZGP(NPZ)/VGM
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVG_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVL_GT
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVL_GT'
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
            INDX = 11
            HDLX = PL(MN,NW) - PL(MP,N) - 0.5D+0*GRVX(NPX)*
     &        (RHOL(MN,NW)*DXGF(N)+RHOL(MP,N)*DXGF(NW))
            IF( M.EQ.1 ) HDL = HDLX
!
!---        Permeability reduction factor
!
            PERM_WX = PERMRF(MN,NW)*PERM(1,IZ(NW))
            PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
            INDX = 11
            KLM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDL,INDX)
            INDX = 8
            RKLM = DIFMN(RKL(1,MN,NW),RKL(1,MP,N),DXGF(NW),
     &        DXGF(N),HDL,INDX)
            IF( PERM_WX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 5
            VLM = DIFMN(VISL(MN,NW),VISL(MP,N),DXGF(NW),DXGF(N),
     &        HDL,INDX)
            UL(M,NPX) = KLM*RKLM*HDLX/DXGP(NPX)/VLM
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
            INDX = 11
            HDLY = PL(MN,NS) - PL(MP,N) - 0.5D+0*GRVY(NPY)*
     &        (RHOL(MN,NS)*DYGF(N)+RHOL(MP,N)*DYGF(NS))
            IF( M.EQ.1 ) HDL = HDLY
!
!---        Permeability reduction factor
!
            PERM_SX = PERMRF(MN,NS)*PERM(2,IZ(NS))
            PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
            INDX = 11
            KLM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDL,INDX)
            INDX = 8
            RKLM = DIFMN(RKL(2,MN,NS),RKL(2,MP,N),DYGF(NS),
     &        DYGF(N),HDL,INDX)
            IF( PERM_SX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 5
            VLM = DIFMN(VISL(MN,NS),VISL(MP,N),DYGF(NS),DYGF(N),
     &        HDL,INDX)
            VL(M,NPY) = KLM*RKLM*HDLY/DYGP(NPY)/VLM/RP(I)
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
            INDX = 11
            HDLZ = PL(MN,NB) - PL(MP,N) - 0.5D+0*GRVZ(NPZ)*
     &        (RHOL(MN,NB)*DZGF(N)+RHOL(MP,N)*DZGF(NB))
            IF( M.EQ.1 ) HDL = HDLZ
!
!---        Permeability reduction factor
!
            PERM_BX = PERMRF(MN,NB)*PERM(3,IZ(NB))
            PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
            INDX = 11
            KLM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDL,INDX)
            INDX = 8
            RKLM = DIFMN(RKL(3,MN,NB),RKL(3,MP,N),DZGF(NB),
     &        DZGF(N),HDL,INDX)
            IF( PERM_BX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 5
            VLM = DIFMN(VISL(MN,NB),VISL(MP,N),DZGF(NB),DZGF(N),
     &        HDL,INDX)
            WL(M,NPZ) = KLM*RKLM*HDLZ/DZGP(NPZ)/VLM
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVL_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGB_GT( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      SUB_LOG(ISUB_LOG) = '/DRCVGB_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPZ = NSZ(N)
      IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HZ = PGB(MP,NB) - PG(MP,N)
     &    - 5.D-1*GRVZ(NPZ)*DZGF(N)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
        WG(M,NPZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          WG(M,NPZ) = MIN( 0.D+0,WG(M,NPZ) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGS_GT( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      SUB_LOG(ISUB_LOG) = '/DRCVGS_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPY = NSY(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HY = PGB(MP,NB) - PG(MP,N)
     &    - 5.D-1*GRVY(NPY)*DYGF(N)*RP(I)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
        VG(M,NPY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          VG(M,NPY) = MIN( 0.D+0,VG(M,NPY) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGW_GT( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      SUB_LOG(ISUB_LOG) = '/DRCVGW_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPX = NSX(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HX = PGB(MP,NB) - PG(MP,N)
     &    - 5.D-1*GRVX(NPX)*DXGF(N)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
        UG(M,NPX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          UG(M,NPX) = MIN( 0.D+0,UG(M,NPX) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGE_GT( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      SUB_LOG(ISUB_LOG) = '/DRCVGE_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQX = NSX(N)+1
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
      DO M = 1,ISVF
        MN = MNEG(M)
        HX = PG(MN,N) - PGB(MN,NB)
     &    - 5.D-1*GRVX(NQX)*DXGF(N)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(1,IZ(N))
        UG(M,NQX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          UG(M,NQX) = MAX( 0.D+0,UG(M,NQX) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGN_GT( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      SUB_LOG(ISUB_LOG) = '/DRCVGN_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQY = NSY(N)+IFLD
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
      DO M = 1,ISVF
        MN = MNEG(M)
        HY = PG(MN,N) - PGB(MN,NB)
     &    - 5.D-1*GRVY(NQY)*DYGF(N)*RP(I)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(2,IZ(N))
        VG(M,NQY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          VG(M,NQY) = MAX( 0.D+0,VG(M,NQY) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGT_GT( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      SUB_LOG(ISUB_LOG) = '/DRCVGT_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
      DO M = 1,ISVF
        MN = MNEG(M)
        HZ = PG(MN,N) - PGB(MN,NB)
     &    - 5.D-1*GRVZ(NQZ)*DZGF(N)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(3,IZ(N))
        WG(M,NQZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          WG(M,NQZ) = MAX( 0.D+0,WG(M,NQZ) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLB_GT( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLB_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPZ = NSZ(N)
      IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        HZ = PLB(MP,NB) - PL(MP,N)
     &    -5.D-1*GRVZ(NPZ)*DZGF(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKLB(3,MP,NB),RKL(3,MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
        WL(M,NPZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          WL(M,NPZ) = MIN( 0.D+0,WL(M,NPZ) )
!
!---    Aqueous seepage boundary condition  ---
!
        ELSEIF( MOD(IBCT(2,NB),100).EQ.13 ) THEN 
          PLBX = PL(MP,N) + RHOL(MP,N)*5.D-1*GRVZ(NPZ)*DZGF(N)
          PBX = MAX( PGB(MP,NB),PLB(MP,NB) )
          IF( PLBX.GT.PBX ) THEN
            WL(M,NSZ(N)) = MIN( 0.D+0,WL(M,NSZ(N)) )
          ELSE
            WL(M,NSZ(N)) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLS_GT( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLS_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPY = NSY(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        HY = PLB(MP,NB) - PL(MP,N)
     &    -5.D-1*GRVY(NPY)*DYGF(N)*RP(I)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKLB(2,MP,NB),RKL(2,MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
        VL(M,NPY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          VL(M,NPY) = MIN( 0.D+0,VL(M,NPY) )
!
!---    Aqueous seepage boundary condition  ---
!
        ELSEIF( MOD(IBCT(2,NB),100).EQ.13 ) THEN 
          PLBX = PL(MP,N) + RHOL(MP,N)*5.D-1*GRVY(NPY)*DYGF(N)*RP(I)
          PBX = MAX( PGB(MP,NB),PLB(MP,NB) )
          IF( PLBX.GT.PBX ) THEN
            VL(M,NPY) = MIN( 0.D+0,VL(M,NPY) )
          ELSE
            VL(M,NPY) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLW_GT( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLW_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPX = NSX(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        HX = PLB(MP,NB) - PL(MP,N)
     &    -5.D-1*GRVX(NPX)*DXGF(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKLB(1,MP,NB),RKL(1,MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
        UL(M,NPX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          UL(M,NPX) = MIN( 0.D+0,UL(M,NPX) )
!
!---    Aqueous seepage boundary condition  ---
!
        ELSEIF( MOD(IBCT(2,NB),100).EQ.13 ) THEN 
          PLBX = PL(MP,N) + RHOL(MP,N)*5.D-1*GRVX(NPX)*DXGF(N)
          PBX = MAX( PGB(MP,NB),PLB(MP,NB) )
          IF( PLBX.GT.PBX ) THEN
            UL(M,NPX) = MIN( 0.D+0,UL(M,NPX) )
          ELSE
            UL(M,NPX) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLE_GT( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLE_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQX = NSX(N)+1
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
      DO M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        HX = PL(MN,N) - PLB(MN,NB)
     &    -5.D-1*GRVX(NQX)*DXGF(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKL(1,MN,N),RKLB(1,MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(1,IZ(N))
        UL(M,NQX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          UL(M,NQX) = MAX( 0.D+0,UL(M,NQX) )
!
!---    Aqueous seepage boundary condition  ---
!
        ELSEIF( MOD(IBCT(2,NB),100).EQ.13 ) THEN 
          PLBX = PL(MN,N) - RHOL(MN,N)*5.D-1*GRVX(NQX)*DXGF(N)
          PBX = MAX( PGB(MN,NB),PLB(MN,NB) )
          IF( PLBX.GT.PBX ) THEN
            UL(M,NQX) = MAX( 0.D+0,UL(M,NQX) )
          ELSE
            UL(M,NQX) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLN_GT( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLN_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQY = NSY(N)+IFLD
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
      DO M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        HY = PL(MN,N) - PLB(MN,NB)
     &    -5.D-1*GRVY(NQY)*DYGF(N)*RP(I)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKL(2,MN,N),RKLB(2,MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(2,IZ(N))
        VL(M,NQY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          VL(M,NQY) = MAX( 0.D+0,VL(M,NQY))
!
!---    Aqueous seepage boundary condition  ---
!
        ELSEIF( MOD(IBCT(2,NB),100).EQ.13 ) THEN 
          PLBX = PL(MN,N) - RHOL(MN,N)*5.D-1*GRVY(NQY)*DYGF(N)*RP(I)
          PBX = MAX( PGB(MN,NB),PLB(MN,NB) )
          IF( PLBX.GT.PBX ) THEN
            VL(M,NQY) = MAX( 0.D+0,VL(M,NQY))
          ELSE
            VL(M,NQY) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLT_GT( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLT_GT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
      DO M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        HZ = PL(MN,N) - PLB(MN,NB)
     &    -5.D-1*GRVZ(NQZ)*DZGF(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKL(3,MN,N),RKLB(3,MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(3,IZ(N))
        WL(M,NQZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( MOD(IBCT(2,NB),100).EQ.9 ) THEN
          WL(M,NQZ) = MAX( 0.D+0,WL(M,NQZ) )
!
!---    Aqueous seepage boundary condition  ---
!
        ELSEIF( MOD(IBCT(2,NB),100).EQ.13 ) THEN 
          PLBX = PL(MN,N) - RHOL(MN,N)*5.D-1*GRVZ(NQZ)*DZGF(N)
          PBX = MAX( PGB(MN,NB),PLB(MN,NB) )
          IF( PLBX.GT.PBX ) THEN
            WL(M,NQZ) = MAX( 0.D+0,WL(M,NQZ) )
          ELSE
            WL(M,NQZ) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THD_GT
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
!     Compute the contribution to the energy flux by thermal conduction
!     for nonboundary node faces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THD_GT'
!
!---  X-direction thermal conduction, excluding boundaries
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
          DTK = T(2,NW)-T(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN          
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZ(N)) + 
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKG(MP,N)*SG(MP,N))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
              TKP = THKS(1,IZ(N)) + 
     &          SQRT(SL(MP,N))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
            ENDIF
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(NW)).EQ.2 ) THEN          
              TKW = MAX(1.D+0-PORD(MN,NW),0.D+0)*THKS(1,IZ(NW)) +
     &          PORD(MN,NW)*(THKL(MN,NW)*SL(MN,NW) +
     &          THKG(MN,NW)*SG(MN,NW))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(NW)).EQ.4 ) THEN
              TKW = THKS(1,IZ(NW)) + 
     &          SQRT(SL(MN,NW))*(THKS(4,IZ(NW))-THKS(1,IZ(NW)))           
            ENDIF
            INDX = 1
            TK = DIFMN( TKW,TKP,DXGF(NW),DXGF(N),DTK,INDX )
            UQ(M,NPX) = TK*(T(MN,NW)-T(MP,N))/DXGP(NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction thermal conduction, excluding boundaries
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
          DTK = T(2,NS)-T(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN          
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZ(N)) +
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKG(MP,N)*SG(MP,N))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
              TKP = THKS(2,IZ(N)) + 
     &          SQRT(SL(MP,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
            ENDIF
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(NS)).EQ.2 ) THEN          
              TKS = MAX(1.D+0-PORD(MN,NS),0.D+0)*THKS(2,IZ(NS)) +
     &          PORD(MN,NS)*(THKL(MN,NS)*SL(MN,NS) + 
     &          THKG(MN,NS)*SG(MN,NS))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(NS)).EQ.4 ) THEN
              TKS = THKS(2,IZ(NS)) + 
     &          SQRT(SL(MN,NS))*(THKS(5,IZ(NS))-THKS(2,IZ(NS)))           
            ENDIF     
            INDX = 1
            TK = DIFMN( TKS,TKP,DYGF(NS),DYGF(N),DTK,INDX )
            VQ(M,NPY) = TK*(T(MN,NS)-T(MP,N))/(DYGP(NPY)*RP(I))
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction thermal conduction, excluding boundaries
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
          DTK = T(2,NB)-T(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN          
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZ(N)) +
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKG(MP,N)*SG(MP,N))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
              TKP = THKS(3,IZ(N)) + 
     &          SQRT(SL(MP,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
            ENDIF
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(NB)).EQ.2 ) THEN          
              TKB = MAX(1.D+0-PORD(MN,NB),0.D+0)*THKS(3,IZ(NB)) +
     &          PORD(MN,NB)*(THKL(MN,NB)*SL(MN,NB) + 
     &          THKG(MN,NB)*SG(MN,NB))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(NB)).EQ.4 ) THEN
              TKB = THKS(3,IZ(NB)) + 
     &          SQRT(SL(MN,NB))*(THKS(6,IZ(NB))-THKS(3,IZ(NB)))           
            ENDIF     
            INDX = 1
            TK = DIFMN( TKB,TKP,DZGF(NB),DZGF(N),DTK,INDX )
            WQ(M,NPZ) = TK*(T(MN,NB)-T(MP,N))/DZGP(NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THD_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAG_GT
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
!     Compute the contribution to the energy flux by gas advection
!     for nonboundary node faces.
!     Donor cell interfacial averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
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
      SUB_LOG(ISUB_LOG) = '/THAG_GT'
!
!---  X-direction gas advection, excluding boundaries
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
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HG(MP,N)*RHOG(MP,N)
            HW = HG(MN,NW)*RHOG(MN,NW)
            UQ(M,NPX) = UQ(M,NPX) + HW*MAX(UG(M,NPX),ZERO)
     &        -HP*MAX(-UG(M,NPX),ZERO)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction gas advection, excluding boundaries
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
            HP = HG(MP,N)*RHOG(MP,N)
            HS = HG(MN,NS)*RHOG(MN,NS)
            VQ(M,NPY) = VQ(M,NPY) + HS*MAX(VG(M,NPY),ZERO)
     &        -HP*MAX(-VG(M,NPY),ZERO)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction gas advection, excluding boundaries
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
            HP = HG(MP,N)*RHOG(MP,N)
            HB = HG(MN,NB)*RHOG(MN,NB)
            WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WG(M,NPZ),ZERO)
     &        -HP*MAX(-WG(M,NPZ),ZERO)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAG_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAL_GT
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
!     Compute the contribution to the energy flux by aqueous advection
!     for nonboundary node faces.
!     Donor cell interfacial averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
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
      SUB_LOG(ISUB_LOG) = '/THAL_GT'
!
!---  X-direction aqueous advection, excluding boundaries
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
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HL(MP,N)*RHOL(MP,N)
            HW = HL(MN,NW)*RHOL(MN,NW)
            UQ(M,NPX) = UQ(M,NPX) + HW*MAX(UL(M,NPX),ZERO)
     &        -HP*MAX(-UL(M,NPX),ZERO)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction aqueous advection, excluding boundaries
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
            HP = HL(MP,N)*RHOL(MP,N)
            HS = HL(MN,NS)*RHOL(MN,NS)
            VQ(M,NPY) = VQ(M,NPY) + HS*MAX(VL(M,NPY),ZERO)
     &        -HP*MAX(-VL(M,NPY),ZERO)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction aqueous advection, excluding boundaries
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
            HP = HL(MP,N)*RHOL(MP,N)
            HB = HL(MN,NB)*RHOL(MN,NB)
            WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WL(M,NPZ),ZERO)
     &        -HP*MAX(-WL(M,NPZ),ZERO)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAL_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDG_GT
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
!     Compute the contribution to the energy flux by water vapor,
!     voc vapor, and air molecular diffusion for nonboundary node faces.
!     Donor cell interfacial averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
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
      SUB_LOG(ISUB_LOG) = '/THDG_GT'
!
!---  X-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
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
            MP = MPOS(M)
            MN = MNEG(M)
            UDGAX = -UDGW(M,NPX)
            UQ(M,NPX) = UQ(M,NPX) 
     &        + HGW(MN,NW)*MAX(UDGW(M,NPX),ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-UDGW(M,NPX),ZERO)*WTMW
     &        + HGA(MN,NW)*MAX(UDGA(M,NPX),ZERO)*WTMA
     &        - HGA(MP,N)*MAX(-UDGAX,ZERO)*WTMA
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
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
            VDGAX = -VDGW(M,NPY)
            VQ(M,NPY) = VQ(M,NPY) 
     &        + HGW(MN,NS)*MAX(VDGW(M,NPY),ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-VDGW(M,NPY),ZERO)*WTMW
     &        + HGA(MN,NS)*MAX(VDGAX,ZERO)*WTMA
     &        - HGA(MP,N)*MAX(-VDGAX,ZERO)*WTMA
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
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
            WDGAX = -WDGW(M,NPZ)
            WQ(M,NPZ) = WQ(M,NPZ) 
     &        + HGW(MN,NB)*MAX(WDGW(M,NPZ),ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-WDGW(M,NPZ),ZERO)*WTMW
     &        + HGA(MN,NB)*MAX(WDGAX,ZERO)*WTMA
     &        - HGA(MP,N)*MAX(-WDGAX,ZERO)*WTMA
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDG_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGB_GT( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a bottom boundary
!     flux_gt.F 1336 2020-06-23 14:31:14Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THAGB_GT'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WG(M,NPZ),ZERO)
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WG(M,NPZ),ZERO)
     &      -HP*MAX(-WG(M,NPZ),ZERO)
        ENDIF 
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGE_GT( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THAGE_GT'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UG(M,NQX),ZERO)
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UG(M,NQX),ZERO)
     &      -HB*MAX(-UG(M,NQX),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGN_GT( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THAGN_GT'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VG(M,NQY),ZERO)
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VG(M,NQY),ZERO)
     &      -HB*MAX(-VG(M,NQY),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGS_GT( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THAGS_GT'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VG(M,NPY),ZERO)
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VG(M,NPY),ZERO)
     &      -HP*MAX(-VG(M,NPY),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGT_GT( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THAGT_GT'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WG(M,NQZ),ZERO)
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WG(M,NQZ),ZERO)
     &      - HB*MAX(-WG(M,NQZ),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGW_GT( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THAGW_GT'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UG(M,NPX),ZERO)
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UG(M,NPX),ZERO)
     &      -HP*MAX(-UG(M,NPX),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALB_GT( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a bottom boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THALB_GT'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WL(M,NPZ),ZERO)
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WL(M,NPZ),ZERO)
     &      -HP*MAX(-WL(M,NPZ),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALE_GT( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THALE_GT'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UL(M,NQX),ZERO)
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UL(M,NQX),ZERO)
     &      - HB*MAX(-UL(M,NQX),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALN_GT( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THALN_GT'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VL(M,NQY),ZERO)
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VL(M,NQY),ZERO)
     &      -HB*MAX(-VL(M,NQY),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALS_GT( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THALS_GT'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VL(M,NPY),ZERO)
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VL(M,NPY),ZERO)
     &      -HP*MAX(-VL(M,NPY),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALT_GT( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THALT_GT'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WL(M,NQZ),ZERO)
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WL(M,NQZ),ZERO)
     &      - HB*MAX(-WL(M,NQZ),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALW_GT( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THALW_GT'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
!
!---    Energy outflow  ---
!
        IF( IBCT(1,NB).EQ.6 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UL(M,NPX),ZERO)
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UL(M,NPX),ZERO)
     &      - HP*MAX(-UL(M,NPX),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGB_GT( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a bottom boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGB_GT'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
        WDGAX = -WDGW(M,NPZ)
        WQ(M,NPZ) = WQ(M,NPZ) 
     &    + HGWB(MP,NB)*MAX(WDGW(M,NPZ),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-WDGW(M,NPZ),ZERO)*WTMW
     &    + HGAB(MP,NB)*MAX(WDGAX,ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-WDGAX,ZERO)*WTMA
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGE_GT( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGE_GT'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        UDGAX = -UDGW(M,NQX)
        UQ(M,NQX) = UQ(M,NQX) 
     &    - HGWB(MN,NB)*MAX(-UDGW(M,NQX),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(UDGW(M,NQX),ZERO)*WTMW
     &    - HGAB(MN,NB)*MAX(-UDGAX,ZERO)*WTMA
     &    + HGA(MN,N)*MAX(UDGAX,ZERO)*WTMA
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGN_GT( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGN_GT'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        VDGAX = -VDGW(M,NQY)
        VQ(M,NQY) = VQ(M,NQY) 
     &    - HGWB(MN,NB)*MAX(-VDGW(M,NQY),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(VDGW(M,NQY),ZERO)*WTMW
     &    - HGAB(MN,NB)*MAX(-VDGAX,ZERO)*WTMA
     &    + HGA(MN,N)*MAX(VDGAX,ZERO)*WTMA
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGS_GT( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGS_GT'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        VDGAX = -VDGW(M,NPY)
        VQ(M,NPY) = VQ(M,NPY) 
     &    + HGWB(MP,NB)*MAX(VDGW(M,NPY),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-VDGW(M,NPY),ZERO)*WTMW
     &    + HGAB(MP,NB)*MAX(VDGAX,ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-VDGAX,ZERO)*WTMA
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGT_GT( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGT_GT'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        WDGAX = -WDGW(M,NQZ)
        WQ(M,NQZ) = WQ(M,NQZ) 
     &    - HGWB(MN,NB)*MAX(-WDGW(M,NQZ),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(WDGW(M,NQZ),ZERO)*WTMW
     &    - HGAB(MN,NB)*MAX(-WDGAX,ZERO)*WTMA
     &    + HGA(MN,N)*MAX(WDGAX,ZERO)*WTMA
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGW_GT( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGW_GT'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        UDGAX = -UDGW(M,NPX)
        UQ(M,NPX) = UQ(M,NPX) 
     &    + HGWB(MP,NB)*MAX(UDGW(M,NPX),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-UDGW(M,NPX),ZERO)*WTMW
     &    + HGAB(MP,NB)*MAX(UDGAX,ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-UDGAX,ZERO)*WTMA
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDB_GT( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a bottom boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVT
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDB_GT'
      K = KD(N)
      DTK = TB(2,NB)-T(2,N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N) 
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          TKB = MAX(1.D+0-PORDB(MP,NB),0.D+0)*THKS(3,IZ(N))
     &      + PORDB(MP,NB)*(THKLB(MP,NB)*SLB(MP,NB)
     &      + THKNB(MP,NB)*SNB(MP,NB) + THKGB(MP,NB)*SGB(MP,NB))
          INDX = 1
          NPZ = NSZ(N)
          IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
          TK = DIFMN( TKB,TKP,DZGF(N),DZGF(N),DTK,INDX )
          WQ(M,NPZ) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DZGF(N))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(3,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))     
          TKB = THKS(3,IZ(N))
     &      + SQRT(SLB(MP,NB))*(THKS(6,IZ(N))-THKS(3,IZ(N)))     
          INDX = 1
          NPZ = NSZ(N)
          IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
          TK = DIFMN( TKB,TKP,DZGF(N),DZGF(N),DTK,INDX )
          WQ(M,NPZ) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DZGF(N))
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDE_GT( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVT
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDE_GT'
      I = ID(N)
      DTK = T(2,N)-TB(2,NB)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(1,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N) 
     &      + THKN(MN,N)*SN(MN,N) + THKG(MN,N)*SG(MN,N))
          TKB = MAX(1.D+0-PORDB(MN,NB),0.D+0)*THKS(1,IZ(N))
     &      + PORDB(MN,NB)*(THKLB(MN,NB)*SLB(MN,NB)
     &      + THKNB(MN,NB)*SNB(MN,NB) + THKGB(MN,NB)*SGB(MN,NB))
          INDX = 1
          NQX = NSX(N)+1
          IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
          TK = DIFMN( TKP,TKB,DXGF(N),DXGF(N),DTK,INDX )
          UQ(M,NQX) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DXGF(N))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(1,IZ(N)) 
     &      + SQRT(SL(MN,N))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          TKB = THKS(1,IZ(N))
     &      + SQRT(SLB(MN,NB))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          INDX = 1
          NQX = NSX(N)+1
          IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
          TK = DIFMN( TKP,TKB,DXGF(N),DXGF(N),DTK,INDX )
          UQ(M,NQX) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DXGF(N))
        ENDDO     
      ENDIF 
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDN_GT( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVT
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDN_GT'
      I = ID(N)
      J = JD(N)
      DTK = T(2,N)-TB(2,NB)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(2,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N)
     &      + THKN(MN,N)*SN(MN,N) + THKG(MN,N)*SG(MN,N))
          TKB = MAX(1.D+0-PORDB(MN,NB),0.D+0)*THKS(2,IZ(N))
     &      + PORDB(MN,NB)*(THKLB(MN,NB)*SLB(MN,NB)
     &      + THKNB(MN,NB)*SNB(MN,NB) + THKGB(MN,NB)*SGB(MN,NB))
          INDX = 1
          NQY = NSY(N)+IFLD
          IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
          TK = DIFMN( TKP,TKB,DYGF(N),DYGF(N),DTK,INDX )
          VQ(M,NQY) = TK*(T(MN,N)-TB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(2,IZ(N)) 
     &        + SQRT(SL(MN,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          TKB = THKS(2,IZ(N))
     &        + SQRT(SLB(MN,NB))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          INDX = 1
          NQY = NSY(N)+IFLD
          IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
          TK = DIFMN( TKP,TKB,DYGF(N),DYGF(N),DTK,INDX )
          VQ(M,NQY) = TK*(T(MN,N)-TB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDS_GT( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVT
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDS_GT'
      I = ID(N)
      J = JD(N)
      DTK = TB(2,NB)-T(2,N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N)
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          TKB = MAX(1.D+0-PORDB(MP,NB),0.D+0)*THKS(2,IZ(N))
     &      + PORDB(MP,NB)*(THKLB(MP,NB)*SLB(MP,NB)
     &      + THKNB(MP,NB)*SNB(MP,NB) + THKGB(MP,NB)*SGB(MP,NB))
        INDX = 1
        NPY = NSY(N)
        TK = DIFMN( TKB,TKP,DYGF(N),DYGF(N),DTK,INDX )
        VQ(M,NPY) = TK*(TB(MP,NB)-T(MP,N))/((5.D-1*DYGF(N))*RP(I))
      ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(2,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          TKB = THKS(2,IZ(N))
     &      + SQRT(SL(MP,NB))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          INDX = 1
          NPY = NSY(N)
          TK = DIFMN( TKB,TKP,DYGF(N),DYGF(N),DTK,INDX )
          VQ(M,NPY) = TK*(TB(MP,NB)-T(MP,N))/((5.D-1*DYGF(N))*RP(I))
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDT_GT( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVT
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDT_GT'
      K = KD(N)
      DTK = T(2,N)-TB(2,NB)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(3,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N)
     &      + THKN(MN,N)*SN(MN,N)  + THKG(MN,N)*SG(MN,N))
          TKB = MAX(1.D+0-PORDB(MN,NB),0.D+0)*THKS(3,IZ(N))
     &      + PORDB(MN,NB)*(THKLB(MN,NB)*SLB(MN,NB)
     &      + THKNB(MN,NB)*SNB(MN,NB) + THKGB(MN,NB)*SGB(MN,NB))
          INDX = 1
          NQZ = NSZ(N)+IJFLD
          IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
          TK = DIFMN( TKP,TKB,DZGF(N),DZGF(N),DTK,INDX )
          WQ(M,NQZ) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DZGF(N))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(3,IZ(N)) 
     &      + SQRT(SL(MN,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
          TKB = THKS(3,IZ(N))
     &      + SQRT(SL(MN,NB))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
          INDX = 1
          NQZ = NSZ(N)+IJFLD
          IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
          TK = DIFMN( TKP,TKB,DZGF(N),DZGF(N),DTK,INDX )
          WQ(M,NQZ) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DZGF(N))
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDW_GT( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVT
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDW_GT'
      I = ID(N)
      DTK = TB(2,NB)-T(2,N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N)
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          TKB = MAX(1.D+0-PORDB(MP,NB),0.D+0)*THKS(1,IZ(N))
     &      + PORDB(MP,NB)*(THKLB(MP,NB)*SLB(MP,NB)
     &      + THKNB(MP,NB)*SNB(MP,NB) + THKGB(MP,NB)*SGB(MP,NB))
          INDX = 1
          NPX = NSX(N)
          TK = DIFMN( TKB,TKP,DXGF(N),DXGF(N),DTK,INDX )
          UQ(M,NPX) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DXGF(N))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(1,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          TKB = THKS(1,IZ(N))
     &      + SQRT(SL(MP,NB))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          INDX = 1
          NPX = NSX(N)
          TK = DIFMN( TKB,TKP,DXGF(N),DXGF(N),DTK,INDX )
          UQ(M,NPX) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DXGF(N))
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDW_GT group  ---
!
      RETURN
      END


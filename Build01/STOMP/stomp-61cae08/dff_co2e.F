!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGW
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
!     Written by MD White, PNNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFGW'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
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
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DXMGW = XMGW(2,NW)-XMGW(2,N)
          DO 100 M = 1,ISVF
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
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
          NPY = NSY(N)
          DXMGW = XMGW(2,NS)-XMGW(2,N)
          DO 300 M = 1,ISVF
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
          IF( NB.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
          NPZ = NSZ(N)
          DXMGW = XMGW(2,NB)-XMGW(2,N)
          DO 500 M = 1,ISVF
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
  500     CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
!
!---  X-direction vapor mole diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 1200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 1200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 1200
          NPX = NSX(N)
          DXMGW = XMGW(2,NW)*RHOMG(2,NW)-XMGW(2,N)*RHOMG(2,N)
          DO 1100 M = 1,ISVF
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
 1100     CONTINUE
 1200   CONTINUE
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 1400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 1400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 1400
          NPY = NSY(N)
          DXMGW = XMGW(2,NS)*RHOMG(2,NS)-XMGW(2,N)*RHOMG(2,N)
          DO 1300 M = 1,ISVF
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
 1300     CONTINUE
 1400   CONTINUE
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 1600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 1600
          NB = ICM(1,1,N)
          IF( NB.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 1600
          NPZ = NSZ(N)
          DXMGW = XMGW(2,NB)*RHOMG(2,NB)-XMGW(2,N)*RHOMG(2,N)
          DO 1500 M = 1,ISVF
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
 1500     CONTINUE
 1600   CONTINUE
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGW group  ---
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
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
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
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DXLA = (XMLA(2,NW)-XMLA(2,N))
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &        *RHOML(MP,N)
            DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*DFLA(MN,NW)
     &        *RHOML(MN,NW)
            INDX = 14
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLA,INDX)
            UDLA(M,NPX) = DFM*(XMLA(MN,NW)-XMLA(MP,N))/DXGP(NPX)
  100     CONTINUE
  200   CONTINUE
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
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
          NPY = NSY(N)
          DXLA = (XMLA(2,NS)-XMLA(2,N))
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &        *RHOML(MP,N)
            DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*DFLA(MN,NS)
     &        *RHOML(MN,NS)
            INDX = 14
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLA,INDX )
            VDLA(M,NPY) = DFM*(XMLA(MN,NS)-XMLA(MP,N))/DYGP(NPY)/RP(I)
  300     CONTINUE
  400   CONTINUE
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
          IF( NB.EQ.0 ) EXIT
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
          NPZ = NSZ(N)
          DXLA = (XMLA(2,NB)-XMLA(2,N))
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &        *RHOML(MP,N)
            DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*DFLA(MN,NB)
     &        *RHOML(MN,NB)
            INDX = 14
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLA,INDX)
            WDLA(M,NPZ) = DFM*(XMLA(MN,NB)-XMLA(MP,N))/DZGP(NPZ)
  500     CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
!
!---  X-direction molar diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 1200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 1200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 1200
          NPX = NSX(N)
          DXLA = (XMLA(2,NW)*RHOML(2,NW)-XMLA(2,N)*RHOML(2,N))
          DO 1100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
            DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*DFLA(MN,NW)
            INDX = 14
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLA,INDX)
            UDLA(M,NPX) = DFM*(XMLA(MN,NW)*RHOML(MN,NW)
     &       - XMLA(MP,N)*RHOML(MP,N))/DXGP(NPX)
 1100     CONTINUE
 1200   CONTINUE
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 1400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 1400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 1400
          NPY = NSY(N)
          DXLA = (XMLA(2,NS)*RHOML(2,NS)-XMLA(2,N)*RHOML(2,N))
          DO 1300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
            DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*DFLA(MN,NS)
            INDX = 14
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLA,INDX )
            VDLA(M,NPY) = DFM*(XMLA(MN,NS)*RHOML(MN,NS)
     &        - XMLA(MP,N)*RHOML(MP,N))/DYGP(NPY)/RP(I)
 1300     CONTINUE
 1400   CONTINUE
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 1600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 1600
          NB = ICM(1,1,N)
          IF( NB.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 1600
          NPZ = NSZ(N)
          DXLA = (XMLA(2,NB)*RHOML(2,NB)-XMLA(2,N)*RHOML(2,N))
          DO 1500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
            DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*DFLA(MN,NB)
            INDX = 14
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLA,INDX)
            WDLA(M,NPZ) = DFM*(XMLA(MN,NB)*RHOML(MN,NB)
     &        - XMLA(MP,N)*RHOML(MP,N))/DZGP(NPZ)
 1500     CONTINUE
 1600   CONTINUE
      ENDIF
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
!     Written by MD White, PNNL, 29 May 2002.
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
          IF( NW.EQ.0 ) CYCLE
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
          IF( NS.EQ.0 ) CYCLE
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
          IF( NB.EQ.0 ) CYCLE
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
      SUBROUTINE DFFGWB( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWB'
      NPZ = NSZ(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      K = KD(N)
      DXMGW = XMGWB(2,NB)-XMGW(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &   *DFGW(MP,N)*RHOMG(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)
     &   *DFGWB(MP,NB)*RHOMGB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGW,INDX )
       WDGW(M,NPZ) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))
     &   /DZGF(N)
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( WG(1,NPZ).LT.-EPSL ) THEN
           WGW(M,NPZ) = WG(M,NPZ)*FGW
         ELSE
           WGW(M,NPZ) = 0.D+0
         ENDIF
       ELSE
         WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      K = KD(N)
      DXMGW = XMGWB(2,NB)*RHOMGB(2,NB) - XMGW(2,N)*RHOMG(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGW(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGW,INDX )
       WDGW(M,NPZ) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB)
     &   - XMGW(MP,N)*RHOMG(MP,N))/DZGF(N)
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( WG(1,NPZ).LT.-EPSL ) THEN
           WGW(M,NPZ) = WG(M,NPZ)*FGW
         ELSE
           WGW(M,NPZ) = 0.D+0
         ENDIF
       ELSE
         WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
       ENDIF
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWE( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWE'
      NQX = NSX(N)+1
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      I = ID(N)
      DXMGW = XMGW(2,N)-XMGWB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &   *DFGW(MN,N)*RHOMG(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*
     &   DFGWB(MN,NB)*RHOMGB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGW,INDX )
       UDGW(M,NQX) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))
     &   /DXGF(N)
       FGWP = XGW(MN,N)*RHOG(MN,N)
       FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
       INDX = 3
       FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( UG(1,NQX).GT.EPSL ) THEN
           UGW(M,NQX) = UG(M,NQX)*FGW
         ELSE
           UGW(M,NQX) = 0.D+0
         ENDIF
       ELSE
         UGW(M,NQX) = UG(M,NQX)*FGW + WTMW*UDGW(M,NQX)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGW(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGW,INDX )
       UDGW(M,NQX) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N)
     &   - XMGWB(MN,NB)*RHOMGB(MN,NB))/DXGF(N)
       FGWP = XGW(MN,N)*RHOG(MN,N)
       FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
       INDX = 3
       FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( UG(1,NQX).GT.EPSL ) THEN
           UGW(M,NQX) = UG(M,NQX)*FGW
         ELSE
           UGW(M,NQX) = 0.D+0
         ENDIF
       ELSE
         UGW(M,NQX) = UG(M,NQX)*FGW + WTMW*UDGW(M,NQX)
       ENDIF
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWN( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/DFFGWN'
      NQY = NSY(N)+IFLD
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      I = ID(N)
      J = JD(N)
      DXMGW = XMGW(2,N)-XMGWB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &   *DFGW(MN,N)*RHOMG(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*
     &   DFGWB(MN,NB)*RHOMGB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGW,INDX )
       VDGW(M,NQY) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))
     &   /(DYGF(N)*RP(I))
       FGWP = XGW(MN,N)*RHOG(MN,N)
       FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
       INDX = 3
       FGW = DIFMN( FGWP,FGWB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( VG(1,NQY).GT.EPSL ) THEN
           VGW(M,NQY) = VG(M,NQY)*FGW
         ELSE
           VGW(M,NQY) = 0.D+0
         ENDIF
       ELSE
         VGW(M,NQY) = VG(M,NQY)*FGW + WTMO*VDGO(M,NQY)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      J = JD(N)
      DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGW(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGW,INDX )
       VDGW(M,NQY) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N)
     &   - XMGWB(MN,NB)*RHOMGB(MN,NB))/(DYGF(N)*RP(I))
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWS( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWS'
      NPY = NSY(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      I = ID(N)
      J = JD(N)
      DXMGW = XMGWB(2,NB)-XMGW(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &   *DFGW(MP,N)*RHOMG(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*
     &   DFGWB(MP,NB)*RHOMGB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGW,INDX )
       VDGW(M,NPY) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))/
     &   (DYGF(N)*RP(I))
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( VG(1,NPY).LT.-EPSL ) THEN
           VGW(M,NPY) = VG(M,NPY)*FGW
         ELSE
           VGW(M,NPY) = 0.D+0
         ENDIF
       ELSE
         VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      J = JD(N)
      DXMGW = XMGWB(2,NB)*RHOMGB(2,NB) - XMGW(2,N)*RHOMG(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGW(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGW,INDX )
       VDGW(M,NPY) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB)
     &   - XMGW(MP,N)*RHOMG(MP,N))/(DYGF(N)*RP(I))
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( VG(1,NPY).LT.-EPSL ) THEN
           VGW(M,NPY) = VG(M,NPY)*FGW
         ELSE
           VGW(M,NPY) = 0.D+0
         ENDIF
       ELSE
         VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
       ENDIF
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWT( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWT'
      NQZ = NSZ(N)+IJFLD
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      K = KD(N)
      DXMGW = XMGW(2,N)-XMGWB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &   *DFGW(MN,N)*RHOMG(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*
     &   DFGWB(MN,NB)*RHOMGB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGW,INDX )
       WDGW(M,NQZ) = 2.D+0*DFM*(XMGW(MN,N)-
     &   XMGWB(MN,NB))/DZGF(N)
       FGWP = XGW(MN,N)*RHOG(MN,N)
       FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
       INDX = 3
       FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( WG(1,NQZ).GT.EPSL ) THEN
           WGW(M,NQZ) = WG(M,NQZ)*FGW
         ELSE
           WGW(M,NQZ) = 0.D+0
         ENDIF
       ELSE
         WGW(M,NQZ) = WG(M,NQZ)*FGW + WTMW*WDGW(M,NQZ)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      K = KD(N)
      DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGW(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGW,INDX )
       WDGW(M,NQZ) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N)
     &   - XMGWB(MN,NB)*RHOMGB(MN,NB))/DZGF(N)
       FGWP = XGW(MN,N)*RHOG(MN,N)
       FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
       INDX = 3
       FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( WG(1,NQZ).GT.EPSL ) THEN
           WGW(M,NQZ) = WG(M,NQZ)*FGW
         ELSE
           WGW(M,NQZ) = 0.D+0
         ENDIF
       ELSE
         WGW(M,NQZ) = WG(M,NQZ)*FGW + WTMW*WDGW(M,NQZ)
       ENDIF
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWW( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWW'
      NPX = NSX(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      I = ID(N)
      DXMGW = XMGWB(2,NB)-XMGW(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &   *DFGW(MP,N)*RHOMG(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*
     &   DFGWB(MP,NB)*RHOMGB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGW,INDX )
       UDGW(M,NPX) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))/DXGF(N)
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( UG(1,NPX).LT.-EPSL ) THEN
           UGW(M,NPX) = UG(M,NPX)*FGW
         ELSE
           UGW(M,NPX) = 0.D+0
         ENDIF
       ELSE
         UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      DXMGW = XMGWB(2,NB)*RHOMGB(2,NB)-XMGW(2,N)*RHOMG(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGW(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGW,INDX )
       UDGW(M,NPX) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB)
     &   - XMGW(MP,N)*RHOMG(MP,N))/DXGF(N)
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( UG(1,NPX).LT.-EPSL ) THEN
           UGW(M,NPX) = UG(M,NPX)*FGW
         ELSE
           UGW(M,NPX) = 0.D+0
         ENDIF
       ELSE
         UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
       ENDIF
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWW group  ---
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
!     Compute dissolved air mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFLAB'
      NPZ = NSZ(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      K = KD(N)
      DXMLA = XMLAB(2,NB)-XMLA(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLA(M,NPZ) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/DZGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      K = KD(N)
      DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLA(M,NPZ) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB)
     &   - XMLA(MP,N)*RHOML(MP,N))/DZGF(N)
  200 CONTINUE
      ENDIF
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
!     Compute dissolved air mole diffusion rates on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFLAE'
      NQX = NSX(N)+1
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      DXMLA = XMLA(2,N)-XMLAB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &  *RHOMLB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLA(M,NQX) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))/DXGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLA(M,NQX) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N)
     &   - XMLAB(MN,NB)*RHOMLB(MN,NB))/DXGF(N)
  200 CONTINUE
      ENDIF
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
!     Compute dissolved air mole diffusion rates on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFLAN'
      NQY = NSY(N)+IFLD
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      J = JD(N)
      DXMLA = XMLA(2,N)-XMLAB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &   *RHOMLB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLA(M,NQY) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))
     &   /(DYGF(N)*RP(I))
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      J = JD(N)
      DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLA(M,NQY) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N)
     &   - XMLAB(MN,NB)*RHOMLB(MN,NB))/(DYGF(N)*RP(I))
  200 CONTINUE
      ENDIF
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
!     Compute dissolved air mole diffusion rates on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFLAS'
      NPY = NSY(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      J = JD(N)
      DXMLA = XMLAB(2,NB)-XMLA(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLA(M,NPY) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/
     &   (DYGF(N)*RP(I))
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      J = JD(N)
      DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLA(M,NPY) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB)
     &   -XMLA(MP,N)*RHOML(MP,N))/(DYGF(N)*RP(I))
  200 CONTINUE
      ENDIF
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
!     Compute dissolved air mole diffusion rates on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFLAT'
      NQZ = NSZ(N)+IJFLD
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      K = KD(N)
      DXMLA = XMLA(2,N)-XMLAB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &   *RHOMLB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLA(M,NQZ) = 2.D+0*DFM*(XMLA(MN,N)-
     &   XMLAB(MN,NB))/DZGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      K = KD(N)
      DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLA(M,NQZ) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N)
     &   - XMLAB(MN,NB)*RHOMLB(MN,NB))/DZGF(N)
  200 CONTINUE
      ENDIF
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
!     Compute dissolved air mole diffusion rates on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/DFFLAW'
      NPX = NSX(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      DXMLA = XMLAB(2,NB)-XMLA(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLA(M,NPX) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/DXGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLA(M,NPX) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB)
     &   - XMLA(MP,N)*RHOML(MP,N))/DXGF(N)
  200 CONTINUE
      ENDIF
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
!     Written by M.D. White, PNNL, 29 May 2002.
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
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
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
!     Written by M.D. White, PNNL, 29 May 2002.
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
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
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
!     Written by M.D. White, PNNL, 29 May 2002.
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
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
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
!     Written by M.D. White, PNNL, 29 May 2002.
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
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
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
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
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
!     Written by M.D. White, PNNL, 29 May 2002.
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
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
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
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
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
!     Written by M.D. White, PNNL, 29 May 2002.
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
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
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
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
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
!---  End of DFFLST group  ---
!
      RETURN
      END

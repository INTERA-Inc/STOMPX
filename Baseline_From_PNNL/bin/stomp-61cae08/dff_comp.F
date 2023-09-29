!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGC(IGC)
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
!     Compute diffusive air or CO2 gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
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
      SUB_LOG(ISUB_LOG) = '/DFFGC'
!
!---  X-direction vapor mole diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 K = 1,KFLD
          DO 200 J = 1,JFLD
            DO 200 I = 2,IFLD
              N = ND(I,J,K)
              NW = N-1
              IF( IXP(N).LE.0 .OR. IXP(N-1).LE.0 ) GOTO 200
              NPX = NSX(N)
              DXGC = XMGC(IGC,2,NW)-XMGC(IGC,2,N)
              DO 100 M = 1,ISVF
                MN = MNEG(M)
                MP = MPOS(M)
                DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &            *DFGC(IGC,MP,N)*RHOMG(MP,N)
                DFW = TORG(MN,NW)*PORD(MN,NW)*(SG(MN,NW)-SGT(MN,NW))
     &            *DFGC(IGC,MN,NW)*RHOMG(MN,NW)
                INDX = 13
                DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXGC,INDX )
                UDGC(IGC,M,NPX) = DFM*(XMGC(IGC,MN,NW)-XMGC(IGC,MP,N))
     &            /DXGP(NPX)
  100         CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 400 I = 1,IFLD
          DO 400 K = 1,KFLD
            DO 400 J = 2,JFLD
              N = ND(I,J,K)
              NS = N-IFLD
              IF( IXP(N).LE.0 .OR. IXP(N-IFLD).LE.0 ) GOTO 400
              NPY = NSY(N)
              DXGC = XMGC(IGC,2,NS)-XMGC(IGC,2,N)
              DO 300 M = 1,ISVF
                MP = MPOS(M)
                MN = MNEG(M)
                DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &            *DFGC(IGC,MP,N)*RHOMG(MP,N)
                DFS = TORG(MN,NS)*PORD(MN,NS)*(SG(MN,NS)-SGT(MN,NS))
     &            *DFGC(IGC,MN,NS)*RHOMG(MN,NS)
                INDX = 13
                DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXGC,INDX )
                VDGC(IGC,M,NPY) = DFM*(XMGC(IGC,MN,NS)-XMGC(IGC,MP,N))
     &            /DYGP(NPY)/RP(I)
  300         CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 600 J = 1,JFLD
          DO 600 I = 1,IFLD
            DO 600 K = 2,KFLD
              N = ND(I,J,K)
              NB = N-IJFLD
              IF( IXP(N).LE.0 .OR. IXP(N-IJFLD).LE.0 ) GOTO 600
              NPZ = NSZ(N)
              DXGC = XMGC(IGC,2,NB)-XMGC(IGC,2,N)
              DO 500 M = 1,ISVF
                MP = MPOS(M)
                MN = MNEG(M)
                DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &            *DFGC(IGC,MP,N)*RHOMG(MP,N)
                DFB = TORG(MN,NB)*PORD(MN,NB)*(SG(MN,NB)-SGT(MN,NB))
     &            *DFGC(IGC,MN,NB)*RHOMG(MN,NB)
                INDX = 13
                DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXGC,INDX )
                WDGC(IGC,M,NPZ) = DFM*(XMGC(IGC,MN,NB)-XMGC(IGC,MP,N))
     &            /DZGP(NPZ)
  500         CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGC group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGW_COMP
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFGW_COMP'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
!
!---  X-direction vapor mole diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 K = 1,KFLD
        DO 200 J = 1,JFLD
        DO 200 I = 2,IFLD
          N = ND(I,J,K)
          NW = N-1
          IF( IXP(N).LE.0 .OR. IXP(N-1).LE.0 ) GOTO 200
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
        DO 400 I = 1,IFLD
        DO 400 K = 1,KFLD
        DO 400 J = 2,JFLD
          N = ND(I,J,K)
          NS = N-IFLD
          IF( IXP(N).LE.0 .OR. IXP(N-IFLD).LE.0 ) GOTO 400
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
        DO 600 J = 1,JFLD
        DO 600 I = 1,IFLD
        DO 600 K = 2,KFLD
          N = ND(I,J,K)
          NB = N-IJFLD
          IF( IXP(N).LE.0 .OR. IXP(N-IJFLD).LE.0 ) GOTO 600
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
        DO 1200 K = 1,KFLD
        DO 1200 J = 1,JFLD
        DO 1200 I = 2,IFLD
          N = ND(I,J,K)
          NW = N-1
          IF( IXP(N).LE.0 .OR. IXP(N-1).LE.0 ) GOTO 1200
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
        DO 1400 I = 1,IFLD
        DO 1400 K = 1,KFLD
        DO 1400 J = 2,JFLD
          N = ND(I,J,K)
          NS = N-IFLD
          IF( IXP(N).LE.0 .OR. IXP(N-IFLD).LE.0 ) GOTO 1400
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
        DO 1600 J = 1,JFLD
        DO 1600 I = 1,IFLD
        DO 1600 K = 2,KFLD
          N = ND(I,J,K)
          NB = N-IJFLD
          IF( IXP(N).LE.0 .OR. IXP(N-IJFLD).LE.0 ) GOTO 1600
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
!---  End of DFFGW_COMP group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLC(IGC)
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
      USE FDVP
      USE FDVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFLC_COMP'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
!
!---  X-direction molar diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
        DO 200 J = 1,JFLD
          DO 200 I = 2,IFLD
            N = ND(I,J,K)
            NW = N-1
            IF( IXP(N).LE.0 .OR. IXP(N-1).LE.0 ) GOTO 200
            NPX = NSX(N)
            DXLC = (XMLC(IGC,2,NW)-XMLC(IGC,2,N))
            DO 100 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
     &          *RHOML(MP,N)
              DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*DFLC(IGC,MN,NW)
     &          *RHOML(MN,NW)
              INDX = 14
              DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLC,INDX)
              UDLC(IGC,M,NPX) = DFM*(XMLC(IGC,MN,NW)-XMLC(IGC,MP,N))
     &             /DXGP(NPX)
  100       CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
      DO 400 I = 1,IFLD
        DO 400 K = 1,KFLD
          DO 400 J = 2,JFLD
            N = ND(I,J,K)
            NS = N-IFLD
            IF( IXP(N).LE.0 .OR. IXP(N-IFLD).LE.0 ) GOTO 400
            NPY = NSY(N)
            DXLC = (XMLC(IGC,2,NS)-XMLC(IGC,2,N))
            DO 300 M = 1,ISVF
              MP = MPOS(M)
              MN = MNEG(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
     &          *RHOML(MP,N)
              DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*DFLC(IGC,MN,NS)
     &          *RHOML(MN,NS)
              INDX = 14
              DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLC,INDX )
              VDLC(IGC,M,NPY) = DFM*(XMLC(IGC,MN,NS)-XMLC(IGC,MP,N))
     &           /DYGP(NPY)/RP(I)
  300       CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
      DO 600 J = 1,JFLD
        DO 600 I = 1,IFLD
          DO 600 K = 2,KFLD
            N = ND(I,J,K)
            NB = N-IJFLD
            IF( IXP(N).LE.0 .OR. IXP(N-IJFLD).LE.0 ) GOTO 600
            NPZ = NSZ(N)
            DXLC = (XMLC(IGC,2,NB)-XMLC(IGC,2,N))
            DO 500 M = 1,ISVF
              MP = MPOS(M)
              MN = MNEG(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
     &          *RHOML(MP,N)
              DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*DFLC(IGC,MN,NB)
     &          *RHOML(MN,NB)
              INDX = 14
              DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLC,INDX)
              WDLC(IGC,M,NPZ) = DFM*(XMLC(IGC,MN,NB)-XMLC(IGC,MP,N))
     &          /DZGP(NPZ)
  500       CONTINUE
  600 CONTINUE
      ENDIF
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
!
!---  X-direction molar diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
      DO 1200 K = 1,KFLD
        DO 1200 J = 1,JFLD
          DO 1200 I = 2,IFLD
            N = ND(I,J,K)
            NW = N-1
            IF( IXP(N).LE.0 .OR. IXP(N-1).LE.0 ) GOTO 1200
            NPX = NSX(N)
            DXLC = (XMLC(IGC,2,NW)*RHOML(2,NW)-XMLC(IGC,2,N)*RHOML(2,N))
            DO 1100 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
              DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*DFLC(IGC,MN,NW)
              INDX = 14
              DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLC,INDX)
              UDLC(IGC,M,NPX) = DFM*(XMLC(IGC,MN,NW)*RHOML(MN,NW)
     &         - XMLC(IGC,MP,N)*RHOML(MP,N))/DXGP(NPX)
 1100       CONTINUE
 1200 CONTINUE
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
      DO 1400 I = 1,IFLD
        DO 1400 K = 1,KFLD
          DO 1400 J = 2,JFLD
            N = ND(I,J,K)
            NS = N-IFLD
            IF( IXP(N).LE.0 .OR. IXP(N-IFLD).LE.0 ) GOTO 1400
            NPY = NSY(N)
            DXLC = (XMLC(IGC,2,NS)*RHOML(2,NS)-XMLC(IGC,2,N)*RHOML(2,N))
            DO 1300 M = 1,ISVF
              MP = MPOS(M)
              MN = MNEG(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
              DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*DFLC(IGC,MN,NS)
              INDX = 14
              DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLC,INDX )
              VDLC(IGC,M,NPY) = DFM*(XMLC(IGC,MN,NS)*RHOML(MN,NS)
     &          - XMLC(IGC,MP,N)*RHOML(MP,N))/DYGP(NPY)/RP(I)
 1300       CONTINUE
 1400 CONTINUE
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
      DO 1600 J = 1,JFLD
        DO 1600 I = 1,IFLD
          DO 1600 K = 2,KFLD
            N = ND(I,J,K)
            NB = N-IJFLD
            IF( IXP(N).LE.0 .OR. IXP(N-IJFLD).LE.0 ) GOTO 1600
            NPZ = NSZ(N)
            DXLC = (XMLC(IGC,2,NB)*RHOML(2,NB)-XMLC(IGC,2,N)*RHOML(2,N))
            DO 1500 M = 1,ISVF
              MP = MPOS(M)
              MN = MNEG(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
              DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*DFLC(IGC,MN,NB)
              INDX = 14
              DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLC,INDX)
              WDLC(IGC,M,NPZ) = DFM*(XMLC(IGC,MN,NB)*RHOML(MN,NB)
     &          - XMLC(IGC,MP,N)*RHOML(MP,N))/DZGP(NPZ)
 1500       CONTINUE
 1600 CONTINUE
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLC_COMP group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLS_COMP
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFLS_COMP'
!
!---  X-direction Darcy velocities, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
        DO 200 J = 1,JFLD
          DO 200 I = 2,IFLD
            N = ND(I,J,K)
            NW = N-1
            IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 ) GOTO 200
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
     &            EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
                DFCLW = SDCLS(1,IZ(NW))*SDCLS(2,IZ(NW))*
     &            EXP(SL(MN,NW)*PORD(MN,NW)*SDCLS(3,IZ(NW)))
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
              DDLW = (DFCLW+DPLW)/DXGP(NPX)
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
            IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 ) GOTO 400
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
     &            EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
                DFCLS = SDCLS(1,IZ(NS))*SDCLS(2,IZ(NS))*
     &            EXP(SL(MN,NS)*PORD(MN,NS)*SDCLS(3,IZ(NS)))
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
                CALL ADVS( PORD,SL,UL,VL,WL,ULS,VLS,WLS,N,M )
                ULS = ULS*ULS
                VLS = VLS*VLS
                WLS = WLS*WLS
                ZVS = SQRT(ULS+VLS+WLS)
                INDX = 17
                DPLS = DIFMN(DPLGS(IZ(NS)),DPLGS(IZ(N)),
     &            DYGF(NS),DYGF(N),VL(1,NPY),INDX)
                DPTS = DIFMN(DPTRS(IZ(NS)),DPTRS(IZ(N)),
     &            DYGF(NS),DYGF(N),VL(1,NPY),INDX)
                DPLS = (DPLS*VLS + DPTS*(WLS+ULS))/(ZVS+SMALL)
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
            IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 ) GOTO 600
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
     &            EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
                DFCLB = SDCLS(1,IZ(NB))*SDCLS(2,IZ(NB))*
     &            EXP(SL(MN,NB)*PORD(MN,NB)*SDCLS(3,IZ(NB)))
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
                CALL ADVB( PORD,SL,UL,VL,WL,ULB,VLB,WLB,N,M )
                ULB = ULB*ULB
                VLB = VLB*VLB
                WLB = WLB*WLB
                ZVB = SQRT(ULB+VLB+WLB)
                INDX = 17
                DPLB = DIFMN(DPLGS(IZ(NB)),DPLGS(IZ(N)),
     &            DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
                DPTB = DIFMN(DPTRS(IZ(NB)),DPTRS(IZ(N)),
     &            DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
                DPLB = (DPLB*WLB + DPTB*(ULB+VLB))/(ZVB+SMALL)
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
!---  End of DFFLS_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCB( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFGCB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMGC = XMGCB(IGC,2,NB)-XMGC(IGC,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGC(IGC,MP,N)
     &    *RHOMG(MP,N)
        DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGCB(IGC,MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGC,INDX )
        WDGC(IGC,M,NPZ) = 2.D+0*DFM*(XMGCB(IGC,MP,NB)-XMGC(IGC,MP,N))
     &    /DZGF(N)
        FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
        FGCB = XGCB(IGC,MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGC = DIFMN( FGCB,FGCP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.19 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGC(IGC,M,NPZ) = WG(M,NPZ)*FGC
          ELSE
            WGC(IGC,M,NPZ) = 0.D+0
          ENDIF
        ELSE
          WGC(IGC,M,NPZ) = WG(M,NPZ)*FGC + GCPP(1,IGC)*WDGC(IGC,M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCE( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFGCE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMGC = XMGC(IGC,2,N)-XMGCB(IGC,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGC(IGC,MN,N)
     &    *RHOMG(MN,N)
        DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGCB(IGC,MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGC,INDX )
        UDGC(IGC,M,NQX) = 2.D+0*DFM*(XMGC(IGC,MN,N)-XMGCB(IGC,MN,NB))
     &     /DXGF(N)
        FGCP = XGC(IGC,MN,N)*RHOG(MN,N)
        FGCB = XGCB(IGC,MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGC = DIFMN( FGCP,FGCB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.19 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGC(IGC,M,NQX) = UG(M,NQX)*FGC
          ELSE
            UGC(IGC,M,NQX) = 0.D+0
          ENDIF
        ELSE
          UGC(IGC,M,NQX) = UG(M,NQX)*FGC + GCPP(1,IGC)*UDGC(IGC,M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCN( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFGCN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGC = XMGC(IGC,2,N)-XMGCB(IGC,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGC(IGC,MN,N)
     &    *RHOMG(MN,N)
        DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGCB(IGC,MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGC,INDX )
        VDGC(IGC,M,NQY) = 2.D+0*DFM*(XMGC(IGC,MN,N)-XMGCB(IGC,MN,NB))
     &    /(DYGF(N)*RP(I))
        FGCP = XGC(IGC,MN,N)*RHOG(MN,N)
        FGCB = XGCB(IGC,MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGC = DIFMN( FGCP,FGCB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.19 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGC(IGC,M,NQY) = VG(M,NQY)*FGC
          ELSE
            VGC(IGC,M,NQY) = 0.D+0
          ENDIF
        ELSE
          VGC(IGC,M,NQY) = VG(M,NQY)*FGC + GCPP(1,IGC)*VDGC(IGC,M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCS( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFGCS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGC = XMGCB(IGC,2,NB)-XMGC(IGC,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGC(IGC,MP,N)
     &    *RHOMG(MP,N)
        DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGCB(IGC,MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGC,INDX )
        VDGC(IGC,M,NPY) = 2.D+0*DFM*(XMGCB(IGC,MP,NB)-XMGC(IGC,MP,N))
     &    /(DYGF(N)*RP(I))
        FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
        FGCB = XGCB(IGC,MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGC = DIFMN( FGCB,FGCP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.19 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGC(IGC,M,NPY) = VG(M,NPY)*FGC
          ELSE
            VGC(IGC,M,NPY) = 0.D+0
          ENDIF
        ELSE
          VGC(IGC,M,NPY) = VG(M,NPY)*FGC + GCPP(1,IGC)*VDGC(IGC,M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCT( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFGCT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGC = XMGC(IGC,2,N)-XMGCB(IGC,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGC(IGC,MN,N)
     &    *RHOMG(MN,N)
        DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGCB(IGC,MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGC,INDX )
        WDGC(IGC,M,NQZ) = 2.D+0*DFM*(XMGC(IGC,MN,N)-XMGCB(IGC,MN,NB))
     &   /DZGF(N)
        FGCP = XGC(IGC,MN,N)*RHOG(MN,N)
        FGCB = XGCB(IGC,MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGC = DIFMN( FGCP,FGCB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.19 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGC(IGC,M,NQZ) = WG(M,NQZ)*FGC
          ELSE
            WGC(IGC,M,NQZ) = 0.D+0
          ENDIF
        ELSE
          WGC(IGC,M,NQZ) = WG(M,NQZ)*FGC + GCPP(1,IGC)*WDGC(IGC,M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCW( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFGCW'
      I = ID(N)
      NPX = NSX(N)
      DXMGA = XMGCB(IGC,2,NB)-XMGC(IGC,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGC(IGC,MP,N)
     &    *RHOMG(MP,N)
        DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGCB(IGC,MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGA,INDX )
        UDGC(IGC,M,NPX) = 2.D+0*DFM*(XMGCB(IGC,MP,NB)-XMGC(IGC,MP,N))
     &    /DXGF(N)
        FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
        FGCB = XGCB(IGC,MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGC = DIFMN( FGCB,FGCP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.19 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGC(IGC,M,NPX) = UG(M,NPX)*FGC
          ELSE
            UGC(IGC,M,NPX) = 0.D+0
          ENDIF
        ELSE
          UGC(IGC,M,NPX) = UG(M,NPX)*FGC + GCPP(1,IGC)*UDGC(IGC,M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCW group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWB_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWB_COMP'
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
       IF( IBCT(3,NB).EQ.19 ) THEN
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
       IF( IBCT(3,NB).EQ.19 ) THEN
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
!---  End of DFFGWB_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWE_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWE_COMP'
      NQX = NSX(N)+1
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
       IF( IBCT(3,NB).EQ.19 ) THEN
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
       IF( IBCT(3,NB).EQ.19 ) THEN
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
!---  End of DFFGWE_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWN_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWN_COMP'
      NQY = NSY(N)+IFLD
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
       IF( IBCT(3,NB).EQ.19 ) THEN
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
!---  End of DFFGWN_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWS_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWS_COMP'
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
       IF( IBCT(3,NB).EQ.19 ) THEN
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
       IF( IBCT(3,NB).EQ.19 ) THEN
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
!---  End of DFFGWS_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWT_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWT_COMP'
      NQZ = NSZ(N)+IJFLD
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
       IF( IBCT(3,NB).EQ.19 ) THEN
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
       IF( IBCT(3,NB).EQ.19 ) THEN
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
!---  End of DFFGWT_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWW_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFGWW_COMP'
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
       UDGW(M,NPX) = DFM*(XMGWB(MP,NB)-XMGW(MP,N))/DXGF(N)
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(3,NB).EQ.19 ) THEN
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
       UDGW(M,NPX) = DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB)
     &   - XMGW(MP,N)*RHOMG(MP,N))/DXGF(N)
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(3,NB).EQ.19 ) THEN
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
!---  End of DFFGWW_COMP group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLCB( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFLCB'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      K = KD(N)
      DXMLA = XMLCB(IGC,2,NB)-XMLC(IGC,2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLCB(IGC,MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLC(IGC,M,NSZ(N)) = 2.D+0*DFM*(XMLCB(IGC,MP,NB)-XMLC(IGC,MP,N))
     &  /DZGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      K = KD(N)
      DXMLA = XMLCB(IGC,2,NB)*RHOMLB(2,NB)-XMLC(IGC,2,N)*RHOML(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLCB(IGC,MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLC(IGC,M,NSZ(N)) = 2.D+0*DFM*(XMLCB(IGC,MP,NB)*RHOMLB(MP,NB)
     &   - XMLC(IGC,MP,N)*RHOML(MP,N))/DZGF(N)
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLCB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLCE( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFLCE'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      DXMLA = XMLC(IGC,2,N)-XMLCB(IGC,2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLC(IGC,MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLCB(IGC,MN,NB)
     &  *RHOMLB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLC(IGC,M,NSX(N)+1) = 2.D+0*DFM*(XMLC(IGC,MN,N)-XMLCB(IGC,MN,
     &  NB))/DXGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      DXMLA = XMLC(IGC,2,N)*RHOML(2,N)-XMLCB(IGC,2,NB)*RHOMLB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLC(IGC,MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLCB(IGC,MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLC(IGC,M,NSX(N)+1) = 2.D+0*DFM*(XMLC(IGC,MN,N)*RHOML(MN,N)
     &   - XMLCB(IGC,MN,NB)*RHOMLB(MN,NB))/DXGF(N)
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLCE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLCN( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFLCN'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      J = JD(N)
      DXMLA = XMLC(IGC,2,N)-XMLCB(IGC,2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLC(IGC,MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLCB(IGC,MN,NB)
     &   *RHOMLB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLC(IGC,M,NSY(N)+IFLD) = 2.D+0*DFM*(XMLC(IGC,MN,N)-
     &   XMLCB(IGC,MN,NB))/(DYGF(N)*RP(I))
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      J = JD(N)
      DXMLA = XMLC(IGC,2,N)*RHOML(2,N)-XMLCB(IGC,2,NB)*RHOMLB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLC(IGC,MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLCB(IGC,MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLC(IGC,M,NSY(N)+IFLD) = 2.D+0*DFM*(XMLC(IGC,MN,N)*RHOML(MN,N)
     &   - XMLCB(IGC,MN,NB)*RHOMLB(MN,NB))/(DYGF(N)*RP(I))
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLCN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLCS( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFLCS'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      J = JD(N)
      DXMLA = XMLCB(IGC,2,NB)-XMLC(IGC,2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLCB(IGC,MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLC(IGC,M,NSY(N)) = 2.D+0*DFM*(XMLCB(IGC,MP,NB)-XMLC(IGC,MP,N))/
     &   (DYGF(N)*RP(I))
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      J = JD(N)
      DXMLA = XMLCB(IGC,2,NB)*RHOMLB(2,NB)-XMLC(IGC,2,N)*RHOML(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLCB(IGC,MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLC(IGC,M,NSY(N)) = 2.D+0*DFM*(XMLCB(IGC,MP,NB)*RHOMLB(MP,NB)
     &   -XMLC(IGC,MP,N)*RHOML(MP,N))/(DYGF(N)*RP(I))
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLCS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLCT( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFLCT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      K = KD(N)
      DXMLA = XMLC(IGC,2,N)-XMLCB(IGC,2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLC(IGC,MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLCB(IGC,MN,NB)
     &   *RHOMLB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLC(IGC,M,NSZ(N)+IJFLD) = 2.D+0*DFM*(XMLC(IGC,MN,N)-
     &   XMLCB(IGC,MN,NB))/DZGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      K = KD(N)
      DXMLA = XMLC(IGC,2,N)*RHOML(2,N)-XMLCB(IGC,2,NB)*RHOMLB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLC(IGC,MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLCB(IGC,MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLC(IGC,M,NSZ(N)+IJFLD) = 2.D+0*DFM*(XMLC(IGC,MN,N)*RHOML(MN,N)
     &   - XMLCB(IGC,MN,NB)*RHOMLB(MN,NB))/DZGF(N)
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLCT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLCW( N,NB,IGC )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE BCVP
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFLCW'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      DXMLA = XMLCB(IGC,2,NB)-XMLC(IGC,2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLCB(IGC,MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLC(IGC,M,NSX(N)) = 2.D+0*DFM*(XMLCB(IGC,MP,NB)-
     &    XMLC(IGC,MP,N))/DXGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      DXMLA = XMLCB(IGC,2,NB)*RHOMLB(2,NB)-XMLC(IGC,2,N)*RHOML(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLC(IGC,MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLCB(IGC,MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLC(IGC,M,NSX(N)) = 2.D+0*DFM*(XMLCB(IGC,MP,NB)*RHOMLB(MP,NB)
     &   - XMLC(IGC,MP,N)*RHOML(MP,N))/DXGF(N)
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLCW group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSB_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSB_COMP'
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
          ULB = UBX*UBX
          VLB = VBX*VBX
          WLB = WBX*WBX
          ZLB = SQRT(ULB+VLB+WLB)
          DPLB = (DPLGS(IZN)*WLB + DPTRS(IZN)*(ULB+VLB))/(ZLB+SMALL)
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
!---  End of DFFLSB_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSS_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSS_COMP'
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
          ULS = USX*USX
          VLS = VSX*VSX
          WLS = WSX*WSX
          ZLS = SQRT(ULS+VLS+WLS)
          DPLS = (DPLGS(IZN)*VLS + DPTRS(IZN)*(ULS+WLS))/(ZLS+SMALL)
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
!---  End of DFFLSS_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSW_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSW_COMP'
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
!---  End of DFFLSW_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSE_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSE_COMP'
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
          ULE = UEX*UEX
          VLE = VEX*VEX
          WLE = WEX*WEX
          ZLE = SQRT(ULE+VLE+WLE)
          DPLE = (DPLGS(IZN)*ULE + DPTRS(IZN)*(WLE+VLE))/(ZLE+SMALL)
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
!---  End of DFFLSE_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSN_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSN_COMP'
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
!---  End of DFFLSN_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLST_COMP( N,NB )
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
!    Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      SUB_LOG(ISUB_LOG) = '/DFFLST_COMP'
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
          ULT = UTX*UTX
          VLT = VTX*VTX
          WLT = WTX*WTX
          ZLT = SQRT(ULT+VLT+WLT)
          DPLT = (DPLGS(IZN)*WLT + DPTRS(IZN)*(ULT+VLT))/(ZLT+SMALL)
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
!---  End of DFFLST_COMP group  ---
!
      RETURN
      END

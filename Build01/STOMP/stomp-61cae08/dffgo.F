!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGO
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
!     Compute diffusive oil gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, October 22, 1998.
!     dffgo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DFFGO'
!
!---  X-direction vapor mole diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 K = 1,KFLD
          DO 200 J = 1,JFLD
            DO 200 I = 2,IFLD
              N = ND(I,J,K)
              NW = N-1
              IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &          INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
              NPX = NSX(N)
              DXGO = XMGO(2,NW)-XMGO(2,N)
              DO 100 M = 1,ISVF
                MN = MNEG(M)
                MP = MPOS(M)
                DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &            DFGO(MP,N)*RHOMG(MP,N)
                DFW = TORG(MN,NW)*PORD(MN,NW)*(SG(MN,NW)-SGT(MN,NW))*
     &            DFGO(MN,NW)*RHOMG(MN,NW)
                INDX = 13
                DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXGO,INDX )
                UDGO(M,NPX) = DFM*(XMGO(MN,NW)-XMGO(MP,N))
     &            /DXGP(NPX)
                FGOP = XGO(MP,N)*RHOG(MP,N)
                FGOW = XGO(MN,NW)*RHOG(MN,NW)
                INDX = 3
                FGO = DIFMN( FGOW,FGOP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
                UGO(M,NPX) = UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
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
              IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &          INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
              NPY = NSY(N)
              DXGO = XMGO(2,NS)-XMGO(2,N)
              DO 300 M = 1,ISVF
                MP = MPOS(M)
                MN = MNEG(M)
                DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &            DFGO(MP,N)*RHOMG(MP,N)
                DFS = TORG(MN,NS)*PORD(MN,NS)*(SG(MN,NS)-SGT(MN,NS))*
     &            DFGO(MN,NS)*RHOMG(MN,NS)
                INDX = 13
                DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXGO,INDX )
                VDGO(M,NPY) = DFM*(XMGO(MN,NS)-XMGO(MP,N))
     &            /DYGP(NPY)/RP(I)
                FGOP = XGO(MP,N)*RHOG(MP,N)
                FGOS = XGO(MN,NS)*RHOG(MN,NS)
                INDX = 3
                FGO = DIFMN( FGOS,FGOP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
                VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
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
              IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &          INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
              NPZ = NSZ(N)
              DXGO = XMGO(2,NB)-XMGO(2,N)
              DO 500 M = 1,ISVF
                MP = MPOS(M)
                MN = MNEG(M)
                DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &            DFGO(MP,N)*RHOMG(MP,N)
                DFB = TORG(MN,NB)*PORD(MN,NB)*(SG(MN,NB)-SGT(MN,NB))*
     &            DFGO(MN,NB)*RHOMG(MN,NB)
                INDX = 13
                DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXGO,INDX )
                WDGO(M,NPZ) = DFM*(XMGO(MN,NB)-XMGO(MP,N))
     &            /DZGP(NPZ)
                FGOP = XGO(MP,N)*RHOG(MP,N)
                FGOB = XGO(MN,NB)*RHOG(MN,NB)
                INDX = 3
                FGO = DIFMN( FGOB,FGOP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
                WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
  500         CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGO group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNA
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
!     Diffusive CO2 nonaqueous-liquid fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
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
      SUB_LOG(ISUB_LOG) = '/DFFNA'
!
!---  X-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNA = XMNA(2,NW) - XMNA(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNA(MP,N)*RHOMN(MP,N)
            DFW = TORN(MN,NW)*PORD(MN,NW)*SN(MN,NW)*
     &        DFNA(MN,NW)*RHOMN(MN,NW)
            INDX = 20
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMNA,INDX )
            UDNA(M,NPX) = DFM*(XMNA(MN,NW)-XMNA(MP,N))/DXGP(NPX)
            FNAP = XNA(MP,N)*RHON(MP,N)
            FNAW = XNA(MN,NW)*RHON(MN,NW)
            INDX = 4
            FNA = DIFMN( FNAW,FNAP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UNA(M,NPX) = UN(M,NPX)*FNA + WTMA*UDNA(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNA = XMNA(2,NS) - XMNA(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNA(MP,N)*RHOMN(MP,N)
            DFS = TORN(MN,NS)*PORD(MN,NS)*SN(MN,NS)*
     &        DFNA(MN,NS)*RHOMN(MN,NS)
            INDX = 20
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMNA,INDX )
            VDNA(M,NPY) = DFM*(XMNA(MN,NS)-XMNA(MP,N))/(DYGP(NPY)*RP(I))
            FNAP = XNA(MP,N)*RHON(MP,N)
            FNAS = XNA(MN,NS)*RHON(MN,NS)
            INDX = 4
            FNA = DIFMN( FNAS,FNAP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VNA(M,NPY) = VN(M,NPY)*FNA + WTMA*VDNA(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNA = XMNA(2,NB) - XMNA(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNA(MP,N)*RHOMN(MP,N)
            DFB = TORN(MN,NB)*PORD(MN,NB)*SN(MN,NB)*
     &        DFNA(MN,NB)*RHOMN(MN,NB)
            INDX = 20
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMNA,INDX )
            WDNA(M,NPZ) = DFM*(XMNA(MN,NB)-XMNA(MP,N))/DZGP(NPZ)
            FNAP = XNA(MP,N)*RHON(MP,N)
            FNAB = XNA(MN,NB)*RHON(MN,NB)
            INDX = 4
            FNA = DIFMN( FNAB,FNAP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WNA(M,NPZ) = WN(M,NPZ)*FNA + WTMA*WDNA(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNA group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNN
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
!     Diffusive N2 nonaqueous-liquid fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
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
      SUB_LOG(ISUB_LOG) = '/DFFNN'
!
!---  X-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNN = XMNN(2,NW) - XMNN(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNN(MP,N)*RHOMN(MP,N)
            DFW = TORN(MN,NW)*PORD(MN,NW)*SN(MN,NW)*
     &        DFNN(MN,NW)*RHOMN(MN,NW)
            INDX = 20
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMNN,INDX )
            UDNN(M,NPX) = DFM*(XMNN(MN,NW)-XMNN(MP,N))/DXGP(NPX)
            FNNP = XNN(MP,N)*RHON(MP,N)
            FNNW = XNN(MN,NW)*RHON(MN,NW)
            INDX = 4
            FNN = DIFMN( FNNW,FNNP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UNN(M,NPX) = UN(M,NPX)*FNN + WTMN*UDNN(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNN = XMNN(2,NS) - XMNN(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNN(MP,N)*RHOMN(MP,N)
            DFS = TORN(MN,NS)*PORD(MN,NS)*SN(MN,NS)*
     &        DFNN(MN,NS)*RHOMN(MN,NS)
            INDX = 20
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMNN,INDX )
            VDNN(M,NPY) = DFM*(XMNN(MN,NS)-XMNN(MP,N))/(DYGP(NPY)*RP(I))
            FNNP = XNN(MP,N)*RHON(MP,N)
            FNNS = XNN(MN,NS)*RHON(MN,NS)
            INDX = 4
            FNN = DIFMN( FNNS,FNNP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VNN(M,NPY) = VN(M,NPY)*FNN + WTMN*VDNN(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNN = XMNN(2,NB) - XMNN(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNN(MP,N)*RHOMN(MP,N)
            DFB = TORN(MN,NB)*PORD(MN,NB)*SN(MN,NB)*
     &        DFNN(MN,NB)*RHOMN(MN,NB)
            INDX = 20
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMNN,INDX )
            WDNN(M,NPZ) = DFM*(XMNN(MN,NB)-XMNN(MP,N))/DZGP(NPZ)
            FNNP = XNN(MP,N)*RHON(MP,N)
            FNNB = XNN(MN,NB)*RHON(MN,NB)
            INDX = 4
            FNN = DIFMN( FNNB,FNNP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WNN(M,NPZ) = WN(M,NPZ)*FNN + WTMN*WDNN(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNO
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
!     Diffusive CH4 nonaqueous-liquid fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXN
      USE FLUXD
      USE FDVP
      USE FDVN
      USE FDVA
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
      SUB_LOG(ISUB_LOG) = '/DFFNO'
!
!---  X-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNO = XMNO(2,NW) - XMNO(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNO(MP,N)*RHOMN(MP,N)
            DFW = TORN(MN,NW)*PORD(MN,NW)*SN(MN,NW)*
     &        DFNO(MN,NW)*RHOMN(MN,NW)
            INDX = 20
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMNO,INDX )
            UDNO(M,NPX) = DFM*(XMNO(MN,NW)-XMNO(MP,N))/DXGP(NPX)
            FNOP = XNO(MP,N)*RHON(MP,N)
            FNOW = XNO(MN,NW)*RHON(MN,NW)
            INDX = 4
            FNO = DIFMN( FNOW,FNOP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UNO(M,NPX) = UN(M,NPX)*FNO + WTMO*UDNO(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNO = XMNO(2,NS) - XMNO(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNO(MP,N)*RHOMN(MP,N)
            DFS = TORN(MN,NS)*PORD(MN,NS)*SN(MN,NS)*
     &        DFNO(MN,NS)*RHOMN(MN,NS)
            INDX = 20
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMNO,INDX )
            VDNO(M,NPY) = DFM*(XMNO(MN,NS)-XMNO(MP,N))/(DYGP(NPY)*RP(I))
            FNOP = XNO(MP,N)*RHON(MP,N)
            FNOS = XNO(MN,NS)*RHON(MN,NS)
            INDX = 4
            FNO = DIFMN( FNOS,FNOP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VNO(M,NPY) = VN(M,NPY)*FNO + WTMO*VDNO(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNO = XMNO(2,NB) - XMNO(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNO(MP,N)*RHOMN(MP,N)
            DFB = TORN(MN,NB)*PORD(MN,NB)*SN(MN,NB)*
     &        DFNO(MN,NB)*RHOMN(MN,NB)
            INDX = 20
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMNO,INDX )
            WDNO(M,NPZ) = DFM*(XMNO(MN,NB)-XMNO(MP,N))/DZGP(NPZ)
            FNOP = XNO(MP,N)*RHON(MP,N)
            FNOB = XNO(MN,NB)*RHON(MN,NB)
            INDX = 4
            FNO = DIFMN( FNOB,FNOP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WNO(M,NPZ) = WN(M,NPZ)*FNO + WTMO*WDNO(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNO group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNW
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
!     Diffusive water nonaqueous-liquid fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
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
      SUB_LOG(ISUB_LOG) = '/DFFNW'
!
!---  X-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNW = XMNW(2,NW) - XMNW(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNW(MP,N)*RHOMN(MP,N)
            DFW = TORN(MN,NW)*PORD(MN,NW)*SN(MN,NW)*
     &        DFNW(MN,NW)*RHOMN(MN,NW)
            INDX = 20
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMNW,INDX )
            UDNW(M,NPX) = DFM*(XMNW(MN,NW)-XMNW(MP,N))/DXGP(NPX)
            FNWP = XGW(MP,N)*RHON(MP,N)
            FNWW = XGW(MN,NW)*RHON(MN,NW)
            INDX = 4
            FNW = DIFMN( FNWW,FNWP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UNW(M,NPX) = UN(M,NPX)*FNW + WTMW*UDNW(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNW = XMNW(2,NS) - XMNW(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNW(MP,N)*RHOMN(MP,N)
            DFS = TORN(MN,NS)*PORD(MN,NS)*SN(MN,NS)*
     &        DFNW(MN,NS)*RHOMN(MN,NS)
            INDX = 20
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMNW,INDX )
            VDNW(M,NPY) = DFM*(XMNW(MN,NS)-XMNW(MP,N))/(DYGP(NPY)*RP(I))
            FNWP = XGW(MP,N)*RHON(MP,N)
            FNWS = XGW(MN,NS)*RHON(MN,NS)
            INDX = 4
            FNW = DIFMN( FNWS,FNWP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VNW(M,NPY) = VN(M,NPY)*FNW + WTMW*VDNW(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNW = XMNW(2,NB) - XMNW(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNW(MP,N)*RHOMN(MP,N)
            DFB = TORN(MN,NB)*PORD(MN,NB)*SN(MN,NB)*
     &        DFNW(MN,NB)*RHOMN(MN,NB)
            INDX = 20
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMNW,INDX )
            WDNW(M,NPZ) = DFM*(XMNW(MN,NB)-XMNW(MP,N))/DZGP(NPZ)
            FNWP = XGW(MP,N)*RHON(MP,N)
            FNWB = XGW(MN,NB)*RHON(MN,NB)
            INDX = 4
            FNW = DIFMN( FNWB,FNWP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WNW(M,NPZ) = WN(M,NPZ)*FNW + WTMW*WDNW(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAB( N,NB )
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
!     Diffusive CO2 nonaqueous-liquid fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFNAB'
      K = KD(N)
      NPZ = NSZ(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNAP = XNA(MP,N)*RHON(MP,N)
        FNAB = XNAB(MP,NB)*RHONB(MP,NB)
        INDX = 4
        FNA = DIFMN( FNAB,FNAP,DZGF(N),DZGF(N),WN(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WN(1,NPZ).LT.-EPSL ) THEN
            WNA(M,NPZ) = WN(M,NPZ)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WN(1,NPZ).GT.EPSL ) THEN
            WNA(M,NPZ) = WN(M,NPZ)*FNA
          ENDIF
        ELSE
          WDNA(M,NPZ) = -(WDNN(M,NPZ)+WDNO(M,NPZ)+WDNW(M,NPZ))
          WNA(M,NPZ) = WN(M,NPZ)*FNA + WTMA*WDNA(M,NPZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAE( N,NB )
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
!     Diffusive CO2 nonaqueous-liquid fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFNAE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMNA = XMNA(2,N) - XMNAB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNAP = XNA(MN,N)*RHON(MN,N)
        FNAB = XNAB(MN,NB)*RHONB(MN,NB)
        INDX = 4
        FNA = DIFMN( FNAP,FNAB,DXGF(N),DXGF(N),UN(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UN(1,NQX).GT.EPSL ) THEN
            UNA(M,NQX) = UN(M,NQX)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UN(1,NQX).LT.-EPSL ) THEN
            UNA(M,NQX) = UN(M,NQX)*FNA
          ENDIF
        ELSE
          UDNA(M,NQX) = -(UDNN(M,NQX)+UDNO(M,NQX)+UDNW(M,NQX))
          UNA(M,NQX) = UN(M,NQX)*FNA + WTMA*UDNA(M,NQX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAN( N,NB )
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
!     Diffusive CO2 nonaqueous-liquid fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFNAN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMNA = XMNA(2,N) - XMNAB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNAP = XNA(MN,N)*RHON(MN,N)
        FNAB = XNAB(MN,NB)*RHONB(MN,NB)
        INDX = 4
        FNA = DIFMN( FNAP,FNAB,DYGF(N),DYGF(N),VN(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VN(1,NQY).GT.EPSL ) THEN
            VNA(M,NQY) = VN(M,NQY)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VN(1,NQY).LT.-EPSL ) THEN
            VNA(M,NQY) = VN(M,NQY)*FNA
          ENDIF
        ELSE
          VDNA(M,NQY) = -(VDNN(M,NQY)+VDNO(M,NQY)+VDNW(M,NQY))
          VNA(M,NQY) = VN(M,NQY)*FNA + WTMN*VDNA(M,NQY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAS( N,NB )
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
!     Diffusive CO2 nonaqueous-liquid fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFNAS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMNA = XMNAB(2,NB) - XMNA(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNAP = XNA(MP,N)*RHON(MP,N)
        FNAB = XNAB(MP,NB)*RHONB(MP,NB)
        INDX = 4
        FNA = DIFMN( FNAB,FNAP,DYGF(N),DYGF(N),VN(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VN(1,NPY).LT.-EPSL ) THEN
            VNA(M,NPY) = VN(M,NPY)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VN(1,NPY).GT.EPSL ) THEN
            VNA(M,NPY) = VN(M,NPY)*FNA
          ENDIF
        ELSE
          VDNA(M,NPY) = -(VDNN(M,NPY)+VDNO(M,NPY)+VDNW(M,NPY))
          VNA(M,NPY) = VN(M,NPY)*FNA + WTMA*VDNA(M,NPY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAT( N,NB )
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
!     Diffusive CO2 nonaqueous-liquid fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFNAT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMNA = XMNA(2,N)-XMNAB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNAP = XNA(MN,N)*RHON(MN,N)
        FNAB = XNAB(MN,NB)*RHONB(MN,NB)
        INDX = 4
        FNA = DIFMN( FNAP,FNAB,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WN(1,NQZ).GT.EPSL ) THEN
            WNA(M,NQZ) = WN(M,NQZ)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WN(1,NQZ).LT.-EPSL ) THEN
            WNA(M,NQZ) = WN(M,NQZ)*FNA
          ENDIF
        ELSE
          WDNA(M,NQZ) = -(WDNN(M,NQZ)+WDNO(M,NQZ)+WDNW(M,NQZ))
          WNA(M,NQZ) = WN(M,NQZ)*FNA + WTMA*WDNA(M,NQZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAW( N,NB )
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
!     Diffusive CO2 nonaqueous-liquid fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFNAW'
      I = ID(N)
      NPX = NSX(N)
      DXMNA = XMNAB(2,NB)-XMNA(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNAP = XNA(MP,N)*RHON(MP,N)
        FNAB = XNAB(MP,NB)*RHONB(MP,NB)
        INDX = 4
        FNA = DIFMN( FNAB,FNAP,DXGF(N),DXGF(N),UN(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UN(1,NPX).LT.-EPSL ) THEN
            UNA(M,NPX) = UN(M,NPX)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UN(1,NPX).GT.EPSL ) THEN
            UNA(M,NPX) = UN(M,NPX)*FNA
          ENDIF
        ELSE
          UDNA(M,NPX) = -(UDNN(M,NPX)+UDNO(M,NPX)+UDNW(M,NPX))
          UNA(M,NPX) = UN(M,NPX)*FNA + WTMA*UDNA(M,NPX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNNB( N,NB )
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
!     Diffusive N2 nonaqueous-liquid fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNNB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMNN = XMNNB(2,NB) - XMNN(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNNP = XNN(MP,N)*RHON(MP,N)
        FNNB = XNNB(MP,NB)*RHONB(MP,NB)
        INDX = 4
        FNN = DIFMN( FNNB,FNNP,DZGF(N),DZGF(N),WN(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WN(1,NPZ).LT.-EPSL ) THEN
            WNN(M,NPZ) = WN(M,NPZ)*FNN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WN(1,NPZ).GT.EPSL ) THEN
            WNN(M,NPZ) = WN(M,NPZ)*FNN
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNN(MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNNB(MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMNN,INDX )
          WDNN(M,NPZ) = 2.D+0*DFM*(XMNNB(MP,NB)-XMNN(MP,N))/DZGF(N)
          WNN(M,NPZ) = WN(M,NPZ)*FNN + WTMN*WDNN(M,NPZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNNB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNNE( N,NB )
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
!     Diffusive N2 nonaqueous-liquid fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNNE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMNN = XMNN(2,N) - XMNNB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNNP = XNN(MN,N)*RHON(MN,N)
        FNNB = XNNB(MN,NB)*RHONB(MN,NB)
        INDX = 4
        FNN = DIFMN( FNNP,FNNB,DXGF(N),DXGF(N),UN(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UN(1,NQX).GT.EPSL ) THEN
            UNN(M,NQX) = UN(M,NQX)*FNN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UN(1,NQX).LT.-EPSL ) THEN
            UNN(M,NQX) = UN(M,NQX)*FNN
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNN(MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNNB(MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMNN,INDX )
          UDNN(M,NQX) = 2.D+0*DFM*(XMNN(MN,N)-XMNNB(MN,NB))/DXGF(N)
          UNN(M,NQX) = UN(M,NQX)*FNN + WTMN*UDNN(M,NQX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNNE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNNN( N,NB )
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
!     Diffusive N2 nonaqueous-liquid fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNNN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMNN = XMNN(2,N) - XMNNB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNNP = XNN(MN,N)*RHON(MN,N)
        FNNB = XNNB(MN,NB)*RHONB(MN,NB)
        INDX = 4
        FNN = DIFMN( FNNP,FNNB,DYGF(N),DYGF(N),VN(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VN(1,NQY).GT.EPSL ) THEN
            VNN(M,NQY) = VN(M,NQY)*FNN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VN(1,NQY).LT.-EPSL ) THEN
            VNN(M,NQY) = VN(M,NQY)*FNN
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNN(MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNNB(MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMNN,INDX )
          VDNN(M,NQY) = 2.D+0*DFM*(XMNN(MN,N)-XMNNB(MN,NB))
     &      /(DYGF(N)*RP(I))
          VNN(M,NQY) = VN(M,NQY)*FNN + WTMN*VDNN(M,NQY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNNN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNNS( N,NB )
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
!     Diffusive N2 nonaqueous-liquid fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNNS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMNN = XMNNB(2,NB) - XMNN(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNNP = XNN(MP,N)*RHON(MP,N)
        FNNB = XNNB(MP,NB)*RHONB(MP,NB)
        INDX = 4
        FNN = DIFMN( FNNB,FNNP,DYGF(N),DYGF(N),VN(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VN(1,NPY).LT.-EPSL ) THEN
            VNN(M,NPY) = VN(M,NPY)*FNN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VN(1,NPY).GT.EPSL ) THEN
            VNN(M,NPY) = VN(M,NPY)*FNN
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNN(MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNNB(MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMNN,INDX )
          VDNN(M,NPY) = 2.D+0*DFM*(XMNNB(MP,NB)-XMNN(MP,N))
     &      /(DYGF(N)*RP(I))
          VNN(M,NPY) = VN(M,NPY)*FNN + WTMN*VDNN(M,NPY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNNS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNNT( N,NB )
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
!     Diffusive N2 nonaqueous-liquid fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNNT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMNN = XMNN(2,N)-XMNNB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNNP = XNN(MN,N)*RHON(MN,N)
        FNNB = XNNB(MN,NB)*RHONB(MN,NB)
        INDX = 4
        FNN = DIFMN( FNNP,FNNB,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WN(1,NQZ).GT.EPSL ) THEN
            WNN(M,NQZ) = WN(M,NQZ)*FNN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WN(1,NQZ).LT.-EPSL ) THEN
            WNN(M,NQZ) = WN(M,NQZ)*FNN
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNN(MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNNB(MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMNN,INDX )
          WDNN(M,NQZ) = 2.D+0*DFM*(XMNN(MN,N)-XMNNB(MN,NB))/DZGF(N)
          WNN(M,NQZ) = WN(M,NQZ)*FNN + WTMN*WDNN(M,NQZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNNT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNNW( N,NB )
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
!     Diffusive CO2 nonaqueous-liquid fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNNW'
      I = ID(N)
      NPX = NSX(N)
      DXMNN = XMNNB(2,NB)-XMNN(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNNP = XNN(MP,N)*RHON(MP,N)
        FNNB = XNNB(MP,NB)*RHONB(MP,NB)
        INDX = 4
        FNN = DIFMN( FNNB,FNNP,DXGF(N),DXGF(N),UN(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UN(1,NPX).LT.-EPSL ) THEN
            UNN(M,NPX) = UN(M,NPX)*FNN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UN(1,NPX).GT.EPSL ) THEN
            UNN(M,NPX) = UN(M,NPX)*FNN
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNN(MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNNB(MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMNN,INDX )
          UDNN(M,NPX) = 2.D+0*DFM*(XMNNB(MP,NB)-XMNN(MP,N))/DXGF(N)
          UNN(M,NPX) = UN(M,NPX)*FNN + WTMN*UDNN(M,NPX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNNW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNOB( N,NB )
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
!     Diffusive CH4 nonaqueous-liquid fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXN
      USE FLUXD
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNOB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMNO = XMNOB(2,NB) - XMNO(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNOP = XNO(MP,N)*RHON(MP,N)
        FNOB = XNOB(MP,NB)*RHONB(MP,NB)
        INDX = 4
        FNO = DIFMN( FNOB,FNOP,DZGF(N),DZGF(N),WN(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WN(1,NPZ).LT.-EPSL ) THEN
            WNO(M,NPZ) = WN(M,NPZ)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WN(1,NPZ).GT.EPSL ) THEN
            WNO(M,NPZ) = WN(M,NPZ)*FNO
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNO(MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNOB(MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMNO,INDX )
          WDNO(M,NPZ) = 2.D+0*DFM*(XMNOB(MP,NB)-XMNO(MP,N))/DZGF(N)
          WNO(M,NPZ) = WN(M,NPZ)*FNO + WTMO*WDNO(M,NPZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNOE( N,NB )
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
!     Diffusive CH4 nonaqueous-liquid fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXN
      USE FLUXD
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNOE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMNO = XMNO(2,N) - XMNOB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNOP = XNO(MN,N)*RHON(MN,N)
        FNOB = XNOB(MN,NB)*RHONB(MN,NB)
        INDX = 4
        FNO = DIFMN( FNOP,FNOB,DXGF(N),DXGF(N),UN(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UN(1,NQX).GT.EPSL ) THEN
            UNO(M,NQX) = UN(M,NQX)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UN(1,NQX).LT.-EPSL ) THEN
            UNO(M,NQX) = UN(M,NQX)*FNO
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNO(MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNOB(MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMNO,INDX )
          UDNO(M,NQX) = 2.D+0*DFM*(XMNO(MN,N)-XMNOB(MN,NB))/DXGF(N)
          UNO(M,NQX) = UN(M,NQX)*FNO + WTMO*UDNO(M,NQX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNON( N,NB )
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
!     Diffusive CH4 nonaqueous-liquid fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXN
      USE FLUXD
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNON'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMNO = XMNO(2,N) - XMNOB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNOP = XNO(MN,N)*RHON(MN,N)
        FNOB = XNOB(MN,NB)*RHONB(MN,NB)
        INDX = 4
        FNO = DIFMN( FNOP,FNOB,DYGF(N),DYGF(N),VN(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VN(1,NQY).GT.EPSL ) THEN
            VNO(M,NQY) = VN(M,NQY)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VN(1,NQY).LT.-EPSL ) THEN
            VNO(M,NQY) = VN(M,NQY)*FNO
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNO(MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNOB(MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMNO,INDX )
          VDNO(M,NQY) = 2.D+0*DFM*(XMNO(MN,N)-XMNOB(MN,NB))
     &      /(DYGF(N)*RP(I))
          VNO(M,NQY) = VN(M,NQY)*FNO + WTMO*VDNO(M,NQY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNOS( N,NB )
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
!     Diffusive CH4 nonaqueous-liquid fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXN
      USE FLUXD
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNOS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMNO = XMNOB(2,NB) - XMNO(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNOP = XNO(MP,N)*RHON(MP,N)
        FNOB = XNOB(MP,NB)*RHONB(MP,NB)
        INDX = 4
        FNO = DIFMN( FNOB,FNOP,DYGF(N),DYGF(N),VN(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VN(1,NPY).LT.-EPSL ) THEN
            VNO(M,NPY) = VN(M,NPY)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VN(1,NPY).GT.EPSL ) THEN
            VNO(M,NPY) = VN(M,NPY)*FNO
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNO(MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNOB(MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMNO,INDX )
          VDNO(M,NPY) = 2.D+0*DFM*(XMNOB(MP,NB)-XMNO(MP,N))
     &      /(DYGF(N)*RP(I))
          VNO(M,NPY) = VN(M,NPY)*FNO + WTMO*VDNO(M,NPY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNOT( N,NB )
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
!     Diffusive CH4 nonaqueous-liquid fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXN
      USE FLUXD
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNOT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMNO = XMNO(2,N)-XMNOB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNOP = XNO(MN,N)*RHON(MN,N)
        FNOB = XNOB(MN,NB)*RHONB(MN,NB)
        INDX = 4
        FNO = DIFMN( FNOP,FNOB,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WN(1,NQZ).GT.EPSL ) THEN
            WNO(M,NQZ) = WN(M,NQZ)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WN(1,NQZ).LT.-EPSL ) THEN
            WNO(M,NQZ) = WN(M,NQZ)*FNO
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNO(MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNOB(MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMNO,INDX )
          WDNO(M,NQZ) = 2.D+0*DFM*(XMNO(MN,N)-XMNOB(MN,NB))/DZGF(N)
          WNO(M,NQZ) = WN(M,NQZ)*FNO + WTMO*WDNO(M,NQZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNOW( N,NB )
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
!     Diffusive CH4 nonaqueous-liquid fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXN
      USE FLUXD
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNOW'
      I = ID(N)
      NPX = NSX(N)
      DXMNO = XMNOB(2,NB)-XMNO(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNOP = XNO(MP,N)*RHON(MP,N)
        FNOB = XNOB(MP,NB)*RHONB(MP,NB)
        INDX = 4
        FNO = DIFMN( FNOB,FNOP,DXGF(N),DXGF(N),UN(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UN(1,NPX).LT.-EPSL ) THEN
            UNO(M,NPX) = UN(M,NPX)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UN(1,NPX).GT.EPSL ) THEN
            UNO(M,NPX) = UN(M,NPX)*FNO
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNO(MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNOB(MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMNO,INDX )
          UDNO(M,NPX) = 2.D+0*DFM*(XMNOB(MP,NB)-XMNO(MP,N))/DXGF(N)
          UNO(M,NPX) = UN(M,NPX)*FNO + WTMO*UDNO(M,NPX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNOW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWB( N,NB )
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
!     Diffusive water nonaqueous-liquid fluxes on bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNWB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMNW = XMNWB(2,NB) - XMNW(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNWB = XGWB(MP,NB)*RHONB(MP,NB)
        FNWP = XGW(MP,N)*RHON(MP,N)
        INDX = 4
        FNW = DIFMN( FNWB,FNWP,DZGF(N),DZGF(N),WN(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WN(1,NPZ).LT.-EPSL ) THEN
            WNW(M,NPZ) = WN(M,NPZ)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WN(1,NPZ).GT.EPSL ) THEN
            WNW(M,NPZ) = WN(M,NPZ)*FNW
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNW(MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNWB(MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMNW,INDX )
          WDNW(M,NPZ) = 2.D+0*DFM*(XMNWB(MP,NB)-XMNW(MP,N))/DZGF(N)
          WNW(M,NPZ) = WN(M,NPZ)*FNW + WTMW*WDNW(M,NPZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWE( N,NB )
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
!     Diffusive water nonaqueous-liquid fluxes on east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNWE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMNW = XMNW(2,N) - XMNWB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNWB = XGWB(MN,NB)*RHONB(MN,NB)
        FNWP = XGW(MN,N)*RHON(MN,N)
        INDX = 4
        FNW = DIFMN( FNWP,FNWB,DXGF(N),DXGF(N),UN(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UN(1,NQX).GT.EPSL ) THEN
            UNW(M,NQX) = UN(M,NQX)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UN(1,NQX).LT.-EPSL ) THEN
            UNW(M,NQX) = UN(M,NQX)*FNW
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNW(MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNWB(MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMNW,INDX )
          UDNW(M,NQX) = 2.D+0*DFM*(XMNW(MN,N)-XMNWB(MN,NB))/DXGF(N)
          UNW(M,NQX) = UN(M,NQX)*FNW + WTMW*UDNW(M,NQX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWN( N,NB )
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
!     Diffusive water nonaqueous-liquid fluxes on north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNWN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMNW = XMNW(2,N) - XMNWB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNWB = XGWB(MN,NB)*RHONB(MN,NB)
        FNWP = XGW(MN,N)*RHON(MN,N)
        INDX = 4
        FNW = DIFMN( FNWP,FNWB,DYGF(N),DYGF(N),VN(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VN(1,NQY).GT.EPSL ) THEN
            VNW(M,NQY) = VN(M,NQY)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VN(1,NQY).LT.-EPSL ) THEN
            VNW(M,NQY) = VN(M,NQY)*FNW
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNW(MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNWB(MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMNW,INDX )
          VDNW(M,NQY) = 2.D+0*DFM*(XMNW(MN,N)-XMNWB(MN,NB))
     &      /(DYGF(N)*RP(I))
          VNW(M,NQY) = VN(M,NQY)*FNW + WTMW*VDNW(M,NQY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWS( N,NB )
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
!     Diffusive water nonaqueous-liquid fluxes on south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNWS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMNW = XMNWB(2,NB) - XMNW(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNWB = XGWB(MP,NB)*RHONB(MP,NB)
        FNWP = XGW(MP,N)*RHON(MP,N)
        INDX = 4
        FNW = DIFMN( FNWB,FNWP,DYGF(N),DYGF(N),VN(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VN(1,NPY).LT.-EPSL ) THEN
            VNW(M,NPY) = VN(M,NPY)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VN(1,NPY).GT.EPSL ) THEN
            VNW(M,NPY) = VN(M,NPY)*FNW
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNW(MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNWB(MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMNW,INDX )
          VDNW(M,NPY) = 2.D+0*DFM*(XMNWB(MP,NB)-XMNW(MP,N))
     &      /(DYGF(N)*RP(I))
          VNW(M,NPY) = VN(M,NPY)*FNW + WTMW*VDNW(M,NPY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWT( N,NB )
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
!     Diffusive water nonaqueous-liquid fluxes on top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNWT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMNW = XMNW(2,N)-XMNWB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FNWB = XGWB(MN,NB)*RHONB(MN,NB)
        FNWP = XGW(MN,N)*RHON(MN,N)
        INDX = 4
        FNW = DIFMN( FNWP,FNWB,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WN(1,NQZ).GT.EPSL ) THEN
            WNW(M,NQZ) = WN(M,NQZ)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WN(1,NQZ).LT.-EPSL ) THEN
            WNW(M,NQZ) = WN(M,NQZ)*FNW
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNW(MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNWB(MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMNW,INDX )
          WDNW(M,NQZ) = 2.D+0*DFM*(XMNW(MN,N)-XMNWB(MN,NB))/DZGF(N)
          WNW(M,NQZ) = WN(M,NQZ)*FNW + WTMW*WDNW(M,NQZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWW( N,NB )
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
!     Diffusive water nonaqueous-liquid fluxes on west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 October 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNWW'
      I = ID(N)
      NPX = NSX(N)
      DXMNW = XMNWB(2,NB)-XMNW(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FNWB = XGWB(MP,NB)*RHONB(MP,NB)
        FNWP = XGW(MP,N)*RHON(MP,N)
        INDX = 4
        FNW = DIFMN( FNWB,FNWP,DXGF(N),DXGF(N),UN(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UN(1,NPX).LT.-EPSL ) THEN
            UNW(M,NPX) = UN(M,NPX)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UN(1,NPX).GT.EPSL ) THEN
            UNW(M,NPX) = UN(M,NPX)*FNW
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNW(MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNWB(MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMNW,INDX )
          UDNW(M,NPX) = 2.D+0*DFM*(XMNWB(MP,NB)-XMNW(MP,N))/DXGF(N)
          UNW(M,NPX) = UN(M,NPX)*FNW + WTMW*UDNW(M,NPX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWW group  ---
!
      RETURN
      END


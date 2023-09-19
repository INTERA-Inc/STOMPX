!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDG
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, June 30, 1994.
!     thdg.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/THDG'
!
!---  X-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
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
              MP = MPOS(M)
              MN = MNEG(M)
              UDGAX = -UDGW(M,NPX)-UDGO(M,NPX)
              UQ(M,NPX) = UQ(M,NPX) 
     &          + HGW(MN,NW)*MAX(UDGW(M,NPX),ZERO)*WTMW
     &          - HGW(MP,N)*MAX(-UDGW(M,NPX),ZERO)*WTMW
     &          + HGA(MN,NW)*MAX(UDGAX,ZERO)*WTMA
     &          - HGA(MP,N)*MAX(-UDGAX,ZERO)*WTMA
     &          + HGO(MN,NW)*MAX(UDGO(M,NPX),ZERO)*WTMO
     &          - HGO(MP,N)*MAX(-UDGO(M,NPX),ZERO)*WTMO
  100       CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
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
              MP = MPOS(M)
              MN = MNEG(M)
              VDGAX = -VDGW(M,NPY)-VDGO(M,NPY)
              VQ(M,NPY) = VQ(M,NPY) 
     &          + HGW(MN,NS)*MAX(VDGW(M,NPY),ZERO)*WTMW
     &          - HGW(MP,N)*MAX(-VDGW(M,NPY),ZERO)*WTMW
     &          + HGA(MN,NS)*MAX(VDGAX,ZERO)*WTMA
     &          - HGA(MP,N)*MAX(-VDGAX,ZERO)*WTMA
     &          + HGO(MN,NS)*MAX(VDGO(M,NPY),ZERO)*WTMO
     &          - HGO(MP,N)*MAX(-VDGO(M,NPY),ZERO)*WTMO
  300       CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
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
              MP = MPOS(M)
              MN = MNEG(M)
              WDGAX = -WDGW(M,NPZ)-WDGO(M,NPZ)
              WQ(M,NPZ) = WQ(M,NPZ) 
     &          + HGW(MN,NB)*MAX(WDGW(M,NPZ),ZERO)*WTMW
     &          - HGW(MP,N)*MAX(-WDGW(M,NPZ),ZERO)*WTMW
     &          + HGA(MN,NB)*MAX(WDGAX,ZERO)*WTMA
     &          - HGA(MP,N)*MAX(-WDGAX,ZERO)*WTMA
     &          + HGO(MN,NB)*MAX(WDGO(M,NPZ),ZERO)*WTMO
     &          - HGO(MP,N)*MAX(-WDGO(M,NPZ),ZERO)*WTMO
  500       CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDG group  ---
!
      RETURN
      END

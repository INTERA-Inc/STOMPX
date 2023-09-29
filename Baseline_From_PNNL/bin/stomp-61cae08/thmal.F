!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMAL
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, April 8, 1994.
!     thmal.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/THMAL'
!
!---  X-direction aqueous advection, excluding boundaries
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
              HP = HL(MP,N)*RHOL(MP,N)
              HW = HL(MN,NW)*RHOL(MN,NW)
              UQ(M,NPX) = UQ(M,NPX) + HW*MAX(UL(M,NPX),ZERO)
     &          -HP*MAX(-UL(M,NPX),ZERO)
  100       CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction aqueous advection, excluding boundaries
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
              HP = HL(MP,N)*RHOL(MP,N)
              HS = HL(MN,NS)*RHOL(MN,NS)
              VQ(M,NPY) = VQ(M,NPY) + HS*MAX(VL(M,NPY),ZERO)
     &          -HP*MAX(-VL(M,NPY),ZERO)
  300       CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction aqueous advection, excluding boundaries
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
              HP = HL(MP,N)*RHOL(MP,N)
              HB = HL(MN,NB)*RHOL(MN,NB)
              WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WL(M,NPZ),ZERO)
     &          -HP*MAX(-WL(M,NPZ),ZERO)
  500       CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMAL group  ---
!
      RETURN
      END

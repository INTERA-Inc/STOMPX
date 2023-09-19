!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVN
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
!     Compute the NAPL Darcy flux from pressure gradients
!     and gravitational body forces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, January, 1992.
!     Last Modified by MD White, PNNL, December 13, 1995.
!     Last Modified by MD White, PNNL, 29 May 2002.
!     drcvn.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 KNM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVN'
!
!---  X-direction Darcy velocities, excluding boundaries
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
              MN = MNEG(M)
              MP = MPOS(M)
              HDNX = PN(MN,NW)-PN(MP,N) - 0.5D+0*GRVX(NPX)*
     &          (RHON(MN,NW)*DXGF(N)+RHON(MP,N)*DXGF(NW))
              IF( M.EQ.1 ) HDN = HDNX
!
!---          Permeability reduction factor
!
              PERM_WX = PERMRF(MN,NW)*PERM(1,IZ(NW))
              PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
              INDX = 11
              KNM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDN,INDX)
              IF( PERM_WX/EPSL.LT.EPSL ) KNM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KNM = 0.D+0
              INDX = 10
              RKNM = DIFMN(RKN(MN,NW),RKN(MP,N),DXGF(NW),DXGF(N),
     &          HDN,INDX)
              INDX = 7
              VNM = DIFMN(VISN(MN,NW),VISN(MP,N),DXGF(NW),DXGF(N),
     &          HDN,INDX)
              UN(M,NPX) = KNM*RKNM*HDNX/DXGP(NPX)/VNM
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
            IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &        INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
            NPY = NSY(N)
            DO 300 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              HDNY = PN(MN,NS)-PN(MP,N) - 0.5D+0*GRVY(NPY)*
     &          (RHON(MN,NS)*DYGF(N)+RHON(MP,N)*DYGF(NS))
              IF( M.EQ.1 ) HDN = HDNY
!
!---          Permeability reduction factor
!
              PERM_SX = PERMRF(MN,NS)*PERM(2,IZ(NS))
              PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
              INDX = 11
              KNM = DIFMN( PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDN,INDX)
              IF( PERM_SX/EPSL.LT.EPSL ) KNM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KNM = 0.D+0
              INDX = 10
              RKNM = DIFMN(RKN(MN,NS),RKN(MP,N),DYGF(NS),DYGF(N),
     &          HDN,INDX)
              INDX = 7
              VNM = DIFMN(VISN(MN,NS),VISN(MP,N),DYGF(NS),DYGF(N),
     &          HDN,INDX)
              VN(M,NPY) = KNM*RKNM*HDNY/DYGP(NPY)/VNM/RP(I)
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
            IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &        INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
            NPZ = NSZ(N)
            DO 500 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              HDNZ = PN(MN,NB)-PN(MP,N) - 0.5D+0*GRVZ(NPZ)*
     &          (RHON(MN,NB)*DZGF(N)+RHON(MP,N)*DZGF(NB))
              IF( M.EQ.1 ) HDN = HDNZ
!
!---          Permeability reduction factor
!
              PERM_BX = PERMRF(MN,NB)*PERM(3,IZ(NB))
              PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
              INDX = 11
              KNM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDN,INDX)
              IF( PERM_BX/EPSL.LT.EPSL ) KNM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KNM = 0.D+0
              INDX = 10
              RKNM = DIFMN(RKN(MN,NB),RKN(MP,N),DZGF(NB),DZGF(N),
     &          HDN,INDX)
              INDX = 7
              VNM = DIFMN(VISN(MN,NB),VISN(MP,N),DZGF(NB),DZGF(N),
     &          HDN,INDX)
              WN(M,NPZ) = KNM*RKNM*HDNZ/DZGP(NPZ)/VNM
  500       CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVN group  ---
!
      RETURN
      END

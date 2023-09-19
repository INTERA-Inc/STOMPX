!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVG
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
!     drcvg.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, January, 1992.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      SUB_LOG(ISUB_LOG) = '/DRCVG'
!
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
              HDGX = PG(MN,NW)-PG(MP,N) - 0.5D+0*GRVX(NPX)*
     &          (RHOG(MN,NW)*DXGF(N)+RHOG(MP,N)*DXGF(NW))
              IF( M.EQ.1 ) HDG = HDGX
!
!---          Permeability reduction factor
!
              PERM_WX = PERMRF(MN,NW)*PERM(1,IZ(NW))
              PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
              INDX = 11
              KGM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDG,INDX)
              IF( PERM_WX/EPSL.LT.EPSL ) KGM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
              INDX = 9
              RKGM = DIFMN(RKG(MN,NW),RKG(MP,N),DXGF(NW),DXGF(N),
     &          HDG,INDX)
              INDX = 6
              VGM = DIFMN(VISG(MN,NW),VISG(MP,N),DXGF(NW),DXGF(N),
     &          HDG,INDX)
              UG(M,NPX) = KGM*RKGM*HDGX/DXGP(NPX)/VGM
  100       CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      IF( IJFLD.GT.1 ) THEN
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
              HDGY = PG(MN,NS)-PG(MP,N) - 0.5D+0*GRVY(NPY)*
     &          (RHOG(MN,NS)*DYGF(N)+RHOG(MP,N)*DYGF(NS))
              IF( M.EQ.1 ) HDG = HDGY
!
!---          Permeability reduction factor
!
              PERM_SX = PERMRF(MN,NS)*PERM(2,IZ(NS))
              PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
              INDX = 11
              KGM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDG,INDX)
              IF( PERM_SX/EPSL.LT.EPSL ) KGM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
              INDX = 9
              RKGM = DIFMN(RKG(MN,NS),RKG(MP,N),DYGF(NS),DYGF(N),
     &          HDG,INDX)
              INDX = 6
              VGM = DIFMN(VISG(MN,NS),VISG(MP,N),DYGF(NS),DYGF(N),
     &          HDG,INDX)
              VG(M,NPY) = KGM*RKGM*HDGY/DYGP(NPY)/VGM/RP(I)
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
              HDGZ = PG(MN,NB)-PG(MP,N) - 0.5D+0*GRVZ(NPZ)*
     &          (RHOG(MN,NB)*DZGF(N)+RHOG(MP,N)*DZGF(NB))
              IF( M.EQ.1 ) HDG = HDGZ
!
!---          Permeability reduction factor
!
              PERM_BX = PERMRF(MN,NB)*PERM(3,IZ(NB))
              PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
              INDX = 11
              KGM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDG,INDX)
              IF( PERM_BX/EPSL.LT.EPSL ) KGM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
              INDX = 9
              RKGM = DIFMN(RKG(MN,NB),RKG(MP,N),DZGF(NB),DZGF(N),
     &          HDG,INDX)
              INDX = 6
              VGM = DIFMN(VISG(MN,NB),VISG(MP,N),DZGF(NB),DZGF(N),
     &          HDG,INDX)
              WG(M,NPZ) = KGM*RKGM*HDGZ/DZGP(NPZ)/VGM
  500       CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVG group  ---
!
      RETURN
      END

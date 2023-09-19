!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVL
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
!     Written by MD White, PNNL, January, 1992.
!     Last Modified by MD White, PNNL, December 13, 1995.
!     Last Modified by MD White, PNNL, 16 January 2002.
!     Last Modified by MD White, PNNL, 29 May 2002.
!     drcvl.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE FDVS
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
      SUB_LOG(ISUB_LOG) = '/DRCVL'
!
!---  X-direction Darcy velocities, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) GOTO 200
          J = JD(N)
          K = KD(N)
            NW = N-1
            IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &        INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
            NPX = NSX(N)
            OECX = 0.D+0
            DO 100 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              INDX = 11
              IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &          OECX = DIFMN(OEC(MN,NW),OEC(MP,N),
     &          DXGF(NW),DXGF(N),ZERO,INDX)
              DPOS = OECX*(POSM(MN,NW)-POSM(MP,N))
              HDLX = PL(MN,NW)-PL(MP,N)+DPOS - 0.5D+0*GRVX(NPX)*
     &          (RHOL(MN,NW)*DXGF(N)+RHOL(MP,N)*DXGF(NW))
              IF( M.EQ.1 ) HDL = HDLX
!
!---          Permeability reduction factor
!
              PERM_WX = PERMRF(MN,NW)*PERM(1,IZ(NW))
              PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
              INDX = 11
              KLM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDL,INDX)
              INDX = 8
              RKLM = DIFMN(RKL(1,MN,NW),RKL(1,MP,N),DXGF(NW),
     &          DXGF(N),HDL,INDX)
              IF( PERM_WX/EPSL.LT.EPSL ) KLM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
              INDX = 5
              VLM = DIFMN(VISL(MN,NW),VISL(MP,N),DXGF(NW),DXGF(N),
     &          HDL,INDX)
              UL(M,NPX) = KLM*RKLM*HDLX/DXGP(NPX)/VLM
  100       CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 400 N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 400
          K = KD(N)
            NS = N-IFLD
            IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &        INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
            NPY = NSY(N)
            OECY = 0.D+0
            DO 300 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              INDX = 11
              IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &          OECY = DIFMN(OEC(MN,NS),OEC(MP,N),
     &          DYGF(NS),DYGF(N),ZERO,INDX)
              DPOS = OECY*(POSM(MN,NS)-POSM(MP,N))
              HDLY = PL(MN,NS)-PL(MP,N)+DPOS - 0.5D+0*GRVY(NPY)*
     &          (RHOL(MN,NS)*DYGF(N)+RHOL(MP,N)*DYGF(NS))
              IF( M.EQ.1 ) HDL = HDLY
!
!---          Permeability reduction factor
!
              PERM_SX = PERMRF(MN,NS)*PERM(2,IZ(NS))
              PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
              INDX = 11
              KLM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDL,INDX)
              INDX = 8
              RKLM = DIFMN(RKL(2,MN,NS),RKL(2,MP,N),DYGF(NS),
     &          DYGF(N),HDL,INDX)
              IF( PERM_SX/EPSL.LT.EPSL ) KLM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
              INDX = 5
              VLM = DIFMN(VISL(MN,NS),VISL(MP,N),DYGF(NS),DYGF(N),
     &          HDL,INDX)
              VL(M,NPY) = KLM*RKLM*HDLY/DYGP(NPY)/VLM/RP(I)
  300       CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction Darcy velocities, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 600 N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 600
            NB = N-IJFLD
            IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &        INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
            NPZ = NSZ(N)
            OECZ = 0.D+0
            DO 500 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              INDX = 11
              IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &          OECZ = DIFMN(OEC(MN,NB),OEC(MP,N),
     &          DZGF(NB),DZGF(N),ZERO,INDX)
              DPOS = OECZ*(POSM(MN,NB)-POSM(MP,N))
              HDLZ = PL(MN,NB)-PL(MP,N)+DPOS - 0.5D+0*GRVZ(NPZ)*
     &          (RHOL(MN,NB)*DZGF(N)+RHOL(MP,N)*DZGF(NB))
              IF( M.EQ.1 ) HDL = HDLZ
!
!---          Permeability reduction factor
!
              PERM_BX = PERMRF(MN,NB)*PERM(3,IZ(NB))
              PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
              INDX = 11
              KLM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDL,INDX)
              INDX = 8
              RKLM = DIFMN(RKL(3,MN,NB),RKL(3,MP,N),DZGF(NB),
     &          DZGF(N),HDL,INDX)
              IF( PERM_BX/EPSL.LT.EPSL ) KLM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
              INDX = 5
              VLM = DIFMN(VISL(MN,NB),VISL(MP,N),DZGF(NB),DZGF(N),
     &          HDL,INDX)
              WL(M,NPZ) = KLM*RKLM*HDLZ/DZGP(NPZ)/VLM
  500       CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVL group  ---
!
      RETURN
      END


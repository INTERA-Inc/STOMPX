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
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NW).LT.0 ) MN = 2
            HDGX = PG(MN,NW)-PG(MP,N)-0.5D+0*GRVX(NPX)
     &        *(RHOG(MN,NW)*DXGF(N)+RHOG(MP,N)*DXGF(NW))
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
            UG(M,NPX) = KGM*RKGM*HDGX/(DXGP(NPX)*VGM)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      IF( IJFLD.GT.1 ) THEN
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
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            HDGY = PG(MN,NS)-PG(MP,N)-0.5D+0*GRVY(NPY)
     &        *(RHOG(MN,NS)*DYGF(N)+RHOG(MP,N)*DYGF(NS))
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
            VG(M,NPY) = KGM*RKGM*HDGY/(DYGP(NPY)*VGM*RP(I))
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
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            HDGZ = PG(MN,NB)-PG(MP,N)-0.5D+0*GRVZ(NPZ)
     &        *(RHOG(MN,NB)*DZGF(N)+RHOG(MP,N)*DZGF(NB))
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
            WG(M,NPZ) = KGM*RKGM*HDGZ/(DZGP(NPZ)*VGM)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVG group  ---
!
      RETURN
      END

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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVL'
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
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NW).LT.0 ) MN = 2
            HDLX = PL(MN,NW)-PL(MP,N)-0.5D+0*GRVX(NPX)
     &        *(RHOL(MN,NW)*DXGF(N)+RHOL(MP,N)*DXGF(NW))
            IF( M.EQ.1 ) HDL = HDLX
!
!---        Permeability reduction factor
!
            PERM_WX = PERMRF(MN,NW)*PERM(1,IZ(NW))
            PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
            INDX = 11
            KLM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDL,INDX)
            IF( PERM_WX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 8
            RKLM = DIFMN(RKL(1,MN,NW),RKL(1,MP,N),DXGF(NW),
     &          DXGF(N),HDL,INDX)
            INDX = 5
            VISL_W = VISL(MN,NW)
            VISL_P = VISL(MP,N)
            VLM = DIFMN(VISL_W,VISL_P,DXGF(NW),DXGF(N),HDL,INDX)
            UL(M,NPX) = KLM*RKLM*HDLX/(DXGP(NPX)*VLM)
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
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            HDLY = PL(MN,NS)-PL(MP,N)-0.5D+0*GRVY(NPY)
     &        *(RHOL(MN,NS)*DYGF(N)+RHOL(MP,N)*DYGF(NS))
            IF( M.EQ.1 ) HDL = HDLY
!
!---        Permeability reduction factor
!
            PERM_SX = PERMRF(MN,NS)*PERM(2,IZ(NS))
            PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
            INDX = 11
            KLM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDL,INDX)
            IF( PERM_SX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 8
            RKLM = DIFMN(RKL(2,MN,NS),RKL(2,MP,N),DYGF(NS),
     &        DYGF(N),HDL,INDX)
            INDX = 5
            VISL_S = VISL(MN,NS)
            VISL_P = VISL(MP,N)
            VLM = DIFMN(VISL_S,VISL_P,DYGF(NS),DYGF(N),HDL,INDX)
            VL(M,NPY) = KLM*RKLM*HDLY/(DYGP(NPY)*VLM*RP(I))
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
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            HDLZ = PL(MN,NB)-PL(MP,N)-0.5D+0*GRVZ(NPZ)
     &        *(RHOL(MN,NB)*DZGF(N)+RHOL(MP,N)*DZGF(NB))
            IF( M.EQ.1 ) HDL = HDLZ
!
!---        Permeability reduction factor
!
            PERM_BX = PERMRF(MN,NB)*PERM(3,IZ(NB))
            PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
            INDX = 11
            KLM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDL,INDX)
            IF( PERM_BX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 8
            RKLM = DIFMN(RKL(3,MN,NB),RKL(3,MP,N),DZGF(NB),
     &        DZGF(N),HDL,INDX)
            INDX = 5
            VISL_B = VISL(MN,NB)
            VISL_P = VISL(MP,N)
            VLM = DIFMN(VISL_B,VISL_P,DZGF(NB),DZGF(N),HDL,INDX)
            WL(M,NPZ) = KLM*RKLM*HDLZ/(DZGP(NPZ)*VLM)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVL group  ---
!
      RETURN
      END

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
!     Written by MD White, PNNL, 10 May 2005.
!     Last modified by MD White, PNNL, 10 May 2005.
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
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NW).LT.0 ) MN = 2
            HDNX = PN(MN,NW)-PN(MP,N)-0.5D+0*GRVX(NPX)
     &        *(RHON(MN,NW)*DXGF(N)+RHON(MP,N)*DXGF(NW))
            IF( M.EQ.1 ) HDN = HDNX
!
!---        Permeability reduction factor
!
            PERM_WX = PERMRF(MN,NW)*PERM(1,IZ(NW))
            PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
            INDX = 11
            KNM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDN,INDX)
            IF( PERM_WX/EPSL.LT.EPSL ) KNM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KNM = 0.D+0
            INDX = 10
            RKNM = DIFMN(RKN(MN,NW),RKN(MP,N),DXGF(NW),DXGF(N),
     &        HDN,INDX)
            INDX = 7
            VNM = DIFMN(VISN(MN,NW),VISN(MP,N),DXGF(NW),DXGF(N),
     &        HDN,INDX)
            UN(M,NPX) = KNM*RKNM*HDNX/(DXGP(NPX)*VNM)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      IF( IJFLD.GT.1 ) THEN
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
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            HDNY = PN(MN,NS)-PN(MP,N)-0.5D+0*GRVY(NPY)
     &        *(RHON(MN,NS)*DYGF(N)+RHON(MP,N)*DYGF(NS))
            IF( M.EQ.1 ) HDN = HDNY
!
!---        Permeability reduction factor
!
            PERM_SX = PERMRF(MN,NS)*PERM(2,IZ(NS))
            PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
            INDX = 11
            KNM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDN,INDX)
            IF( PERM_SX/EPSL.LT.EPSL ) KNM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KNM = 0.D+0
            INDX = 10
            RKNM = DIFMN(RKN(MN,NS),RKN(MP,N),DYGF(NS),DYGF(N),
     &        HDN,INDX)
            INDX = 7
            VNM = DIFMN(VISN(MN,NS),VISN(MP,N),DYGF(NS),DYGF(N),
     &        HDN,INDX)
            VN(M,NPY) = KNM*RKNM*HDNY/(DYGP(NPY)*VNM*RP(I))
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
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            HDNZ = PN(MN,NB)-PN(MP,N)-0.5D+0*GRVZ(NPZ)
     &        *(RHON(MN,NB)*DZGF(N)+RHON(MP,N)*DZGF(NB))
            IF( M.EQ.1 ) HDN = HDNZ
!
!---        Permeability reduction factor
!
            PERM_BX = PERMRF(MN,NB)*PERM(3,IZ(NB))
            PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
            INDX = 11
            KNM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDN,INDX)
            IF( PERM_BX/EPSL.LT.EPSL ) KNM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KNM = 0.D+0
            INDX = 10
            RKNM = DIFMN(RKN(MN,NB),RKN(MP,N),DZGF(NB),DZGF(N),
     &        HDN,INDX)
            INDX = 7
            VNM = DIFMN(VISN(MN,NB),VISN(MP,N),DZGF(NB),DZGF(N),
     &        HDN,INDX)
            WN(M,NPZ) = KNM*RKNM*HDNZ/(DZGP(NPZ)*VNM)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGB( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGB'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPZ = NSZ(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HZ = PGB(MP,NB)-PG(MP,N)-5.D-1*GRVZ(NPZ)
     &    *DZGF(N)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
        WG(M,NPZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          WG(M,NPZ) = MIN( 0.D+0,WG(M,NPZ) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          WG(M,NPZ) = MAX( 0.D+0,WG(M,NPZ) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGS( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGS'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPY = NSY(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HY = PGB(MP,NB)-PG(MP,N)-5.D-1*GRVY(NPY)
     &    *DYGF(N)*RP(I)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
        VG(M,NPY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Dirchlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          VG(M,NPY) = MIN( 0.D+0,VG(M,NPY) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          VG(M,NPY) = MAX( 0.D+0,VG(M,NPY) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGW( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGW'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPX = NSX(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HX = PGB(MP,NB)-PG(MP,N)-5.D-1*GRVX(NPX)
     &    *DXGF(N)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
        UG(M,NPX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          UG(M,NPX) = MIN( 0.D+0,UG(M,NPX) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          UG(M,NPX) = MAX( 0.D+0,UG(M,NPX) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGE( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGE'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQX = NSX(N)+1
      DO M = 1,ISVF
        MN = MNEG(M)
        HX = PG(MN,N)-PGB(MN,NB)-5.D-1*GRVX(NQX)
     &    *DXGF(N)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(1,IZ(N))
        UG(M,NQX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          UG(M,NQX) = MAX( 0.D+0,UG(M,NQX) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          UG(M,NQX) = MIN( 0.D+0,UG(M,NQX) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGN( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGN'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQY = NSY(N)+IFLD
      DO M = 1,ISVF
        MN = MNEG(M)
        HY = PG(MN,N)-PGB(MN,NB)-5.D-1*GRVY(NQY)
     &    *DYGF(N)*RP(I)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(2,IZ(N))
        VG(M,NQY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          VG(M,NQY) = MAX( 0.D+0,VG(M,NQY))
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          VG(M,NQY) = MIN( 0.D+0,VG(M,NQY))
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGT( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DO M = 1,ISVF
        MN = MNEG(M)
        HZ = PG(MN,N)-PGB(MN,NB)-5.D-1*GRVZ(NQZ)
     &    *DZGF(N)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(3,IZ(N))
        WG(M,NQZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          WG(M,NQZ) = MAX( 0.D+0,WG(M,NQZ) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          WG(M,NQZ) = MIN( 0.D+0,WG(M,NQZ) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLB( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLB'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPZ = NSZ(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HZ = PLB(MP,NB)-PL(MP,N)-5.D-1*GRVZ(NPZ)
     &    *DZGF(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKLB(3,MP,NB),RKL(3,MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
        WL(M,NPZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          WL(M,NPZ) = MIN( 0.D+0,WL(M,NPZ) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          WL(M,NPZ) = MAX( 0.D+0,WL(M,NPZ) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLS( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLS'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPY = NSY(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HY = PLB(MP,NB)-PL(MP,N)-5.D-1*GRVY(NPY)
     &    *DYGF(N)*RP(I)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKLB(2,MP,NB),RKL(2,MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
        VL(M,NPY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          VL(M,NPY) = MIN( 0.D+0,VL(M,NPY) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          VL(M,NPY) = MAX( 0.D+0,VL(M,NPY) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLW( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLW'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPX = NSX(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HX = PLB(MP,NB)-PL(MP,N)-5.D-1*GRVX(NPX)
     &    *DXGF(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKLB(1,MP,NB),RKL(1,MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
        UL(M,NPX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          UL(M,NPX) = MIN( 0.D+0,UL(M,NPX) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          UL(M,NPX) = MAX( 0.D+0,UL(M,NPX) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLE( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLE'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQX = NSX(N)+1
      DO M = 1,ISVF
        MN = MNEG(M)
        HX = PL(MN,N)-PLB(MN,NB)-5.D-1*GRVX(NQX)
     &    *DXGF(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKL(1,MN,N),RKLB(1,MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(1,IZ(N))
        UL(M,NQX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          UL(M,NQX) = MAX( 0.D+0,UL(M,NQX) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          UL(M,NQX) = MIN( 0.D+0,UL(M,NQX) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLN( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLN'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQY = NSY(N)+IFLD
      DO M = 1,ISVF
        MN = MNEG(M)
        HY = PL(MN,N)-PLB(MN,NB)-5.D-1*GRVY(NQY)
     &    *DYGF(N)*RP(I)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKL(2,MN,N),RKLB(2,MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(2,IZ(N))
        VL(M,NQY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          VL(M,NQY) = MAX( 0.D+0,VL(M,NQY))
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          VL(M,NQY) = MIN( 0.D+0,VL(M,NQY))
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLT( N,NB )
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
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DO M = 1,ISVF
        MN = MNEG(M)
        HZ = PL(MN,N)-PLB(MN,NB)-5.D-1*GRVZ(NQZ)
     &    *DZGF(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKL(3,MN,N),RKLB(3,MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(3,IZ(N))
        WL(M,NQZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          WL(M,NQZ) = MAX( 0.D+0,WL(M,NQZ) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          WL(M,NQZ) = MIN( 0.D+0,WL(M,NQZ) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLT group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVNB( N,NB )
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
!     Written by MD White, PNNL, 10 May 2005.
!     Last modified by MD White, PNNL, 10 May 2005.
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
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DRCVNB'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPZ = NSZ(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HZ = PNB(MP,NB) - PN(MP,N)
     &    -5.D-1*GRVZ(NPZ)*DZGF(N)*RHONB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 10
        RKM = DIFMN(RKNB(MP,NB),RKN(MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 7
        VM = DIFMN(VISNB(MP,NB),VISN(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
        WN(M,NPZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          WN(M,NPZ) = MIN( 0.D+0,WN(M,NPZ) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          WN(M,NPZ) = MAX( 0.D+0,WN(M,NPZ) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVNB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVNS( N,NB )
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
!     Written by MD White, PNNL, 10 May 2005.
!     Last modified by MD White, PNNL, 10 May 2005.
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
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DRCVNS'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPY = NSY(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HY = PNB(MP,NB) - PN(MP,N)
     &    -5.D-1*GRVY(NPY)*DYGF(N)*RP(I)*RHONB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 10
        RKM = DIFMN(RKNB(MP,NB),RKN(MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 7
        VM = DIFMN(VISNB(MP,NB),VISN(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
        VN(M,NPY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          VN(M,NPY) = MIN( 0.D+0,VN(M,NPY) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          VN(M,NPY) = MAX( 0.D+0,VN(M,NPY) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVNS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVNW( N,NB )
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
!     Written by MD White, PNNL, 10 May 2005.
!     Last modified by MD White, PNNL, 10 May 2005.
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
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DRCVNW'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPX = NSX(N)
      DO M = 1,ISVF
        MP = MPOS(M)
        HX = PNB(MP,NB) - PN(MP,N)
     &    -5.D-1*GRVX(NPX)*DXGF(N)*RHONB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 10
        RKM = DIFMN(RKNB(MP,NB),RKN(MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 7
        VM = DIFMN(VISNB(MP,NB),VISN(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
        UN(M,NPX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          UN(M,NPX) = MIN( 0.D+0,UN(M,NPX) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          UN(M,NPX) = MAX( 0.D+0,UN(M,NPX) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVNW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVNE( N,NB )
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
!     Written by MD White, PNNL, 10 May 2005.
!     Last modified by MD White, PNNL, 10 May 2005.
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
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DRCVNE'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQX = NSX(N)+1
      DO M = 1,ISVF
        MN = MNEG(M)
        HX = PN(MN,N) - PNB(MN,NB)
     &    -5.D-1*GRVX(NQX)*DXGF(N)*RHONB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 10
        RKM = DIFMN(RKN(MN,N),RKNB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 7
        VM = DIFMN(VISN(MN,N),VISNB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(1,IZ(N))
        UN(M,NQX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          UN(M,NQX) = MAX( 0.D+0,UN(M,NQX) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          UN(M,NQX) = MIN( 0.D+0,UN(M,NQX) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVNE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVNN( N,NB )
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
!     Written by MD White, PNNL, 10 May 2005.
!     Last modified by MD White, PNNL, 10 May 2005.
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
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DRCVNN'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQY = NSY(N)+IFLD
      DO M = 1,ISVF
        MN = MNEG(M)
        HY = PN(MN,N) - PNB(MN,NB)
     &    -5.D-1*GRVY(NQY)*DYGF(N)*RP(I)*RHONB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 10
        RKM = DIFMN(RKN(MN,N),RKNB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 7
        VM = DIFMN(VISN(MN,N),VISNB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(2,IZ(N))
        VN(M,NQY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          VN(M,NQY) = MAX( 0.D+0,VN(M,NQY))
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          VN(M,NQY) = MIN( 0.D+0,VN(M,NQY))
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVNN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVNT( N,NB )
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
!     Written by MD White, PNNL, 10 May 2005.
!     Last modified by MD White, PNNL, 10 May 2005.
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
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DRCVNT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DO M = 1,ISVF
        MN = MNEG(M)
        HZ = PN(MN,N) - PNB(MN,NB)
     &    -5.D-1*GRVZ(NQZ)*DZGF(N)*RHONB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 10
        RKM = DIFMN(RKN(MN,N),RKNB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 7
        VM = DIFMN(VISN(MN,N),VISNB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(3,IZ(N))
        WN(M,NQZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Dirichlet-outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          WN(M,NQZ) = MAX( 0.D+0,WN(M,NQZ) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          WN(M,NQZ) = MIN( 0.D+0,WN(M,NQZ) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVNT group  ---
!
      RETURN
      END


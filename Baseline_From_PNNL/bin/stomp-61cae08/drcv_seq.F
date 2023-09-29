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
!     Written by MD White, PNNL, 18 July 2011.
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
            HDGX = PSO(MN,NW) - PSO(MP,N) - 0.5D+0*GRVX(NPX)*
     &          (RHOG(MN,NW)*DXGF(N)+RHOG(MP,N)*DXGF(NW))
!            HDGX = PSO(MN,NW) - PSO(MP,N) - (ZP(N)-ZP(NW))*GRAVZ*
!     &          (RHOG(MN,NW)*DXGF(N)+RHOG(MP,N)*DXGF(NW))/
!     &          (DXGF(N)+DXGF(NW))
            IF( M.EQ.1 ) HDG = HDGX
!
!---        Permeability reduction factor
!
            PERM_WX = PERMRF(MN,NW)*PERMV(1,NW)
            PERM_PX = PERMRF(MP,N)*PERMV(1,N)
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
            UG(M,NPX) = KGM*RKGM*HDGX/DXGP(NPX)/VGM
  100     CONTINUE
  200 CONTINUE
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
            HDGY = PSO(MN,NS) - PSO(MP,N) - 0.5D+0*GRVY(NPY)*
     &          (RHOG(MN,NS)*DYGF(N)+RHOG(MP,N)*DYGF(NS))
!            HDGY = PSO(MN,NS) - PSO(MP,N) - (ZP(N)-ZP(NS))*GRAVZ*
!     &          (RHOG(MN,NS)*DYGF(N)+RHOG(MP,N)*DYGF(NS))/
!     &          (DYGF(N)+DYGF(NS))
            IF( M.EQ.1 ) HDG = HDGY
!
!---        Permeability reduction factor
!
            PERM_SX = PERMRF(MN,NS)*PERMV(2,NS)
            PERM_PX = PERMRF(MP,N)*PERMV(2,N)
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
            VG(M,NPY) = KGM*RKGM*HDGY/DYGP(NPY)/VGM/RP(I)
  300     CONTINUE
  400 CONTINUE
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
            HDGZ = PSO(MN,NB) - PSO(MP,N) - 0.5D+0*GRVZ(NPZ)*
     &          (RHOG(MN,NB)*DZGF(N)+RHOG(MP,N)*DZGF(NB))
!            HDGZ = PSO(MN,NB) - PSO(MP,N) - (ZP(N)-ZP(NB))*GRAVZ*
!     &          (RHOG(MN,NB)*DZGF(N)+RHOG(MP,N)*DZGF(NB))/
!     &          (DZGF(N)+DZGF(NB))
            IF( M.EQ.1 ) HDG = HDGZ
!
!---        Permeability reduction factor
!
            PERM_BX = PERMRF(MN,NB)*PERMV(3,NB)
            PERM_PX = PERMRF(MP,N)*PERMV(3,N)
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
            WG(M,NPZ) = KGM*RKGM*HDGZ/DZGP(NPZ)/VGM
  500     CONTINUE
  600 CONTINUE
      ENDIF
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
!     Written by MD White, PNNL, 18 July 2011.
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
            INDX = 11
            HDLX = PL(MN,NW) - PL(MP,N) - 0.5D+0*GRVX(NPX)*
     &        (RHOL(MN,NW)*DXGF(N)+RHOL(MP,N)*DXGF(NW))
!            HDLX = PL(MN,NW) - PL(MP,N) - (ZP(N)-ZP(NW))*GRAVZ*
!     &          (RHOL(MN,NW)*DXGF(N)+RHOL(MP,N)*DXGF(NW))/
!     &          (DXGF(N)+DXGF(NW))
            IF( M.EQ.1 ) HDL = HDLX
!
!---        Permeability reduction factor
!
            PERM_WX = PERMRF(MN,NW)*PERMV(1,NW)
            PERM_PX = PERMRF(MP,N)*PERMV(1,N)
            INDX = 11
            KLM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDL,INDX)
            INDX = 8
            RKLM = DIFMN(RKL(1,MN,NW),RKL(1,MP,N),DXGF(NW),
     &        DXGF(N),HDL,INDX)
            IF( PERM_WX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 5
            VLM = DIFMN(VISL(MN,NW),VISL(MP,N),DXGF(NW),DXGF(N),
     &        HDL,INDX)
            UL(M,NPX) = KLM*RKLM*HDLX/DXGP(NPX)/VLM
  100     CONTINUE
  200 CONTINUE
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
            INDX = 11
            HDLY = PL(MN,NS) - PL(MP,N) - 0.5D+0*GRVY(NPY)*
     &        (RHOL(MN,NS)*DYGF(N)+RHOL(MP,N)*DYGF(NS))
!            HDLY = PL(MN,NS) - PL(MP,N) - (ZP(N)-ZP(NS))*GRAVZ*
!     &          (RHOL(MN,NS)*DYGF(N)+RHOL(MP,N)*DYGF(NS))/
!     &          (DYGF(N)+DYGF(NS))
            IF( M.EQ.1 ) HDL = HDLY
!
!---        Permeability reduction factor
!
            PERM_SX = PERMRF(MN,NS)*PERMV(2,NS)
            PERM_PX = PERMRF(MP,N)*PERMV(2,N)
            INDX = 11
            KLM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDL,INDX)
            INDX = 8
            RKLM = DIFMN(RKL(2,MN,NS),RKL(2,MP,N),DYGF(NS),
     &        DYGF(N),HDL,INDX)
            IF( PERM_SX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 5
            VLM = DIFMN(VISL(MN,NS),VISL(MP,N),DYGF(NS),DYGF(N),
     &        HDL,INDX)
            VL(M,NPY) = KLM*RKLM*HDLY/DYGP(NPY)/VLM/RP(I)
  300     CONTINUE
  400 CONTINUE
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
            INDX = 11
            HDLZ = PL(MN,NB) - PL(MP,N) - 0.5D+0*GRVZ(NPZ)*
     &        (RHOL(MN,NB)*DZGF(N)+RHOL(MP,N)*DZGF(NB))
!            HDLZ = PL(MN,NB) - PL(MP,N) - (ZP(N)-ZP(NB))*GRAVZ*
!     &          (RHOL(MN,NB)*DZGF(N)+RHOL(MP,N)*DZGF(NB))/
!     &          (DZGF(N)+DZGF(NB))
            IF( M.EQ.1 ) HDL = HDLZ
!
!---        Permeability reduction factor
!
            PERM_BX = PERMRF(MN,NB)*PERMV(3,NB)
            PERM_PX = PERMRF(MP,N)*PERMV(3,N)
            INDX = 11
            KLM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDL,INDX)
            INDX = 8
            RKLM = DIFMN(RKL(3,MN,NB),RKL(3,MP,N),DZGF(NB),
     &        DZGF(N),HDL,INDX)
            IF( PERM_BX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 5
            VLM = DIFMN(VISL(MN,NB),VISL(MP,N),DZGF(NB),DZGF(N),
     &        HDL,INDX)
            WL(M,NPZ) = KLM*RKLM*HDLZ/DZGP(NPZ)/VLM
  500     CONTINUE
  600 CONTINUE
      ENDIF
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 July 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
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
!
!---  X-direction Darcy velocities, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDNX = PN(MN,NW) - PN(MP,N) - 0.5D+0*GRVX(NPX)*
     &          (RHON(MN,NW)*DXGF(N)+RHON(MP,N)*DXGF(NW))
            IF( M.EQ.1 ) HDN = HDNX
!
!---        Permeability reduction factor
!
            PERM_WX = PERMRF(MN,NW)*PERMV(1,NW)
            PERM_PX = PERMRF(MP,N)*PERMV(1,N)
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
            UN(M,NPX) = KNM*RKNM*HDNX/DXGP(NPX)/VNM
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = ICM(1,2,N)
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDNY = PN(MN,NS) - PN(MP,N) - 0.5D+0*GRVY(NPY)*
     &          (RHON(MN,NS)*DYGF(N)+RHON(MP,N)*DYGF(NS))
            IF( M.EQ.1 ) HDN = HDNY
!
!---        Permeability reduction factor
!
            PERM_SX = PERMRF(MN,NS)*PERMV(2,NS)
            PERM_PX = PERMRF(MP,N)*PERMV(2,N)
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
            VN(M,NPY) = KNM*RKNM*HDNY/DYGP(NPY)/VNM/RP(I)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction Darcy velocities, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = ICM(1,1,N)
          IF( NB.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDNZ = PN(MN,NB) - PN(MP,N) - 0.5D+0*GRVZ(NPZ)*
     &          (RHON(MN,NB)*DZGF(N)+RHON(MP,N)*DZGF(NB))
            IF( M.EQ.1 ) HDN = HDNZ
!
!---        Permeability reduction factor
!
            PERM_BX = PERMRF(MN,NB)*PERMV(3,NB)
            PERM_PX = PERMRF(MP,N)*PERMV(3,N)
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
            WN(M,NPZ) = KNM*RKNM*HDNZ/DZGP(NPZ)/VNM
          ENDDO
        ENDDO
      ENDIF
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
!     drcv_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        HZ = PSOB(MP,NB) - PSO(MP,N)
     &    - 5.D-1*GRVZ(NPZ)*DZGF(N)*RHOGB(MP,NB)
!        HZ = PSOB(MP,NB) - PSO(MP,N) - (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERMV(3,N)
        WG(M,NPZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          WG(M,NPZ) = MIN( 0.D+0,WG(M,NPZ) )
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        HY = PSOB(MP,NB) - PSO(MP,N)
     &    - 5.D-1*GRVY(NPY)*DYGF(N)*RP(I)*RHOGB(MP,NB)
!        HY = PSOB(MP,NB) - PSO(MP,N) - (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERMV(2,N)
        VG(M,NPY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          VG(M,NPY) = MIN( 0.D+0,VG(M,NPY) )
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        HX = PSOB(MP,NB) - PSO(MP,N)
     &    - 5.D-1*GRVX(NPX)*DXGF(N)*RHOGB(MP,NB)
!        HX = PSOB(MP,NB) - PSO(MP,N) - (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERMV(1,N)
        UG(M,NPX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          UG(M,NPX) = MIN( 0.D+0,UG(M,NPX) )
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        HX = PSO(MN,N) - PSOB(MN,NB)
     &    - 5.D-1*GRVX(NQX)*DXGF(N)*RHOGB(MN,NB)
!        HX = PSO(MN,N) - PSOB(MN,NB) + (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERMV(1,N)
        UG(M,NQX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          UG(M,NQX) = MAX( 0.D+0,UG(M,NQX) )
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        HY = PSO(MN,N) - PSOB(MN,NB)
     &    - 5.D-1*GRVY(NQY)*DYGF(N)*RP(I)*RHOGB(MN,NB)
!        HY = PSO(MN,N) - PSOB(MN,NB) + (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERMV(2,N)
        VG(M,NQY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          VG(M,NQY) = MAX( 0.D+0,VG(M,NQY) )
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        HZ = PSO(MN,N) - PSOB(MN,NB)
     &    - 5.D-1*GRVZ(NQZ)*DZGF(N)*RHOGB(MN,NB)
!        HZ = PSO(MN,N) - PSOB(MN,NB) + (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERMV(3,N)
        WG(M,NQZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          WG(M,NQZ) = MAX( 0.D+0,WG(M,NQZ) )
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
!     drcv_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        HZ = PLB(MP,NB) - PL(MP,N)
     &    -5.D-1*GRVZ(NPZ)*DZGF(N)*RHOLB(MP,NB)
!        HZ = PLB(MP,NB) - PL(MP,N) - (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKLB(3,MP,NB),RKL(3,MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERMV(3,N)
        WL(M,NPZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          WL(M,NPZ) = MIN( 0.D+0,WL(M,NPZ) )
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        HY = PLB(MP,NB) - PL(MP,N)
     &    -5.D-1*GRVY(NPY)*DYGF(N)*RP(I)*RHOLB(MP,NB)
!        HY = PLB(MP,NB) - PL(MP,N) - (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKLB(2,MP,NB),RKL(2,MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERMV(2,N)
        VL(M,NPY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          VL(M,NPY) = MIN( 0.D+0,VL(M,NPY) )
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        HX = PLB(MP,NB) - PL(MP,N)
     &    -5.D-1*GRVX(NPX)*DXGF(N)*RHOLB(MP,NB)
!        HX = PLB(MP,NB) - PL(MP,N) - (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKLB(1,MP,NB),RKL(1,MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERMV(1,N)
        UL(M,NPX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          UL(M,NPX) = MIN( 0.D+0,UL(M,NPX) )
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        HX = PL(MN,N) - PLB(MN,NB)
     &    -5.D-1*GRVX(NQX)*DXGF(N)*RHOLB(MN,NB)
!        HX = PL(MN,N) - PLB(MN,NB) + (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKL(1,MN,N),RKLB(1,MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERMV(1,N)
        UL(M,NQX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          UL(M,NQX) = MAX( 0.D+0,UL(M,NQX) )
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        HY = PL(MN,N) - PLB(MN,NB)
     &    -5.D-1*GRVY(NQY)*DYGF(N)*RP(I)*RHOLB(MN,NB)
!        HY = PL(MN,N) - PLB(MN,NB) + (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKL(2,MN,N),RKLB(2,MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERMV(2,N)
        VL(M,NQY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          VL(M,NQY) = MAX( 0.D+0,VL(M,NQY))
        ENDIF
  100 CONTINUE
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        HZ = PL(MN,N) - PLB(MN,NB)
     &    -5.D-1*GRVZ(NQZ)*DZGF(N)*RHOLB(MN,NB)
!        HZ = PL(MN,N) - PLB(MN,NB) + (ZP(N)-ZPBC(NB))*GRAVZ*
!     &    RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKL(3,MN,N),RKLB(3,MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERMV(3,N)
        WL(M,NQZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          WL(M,NQZ) = MAX( 0.D+0,WL(M,NQZ) )
        ENDIF
  100 CONTINUE
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a bottom boundary.
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
        PERM_PX = PERMRF(MP,N)*PERMV(3,N)
        WN(M,NPZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(2,NB)).EQ.7 ) THEN
          WN(M,NPZ) = MIN( 0.D+0,WN(M,NPZ) )
        ENDIF
      ENDDO
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
        PERM_PX = PERMRF(MP,N)*PERMV(2,N)
        VN(M,NPY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(2,NB)).EQ.7 ) THEN
          VN(M,NPY) = MIN( 0.D+0,VN(M,NPY) )
        ENDIF
      ENDDO
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
        PERM_PX = PERMRF(MP,N)*PERMV(1,N)
        UN(M,NPX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(2,NB)).EQ.7 ) THEN
          UN(M,NPX) = MIN( 0.D+0,UN(M,NPX) )
        ENDIF
      ENDDO
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
      DO M = 1,ISVF
        MN = MNEG(M)
        HX = PN(MN,N) - PNB(MN,NB)
     &    -5.D-1*GRVX(NQX)*DXGF(N)*RHONB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 10
        RKM = DIFMN(RKN(MN,N),RKNB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 7
        VM = DIFMN(VISN(MN,N),VISNB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERMV(1,N)
        UN(M,NQX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(2,NB)).EQ.7 ) THEN
          UN(M,NQX) = MAX( 0.D+0,UN(M,NQX) )
        ENDIF
      ENDDO
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
      DO M = 1,ISVF
        MN = MNEG(M)
        HY = PN(MN,N) - PNB(MN,NB)
     &    -5.D-1*GRVY(NQY)*DYGF(N)*RP(I)*RHONB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 10
        RKM = DIFMN(RKN(MN,N),RKNB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 7
        VM = DIFMN(VISN(MN,N),VISNB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERMV(2,N)
        VN(M,NQY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(2,NB)).EQ.7 ) THEN
          VN(M,NQY) = MAX( 0.D+0,VN(M,NQY))
        ENDIF
      ENDDO
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
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
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
      DO M = 1,ISVF
        MN = MNEG(M)
        HZ = PN(MN,N) - PNB(MN,NB)
     &    -5.D-1*GRVZ(NQZ)*DZGF(N)*RHONB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 10
        RKM = DIFMN(RKN(MN,N),RKNB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 7
        VM = DIFMN(VISN(MN,N),VISNB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERMV(3,N)
        WN(M,NQZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(2,NB)).EQ.7 ) THEN
          WN(M,NQZ) = MAX( 0.D+0,WN(M,NQZ) )
        ENDIF
      ENDDO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVNT group  ---
!
      RETURN
      END

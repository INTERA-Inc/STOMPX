!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLW_W
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
!     Water Mode (STOMPX-W)
!
!     Diffusive water aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 January 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
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
      SUB_LOG(ISUB_LOG) = '/DFFLW_W'
!
!---  Aqueous diffusive flux, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWW = XLW(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,1,N),INDX )
            ULW(M,1,N) = UL(M,1,N)*FLW
          ENDDO
          DO M = 1,ISVF
            ULW(M,2,NW) = ULW(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            ULW(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWS = XLW(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLW = DIFMN( FLWS,FLWP,DYGF(NS),DYGF(N),VL(1,1,N),INDX )
            VLW(M,1,N) = VL(M,1,N)*FLW
          ENDDO
          DO M = 1,ISVF
            VLW(M,2,NS) = VLW(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            VLW(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWB = XLW(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL(1,1,N),INDX )
            WLW(M,1,N) = WL(M,1,N)*FLW
          ENDDO
          DO M = 1,ISVF
            WLW(M,2,NB) = WLW(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            WLW(M,1,N) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLW_W group  ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVL_W
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
!     Water Mode (STOMPX-W)
!
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 3 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
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
      SUB_LOG(ISUB_LOG) = '/DRCVL_W'
!
!---  Aqueous Darcy velocities, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          OECX = 0.D+0
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            INDX = 11
            HDLX = PL(MN,NW)-PL(MP,N)-0.5D+0*GRVX(1,N)
     &        *(RHOL(MN,NW)*DXGF(N)+RHOL(MP,N)*DXGF(NW))
            IF( M.EQ.1 ) HDL = HDLX
!
!---        Permeability reduction factor
!
            PERM_WX = PERMRF(MN,NW)*PERM(1,NW)
            PERM_PX = PERMRF(MP,N)*PERM(1,N)
            INDX = 11
            KLM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDL,INDX)
            IF( PERM_WX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 8
            RKLM = DIFMN(RKL(MN,NW),RKL(MP,N),DXGF(NW),DXGF(N),
     &        HDL,INDX)
            INDX = 5
            VLM = DIFMN(VISL(MN,NW),VISL(MP,N),DXGF(NW),DXGF(N),
     &        HDL,INDX)
            UL(M,1,N) = KLM*RKLM*HDLX/(DXGP(1,N)*VLM)
          ENDDO
          DO M = 1,ISVF
            UL(M,2,NW) = UL(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            UL(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          OECY = 0.D+0
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            INDX = 11
            HDLY = PL(MN,NS)-PL(MP,N)-0.5D+0*GRVY(1,N)
     &        *(RHOL(MN,NS)*DYGF(N)+RHOL(MP,N)*DYGF(NS))
            IF( M.EQ.1 ) HDL = HDLY
!
!---        Permeability reduction factor
!
            PERM_SX = PERMRF(MN,NS)*PERM(2,NS)
            PERM_PX = PERMRF(MP,N)*PERM(2,N)
            INDX = 11
            KLM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDL,INDX)
            IF( PERM_SX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 8
            RKLM = DIFMN(RKL(MN,NS),RKL(MP,N),DYGF(NS),DYGF(N),HDL,INDX)
            INDX = 5
            VLM = DIFMN(VISL(MN,NS),VISL(MP,N),DYGF(NS),DYGF(N),
     &        HDL,INDX)
            VL(M,1,N) = KLM*RKLM*HDLY/(DYGP(1,N)*VLM*RP(N))
          ENDDO
          DO M = 1,ISVF
            VL(M,2,NS) = VL(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            VL(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          OECZ = 0.D+0
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            INDX = 11
            HDLZ = PL(MN,NB)-PL(MP,N)-0.5D+0*GRVZ(1,N)
     &        *(RHOL(MN,NB)*DZGF(N)+RHOL(MP,N)*DZGF(NB))
            IF( M.EQ.1 ) HDL = HDLZ
!
!---        Permeability reduction factor
!
            PERM_BX = PERMRF(MN,NB)*PERM(3,NB)
            PERM_PX = PERMRF(MP,N)*PERM(3,N)
            INDX = 11
            KLM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDL,INDX)
            IF( PERM_BX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 8
            RKLM = DIFMN(RKL(MN,NB),RKL(MP,N),DZGF(NB),DZGF(N),
     &        HDL,INDX)
            INDX = 5
            VLM = DIFMN(VISL(MN,NB),VISL(MP,N),DZGF(NB),DZGF(N),
     &        HDL,INDX)
            WL(M,1,N) = KLM*RKLM*HDLZ/(DZGP(1,N)*VLM)
          ENDDO
          DO M = 1,ISVF
            WL(M,2,NB) = WL(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            WL(M,1,N) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVL_W group  ---
!
      RETURN
      END



!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGB( N,NB )
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
!     VOC vapor, and air molecular diffusion on a bottom boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, April 8, 1994.
!     thdgb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGB'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        WDGAX = -WDGW(M,NPZ)-WDGO(M,NPZ)
        WQ(M,NPZ) = WQ(M,NPZ) + HGAB(MP,NB)*MAX(WDGAX,ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-WDGAX,ZERO)*WTMA
     &    + HGWB(MP,NB)*MAX(WDGW(M,NPZ),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-WDGW(M,NPZ),ZERO)*WTMW
     &    + HGOB(MP,NB)*MAX(WDGO(M,NPZ),ZERO)*WTMO
     &    - HGO(MP,N)*MAX(-WDGO(M,NPZ),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGE( N,NB )
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
!     VOC vapor, and air molecular diffusion on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGE'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        UDGAX = -UDGW(M,NQX)-UDGO(M,NQX)
        UQ(M,NQX) = UQ(M,NQX) - HGAB(MN,NB)*MAX(-UDGAX,ZERO)*WTMA
     &    + HGA(MN,N)*MAX(UDGAX,ZERO)*WTMA
     &    - HGWB(MN,NB)*MAX(-UDGW(M,NQX),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(UDGW(M,NQX),ZERO)*WTMW
     &    - HGOB(MN,NB)*MAX(-UDGO(M,NQX),ZERO)*WTMO
     &    + HGO(MN,N)*MAX(UDGO(M,NQX),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGN( N,NB )
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
!     VOC vapor, and air molecular diffusion on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGN'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        VDGAX = -VDGW(M,NQY)-VDGO(M,NQY)
        VQ(M,NQY) = VQ(M,NQY) - HGAB(MN,NB)*MAX(-VDGAX,ZERO)*WTMA
     &    + HGA(MN,N)*MAX(VDGAX,ZERO)*WTMA
     &    - HGWB(MN,NB)*MAX(-VDGW(M,NQY),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(VDGW(M,NQY),ZERO)*WTMW
     &    - HGOB(MN,NB)*MAX(-VDGO(M,NQY),ZERO)*WTMO
     &    + HGO(MN,N)*MAX(VDGO(M,NQY),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGS( N,NB )
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
!     VOC vapor, and air molecular diffusion on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGS'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        VDGAX = -VDGW(M,NPY)-VDGO(M,NPY)
        VQ(M,NPY) = VQ(M,NPY) + HGAB(MP,NB)*MAX(VDGAX,ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-VDGAX,ZERO)*WTMA
     &    + HGWB(MP,NB)*MAX(VDGW(M,NPY),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-VDGW(M,NPY),ZERO)*WTMW
     &    + HGOB(MP,NB)*MAX(VDGO(M,NPY),ZERO)*WTMO
     &    - HGO(MP,N)*MAX(-VDGO(M,NPY),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGT( N,NB )
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
!     VOC vapor, and air molecular diffusion on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGT'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        WDGAX = -WDGW(M,NQZ)-WDGO(M,NQZ)
        WQ(M,NQZ) = WQ(M,NQZ) - HGAB(MN,NB)*MAX(-WDGAX,ZERO)*WTMA
     &    + HGA(MN,N)*MAX(WDGAX,ZERO)*WTMA
     &    - HGWB(MN,NB)*MAX(-WDGW(M,NQZ),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(WDGW(M,NQZ),ZERO)*WTMW
     &    - HGOB(MN,NB)*MAX(-WDGO(M,NQZ),ZERO)*WTMO
     &    + HGO(MN,N)*MAX(WDGO(M,NQZ),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGW( N,NB )
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
!     VOC vapor, and air molecular diffusion on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
!     thdgb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGW'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        UDGAX = -UDGW(M,NPX)-UDGO(M,NPX)
        UQ(M,NPX) = UQ(M,NPX) + HGAB(MP,NB)*MAX(UDGAX,ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-UDGAX,ZERO)*WTMA
     &    + HGWB(MP,NB)*MAX(UDGW(M,NPX),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-UDGW(M,NPX),ZERO)*WTMW
     &    + HGOB(MP,NB)*MAX(UDGO(M,NPX),ZERO)*WTMO
     &    - HGO(MP,N)*MAX(-UDGO(M,NPX),ZERO)*WMTO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGW group  ---
!
      RETURN
      END

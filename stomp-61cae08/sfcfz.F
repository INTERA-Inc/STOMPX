!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFCFZ
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
!     Zero surfactant transport flux for interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1997.
!     Last Modified by MD White on January 20, 1997.
!     sfcfz.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFCFZ'
!
!---  X-direction surfactant flux
!
      DO 100 K = 1,KFLD
        DO 100 J = 1,JFLD
          DO 100 I = 2,IFLD
            N = ND(I,J,K)
            NPX = NSX(N)
            DO 100 M = 1,ISVF
              US(M,NPX) = 0.D+0
  100 CONTINUE
!
!---  Y-direction surfactant flux
!
      DO 200 I = 1,IFLD
        DO 200 K = 1,KFLD
          DO 200 J = 2,JFLD
            N = ND(I,J,K)
            NPY = NSY(N)
            DO 200 M = 1,ISVF
              VS(M,NPY) = 0.D+0
  200 CONTINUE
!
!---  Z-direction surfactant flux
!
      DO 300 J = 1,JFLD
        DO 300 I = 1,IFLD
          DO 300 K = 2,KFLD
            N = ND(I,J,K)
            NPZ = NSZ(N)
            DO 300 M = 1,ISVF
              WS(M,NPZ) = 0.D+0
  300 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFCFZ group
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXZ( NSL )
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
!     Zero solute transport flux.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1993.
!     Last Modified by MD White, Battelle, PNL, July 15, 1993.
!     sfxz.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFXZ'
!
!---  X-direction solute flux
!
      DO N = 1,NFLD
        I = ID(N)
        UC(NSX(N),NSL) = 0.D+0
        IF( I.EQ.IFLD ) THEN
          UC(NSX(N)+1,NSL) = 0.D+0
        ENDIF
      ENDDO
!
!---  Y-direction solute flux
!
      DO N = 1,NFLD
        J = JD(N)
        VC(NSY(N),NSL) = 0.D+0
        IF( J.EQ.JFLD ) THEN
          VC(NSY(N)+IFLD,NSL) = 0.D+0
        ENDIF
      ENDDO
!
!---  Z-direction solute flux
!
      DO N = 1,NFLD
        K = KD(N)
        WC(NSZ(N),NSL) = 0.D+0
        IF( K.EQ.KFLD ) THEN
          WC(NSZ(N)+IJFLD,NSL) = 0.D+0
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFXZ group
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE GASVIS( XMW,XMO,XMA,VISW,VISO,VISA,VISGX)
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
!     Calculate the gas viscosity using an extension of the
!     Chapman-Enskog theory for multicomponent gas mixtures at low
!     density. pp. 528-530.
!
!     Hirschfelder, J.O., C.F. Curtiss, and R.B. Bird. 1954. Molecular
!     Theory of Gases and Liquids. John Wiley & Sons.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     gasvis.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/GASVIS'
      PHIAW = ((1.D+0+SQRT(VISA/VISW)*SQRT(SQRT(WTMW/WTMA)))**2)
     &  /SQRT(8.D+0*(1.D+0 + WTMA/WTMW))
      PHIWA = ((1.D+0+SQRT(VISW/VISA)*SQRT(SQRT(WTMA/WTMW)))**2)
     &  /SQRT(8.D+0*(1.D+0 + WTMW/WTMA))
      IF( WTMO.GT.EPSL ) THEN
        PHIOW = ((1.D+0+SQRT(VISO/VISW)*SQRT(SQRT(WTMW/WTMO)))**2)
     &    /SQRT(8.D+0*(1.D+0 + WTMO/WTMW))
        PHIWO = ((1.D+0+SQRT(VISW/VISO)*SQRT(SQRT(WTMO/WTMW)))**2)
     &   /SQRT(8.D+0*(1.D+0 + WTMW/WTMO))
        PHIAO = ((1.D+0+SQRT(VISA/VISO)*SQRT(SQRT(WTMO/WTMA)))**2)
     &    /SQRT(8.D+0*(1.D+0 + WTMA/WTMO))
        PHIOA = ((1.D+0+SQRT(VISO/VISA)*SQRT(SQRT(WTMA/WTMO)))**2)
     &    /SQRT(8.D+0*(1.D+0 + WTMO/WTMA))
        CHIW = XMW + XMA*PHIWA + XMO*PHIWO
        CHIO = XMW*PHIOW + XMA*PHIOA + XMO
        CHIA = XMW*PHIAW + XMA + XMO*PHIAO
        VISGX = XMW*VISW/CHIW + XMO*VISO/CHIO + XMA*VISA/CHIA
      ELSE
        CHIW = XMW + XMA*PHIWA
        CHIA = XMW*PHIAW + XMA
        VISGX = XMW*VISW/CHIW + XMA*VISA/CHIA
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of GASVIS group  ---
!
      RETURN
      END

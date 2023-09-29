!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFCSCF( XLOX,XLSX,RHOX,BAWX,BOWX,BAOX )
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
!     Interfacial tension scaling factors for surfactant solutions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1997.
!     Last Modified by MD White on January 24, 1997.
!     sfcscf.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XMSS(LCMP),XMOL(LCMP)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFCSCF'
!
!---  Compute the aqueous concentration of surfactant  ---
!
      CLSX = XLSX*RHOX
!
!---  Compute phase-pair scaling factors as a function of surfactant
!     concentration assuming a constant air-oil interfacial surface
!     tension  ---
!
      IF( CLSX.GE.SFCSF(1) ) THEN
       BOWX = (BOWX*SFCSF(3))/SFCSF(4)
       BAWX = 1.D+0/((1.D+0/(BOWX+SMALL))+(1.D+0/(BAOX+SMALL)))
      ELSE
       BOWX = BOWX*SFCSF(3)/(CLSX*(SFCSF(4)-SFCSF(3))/SFCSF(1)+SFCSF(3))
       BAWX = 1.D+0/((1.D+0/(BOWX+SMALL))+(1.D+0/(BAOX+SMALL)))
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFCSCF group  ---
!
      RETURN
      END


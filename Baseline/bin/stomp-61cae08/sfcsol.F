!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFCSOL( XLSX,TX,PX,PSOX,HCOWX,XLOX )
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
!     Dissolved VOC solubility in surfactant solution.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1997.
!     Last Modified by MD White on January 21, 1997.
!     sfcsol.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 XMSS(3),XMOL(3),XLOMX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFCSOL'
!
!---  Compute the oil solubility for pure water  ---
!
      XMOL(2) = MIN( PSOX/HCOWX,1.D+0 )
      XMOL(1) = 1.D+0-XMOL(2)
      NCMP = 2
      INDX = 1
      CALL FRAC( NCMP,XMSS,XMOL,INDX )
      CALL WATLQD( TX,PX,RHOX )
      CALL MXLQD( NCMP,XMOL,TX,PX,RHOX )
      CLOX = RHOX*XMSS(2)
!
!---  Estimate the oil solubility w/ surfactant  ---
!
      XLOX = (MAX( (XLSX*RHOX-SFCSF(1)),ZERO )*SFCSF(2)+CLOX)/RHOX
!
!---  Iterate to compute the oil solubility w/ surfactant ---
!
      NCMP = 3
      INDX = 0
      DXLOX = 1.D-6
      XMSS(3) = XLSX
      INDX = 0
 100  CONTINUE
      XMSS(2) = XLOX
      XMSS(1) = MAX( 1.D+0-XMSS(3)-XMSS(2),ZERO )
      CALL FRAC( NCMP,XMSS,XMOL,INDX )
      CALL WATLQD( TX,PX,RHOX )
      CALL MXLQD( NCMP,XMOL,TX,PX,RHOX )
      XLOMX(1) = (MAX((XMSS(3)*RHOX-SFCSF(1)),ZERO)*SFCSF(2)+CLOX)/RHOX
      XMSS(2) = XLOX+DXLOX
      XMSS(1) = MAX( 1.D+0-XMSS(3)-XMSS(2),ZERO )
      CALL FRAC( NCMP,XMSS,XMOL,INDX )
      CALL WATLQD( TX,PX,RHOX )
      CALL MXLQD( NCMP,XMOL,TX,PX,RHOX )
      XLOMX(2) = (MAX((XMSS(3)*RHOX-SFCSF(1)),ZERO)*SFCSF(2)+CLOX)/RHOX
      FX = XLOX - XLOMX(1)
      DFX = 1 - (XLOMX(2)-XLOMX(1))/DXLOX
      DXLO = -FX/DFX
      XLOX = XLOX+DXLO
      IF( ABS(DXLO).GT.1.D-8 ) GOTO 100
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFCSOL group  ---
!
      RETURN
      END


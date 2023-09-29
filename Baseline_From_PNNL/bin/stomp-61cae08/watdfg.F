!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATDFG( XMW,XMO,XMA,DFGAW,DFGOW,DFGWX )
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
!     Calculates water diffusion coefficient by the Wilke
!     method. pp. 34.
!
!     Falta, R.W., K. Pruess, I. Javandel, and P.A. Witherspoon. 1990.
!     Numerical Modeling of Steam Injection for the Removal of
!     Nonaqueous Phase Liquids from the Subsurface: 1 Numerical
!     Formulation. LBL-29615, Lawrence Berkeley Laboratory.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 14, 1992.
!     watdfg.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATDFG'
      IF( (1.D+0-XMW) .LE. 1.D-20 ) THEN
        DFGWX = 2.D+0*DFGAW*DFGOW/(DFGAW+DFGOW)
      ELSE
        DFGWX = (1.D+0-XMW)/((XMA/DFGAW) + (XMO/DFGOW))
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATDFG group  ---
!
      RETURN
      END

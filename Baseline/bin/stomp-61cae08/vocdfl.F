!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VOCDFL( TX,VISLX,DFLOX )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Calculates binary liquid diffusion coefficient for VOC in liquid
!     water, using the Wilke and Chang estimation method; pp. 598.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, June, 1993.
!     Last Modified by MD White, Battelle, PNL, June 2, 1993.
!     vocdfl.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/VOCDFL'
      VBO = 2.85D-1*VCRO**1.048D+0
      DFLOX = 7.4D-15*((2.6D+0*WTMW)**5.D-1)*(TX+TABS)/
     &  (VISLX*(VBO**6.D-1))
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VOCDFL group  ---
!
      RETURN
      END

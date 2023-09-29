!-------------------------Subroutine-----------------------------------!
!
      SUBROUTINE BALANC
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
!            Copyright © Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Compute global balance of mass/activity of the dilute solute
!     identified in the call argument and write the result to the
!     balance reporting file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, Battelle PNNL, May 2000.
!     Last Modified by WE Nichols, Battelle PNNL, June 12, 2003.
!
!     SAC Modification by WE Nichols, PNNL, June 12, 2003.
!     balanc.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BALANC'
!
!---  Loop through grid domain element-by-element;
!     accumulate dilute solute mass  ---
!
      DO 200 NS = 1,NSOLU
        GLBL = 0.
        DO 100 N = 1,NFLD
          GLBL = GLBL +C(N,NS)*VOL(N)
  100   CONTINUE
!
!---  Report global balance of mass/activity to 'balance' file  ---
!
        NCH = LEN_TRIM(SOLUT(NS))
        WRITE(UNIT=IBL,FMT=901) TM,',s,',SOLUT(NS)(1:NCH),',',GLBL,','
!
  200 CONTINUE
!
!---  Format statements  ---
!
  901 FORMAT(1PE19.12,3A,1PE19.12,A)
!
!---  End of BALANC group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END


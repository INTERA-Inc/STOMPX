!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SACREM
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
!     Computes mass or activity remediated in a SAC remediation
!     event and call another subroutine to report this quantity
!     to a remediation transfer file.
!
!     Concentration of the dilute solute is reduced in nodes subject
!     to remediation in accordance with the remedation description.
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, Battelle PNNL, May 2000.
!     Last Modified by WE Nichols, Battelle PNNL, July 14, 2004.
!
!     SAC Modification by WE Nichols, PNNL, July 14, 2004.
!     sacrem.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
#ifdef sac
      USE SACIN
#endif
#ifdef sac
      USE SACEV
#endif
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
      SUB_LOG(ISUB_LOG) = '/SACREM'
!
!---  Determine current remediation event index  ---
!
      DO 400 NRE = 1,INMREM
        IF( ABS(DRMTIM(NRE)-TM) .LT. 1.D-6 ) THEN

!---  Initialize remediation accumulator  ---
!
          REMEDQ = 0.
!
!---  Loop over grid domain  ---
!
          DO 300 K = 1,KFLD
            DO 200 J = 1,JFLD
              DO 100 I = 1,IFLD
                N = ND(I,J,K)
                NK = ND(I,J,KFLD)
!
!---  Check if current node centroid is within remediated depth  ---
!
                IF( (ZP(N) .GE. (ZE(5,NK) -DRMDLO(NRE))) .AND.
     &              (ZP(N) .LE. (ZE(5,NK) -DRMDUP(NRE))) ) THEN
!
!---  Accumulate the remediated mass/activity and
!     change the solute concentration to reflect remediation  ---
!
                  REMEDQ = REMEDQ +C(N,IRMSOL(NRE))*VOL(N)*DRMEFF(NRE)
                  C(N,IRMSOL(NRE)) = C(N,IRMSOL(NRE))*(1.-DRMEFF(NRE))
                ENDIF
  100         CONTINUE
  200       CONTINUE
  300     CONTINUE
!
!---  Write record to the remediation transfer file  ---
!
          NCH1 = LEN_TRIM(CRMSIT(NRE))
          NCH2 = LEN_TRIM(SOLUT(IRMSOL(NRE)))
          NCH3 = LEN_TRIM(CRMWFM(NRE))
          WRITE (UNIT=IVD,FMT=901) TM,',s,',CRMSIT(NRE)(1:NCH1),',',
     &      SOLUT(IRMSOL(NRE))(1:NCH2),',',REMEDQ,',,',
     &      CRMWFM(NRE)(1:NCH3),','
        ENDIF
!
  400 CONTINUE
!
!---  Format statements  ---
!
  901 FORMAT(1PE19.12,5A,1PE19.12,3A)
!
!---  End of SACREM group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

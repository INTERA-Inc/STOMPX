!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SACOUT
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
!     This subroutine is used whenever STOMP ends a simulation before
!     the simulation time is over; i.e., when a time step maximum is
!     exceeeded, a run time maximum is exceeded, or a time step
!     reduction limit (nonconvergence) is exceeded.  This subroutine
!     finishes writing the 'balance', 'release', and 'remediate' files.
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, Battelle PNNL, January 2002.
!     Last Modified by WE Nichols, Battelle PNNL, July 14, 2004.
!
!     SAC Modification by WE Nichols, PNNL, July 14, 2004.
!     sacout.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/SACOUT'
!
!---  Record constant mass/activity balances to 'balance' file  ---
!
      DO 400 NBL = 1,INMBLT
!
!---  Check if any of the mass balance times are after early
!     termination time for this simulaton  ---
!
        IF( DBLTIM(NBL) .GE. TM ) THEN
          IF( NBL .GT. 1 ) THEN
            IF( DBLTIM(NBL-1) .GE. TM ) THEN
              TMLAST = DBLTIM(NBL-1)
            ELSE
              TMLAST = TM
            ENDIF
          ELSE
            TMLAST = TM
          ENDIF
!
!---  Loop through grid domain element-by-element;
!     accumulate dilute solute mass  ---
!
          DO 200 NS = 1,NSOLU
            GLBL = 0.
            DO 100 N = 1,NFLD
              GLBL = GLBL +C(N,NS)*VOL(N)
  100       CONTINUE
!
!---  Decay-correct the global mass balance for radioactives  ---
!
            IF( HLF(NS) .GT. 0 ) THEN
              GLBL = GLBL*EXP(-(LOG(2.0)/HLF(NS))*(DBLTIM(NBL)-TMLAST))
            ENDIF
!
!---  Report global balance of mass/activity to 'balance' file  ---
!
            NCH = LEN_TRIM(SOLUT(NS))
            WRITE(UNIT=IBL,FMT=901) TM,',s,',SOLUT(NS)(1:NCH),',',
     &        GLBL,','
  901       FORMAT(1PE19.12,3A,1PE19.12,A)
!
  200     CONTINUE
        ENDIF
  400 CONTINUE
  401 FORMAT(1PE19.12,3A,1PE19.12,A)
      CLOSE(IBL)
!
!---  Record zero additional mass/activity to 'release' file  ---
!
      WRITE(UNIT=IRL,FMT=501) TM,',s,',0.0,',m,',0.0,','
      WRITE(UNIT=IRL,FMT=501) 1.5E11,',s,',0.0,',m,',0.0,','
  501 FORMAT(1PE19.12,A,1PE19.12,A,1PE19.12,A)
      CLOSE(UNIT=IRL)
!
!---  Record zero additional mass/activity to 'airborne' file  ---
!
      IF( ICODE .EQ. 2 ) THEN
        WRITE(UNIT=IRG,FMT=502) TM,',s,',0.0,',m,',0.0,','
        WRITE(UNIT=IRG,FMT=502) 1.5E11,',s,',0.0,',m,',0.0,','
  502   FORMAT(1PE19.12,A,1PE19.12,A,1PE19.12,A)
        CLOSE(UNIT=IRG)
      ENDIF
!
!---  Record zero additional quantities remediated to 
!     'remediate' file  ---
!
      DO 600 NRE = 1,INMREM
        NCH1 = LEN_TRIM(CRMSIT(NRE))
        NCH2 = LEN_TRIM(SOLUT(IRMSOL(NRE)))
        NCH3 = LEN_TRIM(CRMWFM(NRE))
        IF( DRMTIM(NRE) .GE. TM ) THEN
          WRITE (UNIT=IVD,FMT=601) DRMTIM(NRE),',s,',
     &      CRMSIT(NRE)(1:NCH1),',',
     &      SOLUT(IRMSOL(NRE))(1:NCH2),',',0.0,',,',
     &      CRMWFM(NRE)(1:NCH3),','
        ENDIF
  600 CONTINUE
  601 FORMAT(1PE19.12,5A,1PE19.12,3A)
      CLOSE(UNIT=IVD)
!
!---  End of SACOUT group  ---
!
  900 ISUB_LOG = ISUB_LOG-1
      RETURN
      END

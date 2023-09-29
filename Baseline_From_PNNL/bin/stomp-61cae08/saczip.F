!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SACZIP( IDUN )
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
!     After reading input file, this subroutine is called to check if
!     there is effectively zero analyte source for the simulation, and
!     if so, writes a single zero-release to the 'release' file and
!     then terminates the program without full transport simulation to
!     improve overall SAC execution efficiency.
!
!     Call Argument IDUN:
!       * Set to zero if nonzero analyte source is found;
!       * Set to one if no analyte source is found, signaling STOMP
!         to skip over flow and transport solution.
!
!----------------------SAC Rev. 0 Limitations--------------------------!
!
!     SACZIP is only for detecting mass/analyte sources (not rates),
!     and only in Operation Mode 1 (Water).
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, Battelle PNNL, November 2000.
!     Last Modified by WE Nichols, Battelle PNNL, July 14, 2004.
!
!     SAC Modification by WE Nichols, PNNL, July 14, 2004.
!     saczip.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
#ifdef sac
      USE SACIN
#endif
#ifdef sac
      USE SACEV
#endif
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
      SUB_LOG(ISUB_LOG) = '/SACZIP'
!
!---  By default, this routine will not direct skip over solution  ---
!
      IDUN = 0
!
!---  Only force termination for a transport solution  ---
!
      IF( IEQC .GT. 1 ) THEN
!
!---  Check for presence of a nonzero aqueous mass solute source  ---
!
        DO 300 J = 3,(3+NSOLU)
          DO 200 L = 1,NSR
            IF( ISRT(L) .LT. 0 ) THEN
              DO 100 K = 1,ISRM(L)
                IF( (SRC(J,K,L) .GT. 1.0D-20) ) GOTO 900
  100         CONTINUE
            ENDIF
  200     CONTINUE
  300   CONTINUE
!
!---  Record zero mass/activity balances to 'balance' file  ---
!
        DO 400 NBL = 1,INMBLT
          DO 350 NS = 1,NSOLU
            NCH = LEN_TRIM(SOLUT(NS))
            WRITE(UNIT=IBL,FMT=401) DBLTIM(NBL),',s,',SOLUT(NS)(1:NCH),
     &        ',',0.0,','
  350     CONTINUE
  400   CONTINUE
  401   FORMAT(1PE19.12,3A,1PE19.12,A)
        CLOSE(IBL)
!
!---  Record zero mass/activity to 'release' file  ---
!
        WRITE(UNIT=IRL,FMT=501) 0.0E00,',s,',0.0,',m,',0.0,','
        WRITE(UNIT=IRL,FMT=501) 1.5E11,',s,',0.0,',m,',0.0,','
  501   FORMAT(1PE19.12,A,1PE19.12,A,1PE19.12,A)
        CLOSE(UNIT=IRL)
!
!---  Record zero mass/activity to 'airborne' file  ---
!
        IF( ICODE .EQ. 2 ) THEN
          WRITE(UNIT=IRG,FMT=501) 0.0E00,',s,',0.0,',m,',0.0,','
          WRITE(UNIT=IRG,FMT=501) 1.5E11,',s,',0.0,',m,',0.0,','
  502     FORMAT(1PE19.12,A,1PE19.12,A,1PE19.12,A)
          CLOSE(UNIT=IRG)
        ENDIF
!
!---  Record zero quantities remediated to 'remediate' file  ---
!
        DO 600 NRE = 1,INMREM
          NCH1 = LEN_TRIM(CRMSIT(NRE))
          NCH2 = LEN_TRIM(SOLUT(IRMSOL(NRE)))
          NCH3 = LEN_TRIM(CRMWFM(NRE))
          WRITE (UNIT=IVD,FMT=601) DRMTIM(NRE),',s,',
     &      CRMSIT(NRE)(1:NCH1),',',SOLUT(IRMSOL(NRE))(1:NCH2),',',
     &      0.0,',,',CRMWFM(NRE)(1:NCH3),','
  600   CONTINUE
  601   FORMAT(1PE19.12,5A,1PE19.12,3A)
        CLOSE(UNIT=IVD)
!
!---  Set IDUN to indicate no need to simulate flow and transport  ---
!
        IDUN = 1
        WRITE(IWR,'(/,A)')  'Zero analyte source to simulate; skipping-'
        WRITE(ISC,'(/,A)') ' Zero analyte source to simulate; skipping-'
      ENDIF
!
!---  End of SACZIP group  ---
!
  900 ISUB_LOG = ISUB_LOG-1
      RETURN
      END

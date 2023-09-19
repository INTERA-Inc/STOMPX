!-------------------------Subroutine-----------------------------------!
!
      SUBROUTINE RDSREM
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
!     Read input file for SAC remediation information.
!
!     This information is only of use when running STOMP within the
!     the System Assessment Capability (SAC) framework.
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, Battelle PNNL, May 2000.
!     Last Modified by WE Nichols, Battelle PNNL, July 14, 2004.
!
!     SAC Modification by WE Nichols, PNNL, July 14, 2004.
!     rdsrem.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      LOGICAL FEXIST
      CHARACTER*512 ADUM,CHDUM,UNTS
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSREM'
!
!---  Write card information to ouput file  ---
!
      CARD = 'SAC Remediation Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read the number of SAC remediation events  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of SAC Remediation Events: '
      ISTART = 1
      CALL RDINT(ISTART,ICOMMA,CHDUM,INMREM)
      WRITE (IWR,FMT='(/A,I5)') VARB(1:IVR),INMREM
!
!---  Check for parameter exceedance  ---
!
      IF( INMREM .GT. LREM ) THEN
        INDX = 5
        CHMSG = 'Number of SAC Remediation Events > Parameter LREM'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Loop Over Number of Remediation Events  ---
!
      DO 100 M = 1,INMREM
!
!---  Write SAC Remediation Event Header  ---
!
      VARB = 'SAC Remediation Event # '
      WRITE (IWR,FMT='(/2X,A,I3)') VARB(1:IVR),M
!
!---  Read SAC Remediation Event Time  ---
!
        VARB = 'Remediation Time'
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        IDFLT = 1
        ISTART = 1
        UNTS = 's'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DRMTIM(M))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',DRMTIM(M)
        INDX = 0
        IUNS = 1
        CALL RDUNIT(UNTS,DRMTIM(M),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',DRMTIM(M),', s)'
!
!---  Read Solute Name  ---
!
        VARB = 'Remediated Solute'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---  Search Known Solutes for Matching Name  ---
!
        DO 50 NSL = 1,NSOLU
          IF( ADUM.EQ.SOLUT(NSL) ) GOTO 60
   50   CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute Name: '//ADUM
          CALL WRMSGS( INDX )
   60   CONTINUE
        IRMSOL(M) = NSL
        WRITE(IWR,'(4X,3A)') VARB(1:IVR),': ',SOLUT(NSL)(1:NCH)
!
!---  Read SAC Remediation Top of Depth Range  ---
!
        VARB = 'Remediation Upper Depth'
        UNTS = 'm'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DRMDUP(M))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',DRMDUP(M)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DRMDUP(M),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',DRMDUP(M),', m)'
!
!---  Read SAC Remediation Depth  ---
!
        VARB = 'Remediation Lower Depth'
        UNTS = 'm'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DRMDLO(M))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',DRMDLO(M)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DRMDLO(M),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',DRMDLO(M),', m)'
!
!---  Read SAC Remediation Fraction (efficiency) and
!     Check for realistic remediation efficencies  ---
!
        VARB = 'Remediation Efficiency'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DRMEFF(M))
        IF( DRMEFF(M) .LT. 0.0 ) DRMEFF(M) = 0.0
        IF( DRMEFF(M) .GT. 1.0 ) DRMEFF(M) = 1.0
        WRITE (IWR,FMT='(4X,2A,1PE11.4)') VARB(1:22),': ',DRMEFF(M)
!
!---  Truncate remediation efficencies greater than 100 percent  ---
!
        IF( DRMEFF(M) .GT. 1.0 ) DRMEFF(M) = 1.0
!
!---  Read SAC Destination Site Name for Remediated Mass  ---
!
        VARB = 'SAC Destination Site'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,CRMSIT(M))
        WRITE (IWR,FMT='(4X,3A)') VARB(1:IVR),': ',CRMSIT(M)(1:NCH)
!
!---  Read optional SAC Destination Site waste form (if present)  ---
!
        IF( INDEX(CHDUM((ICOMMA+1):LEN(CHDUM)),',') .GT. 0 ) THEN
          VARB = 'SAC Destination Waste Form'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,CRMWFM(M))
          IF( CRMWFM(M) .EQ. 'null') CRMWFM(M) = 'soil'          
          WRITE (IWR,FMT='(4X,3A)') VARB(1:IVR),':',CRMWFM(M)(1:NCH)
        ELSE
          CRMWFM(M) = 'SOIL'
        ENDIF
  100 CONTINUE
!
!---  Delete any existing 'remediate' file  ---
!
        INQUIRE(FILE=FNVD,EXIST=FEXIST)
        IF( FEXIST ) THEN
          OPEN(UNIT=97,FILE=FNVD,STATUS='unknown')
          CLOSE(UNIT=97,STATUS='delete')
        ENDIF
!
!---  If transport equation is being solved, start 'remediate' file  ---
!
      IF( IEQC .NE. 0 ) THEN
        OPEN(UNIT=IVD, FILE=FNVD, FORM='FORMATTED')
        WRITE (IVD,FMT='(A)') '~STOMP Remediation Transfers Report'
!
!---  Reprint SAC header information, if present in 'input' file  ---
!
        REWIND(UNIT=IRD)
  200   READ(UNIT=IRD,FMT='(A)',END=210) CHDUM
        ADUM = CHDUM
        CALL LCASE( ADUM )
        IF( ADUM(1:4) .EQ. '#sac' ) THEN
          WRITE(UNIT=IVD,FMT='(A)') TRIM(CHDUM)
        ENDIF
        GOTO 200
  210   CONTINUE
!
!---  Write number of remediation transfers to be recorded  ---
!
        WRITE(UNIT=IVD,FMT='(I0,A)') INMREM, ','
!
      ENDIF
!
!---  End of RDSREM group ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

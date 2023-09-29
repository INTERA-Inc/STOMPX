!-------------------------Subroutine-----------------------------------!
!
      SUBROUTINE RDSREL
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
!     Read input file for SAC release plane information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, Battelle PNNL, May 2000.
!     Last Modified by WE Nichols, Battelle PNNL, September 8, 2003.
!
!     SAC Modification by WE Nichols, PNNL, June 30, 2003.
!     rdsrel.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
      SUB_LOG(ISUB_LOG) = '/RDSREL'
!
!---  Write card information to ouput file  ---
!
      CARD = 'SAC Release Plane'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read the number of SAC release plane time-dependent changes  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of SAC Release Plane Time Change, K-Node Indices: '
      ISTART = 1
      CALL RDINT(ISTART,ICOMMA,CHDUM,INMREL)
      WRITE (IWR,FMT='(/A45,I5)') VARB,INMREL
!
!---  Check for parameter exceedance  ---
!
      IF( INMREL .GT. LREL ) THEN
        INDX = 5
        CHMSG = 'Number of SAC Release Plane Changes > Parameter LREL'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Read SAC release plane times and K-node index  ---
!
      DO 100 M = 1,INMREL
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        IDFLT = 1
        ISTART = 1
        UNTS = 's'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DRLTIM(M))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        INDX = 0
        IUNS = 1
        CALL RDUNIT(UNTS,DRLTIM(M),INDX)
        CALL RDINT(ISTART,ICOMMA,CHDUM,IRLK(M))
        WRITE (IWR,FMT='(2X,1PE11.4,1X,A,1X,A,I4)')
     &    DRLTIM(M),UNTS(1:NCH),' at K-Node Index = ',IRLK(M)
  100 CONTINUE
!
!---  Initialize the SAC release plane  ---
!
      IF( DRLTIM(1) .LT. 1. ) THEN
        KEXIT = IRLK(1)
      ELSE
        KEXIT = 1
      ENDIF
!
!---  Delete any existing 'release', 'recharge', or 'airborne files  ---
!
      INQUIRE(FILE=FNRL,EXIST=FEXIST)
      IF( FEXIST ) THEN
        OPEN(UNIT=97,FILE=FNRL,STATUS='unknown')
        CLOSE(UNIT=97,STATUS='delete')
      ENDIF
      INQUIRE(FILE=FNRC,EXIST=FEXIST)
      IF( FEXIST ) THEN
        OPEN(UNIT=97,FILE=FNRC,STATUS='unknown')
        CLOSE(UNIT=97,STATUS='delete')
      ENDIF
      INQUIRE(FILE=FNRG,EXIST=FEXIST)
      IF( FEXIST ) THEN
        OPEN(UNIT=97,FILE=FNRG,STATUS='unknown')
        CLOSE(UNIT=97,STATUS='delete')
      ENDIF
!
!---  If transport equation is being solved, open 'release'
!     and 'recharge' files, plus 'airborne' file if this
!     is STOMP Mode 2 only  ---
!
      IF( IEQC .NE. 0 ) THEN
        OPEN(UNIT=IRL, FILE=FNRL, FORM='FORMATTED')
        OPEN(UNIT=IRC, FILE=FNRC, FORM='FORMATTED')
        WRITE(*,*) ' ICODE=',ICODE
        IF( ICODE .EQ. 2 ) OPEN(UNIT=IRG, FILE=FNRG, FORM='FORMATTED')
!
!---  Write header to 'release' file  ---
!
        WRITE(UNIT=IRL,FMT='(A)')
     &    '~STOMP SAC Mass/Activity Releases to Aquifer'
        WRITE(UNIT=IRC,FMT='(A)')
     &    '~STOMP SAC Artifical Recharge to Aquifer'
        IF( ICODE .EQ. 2 ) WRITE(UNIT=IRG,FMT='(A)')
     &    '~STOMP SAC Mass/Activity Releases to Atmosphere'
!
!---  Reprint SAC header information, if present in 'input' file  ---
!
        REWIND(UNIT=IRD)
  200   READ(UNIT=IRD,FMT='(A)',END=210) CHDUM
        ADUM = CHDUM
        CALL LCASE( ADUM )
        IF( ADUM(1:4) .EQ. '#sac' ) THEN
          WRITE(UNIT=IRL,FMT='(A)') CHDUM
          WRITE(UNIT=IRC,FMT='(A)') CHDUM
          IF( ICODE .EQ. 2 ) WRITE(UNIT=IRG,FMT='(A)') CHDUM
        ENDIF
        GOTO 200
  210   CONTINUE
      ENDIF
!
!---  End of RDSREL group ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

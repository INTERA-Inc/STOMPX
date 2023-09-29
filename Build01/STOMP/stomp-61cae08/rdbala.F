!-------------------------Subroutine-----------------------------------!
!
      SUBROUTINE RDBALA
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
!     Read input file for global balance check times information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, Battelle PNNL, May 2000.
!     Last Modified by WE Nichols, Battelle PNNL, September 15, 2003.
!
!     SAC Modification by WE Nichols, Battelle PNNL, June 12, 2003.
!     rdbala.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0'
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
      CHARACTER*12 FORM1
      CHARACTER*512 ADUM,CHDUM,UNTS
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDBLNC'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Mass/Activity Balance Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read the number of mass balance times  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Global Balance Check Times: '
      ISTART = 1
      CALL RDINT(ISTART,ICOMMA,CHDUM,INMBLT)
      WRITE (IWR,FMT='(/A45,I5)') VARB,INMBLT
!
!---  Check for parameter exceedance  ---
!
      IF( INMBLT .GT. LBAL ) THEN
        INDX = 5
        CHMSG = 'Number of Balance Checks > Parameter LBAL'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Read Mass Balance Times  ---
!
      DO 100 M = 1,INMBLT
        VARB = 'Check Time'
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        IDFLT = 1
        ISTART = 1
        UNTS = 's'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DBLTIM(M))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',DBLTIM(M)
        INDX = 0
        IUNS = 1
        CALL RDUNIT(UNTS,DBLTIM(M),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',DBLTIM(M),', s)'
  100 CONTINUE
!
!---  Open 'balance' reporting file and write header information  ---
!
      INQUIRE(FILE=FNBL,EXIST=FEXIST)
      IF( FEXIST ) THEN
        OPEN(UNIT=97,FILE=FNBL,STATUS='unknown')
        CLOSE(UNIT=97,STATUS='delete')
      ENDIF
      OPEN(UNIT=IBL,FILE=FNBL,FORM='FORMATTED')
!
!---  Initialize format specifer for internal write  ---
!
      FORM1 = '(I1,A1)'
!
!---  Write header to 'balance' file  ---
!
      WRITE(UNIT=IBL,FMT='(A)')
     &  '~STOMP Global Mass/Activity Balance Report'
!
!---  Reprint SAC header information, if present in 'input' file  ---
!
      REWIND(IRD)
  200 READ(UNIT=IRD,FMT='(A)',END=210) CHDUM
      ADUM = CHDUM
      CALL LCASE( ADUM )
      IF( ADUM(1:4) .EQ. '#sac' ) WRITE(UNIT=IBL,FMT='(A)') TRIM(CHDUM)
      GOTO 200
  210 CONTINUE
!
!---  Write number of balance checks to be reported  ---
!
      WRITE(UNIT=IBL,FMT='(I0,A)') INMBLT, ','
!
!---  End of RDBALA group ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END


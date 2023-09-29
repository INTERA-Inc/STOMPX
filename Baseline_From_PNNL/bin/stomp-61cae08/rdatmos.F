!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDATMOS
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
!     Read input file for atmospheric conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 June 2002.
!     Last Modified by MD White, PNNL, 6 June 2002.
!     rdatmos.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PLT_ATM
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,CDUM,FDUM,FMDUM,UNTS
      CHARACTER*512 CHDUM
      LOGICAL FCHK
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDATMOS'
!
!---  Read atmospheric data start time  ---
!
      CARD = 'Atmospheric Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Atmospheric Start Time: Month'
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      VARB = 'Atmospheric Start Time: Day'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IDAY)
      IF( INDEX(ADUM(1:),'jan').NE.0 .OR.
     &  INDEX(ADUM(1:),'01').NE.0 ) THEN
        ATMST = ATMST + 0.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'feb').NE.0 .OR.
     &  INDEX(ADUM(1:),'02').NE.0 ) THEN
        ATMST = ATMST + 31.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.28 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'mar').NE.0 .OR.
     &  INDEX(ADUM(1:),'03').NE.0 ) THEN
        ATMST = ATMST + 59.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'apr').NE.0 .OR.
     &  INDEX(ADUM(1:),'04').NE.0 ) THEN
        ATMST = ATMST + 90.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.30 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'may').NE.0 .OR.
     &  INDEX(ADUM(1:),'05').NE.0 ) THEN
        ATMST = ATMST + 120.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'jun').NE.0 .OR.
     &  INDEX(ADUM(1:),'06').NE.0 ) THEN
        ATMST = ATMST + 151.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.30 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'jul').NE.0 .OR.
     &  INDEX(ADUM(1:),'07').NE.0 ) THEN
        ATMST = ATMST + 181.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'aug').NE.0 .OR.
     &  INDEX(ADUM(1:),'08').NE.0 ) THEN
        ATMST = ATMST + 212.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'sep').NE.0 .OR.
     &  INDEX(ADUM(1:),'09').NE.0 ) THEN
        ATMST = ATMST + 243.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.30 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'oct').NE.0 .OR.
     &  INDEX(ADUM(1:),'10').NE.0 ) THEN
        ATMST = ATMST + 273.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'nov').NE.0 .OR.
     &  INDEX(ADUM(1:),'11').NE.0 ) THEN
        ATMST = ATMST + 304.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.30 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'dec').NE.0 .OR.
     &  INDEX(ADUM(1:),'12').NE.0 ) THEN
        ATMST = ATMST + 334.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Atmospheric Start Month: '
     &    // ADUM(1:NCHA)
        CALL WRMSGS( INDX )
      ENDIF
      VARB = 'Atmospheric Start Time: Year'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IYEAR)
      IF( IYEAR.GE.1900 ) THEN
        REALX = REAL(IYEAR-1900)
        ATMST = ATMST + REALX*365.D+0*8.640D+4
      ELSE
        INDX = 4
        IMSG = IYEAR
        CHMSG = 'Out-of-Range Atmospheric Start Year: '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Correct for leap years  ---
!
      DO 5 M = 1901,IYEAR
        IF( MOD(M,1000).EQ.0 ) THEN
          IF( MOD(M,400).EQ.0 ) THEN
            ATMST = ATMST + 8.64D+4
          ENDIF
        ELSEIF( MOD(M,4).EQ.0 ) THEN
          ATMST = ATMST + 8.64D+4
        ENDIF
    5 CONTINUE
!
!---  Read atmospheric start time in military format  ---
!
      VARB = 'Atmospheric Start Time: Time'
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      READ( ADUM(1:2),'(I2)' ) IHOUR
      IF( IHOUR.GE.0 .OR. IHOUR.LT.24 ) THEN
        REALX = REAL(IHOUR)
        ATMST = ATMST + REALX*3600.D+0
      ELSE
        INDX = 4
        IMSG = IHOUR
        CHMSG = 'Out-of-Range Atmospheric Start Hour: '
        CALL WRMSGS( INDX )
      ENDIF
      READ( ADUM(4:5),'(I2)' ) IMIN
      IF( IMIN.GE.0 .OR. IMIN.LT.60 ) THEN
        REALX = REAL(IMIN)
        ATMST = ATMST + REALX*60.D+0
      ELSE
        INDX = 4
        IMSG = IMIN
        CHMSG = 'Out-of-Range Atmospheric Start Minute: '
        CALL WRMSGS( INDX )
      ENDIF
      READ( ADUM(7:8),'(I2)' ) ISEC
      IF( ISEC.GE.0 .OR. ISEC.LT.60 ) THEN
        REALX = REAL(ISEC)
        ATMST = ATMST + REALX
      ELSE
        INDX = 4
        IMSG = ISEC
        CHMSG = 'Out-of-Range Atmospheric Start Second: '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Read measurement height for wind speed  ---
!
      VARB = 'Wind Speed Measurement Height'
      UNTS = 'm'
      IUNM = 1
      IDFLT = 1
      ATMC(1) = 1.D+0
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(1))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',ATMC(1)
      INDX = 0
      CALL RDUNIT(UNTS,ATMC(1),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(1),', m)'
!
!---  Read measurement height for air temperature and
!     relative humidity  ---
!
      VARB = 'Air Temperature/Relative Humidity Measurement Height'
      UNTS = 'm'
      IUNM = 1
      IDFLT = 1
      ATMC(2) = 1.D+0
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(2))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',ATMC(2)
      INDX = 0
      CALL RDUNIT(UNTS,ATMC(2),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(2),', m)'
!
!---  Read local longitude  ---
!
      VARB = 'Local Longitude'
      UNTS = 'deg'
      IDFLT = 1
      ATMC(3) = 0.D+0
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(3))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',ATMC(3)
      INDX = 0
      CALL RDUNIT(UNTS,ATMC(3),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(3),', radians)'
!
!---  Read local latitude  ---
!
      VARB = 'Local Latitude'
      UNTS = 'deg'
      IDFLT = 1
      ATMC(4) = 0.D+0
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(4))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',ATMC(4)
      INDX = 0
      CALL RDUNIT(UNTS,ATMC(4),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(4),', radians)'
!
!---  Read local meridian, where meridians for the United States
!     are as follows:
!     Eastern Standard Time = 75 deg
!     Central Standard Time = 90 deg
!     Mountain Standard Time = 105 deg
!     Pacific Standard Time = 120 deg  ---
!
      VARB = 'Local Meridian'
      UNTS = 'deg'
      IDFLT = 1
      ATMC(5) = 120.D+0
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(5))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',ATMC(5)
      INDX = 0
      CALL RDUNIT(UNTS,ATMC(5),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(5),', radians)'
!
!---  Read roughness height for momentum transport  ---
!
      CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
      IF( INDX.EQ.1 ) THEN
        VARB = 'Roughness Height for Momentum Transport'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        ATMC(6) = 1.3D-3
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(6))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',ATMC(6)
        INDX = 0
        CALL RDUNIT(UNTS,ATMC(6),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(6),', m)'
      ENDIF
!
!---  Read roughness height for heat and mass transport  ---
!
      CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
      IF( INDX.EQ.1 ) THEN
        VARB = 'Roughness Height for Heat and Mass Transport'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        ATMC(7) = 1.3D-3
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(7))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',ATMC(7)
        INDX = 0
        CALL RDUNIT(UNTS,ATMC(7),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(7),', m)'
      ENDIF
!
!---  Read number of atmospheric condition times  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
!      VARB = 'Input Option [File, Integer]'
!      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
!!
!!---  Read atmospheric conditions from an external file  ---
!!
!      IF( INDEX(ADUM,'file').NE.0 ) THEN
!        VARB = 'External Atmospheric Conditions File Name'
!        CALL RDCHR(ISTART,ICOMMA,NCHF,CHDUM,FDUM)
!        INQUIRE( FILE=FDUM(1:NCHF), FORM=CDUM, EXIST=FCHK )
!        IF( .NOT.FCHK ) THEN
!          INDX = 4
!          CHMSG = 'Atmospheric conditions file does not exist: '
!     &      // FDUM(1:NCHF)
!          CALL WRMSGS( INDX )
!        ELSEIF( CDUM.EQ.'UNFORMATTED' ) THEN
!          INDX = 4
!          CHMSG = 'Atmospheric conditions file is unformatted: '
!     &      // FDUM(1:NCHF)
!          CALL WRMSGS( INDX )
!        END IF
!        OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD', FORM='FORMATTED')
!        WRITE(IWR,'(/,2A)') 'Atmospheric Conditions File: ',FDUM(1:NCHF)
!        READ(26,*,END=10) NATM_T
!        IF( NATM_T.LE.-3 ) THEN
!          IATM_C = 1
!          NATM_T = ABS(NATM_T)
!          WRITE(IWR,'(A)') 'Cyclic Atmospheric Conditions'
!        ELSEIF( NATM_T.GE.1 ) THEN
!          IATM_C = 0
!          WRITE(IWR,'(A)') 'Noncyclic Atmospheric Conditions'
!        ELSEIF( NATM_T.EQ.0 ) THEN
!          INDX = 2
!          CHMSG = 'No Atmospheric Condition Times'
!          CALL WRMSGS( INDX )
!        ELSE
!          INDX = 4
!          CHMSG = 'Number of Cyclic Atmospheric Conditions Times < 3'
!          CALL WRMSGS( INDX )
!        ENDIF
!        GOTO 20
!   10   CONTINUE
!        INDX = 4
!        CHMSG = 'Empty Atmospheric Conditions File'
!        CALL WRMSGS( INDX )
!   20   CONTINUE
!        NTM = 0
!   30   CONTINUE
!        NTM = NTM+1
!        IF( NATM_T.GT.LATM ) THEN
!          INDX = 5
!          CHMSG = 'Number of Atmospheric Condition Times > LATM'
!          CALL WRMSGS( INDX )
!        ENDIF
!        READ(26,*,END=40) (ATMOS(NTM,M),M=1,6)
!        IF( NTM.GE.NATM_T ) GOTO 50
!        GOTO 30
!   40   CONTINUE
!        INDX = 4
!        CHMSG = 'Atmospheric Conditions File Missing Lines'
!        CALL WRMSGS( INDX )
!   50   CONTINUE
!!
!!---  Read atmospheric conditions from input file  ---
!!
!      ELSE
!        ISTART = 1
        VARB = 'Number of Atmospheric Condition Times'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NATM_T)
        IF( NATM_T.LE.-3 ) THEN
          IATM_C = 1
          NATM_T = ABS(NATM_T)
          WRITE(IWR,'(A)') 'Cyclic Atmospheric Conditions'
        ELSEIF( NATM_T.GE.1 ) THEN
          IATM_C = 0
          WRITE(IWR,'(A)') 'Noncyclic Atmospheric Conditions'
        ELSEIF( NATM_T.EQ.0 ) THEN
          INDX = 2
          CHMSG = 'No Atmospheric Condition Times'
          CALL WRMSGS( INDX )
        ELSE
          INDX = 4
          CHMSG = 'Number of Cyclic Atmospheric Conditions Times < 3'
          CALL WRMSGS( INDX )
        ENDIF
        IF( NATM_T.GT.LATM ) THEN
          INDX = 5
          CHMSG = 'Number of Atmospheric Condition Times > LATM'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read number of atmospheric condition times and variables  ---
!
        ACTMO = -SMALL
        WRITE(IWR,'(A)') 'Atmospheric Condition Times and Variables:'
        DO 100 NTM = 1,NATM_T
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
!
!---      Check for external atmospheric condition time file  ---
!
          IF( NTM.EQ.1 ) CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,CDUM)
          IF( INDEX(CDUM(1:),'file').NE.0 ) THEN
            CDUM = ''
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            NCH = INDEX(FDUM,'  ')-1
            IRD = 32
!
!---        Check for external file  ---
!
            INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
            IF( .NOT.FCHK ) THEN
              INDX = 4
              CHMSG = 'Missing Atmospheric Condition File: ' // 
     &          FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ELSEIF( FDUM.EQ.'unformatted' ) THEN
              INDX = 4
              CHMSG = 'Atmospheric Condition File Format: ' // 
     &          FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
            OPEN(UNIT=32,FILE=FDUM(1:NCH),STATUS='OLD',FORM='FORMATTED')
            WRITE(IWR,'(/,2A)') 'Atmospheric Condition Time File: ',
     &        FDUM(1:NCH)
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
          ENDIF
!
!---      Read atmospheric condition time  ---
!
          ISTART = 1
          VARB = 'Atmospheric Condition Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &      ', ',UNTS(1:NCH),': ',ATMOS(NTM,1)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,ATMOS(NTM,1),INDX)
          IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &      ATMOS(NTM,1),', C)'
!
!---      Read atmospheric condition temperature  ---
!
          VARB = 'Atmospheric Condition Temperature'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,2))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &      ', ',UNTS(1:NCH),': ',ATMOS(NTM,2)
          INDX = 0
          IUNK = 1
          CALL RDUNIT(UNTS,ATMOS(NTM,2),INDX)
          IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &      ATMOS(NTM,2),', C)'
!
!---      Read atmospheric condition pressure  ---
!
          VARB = 'Atmospheric Condition Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,3))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &      ', ',UNTS(1:NCH),': ',ATMOS(NTM,3)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,ATMOS(NTM,3),INDX)
          IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &      ATMOS(NTM,3),', Pa)'
!
!---      Read atmospheric condition water-vapor relative
!         humidity  ---
!
          VARB = 'Atmospheric Condition Water-Vapor Relative Humidity'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,4))
          IF( IRD.EQ.21 ) WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),
     &      ': ',ATMOS(NTM,4)
!
!---      Read atmospheric condition net solar radiation  ---
!
          VARB = 'Atmospheric Condition Net Solar Radiation'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,5))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &      ', ',UNTS(1:NCH),': ',ATMOS(NTM,5)
          INDX = 0
          IUNKG = 1
          IUNS = -3
          CALL RDUNIT(UNTS,ATMOS(NTM,5),INDX)
          IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &      ATMOS(NTM,5),', W/m^2)'
!
!---      Read atmospheric condition wind speed  ---
!
          VARB = 'Atmospheric Condition Wind Speed'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,6))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &      ', ',UNTS(1:NCH),': ',ATMOS(NTM,6)
          INDX = 0
          IUNM = 1
          IUNS = -1
          CALL RDUNIT(UNTS,ATMOS(NTM,6),INDX)
          IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &      ATMOS(NTM,6),', m/s)'
  100   CONTINUE
!
!---    Close boundary condition time file  ---
!
        IF( IRD.EQ.32 ) THEN
          IRD = 21
          CLOSE(32)
        ENDIF
!      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDATMOS group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTHER
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
!     Read input file for rock/soil thermal information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, December 1992.
!     Last Modified by MD White, PNNL, December 8, 1995.
!     rdther.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
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
      CHARACTER*64 ADUM,UNTS
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDTHER'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Rock/Soil Thermal Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Loop over the rock/soil thermal information lines  ---
!
      N = 0
      ISGRP = 0
   10 CONTINUE
      IJK = 0
      IF( N.GE.NROCK ) GOTO 500
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Rock/Soil Name'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---  IJK, KIJ, or JKI indexing  ---
!
      IF( INDEX(ADUM(1:),'indexing').NE.0 ) THEN
        IF( INDEX(ROCK(1)(1:),'indexing').EQ.0 ) THEN
          INDX = 4
          CHMSG = 'Indexing Option Not Declared ' // 
     &      'in Rock/Soil Zonation Card'
          CALL WRMSGS( INDX )
        ENDIF
        IF( INDEX(ADUM,'ijk').NE.0 ) THEN
          IJK = 1
        ELSEIF( INDEX(ADUM,'jki').NE.0 ) THEN
          IJK = 2
        ELSEIF( INDEX(ADUM,'kij').NE.0 ) THEN
          IJK = 3
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Indexing Option' // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        GOTO 220
      ENDIF
!
!---  Search known rock types for a matching type ---
!
      DO 100 M = 1, NROCK
        IF( ADUM .EQ. ROCK(M)) THEN
          IROCK = M
          GOTO 200
        ENDIF
  100 CONTINUE
!
!---  Search known scaling groups for a matching type ---
!
      IF( ISLC(19).EQ.1 ) THEN
        DO 110 M = 1,NSCALE
           IF( ADUM.EQ.SCALNM(M) ) THEN
              ISGRP = M
              IROCK = 1
              GOTO 200
           ENDIF
  110   CONTINUE
        INDX = 2
        CHMSG = 'Unrecognized Rock/Soil Type or Scaling Group: '
     &    // ADUM(1:NCH)
        CALL WRMSGS( INDX )
        GOTO 10
      ENDIF
      INDX = 2
      CHMSG = 'Unrecognized Rock/Soil Type: ' // ADUM(1:NCH)
      CALL WRMSGS( INDX )
      GOTO 10
  200 CONTINUE
!
!---  Loop over rock/soils within scaling group  ---
!
      IF( ISLC(19).EQ.1 .AND. ISGRP.NE.0 ) THEN
        DO 202 M = IROCK,NROCK
          IF( ISCALE(M).EQ.ISGRP ) THEN
            IROCK = M
            GOTO 204
          ENDIF
  202   CONTINUE
      ENDIF
  204 CONTINUE
!
!---  Write rock/soil name  ---
!
      WRITE (IWR,'(/,2A)') 'Rock/Soil Name: ',ROCK(IROCK)
      N = N + 1
  220 CONTINUE
!
!---  Effective thermal conductivity function  ---
!
      VARB = 'Thermal Conductivity Function'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
      IF( INDEX(ADUM(1:),'constant').NE.0 ) THEN
        ITHKX = 1
        WRITE(IWR,'(A)') 'Constant Thermal Conductivity Model'
      ELSEIF( INDEX(ADUM(1:),'parallel').NE.0 ) THEN
        ITHKX = 2
        WRITE(IWR,'(A)') 'Parallel Thermal Conductivity Model'
      ELSEIF( INDEX(ADUM(1:),'linear').NE.0 ) THEN
        ITHKX = 3
        WRITE(IWR,'(A)') 'Linear Thermal Conductivity Model'
      ELSEIF( INDEX(ADUM(1:),'somerton').NE.0 ) THEN
        ITHKX = 4
        WRITE(IWR,'(A)') 'Somerton Thermal Conductivity Model'
      ELSEIF( INDEX(ADUM(1:),'campbell').NE.0 ) THEN
        ITHKX = 5
        WRITE(IWR,'(A)') 'Campbell Thermal Conductivity Model'
      ELSEIF( INDEX(ADUM(1:),'jame and norium').NE.0 ) THEN
        ITHKX = 6
        WRITE(IWR,'(A)') 'Jame and Norium Thermal Conductivity Model'
        GOTO 400
      ELSEIF( INDEX(ADUM(1:),'cass').NE.0 ) THEN
        WRITE(IWR,'(A)') 'Cass Thermal Conductivity Model'
        ITHKX = 7
      ELSEIF( INDEX(ADUM(1:),'hydrate composite').NE.0 .OR.
     &  INDEX(ADUM(1:),'kneafsey').NE.0 ) THEN
        WRITE(IWR,'(A)') 'Kneafsey Hydrate Composite ' // 
     &    'Thermal Conductivity Model'
        ITHKX = 8
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Thermal Conductivity Function: '
     &    // ADUM(1:NCH)
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Load model indices  ---
!
      IF( IJK.GT.0 ) THEN
        DO 230 NC = 1,NFLD
          ITHK(IZ(NC)) = ITHKX
  230   CONTINUE
        N = N + NFLD
      ELSE
        ITHK(IROCK) = ITHKX
      ENDIF
!
!---  Hydrate codes  ---
!
      IF( IOM.GE.36 .AND. IOM.LE.39 ) THEN
        IF( ITHKX.NE.2 .AND. ITHKX.NE.4 .AND. ITHKX.NE.8 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Thermal Conductivity Function: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  STOMP-CO2 and STOMP-CO2E  ---
!
      IF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
        IF( ITHKX.NE.2 .AND. ITHKX.NE.4 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Thermal Conductivity Function: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  STOMP-SEQ  ---
!
      IF( IOM.EQ.34 ) THEN
        IF( ITHKX.NE.2 .AND. ITHKX.NE.4 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Thermal Conductivity Function: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  STOMP-COMP   ---
!
      IF( IOM.EQ.40 ) THEN
        IF( ITHKX.NE.2 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Thermal Conductivity Function: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  STOMP-EOR and STOMP-GT  ---
!
      IF( IOM.EQ.43 .OR. IOM.EQ.3 ) THEN
        IF( ITHKX.NE.2 .AND. ITHKX.NE.4) THEN
          INDX = 4
          CHMSG = 'Unrecognized Thermal Conductivity Function: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Hydrate composite thermal conductivity model  ---
!
      IF( ITHKX.EQ.8 ) THEN
        IDFLT = 1
        VARB = 'Dry Thermal Conductivity'
        UNTS = 'w/m k'
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        INDX = 1
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(1,IROCK))
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',THKS(1,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,THKS(1,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(1,IROCK),', W/m K)'
        ENDIF
        IDFLT = 1
        VARB = 'Wet Thermal Conductivity'
        UNTS = 'w/m k'
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        INDX = 2
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(2,IROCK))
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',THKS(2,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,THKS(2,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(2,IROCK),', W/m K)'
        ENDIF
        GOTO 400
!
!---  Campbell effective thermal conductivity model  ---
!
      ELSEIF( ITHKX.EQ.5 ) THEN
!
!---    Load default values  ---
!
        IF( IJK.GT.0 ) THEN
          DO 300 N = 1,NFLD
            THKS(1,IZ(N)) = 7.34D-1
            THKS(2,IZ(N)) = 1.45D+0
            THKS(3,IZ(N)) = 2.01D+0
            THKS(4,IZ(N)) = 2.04D-1
            THKS(5,IZ(N)) = 4.D+0
  300     CONTINUE
        ELSE
          THKS(1,IROCK) = 7.34D-1
          THKS(2,IROCK) = 1.45D+0
          THKS(3,IROCK) = 2.01D+0
          THKS(4,IROCK) = 2.04D-1
          THKS(5,IROCK) = 4.D+0
        ENDIF
        IDFLT = 1
        VARB = 'Parameter a'
        UNTS = 'w/m k'
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        INDX = 1
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(1,IROCK))
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',THKS(1,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,THKS(1,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(1,IROCK),', W/m K)'
        ENDIF
        IDFLT = 1
        VARB = 'Parameter b'
        UNTS = 'w/m k'
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        INDX = 2
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(2,IROCK))
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',THKS(2,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,THKS(2,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(2,IROCK),', W/m K)'
        ENDIF
        IDFLT = 1
        VARB = 'Parameter c'
        UNTS = 'null'
        INDX = 3
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(3,IROCK))
          WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',THKS(3,IROCK)
        ENDIF
        IDFLT = 1
        VARB = 'Parameter d'
        UNTS = 'w/m k'
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        INDX = 4
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(4,IROCK))
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',THKS(4,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,THKS(4,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(4,IROCK),', W/m K)'
        ENDIF
        IDFLT = 1
        VARB = 'Parameter e'
        UNTS = 'null'
        INDX = 5
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(5,IROCK))
          WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',THKS(5,IROCK)
        ENDIF
        GOTO 400
      ELSEIF( ITHKX.EQ.7 ) THEN
!
!---    Load default values  ---
!
        IF( IJK.GT.0 ) THEN
          DO 310 N = 1,NFLD
            THKS(1,IZ(N)) = 0.6D+0
            THKS(2,IZ(N)) = 0.7D+0
            THKS(3,IZ(N)) = 8.0D+0
            THKS(4,IZ(N)) = 0.26D+0
            THKS(5,IZ(N)) = 3.0D+0
  310     CONTINUE
        ELSE
          THKS(1,IROCK) = 0.6D+0
          THKS(2,IROCK) = 0.7D+0
          THKS(3,IROCK) = 8.0D+0
          THKS(4,IROCK) = 0.26D+0
          THKS(5,IROCK) = 3.0D+0
        ENDIF
        IDFLT = 1
        VARB = 'Parameter a'
        UNTS = 'w/m k'
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        INDX = 1
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(1,IROCK))
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',THKS(1,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,THKS(1,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(1,IROCK),', W/m K)'
        ENDIF
        IDFLT = 1
        VARB = 'Parameter b'
        UNTS = 'w/m k'
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        INDX = 2
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(2,IROCK))
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',THKS(2,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,THKS(2,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(2,IROCK),', W/m K)'
        ENDIF
        IDFLT = 1
        VARB = 'Parameter c'
        UNTS = 'null'
        INDX = 3
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(3,IROCK))
          WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',THKS(3,IROCK)
        ENDIF
        IDFLT = 1
        VARB = 'Parameter d'
        UNTS = 'w/m k'
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        INDX = 4
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(4,IROCK))
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',THKS(4,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,THKS(4,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(4,IROCK),', W/m K)'
        ENDIF
        IDFLT = 1
        VARB = 'Parameter e'
        UNTS = 'null'
        INDX = 5
        LNDX = 9
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(5,IROCK))
          WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',THKS(5,IROCK)
        ENDIF
        GOTO 400
      ENDIF
      IDFLT = 1
      VARB = 'X-Direction Thermal Conductivity'
      UNTS = 'w/m k'
      IUNM = 1
      IUNKG = 1
      IUNS = -3
      IUNK = -1
      INDX = 1
      LNDX = 9
      IF( IJK.GT.0 ) THEN
        CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
      ELSE
        IDFLT = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(1,IROCK))
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',THKS(1,IROCK)
        INDX = 0
        CALL RDUNIT(UNTS,THKS(1,IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(1,IROCK),', W/m K)'
        IF( IFLD.GT.1 .AND. THKS(1,IROCK)/EPSL.LT.EPSL ) THEN
          INDX = 4
          CHMSG = 'Zero or Negative ' // VARB(1:IVR) //
     &      ' for ' // ROCK(IROCK)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IDFLT = 1
      VARB = 'Y-Direction Thermal Conductivity'
      UNTS = 'w/m k'
      IUNM = 1
      IUNKG = 1
      IUNS = -3
      IUNK = -1
      INDX = 2
      LNDX = 9
      IF( IJK.GT.0 ) THEN
        CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
      ELSE
        IDFLT = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(2,IROCK))
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',THKS(2,IROCK)
        INDX = 0
        CALL RDUNIT(UNTS,THKS(2,IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(2,IROCK),', W/m K)'
        IF( JFLD.GT.1 .AND. THKS(2,IROCK)/EPSL.LT.EPSL ) THEN
          INDX = 4
          CHMSG = 'Zero or Negative ' // VARB(1:IVR) //
     &      ' for ' // ROCK(IROCK)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IDFLT = 1
      VARB = 'Z-Direction Thermal Conductivity'
      UNTS = 'w/m k'
      IUNM = 1
      IUNKG = 1
      IUNS = -3
      IUNK = -1
      INDX = 3
      LNDX = 9
      IF( IJK.GT.0 ) THEN
        CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
      ELSE
        IDFLT = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(3,IROCK))
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',THKS(3,IROCK)
        INDX = 0
        CALL RDUNIT(UNTS,THKS(3,IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(3,IROCK),', W/m K)'
        IF( KFLD.GT.1 .AND. THKS(3,IROCK)/EPSL.LT.EPSL ) THEN
          INDX = 4
          CHMSG = 'Zero or Negative ' // VARB(1:IVR) //
     &      ' for ' // ROCK(IROCK)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IF( ITHKX.LE.2 ) GOTO 400
      IF( IEQW.EQ.0 ) GOTO 390
      IDFLT = 1
      VARB = 'X-Direction Water Saturated Thermal Conductivity'
      UNTS = 'w/m k'
      IUNM = 1
      IUNKG = 1
      IUNS = -3
      IUNK = -1
      INDX = 4
      LNDX = 9
      IF( IJK.GT.0 ) THEN
        CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
      ELSE
        IDFLT = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(4,IROCK))
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',THKS(4,IROCK)
        INDX = 0
        CALL RDUNIT(UNTS,THKS(4,IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(4,IROCK),', W/m K)'
        IF( IFLD.GT.1 .AND. THKS(4,IROCK)/EPSL.LT.EPSL ) THEN
          INDX = 4
          CHMSG = 'Zero or Negative ' // VARB(1:IVR) //
     &      ' for ' // ROCK(IROCK)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IDFLT = 1
      VARB = 'Y-Direction Water Saturated Thermal Conductivity'
      UNTS = 'w/m k'
      IUNM = 1
      IUNKG = 1
      IUNS = -3
      IUNK = -1
      INDX = 5
      LNDX = 9
      IF( IJK.GT.0 ) THEN
        CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
      ELSE
        IDFLT = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(5,IROCK))
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',THKS(5,IROCK)
        INDX = 0
        CALL RDUNIT(UNTS,THKS(5,IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(5,IROCK),', W/m K)'
        IF( JFLD.GT.1 .AND. THKS(5,IROCK)/EPSL.LT.EPSL ) THEN
          INDX = 4
          CHMSG = 'Zero or Negative ' // VARB(1:IVR) //
     &      ' for ' // ROCK(IROCK)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IDFLT = 1
      VARB = 'Z-Direction Water Saturated Thermal Conductivity'
      UNTS = 'w/m k'
      IUNM = 1
      IUNKG = 1
      IUNS = -3
      IUNK = -1
      INDX = 6
      LNDX = 9
      IF( IJK.GT.0 ) THEN
        CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
      ELSE
        IDFLT = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(6,IROCK))
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',THKS(6,IROCK)
        INDX = 0
        CALL RDUNIT(UNTS,THKS(6,IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(6,IROCK),', W/m K)'
        IF( KFLD.GT.1 .AND. THKS(6,IROCK)/EPSL.LT.EPSL ) THEN
          INDX = 4
          CHMSG = 'Zero or Negative ' // VARB(1:IVR) //
     &      ' for ' // ROCK(IROCK)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
  390 CONTINUE
      IF( IEQO.EQ.0 ) GOTO 400
      IDFLT = 1
      VARB = 'X-Direction NAPL Saturated Thermal Conductivity'
      UNTS = 'w/m k'
      IUNM = 1
      IUNKG = 1
      IUNS = -3
      IUNK = -1
      INDX = 7
      LNDX = 9
      IF( IJK.GT.0 ) THEN
        CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
      ELSE
        IDFLT = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(7,IROCK))
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',THKS(7,IROCK)
        INDX = 0
        CALL RDUNIT(UNTS,THKS(7,IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(7,IROCK),', W/m K)'
        IF( IFLD.GT.1 .AND. THKS(7,IROCK)/EPSL.LT.EPSL ) THEN
          INDX = 4
          CHMSG = 'Zero or Negative ' // VARB(1:IVR) //
     &      ' for ' // ROCK(IROCK)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IDFLT = 1
      VARB = 'Y-Direction NAPL Saturated Thermal Conductivity'
      UNTS = 'w/m k'
      IUNM = 1
      IUNKG = 1
      IUNS = -3
      IUNK = -1
      INDX = 8
      LNDX = 9
      IF( IJK.GT.0 ) THEN
        CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
      ELSE
        IDFLT = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(8,IROCK))
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',THKS(8,IROCK)
        INDX = 0
        CALL RDUNIT(UNTS,THKS(8,IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(8,IROCK),', W/m K)'
        IF( JFLD.GT.1 .AND. THKS(8,IROCK)/EPSL.LT.EPSL ) THEN
          INDX = 4
          CHMSG = 'Zero or Negative ' // VARB(1:IVR) //
     &      ' for ' // ROCK(IROCK)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IDFLT = 1
      VARB = 'Z-Direction NAPL Saturated Thermal Conductivity'
      UNTS = 'w/m k'
      IUNM = 1
      IUNKG = 1
      IUNS = -3
      IUNK = -1
      INDX = 9
      LNDX = 9
      IF( IJK.GT.0 ) THEN
        CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,THKS,INDX,LNDX )
      ELSE
        IDFLT = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,THKS(9,IROCK))
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',THKS(9,IROCK)
        INDX = 0
        CALL RDUNIT(UNTS,THKS(9,IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',THKS(9,IROCK),', W/m K)'
        IF( KFLD.GT.1 .AND. THKS(9,IROCK)/EPSL.LT.EPSL ) THEN
          INDX = 4
          CHMSG = 'Zero or Negative ' // VARB(1:IVR) //
     &      ' for ' // ROCK(IROCK)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
  400 CONTINUE
!
!---  Rock/soil grain specific heat  ---
!
      IDFLT = 1
      VARB = 'Specific Heat'
      UNTS = 'j/kg k'
      IUNM = 2
      IUNS = -2
      IUNK = -1
      IF( IJK.GT.0 ) THEN
        CALL RDIJK( ISTART,IJK,CHDUM,UNTS,CPS )
      ELSE
        IDFLT = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,CPS(IROCK))
        IDFLT = 1
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',CPS(IROCK)
        INDX = 0
        CALL RDUNIT(UNTS,CPS(IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',CPS(IROCK),', J/kg K)'
        IF( CPS(IROCK)/EPSL.LT.EPSL ) THEN
          INDX = 4
          CHMSG = 'Zero or Negative ' // VARB(1:IVR) //
     &      ' for ' // ROCK(IROCK)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Thermal properties card extensions  ---
!
      VARB = 'Thermal Properties Card Extensions'
      CALL CHKCHR(ISTART,ICOMMA,CHDUM,INDX)
      IF( INDX.EQ.1 ) THEN
        VARB = 'Water-Vapor Enhanced Diffusion'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---    Water-vapor enhanced diffusion  ---
!
        IF( INDEX(ADUM(1:),'enhanc').NE.0 ) THEN
          ISLC(2) = 4
          WRITE(IWR,'(/,A)') VARB(1:IVR)
          DFEF(1,IROCK) = 9.5D+0
          DFEF(2,IROCK) = 2.0D+0
          DFEF(3,IROCK) = 8.0D+0
          DFEF(4,IROCK) = 0.5D+0
          DFEF(5,IROCK) = 3.0D+0
          IDFLT = 1
          VARB = 'Parameter a'
          UNTS = 'null'
          INDX = 1
          LNDX = 5
          IF( IJK.GT.0 ) THEN
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,DFEF,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DFEF(1,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        DFEF(1,IROCK)
          ENDIF
          IDFLT = 1
          VARB = 'Parameter b'
          UNTS = 'null'
          INDX = 2
          LNDX = 5
          IF( IJK.GT.0 ) THEN
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,DFEF,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DFEF(2,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        DFEF(2,IROCK)
          ENDIF
          IDFLT = 1
          VARB = 'Parameter c'
          UNTS = 'null'
          INDX = 3
          LNDX = 5
          IF( IJK.GT.0 ) THEN
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,DFEF,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DFEF(3,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        DFEF(3,IROCK)
          ENDIF
          IDFLT = 1
          VARB = 'Parameter d'
          UNTS = 'null'
          INDX = 4
          LNDX = 5
          IF( IJK.GT.0 ) THEN
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,DFEF,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DFEF(4,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        DFEF(4,IROCK)
          ENDIF
          IDFLT = 1
          VARB = 'Parameter e'
          UNTS = 'null'
          INDX = 5
          LNDX = 5
          IF( IJK.GT.0 ) THEN
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,DFEF,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DFEF(5,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        DFEF(5,IROCK)
          ENDIF
        ENDIF
      ENDIF
!
!---  Loop over remaining rock/soils within scaling group  ---
!
      IF( ISLC(19).EQ.1 .AND. IROCK.LT.NROCK ) THEN
        DO 490 M = IROCK+1,NROCK
          IF( ISCALE(M).EQ.ISGRP ) THEN
            N = N+1
            ITHK(M) = ITHK(IROCK)
            CPS(M) = CPS(IROCK)
            DO 480 L = 1,9
              THKS(L,M) = THKS(L,IROCK)
  480       CONTINUE
            DO 482 L = 1,5
              DFEF(L,M) = DFEF(L,IROCK)
  482       CONTINUE
          ENDIF
  490   CONTINUE
      ENDIF
!
!---  Read next rock/soil type or scaling group  ---
!
      IF( N.LT.NROCK ) WRITE(IWR,'(/)')
      GOTO 10
 500  CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTHER group ---
!
      RETURN
      END

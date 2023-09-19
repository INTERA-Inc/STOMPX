!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDNRP
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
!     Read input file for rock/soil NAPL relative permeability
!     function information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, Battelle, PNL, December 1992.
!     rdnrp.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TABL
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
      REAL*8 RPNCX(LRPNC)
      CHARACTER*4 FORM
      CHARACTER*64 ADUM,RDUM,UNTS
      CHARACTER*512 CHDUM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVAR
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM /'(I9)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDNRP'
      ALLOCATE( IVAR(1:LRC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IVAR'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write card information to ouput file  ---
!
      CARD = 'Rock/Soil NAPL Relative Permeability Function Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Loop over the rock/soil NAPL relative permeability
!     information lines  ---
!
      NR = 0
      IJK = 0
      ISGRP = 0
   10 CONTINUE
      DO L = 1,LRPNC
        RPNCX(L) = 0.D+0
      ENDDO
      IF( NR.GE.NROCK .OR. IJK.GT.0 ) GOTO 500
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Rock/Soil Name'
!
!---  Rock/Soil option for IJK Indexing  ---
!
      IF( IJK.LT.0 ) THEN
        RDUM = 'Rock/Soil #'
        NR = NR + 1
        ICX = ICOUNT(NR)
        WRITE( FORM(3:3),'(I1)' ) ICX
        NCH = 12 + ICX - 1
        WRITE( RDUM(12:NCH),FORM ) NR
        WRITE (IWR,'(A)') RDUM(1:NCH)
        GOTO 220
!
!---  Read rock/soil name  ---
!
      ELSE
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,RDUM)
      ENDIF
!
!---  Check for a pair of delimiting slashes in the rock/soil name,
!     indicating a pattern of rock/soil types  ---
!
      KBS = 0
      IBS = INDEX( RDUM(1:),'/' )
      IF( IBS.GT.0 ) THEN
        IBS = IBS + 1
        JBS = INDEX( RDUM(IBS:),'/')
        IF( JBS.GT.0 ) THEN
          JBS = IBS + JBS - 2
          KBS = 1
          ISBS = ISTART
        ENDIF
      ENDIF
      IROCK = 1
   20 CONTINUE
!
!---  IJK, KIJ, or JKI indexing  ---
!
      IF( INDEX(RDUM(1:),'indexing').NE.0 ) THEN
        IF( INDEX(ROCK(1)(1:),'indexing').EQ.0 ) THEN
          INDX = 4
          CHMSG = 'Indexing Option Not Declared ' // 
     &      'in Rock/Soil Zonation Card'
          CALL WRMSGS( INDX )
        ENDIF
        IF( INDEX(RDUM,'ijk').NE.0 ) THEN
          IJK = 1
        ELSEIF( INDEX(RDUM,'jki').NE.0 ) THEN
          IJK = 2
        ELSEIF( INDEX(RDUM,'kij').NE.0 ) THEN
          IJK = 3
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Indexing Option' // RDUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        IROCK = 1
        GOTO 220
      ENDIF
!
!---  Search known rock types for a matching type ---
!
      DO 100 M = IROCK,NROCK
        IF( KBS.EQ.1 ) THEN
          IF( INDEX( ROCK(M)(1:),RDUM(IBS:JBS) ).GT.0 ) THEN
            IROCK = M
            GOTO 200
          ENDIF
        ELSE
          IF( RDUM.EQ.ROCK(M) ) THEN
            IROCK = M
            GOTO 200
          ENDIF
        ENDIF
  100 CONTINUE
!
!---  Search known scaling groups for a matching type ---
!
      IF( ISLC(19).EQ.1 ) THEN
        DO 110 M = 1,NSCALE
           IF( RDUM.EQ.SCALNM(M) ) THEN
              ISGRP = M
              IROCK = 1
              GOTO 200
           ENDIF
  110   CONTINUE
        INDX = 2
        CHMSG = 'Unrecognized Rock/Soil Type or Scaling Group: '
     &    // RDUM(1:NCH)
        CALL WRMSGS( INDX )
        GOTO 10
      ENDIF
      INDX = 2
      CHMSG = 'Unrecognized Rock/Soil Type: ' // RDUM(1:NCH)
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
!
!---  Read gas relative permeability pressure function  ---
!
      NR = NR + 1
  220 CONTINUE
!
!---  Rock/Soil option for IJK Indexing, dissociate from
!     saturation function type  ---
!
      IF( IJK.LT.0 ) THEN
        ISCHRX = 0
      ELSE
        ISCHRX = MOD( ISCHR(IROCK),1000 )
      ENDIF
!
!---  Nonhysteretic saturation functions  ---
!
      IF( ISCHRX.LT.20 .OR. ISCHRX.GT.30 ) THEN
        VARB = 'NAPL Relative Permeability Function'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---    Tabular (relative permeability versus NAPL saturation)  ---
!
        IF( INDEX(ADUM(1:),'tabular').NE.0 ) THEN
          IF( INDEX( ADUM(1:),'head' ).NE.0 ) THEN
            IF( INDEX( ADUM(1:),'log' ).NE.0 ) THEN
              WRITE(IWR,'(A)') 'Tabular NAPL Relative Permeability '
     &          // 'Versus Log Capillary Head Function'
              IRPGX = 14
            ELSE
              WRITE(IWR,'(A)') 'Tabular NAPL Relative Permeability '
     &          // 'Versus Capillary Head Function'
                IRPGX = 12
            ENDIF
          ELSE
            WRITE(IWR,'(A)') 'Tabular NAPL Relative Permeability '
     &        // 'Versus NAPL Saturation Function'
            IRPGX = 10
          ENDIF
          IF( INDEX( ADUM(1:),'spline' ).NE.0 ) THEN
            IRPGX = IRPGX + 1
            WRITE(IWR,'(A)') 'Cubic Spline Interpolation'
          ELSE
            WRITE(IWR,'(A)') 'Linear Interpolation'
          ENDIF
!
!---      IJK Indexing  ---
!
          IF( IJK.GT.0 ) THEN
            VARB = 'Number of Tables'
            CALL RDINT(ISTART,ICOMMA,CHDUM,NTABLX)
            IF( NTABLX.LT.1 ) THEN
              INDX = 4
              CHMSG = 'Invalid Number of NAPL Relative ' // 
     &          'Permeability Tables'
              CALL WRMSGS( INDX )
            ENDIF
            CALL RDIJKI( ISTART,IJK,CHDUM,IVAR )
!
!---        Loop over NAPL relative permeability function tables  ---
!
            DO 280 NTX = 1,NTABLX
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              ISTART = 1
              VARB = 'Number of Table Entries'
              CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
              WRITE(IWR,'(2A,I6)') VARB(1:IVR),': ',NLIN
!
!---          Loop over lines in NAPL relative permeability 
!             function tables  ---
!
              NTBLX = NTBL+1
              DO 270 NL = 1,NLIN
                NTBL = NTBL + 1
                IF( NTBL.GT.LTBL ) THEN
                  INDX = 5
                  CHMSG = 'Number of Table Values > Parameter LTBL'
                  CALL WRMSGS( INDX )
                ENDIF
                ISTART = 1
                CALL RDINPL( CHDUM )
                CALL LCASE( CHDUM )
                IF( INDEX( ADUM(1:),'head' ).NE.0 ) THEN
                  IF( INDEX( ADUM(1:),'log' ).NE.0 ) THEN
                    VARB = 'Log Capillary Head'
                  ELSE
                    VARB = 'Capillary Head'
                  ENDIF
                ELSE
                  VARB = 'NAPL Saturation'
                ENDIF
!
!---            Correct table values for capillary-head units  ---
!
                IF( IRPGX.GE.12 .AND. IRPGX.LE.15 ) THEN
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLX(NTBL))
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &              UNTS(1:NCH),': ',TBLX(NTBL)
                  INDX = 0
                  IUNM = 1
                  VARX = 1.D+0
                  CALL RDUNIT(UNTS,VARX,INDX)
                  IF( IRPGX.GE.14 .AND. IRPGX.LE.15 ) THEN
                    TBLX(NTBL) = LOG( TBLX(NTBL)*VARX )
                  ELSE
                    TBLX(NTBL) = TBLX(NTBL)*VARX
                  ENDIF
                  WRITE(IWR,'(A,1PE11.4,A)') ' (',TBLX(NTBL),', m)'
                ELSE
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLX(NTBL))
                  WRITE(IWR,'(4X,A,1PE11.4)') VARB,TBLX(NTBL)
                ENDIF
                VARB = 'NAPL Relative Permeability'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
                WRITE(IWR,'(4X,A,1PE11.4)') VARB,TBLY(NTBL)
                IF( NL.EQ.2 ) THEN
                  IF( TBLX(NTBL-1).LT.TBLX(NTBL) ) THEN
                    ITDX = 1
                  ELSEIF( TBLX(NTBL-1).GT.TBLX(NTBL) ) THEN
                    ITDX = -1
                  ELSE
                    INDX = 4
                    CHMSG = 'Invalid NAPL Relative ' // 
     &                'Permeability Table'
                    CALL WRMSGS( INDX )
                  ENDIF
                ELSEIF( NL.GT.2 ) THEN
                  IF( (ITDX.EQ.1 .AND. TBLX(NTBL).LE.TBLX(NTBL-1)) .OR.
     &              (ITDX.EQ.-1 .AND. TBLX(NTBL).GE.TBLX(NTBL-1)) ) THEN
                    INDX = 4
                    CHMSG = 'Invalid NAPL Relative ' // 
     &                'Permeability Table'
                    CALL WRMSGS( INDX )
                  ENDIF
                ENDIF
  270         CONTINUE
!
!---          Build cubic splines  ---
!
              IF( IRPGX.EQ.11 .OR. IRPGX.EQ.13 .OR. IRPGX.EQ.15 ) THEN
                CALL SPLINY( NTBLX,NTBL )
              ENDIF
!
!---          Correlate table numbers with nodes  ---
!
              DO 272 N = 1,NFLD
                IF( IVAR(N).EQ.NTX ) THEN
                  IRNTBL(1,N) = NTBLX
                  IRNTBL(2,N) = NTBL
                ELSEIF( IVAR(N).LT.1 .OR. IVAR(N).GT.NTABLX ) THEN
                  INDX = 4
                  CHMSG = 'Invalid NAPL Relative Permeability ' //
     &              'Table Number'
                  CALL WRMSGS( INDX )
                ENDIF
 272          CONTINUE
 280        CONTINUE
!
!---      Rock/soil zonation  ---
!
          ELSE
            VARB = 'Number of Tabular Entries'
            CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
            WRITE(IWR,'(2X,2A,I6)') VARB,': ',NLIN
            IF( NLIN.LT.2 ) THEN
              INDX = 4
              CHMSG = 'Invalid NAPL Relative Permeability Table'
              CALL WRMSGS( INDX )
            ENDIF
            IRNTBL(1,IROCK) = NTBL + 1
            DO 300 NL = 1,NLIN
              NTBL = NTBL + 1
              IF( NTBL.GT.LTBL ) THEN
                INDX = 5
                CHMSG = 'Number of Tables Values > Parameter LTBL'
                CALL WRMSGS( INDX )
              ENDIF
              ISTART = 1
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              IF( INDEX( ADUM(1:),'head' ).NE.0 ) THEN
                IF( INDEX( ADUM(1:),'log' ).NE.0 ) THEN
                  VARB = 'Log Capillary Head'
                ELSE
                  VARB = 'Capillary Head'
                ENDIF
              ELSE
                VARB = 'NAPL Saturation'
              ENDIF
!
!---          Correct table values for capillary-head units  ---
!
              IF( IRPGX.GE.12 .AND. IRPGX.LE.15 ) THEN
                CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLX(NTBL))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',TBLX(NTBL)
                INDX = 0
                IUNM = 1
                VARX = 1.D+0
                CALL RDUNIT(UNTS,VARX,INDX)
                IF( IRPGX.GE.14 .AND. IRPGX.LE.15 ) THEN
                  TBLX(NTBL) = LOG( TBLX(NTBL)*VARX )
                ELSE
                  TBLX(NTBL) = TBLX(NTBL)*VARX
                ENDIF
                WRITE(IWR,'(A,1PE11.4,A)') ' (',TBLX(NTBL),', m)'
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLX(NTBL))
                WRITE(IWR,'(4X,A,1PE11.4)') VARB,TBLX(NTBL)
              ENDIF
              VARB = 'NAPL Relative Permeability'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
              WRITE(IWR,'(4X,A,1PE11.4)') VARB,TBLY(NTBL)
              IF( NL.EQ.2 ) THEN
                IF( TBLX(NTBL-1).LT.TBLX(NTBL) ) THEN
                  ITDX = 1
                ELSEIF( TBLX(NTBL-1).GT.TBLX(NTBL) ) THEN
                  ITDX = -1
                ELSE
                  INDX = 4
                  CHMSG = 'Invalid NAPL Relative ' // 
     &              'Permeability Table'
                  CALL WRMSGS( INDX )
                ENDIF
              ELSEIF( NL.GT.2 ) THEN
                IF( (ITDX.EQ.1 .AND. TBLX(NTBL).LE.TBLX(NTBL-1)) .OR.
     &            (ITDX.EQ.-1 .AND. TBLX(NTBL).GE.TBLX(NTBL-1)) ) THEN
                  INDX = 4
                  CHMSG = 'Invalid NAPL Relative ' // 
     &              'Permeability Table'
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
  300       CONTINUE
            IRNTBL(2,IROCK) = NTBL
            IF( IRPGX.EQ.11 .OR. IRPGX.EQ.13 .OR. IRPGX.EQ.15 ) THEN
              CALL SPLINY( IRNTBL(1,IROCK),IRNTBL(2,IROCK) )
            ENDIF
          ENDIF
          GOTO 460
        ENDIF
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHRX.EQ.1 .OR. ISCHRX.EQ.11 .OR.
     &    ISCHRX.EQ.13 .OR.
     &    ISCHRX.EQ.15 .OR. ISCHRX.EQ.17 .OR.
     &    (ISCHRX.GE.31 .AND. ISCHRX.LE.34) ) THEN
          IF( INDEX(ADUM(1:),'constant') .NE. 0 ) THEN
            IRPNX = 0
            VARB = 'NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem') .NE. 0 ) THEN
            IRPNX = 1
            IDFLT = 1
            VARB = 'Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 305 N = 1,NFLD
                RPNC(1,IZ(N)) = 1.D+0-1.D+0/SCHR(3,IZ(N))
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  305         CONTINUE
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
              DO 310 N = 1,NFLD
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  310         CONTINUE
            ELSE
              RPNCX(1) = 1.D+0-1.D+0/SCHR(3,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 1/n'
              WRITE(IWR,'(4X,A,1PE11.4)') 'm Parameter: ',
     &          RPNCX(1)
              IF( RPNCX(1).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine') .NE. 0 ) THEN
            IRPNX = 2
            IDFLT = 1
            VARB = 'Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 315 N = 1,NFLD
                RPNC(1,IZ(N)) = 1.D+0-2.D+0/SCHR(3,IZ(N))
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  315         CONTINUE
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
              DO 320 N = 1,NFLD
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  320         CONTINUE
            ELSE
              RPNCX(1) = 1.D+0-2.D+0/SCHR(3,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 2/n'
              WRITE(IWR,'(4X,A,1PE11.4)') 'm Parameter: ',
     &          RPNCX(1)
              IF( RPNCX(1).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'fatt and klikoff') .NE. 0 ) THEN
            IRPNX = 4
            WRITE (IWR,'(2X,A)') 'Fatt and Klikoff Relative ' //
     &        'Permeability Function'
          ELSEIF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            IRPNX = 9
            WRITE (IWR,'(2X,A)')'van Genuchten NAPL Relative ' // 
     &        'Permeability Model'
            IF( IJK.GT.0 ) THEN
              DO N = 1,NFLD
                RPNC(1,IZ(N)) = SCHR(4,IZ(N))
              ENDDO
            ELSE
              RPNCX(1) = SCHR(4,IROCK)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Irreducible NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'stone') .NE. 0 ) THEN
            IRPNX = 5
            WRITE (IWR,'(2X,A)')'Stone NAPL Relative ' //
     &        'Permeability Model'
            VARB = 'Stone (Slr)'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
            VARB = 'Stone (Snr)'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Stone (n)'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'variable corey').NE.0 ) THEN
            IRPNX = 20
            WRITE (IWR,'(2X,A)')'Variable Corey Relative ' // 
     &        'Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(1)
            ENDIF
            VARB = 'Irreducible NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(3)
            ENDIF
            VARB = 'Endpoint NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(2)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'extended').NE.0 .AND.
     &      INDEX(ADUM(1:),'power').NE.0 .AND.
     &      INDEX(ADUM(1:),'law').NE.0 ) THEN
            IRPNX = 21
            WRITE (IWR,'(2X,A)')'Extended Power Law Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(1)
            ENDIF
            VARB = 'Exponent NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(3)
            ENDIF
            VARB = 'Residual NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'free corey') .NE. 0 ) THEN
            IRPNX = 7
            WRITE (IWR,'(2X,A)')'Free Corey Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
            VARB = 'Exponent NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
            VARB = 'Residual NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'classical corey').NE.0 ) THEN
            IRPNX = 17
            WRITE (IWR,'(2X,A)')'Classical Corey Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
            VARB = 'Exponent NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
            VARB = 'Residual NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'corey') .NE. 0 ) THEN
            IRPNX = 3
            WRITE (IWR,'(2X,A)')'Corey Relative Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO N = 1,NFLD
                RPNC(1,IZ(N)) = SCHR(4,IZ(N))
              ENDDO
            ELSE
              RPNCX(1) = SCHR(4,IROCK)
            ENDIF
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
            VARB = 'Irreducible NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO N = 1,NFLD
                RPNC(2,IZ(N)) = 1.D+0
              ENDDO
            ELSE
              RPNCX(2) = 1.D+0
            ENDIF
            CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
            IF( INDX.EQ.1 ) THEN
              IF( IJK.GT.0 ) THEN
                INDX = 2
                LNDX = LRPNC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPNCX(2)
              ENDIF
            ENDIF
          ELSE
            INDX = 4
            NCH = INDEX( ADUM(1:),'  ' )-1
            CHMSG = 'Unrecognized NAPL Relative Perm. Function: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHRX.EQ.2 .OR. ISCHRX.EQ.12 .OR.
     &    ISCHRX.EQ.14 .OR.
     &    ISCHRX.EQ.16 .OR. ISCHRX.EQ.18 .OR.
     &    (ISCHRX.GE.35 .AND. ISCHRX.LE.38) ) THEN
          IF( INDEX(ADUM(1:),'constant') .NE. 0 ) THEN
            IRPNX = 0
            VARB = 'NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem') .NE. 0 ) THEN
            IRPNX = 1
            IDFLT = 1
            VARB = 'Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 325 N = 1,NFLD
                RPNC(1,IZ(N)) = SCHR(3,IZ(N))
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative ''lambda'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  325         CONTINUE
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              RPNCX(1) = SCHR(3,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPNCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine') .NE. 0 ) THEN
            IRPNX = 2
            IDFLT = 1
            VARB = 'Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 335 N = 1,NFLD
                RPNC(1,IZ(N)) = SCHR(3,IZ(N))
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative ''lambda'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  335         CONTINUE
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              RPNCX(1) = SCHR(3,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPNCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'fatt and klikoff') .NE. 0 ) THEN
            IRPNX = 4
            WRITE (IWR,'(2X,A)') 'Fatt and Klikoff Relative' //
     &        ' Permeability Function'
          ELSEIF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            IRPNX = 9
            WRITE (IWR,'(2X,A)')'van Genuchten NAPL Relative ' // 
     &        'Permeability Model'
            IF( IJK.GT.0 ) THEN
              DO N = 1,NFLD
                RPNC(1,IZ(N)) = SCHR(4,IZ(N))
              ENDDO
            ELSE
              RPNCX(1) = SCHR(4,IROCK)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Irreducible NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'stone') .NE. 0 ) THEN
            IRPNX = 5
            WRITE (IWR,'(2X,A)')'Stone NAPL Relative ' //
     &        'Permeability Model'
            VARB = 'Stone (Slr)'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
            VARB = 'Stone (Snr)'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Stone (n)'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'classical corey').NE.0 ) THEN
            IRPNX = 17
            WRITE (IWR,'(2X,A)')'Classical Corey Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
            VARB = 'Exponent NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
            VARB = 'Residual NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'variable corey').NE.0 ) THEN
            IRPNX = 20
            WRITE (IWR,'(2X,A)')'Variable Corey Relative ' // 
     &        'Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(1)
            ENDIF
            VARB = 'Irreducible NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(3)
            ENDIF
            VARB = 'Endpoint NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(2)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'extended').NE.0 .AND.
     &      INDEX(ADUM(1:),'power').NE.0 .AND.
     &      INDEX(ADUM(1:),'law').NE.0 ) THEN
            IRPNX = 21
            WRITE (IWR,'(2X,A)')'Extended Power Law Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(1)
            ENDIF
            VARB = 'Exponent NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(3)
            ENDIF
            VARB = 'Residual NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'free corey') .NE. 0 ) THEN
            IRPNX = 7
            WRITE (IWR,'(2X,A)')'Free Corey Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
            VARB = 'Exponent NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
            VARB = 'Residual NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'corey') .NE. 0 ) THEN
            IRPNX = 3
            WRITE (IWR,'(2X,A)')'Corey Relative Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO N = 1,NFLD
                RPNC(1,IZ(N)) = SCHR(4,IZ(N))
              ENDDO
            ELSE
              RPNCX(1) = SCHR(4,IROCK)
            ENDIF
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
            VARB = 'Irreducible NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO N = 1,NFLD
                RPNC(2,IZ(N)) = 1.D+0
              ENDDO
            ELSE
              RPNCX(2) = 1.D+0
            ENDIF
            CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
            IF( INDX.EQ.1 ) THEN
              IF( IJK.GT.0 ) THEN
                INDX = 2
                LNDX = LRPNC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPNCX(2)
              ENDIF
            ENDIF
          ELSE
            INDX = 4
            NCH = INDEX( ADUM(1:),'  ' )-1
            CHMSG = 'Unrecognized NAPL Relative Perm. Function: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Dual Porosity van Genuchten function  ---
!
        ELSEIF( ISCHRX.EQ.3 ) THEN
          IF( INDEX(ADUM(1:),'constant').NE.0 ) THEN
            IRPNX = 0
            VARB = 'Matrix NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Fracture NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem').NE.0 ) THEN
            IRPNX = 1
            IDFLT = 1
            VARB = 'Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 340 N = 1,NFLD
                RPNC(2,IZ(N)) = 1.D+0-1.D+0/SCHR(3,IZ(N))
                IF( RPNC(2,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Matrix van Genuchten' // 
     &              ' ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  340         CONTINUE
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
              DO 342 N = 1,NFLD
                IF( RPNC(2,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Matrix van Genuchten' // 
     &              ' ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  342         CONTINUE
            ELSE
              RPNCX(2) = 1.D+0-1.D+0/SCHR(3,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 1/n'
              WRITE(IWR,'(4X,A,1PE11.4)') 'Matrix m Parameter: ',
     &          RPNCX(2)
              IF( RPNCX(2).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative Matrix van Genuchten' // 
     &            ' ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
            IDFLT = 1
            VARB = 'Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 344 N = 1,NFLD
                RPNC(1,IZ(N)) = 1.D+0-1.D+0/SCHR(6,IZ(N))
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Fracture van Genuchten' // 
     &              ' ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  344         CONTINUE
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
              DO 346 N = 1,NFLD
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Fracture van Genuchten' // 
     &              ' ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  346         CONTINUE
            ELSE
              RPNCX(1) = 1.D+0-1.D+0/SCHR(6,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 1/n'
              WRITE(IWR,'(4X,A,1PE11.4)') 'Fracture m Parameter: ',
     &          RPNCX(1)
              IF( RPNCX(1).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative Fracture van Genuchten' // 
     &            ' ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 ) THEN
            IRPNX = 2
            IDFLT = 1
            VARB = 'Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 350 N = 1,NFLD
                RPNC(2,IZ(N)) = 1.D+0-2.D+0/SCHR(3,IZ(N))
                IF( RPNC(2,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Matrix van Genuchten' // 
     &              ' ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  350         CONTINUE
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
              DO 352 N = 1,NFLD
                IF( RPNC(2,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Matrix van Genuchten' // 
     &              ' ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  352         CONTINUE
            ELSE
              RPNCX(2) = 1.D+0-2.D+0/SCHR(3,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 2/n'
              WRITE(IWR,'(4X,A,1PE11.4)') 'Matrix m Parameter: ',
     &          RPNCX(2)
              IF( RPNCX(2).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative Matrix van Genuchten' // 
     &            ' ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
            IDFLT = 1
            VARB = 'Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 354 N = 1,NFLD
                RPNC(1,IZ(N)) = 1.D+0-2.D+0/SCHR(6,IZ(N))
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Fracture van Genuchten' // 
     &              ' ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  354         CONTINUE
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
              DO 356 N = 1,NFLD
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Fracture van Genuchten' // 
     &              ' ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  356         CONTINUE
            ELSE
              RPNCX(1) = 1.D+0-2.D+0/SCHR(6,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 2/n'
              WRITE(IWR,'(4X,A,1PE11.4)') 'Fracture m Parameter: ',
     &          RPNCX(1)
              IF( RPNCX(1).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative Fracture van Genuchten' // 
     &            ' ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSE
            INDX = 4
            NCH = INDEX( ADUM(1:),'  ' )-1
            CHMSG = 'Unrecognized Relative Perm. Function: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Dual Porosity Brooks and Corey function  ---
!
        ELSEIF( ISCHRX.EQ.4 ) THEN
          IF( INDEX(ADUM(1:),'constant').NE.0 ) THEN
            VARB = 'Matrix NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Fracture NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem').NE.0 ) THEN
            IRPNX = 1
            IDFLT = 1
            VARB = 'Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 360 N = 1,NFLD
                RPNC(2,IZ(N)) = SCHR(3,IZ(N))
                IF( RPNC(2,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Matrix ''lambda'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  360         CONTINUE
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              RPNCX(2) = SCHR(3,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              IF( RPNCX(2).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative Matrix ''lambda'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Matrix ''lambda'' ' // 
     &          'Parameter: ',RPNCX(2)
            ENDIF
            IDFLT = 1
            VARB = 'Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 362 N = 1,NFLD
                RPNC(1,IZ(N)) = SCHR(6,IZ(N))
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Fracture ''lambda'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  362         CONTINUE
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              RPNCX(1) = SCHR(6,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              IF( RPNCX(1).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative Fracture ''lambda'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Fracture ''lambda'' ' // 
     &          'Parameter: ',RPNCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 ) THEN
            IRPNX = 2
            IDFLT = 1
            VARB = 'Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 370 N = 1,NFLD
                RPNC(2,IZ(N)) = SCHR(3,IZ(N))
                IF( RPNC(2,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Matrix ''lambda'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  370         CONTINUE
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              RPNCX(2) = SCHR(3,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              IF( RPNC(2,IZ(N)).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative Matrix ''lambda'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Matrix ''lambda'' ' // 
     &          'Parameter: ',
     &          RPNCX(2)
            ENDIF
            IDFLT = 1
            VARB = 'Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              DO 372 N = 1,NFLD
                RPNC(1,IZ(N)) = SCHR(6,IZ(N))
                IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative Fracture ''lambda'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  372         CONTINUE
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              RPNCX(1) = SCHR(6,IROCK)
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              IF( RPNC(1,IZ(N)).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative Fracture ''lambda'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Fracture ''lambda'' ' // 
     &          'Parameter: ',
     &          RPNCX(1)
            ENDIF
          ELSE
            INDX = 4
            NCH = INDEX( ADUM(1:),'  ' )-1
            CHMSG = 'Unrecognized Relative Perm. Function: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'classical corey').NE.0 ) THEN
          IRPNX = 17
          WRITE (IWR,'(2X,A)')'Classical Corey Relative ' //
     &      'Permeability Model'
          VARB = 'Endpoint NAPL Relative Permeability'
          IF( IJK.GT.0 ) THEN
            INDX = 1
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        RPNCX(1)
          ENDIF
          VARB = 'Exponent NAPL Relative Permeability'
          IF( IJK.GT.0 ) THEN
            INDX = 2
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        RPNCX(2)
          ENDIF
          VARB = 'Residual Aqueous Saturation'
          IF( IJK.GT.0 ) THEN
            INDX = 3
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        RPNCX(3)
          ENDIF
          VARB = 'Residual NAPL Saturation'
          IF( IJK.GT.0 ) THEN
            INDX = 4
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        RPNCX(4)
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'variable corey').NE.0 ) THEN
          IRPNX = 20
          WRITE (IWR,'(2X,A)')'Variable Corey Relative ' // 
     &      'Permeability Model'
          VARB = 'Irreducible Aqueous Saturation'
          IF( IJK.GT.0 ) THEN
            INDX = 1
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(1)
          ENDIF
          VARB = 'Irreducible NAPL Saturation'
          IF( IJK.GT.0 ) THEN
            INDX = 3
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(3)
          ENDIF
          VARB = 'Endpoint NAPL Relative Permeability'
          IF( IJK.GT.0 ) THEN
            INDX = 2
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(2)
          ENDIF
          VARB = 'Exponent'
          IF( IJK.GT.0 ) THEN
            INDX = 4
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(4)
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'extended').NE.0 .AND.
     &    INDEX(ADUM(1:),'power').NE.0 .AND.
     &    INDEX(ADUM(1:),'law').NE.0 ) THEN
          IRPNX = 21
          WRITE (IWR,'(2X,A)')'Extended Power Law Relative ' //
     &      'Permeability Model'
          VARB = 'Endpoint NAPL Relative Permeability'
          IF( IJK.GT.0 ) THEN
            INDX = 1
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(1)
          ENDIF
          VARB = 'Exponent NAPL Relative Permeability'
          IF( IJK.GT.0 ) THEN
            INDX = 2
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(2)
          ENDIF
          VARB = 'Residual Aqueous Saturation'
          IF( IJK.GT.0 ) THEN
            INDX = 3
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(3)
          ENDIF
          VARB = 'Residual NAPL Saturation'
          IF( IJK.GT.0 ) THEN
            INDX = 4
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPNCX(4)
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'free corey') .NE. 0 ) THEN
          IRPGX = 7
          WRITE (IWR,'(2X,A)')'Free Corey Relative ' //
     &      'Permeability Model'
          VARB = 'Endpoint NAPL Relative Permeability'
          IF( IJK.GT.0 ) THEN
            INDX = 1
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        RPNCX(1)
          ENDIF
          VARB = 'Exponent NAPL Relative Permeability'
          IF( IJK.GT.0 ) THEN
            INDX = 2
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        RPNCX(2)
          ENDIF
          VARB = 'Residual Aqueous Saturation'
          IF( IJK.GT.0 ) THEN
            INDX = 3
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        RPNCX(3)
          ENDIF
          VARB = 'Residual Gas Saturation'
          IF( IJK.GT.0 ) THEN
            INDX = 4
            LNDX = LRPNC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(4))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        RPNCX(4)
          ENDIF
        ELSE
          IF( INDEX(ADUM(1:),'constant') .NE. 0 ) THEN
            IRPNX = 0
            VARB = 'NAPL Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'corey') .NE. 0 ) THEN
            IRPNX = 3
            WRITE (IWR,'(2X,A)')'Corey Relative Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO N = 1,NFLD
                RPNC(1,IZ(N)) = SCHR(4,IZ(N))
              ENDDO
            ELSE
              RPNCX(1) = SCHR(4,IROCK)
            ENDIF
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(1)
            ENDIF
            VARB = 'Irreducible NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO N = 1,NFLD
                RPNC(2,IZ(N)) = 1.D+0
              ENDDO
            ELSE
              RPNCX(2) = 1.D+0
            ENDIF
            CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
            IF( INDX.EQ.1 ) THEN
              IF( IJK.GT.0 ) THEN
                INDX = 2
                LNDX = LRPNC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPNCX(2)
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'fatt and klikoff') .NE. 0 ) THEN
            IRPNX = 4
            WRITE (IWR,'(2X,A)') 'Fatt and Klikoff Relative' //
     &        ' Permeability Function'
          ELSEIF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            IRPNX = 9
            WRITE (IWR,'(2X,A)')'van Genuchten NAPL Relative ' // 
     &        'Permeability Model'
            IF( IJK.GT.0 ) THEN
              DO N = 1,NFLD
                RPNC(1,IZ(N)) = SCHR(4,IZ(N))
              ENDDO
            ELSE
              RPNCX(1) = SCHR(4,IROCK)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(2)
            ENDIF
            VARB = 'Irreducible NAPL Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPNC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPNC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPNCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPNCX(3)
            ENDIF
          ELSE
            INDX = 4
            NCH = INDEX( ADUM(1:),'  ' )-1
            CHMSG = 'Unrecognized NAPL Relative Perm. Function: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Hysteretic saturation functions  ---
!
      ELSE
        VARB = 'Porosity Distribution Model'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX(ADUM(1:), 'mualem') .NE. 0 ) THEN
          IRPNX = 1
          WRITE(IWR,'(2X,A)') 'Mualem Porosity Distribution Model'
        ELSEIF( INDEX(ADUM(1:),'burdine') .NE. 0 ) THEN
          IRPNX = 2
          WRITE(IWR,'(2X,A)') 'Burdine Porosity Distribution Model'
        ELSE
          INDX = 4
          NCH = INDEX( ADUM(1:),'  ' )-1
          CHMSG = 'Unrecognized NAPL Relative Perm. Function: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Translate NAPL relative permeability type for IJK indexing  ---
!
  460 CONTINUE
      IF( IJK.GT.0 ) THEN
        DO 462 N = 1,NFLD
          IRPN(IZ(N)) = IRPNX
  462   CONTINUE
!
!---  For IJK indexing with the rock/soil option, correlate rock/soil 
!     numbers with nodes for NAPL relative permeability type and
!     parameters  ---
!
      ELSEIF( IJK.LT.0 ) THEN
        DO 472 N = 1,NFLD
          IF( IVAR(N).EQ.NR ) THEN
            DO 470 L = 1,LRPNC
              RPNC(L,N) = RPNCX(L)
  470       CONTINUE
            IRPN(N) = IRPNX
          ENDIF
  472   CONTINUE
!
!---  For rock/soil zonation input translate NAPL relative permeability 
!     type and parameters  ---
!
      ELSE
        IRPN(IROCK) = IRPNX
        DO 480 L = 1,LRPNC
          RPNC(L,IROCK) = RPNCX(L)
  480   CONTINUE
      ENDIF
!
!---  Loop over remaining rock/soils within scaling group  ---
!
      IF( ISLC(19).EQ.1 .AND. IROCK.LT.NROCK ) THEN
        DO 494 M = IROCK+1,NROCK
          IF( ISCALE(M).EQ.ISGRP ) THEN
            NR = NR + 1
            IRPN(M) = IRPN(IROCK)
            DO 490 L = 1,LRPNC
              RPNC(L,M) = RPNC(L,IROCK)
  490       CONTINUE
            DO 492 L = 1,2
              IRNTBL(L,M) = IRNTBL(L,IROCK)
  492       CONTINUE
          ENDIF
  494   CONTINUE
      ENDIF
!
!---  Read next rock/soil type or scaling group  ---
!
      IF( NR.LT.NROCK ) WRITE(IWR,'(/)')
!
!---  Continue reading rock/soil type names for a pattern match  ---
!
      IF( KBS.EQ.1 .AND. IROCK.LT.NROCK ) THEN
        IROCK = IROCK + 1
        ISTART = ISBS
        GOTO 20
      ENDIF
      GOTO 10
 500  CONTINUE
      IF( ALLOCATED(IVAR) ) THEN
      DEALLOCATE( IVAR,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IVAR'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDNRP group ---
!
      RETURN
      END

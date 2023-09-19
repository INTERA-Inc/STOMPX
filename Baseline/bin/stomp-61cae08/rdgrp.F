!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDGRP
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
!     Read input file for rock/soil gas relative permeability
!     function information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, Battelle, PNL, December 1992.
!     rdgrp.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 RPGCX(LRPGC)
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
      SUB_LOG(ISUB_LOG) = '/RDGRP'
      ALLOCATE( IVAR(1:LRC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IVAR'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write card information to ouput file  ---
!
      CARD = 'Rock/Soil Gas Relative Permeability Function Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Loop over the rock/soil gas relative permeability
!     information lines  ---
!
      NR = 0
      IJK = 0
      ISGRP = 0
   10 CONTINUE
      IF( NR.GE.NROCK .OR. IJK.GT.0 ) GOTO 500
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Rock/Soil Name'
   12 CONTINUE
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
        VARB = 'Gas Relative Permeability Function'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---    Rock/Soil Zonation Option for IJK Indexing  ---
!
        IF( IJK.GT.0 .AND. INDEX(ADUM(1:),'rock/soil').NE.0 ) THEN
          VARB = 'Number of Rock/Soil Entries'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NROCK)
          CALL RDIJKI( ISTART,IJK,CHDUM,IVAR )
          IJK = -IJK
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          GOTO 12
        ENDIF
!
!---    Tabular (relative permeability versus liquid saturation)  ---
!
        IF( INDEX(ADUM(1:),'tabular').NE.0 ) THEN
          IF( INDEX( ADUM(1:),'head' ).NE.0 ) THEN
            IF( INDEX( ADUM(1:),'log' ).NE.0 ) THEN
              WRITE(IWR,'(A)') 'Tabular Gas Relative Permeability '
     &          // 'Versus Log Capillary Head Function'
              IRPGX = 14
            ELSE
              WRITE(IWR,'(A)') 'Tabular Gas Relative Permeability '
     &          // 'Versus Capillary Head Function'
              IRPGX = 12
            ENDIF
          ELSE
            WRITE(IWR,'(A)') 'Tabular Gas Relative Permeability '
     &        // 'Versus Gas Saturation Function'
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
              CHMSG = 'Invalid Number of Gas Relative ' // 
     &          'Permeability Tables'
              CALL WRMSGS( INDX )
            ENDIF
            CALL RDIJKI( ISTART,IJK,CHDUM,IVAR )
!
!---        Loop over gas relative permeability function tables  ---
!
            DO 280 NTX = 1,NTABLX
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              ISTART = 1
              VARB = 'Number of Table Entries'
              CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
              WRITE(IWR,'(2A,I6)') VARB(1:IVR),': ',NLIN
!
!---          Loop over lines in gas relative permeability 
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
                  VARB = 'Gas Saturation'
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
                VARB = 'Gas Relative Permeability'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
                WRITE(IWR,'(4X,A,1PE11.4)') VARB,TBLY(NTBL)
                IF( NL.EQ.2 ) THEN
                  IF( TBLX(NTBL-1).LT.TBLX(NTBL) ) THEN
                    ITDX = 1
                  ELSEIF( TBLX(NTBL-1).GT.TBLX(NTBL) ) THEN
                    ITDX = -1
                  ELSE
                    INDX = 4
                    CHMSG = 'Invalid Gas Relative ' // 
     &                'Permeability Table'
                    CALL WRMSGS( INDX )
                  ENDIF
                ELSEIF( NL.GT.2 ) THEN
                  IF( (ITDX.EQ.1 .AND. TBLX(NTBL).LE.TBLX(NTBL-1)) .OR.
     &              (ITDX.EQ.-1 .AND. TBLX(NTBL).GE.TBLX(NTBL-1)) ) THEN
                    INDX = 4
                    CHMSG = 'Invalid Gas Relative ' // 
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
                  IRGTBL(1,N) = NTBLX
                  IRGTBL(2,N) = NTBL
                ELSEIF( IVAR(N).LT.1 .OR. IVAR(N).GT.NTABLX ) THEN
                  INDX = 4
                  CHMSG = 'Invalid Gas Relative Permeability ' //
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
              CHMSG = 'Invalid Gas Relative Permeability Table'
              CALL WRMSGS( INDX )
            ENDIF
            IRGTBL(1,IROCK) = NTBL + 1
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
                VARB = 'Gas Saturation'
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
              VARB = 'Gas Relative Permeability'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
              WRITE(IWR,'(4X,A,1PE11.4)') VARB,TBLY(NTBL)
              IF( NL.EQ.2 ) THEN
                IF( TBLX(NTBL-1).LT.TBLX(NTBL) ) THEN
                  ITDX = 1
                ELSEIF( TBLX(NTBL-1).GT.TBLX(NTBL) ) THEN
                  ITDX = -1
                ELSE
                  INDX = 4
                  CHMSG = 'Invalid Gas Relative ' // 
     &              'Permeability Table'
                  CALL WRMSGS( INDX )
                ENDIF
              ELSEIF( NL.GT.2 ) THEN
                IF( (ITDX.EQ.1 .AND. TBLX(NTBL).LE.TBLX(NTBL-1)) .OR.
     &            (ITDX.EQ.-1 .AND. TBLX(NTBL).GE.TBLX(NTBL-1)) ) THEN
                  INDX = 4
                  CHMSG = 'Invalid Gas Relative ' // 
     &              'Permeability Table'
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
  300       CONTINUE
            IRGTBL(2,IROCK) = NTBL
            IF( IRPGX.EQ.11 .OR. IRPGX.EQ.13 .OR. IRPGX.EQ.15 ) THEN
              CALL SPLINY( IRGTBL(1,IROCK),IRGTBL(2,IROCK) )
            ENDIF
          ENDIF
          GOTO 460
        ENDIF
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHRX.EQ.1 .OR. ISCHRX.EQ.8 .OR.  ISCHRX.EQ.11 .OR.
     &    ISCHRX.EQ.13 .OR. ISCHRX.EQ.15 .OR.
     &    ISCHRX.EQ.17 .OR. ISCHRX.EQ.101 .OR.
     &    (ISCHRX.GE.31 .AND. ISCHRX.LE.34) ) THEN
          IF( INDEX(ADUM(1:),'constant').NE.0 ) THEN
            IRPGX = 0
            VARB = 'Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 1
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 308 N = 1,NFLD
                RPGC(3,IZ(N)) = 1.D+0-1.D+0/SCHR(3,IZ(N))
  308         CONTINUE
            ELSE
              RPGCX(3) = 1.D+0-1.D+0/SCHR(3,IROCK)
            ENDIF
            VARB = 'Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 310 N = 1,NFLD
                IF( RPGC(3,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  310         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 1/n'
              WRITE(IWR,'(4X,A,1PE11.4)')'m Parameter: ',RPGCX(3)
              IF( RPGCX(3).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0  ) THEN
            IRPGX = 2
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 312 N = 1,NFLD
                RPGC(3,IZ(N)) = 1.D+0-2.D+0/SCHR(3,IZ(N))
  312         CONTINUE
            ELSE
              RPGCX(3) = 1.D+0-2.D+0/SCHR(3,IROCK)
            ENDIF
            VARB = 'Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 314 N = 1,NFLD
                IF( RPGC(3,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  314         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 2/n'
              WRITE(IWR,'(4X,A,1PE11.4)')'m Parameter: ',RPGCX(3)
              IF( RPGCX(3).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'sandia').NE.0 ) THEN
            IRPGX = 8
            WRITE (IWR,'(2X,A)') 'Sandia Relative Permeability Model'
          ELSEIF( INDEX(ADUM(1:),'variable corey').NE.0 ) THEN
            IRPGX = 20
            WRITE (IWR,'(2X,A)')'Variable Corey Relative ' // 
     &        'Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(1)
            ENDIF
            VARB = 'Irreducible Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(3)
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(2)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'extended').NE.0 .AND.
     &      INDEX(ADUM(1:),'power').NE.0 .AND.
     &      INDEX(ADUM(1:),'law').NE.0 ) THEN
            IRPGX = 21
            WRITE (IWR,'(2X,A)')'Extended Power Law Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(1)
            ENDIF
            VARB = 'Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(3)
            ENDIF
            VARB = 'Residual Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'free corey').NE.0 ) THEN
            IRPGX = 7
            WRITE (IWR,'(2X,A)')'Free Corey Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Residual Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'classical corey').NE.0 ) THEN
            IRPGX = 17
            WRITE (IWR,'(2X,A)')'Classical Corey Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Residual Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'doughty').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'drainage').NE.0 .AND. 
     &        INDEX(ADUM(1:),'imbibition').NE.0 ) THEN
              IRPGX = 19
              WRITE (IWR,'(2X,A)')'Doughty Drainage-Imbibition ' //
     &          'Relative Permeability Model'
           ELSE
              IRPGX = 18
              WRITE (IWR,'(2X,A)')'Doughty Relative ' //
     &          'Permeability Model'
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = '''Gamma'' Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 322 N = 1,NFLD
                IF( RPGC(2,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative ''Gamma'' Exponent Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  322         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
              IF( RPGCX(2).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative ''Gamma'' Exponent Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
            IF( IRPGX.EQ.19 ) THEN
              VARB = 'Drainage Residual Aqueous Saturation'
            ELSE
              VARB = 'Residual Aqueous Saturation'
            ENDIF
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = '''m'' Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 324 N = 1,NFLD
                IF( RPGC(4,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative ''m'' Exponent Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  324         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(4)
              IF( RPGCX(4).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative ''m'' Exponent Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
            IF( IRPGX.EQ.19 ) THEN
              VARB = 'Imbibition Residual Aqueous Saturation'
              IF( IJK.GT.0 ) THEN
                INDX = 5
                LNDX = LRPGC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(5))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPGCX(5)
              ENDIF
              VARB = 'Imbibition Actual Trapped Gas Saturation'
              IF( IJK.GT.0 ) THEN
                INDX = 6
                LNDX = LRPGC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(6))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPGCX(6)
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'corey').NE.0 ) THEN
            IRPGX = 3
            WRITE (IWR,'(2X,A)')'Corey Relative Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 328 N = 1,NFLD
                RPGC(1,IZ(N)) = SCHR(4,IZ(N))
  328         CONTINUE
            ELSE
              RPGCX(1) = SCHR(4,IROCK)
            ENDIF
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Irreducible Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 330 N = 1,NFLD
                RPGC(2,IZ(N)) = 1.D+0
  330         CONTINUE
            ELSE
              RPGCX(2) = 1.D+0
            ENDIF
            CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
            IF( INDX.EQ.1 ) THEN
              IF( IJK.GT.0 ) THEN
                INDX = 2
                LNDX = LRPGC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPGCX(2)
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'fatt and klikoff').NE.0 ) THEN
            IRPGX = 4
            WRITE (IWR,'(2X,2A)') 'Fatt and Klikoff Relative ',
     &        'Permeability Function'
          ELSEIF( INDEX(ADUM(1:),'stone').NE.0 ) THEN
            IRPGX = 5
            WRITE (IWR,'(2X,A)')'Stone Gas Relative ' //
     &        'Permeability Model'
            VARB = 'Stone (Slr)'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Stone (Sgr)'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Stone (n)'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'exponential').NE.0 ) THEN
            IRPGX = 6
            VARB = 'Gas Relative Permeability Function Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            IRPGX = 9
            WRITE (IWR,'(2X,A)')'van Genuchten Gas Relative ' // 
     &        'Permeability Model'
            IF( IJK.GT.0 ) THEN
              DO 332 N = 1,NFLD
                RPGC(1,IZ(N)) = SCHR(4,IZ(N))
  332         CONTINUE
            ELSE
              RPGCX(1) = SCHR(4,IROCK)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Irreducible Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSE
            INDX = 4
            NCH = INDEX( ADUM(1:),'  ' )-1
            CHMSG = 'Unrecognized Relative Perm. Function: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHRX.EQ.2 .OR. ISCHRX.EQ.6 .OR. ISCHRX.EQ.12 .OR.
     &    ISCHRX.EQ.14 .OR. ISCHRX.EQ.16 .OR.
     &    ISCHRX.EQ.18 .OR. ISCHRX.EQ.102 .OR.
     &    (ISCHRX.GE.35 .AND. ISCHRX.LE.38) ) THEN
          IF( INDEX(ADUM(1:),'constant').NE.0 ) THEN
            IRPGX = 0
            VARB = 'Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 1
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 338 N = 1,NFLD
                RPGC(3,IZ(N)) = SCHR(3,IZ(N))
  338         CONTINUE
            ELSE
            RPGCX(3) = SCHR(3,IROCK)
            ENDIF
            VARB = 'Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 2
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 340 N = 1,NFLD
                RPGC(3,IZ(N)) = SCHR(3,IZ(N))
  340         CONTINUE
            ELSE
              RPGCX(3) = SCHR(3,IROCK)
            ENDIF
            VARB = 'Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'variable corey').NE.0 ) THEN
            IRPGX = 20
            WRITE (IWR,'(2X,A)')'Variable Corey Relative ' // 
     &        'Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(1)
            ENDIF
            VARB = 'Irreducible Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(3)
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(2)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'extended').NE.0 .AND.
     &      INDEX(ADUM(1:),'power').NE.0 .AND.
     &      INDEX(ADUM(1:),'law').NE.0 ) THEN
            IRPGX = 21
            WRITE (IWR,'(2X,A)')'Extended Power Law Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(1)
            ENDIF
            VARB = 'Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(3)
            ENDIF
            VARB = 'Residual Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'free corey').NE.0 ) THEN
            IRPGX = 7
            WRITE (IWR,'(2X,A)')'Free Corey Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Residual Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'classical corey').NE.0 ) THEN
            IRPGX = 17
            WRITE (IWR,'(2X,A)')'Classical Corey Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Residual Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'doughty').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'drainage').NE.0 .AND. 
     &        INDEX(ADUM(1:),'imbibition').NE.0 ) THEN
              IRPGX = 19
              WRITE (IWR,'(2X,A)')'Doughty Drainage-Imbibition ' //
     &          'Relative Permeability Model'
           ELSE
              IRPGX = 18
              WRITE (IWR,'(2X,A)')'Doughty Relative ' //
     &          'Permeability Model'
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = '''Gamma'' Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 342 N = 1,NFLD
                IF( RPGC(2,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative ''Gamma'' Exponent Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  342         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
              IF( RPGCX(2).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative ''Gamma'' Exponent Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
            IF( IRPGX.EQ.19 ) THEN
              VARB = 'Drainage Residual Aqueous Saturation'
            ELSE
              VARB = 'Residual Aqueous Saturation'
            ENDIF
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = '''m'' Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 344 N = 1,NFLD
                IF( RPGC(4,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative ''m'' Exponent Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  344         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(4)
              IF( RPGCX(4).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative ''m'' Exponent Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
            IF( IRPGX.EQ.19 ) THEN
              VARB = 'Imbibition Residual Aqueous Saturation'
              IF( IJK.GT.0 ) THEN
                INDX = 5
                LNDX = LRPGC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(5))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPGCX(5)
              ENDIF
              VARB = 'Imbibition Actual Trapped Gas Saturation'
              IF( IJK.GT.0 ) THEN
                INDX = 6
                LNDX = LRPGC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(6))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPGCX(6)
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'corey').NE.0 ) THEN
            IRPGX = 3
            WRITE (IWR,'(2X,A)')'Corey Relative Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 348 N = 1,NFLD
                RPGC(1,IZ(N)) = SCHR(4,IZ(N))
  348         CONTINUE
            ELSE
              RPGCX(1) = SCHR(4,IROCK)
            ENDIF
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Irreducible Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 350 N = 1,NFLD
                RPGC(2,IZ(N)) = 1.D+0
  350         CONTINUE
            ELSE
              RPGCX(2) = 1.D+0
            ENDIF
            CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
            IF( INDX.EQ.1 ) THEN
              IF( IJK.GT.0 ) THEN
                INDX = 2
                LNDX = LRPGC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPGCX(2)
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'fatt and klikoff').NE.0 ) THEN
            IRPGX = 4
            WRITE (IWR,'(2X,2A)') 'Fatt and Klikoff Relative ',
     &        'Permeability Function'
          ELSEIF( INDEX(ADUM(1:),'stone').NE.0 ) THEN
            IRPGX = 5
            WRITE (IWR,'(2X,A)')'Stone Gas Relative ' //
     &        'Permeability Model'
            VARB = 'Stone (Slr)'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Stone (Sgr)'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Stone (n)'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'exponential').NE.0 ) THEN
            IRPGX = 6
            VARB = 'Gas Relative Permeability Function Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'sandia').NE.0 ) THEN
            IRPGX = 8
            WRITE (IWR,'(2X,A)') 'Sandia Relative Permeability Model'
          ELSEIF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            IRPGX = 9
            WRITE (IWR,'(2X,A)')'van Genuchten Gas Relative ' // 
     &        'Permeability Model'
            IF( IJK.GT.0 ) THEN
              DO 352 N = 1,NFLD
                RPGC(1,IZ(N)) = SCHR(4,IZ(N))
  352         CONTINUE
            ELSE
              RPGCX(1) = SCHR(4,IROCK)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Irreducible Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSE
            INDX = 4
            NCH = INDEX( ADUM(1:),'  ' )-1
            CHMSG = 'Unrecognized Relative Perm. Function: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Dual Porosity van Genuchten saturation function  ---
!
        ELSEIF( ISCHRX.EQ.3 ) THEN
          IF( INDEX(ADUM(1:),'constant').NE.0 ) THEN
            IRPGX = 0
            VARB = 'Matrix Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Fracture Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem').NE.0 .AND.
     &      INDEX(ADUM(1:),'linear fracture').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 31
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 360 N = 1,NFLD
                RPGC(3,IZ(N)) = 1.D+0-1.D+0/SCHR(3,IZ(N))
  360         CONTINUE
            ELSE
            RPGCX(3) = 1.D+0-1.D+0/SCHR(3,IROCK)
            ENDIF
            VARB = 'Matrix Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 362 N = 1,NFLD
                IF( RPGC(3,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  362         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 1/n'
              WRITE(IWR,'(4X,A,1PE11.4)')'m Parameter: ',RPGCX(3)
              IF( RPGCX(3).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 .AND.
     &      INDEX(ADUM(1:),'linear fracture').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 32
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 364 N = 1,NFLD
                RPGC(3,IZ(N)) = 1.D+0-2.D+0/SCHR(3,IZ(N))
  364         CONTINUE
            ELSE
            RPGCX(3) = 1.D+0-2.D+0/SCHR(3,IROCK)
            ENDIF
            VARB = 'Matrix Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 366 N = 1,NFLD
                IF( RPGC(3,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  366         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 2/n'
              WRITE(IWR,'(4X,A,1PE11.4)')'m Parameter: ',RPGCX(3)
              IF( RPGCX(3).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 1
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 368 N = 1,NFLD
                RPGC(3,IZ(N)) = 1.D+0-1.D+0/SCHR(3,IZ(N))
  368         CONTINUE
            ELSE
            RPGCX(3) = 1.D+0-1.D+0/SCHR(3,IROCK)
            ENDIF
            VARB = 'Matrix Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 370 N = 1,NFLD
                IF( RPGC(3,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  370         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 1/n'
              WRITE(IWR,'(4X,A,1PE11.4)')'m Parameter: ',RPGCX(3)
              IF( RPGCX(3).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 378 N = 1,NFLD
                RPGC(1,IZ(N)) = 1.D+0-1.D+0/SCHR(6,IZ(N))
  378         CONTINUE
            ELSE
            RPGCX(1) = 1.D+0-1.D+0/SCHR(6,IROCK)
            ENDIF
            VARB = 'Fracture Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 380 N = 1,NFLD
                IF( RPGC(3,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  380         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 1/n'
              WRITE(IWR,'(4X,A,1PE11.4)')'m Parameter: ',RPGCX(1)
              IF( RPGCX(3).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 2
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 388 N = 1,NFLD
                RPGC(3,IZ(N)) = 1.D+0-2.D+0/SCHR(3,IZ(N))
  388         CONTINUE
            ELSE
              RPGCX(3) = 1.D+0-2.D+0/SCHR(3,IROCK)
            ENDIF
            VARB = 'Matrix Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 390 N = 1,NFLD
                IF( RPGC(3,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  390         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 2/n'
              WRITE(IWR,'(4X,A,1PE11.4)')'m Parameter: ',RPGCX(3)
              IF( RPGCX(3).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 398 N = 1,NFLD
                RPGC(1,IZ(N)) = 1.D+0-2.D+0/SCHR(6,IZ(N))
  398         CONTINUE
            ELSE
              RPGCX(1) = 1.D+0-2.D+0/SCHR(6,IROCK)
            ENDIF
            VARB = 'Fracture Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 400 N = 1,NFLD
                IF( RPGC(3,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  400         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 2/n'
              WRITE(IWR,'(4X,A,1PE11.4)')'m Parameter: ',RPGCX(3)
              IF( RPGCX(3).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
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
            IRPGX = 0
            VARB = 'Matrix Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Fracture Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem').NE.0 .AND.
     &      INDEX(ADUM(1:),'linear fracture').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 31
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 408 N = 1,NFLD
                RPGC(3,IZ(N)) = SCHR(3,IZ(N))
  408         CONTINUE
            ELSE
            RPGCX(3) = SCHR(3,IROCK)
            ENDIF
            VARB = 'Matrix Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 .AND.
     &      INDEX(ADUM(1:),'linear fracture').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 32
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 410 N = 1,NFLD
                RPGC(3,IZ(N)) = SCHR(3,IZ(N))
  410         CONTINUE
            ELSE
              RPGCX(3) = SCHR(3,IROCK)
            ENDIF
            VARB = 'Matrix Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 1
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 412 N = 1,NFLD
                RPGC(3,IZ(N)) = SCHR(3,IZ(N))
  412         CONTINUE
            ELSE
              RPGCX(3) = SCHR(3,IROCK)
            ENDIF
            VARB = 'Matrix Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPGCX(3)
            ENDIF
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 416 N = 1,NFLD
                RPGC(1,IZ(N)) = SCHR(6,IZ(N))
  416         CONTINUE
            ELSE
            RPGCX(1) = SCHR(6,IROCK)
            ENDIF
            VARB = 'Fracture Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPGCX(1)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').EQ.0 .AND.
     &      INDEX(ADUM(1:),'bc').EQ.0 ) THEN
            IRPGX = 2
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 428 N = 1,NFLD
                RPGC(3,IZ(N)) = SCHR(3,IZ(N))
  428         CONTINUE
            ELSE
            RPGCX(3) = SCHR(3,IROCK)
            ENDIF
            VARB = 'Matrix Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPGCX(3)
            ENDIF
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 430 N = 1,NFLD
                RPGC(1,IZ(N)) = SCHR(6,IZ(N))
  430         CONTINUE
            ELSE
              RPGCX(1) = SCHR(6,IROCK)
            ENDIF
            VARB = 'Fracture Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPGCX(1)
            ENDIF
          ELSE
            INDX = 4
            NCH = INDEX( ADUM(1:),'  ' )-1
            CHMSG = 'Unrecognized Relative Perm. Function: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
!
!---  Unknown saturation function  ---
!
        ELSE
          IF( INDEX(ADUM(1:),'constant').NE.0 ) THEN
            IRPGX = 0
            VARB = 'Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').NE.0 ) THEN
            IRPGX = 41
            IDFLT = 0
            VARB = 'Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO N = 1,NFLD
                IF( RPGC(3,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
              ENDDO
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 1/n'
              WRITE(IWR,'(4X,A,1PE11.4)')'m Parameter: ',RPGCX(3)
              IF( RPGCX(3).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mualem').NE.0 .AND.
     &      INDEX(ADUM(1:),'bc').NE.0 ) THEN
            IRPGX = 42
            IDFLT = 0
            VARB = 'Mualem Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 .AND.
     &      INDEX(ADUM(1:),'vg').NE.0  ) THEN
            IRPGX = 43
            IDFLT = 0
            VARB = 'Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO N = 1,NFLD
                IF( RPGC(3,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative van Genuchten ''m'' Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
              ENDDO
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A)') 'Default Value: m = 1 - 2/n'
              WRITE(IWR,'(4X,A,1PE11.4)')'m Parameter: ',RPGCX(3)
              IF( RPGCX(3).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative van Genuchten ''m'' Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 .AND.
     &      INDEX(ADUM(1:),'bc').NE.0 ) THEN
            IRPGX = 44
            IDFLT = 0
            VARB = 'Burdine Porosity Distribution Model'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE(IWR,'(2X,A)') VARB(1:IVR)
              WRITE(IWR,'(4X,A,1PE11.4)') 'Lambda Parameter: ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'variable corey').NE.0 ) THEN
            IRPGX = 20
            WRITE (IWR,'(2X,A)')'Variable Corey Relative ' // 
     &        'Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(1)
            ENDIF
            VARB = 'Irreducible Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(3)
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(2)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'extended').NE.0 .AND.
     &      INDEX(ADUM(1:),'power').NE.0 .AND.
     &      INDEX(ADUM(1:),'law').NE.0 ) THEN
            IRPGX = 21
            WRITE (IWR,'(2X,A)')'Extended Power Law Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(1)
            ENDIF
            VARB = 'Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(3)
            ENDIF
            VARB = 'Residual Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'free corey').NE.0 ) THEN
            IRPGX = 7
            WRITE (IWR,'(2X,A)')'Free Corey Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Residual Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'classical corey').NE.0 ) THEN
            IRPGX = 17
            WRITE (IWR,'(2X,A)')'Classical Corey Relative ' //
     &        'Permeability Model'
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Residual Aqueous Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Residual Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(4)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'doughty').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'drainage').NE.0 .AND. 
     &        INDEX(ADUM(1:),'imbibition').NE.0 ) THEN
              IRPGX = 19
              WRITE (IWR,'(2X,A)')'Doughty Drainage-Imbibition ' //
     &          'Relative Permeability Model'
           ELSE
              IRPGX = 18
              WRITE (IWR,'(2X,A)')'Doughty Relative ' //
     &          'Permeability Model'
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = '''Gamma'' Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 432 N = 1,NFLD
                IF( RPGC(2,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative ''Gamma'' Exponent Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  432         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
              IF( RPGCX(2).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative ''Gamma'' Exponent Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
            IF( IRPGX.EQ.19 ) THEN
              VARB = 'Drainage Residual Aqueous Saturation'
            ELSE
              VARB = 'Residual Aqueous Saturation'
            ENDIF
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = '''m'' Exponent Gas Relative Permeability'
            IF( IJK.GT.0 ) THEN
              INDX = 4
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              DO 434 N = 1,NFLD
                IF( RPGC(4,IZ(N)).LT.EPSL ) THEN
                  INDX = 4
                  CHMSG = 'Negative ''m'' Exponent Parameter'
                  CALL WRMSGS( INDX )
                ENDIF
  434         CONTINUE
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(4))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(4)
              IF( RPGCX(4).LT.EPSL ) THEN
                INDX = 4
                CHMSG = 'Negative ''m'' Exponent Parameter'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
            IF( IRPGX.EQ.19 ) THEN
              VARB = 'Imbibition Residual Aqueous Saturation'
              IF( IJK.GT.0 ) THEN
                INDX = 5
                LNDX = LRPGC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(5))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPGCX(5)
              ENDIF
              VARB = 'Imbibition Actual Trapped Gas Saturation'
              IF( IJK.GT.0 ) THEN
                INDX = 6
                LNDX = LRPGC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(6))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPGCX(6)
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'corey').NE.0 ) THEN
            IRPGX = 3
            WRITE (IWR,'(2X,A)')'Corey Relative Permeability Model'
            VARB = 'Irreducible Aqueous Saturation'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 438 N = 1,NFLD
                RPGC(1,IZ(N)) = SCHR(4,IZ(N))
  438         CONTINUE
            ELSE
              RPGCX(1) = SCHR(4,IROCK)
            ENDIF
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Irreducible Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
            VARB = 'Endpoint Gas Relative Permeability'
            IDFLT = 1
            IF( IJK.GT.0 ) THEN
              DO 440 N = 1,NFLD
                RPGC(2,IZ(N)) = 1.D+0
  440         CONTINUE
            ELSE
              RPGCX(2) = 1.D+0
            ENDIF
            CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
            IF( INDX.EQ.1 ) THEN
              IF( IJK.GT.0 ) THEN
                INDX = 2
                LNDX = LRPGC
                UNTS = 'null'
                CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
              ELSE
                CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
                WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &            RPGCX(2)
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'fatt and klikoff').NE.0 ) THEN
            IRPGX = 4
            WRITE (IWR,'(2X,A)') 'Fatt and Klikoff Relative ',
     &        'Permeability Function'
          ELSEIF( INDEX(ADUM(1:),'stone').NE.0 ) THEN
            IRPGX = 5
            WRITE (IWR,'(2X,A)')'Stone Gas Relative ' //
     &        'Permeability Model'
            VARB = 'Stone (Slr)'
            IF( IJK.GT.0 ) THEN
              INDX = 1
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(1))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(1)
            ENDIF
            VARB = 'Stone (Sgr)'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Stone (n)'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'exponential').NE.0 ) THEN
            IRPGX = 6
            VARB = 'Gas Relative Permeability Function Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'sandia').NE.0 ) THEN
            IRPGX = 8
            WRITE (IWR,'(2X,A)') 'Sandia Relative Permeability Model'
          ELSEIF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            IRPGX = 9
            WRITE (IWR,'(2X,A)')'van Genuchten Gas Relative ' // 
     &        'Permeability Model'
            IF( IJK.GT.0 ) THEN
              DO 442 N = 1,NFLD
                RPGC(1,IZ(N)) = SCHR(4,IZ(N))
  442         CONTINUE
            ELSE
              RPGCX(1) = SCHR(4,IROCK)
            ENDIF
            VARB = 'Exponent'
            IF( IJK.GT.0 ) THEN
              INDX = 2
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(2))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(2)
            ENDIF
            VARB = 'Irreducible Gas Saturation'
            IF( IJK.GT.0 ) THEN
              INDX = 3
              LNDX = LRPGC
              UNTS = 'null'
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
            ELSE
              CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(3))
              WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &          RPGCX(3)
            ENDIF
          ELSE
            INDX = 4
            NCH = INDEX( ADUM(1:),'  ' )-1
            CHMSG = 'Unrecognized Relative Perm. Function: '
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
        IF( INDEX(ADUM(1:), 'mualem').NE.0 ) THEN
          IRPGX = 1
          WRITE(IWR,'(2X,A)') 'Mualem Porosity Distribution Model'
        ELSEIF( INDEX(ADUM(1:),'burdine').NE.0 ) THEN
          IRPGX = 2
          WRITE(IWR,'(2X,A)') 'Burdine Porosity Distribution Model'
        ELSE
          INDX = 4
          NCH = INDEX( ADUM(1:),'  ' )-1
          CHMSG = 'Unrecognized Relative Perm. Function: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Klinkenberg effect  ---
!
      CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
      IF( INDX.EQ.1 ) THEN
        VARB = 'Gas Relative Permeability Extension'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
        IF( INDEX(ADUM(1:), 'klink').NE.0 ) THEN
          IRPGX = IRPGX + 100
          WRITE(IWR,'(2X,A)') 'Klinkenberg Effect Model'
          VARB = 'Klinkenberg Scaling Parameter (C1)'
          IF( IJK.GT.0 ) THEN
            INDX = 5
            LNDX = LRPGC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(5))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        RPGCX(5)
          ENDIF
          VARB = 'Klinkenberg Expontential Parameter (C2)'
          IF( IJK.GT.0 ) THEN
            INDX = 6
            LNDX = LRPGC
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,RPGC,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RPGCX(6))
            WRITE (IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        RPGCX(6)
          ENDIF
          VARB = 'Klinkenberg Pressure Units'
          CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,UNTS )
          VAR = 1.D+0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          INDX = 0
          CALL RDUNIT( UNTS,VAR,INDX )
          IF( IJK.GT.0 ) THEN
            DO 450 N = 1,NFLD
              RPGC(5,N) = RPGC(5,N)*(VAR**(1.D+0-RPGC(6,N)))
  450       CONTINUE
          ELSE
            RPGCX(5) = RPGCX(5)*(VAR**(1.D+0-RPGCX(6)))
          ENDIF
!        ELSE
!          INDX = 4
!          NCH = INDEX( ADUM(1:),'  ' )-1
!          CHMSG = 'Unrecognized Gas Relative Permeability Extension: '
!     &      // ADUM(1:NCH)
!          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Translate gas relative permeability type for IJK indexing  ---
!
  460 CONTINUE
      IF( IJK.GT.0 ) THEN
        DO 462 N = 1,NFLD
          IRPG(IZ(N)) = IRPGX
  462   CONTINUE
!
!---  For IJK indexing with the rock/soil option, correlate rock/soil 
!     numbers with nodes for gas relative permeability type and
!     parameters  ---
!
      ELSEIF( IJK.LT.0 ) THEN
        DO 472 N = 1,NFLD
          IF( IVAR(N).EQ.NR ) THEN
            DO 470 L = 1,LRPGC
              RPGC(L,N) = RPGCX(L)
  470       CONTINUE
            IRPG(N) = IRPGX
          ENDIF
  472   CONTINUE
!
!---  For rock/soil zonation input translate gas relative permeability 
!     type and parameters  ---
!
      ELSE
        IRPG(IROCK) = IRPGX
        DO 480 L = 1,LRPGC
          RPGC(L,IROCK) = RPGCX(L)
  480   CONTINUE
      ENDIF
!
!---  Loop over remaining rock/soils within scaling group  ---
!
      IF( ISLC(19).EQ.1 .AND. IROCK.LT.NROCK ) THEN
        DO 494 M = IROCK+1,NROCK
          IF( ISCALE(M).EQ.ISGRP ) THEN
            NR = NR + 1
            IRPG(M) = IRPG(IROCK)
            DO 490 L = 1,LRPGC
              RPGC(L,M) = RPGC(L,IROCK)
  490       CONTINUE
            DO 492 L = 1,2
              IRGTBL(L,M) = IRGTBL(L,IROCK)
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
!---  End of RDGRP group ---
!
      RETURN
      END

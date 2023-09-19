!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIJK( ISTART,IJK,CHDUM,UNTS,VAR )
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
!     Read dimensionless property values from an external file for
!     IJK, JKI, or KIJ indexed inputs
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 March 2002.
!     Last Modified by MD White, PNNL, 19 March 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,CDUM,FDUM,UNTS,UNTSX
      CHARACTER*132 VARBX
      CHARACTER*512 CHDUM,CHDUMX
      REAL*8 VAR(*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARL
      LOGICAL FCHK,FBIN
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIJK'
      ALLOCATE( VARL(1:LRC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: VARL'
        CALL WRMSGS( INDX )
      ENDIF
      ISTX = ISTART
      VARBX = VARB
      IDFLT = 0
      VARB = 'External File Name'
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      FBIN = .FALSE.
      IF( INDEX(ADUM(1:),'binary').NE.0 .OR.
     &  INDEX(ADUM(1:),'bfile').NE.0 .OR.
     &  INDEX(ADUM(1:),'b_file').NE.0 ) FBIN = .TRUE.
!
!---  External-file input  ---
!
      IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
        ICOLON = INDEX(ADUM(1:),':') + 1
        FDUM = ADUM(ICOLON:NCHA)
        NCHF = NCHA-ICOLON+1
        INQUIRE( FILE=FDUM(1:NCHF), FORM=CDUM, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          INDX = 4
          CHMSG = 'IJK Indexing file does not exist: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'UNFORMATTED' .AND. (.NOT. FBIN) ) THEN
          INDX = 4
          CHMSG = 'Formatted IJK Indexing file is unformatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'FORMATTED' .AND. FBIN ) THEN
          INDX = 4
          CHMSG = 'Unformatted IJK Indexing file is formatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        END IF
        IF( FBIN ) THEN
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='UNFORMATTED')
        ELSE
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='FORMATTED')
        ENDIF
        WRITE(IWR,'(/,2A)') 'IJK Indexing File: ',FDUM(1:NCHF)
!
!---    Check for units  ---
!
        VARB = VARBX
        VARX = 1.D+0
        NCHU = INDEX( UNTS(1:),'  ' ) - 1
        IF( UNTS(1:NCHU).NE.'null' ) THEN
          IDFLT = 1
          UNTSX = UNTS
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTSX)
          WRITE(IWR,'(2X,3A)') VARB(1:IVR),', ',UNTSX(1:NCH)
          INDX = 0
          CALL RDUNIT(UNTSX,VARX,INDX)
        ENDIF
        IF( FBIN ) THEN
          READ(26) (VARL(IROCK),IROCK=1,NFLD)
        ELSE
          READ(26,'(A)') CHDUMX
          IF( INDEX( CHDUMX(1:),'@' ).NE.0 ) THEN
            REWIND(26)
            NC = 0
            DO
              ISTARTX = 1
              READ(26,'(A)',END=90) CHDUMX
              IATX = INDEX( CHDUMX(1:), '@' )
              ICMX = INDEX( CHDUMX(1:), '  ' )
              CHDUMX(IATX:IATX) = ','
              CHDUMX(ICMX:ICMX) = ','
              CALL RDINT(ISTARTX,ICOMMAX,CHDUMX,IATX )
              MC = NC + 1
              NC = MIN( NC+IATX,NFLD )
              CALL RDDPR(ISTARTX,ICOMMAX,CHDUMX,VARLX )
              DO IROCK = MC,NC
                VARL(IROCK) = VARLX
              ENDDO
              IF( NC.GE.NFLD ) EXIT
            ENDDO
   90       IF( NC.LT.NFLD ) THEN
              INDX = 4
              CHMSG = 'End of IJK Indexing File Encountered: ' // 
     &          FDUM(1:NCHF)
              CALL WRMSGS( INDX )
            ENDIF
          ELSE
            REWIND(26)
            READ(26,*) (VARL(IROCK),IROCK=1,NFLD)
          ENDIF
        ENDIF
!
!---    IJK indexing  ---
!
        IF( IJK.EQ.1 ) THEN
          NC = 0
          DO K = 1,KFLD
          DO J = 1,JFLD
          DO I = 1,IFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            VAR(IROCK) = VARL(NC)*VARX
          ENDDO
          ENDDO
          ENDDO
!
!---    JKI indexing  ---
!
        ELSEIF( IJK.EQ.2 ) THEN
          NC = 0
          DO I = 1,IFLD
          DO K = 1,KFLD
          DO J = 1,JFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            VAR(IROCK) = VARL(NC)*VARX
          ENDDO
          ENDDO
          ENDDO
!
!---    KIJ indexing  ---
!
        ELSEIF( IJK.EQ.3 ) THEN
          NC = 0
          DO J = 1,JFLD
          DO I = 1,IFLD
          DO K = 1,KFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            VAR(IROCK) = VARL(NC)*VARX
          ENDDO
          ENDDO
          ENDDO
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Indexing Option' // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        CLOSE(UNIT=26)
!
!---  Input-file input  ---
!
      ELSE
        ISTART = ISTX
        VARB = VARBX
        IDFLT = 1
        VARX = 0.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VARX)
        IDFLTDX = IDFLTD
!
!---    Check for units  ---
!
        NCHU = INDEX( UNTS(1:),'  ' ) - 1
        IF( UNTS(1:NCHU).NE.'null' ) THEN
          IDFLT = 1
          UNTSX = UNTS
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTSX)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTSX(1:NCH),
     &      ': ',VARX
          INDX = 0
          CALL RDUNIT(UNTSX,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,3A)') ' (',VARX,', ',UNTS(1:NCHU),')'
        ELSE
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VARX
        ENDIF
        IF( IDFLTDX.EQ.0 ) THEN
          DO 400 N = 1,NFLD
            VAR(IZ(N)) = VARX
 400      CONTINUE
        ENDIF
      ENDIF
      IF( ALLOCATED(VARL) ) THEN
      DEALLOCATE( VARL,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: VARL'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
!
!---  Reset subroutine name string  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIJK group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIJKC( ISTART,IJK,CHDUM,UNTS,VAR,INDC,LNDX )
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
!     Read dimensionless property values from an external file for
!     IJK indexed inputs with double indices independent of the number
!     of rock types
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 March 2002.
!     Last Modified by MD White, PNNL, 19 March 2002.
!     Last Modified by VL Freedman, PNNL, 4 May 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,CDUM,FDUM,UNTS,UNTSX
      CHARACTER*132 VARBX
      CHARACTER*512 CHDUM,CHDUMX
      REAL*8 VAR(LNDX,*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARL
      LOGICAL FCHK,FBIN
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIJKC'
      ALLOCATE( VARL(1:LRC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: VARL'
        CALL WRMSGS( INDX )
      ENDIF
      ISTX = ISTART
      VARBX = VARB
      IDFLT = 0
      VARB = 'External File Name'
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      FBIN = .FALSE.
!
!---  External-file input  ---
!
      IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
        ICOLON = INDEX(ADUM(1:),':') + 1
        FDUM = ADUM(ICOLON:NCHA)
        NCHF = NCHA-ICOLON+1
        INQUIRE( FILE=FDUM(1:NCHF), FORM=CDUM, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          INDX = 4
          CHMSG = 'IJK Indexing file does not exist: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'UNFORMATTED' .AND. (.NOT. FBIN) ) THEN
          INDX = 4
          CHMSG = 'Formatted IJK Indexing file is unformatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'FORMATTED' .AND. FBIN ) THEN
          INDX = 4
          CHMSG = 'Unformatted IJK Indexing file is formatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        END IF
        IF( FBIN ) THEN
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='UNFORMATTED')
        ELSE
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='FORMATTED')
        ENDIF
        WRITE(IWR,'(/,2A)') 'IJK Indexing File: ',FDUM(1:NCHF)
!
!---    Check for units  ---
!
        VARB = VARBX
        VARX = 1.D+0
        NCHU = INDEX( UNTS(1:),'  ' ) - 1
        IF( UNTS(1:NCHU).NE.'null' ) THEN
          IDFLT = 1
          UNTSX = UNTS
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTSX)
          WRITE(IWR,'(2X,3A)') VARB(1:IVR),', ',UNTSX(1:NCH)
          INDX = 0
          CALL RDUNIT(UNTSX,VARX,INDX)
        ENDIF
        IF( FBIN ) THEN
          READ(26) (VAR(INDC,IROCK),IROCK=1,NFLD)
        ELSE
          READ(26,'(A)') CHDUMX
          IF( INDEX( CHDUMX(1:),'@' ).NE.0 ) THEN
            REWIND(26)
            NC = 0
            DO
              ISTARTX = 1
              READ(26,'(A)',END=90) CHDUMX
              IATX = INDEX( CHDUMX(1:), '@' )
              ICMX = INDEX( CHDUMX(1:), '  ' )
              CHDUMX(IATX:IATX) = ','
              CHDUMX(ICMX:ICMX) = ','
              CALL RDINT(ISTARTX,ICOMMAX,CHDUMX,IATX )
              MC = NC + 1
              NC = MIN( NC+IATX,NFLD )
              CALL RDDPR(ISTARTX,ICOMMAX,CHDUMX,VARLX )
              DO IROCK = MC,NC
                VARL(IROCK) = VARLX
              ENDDO
              IF( NC.GE.NFLD ) EXIT
            ENDDO
   90       IF( NC.LT.NFLD ) THEN
              INDX = 4
              CHMSG = 'End of IJK Indexing File Encountered: ' // 
     &          FDUM(1:NCHF)
              CALL WRMSGS( INDX )
            ENDIF
          ELSE
            REWIND(26)
            READ(26,*) (VARL(IROCK),IROCK=1,NFLD)
          ENDIF
        ENDIF
!
!---    IJK indexing  ---
!
        IF( IJK.EQ.1 ) THEN
          NC = 0
          DO K = 1,KFLD
          DO J = 1,JFLD
          DO I = 1,IFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            VAR(INDC,IROCK) = VARL(NC)*VARX
          ENDDO
          ENDDO
          ENDDO
!
!---    JKI indexing  ---
!
        ELSEIF( IJK.EQ.2 ) THEN
          NC = 0
          DO I = 1,IFLD
          DO K = 1,KFLD
          DO J = 1,JFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            VAR(INDC,IROCK) = VARL(NC)*VARX
          ENDDO
          ENDDO
          ENDDO
!
!---    KIJ indexing  ---
!
        ELSEIF( IJK.EQ.3 ) THEN
          NC = 0
          DO J = 1,JFLD
          DO I = 1,IFLD
          DO K = 1,KFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            VAR(INDC,IROCK) = VARL(NC)*VARX
          ENDDO
          ENDDO
          ENDDO
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Indexing Option' // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        CLOSE(UNIT=26)
!
!---  Input-file input  ---
!
      ELSE
        ISTART = ISTX
        VARB = VARBX
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VARX)
        IDFLTDX = IDFLTD
!
!---    Check for units  ---
!
        NCHU = INDEX( UNTS(1:),'  ' ) - 1
        IF( UNTS(1:NCHU).NE.'null' ) THEN
          IDFLT = 1
          UNTSX = UNTS
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTSX)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTSX(1:NCH),
     &      ': ',VARX
          INDX = 0
          CALL RDUNIT(UNTSX,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,3A)') ' (',VARX,', ',UNTS(1:NCHU),')'
        ELSE
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VARX
        ENDIF
        IF( IDFLTDX.EQ.0 ) THEN
          DO 400 N = 1,NFLD
            VAR(INDC,IZ(N)) = VARX
  400     CONTINUE
        ENDIF
      ENDIF
      IF( ALLOCATED(VARL) ) THEN
      DEALLOCATE( VARL,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: VARL'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
!
!---  Reset subroutine name string  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIJKC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIJKD( ISTART,IJK,CHDUM,UNTS,VAR,INDC,LNDX )
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
!     Read dimensionless property values from an external file for
!     IJK, JKI, or KIJ indexed inputs with double indices
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 March 2002.
!     Last Modified by MD White, PNNL, 19 March 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,CDUM,FDUM,UNTS,UNTSX
      CHARACTER*132 VARBX
      CHARACTER*512 CHDUM,CHDUMX
      REAL*8 VAR(LNDX,*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARL
      LOGICAL FCHK,FBIN
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIJKD'
      ALLOCATE( VARL(1:LRC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: VARL'
        CALL WRMSGS( INDX )
      ENDIF
      ISTX = ISTART
      VARBX = VARB
      IDFLT = 0
      VARB = 'External File Name'
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      FBIN = .FALSE.
!
!---  External-file input  ---
!
      IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
        ICOLON = INDEX(ADUM(1:),':') + 1
        FDUM = ADUM(ICOLON:NCHA)
        NCHF = NCHA-ICOLON+1
        INQUIRE( FILE=FDUM(1:NCHF), FORM=CDUM, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          INDX = 4
          CHMSG = 'IJK Indexing file does not exist: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'UNFORMATTED' .AND. (.NOT. FBIN) ) THEN
          INDX = 4
          CHMSG = 'Formatted IJK Indexing file is unformatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'FORMATTED' .AND. FBIN ) THEN
          INDX = 4
          CHMSG = 'Unformatted IJK Indexing file is formatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        END IF
        IF( FBIN ) THEN
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='UNFORMATTED')
        ELSE
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='FORMATTED')
        ENDIF
        WRITE(IWR,'(/,2A)') 'IJK Indexing File: ',FDUM(1:NCHF)
!
!---    Check for units  ---
!
        VARB = VARBX
        VARX = 1.D+0
        NCHU = INDEX( UNTS(1:),'  ' ) - 1
        IF( UNTS(1:NCHU).NE.'null' ) THEN
          IDFLT = 1
          UNTSX = UNTS
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTSX)
          WRITE(IWR,'(2X,3A)') VARB(1:IVR),', ',UNTSX(1:NCH)
          INDX = 0
          CALL RDUNIT(UNTSX,VARX,INDX)
        ENDIF
        IF( FBIN ) THEN
          READ(26) (VARL(IROCK),IROCK=1,NFLD)
        ELSE
          READ(26,'(A)') CHDUMX
          IF( INDEX( CHDUMX(1:),'@' ).NE.0 ) THEN
            REWIND(26)
            NC = 0
            DO
              ISTARTX = 1
              READ(26,'(A)',END=90) CHDUMX
              IATX = INDEX( CHDUMX(1:), '@' )
              ICMX = INDEX( CHDUMX(1:), '  ' )
              CHDUMX(IATX:IATX) = ','
              CHDUMX(ICMX:ICMX) = ','
              CALL RDINT(ISTARTX,ICOMMAX,CHDUMX,IATX )
              MC = NC + 1
              NC = MIN( NC+IATX,NFLD )
              CALL RDDPR(ISTARTX,ICOMMAX,CHDUMX,VARLX )
              DO IROCK = MC,NC
                VARL(IROCK) = VARLX
              ENDDO
              IF( NC.GE.NFLD ) EXIT
            ENDDO
   90       IF( NC.LT.NFLD ) THEN
              INDX = 4
              CHMSG = 'End of IJK Indexing File Encountered: ' // 
     &          FDUM(1:NCHF)
              CALL WRMSGS( INDX )
            ENDIF
          ELSE
            REWIND(26)
            READ(26,*) (VARL(IROCK),IROCK=1,NFLD)
          ENDIF
        ENDIF
!
!---    IJK indexing  ---
!
        IF( IJK.EQ.1 ) THEN
          NC = 0
          DO K = 1,KFLD
          DO J = 1,JFLD
          DO I = 1,IFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            VAR(INDC,IROCK) = VARL(NC)*VARX
          ENDDO
          ENDDO
          ENDDO
!
!---    JKI indexing  ---
!
        ELSEIF( IJK.EQ.2 ) THEN
          NC = 0
          DO I = 1,IFLD
          DO K = 1,KFLD
          DO J = 1,JFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            VAR(INDC,IROCK) = VARL(NC)*VARX
          ENDDO
          ENDDO
          ENDDO
!
!---    KIJ indexing  ---
!
        ELSEIF( IJK.EQ.3 ) THEN
          NC = 0
          DO J = 1,JFLD
          DO I = 1,IFLD
          DO K = 1,KFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            VAR(INDC,IROCK) = VARL(NC)*VARX
          ENDDO
          ENDDO
          ENDDO
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Indexing Option' // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        CLOSE(UNIT=26)
!
!---  Input-file input  ---
!
      ELSE
        ISTART = ISTX
        VARB = VARBX
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VARX)
        IDFLTDX = IDFLTD
!
!---    Check for units  ---
!
        NCHU = INDEX( UNTS(1:),'  ' ) - 1
        IF( UNTS(1:NCHU).NE.'null' ) THEN
          IDFLT = 1
          UNTSX = UNTS
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTSX)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTSX(1:NCH),
     &      ': ',VARX
          INDX = 0
          CALL RDUNIT(UNTSX,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,3A)') ' (',VARX,', ',UNTS(1:NCHU),')'
        ELSE
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VARX
        ENDIF
        IF( IDFLTDX.EQ.0 ) THEN
          DO 400 N = 1,NFLD
            VAR(INDC,IZ(N)) = VARX
  400     CONTINUE
        ENDIF
      ENDIF
      IF( ALLOCATED(VARL) ) THEN
      DEALLOCATE( VARL,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: VARL'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
!
!---  Reset subroutine name string  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIJKD group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIJKI( ISTART,IJK,CHDUM,IVAR )
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
!     Read table numbers from an external file for
!     IJK, JKI, or KIJ indexed inputs.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 March 2013.
!     Last Modified by MD White, PNNL, 18 March 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,CDUM,FDUM
      CHARACTER*132 VARBX
      CHARACTER*512 CHDUM,CHDUMX
      INTEGER IVAR(*)
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVARL
      LOGICAL FCHK,FBIN
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIJKI'
      ALLOCATE( IVARL(1:LRC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IVARL'
        CALL WRMSGS( INDX )
      ENDIF
      ISTX = ISTART
      VARBX = VARB
      IDFLT = 0
      VARB = 'External File Name'
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      FBIN = .FALSE.
      IF( INDEX(ADUM(1:),'binary').NE.0 .OR.
     &  INDEX(ADUM(1:),'bfile').NE.0 .OR.
     &  INDEX(ADUM(1:),'b_file').NE.0 ) FBIN = .TRUE.
!
!---  External-file input  ---
!
      IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
        ICOLON = INDEX(ADUM(1:),':') + 1
        FDUM = ADUM(ICOLON:NCHA)
        NCHF = NCHA-ICOLON+1
        INQUIRE( FILE=FDUM(1:NCHF), FORM=CDUM, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          INDX = 4
          CHMSG = 'IJK Indexing file does not exist: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'UNFORMATTED' .AND. (.NOT. FBIN) ) THEN
          INDX = 4
          CHMSG = 'Formatted IJK Indexing file is unformatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'FORMATTED' .AND. FBIN ) THEN
          INDX = 4
          CHMSG = 'Unformatted IJK Indexing file is formatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        END IF
        IF( FBIN ) THEN
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='UNFORMATTED')
        ELSE
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='FORMATTED')
        ENDIF
        WRITE(IWR,'(/,2A)') 'IJK Indexing File: ',FDUM(1:NCHF)
!
!---    Read file for table numbers  ---
!
        IF( FBIN ) THEN
          READ(26) (IVARL(IROCK),IROCK=1,NFLD)
        ELSE
          READ(26,'(A)') CHDUMX
          IF( INDEX( CHDUMX(1:),'@' ).NE.0 ) THEN
            REWIND(26)
            NC = 0
            DO
              ISTARTX = 1
              READ(26,'(A)',END=90) CHDUMX
              IATX = INDEX( CHDUMX(1:), '@' )
              ICMX = INDEX( CHDUMX(1:), '  ' )
              CHDUMX(IATX:IATX) = ','
              CHDUMX(ICMX:ICMX) = ','
              CALL RDINT(ISTARTX,ICOMMAX,CHDUMX,IATX )
              MC = NC + 1
              NC = MIN( NC+IATX,NFLD )
              CALL RDINT(ISTARTX,ICOMMAX,CHDUMX,IVARLX )
              DO IROCK = MC,NC
                IVARL(IROCK) = IVARLX
              ENDDO
              IF( NC.GE.NFLD ) EXIT
            ENDDO
   90       IF( NC.LT.NFLD ) THEN
              INDX = 4
              CHMSG = 'End of IJK Indexing File Encountered: ' // 
     &          FDUM(1:NCHF)
              CALL WRMSGS( INDX )
            ENDIF
          ELSE
            REWIND(26)
            READ(26,*) (IVARL(IROCK),IROCK=1,NFLD)
          ENDIF
        ENDIF
!
!---    IJK indexing  ---
!
        IF( IJK.EQ.1 ) THEN
          NC = 0
          DO K = 1,KFLD
          DO J = 1,JFLD
          DO I = 1,IFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            IVAR(IROCK) = IVARL(NC)
          ENDDO
          ENDDO
          ENDDO
!
!---    JKI indexing  ---
!
        ELSEIF( IJK.EQ.2 ) THEN
          NC = 0
          DO I = 1,IFLD
          DO K = 1,KFLD
          DO J = 1,JFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            IVAR(IROCK) = IVARL(NC)
          ENDDO
          ENDDO
          ENDDO
!
!---    KIJ indexing  ---
!
        ELSEIF( IJK.EQ.3 ) THEN
          NC = 0
          DO J = 1,JFLD
          DO I = 1,IFLD
          DO K = 1,KFLD
            NC = NC + 1
            IROCK = (K-1)*IJFLD + (J-1)*IFLD + I
            IVAR(IROCK) = IVARL(NC)
          ENDDO
          ENDDO
          ENDDO
        ELSE
          INDX = 7
          CHMSG = 'Unrecognized Indexing Option'
          IMSG = IJK
          CALL WRMSGS( INDX )
        ENDIF
        CLOSE(UNIT=26)
!
!---  Input-file input  ---
!
      ELSE
        ISTART = ISTX
        CALL RDINT(ISTART,ICOMMA,CHDUM,IVARX)
        DO 400 N = 1,NFLD
          IVAR(IZ(N)) = IVARX
 400    CONTINUE
      ENDIF
      IF( ALLOCATED(IVARL) ) THEN
      DEALLOCATE( IVARL,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IVARL'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
!
!---  Reset subroutine name string  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIJKI group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDCM( ISTART,CHDUM,UNTS,VAR )
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
!     Read dimensionless property values from an external file for
!     IJK, JKI, or KIJ indexed inputs with double indices
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 March 2002.
!     Last Modified by MD White, PNNL, 19 March 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,CDUM,FDUM,UNTS,UNTSX
      CHARACTER*132 VARBX
      CHARACTER*512 CHDUM
      REAL*8 VAR(*)
      LOGICAL FCHK,FBIN
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDCM'
      ISTX = ISTART
      VARBX = VARB
      IDFLT = 0
      VARB = 'External File Name'
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      FBIN = .FALSE.
!
!---  External-file input  ---
!
      IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
        ICOLON = INDEX(ADUM(1:),':') + 1
        FDUM = ADUM(ICOLON:NCHA)
        NCHF = NCHA-ICOLON+1
        INQUIRE( FILE=FDUM(1:NCHF), FORM=CDUM, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          INDX = 4
          CHMSG = 'IJK Indexing file does not exist: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'UNFORMATTED' .AND. (.NOT. FBIN) ) THEN
          INDX = 4
          CHMSG = 'Formatted IJK Indexing file is unformatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'FORMATTED' .AND. FBIN ) THEN
          INDX = 4
          CHMSG = 'Unformatted IJK Indexing file is formatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        END IF
        IF( FBIN ) THEN
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='UNFORMATTED')
        ELSE
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='FORMATTED')
        ENDIF
        WRITE(IWR,'(/,2A)') 'IJK Indexing File: ',FDUM(1:NCHF)
!
!---    Check for units  ---
!
        VARB = VARBX
        VARX = 1.D+0
        NCHU = INDEX( UNTS(1:),'  ' ) - 1
        IF( UNTS(1:NCHU).NE.'null' ) THEN
          IDFLT = 1
          UNTSX = UNTS
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTSX)
          WRITE(IWR,'(2X,3A)') VARB(1:IVR),', ',UNTSX(1:NCH)
          INDX = 0
          CALL RDUNIT(UNTSX,VARX,INDX)
        ENDIF
        IF( FBIN ) THEN
          READ(26) (VAR(IROCK),IROCK=1,NFLD)
        ELSE
          READ(26,*) (VAR(IROCK),IROCK=1,NFLD)
        ENDIF
        CLOSE(UNIT=26)
!
!---  Input-file input  ---
!
      ELSE
        ISTART = ISTX
        VARB = VARBX
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VARX)
        IDFLTDX = IDFLTD
!
!---    Check for units  ---
!
        NCHU = INDEX( UNTS(1:),'  ' ) - 1
        IF( UNTS(1:NCHU).NE.'null' ) THEN
          IDFLT = 1
          UNTSX = UNTS
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTSX)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTSX(1:NCH),
     &      ': ',VARX
          INDX = 0
          CALL RDUNIT(UNTSX,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,3A)') ' (',VARX,', ',UNTS(1:NCHU),')'
        ELSE
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VARX
        ENDIF
        IF( IDFLTDX.EQ.0 ) THEN
          DO 400 N = 1,NFLD
            VAR(N) = VARX
  400     CONTINUE
        ENDIF
      ENDIF
!
!---  Reset subroutine name string  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDCM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIJKT( ISTART,IJK,CHDUM,UNTS,TBL,ITBL,NL,NT,ILOG )
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
!     Read dimensionless property values from an external file for
!     IJK, JKI, or KIJ indexed inputs for tabular input.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 July 2002.
!     Last Modified by MD White, PNNL, 29 July 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,CDUM,FDUM,UNTS,UNTSX
      CHARACTER*132 VARBX
      CHARACTER*512 CHDUM
      REAL*8 VAR(LTBL),TBL(*)
      INTEGER ITBL(2,*)
      LOGICAL FCHK,FBIN
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIJKT'
      ISTX = ISTART
      VARBX = VARB
      VARB = 'External File Name'
      IDFLT = 0
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      FBIN = .FALSE.
      IF( INDEX(ADUM(1:),'binary').NE.0 .OR.
     &  INDEX(ADUM(1:),'bfile').NE.0 .OR.
     &  INDEX(ADUM(1:),'b_file').NE.0 ) FBIN = .TRUE.
!
!---  External-file input  ---
!
      IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
        ICOLON = INDEX(ADUM(1:),':') + 1
        FDUM = ADUM(ICOLON:NCHA)
        NCHF = NCHA-ICOLON+1
        INQUIRE( FILE=FDUM(1:NCHF), FORM=CDUM, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          INDX = 4
          CHMSG = 'IJK Indexing file does not exist: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'UNFORMATTED' .AND. (.NOT. FBIN) ) THEN
          INDX = 4
          CHMSG = 'Formatted IJK Indexing file is unformatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        ELSEIF( CDUM.EQ.'FORMATTED' .AND. FBIN ) THEN
          INDX = 4
          CHMSG = 'Unformatted IJK Indexing file is formatted: '
     &      // FDUM(1:NCHF)
          CALL WRMSGS( INDX )
        END IF
        IF( FBIN ) THEN
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='UNFORMATTED')
        ELSE
          OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &      FORM='FORMATTED')
        ENDIF
        WRITE(IWR,'(/,2A)') 'IJK Indexing File: ',FDUM(1:NCHF)
!
!---    Check for units  ---
!
        VARB = VARBX
        VARX = 1.D+0
        NCHU = INDEX( UNTS(1:),'  ' ) - 1
        IF( UNTS(1:NCHU).NE.'null' ) THEN
          IDFLT = 1
          UNTSX = UNTS
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTSX)
          WRITE(IWR,'(2X,3A)') VARB(1:IVR),', ',UNTSX(1:NCH)
          INDX = 0
          CALL RDUNIT(UNTSX,VARX,INDX)
        ENDIF
        NTBV = NFLD*NL
        IF( FBIN ) THEN
          READ(26) (VAR(N),N=1,NTBV)
        ELSE
          READ(26,*) (VAR(N),N=1,NTBV)
        ENDIF
        NV = 0
!
!---    IJK indexing  ---
!
        IF( IJK.EQ.1 ) THEN
          DO K = 1,KFLD
          DO J = 1,JFLD
          DO I = 1,IFLD
            N = (K-1)*IJFLD + (J-1)*IFLD + I
            ITBL(1,N) = NT + 1
            DO 90 L = 1,NL
              NT = NT + 1
              NV = NV + 1
              IF( ILOG.EQ.1 ) THEN
                TBL(NT) = LOG( EXP(VAR(NV))*VARX )
              ELSE
                TBL(NT) = VAR(NV)*VARX
              ENDIF
   90       CONTINUE
            ITBL(2,N) = NT
          ENDDO
          ENDDO
          ENDDO
!
!---    JKI indexing  ---
!
        ELSEIF( IJK.EQ.2 ) THEN
          DO I = 1,IFLD
          DO K = 1,KFLD
          DO J = 1,JFLD
            N = (K-1)*IJFLD + (J-1)*IFLD + I
            ITBL(1,N) = NT + 1
            DO 190 L = 1,NL
              NT = NT + 1
              NV = NV + 1
              IF( ILOG.EQ.1 ) THEN
                TBL(NT) = LOG( EXP(VAR(NV))*VARX )
              ELSE
                TBL(NT) = VAR(NV)*VARX
              ENDIF
  190       CONTINUE
            ITBL(2,N) = NT
          ENDDO
          ENDDO
          ENDDO
!
!---    KIJ indexing  ---
!
        ELSEIF( IJK.EQ.3 ) THEN
          DO J = 1,JFLD
          DO I = 1,IFLD
          DO K = 1,KFLD
            N = (K-1)*IJFLD + (J-1)*IFLD + I
            ITBL(1,N) = NT + 1
            DO 290 L = 1,NL
              NT = NT + 1
              NV = NV + 1
              IF( ILOG.EQ.1 ) THEN
                TBL(NT) = LOG( EXP(VAR(NV))*VARX )
              ELSE
                TBL(NT) = VAR(NV)*VARX
              ENDIF
  290       CONTINUE
            ITBL(2,N) = NT
          ENDDO
          ENDDO
          ENDDO
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Indexing Option' // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        CLOSE(UNIT=26)
!
!---  Input-file input  ---
!
      ELSE
        INDX = 4
        CHMSG = 'Input-file input not an option for tabular input'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine name string  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIJKT group  ---
!
      RETURN
      END


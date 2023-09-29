!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FRAC_RDROCK
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
!     Read input file for rock/soil zonation information for fractures.
!     Overlapping rock/soil zone definitions are sequence dependent.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Ramesh Sarathi, PNNL, 6 March 2020
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
      USE GEOM_FRC
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM
      CHARACTER*512 CHDUM
      CHARACTER*17 FORM2
      CHARACTER*4 FORM1
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM1,FORM2
      DATA FORM2 /'Soil 000000000'/
      DATA FORM1 /'(I )'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FRAC_RDROCK'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Fracture/Fault Zonation Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read the number of rock/soil zonation information lines  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ! NROCK = 0  !rss_add do not reset since we are appending to the list
      ! NROCK2 = 0
      VARB = 'Input Option [File, Zonation_File, Integer, Indexing]'
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)

!
!---  Read rock/soil zonation information from the input file  ---
!
      ISTART = 1
      VARB = 'Number of Fracture/Fault Zonation Lines'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO NL = 1, NLIN
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        ADUM(1:) = ' '
        VARB = 'Rock/Soil Name: '
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        DO 210 M = 1,NROCK
          IF( ROCK(M).EQ.ADUM ) THEN
            IROCK = M
            GOTO 220
          ENDIF
  210   CONTINUE
        NROCK = NROCK+1
        IF( NROCK.GT.LRC ) THEN
          INDX = 5
          CHMSG = 'Number of Rock/Soil Types > Parameter LRC'
          CALL WRMSGS( INDX )
        ENDIF
        ROCK(NROCK) = ADUM
        IROCK = NROCK
  220   CONTINUE
!
!---    Read rock/soil domain  ---
!
        VARB = 'Rock/Soil Domain Index: '
        CALL RDINT(ISTART,ICOMMA,CHDUM,NFX)  ! fracture/fault number
        CALL RDINT(ISTART,ICOMMA,CHDUM,I1)   ! local triangle starting index
        CALL RDINT(ISTART,ICOMMA,CHDUM,I2)   ! local triangle ending index
        I1 = MAX( 1,I1 )
        I1 = MIN( NTP_FRC(NFX),I1,I2 )       ! NTP_FRC is the number of triangles in the fracture
        I2 = MAX( 1,I1,I2 )
        I2 = MIN( NTP_FRC(NFX),I2 )
        WRITE(IWR,'(/,2A)') 'Rock/Soil Name: ',ROCK(IROCK)
        WRITE(IWR,'(2X,A)') 'Fracture Rock/Soil Domain:'
        WRITE(IWR,'(4X,A,I6)') 'NFX = ',NFX 
        WRITE(IWR,'(4X,A,I6,A,I6)') 'I = ',I1,' to ',I2
!
!---    Read scaling group associations  ---
!
        IF( ISLC(19).EQ.1 ) THEN
          ADUM(1:) = ' '
          VARB = 'Scaling Group Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
          DO 230 M = 1,NSCALE
            IF( SCALNM(M).EQ.ADUM ) THEN
              IF( ISCALE(IROCK).NE.0 .AND. ISCALE(IROCK).NE.M ) THEN
                INDX = 4
                CHMSG = 'Rock/Soil Name Associated with ' //
     &            'Multiple Scaling Groups: ' // ROCK(IROCK)
                CALL WRMSGS( INDX )
              ENDIF
              ISCALE(IROCK) = M
              GOTO 240
            ENDIF
  230     CONTINUE
          NSCALE = NSCALE+1
          IF( NSCALE.GT.LRC ) THEN
            INDX = 5
            CHMSG = 'Number of Scaling Groups > Parameter LRC'
            CALL WRMSGS( INDX )
          ENDIF
          SCALNM(NSCALE) = ADUM
          ISCALE(IROCK) = NSCALE
  240     CONTINUE
        ENDIF
!
!---    Label triangles with rock/soil types  ---
!
        DO I = I1,I2
          N = I + IP_FRC(1,NFX) - 1 !rss_ask - IP_FRC is the global triangle index offset
          IZ_FRC(N) = IROCK
        END DO
        IF( NL.LT.NLIN ) WRITE(IWR,'(/)')
      ENDDO
!
!---  Check for fault triangles without rock/soil types  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
          IF( IZ_FRC(NTX).LE.0 .OR. IZ_FRC(NTX).GT.NROCK ) THEN
            INDX = 7
            CHMSG = 'Undefined Fault Triangle Rock/Soil Type: ' //
     &        '@ Fault Triangle'
            IMSG = NTX
            CALL WRMSGS( INDX )
          ENDIF
        ENDDO
      ENDDO
!
!---  Check for conflicts between scaling group and rock/soil names  ---
!
      DO 402 M = 1,NROCK
        DO 400 L = 1,NSCALE
          IF( SCALNM(L).EQ.ROCK(M) )THEN
            INDX = 4
            CHMSG = 'Scaling Group Name Equals Rock/Soil Name: ' //
     &        SCALNM(L)
            CALL WRMSGS( INDX )
          ENDIF
  400   CONTINUE
  402 CONTINUE
  500 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FRAC_RDROCK group  ---
!
      RETURN
      END

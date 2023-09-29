!----------------------Program-----------------------------------------!
!
      PROGRAM SURFACE_BIN_TO_TXT
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
!----------------------Description-------------------------------------!
!
!     STOMPX: Subsurface Transport Over Multiple Phases Extensible
!
!     This utility routine converts surface.bin binary 
!     formatted files to a text surface file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 4 March 2022.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SF
      CHARACTER(64) :: CDUM,SF_BIN_FILE,SF_TXT_FILE
      CHARACTER(64) :: UNTM
      CHARACTER(64), DIMENSION(:,:), ALLOCATABLE :: CHSF,UNSF
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: GCNM
      CHARACTER(4) FORM1
      CHARACTER(5) FORM2
      CHARACTER(11) FORM3
      CHARACTER(9) FORM4
      CHARACTER(14) FORM5
      CHARACTER(6) FORM6
      CHARACTER(7) FORM7
      CHARACTER(9) FORM8
      CHARACTER*3 CHX
      LOGICAL FCHK
      INTEGER :: NBYTR = 8
      INTEGER :: NBYTI = 4
      INTEGER :: NBYTC = 1
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFT,ISFGC
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
      DATA FORM2 /'(A,$)'/
      DATA FORM3 /'(5X,A,2X,$)'/
      DATA FORM4 /'(2X,3A,$)'/
      DATA FORM5 /'(1PE13.6,A,$)'/
      DATA FORM6 /'(5X,A)'/
      DATA FORM7 /'(2X,3A)'/
      DATA FORM8 /'(1PE13.6)'/
!
!----------------------Executable Lines--------------------------------!
!
      LNGC = 1
      NCLA = IARGC()
      IHSF = 0
      IF( NCLA.EQ.0 ) THEN
        WRITE(6,'(/,A,/)')'Enter the binary surface filename: '
        READ(5,*) SF_BIN_FILE
        CALL LCASE( SF_BIN_FILE )
        CDUM = SF_BIN_FILE
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts a binary formatted ' // 
     &      'surface file created from executing STOMPX '
          WRITE(6,'(A)') '(e.g., surface.bin) to a text formatted ' //
     &      'surface file (e.g., surface). The script takes two '
          WRITE(6,'(A)') 'arguments, the name of the binary surface ' //
     &      'file and the name of the text surface file. If the '
          WRITE(6,'(A)') 'arguments are missing, the user will be ' // 
     &      'prompted for the filenames.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'surfacebinTo_txt.exe surface.bin surface'
          WRITE(6,'(A)') 'surfacebinTo_txt.exe surface_1.bin ' // 
     &      'surface_1.txt'
          STOP
        ENDIF
!
!---    No surface file, exit program with error message  ---
!
        SF_BIN_FILE = ADJUSTL(SF_BIN_FILE)
        INQUIRE( FILE=TRIM(SF_BIN_FILE), EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          WRITE(6,'(3A)') 'ERROR: Nonexistent "',TRIM(SF_BIN_FILE),
     &      '" file.'
          STOP
        ELSE
          OPEN(UNIT=55, FILE=TRIM(SF_BIN_FILE), STATUS='OLD',
     &      FORM='UNFORMATTED', ACCESS='STREAM')
        ENDIF
        WRITE(6,'(/,A,/)')'Enter the text "csv" surface filename: '
        READ(5,*) SF_TXT_FILE
        CALL LCASE( SF_TXT_FILE )
!  
!---    Open surface.csv file
!
        SF_TXT_FILE = ADJUSTL(SF_TXT_FILE)
        WRITE(6,'(2A)') 'Creating text surface file: ',TRIM(SF_TXT_FILE)
        OPEN(UNIT=56, FILE=TRIM(SF_TXT_FILE))
        CLOSE(UNIT=56,STATUS='DELETE')
        OPEN(UNIT=56, FILE=TRIM(SF_TXT_FILE), STATUS='NEW', 
     &    FORM='FORMATTED')
      ELSEIF( NCLA.EQ.1 .OR. NCLA.GT.2 ) THEN
        CALL GETARG( NA,CDUM )
        CALL LCASE( CDUM )
        IF( INDEX(CDUM(1:),'help').EQ.0 ) THEN
          WRITE(6,'(A)') 'ERROR: Incorrect Number of Arguments'
        ENDIF
        WRITE(6,'(/,A)') 'This script converts a binary formatted ' // 
     &    'surface file created from executing STOMPX '
        WRITE(6,'(A)') '(e.g., surface.bin) to a text formatted ' //
     &    'surface file (e.g., surface). The script takes two ' 
        WRITE(6,'(A)') 'arguments, the name of the binary surface ' //
     &    'file and the name of the text surface file. If the '
        WRITE(6,'(A)') 'arguments are missing, the user will be ' // 
     &    'prompted for the filenames.'
        WRITE(6,'(A,/)') 'Example command line entries'
        WRITE(6,'(A)') 'surfacebinTo_txt.exe surface.bin surface'
        WRITE(6,'(A)') 'surfacebinTo_txt.exe surface_1.bin ' // 
     &    'surface_1.txt'
        STOP
      ELSEIF( NCLA.EQ.2 ) THEN
        CALL GETARG( 1,SF_BIN_FILE )
        CALL LCASE( SF_BIN_FILE )
        SF_BIN_FILE = ADJUSTL(SF_BIN_FILE)
        INQUIRE( FILE=TRIM(SF_BIN_FILE), EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          WRITE(6,'(3A)') 'ERROR: Nonexistent "',TRIM(SF_BIN_FILE),
     &      '" file.'
          STOP
        ELSE
          OPEN(UNIT=55, FILE=TRIM(SF_BIN_FILE), STATUS='OLD',
     &      FORM='UNFORMATTED', ACCESS='STREAM')
        ENDIF
        CALL GETARG( 2,SF_TXT_FILE )
        CALL LCASE( SF_TXT_FILE )
        SF_TXT_FILE = ADJUSTL(SF_TXT_FILE)
        WRITE(6,'(2A)') 'Creating text surface file: ',TRIM(SF_TXT_FILE)
        OPEN(UNIT=56, FILE=TRIM(SF_TXT_FILE))
        CLOSE(UNIT=56,STATUS='DELETE')
        OPEN(UNIT=56, FILE=TRIM(SF_TXT_FILE), STATUS='NEW', 
     &    FORM='FORMATTED')
      ENDIF
!
!---  Write banner to screen  ---
!
      WRITE(6,'(A,//)')' Welcome to ...'
      WRITE(6,'(A,/)') '                       SURFACE_BIN_TO_TXT'
      WRITE(6,'(A)')   ' A surface flux file converter routine for ' //
     &  'the STOMPX simulator.  This'
      WRITE(6,'(A)')   ' routine, developed by the Pacific ' //
     &  'Northwest National Laboratory, reads'
      WRITE(6,'(A)')   ' a STOMPX surface.bin flux ' //
     &  'file and converts the data in the surface '
      WRITE(6,'(A,/)')   ' section to a text formatted file ' //
     &  'defined by the user (e.g., surface.dat).'
      WRITE(6,'(A)')   ' For support:  Tel: 509.372.6070'
      WRITE(6,'(A,/)') '               E-mail:  mark.white@pnnl.gov'
!
!---  Initialize binary file position  ---
!
      IPOS = 1
!
!---  Read operational mode index  ---
!
      READ(55, POS=IPOS, IOSTAT=ICHK) IOM
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'IOM: End of file encountered during ' //
     &    'binary surface file read: ',TRIM(SF_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'IOM: Read error encountered during ' //
     &    'binary surface file read: ',TRIM(SF_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NBYTI
!
!---  Operational mode dependent component names  ---
!
      IF( IOM.EQ.32 .OR. IOM.EQ.33 .OR. IOM.EQ.37 .OR. IOM.EQ.38
     &  .OR. IOM.EQ.39 ) THEN
        CHX = 'CO2'
      ELSE
        CHX = 'Air'
      ENDIF
!
!---  Read number of surfaces  ---
!
      READ(55, POS=IPOS, IOSTAT=ICHK) ISFGP
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'ISFGP: End of file encountered during ' //
     &    'binary surface file read: ',TRIM(SF_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'ISFGP: Read error encountered during ' //
     &    'binary surface file read: ',TRIM(SF_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NBYTI
!  
!---  Read time units
!
      READ(55, POS=IPOS, IOSTAT=ICHK) UNTM
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'UNTM: End of file encountered during ' //
     &    'binary surface file read: ',TRIM(SF_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'UNTM: Read error encountered during ' //
     &    'binary surface file read: ',TRIM(SF_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + 64*NBYTC
!
!---  Write header  ---
!
      WRITE(56,'(A,//)')' Welcome to ...'
      WRITE(56,'(A)')   '                           STOMP'
      WRITE(56,'(A,//)')'        Subsurface Transport Over '
     &   // 'Multiple Phases'
      WRITE(56,'(A)')   ' This file was produced by STOMP, '
     &   // 'a numerical simulator'
      WRITE(56,'(A)')   ' developed by the Pacific Northwest '
     &   // 'Laboratory, with'
      WRITE(56,'(A)')   ' support from the VOC-Arid '
     &   // 'Integrated Demonstration Project,'
      WRITE(56,'(A)')   ' Office of Technology Development, '
     &   // 'U.S. Department of Energy.'
      WRITE(56,'(A)')   ' Results from this version of STOMP '
     &   // 'should not be used for'
      WRITE(56,'(A,/)') ' license related applications.'
      WRITE(56,'(A,/)') ' For inquiries or assistance:  '
     &   // 'Call (509) 372-6070'
      WRITE(56,'(A,//)')'                      ---  SURFACE '
     &  // ' ---'
      WRITE(56,'(2A)') 'Version: 1087'
      WRITE(56,'(2A)') 'Date: ','Date system call inactive.'
      WRITE(56,'(2A)') 'Time: ','Time system call inactive.'
      WRITE(56,'(A,I4)') 'Number of Surfaces: ',ISFGP
      WRITE(56,'(A,I6,/)') 'Surface Variables: ',ISFGP
!  
!---  Allocate memory
!
      ALLOCATE( ISFT(1:ISFGP),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: ISFT'
        STOP
      ENDIF
      ALLOCATE( CHSF(1:2,1:ISFGP),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: CHSF'
        STOP
      ENDIF
      ALLOCATE( UNSF(1:2,1:ISFGP),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: CHSF'
        STOP
      ENDIF
      ALLOCATE( SF(1:2,1:ISFGP),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: SF'
        STOP
      ENDIF
      ALLOCATE( ISFGC(1:ISFGP),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: ISFGC'
        STOP
      ENDIF
      ALLOCATE( GCNM(1:LNGC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: GCNM'
        STOP
      ENDIF
!  
!---  Loop over surfaces
!
      DO NS = 1,ISFGP
!
!---    Read surface type  ---
!
        READ(55, POS=IPOS, IOSTAT=ICHK) ISFT(NS)
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'ISFGP: End of file encountered during ' //
     &      'binary surface file read: ',TRIM(SF_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'ISFGP: Read error encountered during ' //
     &      'binary surface file read: ',TRIM(SF_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read surface flux header character strings
!
        READ(55, POS=IPOS, IOSTAT=ICHK) (CHSF(M,NS),M=1,2)
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'CHSF: End of file encountered during ' //
     &      'binary surface file read: ',TRIM(SF_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'CHSF: Read error encountered during ' //
     &      'binary surface file read: ',TRIM(SF_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + 2*64*NBYTC
!  
!---    Read surface flux unit character strings
!
        READ(55, POS=IPOS, IOSTAT=ICHK) (UNSF(M,NS),M=1,2)
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'UNSF: End of file encountered during ' //
     &      'binary surface file read: ',TRIM(SF_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'UNSF: Read error encountered during ' //
     &      'binary surface file read: ',TRIM(SF_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + 2*64*NBYTC
      ENDDO
!  
!---  Write list of outputs  ---
!
      WRITE(56,'(3A)') 'Time,',TRIM(UNTM),','
      DO NS = 1,ISFGP
        IF( ISFT(NS).EQ.1 ) THEN
          WRITE(56,'(6A)') 'Heat Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.2 ) THEN
          WRITE(56,'(6A)') 'Aqueous Volumetric Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.3 ) THEN
          WRITE(56,'(6A)') 'Gas Volumetric Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.4 ) THEN
          WRITE(56,'(6A)') 'Nonaqueous-Liquid Volumetric Flux'
     &      ,',',TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.5 ) THEN
          WRITE(56,'(6A)') 'Aqueous Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.6 ) THEN
          WRITE(56,'(6A)') 'Gas Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.7 ) THEN
          WRITE(56,'(6A)') 'Nonaqueous-Liquid Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.8 ) THEN
          WRITE(56,'(6A)') 'Salt Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.9 ) THEN
          WRITE(56,'(6A)') 'Dissolved-Oil Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.10 ) THEN
          WRITE(56,'(6A)') 'Condensate Water Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.11 ) THEN
          IF( IOM.GE.36 .AND. IOM.LE.39 ) THEN
            WRITE(56,'(6A)') 'Gas CH4 Mass Flux',',',
     &        TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
          ELSEIF( IOM.EQ.35 ) THEN
            WRITE(56,'(6A)') 'Gas CO2 Mass Flux',',',
     &        TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
          ELSE
            WRITE(56,'(6A)') 'Gas Oil Mass Flux',',',
     &        TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
          ENDIF
        ELSEIF( ISFT(NS).EQ.12 ) THEN
          IF( IOM.GE.36 .AND. IOM.LE.39 ) THEN
            WRITE(56,'(6A)') 'Aqueous CH4 Mass Flux',',',
     &        TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
          ELSEIF( IOM.EQ.35 ) THEN
            WRITE(56,'(6A)') 'Aqueous CO2 Mass Flux',',',
     &        TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
          ELSE
            WRITE(56,'(6A)') 'Aqueous Oil Mass Flux',',',
     &        TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
          ENDIF
        ELSEIF( ISFT(NS).EQ.13 ) THEN
          IF( IOM.GE.36 .AND. IOM.LE.39 ) THEN
            WRITE(56,'(6A)') 'Total CH4 Mass Flux',',',
     &        TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
          ELSEIF( IOM.EQ.35 ) THEN
            WRITE(56,'(6A)') 'Total CO2 Mass Flux',',',
     &        TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
          ELSE
            WRITE(56,'(6A)') 'Total Oil Mass Flux',',',
     &        TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
          ENDIF
        ELSEIF( ISFT(NS).EQ.20 ) THEN
          WRITE(56,'(6A)') 'Gas-Advective Heat Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.21 ) THEN
          WRITE(56,'(6A)') 'Gas-Advective Water-Mass Flux',
     &      ',',TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.22 ) THEN
          WRITE(56,'(7A)') 'Gas-Advective ',CHX,'-Mass Flux',
     &      ',',TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.25 ) THEN
          WRITE(56,'(6A)') 'Gas-Diffusive Heat Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.26 ) THEN
          WRITE(56,'(6A)') 'Gas-Diffusive Water-Mass Flux',
     &      ',',TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.27 ) THEN
          WRITE(56,'(7A)') 'Gas-Diffusive ',CHX,'-Mass Flux',
     &      ',',TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.28 ) THEN
          WRITE(56,'(6A)') 'Gas CO2 Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.29 ) THEN
          WRITE(56,'(6A)') 'Aqueous CO2 Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.30 ) THEN
          WRITE(56,'(8A)') 'Total ',CHX,'-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.31 ) THEN
          WRITE(56,'(6A)') 'Gas-Advective Oil-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.32 ) THEN
          WRITE(56,'(6A)') 'Gas-Diffusive Oil-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.33 ) THEN
          WRITE(56,'(6A)') 'Gas-Total Oil-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.34 ) THEN
          WRITE(56,'(5A)') 'Surface Actual Evaporation,',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.35 ) THEN
          WRITE(56,'(5A)') 'Surface Potential Evaporation,',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.36 ) THEN
          WRITE(56,'(5A)') 'Surface Actual Transpiration,',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.37 ) THEN
          WRITE(56,'(5A)') 'Surface Potential Transpiration,',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.38 ) THEN
          WRITE(56,'(5A)') 'Surface Net Total Radiaion,',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.39 ) THEN
          WRITE(56,'(5A)') 'Surface Net Short-Wave Radiaion,',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.40 ) THEN
          WRITE(56,'(5A)') 'Surface Net Long-Wave Radiaion,',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.41 ) THEN
          WRITE(56,'(5A)') 'Surface Water-Mass Balance,',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.42 ) THEN
          WRITE(56,'(5A)') 'Surface Rain-Water Runoff,',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.43 ) THEN
          WRITE(56,'(6A)') 'Aqueous Water Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.44 ) THEN
          WRITE(56,'(6A)') 'Gas Water Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.45 ) THEN
          WRITE(56,'(6A)') 'Total Water Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.46 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Aqueous-Advective ' // 
     &      GCNM(IGC)(1:NCH0) // '-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.47 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Aqueous-Diffusive ' // 
     &      GCNM(IGC)(1:NCH0) // '-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.48 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Gas-Advective ' // 
     &      GCNM(IGC)(1:NCH0) // '-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.49 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Gas-Diffusive ' // 
     &      GCNM(IGC)(1:NCH0) // '-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.50 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Total-Advective ' // 
     &      GCNM(IGC)(1:NCH0) // '-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.51 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Total-Diffusive ' // 
     &      GCNM(IGC)(1:NCH0) // '-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.52 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Nonaqueous-Liquid-Advective ' // 
     &      GCNM(IGC)(1:NCH0) // '-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.53 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Nonaqueous-Liquid-Diffusive ' // 
     &      GCNM(IGC)(1:NCH0) // '-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.54 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Total-Advective-Diffusive ' // 
     &      GCNM(IGC)(1:NCH0) // '-Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.55 ) THEN
          WRITE(56,'(6A)') 'Liquid-CO2 CO2 Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.56 ) THEN
          WRITE(56,'(6A)') 'Liquid-CO2 CH4 Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.57 ) THEN
          WRITE(56,'(6A)') 'Liquid-CO2 H20 Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.58 ) THEN
          WRITE(56,'(6A)') 'Total N2 Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.59 ) THEN
          WRITE(56,'(6A)') 'Aqueous N2 Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.60 ) THEN
          WRITE(56,'(6A)') 'Gas N2 Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.61 ) THEN
          WRITE(56,'(6A)') 'Nonaqueous-Liquid N2 Mass Flux',
     &      ',',TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.62 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Gas Component ' // 
     &      GCNM(IGC)(1:NCH0) // ' Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.63 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Aqueous Component ' // 
     &      GCNM(IGC)(1:NCH0) // ' Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.64 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Nonaqueous-Liquid Component ' // 
     &      GCNM(IGC)(1:NCH0) // ' Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.65 ) THEN
          IGC = ISFGC(NS)
          NCH0 = INDEX(GCNM(IGC),'  ')-1
          WRITE(56,'(6A)') 'Total Component ' // 
     &      GCNM(IGC)(1:NCH0) // ' Mass Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.66 ) THEN
          WRITE(56,'(6A)') 'Aqueous Volumetric ' // 
     &      'Fracture-Matrix Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ELSEIF( ISFT(NS).EQ.67 ) THEN
          WRITE(56,'(6A)') 'Aqueous Volumetric ' // 
     &      'Borehole Flux',',',
     &      TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
!        ELSEIF( ISFT(NS).GT.100 .AND. ISFT(NS).LE.(100+NSOLU) ) THEN
!          NSL = ISFT(NS)-100
!          NCH = INDEX(SOLUT(NSL),'  ')-1
!          WRITE(56,'(7A)') 'Solute Flux (',SOLUT(NSL)(1:NCH),
!     &      '),',TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
!        ELSEIF( ISFT(NS).GT.(100+NSOLU) .AND. 
!     &    ISFT(NS).LE.(100+NSOLU+NEQC) ) THEN
!          NSL = ISFT(NS)-100
!          NCH = INDEX(SOLUT(NSL),'  ')-1
!          WRITE(56,'(7A)') 'Conservation Component ' // 
!     &      'Species Flux (',SOLUT(NSL)(1:NCH),
!     &      '),',TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
!        ELSEIF( ISFT(NS).GT.(100+NSOLU+NEQC) .AND. 
!     &    ISFT(NS).LE.(100+NSOLU+NEQC+NEQK) ) THEN
!          NSL = ISFT(NS)-100
!          NCH = INDEX(SOLUT(NSL),'  ')-1
!          WRITE(56,'(7A)') 'Kinetic Component ' // 
!     &      'Species Flux (',SOLUT(NSL)(1:NCH),
!     &      '),',TRIM(UNSF(1,NS)),',',TRIM(UNSF(2,NS)),','
        ENDIF
      ENDDO
      WRITE(56,'(/)')
!
!---  Loop over time steps in simulation
!
      DO
!
!---    Write header every ten time steps  ---
!
        IF( IHSF.EQ.0 .OR. IHSF.EQ.10 ) THEN
          IHSF = 0
          WRITE(56,FORM2) '   Time     '
          DO NS = 1,ISFGP
            IF( NS.EQ.ISFGP ) THEN
              WRITE(56,FORM3) TRIM(CHSF(1,NS))
              WRITE(56,FORM6) TRIM(CHSF(2,NS))
            ELSE
              WRITE(56,FORM3) TRIM(CHSF(1,NS))
              WRITE(56,FORM3) TRIM(CHSF(2,NS))
            ENDIF
          ENDDO
          WRITE(56,FORM4) ' [',UNTM(1:8),'] '
          DO NS = 1,ISFGP
            IF( NS.EQ.ISFGP ) THEN
              WRITE(56,FORM4) ' [',UNSF(1,NS)(1:8),'] '
              WRITE(56,FORM7) ' [',UNSF(2,NS)(1:8),'] '
            ELSE
              WRITE(56,FORM4) ' [',UNSF(1,NS)(1:8),'] '
              WRITE(56,FORM4) ' [',UNSF(2,NS)(1:8),'] '
            ENDIF
          ENDDO
        ENDIF
!  
!---    Read time  ---
!
        READ(55, POS=IPOS, IOSTAT=ICHK) TM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'TM: End of file encountered during ' //
     &      'binary surface file: ',TRIM(SF_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'TM: Read error encountered during ' //
     &      'binary surface file: ',TRIM(SF_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTR
!  
!---    End of file  ---
!
        IF( TM.LT.-8.999D+20 ) EXIT 
!  
!---    Write time  ---
!
        WRITE(56,FORM5) TM,' '
!  
!---    Loop over surfaces
!
        DO NS = 1,ISFGP
          DO M = 1,2
            READ(55, POS=IPOS, IOSTAT=ICHK) VARX
            SF(M,NS) = VARX
            IF( ICHK.GT.0 ) THEN
              WRITE(6,'(2A)') 'SF: End of file encountered during ' //
     &          'binary surface file: ',TRIM(SF_BIN_FILE)
              STOP
            ELSEIF( ICHK.LT.0 ) THEN
              WRITE(6,'(2A)') 'SF: Read error encountered during ' //
     &          'binary surface file: ',TRIM(SF_BIN_FILE)
              STOP
            ENDIF
            IPOS = IPOS + NBYTR
!  
!---        Surface flux rate and integral  ---
!
            IF( NS.EQ.ISFGP .AND. M.EQ.2 ) THEN
              WRITE(56,FORM8) SF(M,NS)
            ELSE
              WRITE(56,FORM5) SF(M,NS),' '
            ENDIF
          ENDDO
        ENDDO
        IHSF = IHSF + 1
      ENDDO
!  
!---  Deallocate memory
!
      DEALLOCATE( ISFT,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Deallocation Error: ISFT'
        STOP
      ENDIF
      DEALLOCATE( CHSF,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Deallocation Error: CHSF'
        STOP
      ENDIF
      DEALLOCATE( UNSF,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Deallocation Error: CHSF'
        STOP
      ENDIF
      DEALLOCATE( SF,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Deallocation Error: SF'
        STOP
      ENDIF
!
!---  Close the binary and text surface files  ---
!
      CLOSE( UNIT=55 )
      CLOSE( UNIT=56 )
      STOP
!
!---  End of SURFACE_BIN_TO_TXT program  ---
!
      END

!----------------------Function----------------------------------------!
!
      FUNCTION ICOUNT( I )
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
!     Count the number of digits in an integer variable.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNL, March, 1993.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IC = ABS(I)
      ICOUNT = 0
   10 CONTINUE
      ICOUNT = ICOUNT + 1
      IC = IC/10
      IF( IC.GT.0 ) GOTO 10
      IF( I.LT.0 ) ICOUNT = ICOUNT + 1
!
!---  End of ICOUNT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LCASE( CHDUM )
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
!     Convert all upper-case characters in a variable-length string
!     variable to lower case.  This subroutine does not disturb
!     non-alphabetic characters; only captial letters
!     (ASCII 65 through 90) are modified.
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, Battelle, March, 1991.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*(*) CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      DO N = 1,LEN(CHDUM)
        M = ICHAR(CHDUM(N:N))
        IF( M.GE.65 .AND. M.LE.90 ) THEN
          M = M + 32
          CHDUM(N:N) = CHAR(M)
        ENDIF
      ENDDO
!
!---  End of LCASE group
!
      RETURN
      END

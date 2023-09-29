!----------------------Program-----------------------------------------!
!
      PROGRAM SURFACE_BIN_TO_CSV
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
!     formatted files to csv format for plotting.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 March 2022.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SF
      CHARACTER(64) :: CDUM,SF_BIN_FILE,OUTFILE
      CHARACTER(64) :: UNTM
      CHARACTER(64), DIMENSION(:,:), ALLOCATABLE :: CHSF,UNSF
      CHARACTER(4) FORM1
      LOGICAL FCHK
      INTEGER :: NBYTR = 8
      INTEGER :: NBYTI = 4
      INTEGER :: NBYTC = 1
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFT
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
      NCLA = IARGC()
      IF( NCLA.EQ.0 ) THEN
        WRITE(6,'(/,A,/)')'Enter the binary surface filename: '
        READ(5,*) SF_BIN_FILE
        CALL LCASE( SF_BIN_FILE )
        CDUM = SF_BIN_FILE
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts binary formatted ' // 
     &      'surface.bin files created from executing STOMPX to a '
          WRITE(6,'(A)') '"csv" (comma separated variables) ' //
     &      'formatted file. The script takes two arguments, the'
          WRITE(6,'(A)') 'name of the binary surface file and. ' //
     &      'the name of the text "csv" file, which should have a '
          WRITE(6,'(A)') '"csv" extension (e.g., surface.csv). ' //
     &      'If the arguments are missing, the user will be prompted  '
          WRITE(6,'(A)') 'for the filenames. The script determines ' //
     &      'the variables in the surface file, and prompts the user '
          WRITE(6,'(A)') 'for a list of variables to be included ' //
     &      'in the surface.csv file. The variable order will be that '
          WRITE(6,'(A,/)') 'specified by the user, which allows the ' //
     &      'variables to be reordered from the order in the surface '
          WRITE(6,'(A,/)') 'file. Generated "csv" files will ' //
     &      'contain one header line, and can be read by a number of '
          WRITE(6,'(A,/)') 'plotting packages for making plots of ' //
     &      'variables versus time. '
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'surfacebinTo_csv.exe surface.bin surface.csv'
          WRITE(6,'(A)') 'surfacebinTo_csv.exe surface_1.bin ' // 
     &      'surface_1.csv'
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
        READ(5,*) OUTFILE
        CALL LCASE( OUTFILE )
!  
!---    Open surface.csv file
!
        OUTFILE = ADJUSTL(OUTFILE)
        WRITE(6,'(2A)') 'Creating surface "csv" file: ',TRIM(OUTFILE)
        OPEN(UNIT=56, FILE=TRIM(OUTFILE))
        CLOSE(UNIT=56,STATUS='DELETE')
        OPEN(UNIT=56, FILE=TRIM(OUTFILE), STATUS='NEW', 
     &    FORM='FORMATTED')
      ELSEIF( NCLA.EQ.1 .OR. NCLA.GT.2 ) THEN
        CALL GETARG( NA,CDUM )
        CALL LCASE( CDUM )
        IF( INDEX(CDUM(1:),'help').EQ.0 ) THEN
          WRITE(6,'(A)') 'ERROR: Incorrect Number of Arguments'
        ENDIF
        WRITE(6,'(/,A)') 'This script converts binary formatted ' // 
     &    'surface.bin files created from executing STOMPX to a '
        WRITE(6,'(A)') '"csv" (comma separated variables) ' //
     &    'formatted file. The script takes two arguments, the'
        WRITE(6,'(A)') 'name of the binary surface file and. ' //
     &    'the name of the text "csv" file, which should have a '
        WRITE(6,'(A)') '"csv" extension (e.g., surface.csv). ' //
     &    'If the arguments are missing, the user will be prompted  '
        WRITE(6,'(A)') 'for the filenames. The script determines ' //
     &    'the variables in the surface file, and prompts the user '
        WRITE(6,'(A)') 'for a list of variables to be included ' //
     &    'in the surface.csv file. The variable order will be that '
        WRITE(6,'(A,/)') 'specified by the user, which allows the ' //
     &    'variables to be reordered from the order in the surface '
        WRITE(6,'(A,/)') 'file. Generated "csv" files will ' //
     &    'contain one header line, and can be read by a number of '
        WRITE(6,'(A,/)') 'plotting packages for making plots of ' //
     &    'variables versus time. '
        WRITE(6,'(A,/)') 'Example command line entries'
        WRITE(6,'(A)') 'surfacebinTo_csv.exe surface.bin surface.csv'
        WRITE(6,'(A)') 'surfacebinTo_csv.exe surface_1.bin ' // 
     &    'surface_1.csv'
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
        CALL GETARG( 2,OUTFILE )
        CALL LCASE( OUTFILE )
        OUTFILE = ADJUSTL(OUTFILE)
        WRITE(6,'(2A)') 'Creating surface "csv" file: ',TRIM(OUTFILE)
        OPEN(UNIT=56, FILE=TRIM(OUTFILE))
        CLOSE(UNIT=56,STATUS='DELETE')
        OPEN(UNIT=56, FILE=TRIM(OUTFILE), STATUS='NEW', 
     &    FORM='FORMATTED')
      ENDIF
!
!---  Write banner to screen  ---
!
      WRITE(6,'(A,//)')' Welcome to ...'
      WRITE(6,'(A,/)') '                       SURFACE_BIN_TO_CSV'
      WRITE(6,'(A)')   ' A surface flux file converter routine for ' //
     &  'the STOMPX simulator.  This'
      WRITE(6,'(A)')   ' routine, developed by the Pacific ' //
     &  'Northwest National Laboratory, reads'
      WRITE(6,'(A)')   ' a STOMPX surface.bin flux ' //
     &  'file and converts the data in the surface '
      WRITE(6,'(A,/)')   ' section to a csv formatted file ' //
     &  'defined by the user (e.g., surface.csv).'
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
      WRITE(56,'(A,I6,/)') 'Surface Variables: ',NSF
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
!---  Loop over time steps in simulation
!
      DO
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
          ENDDO
        ENDDO
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
!---  Close the binary surface.bin file and text surface.csv file  ---
!
      CLOSE( UNIT=55 )
      CLOSE( UNIT=56 )
      STOP
!
!---  End of SURFACE_BIN_TO_CSV program  ---
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

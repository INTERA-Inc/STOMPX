!----------------------Program-----------------------------------------!
!
      PROGRAM RESTART_BIN_TO_TXT
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
!     This utility routine converts one or more restart_xxx.bin binary 
!     formatted files to text (ASCII) formatted restart.xxx files
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 16 February 2022.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: VAR
      CHARACTER(128) :: CDUM,RST_BIN_FILE,RST_TXT_FILE
      CHARACTER*20 FORM5
      CHARACTER*38 FORM6
      CHARACTER*39 FORM7
      CHARACTER*40 FORM14
      CHARACTER(64) :: CH_VRSN
      LOGICAL FCHK
      INTEGER :: NBYTR = 8
      INTEGER :: NBYTI = 4
      INTEGER :: NBYTC = 1
      INTEGER, DIMENSION(:), ALLOCATABLE :: NFCGC
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ND
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM5 /'(10(1PE22.15,1X),I3)'/
      DATA FORM6 /'(10(1PE22.15,1X),I3,1X,1(1PE22.15,1X))'/
      DATA FORM7 /'(10(1PE22.15,1X),I3,1X,10(1PE22.15,1X))'/
      DATA FORM14 /'(10(1PE22.15,1X),I3,1X,100(1PE22.15,1X))'/
!
!----------------------Executable Lines--------------------------------!
!
      NCLA = IARGC()
      IF( NCLA.EQ.0 ) THEN
        WRITE(6,'(/,A,/)')'Enter binary restart file name or "help": '
        READ(5,*) RST_BIN_FILE
        CALL LCASE( RST_BIN_FILE )
        CDUM = RST_BIN_FILE
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts binary formatted ' // 
     &      'restart_xxx.bin files created from executing STOMPX'
          WRITE(6,'(A)') 'to conventional text (ASCII) formatted' //
     &      ' restart.xxx files created from executing STOMP.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'restartbinTo.exe restart_0000007.bin'
          WRITE(6,'(A,/)') 'restartbinTo.exe restart_*.bin '
          STOP
        ENDIF
        WRITE(6,'(/,A,/)')'Enter a single binary restart filename ' // 
     &    '(e.g., restart_0000007.bin): '
        READ(5,*) RST_BIN_FILE
        NCLA = 1
      ELSE
        CALL GETARG( 1,RST_BIN_FILE )
        CALL LCASE( RST_BIN_FILE )
        CDUM = RST_BIN_FILE
        CALL LCASE( CDUM )
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts binary formatted ' // 
     &      'restart_xxx.bin files created from executing STOMPX'
          WRITE(6,'(A)') 'to conventional text (ASCII) formatted' //
     &      ' restart.xxx files created from executing STOMP.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'restartbinTo.exe restart_0000007.bin'
          WRITE(6,'(A,/)') 'restartbinTo.exe restart_*.bin '
          STOP
        ENDIF
      ENDIF
!
!---  Write banner to screen  ---
!       
      WRITE(6,'(A,//)')' Welcome to ...'
      WRITE(6,'(A,/)') '                       RESTART_BIN_TO_TXT'
      WRITE(6,'(A)')   ' A binary restart file converter ' //
     &  'routine for the STOMPX simulator.'
      WRITE(6,'(A)')   ' This routine, developed by the Pacific ' //
     &  'Northwest National Laboratory,'
      WRITE(6,'(A)')   ' reads one or more STOMPX binary formatted ' //
     &  ' restart_xxx.bin files to '
      WRITE(6,'(A,/)')   ' text (ASCII) formatted restart.xxx files.'
      WRITE(6,'(A)')   ' For support:  Tel: 509.372.6070'
      WRITE(6,'(A,/)') '               E-mail:  mark.white@pnnl.gov'
!
!---  Loop over STOMPX restart_xxx.bin files ---
!
      DO NA = 1,NCLA
        IF( NA.GE.2 ) CALL GETARG( NA,RST_BIN_FILE )
        RST_TXT_FILE = 'restart.'
        RST_TXT_FILE(9:15) = RST_BIN_FILE(9:15)
!
!---    No restart file, exit program with error message  ---
!
        INQUIRE( FILE=TRIM(RST_BIN_FILE), EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          WRITE(6,'(A)') 'ERROR: Nonexistent restart file: ',
     &      TRIM(RST_BIN_FILE)
          STOP
        ELSE
          OPEN(UNIT=55, FILE=TRIM(RST_BIN_FILE), STATUS='OLD', 
     &      FORM='UNFORMATTED', ACCESS='STREAM')
          OPEN(UNIT=56, FILE=TRIM(RST_TXT_FILE), STATUS='UNKNOWN', 
     &      FORM='FORMATTED')
          CLOSE(UNIT=56, STATUS='DELETE')
          OPEN(UNIT=56, FILE=TRIM(RST_TXT_FILE), STATUS='NEW', 
     &      FORM='FORMATTED')
        ENDIF
        WRITE(6,'(4A)') 'Converting binary restart file: ',
     &    TRIM(RST_BIN_FILE),' to text (ASCII) restart file: ',
     &    TRIM(RST_TXT_FILE)
!
!---    Initialize binary file position  ---
!
        IPOS = 1
!  
!---    Read maximum number of Newton-Raphson iterations
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NRIMX
        IF( ICHK.GT.0 ) THEN
          NCH = INDEX( RST_BIN_FILE(1:),'  ' ) - 1
          WRITE(6,'(2A)') 'NRIMX: End of file encountered during ' //
     &      'restart file read: ',RST_BIN_FILE(1:NCH)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NRIMX: Read error encountered during ' //
     &      'restart file read: ',RST_BIN_FILE(1:NCH)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read time step
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NSTEP
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NSTEP: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NSTEP: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of field grid cells
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of x-direction field grid cells
!
        READ(55, POS=IPOS, IOSTAT=ICHK) IFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'IFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'IFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of y-direction field grid cells
!
        READ(55, POS=IPOS, IOSTAT=ICHK) JFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'JFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'JFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of z-direction field grid cells
!
        READ(55, POS=IPOS, IOSTAT=ICHK) KFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'KFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'KFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of field processors
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NPFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of x-direction field processors
!
        READ(55, POS=IPOS, IOSTAT=ICHK) IPFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'IFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'IFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of y-direction field processors
!
        READ(55, POS=IPOS, IOSTAT=ICHK) JPFLD
        IF( ICHK.GT.0 ) THEN
          NCH = INDEX( RST_BIN_FILE(1:),'  ' ) - 1
          WRITE(6,'(2A)') 'JFLD: End of file encountered during ' //
     &      'plot file read: ',RST_BIN_FILE(1:NCH)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'JFLD: Read error encountered during plot ' //
     &      'file read: ',RST_BIN_FILE(1:NCH)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of z-direction field processors
!
        READ(55, POS=IPOS, IOSTAT=ICHK) KPFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'KFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'KFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!
!---    Plane node counts  ---
!
        IJFLD = IFLD*JFLD
        JKFLD = JFLD*KFLD
        KIFLD = KFLD*IFLD
!  
!---    Read number of processors
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NP
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NP: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NP: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Allocate memory for number of grid cells
!
        ALLOCATE( NFCGC(1:NP),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: NFCGC'
          STOP
        ENDIF
!  
!---    Loop over number of processors, reading number of nodes per
!       processor, including ghost cells
!
        NFCGC_G = 0
        DO MP = 1,NP
          READ(55, POS=IPOS, IOSTAT=ICHK) NFCGC(MP)
          NFCGC_G = NFCGC_G + NFCGC(MP)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'NFCGC: End of file encountered during ' //
     &        'restart file read: ',TRIM(RST_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'NFCGC: Read error encountered during ' //
     &        'restart file read: ',TRIM(RST_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NBYTI
        ENDDO
!  
!---    Read number of solutes
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NSOLU
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NSOLU: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NSOLU: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read operational mode index
!
        READ(55, POS=IPOS, IOSTAT=ICHK) IOM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'IOM: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'IOM: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read soil freezing solution control index
!
        READ(55, POS=IPOS, IOSTAT=ICHK) ISLC5X
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'ISLC(5): End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'ISLC(5): Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of reactive species
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NSPR
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NSPR: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NSPR: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of reactive solid species
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NSPS
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NSPS: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NSPS: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!
!---    Read time  ---
!
        READ(55, POS=IPOS, IOSTAT=ICHK) TM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'TM: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'TM: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTR
!  
!---    Read time step
!
        READ(55, POS=IPOS, IOSTAT=ICHK) DT
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'DT: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'DT: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTR
!  
!---    Read maximum time step
!
        READ(55, POS=IPOS, IOSTAT=ICHK) DTMX
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'DTMX: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'DTMX: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTR
!  
!---    Read minimum time step
!
        READ(55, POS=IPOS, IOSTAT=ICHK) DTMN
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'DTMN: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'DTMN: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTR
!  
!---    Read time step acceleration factor
!
        READ(55, POS=IPOS, IOSTAT=ICHK) DTAF
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'DTAF: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'DTAF: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTR
!  
!---    Read maximum residual
!
        READ(55, POS=IPOS, IOSTAT=ICHK) RSDMX
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'RSDMX: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'RSDMX: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTR
!  
!---    Read time step cut factor
!
        READ(55, POS=IPOS, IOSTAT=ICHK) DTCF
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'DTCF: End of file encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'DTCF: Read error encountered during ' //
     &      'restart file read: ',TRIM(RST_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTR
!
!---    Write header to text (ASCII) restart.xxx file  ---
!
        WRITE(56,'(A,//)')' Welcome to ...'
        WRITE(56,'(A)')   '                           STOMP'
        WRITE(56,'(A,//)')'        Subsurface Transport Over ' // 
     &    'Multiple Phases'
        WRITE(56,'(A)')   ' This file was produced by STOMP, a ' // 
     &    'numerical simulator'
        WRITE(56,'(A)')   ' developed by the Pacific Northwest ' // 
     &    'Laboratory, with'
        WRITE(56,'(A)')   ' support from the VOC-Arid Integrated ' // 
     &    'Demonstration Project,'
        WRITE(56,'(A)')   ' Office of Technology Development, U.S. ' // 
     &    'Department of Energy.'
        WRITE(56,'(A)')   ' Results from this version of STOMP ' // 
     &   'should not be used for'
        WRITE(56,'(A,/)') ' license related applications.'
        WRITE(56,'(A,/)') ' For inquiries or assistance:  ' // 
     &  'Call (509) 372-6070'
        WRITE(56,'(A,//)')'                       ---  RESTART  ---'
        CH_VRSN = '1087'
        WRITE(56,'(2A)') 'Version: ',CH_VRSN
        WRITE(56,'(A)') 'Date: '
        WRITE(56,'(A)') 'Time: '
!
!---    Write timing data, field data by node numbers  ---
!
        IF( NSPR.GT.0 ) THEN
          WRITE(56,'(/,7(1PE22.15),10(I9))') TM,DT,DTMX,DTMN,DTAF,
     &      RSDMX,DTCF,NRIMX,NSTEP,NFLD,NSOLU,IOM,ISLC5X,NSPR,NSPS,
     &      NFBN,NWN_LW
        ELSE
          WRITE(56,'(/,7(1PE22.15),8(I9))') TM,DT,DTMX,DTMN,DTAF,
     &      RSDMX,DTCF,NRIMX,NSTEP,NFLD,NSOLU,IOM,ISLC5X,NFBN,NWN_LW
        ENDIF
!  
!---    STOMPX-CO2 to STOMP-CO2
!
        IF( IOM.EQ.32 ) THEN
!  
!---      Allocate memory for restart variables
!
          NRSV = 14
          NVAR = NRSV+NSOLU+NSPR+NSPS
          ALLOCATE( VAR(1:NVAR),STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            WRITE(6,'(A)') 'Allocation Error: VAR'
            STOP
          ENDIF
!  
!---      Set output formatting
!
          IF( (NSOLU+NSPR+NSPS).LE.0 ) THEN
            WRITE( FORM5(2:3),'(I2)' ) NRSV-1
          ELSEIF( (NSOLU+NSPR+NSPS).LT.10 ) THEN
            WRITE( FORM6(2:3),'(I2)' ) NRSV-1
            WRITE( FORM6(24:24),'(I1)' ) NSOLU+NSPR+NSPS
          ELSEIF( (NSOLU+NSPR+NSPS).LT.100 ) THEN
            WRITE( FORM7(2:3),'(I2)' ) NRSV-1
            WRITE( FORM7(24:25),'(I2)' ) NSOLU+NSPR+NSPS
          ELSE
            WRITE( FORM14(2:3),'(I2)' ) NRSV-1
            WRITE( FORM14(23:25),'(I3)' ) NSOLU+NSPR+NSPS
          ENDIF
        ENDIF
!
!---    Define the grid pointers  ---
!
        ALLOCATE( ND(1:IFLD,1:JFLD,1:KFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: ND'
          STOP
        ENDIF
        DO K = 1,KFLD
          DO J = 1,JFLD
            DO I = 1,IFLD
              ND(I,J,K) = (K-1)*IJFLD + (J-1)*IFLD + I
            ENDDO
          ENDDO
        ENDDO
!
!---    Loop over nodes in i,j,k order  ---
!
        NC = 0
        IL = IFLD/IPFLD
        JL = JFLD/JPFLD
        KL = KFLD/KPFLD
        DO KP = 1,KPFLD
          IF( KP.EQ.1 ) THEN
            KD1X = 1
            KD2X = KL + 1
          ELSEIF( KP.EQ.KPFLD ) THEN
            KD1X = (KP-1)*KL
            KD2X = KFLD
          ELSE
            KD1X = (KP-1)*KL
            KD2X = KP*KL + 1
          ENDIF
          IF( KPFLD.EQ.1 ) THEN
            KD1X = 1
            KD2X = 1
          ENDIF
          DO K = KD1X,KD2X
            DO JP = 1,JPFLD
              IF( JP.EQ.1 ) THEN
                JD1X = 1
                JD2X = JL + 1
              ELSEIF( JP.EQ.JPFLD ) THEN
                JD1X = (JP-1)*JL
                JD2X = JFLD
              ELSE
                JD1X = (JP-1)*JL
                JD2X = JP*JL + 1
              ENDIF
              IF( JPFLD.EQ.1 ) THEN
                JD1X = 1
                JD2X = 1
              ENDIF
              DO J = JD1X,JD2X
                DO IP = 1,IPFLD
                  IF( IP.EQ.1 ) THEN
                    ID1X = 1
                    ID2X = IL + 1
                  ELSEIF( IP.EQ.IPFLD ) THEN
                    ID1X = (IP-1)*IL
                    ID2X = IFLD
                  ELSE
                    ID1X = (IP-1)*IL
                    ID2X = IP*IL + 1
                  ENDIF
                  IF( IPFLD.EQ.1 ) THEN
                    ID1X = 1
                    ID2X = 1
                  ENDIF
                  DO I = ID1X,ID2X
!
!---                Count field and ghost cells  ---
!
                    NC = NC + 1
                    IPOSX = IPOS + (NC-1)*NBYTR
!  
!---                Read restart variables
!
                    DO M = 1,NVAR
                      READ(55, POS=IPOSX, IOSTAT=ICHK) VAR(M)
                      IPOSX = IPOSX + NFCGC_G*NBYTR
                    ENDDO
!
!---                Bottom ghost cell  ---
!
                    IF( K.EQ.KD1X .AND. KP.GT.1 ) CYCLE
!
!---                South ghost cell  ---
!
                    IF( J.EQ.JD1X .AND. JP.GT.1 ) CYCLE
!
!---                West ghost cell  ---
!
                    IF( I.EQ.ID1X .AND. IP.GT.1 ) CYCLE
!
!---                East ghost cell  ---
!
                    IF( I.EQ.ID2X .AND. IP.LT.IPFLD ) CYCLE
!
!---                North ghost cell  ---
!
                    IF( J.EQ.JD2X .AND. JP.LT.JPFLD ) CYCLE
!
!---                Top ghost cell  ---
!
                    IF( K.EQ.KD2X .AND. KP.LT.KPFLD ) CYCLE
!
!---                0 solutes and zero reactive species  ---
!
                    N = ND(I,J,K)
                    IF( (NSOLU+NSPR+NSPS).LE.0 ) THEN
                      WRITE(56,FORM5) (VAR(M),M=1,NRSV-1),
     &                  INT(VAR(NRSV))
!
!---                <10 solutes and zero reactive species  ---
!
                    ELSEIF( (NSOLU+NSPR+NSPS).LT.10 ) THEN
                      WRITE(56,FORM6) (VAR(M),M=1,NRSV-1),
     &                  INT(VAR(NRSV)),
     &                  (VAR(M),M=NRSV+1,NRSV+NSOLU+NSPR+NSPS)
!
!---                <100 solutes and zero reactive species  ---
!
                    ELSEIF( (NSOLU+NSPR+NSPS).LT.100 ) THEN
                      WRITE(56,FORM7) (VAR(M),M=1,NRSV-1),
     &                  INT(VAR(NRSV)),
     &                  (VAR(M),M=NRSV+1,NRSV+NSOLU+NSPR+NSPS)
!
!---                >=100 solutes and zero reactive species  ---
!
                    ELSE
                      WRITE(56,FORM14) (VAR(M),M=1,NRSV-1),
     &                  INT(VAR(NRSV)),
     &                  (VAR(M),M=NRSV+1,NRSV+NSOLU+NSPR+NSPS)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!  
!---    Deallocate memory
!
        DEALLOCATE( NFCGC,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: NFCGC'
          STOP
        ENDIF
        DEALLOCATE( VAR,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: VAR'
          STOP
        ENDIF
        DEALLOCATE( ND,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: ND'
          STOP
        ENDIF
!
!---    Close the binary restart_xxx.bin file  ---
!
        CLOSE( UNIT=55 )
!
!---    Close the text (ASCII) restart.xxx file  ---
!
        CLOSE( UNIT=56 )
      ENDDO
      STOP
!
!---  End of RESTART_BIN_TO_TXT program  ---
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
!     lcase.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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

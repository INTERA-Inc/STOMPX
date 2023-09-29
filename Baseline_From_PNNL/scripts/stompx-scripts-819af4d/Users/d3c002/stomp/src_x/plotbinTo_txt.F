!----------------------Program-----------------------------------------!
!
      PROGRAM PLOT_BIN_TO_TXT
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
!     This utility routine converts one or more plot_xxx.bin binary 
!     formatted files to text (ASCII) formatted plot.xxx files
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 November 2021.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: VER,VERX
      REAL*8, DIMENSION(:), ALLOCATABLE :: VAR,VARX
      CHARACTER(128) :: CDUM,PLT_BIN_FILE,PLT_TXT_FILE
      CHARACTER(64) :: UNTM,UNLN,PLTVNM
      CHARACTER*4 FORM1
      CHARACTER*16 FORM2
      CHARACTER*11 FORM3
      LOGICAL FCHK
      INTEGER :: NBYTR = 8
      INTEGER :: NBYTI = 4
      INTEGER :: NBYTC = 1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IXP,IXPX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ID_MPI,JD_MPI,KD_MPI
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ND
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 / '(I )' /
      DATA FORM2 / '(10(1PE16.9,1X))' /
      DATA FORM3 / '(10(I2,1X))' /
!
!----------------------Executable Lines--------------------------------!
!
      NCLA = IARGC()
      IF( NCLA.EQ.0 ) THEN
        WRITE(6,'(/,A,/)')'Enter binary plot file name or "help": '
        READ(5,*) PLT_BIN_FILE
        CALL LCASE( PLT_BIN_FILE )
        CDUM = PLT_BIN_FILE
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts binary formatted ' // 
     &      'plot_xxx.bin files created from executing STOMPX'
          WRITE(6,'(A)') 'to conventional text (ASCII) formatted' //
     &      ' plot.xxx files created from executing STOMP.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'plotbinTo.exe plot_0000007.bin'
          WRITE(6,'(A,/)') 'plotbinTo.exe plot_*.bin '
          STOP
        ENDIF
        WRITE(6,'(/,A,/)')'Enter a single binary plot filename ' // 
     &    '(e.g., plot_0000007.bin): '
        READ(5,*) PLT_BIN_FILE
        NCLA = 1
      ELSE
        CALL GETARG( 1,PLT_BIN_FILE )
        CALL LCASE( PLT_BIN_FILE )
        CDUM = PLT_BIN_FILE
        CALL LCASE( CDUM )
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts binary formatted ' // 
     &      'plot_xxx.bin files created from executing STOMPX'
          WRITE(6,'(A)') 'to conventional text (ASCII) formatted' //
     &      ' plot.xxx files created from executing STOMP.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'plotbinTo.exe plot_0000007.bin'
          WRITE(6,'(A,/)') 'plotbinTo.exe plot_*.bin '
          STOP
        ENDIF
      ENDIF
!
!---  Write banner to screen  ---
!       
      WRITE(6,'(A,//)')' Welcome to ...'
      WRITE(6,'(A,/)') '                       PLOT_BIN_TO_TXT'
      WRITE(6,'(A)')   ' A binary plot file converter ' //
     &  'routine for the STOMPX simulator.'
      WRITE(6,'(A)')   ' This routine, developed by the Pacific ' //
     &  'Northwest National Laboratory,'
      WRITE(6,'(A)')   ' reads one or more STOMPX binary formatted ' //
     &  ' plot_xxx.bin files to '
      WRITE(6,'(A,/)')   ' text (ASCII) formatted plot.xxx files.'
      WRITE(6,'(A)')   ' For support:  Tel: 509.372.6070'
      WRITE(6,'(A,/)') '               E-mail:  mark.white@pnnl.gov'
!
!---  No plot_xgrid.bin file, exit program with error message  ---
!
      INQUIRE( FILE='plot_xgrid.bin', EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        WRITE(6,'(A)') 'ERROR: Nonexistent plot_xgrid.bin file'
        STOP
      ELSE
        OPEN(UNIT=51, FILE='plot_xgrid.bin', STATUS='OLD', 
     &    FORM='UNFORMATTED', ACCESS='STREAM')
      ENDIF
!
!---  No plot_ygrid.bin file, exit program with error message  ---
!
      INQUIRE( FILE='plot_ygrid.bin', EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        WRITE(6,'(A)') 'ERROR: Nonexistent plot_ygrid.bin file'
        STOP
      ELSE
        OPEN(UNIT=52, FILE='plot_ygrid.bin', STATUS='OLD', 
     &    FORM='UNFORMATTED', ACCESS='STREAM')
      ENDIF
!
!---  No plot_zgrid.bin file, exit program with error message  ---
!
      INQUIRE( FILE='plot_zgrid.bin', EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        WRITE(6,'(A)') 'ERROR: Nonexistent plot_zgrid.bin file'
        STOP
      ELSE
        OPEN(UNIT=53, FILE='plot_zgrid.bin', STATUS='OLD', 
     &    FORM='UNFORMATTED', ACCESS='STREAM')
      ENDIF
!
!---  Loop over STOMPX plot_xxx.bin files ---
!
      DO NA = 1,NCLA
        IF( NA.GE.2 ) CALL GETARG( NA,PLT_BIN_FILE )
        PLT_TXT_FILE = 'plot.'
        PLT_TXT_FILE(6:12) = PLT_BIN_FILE(6:12)
!
!---    No plot file, exit program with error message  ---
!
        INQUIRE( FILE=TRIM(PLT_BIN_FILE), EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          WRITE(6,'(A)') 'ERROR: Nonexistent plot file: ',
     &      TRIM(PLT_BIN_FILE)
          STOP
        ELSE
          OPEN(UNIT=55, FILE=TRIM(PLT_BIN_FILE), STATUS='OLD', 
     &      FORM='UNFORMATTED', ACCESS='STREAM')
          OPEN(UNIT=56, FILE=TRIM(PLT_TXT_FILE), STATUS='UNKNOWN', 
     &      FORM='FORMATTED')
          CLOSE(UNIT=56, STATUS='DELETE')
          OPEN(UNIT=56, FILE=TRIM(PLT_TXT_FILE), STATUS='NEW', 
     &      FORM='FORMATTED')
        ENDIF
        WRITE(6,'(4A)') 'Converting binary plot file: ',
     &    TRIM(PLT_BIN_FILE),' to text (ASCII) plot file: ',
     &    TRIM(PLT_TXT_FILE)
!
!---    Initialize binary file positions  ---
!
        IPOS = 1
        IPOSX = 1
        IPOSY = 1
        IPOSZ = 1
!
!---    Read time step  ---
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NSTEP
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NSTEP: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NSTEP: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of field grid cells
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of x-direction field grid cells
!
        READ(55, POS=IPOS, IOSTAT=ICHK) IFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'IFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'IFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of y-direction field grid cells
!
        READ(55, POS=IPOS, IOSTAT=ICHK) JFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'JFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'JFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of z-direction field grid cells
!
        READ(55, POS=IPOS, IOSTAT=ICHK) KFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'KFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'KFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of field processors
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NPFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of x-direction field processors
!
        READ(55, POS=IPOS, IOSTAT=ICHK) IPFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'IFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'IFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of y-direction field processors
!
        READ(55, POS=IPOS, IOSTAT=ICHK) JPFLD
        IF( ICHK.GT.0 ) THEN
          NCH = INDEX( PLT_BIN_FILE(1:),'  ' ) - 1
          WRITE(6,'(2A)') 'JFLD: End of file encountered during ' //
     &      'plot file read: ',PLT_BIN_FILE(1:NCH)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'JFLD: Read error encountered during plot ' //
     &      'file read: ',PLT_BIN_FILE(1:NCH)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of z-direction field processors
!
        READ(55, POS=IPOS, IOSTAT=ICHK) KPFLD
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'KFLD: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'KFLD: Read error encountered during plot ' //
     &      'file read: ',TRIM(PLT_BIN_FILE)
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
!---    Read number of inactive field grid cells
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NXP
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NXP: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NXP: Read error encountered during plot ' //
     &      'file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of plot variables
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NVPLOT
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NVPLOT: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NVPLOT: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read number of plot file significant digits
!
        READ(55, POS=IPOS, IOSTAT=ICHK) ISGNP
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'ISGNP: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'ISGNP: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
        WRITE(FORM2(8:9),'(I2)') MAX( 10,ISGNP+6 )
        WRITE(FORM2(11:11),'(I1)') MIN( 9,ISGNP-1 )
!  
!---    Read index for grid type
!
        READ(55, POS=IPOS, IOSTAT=ICHK) ICS
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'ICS: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'ICS: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read index for grid option
!
        READ(55, POS=IPOS, IOSTAT=ICHK) ISLC63
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'ISLC63: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'ISLC63: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!  
!---    Read simulation time, s
!
        READ(55, POS=IPOS, IOSTAT=ICHK) TM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'TM: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'TM: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTR
!  
!---    Read time units
!
        READ(55, POS=IPOS, IOSTAT=ICHK) UNTM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'UNTM: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'UNTM: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + 64*NBYTC
!  
!---    Read length units
!
        READ(55, POS=IPOS, IOSTAT=ICHK) UNLN
        IDB = INDEX(UNLN,'  ')-1
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'UNTM: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'UNTM: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + 64*NBYTC
!  
!---    Read number of vertices
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NVTC
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NVTC: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NVTC: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!
!---    Write header to text (ASCII) plot.xxx file  ---
!
        WRITE(56,'(A,//)')' Welcome to ...'
        WRITE(56,'(A)')   '                           STOMP'
        WRITE(56,'(A,//)')'        Subsurface Transport Over ' // 
     &    'Multiple Phases'
        WRITE(56,'(A)')   ' This file was produced by STOMP, ' //
     &    'a numerical simulator'
        WRITE(56,'(A)')   ' developed by the Pacific Northwest ' //
     &    'Laboratory, with'
        WRITE(56,'(A)')   ' support from the VOC-Arid Integrated ' //
     &    'Demonstration Project,'
        WRITE(56,'(A)')   ' Office of Technology Development, ' //
     &    'U.S. Department of Energy.'
        WRITE(56,'(A)')   ' Results from this version of STOMP ' //
     &    'should not be used for'
        WRITE(56,'(A,/)') ' license related applications.'
        WRITE(56,'(A,/)') ' For inquiries or assistance:  ' //
     &    'Call (509) 372-6070'
        WRITE(56,'(A,//)')'                        ---  PLOT  ---'

        WRITE(56,'(2A)') ' Version: STOMPX to STOMP Conversion'
        WRITE(56,'(2A)') ' Date: ','Date system call inactive'
        WRITE(56,'(2A,//)') ' Time: ','Time system call inactive'
        WRITE(56,'(/,A)') '--- Field Variable Data by Node Numbers'
        WRITE(56,'(A)') '    Flux Variable Data by Surface Numbers  ---'
        WRITE(56,'(/,A,I9)') 'Number of Time Steps = ',NSTEP
        WRITE(56,'(2A)') 'Time Output Units = ',TRIM(UNTM)
        TMIN = TM/60.
        THR = TM/3600.
        TDAY = THR/24.
        TWK = TDAY/7.
        TYR = TDAY/365.25D+0
        WRITE(56,'(A,6(1PE13.6,A))') 'Time = ',TM,',s ',TMIN,',min ',
     &    THR,',h ',TDAY,',day ',TWK,',wk ',TYR,',yr '
        WRITE(56,'(/,A,I9)') 'Number of Plot Variables = ',NVPLOT
        WRITE(56,'(A,I9)') 'Number of X or R-Direction Nodes = ',IFLD
        WRITE(56,'(A,I9)') 'Number of Y or Theta-Direction Nodes = ',
     &    JFLD
        WRITE(56,'(A,I9)') 'Number of Z-Direction Nodes = ',KFLD
        WRITE(56,'(A,I9)') 'Number of Field Nodes = ',NFLD
        WRITE(56,'(A,I9)') 'Number of Active Nodes = ',NFLD-NXP
        WRITE(56,'(A,I9)') 'Number of Vertices = ',NVTC
!  
!---    Allocate memory for plot variables and vertices
!
        ALLOCATE( VAR(1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VAR'
          STOP
        ENDIF
        ALLOCATE( VARX(1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VARX'
          STOP
        ENDIF
        ALLOCATE( IXP(1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: IXP'
          STOP
        ENDIF
        ALLOCATE( IXPX(1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: IXPX'
          STOP
        ENDIF
        ALLOCATE( VER(1:NVTC,1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VER'
          STOP
        ENDIF
        ALLOCATE( VERX(1:NVTC,1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VERX'
          STOP
        ENDIF
        ALLOCATE( ID_MPI(1:2,1:IPFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: ID_MPI'
          STOP
        ENDIF
        ALLOCATE( JD_MPI(1:2,1:JPFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: JD_MPI'
          STOP
        ENDIF
        ALLOCATE( KD_MPI(1:2,1:KPFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: KD_MPI'
          STOP
        ENDIF
        ALLOCATE( ND(1:IFLD,1:JFLD,1:KFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: ND'
          STOP
        ENDIF
!
!---    Define the grid pointers  ---
!
        DO K = 1,KFLD
          DO J = 1,JFLD
            DO I = 1,IFLD
              ND(I,J,K) = (K-1)*IJFLD + (J-1)*IFLD + I
            ENDDO
          ENDDO
        ENDDO
!
!---    Fill processor indice limits  ---
!
        IL = IFLD/IPFLD
        JL = JFLD/JPFLD
        KL = KFLD/KPFLD
        DO IP = 1,IPFLD
          IF( IP.EQ.1 ) THEN
            ID_MPI(1,IP) = 1
            ID_MPI(2,IP) = IL
          ELSEIF( IP.EQ.IPFLD ) THEN
            ID_MPI(1,IP) = (IP-1)*IL + 1
            ID_MPI(2,IP) = IFLD
          ELSE
            ID_MPI(1,IP) = (IP-1)*IL + 1
            ID_MPI(2,IP) = IP*IL
          ENDIF
        ENDDO
        DO JP = 1,JPFLD
          IF( JP.EQ.1 ) THEN
            JD_MPI(1,JP) = 1
            JD_MPI(2,JP) = JL
          ELSEIF( JP.EQ.JPFLD ) THEN
            JD_MPI(1,JP) = (JP-1)*JL + 1
            JD_MPI(2,JP) = JFLD
          ELSE
            JD_MPI(1,JP) = (JP-1)*JL + 1
            JD_MPI(2,JP) = JP*JL
          ENDIF
        ENDDO
        DO KP = 1,KPFLD
          IF( KP.EQ.1 ) THEN
            KD_MPI(1,KP) = 1
            KD_MPI(2,KP) = KL
          ELSEIF( KP.EQ.KPFLD ) THEN
            KD_MPI(1,KP) = (KP-1)*KL + 1
            KD_MPI(2,KP) = KFLD
          ELSE
            KD_MPI(1,KP) = (KP-1)*KL + 1
            KD_MPI(2,KP) = KP*KL
          ENDIF
        ENDDO
        IF( IPFLD.EQ.1 ) THEN
          ID_MPI(1,IPFLD) = 1
          ID_MPI(2,IPFLD) = IFLD
        ENDIF
        IF( JPFLD.EQ.1 ) THEN
          JD_MPI(1,JPFLD) = 1
          JD_MPI(2,JPFLD) = JFLD
        ENDIF
        IF( KPFLD.EQ.1 ) THEN
          KD_MPI(1,KPFLD) = 1
          KD_MPI(2,KPFLD) = KFLD
        ENDIF
!  
!---    Read x vertices
!
        READ(51, POS=IPOSX, IOSTAT=ICHK) ((VER(M,N),M=1,NVTC),N=1,NFLD)
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'XE: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'XE: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOSX = IPOSX + NVTC*NBYTR*NFLD
!
!---    Write x-direction vertices
!
        WRITE(56,'(/,2A)') 'X-Direction Nodal Vertices, ',TRIM(UNLN)
        NC = 0
        DO KP = 1,KPFLD
          DO JP = 1,JPFLD
            DO IP = 1,IPFLD
              DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                  DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                    NC = NC + 1
                    DO M = 1,NVTC
                      VERX(M,ND(I,J,K)) = VER(M,NC)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO               
        DO N = 1,NFLD
          WRITE(56,'(8(1PE16.9,1X))') (VERX(M,N),M=1,NVTC)
        ENDDO
!  
!---    Read y vertices
!
        READ(52, POS=IPOSY, IOSTAT=ICHK) ((VER(M,N),M=1,NVTC),N=1,NFLD)
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'YE: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'YE: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOSY = IPOSY + NVTC*NBYTR*NFLD
!
!---    Write Y-Direction vertices
!
        WRITE(56,'(/,2A)') 'Y-Direction Nodal Vertices, ',TRIM(UNLN)
        NC = 0
        DO KP = 1,KPFLD
          DO JP = 1,JPFLD
            DO IP = 1,IPFLD
              DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                  DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                    NC = NC + 1
                    DO M = 1,NVTC
                      VERX(M,ND(I,J,K)) = VER(M,NC)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO               
        DO N = 1,NFLD
          WRITE(56,'(8(1PE16.9,1X))') (VERX(M,N),M=1,NVTC)
        ENDDO
!  
!---    Read z vertices
!
        READ(53, POS=IPOSZ, IOSTAT=ICHK) ((VER(M,N),M=1,NVTC),N=1,NFLD)
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'ZE: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'ZE: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOSZ = IPOSZ + NVTC*NBYTR*NFLD
!
!---    Write Z-Direction vertices
!
        WRITE(56,'(/,2A)') 'Z-Direction Nodal Vertices, ',TRIM(UNLN)
        NC = 0
        DO KP = 1,KPFLD
          DO JP = 1,JPFLD
            DO IP = 1,IPFLD
              DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                  DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                    NC = NC + 1
                    DO M = 1,NVTC
                      VERX(M,ND(I,J,K)) = VER(M,NC)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO               
        DO N = 1,NFLD
          WRITE(56,'(8(1PE16.9,1X))') (VERX(M,N),M=1,NVTC)
        ENDDO
!  
!---    Deallocate memory
!
        DEALLOCATE( VER,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: VER'
          STOP
        ENDIF
        DEALLOCATE( VERX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: VERX'
          STOP
        ENDIF
!  
!---    Read grid cell volume
!
        READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(N),N=1,NFLD)
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'VOL: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'VOL: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTR*NFLD
!
!---    Write grid cell volume  ---
!
        WRITE(56,'(/,3A)') 'Node Volume, ',TRIM(UNLN),'^3'
        NC = 0
        DO KP = 1,KPFLD
          DO JP = 1,JPFLD
            DO IP = 1,IPFLD
              DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                  DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                    NC = NC + 1
                    VARX(ND(I,J,K)) = VAR(NC)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO               
        WRITE(56,FORM2) (VARX(N),N=1,NFLD)
!  
!---    Read active grid cell number
!
        READ(55, POS=IPOS, IOSTAT=ICHK) (IXPX(N),N=1,NFLD)
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'IXP: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'IXP: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + 4*NFLD
!
!---    Write node map  ---
!
        WRITE(FORM3(6:6),'(I1)') ICOUNT(NFLD)+1
        WRITE(56,'(/,A)') 'Node Map'
        NC = 0
        DO KP = 1,KPFLD
          DO JP = 1,JPFLD
            DO IP = 1,IPFLD
              DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                  DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                    NC = NC + 1
                    IF( IXP(NC).EQ.0 ) THEN
                      IXP(ND(I,J,K)) = 0
                    ELSE
                      IXP(ND(I,J,K)) = 1
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        WRITE(56,FORM3) (IXP(N),N=1,NFLD)
!  
!---    Deallocate memory
!
        DEALLOCATE( IXP,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: IXP'
          STOP
        ENDIF
        DEALLOCATE( IXPX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: IXPX'
          STOP
        ENDIF
!  
!---    Plot file variables
!
        DO NV = 1,NVPLOT
!  
!---      Read plot variable name and units
!
          READ(55, POS=IPOS, IOSTAT=ICHK) PLTVNM
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'PLTVNM: End of file encountered during ' //
     &        'plot file read: ',TRIM(PLT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'PLTVNM: Read error encountered during ' //
     &        'plot file read: ',TRIM(PLT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + 64*NBYTC
!  
!---      Write plot variable name and units
!
          WRITE(56,'(/,2A)') TRIM(PLTVNM)
!  
!---      Read plot variable
!
          READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(N),N=1,NFLD)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(3A)') TRIM(PLTVNM),': End of file encountered ' //
     &        'during plot file read: ',TRIM(PLT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(3A)') TRIM(PLTVNM),': Read error encountered ' //
     &        'during plot file read: ',TRIM(PLT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NBYTR*NFLD
          NC = 0
          DO KP = 1,KPFLD
            DO JP = 1,JPFLD
              DO IP = 1,IPFLD
                DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                  DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                    DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                      NC = NC + 1
                      VARX(ND(I,J,K)) = VAR(NC)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO               
!  
!---      Write plot variable
!
          WRITE(56,FORM2) (VARX(N),N=1,NFLD)
        ENDDO
!  
!---    Deallocate memory
!
        DEALLOCATE( VAR,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: VAR'
          STOP
        ENDIF
        DEALLOCATE( VARX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: VARX'
          STOP
        ENDIF
!
!---    Close the binary plot_xxx.bin file  ---
!
        CLOSE( UNIT=55 )
!
!---    Close the text (ASCII) plot.xxx file  ---
!
        CLOSE( UNIT=56 )
      ENDDO
      STOP
!
!---  End of PLOT_BIN_TO_TXT program  ---
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

!----------------------Program-----------------------------------------!
!
      PROGRAM PLOT_BIN_TO_CSV
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
!     formatted files to csv format for plotting.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 1 December 2021.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XE,YE,ZE
      REAL*8, DIMENSION(:), ALLOCATABLE :: XP,YP,ZP,VOL,VARX
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: VAR
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: VARNM
      CHARACTER(64) :: CDUM,PLT_BIN_FILE,OUTFILE
      CHARACTER(64) :: UNTM,UNLN
      CHARACTER(4) FORM1
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
      DATA FORM1 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
      NCLA = IARGC()
      IFILE = 0
      IF( NCLA.EQ.0 ) THEN
        WRITE(6,'(/,A,/)')'Enter a single plot filename: '
        READ(5,*) PLT_BIN_FILE
        CALL LCASE( PLT_BIN_FILE )
        CDUM = PLT_BIN_FILE
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts binary formatted ' // 
     &      'plot_xxx.bin files created from executing STOMPX to a '
          WRITE(6,'(A)') '"csv" (comma separated variables) ' //
     &      'formatted file. The script takes one or more argument(s), '
          WRITE(6,'(A)') 'the names of plot files to convert. ' //
     &      'For each plot file in the argument list, a new file '
          WRITE(6,'(A)') 'will be created that contains comma ' //
     &      'separated columns of x, y, z node centroids followed '
          WRITE(6,'(A)') 'by all variables in the plot files at ' //
     &      'the node centroids. The output file names will be '
          WRITE(6,'(A)') 'plot-1234.csv where 1234 was the time ' //
     &      'step extension of the plot files in the argument list. '
          WRITE(6,'(A,/)') 'If an argument is missing the user will ' //
     &      'be prompted for a single plot file. '
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'plotbinTo_tec.exe plot_*.bin'
          WRITE(6,'(A)') 'plotbinTo_tec.exe plot_0000007.bin'
          WRITE(6,'(A,/)') 'plotbinTo_tec.exe ' // 
     &      'plot_0000007.bin plot_0000057.bin'
          STOP
        ENDIF
        NCLA = 1
        IFILE = 1
      ENDIF
      DO NA = 1,NCLA
        IF( IFILE.EQ.0 ) THEN
          CALL GETARG( NA,CDUM )
          CALL LCASE( CDUM )
          IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
            WRITE(6,'(/,A)') 'This script converts binary formatted ' // 
     &        'plot_xxx.bin files created from executing STOMPX to a '
            WRITE(6,'(A)') '"csv" (comma separated variables) ' //
     &        'formatted file. The script takes one or more '
            WRITE(6,'(A)') 'argument(s), the names of plot files to ' //
     &        'convert. For each plot file in the argument list, a '
            WRITE(6,'(A)') 'new file will be created that contains ' //
     &        'comma separated columns of x, y, z node centroids '
            WRITE(6,'(A)') 'followed by all variables in the plot ' //
     &        'files at the node centroids. The output file names '
            WRITE(6,'(A)') 'will be plot-1234.csv where 1234 was the' //
     &        ' time step extension of the plot files in the argument '
            WRITE(6,'(A,/)') 'list. If an argument is missing the ' //
     &        'user will be prompted for a single plot file. '
            WRITE(6,'(A,/)') 'Example command line entries'
            WRITE(6,'(A)') 'plot_to_csv.exe plot.*'
            WRITE(6,'(A)') 'plot_to_csv.exe plot.002'
            WRITE(6,'(A,/)') 'plot_to_csv.exe plot.01 plot.02 plot.10'
            STOP
          ELSE
            PLT_BIN_FILE = CDUM
          ENDIF
        ENDIF
!
!---    No plot_xgrid.bin file, exit program with error message  ---
!
        INQUIRE( FILE='plot_xgrid.bin', EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          WRITE(6,'(A)') 'ERROR: Nonexistent plot_xgrid.bin file'
          STOP
        ELSE
          OPEN(UNIT=51, FILE='plot_xgrid.bin', STATUS='OLD', 
     &      FORM='UNFORMATTED', ACCESS='STREAM')
        ENDIF
!
!---    No plot_ygrid.bin file, exit program with error message  ---
!
        INQUIRE( FILE='plot_ygrid.bin', EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          WRITE(6,'(A)') 'ERROR: Nonexistent plot_ygrid.bin file'
          STOP
        ELSE
          OPEN(UNIT=52, FILE='plot_ygrid.bin', STATUS='OLD', 
     &      FORM='UNFORMATTED', ACCESS='STREAM')
        ENDIF
!
!---    No plot_zgrid.bin file, exit program with error message  ---
!
        INQUIRE( FILE='plot_zgrid.bin', EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          WRITE(6,'(A)') 'ERROR: Nonexistent plot_zgrid.bin file'
          STOP
        ELSE
          OPEN(UNIT=53, FILE='plot_zgrid.bin', STATUS='OLD', 
     &      FORM='UNFORMATTED', ACCESS='STREAM')
        ENDIF
!
!---    No plot file, exit program with error message  ---
!
        INQUIRE( FILE=PLT_BIN_FILE, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          WRITE(6,'(A)') 'ERROR: Nonexistent plot file: ',
     &      TRIM(PLT_BIN_FILE)
          STOP
        ELSE
          OPEN(UNIT=55, FILE=PLT_BIN_FILE, STATUS='OLD', 
     &      FORM='UNFORMATTED', ACCESS='STREAM')
        ENDIF
!
!---    Write banner to screen  ---
!       
        IF( NA.EQ.1 ) THEN
          WRITE(6,'(A,//)')' Welcome to ...'
          WRITE(6,'(A,/)') '                       PLOT_BIN_TO_CSV'
          WRITE(6,'(A)')   ' An output file converter routine for ' //
     &      'the STOMP simulator.  This'
          WRITE(6,'(A)')   ' routine, developed by the Pacific ' //
     &      'Northwest National Laboratory,'
          WRITE(6,'(A)')   ' reads one or more STOMPX plot file(s) ' //
     &      ' and converts the data to one '
          WRITE(6,'(A)')   ' or more plot files with an extension ' //
     &      'of ".csv". These files contain '
          WRITE(6,'(A)')   ' columns of user-requested data at ' //
     &      'every grid cell, starting with the'
          WRITE(6,'(A,/)')   ' x, y, and z centroid coordinates and ' //
     &      'volume of the grid cell.'
          WRITE(6,'(A)')   ' For support:  Tel: 509.372.6070'
          WRITE(6,'(A,/)') '               E-mail:  mark.white@pnnl.gov'
        ENDIF
!
!---    Open plot-1234.csv file  ---
!
        OUTFILE = PLT_BIN_FILE
        NDOT = INDEX(OUTFILE(1:),'.bin')
        OUTFILE(NDOT:NDOT+3) = '.csv'
        OUTFILE = TRIM(ADJUSTL(OUTFILE))
        WRITE(6,'(2A)') 'Creating file: ',TRIM(OUTFILE)
        OPEN(UNIT=56, FILE=OUTFILE)
        CLOSE(UNIT=56,STATUS='DELETE')
        OPEN(UNIT=56, FILE=OUTFILE, STATUS='NEW', FORM='FORMATTED')
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
!---    Allocate memory
!
        ALLOCATE( VARNM(1:NVPLOT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VARNM'
          STOP
        ENDIF
        ALLOCATE( VAR(1:NFLD,1:NVPLOT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VAR'
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
        ALLOCATE( XE(1:NVTC,1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: XE'
          STOP
        ENDIF
        ALLOCATE( YE(1:NVTC,1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VER'
          STOP
        ENDIF
        ALLOCATE( ZE(1:NVTC,1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VER'
          STOP
        ENDIF
        ALLOCATE( XP(1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: XP'
          STOP
        ENDIF
        ALLOCATE( YP(1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: YP'
          STOP
        ENDIF
        ALLOCATE( ZP(1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: ZP'
          STOP
        ENDIF
        ALLOCATE( VOL(1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VOL'
          STOP
        ENDIF
        ALLOCATE( VARX(1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VARX'
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
        DO N = 1,NFLD
          DO NV = 1,NVTC
            XE(NV,N) = 0.D+0
            YE(NV,N) = 0.D+0
            ZE(NV,N) = 0.D+0
          ENDDO
          XP(N) = 0.D+0
          YP(N) = 0.D+0
          ZP(N) = 0.D+0
          VOL(N) = 0.D+0
        ENDDO
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
        READ(51, POS=IPOSX, IOSTAT=ICHK) ((XE(M,N),M=1,NVTC),N=1,NFLD)
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
        NC = 0
        DO KP = 1,KPFLD
          DO JP = 1,JPFLD
            DO IP = 1,IPFLD
              DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                  DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                    NC = NC + 1
                    XP(ND(I,J,K)) = 0.D+0
                    DO M = 1,NVTC
                      XP(ND(I,J,K)) = XP(ND(I,J,K)) + XE(M,NC)
                    ENDDO
                    XP(ND(I,J,K)) = XP(ND(I,J,K))/REAL(NVTC)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO               
!  
!---    Read y vertices
!
        READ(52, POS=IPOSY, IOSTAT=ICHK) ((YE(M,N),M=1,NVTC),N=1,NFLD)
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
        NC = 0
        DO KP = 1,KPFLD
          DO JP = 1,JPFLD
            DO IP = 1,IPFLD
              DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                  DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                    NC = NC + 1
                    YP(ND(I,J,K)) = 0.D+0
                    DO M = 1,NVTC
                      YP(ND(I,J,K)) = YP(ND(I,J,K)) + YE(M,NC)
                    ENDDO
                    YP(ND(I,J,K)) = YP(ND(I,J,K))/REAL(NVTC)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO               
!  
!---    Read z vertices
!
        READ(53, POS=IPOSZ, IOSTAT=ICHK) ((ZE(M,N),M=1,NVTC),N=1,NFLD)
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
        NC = 0
        DO KP = 1,KPFLD
          DO JP = 1,JPFLD
            DO IP = 1,IPFLD
              DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                  DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                    NC = NC + 1
                    YP(ND(I,J,K)) = 0.D+0
                    DO M = 1,NVTC
                      ZP(ND(I,J,K)) = ZP(ND(I,J,K)) + ZE(M,NC)
                    ENDDO
                    ZP(ND(I,J,K)) = ZP(ND(I,J,K))/REAL(NVTC)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO               
!  
!---    Read grid cell volume
!
        READ(55, POS=IPOS, IOSTAT=ICHK) (VARX(N),N=1,NFLD)
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
        NC = 0
        DO KP = 1,KPFLD
          DO JP = 1,JPFLD
            DO IP = 1,IPFLD
              DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                  DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                    NC = NC + 1
                    VOL(ND(I,J,K)) = VARX(NC)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO               
!  
!---    Read active grid cell number
!
        READ(55, POS=IPOS, IOSTAT=ICHK) (IXPX(N),N=1,NFLD)
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'IXPX: End of file encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'IXPX: Read error encountered during ' //
     &      'plot file read: ',TRIM(PLT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + 4*NFLD
        NC = 0
        DO KP = 1,KPFLD
          DO JP = 1,JPFLD
            DO IP = 1,IPFLD
              DO K = KD_MPI(1,KP),KD_MPI(2,KP)
                DO J = JD_MPI(1,JP),JD_MPI(2,JP)
                  DO I = ID_MPI(1,IP),ID_MPI(2,IP)
                    NC = NC + 1
                    IF( IXPX(NC).EQ.0 ) THEN
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
!  
!---    Plot file variables
!
        DO NV = 1,NVPLOT
!  
!---      Read plot variable name and units
!
          READ(55, POS=IPOS, IOSTAT=ICHK) VARNM(NV)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'VARNM: End of file encountered during ' //
     &        'plot file read: ',TRIM(PLT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'VARNM: Read error encountered during ' //
     &        'plot file read: ',TRIM(PLT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + 64*NBYTC
!  
!---      Read plot variable
!
          READ(55, POS=IPOS, IOSTAT=ICHK) (VARX(N),N=1,NFLD)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(3A)') TRIM(VARNM(NV)),': End of file ' //
     &        'encountered during plot file read: ',TRIM(PLT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(3A)') TRIM(VARNM(NV)),': Read error ' //
     &        'encountered during plot file read: ',TRIM(PLT_BIN_FILE)
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
                      VAR(ND(I,J,K),NV) = VARX(NC)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO               
        ENDDO
!
!---    Convert time to UNTM units  ---
!
        IF( INDEX(UNTM(1:),'s').NE.0 ) THEN
          TMX = TM
        ELSEIF( INDEX(UNTM(1:),'sec').NE.0 ) THEN
          TMX = TM
        ELSEIF( INDEX(UNTM(1:),'min').NE.0 ) THEN
          TMX = TM/6.D+1
        ELSEIF( INDEX(UNTM(1:),'hour').NE.0 ) THEN
          TMX = TM/3.6D+3
        ELSEIF( INDEX(UNTM(1:),'hr').NE.0 ) THEN
          TMX = TM/3.6D+3
        ELSEIF( INDEX(UNTM(1:),'day').NE.0 ) THEN
          TMX = TM/3.6D+3/2.4D+1
        ELSEIF( INDEX(UNTM(1:),'d').NE.0 ) THEN
          TMX = TM/3.6D+3/2.4D+1
        ELSEIF( INDEX(UNTM(1:),'week').NE.0 ) THEN
          TMX = TM/3.6D+3/2.4D+1/7.D+0
        ELSEIF( INDEX(UNTM(1:),'wk').NE.0 ) THEN
          TMX = TM/3.6D+3/2.4D+1/7.D+0
        ELSEIF( INDEX(UNTM(1:),'year').NE.0 ) THEN
          TMX = TM/3.6D+3/2.4D+1/3.6525D+2
        ELSEIF( INDEX(UNTM(1:),'yr').NE.0 ) THEN
          TMX = TM/3.6D+3/2.4D+1/3.6525D+2
        ENDIF
!
!---    Write header line to plot-1234.csv file  ---
!
        IF( IFLD.GT.1 ) WRITE(56,'(3A,$)') 
     &   'X-Direction Grid-Cell Centroid (',TRIM(UNLN),'),'
        IF( JFLD.GT.1 ) WRITE(56,'(3A,$)') 
     &   'Y-Direction Grid-Cell Centroid (',TRIM(UNLN),'),'
        IF( KFLD.GT.1 ) WRITE(56,'(3A,$)') 
     &   'Z-Direction Grid-Cell Centroid (',TRIM(UNLN),'),'
        WRITE(56,'(3A,$)') 
     &   'Grid-Cell Volume (',TRIM(UNLN),'^3),'
        DO NV = 1,NVPLOT
          NCM = INDEX(VARNM(NV)(1:),', ')
          IF( NCM.NE.0 ) THEN
            VARNM(NV)(NCM:NCM+1) = ' ('
            NCH = INDEX(VARNM(NV)(1:),'  ')
            VARNM(NV)(NCH:NCH) = ')'
          ENDIF
          WRITE(56,'(2A,$)') TRIM(VARNM(NV)),','
        ENDDO
        WRITE(56,'(A)') ' '
!
!---    Write variable data to plot-1234.csv file  ---
!
        DO N = 1,NFLD
          IF( IFLD.GT.1 ) WRITE(56,'(1PE14.7,A,$)') XP(N),','
          IF( JFLD.GT.1 ) WRITE(56,'(1PE14.7,A,$)') YP(N),','
          IF( KFLD.GT.1 ) WRITE(56,'(1PE14.7,A,$)') ZP(N),','
          WRITE(56,'(1PE14.7,A,$)') VOL(N),','
          DO NV = 1,NVPLOT
            WRITE(56,'(1PE14.7,A,$)') VAR(N,NV),','
          ENDDO
          WRITE(56,'(A)') ' '
        ENDDO
!  
!---    Deallocate memory
!
        DEALLOCATE( VARNM,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: VARNM'
          STOP
        ENDIF
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
        DEALLOCATE( XE,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: XE'
          STOP
        ENDIF
        DEALLOCATE( YE,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: YE'
          STOP
        ENDIF
        DEALLOCATE( ZE,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: ZE'
          STOP
        ENDIF
        DEALLOCATE( XP,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: XP'
          STOP
        ENDIF
        DEALLOCATE( YP,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: YP'
          STOP
        ENDIF
        DEALLOCATE( ZP,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: ZP'
          STOP
        ENDIF
        DEALLOCATE( VOL,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: VOL'
          STOP
        ENDIF
!
!---    Close the binary plot_xxx.bin file  ---
!
        CLOSE( UNIT=55 )
        CLOSE( UNIT=56 )
      ENDDO
      STOP
!
!---  End of PLOT_BIN_TO_CSV program  ---
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

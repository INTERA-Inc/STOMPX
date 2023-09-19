!----------------------Program-----------------------------------------!
!
      PROGRAM PLOT_FRCTO_TEC
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
!     STOMP: Subsurface Transport Over Multiple Phases
!
!     This utility routine converts one or more plot_frc files to 
!     Tecplot data format for plotting.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 March 2021.
!     plot_frcTo_tec.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(10) :: VARX
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XE,YE,ZE,VAR
      REAL*8, DIMENSION(:), ALLOCATABLE :: XP,YP,ZP,AREA,APER
      INTEGER, DIMENSION(:), ALLOCATABLE :: NMAP
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NCMAP
      CHARACTER(32), DIMENSION(5) :: UNITS
      CHARACTER(32), DIMENSION(6) :: TIME
      CHARACTER(128), DIMENSION(:), ALLOCATABLE :: VARNM
      CHARACTER(128) :: CDUM,PLTFILE,PLTTYPE,OUTFILE,CNTFILE
      CHARACTER(16) :: UNTM
      CHARACTER(6) :: FORM1
      CHARACTER(512) :: CHDUM
      LOGICAL FCHK
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1,$)'/
!
!----------------------Executable Lines--------------------------------!
!
      DO M = 1,5
        UNITS(M) = ''
      ENDDO
      CNTFILE = 'connect_frc'
      NCLA = IARGC()
      IFILE = 0
      IF( NCLA.EQ.0 ) THEN
        WRITE(6,'(/,A,/)')'Enter "contour" or "help": '
        READ(5,*) PLTTYPE
        CALL LCASE( PLTTYPE )
        CDUM = PLTTYPE
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts outputs in the ' //
     &      'data section of one or more STOMP plot_frc file(s) to a '
          WRITE(6,'(A)') 'Tecplot data formatted file with an ".dat"' //
     &      ' extension. The script takes three or more argument(s): '
          WRITE(6,'(A)') '  Plotting Option: "contour"'
          WRITE(6,'(A)') '  Output File Name: (e.g. plot_frcs.dat)'
          WRITE(6,'(A)') '  Plot File(s): (e.g. plot_frc.*)'
          WRITE(6,'(A)') 'Currently only the "contour" plotting ' //
     &       'option is available. The output file is created by this '
          WRITE(6,'(A)') 'script, which will become the input ' //
     &      'file to Tecplot. The data in each plot file, will be '
          WRITE(6,'(A,/)') 'stored in the output file, as a ' //
     &      'separate zone, titled by the time stamp in the plot file.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'plot_frcTo_tec.exe contour plot_frcs.dat ' // 
     &      'plot_frc.*'
          WRITE(6,'(A)') 'plot_frcTo_tec.exe contour plot_frcs.dat ' // 
     &       'plot_frc.*'
          WRITE(6,'(A)') 'plot_frcTo_tec.exe contour pressure_frc.dat '
     &      // 'plot_frc.002'
          WRITE(6,'(A,/)') 'plot_frcTo_tec.exe contour run1_frc.dat ' // 
     &      'plot_frc.01 plot_frc.02 plot_frc.10'
          STOP
        ENDIF
        WRITE(6,'(/,A,/)')'Enter the output filename: '
        READ(5,*) OUTFILE
        WRITE(6,'(/,A,/)')'Enter a single plot filename ' // 
     &    '(e.g., plot_frc.001): '
        READ(5,*) PLTFILE
        NCLA = 3
        IFILE = 1
      ELSE
        CALL GETARG( 1,PLTTYPE )
        CALL LCASE( PLTTYPE )
        CALL GETARG( 2,OUTFILE )
        CDUM = PLTTYPE
        CALL LCASE( CDUM )
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts outputs in the ' //
     &      'data section of one or more STOMP plot_frc file(s) to a '
          WRITE(6,'(A)') 'Tecplot data formatted file with an ".dat"' //
     &      ' extension. The script takes three or more argument(s): '
          WRITE(6,'(A)') '  Plotting Option: "contour"'
          WRITE(6,'(A)') '  Output File Name: (e.g. plot_frcs.dat)'
          WRITE(6,'(A)') '  Plot File(s): (e.g. plot_frc.*)'
          WRITE(6,'(A)') 'Currently only the "contour" plotting ' //
     &       'option is available. The output file is created by this '
          WRITE(6,'(A)') 'script, which will become the input ' //
     &      'file to Tecplot. The data in each plot file, will be '
          WRITE(6,'(A,/)') 'stored in the output file, as a ' //
     &      'separate zone, titled by the time stamp in the plot file.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'plot_frcTo_tec.exe contour plot_frcs.dat ' // 
     &      'plot_frc.*'
          WRITE(6,'(A)') 'plot_frcTo_tec.exe contour plot_frcs.dat ' // 
     &      'plot_frc.*'
          WRITE(6,'(A)') 'plot_frcTo_tec.exe contour pressure_frc.dat '
     &      // 'plot_frc.002'
          WRITE(6,'(A,/)') 'plot_frcTo_tec.exe contour run1_frc.dat ' // 
     &      'plot_frc.01 plot_frc.02 plot_frc.10'
          STOP
        ENDIF
      ENDIF
!
!---  Check plotting type  ---
!
      IF( INDEX(PLTTYPE(1:),'contour').NE.1 ) THEN
        NCH = INDEX(PLTTYPE(1:),' ')-1
        WRITE(6,'(A)') 'Unrecognized Plot Type: ',PLTTYPE(1:NCH)
        STOP
      ENDIF
!
!---  Open output file  ---
!
      OUTFILE = TRIM(ADJUSTL(OUTFILE))
      WRITE(6,'(2A)') 'Creating file: ',TRIM(OUTFILE)
      OPEN(UNIT=56, FILE=OUTFILE)
      CLOSE(UNIT=56,STATUS='DELETE')
      OPEN(UNIT=56, FILE=OUTFILE, STATUS='NEW', FORM='FORMATTED')
!
!---  No connect file, exit program with error message or open
!     connect file  ---
!
      INQUIRE( FILE=CNTFILE, EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        WRITE(6,'(A)') 'ERROR: Nonexistent "connect_frc" file.'
        STOP
      ELSE
        OPEN(UNIT=45, FILE=CNTFILE, STATUS='OLD', FORM='FORMATTED')
      ENDIF
!
!---  Write banner to screen  ---
!       
      WRITE(6,'(A,//)')' Welcome to ...'
      WRITE(6,'(A,/)') '                       PLOT_FRCTO_TEC'
      WRITE(6,'(A)')   ' An output file converter routine for ' //
     &  'the STOMP simulator.  This'
      WRITE(6,'(A)')   ' routine, developed by the Pacific ' //
     &  'Northwest National Laboratory,'
      WRITE(6,'(A)')   ' reads one or more STOMP plot_frc file(s) ' //
     &  ' and converts the data to one '
      WRITE(6,'(A)')   ' Tecplot data file with an extension ' //
     &  'of ".dat". This file contains the '
      WRITE(6,'(A)')   'user-requested data at every grid cell' //
     &  'in Tecplot format for visualization with Tecplot. '
      WRITE(6,'(A)')   'Currently fracture/fault domains can only ' //
     &  'be visualized as contour plots.'
      WRITE(6,'(A)')   ' For support:  Tel: 509.372.6070'
      WRITE(6,'(A,/)') '               E-mail:  mark.white@pnnl.gov'
!
!---  Loop over STOMP plot_frc files ---
!
      DO NA = 3,NCLA
        IF( IFILE.EQ.0 ) THEN
          CALL GETARG( NA,PLTFILE )
        ENDIF
!
!---    No plot_frc file, exit program with error message  ---
!
        INQUIRE( FILE=PLTFILE, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          NCH = INDEX(PLTFILE(1:),'  ')-1
          WRITE(6,'(A)') 'ERROR: Nonexistent plot_frc file: ',
     &      PLTFILE(1:NCH)
          STOP
        ELSE
          OPEN(UNIT=55, FILE=PLTFILE, STATUS='OLD', FORM='FORMATTED')
        ENDIF
        NCH = INDEX(PLTFILE(1:),'  ')-1
        WRITE(6,'(2A)') 'Processing plot_frc file: ',PLTFILE(1:NCH)
!
!---    Search for coordinate dimensions  ---
!
        DO
          READ(55,'(A)', IOSTAT=ICHK) CHDUM
          IF( ICHK.GT.0 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'End of file encountered during plot_frc '
     &         // 'file read: ',PLTFILE(1:NCH)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Read error encountered during plot_frc ' //
     &        'file read: ',PLTFILE(1:NCH)
            STOP
          ENDIF
          UNTMX = INDEX(CHDUM(1:),'Time Output Units =')
          TIMEX = INDEX(CHDUM(1:),'Time =')
          NVPTX = INDEX(CHDUM(1:),'Number of Plot Variables =')
          NFTX = INDEX(CHDUM(1:),'Number of Fracture Triangles =')
          NAFTX = INDEX(CHDUM(1:),'Number of Active Fracture Triang')
          IF( UNTMX.NE.0 ) THEN
            UNTMX = INDEX(CHDUM(1:),'=')
            I1 = UNTMX+2
            I2 = UNTMX+16
            READ(CHDUM(I1:I2),'(A)') UNTM
            IF( INDEX(UNTM(1:1),'s').NE.0 ) THEN
              IUNTM = 1
            ELSEIF( INDEX(UNTM(1:1),'m').NE.0 ) THEN
              IUNTM = 2
            ELSEIF( INDEX(UNTM(1:1),'h').NE.0 ) THEN
              IUNTM = 3
            ELSEIF( INDEX(UNTM(1:1),'d').NE.0 ) THEN
              IUNTM = 4
            ELSEIF( INDEX(UNTM(1:1),'w').NE.0 ) THEN
              IUNTM = 5
            ELSEIF( INDEX(UNTM(1:1),'y').NE.0 ) THEN
              IUNTM = 6
            ENDIF
          ENDIF
          IF( TIMEX.NE.0 ) THEN
            TIMEX = INDEX(CHDUM(1:),'=')
            I1 = TIMEX+3
            DO M = 1,6
              I2 = I1 + 11
              TIME(M) = CHDUM(I1:I2)
              I1 = INDEX(CHDUM(I2:),'  ') + I2 + 1
            ENDDO
          ENDIF
          IF( NVPTX.NE.0 ) THEN
            NVPTX = INDEX(CHDUM(1:),'=')
            I1 = NVPTX+1
            I2 = NVPTX+10
            READ(CHDUM(I1:I2),'(I10)') NVPLOT
          ENDIF
          IF( NFTX.NE.0 ) THEN
            NFTX = INDEX(CHDUM(1:),'=')
            N1 = NFTX+1
            N2 = NFTX+10
            READ(CHDUM(N1:N2),'(I10)') NFT
          ENDIF
          IF( NAFTX.NE.0 ) THEN
            NAFTX = INDEX(CHDUM(1:),'=')
            N1 = NAFTX+1
            N2 = NAFTX+10
            READ(CHDUM(N1:N2),'(I10)') NAFT
            EXIT
          ENDIF
        ENDDO  
        ITYPE = 0
        IF( NA.GT.3 ) THEN
          IF( NFT.NE.NFT1 .OR. NAFT.NE.NAFT1 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'Plot_frc file not compatible with ' //
     &        'previous plot_frc files: ',PLTFILE(1:NCH)
            STOP
          ENDIF
        ELSE
          NFT1 = NFT
          NAFT1 = NAFT
        ENDIF     
!  
!---    Allocate memory
!
        ALLOCATE( VARNM(1:NVPLOT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VARNM'
          STOP
        ENDIF
        ALLOCATE( VAR(1:NFT,1:NVPLOT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: VAR'
          STOP
        ENDIF
        ALLOCATE( XE(1:3,1:NFT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: XE'
          STOP
        ENDIF
        ALLOCATE( YE(1:3,1:NFT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: YE'
          STOP
        ENDIF
        ALLOCATE( ZE(1:3,1:NFT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: ZE'
          STOP
        ENDIF
        ALLOCATE( XP(1:NFT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: XP'
          STOP
        ENDIF
        ALLOCATE( YP(1:NFT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: YP'
          STOP
        ENDIF
        ALLOCATE( ZP(1:NFT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: ZP'
          STOP
        ENDIF
        ALLOCATE( AREA(1:NFT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: AREA'
          STOP
        ENDIF
        ALLOCATE( APER(1:NFT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: APER'
          STOP
        ENDIF
        ALLOCATE( NMAP(1:NFT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: NMAP'
          STOP
        ENDIF
        IF( NA.EQ.3 ) THEN
          ALLOCATE( NCMAP(1:3,1:NFT),STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            WRITE(6,'(A)') 'Allocation Error: NCMAP'
            STOP
          ENDIF
          DO N = 1,NFT
            DO NV = 1,3
              NCMAP(NV,N) = 0
            ENDDO
          ENDDO
        ENDIF
        DO N = 1,NFT
          DO NV = 1,3
            XE(NV,N) = 0.D+0
            YE(NV,N) = 0.D+0
            ZE(NV,N) = 0.D+0
          ENDDO
          XP(N) = 0.D+0
          YP(N) = 0.D+0
          ZP(N) = 0.D+0
          AREA(N) = 0.D+0
          APER(N) = 0.D+0
        ENDDO
!
!---    Search for x vertices  ---
!
        DO
          READ(55,'(A)', IOSTAT=ICHK) CHDUM
          IF( ICHK.GT.0 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &        'file read: X-Vertices: ',PLTFILE(1:NCH)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &        'file read: X-Vertices: ',PLTFILE(1:NCH)
            STOP
          ENDIF
          IXDNV = INDEX(CHDUM(1:),'X-Direction Fracture Triangle')
          IF( IXDNV.NE.0 ) THEN
            I1 = INDEX(CHDUM(1:),',')+1
            I2 = INDEX(CHDUM(1:),'   ')-1
            UNITS(1) = CHDUM(I1:I2)
            EXIT
          ENDIF
        ENDDO
        DO N = 1,NFT
          READ(55,*) (XE(M,N),M=1,3)
          XP(N) = 0.D+0
          DO M = 1,3
            XP(N) = XP(N) + XE(M,N)
          ENDDO
          XP(N) = XP(N)/3.D+0
        ENDDO
!
!---    Search for y vertices  ---
!
        DO
          READ(55,'(A)', IOSTAT=ICHK) CHDUM
          IF( ICHK.GT.0 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &        'file read: Y-Vertices: ',PLTFILE(1:NCH)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &        'file read: Y-Vertices: ',PLTFILE(1:NCH)
            STOP
          ENDIF
          IXDNV = INDEX(CHDUM(1:),'Y-Direction Fracture Triangle')
          IF( IXDNV.NE.0 ) THEN
            I1 = INDEX(CHDUM(1:),',')+1
            I2 = INDEX(CHDUM(1:),'   ')-1
            UNITS(2) = CHDUM(I1:I2)
            EXIT
          ENDIF
        ENDDO
        DO N = 1,NFT
          READ(55,*) (YE(M,N),M=1,3)
          YP(N) = 0.D+0
          DO M = 1,3
            YP(N) = YP(N) + YE(M,N)
          ENDDO
          YP(N) = YP(N)/3.D+0
        ENDDO
!
!---    Search for z vertices  ---
!
        DO
          READ(55,'(A)', IOSTAT=ICHK) CHDUM
          IF( ICHK.GT.0 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &        'file read: Z-Vertices: ',PLTFILE(1:NCH)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &        'file read: Z-Vertices: ',PLTFILE(1:NCH)
            STOP
          ENDIF
          IXDNV = INDEX(CHDUM(1:),'Z-Direction Fracture Triangle')
          IF( IXDNV.NE.0 ) THEN
            I1 = INDEX(CHDUM(1:),',')+1
            I2 = INDEX(CHDUM(1:),'   ')-1
            UNITS(3) = CHDUM(I1:I2)
            EXIT
          ENDIF
        ENDDO
        DO N = 1,NFT
          READ(55,*) (ZE(M,N),M=1,3)
          ZP(N) = 0.D+0
          DO M = 1,3
            ZP(N) = ZP(N) + ZE(M,N)
          ENDDO
          ZP(N) = ZP(N)/3.D+0
        ENDDO
!
!---    Search for fracture triangle map  ---
!
        DO
          READ(55,'(A)', IOSTAT=ICHK) CHDUM
          IF( ICHK.GT.0 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &        'file read: Fracture Triangle Map: ',PLTFILE(1:NCH)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &        'file read: Fracture Triangle Map: ',PLTFILE(1:NCH)
            STOP
          ENDIF
          IXDNV = INDEX(CHDUM(1:),'Fracture Triangle Map')
          IF( IXDNV.NE.0 ) EXIT
        ENDDO
        DO N = 1,NFT,10
          READ(55,*) (NMAP(M),M=N,MIN(N+9,NFT))
        ENDDO
!
!---    Read connection map  ---
!
        IF( NA.EQ.3 ) THEN
          DO N = 1,NFT
            IF( NMAP(N).EQ.0 ) CYCLE
            READ(45,*) (NCMAP(NV,N),NV=1,3)
          ENDDO
!
!---      Close the connect file  ---
!
          CLOSE( UNIT=45 )
        ENDIF
!
!---    Search for fracture triangle area  ---
!
        DO
          READ(55,'(A)', IOSTAT=ICHK) CHDUM
          IF( ICHK.GT.0 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &        'file read: Fracture Triangle Area: ',PLTFILE(1:NCH)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &        'file read: Fracture Triangle Area: ',PLTFILE(1:NCH)
            STOP
          ENDIF
          IXDNV = INDEX(CHDUM(1:),'Fracture Triangle Area')
          IF( IXDNV.NE.0 ) THEN
            I1 = INDEX(CHDUM(1:),',')+1
            I2 = INDEX(CHDUM(1:),'   ')-1
            UNITS(4) = CHDUM(I1:I2)
            EXIT
          ENDIF
        ENDDO
        DO N = 1,NFT,10
          READ(55,*) (AREA(M),M=N,MIN(N+9,NFT))
        ENDDO
!
!---    Search for fracture triangle aperture  ---
!
        DO
          READ(55,'(A)', IOSTAT=ICHK) CHDUM
          IF( ICHK.GT.0 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &        'file read: Fracture Triangle Aperture: ',PLTFILE(1:NCH)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &        'file read: Fracture Triangle Aperture: ',PLTFILE(1:NCH)
            STOP
          ENDIF
          IXDNV = INDEX(CHDUM(1:),'Fracture Triangle Aperture')
          IF( IXDNV.NE.0 ) THEN
            I1 = INDEX(CHDUM(1:),',')+1
            I2 = INDEX(CHDUM(1:),'   ')-1
            UNITS(5) = CHDUM(I1:I2)
            EXIT
          ENDIF
        ENDDO
        DO N = 1,NFT,10
          READ(55,*) (APER(M),M=N,MIN(N+9,NFT))
        ENDDO
!
!---    Search for variable data  ---
!
        DO NV = 1,NVPLOT
          DO
            READ(55,'(A)', IOSTAT=ICHK) CHDUM
            IF( ICHK.GT.0 ) THEN
              NCH = INDEX( PLTFILE(1:),'  ' ) - 1
              WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &          'file read: Variable Data: ',PLTFILE(1:NCH)
              STOP
            ELSEIF( ICHK.LT.0 ) THEN
              WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &          'file read: Variable Data: ',PLTFILE(1:NCH)
              STOP
            ENDIF
            IF( (ICHAR(CHDUM(1:1)).GE.97 .AND. ICHAR(CHDUM(1:1)).LE.122)
     &        .OR. (ICHAR(CHDUM(1:1)).GE.65 .AND.
     &        ICHAR(CHDUM(1:1)).LE.90) ) THEN
              NCH = INDEX(CHDUM(1:),'   ')-1
              VARNM(NV) = CHDUM(1:NCH)
              EXIT
            ENDIF
          ENDDO
          DO N = 1,NFT,10
            READ(55,*) (VAR(M,NV),M=N,MIN(N+9,NFT))
          ENDDO
        ENDDO
!
!---    Write file version and identifier to the Tecplot data file  ---
!
        IF( NA.EQ.3 ) THEN
!
!---      Contour plotting type  ---
!
          IF( ITYPE.EQ.0 ) THEN
            WRITE(56,'(A)') 'TITLE = "STOMP Fracture/Fault Plot ' // 
     &        'File Visualization"'
            WRITE(56,'(A,$)') 'VARIABLES = '
            NCH = INDEX(UNITS(1)(1:),'  ') - 1
            WRITE(56,'(3A,$)') '"X, ',UNITS(1)(1:NCH),'" '
            NCH = INDEX(UNITS(2)(1:),'  ') - 1
            WRITE(56,'(3A,$)') '"Y, ',UNITS(2)(1:NCH),'" '
            NCH = INDEX(UNITS(3)(1:),'  ') - 1
            WRITE(56,'(3A,$)') '"Z, ',UNITS(3)(1:NCH),'" '
            NCH = INDEX(UNITS(4)(1:),'  ') - 1
            WRITE(56,'(3A,$)') '"Fracture/Fault Triangle Area, ',
     &        UNITS(4)(1:NCH),'" '
            NCH = INDEX(UNITS(5)(1:),'  ') - 1
            WRITE(56,'(3A,$)') '"Fracture/Fault Triangle ' // 
     &        'Aperture/Width, ',UNITS(5)(1:NCH),'" '
            DO NV = 1,NVPLOT
              NCH = INDEX(VARNM(NV)(1:),'  ') - 1
              WRITE(56,'(3A,$)') '"',VARNM(NV)(1:NCH),'" '
            ENDDO
            WRITE(56,'(A)') ''
            NCH1 = INDEX( TIME(IUNTM)(1:),'  ' ) - 1
            NCH2 = INDEX( UNTM(1:),'  ' ) - 1
            WRITE(56,'(5A,$)') 'ZONE T = "',TIME(IUNTM)(1:NCH1),', ',
     &        UNTM(1:NCH2),'"'         
            WRITE(56,'(2A,$)') ', STRANDID = 1, SOLUTIONTIME = ',
     &        TIME(IUNTM)(1:NCH1)
            WRITE(56,'(A,I9,$)') ', NODES = ',NFT*3
            WRITE(56,'(A,I9,$)') ', ELEMENTS = ',NAFT
            WRITE(56,'(A)') ', DATAPACKING = BLOCK, ' // 
     &          'ZONETYPE = FETRIANGLE'
            WRITE(56,'(A,$)') 'VARLOCATION = (['
            DO N = 1,NVPLOT+2
              IF( N.GT.1 ) WRITE(56,'(A,$)') ','
              M = N + 3
              NCT = ICOUNT(M)
              WRITE(FORM1(3:3),'(I1)') NCT
              WRITE(56,FORM1) M
            ENDDO
            WRITE(56,'(A)') '] = CELLCENTERED)'
            DO N = 1,NFT
              WRITE(56,'(3(1PE14.7,1X))') (XE(NV,N),NV=1,3)
            ENDDO
            DO N = 1,NFT
              WRITE(56,'(3(1PE14.7,1X))') (YE(NV,N),NV=1,3)
            ENDDO
            DO N = 1,NFT
              WRITE(56,'(3(1PE14.7,1X))') (ZE(NV,N),NV=1,3)
            ENDDO
            NC = 0
            DO N = 1,NFT
              IF( NMAP(N).EQ.0 ) CYCLE
              NC = NC + 1
              VARX(NC) = AREA(N)
              IF( NC.EQ.10 .OR. N.EQ.NFT ) THEN
                WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                NC = 0
              ENDIF
            ENDDO
            IF( NC.NE.0 ) THEN
              WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
              NC = 0
            ENDIF
            DO N = 1,NFT
              IF( NMAP(N).EQ.0 ) CYCLE
              NC = NC + 1
              VARX(NC) = APER(N)
              IF( NC.EQ.10 .OR. N.EQ.NFT ) THEN
                WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                NC = 0
              ENDIF
            ENDDO
            IF( NC.NE.0 ) THEN
              WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
              NC = 0
            ENDIF
            DO NV = 1,NVPLOT
              NC = 0
              DO N = 1,NFT
                IF( NMAP(N).EQ.0 ) CYCLE
                NC = NC + 1
                VARX(NC) = VAR(N,NV)
                IF( NC.EQ.10 .OR. N.EQ.NFT ) THEN
                  WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                  NC = 0
                ENDIF
              ENDDO
              IF( NC.NE.0 ) THEN
                WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                NC = 0
              ENDIF
            ENDDO
            DO N = 1,NFT
              IF( NMAP(N).EQ.0 ) CYCLE
              WRITE(56,'(8(I9,1X))') (NCMAP(NV,N),NV=1,3)
            ENDDO
          ENDIF       
        ELSE
!
!---      Contour plotting type  ---
!
          IF( ITYPE.EQ.0 ) THEN
            NCH1 = INDEX( TIME(IUNTM)(1:),'  ' ) - 1
            NCH2 = INDEX( UNTM(1:),'  ' ) - 1
            WRITE(56,'(5A,$)') 'ZONE T = "',TIME(IUNTM)(1:NCH1),', ',
     &        UNTM(1:NCH2),'"'         
            WRITE(56,'(2A,$)') ', STRANDID = 1, SOLUTIONTIME = ',
     &        TIME(IUNTM)(1:NCH1)
            WRITE(56,'(A,I9,$)') ', NODES = ',NFT*3
            WRITE(56,'(A,I9,$)') ', ELEMENTS = ',NAFT
            WRITE(56,'(A,$)') ', DATAPACKING = BLOCK, ' // 
     &        'ZONETYPE = FETRIANGLE,'
            WRITE(56,'(A)')  ' VARSHARELIST = ([1,2,3]=1), ' // 
     &        'CONNECTIVITYSHAREZONE = 1'
            WRITE(56,'(A,$)') 'VARLOCATION = (['
            DO N = 1,NVPLOT+2
              IF( N.GT.1 ) WRITE(56,'(A,$)') ','
              M = N + 3
              NCT = ICOUNT(M)
              WRITE(FORM1(3:3),'(I1)') NCT
              WRITE(56,FORM1) M
            ENDDO
            WRITE(56,'(A)') '] = CELLCENTERED)'
            NC = 0
            DO N = 1,NFT
              IF( NMAP(N).EQ.0 ) CYCLE
              NC = NC + 1
              VARX(NC) = AREA(N)
              IF( NC.EQ.10 .OR. N.EQ.NFT ) THEN
                WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                NC = 0
              ENDIF
            ENDDO
            DO N = 1,NFT
              IF( NMAP(N).EQ.0 ) CYCLE
              NC = NC + 1
              VARX(NC) = APER(N)
              IF( NC.EQ.10 .OR. N.EQ.NFT ) THEN
                WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                NC = 0
              ENDIF
            ENDDO
            DO NV = 1,NVPLOT
              NC = 0
              DO N = 1,NFT
                IF( NMAP(N).EQ.0 ) CYCLE
                NC = NC + 1
                VARX(NC) = VAR(N,NV)
                IF( NC.EQ.10 .OR. N.EQ.NFT ) THEN
                  WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                  NC = 0
                ENDIF
              ENDDO
            ENDDO
          ENDIF       
        ENDIF
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
        DEALLOCATE( AREA,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: AREA'
          STOP
        ENDIF
        DEALLOCATE( APER,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: APER'
          STOP
        ENDIF
        DEALLOCATE( NMAP,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: NMAP'
          STOP
        ENDIF
!
!---    Close the STOMP plot file  ---
!
        CLOSE( UNIT=55 )
      ENDDO
!
!---  Close the Tecplot data file  ---
!
      CLOSE( UNIT=56 )
      STOP
!
!---  End of PLOT_FRCTO_TEC program  ---
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
!     Written by MD White, Battelle, PNL, March, 1993.
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

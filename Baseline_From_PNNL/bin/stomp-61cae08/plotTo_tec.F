!----------------------Program-----------------------------------------!
!
      PROGRAM PLOT_TO_TEC
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
!     This utility routine converts one or more plot files to tecplot
!     data format for plotting.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 January 2021.
!     plot_to_tec.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8, DIMENSION(:), ALLOCATABLE :: XP,YP,ZP,VOL
      INTEGER, DIMENSION(:), ALLOCATABLE :: NMAP
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NCMAP
      CHARACTER(32), DIMENSION(4) :: UNITS
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
      DO M = 1,4
        UNITS(M) = ''
      ENDDO
      CNTFILE = 'connect'
      NCLA = IARGC()
      IFILE = 0
      IF( NCLA.EQ.0 ) THEN
        WRITE(6,'(/,A,/)')'Enter "line" or "contour" or "help": '
        READ(5,*) PLTTYPE
        CALL LCASE( PLTTYPE )
        CDUM = PLTTYPE
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts outputs in the ' //
     &      'data section of one or more STOMP plot file(s) to a '
          WRITE(6,'(A)') 'tecplot data formatted file with an ".dat"' //
     &      ' extension. The script takes three or more argument(s): '
          WRITE(6,'(A)') '  Plotting Option: "line" or "contour"'
          WRITE(6,'(A)') '  Output File Name: (e.g. plots.dat)'
          WRITE(6,'(A)') '  Plot File(s): (e.g. plot.* or plot.0001)'
          WRITE(6,'(A)') 'The "line" plotting option is only ' //
     &       'allowed for one-dimensional domains, and is used '
          WRITE(6,'(A)') 'for making line plots. The "contour"' //
     &       'plotting option is for one-, two-, and three-dimensional '
          WRITE(6,'(A)') 'domains, and is used to make contour ' //
     &      'plots. The output file is a new file created by this '
          WRITE(6,'(A)') 'script, which will become the input ' //
     &      'file to tecplot. The data in each plot file, will be '
          WRITE(6,'(A,/)') 'stored in the output file, as a ' //
     &      'separate zone, titled by the time stamp in the plot file.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'plotTo_tec.exe line plots.dat plot.*'
          WRITE(6,'(A)') 'plotTo_tec.exe contour plots.dat plot.*'
          WRITE(6,'(A)') 'plotTo_tec.exe line pressure.dat plot.002'
          WRITE(6,'(A,/)') 'plotTo_tec.exe contour run1.dat ' // 
     &      'plot.01 plot.02 plot.10'
          STOP
        ENDIF
        WRITE(6,'(/,A,/)')'Enter the output filename: '
        READ(5,*) OUTFILE
        WRITE(6,'(/,A,/)')'Enter a single plot filename ' // 
     &    '(e.g., plot.001): '
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
     &      'data section of one or more STOMP plot file(s) to a '
          WRITE(6,'(A)') 'tecplot data formatted file with an ".dat"' //
     &      ' extension. The script takes three or more argument(s): '
          WRITE(6,'(A)') '  Plotting Option: "line" or "contour"'
          WRITE(6,'(A)') '  Output File Name: (e.g. plots.dat)'
          WRITE(6,'(A)') '  Plot File(s): (e.g. plot.* or plot.0001)'
          WRITE(6,'(A)') 'The "line" plotting option is only ' //
     &       'allowed for one-dimensional domains, and is used '
          WRITE(6,'(A)') 'for making line plots. The "contour"' //
     &       'plotting option is for one-, two-, and three-dimensional '
          WRITE(6,'(A)') 'domains, and is used to make contour ' //
     &      'plots. The output file is a new file created by this '
          WRITE(6,'(A)') 'script, which will become the input ' //
     &      'file to tecplot. The data in each plot file, will be '
          WRITE(6,'(A,/)') 'stored in the output file, as a ' //
     &      'separate zone, titled by the time stamp in the plot file.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'plotTo_tec.exe line plots.dat plot.*'
          WRITE(6,'(A)') 'plotTo_tec.exe contour plots.dat plot.*'
          WRITE(6,'(A)') 'plotTo_tec.exe line pressure.dat plot.002'
          WRITE(6,'(A,/)') 'plotTo_tec.exe contour run1.dat ' // 
     &      'plot.01 plot.02 plot.10'
          STOP
        ENDIF
      ENDIF
!
!---  Check plotting type  ---
!
      IF( INDEX(PLTTYPE(1:),'line').NE.1 .AND. 
     &  INDEX(PLTTYPE(1:),'contour').NE.1 ) THEN
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
        WRITE(6,'(A)') 'ERROR: Nonexistent "connect" file.'
        STOP
      ELSE
        OPEN(UNIT=45, FILE=CNTFILE, STATUS='OLD', FORM='FORMATTED')
      ENDIF
!
!---  Write banner to screen  ---
!       
      WRITE(6,'(A,//)')' Welcome to ...'
      WRITE(6,'(A,/)') '                       PLOT_TO_TEC'
      WRITE(6,'(A)')   ' An output file converter routine for ' //
     &  'the STOMP simulator.  This'
      WRITE(6,'(A)')   ' routine, developed by the Pacific ' //
     &  'Northwest National Laboratory,'
      WRITE(6,'(A)')   ' reads one or more STOMP plot file(s) ' //
     &  ' and converts the data to one '
      WRITE(6,'(A)')   ' tecplot data file with an extension ' //
     &  'of ".dat". This file contains the '
      WRITE(6,'(A)')   'user-requested data at every grid cell' //
     &  'in tecplot format for visualization with tecplot. '
      WRITE(6,'(A)')   'One-dimensional domains can be ' //
     &  'visualized as line plots or contour plots. Two- and '
      WRITE(6,'(A,/)')   'three-dimensional domains are only' //
     &  'visualized as contour plots. '
      WRITE(6,'(A)')   ' For support:  Tel: 509.372.6070'
      WRITE(6,'(A,/)') '               E-mail:  mark.white@pnnl.gov'
!
!---  Loop over STOMP plot files ---
!
      DO NA = 3,NCLA
        IF( IFILE.EQ.0 ) THEN
          CALL GETARG( NA,PLTFILE )
        ENDIF
!
!---    No plot file, exit program with error message  ---
!
        INQUIRE( FILE=PLTFILE, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          NCH = INDEX(PLTFILE(1:),'  ')-1
          WRITE(6,'(A)') 'ERROR: Nonexistent plot file: ',PLTFILE(1:NCH)
          STOP
        ELSE
          OPEN(UNIT=55, FILE=PLTFILE, STATUS='OLD', FORM='FORMATTED')
        ENDIF
        NCH = INDEX(PLTFILE(1:),'  ')-1
        WRITE(6,'(2A)') 'Processing plot file: ',PLTFILE(1:NCH)
!
!---    Search for coordinate dimensions  ---
!
        DO
          READ(55,'(A)', IOSTAT=ICHK) CHDUM
          IF( ICHK.GT.0 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &        'file read: ',PLTFILE(1:NCH)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &        'file read: ',PLTFILE(1:NCH)
            STOP
          ENDIF
          UNTMX = INDEX(CHDUM(1:),'Time Output Units =')
          TIMEX = INDEX(CHDUM(1:),'Time =')
          NVPTX = INDEX(CHDUM(1:),'Number of Plot Variables =')
          IFLDX = INDEX(CHDUM(1:),'X or R-Direction Nodes =')
          JFLDX = INDEX(CHDUM(1:),'Y or Theta-Direction Nodes =')
          KFLDX = INDEX(CHDUM(1:),'Z-Direction Nodes =')
          NACTX = INDEX(CHDUM(1:),'Number of Active Nodes =')
          NVTCX = INDEX(CHDUM(1:),'Number of Vertices =')
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
          IF( IFLDX.NE.0 ) THEN
            IFLDX = INDEX(CHDUM(1:),'=')
            I1 = IFLDX+1
            I2 = IFLDX+10
            READ(CHDUM(I1:I2),'(I10)') IFLD
          ENDIF
          IF( JFLDX.NE.0 ) THEN
            JFLDX = INDEX(CHDUM(1:),'=')
            J1 = JFLDX+1
            J2 = JFLDX+10
            READ(CHDUM(J1:J2),'(I10)') JFLD
          ENDIF
          IF( KFLDX.NE.0 ) THEN
            KFLDX = INDEX(CHDUM(1:),'=')
            K1 = KFLDX+1
            K2 = KFLDX+10
            READ(CHDUM(K1:K2),'(I10)') KFLD
          ENDIF
          IF( NACTX.NE.0 ) THEN
            NACTX = INDEX(CHDUM(1:),'=')
            N1 = NACTX+1
            N2 = NACTX+10
            READ(CHDUM(N1:N2),'(I10)') NACT
          ENDIF
          IF( NVTCX.NE.0 ) THEN
            NVTCX = INDEX(CHDUM(1:),'=')
            N1 = NVTCX+1
            N2 = NVTCX+10
            READ(CHDUM(N1:N2),'(I10)') NVTC
            EXIT
          ENDIF
        ENDDO  
        NFLD = IFLD*JFLD*KFLD
        IF( INDEX(PLTTYPE(1:),'line').NE.0 ) THEN
          IF( NFLD.EQ.IFLD ) THEN
            ITYPE = 1
          ELSEIF( NFLD.EQ.JFLD ) THEN
            ITYPE = 2
          ELSEIF( NFLD.EQ.KFLD ) THEN
            ITYPE = 3
          ELSE
            WRITE(6,'(2A)') 'Line plotting option declared for ' //
     &        'a two- or three-dimensional domain.'
            STOP
          ENDIF
        ELSE
          ITYPE = 0
        ENDIF
        IF( NA.GT.3 ) THEN
          IF( IFLD.NE.IFLD1 .OR. JFLD.NE.JFLD1 .OR. KFLD.NE.KFLD1 .OR.
     &      NACT.NE.NACT1 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'Plot file not compatible with previous ' //
     &        'plot files: ',PLTFILE(1:NCH)
            STOP
          ENDIF
        ELSE
          IFLD1 = IFLD
          JFLD1 = JFLD
          KFLD1 = KFLD
          NACT1 = NACT
        ENDIF     
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
        ALLOCATE( XE(1:NVTC,1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: XE'
          STOP
        ENDIF
        ALLOCATE( YE(1:NVTC,1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: YE'
          STOP
        ENDIF
        ALLOCATE( ZE(1:NVTC,1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: ZE'
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
        ALLOCATE( NMAP(1:NFLD),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: NMAP'
          STOP
        ENDIF
        IF( NA.EQ.3 ) THEN
          ALLOCATE( NCMAP(1:NVTC,1:NFLD),STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            WRITE(6,'(A)') 'Allocation Error: NCMAP'
            STOP
          ENDIF
          DO N = 1,NFLD
            DO NV = 1,NVTC
              NCMAP(NV,N) = 0
            ENDDO
          ENDDO
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
!---    Search for x vertices  ---
!
        IF( IFLD.GT.1 .OR. NVTC.EQ.8 ) THEN
          DO
            READ(55,'(A)', IOSTAT=ICHK) CHDUM
            IF( ICHK.GT.0 ) THEN
              NCH = INDEX( PLTFILE(1:),'  ' ) - 1
              WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &          'file read: X-Vertices: ',PLTFILE(1:NCH)
              STOP
            ELSEIF( ICHK.LT.0 ) THEN
              WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &          'file read: X-Vertices: ',PLTFILE(1:NCH)
              STOP
            ENDIF
            IXDNV = INDEX(CHDUM(1:),'X-Direction Nodal Vertices')
            IF( IXDNV.NE.0 ) THEN
              I1 = INDEX(CHDUM(1:),',')+1
              I2 = INDEX(CHDUM(1:),'   ')-1
              UNITS(1) = CHDUM(I1:I2)
              EXIT
            ENDIF
          ENDDO
          DO N = 1,NFLD
            READ(55,*) (XE(M,N),M=1,NVTC)
            XP(N) = 0.D+0
            DO M = 1,NVTC
              XP(N) = XP(N) + XE(M,N)
            ENDDO
            XP(N) = XP(N)/REAL(NVTC)
          ENDDO
        ENDIF
!
!---    Search for y vertices  ---
!
        IF( JFLD.GT.1 .OR. NVTC.EQ.8 ) THEN
          DO
            READ(55,'(A)', IOSTAT=ICHK) CHDUM
            IF( ICHK.GT.0 ) THEN
              NCH = INDEX( PLTFILE(1:),'  ' ) - 1
              WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &          'file read: Y-Vertices: ',PLTFILE(1:NCH)
              STOP
            ELSEIF( ICHK.LT.0 ) THEN
              WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &          'file read: Y-Vertices: ',PLTFILE(1:NCH)
              STOP
            ENDIF
            IXDNV = INDEX(CHDUM(1:),'Y-Direction Nodal Vertices')
            IF( IXDNV.NE.0 ) THEN
              I1 = INDEX(CHDUM(1:),',')+1
              I2 = INDEX(CHDUM(1:),'   ')-1
              UNITS(2) = CHDUM(I1:I2)
              EXIT
            ENDIF
          ENDDO
          DO N = 1,NFLD
            READ(55,*) (YE(M,N),M=1,NVTC)
            YP(N) = 0.D+0
            DO M = 1,NVTC
              YP(N) = YP(N) + YE(M,N)
            ENDDO
            YP(N) = YP(N)/REAL(NVTC)
          ENDDO
        ENDIF
!
!---    Search for z vertices  ---
!
        IF( KFLD.GT.1 .OR. NVTC.EQ.8 ) THEN
          DO
            READ(55,'(A)', IOSTAT=ICHK) CHDUM
            IF( ICHK.GT.0 ) THEN
              NCH = INDEX( PLTFILE(1:),'  ' ) - 1
              WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &          'file read: Z-Vertices: ',PLTFILE(1:NCH)
              STOP
            ELSEIF( ICHK.LT.0 ) THEN
              WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &          'file read: Z-Vertices: ',PLTFILE(1:NCH)
              STOP
            ENDIF
            IXDNV = INDEX(CHDUM(1:),'Z-Direction Nodal Vertices')
            IF( IXDNV.NE.0 ) THEN
              I1 = INDEX(CHDUM(1:),',')+1
              I2 = INDEX(CHDUM(1:),'   ')-1
              UNITS(3) = CHDUM(I1:I2)
              EXIT
            ENDIF
          ENDDO
          DO N = 1,NFLD
            READ(55,*) (ZE(M,N),M=1,NVTC)
            ZP(N) = 0.D+0
            DO M = 1,NVTC
              ZP(N) = ZP(N) + ZE(M,N)
            ENDDO
            ZP(N) = ZP(N)/REAL(NVTC)
          ENDDO
        ENDIF
!
!---    Search for node volume  ---
!
        DO
          READ(55,'(A)', IOSTAT=ICHK) CHDUM
          IF( ICHK.GT.0 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &        'file read: Node Volume: ',PLTFILE(1:NCH)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &        'file read: Node Volume: ',PLTFILE(1:NCH)
            STOP
          ENDIF
          IXDNV = INDEX(CHDUM(1:),'Node Volume')
          IF( IXDNV.NE.0 ) THEN
            I1 = INDEX(CHDUM(1:),',')+1
            I2 = INDEX(CHDUM(1:),'   ')-1
            UNITS(4) = CHDUM(I1:I2)
            EXIT
          ENDIF
        ENDDO
        DO N = 1,NFLD,10
          READ(55,*) (VOL(M),M=N,MIN(N+9,NFLD))
        ENDDO
!
!---    Search for node map  ---
!
        DO
          READ(55,'(A)', IOSTAT=ICHK) CHDUM
          IF( ICHK.GT.0 ) THEN
            NCH = INDEX( PLTFILE(1:),'  ' ) - 1
            WRITE(6,'(2A)') 'End of file encountered during plot ' //
     &        'file read: Node Map: ',PLTFILE(1:NCH)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &        'file read: Node Map: ',PLTFILE(1:NCH)
            STOP
          ENDIF
          IXDNV = INDEX(CHDUM(1:),'Node Map')
          IF( IXDNV.NE.0 ) EXIT
        ENDDO
        DO N = 1,NFLD,10
          READ(55,*) (NMAP(M),M=N,MIN(N+9,NFLD))
        ENDDO
!
!---    Read connection map  ---
!
        IF( NA.EQ.3 ) THEN
          DO N = 1,NFLD
            IF( NMAP(N).EQ.0 ) CYCLE
            READ(45,*) (NCMAP(NV,N),NV=1,NVTC)
          ENDDO
!
!---      Close the connect file  ---
!
          CLOSE( UNIT=45 )
        ENDIF
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
          DO N = 1,NFLD,10
            READ(55,*) (VAR(M,NV),M=N,MIN(N+9,NFLD))
          ENDDO
        ENDDO
!
!---    Write file version and identifier to the tecplot data file  ---
!
        IF( NA.EQ.3 ) THEN
!
!---      Contour plotting type  ---
!
          IF( ITYPE.EQ.0 ) THEN
            WRITE(56,'(A)') 'TITLE = "STOMP Plot File Visualization"'
            WRITE(56,'(A,$)') 'VARIABLES = '
            NCH = INDEX(UNITS(1)(1:),'  ') - 1
            IF( IFLD.GT.1 .OR. NVTC.EQ.8 ) WRITE(56,'(3A,$)') '"X, ',
     &        UNITS(1)(1:NCH),'" '
            NCH = INDEX(UNITS(2)(1:),'  ') - 1
            IF( JFLD.GT.1 .OR. NVTC.EQ.8 ) WRITE(56,'(3A,$)') '"Y, ',
     &        UNITS(2)(1:NCH),'" '
            NCH = INDEX(UNITS(3)(1:),'  ') - 1
            IF( KFLD.GT.1 .OR. NVTC.EQ.8 ) WRITE(56,'(3A,$)') '"Z, ',
     &        UNITS(3)(1:NCH),'" '
            NCH = INDEX(UNITS(4)(1:),'  ') - 1
            WRITE(56,'(3A,$)') '"Volume, ',UNITS(4)(1:NCH),'" '
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
            WRITE(56,'(A,I9,$)') ', NODES = ',NFLD*NVTC
            WRITE(56,'(A,I9,$)') ', ELEMENTS = ',NACT
            IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1) 
     &        .OR. (NVTC.EQ.8) ) THEN
              WRITE(56,'(A)') ', DATAPACKING = BLOCK, ' // 
     &          'ZONETYPE = FEBRICK'
            ELSE
              WRITE(56,'(A)') ', DATAPACKING = BLOCK, ' // 
     &          'ZONETYPE = FEQUADRILATERAL'
            ENDIF
            WRITE(56,'(A,$)') 'VARLOCATION = (['
            IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1)
     &         .OR. (NVTC.EQ.8) ) THEN
              DO N = 1,NVPLOT+1
                IF( N.GT.1 ) WRITE(56,'(A,$)') ','
                M = N + 3
                NCT = ICOUNT(M)
                WRITE(FORM1(3:3),'(I1)') NCT
                WRITE(56,FORM1) M
              ENDDO
            ELSE
              DO N = 1,NVPLOT+1
                IF( N.GT.1 ) WRITE(56,'(A,$)') ','
                M = N + 2
                NCT = ICOUNT(M)
                WRITE(FORM1(3:3),'(I1)') NCT
                WRITE(56,FORM1) M
              ENDDO
            ENDIF
            WRITE(56,'(A)') '] = CELLCENTERED)'
            IF( IFLD.GT.1 .OR. NVTC.EQ.8 ) THEN
              DO N = 1,NFLD
                WRITE(56,'(8(1PE14.7,1X))') (XE(NV,N),NV=1,NVTC)
              ENDDO
            ENDIF
            IF( JFLD.GT.1 .OR. NVTC.EQ.8 ) THEN
              DO N = 1,NFLD
                WRITE(56,'(8(1PE14.7,1X))') (YE(NV,N),NV=1,NVTC)
              ENDDO
            ENDIF
            IF( KFLD.GT.1 .OR. NVTC.EQ.8 ) THEN
              DO N = 1,NFLD
                WRITE(56,'(8(1PE14.7,1X))') (ZE(NV,N),NV=1,NVTC)
              ENDDO
            ENDIF
            NC = 0
            DO N = 1,NFLD
              IF( NMAP(N).EQ.0 ) CYCLE
              NC = NC + 1
              VARX(NC) = VOL(N)
              IF( NC.EQ.10 .OR. N.EQ.NFLD ) THEN
                WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                NC = 0
              ENDIF
            ENDDO
            IF( NC.NE.0 ) WRITE(56,'(10(1PE14.7,1X))')(VARX(M),M=1,NC)
            DO NV = 1,NVPLOT
              NC = 0
              DO N = 1,NFLD
                IF( NMAP(N).EQ.0 ) CYCLE
                NC = NC + 1
                VARX(NC) = VAR(N,NV)
                IF( NC.EQ.10 .OR. N.EQ.NFLD ) THEN
                  WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                  NC = 0
                ENDIF
              ENDDO
              IF( NC.NE.0 ) WRITE(56,'(10(1PE14.7,1X))')(VARX(M),M=1,NC)
            ENDDO
            DO N = 1,NFLD
              IF( NMAP(N).EQ.0 ) CYCLE
              WRITE(56,'(8(I9,1X))') (NCMAP(NV,N),NV=1,NVTC)
            ENDDO
!
!---      Line plotting type  ---
!
          ELSE
            WRITE(56,'(A)') 'TITLE = "STOMP Plot File Visualization"'
            WRITE(56,'(A,$)') 'VARIABLES = '
            IF( ITYPE.EQ.1 ) THEN
              NCH = INDEX(UNITS(1)(1:),'  ') - 1
              WRITE(56,'(3A,$)') '"X, ',UNITS(1)(1:NCH),'" '
            ELSEIF( ITYPE.EQ.2 ) THEN
              NCH = INDEX(UNITS(2)(1:),'  ') - 1
              WRITE(56,'(3A,$)') '"Y, ',UNITS(2)(1:NCH),'" '
            ELSEIF( ITYPE.EQ.3 ) THEN
              NCH = INDEX(UNITS(3)(1:),'  ') - 1
              WRITE(56,'(3A,$)') '"Z, ',UNITS(3)(1:NCH),'" '
            ENDIF
            NCH = INDEX(UNITS(4)(1:),'  ') - 1
            WRITE(56,'(3A,$)') '"Volume, ',UNITS(4)(1:NCH),'" '
            DO NV = 1,NVPLOT
              NCH = INDEX(VARNM(NV)(1:),'  ') - 1
              WRITE(56,'(3A,$)') '"',VARNM(NV)(1:NCH),'" '
            ENDDO
            WRITE(56,'(A)') ''
            NCH1 = INDEX( TIME(IUNTM)(1:),'  ' ) - 1
            NCH2 = INDEX( UNTM(1:),'  ' ) - 1
            WRITE(56,'(5A)') 'ZONE T = "',TIME(IUNTM)(1:NCH1),', ',
     &        UNTM(1:NCH2),'", F=POINT'  
            DO N = 1,NFLD
              IF( NMAP(N).EQ.0 ) CYCLE
              DVAR = 0.D+0
              IF( ITYPE.EQ.1 ) THEN
                DO NV = 1,8
                  DVAR = DVAR + XE(NV,N)
                ENDDO
              ELSEIF( ITYPE.EQ.2 ) THEN
                DO NV = 1,8
                  DVAR = DVAR + YE(NV,N)
                ENDDO
              ELSEIF( ITYPE.EQ.3 ) THEN
                DO NV = 1,8
                  DVAR = DVAR + ZE(NV,N)
                ENDDO
              ENDIF
              DVAR = DVAR/8.D+0
              WRITE(56,'(1PE15.7,1X,$)') DVAR
              WRITE(56,'(1PE15.7,1X,$)') VOL(N)
              DO NV = 1,NVPLOT
                IF( NV.EQ.NVPLOT ) THEN
                  WRITE(56,'(1PE15.7)') VAR(N,NV)
                ELSE
                  WRITE(56,'(1PE15.7,1X,$)') VAR(N,NV)
                ENDIF
              ENDDO
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
            WRITE(56,'(A,I9,$)') ', NODES = ',NFLD*NVTC
            WRITE(56,'(A,I9,$)') ', ELEMENTS = ',NACT
            IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1)
     &         .OR. (NVTC.EQ.8) ) THEN
              WRITE(56,'(A,$)') ', DATAPACKING = BLOCK, ' // 
     &          'ZONETYPE = FEBRICK,'
            ELSE
              WRITE(56,'(A,$)') ', DATAPACKING = BLOCK, ' // 
     &          'ZONETYPE = FEQUADRILATERAL,'
            ENDIF
             IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1)
     &          .OR. (NVTC.EQ.8) ) THEN
              WRITE(56,'(A)')  ' VARSHARELIST = ([1,2,3]=1), ' // 
     &          'CONNECTIVITYSHAREZONE = 1'
            ELSE
              WRITE(56,'(A)')  ' VARSHARELIST = ([1,2]=1), ' // 
     &          'CONNECTIVITYSHAREZONE = 1'
            ENDIF
            WRITE(56,'(A,$)') 'VARLOCATION = (['
            IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1)
     &         .OR. (NVTC.EQ.8) ) THEN
              DO N = 1,NVPLOT+1
                IF( N.GT.1 ) WRITE(56,'(A,$)') ','
                M = N + 3
                NCT = ICOUNT(M)
                WRITE(FORM1(3:3),'(I1)') NCT
                WRITE(56,FORM1) M
              ENDDO
            ELSE
              DO N = 1,NVPLOT+1
                IF( N.GT.1 ) WRITE(56,'(A,$)') ','
                M = N + 2
                NCT = ICOUNT(M)
                WRITE(FORM1(3:3),'(I1)') NCT
                WRITE(56,FORM1) M
              ENDDO
            ENDIF
            WRITE(56,'(A)') '] = CELLCENTERED)'
            NC = 0
            DO N = 1,NFLD
              IF( NMAP(N).EQ.0 ) CYCLE
              NC = NC + 1
              VARX(NC) = VOL(N)
              IF( NC.EQ.10 .OR. N.EQ.NFLD ) THEN
                WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                NC = 0
              ENDIF
            ENDDO
            IF( NC.NE.0 ) WRITE(56,'(10(1PE14.7,1X))')(VARX(M),M=1,NC)
            DO NV = 1,NVPLOT
              NC = 0
              DO N = 1,NFLD
                IF( NMAP(N).EQ.0 ) CYCLE
                NC = NC + 1
                VARX(NC) = VAR(N,NV)
                IF( NC.EQ.10 .OR. N.EQ.NFLD ) THEN
                  WRITE(56,'(10(1PE14.7,1X))') (VARX(M),M=1,NC)
                  NC = 0
                ENDIF
              ENDDO
              IF( NC.NE.0 ) WRITE(56,'(10(1PE14.7,1X))')(VARX(M),M=1,NC)
            ENDDO
!
!---      Line plotting type  ---
!
          ELSE
            WRITE(56,'(A,$)') 'VARIABLES = '
            IF( ITYPE.EQ.1 ) THEN
              NCH = INDEX(UNITS(1)(1:),'  ') - 1
              WRITE(56,'(3A,$)') '"X, ',UNITS(1)(1:NCH),'" '
            ELSEIF( ITYPE.EQ.2 ) THEN
              NCH = INDEX(UNITS(2)(1:),'  ') - 1
              WRITE(56,'(3A,$)') '"Y, ',UNITS(2)(1:NCH),'" '
            ELSEIF( ITYPE.EQ.3 ) THEN
              NCH = INDEX(UNITS(3)(1:),'  ') - 1
              WRITE(56,'(3A,$)') '"Z, ',UNITS(3)(1:NCH),'" '
            ENDIF
            NCH = INDEX(UNITS(4)(1:),'  ') - 1
            WRITE(56,'(3A,$)') '"Volume, ',UNITS(4)(1:NCH),'" '
            DO NV = 1,NVPLOT
              NCH = INDEX(VARNM(NV)(1:),'  ') - 1
              WRITE(56,'(3A,$)') '"',VARNM(NV)(1:NCH),'" '
            ENDDO
            WRITE(56,'(A)') ''
            NCH1 = INDEX( TIME(IUNTM)(1:),'  ' ) - 1
            NCH2 = INDEX( UNTM(1:),'  ' ) - 1
            WRITE(56,'(5A)') 'ZONE T = "',TIME(IUNTM)(1:NCH1),', ',
     &        UNTM(1:NCH2),'", F=POINT'  
            DO N = 1,NFLD
              IF( NMAP(N).EQ.0 ) CYCLE
              DVAR = 0.D+0
              IF( ITYPE.EQ.1 ) THEN
                DO NV = 1,8
                  DVAR = DVAR + XE(NV,N)
                ENDDO
              ELSEIF( ITYPE.EQ.2 ) THEN
                DO NV = 1,8
                  DVAR = DVAR + YE(NV,N)
                ENDDO
              ELSEIF( ITYPE.EQ.3 ) THEN
                DO NV = 1,8
                  DVAR = DVAR + ZE(NV,N)
                ENDDO
              ENDIF
              DVAR = DVAR/8.D+0
              WRITE(56,'(1PE15.7,1X,$)') DVAR
              WRITE(56,'(1PE15.7,1X,$)') VOL(N)
              DO NV = 1,NVPLOT
                IF( NV.EQ.NVPLOT ) THEN
                  WRITE(56,'(1PE15.7)') VAR(N,NV)
                ELSE
                  WRITE(56,'(1PE15.7,1X,$)') VAR(N,NV)
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
        DEALLOCATE( VOL,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Deallocation Error: VOL'
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
!---  Close the tecplot data file  ---
!
      CLOSE( UNIT=56 )
      STOP
!
!---  End of PLOT_TO_VTK program  ---
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

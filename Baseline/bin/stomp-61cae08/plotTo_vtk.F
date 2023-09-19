!----------------------Program-----------------------------------------!
!
      PROGRAM PLOT_TO_VTK
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
!     This utility routine converts one or more plot files to VTK format
!     for plotting.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 May 2022.
!     plotTo_vtk.F
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
      INTEGER, DIMENSION(10) :: IVARX
      CHARACTER(32), DIMENSION(4) :: UNITS
      CHARACTER(128), DIMENSION(:), ALLOCATABLE :: VARNM
      CHARACTER(128) :: CDUM,PLTFILE,OUTFILE,CNTFILE
      CHARACTER(4) :: FORM1
      CHARACTER(512) :: CHDUM
      LOGICAL FCHK
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
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
        WRITE(6,'(/,A,/)')'Enter a single plot filename: '
        READ(5,*) PLTFILE
        CALL LCASE( PLTFILE )
        CDUM = PLTFILE
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts outputs in the ' //
     &      'data section of one or more STOMP plot file(s) to a '
          WRITE(6,'(A)') 'VTK (Visualization Toolkit) ' //
     &      'formatted file. The script takes one or more argument(s), '
          WRITE(6,'(A)') 'the names of plot files to convert. ' //
     &      'For each plot file in the argument list, a new file '
          WRITE(6,'(A)') 'will be created that contains a VTK ' //
     &      'legacy serial formatted file that contains the five '
          WRITE(6,'(A)') 'basic parts: 1) file version and ' //
     &      'identifer, 2) header, 3) file format, 4) dataset '
          WRITE(6,'(A)') 'strucuture, and 5) dataset attributes. ' //
     &      'The output file names will be plot-1234.vtk, where 1234 '
          WRITE(6,'(A)') 'was the time step extension of the plot ' //
     &      'files in the argument list. If an argument is  missing '
          WRITE(6,'(A,/)') 'the user will be prompted for a single ' //
     &      'plot file. '
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'plotTo_vtk.exe plot.*'
          WRITE(6,'(A)') 'plotTo_vtk.exe plot.002'
          WRITE(6,'(A,/)') 'plotTo_vtk.exe plot.01 plot.02 plot.10'
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
            WRITE(6,'(/,A)') 'This script converts outputs in the ' //
     &        'data section of one or more STOMP plot file(s) to a '
            WRITE(6,'(A)') 'VTK (Visualization Toolkit) formatted ' //
     &        'file. The script takes one or more argument(s), '
            WRITE(6,'(A)') 'the names of plot files to convert. ' //
     &        'For each plot file in the argument list, a new file '
            WRITE(6,'(A)') 'will be created that contains a VTK ' //
     &        'legacy serial formatted file that contains the five '
            WRITE(6,'(A)') 'basic parts: 1) file version and ' //
     &        'identifer, 2) header, 3) file format, 4) dataset '
            WRITE(6,'(A)') 'strucuture, and 5) dataset attributes. ' //
     &        'The output file names will be plot-1234.vtk, where 1234 '
            WRITE(6,'(A)') 'was the time step extension of the plot ' //
     &        'files in the argument list. If an argument is  missing '
            WRITE(6,'(A,/)') 'the user will be prompted for a ' //
     &        'single plot file. '
            WRITE(6,'(A,/)') 'Example command line entries'
            WRITE(6,'(A)') 'plotTo_vtk.exe plot.*'
            WRITE(6,'(A)') 'plotTo_vtk.exe plot.002'
            WRITE(6,'(A,/)') 'plotTo_vtk.exe plot.01 plot.02 plot.10'
            STOP
          ELSE
            PLTFILE = CDUM
          ENDIF
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
!
!---    No connect file, exit program with error message or open
!       connect file  ---
!
        INQUIRE( FILE=CNTFILE, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          WRITE(6,'(A)') 'ERROR: Nonexistent "connect" file.'
          STOP
        ELSE
          OPEN(UNIT=45, FILE=CNTFILE, STATUS='OLD', FORM='FORMATTED')
        ENDIF
!
!---    Write banner to screen  ---
!       
        IF( NA.EQ.1 ) THEN
          WRITE(6,'(A,//)')' Welcome to ...'
          WRITE(6,'(A,/)') '                       PLOT_TO_VTK'
          WRITE(6,'(A)')   ' An output file converter routine for ' //
     &      'the STOMP simulator.  This'
          WRITE(6,'(A)')   ' routine, developed by the Pacific ' //
     &      'Northwest National Laboratory,'
          WRITE(6,'(A)')   ' reads one or more STOMP plot file(s) ' //
     &      ' and converts the data to one '
          WRITE(6,'(A)')   ' or more plot files with an extension ' //
     &      'of ".vtk". These files contain '
          WRITE(6,'(A)')   'columns of user-requested data at ' //
     &      'every grid cell, starting with the'
          WRITE(6,'(A,/)')   'x, y, and z centroid coordinates and ' //
     &      'volume of the grid cell.'
          WRITE(6,'(A)')   ' For support:  Tel: 509.372.6070'
          WRITE(6,'(A,/)') '               E-mail:  mark.white@pnnl.gov'
        ENDIF
!
!---    Open plot-1234.vtk file  ---
!
        OUTFILE = PLTFILE
        NDOT = INDEX(OUTFILE(1:),'.')
        OUTFILE(NDOT:NDOT) = '-'
        NCH = INDEX(OUTFILE(1:),'  ')
        OUTFILE(NCH:NCH+3) = '.vtk'
        OUTFILE = TRIM(ADJUSTL(OUTFILE))
        WRITE(6,'(2A)') 'Creating file: ',TRIM(OUTFILE)
        OPEN(UNIT=56, FILE=OUTFILE)
        CLOSE(UNIT=56,STATUS='DELETE')
        OPEN(UNIT=56, FILE=OUTFILE, STATUS='NEW', FORM='FORMATTED')
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
          NVPTX = INDEX(CHDUM(1:),'Number of Plot Variables =')
          IFLDX = INDEX(CHDUM(1:),'X or R-Direction Nodes =')
          JFLDX = INDEX(CHDUM(1:),'Y or Theta-Direction Nodes =')
          KFLDX = INDEX(CHDUM(1:),'Z-Direction Nodes =')
          NACTX = INDEX(CHDUM(1:),'Number of Active Nodes =')
          NVTCX = INDEX(CHDUM(1:),'Number of Vertices =')
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
        IF( NA.GT.1 ) THEN
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
        IF( NA.EQ.1 ) THEN
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
            I1 = INDEX(CHDUM(1:),',')+2
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
        IF( NA.EQ.1 ) THEN
          DO N = 1,NFLD
            IF( NMAP(N).EQ.0 ) CYCLE
            READ(45,*) (NCMAP(NV,N),NV=1,NVTC)
            DO NV = 1,NVTC
              NCMAP(NV,N) = NCMAP(NV,N)-1
            ENDDO
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
     &          'file read: Node Map: ',PLTFILE(1:NCH)
              STOP
            ELSEIF( ICHK.LT.0 ) THEN
              WRITE(6,'(2A)') 'Read error encountered during plot ' //
     &          'file read: Node Map: ',PLTFILE(1:NCH)
              STOP
            ENDIF
            IF( (ICHAR(CHDUM(1:1)).GE.97 .AND. ICHAR(CHDUM(1:1)).LE.122)
     &        .OR. (ICHAR(CHDUM(1:1)).GE.65 .AND.
     &        ICHAR(CHDUM(1:1)).LE.90) ) THEN
              NCH = INDEX(CHDUM(1:),'   ')-1
              VARNM(NV) = CHDUM(1:NCH)
              DO NC = 1,NCH
                IF( INDEX(VARNM(NV)(NC:NC),' ').NE.0 ) THEN
                  VARNM(NV)(NC:NC) = '_'
                ENDIF
                IF( INDEX(VARNM(NV)(NC:NC+1),', ').NE.0 ) THEN
                  VARNM(NV)(NC:NC+1) = '_('
                  VARNM(NV)(NCH+1:NCH+1) = ')'
                ENDIF
              ENDDO
              EXIT
            ENDIF
          ENDDO
          DO N = 1,NFLD,10
            READ(55,*) (VAR(M,NV),M=N,MIN(N+9,NFLD))
          ENDDO
        ENDDO
!
!---    Write header line to plot-1234.vtk file  ---
!
        WRITE(56,'(A)') '# vtk DataFile Version 3.0'
        WRITE(56,'(A)') 'STOMP Visualization'
        WRITE(56,'(A)') 'ASCII'
        WRITE(56,'(A)') 'DATASET UNSTRUCTURED_GRID'
        NPOINTS = NFLD*NVTC
!
!---    VTK_HEXAHEDRON  ---
!
        IF( (IFLD.GT.1 .AND. JFLD.GT.1 .AND. KFLD.GT.1) 
     &    .OR. (NVTC.EQ.8) ) THEN
!
!---      Write x,y,z vertices to plot-1234.vtk file  ---
!
          WRITE(56,'(/,A,I9,A)') 'POINTS ',NPOINTS,' float'
          DO N = 1,NFLD
            DO M = 1,NVTC
              WRITE(56,'(3(1PE14.7,1X))') XE(M,N),YE(M,N),ZE(M,N)
            ENDDO
          ENDDO
!
!---      Write vertice map for active nodes to plot-1234.vtk file  ---
!
          WRITE(56,'(/,A,I9,1X,I9)') 'CELLS ',NACT,(NVTC+1)*NACT
          DO N = 1,NFLD
            IF( NMAP(N).EQ.0 ) CYCLE
            WRITE(56,'(5(I9,1X))') NVTC,(NCMAP(NV,N),NV=1,NVTC)
          ENDDO
!
!---      Write cell-type map for active nodes to plot-1234.vtk file  ---
!
          WRITE(56,'(/,A,I9)') 'CELL_TYPES ',NACT
          NC = 0
          DO N = 1,NFLD
            IF( NMAP(N).EQ.0 ) CYCLE
            NC = NC + 1
            IVARX(NC) = 12
            IF( NC.EQ.10 .OR. N.EQ.NFLD ) THEN
              WRITE(56,'(10X,10(I2,1X))') (IVARX(M),M=1,NC)
              NC = 0
            ENDIF
          ENDDO
!
!---      Write node volume for active nodes to plot-1234.vtk file  ---
!
          WRITE(56,'(/,A,2(I9,1X))') 'CELL_DATA ',NACT
          NCH = INDEX(UNITS(4)(1:),'  ') - 1
          WRITE(56,'(3A)') 'SCALARS Node_Volume_(',UNITS(4)(1:NCH),
     &      ') float 1'
          WRITE(56,'(A)') 'LOOKUP_TABLE default'
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
!
!---      Write plot file variables for active nodes to 
!         plot-1234.vtk file  ---
!
          DO NV = 1,NVPLOT
            NCH = INDEX(VARNM(NV)(1:),'  ') - 1
            WRITE(56,'(/,3A)') 'SCALARS ',VARNM(NV)(1:NCH),' float 1'
            WRITE(56,'(A)') 'LOOKUP_TABLE default'
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
          ENDDO
!
!---    VTK_QUAD   ---
!
        ELSE
!
!---      Write x,y,z vertices to plot-1234.vtk file  ---
!
          WRITE(56,'(/,A,I9,A)') 'POINTS ',NPOINTS,' float'
          DO N = 1,NFLD
            DO M = 1,NVTC
              WRITE(56,'(3(1PE14.7,1X))') XE(M,N),YE(M,N),ZE(M,N)
            ENDDO
          ENDDO
!
!---      Write vertice map for active nodes to plot-1234.vtk file  ---
!
          WRITE(56,'(/,A,I9,1X,I9)') 'CELLS ',NACT,(NVTC+1)*NACT
          DO N = 1,NFLD
            IF( NMAP(N).EQ.0 ) CYCLE
            WRITE(56,'(9(I9,1X))') NVTC,(NCMAP(NV,N),NV=1,NVTC)
          ENDDO
!
!---      Write cell-type map for active nodes to plot-1234.vtk file  ---
!
          WRITE(56,'(/,A,I9)') 'CELL_TYPES ',NACT
          NC = 0
          DO N = 1,NFLD
            IF( NMAP(N).EQ.0 ) CYCLE
            NC = NC + 1
            IVARX(NC) = 9
            IF( NC.EQ.10 .OR. N.EQ.NFLD ) THEN
              WRITE(56,'(10X,10(I2,1X))') (IVARX(M),M=1,NC)
              NC = 0
            ENDIF
          ENDDO
!
!---      Write node volume for active nodes to plot-1234.vtk file  ---
!
          WRITE(56,'(/,A,2(I9,1X))') 'CELL_DATA ',NACT
          NCH = INDEX(UNITS(4)(1:),'  ') - 1
          WRITE(56,'(/,3A)') 'SCALARS Node_Volume_(',UNITS(4)(1:NCH),
     &      ') float 1'
          WRITE(56,'(A)') 'LOOKUP_TABLE default'
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
!
!---      Write plot file variables for active nodes to 
!         plot-1234.vtk file  ---
!
          DO NV = 1,NVPLOT
            NCH = INDEX(VARNM(NV)(1:),'  ') - 1
            WRITE(56,'(3A)') 'SCALARS ',VARNM(NV)(1:NCH),' float 1'
            WRITE(56,'(A)') 'LOOKUP_TABLE default'
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
          ENDDO
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
!---    Close the plot.1234 and plot-1234.csv files  ---
!
        CLOSE( UNIT=55 )
        CLOSE( UNIT=56 )
      ENDDO
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

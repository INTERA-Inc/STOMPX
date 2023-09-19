!----------------------Program-----------------------------------------!
!
      PROGRAM OUTPUT_TO_TECPLOT
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
!     This utility routine converts output files to Tecplot format
!     for plotting.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 February 2021.
!     outputTo_tec.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: VAR
      INTEGER, DIMENSION(:), ALLOCATABLE :: N_GC,N_FT,N_BH
      INTEGER, DIMENSION(:), ALLOCATABLE :: NVPK,NNPK
      CHARACTER(132), DIMENSION(:), ALLOCATABLE :: CHREF,SOLSPNM
      CHARACTER(132), DIMENSION(:), ALLOCATABLE :: CH_GC,CH_FT,CH_BH
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: CH_HD
      CHARACTER(512) :: CHDUM
      CHARACTER(1600) :: CHREC
      CHARACTER(132) :: CDUM,OUTFILE
      CHARACTER(32) :: HDUM
      CHARACTER(4) :: FORM1
      CHARACTER*14 FORM2
      LOGICAL FCHK
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
      DATA FORM2 /'(1X,1PE10.3,$)'/
!
!----------------------Executable Lines--------------------------------!
!
      IGEOMECH = 0
!
!---  No output file, exit program with error message  ---
!
      INQUIRE( FILE='output', EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        WRITE(6,'(A)') 'ERROR: Nonexistent "output" file.'
        STOP
      ELSE
        OPEN(UNIT=55, FILE='output', STATUS='OLD', FORM='FORMATTED')
      ENDIF
!
!---  Get inputs from command line or by prompting the user
!
      NCLA = IARGC()
      IF( NCLA.EQ.1 ) THEN
        CALL GETARG( 1,CDUM )
        CALL LCASE( CDUM )
        IF( INDEX(CDUM(1:),'-help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts data in the ' //
     &      'reference-node section of a STOMP output file to an '
          WRITE(6,'(A)') 'Tecplot "Point" formatted ' //
     &      'file. The script takes one argument, the name of '
          WRITE(6,'(A)') 'Tecplot file, which should have a "dat" ' //
     &      'extension (e.g., output.dat). If the argument is missing '
          WRITE(6,'(A)') 'the user will be prompted for the ' //
     &      'filename. The script determines the variables in the '
          WRITE(6,'(A)') 'output file, and prompts the user for ' //
     &      'a list of variables to be included in the output.dat '
          WRITE(6,'(A)') 'file. The variable order will be that ' //
     &      'specified by the user, which allows the variables to be '
          WRITE(6,'(A)') 'reordered from the order in the output ' //
     &      'file. Likewise, the script determines the nodes in the '
          WRITE(6,'(A)') 'output file and prompts the user for ' //
     &      'a list of nodes to be included in the output.dat file.'
          WRITE(6,'(A)') 'The node order will be that specified by ' //
     &      'the user, which allows the nodes to be reordered from the '
          WRITE(6,'(A,/)') 'order in the output file.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'output_to_tec.exe output.dat'
          WRITE(6,'(A)') 'output_to_tec.exe pressure.dat'
          WRITE(6,'(A)') 'output_to_tec.exe temperature.dat'
          WRITE(6,'(A,/)') 'output_to_tec.exe aqueous_saturation.dat'
          STOP
        ELSE
          OUTFILE = CDUM
        ENDIF
      ELSE
        WRITE(6,'(/,A,/)')'Enter the Tecplot filename ' // 
     &    '(e.g., output.dat, pressure.dat): '
        READ(5,*) OUTFILE
        CDUM = OUTFILE
        CALL LCASE( CDUM )
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts data in the ' //
     &      'reference-node section of a STOMP output file to an '
          WRITE(6,'(A)') 'Tecplot "Point" formatted ' //
     &      'file. The script takes one argument, the name of '
          WRITE(6,'(A)') 'Tecplot file, which should have a "dat" ' //
     &      'extension (e.g., output.dat). If the argument is missing '
          WRITE(6,'(A)') 'the user will be prompted for the ' //
     &      'filename. The script determines the variables in the '
          WRITE(6,'(A)') 'output file, and prompts the user for ' //
     &      'a list of variables to be included in the output.dat '
          WRITE(6,'(A)') 'file. The variable order will be that ' //
     &      'specified by the user, which allows the variables to be '
          WRITE(6,'(A)') 'reordered from the order in the output ' //
     &      'file. Likewise, the script determines the nodes in the '
          WRITE(6,'(A)') 'output file and prompts the user for ' //
     &      'a list of nodes to be included in the output.dat file.'
          WRITE(6,'(A)') 'The node order will be that specified by ' //
     &      'the user, which allows the nodes to be reordered from the '
          WRITE(6,'(A,/)') 'order in the output file.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'output_to_tec.exe output.dat'
          WRITE(6,'(A)') 'output_to_tec.exe pressure.dat'
          WRITE(6,'(A)') 'output_to_tec.exe temperature.dat'
          WRITE(6,'(A,/)') 'output_to_tec.exe aqueous_saturation.dat'
          STOP
        ENDIF
      ENDIF
!  
!---  Open output.dat file
!
      OUTFILE = TRIM(ADJUSTL(OUTFILE))
      WRITE(6,'(2A)') 'Creating output.dat ',TRIM(OUTFILE)
      OPEN(UNIT=56, FILE=OUTFILE)
      CLOSE(UNIT=56,STATUS='DELETE')
      OPEN(UNIT=56, FILE=OUTFILE, STATUS='NEW', FORM='FORMATTED')
!
!---  Write banner to screen  ---
!
      WRITE(6,'(A,//)')' Welcome to ...'
      WRITE(6,'(A,/)') '                       OUTPUT_TO_TECPLOT'
      WRITE(6,'(A)')   ' An output file converter routine for ' //
     &  'the STOMP simulator.  This'
      WRITE(6,'(A)')   ' routine, developed by the Pacific ' //
     &  'Northwest National Laboratory,'
      WRITE(6,'(A)')   ' reads a STOMP output ' //
     &  'file and converts the data in the reference'
      WRITE(6,'(A,/)')   ' node section to an Tecplot formatted ' //
     &  'file defined by the user '
      WRITE(6,'(A,/)')   ' (e.g., output.dat, pressure.dat).'
      WRITE(6,'(A)')   ' For support:  Tel: 509.372.6070'
      WRITE(6,'(A,/)') '               E-mail:  mark.white@pnnl.gov'
!
!---  Search for coordinate dimensions  ---
!
      DO
        READ(55,'(A)', IOSTAT=ICHK) CHDUM
        IF( INDEX(CHDUM(1:),'w/ Geomechanics').NE.0 ) IGEOMECH = 1
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: Number of Nodes'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: Number of Nodes'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Number of I-indexed Nodes:')
        JNDX = INDEX(CHDUM(1:),'Number of J-indexed Nodes:')
        KNDX = INDEX(CHDUM(1:),'Number of K-indexed Nodes:')
        IF( INDX.NE.0 ) THEN
          INDX = INDEX(CHDUM(1:),':')
          I1 = INDX+1
          I2 = INDX+12
          READ(CHDUM(I1:I2),'(I12)') IFLD
        ENDIF
        IF( JNDX.NE.0 ) THEN
          JNDX = INDEX(CHDUM(1:),':')
          J1 = JNDX+1
          J2 = JNDX+12
          READ(CHDUM(J1:J2),'(I12)') JFLD
        ENDIF
        IF( KNDX.NE.0 ) THEN
          KNDX = INDEX(CHDUM(1:),':')
          K1 = KNDX+1
          K2 = KNDX+12
          READ(CHDUM(K1:K2),'(I12)') KFLD
          EXIT
        ENDIF
      ENDDO
      NFLD = IFLD*JFLD*KFLD
!
!---  Search for reference node information  ---
!
      NREF = 0
      DO
        READ(55,'(A)', IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: Reference Node No. and Indices'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: Reference Node No. and Indices'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Fracture Nos.')
        IF( INDX.NE.0 ) THEN
          BACKSPACE( UNIT=55 )
          EXIT
        ENDIF
        INDX = INDEX(CHDUM(1:),'Borehole Nos.')
        IF( INDX.NE.0 ) THEN
          BACKSPACE( UNIT=55 )
          EXIT
        ENDIF
        INDX = INDEX(CHDUM(1:),'Reference Node Screen')
        IF( INDX.NE.0 ) THEN
          BACKSPACE( UNIT=55 )
          EXIT
        ENDIF
        INDX = INDEX(CHDUM(1:),'Reference Node No. and Indices')
        IF( INDX.NE.0 ) THEN
          INDX = INDEX(CHDUM(1:),':')
          I1 = INDX+2
          I2 = INDX+7
          READ(CHDUM(I1:I2),'(I6)') NREF
          EXIT
        ENDIF
      ENDDO
      ALLOCATE( N_GC(1:NREF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: N_GC'
        STOP
      ENDIF
      ALLOCATE( CH_GC(1:100),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: CH_GC'
        STOP
      ENDIF
      DO N = 1,NREF
        N_GC(N) = 0
        CH_GC(N) = ''
      ENDDO
      DO NC = 1,NREF
        READ(55,'(A)', IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: Reference Node No.'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: Reference Node No.'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Reference Node No.')
        IF( INDX.NE.0 ) THEN
          I1 = INDEX(CHDUM(1:),'.') + 1
          I2 = INDEX(CHDUM(1:),'I =') - 1
          READ(CHDUM(I1:I2),*) N_GC(NC)
          NCH = INDEX(CHDUM(1:),'   ')-1
          CH_GC(NC) = CHDUM(3:NCH)
        ELSE
          EXIT
        ENDIF
      ENDDO
      NC_GC = NREF
!
!---  Search for fracture triangle information  ---
!
      NREF_FRC = 0
      DO
        READ(55,'(A)', IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Borehole Nos.')
        IF( INDX.NE.0 ) THEN
          BACKSPACE( UNIT=55 )
          EXIT
        ENDIF
        INDX = INDEX(CHDUM(1:),'Reference Node Screen')
        IF( INDX.NE.0 ) THEN
          BACKSPACE( UNIT=55 )
          EXIT
        ENDIF
        INDX = INDEX(CHDUM(1:),'Fracture Nos., Local Triangle Nos.')
        IF( INDX.NE.0 ) THEN
          INDX = INDEX(CHDUM(1:),':')
          I1 = INDX+2
          I2 = INDX+7
          READ(CHDUM(I1:I2),'(I6)') NREF_FRC
          EXIT
        ENDIF
      ENDDO
      IF( NREF_FRC.GT.0 ) THEN
        ALLOCATE( N_FT(1:NREF_FRC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: N_FT'
          STOP
        ENDIF
        ALLOCATE( CH_FT(1:NREF_FRC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: CH_FT'
          STOP
        ENDIF
        DO N = 1,NREF_FRC
          N_FT(N) = 0
          CH_FT(N) = ''
        ENDDO
      ENDIF
      DO NC = 1,NREF_FRC
        READ(55,'(A)', IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: Global Triangle No. ='
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: Global Triangle No. ='
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Global Triangle No. =')
        IF( INDX.NE.0 ) THEN
          INDX = INDX + 20
          I1 = INDX+1
          I2 = INDEX(CHDUM(I1:),'  ') + I1
          READ(CHDUM(I1:I2),*) N_FT(NC)
          NCH = INDEX(CHDUM(1:),'   ')-1
          CH_FT(NC) = CHDUM(3:NCH)
        ELSE
          EXIT
        ENDIF
      ENDDO
      NC_FT = NREF_FRC
!
!---  Search for borehole node information  ---
!
      NREF_BH = 0
      DO
        READ(55,'(A)', IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Reference Node Screen')
        IF( INDX.NE.0 ) THEN
          BACKSPACE( UNIT=55 )
          EXIT
        ENDIF
        INDX = INDEX(CHDUM(1:),'Borehole Nos., Local Borehole Nos.')
        IF( INDX.NE.0 ) THEN
          INDX = INDEX(CHDUM(1:),':')
          I1 = INDX+2
          I2 = INDX+7
          READ(CHDUM(I1:I2),'(I6)') NREF_BH
          EXIT
        ENDIF
      ENDDO
      IF( NREF_BH.GT.0 ) THEN
        ALLOCATE( N_BH(1:NREF_BH),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: N_BH'
          STOP
        ENDIF
        ALLOCATE( CH_BH(1:NREF_BH),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: CH_BH'
          STOP
        ENDIF
        DO N = 1,NREF_BH
          N_BH(N) = 0
          CH_BH(N) = ''
        ENDDO
      ENDIF
      DO NC = 1,NREF_BH
        READ(55,'(A)', IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: Global Borehole Node No. ='
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: Global Borehole Node No. ='
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Reference Node Screen')
        IF( INDX.NE.0 ) THEN
          BACKSPACE( UNIT=55 )
          EXIT
        ENDIF
        INDX = INDEX(CHDUM(1:),'Global Borehole Node No. =')
        IF( INDX.NE.0 ) THEN
          INDX = INDX + 25
          I1 = INDX+1
          I2 = INDEX(CHDUM(I1:),'  ') + I1
          READ(CHDUM(I1:I2),*) N_BH(NC)
          NCH = INDEX(CHDUM(1:),'   ')-1
          CH_BH(NC) = CHDUM(3:NCH)
        ELSE
          EXIT
        ENDIF
      ENDDO
      NC_BH = NREF_BH
      IF( NC_GC+NC_FT+NC_BH.EQ.0 ) THEN
        WRITE(6,'(A)') 'No Reference Node, Fracture Triangle, or ' //
     &    'Borehole Node Information Found.'
        STOP
      ENDIF
!
!---  Find the number of output file significant digits  ---
!
      DO
        READ(55,'(A)',IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: Output File Significant Digits'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: Output File Significant Digits'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Output File Significant Digits')
        IF( INDX.NE.0 ) THEN
          EXIT
        ENDIF
      ENDDO
      INDX = INDEX(CHDUM(1:),':')
      I1 = INDX+1
      I2 = INDX+3
      READ(CHDUM(I1:I2),'(I3)') NOFSD
      WRITE(FORM2(8:9),'(I2)') MAX( 10,NOFSD+6 )
      WRITE(FORM2(11:11),'(I1)') MIN( 9,NOFSD-1 )
!
!---  Search for reference node variables  ---
!
      NVREF = 0
      DO
        READ(55,'(A)',IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: Reference Node Variables'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: Reference Node Variables'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Reference Node Variables')
        IF( INDX.NE.0 ) THEN
          INDX = INDEX(CHDUM(1:),':')
          I1 = INDX+2
          I2 = INDX+7
          READ(CHDUM(I1:I2),'(I6)') NVREF
          EXIT
        ENDIF
      ENDDO
!  
!---  Allocate memory
!
      IF( NVREF.GT.0 ) THEN
        ALLOCATE( CHREF(1:NVREF),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: CHREF'
          STOP
        ENDIF
        ALLOCATE( SOLSPNM(1:NVREF),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          WRITE(6,'(A)') 'Allocation Error: SOLSPNM'
          STOP
        ENDIF
        DO N = 1,NVREF
          CHREF(N) = ''
          SOLSPNM(N) = ''
        ENDDO
      ENDIF
      DO NC = 1,NVREF
        READ(55,'(A)', IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: Variable List'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: Variable List'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'   ')
        IF( INDX.NE.1 ) THEN
          INDX = INDEX(CHDUM(3:),'   ')
          I1 = 3
          I2 = INDX+2
          READ(CHDUM(I1:I2),'(A)') CHREF(NC)
!
!---      Species and solute names
!
          IF( INDEX(CHREF(NC)(1:),'species').NE.0 .OR. 
     &      INDEX(CHREF(NC)(1:),'solute').NE.0 ) THEN
            I1 = INDEX(CHREF(NC)(1:),'  ') + 2
            I2 = INDEX(CHREF(NC)(I1:),' S') + I1 -1
            SOLSPNM(NC) = CHREF(NC)(I1:I2)
          ENDIF
        ELSE
          EXIT
        ENDIF
      ENDDO
      NC_RNV = NVREF
!
!---  List variables and collect user input  ---
!
      ALLOCATE( NVPK(1:NC_RNV),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: NVPK'
        STOP
      ENDIF
      IF( NC_RNV.EQ.0 ) THEN
        WRITE(6,'(A)') 'No Reference Node Variables Found'
        STOP
      ELSE
        DO N = 1,NC_RNV
          NVPK(N) = 0
        ENDDO
        WRITE(6,'(A)') 'The STOMP output file, "output", contains the'
        WRITE(6,'(A,/)') 'following reference-node variables:'
        DO N = 1,NC_RNV
          NCH = INDEX(CHREF(N),'   ')-1
          WRITE(6,'(I3,2A)') N,' -- ',CHREF(N)(1:NCH)
        ENDDO
        WRITE(6,'(A,/)') '  a -- all reference-node variables'
        WRITE(6,'(A)') 'Enter the reference-node variables to include'
        WRITE(6,'(A)') 'in the "output.dat" file, by entering a'
        WRITE(6,'(A)') 'string of comma or space delimited indices or'
        WRITE(6,'(A,/)')'entering "a" for all reference node variables.'
        READ(5,'(A)') CHDUM
        NCH = 512
        I1 = 1
        NC = 0
        DO
          IS = INDEX(CHDUM(I1:I1),' ')
          IC = INDEX(CHDUM(I1:I1),',')
          IF( IS.NE.0 .OR. IC.NE.0 ) THEN
            I1 = I1 + 1
            IF( I1.GE.NCH ) EXIT
            CYCLE
          ENDIF
          IF( INDEX(CHDUM(I1:I1),'a').NE.0 .OR. 
     &      INDEX(CHDUM(I1:I1),'A').NE.0 ) THEN
            DO N = 1,NC_RNV
              NVPK(N) = N
            ENDDO
            EXIT
          ENDIF
          IF( ICHAR(CHDUM(I1:I1)).LT.48 .OR. 
     &      ICHAR(CHDUM(I1:I1)).GT.57 ) THEN
            WRITE(6,'(2A)') 'Unrecognized entry: ',CHDUM(I1:I1)
            STOP
          ENDIF
          IS = INDEX(CHDUM(I1:),' ')
          IC = INDEX(CHDUM(I1:),',')
          IF( IS.NE.0 .AND. IC.NE.0 ) THEN
            I2 = I1 + MIN(IS,IC) - 1
          ELSEIF( IC.NE.0 ) THEN
            I2 = I1 + IC - 1
          ELSEIF( IS.NE.0 ) THEN
            I2 = I1 + IS - 1
          ELSE
            I2 = NCH
          ENDIF
          NC = NC + 1
          READ(CHDUM(I1:I2),*,IOSTAT=ICHK) NVPK(NC)
          IF( ICHK.NE.0 ) THEN
            WRITE(6,'(2A)') 'Unrecognized entry: ',CHDUM(I1:I2)
            STOP
          ENDIF
          IF( NVPK(NC).LT.1 .OR. NVPK(NC).GT.NC_RNV ) THEN
            WRITE(6,'(A,I3)') 'Out-of-range selection: ',NVPK(NC)
            STOP
          ENDIF
          DO N = 1,NC-1
            IF( NVPK(NC).EQ.NVPK(N) ) THEN
              WRITE(6,'(A,I3)') 'Duplicate selection: ',NVPK(NC)
              STOP
            ENDIF
          ENDDO
          IF( I2.GE.NCH ) EXIT
          I1 = I2+1
        ENDDO
      ENDIF
!
!---  List grid nodes, fracture triangles, and borehole nodes
!     and collect user input  ---
!
      NC_TTL = NC_GC+NC_FT+NC_BH
      ALLOCATE( NNPK(1:NC_TTL),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: NNPK'
        STOP
      ENDIF
      IF( NC_RNV.EQ.0 ) THEN
        WRITE(6,'(A)') 'No Reference Node Variables Found'
        STOP
      ELSE
        DO N = 1,NC_TTL
         NNPK(N) = 0
        ENDDO
        WRITE(6,'(A)') 'The STOMP output file, "output", contains the'
        WRITE(6,'(A)')'following grid-cell nodes, fracture triangles,'
        WRITE(6,'(A,/)')'and borehole nodes:'
        NC = 0
        DO N = 1,NC_GC
          NCH = INDEX(CH_GC(N)(1:),'   ')-1
          NC = NC + 1
          WRITE(6,'(I3,2A)') NC,' -- ',CH_GC(N)(1:NCH)
        ENDDO
        DO N = 1,NC_FT
          NCH = INDEX(CH_FT(N)(1:),'   ')-1
          NC = NC + 1
          WRITE(6,'(I3,2A)') NC,' -- ',CH_FT(N)(1:NCH)
        ENDDO
        DO N = 1,NC_BH
          NCH = INDEX(CH_BH(N)(1:),'   ')-1
          NC = NC + 1
          WRITE(6,'(I3,2A)') NC,' -- ',CH_BH(N)(1:NCH)
        ENDDO
        WRITE(6,'(A,/)') '  a -- all grid-cell nodes, fracture ' // 
     &    'triangles, and borehole nodes'
        WRITE(6,'(A)') 'Enter the reference grid-cell nodes, ' //
     &    'fracture triangles, and borehole nodes to include'
        WRITE(6,'(A)') 'in the "output.dat" file, by entering a' //
     &    'string of comma or space delimited indices or'
        WRITE(6,'(A,/)')'entering "a" for all reference grid-cell ' //
     &    'nodes, fracture triangles, and borehole nodes'
        READ(5,'(A)') CHDUM
        NCH = 512
        I1 = 1
        NC = 0
        DO
          IS = INDEX(CHDUM(I1:I1),' ')
          IC = INDEX(CHDUM(I1:I1),',')
          IF( IS.NE.0 .OR. IC.NE.0 ) THEN
            I1 = I1 + 1
            IF( I1.GE.NCH ) EXIT
            CYCLE
          ENDIF
          IF( INDEX(CHDUM(I1:I1),'a').NE.0 .OR. 
     &      INDEX(CHDUM(I1:I1),'A').NE.0 ) THEN
            DO N = 1,NC_TTL
              NNPK(N) = N
            ENDDO
            EXIT
          ENDIF
          IF( ICHAR(CHDUM(I1:I1)).LT.48 .OR. 
     &      ICHAR(CHDUM(I1:I1)).GT.57 ) THEN
            WRITE(6,'(2A)') 'Unrecognized entry: ',CHDUM(I1:I1)
            STOP
          ENDIF
          IS = INDEX(CHDUM(I1:),' ')
          IC = INDEX(CHDUM(I1:),',')
          IF( IS.NE.0 .AND. IC.NE.0 ) THEN
            I2 = I1 + MIN(IS,IC) - 1
          ELSEIF( IC.NE.0 ) THEN
            I2 = I1 + IC - 1
          ELSEIF( IS.NE.0 ) THEN
            I2 = I1 + IS - 1
          ELSE
            I2 = NCH
          ENDIF
          NC = NC + 1
          READ(CHDUM(I1:I2),*,IOSTAT=ICHK) NNPK(NC)
          IF( ICHK.NE.0 ) THEN
            WRITE(6,'(2A)') 'Unrecognized entry: ',CHDUM(I1:I2)
            STOP
          ENDIF
          IF( NNPK(NC).LT.1 .OR. NNPK(NC).GT.NC_TTL ) THEN
            WRITE(6,'(A,I3)') 'Out-of-range selection: ',NNPK(NC)
            STOP
          ENDIF
          DO N = 1,NC-1
            IF( NNPK(NC).EQ.NNPK(N) ) THEN
              WRITE(6,'(A,I3)') 'Duplicate selection: ',NNPK(NC)
              STOP
            ENDIF
          ENDDO
          IF( I2.GE.NCH ) EXIT
          I1 = I2+1
        ENDDO
      ENDIF
!
!---  Find the start of the reference node output record  ---
!
      DO
        READ(55,'(A)',IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: ---  Reference Node Output Record  ---'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: ---  Reference Node Output Record  ---'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'---  Reference Node Output Record  ---')
        IF( INDX.NE.0 ) THEN
          EXIT
        ENDIF
      ENDDO
!
!---  Find the header line for the reference node output  ---
!
      DO
        READ(55,'(A)',IOSTAT=ICHK) CHREC
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: Header line'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: Header line'
          STOP
        ENDIF
        INDX = INDEX(CHREC(1:),'Step')
        IF( INDX.NE.0 ) THEN
          EXIT
        ENDIF
      ENDDO
      ALLOCATE( CH_HD(1:NVREF+5),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: CH_HD'
        STOP
      ENDIF
      DO N = 1,NVREF+5
        CH_HD(N) = ''
      ENDDO
      NCH = 1600
      I1 = 1
      NC = 0
      DO
        IS = INDEX(CHREC(I1:I1),' ')
        IF( IS.NE.0 ) THEN
          I1 = I1 + 1
          IF( I1.GE.NCH ) EXIT
          CYCLE
        ENDIF
        IS = INDEX(CHREC(I1:),' ')
        IF( IS.NE.0 ) THEN
          I2 = I1 + IS - 1
        ELSE
          I2 = NCH
        ENDIF
        NC = NC + 1
        READ(CHREC(I1:I2),*,IOSTAT=ICHK) CH_HD(NC)
        IF( ICHK.NE.0 ) THEN
          WRITE(6,'(2A)') 'Header Line Read Error: ',CHREC(I1:I2)
          STOP
        ENDIF
        IF( NC.GT.5 ) THEN
          J1 = 1
          J2 = INDEX(SOLSPNM(NC-5)(1:),'   ')-1
          IF( J2.GE.1 ) THEN
            JCH = INDEX(CH_HD(NC)(1:),'  ')
            CH_HD(NC)(JCH:JCH) = '-'
            CH_HD(NC)(JCH+1:JCH+J2) = SOLSPNM(NC-5)(J1:J2)
          ENDIF
        ENDIF
        IF( I2.GE.NCH ) EXIT
        I1 = I2+1
      ENDDO
      NHD = NC
!
!---  Write header line to output.dat file  ---
!
      WRITE(56,'(A)') 'TITLE = "STOMP Output File Visualization"'
      WRITE(56,'(A,$)') 'VARIABLES = '
      WRITE(56,'(A,$)') '"TM" '
      DO N = 1,NC_TTL
        IF( NNPK(N).EQ.0 ) EXIT
        DO M = 1,NC_RNV
          IF( NVPK(M).EQ.0 ) EXIT
          NV = NVPK(M) + 5
          NCH = INDEX(CH_HD(NV)(1:),' ')-1
          HDUM(1:NCH+1) = '"' // CH_HD(NV)(1:NCH)
          NCH = NCH + 2
          HDUM(NCH:NCH) = '_'
          NCH = NCH + 1
          IF( NNPK(N).LE.NC_GC ) THEN
            NN = NNPK(N)
            NCT = ICOUNT(N_GC(NN))
            WRITE(FORM1(3:3),'(I1)') NCT
            WRITE(HDUM(NCH:NCH+NCT-1),FORM1) N_GC(NN)
            WRITE(56,'(2A,$)') HDUM(1:NCH+NCT-1),'" '
          ELSEIF( NNPK(N).LE.NC_GC+NC_FT ) THEN
            HDUM(NCH:NCH+1) = '"T'
            NCH = NCH + 2
            NT = NNPK(N)-NC_GC
            NCT = ICOUNT(N_FT(NT))
            WRITE(FORM1(3:3),'(I1)') NCT
            WRITE(HDUM(NCH:NCH+NCT-1),FORM1) N_FT(NT)
            WRITE(56,'(2A,$)') HDUM(1:NCH+NCT-1),'" '
          ELSE
            HDUM(NCH:NCH) = '"B'
            NCH = NCH + 2
            NB = NNPK(N)-NC_GC-NC_FT
            NCT = ICOUNT(N_BH(NB))
            WRITE(FORM1(3:3),'(I1)') NCT
            WRITE(HDUM(NCH:NCH+NCT-1),FORM1) N_BH(NB)
            WRITE(56,'(2A,$)') HDUM(1:NCH+NCT-1),'" '
          ENDIF
        ENDDO
      ENDDO
      WRITE(56,'(A)') ''
      WRITE(56,'(A)') 'ZONE T="Reference Node Variables" F=POINT'
!  
!---  Allocate memory
!
      ALLOCATE( VAR(1:NC_TTL,1:NC_RNV),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: VAR'
        STOP
      ENDIF
!
!---  Search output file for reference node output records  ---
!
      DO
        READ(55,'(A)',IOSTAT=ICHK) CHREC
        IF( INDEX(CHREC(1:),'NOTE: Simulation Stopped').NE.0 ) EXIT
        IF( INDEX(CHREC(1:),'End of STOMP Simulation').NE.0 ) EXIT
        IF( INDEX(CHREC(1:),'Configuration Version Record').NE.0 ) EXIT
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during output ' //
     &      'file read: Reference node output records'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during output ' //
     &      'file read: Reference node output records'
          STOP
        ENDIF
!
!---    Check first entry of record for the time step  ---
!
        IFIND = 0
        DO I1 = 1,1600
          IF( ICHAR(CHREC(I1:I1)).EQ.32 ) CYCLE
          IF( ICHAR(CHREC(I1:I1)).GE.48 .AND. 
     &      ICHAR(CHREC(I1:I1)).LE.57 ) IFIND = 1
          EXIT
        ENDDO
!
!---    First entry is a time step, start sequence  ---
!
        IF( IFIND.EQ.1 ) THEN
          I1 = INDEX(CHREC(1:),'.')-2
          I2 = INDEX(CHREC(1:),'[')
          I3 = INDEX(CHREC(1:),'/')
          I4 = INDEX(CHREC(1:),']')
!
!---      Courant limit option  ---
!
          IF( I2.NE.0 ) THEN
            READ(CHREC(I1:I2-1),*) TM,TSTEP
            READ(CHREC(I4+1:),*) (VAR(1,J),J=1,NC_RNV)
!
!---      Geomechanics option  ---
!
          ELSEIF( IGEOMECH.EQ.1 ) THEN
            READ(CHREC(I1:),*) TM,TSTEP,NITER,MLP,(VAR(1,J),J=1,NC_RNV)
          ELSE
            READ(CHREC(I1:),*) TM,TSTEP,NITER,(VAR(1,J),J=1,NC_RNV)
          ENDIF
          DO I = 2,NC_TTL
            READ(55,'(A)',IOSTAT=ICHK) CHREC
            IF( ICHK.GT.0 ) THEN
              WRITE(6,'(A)') 'End of file encountered during output ' //
     &          'file read: Reference node output records'
              STOP
            ELSEIF( ICHK.LT.0 ) THEN
              WRITE(6,'(A)') 'Read error encountered during output ' //
     &          'file read: Reference node output records'
              STOP
            ENDIF
            I1 = INDEX(CHREC(1:),'.')-2
            I2 = INDEX(CHREC(1:),'[')
            I3 = INDEX(CHREC(1:),'/')
            I4 = INDEX(CHREC(1:),']')
!
!---        Courant limit option  ---
!
            IF( I2.NE.0 ) THEN
              READ(CHREC(I1:I2-1),*) TM,TSTEP
              READ(CHREC(I4+1:),*) (VAR(I,J),J=1,NC_RNV)
!
!---        Geomechanics option  ---
!
            ELSEIF( IGEOMECH.EQ.1 ) THEN
              READ(CHREC(I1:),*) TM,TSTEP,NITER,MLP,
     &          (VAR(I,J),J=1,NC_RNV)
            ELSE
              READ(CHREC(I1:),*) TM,TSTEP,NITER,(VAR(I,J),J=1,NC_RNV)
            ENDIF
          ENDDO
!
!---      Write variable record line to output.dat file  ---
!
          WRITE(56,'(1PE12.5,$)') TM
          DO N = 1,NC_TTL
            IF( NNPK(N).EQ.0 ) EXIT
            NN = NNPK(N)
            DO M = 1,NC_RNV
              IF( NVPK(M).EQ.0 ) EXIT
              NV = NVPK(M)
              WRITE(56,FORM2) VAR(NN,NV)
            ENDDO
          ENDDO
          WRITE(56,'(A)') ' '
        ENDIF        
      ENDDO
!
!---  Close the output and output.dat files  ---
!
      CLOSE( UNIT=55 )
      CLOSE( UNIT=56 )
      STOP
!
!---  End of OUTPUT_TO_TECPLOT program  ---
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

      

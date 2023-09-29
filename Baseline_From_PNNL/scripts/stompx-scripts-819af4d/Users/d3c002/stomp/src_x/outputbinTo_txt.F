!----------------------Program-----------------------------------------!
!
      PROGRAM OUTPUT_BIN_TO_TXT
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
!     This utility routine converts an output.bin binary 
!     formatted file to text (ASCII) formatted output file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 December 2021.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: VAR
      CHARACTER(128) :: CDUM,OUT_BIN_FILE,OUT_TXT_FILE
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: CHREF,UNREF
      CHARACTER*12 CHREFX
      CHARACTER*7  FORM1
      CHARACTER*27 FORM2
      CHARACTER*9  FORM3
      CHARACTER*12 FORM4
      CHARACTER*28 FORM5
      CHARACTER*14 FORM6
      CHARACTER*16 FORM7
      CHARACTER*14 FORM8
      CHARACTER*5  FORM9
      CHARACTER*28 FORM10
      CHARACTER*43 FORM11
      CHARACTER*39 FORM12
      CHARACTER*49 FORM13
      CHARACTER*12 FORM14
      CHARACTER*40 FORM15
      CHARACTER*14 FORM16
      CHARACTER*4  FORM17
      CHARACTER*14 FORM18
      CHARACTER*5  FORM19
      CHARACTER*10 FORM20
      CHARACTER*6 FORM21
      CHARACTER*42 FORM22
      CHARACTER*33 FORM23
      CHARACTER*49 FORM24
      CHARACTER*22 FORM25
      CHARACTER*38 FORM26
      CHARACTER*48 FORM27
      CHARACTER*42 FORM28
      CHARACTER*20 FORM29
      CHARACTER*64 SPACES
      CHARACTER*64 TXTX,STRX,UNTM
      CHARACTER*4096 N_IWR
      CHARACTER*64 PH_CND_CO2(10)
      CHARACTER*128 PH_CND_HYD(43)
      LOGICAL FCHK
      INTEGER :: NBYTR = 8
      INTEGER :: NBYTI = 4
      INTEGER :: NBYTC = 1
      INTEGER, DIMENSION(:), ALLOCATABLE :: NDREF
      INTEGER, DIMENSION(16) :: IVAR
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(/,A,$)'/
      DATA FORM2 /'(3X,A,6X,A,5X,A,7X,A,2X,$)'/
      DATA FORM3 /'(3X,2A,$)'/
      DATA FORM4 /'(4X,A4,A,$)'/
      DATA FORM5 /'(A,I3,A,I3,A,I3,A,I8,A,$)'/
      DATA FORM6 /'(1X,1PE10.3,$)'/
      DATA FORM7 /'(20X,3A,5X,4A,$)'/
      DATA FORM8 /'(1X,A,A7,2A,$)'/
      DATA FORM9 /'(A,$)'/
      DATA FORM10 /'(3X,A,6X,A,5X,A,7X,A,3X,A,$)'/
      DATA FORM11 /'(1X,I6,1X,I9,1X,1PE12.5,1X,1PE12.5,1X,I2,$)'/
      DATA FORM12 /'(1X,I6,1X,I9,1X,1PE12.5,1X,1PE12.5,A,$)'/
      DATA FORM13 /'(1X,I6,1X,I9,1X,1PE12.5,1X,1PE12.5,1X,I2,1X,I6,$)'/
      DATA FORM14 /'(4X,A4,A,$)'/
      DATA FORM15 /'(3X,A,6X,A,5X,A,7X,A,3X,A,1X,A,$)'/
      DATA FORM16 /'(1X,1PE10.3,$)'/
      DATA FORM17 /'(1X)'/
      DATA FORM18 /'(1X,A,A7,2A,$)'/
      DATA FORM19 /'(A,$)'/
      DATA FORM20 /'(10X,2A,$)'/
      DATA FORM21 /'(I1)'/
      DATA FORM22 /'(A,I3,A,I3,A,I3,A,I3,A,I3,A,I3,A,I8,A,$)'/
      DATA FORM23 /'(3X,A,6X,A,5X,A,7X,A,4X,A,2X,A,$)'/
      DATA FORM24 /'(1X,I6,1X,I9,1X,1PE12.5,1X,1PE12.5,1X,I3,2X,I3,$)'/
      DATA FORM25 /'(A,I4,A,I4,A,I8,A,$)'/
      DATA FORM26 /'(1X,I6,1X,A,1X,1PE12.5,1X,1PE12.5,A,$)'/
      DATA FORM27 /'(1X,I6,1X,A,1X,1PE12.5,1X,1PE12.5,1X,I2,2X,I2,$)'/
      DATA FORM28 /'(1X,I6,1X,A,1X,1PE12.5,1X,1PE12.5,1X,I2,$)'/
      DATA FORM29 /'(A,I4,A,I8,A,I8,A,$)'/
      DATA SPACES /'                                                    
     &      '/
      DATA PH_CND_CO2 /'Saturated w/ Aqueous CO2',
     &   'Unsaturated Subcritical Pressure Gas', 
     &   'Saturated w/ Trapped Subcritical Pressure Gas',
     &   'Unsaturated Subcritical Pressure Liquid', 
     &   'Saturated w/ Trapped Subcritical Presssure Liquid',
     &   'Unsaturated Supercritical Pressure Gas/Liquid', 
     &   'Saturated w/ Trapped Supercritical Pressure Gas/Liquid',
     &   'Fully Unsaturated Subcritical Pressure Gas',
     &   'Fully Unsaturated Subcritical Pressure Liquid',
     &   'Fully Unsaturated Supercritical Pressure Gas/Liquid'/
      DATA PH_CND_HYD /9*'','SH = 0, SG = 0, SN = 0',10*'',
     &  'SH = 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN )',
     &  'SH = 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN )',
     &  'SH = 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO )',7*'',
     &  'SH > 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN ), YMHGA >= MAX(YMH
     &GO,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN ), YMHGA >= MAX(YMH
     &GO,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO ), YMHGA >= MAX(YMH
     &GO,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN ), YMHGO >= MAX(YMH
     &GA,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN ), YMHGO >= MAX(YMH
     &GA,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO ), YMHGO >= MAX(YMH
     &GA,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN ), YMHGN >= MAX(YMH
     &GA,YMHGO )',
     &  'SH > 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN ), YMHGN >= MAX(YMH
     &GA,YMHGO )',
     &  'SH > 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO ), YMHGN >= MAX(YMH
     &GA,YMHGO )','',
     &  'SH > 0, SG + SN = 0, YMHGA >= MAX(YMHGO,YMHGN )',
     &  'SH > 0, SG + SN = 0, YMHGO >= MAX(YMHGA,YMHGN )',
     &  'SH > 0, SG + SN = 0, YMHGN >= MAX(YMHGA,YMHGO )'/
!
!----------------------Executable Lines--------------------------------!
!
      ISKIP = 0
      NCLA = IARGC()
      IF( NCLA.EQ.0 ) THEN
        WRITE(6,'(/,A,/)')'Enter text output file name or "help": '
        READ(5,*) OUT_TXT_FILE
        CALL LCASE( OUT_TXT_FILE )
        CDUM = OUT_TXT_FILE
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts a binary formatted ' // 
     &      'output.bin file created from executing STOMPX'
          WRITE(6,'(A)') 'to a conventional text (ASCII) formatted' //
     &      ' output file created from executing STOMP.'
          WRITE(6,'(A)') 'The created text file only contains' //
     &      ' the reference node section of the output. This file'
          WRITE(6,'(A)') 'should be appended to the output file' //
     &      ' created during the parallel preprocessing stage'
          WRITE(6,'(A)') 'to create a full output created from' //
     &      ' executing STOMP.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'outputbinTo.exe output2'
          WRITE(6,'(A,/)') 'cat output2 >> output '
          STOP
        ENDIF
      ELSE
        CALL GETARG( 1,OUT_TXT_FILE )
        CALL LCASE( OUT_TXT_FILE )
        CDUM = OUT_TXT_FILE
        CALL LCASE( CDUM )
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts a binary formatted ' // 
     &      'output.bin file created from executing STOMPX'
          WRITE(6,'(A)') 'to a conventional text (ASCII) formatted' //
     &      ' output file created from executing STOMP.'
          WRITE(6,'(A)') 'The created text file only contains' //
     &      ' the reference node section of the output. This file'
          WRITE(6,'(A)') 'should be appended to the output file' //
     &      ' created during the parallel preprocessing stage'
          WRITE(6,'(A)') 'to create a full output created from' //
     &      ' executing STOMP.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'outputbinTo.exe output2'
          WRITE(6,'(A,/)') 'cat output2 >> output '
          STOP
        ENDIF
      ENDIF
!
!---  Write banner to screen  ---
!       
      WRITE(6,'(A,//)')' Welcome to ...'
      WRITE(6,'(A,/)') '                       OUTPUT_BIN_TO_TXT'
      WRITE(6,'(A)')   ' A binary output file converter ' //
     &  'routine for the STOMPX simulator.'
      WRITE(6,'(A)')   ' This routine, developed by the Pacific ' //
     &  'Northwest National Laboratory,'
      WRITE(6,'(A)')   ' converts a STOMPX binary formatted ' //
     &  ' output.bin file to '
      WRITE(6,'(A)') ' a text (ASCII) formatted user specified' //
     &  ' file. The created text '
      WRITE(6,'(A)') ' file only contains' //
     &  ' the reference node section of the output. This file'
      WRITE(6,'(A)') ' should be appended to the output file' //
     &  ' created during the parallel '
      WRITE(6,'(A)') ' preprocessing stage to create a full ' //
     &  ' output created from executing STOMP.'
      WRITE(6,'(A,/)') ' For example with the following command: ' //
     & ' cat output2 >> output '
      WRITE(6,'(A)')   ' For support:  Tel: 509.372.6070'
      WRITE(6,'(A,/)') '               E-mail:  mark.white@pnnl.gov'
!
!---  The binary formatted output file created by STOMPX uses stream
!     access, requiring a tracking pointer to access the data. This
!     file only contains the information in the reference node
!     section of the output file and should be appended to the 
!     output file created in the parallel preprocessing step. The file
!     contains header information with the number of reference
!     nodes, reference node variables, and active equations then 
!     for each time step and reference node the time step number, 
!     number of iterations, time,time step, node number, and reference 
!     node variables. Negative time steps indicate failures:
!
!     -1 = Convergence Failure (cut time step)
!          node with maximum residual, phase condition of that node
!          and residuals for each active equation + old time step,
!          new time step
!     -2 = Exceed limit on time step reduction (stop simulation)
!     -3 = Linear system failure (cut time step)
!          old time step, new time step
!     -4 = Excessive primary variable failure (cut time step)
!          lowest global node number with excessive primary variable
!          phase condition of that node, and changes in primary 
!          variables for that node for each equation
!          old time step + new time step
!
      OUT_BIN_FILE = 'output.bin'
!
!---  No output.bin file, exit program with error message  ---
!
      INQUIRE( FILE=TRIM(OUT_BIN_FILE), EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        WRITE(6,'(A)') 'ERROR: Nonexistent output.bin file.'
        STOP
      ELSE
        OPEN(UNIT=55, FILE=TRIM(OUT_BIN_FILE), STATUS='OLD', 
     &    FORM='UNFORMATTED', ACCESS='STREAM')
        OPEN(UNIT=56, FILE=TRIM(OUT_TXT_FILE), STATUS='UNKNOWN', 
     &    FORM='FORMATTED')
        CLOSE(UNIT=56, STATUS='DELETE')
        OPEN(UNIT=56, FILE=TRIM(OUT_TXT_FILE), STATUS='NEW', 
     &    FORM='FORMATTED')
      ENDIF
      WRITE(6,'(4A)') 'Converting binary output file: ',
     &  TRIM(OUT_BIN_FILE),' to text (ASCII) output file: ',
     &  TRIM(OUT_TXT_FILE)
!
!---  Initialize binary file position  ---
!
      IPOS = 1
!
!---  Read number of reference nodes  ---
!
      READ(55, POS=IPOS, IOSTAT=ICHK) NREF
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'NREF: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'NREF: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NBYTI
!  
!---  Allocate memory for reference node variables
!
      ALLOCATE( NDREF(1:NREF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: NDREF'
        STOP
      ENDIF
!  
!---  Read number of reference node variables
!
      READ(55, POS=IPOS, IOSTAT=ICHK) NVREF
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'NVREF: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'NVREF: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NBYTI
!  
!---  Allocate memory for reference node variables and reference
!     node variable header strings
!
      NVARX = MAX( 28,NVREF )
      ALLOCATE( VAR(1:NVARX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: VAR'
        STOP
      ENDIF
      ALLOCATE( CHREF(1:NVREF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: CHREF'
        STOP
      ENDIF
      ALLOCATE( UNREF(1:NVREF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: UNREF'
        STOP
      ENDIF
!  
!---  Read number of active equations
!
      READ(55, POS=IPOS, IOSTAT=ICHK) ISVC
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'ISVC: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'ISVC: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NBYTI
!  
!---  Read number of output file significant digits
!
      READ(55, POS=IPOS, IOSTAT=ICHK) ISGNO
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'ISGNO: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'ISGNO: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NBYTI
!  
!---  Read number of x-direction nodes
!
      READ(55, POS=IPOS, IOSTAT=ICHK) IFLD
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'IFLD: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'IFLD: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NBYTI
!  
!---  Read number of y-direction nodes
!
      READ(55, POS=IPOS, IOSTAT=ICHK) JFLD
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'JFLD: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'JFLD: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NBYTI
!  
!---  Read number of z-direction nodes
!
      READ(55, POS=IPOS, IOSTAT=ICHK) KFLD
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'KFLD: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'KFLD: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NBYTI
!  
!---  Read reference header strings
!
      READ(55, POS=IPOS, IOSTAT=ICHK) (CHREF(M),M=1,NVREF)
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'CHREF: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'CHREF: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NVREF*NBYTC*64
!  
!---  Read reference unit header strings
!
      READ(55, POS=IPOS, IOSTAT=ICHK) (UNREF(M),M=1,NVREF)
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'UNREF: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'UNREF: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NVREF*NBYTC*64
!  
!---  Read time unit string
!
      READ(55, POS=IPOS, IOSTAT=ICHK) UNTM
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'UNTM: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'UNTM: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NBYTC*64
!  
!---  Read reference nodes
!
      READ(55, POS=IPOS, IOSTAT=ICHK) (NDREF(M),M=1,NREF)
      IF( ICHK.GT.0 ) THEN
        WRITE(6,'(2A)') 'NDREF: End of file encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ELSEIF( ICHK.LT.0 ) THEN
        WRITE(6,'(2A)') 'NDREF: Read error encountered during ' //
     &    'output.bin file read: ',TRIM(OUT_BIN_FILE)
        STOP
      ENDIF
      IPOS = IPOS + NREF*NBYTI
!
!---  Compute format spacings according to the number of significant
!     digits requested  ---
!
      IVARX = (ISGNO-2)/2
      WRITE(FORM14(2:2),'(I1)') MAX( 4,ISGNO+3-IVARX )
      WRITE(FORM16(8:9),'(I2)') MAX( 10,ISGNO+6 )
      WRITE(FORM16(11:11),'(I1)') MIN( 9,ISGNO-1 )
      WRITE(FORM18(2:2),'(I1)') MAX( 1,ISGNO-2-IVARX )
!
!---  Write main header line  ---
!
      WRITE(56,'(//,A)') ' ---  Reference Node Output Record  ---'
      ICNO = 10
!  
!---  Loop over time steps in simulation
!
      DO
!
!---    Write subheader lines to the output file ---
!
        IF( ICNO.EQ.10 ) THEN
          ICNO = 0
          WRITE(56,FORM1) 'Reference Node(s)'
          DO M = 1,NREF
            N = NDREF(M)
            K = N/(IFLD*JFLD) + 1
            IF( MOD( N,(IFLD*JFLD) ).EQ.0 ) K = K - 1
            J = N/IFLD + 1
            IF( MOD( N,IFLD ).EQ.0 ) J = J - 1
            I = MOD(N,IFLD)
            IF( MOD( N,IFLD ).EQ.0 ) I = IFLD
            WRITE(FORM5(5:5),'(I1)') ICOUNT(I)
            WRITE(FORM5(10:10),'(I1)') ICOUNT(J)
            WRITE(FORM5(15:15),'(I1)') ICOUNT(K)
            WRITE(FORM5(20:20),'(I1)') ICOUNT(N)
            WRITE(56,FORM5) ' (',I,',',J,',',K,':',N,')'
          ENDDO
          WRITE(56,FORM17)
          IN_IWR = 0
          WRITE(56,FORM10) 'Step','Node','Time','Timestep','Itr'
          INX = 1
          WRITE(N_IWR(INX:INX+46),'(A)') SPACES(1:47)
          INX = INX+47
          DO NV = 1,NVREF
            ISX = ISGNO+7
            CHREFX = '            '
            CHREFX(1:6) = CHREF(NV)
            TXTX = CHREFX
            ITX = INDEX(TXTX(1:),' ') - 1
            CALL CNTRTXT( TXTX,STRX,ITX,ISX )
            WRITE(56,'(A,$)') STRX(1:ISX)
            IN_IWR = 1
          ENDDO
          WRITE(56,FORM17)
          IF( IN_IWR.EQ.1 ) WRITE(56,'(A)') N_IWR(1:INX-1)
          WRITE(56,FORM7) '[',UNTM(1:6),']','[',UNTM(1:6),']',
     &      SPACES(1:6)
          DO NV = 1,NVREF
            ITX = 9
            IF( UNREF(NV).NE.'null' ) THEN
               TXTX = '[' // UNREF(NV)(1:7) // ']'
               CALL CNTRTXT( TXTX,STRX,ITX,ISX )
               WRITE(56,'(A,$)') STRX(1:ISX)
            ELSE
               TXTX = SPACES(1:9)
               CALL CNTRTXT( TXTX,STRX,ITX,ISX )
               WRITE(56,'(A,$)') STRX(1:ISX)
            ENDIF
          ENDDO
          WRITE(56,FORM17)
        ENDIF
!
!---    Read time step  ---
!
        READ(55, POS=IPOS, IOSTAT=ICHK) NSTEP
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(2A)') 'NSTEP: End of file encountered during ' //
     &      'output.bin file read: ',TRIM(OUT_BIN_FILE)
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(2A)') 'NSTEP: Read error encountered during ' //
     &      'output.bin file read: ',TRIM(OUT_BIN_FILE)
          STOP
        ENDIF
        IPOS = IPOS + NBYTI
!
!---    Postive time step number, read normal reference node output  ---
!
        IF( NSTEP.GE.0 ) THEN
!
!---      Read number of iterations  ---
!
          READ(55, POS=IPOS, IOSTAT=ICHK) NITER
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'NSTEP: End of file encountered ' //
     &        'during output.bin file read: ',TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'NSTEP: Read error encountered ' //
     &        'during output.bin file read: ',TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NBYTI
!  
!---      Read time  ---
!
          READ(55, POS=IPOS, IOSTAT=ICHK) TM
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'TM: End of file encountered during ' //
     &        'output.bin file read: ',TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'TM: Read error encountered during ' //
     &        'output.bin file read: ',TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NBYTR
!  
!---      Read time step  ---
!
          READ(55, POS=IPOS, IOSTAT=ICHK) DT
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'TM: End of file encountered during ' //
     &        'output.bin file read: ',TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'TM: Read error encountered during ' //
     &        'output.bin file read: ',TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NBYTR
!  
!---      Loop over reference nodes  ---
!
          DO M = 1,NREF
!
!---        Read node number  ---
!
            READ(55, POS=IPOS, IOSTAT=ICHK) N
            IF( ICHK.GT.0 ) THEN
              WRITE(6,'(2A)') 'NSTEP: End of file encountered ' //
     &          'during output.bin file read: ',TRIM(OUT_BIN_FILE)
              STOP
            ELSEIF( ICHK.LT.0 ) THEN
              WRITE(6,'(2A)') 'NSTEP: Read error encountered ' //
     &          'during output.bin file read: ',TRIM(OUT_BIN_FILE)
              STOP
            ENDIF
            IPOS = IPOS + NBYTI
!  
!---        Read reference node variables  ---
!
            READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(NV),NV=1,NVREF)
            IF( ICHK.GT.0 ) THEN
              WRITE(6,'(2A)') 'TM: End of file encountered during ' //
     &          'output.bin file read: ',TRIM(OUT_BIN_FILE)
              STOP
            ELSEIF( ICHK.LT.0 ) THEN
              WRITE(6,'(2A)') 'TM: Read error encountered during ' //
     &          'output.bin file read: ',TRIM(OUT_BIN_FILE)
              STOP
            ENDIF
            IPOS = IPOS + NVREF*NBYTR
!  
!---        Write line to text formatted output file  ---
!
            WRITE(56,FORM11) NSTEP,N,TM,DT,NITER
            DO NV = 1,NVREF
              WRITE(56,FORM16) VAR(NV)
            ENDDO
            WRITE(56,FORM17)
          ENDDO
          ICNO = ICNO + 1
!
!---    NSTEP = -1, Convergence failure for STOMP-CO2
!       (cut time step)   ---
!
        ELSEIF( NSTEP.EQ.-1 ) THEN
!  
!---      Read convergence failure variables  ---
!
          NVAR = 6
          READ(55, POS=IPOS, IOSTAT=ICHK) (IVAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Integer: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Integer: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTI
          NVAR = 5
          READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Real: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Real: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTR
          IF( ISKIP.EQ.0 ) 
     &      WRITE(56,'(10X,A)') '---  Convergence Failure  ---'
          ISKIP = 0
          WRITE(56,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Equation Maximum Residual = ',VAR(1),
     &      ': Node = ',IVAR(1),
     &      ': Phase Condition = ',TRIM(PH_CND_CO2(IVAR(2)))
          WRITE(56,'(4X,A,1PE11.4,A,I6,2A)')
     &      'CO2 Equation Maximum Residual = ',VAR(2),
     &      ': Node = ',IVAR(3),
     &      ': Phase Condition = ',TRIM(PH_CND_CO2(IVAR(4)))
          IF( IVAR(6).GT.0 ) THEN
            WRITE(56,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Salt Equation Maximum Residual = ',VAR(3),
     &      ': Node = ',IVAR(5),
     &      ': Phase Condition = ',TRIM(PH_CND_CO2(IVAR(6)))
          WRITE(56,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VAR(4),TRIM(UNTM),' to ',
     &      VAR(5),TRIM(UNTM)
          ENDIF
!
!---    NSTEP = -2, Exceed limit on time step reduction 
!         (stop simulation)   ---
!
        ELSEIF( NSTEP.EQ.-2 ) THEN
          EXIT
!
!---    NSTEP = -3, Linear system failure (cut time step)   ---
!
        ELSEIF( NSTEP.EQ.-3 ) THEN
          NVAR = 2
          READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Time Step Real: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Time Step Real: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTR
          IF( ISKIP.EQ.0 ) 
     &      WRITE(56,'(10X,A)') '---  Linear System Solver Failure  ---'
          ISKIP = 0
          WRITE(56,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VAR(1),TRIM(UNTM),' to ',
     &      VAR(2),TRIM(UNTM)
!
!---    NSTEP = -4, Excessive primary variable change failure for 
!       STOMPX-CO2 (cut time step)   ---
!
        ELSEIF( NSTEP.EQ.-4 ) THEN
!  
!---      Read excessive residual variables  ---
!
          NVAR = 2
          READ(55, POS=IPOS, IOSTAT=ICHK) (IVAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Integer: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Integer: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTI
          NVAR = 9
          READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Real: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Real: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTR
          IF( ISKIP.EQ.0 ) 
     &      WRITE(56,'(10X,A)') '---  Excessive Primary Variable ' // 
     &        'Change  ---'
          ISKIP = 0
          WRITE(56,'(4X,A,I6)') 'Node = ',IVAR(1)
          WRITE(56,'(4X,2A)') 'Phase Condition = ',PH_CND_CO2(IVAR(2))
          WRITE(56,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',VAR(1)
          WRITE(56,'(4X,A,1PE12.5)') 'Gas Pressure = ',VAR(2)
          WRITE(56,'(4X,A,1PE12.5)') 'Gas Saturation = ',VAR(3)
          WRITE(56,'(4X,A,1PE12.5)')
     &      'Aqueous-CO2 Mass Fraction = ',VAR(4)
          WRITE(56,'(4X,A,1PE12.5)')'Water Vapor Pressure = ',VAR(5)
          WRITE(56,'(4X,A,1PE12.5,A,I6)')
     &      'Total-Salt Aqu. Mass Fraction = ',VAR(6)
          WRITE(56,'(4X,A,1PE12.5,A,I6)')
     &      'Salt Volumetric Concentration = ',VAR(7)
          WRITE(56,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VAR(8),TRIM(UNTM),' to ',
     &      VAR(9),TRIM(UNTM)
!
!---    NSTEP = -5, Excessive residual for STOMPX-CO2 
!       (cut time step)   ---
!
        ELSEIF( NSTEP.EQ.-5 ) THEN
!  
!---      Read excessive residual variables  ---
!
          NVAR = 6
          READ(55, POS=IPOS, IOSTAT=ICHK) (IVAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Integer: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Integer: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTI
          NVAR = 5
          READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Real: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Real: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTR
          IF( ISKIP.EQ.0 ) 
     &      WRITE(56,'(10X,A)') '---  Excessive Residual  ---'
          ISKIP = 0
          WRITE(56,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Equation Maximum Residual = ',VAR(1),
     &      ': Node = ',IVAR(1),
     &      ': Phase Condition = ',TRIM(PH_CND_CO2(IVAR(2)))
          WRITE(56,'(4X,A,1PE11.4,A,I6,2A)')
     &      'CO2 Equation Maximum Residual = ',VAR(2),
     &      ': Node = ',IVAR(3),
     &      ': Phase Condition = ',TRIM(PH_CND_CO2(IVAR(4)))
          IF( IVAR(6).GT.0 ) THEN
            WRITE(56,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',VAR(3),
     &        ': Node = ',IVAR(5),
     &        ': Phase Condition = ',TRIM(PH_CND_CO2(IVAR(6)))
            WRITE(56,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &        'Time Step Reduced From ',VAR(4),TRIM(UNTM),' to ',
     &        VAR(5),TRIM(UNTM)
          ENDIF
!
!---    NSTEP = -7, Excessive residual with coupled well
!       for STOMPX-CO2   ---
!
        ELSEIF( NSTEP.EQ.-7 ) THEN
!  
!---      Read Excessive residual variables  ---
!
          NVAR = 1
          READ(55, POS=IPOS, IOSTAT=ICHK) (IVAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Integer: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Integer: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTI
          NVAR = 1
          READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Real: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Real: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTR
          WRITE(56,'(10X,A)') '---  Excessive Residual  ---'
          ISKIP = 1
          WRITE(56,'(4X,A,1PE11.4,A,I6)')
     &      'Coupled Well Maximum Residual = ',VAR(1),
     &      ': Coupled Well Number = ',IVAR(1)
!
!---    NSTEP = -8, Convergence failure with coupled well   ---
!
        ELSEIF( NSTEP.EQ.-8 ) THEN
!  
!---      Read convergence failure variables  ---
!
          NVAR = 1
          READ(55, POS=IPOS, IOSTAT=ICHK) (IVAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Integer: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Integer: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTI
          NVAR = 1
          READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Real: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Real: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTR
          WRITE(56,'(10X,A)') '---  Convergence Failure  ---'
          ISKIP = 1
          WRITE(56,'(4X,A,1PE11.4,A,I6)')
     &      'Coupled Well Maximum Residual = ',VAR(1),
     &      ': Coupled Well Number = ',IVAR(1)
!
!---    NSTEP = -9, End of file   ---
!
        ELSEIF( NSTEP.EQ.-9 ) THEN
          EXIT
!
!---    NSTEP = -11, Convergence failure for STOMP-HYDT-KE
!      (cut time step)   ---
!
        ELSEIF( NSTEP.EQ.-11 ) THEN
!  
!---      Read convergence failure variables  ---
!
          NVAR = 18
          READ(55, POS=IPOS, IOSTAT=ICHK) (IVAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Integer: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Integer: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTI
          NVAR = 11
          READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Real: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Convergence Failure Real: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTR
          IF( ISKIP.EQ.0 ) 
     &      WRITE(56,'(10X,A)') '---  Convergence Failure  ---'
          ISKIP = 0
          WRITE(56,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Energy Equation Maximum Residual = ',VAR(1),
     &      ': Node = ',IVAR(1),
     &      ': Phase Condition = ',TRIM(PH_CND_HYD(IVAR(2)))
          WRITE(56,'(4X,A,1PE11.4,A,I6,2A)')
     &      'CO2 Equation Maximum Residual = ',VAR(2),
     &      ': Node = ',IVAR(3),
     &      ': Phase Condition = ',TRIM(PH_CND_HYD(IVAR(4)))
          IF( IVAR(6).GT.0 ) THEN
            WRITE(56,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Salt Equation Maximum Residual = ',VAR(3),
     &      ': Node = ',IVAR(5),
     &      ': Phase Condition = ',TRIM(PH_CND_HYD(IVAR(6)))
            WRITE(56,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &        'Time Step Reduced From ',VAR(4),TRIM(UNTM),' to ',
     &        VAR(5),TRIM(UNTM)
          ENDIF
!
!---    NSTEP = -14, Excessive primary variable change failure for 
!       STOMPX-HYDT-KE (cut time step)   ---
!
        ELSEIF( NSTEP.EQ.-14 ) THEN
!  
!---      Read excessive residual variables  ---
!
          NVAR = 5
          READ(55, POS=IPOS, IOSTAT=ICHK) (IVAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Integer: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Integer: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTI
          NVAR = 28
          READ(55, POS=IPOS, IOSTAT=ICHK) (VAR(NV),NV=1,NVAR)
          IF( ICHK.GT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Real: ' //
     &        'End of file encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ELSEIF( ICHK.LT.0 ) THEN
            WRITE(6,'(2A)') 'Excessive Residual Real: Read ' //
     &        'error encountered during output.bin file read: ',
     &        TRIM(OUT_BIN_FILE)
            STOP
          ENDIF
          IPOS = IPOS + NVAR*NBYTR
          IF( ISKIP.EQ.0 ) 
     &      WRITE(56,'(10X,A)') '---  Excessive Primary Variable ' // 
     &        'Change  ---'
          ISKIP = 0
          WRITE(56,'(4X,A,I6)') 'Node = ',IVAR(1)
          WRITE(56,'(4X,2A)') 'Phase Condition = ',PH_CND_HYD(IVAR(2))
          WRITE(56,'(4X,A,1PE12.5)') 'DPT = ',VAR(1)
          WRITE(56,'(4X,A,1PE12.5)') 'DPW = ',VAR(2)
          WRITE(56,'(4X,A,1PE12.5)') 'DPA = ',VAR(3)
          WRITE(56,'(4X,A,1PE12.5)') 'DPO = ',VAR(4)
          WRITE(56,'(4X,A,1PE12.5)') 'DPN = ',VAR(5)
          WRITE(56,'(4X,A,1PE12.5)') 'DPHA = ',VAR(6)
          WRITE(56,'(4X,A,1PE12.5)') 'DPHO = ',VAR(7)
          WRITE(56,'(4X,A,1PE12.5)') 'DPHN = ',VAR(8)
          WRITE(56,'(4X,A,1PE12.5)') 'Temperature, C = ',VAR(9)
          WRITE(56,'(4X,A,1PE12.5)') 'Pressure, Pa = ',VAR(10)
          WRITE(56,'(4X,A,1PE12.5)') 'Aqueous Saturation = ',VAR(11)
          WRITE(56,'(4X,A,1PE12.5)') 'Ice Saturation = ',VAR(12)
!
!---      Core phase condition  ---
!
          NPHZX = MOD(IVAR(2),100)
!
!---      Phase condition #10 (10 series)  ---
!
!         SH = 0, SG + SN = 0, SL + SI = 1
!
!         energy - temperature, C - T
!         water mass - pressure, Pa (gauge) - PSO
!         mobile CO2 mass - CO2 vapor pressure, Pa absolute - PVA
!         mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!         mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!         hydrate CO2 mass - hydrate CO2 mass, kg - TMHA
!         hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!         hydrate N2 mass - hydrate N2 mass, kg - TMHN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          IF( NPHZX.EQ.10 ) THEN
            IF( IVAR(3).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'CO2 vapor pressure, Pa absolute = ',VAR(13)
            ENDIF
            IF( IVAR(4).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'CH4 vapor pressure, Pa absolute = ',VAR(14)
            ENDIF
            IF( IVAR(5).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'N2 vapor pressure, Pa absolute = ',VAR(15)
            ENDIF
            IF( IVAR(3).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)') 'CO2 Hydrate Mass = ',VAR(16)
            ENDIF
            IF( IVAR(4).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)') 'CH4 Hydrate Mass = ',VAR(17)
            ENDIF
            IF( IVAR(5).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)') 'N2 Hydrate Mass = ',VAR(18)
            ENDIF
!
!---      Phase conditions 21 through 23 (20 series)  ---
!
!         SH = 0, SG + SN > 0, SL + SI < 1
!         ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #21
!         ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #22
!         ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #23
!
!         energy - temperature, C - T
!         water mass - aqueous + ice saturation - POSM
!         maximum mobile mass component - pressure, Pa (gauge) - PSO
!         mobile CO2 mass - total nonaqueous CH4 mole fraction - ZMCA
!         mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!         mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!         hydrate CO2 mass - hydrate CO2 mass, kg - TMHA
!         hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!         hydrate N2 mass - hydrate N2 mass, kg - TMHN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          ELSEIF( NPHZX.GE.21 .AND. NPHZX.LE.23  ) THEN
            IF( IVAR(3).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Total Nonaqueous CO2 Mole Fraction = ',VAR(19)
            ENDIF
            IF( IVAR(4).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Total Nonaqueous CH4 Mole Fraction = ',VAR(20)
            ENDIF
            IF( IVAR(5).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Total Nonaqueous N2 Mole Fraction = ',VAR(21)
            ENDIF
            IF( IVAR(3).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'CO2 Hydrate Mass = ',VAR(16)
            ENDIF
            IF( IVAR(4).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'CH4 Hydrate Mass = ',VAR(17)
            ENDIF
            IF( IVAR(5).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'N2 Hydrate Mass = ',VAR(18)
            ENDIF
!
!---      Phase conditions #31 through #39 (30 series)  ---
!
!         SH > 0, SG + SN > 0, SL + SI < 1
!         ZMCA = MAX( ZMCA,ZMCO,ZMCN )
!           YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #31
!           YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #32
!           YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #33
!         ZMCO = MAX( ZMCA,ZMCO,ZMCN )
!           YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #34
!           YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #35
!           YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #36
!         ZMCN = MAX( ZMCA,ZMCO,ZMCN )
!           YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #37
!           YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #38
!           YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #39
!
!         energy - temperature, C - T
!         water mass - aqueous + ice saturation, POSM
!         maximum mobile mass component - pressure, Pa (gauge) - PSO
!         mobile CO2 mass - total nonaqueous CH4 mole fraction - ZMCA
!         mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!         mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!         maximum hydrate mass component - hydrate saturation - SH
!         hydrate CO2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!         hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!         hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          ELSEIF( NPHZX.GE.31 .AND. NPHZX.LE.39 ) THEN
            IF( IVAR(3).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Total Nonaqueous CO2 Mole Fraction = ',VAR(19)
            ENDIF
            IF( IVAR(4).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Total Nonaqueous CH4 Mole Fraction = ',VAR(20)
            ENDIF
            IF( IVAR(5).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Total Nonaqueous N2 Mole Fraction = ',VAR(21)
            ENDIF
            WRITE(56,'(4X,A,1PE12.5)')
     &        'Hydrate Saturation = ',VAR(22)
            IF( IVAR(3).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Hydrate Gas CO2 Mole Fraction of Formers = ',VAR(23)
            ENDIF
            IF( IVAR(4).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Hydrate Gas CH4 Mole Fraction of Formers = ',VAR(24)
            ENDIF
            IF( IVAR(5).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Hydrate Gas N2 Mole Fraction of Formers = ',VAR(25)
            ENDIF
!
!---      Phase conditions #41 through #43 (40 series)  ---
!
!         SH > 0, SG + SN = 0, SL + SI < 1
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #41
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #42
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #43
!
!         energy - temperature, C - T
!         water mass - pressure, Pa (gauge) - PSO
!         mobile CO2 mass - CO2 vapor pressure, Pa absolute - PVA
!         mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!         mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!         maximum hydrate mass component - hydrate saturation - SH
!         hydrate CO2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!         hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!         hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          ELSEIF( NPHZX.GE.41 .AND. NPHZX.LE.43 ) THEN
            IF( IVAR(3).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'CO2 vapor pressure, Pa absolute = ',VAR(13)
            ENDIF
            IF( IVAR(4).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'CH4 vapor pressure, Pa absolute = ',VAR(14)
            ENDIF
            IF( IVAR(5).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'N2 vapor pressure, Pa absolute = ',VAR(15)
            ENDIF
            WRITE(56,'(4X,A,1PE12.5)')
     &        'Hydrate Saturation = ',VAR(22)
            IF( IVAR(3).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Hydrate Gas CO2 Mole Fraction of Formers = ',VAR(23)
            ENDIF
            IF( IVAR(4).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Hydrate Gas CH4 Mole Fraction of Formers = ',VAR(24)
            ENDIF
            IF( IVAR(5).EQ.0 ) THEN
              WRITE(56,'(4X,A,1PE12.5)')
     &          'Hydrate Gas N2 Mole Fraction of Formers = ',VAR(25)
            ENDIF
          ENDIF
          WRITE(56,'(4X,A,1PE12.5,A,I6)')
     &      'Total-Salt Aqu. Mass Fraction = ',VAR(26)
          WRITE(56,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VAR(27),TRIM(UNTM),' to ',
     &      VAR(28),TRIM(UNTM)
        ENDIF
      ENDDO
      WRITE(56,'(A)') 'NOTE: Simulation Stopped:'
      WRITE(56,'(/,A)') '---  End of STOMP Simulation ---'
!  
!---  Deallocate memory
!
      DEALLOCATE( VAR,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Deallocation Error: VAR'
        STOP
      ENDIF
      DEALLOCATE( NDREF,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Deallocation Error: NDREF'
        STOP
      ENDIF
      DEALLOCATE( CHREF,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Deallocation Error: CHREF'
        STOP
      ENDIF
      DEALLOCATE( UNREF,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Deallocation Error: UNREF'
        STOP
      ENDIF
!
!---  Close the binary output.bin file  ---
!
      CLOSE( UNIT=55 )
!
!---  Close the text (ASCII) output file  ---
!
      CLOSE( UNIT=56 )
      STOP
!
!---  End of OUTPUT_BIN_TO_TXT program  ---
!
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CNTRTXT( TXTX,STRX,ITX,ISX )
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
!     Center text of length ITX in a string of length ISX.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 January 2004.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*(*) TXTX,STRX
!
!----------------------Executable Lines--------------------------------!
!
      IF( ITX.GT.ISX ) THEN
        PRINT *,'Text Length > String Length'
        STOP
      ENDIF
      ISRX = (ISX-ITX)/2
      ISLX = ISX-ITX-ISRX
      DO I = 1,ISLX
        STRX(I:I) = ' '
      ENDDO
      STRX(ISLX+1:ISLX+ITX) = TXTX(1:ITX)
      DO I = ISLX+ITX+1,ISLX+ITX+ISRX
        STRX(I:I) = ' '
      ENDDO
!
!---  End of CNTRTXT group  ---
!
      RETURN
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

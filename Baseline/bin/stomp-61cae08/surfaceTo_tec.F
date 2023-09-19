!----------------------Program-----------------------------------------!
!
      PROGRAM SURFACE_TO_TECPLOT
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
!     This utility routine converts surface files to Tecplot format
!     for plotting.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 February 2021.
!     surfaceTo_tec.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: VAR
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSPK
      CHARACTER(132), DIMENSION(:), ALLOCATABLE :: CHSURF,DHSURF
      CHARACTER(32), DIMENSION(:), ALLOCATABLE :: CH_HD
      CHARACTER(512) :: CHDUM
      CHARACTER(1600) :: CHREC
      CHARACTER(132) :: CDUM,OUTFILE
      CHARACTER(4) :: FORM1
      LOGICAL FCHK
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
!
!---  No surface file, exit program with error message  ---
!
      INQUIRE( FILE='surface', EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        WRITE(6,'(A)') 'ERROR: Nonexistent "surface" file.'
        STOP
      ELSE
        OPEN(UNIT=55, FILE='surface', STATUS='OLD', FORM='FORMATTED')
      ENDIF
!
!---  Get inputs from command line or by prompting the user
!
      NCLA = IARGC()
      IF( NCLA.EQ.1 ) THEN
        CALL GETARG( 1,CDUM )
        CALL LCASE( CDUM )
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts data in the ' //
     &      'surface flux data section of a STOMP surface file to '
          WRITE(6,'(A)') 'an Tecplot "Point" formatted ' //
     &      'file. The script takes one argument, the name of'
          WRITE(6,'(A)') 'Tecplot file, which should have a "dat" ' //
     &      'extension (e.g., surface.dat). If the argument is missing '
          WRITE(6,'(A)') 'the user will be prompted for the ' //
     &      'filename. The script determines the variables in the '
          WRITE(6,'(A)') 'surface file, and prompts the user for ' //
     &      'a list of variables to be included in the surface.dat '
          WRITE(6,'(A)') 'file. The variable order will be that ' //
     &      'specified by the user, which allows the variables to be '
          WRITE(6,'(A,/)') 'reordered from the order in the surface ' //
     &      'file. '
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'surfaceTo_tec.exe surface.dat'
          WRITE(6,'(A)') 'surfaceTo_tec.exe aqueous_flux.dat'
          WRITE(6,'(A)') 'surfaceTo_tec.exe heat_flux.dat'
          WRITE(6,'(A,/)') 'surfaceTo_tec.exe co2_mass_flux.dat'
          STOP
        ELSE
          OUTFILE = CDUM
        ENDIF
      ELSE
        WRITE(6,'(/,A,/)')'Enter the Tecplot filename ' // 
     &    '(e.g., surface.dat): '
        READ(5,*) OUTFILE
        CDUM = OUTFILE
        CALL LCASE( CDUM )
        IF( INDEX(CDUM(1:),'help').NE.0 ) THEN
          WRITE(6,'(/,A)') 'This script converts data in the ' //
     &      'surface flux data section of a STOMP surface file to '
          WRITE(6,'(A)') 'an Tecplot "Point"  formatted ' //
     &      'file. The script takes one argument, the name of'
          WRITE(6,'(A)') 'Tecplot file, which should have a "dat" ' //
     &      'extension (e.g., surface.dat). If the argument is missing '
          WRITE(6,'(A)') 'the user will be prompted for the ' //
     &      'filename. The script determines the variables in the '
          WRITE(6,'(A)') 'surface file, and prompts the user for ' //
     &      'a list of variables to be included in the surface.dat '
          WRITE(6,'(A)') 'file. The variable order will be that ' //
     &      'specified by the user, which allows the variables to be '
          WRITE(6,'(A,/)') 'reordered from the order in the surface ' //
     &      'file.'
          WRITE(6,'(A,/)') 'Example command line entries'
          WRITE(6,'(A)') 'surfaceTo_tec.exe surface.dat'
          WRITE(6,'(A)') 'surfaceTo_tec.exe aqueous_flux.dat'
          WRITE(6,'(A)') 'surfaceTo_tec.exe heat_flux.dat'
          WRITE(6,'(A,/)') 'surfaceTo_tec.exe co2_mass_flux.dat'
          STOP
        ENDIF
      ENDIF
!  
!---  Open surface.dat file
!
      OUTFILE = TRIM(ADJUSTL(OUTFILE))
      WRITE(6,'(2A)') 'Creating surface.dat ',TRIM(OUTFILE)
      OPEN(UNIT=56, FILE=OUTFILE)
      CLOSE(UNIT=56,STATUS='DELETE')
      OPEN(UNIT=56, FILE=OUTFILE, STATUS='NEW', FORM='FORMATTED')
!  
!---  Allocate memory
!
      ALLOCATE( CHSURF(1:100),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: CHSURF'
        STOP
      ENDIF
!
!---  Write banner to screen  ---
!
      WRITE(6,'(A,//)')' Welcome to ...'
      WRITE(6,'(A,/)') '                       SURFACE_TO_IGOR'
      WRITE(6,'(A)')   ' An output file converter routine for ' //
     &  'the STOMP simulator.  This'
      WRITE(6,'(A)')   ' routine, developed by the Pacific ' //
     &  'Northwest National Laboratory,'
      WRITE(6,'(A)')   ' reads a STOMP surface flux ' //
     &  'file and converts the data in the surface output'
      WRITE(6,'(A,/)')   ' section to a Tecplot formatted file ' //
     &  'defined by the user (e.g., surface.dat).'
      WRITE(6,'(A)')   ' For support:  Tel: 509.372.6070'
      WRITE(6,'(A,/)') '               E-mail:  mark.white@pnnl.gov'
!
!---  Search for number of surfaces  ---
!
      DO
        READ(55,'(A)', IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during surface ' //
     &      'file read: Number of Surfaces'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during surface ' //
     &      'file read: Number of Surfaces'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Number of Surfaces:')
        IF( INDX.NE.0 ) THEN
          INDX = INDEX(CHDUM(1:),':')
          I1 = INDX+1
          I2 = INDX+6
          READ(CHDUM(I1:I2),'(I6)') NSURF
          EXIT
        ENDIF
      ENDDO
!
!---  Search for surface variables  ---
!
      DO
        READ(55,'(A)',IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during surface ' //
     &      'file read: Surface Variables'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during surface ' //
     &      'file read: Surface Variables'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'Surface Variables')
        IF( INDX.NE.0 ) THEN
          EXIT
        ENDIF
      ENDDO
      DO N = 1,2
        READ(55,'(A)',IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during surface ' //
     &      'file read: Surface Variables'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during surface ' //
     &      'file read: Surface Variables'
          STOP
        ENDIF
      ENDDO
      NC = 0
      DO 
        READ(55,'(A)', IOSTAT=ICHK) CHDUM
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during surface ' //
     &      'file read: Variable List'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during surface ' //
     &      'file read: Variable List'
          STOP
        ENDIF
        INDX = INDEX(CHDUM(1:),'   ')
        IF( INDX.NE.1 ) THEN
          NC = NC + 1
          INDX = INDEX(CHDUM(1:),'  ')
          I1 = 1
          I2 = INDX-1
          READ(CHDUM(I1:I2),'(A)') CHSURF(NC)
          IF( MOD(NC,100).EQ.0 ) THEN
            DO N = 1,NC
              DHSURF(N) = CHSURF(N)
            ENDDO
            DEALLOCATE( CHSURF,STAT=ISTAT )
            IF( ISTAT.NE.0 ) THEN
              WRITE(6,'(A)') 'Deallocation Error: CHSURF'
              STOP
            ENDIF
            ALLOCATE( CHSURF(1:NC+100),STAT=ISTAT )
            IF( ISTAT.NE.0 ) THEN
              WRITE(6,'(A)') 'Allocation Error: CHSURF'
              STOP
            ENDIF
            DO N = 1,NC
              CHSURF(N) = CHSURF(N)
            ENDDO
            DO N = NC+1,NC+100
              CHSURF(N) = ''
            ENDDO
            DEALLOCATE( DHSURF,STAT=ISTAT )
            IF( ISTAT.NE.0 ) THEN
              WRITE(6,'(A)') 'Deallocation Error: DHSURF'
              STOP
            ENDIF
            ALLOCATE( DHSURF(1:NC+100),STAT=ISTAT )
            IF( ISTAT.NE.0 ) THEN
              WRITE(6,'(A)') 'Allocation Error: DHSURF'
              STOP
            ENDIF
          ENDIF
        ELSE
          EXIT
        ENDIF
      ENDDO
      NC_SFV = NC
      ALLOCATE( CH_HD(1:2*NC_SFV+1),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: CH_HD'
        STOP
      ENDIF
      DO N = 1,2*NC_SFV+1
        CH_HD(N) = ''
      ENDDO
!
!---  List surface variables and collect user input  ---
!
      ALLOCATE( NSPK(1:NC_SFV),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: NSPK'
        STOP
      ENDIF
      IF( NC_SFV.EQ.0 ) THEN
        WRITE(6,'(A)') 'No Surface Variables Found'
        STOP
      ELSE
        DO N = 1,NC_SFV
          NSPK(N) = 0
        ENDDO
        WRITE(6,'(A)') 'The STOMP surface file, "surface", contains the'
        WRITE(6,'(A,/)') 'following surface variables:'
        DO N = 1,NC_SFV
          NCH = INDEX(CHSURF(N),'  ')-1
          WRITE(6,'(I3,2A)') N,' -- ',CHSURF(N)(1:NCH)
        ENDDO
        WRITE(6,'(A,/)') '  a -- all surface flux variables'
        WRITE(6,'(A)') 'Enter the surface flux variables to include'
        WRITE(6,'(A)') 'in the "surface.dat" file, by entering a'
        WRITE(6,'(A)') 'string of comma or space delimited indices or'
        WRITE(6,'(A,/)')'entering "a" for all surface flux variables.'
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
            DO N = 1,NC_SFV
              NSPK(N) = N
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
          READ(CHDUM(I1:I2),*,IOSTAT=ICHK) NSPK(NC)
          IF( ICHK.NE.0 ) THEN
            WRITE(6,'(2A)') 'Unrecognized entry: ',CHDUM(I1:I2)
            STOP
          ENDIF
          IF( NSPK(NC).LT.1 .OR. NSPK(NC).GT.NC_SFV ) THEN
            WRITE(6,'(A,I3)') 'Out-of-range selection: ',NSPK(NC)
            STOP
          ENDIF
          DO N = 1,NC-1
            IF( NSPK(NC).EQ.NSPK(N) ) THEN
              WRITE(6,'(A,I3)') 'Duplicate selection: ',NSPK(NC)
              STOP
            ENDIF
          ENDDO
          IF( I2.GE.NCH ) EXIT
          I1 = I2+1
        ENDDO
      ENDIF
!
!---  Find the header line for the surface output  ---
!
      DO
        READ(55,'(A)',IOSTAT=ICHK) CHREC
        IF( ICHK.GT.0 ) THEN
          WRITE(6,'(A)') 'End of file encountered during surface ' //
     &      'file read: Header line'
          STOP
        ELSEIF( ICHK.LT.0 ) THEN
          WRITE(6,'(A)') 'Read error encountered during surface ' //
     &      'file read: Header line'
          STOP
        ENDIF
        INDX = INDEX(CHREC(1:),'Time')
        IF( INDX.NE.0 ) THEN
          EXIT
        ENDIF
      ENDDO
      NCHREC = 1600
      I1 = INDEX(CHREC(1:),'Time') + 4
      NC = 0
      DO
        IS = INDEX(CHREC(I1:I1),' ')
        IF( IS.NE.0 ) THEN
          I1 = I1 + 1
          IF( I1.GE.NCHREC ) EXIT
          CYCLE
        ENDIF
        IS = INDEX(CHREC(I1:),'(')
        IF( IS.NE.0 ) THEN
          I2 = I1 + IS - 2
        ELSE
          I2 = NCHREC
        ENDIF
        NC = NC + 1
        READ(CHREC(I1:I2),*,IOSTAT=ICHK) CH_HD(NC)
        IF( ICHK.NE.0 ) THEN
          WRITE(6,'(2A)') 'Header Line Read Error: ',CHREC(I1:I2)
          STOP
        ENDIF
        NCH = INDEX(CH_HD(NC),'  ')
        CH_HD(NC)(NCH:NCH) = '-'
        I1 = I2+2
        IS = INDEX(CHREC(I1:),')')
        IF( IS.NE.0 ) THEN
          I2 = I1 + IS - 2
        ELSE
          I2 = NCHREC
        ENDIF
        READ(CHREC(I1:I2),*,IOSTAT=ICHK) NS
        NCT = ICOUNT(NS)
        WRITE(FORM1(3:3),'(I1)') NCT
        NCH = INDEX(CH_HD(NC),'  ')
        WRITE(CH_HD(NC)(NCH:NCH+NCT-1),FORM1) NS
        I2 = I2+1
        IF( I2.GE.NCHREC ) EXIT
        I1 = I2+1
      ENDDO
      NHD = NC
!
!---  Write header line to surface.dat file  ---
!
      WRITE(56,'(A)') 'TITLE = "STOMP Surface File Visualization"'
      WRITE(56,'(A,$)') 'VARIABLES = '
      WRITE(56,'(A,$)') '"TM" '
      DO N = 1,NC_SFV
        IF( NSPK(N).EQ.0 ) EXIT
        NS = (N-1)*2 + 1
        NCH = INDEX(CH_HD(NS)(1:),'  ')-1
        WRITE(56,'(3A,$)') '"',CH_HD(NS)(1:NCH),'" '
        NS = NS + 1
        NCH = INDEX(CH_HD(NS)(1:),'  ')-1
        WRITE(56,'(3A,$)') '"',CH_HD(NS)(1:NCH),'" '
      ENDDO
      WRITE(56,'(A)') ''
!  
!---  Allocate memory
!
      ALLOCATE( VAR(1:NHD+1),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        WRITE(6,'(A)') 'Allocation Error: VAR'
        STOP
      ENDIF
!
!---  Search surface file for surface output records  ---
!
      DO
        READ(55,'(A)',IOSTAT=ICHK) CHREC
        IF( ICHK.NE.0 ) THEN
          EXIT
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
          READ(CHREC(I1:),*) TM,(VAR(J),J=1,NHD)
!
!---      Write variable record line to surface.dat file  ---
!
          WRITE(56,'(1PE12.5,$)') TM
          DO N = 1,NC_SFV
            IF( NSPK(N).EQ.0 ) EXIT
            NV = (NSPK(N)-1)*2 + 1
            WRITE(56,'(1X,1PE12.5,$)') VAR(NV)
            NV = NV + 1
            WRITE(56,'(1X,1PE12.5,$)') VAR(NV)
          ENDDO
          WRITE(56,'(A)') ' '
        ENDIF        
      ENDDO
!
!---  Close the surface and surface.dat files  ---
!
      CLOSE( UNIT=55 )
      CLOSE( UNIT=56 )
      STOP
!
!---  End of SURFACE_TO_TECPLOT program  ---
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

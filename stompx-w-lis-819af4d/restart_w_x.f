!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDRST_W
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
!     Water Mode (STOMP-W)
!
!     Write restart files.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 June 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PROP
      USE HYST
      USE GRID
      USE FLUX
      USE FILES
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARX
      CHARACTER*16 FN
      CHARACTER(32) :: CHMSG
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVARX
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDRST_W'
!
!---  Open restart.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'restart_0000001.bin',
     &  MPI_MODE_RDONLY,MPI_INFO_NULL,IRS,IERR )
!
!---  Set initial offset  ---
!
      IOFFSET = 0
!
!---  Allocate temporary memory for dummy integer array  ---
!
      NVAR = 16+NP
      ALLOCATE( IVARX(1:NVAR),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read integer header data (duplicated across processors)  ---
!
      OFFSET = IOFFSET
      IOFFSET = IOFFSET + NVAR*NBYTI
      CALL MPI_FILE_READ_AT( IRS,OFFSET,IVARX,NVAR,
     &   MPI_INTEGER,STATUS,IERR )
      IF( IVARX(3).NE.NFLD ) THEN
        M_ERR(1) = 'Restart Conflict: Restart File NFLD =  '
        M_ERR(2) = ''
        CALL PATH
        R_ERR = 0.D+0
        I_ERR(1) = IVARX(3)
        I_ERR(2) = 0
        I_ERR(3) = 1
        I_ERR(4) = ID
        CALL MPI_FILE_CLOSE( IRS,IERR )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IVARX(4).NE.IFLD ) THEN
        M_ERR(1) = 'Restart Conflict: Restart File IFLD =  '
        M_ERR(2) = ''
        CALL PATH
        R_ERR = 0.D+0
        I_ERR(1) = IVARX(4)
        I_ERR(2) = 0
        I_ERR(3) = 1
        I_ERR(4) = ID
        CALL MPI_FILE_CLOSE( IRS,IERR )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IVARX(5).NE.JFLD ) THEN
        M_ERR(1) = 'Restart Conflict: Restart File JFLD =  '
        M_ERR(2) = ''
        CALL PATH
        R_ERR = 0.D+0
        I_ERR(1) = IVARX(5)
        I_ERR(2) = 0
        I_ERR(3) = 1
        I_ERR(4) = ID
        CALL MPI_FILE_CLOSE( IRS,IERR )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IVARX(6).NE.KFLD ) THEN
        M_ERR(1) = 'Restart Conflict: Restart File KFLD =  '
        M_ERR(2) = ''
        CALL PATH
        R_ERR = 0.D+0
        I_ERR(1) = IVARX(6)
        I_ERR(2) = 0
        I_ERR(3) = 1
        I_ERR(4) = ID
        CALL MPI_FILE_CLOSE( IRS,IERR )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IVARX(7).NE.NPFLD ) THEN
        M_ERR(1) = 'Restart Conflict: Restart File NPFLD =  '
        M_ERR(2) = ''
        CALL PATH
        R_ERR = 0.D+0
        I_ERR(1) = IVARX(7)
        I_ERR(2) = 0
        I_ERR(3) = 1
        I_ERR(4) = ID
        CALL MPI_FILE_CLOSE( IRS,IERR )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IVARX(8).NE.IPFLD ) THEN
        M_ERR(1) = 'Restart Conflict: Restart File IPFLD =  '
        M_ERR(2) = ''
        CALL PATH
        R_ERR = 0.D+0
        I_ERR(1) = IVARX(8)
        I_ERR(2) = 0
        I_ERR(3) = 1
        I_ERR(4) = ID
        CALL MPI_FILE_CLOSE( IRS,IERR )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IVARX(9).NE.JPFLD ) THEN
        M_ERR(1) = 'Restart Conflict: Restart File JPFLD =  '
        M_ERR(2) = ''
        CALL PATH
        R_ERR = 0.D+0
        I_ERR(1) = IVARX(9)
        I_ERR(2) = 0
        I_ERR(3) = 1
        I_ERR(4) = ID
        CALL MPI_FILE_CLOSE( IRS,IERR )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IVARX(10).NE.JPFLD ) THEN
        M_ERR(1) = 'Restart Conflict: Restart File KPFLD =  '
        M_ERR(2) = ''
        CALL PATH
        R_ERR = 0.D+0
        I_ERR(1) = IVARX(10)
        I_ERR(2) = 0
        I_ERR(3) = 1
        I_ERR(4) = ID
        CALL MPI_FILE_CLOSE( IRS,IERR )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      NRIMX = IVARX(1)
      NSTEP = IVARX(2)
      NFLD = IVARX(3)
      IFLD = IVARX(4)
      JFLD = IVARX(5)
      KFLD = IVARX(6)
      NPFLD = IVARX(7)
      IPFLD = IVARX(8)
      JPFLD = IVARX(9)
      KPFLD = IVARX(10)
      NP = IVARX(11)
      DO MP = 1,NP
        NFCGC(MP) = IVARX(11+MP)
      ENDDO
      NSOLU = IVARX(12+NP)
      IF( IVARX(13+NP).NE.IOM ) THEN
        M_ERR(1) = 'Restart Conflict: Restart File IOM =  '
        M_ERR(2) = ''
        CALL PATH
        R_ERR = 0.D+0
        I_ERR(1) = IVARX(13+NP)
        I_ERR(2) = 0
        I_ERR(3) = 1
        I_ERR(4) = ID
        CALL MPI_FILE_CLOSE( IRS,IERR )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IOM = IVARX(13+NP)
      ISLC(5) = IVARX(14+NP)
      NSPR = IVARX(15+NP)
      NSPS = IVARX(16+NP)
!      PRINT *,'IVARX = ',(IVARX(M),M=1,NVAR),' ID = ',ID
!
!---  Deallocate temporary memory for dummy integer array  ---
!
      DEALLOCATE( IVARX,STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate temporary memory for dummy real array  ---
!
      ALLOCATE( VARX(1:7),STAT=ISTAT )
      CHMSG = 'VARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read real header data (duplicated across processors)  ---
!
      NVAR = 7
      OFFSET = IOFFSET
      IOFFSET = IOFFSET + NVAR*NBYTR
!      PRINT *,'OFFSET = ',OFFSET,' ID = ',ID
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NVAR,
     &   MPI_REAL8,STATUS,IERR )
      TM = VARX(1)
      DT = VARX(2)
      DTMX = VARX(3)
      DTMN = VARX(4)
      DTAF = VARX(5)
      RSDMX = VARX(6)
      DTCF = VARX(7)
!      PRINT *,'VARX = ',(VARX(M),M=1,NVAR),' ID = ',ID
!
!---  Deallocate temporary memory for dummy real array  ---
!
      DEALLOCATE( VARX,STAT=ISTAT )
      CHMSG = 'VARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Set local starting point for local copies of nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFCGC(I)
      ENDDO
!
!---  Allocate local temporary state condition arrays 
!     (including ghost cells)  ---
!
      ALLOCATE( VARX(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'VARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Read local copies of temperature 
!     (including ghost cells)  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        T(2,N) = VARX(N)
!        IF( N.EQ.11 ) PRINT *,'T(2,',ND(N),') = ',T(2,N),' ID = ',ID
      ENDDO
!
!---  Read local copies of aqueous pressure 
!     (including ghost cells)  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        PL(2,N) = VARX(N)
!        IF( N.EQ.11 ) PRINT *,'PL(2,',ND(N),') = ',PL(2,N),' ID = ',ID
      ENDDO
!
!---  Read local copies of gas pressure 
!     (including ghost cells)  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        PG(2,N) = VARX(N)
!        IF( N.EQ.11 ) PRINT *,'PG(2,',ND(N),') = ',PG(2,N),' ID = ',ID
      ENDDO
!
!---  Read local copies of aqueous saturation 
!     (including ghost cells)  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        SL(2,N) = VARX(N)
!        IF( N.EQ.11 ) PRINT *,'SL(2,',ND(N),') = ',SL(2,N),' ID = ',ID
      ENDDO
!
!---  Read local copies of trapped gas saturation 
!     (including ghost cells)  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        SGT(2,N) = VARX(N)
!        IF( N.EQ.11 ) PRINT *,'SGT(2,',ND(N),') = ',SGT(2,N),' ID = ',ID
      ENDDO
!
!---  Read local copies of minimum apparent aqueous saturation 
!     (including ghost cells)  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        ASLMIN(2,N) = VARX(N)
!        IF( N.EQ.11 ) PRINT *,'ASLMIN(2,',ND(N),') = ',ASLMIN(2,N),
!     &    ' ID = ',ID
      ENDDO
!
!---  Read local copies of compressibility reference pressure
!     (including ghost cells)  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        PCMP(N) = VARX(N)
!        IF( N.EQ.11 ) PRINT *,'PCMP(',ND(N),') = ',PCMP(N),' ID = ',ID
      ENDDO
!
!---  Read local copies of compressibility reference temperature
!     (including ghost cells)  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        TCMP(N) = VARX(N)
!        IF( N.EQ.11 ) PRINT *,'TCMP(',ND(N),') = ',TCMP(N),' ID = ',ID
      ENDDO
!
!---  Read local copies of phase condition as a real
!     (including ghost cells)  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        NPHAZ(2,N) = INT( VARX(N) )
!        IF( N.EQ.11 ) PRINT *,'NPHAZ(2,',ND(N),') = ',NPHAZ(2,N),
!     &    ' ID = ',ID
      ENDDO
!
!---  Read local copies of drainage/imbibition index as a real
!     (including ghost cells)  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &  STATUS,IERR)
      DO N = 1,NFCGC(ID+1)
        IPH(2,N) = INT( VARX(N) )
!        IF( N.EQ.11 ) PRINT *,'NPHAZ(2,',ND(N),') = ',NPHAZ(2,N),
!     &    ' ID = ',ID
      ENDDO
!
!---  Loop over number of solutes  ---
!
      DO NSL = 1,NSOLU
!
!---    Read local copies of solute volumetric concentration
!       (including ghost cells)  ---
!
        OFFSET = IOFFSET + NC*NBYTR
        IOFFSET = IOFFSET + NFCGC_G*NBYTR
        CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR)
        DO N = 1,NFCGC(ID+1)
          C(N,NSL) = VARX(N)
!          IF( N.EQ.11 ) PRINT *,'C(',ND(N),',',NSL,') = ',C(N,NSL),
!     &      ' ID = ',ID
        ENDDO
      ENDDO
!---  Loop over number of reactive species  ---
!
      DO NSP = 1,NSPR
!
!---    Read local copies of species volumetric concentration
!       (including ghost cells)  ---
!
        OFFSET = IOFFSET + NC*NBYTR
        IOFFSET = IOFFSET + NFCGC_G*NBYTR
        CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR)
        DO N = 1,NFCGC(ID+1)
          SP_C(N,NSP) = VARX(N)
!          IF( N.EQ.11 ) PRINT *,'SP_C(',ND(N),',',NSP,') = ',
!     &      SP_C(N,NSP),' ID = ',ID
        ENDDO
      ENDDO
!
!
!---  Loop over number of mineral (solid) reactive species  ---
!
      DO NSP = 1,NSPS
!
!---    Read local copies of base mineral species 
!       volumetric concentration (including ghost cells)  ---
!
        OFFSET = IOFFSET + NC*NBYTR
        IOFFSET = IOFFSET + NFCGC_G*NBYTR
        CALL MPI_FILE_READ_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR)
        DO N = 1,NFCGC(ID+1)
          SP_CMN(N,NSP) = VARX(N)
!          IF( N.EQ.11 ) PRINT *,'SP_CMN(',ND(N),',',NSP,') = ',
!     &      SP_CMN(N,NSP),' ID = ',ID
        ENDDO
      ENDDO
!
!---  Deallocate local temporary state condition arrays  ---
!
      DEALLOCATE( VARX,STAT=ISTAT )
      CHMSG = 'VARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Close restart.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRS,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDRST_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRRST_W
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
!     Water Mode (STOMP-W)
!
!     Write restart files.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 June 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PROP
      USE HYST
      USE GRID
      USE FLUX
      USE FILES
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: VARX
      CHARACTER*19 FN
      CHARACTER(32) :: CHMSG
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVARX
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRRST_W'
!
!---  Create a new restart file with number of time steps
!     as the file name extension  ---
!
      FN = 'restart_0000000.bin'
      IF( NSTEP.LE.9 ) THEN
        FN(15:15) = CHAR(48+NSTEP)
      ELSEIF( NSTEP.LE.99 ) THEN
        I1 = INT(NSTEP/10)
        I2 = MOD( NSTEP,10 )
        FN(14:14) = CHAR(48+I1)
        FN(15:15) = CHAR(48+I2)
      ELSEIF( NSTEP.LE.999 ) THEN
        I1 = INT(NSTEP/100)
        I2 = INT( MOD( NSTEP,100 )/10 )
        I3 = MOD( NSTEP,10 )
        FN(13:13) = CHAR(48+I1)
        FN(14:14) = CHAR(48+I2)
        FN(15:15) = CHAR(48+I3)
      ELSEIF( NSTEP.LE.9999 ) THEN
        I1 = INT(NSTEP/1000)
        I2 = INT( MOD( NSTEP,1000 )/100 )
        I3 = INT( MOD( NSTEP,100 )/10 )
        I4 = MOD( NSTEP,10 )
        FN(12:12) = CHAR(48+I1)
        FN(13:13) = CHAR(48+I2)
        FN(14:14) = CHAR(48+I3)
        FN(15:15) = CHAR(48+I4)
      ELSEIF( NSTEP.LE.99999 ) THEN
        I1 = INT(NSTEP/10000)
        I2 = INT( MOD( NSTEP,10000 )/1000 )
        I3 = INT( MOD( NSTEP,1000 )/100 )
        I4 = INT( MOD( NSTEP,100 )/10 )
        I5 = MOD( NSTEP,10 )
        FN(11:11) = CHAR(48+I1)
        FN(12:12) = CHAR(48+I2)
        FN(13:13) = CHAR(48+I3)
        FN(14:14) = CHAR(48+I4)
        FN(15:15) = CHAR(48+I5)
      ELSEIF( NSTEP.LE.999999 ) THEN
        I1 = INT(NSTEP/100000)
        I2 = INT( MOD( NSTEP,100000 )/10000 )
        I3 = INT( MOD( NSTEP,10000 )/1000 )
        I4 = INT( MOD( NSTEP,1000 )/100 )
        I5 = INT( MOD( NSTEP,100 )/10 )
        I6 = MOD( NSTEP,10 )
        FN(10:10) = CHAR(48+I1)
        FN(11:11) = CHAR(48+I2)
        FN(12:12) = CHAR(48+I3)
        FN(13:13) = CHAR(48+I4)
        FN(14:14) = CHAR(48+I5)
        FN(15:15) = CHAR(48+I6)
      ELSEIF( NSTEP.LE.9999999 ) THEN
        I1 = INT(NSTEP/1000000)
        I2 = INT( MOD( NSTEP,1000000 )/100000 )
        I3 = INT( MOD( NSTEP,100000 )/10000 )
        I4 = INT( MOD( NSTEP,10000 )/1000 )
        I5 = INT( MOD( NSTEP,1000 )/100 )
        I6 = INT( MOD( NSTEP,100 )/10 )
        I7 = MOD( NSTEP,10 )
        FN(9:9) = CHAR(48+I1)
        FN(10:10) = CHAR(48+I2)
        FN(11:11) = CHAR(48+I3)
        FN(12:12) = CHAR(48+I4)
        FN(13:13) = CHAR(48+I5)
        FN(14:14) = CHAR(48+I6)
        FN(15:15) = CHAR(48+I7)
      ELSE
        FN = 'restart_xxxxxxx.bin'
      ENDIF
!
!---  Open restart_xxxxxxx.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,FN,
     &  MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,IRS,IERR )
!
!---  Set initial offset  ---
!
      IOFFSET = 0
!
!---  Allocate temporary memory for dummy integer array  ---
!
      NVAR = 16+NP
      ALLOCATE( IVARX(1:NVAR),STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Write integer header data  ---
!
      IVARX(1) = NRIMX
      IVARX(2) = NSTEP
      IVARX(3) = NFLD
      IVARX(4) = IFLD
      IVARX(5) = JFLD
      IVARX(6) = KFLD
      IVARX(7) = NPFLD
      IVARX(8) = IPFLD
      IVARX(9) = JPFLD
      IVARX(10) = KPFLD
      IVARX(11) = NP
      DO MP = 1,NP
       IVARX(11+MP) = NFCGC(MP)
      ENDDO
      IVARX(12+NP) = NSOLU
      IVARX(13+NP) = IOM
      IVARX(14+NP) = ISLC(5)
      IVARX(15+NP) = NSPR
      IVARX(16+NP) = NSPS
      OFFSET = IOFFSET
      IOFFSET = IOFFSET + NVAR*NBYTI
!      PRINT *,'1: OFFSET = ',OFFSET,' IOFFSET = ',IOFFSET,' ID = ',ID
      IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IRS,OFFSET,IVARX,NVAR,
     &   MPI_INTEGER,STATUS,IERR )
!      PRINT *,'IVARX = ',(IVARX(M),M=1,NVAR),' ID = ',ID
!
!---  Deallocate temporary memory for dummy integer array  ---
!
      DEALLOCATE( IVARX,STAT=ISTAT )
      CHMSG = 'IVARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate temporary memory for dummy real array  ---
!
      ALLOCATE( VARX(1:7),STAT=ISTAT )
      CHMSG = 'VARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Write real header data  ---
!
      NVAR = 7
      VARX(1) = TM
      VARX(2) = DT
      VARX(3) = DTMX
      VARX(4) = DTMN
      VARX(5) = DTAF
      VARX(6) = RSDMX
      VARX(7) = DTCF
      OFFSET = IOFFSET
      IOFFSET = IOFFSET + NVAR*NBYTR
!      PRINT *,'2: OFFSET = ',OFFSET,' IOFFSET = ',IOFFSET,' ID = ',ID
      IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NVAR,
     &   MPI_REAL8,STATUS,IERR )
!      PRINT *,'VARX = ',(VARX(M),M=1,NVAR),' ID = ',ID
!
!---  Deallocate temporary memory for dummy real array  ---
!
      DEALLOCATE( VARX,STAT=ISTAT )
      CHMSG = 'VARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Allocate temporary memory for restart file variables  ---
!
      ALLOCATE( VARX(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'VARX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Set local starting point for local copies of nodal variables  ---
!
      NC = 0
      DO I = 1,ID
        NC = NC + NFCGC(I)
      ENDDO
!
!---  Load temporary array with temperature  ---
!
      DO N = 1,NFCGC(ID+1)
        VARX(N) = T(2,N)
      ENDDO
!
!---  Write temperature to restart file  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR )
!
!---  Load temporary array with aqueous pressure  ---
!
      DO N = 1,NFCGC(ID+1)
        VARX(N) = PL(2,N)
      ENDDO
!
!---  Write aqueous pressure to restart file  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR )
!
!---  Load temporary array with gas pressure  ---
!
      DO N = 1,NFCGC(ID+1)
        VARX(N) = PG(2,N)
      ENDDO
!
!---  Write gas pressure to restart file  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR )
!
!---  Load temporary array with aqueous saturation  ---
!
      DO N = 1,NFCGC(ID+1)
        VARX(N) = SL(2,N)
      ENDDO
!
!---  Write aqueous saturation to restart file  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR )
!
!---  Load temporary array with trapped gas saturation  ---
!
      DO N = 1,NFCGC(ID+1)
        VARX(N) = SGT(2,N)
      ENDDO
!
!---  Write trapped gas saturation to restart file  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR )
!
!---  Load temporary array with minimum apparent aqueous saturation  ---
!
      DO N = 1,NFCGC(ID+1)
        VARX(N) = ASLMIN(2,N)
      ENDDO
!
!---  Write minimum apparent aqueous saturation to restart file  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR )
!
!---  Load temporary array with compressibility reference pressure   ---
!
      DO N = 1,NFCGC(ID+1)
        VARX(N) = PCMP(N)
      ENDDO
!
!---  Write compressibility reference pressure to restart file  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR )
!
!---  Load temporary array with compressibility reference 
!     temperature   ---
!
      DO N = 1,NFCGC(ID+1)
        VARX(N) = TCMP(N)
      ENDDO
!
!---  Write compressibility reference temperature to restart file  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR )
!
!---  Load temporary array with phase condition as a real   ---
!
      DO N = 1,NFCGC(ID+1)
        VARX(N) = REAL(NPHAZ(2,N))
      ENDDO
!
!---  Write phase condition as a real to restart file  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR )
!
!---  Load temporary array with draining/imbibition index as a 
!     real   ---
!
      DO N = 1,NFCGC(ID+1)
        VARX(N) = REAL(IPH(2,N))
      ENDDO
!
!---  Write draining/imbibition index as a real to restart file  ---
!
      OFFSET = IOFFSET + NC*NBYTR
      IOFFSET = IOFFSET + NFCGC_G*NBYTR
      CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &    STATUS,IERR )
!
!---  Loop over number of solutes  ---
!
      DO NSL = 1,NSOLU
!
!---    Load temporary array with solute volumetric concentration   ---
!
        DO N = 1,NFCGC(ID+1)
          VARX(N) = C(N,NSL)
        ENDDO
!
!---    Write solute volumetric concentration to restart file  ---
!
        OFFSET = IOFFSET + NC*NBYTR
        IOFFSET = IOFFSET + NFCGC_G*NBYTR
!     &    ' ID = ',ID
        CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &      STATUS,IERR )
      ENDDO
!
!---  Loop over number of reactive species  ---
!
      DO NSP = 1,NSPR
!
!---    Load temporary array with species volumetric concentration   ---
!
        DO N = 1,NFCGC(ID+1)
          VARX(N) = SP_C(N,NSP)
        ENDDO
!
!---    Write species volumetric concentration to restart file  ---
!
        OFFSET = IOFFSET + NC*NBYTR
        IOFFSET = IOFFSET + NFCGC_G*NBYTR
!     &    ' ID = ',ID
        CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &      STATUS,IERR )
      ENDDO
!
!---  Loop over number of mineral (solid) reactive species  ---
!
      DO NSP = 1,NSPS
!
!---    Load temporary array with base mineral species 
!       volumetric concentration   ---
!
        DO N = 1,NFCGC(ID+1)
          VARX(N) = SP_CMN(N,NSP)
        ENDDO
!
!---    Write base mineral species 
!       volumetric concentration to restart file  ---
!
        OFFSET = IOFFSET + NC*NBYTR
        IOFFSET = IOFFSET + NFCGC_G*NBYTR
!     &    ' IOFFSET = ',IOFFSET,' ID = ',ID
        CALL MPI_FILE_WRITE_AT( IRS,OFFSET,VARX,NFCGC(ID+1),MPI_REAL8,
     &      STATUS,IERR )
      ENDDO
!
!---  Deallocate temporary memory for restart file variables  ---
!
      DEALLOCATE( VARX,STAT=ISTAT )
      CHMSG = 'VARX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Close restart_xxxxxxx.bin file  ---
!
      CALL MPI_FILE_CLOSE( IRS,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRRST_W group  ---
!
      RETURN
      END


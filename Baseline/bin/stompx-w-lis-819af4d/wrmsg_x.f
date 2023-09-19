!----------------------Function----------------------------------------!
!
      SUBROUTINE ALLOC_ERROR( CHMSG,ISTAT )
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
!     Allocation error occurs print message and abort simulation
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 6 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GRID
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER(32) :: CHMSG
      CHARACTER(64) :: SPATH
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_ERROR'
      CALL MPI_ALLREDUCE( ISTAT,ISTATX,1,MPI_INTEGER,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
      IF( ISTATX.NE.0 ) THEN
        IF( ID.EQ.0 ) THEN
          SPATH(1:6) = 'Path: '
          IC1 = 7
          DO I = 1,ISUB_LOG
            ICSN = INDEX( SUB_LOG(I)(1:),'  ' ) - 1
            IC2 = IC1 + ICSN - 1
            SPATH(IC1:IC2) = SUB_LOG(I)(1:ICSN)
            IC1 = IC2 + 1
          ENDDO
          PRINT *,'Allocation Error for ',TRIM(CHMSG)
          PRINT *,SPATH(1:IC2)
        ENDIF
        CALL MPI_FINALIZE(IERR)
        STOP
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_ERROR group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      SUBROUTINE CHK_ERROR
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
!     Check for fatal error messages and stop simulation if indicated.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GRID
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER, DIMENSION(1:4) :: NCH
      CHARACTER(32) :: CHMSG
      CHARACTER(64) :: SPATH
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK_ERROR'
!
!---  Check across processors for a processor rank other than
!     the initial value of the number of processors plus one   ---
!
      CALL MPI_ALLREDUCE( I_ERR(4),ID_ERR,1,MPI_INTEGER,MPI_MIN,
     &  MPI_COMM_WORLD,IERR )
!
!---  If an active processor rank is detected then broadcast the
!     error messaging details and print the error message from P0
!     and stop the simulation  ---
!
      IF( ID_ERR.GE.0 .AND. ID_ERR.LE.NP ) THEN
        CALL MPI_BCAST( M_ERR,4*128,MPI_CHAR,ID_ERR,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST( I_ERR,4,MPI_INTEGER,ID_ERR,MPI_COMM_WORLD,IERR )
        CALL MPI_BCAST( R_ERR,1,MPI_REAL8,ID_ERR,MPI_COMM_WORLD,IERR )
        IF( ID.EQ.0 ) THEN
          DO M = 1,4
            NCH(M) = INDEX(M_ERR(M),'   ') - 1
          ENDDO
!
!---      No real output with error message  ---
!
          IF( I_ERR(2).EQ.0 ) THEN
!
!---        No integer output with error message  ---
!
            IF( I_ERR(3).EQ.0 ) THEN
              PRINT *,'Error: ',(M_ERR(I)(1:NCH(I)),I=1,3)
!
!---        Integer output between error string 1 and 2  ---
!
            ELSEIF( I_ERR(3).EQ.1 ) THEN
              PRINT *,'Error: ',M_ERR(1)(1:NCH(1)),I_ERR(1),
     &          (M_ERR(I)(1:NCH(I)),I=2,3)
!
!---        Integer output between error string 2 and 3  ---
!
            ELSEIF( I_ERR(3).EQ.2 ) THEN
              PRINT *,'Error: ',(M_ERR(I)(1:NCH(I)),I=1,2),I_ERR(1),
     &          (M_ERR(I)(1:NCH(I)),I=3,3)
!
!---        Integer output between error string 3 and 4  ---
!
            ELSEIF( I_ERR(3).EQ.3 ) THEN
              PRINT *,'Error: ',(M_ERR(I)(1:NCH(I)),I=1,3),I_ERR(1)
            ENDIF
!
!---      Real output between error string 1 and 2  ---
!
          ELSEIF( I_ERR(2).EQ.1 ) THEN
!
!---        No integer output with error message  ---
!
            IF( I_ERR(3).EQ.0 ) THEN
              PRINT *,'Error: ',M_ERR(1)(1:NCH(1)),R_ERR,
     &          (M_ERR(I)(1:NCH(I)),I=2,3)
!
!---        Integer output between error string 2 and 3  ---
!
            ELSEIF( I_ERR(3).EQ.2 ) THEN
              PRINT *,'Error: ',(M_ERR(I)(1:NCH(I)),I=1,1),R_ERR,
     &          (M_ERR(I)(1:NCH(I)),I=2,2),I_ERR(1),
     &          (M_ERR(I)(1:NCH(I)),I=3,3)
!
!---        Integer output between error string 3 and 4  ---
!
            ELSEIF( I_ERR(3).EQ.3 ) THEN
              PRINT *,'Error: ',M_ERR(1)(1:NCH(1)),R_ERR,
     &          (M_ERR(I)(1:NCH(I)),I=2,3),I_ERR(1)
            ENDIF
!
!---      Real output between error string 2 and 3  ---
!
          ELSEIF( I_ERR(2).EQ.2 ) THEN
!
!---        No integer output with error message  ---
!
            IF( I_ERR(3).EQ.0 ) THEN
              PRINT *,'Error: ',(M_ERR(I)(1:NCH(I)),I=1,2),R_ERR,
     &          (M_ERR(I)(1:NCH(I)),I=3,3)
!
!---        Integer output between error string 1 and 2  ---
!
            ELSEIF( I_ERR(3).EQ.1 ) THEN
              PRINT *,'Error: ',M_ERR(1)(1:NCH(1)),I_ERR(1),
     &          M_ERR(2)(1:NCH(2)),R_ERR,
     &          M_ERR(3)(1:NCH(3))
!
!---        Integer output between error string 3 and 4  ---
!
            ELSEIF( I_ERR(3).EQ.3 ) THEN
              PRINT *,'Error: ',(M_ERR(I)(1:NCH(I)),I=1,2),R_ERR,
     &          M_ERR(3)(1:NCH(3)),I_ERR(1)
            ENDIF
!
!---      Real output between error string 3 and 4  ---
!
          ELSEIF( I_ERR(2).EQ.3 ) THEN
!
!---        No integer output with error message  ---
!
            IF( I_ERR(3).EQ.0 ) THEN
              PRINT *,'Error: ',(M_ERR(I)(1:NCH(I)),I=1,3),R_ERR
!
!---        Integer output between error string 1 and 2  ---
!
            ELSEIF( I_ERR(3).EQ.1 ) THEN
              PRINT *,'Error: ',M_ERR(1)(1:NCH(1)),I_ERR(1),
     &          (M_ERR(I)(1:NCH(I)),I=2,3),R_ERR
!
!---        Integer output between error string 2 and 3  ---
!
            ELSEIF( I_ERR(3).EQ.2 ) THEN
              PRINT *,'Error: ',(M_ERR(I)(1:NCH(I)),I=1,2),I_ERR(1),
     &          M_ERR(3)(1:NCH(3)),R_ERR
            ENDIF
          ENDIF
          PRINT *,M_ERR(4)(1:NCH(4))
          PRINT *,'Processor: ',I_ERR(4)
        ENDIF
        CALL MPI_FINALIZE(IERR)
        STOP
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_ERROR group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      SUBROUTINE DEALLOC_ERROR( CHMSG,ISTAT )
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
!     Deallocation error occurs print message and abort simulation
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 6 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GRID
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER(32) :: CHMSG
      CHARACTER(64) :: SPATH
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DEALLOC_ERROR'
      CALL MPI_ALLREDUCE( ISTAT,ISTATX,1,MPI_INTEGER,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
      IF( ISTATX.NE.0 ) THEN
        IF( ID.EQ.0 ) THEN
          SPATH(1:6) = 'Path: '
          IC1 = 7
          DO I = 1,ISUB_LOG
            ICSN = INDEX( SUB_LOG(I)(1:),'  ' ) - 1
            IC2 = IC1 + ICSN - 1
            SPATH(IC1:IC2) = SUB_LOG(I)(1:ICSN)
            IC1 = IC2 + 1
          ENDDO
          PRINT *,'Deallocation Error for ',TRIM(CHMSG)
          PRINT *,SPATH(1:IC2)
        ENDIF
        CALL MPI_FINALIZE(IERR)
        STOP
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DEALLOC_ERROR group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      SUBROUTINE PATH
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
!     Load message error string with subroutine log path
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 6 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      M_ERR(4)(1:6) = 'Path: '
      IC1 = 7
      DO I = 1,ISUB_LOG
        ICSN = INDEX( SUB_LOG(I)(1:),'  ' ) - 1
        IC2 = IC1 + ICSN - 1
        M_ERR(4)(IC1:IC2) = SUB_LOG(I)(1:ICSN)
        IC1 = IC2 + 1
      ENDDO
!
!---  End of PATH group
!
      RETURN
      END


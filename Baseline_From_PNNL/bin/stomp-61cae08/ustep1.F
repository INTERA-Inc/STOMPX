!----------------------Program-----------------------------------------!
!
      PROGRAM USTEP
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
!     USTEP: UCode Preprocessor for STOMP
!
!     Water Mode
!
!     This program reads a STOMP input file and generates two control
!     files for UCode.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 June 2001.
!     Last Modified by MD White, PNNL, 15 June 2001.
!     ustep1.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = 1
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = 'USTEP1'
      ICODE = 1
!
!---  Intialize variables in common blocks and open files  ---
!
      CALL INTLZ
!
!---  Search input file for grid card  ---
!
   10 READ(IRD,'(A)', END=20) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 10
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'grid').NE.0 ) THEN
        CALL RDGRID
        REWIND(IRD)
        GOTO 100
      ELSE
        GOTO 10
      ENDIF
   20 CONTINUE
      INDX = 18
      CHMSG = 'Missing Grid Card'
      CALL WRMSGS( INDX )
  100 CONTINUE
!
!---  Search input file for inactive nodes card  ---
!
  110 READ(IRD,'(A)', END=120) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 110
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'inactive').NE.0 ) THEN
        CALL RDINAC
        REWIND(IRD)
        GOTO 200
      ELSE
        GOTO 100
      ENDIF
  120 CONTINUE
      INDX = 1
      CHMSG = 'Missing Inactive Nodes Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
  200 CONTINUE
!
!---  Search input file for observed data card --
!
  310 READ(IRD,'(A)', END=320) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 310
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'observed').NE.0 ) THEN
        CALL RDOBDA
        REWIND(IRD)
        GOTO 400
      ELSE
        GOTO 310
      ENDIF
  320 CONTINUE
      INDX = 18
      CHMSG = 'Missing Observed Data Card'
      CALL WRMSGS( INDX )
  400 CONTINUE
!
!---  Search input file for UCode control card --
!
  410 READ(IRD,'(A)', END=420) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 410
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'ucode').NE.0 ) THEN
        CALL RDUCODE
        REWIND(IRD)
        GOTO 500
      ELSE
        GOTO 410
      ENDIF
  420 CONTINUE
      INDX = 18
      CHMSG = 'Missing UCode Control Card'
      CALL WRMSGS( INDX )
  500 CONTINUE
!
!---  Write data to UCode Control Files --
!
      CALL WROBDA
!
!---  End of USTEP program  ---
!
      STOP
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINIGC( DFV,DVAR,DADD,IDOM,MM,IGC )
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
!     Load component variables w/ initial conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 March 2008.
!     Last Modified by MD White, PNNL, 31 March 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER IDOM(6)
      REAL*8 DVAR(*),DFV(LNGC,LSV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINIC'
      IDOM(1) = MAX( 1,IDOM(1) )
      IDOM(1) = MIN( IDOM(1),IDOM(2),IFLD )
      IDOM(2) = MAX( 1,IDOM(1),IDOM(2) )
      IDOM(2) = MIN( IDOM(2),IFLD )
      IDOM(3) = MAX( 1,IDOM(3) )
      IDOM(3) = MIN( IDOM(3),IDOM(4),JFLD )
      IDOM(4) = MAX( 1,IDOM(3),IDOM(4) )
      IDOM(4) = MIN( IDOM(4),JFLD )
      IDOM(5) = MAX( 1,IDOM(5) )
      IDOM(5) = MIN( IDOM(5),IDOM(6),KFLD )
      IDOM(6) = MAX( 1,IDOM(5),IDOM(6) )
      IDOM(6) = MIN( IDOM(6),KFLD )
      WRITE (IWR,'(4X,A)') 'Domain: '
      WRITE (IWR,'(6X,2(A,I6))') 'I = ',IDOM(1),' to ',IDOM(2)
      WRITE (IWR,'(6X,2(A,I6))') 'J = ',IDOM(3),' to ',IDOM(4)
      WRITE (IWR,'(6X,2(A,I6))') 'K = ',IDOM(5),' to ',IDOM(6)
      IF( IDOM(1) .GT. IDOM(2) .OR. IDOM(3) .GT. IDOM(4) .OR.
     &  IDOM(5) .GT. IDOM(6) ) THEN
        NCH = INDEX( VARB(1:),'  ' )-1
        CHMSG = 'Invalid Initial Condition Domain: '//VARB(1:NCH)
        INDX = 4
        CALL WRMSGS( INDX )
      ENDIF
      DO K = IDOM(5), IDOM(6)
      DO J = IDOM(3), IDOM(4)
      DO I = IDOM(1), IDOM(2)
        N = ND(I,J,K)
        M = ND(IDOM(1),IDOM(3),IDOM(5))
        IF( IXP(N).EQ.0 ) CYCLE
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(IGC,MM,NRX) = DVAR(1) + DADD +
     &      (XP(NRX) - XP(M))*DVAR(2) +
     &      (YP(NRX) - YP(M))*DVAR(3) +
     &      (ZP(NRX) - ZP(M))*DVAR(4)
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINIC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINIP( DFV,DVAR,DADD,IFV,IVAR,IDOM )
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
!     Load primary variables w/ initial conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 December 1992.
!     Last Modified by MD White, PNNL, 31 December 1992.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER IDOM(6),IVAR,IFV(*)
      REAL*8 DVAR(*),DFV(*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINIP'
      IDOM(1) = MAX( 1,IDOM(1) )
      IDOM(1) = MIN( IDOM(1),IDOM(2),IFLD )
      IDOM(2) = MAX( 1,IDOM(1),IDOM(2) )
      IDOM(2) = MIN( IDOM(2),IFLD )
      IDOM(3) = MAX( 1,IDOM(3) )
      IDOM(3) = MIN( IDOM(3),IDOM(4),JFLD )
      IDOM(4) = MAX( 1,IDOM(3),IDOM(4) )
      IDOM(4) = MIN( IDOM(4),JFLD )
      IDOM(5) = MAX( 1,IDOM(5) )
      IDOM(5) = MIN( IDOM(5),IDOM(6),KFLD )
      IDOM(6) = MAX( 1,IDOM(5),IDOM(6) )
      IDOM(6) = MIN( IDOM(6),KFLD )
      WRITE(IWR,'(4X,A)') 'Domain: '
      WRITE (IWR, '(6X,2(A,I6))') 'I = ',IDOM(1),  ' to ',IDOM(2)
      WRITE (IWR, '(6X,2(A,I6))') 'J = ',IDOM(3),  ' to ',IDOM(4)
      WRITE (IWR, '(6X,2(A,I6))') 'K = ',IDOM(5),  ' to ',IDOM(6)
      IF( IDOM(1) .GT. IDOM(2) .OR. IDOM(3) .GT. IDOM(4) .OR.
     &  IDOM(5) .GT. IDOM(6) ) THEN
        NCH = INDEX( VARB(1:),'  ' )-1
        CHMSG = 'Invalid Initial Condition Domain: '//VARB(1:NCH)
        INDX = 4
        CALL WRMSGS( INDX )
      ENDIF
!
!---  pH  ---
!
      IF( ABS(DADD-7.D+0)/EPSL.LT.EPSL ) THEN
        DO K = IDOM(5), IDOM(6)
        DO J = IDOM(3), IDOM(4)
        DO I = IDOM(1), IDOM(2)
          N = ND(I,J,K)
          M = ND(IDOM(1),IDOM(3),IDOM(5))
          IF( IXP(N).EQ.0 ) CYCLE
          DO NRX = IBR(4,N),IBR(5,N)
            DFV(NRX) = DVAR(1) +
     &        (XP(NRX) - XP(M))*DVAR(2) +
     &        (YP(NRX) - YP(M))*DVAR(3) +
     &        (ZP(NRX) - ZP(M))*DVAR(4)
            DFV(NRX) = MIN( MAX( DFV(NRX),0.D+0 ), 1.4D+1 )
            DFV(NRX) = 1.D+1**(-DFV(NRX))
            IFV(NRX) = IVAR
          ENDDO
        ENDDO
        ENDDO
        ENDDO
      ELSE
        DO K = IDOM(5), IDOM(6)
        DO J = IDOM(3), IDOM(4)
        DO I = IDOM(1), IDOM(2)
          N = ND(I,J,K)
          M = ND(IDOM(1),IDOM(3),IDOM(5))
          IF( IXP(N).EQ.0 ) CYCLE
          DO NRX = IBR(4,N),IBR(5,N)
            DFV(NRX) = DVAR(1) + DADD +
     &        (XP(NRX) - XP(M))*DVAR(2) +
     &        (YP(NRX) - YP(M))*DVAR(3) +
     &        (ZP(NRX) - ZP(M))*DVAR(4)
            IFV(NRX) = IVAR
          ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINIP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINIR( DFV,DVAR,DADD,IFV,IVAR,IDOM,MM )
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
!     Load secondary variables w/ initial conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 June 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER IDOM(6),IFV(LSU,*)
      REAL*8 DVAR(*),DFV(LSV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINIR'
      IDOM(1) = MAX( 1,IDOM(1) )
      IDOM(1) = MIN( IDOM(1),IDOM(2),IFLD )
      IDOM(2) = MAX( 1,IDOM(1),IDOM(2) )
      IDOM(2) = MIN( IDOM(2),IFLD )
      IDOM(3) = MAX( 1,IDOM(3) )
      IDOM(3) = MIN( IDOM(3),IDOM(4),JFLD )
      IDOM(4) = MAX( 1,IDOM(3),IDOM(4) )
      IDOM(4) = MIN( IDOM(4),JFLD )
      IDOM(5) = MAX( 1,IDOM(5) )
      IDOM(5) = MIN( IDOM(5),IDOM(6),KFLD )
      IDOM(6) = MAX( 1,IDOM(5),IDOM(6) )
      IDOM(6) = MIN( IDOM(6),KFLD )
      WRITE (IWR,'(4X,A)') 'Domain: '
      WRITE (IWR,'(6X,2(A,I6))') 'I = ',IDOM(1),' to ',IDOM(2)
      WRITE (IWR,'(6X,2(A,I6))') 'J = ',IDOM(3),' to ',IDOM(4)
      WRITE (IWR,'(6X,2(A,I6))') 'K = ',IDOM(5),' to ',IDOM(6)
      IF( IDOM(1) .GT. IDOM(2) .OR. IDOM(3) .GT. IDOM(4) .OR.
     &  IDOM(5) .GT. IDOM(6) ) THEN
        NCH = INDEX( VARB(1:),'  ' )-1
        CHMSG = 'Invalid Initial Condition Domain: '//VARB(1:NCH)
        INDX = 4
        CALL WRMSGS( INDX )
      ENDIF
      DO K = IDOM(5), IDOM(6)
      DO J = IDOM(3), IDOM(4)
      DO I = IDOM(1), IDOM(2)
        N = ND(I,J,K)
        M = ND(IDOM(1),IDOM(3),IDOM(5))
        IF( IXP(N).EQ.0 ) CYCLE
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(MM,NRX) = DVAR(1) + DADD +
     &      (XP(NRX) - XP(M))*DVAR(2) +
     &      (YP(NRX) - YP(M))*DVAR(3) +
     &      (ZP(NRX) - ZP(M))*DVAR(4)
          IFV(MM,NRX) = IVAR
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINIR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINIS( DFV,DVAR,DADD,IDOM,MM )
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
!     Load secondary variables w/ initial conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 December 1992.
!     Last Modified by MD White, PNNL, 31 December 1992.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER IDOM(6)
      REAL*8 DVAR(*),DFV(LSV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINIS'
      IDOM(1) = MAX( 1,IDOM(1) )
      IDOM(1) = MIN( IDOM(1),IDOM(2),IFLD )
      IDOM(2) = MAX( 1,IDOM(1),IDOM(2) )
      IDOM(2) = MIN( IDOM(2),IFLD )
      IDOM(3) = MAX( 1,IDOM(3) )
      IDOM(3) = MIN( IDOM(3),IDOM(4),JFLD )
      IDOM(4) = MAX( 1,IDOM(3),IDOM(4) )
      IDOM(4) = MIN( IDOM(4),JFLD )
      IDOM(5) = MAX( 1,IDOM(5) )
      IDOM(5) = MIN( IDOM(5),IDOM(6),KFLD )
      IDOM(6) = MAX( 1,IDOM(5),IDOM(6) )
      IDOM(6) = MIN( IDOM(6),KFLD )
      WRITE (IWR,'(4X,A)') 'Domain: '
      WRITE (IWR,'(6X,2(A,I6))') 'I = ',IDOM(1),' to ',IDOM(2)
      WRITE (IWR,'(6X,2(A,I6))') 'J = ',IDOM(3),' to ',IDOM(4)
      WRITE (IWR,'(6X,2(A,I6))') 'K = ',IDOM(5),' to ',IDOM(6)
      IF( IDOM(1) .GT. IDOM(2) .OR. IDOM(3) .GT. IDOM(4) .OR.
     &  IDOM(5) .GT. IDOM(6) ) THEN
        NCH = INDEX( VARB(1:),'  ' )-1
        CHMSG = 'Invalid Initial Condition Domain: '//VARB(1:NCH)
        INDX = 4
        CALL WRMSGS( INDX )
      ENDIF
      DO K = IDOM(5), IDOM(6)
      DO J = IDOM(3), IDOM(4)
      DO I = IDOM(1), IDOM(2)
        N = ND(I,J,K)
        M = ND(IDOM(1),IDOM(3),IDOM(5))
        IF( IXP(N).EQ.0 ) CYCLE
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(MM,NRX) = DVAR(1) + DADD +
     &      (XP(NRX) - XP(M))*DVAR(2) +
     &      (YP(NRX) - YP(M))*DVAR(3) +
     &      (ZP(NRX) - ZP(M))*DVAR(4)
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINIS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINFGC( DFV,DVAR,DADD,UNTS,MM,IGC )
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
!     Load primary variables w/ initial conditions from an external
!     file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 March 2008.
!     Last Modified by MD White, PNNL, 31 March 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      REAL*8 DVAR(*),DFV(LNGC,LSV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINFC'
      ILUM = IUNM
      ILUS = IUNS
      ILUKG = IUNKG
      ILUK = IUNK
      ILUMOL = IUNMOL
      DO 10 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 10
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(IGC,MM,NRX) = DVAR(1)+DADD
        ENDDO
   10 CONTINUE
   20 CONTINUE
      READ(26,*,END=30) I,J,K,VAR
      IF( I.LT.1 .OR. I.GT.IFLD .OR. J.LT.1 .OR. J.GT.JFLD
     &  .OR. K.LT.1 .OR. K.GT.KFLD ) THEN
        CHMSG = 'Domain Index Out of Range: '
        INDX = 4
        CALL WRMSGS( INDX )
      ENDIF
      N = ND(I,J,K)
      IF( IXP(N).EQ.0 ) GOTO 20
      INDX = 0
      IUNM = ILUM
      IUNS = ILUS
      IUNKG = ILUKG
      IUNK = ILUK
      IUNMOL = ILUMOL
      CALL RDUNIT( UNTS,VAR,INDX )
      DO NRX = IBR(4,N),IBR(5,N)
        DFV(IGC,MM,NRX) = VAR + DADD
      ENDDO
      GOTO 20
   30 CONTINUE
      IUNM = 0
      IUNS = 0
      IUNKG = 0
      IUNK = 0
      IUNMOL = 0
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINFC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINFP( DFV,DVAR,DADD,IFV,IVAR,UNTS )
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
!     Load primary variables w/ initial conditions from an external
!     file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 December 1992.
!     Last Modified by MD White, PNNL, 31 December 1992.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      INTEGER*4 IVAR,IFV(*)
      REAL*8 DVAR(*),DFV(*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINFP'
      ILUM = IUNM
      ILUS = IUNS
      ILUKG = IUNKG
      ILUK = IUNK
      ILUMOL = IUNMOL
!
!---  pH  ---
!
      IF( ABS(DADD-7.D+0)/EPSL.LT.EPSL ) THEN
        DO 10 N = 1,NFLD
          IF( IXP(N).EQ.0 ) GOTO 10
          DO NRX = IBR(4,N),IBR(5,N)
            DVAR(1) = MIN( MAX( DVAR(1),0.D+0 ), 1.4D+1 )
            DFV(NRX) = 1.D+1**(-DVAR(1))
            IFV(NRX) = IVAR
          ENDDO
   10   CONTINUE
   20   CONTINUE
        READ(26,*,END=30) I,J,K,VAR
        IF( I.LT.1 .OR. I.GT.IFLD .OR. J.LT.1 .OR. J.GT.JFLD
     &    .OR. K.LT.1 .OR. K.GT.KFLD ) THEN
          CHMSG = 'Domain Index Out of Range: '
          INDX = 4
          CALL WRMSGS( INDX )
        ENDIF
        N = ND(I,J,K)
        IF( IXP(N).EQ.0 ) GOTO 20
        INDX = 0
        IUNM = ILUM
        IUNS = ILUS
        IUNKG = ILUKG
        IUNK = ILUK
        IUNMOL = ILUMOL
        CALL RDUNIT( UNTS,VAR,INDX )
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(NRX) = 1.D+3*(1.D+1**(-VAR))
          IFV(NRX) = IVAR
        ENDDO
        GOTO 20
   30   CONTINUE
      ELSE
        DO N = 1,NFLD
          IF( IXP(N).EQ.0 ) CYCLE
          DO NRX = IBR(4,N),IBR(5,N)
            DFV(NRX) = DVAR(1)+DADD
            IFV(NRX) = IVAR
          ENDDO
        ENDDO
  120   CONTINUE
        READ(26,*,END=130) I,J,K,VAR
        IF( I.LT.1 .OR. I.GT.IFLD .OR. J.LT.1 .OR. J.GT.JFLD
     &    .OR. K.LT.1 .OR. K.GT.KFLD ) THEN
          CHMSG = 'Domain Index Out of Range: '
          INDX = 4
          CALL WRMSGS( INDX )
        ENDIF
        N = ND(I,J,K)
        IF( IXP(N).EQ.0 ) GOTO 120
        INDX = 0
        IUNM = ILUM
        IUNS = ILUS
        IUNKG = ILUKG
        IUNK = ILUK
        IUNMOL = ILUMOL
        CALL RDUNIT( UNTS,VAR,INDX )
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(NRX) = VAR+DADD
          IFV(NRX) = IVAR
        ENDDO
        GOTO 120
  130   CONTINUE
      ENDIF
      IUNM = 0
      IUNS = 0
      IUNKG = 0
      IUNK = 0
      IUNMOL = 0
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINFP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINFR( DFV,DVAR,DADD,IFV,IVAR,UNTS,MM )
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
!     Load primary variables w/ initial conditions from an external
!     file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 June 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      REAL*8 DVAR(*),DFV(LSV,*)
      INTEGER IFV(LSU,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINFR'
      ILUM = IUNM
      ILUS = IUNS
      ILUKG = IUNKG
      ILUK = IUNK
      ILUMOL = IUNMOL
      DO 10 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 10
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(MM,NRX) = DVAR(1)+DADD
          IFV(MM,NRX) = IVAR
        ENDDO
   10 CONTINUE
   20 CONTINUE
      READ(26,*,END=30) I,J,K,VAR
      IF( I.LT.1 .OR. I.GT.IFLD .OR. J.LT.1 .OR. J.GT.JFLD
     &  .OR. K.LT.1 .OR. K.GT.KFLD ) THEN
        CHMSG = 'Domain Index Out of Range: '
        INDX = 4
        CALL WRMSGS( INDX )
      ENDIF
      N = ND(I,J,K)
      IF( IXP(N).EQ.0 ) GOTO 20
      INDX = 0
      IUNM = ILUM
      IUNS = ILUS
      IUNKG = ILUKG
      IUNK = ILUK
      IUNMOL = ILUMOL
      CALL RDUNIT( UNTS,VAR,INDX )
      DO NRX = IBR(4,N),IBR(5,N)
        DFV(MM,NRX) = VAR + DADD
        IFV(MM,NRX) = IVAR
      ENDDO
      GOTO 20
   30 CONTINUE
      IUNM = 0
      IUNS = 0
      IUNKG = 0
      IUNK = 0
      IUNMOL = 0
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINFR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINFS( DFV,DVAR,DADD,UNTS,MM )
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
!     Load primary variables w/ initial conditions from an external
!     file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 December 1992.
!     Last Modified by MD White, PNNL, 16 June 1997.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      REAL*8 DVAR(*),DFV(LSV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINFS'
      ILUM = IUNM
      ILUS = IUNS
      ILUKG = IUNKG
      ILUK = IUNK
      ILUMOL = IUNMOL
      DO 10 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 10
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(MM,NRX) = DVAR(1)+DADD
        ENDDO
   10 CONTINUE
   20 CONTINUE
      READ(26,*,END=30) I,J,K,VAR
      IF( I.LT.1 .OR. I.GT.IFLD .OR. J.LT.1 .OR. J.GT.JFLD
     &  .OR. K.LT.1 .OR. K.GT.KFLD ) THEN
        CHMSG = 'Domain Index Out of Range: '
        INDX = 4
        CALL WRMSGS( INDX )
      ENDIF
      N = ND(I,J,K)
      IF( IXP(N).EQ.0 ) GOTO 20
      INDX = 0
      IUNM = ILUM
      IUNS = ILUS
      IUNKG = ILUKG
      IUNK = ILUK
      IUNMOL = ILUMOL
      CALL RDUNIT( UNTS,VAR,INDX )
      DO NRX = IBR(4,N),IBR(5,N)
        DFV(MM,NRX) = VAR + DADD
      ENDDO
      GOTO 20
   30 CONTINUE
      IUNM = 0
      IUNS = 0
      IUNKG = 0
      IUNK = 0
      IUNMOL = 0
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINFS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINZGC( DFV,DVAR,DADD,IZN,MM,IGC )
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
!     Load secondary variables w/ initial conditions according to
!     rock/soil zonations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 March 2008.
!     Last Modified by MD White, PNNL, 31 March 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DFV(LNGC,LSV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINZS'
      DO N = 1,NFLD
        IF( IXP(N).EQ.0 ) CYCLE
        IF( IZ(N).EQ.IZN ) THEN
          DO NRX = IBR(4,N),IBR(5,N)
            DFV(IGC,MM,NRX) = DVAR + DADD
          ENDDO
        ENDIF
      ENDDO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINZS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINZP( DFV,DVAR,DADD,IFV,IVAR,IZN )
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
!     Load primary variables w/ initial conditions according to
!     rock/soil zonations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 August 1998.
!     Last Modified by MD White, PNNL, 25 August 1998.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER IVAR,IFV(*)
      REAL*8 DVAR,DFV(*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINZP'
!
!---  pH  ---
!
      IF( ABS(DADD-7.D+0)/EPSL.LT.EPSL ) THEN
        DO N = 1,NFLD
          IF( IXP(N).EQ.0 ) CYCLE
          IF( IZ(N).EQ.IZN ) THEN
            DO NRX = IBR(4,N),IBR(5,N)
              DVAR = MIN( MAX( DVAR,0.D+0 ), 1.4D+1 )
              DFV(NRX) = 1.D+1**(-DVAR)
              IFV(NRX) = IVAR
            ENDDO
          ENDIF
        ENDDO
      ELSE
        DO N = 1,NFLD
          IF( IXP(N).EQ.0 ) CYCLE
          IF( IZ(N).EQ.IZN ) THEN
            DO NRX = IBR(4,N),IBR(5,N)
              DFV(NRX) = DVAR + DADD
              IFV(NRX) = IVAR
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINZP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINZR( DFV,DVAR,DADD,IFV,IVAR,IZN,MM )
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
!     Load secondary variables w/ initial conditions according to
!     rock/soil zonations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 June 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DFV(LSV,*)
      INTEGER IFV(LSU,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINZR'
      DO N = 1,NFLD
        IF( IXP(N).EQ.0 ) CYCLE
        IF( IZ(N).EQ.IZN ) THEN
          DO NRX = IBR(4,N),IBR(5,N)
            DFV(MM,NRX) = DVAR + DADD
            IFV(MM,NRX) = IVAR
          ENDDO
        ENDIF
      ENDDO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINZR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINZS( DFV,DVAR,DADD,IZN,MM )
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
!     Load secondary variables w/ initial conditions according to
!     rock/soil zonations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 August 1998.
!     Last Modified by MD White, PNNL, 25 August 1998.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DFV(LSV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINZS'
      DO N = 1,NFLD
        IF( IXP(N).EQ.0 ) CYCLE
        IF( IZ(N).EQ.IZN ) THEN
          DO NRX = IBR(4,N),IBR(5,N)
            DFV(MM,NRX) = DVAR + DADD
          ENDDO
        ENDIF
      ENDDO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINZS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINBGC( DFV,DADD,UNTS,M,IGC )
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
!     Load primary variables w/ initial conditions from an external
!     binary file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 March 2008.
!     Last Modified by MD White, PNNL, 31 March 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      REAL*8 DFV(LNGC,LSV,*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: SFV
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINBC'
      ALLOCATE( SFV(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      INDX = 0
      VAR = 1.D+0
      CALL RDUNIT( UNTS,VAR,INDX )
      READ(26) (SFV(N),N=1,NFLD)
      DO 10 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 10
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(IGC,M,NRX) = SFV(N)*VAR + DADD
        ENDDO
   10 CONTINUE
      IF( ALLOCATED(SFV) ) THEN
      DEALLOCATE( SFV,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINBC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINBP( DFV,DADD,IFV,IVAR,UNTS )
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
!     Load primary variables w/ initial conditions from an external
!     binary file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 July 2000.
!     Last Modified by MD White, PNNL, 27 July 2000.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      INTEGER*4 IVAR,IFV(*)
      REAL*8 DFV(*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: SFV
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINBP'
      ALLOCATE( SFV(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      INDX = 0
      VAR = 1.D+0
      CALL RDUNIT( UNTS,VAR,INDX )
      READ(26) (SFV(N),N=1,NFLD)
!
!---  pH  ---
!
      IF( ABS(DADD-7.D+0)/EPSL.LT.EPSL ) THEN
        DO 10 N = 1,NFLD
          IF( IXP(N).EQ.0 ) GOTO 10
          DO NRX = IBR(4,N),IBR(5,N)
            IFV(NRX) = IVAR
            SFV(NRX) = MIN( MAX( SFV(N),0.D+0 ), 14.D+0 )
            DFV(NRX) = 1.D+1**(-SFV(N))
          ENDDO
   10  CONTINUE
      ELSE
        DO 20 N = 1,NFLD
          IF( IXP(N).EQ.0 ) GOTO 20
          DO NRX = IBR(4,N),IBR(5,N)
            IFV(NRX) = IVAR
            DFV(NRX) = SFV(N)*VAR + DADD
          ENDDO
   20   CONTINUE
      ENDIF      
      IF( ALLOCATED(SFV) ) THEN
      DEALLOCATE( SFV,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINBP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINBR( DFV,DADD,IFV,IVAR,UNTS,M )
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
!     Load primary variables w/ initial conditions from an external
!     binary file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, 11 June 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      REAL*8 DFV(LSV,*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: SFV
      INTEGER IFV(LSU,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINBR'
      ALLOCATE( SFV(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      INDX = 0
      VAR = 1.D+0
      CALL RDUNIT( UNTS,VAR,INDX )
      READ(26) (SFV(N),N=1,NFLD)
      DO 10 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 10
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(M,NRX) = SFV(N)*VAR + DADD
          IFV(M,NRX) = IVAR
        ENDDO
   10 CONTINUE
      IF( ALLOCATED(SFV) ) THEN
      DEALLOCATE( SFV,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINBR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINBS( DFV,DADD,UNTS,M )
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
!     Load primary variables w/ initial conditions from an external
!     binary file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 July 2000.
!     Last Modified by MD White, PNNL, 27 July 2000.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      REAL*8 DFV(LSV,*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: SFV
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINBS'
      ALLOCATE( SFV(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      INDX = 0
      VAR = 1.D+0
      CALL RDUNIT( UNTS,VAR,INDX )
      READ(26) (SFV(N),N=1,NFLD)
      DO 10 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 10
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(M,NRX) = SFV(N)*VAR + DADD
        ENDDO
   10 CONTINUE
      IF( ALLOCATED(SFV) ) THEN
      DEALLOCATE( SFV,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINBS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINAGC( DFV,DADD,UNTS,M,IGC )
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
!     Load primary variables w/ initial conditions from an external
!     binary file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 March 2008.
!     Last Modified by MD White, PNNL, 31 March 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      REAL*8 DFV(LNGC,LSV,*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: SFV
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINAC'
      ALLOCATE( SFV(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      INDX = 0
      VAR = 1.D+0
      CALL RDUNIT( UNTS,VAR,INDX )
      READ(26,*) (SFV(N),N=1,NFLD)
      DO 10 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 10
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(IGC,M,NRX) = SFV(N)*VAR + DADD
        ENDDO
   10 CONTINUE
      IF( ALLOCATED(SFV) ) THEN
      DEALLOCATE( SFV,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINAC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINAP( DFV,DADD,IFV,IVAR,UNTS )
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
!     Load primary variables w/ initial conditions from an external
!     binary file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 July 2000.
!     Last Modified by MD White, PNNL, 27 July 2000.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      INTEGER*4 IVAR,IFV(*)
      REAL*8 DFV(*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: SFV
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINAP'
      ALLOCATE( SFV(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      INDX = 0
      VAR = 1.D+0
      CALL RDUNIT( UNTS,VAR,INDX )
      READ(26,*) (SFV(N),N=1,NFLD)
!
!---  pH  ---
!
      IF( ABS(DADD-7.D+0)/EPSL.LT.EPSL ) THEN
        DO 10 N = 1,NFLD
          IF( IXP(N).EQ.0 ) GOTO 10
          DO NRX = IBR(4,N),IBR(5,N)
            IFV(NRX) = IVAR
            SFV(NRX) = MIN( MAX( SFV(N),0.D+0 ), 14.D+0 )
            DFV(NRX) = 1.D+1**(-SFV(NRX))
          ENDDO
   10   CONTINUE
      ELSE
        DO 20 N = 1,NFLD
          IF( IXP(N).EQ.0 ) GOTO 20
          DO NRX = IBR(4,N),IBR(5,N)
            IFV(NRX) = IVAR
            DFV(NRX) = SFV(N)*VAR + DADD
          ENDDO
   20   CONTINUE
      ENDIF
      IF( ALLOCATED(SFV) ) THEN
      DEALLOCATE( SFV,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINAP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINAR( DFV,DADD,IFV,IVAR,UNTS,M )
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
!     Load primary variables w/ initial conditions from an external
!     binary file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 June 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      REAL*8 DFV(LSV,*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: SFV
      INTEGER IFV(LSU,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINAR'
      ALLOCATE( SFV(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      INDX = 0
      VAR = 1.D+0
      CALL RDUNIT( UNTS,VAR,INDX )
      READ(26,*) (SFV(N),N=1,NFLD)
      DO 10 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 10
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(M,NRX) = SFV(N)*VAR + DADD
          IFV(M,NRX) = IVAR
        ENDDO
   10 CONTINUE
      IF( ALLOCATED(SFV) ) THEN
      DEALLOCATE( SFV,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINAR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINAS( DFV,DADD,UNTS,M )
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
!     Load primary variables w/ initial conditions from an external
!     binary file.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 July 2000.
!     Last Modified by MD White, PNNL, 27 July 2000.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS
      CHARACTER*512 CHDUM
      REAL*8 DFV(LSV,*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: SFV
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINAS'
      ALLOCATE( SFV(1:LFD),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      INDX = 0
      VAR = 1.D+0
      CALL RDUNIT( UNTS,VAR,INDX )
      READ(26,'(A)') CHDUM
      IF( INDEX( CHDUM(1:),'@' ).NE.0 ) THEN
        REWIND(26)
        NC = 0
        DO
          ISTART = 1
          READ(26,'(A)',END=90) CHDUM
          IATX = INDEX( CHDUM(1:), '@' )
          ICMX = INDEX( CHDUM(1:), '  ' )
          CHDUM(IATX:IATX) = ','
          CHDUM(ICMX:ICMX) = ','
          CALL RDINT(ISTART,ICOMMA,CHDUM,IATX )
          MC = NC + 1
          NC = MIN( NC+IATX,NFLD )
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SFVX )
          DO N = MC,NC
            SFV(N) = SFVX
          ENDDO
          IF( NC.GE.NFLD ) EXIT
        ENDDO
   90   IF( NC.LT.NFLD ) THEN
          INDX = 4
          CHMSG = 'End of Initial Conditions File Encountered '
          CALL WRMSGS( INDX )
        ENDIF
      ELSE
        REWIND(26)
        READ(26,*) (SFV(N),N=1,NFLD)
      ENDIF
      DO 10 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 10
        DO NRX = IBR(4,N),IBR(5,N)
          DFV(M,NRX) = SFV(N)*VAR + DADD
        ENDDO
   10 CONTINUE
      IF( ALLOCATED(SFV) ) THEN
      DEALLOCATE( SFV,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: SFV'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINAS group  ---
!
      RETURN
      END



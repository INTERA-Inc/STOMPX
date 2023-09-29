!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CLOCKS
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
!     Call system dependent subroutines to determine clock time, date,
!     and execution time.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, December 1992.
!     Last Modified by MD White, Battelle, PNL, December 31, 1992.
!     clocks.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
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
      CHARACTER*2 CH2X
      CHARACTER*4 CH4X
      CHARACTER*8 CH8X
      INTEGER IARR(3)
      REAL TARR(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CLOCKS'
#ifdef macintosh
      CALL DATE( IMON,IDAY,IYR )
      IF( IYR .LT. 92 ) THEN
        IYR = IYR + 2000
      ELSE
        IYR = IYR + 1900
      ENDIF
      CALL TIME( NSEC )
      IHR = NSEC/3600
      IMIN = (NSEC - IHR*3600)/60
      ISEC = NSEC - IHR*3600 - IMIN*60
      CLKSEC = NSEC - MOD(CLKSEC,86400.D+0)
      CPUSEC = CLKSEC
      TARR(1) = 0.
      TARR(2) = 0.
      WRITE( CH8X(1:2),'(A2)') IHR
      CH8X(3:3) = ':'
      WRITE( CH8X(4:5),'(A2)') IMIN
      CH8X(6:6) = ':'
      WRITE( CH8X(7:8),'(A2)') ISEC
#endif
#ifdef sun
      CALL IDATE(IARR)
      IMON = IARR(1)
      IDAY = IARR(2)
      IYR = IARR(3)
      CALL ITIME(IARR)
      IHR = IARR(1)
      IMIN = IARR(2)
      ISEC = IARR(3)
      CLKSEC = REAL( TIME() )
      CPUSEC = ETIME( TARR )
      WRITE( CH8X(1:2),'(A2)') IHR
      CH8X(3:3) = ':'
      WRITE( CH8X(4:5),'(A2)') IMIN
      CH8X(6:6) = ':'
      WRITE( CH8X(7:8),'(A2)') ISEC
#endif
#ifdef digital
      CALL IDATE(IARR)
      IMON = IARR(1)
      IDAY = IARR(2)
      IYR = IARR(3)
      CALL ITIME(IARR)
      IHR = IARR(1)
      IMIN = IARR(2)
      ISEC = IARR(3)
      CLKSEC = REAL( TIME() )
      CPUSEC = ETIME( TARR )
      WRITE( CH8X(1:2),'(A2)') IHR
      CH8X(3:3) = ':'
      WRITE( CH8X(4:5),'(A2)') IMIN
      CH8X(6:6) = ':'
      WRITE( CH8X(7:8),'(A2)') ISEC
#endif
#ifdef silcon_graphics
      CALL IDATE(IMON,IDAY,IYR)
      CALL ITIME(IARR)
      IHR = IARR(1)
      IMIN = IARR(2)
      ISEC = IARR(3)
      CLKSEC = SECNDS( ZERO )
      CPUSEC = ETIME( TARR )
      WRITE( CH8X(1:2),'(A2)') IHR
      CH8X(3:3) = ':'
      WRITE( CH8X(4:5),'(A2)') IMIN
      CH8X(6:6) = ':'
      WRITE( CH8X(7:8),'(A2)') ISEC
#endif
#ifdef ibm
      IMON = 1
      IDAY = 1
      IYR = 1992
      IHR = 0
      IMIN = 0
      ISEC = 0
      IARR(1) = 0
      IARR(2) = 0
      IARR(3) = 0
      CLKSEC = 0.D+0
      CPUSEC = 0.D+0
      TARR(1) = 0.
      TARR(2) = 0.
      WRITE( CH8X(1:2),'(A2)') IHR
      CH8X(3:3) = ':'
      WRITE( CH8X(4:5),'(A2)') IMIN
      CH8X(6:6) = ':'
      WRITE( CH8X(7:8),'(A2)') ISEC
#endif
#ifdef stardent
      CALL IDATE(IMON,IDAY,IYR)
      CALL ITIME(IARR)
      IHR = IARR(1)
      IMIN = IARR(2)
      ISEC = IARR(3)
      CLKSEC = SECNDS( ZERO )
      CPUSEC = ETIME( TARR )
      WRITE( CH8X(1:2),'(A2)') IHR
      CH8X(3:3) = ':'
      WRITE( CH8X(4:5),'(A2)') IMIN
      CH8X(6:6) = ':'
      WRITE( CH8X(7:8),'(A2)') ISEC
#endif
#ifdef hp
      CALL IDATE(IARR)
      IMON = IARR(1)
      IDAY = IARR(2)
      IYR = IARR(3)
      CALL ITIME(IARR)
      IHR = IARR(1)
      IMIN = IARR(2)
      ISEC = IARR(3)
      CLKSEC = REAL( TIME() )
      CPUSEC = ETIME( TARR )
      WRITE( CH8X(1:2),'(A2)') IHR
      CH8X(3:3) = ':'
      WRITE( CH8X(4:5),'(A2)') IMIN
      CH8X(6:6) = ':'
      WRITE( CH8X(7:8),'(A2)') ISEC
      TARR(1) = 0.
      TARR(2) = 0.
#endif
#ifdef convex
      CALL IDATE(IARR)
      IMON = IARR(1)
      IDAY = IARR(2)
      IYR = IARR(3)
      CALL ITIME(IARR)
      IHR = IARR(1)
      IMIN = IARR(2)
      ISEC = IARR(3)
      CLKSEC = SECNDS( ZERO )
      CPUSEC = ETIME( TARR )
      TARR(1) = 0.
      TARR(2) = 0.
      WRITE( CH8X(1:2),'(A2)') IHR
      CH8X(3:3) = ':'
      WRITE( CH8X(4:5),'(A2)') IMIN
      CH8X(6:6) = ':'
      WRITE( CH8X(7:8),'(A2)') ISEC
#endif
#ifdef cray
      CALL DATE( VAR )
      WRITE( CH8X,'(A8)' ) VAR
      READ( CH8X(1:2),'(I2)' ) IMON
      READ( CH8X(4:5),'(I2)' ) IDAY
      READ( CH8X(7:8),'(I2)' ) IYR
      IF( IYR .LT. 92 ) THEN
        IYR = IYR + 2000
      ELSE
        IYR = IYR + 1900
      ENDIF
      CALL CLOCK( VAR )
      WRITE( CH8X,'(A8)' ) VAR
      READ( CH8X(1:2),'(I2)' ) IHR
      READ( CH8X(4:5),'(I2)' ) IMIN
      READ( CH8X(7:8),'(I2)' ) ISEC
      IARR(1) = 0
      IARR(2) = 0
      IARR(3) = 0
      CALL SECOND( CPUSEC )
      CLKSEC = RTC()
      TARR(1) = 0.
      TARR(2) = 0.
#endif
!
!---  Load date and time character strings  ---
!
      CHDATE(3:3) = '/'
      CHDATE(6:6) = '/'
      CHTIME(3:3) = ':'
      CHTIME(6:6) = ':'
      WRITE(CH2X,'(I2)') IMON
      READ(CH2X,'(A2)') CHDATE(1:2)
      WRITE(CH2X,'(I2)') IDAY
      READ(CH2X,'(A2)') CHDATE(4:5)
      WRITE(CH4X,'(I4)') IYR
      READ(CH4X,'(A4)') CHDATE(7:10)
      WRITE(CH2X,'(I2)') IHR
      READ(CH2X,'(A2)') CHTIME(1:2)
      WRITE(CH2X,'(I2)') IMIN
      READ(CH2X,'(A2)') CHTIME(4:5)
      WRITE(CH2X,'(I2)') ISEC
      READ(CH2X,'(A2)') CHTIME(7:8)
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CLOCKS group ---
!
      RETURN
      END

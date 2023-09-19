!----------------------Function----------------------------------------!
!
      FUNCTION FSPLNY( FX,ITS,ITE )
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
!     Given arrays TBLX and TBLY containing a tabulated function
!     (i.e., y(n) = f(x(n))), with x(1) < x(2) < ... < x(n), and given
!     the array TBLDDY, which is output from subroutine SPLINE, and
!     given a value of FX this function returns a cubic spline
!     interpolated value.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1994.
!     Last Modified by MD White, Battelle, PNL, September 9, 1994.




!     fsplny.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TABL
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
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FSPLNY'
!
!---  Ascending table order  ---
!
      IF( TBLX(ITE).GT.TBLX(ITS) ) THEN
        IF( FX.LE.TBLX(ITS) ) THEN
          FSPLNY = TBLY(ITS)
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ELSEIF( FX.GE.TBLX(ITE) ) THEN
          FSPLNY = TBLY(ITE)
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
!
!---  Find the right place in the table by means of bisection  ---
!
        KLO = ITS
        KHI = ITE
   10   CONTINUE
        IF( KHI-KLO.GT.1 ) THEN
          K = (KHI+KLO)/2
          IF( TBLX(K).GT.FX ) THEN
            KHI = K
          ELSE
            KLO = K
          ENDIF
          GOTO 10
        ENDIF
        H = TBLX(KHI)-TBLX(KLO)
        IF( ABS(H)/EPSL.LT.EPSL ) THEN
          INDX = 3
          CHMSGX = 'Invalid Table'
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FSPLNY'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
!
!---  Evaluate cubic spline  ---
!
        A = (TBLX(KHI)-FX)/H
        B = (FX-TBLX(KLO))/H
        FSPLNY = A*TBLY(KLO)+B*TBLY(KHI)+
     &    ((A**3-A)*TBLDDY(KLO)+(B**3-B)*TBLDDY(KHI))*(H**2)/6.D+0
!
!---  Descending table order  ---
!
      ELSEIF( TBLX(ITE).LT.TBLX(ITS) ) THEN
        IF( FX.GE.TBLX(ITS) ) THEN
          FSPLNY = TBLY(ITS)
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ELSEIF( FX.LE.TBLX(ITE) ) THEN
          FSPLNY = TBLY(ITE)
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
!
!---  Find the right place in the table by means of bisection  ---
!
        KLO = ITS
        KHI = ITE
   20   CONTINUE
        IF( KHI-KLO.GT.1 ) THEN
          K = (KHI+KLO)/2
          IF( TBLX(K).LT.FX ) THEN
            KHI = K
          ELSE
            KLO = K
          ENDIF
          GOTO 20
        ENDIF
        H = TBLX(KLO)-TBLX(KHI)
        IF( ABS(H)/EPSL.LT.EPSL ) THEN
          INDX = 3
          CHMSGX = 'Invalid Table'
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FSPLNY'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
!
!---  Evaluate cubic spline  ---
!
        A = (TBLX(KLO)-FX)/H
        B = (FX-TBLX(KHI))/H
        FSPLNY = A*TBLY(KHI)+B*TBLY(KLO)+
     &    ((A**3-A)*TBLDDY(KHI)+(B**3-B)*TBLDDY(KLO))*(H**2)/6.D+0
      ELSE
        INDX = 3
        CHMSGX = 'Invalid Table'
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'FSPLNY'
        RLMSGX = 0.D+0
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FSPLNY group
!
      RETURN
      END



!----------------------Function----------------------------------------!
!
      FUNCTION FSPLNX( FY,ITS,ITE )
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
!     Given arrays TBLY and TBLX containing a tabulated function
!     (i.e., x(n) = f(y(n))), with y(1) < y(2) < ... < y(n), and given
!     the array TBLDDX, which is output from subroutine SPLINE, and
!     given a value of FY this function returns a cubic spline
!     interpolated value.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1994.
!     Last Modified by MD White, Battelle, PNL, September 9, 1994.
#ifdef sac
!
!     SAC Modification by WE Nichols, PNNL, June 12, 2003.
#endif
!     fsplnx.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/FSPLNX'
!
!---  Ascending table order  ---
!
      IF( TBLY(ITE).GT.TBLY(ITS) ) THEN
        IF( FY.LE.TBLY(ITS) ) THEN
          FSPLNX = TBLX(ITS)
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ELSEIF( FY.GE.TBLY(ITE) ) THEN
          FSPLNX = TBLX(ITE)
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
          IF( TBLY(K).GT.FY ) THEN
            KHI = K
          ELSE
            KLO = K
          ENDIF
          GOTO 10
        ENDIF
        H = TBLY(KHI)-TBLY(KLO)
        IF( ABS(H)/EPSL.LT.EPSL ) THEN
          INDX = 3
          CHMSGX = 'Invalid Table'
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FSPLNX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
!
!---  Evaluate cubic spline  ---
!
        A = (TBLY(KHI)-FY)/H
        B = (FY-TBLY(KLO))/H
        FSPLNX = A*TBLX(KLO)+B*TBLX(KHI)+
     &    ((A**3-A)*TBLDDX(KLO)+(B**3-B)*TBLDDX(KHI))*(H**2)/6.D+0
!
!---  Descending table order  ---
!
      ELSEIF( TBLY(ITE).LT.TBLY(ITS) ) THEN
        IF( FY.GE.TBLY(ITS) ) THEN
          FSPLNX = TBLX(ITS)
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ELSEIF( FY.LE.TBLY(ITE) ) THEN
          FSPLNX = TBLX(ITE)
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
          IF( TBLY(K).LT.FY ) THEN
            KHI = K
          ELSE
            KLO = K
          ENDIF
          GOTO 20
        ENDIF
        H = TBLY(KLO)-TBLY(KHI)
        IF( ABS(H)/EPSL.LT.EPSL ) THEN
          INDX = 3
          CHMSGX = 'Invalid Table'
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FSPLNX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
!
!---  Evaluate cubic spline  ---
!
        A = (TBLY(KLO)-FY)/H
        B = (FY-TBLY(KHI))/H
        FSPLNX = A*TBLX(KHI)+B*TBLX(KLO)+
     &    ((A**3-A)*TBLDDX(KHI)+(B**3-B)*TBLDDX(KLO))*(H**2)/6.D+0
      ELSE
        INDX = 3
        CHMSGX = 'Invalid Table'
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'FSPLNX'
        RLMSGX = 0.D+0
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FSPLNX group
!
      RETURN
      END


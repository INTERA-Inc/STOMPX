!----------------------Function----------------------------------------!
!
      FUNCTION DIFMN( FDL,FDH,DXL,DXH,FLX,INDX )
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
!     Computes interfacial averages
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNL, March, 1993.
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
      EXTERNAL ERFX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DIFMN'
!
!---  Negative indices for manual configuration  ---
!
      IF( INDX.LT.0 ) THEN
        IDMNX = ABS(INDX)
      ELSE
        IDMNX = IDMN(INDX)
      ENDIF
!
!---  Harmonic mean: default mode  ---
!
      IF( IDMNX.EQ.1 ) THEN
        IF( ABS(FDL*DXL+FDH*DXH)/EPSL.LT.EPSL ) THEN
          DIFMN = 0.D+0
        ELSE
          DIFMN = (FDL*FDH*(DXL+DXH))/(FDL*DXL+FDH*DXH)
        ENDIF
!
!---  Geometric mean  ---
!
      ELSEIF( IDMNX.EQ.2 ) THEN
        IF( ABS(DXL+DXH)/EPSL.LT.EPSL ) THEN
          DIFMN = 0.D+0
        ELSE
          DIFMN = SQRT((FDH**((2.D+0*DXL)/(DXL+DXH)))*
     &      (FDL**((2.D+0*DXH)/(DXL+DXH))))
        ENDIF
!
!---  Arithmetic mean  ---
!
      ELSEIF( IDMNX.EQ.3 ) THEN
        IF( ABS(DXL+DXH)/EPSL.LT.EPSL ) THEN
          DIFMN = 0.D+0
        ELSE
          DIFMN = (FDH*DXL+FDL*DXH)/(DXL+DXH)
        ENDIF
!
!---  Upwind mean  ---
!
      ELSEIF( IDMNX.EQ.4 ) THEN
        IF( FLX.GT.EPSL ) THEN
          DIFMN = FDL
        ELSEIF( FLX.LT.-EPSL ) THEN
          DIFMN = FDH
        ELSEIF( ABS(DXL+DXH)/EPSL.LT.EPSL ) THEN
          DIFMN = 0.D+0
        ELSE
          DIFMN = (FDH*DXL+FDL*DXH)/(DXL+DXH)
        ENDIF
!
!---  Downstream mean  ---
!
      ELSEIF( IDMNX.EQ.5 ) THEN
        IF( FLX.GT.EPSL ) THEN
          DIFMN = 5.D-1*((1.D+0+WFMN(INDX))*FDL+(1.D+0-WFMN(INDX))*FDH)
        ELSEIF( FLX.LT.-EPSL ) THEN
          DIFMN = 5.D-1*((1.D+0+WFMN(INDX))*FDH+(1.D+0-WFMN(INDX))*FDL)
        ELSEIF( ABS(DXL+DXH)/EPSL.LT.EPSL ) THEN
          DIFMN = 0.D+0
        ELSE
          DIFMN = (FDH*DXL+FDL*DXH)/(DXL+DXH)
        ENDIF
!
!---  Moderated upwind mean  ---
!
      ELSEIF( IDMNX.EQ.6 ) THEN
        IF( WFMN(INDX)/EPSL.LT.EPSL ) THEN
          WFMNX=ABS(FLX)
        ELSE
          WFMNX=WFMN(INDX)
        ENDIF
        DIFMN = MAX( ERFX(FLX/WFMNX),0.D+0 )*FDL +
     &    MAX( -ERFX(FLX/WFMNX),0.D+0 )*FDH
!
!---  Nieber downstream mean  ---
!
      ELSEIF( IDMNX.EQ.7 ) THEN
        IF( FLX.GT.EPSL ) THEN
         IF( FDH.LT.EPSL ) THEN
          DIFMN = 5.D-1*((1.D+0+WFMN(INDX))*FDL+(1.D+0-WFMN(INDX))*FDH)
         ELSEIF( ABS(DXL+DXH)/EPSL.LT.EPSL ) THEN
          DIFMN = 0.D+0
         ELSE
          DIFMN = (FDH*DXL+FDL*DXH)/(DXL+DXH)
         ENDIF
        ELSEIF( FLX.LT.-EPSL ) THEN
         IF( FDL.LT.EPSL ) THEN
          DIFMN = 5.D-1*((1.D+0+WFMN(INDX))*FDH+(1.D+0-WFMN(INDX))*FDL)
         ELSEIF( ABS(DXL+DXH)/EPSL.LT.EPSL ) THEN
          DIFMN = 0.D+0
         ELSE
          DIFMN = (FDH*DXL+FDL*DXH)/(DXL+DXH)
         ENDIF
        ELSEIF( ABS(DXL+DXH)/EPSL.LT.EPSL ) THEN
          DIFMN = 0.D+0
        ELSE
          DIFMN = (FDH*DXL+FDL*DXH)/(DXL+DXH)
        ENDIF
!
!---  Downwind  ---
!
      ELSEIF( IDMNX.EQ.8 ) THEN
        IF( FLX.GT.EPSL ) THEN
          DIFMN = FDH
        ELSEIF( FLX.LT.-EPSL ) THEN
          DIFMN = FDL
        ELSEIF( ABS(DXL+DXH)/EPSL.LT.EPSL ) THEN
          DIFMN = 0.D+0
        ELSE
          DIFMN = (FDH*DXL+FDL*DXH)/(DXL+DXH)
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIFMN group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION FLIMIT( R,CRN,INDX )
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
!     Computes interfacial averages
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNL, December 14, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLIMIT'
!
!---  Leonard-TVD flux limiter  ---
!
      IF( INDX.EQ.1 ) THEN
        FLIMIT = MAX( 0.D+0,MIN( 2.D+0,2.D+0*R,
     &    ((2.D+0-CRN+R*(1.D+0+CRN))/3.D+0) ) )
!
!---  Roe's Superbee flux limiter  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
        FLIMIT = MAX( 0.D+0,MIN( 2.D+0*R,1.D+0 ),MIN( R,2.D+0 ) )
!
!---  First-order upwind  ---
!
      ELSEIF( INDX.EQ.3 ) THEN
        FLIMIT = 0.D+0
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLIMIT group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION RLIMIT( V1X,V2X,V3X,V4X )
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
!     Computes interfacial averages
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 06 August 2010.
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RLIMIT'
!
!---  Zero numerator  ---
!
      IF( ABS(V1X)/EPSL.LT.EPSL ) THEN
        RLIMIT = 0.D+0
!
!---  Zero denominator  ---
!
      ELSEIF( ABS(V2X)/EPSL.LT.EPSL ) THEN
        RLIMIT = 1.D+20
!
!---  Non-zero numerator and denominator  ---
!
      ELSE
        RLIMIT = (V1X/V2X)*(V3X/V4X)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RLIMIT group
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
!     Last Modified by M.D. White, Battelle, PNL, April 14, 1994.
!     difmn.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ICOUNT'
      IC = ABS(I)
      ICOUNT = 0
   10 CONTINUE
      ICOUNT = ICOUNT + 1
      IC = IC/10
      IF( IC.GT.0 ) GOTO 10
      IF( I.LT.0 ) ICOUNT = ICOUNT + 1
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ICOUNT group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION DSHIFT( SV,SVP,DPV,DPVP,IEQX,M )
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
!     Written by M.D. White, Battelle, PNL, March, 2000.
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DSHIFT'
      IF( M.EQ.(IEQX+2) ) THEN
        DSHIFT = (SVP*DPV-SV*DPV+SV*DPVP)/(DPVP+SMALL)
      ELSE
        DSHIFT = SVP
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DSHIFT group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION ERFX( XX )
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
!     Error function.
!
!     Numerical Recipes, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 September 2008.
!     Last Modified by M.D. White, PNNL, 22 September 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ERFX'
      IF( XX.LT.0.D+0 ) THEN
        ERFX = -GAMMP(5.D-1,(XX**2))
      ELSE
        ERFX = GAMMP(5.D-1,(XX**2))
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ERFX group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION GAMMLN( XX )
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
!     Gamma function.
!
!     Numerical Recipes, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 September 2008.
!     Last Modified by M.D. White, PNNL, 22 September 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 COFX(6)
!
!----------------------Data Statements---------------------------------!
!
      DATA COFX /76.18009173D+0,-86.50532033D+0,24.01409822D+0,
     &  -1.231739516D+0,0.120858003D-2,-0.536382D-5/
      DATA STPX /2.50662827465D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/GAMMLN'
      XMX = XX - 1.D+0
      GX = XMX + 5.5D+0
      GX = (XMX + 5.D-1)*LOG(GX) - GX
      SERX = 1.D+0
      DO 10 J = 1,6
        XMX = XMX + 1.D+0
        SERX = SERX + COFX(J)/XMX
   10 CONTINUE
      GAMMLN = GX + LOG(STPX*SERX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of GAMMLN group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION GAMMP( AX,XX )
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
!     Incomplete gamma function.
!
!     Numerical Recipes, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 September 2008.
!     Last Modified by M.D. White, PNNL, 22 September 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/GAMMP'
!
!---  Variable out of range  ---
!
      IF( XX.LT.0.D+0 .OR. AX.LE.0.D+0 ) THEN
        M_ERR(1) = 'GAMMP Function Argument Out of Range: XX = '
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = AX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        XX = 5.D-1
        AX = 5.D-1
      ENDIF
      IF( XX.LT.(AX+1.D+0) ) THEN
        CALL GSER( GAMSER,AX,XX,GLN )
        GAMMP = GAMSER
      ELSE
        CALL GCF( GAMMCF,AX,XX,GLN )
        GAMMP = 1.D+0-GAMMCF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of GAMMP group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      SUBROUTINE GSER( GAMSER,AX,XX,GLN )
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
!     Incomplete gamma function P(AX,XX) evaluated by its series 
!     representation as GAMSER.
!
!     Numerical Recipes, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 September 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      EXTERNAL GAMMLN
      PARAMETER (ITMAX=100,EPS=3.D-7)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/GSER'
      GLN = GAMMLN(AX)
      IF( XX.LE.0.D+0 ) THEN
        GAMSER = 0.D+0
      ELSE
        APX = AX
        SUMX = 1.D+0/AX
        DELX = SUMX
        DO 10 N = 1,ITMAX
          APX = APX + 1.D+0
          DELX = DELX*XX/APX
          SUMX = SUMX + DELX
          IF( ABS(DELX).LT.ABS(SUMX)*EPS ) GOTO 20
   10   CONTINUE
!
!---    No convergence  ---
!
        M_ERR(1) = 'GSER Function Failure: AX = '
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = AX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        SUMX = 1.D+0/AX
   20   CONTINUE
        GAMSER = SUMX*EXP(-XX+AX*LOG(XX)-GLN)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of GSER group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      SUBROUTINE GCF( GAMMCF,AX,XX,GLN )
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
!     Incomplete gamma function Q(AX,XX) evaluated by its continued
!     fraction representation as GAMMCF.
!
!     Numerical Recipes, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 September 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      EXTERNAL GAMMLN
      PARAMETER (ITMAX=100,EPS=3.D-7)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/GCF'
      GLN = GAMMLN(AX)
      GOLDX = 0.D+0
      A0X = 1.D+0
      A1X = XX
      B0X = 0.D+0
      B1X = 1.D+0
      FACX = 1.D+0
      DO 10 N = 1,ITMAX
        ANX = REAL(N)
        ANAX = ANX-AX
        A0X = (A1X+A0X*ANAX)*FACX
        B0X = (B1X+B0X*ANAX)*FACX
        ANFX = ANX*FACX
        A1X = XX*A0X+ANFX*A1X
        B1X = XX*B0X+ANFX*B1X
        IF( ABS(A1X)/EPSL.GT.EPSL ) THEN
          FACX = 1.D+0/A1X
          GX = B1X*FACX
          IF( ABS((GX-GOLDX)/GX).LT.EPS ) GOTO 20
          GOLDX = GX
        ENDIF
   10 CONTINUE
!
!---  No convergence  ---
!
      M_ERR(1) = 'GCF Function Failure: AX = '
      IF( N_DB.LT.0 ) THEN
        M_ERR(2) = ' at Boundary Surface: '
      ELSE
        M_ERR(2) = ' at Node: '
      ENDIF
      CALL PATH
      R_ERR = AX
      I_ERR(1) = ABS(N_DB)
      I_ERR(2) = 1
      I_ERR(3) = 2
      I_ERR(4) = ID
   20 CONTINUE
      GAMMCF = EXP(-XX+AX*LOG(XX)-GLN)*GX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of GCF group
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ELC_DEN( RHOLX,CLX,DCFX )
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
!     Compute brine density from brine solution mass fraction and
!     liquid water density.
!
!     Leijnse, A.  1992.  Three-Dimensional Modeling of Coupled Flow
!     and Transport in Porous Media.  Ph.D. Dissertation, Department
!     of Civil Engineering and Geological Sciences, University of
!     Notre Dame, Notre Dame, Indiana.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DCFX(4)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ELC_DEN'
!
!---  Electrolyte aqeuous density function by A. Leijnse  ---
!
      IF( IDF_ELC.EQ.1 ) THEN
        NC = 0
        RHOLSX = RHOLX
        DO
          NC = NC+1
          F = RHOLSX - RHOLX*EXP(DCFX(1)*CLX/RHOLSX)
          DF = 1.D+0 + RHOLX*DCFX(1)*CLX*EXP(DCFX(1)*CLX/RHOLSX)/
     &      (RHOLSX**2)
          DRHOX = -F/DF
          RHOLSX = RHOLSX + DRHOX
          IF( NC.GT.32 ) THEN
            M_ERR(1) = 'Convergence Failure on Electrolyte Density: '
            M_ERR(2) = ' Density = '
            IF( N_DB.LT.0 ) THEN
              M_ERR(3) = ' at Boundary Surface: '
            ELSE
              M_ERR(3) = ' at Node: '
            ENDIF
            CALL PATH
            R_ERR = RHOLSX
            I_ERR(1) = ABS(N_DB)
            I_ERR(2) = 2
            I_ERR(3) = 3
            I_ERR(4) = ID
          ENDIF
          IF( ABS(DRHOX/RHOLSX).LE.1.D-9 ) EXIT
        ENDDO
        RHOLX = RHOLSX
!
!---  Fourth-order electrolyte aqeuous density function  ---
!
      ELSEIF( IDF_ELC.EQ.2 ) THEN
        CLXX = CLX/ELC_DUN
        RHOLX = RHOLX*(DCFX(1) + DCFX(2)*CLXX + DCFX(3)*(CLXX**2) +
     &    DCFX(4)*(CLXX**3))
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ELC_DEN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ELC_VIS( VISLX,CLX,VCFX )
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
!     Compute brine viscosity from brine solution mass fraction and
!     liquid water viscosity.
!
!     Leijnse, A.  1992.  Three-Dimensional Modeling of Coupled Flow
!     and Transport in Porous Media.  Ph.D. Dissertation, Department
!     of Civil Engineering and Geological Sciences, University of
!     Notre Dame, Notre Dame, Indiana.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VCFX(4)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ELC_VIS'
!
!---  Electrolyte aqeuous viscosity function by A. Leijnse  ---
!
      IF( IVF_ELC.EQ.1 ) THEN
        VISLX = VISLX*(VCFX(1) + VCFX(2)*CLX + VCFX(3)*(CLX**2)
     &    + VCFX(4)*(CLX**3))
!
!---  Fourth-order electrolyte aqeuous density function  ---
!
      ELSEIF( IDF_ELC.EQ.2 ) THEN
        CLXX = CLX/ELC_VUN
        VISLX = VISLX*(VCFX(1) + VCFX(2)*CLXX + VCFX(3)*(CLXX**2)
     &    + VCFX(4)*(CLXX**3))
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ELC_VIS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FNHGBL( XSBCX,YSBCX,ZSBCX,XPBCX,YPBCX,ZPBCX,TSX,
     &  TX,PLX )
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
!     Compute the aqueous pressure at a boundary surface, assuming
!     hydrostatic conditions, from a given pressure at a boundary
!     surface. Use twenty density calculations being points. 
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 July 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FNHGBL'
!
!---  x,y,z deltas  ---
!
      XDX = 5.D-2*(XPBCX - XSBCX)
      YDX = 5.D-2*(YPBCX - YSBCX)
      ZDX = 5.D-2*(ZPBCX - ZSBCX)
      DO N = 1,20
        TDX = (TX-TSX)*5.D-2 + TSX
        CALL WATLQD( TDX,PLX,RHOLX )
        PLX = PLX - RHOLX*(XDX*GRAVX + YDX*GRAVY + ZDX*GRAVZ)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FNHGBL group  ---
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION FNTBLX( FY,N1,N2,INDX )
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
!     Linearly interpolate from a table of values.
!
!     INDX = 0 : Table truncation beyond table limits.
!     INDX = 1 : Table extrapolation beyond table limits.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNL, March, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TABL
      USE GRID
      USE SOLTN
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
      SUB_LOG(ISUB_LOG) = '/FNTBLX'
!
!---  Ascending table order  ---
!
      IF( TBLY(N2).GT.TBLY(N1) ) THEN
        IF( FY.LE.TBLY(N1) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLX = TBLX(N1)
          ELSE
            FNTBLX = (FY-TBLY(N1))*(TBLX(N1+1)-TBLX(N1))/
     &        (TBLY(N1+1)-TBLY(N1)) + TBLX(N1)
          ENDIF
        ELSEIF( FY.GE.TBLY(N2) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLX = TBLX(N2)
          ELSE
            FNTBLX = (FY-TBLY(N2))*(TBLX(N2)-TBLX(N2-1))/
     &        (TBLY(N2)-TBLY(N2-1)) + TBLX(N2)
          ENDIF
        ELSE
          DO N = N1+1,N2
            IF( FY.LE.TBLY(N) ) THEN
              FNTBLX = (FY-TBLY(N-1))*(TBLX(N)-TBLX(N-1))/
     &          (TBLY(N)-TBLY(N-1)) + TBLX(N-1)
              EXIT
            ENDIF
          ENDDO
        ENDIF
!
!---  Descending table order  ---
!
      ELSEIF( TBLY(N2) .LT. TBLY(N1) ) THEN
        IF( FY.LE.TBLY(N2) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLX = TBLX(N2)
          ELSE
            FNTBLX = (FY-TBLY(N2))*(TBLX(N2-1)-TBLX(N2))/
     &        (TBLY(N2-1)-TBLY(N2)) + TBLX(N2)
          ENDIF
        ELSEIF( FY.GE.TBLY(N1) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLX = TBLX(N1)
          ELSE
            FNTBLX = (FY-TBLY(N1))*(TBLX(N1)-TBLX(N1+1))/
     &        (TBLY(N1)-TBLY(N1+1)) + TBLX(N1)
          ENDIF
        ELSE
          DO N = N1+1,N2
            IF( FY.GE.TBLY(N) ) THEN
              FNTBLX = (FY-TBLY(N))*(TBLX(N-1)-TBLX(N))/
     &          (TBLY(N-1)-TBLY(N)) + TBLX(N)
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ELSE
        A = 5.D-1
        B = 5.D-1
        FNTBLX = A*TBLX(N1) + B*TBLX(N2)
        M_ERR(1) = 'Invalid Table in FNTBLX '
        CALL PATH
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
        ICNV = 4
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FNTBLX group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION FNTBLY( FX,N1,N2,INDX )
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
!
!----------------------Description-------------------------------------!
!
!     Linearly interpolate from a table of values.
!
!     INDX = 0 : Table truncation beyond table limits.
!     INDX = 1 : Table extrapolation beyond table limits.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNL, March, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TABL
      USE GRID
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
      SUB_LOG(ISUB_LOG) = '/FNTBLY'
!
!---  Ascending table order  ---
!
      IF( TBLX(N2).GT.TBLX(N1) ) THEN
        IF( FX .LE. TBLX(N1) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLY = TBLY(N1)
          ELSE
            FNTBLY = (FX-TBLX(N1))*(TBLY(N1+1)-TBLY(N1))/
     &        (TBLX(N1+1)-TBLX(N1)) + TBLY(N1)
          ENDIF
        ELSEIF( FX .GE. TBLX(N2) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLY = TBLY(N2)
          ELSE
            FNTBLY = (FX-TBLX(N2))*(TBLY(N2)-TBLY(N2-1))/
     &        (TBLX(N2)-TBLX(N2-1)) + TBLY(N2)
          ENDIF
        ELSE
          DO N = N1+1,N2
            IF( FX .LE. TBLX(N) ) THEN
              FNTBLY = (FX-TBLX(N-1))*(TBLY(N)-TBLY(N-1))/
     &          (TBLX(N)-TBLX(N-1)) + TBLY(N-1)
              EXIT
            ENDIF
          ENDDO
        ENDIF
!
!---  Descending table order  ---
!
      ELSEIF( TBLX(N2).LT.TBLX(N1) ) THEN
        IF( FX .LE. TBLX(N2) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLY = TBLY(N2)
          ELSE
            FNTBLY = (FX-TBLX(N2))*(TBLY(N2-1)-TBLY(N2))/
     &        (TBLX(N2-1)-TBLX(N2)) + TBLY(N2)
          ENDIF
        ELSEIF( FX .GE. TBLX(N1) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLY = TBLY(N1)
          ELSE
            FNTBLY = (FX-TBLX(N1))*(TBLY(N1)-TBLY(N1+1))/
     &        (TBLX(N1)-TBLX(N1+1)) + TBLY(N1)
          ENDIF
        ELSE
          DO  N = N1+1,N2
            IF( FX .GE. TBLX(N) ) THEN
              FNTBLY = (FX-TBLX(N))*(TBLY(N-1)-TBLY(N))/
     &          (TBLX(N-1)-TBLX(N)) + TBLY(N)
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ELSE
        A = 5.D-1
        B = 5.D-1
        FNTBLY = A*TBLY(N1) + B*TBLY(N2)
        M_ERR(1) = 'Invalid Table in FNTBLY '
        CALL PATH
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FNTBLY group
!
      RETURN
      END

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
!     Written by M.D. White, Battelle, PNL, September, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TABL
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
!
!---      Reset subroutine string sequence  ---
!
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ELSEIF( FX.GE.TBLX(ITE) ) THEN
          FSPLNY = TBLY(ITE)
!
!---      Reset subroutine string sequence  ---
!
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
          H = 1.D+0
          M_ERR(1) = 'Invalid Table in FSPLNY '
          CALL PATH
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
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
!
!---      Reset subroutine string sequence  ---
!
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ELSEIF( FX.LE.TBLX(ITE) ) THEN
          FSPLNY = TBLY(ITE)
!
!---      Reset subroutine string sequence  ---
!
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
          H = 1.D+0
          M_ERR(1) = 'Invalid Table in FSPLNY '
          CALL PATH
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
        ENDIF
!
!---  Evaluate cubic spline  ---
!
        A = (TBLX(KLO)-FX)/H
        B = (FX-TBLX(KHI))/H
        FSPLNY = A*TBLY(KHI)+B*TBLY(KLO)+
     &    ((A**3-A)*TBLDDY(KHI)+(B**3-B)*TBLDDY(KLO))*(H**2)/6.D+0
      ELSE
        M_ERR(1) = 'Invalid Table in FSPLNY '
        CALL PATH
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
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

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LUDCMP( A,N,NPX,IX,D )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Numerical Recipes, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, August 1, 2000.
!     Last Modified by Mark White, PNNL, August 1, 2000.
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 A(NPX,NPX),VV(NPX)
      INTEGER IX(NPX)
!
!----------------------Executable Lines--------------------------------!
!
      D = 1.D+0
      DO I = 1,N
        AAMAX = 0.D+0
        DO J = 1,N
          IF( ABS(A(I,J)).GT.AAMAX ) AAMAX = ABS(A(I,J))
        ENDDO
        IF( ABS(AAMAX)/EPSL.LT.EPSL ) THEN
          PRINT *,'LUDCMP: Singular Matrix: NPX = ',NPX,' ID = ',ID
          ICUTTS = 1
        ENDIF
        VV(I) = 1.D+0/AAMAX
      ENDDO
      IMAX = 0
      DO J = 1,N
        DO I = 1,J-1
          SUM = A(I,J)
          DO K = 1,I-1
            SUM = SUM - A(I,K)*A(K,J)
          ENDDO
          A(I,J) = SUM
        ENDDO
        AAMAX = 0.D+0
        DO I = J,N
          SUM = A(I,J)
          DO K = 1,J-1
            SUM = SUM - A(I,K)*A(K,J)
          ENDDO
          A(I,J) = SUM
          DUM = VV(I)*ABS(SUM)
          IF( DUM.GE.AAMAX ) THEN
            IMAX = I
            AAMAX = DUM
          ENDIF
        ENDDO
        IF( J.NE.IMAX ) THEN
          DO K = 1,N
            DUM = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = DUM
          ENDDO
          D = -D
          VV(IMAX) = VV(J)
        ENDIF
        IX(J) = IMAX
        IF( ABS(A(J,J))/EPSL.LT.EPSL ) A(J,J) = 1.D-30
        IF( J.NE.N ) THEN
          DUM = 1.D+0/A(J,J)
          DO I = J+1,N
            A(I,J) = A(I,J)*DUM
          ENDDO
        ENDIF
      ENDDO
!
!---  End of LUDCMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LUBKSB( A,N,NPX,IX,B )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Numerical Recipes, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, August 1, 2000.
!     Last Modified by Mark White, PNNL, August 1, 2000.
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 A(NPX,NPX),B(NPX)
      INTEGER IX(NPX)
!
!----------------------Executable Lines--------------------------------!
!
      II = 0
      DO I = 1,N
        IL = IX(I)
        SUM = B(IL)
        B(IL) = B(I)
        IF( II.NE.0 ) THEN
          DO J = II,I-1
            SUM = SUM - A(I,J)*B(J)
          ENDDO
        ELSEIF( ABS(SUM)/EPSL.GT.EPSL ) THEN
          II = I
        ENDIF
        B(I) = SUM
      ENDDO
      DO I = N,1,-1
        SUM = B(I)
        IF( I.LT.N ) THEN
          DO J = I+1,N
            SUM = SUM - A(I,J)*B(J)
          ENDDO
        ENDIF
        B(I) = SUM/A(I,I)
      ENDDO
!
!---  End of LUBKSB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PERM_I( PERMRFX,PORDX )
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
!     Kozeny-Carmen  Porosity-Permeability Relationship
!
!     Original:
!     Carman PC. 1937. Fluid flow through granular beds. Transactions of
!     the Institution of Chemical Engineers, 15:150-156.
!     Reprinted:
!     Carman PC. 1997. Fluid flow through granular beds. Chemical 
!     Engineering Research & Design, 75:S32-S48.
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 11 Nov 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
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
      SUB_LOG(ISUB_LOG) = '/PERM_I'
!
!---  Permeability reduction factor  ---
!
      PERMRFX = PORDX*PORDX*PORDX/((1.D+0 - PORDX)*(1.D+0 - PORDX))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PERM_I group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY_W( N,PX,PREFX,PORDX,PORTX )
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
!     Compute diffusive and total porosities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 June 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE GEO_MECH
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PORSTY_W'
      DPX = PX-PREFX
      ISLC50X = ABS(ISLC(50))
!
!---  Geomechanics simulations  --
!
      IF( ISLC(50).NE.0 .AND. MOD(ISLC50X,10).NE.2 .AND. 
     &  MOD(ISLC50X,10).NE.4 ) THEN
!
!---    Drained bulk modulus  --
!
        BLKDRNX = PROP_GM(1,N)/(3.0D+0*(1.0D+0-2.0D+0*PROP_GM(2,N)))
!
!---    1/N  ---
!
        OONMODX = (PROP_GM(3,N)-POR(1,N))*(1.0D+0-PROP_GM(3,N))
     &    /BLKDRNX
!
!---    Volumetric strain differential, at iterate level k  ---
!
        DEVKX = EPSV_GM(2,N)-EPSV_CMP(N)
!
!---    Pressure differential, at iterate level k  ---
!
        DPKX = P_GM(2,N)-PCMP(N)
!
!---    Diffusive and total mechanical porosity, at iterate level k  ---
!
        PORDMCHKX = POR(1,N) - PROP_GM(3,N)*DEVKX + OONMODX*DPKX
        PORTMCHKX = POR(2,N) - PROP_GM(3,N)*DEVKX + OONMODX*DPKX
!
!---    Pressure differential, at iterate level k+1  ---
!
        DPK1X = PX-P_GM(2,N)
!
!---    Diffusive and total flow porosity, at iterate level k+1  ---
!
        PORDX = PORDMCHKX + ((PROP_GM(3,N)**2)/BLKDRNX+OONMODX)*DPK1X
        PORTX = PORTMCHKX + ((PROP_GM(3,N)**2)/BLKDRNX+OONMODX)*DPK1X
!
!---  No geomechanics  --
!
      ELSE        
!
!---    Pore compressibility w/ fixed bulk volume  ---
!
        IF( ISLC(15).EQ.1 ) THEN
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0(1,N)*EXP(DPX*CMP(1,N))
            PORDX = POR0(2,N)*EXP(DPX*CMP(1,N))
          ELSE
            PORTX = POR(1,N)*EXP(DPX*CMP(1,N))
            PORDX = POR(2,N)*EXP(DPX*CMP(1,N))
          ENDIF
!
!---    Bulk compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.10 ) THEN
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0(1,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR0(1,N)/
     &        POR0(1,N)))
            PORDX = POR0(2,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR0(2,N)/
     &        POR0(2,N)))
          ELSE
            PORTX = POR(1,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR(1,N)/
     &        POR(1,N)))
            PORDX = POR(2,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR(2,N)/
     &        POR(2,N)))
          ENDIF
!
!---    Pore compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.11 ) THEN
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0(1,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR0(1,N)))
            PORDX = POR0(2,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR0(2,N)))
          ELSE
            PORTX = POR(1,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR(1,N)))
            PORDX = POR(2,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR(2,N)))
          ENDIF
!
!---    Bulk compressibility w/ fixed bulk volume  ---
!
        ELSE
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0(1,N) + DPX*CMP(1,N)
            PORDX = POR0(2,N) + DPX*CMP(1,N)
          ELSE
            PORTX = POR(1,N) + DPX*CMP(1,N)
            PORDX = POR(2,N) + DPX*CMP(1,N)
          ENDIF
        ENDIF
      ENDIF
      PORTX = MAX( MIN( PORTX,1.D+0 ),1.D-6 )
      PORDX = MAX( MIN( PORDX,1.D+0 ),1.D-6 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TORTU_W( N,SLX,PORDX,TORLX )
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
!     Compute phase tortuosity.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
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
      SUB_LOG(ISUB_LOG) = '/TORTU_W'
!
!---  Constant phase tortuosity  ---
!
      IF( ITOR(IZN).EQ.1 .OR. ITOR(IZN).EQ.5 ) THEN
        TORLX = TOR(1,IZN)
!
!---  Millington and Quirk tortuosity model  ---
!
      ELSEIF( ITOR(IZN).EQ.2 .OR. ITOR(IZN).EQ.3 ) THEN
        IF( SLX*PORDX.LT.EPSL ) THEN
          TORLX = 0.D+0
        ELSE
          TORLX = (PORDX*(SLX**7))**(THIRD)
        ENDIF
!
!---  Marshal tortuosity model  ---
!
      ELSEIF( ITOR(IZN).EQ.4 ) THEN
        IF( SLX*PORDX.LT.EPSL ) THEN
          TORLX = 0.D+0
        ELSE
          TORLX = SQRT(PORDX*SLX)
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TORTU_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATGSD( TX,PX,RHOX,INDX )
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
!     INDX = 0
!
!     Calculate the water vapor component density
!     with the ideal gas law equation of state.
!
!     INDX = 1
!
!     Calculate the water vapor component density, as a function of
!     temperature and pressure per the steam table equations
!     as given by the 1967 International Formulation Committee:
!     Formulation for Industrial Use.
!
!     Thermodynamic and Transport Properties of Steam.
!     1967. ASME Steam Tables.
!     The American Society of Mechanical Engineers.
!     United Engineering Center, 345 East 47th Street, New York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < T > 364.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June 2022.
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
      REAL*8 BX(31),SBX(5),L(3)
!
!----------------------Data Statements---------------------------------!
!
      DATA BX /1.683599274D+1,2.856067796D+1,-5.438923329D+1,
     &4.330662834D-1,
     &-6.547711697D-1,8.565182058D-2,6.670375918D-2,1.388983801D+0,
     &8.390104328D-2,2.614670893D-2,-3.373439453D-2,4.520918904D-1,
     &1.069036614D-1,-5.975336707D-1,-8.847535804D-2,5.958051609D-1,
     &-5.159303373D-1,2.075021122D-1,1.190610271D-1,-9.867174132D-2,
     &1.683998803D-1,-5.809438001D-2,6.552390126D-3,5.710218649D-4,
     &1.936587558D+2,-1.388522425D+3,4.126607219D+3,-6.508211677D+3,
     &5.745984054D+3,-2.693088365D+3,5.235718623D+2/
      DATA SBX /7.633333333D-1,4.006073948D-1,8.636081627D-2,
     &-8.532322921D-1,3.460208861D-1/
      DATA L/1.574373327D+1,-3.417061978D+1,1.931380707D+1/
      DATA SL1/4.260321148D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATGSD'
      IF( INDX .EQ. 0 ) THEN
        RHOX = PX/((TX+TABS)*RCW)
      ELSEIF( INDX .EQ. 1 ) THEN
        IF( PX .LT. SMALL ) THEN
          RHOX = 0.D+0
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
        TR = (TX+TABS)/TCRW
        PR = PX/PCRW
        CX = EXP(SBX(1)*(1.D+0-TR))
        BL = L(1) + L(2)*TR + L(3)*TR*TR
!----------------------------------------------------------------------!
!     Compute the vapor density.
!----------------------------------------------------------------------!
        CHX = SL1*TR/PR
        CHX = CHX -(BX(7)*CX**13 + BX(8)*CX**3)
        CHX = CHX -2.D+0*PR*(BX(9)*CX**18 + BX(10)*CX**2 + BX(11)*CX)
        CHX = CHX -3.D+0*(PR**2)*(BX(12)*CX**18 + BX(13)*CX**10)
        CHX = CHX -4.D+0*(PR**3)*(BX(14)*CX**25 + BX(15)*CX**14)
        CHX = CHX -5.D+0*(PR**4)*
     &    (BX(16)*CX**32 + BX(17)*CX**28 + BX(18)*CX**24)
        CHX = CHX -4.D+0*(1.D+0/(PR**5))*(BX(19)*CX**12 +
     &    BX(20)*CX**11)/
     &    (((1.D+0/(PR**4))+(SBX(2)*CX**14))**2)
        CHX = CHX -5.D+0*(1.D+0/(PR**6))*(BX(21)*CX**24 +
     &    BX(22)*CX**18)/
     &    (((1.D+0/(PR**5))+(SBX(3)*CX**19))**2)
        CHX = CHX -6.D+0*(1.D+0/(PR**7))*(BX(23)*CX**24 +
     &    BX(24)*CX**14)/
     &    (((1.D+0/(PR**6))+(SBX(4)*CX**54 + SBX(5)*CX**27))**2)
        T1 = 1.1D+1*((PR/BL)**10)
        CHX = CHX + T1*BX(25)
        CHX = CHX + T1*BX(26)*CX
        CHX = CHX + T1*BX(27)*CX**2
        CHX = CHX + T1*BX(28)*CX**3
        CHX = CHX + T1*BX(29)*CX**4
        CHX = CHX + T1*BX(30)*CX**5
        CHX = CHX + T1*BX(31)*CX**6
        RHOX = 1.D+0/(CHX*VCRW*1.D-3/WTMW)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATGSD group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATLQD( TX,PX,RHO )
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
!     Calculate the subcooled or saturated density, as a function of
!     temperature and pressure per the steam table equations
!     as given by the 1967 International Formulation Committee:
!     Formulation for Industrial Use.
!
!     Thermodynamic and Transport Properties of Steam.
!     1967. ASME Steam Tables.
!     The American Society of Mechanical Engineers.
!     United Engineering Center, 345 East 47th Street, New York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < T > 364.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June 2022.
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
      REAL*8 A(23),SA(12),B(10)
!
!----------------------Data Statements---------------------------------!
!
      DATA A /6.824687741D+3,-5.422063673D+2,-2.096666205D+4,
     &3.941286787D+4,-6.733277739D+4,9.902381028D+4,-1.093911774D+5,
     &8.590841667D+4,-4.511168742D+4,1.418138926D+4,-2.017271113D+3,
     &7.982692717D+0,-2.616571843D-2,1.522411790D-3,2.284279054D-2,
     &2.421647003D+2,1.269716088D-10,2.074838328D-7,2.174020350D-8,
     &1.105710498D-9,1.293441934D+1,1.308119072D-5,6.047626338D-14/
      DATA SA /8.438375405D-1,5.362162162D-4,1.720D+0,7.342278489D-2,
     &4.975858870D-2,6.537154300D-1,1.15D-6,1.5108D-5,1.4188D-1,
     &7.002753165D+0,2.995284926D-4,2.040D-1/
      DATA B/9.99667D+2,6.85021D-2,-7.0966D-3,2.76483D-5,-5.4108D-8,
     &5.20175D-7,-7.41396D-9,1.41879D-10,-8.82877D-13,1.92152D-15/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATLQD'
!
!---  Constant aqueous density  ---
!
      IF( ISLC(9).EQ.1 .OR. ISLC(53).EQ.1 ) THEN
        RHO = RHOLI
!
!---  Polynomial formulation  ---
!
      ELSEIF( IOM.NE.3 .AND. IOM.NE.13 ) THEN
        RHO = B(1) + PX*B(6)
        DO M = 1,4
          RHO = RHO + B(M+1)*(TX**M) + PX*B(M+6)*(TX**M)
        ENDDO
!
!---  ASME formulation  ---
!
      ELSE
        TR = MIN( (TX+TABS)/TCRW,1.D+0 )
        PR = PX/PCRW
        YC = 1.0D+0 - SA(1)*TR*TR-SA(2)/(TR**6)
        ZC = YC +
     &    SQRT(MAX( ZERO,(SA(3)*YC*YC-2.D+0*SA(4)*TR+2.D+0*SA(5)*PR)))
        CHI = A(12)*SA(5)/ZC**(2.9412D-1) +A(13) +A(14)*TR +A(15)*TR*TR
     &    +A(16)*(SA(6)-TR)**10 +A(17)/(SA(7)+TR**19)
     &    -(A(18)+2.D+0*A(19)*PR+3.D+0*A(20)*PR*PR)/(SA(8)+TR**11)
     &    -(A(21)*TR**18*(SA(9)+TR*TR)*(-3.D+0/((SA(10)+PR)**4)+SA(11)))
     &    +3.D+0*A(22)*(SA(12)-TR)*PR*PR +4.D+0*A(23)*PR*PR*PR/(TR**20)
        RHO = 1.D+0/(CHI*VCRW*1.D-3/WTMW)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATLQD group  ---
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION WLQDF( TX )
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
!     Water liquid density function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
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
      REAL*8 B(10)
!
!----------------------Data Statements---------------------------------!
!
      DATA B/9.99667D+2,6.85021D-2,-7.0966D-3,2.76483D-5,-5.4108D-8,
     &5.20175D-7,-7.41396D-9,1.41879D-10,-8.82877D-13,1.92152D-15/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WLQDF'
!
!---  Water liquid density function  ---
!
      WLQDF = B(1)+B(2)*TX+B(3)*(TX**2)+B(4)*(TX**3)+B(5)*(TX**4)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WLQDF group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION WLQDG( TX )
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
!     Water liquid density function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
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
      REAL*8 B(10)
!
!----------------------Data Statements---------------------------------!
!
      DATA B/9.99667D+2,6.85021D-2,-7.0966D-3,2.76483D-5,-5.4108D-8,
     &5.20175D-7,-7.41396D-9,1.41879D-10,-8.82877D-13,1.92152D-15/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WLQDG'
!
!---  Water liquid density function  ---
!
      WLQDG = B(6)+B(7)*TX+B(8)*(TX**2)+B(9)*(TX**3)+B(10)*(TX**4)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WLQDG group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATLQV( TX,PX,PSWX,VIS )
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
!     Calculate the subcooled or saturated viscosity, as a function of
!     temperature and pressure per the steam table equations
!     as given by the 1967 International Formulation Committee:
!     Formulation for Industrial Use.
!
!     Thermodynamic and Transport Properties of Steam.
!     1967. ASME Steam Tables.
!     The American Society of Mechanical Engineers.
!     United Engineering Center, 345 East 47th Street, New York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < T > 364.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June 2022.
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
      REAL*8 K(4)
!
!----------------------Data Statements---------------------------------!
!
      DATA K/-2.471D+1, 4.209D+3, 4.527D-02, -3.376D-5/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATLQV'
      TK = MIN( TX+TABS,TCRW )
      DPX = MAX( PX,PSWX ) - PSWX
      PHI = 1.D+0 + 1.0467D+0*(TK-3.05D+2)*DPX/1.D+11
      VIS = PHI*EXP( K(1) + K(2)/TK + K(3)*TK + K(4)*TK*TK )*1.D-3
      IF( ISLC(9).EQ.1 ) VIS = VISLI
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATLQV group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATSP( TX,PX )
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
!     Calculate the saturation pressure of water as a function of
!     temperature per the 1967 International Formulation Committee
!     Formulation for Industrial Use.
!
!     Thermodynamic and Transport Properties of Steam.
!     1967. ASME Steam Tables.
!     The American Society of Mechanical Engineers.
!     United Engineering Center, 345 East 47th Street, New York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < T > 364.0 !
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 June 2022
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
      REAL*8 K(9)
!
!----------------------Data Statements---------------------------------!
!
      DATA K/-7.691234564D+0,-2.608023696D+1,-1.681706546D+2,
     &  6.423285504D+1,-1.189646225D+2,4.167117320D+0,2.097506760D+1,
     &  1.D+9,6.D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATSP'
      TR = MIN( (TX+TABS)/TCRW,1.D+0 )
      TR1 = 1.D+0-TR
      SVM = K(1)*(TR1**1) + K(2)*(TR1**2) + K(3)*(TR1**3)
     &  + K(4)*(TR1**4) + K(5)*(TR1**5)
      PR = EXP((SVM/(TR*(1.D+0+K(6)*TR1+K(7)*TR1*TR1)))
     &  + (TR1/(K(8)*TR1*TR1+K(9))))
      PX = PR*PCRW
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATSP group  ---
!
      RETURN
      END


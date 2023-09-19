!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWB( N,NB,NPZ )
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
!     SKiP -- Compute the water mass conservation equation residuals
!     and partial derivatives for bottom boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, December 1996.
!     Last Modified by MD White, PNNL, December 20, 1996.
!     jcbwb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RWP(LSV)
      REAL*8 VARL(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWB'
      IEQW2 = IEQW+2
      K = KD(N)
      FLWB = XLWB(2,NB)*RHOLB(2,NB)
      FLWP = XLW(2,N)*RHOL(2,N)
      INDX = 2
      CALL DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,NPZ),VARL,INDX )
      FLW = VARL(1)
      RWP(2) = -AFZ(NPZ)*(WL(1,NPZ)*FLW - WTMW*WDLA(1,NPZ)
     &  - WTMW*WDS(1,NPZ)/WTMS )
      IF( IEQW.GT.0 ) THEN
        MP = 2*(IEQW-1)+2
        FLW = VARL(2)*(XLW(IEQW2,N)*RHOL(2,N) + XLW(2,N)*RHOL(IEQW2,N))
     &    + VARL(3)*(XLWB(IEQW2,NB)*RHOLB(2,NB) +
     &    XLWB(2,NB)*RHOLB(IEQW2,NB))
        RWP(IEQW2) = -AFZ(NPZ)*(WL(MP,NPZ)*FLW -WTMW*WDLA(MP,NPZ)
     &    - WTMW*WDS(MP,NPZ)/WTMS)
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RWP,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWS( N,NB,NPY )
!
!-------------------------Notice---------------------------------------!
!
!     This software is being made available for testing and evaluation
!     purposes only.  Because it has not yet been authenticated and
!     cleared for publication, this software is not intended for
!     duplication or distribution to third parties without the
!     permission of the Manager of Software Products at Battelle,
!     Pacific Northwest Division, Richland, Washington, 99352.
!
!----------------------Description-------------------------------------!
!
!     SKiP -- Compute the water mass conservation equation residuals
!     and partial derivatives for south boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, December 1996.
!     Last Modified by MD White, PNNL, December 20, 1996.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWS'
      IEQW2 = IEQW+2
      J = JD(N)
      FLWB = XLWB(2,NB)*RHOLB(2,NB)
      FLWP = XLW(2,N)*RHOL(2,N)
      INDX = 2
      CALL DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,NPY),VARL,INDX )
      FLW = VARL(1)
      RWP(2) = -AFY(NPY)*(VL(1,NPY)*FLW - WTMW*VDLA(1,NPY)
     &  - WTMW*VDS(1,NPY)/WTMS )
      IF( IEQW.GT.0 ) THEN
        MP = 2*(IEQW-1)+2
        FLW = VARL(2)*(XLW(IEQW2,N)*RHOL(2,N) + XLW(2,N)*RHOL(IEQW2,N))
     &    + VARL(3)*(XLWB(IEQW2,NB)*RHOLB(2,NB) +
     &    XLWB(2,NB)*RHOLB(IEQW2,NB))
        RWP(IEQW2) = -AFY(NPY)*(VL(MP,NPY)*FLW -WTMW*VDLA(MP,NPY)
     &    - WTMW*VDS(MP,NPY)/WTMS)
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RWP,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWW( N,NB,NPX )
!
!-------------------------Notice---------------------------------------!
!
!     This software is being made available for testing and evaluation
!     purposes only.  Because it has not yet been authenticated and
!     cleared for publication, this software is not intended for
!     duplication or distribution to third parties without the
!     permission of the Manager of Software Products at Battelle,
!     Pacific Northwest Division, Richland, Washington, 99352.
!
!----------------------Description-------------------------------------!
!
!     SKiP -- Compute the water mass conservation equation residuals
!     and partial derivatives for west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, December 1996.
!     Last Modified by MD White, PNNL, December 20, 1996.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWW'
      IEQW2 = IEQW+2
      I = ID(N)
      FLWB = XLWB(2,NB)*RHOLB(2,NB)
      FLWP = XLW(2,N)*RHOL(2,N)
      INDX = 2
      CALL DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,NPX),VARL,INDX )
      FLW = VARL(1)
      RWP(2) = -AFX(NPX)*(UL(1,NPX)*FLW - WTMW*UDLA(1,NPX)
     &   - WTMW*UDS(1,NPX)/WTMS )
      IF( IEQW.GT.0 ) THEN
        MP = 2*(IEQW-1)+2
        FLW = VARL(2)*(XLW(IEQW2,N)*RHOL(2,N) + XLW(2,N)*RHOL(IEQW2,N))
     &    + VARL(3)*(XLWB(IEQW2,NB)*RHOLB(2,NB) +
     &    XLWB(2,NB)*RHOLB(IEQW2,NB))
        RWP(IEQW2) = -AFX(NPX)*(UL(MP,NPX)*FLW -WTMW*UDLA(MP,NPX)
     &    - WTMW*UDS(MP,NPX)/WTMS)
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RWP,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWE( N,NB,NQX )
!
!-------------------------Notice---------------------------------------!
!
!     This software is being made available for testing and evaluation
!     purposes only.  Because it has not yet been authenticated and
!     cleared for publication, this software is not intended for
!     duplication or distribution to third parties without the
!     permission of the Manager of Software Products at Battelle,
!     Pacific Northwest Division, Richland, Washington, 99352.
!
!----------------------Description-------------------------------------!
!
!     SKiP -- Compute the water mass conservation equation residuals
!     and partial derivatives for east boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, December 1996.
!     Last Modified by MD White, PNNL, December 20, 1996.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWE'
      IEQW2 = IEQW+2
      I = ID(N)
      FLWB = XLWB(2,NB)*RHOLB(2,NB)
      FLWP = XLW(2,N)*RHOL(2,N)
      INDX = 2
      CALL DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,NQX),VARL,INDX )
      FLW = VARL(1)
      RWP(2) = -AFX(NQX)*(UL(1,NQX)*FLW - WTMW*UDLA(1,NQX)
     &   - WTMW*UDS(1,NQX)/WTMS )
      IF( IEQW.GT.0 ) THEN
        MP = 2*(IEQW-1)+2
        FLW = VARL(3)*(XLW(IEQW2,N)*RHOL(2,N) + XLW(2,N)*RHOL(IEQW2,N))
     &    + VARL(2)*(XLWB(IEQW2,NB)*RHOLB(2,NB) +
     &    XLWB(2,NB)*RHOLB(IEQW2,NB))
        RWP(IEQW2) = -AFX(NQX)*(UL(MP,NQX)*FLW -WTMW*UDLA(MP,NQX)
     &    - WTMW*UDS(MP,NQX)/WTMS)
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RWP,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWN( N,NB,NQY )
!
!-------------------------Notice---------------------------------------!
!
!     This software is being made available for testing and evaluation
!     purposes only.  Because it has not yet been authenticated and
!     cleared for publication, this software is not intended for
!     duplication or distribution to third parties without the
!     permission of the Manager of Software Products at Battelle,
!     Pacific Northwest Division, Richland, Washington, 99352.
!
!----------------------Description-------------------------------------!
!
!     SKiP -- Compute the water mass conservation equation residuals
!     and partial derivatives for north boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, December 1996.
!     Last Modified by MD White, PNNL, December 20, 1996.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWN'
      IEQW2 = IEQW+2
      J = JD(N)
      FLWB = XLWB(2,NB)*RHOLB(2,NB)
      FLWP = XLW(2,N)*RHOL(2,N)
      INDX = 2
      CALL DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,NQY),VARL,INDX )
      FLW = VARL(1)
      RWP(2) = -AFY(NQY)*(VL(1,NQY)*FLW - WTMW*VDLA(1,NQY)
     &   - WTMW*VDS(1,NQY)/WTMS )
      IF( IEQW.GT.0 ) THEN
        MP = 2*(IEQW-1)+2
        FLW = VARL(3)*(XLW(IEQW2,N)*RHOL(2,N) + XLW(2,N)*RHOL(IEQW2,N))
     &    + VARL(2)*(XLWB(IEQW2,NB)*RHOLB(2,NB) +
     &    XLWB(2,NB)*RHOLB(IEQW2,NB))
        RWP(IEQW2) = -AFY(NQY)*(VL(MP,NQY)*FLW -WTMW*VDLA(MP,NQY)
     &    - WTMW*VDS(MP,NQY)/WTMS)
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RWP,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWT( N,NB,NQZ )
!
!-------------------------Notice---------------------------------------!
!
!     This software is being made available for testing and evaluation
!     purposes only.  Because it has not yet been authenticated and
!     cleared for publication, this software is not intended for
!     duplication or distribution to third parties without the
!     permission of the Manager of Software Products at Battelle,
!     Pacific Northwest Division, Richland, Washington, 99352.
!
!----------------------Description-------------------------------------!
!
!     SKiP -- Compute the water mass conservation equation residuals
!     and partial derivatives for top boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, December 1996.
!     Last Modified by MD White, PNNL, December 20, 1996.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWT'
      IEQW2 = IEQW+2
      K = KD(N)
      FLWB = XLWB(2,NB)*RHOLB(2,NB)
      FLWP = XLW(2,N)*RHOL(2,N)
      INDX = 2
      CALL DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,NQZ),VARL,INDX )
      FLW = VARL(1)
      RWP(2) = -AFZ(NQZ)*(WL(1,NQZ)*FLW - WTMW*WDLA(1,NQZ)
     &   - WTMW*WDS(1,NQZ)/WTMS )
      IF( IEQW.GT.0 ) THEN
        MP = 2*(IEQW-1)+2
        FLW = VARL(3)*(XLW(IEQW2,N)*RHOL(2,N) + XLW(2,N)*RHOL(IEQW2,N))
     &    + VARL(2)*(XLWB(IEQW2,NB)*RHOLB(2,NB) +
     &    XLWB(2,NB)*RHOLB(IEQW2,NB))
        RWP(IEQW2) = -AFZ(NQZ)*(WL(MP,NQZ)*FLW -WTMW*WDLA(MP,NQZ)
     &    - WTMW*WDS(MP,NQZ)/WTMS)
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RWP,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWT group  ---
!
      RETURN
      END

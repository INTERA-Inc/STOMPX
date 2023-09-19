!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAB( N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!     jcblab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/JCBLAB'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
        RS(M) = -AFZ(NPZ)*(WL(MP,NPZ)*FLA + WDLA(MP,NPZ)*WTMA)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAS( N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/JCBLAS'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
        RS(M) = -AFY(NPY)*(VL(MP,NPY)*FLA + VDLA(MP,NPY)*WTMA)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAW( N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/JCBLAW'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
        RS(M) = -AFX(NPX)*(UL(MP,NPX)*FLA + UDLA(MP,NPX)*WTMA)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAE( N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/JCBLAE'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
        RS(M) = AFX(NQX)*(UL(MN,NQX)*FLA + UDLA(MN,NQX)*WTMA)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAN( N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/JCBLAN'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
        RS(M) = AFY(NQY)*(VL(MN,NQY)*FLA + VDLA(MN,NQY)*WTMA)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAT( N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/JCBLAT'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
        RS(M) = AFZ(NQZ)*(WL(MN,NQZ)*FLA + WDLA(MN,NQZ)*WTMA )
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAT group  ---
!
      RETURN
      END

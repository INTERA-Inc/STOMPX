!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAB( N,NB,NPZ )
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
!     (gas boundary, air equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!     jcbgab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/JCBGAB'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
        RS(M) = -AFZ(NPZ)*(WG(MP,NPZ)*FLA - WTMA*WDGW(MP,NPZ) -
     &    WTMA*WDGO(MP,NPZ))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAS( N,NB,NPY )
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
!     (gas boundary, air equation, south surface)
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
      SUB_LOG(ISUB_LOG) = '/JCBGAS'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
        RS(M) = -AFY(NPY)*(VG(MP,NPY)*FLA - WTMA*VDGW(MP,NPY) -
     &    WTMA*VDGO(MP,NPY))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAW( N,NB,NPX )
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
!     (gas boundary, air equation, west surface)
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
      SUB_LOG(ISUB_LOG) = '/JCBGAW'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
        RS(M) = -AFX(NPX)*(UG(MP,NPX)*FLA - WTMA*UDGW(MP,NPX) -
     &    WTMA*UDGO(MP,NPX))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAE( N,NB,NQX )
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
!     (gas boundary, air equation, east surface)
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
      SUB_LOG(ISUB_LOG) = '/JCBGAE'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
        RS(M) = AFX(NQX)*(UG(MN,NQX)*FLA - WTMA*UDGW(MN,NQX) -
     &    WTMA*UDGO(MN,NQX))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAN( N,NB,NQY )
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
!     (gas boundary, air equation, north surface)
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
      SUB_LOG(ISUB_LOG) = '/JCBGAN'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
        RS(M) = AFY(NQY)*(VG(MN,NQY)*FLA - WTMA*VDGW(MN,NQY) -
     &    WTMA*VDGO(MN,NQY))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAT( N,NB,NQZ )
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
!     (gas boundary, air equation, top surface)
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
      SUB_LOG(ISUB_LOG) = '/JCBGAT'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
        RS(M) = AFZ(NQZ)*(WG(MN,NQZ)*FLA - WTMA*WDGW(MN,NQZ) -
     &    WTMA*WDGO(MN,NQZ))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAT group  ---
!
      RETURN
      END

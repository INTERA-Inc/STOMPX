!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLOB( N,NB,NPZ )
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
!     (aqueous boundary, VOC equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!     Last Modified by MD White, PNNL, 20 March 2003.
!     jcblob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/JCBLOB'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = -AFZ(NPZ)*WLO(MP,NPZ)
        ELSE
          FLOB = XLOB(MB,NB)*RHOLB(MB,NB)
          FLOP = XLO(MB,N)*RHOL(MB,N)
          INDX = 2
          FLO = DIFMN( FLOB,FLOP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
          RS(M) = -AFZ(NPZ)*(WL(MP,NPZ)*FLO + WDLO(MP,NPZ)*WTMO)
          WLO(MP,NPZ) = WL(MP,NPZ)*FLO + WDLO(MP,NPZ)*WTMO
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLOS( N,NB,NPY )
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
!     (aqueous boundary, VOC equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!     Last Modified by MD White, PNNL, 20 March 2003.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/JCBLOS'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = -AFY(NPY)*VLO(MP,NPY)
        ELSE
          FLOB = XLOB(MB,NB)*RHOLB(MB,NB)
          FLOP = XLO(MB,N)*RHOL(MB,N)
          INDX = 2
          FLO = DIFMN( FLOB,FLOP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
          RS(M) = -AFY(NPY)*(VL(MP,NPY)*FLO + VDLO(MP,NPY)*WTMO)
          VLO(MP,NPY) = VL(MP,NPY)*FLO + VDLO(MP,NPY)*WTMO
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLOW( N,NB,NPX )
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
!     (aqueous boundary, VOC equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!     Last Modified by MD White, PNNL, 20 March 2003.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/JCBLOW'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = -AFX(NPX)*ULO(MP,NPX)
        ELSE
          FLOB = XLOB(MB,NB)*RHOLB(MB,NB)
          FLOP = XLO(MB,N)*RHOL(MB,N)
          INDX = 2
          FLO = DIFMN( FLOB,FLOP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
          RS(M) = -AFX(NPX)*(UL(MP,NPX)*FLO + UDLO(MP,NPX)*WTMO)
          ULO(MP,NPX) = UL(MP,NPX)*FLO + UDLO(MP,NPX)*WTMO
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLOW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLOE( N,NB,NQX )
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
!     (aqueous boundary, VOC equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!     Last Modified by MD White, PNNL, 20 March 2003.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/JCBLOE'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = AFX(NQX)*ULO(MN,NQX)
        ELSE
          FLOB = XLOB(MB,NB)*RHOLB(MB,NB)
          FLOP = XLO(MB,N)*RHOL(MB,N)
          INDX = 2
          FLO = DIFMN( FLOP,FLOB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
          RS(M) = AFX(NQX)*(UL(MN,NQX)*FLO + UDLO(MN,NQX)*WTMO)
          ULO(MN,NQX) = UL(MN,NQX)*FLO + UDLO(MN,NQX)*WTMO
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLON( N,NB,NQY )
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
!     (aqueous boundary, VOC equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!     Last Modified by MD White, PNNL, 20 March 2003.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/JCBLON'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = AFY(NQY)*VLO(MN,NQY)
        ELSE
          FLOB = XLOB(MB,NB)*RHOLB(MB,NB)
          FLOP = XLO(MB,N)*RHOL(MB,N)
          INDX = 2
          FLO = DIFMN( FLOP,FLOB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
          RS(M) = AFY(NQY)*(VL(MN,NQY)*FLO + VDLO(MN,NQY)*WTMO)
          VLO(MN,NQY) = VL(MN,NQY)*FLO + VDLO(MN,NQY)*WTMO
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLOT( N,NB,NQZ )
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
!     (aqueous boundary, VOC equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!     Last Modified by MD White, PNNL, 20 March 2003.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/JCBLOT'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = AFZ(NQZ)*WLO(MN,NQZ)
        ELSE
          FLOB = XLOB(MB,NB)*RHOLB(MB,NB)
          FLOP = XLO(MB,N)*RHOL(MB,N)
          INDX = 2
          FLO = DIFMN( FLOP,FLOB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
          RS(M) = AFZ(NQZ)*(WL(MN,NQZ)*FLO + WDLO(MN,NQZ)*WTMO )
          WLO(MN,NQZ) = WL(MN,NQZ)*FLO + WDLO(MN,NQZ)*WTMO
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLOT group  ---
!
      RETURN
      END

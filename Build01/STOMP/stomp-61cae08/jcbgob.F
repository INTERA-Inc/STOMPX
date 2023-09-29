!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGOB( N,NB,NPZ )
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
!     (gas boundary, VOC equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, February 19, 1993.
!     jcbgob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/JCBGOB'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
!
!---    Leonard-TVD and Roe-Superbee schemes  ---
!
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = -AFZ(NPZ)*WGO(MP,NPZ)
        ELSE
          FGOB = XGOB(MB,NB)*RHOGB(MB,NB)
          FGOP = XGO(MB,N)*RHOG(MB,N)
          INDX = 3
          FGO = DIFMN( FGOB,FGOP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( WG(1,NPZ).LT.-EPSL ) THEN
              RS(M) = -AFZ(NPZ)*WG(MP,NPZ)*FGO
              WGO(MP,NPZ) = WG(MP,NPZ)*FGO
            ELSE
              RS(M) = 0.D+0
              WGO(MP,NPZ) = 0.D+0
            ENDIF
          ELSE
            RS(M) = -AFZ(NPZ)*(WG(MP,NPZ)*FGO + WTMO*WDGO(MP,NPZ))
            WGO(MP,NPZ) = WG(MP,NPZ)*FGO + WTMO*WDGO(MP,NPZ)
          ENDIF
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGOB group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGOS( N,NB,NPY )
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
!     (gas boundary, VOC equation, south surface)
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
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/JCBGOS'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
!
!---    Leonard-TVD and Roe-Superbee schemes  ---
!
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = -AFY(NPY)*VGO(MP,NPY)
        ELSE
          FGOB = XGOB(MB,NB)*RHOGB(MB,NB)
          FGOP = XGO(MB,N)*RHOG(MB,N)
          INDX = 3
          FGO = DIFMN( FGOB,FGOP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( VG(1,NPY).LT.-EPSL ) THEN
              RS(M) = -AFY(NPY)*VG(MP,NPY)*FGO
              VGO(MP,NPY) = VG(MP,NPY)*FGO
            ELSE
              RS(M) = 0.D+0
              VGO(MP,NPY) = 0.D+0
            ENDIF
          ELSE
            RS(M) = -AFY(NPY)*(VG(MP,NPY)*FGO + WTMO*VDGO(MP,NPY))
            VGO(MP,NPY) = VG(MP,NPY)*FGO + WTMO*VDGO(MP,NPY)
          ENDIF
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGOS group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGOW( N,NB,NPX )
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
!     (gas boundary, VOC equation, west surface)
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
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/JCBGOW'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
!
!---    Leonard-TVD and Roe-Superbee schemes  ---
!
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = -AFX(NPX)*UGO(MP,NPX)
        ELSE
          FGOB = XGOB(MB,NB)*RHOGB(MB,NB)
          FGOP = XGO(MB,N)*RHOG(MB,N)
          INDX = 3
          FGO = DIFMN( FGOB,FGOP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( UG(1,NPX).LT.-EPSL ) THEN
              RS(M) = -AFX(NPX)*UG(MP,NPX)*FGO
              UGO(MP,NPX) = UG(MP,NPX)*FGO
            ELSE
              RS(M) = 0.D+0
              UGO(MP,NPX) = 0.D+0
            ENDIF
          ELSE
            RS(M) = -AFX(NPX)*(UG(MP,NPX)*FGO + WTMO*UDGO(MP,NPX))
            UGO(MP,NPX) = UG(MP,NPX)*FGO + WTMO*UDGO(MP,NPX)
          ENDIF
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGOW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGOE( N,NB,NQX )
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
!     (gas boundary, VOC equation, east surface)
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
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/JCBGOE'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
!
!---    Leonard-TVD and Roe-Superbee schemes  ---
!
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = AFX(NQX)*UGO(MN,NQX)
        ELSE
          FGOB = XGOB(MB,NB)*RHOGB(MB,NB)
          FGOP = XGO(MB,N)*RHOG(MB,N)
          INDX = 3
          FGO = DIFMN( FGOP,FGOB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( UG(1,NQX).GT.EPSL ) THEN
              RS(M) = AFX(NQX)*UG(MN,NQX)*FGO
              UGO(MN,NQX) = UG(MN,NQX)*FGO
            ELSE
              RS(M) = 0.D+0
              UGO(MN,NQX) = 0.D+0
            ENDIF
          ELSE
            RS(M) = AFX(NQX)*(UG(MN,NQX)*FGO + WTMO*UDGO(MN,NQX))
            UGO(MN,NQX) = UG(MN,NQX)*FGO + WTMO*UDGO(MN,NQX)
          ENDIF
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGON( N,NB,NQY )
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
!     (gas boundary, VOC equation, north surface)
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
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/JCBGON'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
!
!---    Leonard-TVD and Roe-Superbee schemes  ---
!
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = AFY(NQY)*VGO(MN,NQY)
        ELSE
          FGOB = XGOB(MB,NB)*RHOGB(MB,NB)
          FGOP = XGO(MB,N)*RHOG(MB,N)
          INDX = 3
          FGO = DIFMN( FGOP,FGOB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( VG(1,NQY).GT.EPSL ) THEN
              RS(M) = AFY(NQY)*VG(MN,NQY)*FGO
              VGO(MN,NQY) = VG(MN,NQY)*FGO
            ELSE
              RS(M) = 0.D+0
              VGO(MN,NQY) = 0.D+0
            ENDIF
          ELSE
            RS(M) = AFY(NQY)*(VG(MN,NQY)*FGO + WTMO*VDGO(MN,NQY))
            VGO(MN,NQY) = VG(MN,NQY)*FGO + WTMO*VDGO(MN,NQY)
          ENDIF
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGOT( N,NB,NQZ )
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
!     (gas boundary, VOC equation, top surface)
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
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/JCBGOT'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
!
!---    Leonard-TVD and Roe-Superbee schemes  ---
!
        IF( ISLC(12).EQ.1 ) THEN
          RS(M) = AFZ(NQZ)*WGO(MN,NQZ)
        ELSE
          FGOB = XGOB(MB,NB)*RHOGB(MB,NB)
          FGOP = XGO(MB,N)*RHOG(MB,N)
          INDX = 3
          FGO = DIFMN( FGOP,FGOB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( WG(1,NQZ).GT.EPSL ) THEN
              RS(M) = AFZ(NQZ)*WG(MN,NQZ)*FGO
              WGO(MN,NQZ) = WG(MN,NQZ)*FGO
            ELSE
              RS(M) = 0.D+0
              WGO(MN,NQZ) = 0.D+0
            ENDIF
          ELSE
            RS(M) = AFZ(NQZ)*(WG(MN,NQZ)*FGO + WTMO*WDGO(MN,NQZ))
            WGO(MN,NQZ) = WG(MN,NQZ)*FGO + WTMO*WDGO(MN,NQZ)
          ENDIF
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGOT group  ---
!
      RETURN
      END

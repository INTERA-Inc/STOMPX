!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SHDPN( N,DPB,DPS,DPW,DPE,DPN,DPT )
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
!     Calculates hydrodynamic dispersion coefficients for the NAPL
!     phase from phase velocities and user-specified dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by ML Rockhold, Battelle, PNL, March 1993.
!     Last Modified by MD White, Battelle, PNL, July 16, 1993.
!     shdpn.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE FLUXN
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/SHDPN'
      M = 1
      DPB = 0.D+0
      DPS = 0.D+0
      DPW = 0.D+0
      DPE = 0.D+0
      DPN = 0.D+0
      DPT = 0.D+0
      IF( IDISP .EQ. 0 ) GOTO 170

!
!---  Define indexing for grid block surfaces  ---
!
      NW = N-1
      NE = N+1
      NS = N-IFLD
      NN = N+IFLD
      NB = N-IJFLD
      NT = N+IJFLD
      I = ID(N)
      J = JD(N)
      K = KD(N)
!
!---  Bottom face  ---
!
      IF( K .NE. 1 ) THEN
        IF( IXP(NB) .LE. 0 ) GOTO 110
        CALL ADVB( PORD,SN,UN,VN,WN,UNBX,VNBX,WNBX,N,M )
        UNBSQ = UNBX*UNBX
        VNBSQ = VNBX*VNBX
        WNBSQ = WNBX*WNBX
        ZVB = SQRT(UNBSQ+VNBSQ+WNBSQ)
        INDX = 17
        DPLB = DIFMN(DISPL(IZ(NB)),DISPL(IZ(N)),DZGF(NB),DZGF(N),
     &    WNBX,INDX)
        DPTB = DIFMN(DISPT(IZ(NB)),DISPT(IZ(N)),DZGF(NB),DZGF(N),
     &    WNBX,INDX)
        DPB = (DPLB*WNBSQ + DPTB*(VNBSQ+UNBSQ))/(ZVB+SMALL)
      ENDIF
  110 CONTINUE
!
!---  South face  ---
!
      IF( J .NE. 1 ) THEN
        IF( IXP(NS) .LE. 0 ) GOTO 120
        CALL ADVS( PORD,SN,UN,VN,WN,UNSX,VNSX,WNSX,N,M )
        UNSSQ = UNSX*UNSX
        VNSSQ = VNSX*VNSX
        WNSSQ = WNSX*WNSX
        ZVS = SQRT(UNSSQ+VNSSQ+WNSSQ)
        INDX = 17
        DPLS = DIFMN(DISPL(IZ(NS)),DISPL(IZ(N)),DYGF(NS),DYGF(N),
     &    VNSX,INDX)
        DPTS = DIFMN(DISPT(IZ(NS)),DISPT(IZ(N)),DYGF(NS),DYGF(N),
     &    VNSX,INDX)
        DPS = (DPLS*VNSSQ + DPTS*(UNSSQ+WNSSQ))/(ZVS+SMALL)
      ENDIF
  120 CONTINUE
!
!---  West face  ---
!
      IF( I .NE. 1 ) THEN
        IF( IXP(NW) .LE. 0 ) GOTO 130
        CALL ADVW( PORD,SN,UN,VN,WN,UNWX,VNWX,WNWX,N,M )
        UNWSQ = UNWX*UNWX
        VNWSQ = VNWX*VNWX
        WNWSQ = WNWX*WNWX
        ZVW = SQRT(UNWSQ+VNWSQ+WNWSQ)
        INDX = 17
        DPLW = DIFMN(DISPL(IZ(NW)),DISPL(IZ(N)),DXGF(NW),DXGF(N),
     &    UNWX,INDX)
        DPTW = DIFMN(DISPT(IZ(NW)),DISPT(IZ(N)),DXGF(NW),DXGF(N),
     &    UNWX,INDX)
        DPW = (DPLW*UNWSQ + DPTW*(VNWSQ+WNWSQ))/(ZVW+SMALL)
      ENDIF
  130 CONTINUE

!
!---  East face  ---
!
      IF( I .NE. IFLD ) THEN
        IF( IXP(NE) .LE. 0 ) GOTO 140
        CALL ADVE( PORD,SN,UN,VN,WN,UNEX,VNEX,WNEX,N,M )
        UNESQ = UNEX*UNEX
        VNESQ = VNEX*VNEX
        WNESQ = WNEX*WNEX
        ZVE = SQRT(UNESQ+VNESQ+WNESQ)
        INDX = 17
        DPLE = DIFMN(DISPL(IZ(N)),DISPL(IZ(NE)),DXGF(N),DXGF(NE),
     &    UNEX,INDX)
        DPTE = DIFMN(DISPT(IZ(N)),DISPT(IZ(NE)),DXGF(N),DXGF(NE),
     &    UNEX,INDX)
        DPE = (DPLE*UNESQ + DPTE*(VNESQ+WNESQ))/(ZVE+SMALL)
      ENDIF
  140 CONTINUE
!
!---  North face  ---
!
      IF( J .NE. JFLD ) THEN
        IF( IXP(NN) .LE. 0 ) GOTO 150
        CALL ADVN( PORD,SN,UN,VN,WN,UNNX,VNNX,WNNX,N,M )
        UNNSQ = UNNX*UNNX
        VNNSQ = VNNX*VNNX
        WNNSQ = WNNX*WNNX
        ZVN = SQRT(UNNSQ+VNNSQ+WNNSQ)
        INDX = 17
        DPLN = DIFMN(DISPL(IZ(N)),DISPL(IZ(NN)),DYGF(N),DYGF(NN),
     &    VNNX,INDX)
        DPTN = DIFMN(DISPT(IZ(N)),DISPT(IZ(NN)),DYGF(N),DYGF(NN),
     &    VNNX,INDX)
        DPN = (DPLN*VNNSQ + DPTN*(UNNSQ+WNNSQ))/(ZVN+SMALL)
      ENDIF
  150 CONTINUE
!
!---  Top face  ---
!
      IF( K .NE. KFLD ) THEN
        IF( IXP(NT) .LE. 0 ) GOTO 160
        CALL ADVT( PORD,SN,UN,VN,WN,UNTX,VNTX,WNTX,N,M )
        UNTSQ = UNTX*UNTX
        VNTSQ = VNTX*VNTX
        WNTSQ = WNTX*WNTX
        ZVT = SQRT(UNTSQ+VNTSQ+WNTSQ)
        INDX = 17
        DPLT = DIFMN(DISPL(IZ(N)),DISPL(IZ(NT)),DZGF(N),DZGF(NT),
     &    WNTX,INDX)
        DPTT = DIFMN(DISPT(IZ(N)),DISPT(IZ(NT)),DZGF(N),DZGF(NT),
     &    WNTX,INDX)
        DPT = (DPLT*WNTSQ + DPTT*(VNTSQ+UNTSQ))/(ZVT+SMALL)
      ENDIF
  160 CONTINUE
  170 CONTINUE
!
!---  End of SHDPN group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END




!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLO
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
!     Compute dissolved oil fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!     Last Modified by MD White, PNNL, October 27, 1997.
!     dfflo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLO'
!
!---  X-direction molar diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
      DO 200 J = 1,JFLD
      DO 200 I = 2,IFLD
        N = ND(I,J,K)
        NW = N-1
        IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &    INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
        NPX = NSX(N)
        DXMLO = XMLO(2,NW)-XMLO(2,N)
        DO 100 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
          DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFSW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*DFLO(MN,NW)
     &      *RHOML(MN,NW)
          INDX = 15
          DFM = DIFMN( DFSW,DFSP,DXGF(NW),DXGF(N),DXMLO,INDX )
          UDLO(M,NPX) = DFM*(XMLO(MN,NW)-XMLO(MP,N))
     &      /DXGP(NPX)
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOW = XLO(MN,NW)*RHOL(MN,NW)
          INDX = 2
          FLO = DIFMN( FLOW,FLOP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
          ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
  100   CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
      DO 400 I = 1,IFLD
      DO 400 K = 1,KFLD
      DO 400 J = 2,JFLD
        N = ND(I,J,K)
        NS = N-IFLD
        IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &    INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
        NPY = NSY(N)
        DXMLO = XMLO(2,NS)-XMLO(2,N)
        DO 300 M = 1,ISVF
         MN = MNEG(M)
         MP = MPOS(M)
         DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &     *RHOML(MP,N)
         DFSS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*DFLO(MN,NS)
     &     *RHOML(MN,NS)
         INDX = 15
         DFM = DIFMN( DFSS,DFSP,DYGF(NS),DYGF(N),DXMLO,INDX )
         VDLO(M,NPY) = DFM*(XMLO(MN,NS)-XMLO(MP,N))
     &     /DYGP(NPY)/RP(I)
         FLOP = XLO(MP,N)*RHOL(MP,N)
         FLOS = XLO(MN,NS)*RHOL(MN,NS)
         INDX = 2
         FLO = DIFMN( FLOS,FLOP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
         VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
  300   CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
      DO 600 J = 1,JFLD
      DO 600 I = 1,IFLD
      DO 600 K = 2,KFLD
        N = ND(I,J,K)
        NB = N-IJFLD
        IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &    INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
        NPZ = NSZ(N)
        DXMLO = XMLO(2,NB)-XMLO(2,N)
        DO 500 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
          DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFSB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*DFLO(MN,NB)
     &      *RHOML(MN,NB)
          INDX = 15
          DFM = DIFMN( DFSB,DFSP,DZGF(NB),DZGF(N),DXMLO,INDX )
          WDLO(M,NPZ) = DFM*(XMLO(MN,NB)-XMLO(MP,N))
     &      /DZGP(NPZ)
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOB = XLO(MN,NB)*RHOL(MN,NB)
          INDX = 2
          FLO = DIFMN( FLOB,FLOP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
          WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
  500   CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLO group  ---
!
      RETURN
      END

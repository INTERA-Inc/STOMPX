!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS33
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
!     Load the Jacobian matrix for the salt equation aqueous-phase
!     contributions (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July, 1995.
!     Last Modified by MD White, Battelle, PNL, October 19, 1999.
!     jcbs_co2e.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STS(LUK+1),RSP(LUK),RSA(LUK,6),FS(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBS33'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 1000 K = 1,KFLD
        DO 1000 J = 1,JFLD
          DO 1000 I = 1,IFLD
        N = ND(I,J,K)
        IF( IXP(N).EQ.0 ) GOTO 1000
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX + 1
        NQY = NPY + IFLD
        NQZ = NPZ + IJFLD
        DO 10 M = 1,ISVC+1
          MP = M + 1
          STS(M) = (TMS(MP,N)-TMS(1,N))*VOL(N)*DTI
   10   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 30 MD = 1,6
          DO 20 M = 1,ISVF
            FS(M,MD) = 0.D+0
   20     CONTINUE
   30   CONTINUE
!
!---  Bottom surface flux component  ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 110
          DO 100 M = 1,ISVF
            FS(M,1) = -AFZ(NPZ)*WS(M,NPZ)
  100     CONTINUE
        ENDIF
  110   CONTINUE
!
!---  South surface flux component ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 210
          DO 200 M = 1,ISVF
            FS(M,2) = -AFY(NPY)*VS(M,NPY)
  200     CONTINUE
        ENDIF
  210   CONTINUE
!
!---  West surface flux component ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 310
          DO 300 M = 1,ISVF
            FS(M,3) = -AFX(NPX)*US(M,NPX)
  300     CONTINUE
        ENDIF
  310   CONTINUE
!
!---  East surface flux component  ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 410
          DO 400 M = 1,ISVF
            MF = MFLX(M)
            FS(M,4) = AFX(NQX)*US(MF,NQX)
  400     CONTINUE
        ENDIF
  410   CONTINUE
!
!---  North surface flux component  ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 510
          DO 500 M = 1,ISVF
            MF = MFLX(M)
            FS(M,5) = AFY(NQY)*VS(MF,NQY)
  500     CONTINUE
        ENDIF
  510   CONTINUE
!
!---  Top surface flux component  ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 610
          DO 600 M = 1,ISVF
            MF = MFLX(M)
            FS(M,6) = AFZ(NQZ)*WS(MF,NQZ)
  600     CONTINUE
        ENDIF
  610   CONTINUE
!
!---  Salt equation residuals  ---
!
        RSS = STS(1) - SRCS(2,N)
        DO 900 MD = 1,6
          RSS = RSS + FS(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RSP(M) = STS(M+1) - SRCS(M+2,N)
          MM = 2*M
          DO 910 MD = 1,6
            RSP(M) = RSP(M) + FS(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,6
            RSA(M,MD) = RSS - FS(1,MD) + FS(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Jacobian matrix loader  --
!
        CALL JCBBM( RSS,RSP,RSA,N,I,J,K,IEQS )
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS33 group  ---
!
      RETURN
      END

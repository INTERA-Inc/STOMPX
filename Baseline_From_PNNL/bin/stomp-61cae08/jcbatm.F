!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBATM
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
!     Load the Jacobian matrix for the alcohol equation for
!     an aqueous-phase/NAPL ternary mixture system
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, October 19, 1999.
!     jcbatm.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STA(LUK+1),RAP(LUK),RAA(LUK,6),FA(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBATM'
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
!
!---    First-order, forward-difference, time differential  ---
!
        STA1 = PORD(1,N)*(XLA(1,N)*RHOL(1,N)*SL(1,N) +
     &    XNA(1,N)*RHON(1,N)*SN(1,N))
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STA0 = PORD(MP,N)*(XLA(MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XNA(1,N)*RHON(MP,N)*SN(MP,N))
          STA(M) = (STA0-STA1)*DTI*VOL(N)
  100   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 210 MD = 1,6
          DO 200 M = 1,ISVF
            FA(M,MD) = 0.D+0
  200     CONTINUE
  210   CONTINUE
!
!---  Compute surface fluxes  ---
!
!---  Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 310
          DO 300 M = 1,ISVF
            MB = MADJ(M)
            MP = MNOD(M)
            FNAB = XNA(MB,NB)*RHON(MB,NB)
            FNAP = XNA(MP,N)*RHON(MP,N)
            INDX = 4
            FNA = DIFMN( FNAB,FNAP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            FLAB = XLA(MB,NB)*RHOL(MB,NB)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAB,FLAP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            FA(M,1) = -AFZ(NPZ)*(WL(M,NPZ)*FLA + WN(M,NPZ)*FNA +
     &        (WDLA(M,NPZ)+WDGO(M,NPZ))*WTMA)
  300     CONTINUE
        ENDIF
  310   CONTINUE
!
!---  South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 410
          DO 400 M = 1,ISVF
            MS = MADJ(M)
            MP = MNOD(M)
            FNAS = XNA(MS,NS)*RHON(MS,NS)
            FNAP = XNA(MP,N)*RHON(MP,N)
            INDX = 4
            FNA = DIFMN( FNAS,FNAP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            FLAS = XLA(MS,NS)*RHOL(MS,NS)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAS,FLAP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            FA(M,2) = -AFY(NPY)*(VL(M,NPY)*FLA + VN(M,NPY)*FNA +
     &        (VDLA(M,NPY)+VDGO(M,NPY))*WTMA)
  400     CONTINUE
        ENDIF
  410   CONTINUE
!
!---  West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 510
          DO 500 M = 1,ISVF
            MW = MADJ(M)
            MP = MNOD(M)
            FNAW = XNA(MW,NW)*RHON(MW,NW)
            FNAP = XNA(MP,N)*RHON(MP,N)
            INDX = 4
            FNA = DIFMN( FNAW,FNAP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            FLAW = XLA(MW,NW)*RHOL(MW,NW)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAW,FLAP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            FA(M,3) = -AFX(NPX)*(UL(M,NPX)*FLA + UN(M,NPX)*FNA +
     &        (UDLA(M,NPX)+UDGO(M,NPX))*WTMA)
  500     CONTINUE
        ENDIF
  510   CONTINUE
!
!---  East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 610
          DO 600 M = 1,ISVF
            ME = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FNAP = XNA(MP,N)*RHON(MP,N)
            FNAE = XNA(ME,NE)*RHON(ME,NE)
            INDX = 4
            FNA = DIFMN( FNAP,FNAE,DXGF(N),DXGF(NE),UN(1,NQX),INDX )
            FLAE = XLA(ME,NE)*RHOL(ME,NE)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAP,FLAE,DXGF(N),DXGF(NE),UL(1,NQX),INDX )
            FA(M,4) = AFX(NQX)*(UL(MF,NQX)*FLA + UN(MF,NQX)*FNA +
     &        (UDLA(MF,NQX)+UDGO(MF,NQX))*WTMA)
  600     CONTINUE
        ENDIF
  610   CONTINUE
!
!---  North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 710
          DO 700 M = 1,ISVF
            MN = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FNAP = XNA(MP,N)*RHON(MP,N)
            FNAN = XNA(MN,NN)*RHON(MN,NN)
            INDX = 4
            FNA = DIFMN( FNAP,FNAN,DYGF(N),DYGF(NN),VN(1,NQY),INDX )
            FLAN = XLA(MN,NN)*RHOL(MN,NN)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAP,FLAN,DYGF(N),DYGF(NN),VL(1,NQY),INDX )
            FA(M,5) = AFY(NQY)*(VL(MF,NQY)*FLA + VN(MF,NQY)*FNA +
     &        (VDLA(MF,NQY)+VDGO(MF,NQY))*WTMA)
  700     CONTINUE
        ENDIF
  710   CONTINUE
!
!---  Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 810
          DO 800 M = 1,ISVF
            MT = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FNAP = XNA(MP,N)*RHON(MP,N)
            FNAT = XNA(MT,NT)*RHON(MT,NT)
            INDX = 4
            FNA = DIFMN( FNAP,FNAT,DZGF(N),DZGF(NT),WN(1,NQZ),INDX )
            FLAT = XLA(MT,NT)*RHOL(MT,NT)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAP,FLAT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX )
            FA(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLA + WN(MF,NQZ)*FNA +
     &        (WDLA(MF,NQZ)+WDGO(MF,NQZ))*WTMA)
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute VOC equation residuals  ---
!
        RAS = STA(1) - SRCS(2,N)
        DO 900 MD = 1,6
          RAS = RAS + FA(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RAP(M) = STA(M+1) - SRCS(M+2,N)
          MM = 2*M
          DO 910 MD = 1,6
            RAP(M) = RAP(M) + FA(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,6
            RAA(M,MD) = RAS - FA(1,MD) + FA(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBM( RAS,RAP,RAA,N,I,J,K,IEQALC )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBATM group
!
      RETURN
      END

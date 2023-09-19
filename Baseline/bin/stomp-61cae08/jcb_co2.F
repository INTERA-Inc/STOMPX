!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBA32
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
!     Load the Jacobian matrix for the air equation with
!     aqueous-phase and gas-phase contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 May 2002.
!     Last Modified by MD White, PNNL, 29 May 2002.
!     jcb_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/JCBA32'
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
     &      XGA(1,N)*RHOG(1,N)*SG(1,N))
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STA0 = PORD(MP,N)*(XLA(MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XGA(MP,N)*RHOG(MP,N)*SG(MP,N))
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
            FLAB = XLA(MB,NB)*RHOL(MB,NB)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAB,FLAP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            FGAB = XGA(MB,NB)*RHOG(MB,NB)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAB,FGAP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            FA(M,1) = -AFZ(NPZ)*(WL(M,NPZ)*FLA + WG(M,NPZ)*FGA
     &        + WTMA*WDGA(M,NPZ) + WTMA*WDLA(M,NPZ))
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
            FLAS = XLA(MS,NS)*RHOL(MS,NS)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAS,FLAP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            FGAS = XGA(MS,NS)*RHOG(MS,NS)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAS,FGAP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            FA(M,2) = -AFY(NPY)*(VL(M,NPY)*FLA + VG(M,NPY)*FGA
     &        + WTMA*VDGA(M,NPY) + WTMA*VDLA(M,NPY))
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
            FLAW = XLA(MW,NW)*RHOL(MW,NW)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAW,FLAP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            FGAW = XGA(MW,NW)*RHOG(MW,NW)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAW,FGAP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            FA(M,3) = -AFX(NPX)*(UL(M,NPX)*FLA + UG(M,NPX)*FGA
     &        + WTMA*UDGA(M,NPX) + WTMA*UDLA(M,NPX))
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
            FLAE = XLA(ME,NE)*RHOL(ME,NE)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAP,FLAE,DXGF(N),DXGF(NE),UL(1,NQX),INDX )
            FGAE = XGA(ME,NE)*RHOG(ME,NE)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAP,FGAE,DXGF(N),DXGF(NE),UG(1,NQX),INDX )
            FA(M,4) = AFX(NQX)*(UL(MF,NQX)*FLA + UG(MF,NQX)*FGA
     &        + WTMA*UDGA(MF,NQX) + WTMA*UDLA(MF,NQX))
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
            FLAN = XLA(MN,NN)*RHOL(MN,NN)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAP,FLAN,DYGF(N),DYGF(NN),VL(1,NQY),INDX )
            FGAN = XGA(MN,NN)*RHOG(MN,NN)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAP,FGAN,DYGF(N),DYGF(NN),VG(1,NQY),INDX )
            FA(M,5) = AFY(NQY)*(VL(MF,NQY)*FLA + VG(MF,NQY)*FGA
     &       + WTMA*VDGA(MF,NQY) + WTMA*VDLA(MF,NQY))
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
            FLAT = XLA(MT,NT)*RHOL(MT,NT)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAP,FLAT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX )
            FGAT = XGA(MT,NT)*RHOG(MT,NT)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAP,FGAT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX )
            FA(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLA + WG(MF,NQZ)*FGA
     &        + WTMA*WDGA(MF,NQZ) + WTMA*WDLA(MF,NQZ))
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute air equation residuals  ---
!
        RAS = STA(1) - SRCA(2,N)
        DO 900 MD = 1,6
          RAS = RAS + FA(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RAP(M) = STA(M+1) - SRCA(M+2,N)
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
        CALL JCBBM( RAS,RAP,RAA,N,I,J,K,IEQA )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBA32 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS32
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
!     jcb_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/JCBS32'
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
!---  End of JCBS32 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW32
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
!     Load the Jacobian matrix for the water equation with
!     aqueous-phase, gas-phase and salt contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 May 2002.
!     Last Modified by MD White, PNNL, 29 May 2002.
!     jcb_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE FDVP
      USE FDVI
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STWX(LUK+1),RWP(LUK),RWA(LUK,6),FW(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBW32'
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
        STW1 = PORD(1,N)*(XLW(1,N)*RHOL(1,N)*SL(1,N) +
     &      XGW(1,N)*RHOG(1,N)*SG(1,N))
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STW0 = PORD(MP,N)*(XLW(MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XGW(MP,N)*RHOG(MP,N)*SG(MP,N))
          STWX(M) = (STW0-STW1)*DTI*VOL(N)
  100   CONTINUE
!
!---  Freezing conditions  ---
!
        IF( ISLC(5).EQ.1 ) THEN
!
!---      First-order, forward-difference, time differential  ---
!
          DO 104 M = 1,ISVC+1
            MP = M + 1
            STWX(M) = STWX(M) + (PORD(MP,N)*RHOI(MP,N)*SI(MP,N) -
     &        PORD(1,N)*RHOI(1,N)*SI(1,N))*DTI*VOL(N)
  104     CONTINUE
        ENDIF
!
!---  Initialize surface fluxes  ---
!
        DO 210 MD = 1,6
          DO 200 M = 1,ISVF
            FW(M,MD) = 0.D+0
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
            FLWB = XLW(MB,NB)*RHOL(MB,NB)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            FGWB = XGW(MB,NB)*RHOG(MB,NB)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWB,FGWP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            FW(M,1) = -AFZ(NPZ)*(WL(M,NPZ)*FLW + WG(M,NPZ)*FGW
     &        - WTMW*WDGA(M,NPZ) - WTMW*WDLA(M,NPZ)
     &        - WTMW*WDS(M,NPZ)/WTMS)
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
            FLWS = XLW(MS,NS)*RHOL(MS,NS)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWS,FLWP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            FGWS = XGW(MS,NS)*RHOG(MS,NS)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWS,FGWP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            FW(M,2) = -AFY(NPY)*(VL(M,NPY)*FLW + VG(M,NPY)*FGW
     &        - WTMW*VDGA(M,NPY) - WTMW*VDLA(M,NPY)
     &        - WTMW*VDS(M,NPY)/WTMS)
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
            FLWW = XLW(MW,NW)*RHOL(MW,NW)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            FGWW = XGW(MW,NW)*RHOG(MW,NW)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            FW(M,3) = -AFX(NPX)*(UL(M,NPX)*FLW + UG(M,NPX)*FGW
     &        - WTMW*UDGA(M,NPX) - WTMW*UDLA(M,NPX)
     &        - WTMW*UDS(M,NPX)/WTMS)
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
            FLWE = XLW(ME,NE)*RHOL(ME,NE)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWE,DXGF(N),DXGF(NE),UL(1,NQX),INDX )
            FGWE = XGW(ME,NE)*RHOG(ME,NE)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWP,FGWE,DXGF(N),DXGF(NE),UG(1,NQX),INDX )
            FW(M,4) = AFX(NQX)*(UL(MF,NQX)*FLW + UG(MF,NQX)*FGW
     &        - WTMW*UDGA(MF,NQX) - WTMW*UDLA(MF,NQX)
     &        - WTMW*UDS(MF,NQX)/WTMS)
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
            FLWN = XLW(MN,NN)*RHOL(MN,NN)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWN,DYGF(N),DYGF(NN),VL(1,NQY),INDX )
            FGWN = XGW(MN,NN)*RHOG(MN,NN)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWP,FGWN,DYGF(N),DYGF(NN),VG(1,NQY),INDX )
            FW(M,5) = AFY(NQY)*(VL(MF,NQY)*FLW + VG(MF,NQY)*FGW
     &        - WTMW*VDGA(MF,NQY) - WTMW*VDLA(MF,NQY)
     &        - WTMW*VDS(MF,NQY)/WTMS)
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
            FLWT = XLW(MT,NT)*RHOL(MT,NT)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX )
            FGWT = XGW(MT,NT)*RHOG(MT,NT)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWP,FGWT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX )
            FW(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLW + WG(MF,NQZ)*FGW
     &        - WTMW*WDGA(MF,NQZ) - WTMW*WDLA(MF,NQZ)
     &        - WTMW*WDS(MF,NQZ)/WTMS)
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute water equation residuals  ---
!
        RWS = STWX(1) - SRCW(2,N)
        DO 900 MD = 1,6
          RWS = RWS + FW(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RWP(M) = STWX(M+1) - SRCW(M+2,N)
          MM = 2*M
          DO 910 MD = 1,6
            RWP(M) = RWP(M) + FW(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,6
            RWA(M,MD) = RWS - FW(1,MD) + FW(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Load Jacobian Matrix  ---
!
        CALL JCBBM( RWS,RWP,RWA,N,I,J,K,IEQW )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW32 group
!
      RETURN
      END

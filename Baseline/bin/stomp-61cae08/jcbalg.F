!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBALG
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
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, October 19, 1999.
!     jcbalg.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/JCBALG'
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
        NB = N-IJFLD
        NS = N-IFLD
        NW = N-1
        NE = N+1
        NN = N+IFLD
        NT = N+IJFLD
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
!
!---  TVD air transport  ---
!
            IF( ISLC(10).EQ.1 ) THEN
              FA(M,1) = -AFZ(NPZ)*WTMA*(WDLA(M,NPZ)-WDGW(M,NPZ))
              IF( WL(1,NPZ).GE.0.D+0 .AND. K.GT.2 ) THEN
                NBB = NB-IJFLD
                IF( IXP(NBB).GT.0 ) THEN
                 R = ((XLA(MB,NB)*RHOL(MB,NB)-XLA(2,NBB)*RHOL(2,NBB))/
     &             (XLA(MP,N)*RHOL(MP,N)-XLA(MB,NB)*RHOL(MB,NB)+SMALL))*
     &             ((DZGF(N)+DZGF(NB))/(DZGF(NB)+DZGF(NBB)))
                 THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,(2.D+0+R)/3.D+0))
                 DZF = DZGF(NB)/(DZGF(N)+DZGF(NB))
                 WLAB = WL(M,NPZ)*XLA(MB,NB)*RHOL(MB,NB)*
     &             (1.D+0-THETA*DZF)
     &             + WL(M,NPZ)*XLA(MP,N)*RHOL(MP,N)*THETA*DZF
                ELSE
                 WLAB = WL(M,NPZ)*XLA(MB,NB)*RHOL(MB,NB)
                ENDIF
              ELSEIF( WL(1,NPZ).GE.0.D+0 .AND. K.EQ.2 ) THEN
                WLAB = WL(M,NPZ)*XLA(MB,NB)*RHOL(MB,NB)
              ELSEIF( WL(1,NPZ).LT.0.D+0 .AND. K.LT.KFLD ) THEN
                NT = N+IJFLD
                IF( IXP(NT).GT.0 ) THEN
                 R = ((XLA(MP,N)*RHOL(MP,N)-XLA(2,NT)*RHOL(2,NT))/
     &             (XLA(MB,NB)*RHOL(MB,NB)-XLA(MP,N)*RHOL(MP,N)+SMALL))*
     &             ((DZGF(NB)+DZGF(N))/(DZGF(N)+DZGF(NT)))
                  THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,(2.D+0+R)/3.D+0))
                  DZF = DZGF(N)/(DZGF(N)+DZGF(NB))
                  WLAB = WL(M,NPZ)*XLA(MB,NB)*RHOL(MB,NB)*THETA*DZF
     &              + WL(M,NPZ)*XLA(MP,N)*RHOL(MP,N)*(1.D+0-THETA*DZF)
                ELSE
                  WLAB = WL(M,NPZ)*XLA(MP,N)*RHOL(MP,N)
                ENDIF
              ELSEIF( WL(1,NPZ).LT.0.D+0 .AND. K.EQ.KFLD ) THEN
                WLAB = WL(M,NPZ)*XLA(MP,N)*RHOL(MP,N)
              ENDIF
              FA(M,1) = FA(M,1) - AFZ(NPZ)*WLAB
!             IF( WG(1,NPZ).GE.0.D+0 .AND. K.GT.2 ) THEN
!               IF( IXP(NB-IJFLD).GT.0 ) THEN
!                NBB = NB-IJFLD
!                R = ((XGA(MB,NB)*RHOG(MB,NB)-XGA(2,NBB)*RHOG(2,NBB))/
!    &             (XGA(MP,N)*RHOG(MP,N)-XGA(MB,NB)*RHOG(MB,NB)+SMALL))*
!    &             ((DZGF(N)+DZGF(NB))/(DZGF(NB)+DZGF(NBB)))
!                THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,(2.D+0+R)/3.D+0))
!                DZF = DZGF(NB)/(DZGF(N)+DZGF(NB))
!                WGAB = WG(M,NPZ)*XGA(MB,NB)*RHOG(MB,NB)*
!    &             (1.D+0-THETA*DZF)
!    &             + WG(M,NPZ)*XGA(MP,N)*RHOG(MP,N)*THETA*DZF
!               ELSE
!                WGAB = WG(M,NPZ)*XGA(MB,NB)*RHOG(MB,NB)
!               ENDIF
!             ELSEIF( WG(1,NPZ).GE.0.D+0 .AND. K.EQ.2 ) THEN
!               WGAB = WG(M,NPZ)*XGA(MB,NB)*RHOG(MB,NB)
!             ELSEIF( WG(1,NPZ).LT.0.D+0 .AND. K.LT.KFLD ) THEN
!               NT = N+IJFLD
!               IF( IXP(NT).GT.0 ) THEN
!                R = ((XGA(MP,N)*RHOG(MP,N)-XGA(2,NT)*RHOG(2,NT))/
!    &             (XGA(MB,NB)*RHOG(MB,NB)-XGA(MP,N)*RHOG(MP,N)+SMALL))*
!    &             ((DZGF(NB)+DZGF(N))/(DZGF(N)+DZGF(NT)))
!                 THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,(2.D+0+R)/3.D+0))
!                 DZF = DZGF(N)/(DZGF(N)+DZGF(NB))
!                 WGAB = WG(M,NPZ)*XGA(MB,NB)*RHOG(MB,NB)*THETA*DZF
!    &              + WG(M,NPZ)*XGA(MP,N)*RHOG(MP,N)*(1.D+0-THETA*DZF)
!               ELSE
!                 WGAB = WG(M,NPZ)*XGA(MP,N)*RHOG(MP,N)
!               ENDIF
!             ELSEIF( WG(1,NPZ).LT.0.D+0 .AND. K.EQ.KFLD ) THEN
!               WGAB = WG(M,NPZ)*XGA(MP,N)*RHOG(MP,N)
!             ENDIF
              FGAB = XGA(MB,NB)*RHOG(MB,NB)
              FGAP = XGA(MP,N)*RHOG(MP,N)
              INDX = 3
              FGA = DIFMN( FGAB,FGAP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
              WGLAB = WG(M,NPZ)*FGA
              FA(M,1) = FA(M,1) - AFZ(NPZ)*WGAB
!
!---  Upwind air transport  ---
!
            ELSE
              FLAB = XLA(MB,NB)*RHOL(MB,NB)
              FLAP = XLA(MP,N)*RHOL(MP,N)
              INDX = 2
              FLA = DIFMN( FLAB,FLAP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
              FGAB = XGA(MB,NB)*RHOG(MB,NB)
              FGAP = XGA(MP,N)*RHOG(MP,N)
              INDX = 3
              FGA = DIFMN( FGAB,FGAP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
              FA(M,1) = -AFZ(NPZ)*(WL(M,NPZ)*FLA + WG(M,NPZ)*FGA -
     &          WTMA*WDGW(M,NPZ) + WTMA*WDLA(M,NPZ))
            ENDIF
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
            FA(M,2) = -AFY(NPY)*(VL(M,NPY)*FLA + VG(M,NPY)*FGA -
     &        WTMA*VDGW(M,NPY) + WTMA*VDLA(M,NPY))
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
            FA(M,3) = -AFX(NPX)*(UL(M,NPX)*FLA + UG(M,NPX)*FGA -
     &        WTMA*UDGW(M,NPX) + WTMA*UDLA(M,NPX))
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
            FA(M,4) = AFX(NQX)*(UL(MF,NQX)*FLA + UG(MF,NQX)*FGA -
     &        WTMA*UDGW(MF,NQX) + WTMA*UDLA(MF,NQX))
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
            FA(M,5) = AFY(NQY)*(VL(MF,NQY)*FLA + VG(MF,NQY)*FGA -
     &        WTMA*VDGW(MF,NQY) + WTMA*VDLA(MF,NQY))
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
!
!---  TVD air transport  ---
!
            IF( ISLC(10).EQ.1 ) THEN
              FA(M,6) = AFZ(NQZ)*WTMA*(WDLA(MF,NQZ)-WDGW(MF,NQZ))
              IF( WL(1,NQZ).GE.0.D+0 .AND. K.GT.1 ) THEN
                IF( IXP(NB).GT.0 ) THEN
                 R = ((XLA(MP,N)*RHOL(MP,N)-XLA(2,NB)*RHOL(2,NB))/
     &             (XLA(MT,NT)*RHOL(MT,NT)-XLA(MP,N)*RHOL(MP,N)+SMALL))*
     &             ((DZGF(NT)+DZGF(N))/(DZGF(N)+DZGF(NB)))
                 THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,(2.D+0+R)/3.D+0))
                 DZF = DZGF(N)/(DZGF(NT)+DZGF(N))
                 WLAB = WL(MF,NQZ)*XLA(MT,NT)*RHOL(MT,NT)*THETA*DZF
     &             + WL(MF,NQZ)*XLA(MP,N)*RHOL(MP,N)*(1.D+0-THETA*DZF)
                ELSE
                 WLAB = WL(MF,NQZ)*XLA(MP,N)*RHOL(MP,N)
                ENDIF
              ELSEIF( WL(1,NQZ).GE.0.D+0 .AND. K.EQ.1 ) THEN
                WLAB = WL(MF,NQZ)*XLA(MP,N)*RHOL(MP,N)
              ELSEIF( WL(1,NQZ).LT.0.D+0 .AND. K.LT.KFLD-1 ) THEN
                NTT = NT+IJFLD
                IF( IXP(NTT).GT.0 ) THEN
                 R = ((XLA(MT,NT)*RHOL(MT,NT)-XLA(2,NTT)*RHOL(2,NTT))/
     &             (XLA(MP,N)*RHOL(MP,N)-XLA(MT,NT)*RHOL(MT,NT)+SMALL))*
     &             ((DZGF(N)+DZGF(NT))/(DZGF(NT)+DZGF(NTT)))
                 THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,(2.D+0+R)/3.D+0))
                 DZF = DZGF(NT)/(DZGF(NT)+DZGF(N))
                 WLAB = WL(MF,NQZ)*XLA(MT,NT)*RHOL(MT,NT)*
     &             (1.D+0-THETA*DZF)
     &             + WL(MF,NQZ)*XLA(MP,N)*RHOL(MP,N)*THETA*DZF
                ELSE
                 WLAB = WL(MF,NQZ)*XLA(MT,NT)*RHOL(MT,NT)
                ENDIF
              ELSEIF( WL(1,NQZ).LT.0.D+0 .AND. K.EQ.KFLD-1 ) THEN
                WLAB = WL(MF,NQZ)*XLA(MT,NT)*RHOL(MT,NT)
              ENDIF
              FA(M,6) = FA(M,6) + AFZ(NQZ)*WLAB
!             IF( WG(1,NQZ).GE.0.D+0 .AND. K.GT.1 ) THEN
!               IF( IXP(NB).GT.0 ) THEN
!                R = ((XGA(MP,N)*RHOG(MP,N)-XGA(2,NB)*RHOG(2,NB))/
!    &             (XGA(MT,NT)*RHOG(MT,NT)-XGA(MP,N)*RHOG(MP,N)+SMALL))*
!    &             ((DZGF(NT)+DZGF(N))/(DZGF(N)+DZGF(NB)))
!                THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,(2.D+0+R)/3.D+0))
!                DZF = DZGF(N)/(DZGF(NT)+DZGF(N))
!                WGAB = WG(MF,NQZ)*XGA(MT,NT)*RHOG(MT,NT)*THETA*DZF
!    &             + WG(MF,NQZ)*XGA(MP,N)*RHOG(MP,N)*(1.D+0-THETA*DZF)
!               ELSE
!                WGAB = WG(MF,NQZ)*XGA(MP,N)*RHOG(MP,N)
!               ENDIF
!             ELSEIF( WG(1,NQZ).GE.0.D+0 .AND. K.EQ.1 ) THEN
!               WGAB = WG(MF,NQZ)*XGA(MP,N)*RHOG(MP,N)
!             ELSEIF( WG(1,NQZ).LT.0.D+0 .AND. K.LT.KFLD-1 ) THEN
!               NTT = NT+IJFLD
!               IF( IXP(NTT).GT.0 ) THEN
!                R = ((XGA(MT,NT)*RHOG(MT,NT)-XGA(2,NTT)*RHOG(2,NTT))/
!    &             (XGA(MP,N)*RHOG(MP,N)-XGA(MT,NT)*RHOG(MT,NT)+SMALL))*
!    &             ((DZGF(N)+DZGF(NT))/(DZGF(NT)+DZGF(NTT)))
!                THETA = MAX( ZERO,MIN(2.D+0,2.D+0*R,(2.D+0+R)/3.D+0))
!                DZF = DZGF(NT)/(DZGF(NT)+DZGF(N))
!                WGAB = WG(MF,NQZ)*XGA(MT,NT)*RHOG(MT,NT)*
!    &             (1.D+0-THETA*DZF)
!    &             + WG(MF,NQZ)*XGA(MP,N)*RHOG(MP,N)*THETA*DZF
!               ELSE
!                WGAB = WG(MF,NQZ)*XGA(MT,NT)*RHOG(MT,NT)
!               ENDIF
!             ELSEIF( WG(1,NQZ).LT.0.D+0 .AND. K.EQ.KFLD-1 ) THEN
!               WGAB = WG(MF,NQZ)*XGA(MT,NT)*RHOG(MT,NT)
!             ENDIF
              FGAT = XGA(MT,NT)*RHOG(MT,NT)
              FGAP = XGA(MP,N)*RHOG(MP,N)
              INDX = 3
              FGA = DIFMN( FGAP,FGAT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX )
              WGAB = WG(MF,NQZ)*FGA
              FA(M,6) = FA(M,6) + AFZ(NQZ)*WGAB
!
!---  Upwind air transport  ---
!
            ELSE
              FLAT = XLA(MT,NT)*RHOL(MT,NT)
              FLAP = XLA(MP,N)*RHOL(MP,N)
              INDX = 2
              FLA = DIFMN( FLAP,FLAT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX )
              FGAT = XGA(MT,NT)*RHOG(MT,NT)
              FGAP = XGA(MP,N)*RHOG(MP,N)
              INDX = 3
              FGA = DIFMN( FGAP,FGAT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX )
              FA(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLA + WG(MF,NQZ)*FGA -
     &         WTMA*WDGW(MF,NQZ) + WTMA*WDLA(MF,NQZ))
            ENDIF
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
!---  End of JCBALG group
!
      RETURN
      END

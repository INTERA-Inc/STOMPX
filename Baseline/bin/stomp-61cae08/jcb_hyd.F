!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBA37
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
!     Load the Jacobian matrix for the CO2 equation with
!     aqueous, gas, liquid CO2, and hydrate contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
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
      USE FLUXGC
      USE FDVP
      USE FDVH
      USE FDVGC
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
      SUB_LOG(ISUB_LOG) = '/JCBA37'
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
        USLX = PORD(1,N)*XLA(1,N)*RHOL(1,N)*SL(1,N)
        USGX = PORD(1,N)*XGA(1,N)*RHOG(1,N)*SG(1,N)
        USNX = PORD(1,N)*XNC(1,1,N)*RHON(1,N)*SN(1,N)
        USHX = PORD(1,N)*XHA(1,N)*RHOH(1,N)*SH(1,N)
        DO 100 M = 1,ISVC+1
          STA(M) = 0.D+0
          MP = M + 1
          STA(M) = STA(M) + 
     &      (PORD(MP,N)*XLA(MP,N)*RHOL(MP,N)*SL(MP,N)-USLX)
          STA(M) = STA(M) + 
     &      (PORD(MP,N)*XGA(MP,N)*RHOG(MP,N)*SG(MP,N)-USGX)
          STA(M) = STA(M) + 
     &      (PORD(MP,N)*XNC(1,MP,N)*RHON(MP,N)*SN(MP,N)-USNX)
          STA(M) = STA(M) + 
     &      (PORD(MP,N)*XHA(MP,N)*RHOH(MP,N)*SH(MP,N)-USHX)
          STA(M) = STA(M)*DTI*VOL(N)
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
            FA(M,1) = -AFZ(NPZ)*(WLA(M,NPZ)+WGA(M,NPZ)+WNC(1,M,NPZ))
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
            FA(M,2) = -AFY(NPY)*(VLA(M,NPY)+VGA(M,NPY)+VNC(1,M,NPY))
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
            FA(M,3) = -AFX(NPX)*(ULA(M,NPX)+UGA(M,NPX)+UNC(1,M,NPX))
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
            MF = MFLX(M)
            FA(M,4) = AFX(NQX)*(ULA(MF,NQX)+UGA(MF,NQX)+UNC(1,MF,NQX))
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
            MF = MFLX(M)
            FA(M,5) = AFY(NQY)*(VLA(MF,NQY)+VGA(MF,NQY)+VNC(1,MF,NQY))
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
            MF = MFLX(M)
            FA(M,6) = AFZ(NQZ)*(WLA(MF,NQZ)+WGA(MF,NQZ)+WNC(1,MF,NQZ))
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute CO2 equation residuals  ---
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
!---  End of JCBA37 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBO37
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
!     Load the Jacobian matrix for the CH4 equation with
!     aqueous, gas, liquid CO2, and hydrate contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
      USE FLUXD
      USE FDVP
      USE FDVH
      USE FDVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STO(LUK+1),ROP(LUK),ROA(LUK,6),FO(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBO37'
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
        USLX = PORD(1,N)*XLO(1,N)*RHOL(1,N)*SL(1,N)
        USGX = PORD(1,N)*XGO(1,N)*RHOG(1,N)*SG(1,N)
        USNX = PORD(1,N)*XNC(2,1,N)*RHON(1,N)*SN(1,N)
        USHX = PORD(1,N)*XHO(1,N)*RHOH(1,N)*SH(1,N)
        DO 100 M = 1,ISVC+1
          STO(M) = 0.D+0
          MP = M + 1
          STO(M) = STO(M) + 
     &      (PORD(MP,N)*XLO(MP,N)*RHOL(MP,N)*SL(MP,N)-USLX)
          STO(M) = STO(M) + 
     &      (PORD(MP,N)*XGO(MP,N)*RHOG(MP,N)*SG(MP,N)-USGX)
          STO(M) = STO(M) + 
     &      (PORD(MP,N)*XNC(2,MP,N)*RHON(MP,N)*SN(MP,N)-USNX)
          STO(M) = STO(M) + 
     &      (PORD(MP,N)*XHO(MP,N)*RHOH(MP,N)*SH(MP,N)-USHX)
          STO(M) = STO(M)*DTI*VOL(N)
  100   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 210 MD = 1,6
          DO 200 M = 1,ISVF
            FO(M,MD) = 0.D+0
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
            FO(M,1) = -AFZ(NPZ)*(WLO(M,NPZ)+WGO(M,NPZ)+WNC(2,M,NPZ))
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
            FO(M,2) = -AFY(NPY)*(VLO(M,NPY)+VGO(M,NPY)+VNC(2,M,NPY))
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
            FO(M,3) = -AFX(NPX)*(ULO(M,NPX)+UGO(M,NPX)+UNC(2,M,NPX))
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
            MF = MFLX(M)
            FO(M,4) = AFX(NQX)*(ULO(MF,NQX)+UGO(MF,NQX)+UNC(2,MF,NQX))
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
            MF = MFLX(M)
            FO(M,5) = AFY(NQY)*(VLO(MF,NQY)+VGO(MF,NQY)+VNC(2,MF,NQY))
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
            MF = MFLX(M)
            FO(M,6) = AFZ(NQZ)*(WLO(MF,NQZ)+WGO(MF,NQZ)+WNC(2,MF,NQZ))
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute CH4 equation residuals  ---
!
        ROS = STO(1) - SRCO(2,N)
        DO 900 MD = 1,6
          ROS = ROS + FO(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          ROP(M) = STO(M+1) - SRCO(M+2,N)
          MM = 2*M
          DO 910 MD = 1,6
            ROP(M) = ROP(M) + FO(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,6
            ROA(M,MD) = ROS - FO(1,MD) + FO(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Load Jacobian Matrix  ---
!
        CALL JCBBM( ROS,ROP,ROA,N,I,J,K,IEQDO )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBO37 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS37
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
!     Load the Jacobian matrix for the salt equation with
!     aqueous contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
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
      SUB_LOG(ISUB_LOG) = '/JCBS37'
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
        DO 10 M = 1,ISVC+1
          MP = M + 1
          STS(M) = (TMS(MP,N)-TMS(1,N))*VOL(N)*DTI
   10   CONTINUE
!
!---    Initialize surface fluxes  ---
!
        DO 30 MD = 1,6
          DO 20 M = 1,ISVF
            FS(M,MD) = 0.D+0
   20     CONTINUE
   30   CONTINUE
!
!---    Bottom surface flux component  ---
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
!---    South surface flux component ---
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
!---    West surface flux component ---
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
!---    East surface flux component  ---
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
!---    North surface flux component  ---
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
!---    Top surface flux component  ---
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
!---    Salt equation residuals  ---
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
!---  End of JCBS37 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT37
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
!     Load the Jacobian matrix for the energy equation with
!     aqueous, gas, liquid CO2, and hydrate contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVI
      USE FDVH
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STT(LUK+1),RTP(LUK),RTA(LUK,6),FT(LSFV,6)
      REAL*8 UEL(LSV),UEN(LSV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT37'
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
        IZN = IZ(N)
!
!---  First-order, forward-difference, time differential
!     assume rock and entrapped water density increases
!     porportionally with increasing porosity  ---
!
        UEL(1) = HL(1,N) - (PG(1,N)+PATM)/RHOL(1,N)
        UEN(1) = HN(1,N) - (PG(1,N)+PATM)/RHON(1,N)
        USTOX = (1.D+0-PORD(1,N))*T(1,N)
        USPOX = PORD(1,N)*SS(1,N)*RHOSP(1,N)*HSP(1,N)
        USLOX = PORD(1,N)*SL(1,N)*RHOL(1,N)*UEL(1)
        USGOX = PORD(1,N)*SG(1,N)*RHOG(1,N)*UEG(1,N)
        USNOX = PORD(1,N)*SN(1,N)*RHON(1,N)*UEN(1)
        USHOX = PORD(1,N)*SH(1,N)*RHOH(1,N)*HH(1,N)
        USIOX = PORD(1,N)*SI(1,N)*RHOI(1,N)*HI(1,N)
        DO 100 M = 1,ISVC+1
          STT(M) = 0.D+0
          MP = M + 1
          UEL(MP) = HL(MP,N) - (PG(MP,N)+PATM)/RHOL(MP,N)
          UEN(MP) = HN(MP,N) - (PG(MP,N)+PATM)/RHON(MP,N)
          USTX = (1.D+0-PORD(MP,N))*T(MP,N)
          USPX = PORD(MP,N)*SS(MP,N)*RHOSP(MP,N)*HSP(MP,N)
          USLX = PORD(MP,N)*SL(MP,N)*RHOL(MP,N)*UEL(MP)
          USGX = PORD(MP,N)*SG(MP,N)*RHOG(MP,N)*UEG(MP,N)
          USNX = PORD(MP,N)*SN(MP,N)*RHON(MP,N)*UEN(MP)
          USHX = PORD(MP,N)*SH(MP,N)*RHOH(MP,N)*HH(MP,N)
          USIX = PORD(MP,N)*SI(MP,N)*RHOI(MP,N)*HI(MP,N)
          STT(M) = STT(M) + RHOS(IZN)*CPS(IZN)*(USTX-USTOX)
          STT(M) = STT(M) + (USPX-USPOX)
          STT(M) = STT(M) + (USLX-USLOX)
          STT(M) = STT(M) + (USGX-USGOX)
          STT(M) = STT(M) + (USNX-USNOX)
          STT(M) = STT(M) + (USHX-USHOX)
          STT(M) = STT(M) + (USIX-USIOX)
          STT(M) = STT(M)*DTI*VOL(N)
  100   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 210 MD = 1,6
          DO 200 M = 1,ISVF
            FT(M,MD) = 0.D+0
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
            FT(M,1) = -AFZ(NPZ)*WQ(M,NPZ)
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
            FT(M,2) = -AFY(NPY)*VQ(M,NPY)
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
            FT(M,3) = -AFX(NPX)*UQ(M,NPX)
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
            MF = MFLX(M)
            FT(M,4) = AFX(NQX)*UQ(MF,NQX)
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
            MF = MFLX(M)
            FT(M,5) = AFY(NQY)*VQ(MF,NQY)
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
            MF = MFLX(M)
            FT(M,6) = AFZ(NQZ)*WQ(MF,NQZ)
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute air equation residuals  ---
!
        RTS = STT(1) - SRCT(2,N)
        DO 900 MD = 1,6
          RTS = RTS + FT(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RTP(M) = STT(M+1) - SRCT(M+2,N)
          MM = 2*M
          DO 910 MD = 1,6
            RTP(M) = RTP(M) + FT(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,6
            RTA(M,MD) = RTS - FT(1,MD) + FT(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Load Jacobian Matrix  ---
!
        CALL JCBBM( RTS,RTP,RTA,N,I,J,K,IEQT )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT37 group
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW37
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
!     Load the Jacobian matrix for the H2O equation with
!     aqueous, gas, liquid CO2, and hydrate contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
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
      USE FLUXGC
      USE FDVP
      USE FDVI
      USE FDVH
      USE FDVGC
      USE FDVA
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
      SUB_LOG(ISUB_LOG) = '/JCBW37'
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
        USLX = PORD(1,N)*XLW(1,N)*RHOL(1,N)*SL(1,N)
        USGX = PORD(1,N)*XGW(1,N)*RHOG(1,N)*SG(1,N)
        USNX = PORD(1,N)*XNC(3,1,N)*RHON(1,N)*SN(1,N)
        USHX = PORD(1,N)*XHW(1,N)*RHOH(1,N)*SH(1,N)
        USIX = PORD(1,N)*RHOI(1,N)*SI(1,N)
        DO 100 M = 1,ISVC+1
          STWX(M) = 0.D+0
          MP = M + 1
          STWX(M) = STWX(M) + 
     &      (PORD(MP,N)*XLW(MP,N)*RHOL(MP,N)*SL(MP,N)-USLX)
          STWX(M) = STWX(M) + 
     &      (PORD(MP,N)*XGW(MP,N)*RHOG(MP,N)*SG(MP,N)-USGX)
          STWX(M) = STWX(M) + 
     &      (PORD(MP,N)*XNW(MP,N)*RHON(MP,N)*SN(MP,N)-USNX)
          STWX(M) = STWX(M) + 
     &      (PORD(MP,N)*XHW(MP,N)*RHOH(MP,N)*SH(MP,N)-USHX)
          STWX(M) = STWX(M) + 
     &      (PORD(MP,N)*RHOI(MP,N)*SI(MP,N)-USIX)
          STWX(M) = STWX(M)*DTI*VOL(N)
  100   CONTINUE
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
            FW(M,1) = -AFZ(NPZ)*(WLW(M,NPZ)+WGW(M,NPZ)+WNC(3,M,NPZ))
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
            FW(M,2) = -AFY(NPY)*(VLW(M,NPY)+VGW(M,NPY)+VNC(3,M,NPY))
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
            FW(M,3) = -AFX(NPX)*(ULW(M,NPX)+UGW(M,NPX)+UNC(3,M,NPX))
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
            MF = MFLX(M)
            FW(M,4) = AFX(NQX)*(ULW(MF,NQX)+UGW(MF,NQX)+UNC(3,MF,NQX))
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
            MF = MFLX(M)
            FW(M,5) = AFY(NQY)*(VLW(MF,NQY)+VGW(MF,NQY)+VNC(3,MF,NQY))
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
            MF = MFLX(M)
            FW(M,6) = AFZ(NQZ)*(WLW(MF,NQZ)+WGW(MF,NQZ)+WNC(3,MF,NQZ))
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
!---  End of JCBW37 group
!
      RETURN
      END

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
!     (gas boundary, CO2 equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = -AFZ(NPZ)*WGA(MP,NPZ)
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
!     (gas boundary, CO2 equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = -AFY(NPY)*VGA(MP,NPY)
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
!     (gas boundary, CO2 equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = -AFX(NPX)*UGA(MP,NPX)
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
!     (gas boundary, CO2 equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = AFX(NQX)*UGA(MN,NQX)
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
!     (gas boundary, CO2 equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = AFY(NQY)*VGA(MN,NQY)
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
!     (gas boundary, CO2 equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = AFZ(NQZ)*WGA(MN,NQZ)
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
!     (gas boundary, CH4 equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = -AFZ(NPZ)*WGO(MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
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
!     (gas boundary, CH4 equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = -AFY(NPY)*VGO(MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
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
!     (gas boundary, CH4 equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = -AFX(NPX)*UGO(MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
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
!     (gas boundary, CH4 equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = AFX(NQX)*UGO(MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
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
!     (gas boundary, CH4 equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = AFY(NQY)*VGO(MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
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
!     (gas boundary, CH4 equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = AFZ(NQZ)*WGO(MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWB( N,NB,NPZ )
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
!     (gas boundary, H2O equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBGWB'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WGW(MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWS( N,NB,NPY )
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
!     (gas boundary, H2O equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBGWS'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VGW(MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWW( N,NB,NPX )
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
!     (gas boundary, H2O equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBGWW'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*UGW(MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWE( N,NB,NQX )
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
!     (gas boundary, H2O equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBGWE'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*UGW(MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWN( N,NB,NQY )
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
!     (gas boundary, H2O equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBGWN'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VGW(MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWT( N,NB,NQZ )
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
!     (gas boundary, H2O equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBGWT'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WGW(MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWT group  ---
!
      RETURN
      END

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
!     (aqueous boundary, CO2 equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = -AFZ(NPZ)*WLA(MP,NPZ)
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
!     (aqueous boundary, CO2 equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = -AFY(NPY)*VLA(MP,NPY)
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
!     (aqueous boundary, CO2 equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = -AFX(NPX)*ULA(MP,NPX)
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
!     (aqueous boundary, CO2 equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = AFX(NQX)*ULA(MN,NQX)
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
!     (aqueous boundary, CO2 equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = AFY(NQY)*VLA(MN,NQY)
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
!     (aqueous boundary, CO2 equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
        RS(M) = AFZ(NQZ)*WLA(MN,NQZ)
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
!     (aqueous boundary, CH4 equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = -AFZ(NPZ)*WLO(MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
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
!     (aqueous boundary, CH4 equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = -AFY(NPY)*VLO(MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
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
!     (aqueous boundary, CH4 equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = -AFX(NPX)*ULO(MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
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
!     (aqueous boundary, CH4 equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = AFX(NQX)*ULO(MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
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
!     (aqueous boundary, CH4 equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = AFY(NQY)*VLO(MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
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
!     (aqueous boundary, CH4 equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXD
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
        RS(M) = AFZ(NQZ)*WLO(MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWB( N,NB,NPZ )
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
!     (aqueous boundary, H2O equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBLWB'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WLW(MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWS( N,NB,NPY )
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
!     (aqueous boundary, H2O equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBLWS'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VLW(MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWW( N,NB,NPX )
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
!     (aqueous boundary, H2O equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBLWW'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*ULW(MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWE( N,NB,NQX )
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
!     (aqueous boundary, H2O equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBLWE'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*ULW(MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWN( N,NB,NQY )
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
!     (aqueous boundary, H2O equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBLWN'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VLW(MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWT( N,NB,NQZ )
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
!     (aqueous boundary, H2O equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
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
      SUB_LOG(ISUB_LOG) = '/JCBLWT'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WLW(MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNAB( N,NB,NPZ )
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
!     (liquid-CO2 boundary, CO2 equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNAB'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WNC(1,MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNAS( N,NB,NPY )
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
!     (liquid-CO2 boundary, CO2 equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNAS'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VNC(1,MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNAW( N,NB,NPX )
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
!     (liquid-CO2 boundary, CO2 equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNAW'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*UNC(1,MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNAE( N,NB,NQX )
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
!     (liquid-CO2 boundary, CO2 equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNAE'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*UNC(1,MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNAN( N,NB,NQY )
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
!     (liquid-CO2 boundary, CO2 equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNAN'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VNC(1,MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNAT( N,NB,NQZ )
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
!     (liquid-CO2 boundary, CO2 equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNAT'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WNC(1,MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQA )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNAT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNOB( N,NB,NPZ )
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
!     (liquid-CO2 boundary, CH4 equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNOB'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WNC(2,MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNOS( N,NB,NPY )
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
!     (liquid-CO2 boundary, CH4 equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNOS'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VNC(2,MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNOW( N,NB,NPX )
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
!     (liquid-CO2 boundary, CH4 equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNOW'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*UNC(2,MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNOW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNOE( N,NB,NQX )
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
!     (liquid-CO2 boundary, CH4 equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNOE'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*UNC(2,MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNON( N,NB,NQY )
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
!     (liquid-CO2 boundary, CH4 equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNON'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VNC(2,MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNOT( N,NB,NQZ )
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
!     (liquid-CO2 boundary, CH4 equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNOT'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WNC(2,MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQDO )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNWB( N,NB,NPZ )
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
!     (liquid-CO2 boundary, H2O equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNWB'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WNC(3,MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNWS( N,NB,NPY )
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
!     (liquid-CO2 boundary, H2O equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNWS'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VNC(3,MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNWW( N,NB,NPX )
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
!     (liquid-CO2 boundary, H2O equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNWW'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*UNC(3,MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNWW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNWE( N,NB,NQX )
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
!     (liquid-CO2 boundary, H2O equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNWE'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*UNC(3,MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNWN( N,NB,NQY )
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
!     (liquid-CO2 boundary, H2O equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNWN'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VNC(3,MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBNWT( N,NB,NQZ )
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
!     (liquid-CO2 boundary, H2O equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXGC
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
      SUB_LOG(ISUB_LOG) = '/JCBNWT'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WNC(3,MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBNWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSB( N,NB,NPZ )
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
!     (salt boundary, salt equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBSB'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WS(MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSS( N,NB,NPY )
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
!     (salt boundary, salt equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBSS'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VS(MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSW( N,NB,NPX )
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
!     (salt boundary, salt equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBSW'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*US(MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSE( N,NB,NQX )
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
!     (salt boundary, salt equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBSE'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*US(MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSN( N,NB,NQY )
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
!     (salt boundary, salt equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBSN'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VS(MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBST( N,NB,NQZ )
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
!     (salt boundary, salt equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBST'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WS(MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBST group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTB( N,NB,NPZ )
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
!     (thermal boundary, energy equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTB'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WQ(MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTS( N,NB,NPY )
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
!     (thermal boundary, energy equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTS'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VQ(MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTW( N,NB,NPX )
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
!     (thermal boundary, energy equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTW'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*UQ(MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTE( N,NB,NQX )
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
!     (thermal boundary, energy equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTE'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*UQ(MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTN( N,NB,NQY )
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
!     (thermal boundary, energy equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTN'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VQ(MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTT( N,NB,NQZ )
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
!     (thermal boundary, energy equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 September 2004
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTT'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WQ(MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTT group  ---
!
      RETURN
      END

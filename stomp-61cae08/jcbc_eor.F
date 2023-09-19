!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBC_EOR( IGC )
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
!     Load the Jacobian matrix for the petroleum component equation
!     with aqueous, gas, and nonaqueous-liquid contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 March 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
      USE FDVGC
      USE DUAL_POR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STCX(LUK+1),RCPX(LUK),RCAX(LUK,LSTC-1),FCX(LSFV,LSTC-1)
#ifdef dual_porosity
      REAL*8 RCTX(LUK)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBC_EOR'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 1000 N = 1,NFBN
!
!---    Inactive node or block refinement node  ---
!
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 1000
        I = ID(N)
        J = JD(N)
        K = KD(N)
!
!---    First-order, forward-difference, time differential  ---
!
        DO 10 M = 1,ISVC+1
          MP = M + 1
          STCX(M) = DTI*(TMC(IGC,MP,N)-TMC(IGC,1,N))
   10   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 30 MD = 1,LSTC-1
          DO 20 M = 1,ISVF
            FCX(M,MD) = 0.D+0
   20     CONTINUE
   30   CONTINUE
!
!---  Compute surface fluxes  ---
!
        MA = 1
!
!---    Bottom ---
!
        IF( K.EQ.1 .AND. N.LE.NFLD ) GOTO 110
!
!---    Loop over bottom connections  ---
!
        DO NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
!
!---      Single bottom connections  ---
!
          ELSE
            NPZ = NSZ(N)
          ENDIF
          IF( IGC.EQ.1 ) THEN
            DO 100 M = 1,ISVF
              FCX(M,MA) = -AFZ(NPZ)*(WLA(M,NPZ) + WGC(IGC,M,NPZ) + 
     &          WNC(IGC,M,NPZ))
  100       CONTINUE
          ELSE
            DO 102 M = 1,ISVF
              FCX(M,MA) = -AFZ(NPZ)*(WGC(IGC,M,NPZ) +  WNC(IGC,M,NPZ))
  102       CONTINUE
          ENDIF
          MA = MA + 1
        ENDDO
  110   CONTINUE
!
!---    South ---
!
        IF( J.EQ.1 .AND. N.LE.NFLD ) GOTO 210
!
!---    Loop over south connections  ---
!
        DO NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
!
!---      Single south connections  ---
!
          ELSE
            NPY = NSY(N)
          ENDIF
          IF( IGC.EQ.1 ) THEN
            DO 200 M = 1,ISVF
              FCX(M,MA) = -AFY(NPY)*(VLA(M,NPY) + VGC(IGC,M,NPY) + 
     &          VNC(IGC,M,NPY))
  200       CONTINUE
          ELSE
            DO 202 M = 1,ISVF
              FCX(M,MA) = -AFY(NPY)*(VGC(IGC,M,NPY) + VNC(IGC,M,NPY))
  202       CONTINUE
          ENDIF
          MA = MA + 1
        ENDDO
  210   CONTINUE
!
!---    West ---
!
        IF( I.EQ.1 .AND. N.LE.NFLD ) GOTO 310
!
!---    Loop over west connections  ---
!
        DO NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
!
!---      Single west connections  ---
!
          ELSE
            NPX = NSX(N)
          ENDIF
          IF( IGC.EQ.1 ) THEN
            DO 300 M = 1,ISVF
              FCX(M,MA) = -AFX(NPX)*(ULA(M,NPX) + UGC(IGC,M,NPX) + 
     &          UNC(IGC,M,NPX))
  300       CONTINUE
          ELSE
            DO 302 M = 1,ISVF
              FCX(M,MA) = -AFX(NPX)*(UGC(IGC,M,NPX) + UNC(IGC,M,NPX))
  302       CONTINUE
          ENDIF
          MA = MA + 1
        ENDDO
  310   CONTINUE
!
!---    East ---
!
        IF( I.EQ.IFLD .AND. N.LE.NFLD ) GOTO 410
!
!---    Loop over east connections  ---
!
        DO NC = 1,4
          NE = ICM(NC,4,N)
          IF( NE.EQ.0 ) EXIT
!
!---      Multiple west connections for east node  ---
!
          IF( INBS(3,NE).LT.0 ) THEN
            NQX = NSSX(N)
!
!---      Single west connection for east node  ---
!
          ELSE
            NQX = NSX(NE)
          ENDIF
          IF( IGC.EQ.1 ) THEN
            DO 400 M = 1,ISVF
              MF = MFLX(M)
              FCX(M,MA) = AFX(NQX)*(ULA(MF,NQX) + UGC(IGC,MF,NQX) + 
     &          UNC(IGC,MF,NQX))
  400       CONTINUE
          ELSE
            DO 402 M = 1,ISVF
              MF = MFLX(M)
              FCX(M,MA) = AFX(NQX)*(UGC(IGC,MF,NQX) + UNC(IGC,MF,NQX))
  402       CONTINUE
          ENDIF
          MA = MA + 1
        ENDDO
  410   CONTINUE
!
!---    North ---
!
        IF( J.EQ.JFLD .AND. N.LE.NFLD ) GOTO 510
!
!---    Loop over north connections  ---
!
        DO NC = 1,4
          NN = ICM(NC,5,N)
          IF( NN.EQ.0 ) EXIT
!
!---      Multiple south connections for north node  ---
!
          IF( INBS(2,NN).LT.0 ) THEN
            NQY = NSSY(N)
!
!---      Single south connection for north node  ---
!
          ELSE
            NQY = NSY(NN)
          ENDIF
          IF( IGC.EQ.1 ) THEN
            DO 500 M = 1,ISVF
              MF = MFLX(M)
              FCX(M,MA) = AFY(NQY)*(VLA(MF,NQY) + VGC(IGC,MF,NQY) + 
     &          VNC(IGC,MF,NQY))
  500       CONTINUE
          ELSE
            DO 502 M = 1,ISVF
              MF = MFLX(M)
              FCX(M,MA) = AFY(NQY)*(VGC(IGC,MF,NQY) + VNC(IGC,MF,NQY))
  502       CONTINUE
          ENDIF
          MA = MA + 1
        ENDDO
  510   CONTINUE
!
!---    Top ---
!
        IF( K.EQ.KFLD .AND. N.LE.NFLD ) GOTO 610
!
!---    Loop over top connections  ---
!
        DO NC = 1,4
          NT = ICM(NC,6,N)
          IF( NT.EQ.0 ) EXIT
!
!---      Multiple bottom connections for top node  ---
!
          IF( INBS(1,NT).LT.0 ) THEN
            NQZ = NSSZ(N)
!
!---      Single bottom connection for top node  ---
!
          ELSE
            NQZ = NSZ(NT)
          ENDIF
          IF( IGC.EQ.1 ) THEN
            DO 600 M = 1,ISVF
              MF = MFLX(M)
              FCX(M,MA) = AFZ(NQZ)*(WLA(MF,NQZ) + WGC(IGC,MF,NQZ) + 
     &          WNC(IGC,MF,NQZ))
  600       CONTINUE
          ELSE
            DO 602 M = 1,ISVF
              MF = MFLX(M)
              FCX(M,MA) = AFZ(NQZ)*(WGC(IGC,MF,NQZ) + WNC(IGC,MF,NQZ))
  602       CONTINUE
          ENDIF
          MA = MA + 1
        ENDDO
  610   CONTINUE
        MA = MA - 1
!
!---    Compute petroleum component equation residuals  ---
!
        RCSX = STCX(1) - SRCGC(IGC,2,N)
        DO 900 MD = 1,MA
          RCSX = RCSX + FCX(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RCPX(M) = STCX(M+1) - SRCGC(IGC,M+2,N)
          MM = 2*M
          DO 910 MD = 1,MA
            RCPX(M) = RCPX(M) + FCX(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,MA
            RCAX(M,MD) = RCSX - FCX(1,MD) + FCX(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBBM_BR( RCSX,RCPX,RCAX,N,IEQGC(IGC) )
#ifdef dual_porosity
!
!---    Fracture to matrix transfer contributions  ---
!
        RCSX = QGC_FM(IGC,1,N) + QNC_FM(IGC,1,N)
        DO M = 1,ISVC
          MM = 2*M
          RCPX(M) = QGC_FM(IGC,MM,N) + QNC_FM(IGC,MM,N)
          MM = 2*M + 1
          RCTX(M) = QGC_FM(IGC,MM,N) + QNC_FM(IGC,MM,N)
        ENDDO
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBBM_F( RCSX,RCPX,RCTX,N,IEQGC(IGC) )
#endif
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
       ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBC_EOR group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBCB_EOR( IGC,N,NB,NPZ )
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
!     (fluid boundary, petroleum component equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSX(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBCB_EOR'
      NBX = NB
      IF( IGC.EQ.1 ) THEN
        DO 100 M = 1,ISVC+1
          MP = MPOSB(M)
          RSX(M) = -AFZ(NPZ)*(WLA(MP,NPZ) + WGC(IGC,MP,NPZ) + 
     &      WNC(IGC,MP,NPZ))
  100   CONTINUE
      ELSE
        DO 102 M = 1,ISVC+1
          MP = MPOSB(M)
          RSX(M) = -AFZ(NPZ)*(WGC(IGC,MP,NPZ) + WNC(IGC,MP,NPZ))
  102   CONTINUE
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBCB_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBCS_EOR( IGC,N,NB,NPY )
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
!     (gas boundary, petroleum component equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 March 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSX(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBCS_EOR'
      NBX = NB
      IF( IGC.EQ.1 ) THEN
        DO 100 M = 1,ISVC+1
          MP = MPOSB(M)
          RSX(M) = -AFY(NPY)*(VLA(MP,NPY) + VGC(IGC,MP,NPY) + 
     &      VNC(IGC,MP,NPY))
  100   CONTINUE
      ELSE
        DO 102 M = 1,ISVC+1
          MP = MPOSB(M)
          RSX(M) = -AFY(NPY)*(VGC(IGC,MP,NPY) + VNC(IGC,MP,NPY))
  102   CONTINUE
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBCS_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBCW_EOR( IGC,N,NB,NPX )
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
!     (gas boundary, petroleum component equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 March 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSX(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBCW_EOR'
      NBX = NB
      IF( IGC.EQ.1 ) THEN
        DO 100 M = 1,ISVC+1
          MP = MPOSB(M)
          RSX(M) = -AFX(NPX)*(ULA(MP,NPX) + UGC(IGC,MP,NPX) + 
     &      UNC(IGC,MP,NPX))
  100   CONTINUE
      ELSE
        DO 102 M = 1,ISVC+1
          MP = MPOSB(M)
          RSX(M) = -AFX(NPX)*(UGC(IGC,MP,NPX) + UNC(IGC,MP,NPX))
  102   CONTINUE
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBCW_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBCE_EOR( IGC,N,NB,NQX )
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
!     (gas boundary, petroleum component equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 March 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSX(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBCE_EOR'
      NBX = NB
      IF( IGC.EQ.1 ) THEN
        DO 100 M = 1,ISVC+1
          MN = MNEGB(M)
          RSX(M) = AFX(NQX)*(ULA(MN,NQX) + UGC(IGC,MN,NQX) + 
     &      UNC(IGC,MN,NQX))
  100   CONTINUE
      ELSE
        DO 102 M = 1,ISVC+1
          MN = MNEGB(M)
          RSX(M) = AFX(NQX)*(UGC(IGC,MN,NQX) + UNC(IGC,MN,NQX))
  102   CONTINUE
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBCE_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBCN_EOR( IGC,N,NB,NQY )
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
!     (gas boundary, petroleum component equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 March 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSX(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBCN_EOR'
      NBX = NB
      IF( IGC.EQ.1 ) THEN
        DO 100 M = 1,ISVC+1
          MN = MNEGB(M)
          RSX(M) = AFY(NQY)*(VLA(MN,NQY) + VGC(IGC,MN,NQY) + 
     &      VNC(IGC,MN,NQY))
  100   CONTINUE
      ELSE
        DO 102 M = 1,ISVC+1
          MN = MNEGB(M)
          RSX(M) = AFY(NQY)*(VGC(IGC,MN,NQY) + VNC(IGC,MN,NQY))
  102   CONTINUE
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBCN_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBCT_EOR( IGC,N,NB,NQZ )
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
!     (gas boundary, petroleum component equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 March 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSX(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBCT_EOR'
      NBX = NB
      IF( IGC.EQ.1 ) THEN
        DO 100 M = 1,ISVC+1
          MN = MNEGB(M)
          RSX(M) = AFZ(NQZ)*(WLA(MN,NQZ) + WGC(IGC,MN,NQZ) + 
     &      WNC(IGC,MN,NQZ))
  100   CONTINUE
      ELSE
        DO 102 M = 1,ISVC+1
          MN = MNEGB(M)
          RSX(M) = AFZ(NQZ)*(WGC(IGC,MN,NQZ) + WNC(IGC,MN,NQZ))
  102   CONTINUE
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBCT_EOR group  ---
!
      RETURN
      END


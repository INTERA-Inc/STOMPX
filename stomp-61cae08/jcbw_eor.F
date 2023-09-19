!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_EOR
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
!     aqueous, gas, nonaqueous-liquid, ice, and hydrate contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 September 2010
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
      USE FDVP
      USE DUAL_POR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STWX(LUK+1),RWPX(LUK),RWAX(LUK,LSTC-1),FWX(LSFV,LSTC-1)
#ifdef dual_porosity
      REAL*8 RWTX(LUK)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBW_EOR'
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
        USLOX = PORD(1,N)*XLW(1,N)*RHOL(1,N)*SL(1,N)
        USGOX = PORD(1,N)*XGW(1,N)*RHOG(1,N)*SG(1,N)
        DO 10 M = 1,ISVC+1
          MP = M + 1
          USLX = PORD(MP,N)*XLW(MP,N)*RHOL(MP,N)*SL(MP,N)
          USGX = PORD(MP,N)*XGW(MP,N)*RHOG(MP,N)*SG(MP,N)
          DUSLX = (USLX-USLOX)
          DUSGX = (USGX-USGOX)
          STWX(M) = DUSLX + DUSGX
          STWX(M) = STWX(M)*DTI*VOL(N)
   10   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 30 MD = 1,LSTC-1
          DO 20 M = 1,ISVF
            FWX(M,MD) = 0.D+0
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
          DO 100 M = 1,ISVF
            FWX(M,MA) = -AFZ(NPZ)*(WLW(M,NPZ)+WGW(M,NPZ))
  100     CONTINUE
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
          DO 200 M = 1,ISVF
            FWX(M,MA) = -AFY(NPY)*(VLW(M,NPY)+VGW(M,NPY))
  200     CONTINUE
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
          DO 300 M = 1,ISVF
            FWX(M,MA) = -AFX(NPX)*(ULW(M,NPX)+UGW(M,NPX))
  300     CONTINUE
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
          DO 400 M = 1,ISVF
            MF = MFLX(M)
            FWX(M,MA) = AFX(NQX)*(ULW(MF,NQX)+UGW(MF,NQX))
  400     CONTINUE
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
          DO 500 M = 1,ISVF
            MF = MFLX(M)
            FWX(M,MA) = AFY(NQY)*(VLW(MF,NQY)+VGW(MF,NQY))
  500     CONTINUE
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
          DO 600 M = 1,ISVF
            MF = MFLX(M)
            FWX(M,MA) = AFZ(NQZ)*(WLW(MF,NQZ)+WGW(MF,NQZ))
  600     CONTINUE
          MA = MA + 1
        ENDDO
  610   CONTINUE
        MA = MA - 1
!
!---  Compute water equation residuals  ---
!
        RWSX = STWX(1) - SRCW(2,N)
        DO 900 MD = 1,MA
          RWSX = RWSX + FWX(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RWPX(M) = STWX(M+1) - SRCW(M+2,N)
          MM = 2*M
          DO 910 MD = 1,MA
            RWPX(M) = RWPX(M) + FWX(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,MA
            RWAX(M,MD) = RWSX - FWX(1,MD) + FWX(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBBM_BR( RWSX,RWPX,RWAX,N,IEQW )
#ifdef dual_porosity
!
!---    Fracture to matrix transfer contributions  ---
!
        RWSX = QLW_FM(1,N) + QGW_FM(1,N)
        DO M = 1,ISVC
          MM = 2*M
          RWPX(M) = QLW_FM(MM,N) + QGW_FM(MM,N)
          MM = 2*M + 1
          RWTX(M) = QLW_FM(MM,N) + QGW_FM(MM,N)
        ENDDO
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBBM_F( RWSX,RWPX,RWTX,N,IEQW )
#endif
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
       ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_EOR group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBWB_EOR( N,NB,NPZ )
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
!     (fluid flow boundary, water equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBWB_EOR'
      NBX = NB
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RSX(M) = -AFZ(NPZ)*(WGW(MP,NPZ) + WLW(MP,NPZ))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBWB_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBWS_EOR( N,NB,NPY )
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
!     (fluid flow boundary, water equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBWS_EOR'
      NBX = NB
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RSX(M) = -AFY(NPY)*(VGW(MP,NPY) + VLW(MP,NPY))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBWS_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBWW_EOR( N,NB,NPX )
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
!     (fluid flow boundary, water equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBWW_EOR'
      NBX = NB
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RSX(M) = -AFX(NPX)*(UGW(MP,NPX) + ULW(MP,NPX))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBWW_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBWE_EOR( N,NB,NQX )
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
!     (fluid flow boundary, water equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBWE_EOR'
      NBX = NB
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RSX(M) = AFX(NQX)*(UGW(MN,NQX) + ULW(MN,NQX))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBWE_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBWN_EOR( N,NB,NQY )
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
!     (fluid flow boundary, water equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBWN_EOR'
      NBX = NB
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RSX(M) = AFY(NQY)*(VGW(MN,NQY) + VLW(MN,NQY))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBWN_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBWT_EOR( N,NB,NQZ )
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
!     (fluid flow boundary, water equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBWT_EOR'
      NBX = NB
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RSX(M) = AFZ(NQZ)*(WGW(MN,NQZ) + WLW(MN,NQZ))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBWT_EOR group  ---
!
      RETURN
      END



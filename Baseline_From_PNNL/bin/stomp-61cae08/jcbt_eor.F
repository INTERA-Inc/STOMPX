!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_EOR
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
!     Written by MD White, PNNL, 14 September 2010
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
      USE GLB_PAR
      USE FLUXT
      USE FDVT
      USE FDVS
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
      REAL*8 STTX(LUK+1),RTPX(LUK),RTAX(LUK,LSTC-1),FTX(LSFV,LSTC-1)
#ifdef dual_porosity
      REAL*8 RTTX(LUK)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT_EOR'
      DO 1000 N = 1,NFBN
!
!---    Inactive node or block refinement node  ---
!
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 1000
        I = ID(N)
        J = JD(N)
        K = KD(N)
        IZN = IZ(N)
!
!---  First-order, forward-difference, time differential  ---
!
        USTOX = RHOS(IZN)*CPS(IZN)*(1.D+0-PORD(1,N))*T(1,N)
        USLOX = PORD(1,N)*SL(1,N)*RHOL(1,N)*UEL(1,N)
        USGOX = PORD(1,N)*SG(1,N)*RHOG(1,N)*UEG(1,N)
        USNOX = PORD(1,N)*SN(1,N)*RHON(1,N)*UEN(1,N)
        USSOX = PORD(1,N)*SS(1,N)*RHOSP(1,N)*HSP(1,N)
        DO 10 M = 1,ISVC+1
          MP = M + 1
          USTX = RHOS(IZN)*CPS(IZN)*(1.D+0-PORD(MP,N))*T(MP,N)
          USLX = PORD(MP,N)*SL(MP,N)*RHOL(MP,N)*UEL(MP,N)
          USGX = PORD(MP,N)*SG(MP,N)*RHOG(MP,N)*UEG(MP,N)
          USNX = PORD(MP,N)*SN(MP,N)*RHON(MP,N)*UEN(MP,N)
          USSX = PORD(MP,N)*SS(MP,N)*RHOSP(MP,N)*HSP(MP,N)
          DUSTX = (USTX-USTOX)
          DUSLX = (USLX-USLOX)
          DUSGX = (USGX-USGOX)
          DUSNX = (USNX-USNOX)
          DUSSX = (USSX-USSOX)
          STTX(M) = DUSTX + DUSLX + DUSGX + DUSNX + DUSSX
          STTX(M) = STTX(M)*DTI*VOL(N)
   10   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 30 MD = 1,LSTC-1
          DO 20 M = 1,ISVF
            FTX(M,MD) = 0.D+0
   20     CONTINUE
   30  CONTINUE
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
            FTX(M,MA) = -AFZ(NPZ)*WQ(M,NPZ)
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
            FTX(M,MA) = -AFY(NPY)*VQ(M,NPY)
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
            FTX(M,MA) = -AFX(NPX)*UQ(M,NPX)
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
            FTX(M,MA) = AFX(NQX)*UQ(MF,NQX)
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
            FTX(M,MA) = AFY(NQY)*VQ(MF,NQY)
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
            FTX(M,MA) = AFZ(NQZ)*WQ(MF,NQZ)
  600     CONTINUE
          MA = MA + 1
        ENDDO
  610   CONTINUE
        MA = MA - 1
!
!---  Compute air equation residuals  ---
!
        RTSX = STTX(1) - SRCT(2,N)
        DO 900 MD = 1,MA
          RTSX = RTSX + FTX(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RTPX(M) = STTX(M+1) - SRCT(M+2,N)
          MM = 2*M
          DO 910 MD = 1,MA
            RTPX(M) = RTPX(M) + FTX(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,MA
            RTAX(M,MD) = RTSX - FTX(1,MD) + FTX(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Load Jacobian Matrix  ---
!
        CALL JCBBM_BR( RTSX,RTPX,RTAX,N,IEQT )
#ifdef dual_porosity
!
!---    Fracture to matrix transfer contributions  ---
!
        RTSX = QQ_FM(1,N)
        DO M = 1,ISVC
          MM = 2*M
          RTPX(M) = QQ_FM(MM,N)
          MM = 2*M + 1
          RTTX(M) = QQ_FM(MM,N)
        ENDDO
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBBM_F( RTSX,RTPX,RTTX,N,IEQT )
#endif
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
       ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_EOR group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTB_EOR( N,NB,NPZ )
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
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTB_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RSX(M) = -AFZ(NPZ)*WQ(MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTB_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTS_EOR( N,NB,NPY )
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
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTS_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RSX(M) = -AFY(NPY)*VQ(MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTS_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTW_EOR( N,NB,NPX )
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
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTW_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RSX(M) = -AFX(NPX)*UQ(MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTW_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTE_EOR( N,NB,NQX )
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
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTE_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RSX(M) = AFX(NQX)*UQ(MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTE_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTN_EOR( N,NB,NQY )
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
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTN_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RSX(M) = AFY(NQY)*VQ(MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTN_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTT_EOR( N,NB,NQZ )
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
      USE FLUXT
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
      SUB_LOG(ISUB_LOG) = '/JCBTT_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RSX(M) = AFZ(NQZ)*WQ(MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTT_EOR group  ---
!
      RETURN
      END

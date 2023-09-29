!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_EOR
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
      USE FLUXS
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
      REAL*8 STSX(LUK+1),RSPX(LUK),RSAX(LUK,LSTC-1),FSX(LSFV,LSTC-1)
#ifdef dual_porosity
      REAL*8 RSTX(LUK)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBS_EOR'
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
          STSX(M) = (TMS(MP,N)-TMS(1,N))*VOL(N)*DTI
   10   CONTINUE
!
!---    Initialize surface fluxes  ---
!
        DO 30 MD = 1,LSTC-1
          DO 20 M = 1,ISVF
            FSX(M,MD) = 0.D+0
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
            FSX(M,MA) = -AFZ(NPZ)*WS(M,NPZ)
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
            FSX(M,MA) = -AFY(NPY)*VS(M,NPY)
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
            FSX(M,MA) = -AFX(NPX)*US(M,NPX)
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
            FSX(M,MA) = AFX(NQX)*US(MF,NQX)
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
            FSX(M,MA) = AFY(NQY)*VS(MF,NQY)
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
            FSX(M,MA) = AFZ(NQZ)*WS(MF,NQZ)
  600     CONTINUE
          MA = MA + 1
        ENDDO
  610   CONTINUE
        MA = MA - 1
!
!---    Salt equation residuals  ---
!
        RSSX = STSX(1) - SRCS(2,N)
        DO 900 MD = 1,MA
          RSSX = RSSX + FSX(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RSPX(M) = STSX(M+1) - SRCS(M+2,N)
          MM = 2*M
          DO 910 MD = 1,MA
            RSPX(M) = RSPX(M) + FSX(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,MA
            RSAX(M,MD) = RSSX - FSX(1,MD) + FSX(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Jacobian matrix loader  --
!
        CALL JCBBM_BR( RSSX,RSPX,RSAX,N,IEQS )
#ifdef dual_porosity
!
!---    Fracture to matrix transfer contributions  ---
!
        RSSX = QS_FM(1,N)
        DO M = 1,ISVC
          MM = 2*M
          RSPX(M) = QS_FM(MM,N)
          MM = 2*M + 1
          RSTX(M) = QS_FM(MM,N)
        ENDDO
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBBM_F( RSSX,RSPX,RSTX,N,IEQS )
#endif
 1000 CONTINUE
       ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSB_EOR( N,NB,NPZ )
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
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBSB_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RSX(M) = -AFZ(NPZ)*WS(MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSB_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSS_EOR( N,NB,NPY )
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
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBSS_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RSX(M) = -AFY(NPY)*VS(MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSS_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSW_EOR( N,NB,NPX )
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
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBSW_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RSX(M) = -AFX(NPX)*US(MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSW_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSE_EOR( N,NB,NQX )
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
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBSE_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RSX(M) = AFX(NQX)*US(MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSE_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSN_EOR( N,NB,NQY )
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
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBSN_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RSX(M) = AFY(NQY)*VS(MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSN_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBST_EOR( N,NB,NQZ )
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
      USE FLUXS
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
      SUB_LOG(ISUB_LOG) = '/JCBST_EOR'
      NBX = NB
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RSX(M) = AFZ(NQZ)*WS(MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RSX,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBST_EOR group  ---
!
      RETURN
      END


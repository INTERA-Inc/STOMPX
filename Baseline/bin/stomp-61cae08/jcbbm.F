!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBBM( RSS,RSP,RSA,N,I,J,K,MEQ )
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
!     Load the Jacobian matrix (banded matrix solver).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on September 5, 1996.
#ifdef sac
!
!     SAC Modification by WE Nichols, PNNL, June 12, 2003.
#endif
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_CL
      USE SOLTN
      USE JACOB
      USE GRID
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP(LUK),RSA(LUK,6)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBBM'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Node  ---
!
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        DO M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSP(M)-RSS)/DNR(M,N)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 210
          NMD = IXP(NB)
          DO M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,1)-RSS)/DNR(M,NB)
          ENDDO
        ENDIF
  210 CONTINUE
!
!---    South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 310
          NMD = IXP(NS)
          DO M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,2)-RSS)/DNR(M,NS)
          ENDDO
        ENDIF
  310   CONTINUE
!
!---    West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 410
          NMD = IXP(NW)
          DO M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,3)-RSS)/DNR(M,NW)
          ENDDO
        ENDIF
  410   CONTINUE
!
!---    East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 510
          NMD = IXP(NE)
          DO M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,4)-RSS)/DNR(M,NE)
          ENDDO
        ENDIF
  510   CONTINUE
!
!---    North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 610
          NMD = IXP(NN)
          DO M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,5)-RSS)/DNR(M,NN)
          ENDDO
        ENDIF
  610   CONTINUE
!
!---    Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 710
          NMD = IXP(NT)
          DO M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,6)-RSS)/DNR(M,NT)
          ENDDO
        ENDIF
  710   CONTINUE
!
!---  SPLib or Lis or PETSc solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 .OR. ILES.EQ.5 ) THEN
!
!---    Node  ---
!
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        MA = 0
!
!---    Coupled well  ---
!
        IF( IXW(N).NE.0 ) MA = MA + ISVC
        DO M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR(M,N)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
        MA = MA + ISVC
!
!---    Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 2210
          NMD = IXP(NB)
          DO M = 1,ISVC
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,1)-RSS)/DNR(M,NB)
          ENDDO
          MA = MA + ISVC
        ENDIF
 2210   CONTINUE
!
!---    South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 2310
          NMD = IXP(NS)
!
!---      Coupled well  ---
!
          IF( IXW(N).NE.0 ) MA = MA + ISVC
          DO M = 1,ISVC
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,2)-RSS)/DNR(M,NS)
          ENDDO
          MA = MA + ISVC
        ENDIF
 2310   CONTINUE
!
!---    West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 2410
          NMD = IXP(NW)
!
!---      Coupled well  ---
!
          IF( IXW(N).NE.0 ) MA = MA + ISVC
          DO M = 1,ISVC
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,3)-RSS)/DNR(M,NW)
          ENDDO
          MA = MA + ISVC
        ENDIF
 2410   CONTINUE
!
!---    East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 2510
          NMD = IXP(NE)
!
!---      Coupled well  ---
!
          IF( IXW(N).NE.0 ) MA = MA + 2*ISVC
          DO M = 1,ISVC
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,4)-RSS)/DNR(M,NE)
          ENDDO
          MA = MA + ISVC
        ENDIF
 2510   CONTINUE
!
!---    North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 2610
          NMD = IXP(NN)
!
!---      Coupled well  ---
!
          IF( IXW(N).NE.0 ) MA = MA + 2*ISVC
          DO M = 1,ISVC
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,5)-RSS)/DNR(M,NN)
          ENDDO
          MA = MA + ISVC
        ENDIF
 2610   CONTINUE
!
!---    Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 2710
          NMD = IXP(NT)
!
!---      Coupled well  ---
!
          IF( IXW(N).NE.0 ) MA = MA + 2*ISVC
          DO M = 1,ISVC
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,6)-RSS)/DNR(M,NT)
          ENDDO
          MA = MA + ISVC
        ENDIF
 2710   CONTINUE
      ELSE
#ifdef sac
        OPEN(UNIT=98,FILE='stomp.run')
        CLOSE(UNIT=98,STATUS='delete')
        OPEN(UNIT=99,FILE='stomp.fail')
        CLOSE(UNIT=99)
#endif
        INDX = 3
        CHMSGX = 'Unknown Linear Equation Solver'
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'JCBBM'
        RLMSGX = 0.D+0
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBBM group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBBM_BR( RSS,RSP,RSA,N,MEQ )
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
!     Load the Jacobian matrix for the block refinement configuration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 April 2016.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP(LUK),RSA(LUK,LSTC-1)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBBM_BR'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Node  ---
!
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        DO M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSP(M)-RSS)/DNR(M,N)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Loop over surface connections  ---
!
        MA = 1
        DO ISX = 1,6
          DO NCX = 1,4
            NA = ICM(NCX,ISX,N)
            IF( NA.EQ.0 ) EXIT
             NMD = IXP(NA)
             DO M = 1,ISVC
              DNRX = DNR(M,NA)
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,MA)-RSS)/DNRX
            ENDDO
            MA = MA + 1
          ENDDO
        ENDDO
!
!---  SPLib or Lis or PETSc solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 .OR. ILES.EQ.5 ) THEN
!
!---    Node  ---
!
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        MA = 0
        DO M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR(M,N)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Loop over surface connections  ---
!
        MA = 1
        DO ISX = 1,6
          DO NCX = 1,4
            NA = ICM(NCX,ISX,N)
            IF( NA.EQ.0 ) EXIT
            DO M = 1,ISVC
              DNRX = DNR(M,NA)
              MCOL = KLU(MP,M+(MA*ISVC))
              DLU(MCOL) = DLU(MCOL) + (RSA(M,MA)-RSS)/DNRX
            ENDDO
            MA = MA + 1
          ENDDO
        ENDDO
      ELSE
        INDX = 3
        CHMSGX = 'Unknown Linear Equation Solver'
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'JCBBM_BR'
        RLMSGX = 0.D+0
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBBM_BR group
!
      RETURN
      END


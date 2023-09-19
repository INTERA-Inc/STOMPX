!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BATEMAN( CX,SLX,ICX,NPATH,NSOLC,NSOLP,INDX )
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
!     Updates solute concentrations for chain decay via a matrix
!     solution to the Batemann equations
!
!    Yuan, D. and W. Kernan. 2007. "Explicit solutions for exit-only
!    radioactive decay chains." J. Appl. Phys., 101, 094907.  ---
!
!----------------------Authors-----------------------------------------!
!
!    Written by Mark D. White, PNNL, 20 July 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(NSOLC,NSOLC) :: CMX,XPF
      REAL*8, DIMENSION(NSOLC) :: CBX,CBNX,LMDX
      REAL*8, DIMENSION(NSOLC) :: CX
      REAL*8, DIMENSION(NSOLC) :: FRACP
      REAL*8 :: LAMBDAX
      INTEGER, DIMENSION(NSOLC) :: ICX
      INTEGER, DIMENSION(LCDS,LCDP) :: NSOLP
!
!----------------------Data Statements---------------------------------!
!
!      DATA NSOLP /1,2,3,5,6,7,0,1,2,3,4,5,6,7,1,2,3,5,7,0,0,
!     &  1,2,3,4,5,7,0,1,2,4,5,6,7,0,1,2,4,5,7,0,0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BATEMAN'
!
!---  Radioactive decay for solutes outside of decay chain  ---
!
      IF( INDX.EQ.0 ) THEN
        DO I = 1,NSOLC
          ISOLUX = ICX(I)
!
!---      Schroth half-life modification model  ---
!
          IF( PCLN(3,ISOLUX).LT.0.D+0 ) THEN
            IF( SLX.LT.PCLN(1,ISOLUX) ) THEN
              HLFX = 1.D+12
            ELSE
              HLFX = HLF(ISOLUX)*(PCLN(2,ISOLUX) + SLX - 
     &          2.D+0*PCLN(1,ISOLUX))/(SLX - PCLN(1,ISOLUX) + SMALL)
            ENDIF
          ELSE
            HLFX = HLF(ISOLUX)
          ENDIF
          LAMBDAX = LOG(2.D+0)/HLFX
          CX(I) = CX(I)*EXP(-LAMBDAX*DT)
          IF( CX(I).LT.1.D-20 ) CX(I) = 0.D+0
        ENDDO
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      DO I = 1,NSOLC
        ISOLUX = ICX(I)
        DO J = 1,NSOLC
          JSOLUX = ICX(J)
          CMX(I,J) = 0.D+0
          XPF(I,J) = CHDF(ISOLUX,JSOLUX)
        ENDDO
        CBX(I) = CX(I)
!
!---    Schroth half-life modification model  ---
!
        IF( PCLN(3,ISOLUX).LT.0.D+0 ) THEN
          IF( SLX.LT.PCLN(1,ISOLUX) ) THEN
            HLFX = 1.D+12
          ELSE
            HLFX = HLF(ISOLUX)*(PCLN(2,ISOLUX) + SLX - 
     &        2.D+0*PCLN(1,ISOLUX))/(SLX - PCLN(1,ISOLUX) + SMALL)
          ENDIF
        ELSE
          HLFX = HLF(ISOLUX)
        ENDIF
        LMDX(I) = LOG(2.D+0)/HLFX
      ENDDO
!!
!!     Compute C matrix using recursive formula of Yuan and Kernan
!!
!      CMX(1,1) = CBX(1)
!      DO I = 2,NSOLC
!        DO J = 1,I
!          IF( I.EQ.J ) THEN
!            CMX(I,J) = CBX(I)
!            DO K = 1,I-1
!              CMX(I,J) = CMX(I,J) - CMX(I,K)
!            ENDDO
!          ELSE
!            CMX(I,J) = 0.D+0
!            DO K = J,I-1
!              CMX(I,J) = CMX(I,J) + XPF(K,I)*LMDX(K)*CMX(K,J)/
!     &          (LMDX(I)-LMDX(J))
!            ENDDO
!          ENDIF
!        ENDDO
!      ENDDO
!!
!!     Compute new time level radionuclide amounts
!!
!      DO I = 1,NSOLC
!        CBNX(I) = 0.D+0
!        DO J = 1,I
!          CBNX(I) = CBNX(I) + CMX(I,J)*EXP(-LMDX(J)*DT)
!        ENDDO
!        IF( CBNX(I).LT.1.D-20 ) CBNX(I) = 0.D+0
!        CX(I) = CBNX(I)
!      ENDDO
!
!     Wikipedia Bateman equation
!
      DO N = 1,NSOLC
        CBNX(N) = 0.D+0
        FRACP(N) = 0.D+0
      ENDDO
      DO NP = 1,NPATH
      DO N = 1,NSOLC
        NX = NSOLP(N,NP)
        IF( NX.EQ.0 ) CYCLE
        VARX = 0.D+0
        DO I = 1,N
          IX = NSOLP(I,NP)
          IF( IX.EQ.0 ) CYCLE
          VAR1X = 1.D+0
          DO J = I,N-1
            JX = NSOLP(J,NP)
            IF( JX.EQ.0 ) CYCLE
            KX = NSOLP(J+1,NP)
            VAR1X = VAR1X*XPF(JX,KX)*LMDX(JX)
          ENDDO
          VAR2X = 0.D+0
          DO J = I,N
            JX = NSOLP(J,NP)
            IF( JX.EQ.0 ) CYCLE
            VAR3X = 1.D+0
            DO K = I,N
              KX = NSOLP(K,NP)
              IF( KX.EQ.0 ) CYCLE
              IF( KX.NE.JX ) THEN
                VAR3X = VAR3X*(LMDX(KX)-LMDX(JX))
              ENDIF
            ENDDO
            VAR2X = VAR2X + EXP(-LMDX(JX)*DT)/VAR3X
          ENDDO
          VARX = VARX + CX(IX)*VAR1X*VAR2X
        ENDDO
        CBNX(NX) = CBNX(NX) + VARX
        FRACP(NX) = FRACP(NX) + 1.D+0
      ENDDO
      ENDDO
      DO N = 1,NSOLC
        CBNX(N) = CBNX(N)/FRACP(N)
        IF( CBNX(N).LT.1.D-20 ) CBNX(N) = 0.D+0
      ENDDO
      DO N = 1,NSOLC
        CX(N) = CBNX(N)
      ENDDO
!
!---  End of BATEMAN group  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHAIN_DECAY
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
!     Alter solute concentrations for chain decay.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 21 July 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FDVP_FRC
      USE FDVP_BH
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
      REAL*8, DIMENSION(LCDS) :: CX
      REAL*8 :: LAMBDAX
      INTEGER, DIMENSION(LCDS) :: ICX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHAIN_DECAY'
!
!---  Loop over matrix grid cells  ---
!
      DO N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
!
!---    Zero chain-decay series  ---
!
        IF( NBCDS.EQ.0 ) THEN
          DO NSL = 1,NSOLU
!
!---        Schroth half-life modification model  ---
!
            IF( PCLN(3,NSL).LT.0.D+0 ) THEN
              IF( SL(2,N).LT.PCLN(1,NSL) ) THEN
                HLFX = 1.D+12
              ELSE
                HLFX = HLF(NSL)*(PCLN(2,NSL) + SL(2,N) - 
     &            2.D+0*PCLN(1,NSL))/(SL(2,N) - PCLN(1,NSL) + SMALL)
              ENDIF
            ELSE
              HLFX = HLF(NSL)
            ENDIF
            LAMBDAX = LOG(2.D+0)/HLFX
            C(N,NSL) = C(N,NSL)*EXP(-LAMBDAX*DT)
            IF( C(N,NSL).LT.1.D-20 ) C(N,NSL) = 0.D+0
          ENDDO
!
!---    Non-zero number of chain-decay series  ---
!
        ELSE
!
!---    Loop over the number of chain-decay series  ---
!
          MC = 0
          DO NCS = 1,NBCDS
            MC = MC + 1
            MC0 = MC
            NSOLC = IBCDS(MC)
            DO NSL = 1,NSOLC
              MC = MC + 1
              ICX(NSL) = IBCDS(MC)
              CX(NSL) = C(N,IBCDS(MC))
            ENDDO
!
!---        Radioactive chain decay  ---
!
            INDX = 1
!
!---        Radioactive decay of individual solutes  ---
!
            IF( NCS.EQ.NBCDS ) INDX = 0
            CALL BATEMAN( CX,SL(2,N),ICX,NBCDP(NCS),NSOLC,
     &        IBCDP(1,1,NCS),INDX )
            MC = MC0
            DO NSL = 1,NSOLC
              MC = MC + 1
              C(N,IBCDS(MC)) = CX(NSL)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!
!---  Fracture flow and transport  ---
!
      IF( ISLC(17).NE.0 .AND. (ISLC(74).EQ.1 .OR. ISLC(74).EQ.3) ) THEN
!
!---    Loop over fractures  ---
!
        DO NFX = 1,NF_FRC
!
!---      Loop over fracture triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---        Zero chain-decay series  ---
!
            IF( NBCDS.EQ.0 ) THEN
              DO NSL = 1,NSOLU
!
!---            Schroth half-life modification model  ---
!
                IF( PCLN(3,NSL).LT.0.D+0 ) THEN
                  IF( SL_FRC(2,NTX).LT.PCLN(1,NSL) ) THEN
                    HLFX = 1.D+12
                  ELSE
                    HLFX = HLF(NSL)*(PCLN(2,NSL) + SL_FRC(2,NTX) - 
     &                2.D+0*PCLN(1,NSL))/(SL_FRC(2,NTX) - PCLN(1,NSL)
     &                + SMALL)
                  ENDIF
                ELSE
                  HLFX = HLF(NSL)
                ENDIF
                LAMBDAX = LOG(2.D+0)/HLFX
                C_FRC(NTX,NSL) = C_FRC(NTX,NSL)*EXP(-LAMBDAX*DT)
                IF( C_FRC(NTX,NSL).LT.1.D-20 ) C_FRC(NTX,NSL) = 0.D+0
              ENDDO
            ELSE
!
!---          Loop over the number of chain-decay series  ---
!
              MC = 0
              DO NCS = 1,NBCDS
                MC = MC + 1
                MC0 = MC
                NSOLC = IBCDS(MC)
                DO NSL = 1,NSOLC
                  MC = MC + 1
                  ICX(NSL) = IBCDS(MC)
                  CX(NSL) = C_FRC(NTX,IBCDS(MC))
                ENDDO
!
!---            Radioactive chain decay  ---
!
                INDX = 1
!
!---            Radioactive decay of individual solutes  ---
!
                IF( NCS.EQ.NBCDS ) INDX = 0
                CALL BATEMAN( CX,SL_FRC(2,NTX),ICX,NBCDP(NCS),NSOLC,
     &            IBCDP(1,1,NCS),INDX )
                MC = MC0
                DO NSL = 1,NSOLC
                  MC = MC + 1
                  C_FRC(NTX,IBCDS(MC)) = CX(NSL)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
!---  Borehole flow and transport  ---
!
      IF( ISLC(17).NE.0 .AND. (ISLC(74).EQ.2 .OR. ISLC(74).EQ.3) ) THEN
!
!---    Loop over boreholes  ---
!
        DO NBH = 1,N_BH
!
!---      Loop over borehole nodes  ---
!
          DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---        Zero chain-decay series  ---
!
            IF( NBCDS.EQ.0 ) THEN
              DO NSL = 1,NSOLU
!
!---            Schroth half-life modification model  ---
!
                IF( PCLN(3,NSL).LT.0.D+0 ) THEN
                  IF( SL_BH(2,NBN).LT.PCLN(1,NSL) ) THEN
                    HLFX = 1.D+12
                  ELSE
                    HLFX = HLF(NSL)*(PCLN(2,NSL) + SL_BH(2,NBN) - 
     &                2.D+0*PCLN(1,NSL))/(SL_BH(2,NBN) - PCLN(1,NSL) + 
     &                SMALL)
                  ENDIF
                ELSE
                  HLFX = HLF(NSL)
                ENDIF
                LAMBDAX = LOG(2.D+0)/HLFX
                C_BH(NBN,NSL) = C_BH(NBN,NSL)*EXP(-LAMBDAX*DT)
                IF( C_BH(NBN,NSL).LT.1.D-20 ) C_BH(NBN,NSL) = 0.D+0
              ENDDO
            ELSE
!
!---          Loop over the number of chain-decay series  ---
!
              MC = 0
              DO NCS = 1,NBCDS
                MC = MC + 1
                MC0 = MC
                NSOLC = IBCDS(MC)
                DO NSL = 1,NSOLC
                  MC = MC + 1
                  ICX(NSL) = IBCDS(MC)
                  CX(NSL) = C_BH(NBN,IBCDS(MC))
                ENDDO
!
!---            Radioactive chain decay  ---
!
                INDX = 1
!
!---            Radioactive decay of individual solutes  ---
!
                IF( NCS.EQ.NBCDS ) INDX = 0
                CALL BATEMAN( CX,SL_BH(2,NBN),ICX,NBCDP(NCS),NSOLC,
     &            IBCDP(1,1,NCS),INDX )
                MC = MC0
                DO NSL = 1,NSOLC
                  MC = MC + 1
                  C_BH(NBN,IBCDS(MC)) = CX(NSL)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
!---  Loop over number of solutes  ---
!
      DO NSL = 1,NSOLU
!
!---    Load old sub-time-step concentrations  ---
!
        IF( ISLC(17).NE.0 ) CALL UPDTCO( NSL )
!
!---    Load old sub-time-step concentrations for fracture flow  ---
!
        IF( ISLC(17).NE.0 .AND. (ISLC(74).EQ.1 .OR. ISLC(74).EQ.3) )
     &    CALL UPDTCO_FRC( NSL )
!
!---    Load old sub-time-step concentrations for borehole flow  ---
!
        IF( ISLC(17).NE.0 .AND. (ISLC(74).EQ.2 .OR. ISLC(74).EQ.3) )
     &    CALL UPDTCO_BH( NSL )
      ENDDO
!
!---  End of CHAIN_DECAY group  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTC(NSL)
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
!     Updates concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by ML Rockhold, Battelle, PNL, May 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE JACOB
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDTC'
!
!---  Put the B array into the A array, skipping masked nodes  ---
!
      DO N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
        NMD = IXP(N)
        C(N,NSL) = BLU(NMD)
        IF( ABS(C(N,NSL)).LT.EPSL ) C(N,NSL) = 0.D+0
      ENDDO
!
!---  End of UPDTC group  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTCO(NSL)
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
!     Load old time step concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by ML Rockhold, Battelle, PNL, May 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDTCO'
      DO N = 1,NFLD
        IF( IXP(N).LE.0 ) CYCLE
        CO(N,NSL) = C(N,NSL)
      ENDDO
!
!---  End of UPDTCO group  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTC_FRC( NSL )
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
!     Update concentrations for fractures
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 July 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE SOLTN
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDTC_FRC'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fracture triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
          MTX = NFLD - NXP + NTX
          C_FRC(NTX,NSL) = BLU(MTX)
          IF( ABS(C_FRC(NTX,NSL)).LT.EPSL ) C_FRC(NTX,NSL) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDTC_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTC_BH( NSL )
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
!     Update concentrations for borehole nodes
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE SOLTN
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDTC_BH'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
          MBN = NFLD - NXP + NT_FRC + NBN 
          C_BH(NBN,NSL) = BLU(MBN)
          IF( ABS(C_BH(NBN,NSL)).LT.EPSL ) C_BH(NBN,NSL) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDTC_BH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTCO_FRC( NSL )
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
!     Load old time step concentrations for fractures
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE SOLTN
      USE GEOM_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDTCO_FRC'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fracture triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
          CO_FRC(NTX,NSL) = C_FRC(NTX,NSL)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDTCO_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTCO_BH( NSL )
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
!     Load old time step concentrations for fractures
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 July 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE SOLTN
      USE PARM_BH
      USE GEOM_BH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDTCO_BH'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
          CO_BH(NBN,NSL) = C_BH(NBN,NSL)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDTCO_BH group  ---
!
      RETURN
      END

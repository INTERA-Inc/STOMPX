!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TPORT_CO2( NSL )
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
!     STOMPX-CO2
!
!     Solute/Reactive Species Transport Shell.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 January 2022
!

!----------------------LIS Modules-------------------------------------!
!
      USE LIS_STOMP
!







!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE GRID
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/TPORT_CO2'
!
!---  Zero Jacobian matrix ---
!
      CALL JCBTZ_CO2
!
!---  Compute solute sources ---
!
      CALL SORT_CO2( NSL )
!
!---  Compute solute sources from injection and production wells ---
!
      CALL SORT_COUP_WELL( NSL )
!
!---  Zero solute transport fluxes  ---
!
      CALL SFXZ( NSL )
!
!---  Load Jacobian matrix (aqueous-phase transport)  ---
!
      CALL SJCBL( NSL )
!
!---  Load Jacobian matrix (gas-phase transport)  ---
!
      CALL SJCBG( NSL )
!
!---  Modify Jacobian matrix for boundary conditions ---
!
      CALL SBND_CO2( NSL )
!
!---  Set values of the Jacobian matrix  ---
!
      CALL JCBT_SV
!
!---  Solve the linear system A x = b for solute transport  ---
!

      INDX = 1
      CALL STOMP_LIS_SOLVE( T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,
     &    NUKTO(ID+1),NUKTL(ID+1),INDX )






      IF( ICNV.EQ.4 ) RETURN
!
!---  Check for fatal execution errors and stop simulation
!     if detected  ---
!
      CALL CHK_ERROR
!
!---  Update solute concentrations ---
!
      CALL UPDTC( NSL )

!
!---  Update solute concentrations on ghost cells  ---
!
      CALL UPDTC_GC( NSL )

!
!---  Compute solute aqueous-phase fluxes (interior nodes)  ---
!
      CALL SFXL( NSL )
!
!---  Compute solute gas-phase fluxes (interior nodes)  ---
!
      CALL SFXG( NSL )
!
!---  Compute solute aqueous and gas fluxes (boundary surfaces)  ---
!
      CALL SFXB_CO2( NSL )
!
!---  Integrate solute sources  ---
!
      CALL SORIT_CO2( NSL )
!
!---  Load old sub-time-step concentrations  ---
!
      IF( ISLC(17).NE.0 ) CALL UPDTCO( NSL)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TPORT_CO2 group  ---
!
      RETURN
      END

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
!    1 decay chains." J. Appl. Phys., 101, 094907.  ---
!
!----------------------Authors-----------------------------------------!
!
!    Written by Mark D. White, PNNL, 27 January 2022.
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
          LAMBDAX = LOG(2.D+0)/HLF(ISOLUX)
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
        LMDX(I) = LOG(2.D+0)/HLF(ISOLUX)
      ENDDO
!
!---  Wikipedia Bateman equation
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
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BATEMAN group  ---
!
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
!     Written by M.D. White, PNNL, 27 January 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
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
      REAL*8, DIMENSION(LCDS) :: CX
      REAL*8 :: LAMBDAX
      INTEGER, DIMENSION(LCDS) :: ICX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHAIN_DECAY'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = ND(N)
!
!---    Zero chain-decay series  ---
!
        IF( NBCDS.EQ.0 ) THEN
          DO NSL = 1,NSOLU
            LAMBDAX = LOG(2.D+0)/HLF(NSL)
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
!---  Loop over number of solutes  ---
!
      DO NSL = 1,NSOLU
!
!---    Load old sub-time-step concentrations  ---
!
        IF( ISLC(17).NE.0 ) CALL UPDTCO( NSL )
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHAIN_DECAY group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC_CO2
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
!     STOMPX-CO2
!
!     Compute initial solute concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 27 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE GRID
      USE FDVP
      USE CONST
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( IEQC.EQ.0 .AND. ISLC(40).EQ.0 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CISC_CO2'
      DO NSL = 1,NSOLU
!
!---    Loop over all nodes, skipping inactive nodes  ---
!
        DO N = 1,NFCGC(ID+1)
          IF( IXP(N).EQ.0 ) CYCLE
          N_DB = ND(N)
!
!---      Load old-time-step concentrations  ---
!
          CO(N,NSL) = C(N,NSL)
        ENDDO
!
!---    Assign boundary solute concentrations for initial condition
!       type boundary conditions  ---
!
        DO NB = 1,NBC(ID+1)
          IF( IBCT(NSL+LUK+LPH,NB).EQ.12 ) THEN
            N = IBCN(NB)
            CBO(NB,NSL) = C(N,NSL)
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CISC_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CRN_LIM( NSL )
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
!     STOMPX-CO2
!
!     Compute a time step that globally satisfies the Courant
!     number limit.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 January 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE HYST
      USE GRID
      USE FLUX
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CRN_LIM'
!
!---  Compute the maximum Courant number ---
!
      CRNLX = 0.D+0
      CRNGX = 0.D+0
!
!---  Courant number control ---
!
      IF( ISLC(17).EQ.1 ) THEN
!
!---    Loop over local nodes  ---
!
        DO N = 1,NFCGC(ID+1)
!
!---      Skip for inactive nodes or ghost cells  ---
!
          IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
          CRNLX = MAX( CRNTL(N),CRNLX )
          CRNGX = MAX( CRNTG(N),CRNGX )
        ENDDO
!
!---  Special Vadose Zone Courant Number Control ---
!
      ELSEIF( ISLC(17).EQ.2 ) THEN
!
!---    Loop over local nodes  ---
!
        DO N = 1,NFCGC(ID+1)
!
!---      Skip for inactive nodes or ghost cells  ---
!
          IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
          IF( 1.D+0-SL(2,N).LT.EPSL ) CYCLE
          CLX = C(N,NSL)*YL(N,NSL)/(SL(2,N)*PORD(2,N)+EPSL)
          IF( CLX.GT.CCL_CRN(NSL) ) CRNLX = MAX( CRNTL(N),CRNLX )
        ENDDO
!
!---  Aqueous-Only Courant Number Control ---
!
      ELSEIF( ISLC(17).EQ.3 ) THEN
!
!---    Loop over local nodes  ---
!
        DO N = 1,NFCGC(ID+1)
!
!---      Skip for inactive nodes or ghost cells  ---
!
          IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
          CRNLX = MAX( CRNTL(N),CRNLX )
        ENDDO
      ENDIF
      CRNMLX = MAX( CRNLX,CRNGX )
!
!---  Maximum Courant number  ---
!
      CALL MPI_ALLREDUCE( CRNMLX,CRNMX,ISVC,MPI_REAL8,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
      DT_CRN = DT
      DTI_CRN = DTI
      TM_CRN = TM
      TM = TM-DT
      IF( CRNMX.GT.CRNTMXT ) THEN
        N_CRN(NSL) = INT(CRNMX/CRNTMXT) + 1
        REALX = REAL(N_CRN(NSL))
        DT = DT_CRN/REALX
        DTI = 1.D+0/(DT+EPSL)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CRN_LIM group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CRNTNB
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
!     STOMPX-CO2
!
!     Compute the local maximum Courant numbers.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 6 January 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE HYST
      USE GRID
      USE FLUX
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
      REAL*8 KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CRNTNB'
!
!---  Loop over local nodes  ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    Aqueous Phase Courant Number  ---
!
        CRNTL(N) = MAX( ABS(UL(1,1,N)*DT/DXGP(1,N)),
     &    ABS(UL(1,2,N)*DT/DXGP(2,N)),
     &    ABS(VL(1,1,N)*DT/(DYGP(1,N)*RP(N))),
     &    ABS(VL(1,2,N)*DT/(DYGP(2,N)*RP(N))),
     &    ABS(WL(1,1,N)*DT/DZGP(1,N)),
     &    ABS(WL(1,2,N)*DT/DZGP(2,N)) )
        IF( SL(2,N)*PORD(2,N).GT.EPSL ) THEN
          CRNTL(N) = CRNTL(N)/(SL(2,N)*PORD(2,N))
        ELSE
          CRNTL(N) = 0.D+0
        ENDIF
!
!---    Gas Phase Courant Number  ---
!
        CRNTG(N) = MAX( ABS(UG(1,1,N)*DT/DXGP(1,N)),
     &    ABS(UG(1,2,N)*DT/DXGP(2,N)),
     &    ABS(VG(1,1,N)*DT/(DYGP(1,N)*RP(N))),
     &    ABS(VG(1,2,N)*DT/(DYGP(2,N)*RP(N))),
     &    ABS(WG(1,1,N)*DT/DZGP(1,N)),
     &    ABS(WG(1,2,N)*DT/DZGP(2,N)) )
        IF( (SG(2,N)-SGT(2,N))*PORD(2,N).GT.EPSL ) THEN
          CRNTG(N) = CRNTG(N)/((SG(2,N)-SGT(2,N))*PORD(2,N))
        ELSE
          CRNTG(N) = 0.D+0
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CRNTNB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_SV
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
!     Set values of the transport Jacobian matrix
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 March 2022
!

!----------------------LIS Modules-------------------------------------!
!
      USE LIS_STOMP
!







!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
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
      REAL*8, DIMENSION(NNZTR) :: BUFFER
      INTEGER, DIMENSION(NNZTR) :: ICOL
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT_SV'

!
!---  Lis solver  ---
!
      CALL lis_matrix_set_csr_f( NNZC,NLUC,MLUC,DLU,T_MAT,IERR )
      IF( IERR.NE.0 ) THEN
        I_ERR(1) = 0
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
        M_ERR(1) = 'Solute Transport Lis Matrix Set CSR'
        RETURN
      ENDIF

!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_SV group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTZ_CO2
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
!     Zero the transport Jacobian matrix.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 March 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTZ_CO2'
!
!---  Deallocate problem/solution vector array  ---
!
      IF( ALLOCATED(BLU) ) THEN
        DEALLOCATE( BLU,STAT=ISTAT )
        CHMSG = 'BLU'
        CALL DEALLOC_ERROR( CHMSG,ISTAT )
      ENDIF
!
!---  Allocate problem/solution vector array for solute transport  ---
!
      ALLOCATE( BLU(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'BLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Zero problem vector array for solute transport  ---
!
      DO N = 1,NFCGC(ID+1)
        BLU(N) = 0.D+0
      ENDDO
!
!---  Deallocate Jacobian matrix array  ---
!
      IF( ALLOCATED(DLU) ) THEN
        DEALLOCATE( DLU,STAT=ISTAT )
        CHMSG = 'DLU'
        CALL DEALLOC_ERROR( CHMSG,ISTAT )
      ENDIF
!
!---  Allocate Jacobian matrix array for solute transport  ---
!
      ALLOCATE( DLU(1:NNZC),STAT=ISTAT )
      CHMSG = 'DLU'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Zero Jacobian matrix array for solute transport  ---
!
      DO N = 1,NNZC
        DLU(N) = 0.D+0
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTZ_CO2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SBND_CO2( NSL )
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
!     STOMPX-CO2
!
!     Modify the Jacobian matrix for the solute transport equation
!     to incorporate boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 January 2022

!----------------------LIS Modules-------------------------------------!
!
      USE LIS_STOMP
!







!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PROP
      USE JACOB
      USE HYST
      USE GRID
      USE FLUX
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)






!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 BCX(LSPBC+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SBND_CO2'
!
!---  Loop over boundary conditions  ---
!
      DO NB = 1,NBC(ID+1)
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) CYCLE
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) CYCLE
!
!---    Solute transport  ---
!
        IF( NSL.LE.NSOLU ) THEN
          IBCTX = IBCT(NSL+LUK,NB)
!
!---    Reactive species transport  ---
!
        ELSE
          IBCTX = IBCT(NSOLU+LUK+1,NB)
        ENDIF
!
!---    Zero flux boundary condition  ---
!
        IF( IBCTX.EQ.3 ) CYCLE
!
!---    Single boundary condition time  ---
!
        IF( IBCM(NB).EQ.1 ) THEN
!
!---      Solute transport  ---
!
          IF( NSL.LE.NSOLU ) THEN
            BCX(1) = BC(NSL+LBCU,1,MB)
            IF( IBCTX.EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---      Reactive species transport  ---
!
          ELSE
            BCX(1) = 0.D+0
            DO NSPX = 1,IBCSP(1,NB)
              NSP = IBCSP(NSPX+1,NB)
              MX = NSOLU+LBCU+NSPX
              BCX(NSPX+1) = BC(MX,1,MB)
!
!---          Aqueous species ---
!
              IF( NSP.LE.NSPL ) THEN
                IF( IBCT(NSOLU+LUK+1,NB).EQ.12 ) 
     &            BCX(NSPX+1) = SP_CBO(NB,NSP)
!
!---          Gas species ---
!
              ELSE
                IF( IBCT(NSOLU+LUK+2,NB).EQ.12 ) 
     &            BCX(NSPX+1) = SP_CBO(NB,NSP)
              ENDIF
            ENDDO
          ENDIF
!
!---    Multiple boundary condition times  ---
!
        ELSE
          IFIND = 0
          DO M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
              TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
              DTBC = MIN( BC(1,M,MB)-TMZ,DT )
              TFBC = (TMZ-5.D-1*DTBC-BC(1,M-1,MB))/TDBC
!
!---          Solute transport  ---
!
              IF( NSL.LE.NSOLU ) THEN
                BCX(1) = BC(NSL+LBCU,M-1,MB) +
     &            TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
                IF( IBCT(NSL+LUK,NB).EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---          Reactive species transport  ---
!
              ELSE
                BCX(1) = 0.D+0
                DO NSPX = 1,IBCSP(1,NB)
                  NSP = IBCSP(NSPX+1,NB)
                  MX = NSOLU+LBCU+NSPX
                  BCX(NSPX+1) = BC(MX,M-1,MB) +
     &              TFBC*(BC(MX,M,MB)-BC(MX,M-1,MB))
!
!---              Aqueous species ---
!
                  IF( NSP.LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.12 ) 
     &                BCX(NSPX+1) = SP_CBO(NB,NSP)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.12 ) 
     &                BCX(NSPX+1) = SP_CBO(NB,NSP)
                  ENDIF
                ENDDO
              ENDIF
              IFIND = 1
              EXIT
            ENDIF
          ENDDO
          IF( IFIND.EQ.0 ) CYCLE
        ENDIF
        N = IBCN(NB)
        N_DB = -NB
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
        MP = IXP(N)

        MA = 1
        MCOL = KLUC(MP-IEQC_OFFSET,MA)
        MA = MA + 1






!
!---  Diffusion coefficients at node adjacent to boundary  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SMDLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
        DLP = TORL(2,N)*SL(2,N)*PORD(2,N)*SMDLP
        PCOR = (PG(2,N)+PATM)/PATM
        SMDGP = SMDG(NSL)*(TCOR**1.75)/PCOR
        DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SMDGP
!
!---  Phase fraction factors at node adjacent to boundary  ---
!
        XVLP = SL(2,N)*PORD(2,N)
        FCLP = 0.D+0
        IF( XVLP.GT.SMALL ) FCLP = YL(N,NSL)/XVLP
        XVGP = SG(2,N)*PORD(2,N)
        FCGP = 0.D+0
        IF( XVGP.GT.SMALL ) FCGP = YG(N,NSL)/XVGP
!
!---    Solute transport only, skip calculations for reactive
!       species transport  ---
!
        XVLB = SLB(2,NB)*PORDB(2,NB)
        XVGB = SGB(2,NB)*PORDB(2,NB)
        IF( NSL.LE.NSOLU ) THEN
!
!---      Phase fraction factors at boundary  ---
!
          IF( IPCL(NSL).EQ.2 ) THEN
            XVSB = RHOS(N)*PCSL(1,N,NSL)*(1.D+0-PORTB(2,NB))
     &        *SLB(2,NB)
          ELSE
            XVSB = RHOS(N)*PCSL(1,N,NSL)*(1.D+0-PORTB(2,NB))
          ENDIF
!
!---      Constant gas-aqueous partition coefficient  ---
!
          IF( IPCGL(NSL).EQ.0 ) THEN
            PCGLX = PCGL(1,NSL)
!
!---      Temperature dependent gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.1 ) THEN
            TK = TB(2,NB)+TABS
            PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &        + PCGL(3,NSL)*LOG(TK)
     &        + PCGL(4,NSL)*TK + PCGL(5,NSL)*TK**2 )
!
!---      Water-vapor 1 gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.2 ) THEN
            PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
          ENDIF
          PCGLX = MAX( PCGLX,1.D-20 )
          PCGLX = MIN( PCGLX,1.D+20 )
!
!---      Phase-volumetric concentration ratios  ---
!
          FCL = 1.D+0/(XVSB + XVLB + XVGB*PCGLX)
          FCG = 1.D+0/((XVSB + XVLB)/PCGLX + XVGB)
!
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = XVLB*FCL
          YGB(NB,NSL) = XVGB*FCG
!
!---      Convert boundary phase concentrations to
!         volumetric concentrations  ---
!
          IF( IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.14 .OR.
     &      IBCT(NSL+LUK,NB).EQ.23 ) THEN
            BCX(1) = BCX(1)/(FCL+SMALL)
          ELSEIF( IBCT(NSL+LUK,NB).EQ.9 .OR.
     &      IBCT(NSL+LUK,NB).EQ.15 .OR.
     &      IBCT(NSL+LUK,NB).EQ.43 ) THEN
            BCX(1) = BCX(1)/(FCG+SMALL)
          ENDIF
          CB(NB,NSL) = BCX(1)
        ELSE
!
!---      Convert species concentrations to total-component
!         concentrations  ---
!
          IF( NSL.LE.NSOLU+NEQC ) THEN
            NEQ = NSL-NSOLU
            YSPLX = 0.D+0
            YSPGX = 0.D+0
            DO NSP = 1,IEQ_C(1,NEQ)
              DO NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_C(NSP+1,NEQ) ) THEN
!
!---              Aqueous species ---
!
                  IF( IEQ_C(NSP+1,NEQ).LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+1,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPLX = YSPLX + EQ_C(NSP,NEQ)*BCX(NSPX+1)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+2,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPGX = YSPGX + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                  ENDIF                    
                  BCX(1) = BCX(1) + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
              ENDDO
            ENDDO
!
!---        Linked aqueous CO2   ---
!
            IF( ISPLK(6).EQ.NSL ) BCX(1) = 1.D+3*XLAB(2,NB)*
     &          RHOLB(2,N)*SLB(2,NB)*PORDB(2,NB)/WTMA
!
!---        Convert species concentrations to total-kinetic
!           concentrations  ---
!
          ELSEIF( NSL.LE.NSOLU+NEQC+NEQK ) THEN
            NEQ = NSL-NSOLU-NEQC
            YSPLX = 0.D+0
            DO NSP = 1,IEQ_K(1,NEQ)
              DO NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_K(NSP+1,NEQ) ) THEN
!
!---              Aqueous species ---
!
                  IF( IEQ_K(NSP+1,NEQ).LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+1,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPLX = YSPLX + EQ_K(NSP,NEQ)*BCX(NSPX+1)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+2,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                  ENDIF                    
                  BCX(1) = BCX(1) + EQ_K(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
              ENDDO
            ENDDO
!
!---        Linked aqueous CO2   ---
!
            IF( ISPLK(6).EQ.NSL ) BCX(1) = 1.D+3*XLAB(2,NB)*
     &          RHOLB(2,N)*SLB(2,NB)*PORDB(2,NB)/WTMA
          ENDIF
          IF( ABS(BCX(1))/EPSL.LT.EPSL ) THEN
            YSPLX = 0.D+0
          ELSE
            YSPLX = YSPLX/BCX(1)
          ENDIF
!
!---      Phase-volumetric concentration ratios  ---
!
          YLBX = MAX( MIN( 1.D+0,YSPLX ),0.D+0 )
          YGBX = 1.D+0-YLBX
          FCL = 0.D+0
          IF( XVLB/EPSL.GT.EPSL ) FCL = YLBX/XVLB
          FCG = 0.D+0
          IF( XVGB/EPSL.GT.EPSL ) FCG = YGBX/XVGB
!
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = YLBX
          YGB(NB,NSL) = YGBX
          CB(NB,NSL) = BCX(1)
        ENDIF
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULBX = 5.D-1*(UL(1,1,N)+UL(1,2,N))
            VLBX = 5.D-1*(VL(1,1,N)+VL(1,2,N))
            WLBX = WL(1,1,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WLBX,INDX)
              ULBX = ULBX/VMCX
              VLBX = VLBX/VMCX
              WLBX = WLBX/VMCX
            ENDIF
            ULBSQ = ULBX*ULBX
            VLBSQ = VLBX*VLBX
            WLBSQ = WLBX*WLBX
            ZLVB = SQRT(ULBSQ+VLBSQ+WLBSQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLB = (DISPL(N)*WLBSQ + DISPT(N)*(ULBSQ+VLBSQ))/
     &        (ZLVB+SMALL)
            UGBX = 5.D-1*(UG(1,1,N)+UG(1,2,N))
            VGBX = 5.D-1*(VG(1,1,N)+VG(1,2,N))
            WGBX = WG(1,1,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WGBX,INDX)
              UGBX = UGBX/VMCX
              VGBX = VGBX/VMCX
              WGBX = WGBX/VMCX
            ENDIF
            UGBSQ = UGBX*UGBX
            VGBSQ = VGBX*VGBX
            WGBSQ = WGBX*WGBX
            ZGVB = SQRT(UGBSQ+VGBSQ+WGBSQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGB = (DISPL(N)*WGBSQ + DISPT(N)*(UGBSQ+VGBSQ))/
     &        (ZGVB+SMALL)
          ELSE
            DPLB = 0.D+0
            DPGB = 0.D+0
          ENDIF
          FLB = AFZ(1,N)*WL(1,1,N)
          FGB = AFZ(1,N)*WG(1,1,N)
          CRLB = ABS( WL(1,1,N) )*DT/(DZGF(N)*XVLB+SMALL)
          CRGB = ABS( WG(1,1,N) )*DT/(DZGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLZ = DIFMN(DLB,DLP,DZGF(N),DZGF(N),WL(1,1,N),INDX)
            DLZ = AFZ(1,N)*(DLZ+DPLB)/(5.D-1*DZGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGZ = DIFMN(DGB,DGP,DZGF(N),DZGF(N),WG(1,1,N),INDX)
            DGZ = AFZ(1,N)*(DGZ+DPGB)/(5.D-1*DZGF(N))
            ALB = MAX( FLB,ZERO ) +
     &        DLZ*MAX((ONE-(TENTH*ABS(FLB)/(DLZ+SMALL)))**5,ZERO)
            AGB = MAX( FGB,ZERO ) +
     &        DGZ*MAX((ONE-(TENTH*ABS(FGB)/(DGZ+SMALL)))**5,ZERO)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLB = MIN( FLB,0.D+0 )
            FGB = MIN( FGB,0.D+0 )
            ALB = MAX( FLB,ZERO )
            AGB = MAX( FGB,ZERO )
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLB = MAX( FLB,0.D+0 )
            FGB = MAX( FGB,0.D+0 )
            ALB = MAX( FLB,ZERO )
            AGB = MAX( FGB,ZERO )
          ENDIF
          AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP
          AB = ALB*FCL + AGB*FCG
!
!---      Solution vector  ---
!

          BUFFER = AB*BCX(1)
          CALL lis_vector_set_value_f( 1,MP,BUFFER,
     &      T_RHS_VEC,IERR )






!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULSX = 5.D-1*(UL(1,1,N)+UL(1,2,N))
            VLSX = VL(1,1,N)
            WLSX = 5.D-1*(WL(1,1,N)+WL(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DYGF(N),DYGF(N),VLBX,INDX)
              ULSX = ULSX/VMCX
              VLSX = VLSX/VMCX
              WLSX = WLSX/VMCX
            ENDIF
            ULSSQ = ULSX*ULSX
            VLSSQ = VLSX*VLSX
            WLSSQ = WLSX*WLSX
            ZLVS = SQRT(ULSSQ+VLSSQ+WLSSQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLS = (DISPL(N)*VLSSQ + DISPT(N)*(ULSSQ+WLSSQ))/
     &        (ZLVS+SMALL)
            UGSX = 5.D-1*(UG(1,1,N)+UG(1,2,N))
            VGSX = VG(1,1,N)
            WGSX = 5.D-1*(WG(1,1,N)+WG(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DYGF(N),DYGF(N),VGSX,INDX)
              UGSX = UGSX/VMCX
              VGSX = VGSX/VMCX
              WGSX = WGSX/VMCX
            ENDIF
            UGSSQ = UGSX*UGSX
            VGSSQ = VGSX*VGSX
            WGSSQ = WGSX*WGSX
            ZGVS = SQRT(UGSSQ+VGSSQ+WGSSQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGS = (DISPL(N)*VGSSQ + DISPT(N)*(WGSSQ+VGSSQ))/
     &        (ZGVS+SMALL)
          ELSE
            DPLS = 0.D+0
            DPGS = 0.D+0
          ENDIF
          FLS = AFY(1,N)*VL(1,1,N)
          FGS = AFY(1,N)*VG(1,1,N)
          CRLS = ABS( VL(1,1,N) )*DT/(RP(N)*DYGF(N)*XVLB+SMALL)
          CRGS = ABS( VG(1,1,N) )*DT/(RP(N)*DYGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLY = DIFMN(DLB,DLP,DYGF(N),DYGF(N),VL(1,1,N),INDX)
            DLY = AFY(1,N)*(DLY+DPLS)/RP(N)/(5.D-1*DYGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGY = DIFMN(DGB,DGP,DYGF(N),DYGF(N),VG(1,1,N),INDX)
            DGY = AFY(1,N)*(DGY+DPGS)/RP(N)/(5.D-1*DYGF(N))
            ALS = MAX( FLS,ZERO ) +
     &        DLY*MAX((ONE-(TENTH*ABS(FLS)/(DLY+SMALL)))**5,ZERO)
            AGS = MAX( FGS,ZERO ) +
     &        DGY*MAX((ONE-(TENTH*ABS(FGS)/(DGY+SMALL)))**5,ZERO)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLS = MIN( FLS,0.D+0 )
            FGS = MIN( FGS,0.D+0 )
            ALS = MAX( FLS,ZERO )
            AGS = MAX( FGS,ZERO )
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLS = MAX( FLS,0.D+0 )
            FGS = MAX( FGS,0.D+0 )
            ALS = MAX( FLS,ZERO )
            AGS = MAX( FGS,ZERO )
          ENDIF
          AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP
          AS = ALS*FCL + AGS*FCG
!
!---      Solution vector  ---
!

          BUFFER = AS*BCX(1)
          CALL lis_vector_set_value_f( 1,MP,BUFFER,
     &      T_RHS_VEC,IERR )






!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULWX = UL(1,1,N)
            VLWX = 5.D-1*(VL(1,1,N)+VL(1,2,N))
            WLWX = 5.D-1*(WL(1,1,N)+WL(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DXGF(N),DXGF(N),ULBX,INDX)
              ULWX = ULWX/VMCX
              VLWX = VLWX/VMCX
              WLWX = WLWX/VMCX
            ENDIF
            ULWSQ = ULWX*ULWX
            VLWSQ = VLWX*VLWX
            WLWSQ = WLWX*WLWX
            ZLVW = SQRT(ULWSQ+VLWSQ+WLWSQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLW = (DISPL(N)*VLWSQ + DISPT(N)*(ULWSQ+WLWSQ))/
     &        (ZLVW+SMALL)
            UGWX = UG(1,1,N)
            VGWX = 5.D-1*(VG(1,1,N)+VG(1,2,N))
            WGWX = 5.D-1*(WG(1,1,N)+WG(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DXGF(N),DXGF(N),UGWX,INDX)
              UGWX = UGWX/VMCX
              VGWX = VGWX/VMCX
              WGWX = WGWX/VMCX
            ENDIF
            UGWSQ = UGWX*UGWX
            VGWSQ = VGWX*VGWX
            WGWSQ = WGWX*WGWX
            ZGVW = SQRT(UGWSQ+VGWSQ+WGWSQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGW = (DISPL(N)*VGWSQ + DISPT(N)*(WGWSQ+VGWSQ))/
     &        (ZGVW+SMALL)
          ELSE
            DPLW = 0.D+0
            DPGW = 0.D+0
          ENDIF
          FLW = AFX(1,N)*UL(1,1,N)
          FGW = AFX(1,N)*UG(1,1,N)
          CRLW = ABS( UL(1,1,N) )*DT/(DXGF(N)*XVLB+SMALL)
          CRGW = ABS( UG(1,1,N) )*DT/(DXGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLX = DIFMN(DLB,DLP,DXGF(N),DXGF(N),UL(1,1,N),INDX)
            DLX = AFX(1,N)*(DLX+DPLW)/(5.D-1*DXGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGX = DIFMN(DGB,DGP,DXGF(N),DXGF(N),UG(1,1,N),INDX)
            DGX = AFX(1,N)*(DGX+DPGW)/(5.D-1*DXGF(N))
            ALW = MAX(FLW,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLW)/(DLX+SMALL)))**5,ZERO)
            AGW = MAX(FGW,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGW)/(DGX+SMALL)))**5,ZERO)
            AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
            AW = ALW*FCL + AGW*FCG
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLW = MIN( FLW,0.D+0 )
            FGW = MIN( FGW,0.D+0 )
            ALW = MAX(FLW,ZERO)
            AGW = MAX(FGW,ZERO)
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLW = MAX( FLW,0.D+0 )
            FGW = MAX( FGW,0.D+0 )
              ALW = MAX(FLW,ZERO)
              AGW = MAX(FGW,ZERO)
          ENDIF
          AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
          AW = ALW*FCL + AGW*FCG
!
!---      Solution vector  ---
!

          BUFFER = AW*BCX(1)
          CALL lis_vector_set_value_f( 1,MP,BUFFER,
     &      T_RHS_VEC,IERR )






!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULEX = UL(1,2,N)
            VLEX = 5.D-1*(VL(1,1,N)+VL(1,2,N))
            WLEX = 5.D-1*(WL(1,1,N)+WL(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DXGF(N),DXGF(N),ULBX,INDX)
              ULEX = ULEX/VMCX
              VLEX = VLEX/VMCX
              WLEX = WLEX/VMCX
            ENDIF
            ULESQ = ULEX*ULEX
            VLESQ = VLEX*VLEX
            WLESQ = WLEX*WLEX
            ZLVE = SQRT(ULESQ+VLESQ+WLESQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLE = (DISPL(N)*VLESQ + DISPT(N)*(ULESQ+WLESQ))/
     &        (ZLVE+SMALL)
            UGEX = UG(1,2,N)
            VGEX = 5.D-1*(VG(1,1,N)+VG(1,2,N))
            WGEX = 5.D-1*(WG(1,1,N)+WG(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DXGF(N),DXGF(N),UGEX,INDX)
              UGEX = UGEX/VMCX
              VGEX = VGEX/VMCX
              WGEX = WGEX/VMCX
            ENDIF
            UGESQ = UGEX*UGEX
            VGESQ = VGEX*VGEX
            WGESQ = WGEX*WGEX
            ZGVE = SQRT(UGESQ+VGESQ+WGESQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGE = (DISPL(N)*VGESQ + DISPT(N)*(WGESQ+VGESQ))/
     &        (ZGVE+SMALL)
          ELSE
            DPLE = 0.D+0
            DPGE = 0.D+0
          ENDIF
          FLE = AFX(2,N)*UL(1,2,N)
          FGE = AFX(2,N)*UG(1,2,N)
          CRLE = ABS( UL(1,2,N) )*DT/(DXGF(N)*XVLB+SMALL)
          CRGE = ABS( UG(1,2,N) )*DT/(DXGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLX = DIFMN(DLP,DLB,DXGF(N),DXGF(N),UL(1,2,N),INDX)
            DLX = AFX(2,N)*(DLX+DPLE)/(5.D-1*DXGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGX = DIFMN(DGP,DGB,DXGF(N),DXGF(N),UG(1,2,N),INDX)
            DGX = AFX(2,N)*(DGX+DPGE)/(5.D-1*DXGF(N))
            ALE = MAX( -FLE,ZERO ) +
     &        DLX*MAX((ONE-(TENTH*ABS(FLE)/(DLX+SMALL)))**5,ZERO)
            AGE = MAX( -FGE,ZERO ) +
     &        DGX*MAX((ONE-(TENTH*ABS(FGE)/(DGX+SMALL)))**5,ZERO)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLE = MAX( FLE,0.D+0 )
            FGE = MAX( FGE,0.D+0 )
            ALE = MAX( -FLE,ZERO )
            AGE = MAX( -FGE,ZERO )
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLE = MIN( FLE,0.D+0 )
            FGE = MIN( FGE,0.D+0 )
            ALE = MAX( -FLE,ZERO )
            AGE = MAX( -FGE,ZERO )
          ENDIF
          AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP
          AE = ALE*FCL + AGE*FCG
!
!---      Solution vector  ---
!

          BUFFER = AE*BCX(1)
          CALL lis_vector_set_value_f( 1,MP,BUFFER,
     &      T_RHS_VEC,IERR )






!
!---    North boundary
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULNX = 5.D-1*(UL(1,1,N)+UL(1,2,N))
            VLNX = VL(1,2,N)
            WLNX = 5.D-1*(WL(1,1,N)+WL(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DYGF(N),DYGF(N),VLBX,INDX)
              ULNX = ULNX/VMCX
              VLNX = VLNX/VMCX
              WLNX = WLNX/VMCX
            ENDIF
            ULNSQ = ULNX*ULNX
            VLNSQ = VLNX*VLNX
            WLNSQ = WLNX*WLNX
            ZLVN = SQRT(ULNSQ+VLNSQ+WLNSQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLN = (DISPL(N)*VLNSQ + DISPT(N)*(ULNSQ+WLNSQ))/
     &        (ZLVN+SMALL)
            UGNX = 5.D-1*(UG(1,1,N)+UG(1,2,N))
            VGNX = VG(1,2,N)
            WGNX = 5.D-1*(WG(1,1,N)+WG(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DYGF(N),DYGF(N),VGNX,INDX)
              UGNX = UGNX/VMCX
              VGNX = VGNX/VMCX
              WGNX = WGNX/VMCX
            ENDIF
            UGNSQ = UGNX*UGNX
            VGNSQ = VGNX*VGNX
            WGNSQ = WGNX*WGNX
            ZGVN = SQRT(UGNSQ+VGNSQ+WGNSQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGN = (DISPL(N)*VGNSQ + DISPT(N)*(WGNSQ+VGNSQ))/
     &        (ZGVN+SMALL)
          ELSE
            DPLN = 0.D+0
            DPGN = 0.D+0
          ENDIF
          FLN = AFY(2,N)*VL(1,2,N)
          FGN = AFY(2,N)*VG(1,2,N)
          CRLN = ABS( VL(1,2,N) )*DT/(RP(N)*DYGF(N)*XVLB+SMALL)
          CRGN = ABS( VG(1,2,N) )*DT/(RP(N)*DYGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLY = DIFMN(DLP,DLB,DYGF(N),DYGF(N),VL(1,2,N),INDX)
            DLY = AFY(2,N)*(DLY+DPLN)/RP(N)/(5.D-1*DYGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGY = DIFMN(DGP,DGB,DYGF(N),DYGF(N),VG(1,2,N),INDX)
            DGY = AFY(2,N)*(DGY+DPGN)/RP(N)/(5.D-1*DYGF(N))
            ALN = MAX( -FLN,ZERO ) +
     &        DLY*MAX((ONE-(TENTH*ABS(FLN)/(DLY+SMALL)))**5,ZERO)
            AGN = MAX( -FGN,ZERO ) +
     &        DGY*MAX((ONE-(TENTH*ABS(FGN)/(DGY+SMALL)))**5,ZERO)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLN = MAX( FLN,0.D+0 )
            FGN = MAX( FGN,0.D+0 )
            ALN = MAX( -FLN,ZERO )
            AGN = MAX( -FGN,ZERO )
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLN = MIN( FLN,0.D+0 )
            FGN = MIN( FGN,0.D+0 )
            ALN = MAX( -FLN,ZERO )
            AGN = MAX( -FGN,ZERO )
          ENDIF
          AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP
          AN = ALN*FCL + AGN*FCG
!
!---      Solution vector  ---
!

          BUFFER = AN*BCX(1)
          CALL lis_vector_set_value_f( 1,MP,BUFFER,
     &      T_RHS_VEC,IERR )






!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULTX = 5.D-1*(UL(1,1,N)+UL(1,2,N))
            VLTX = 5.D-1*(VL(1,1,N)+VL(1,2,N))
            WLTX = WL(1,2,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WLTX,INDX)
              ULTX = ULTX/VMCX
              VLTX = VLTX/VMCX
              WLTX = WLTX/VMCX
            ENDIF
            ULTSQ = ULTX*ULTX
            VLTSQ = VLTX*VLTX
            WLTSQ = WLTX*WLTX
            ZLVT = SQRT(ULTSQ+VLTSQ+WLTSQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLT = (DISPL(N)*WLTSQ + DISPT(N)*(ULTSQ+VLTSQ))/
     &        (ZLVT+SMALL)
            UGTX = 5.D-1*(UG(1,1,N)+UG(1,2,N))
            VGTX = 5.D-1*(VG(1,1,N)+VG(1,2,N))
            WGTX = WG(1,2,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WGTX,INDX)
              UGTX = UGTX/VMCX
              VGTX = VGTX/VMCX
              WGTX = WGTX/VMCX
            ENDIF
            UGTSQ = UGTX*UGTX
            VGTSQ = VGTX*VGTX
            WGTSQ = WGTX*WGTX
            ZGVT = SQRT(UGTSQ+VGTSQ+WGTSQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGT = (DISPL(N)*WGTSQ + DISPT(N)*(UGTSQ+VGTSQ))/
     &        (ZGVT+SMALL)
          ELSE
            DPLT = 0.D+0
            DPGT = 0.D+0
          ENDIF
          FLT = AFZ(2,N)*WL(1,2,N)
          FGT = AFZ(2,N)*WG(1,2,N)
          CRLT = ABS( WL(1,2,N) )*DT/(DZGF(N)*XVLB+SMALL)
          CRGT = ABS( WG(1,2,N) )*DT/(DZGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLZ = DIFMN(DLP,DLB,DZGF(N),DZGF(N),WL(1,2,N),INDX)
            DLZ = AFZ(2,N)*(DLZ+DPLT)/(5.D-1*DZGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGZ = DIFMN(DGP,DGB,DZGF(N),DZGF(N),WG(1,2,N),INDX)
            DGZ = AFZ(2,N)*(DGZ+DPGT)/(5.D-1*DZGF(N))
            ALT = MAX( -FLT,ZERO ) +
     &        DLZ*MAX((ONE-(TENTH*ABS(FLT)/(DLZ+SMALL)))**5,ZERO)
            AGT = MAX( -FGT,ZERO ) +
     &        DGZ*MAX((ONE-(TENTH*ABS(FGT)/(DGZ+SMALL)))**5,ZERO)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLT = MAX( FLT,0.D+0 )
            FGT = MAX( FGT,0.D+0 )
            ALT = MAX( -FLT,ZERO )
            AGT = MAX( -FGT,ZERO )
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLT = MIN( FLT,0.D+0 )
            FGT = MIN( FGT,0.D+0 )
            ALT = MAX( -FLT,ZERO )
            AGT = MAX( -FGT,ZERO )
          ENDIF
          AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP
          AT = ALT*FCL + AGT*FCG
!
!---      Solution vector  ---
!

          BUFFER = AT*BCX(1)
          CALL lis_vector_set_value_f( 1,MP,BUFFER,
     &      T_RHS_VEC,IERR )






        ENDIF

        DLU(MCOL) = DLU(MCOL) + AP




      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SBND_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXB_CO2( NSL )
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
!     STOMPX-CO2
!
!     Compute solute aqueous-phase fluxes on boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PROP
      USE HYST
      USE GRID
      USE FLUX
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 BCX(LSPBC+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFXB_CO2'
!
!---  Loop over boundary conditions  ---
!
      DO NB = 1,NBC(ID+1)
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) CYCLE
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) CYCLE
!
!---    Solute transport  ---
!
        IF( NSL.LE.NSOLU ) THEN
          IBCTX = IBCT(NSL+LUK,NB)
!
!---    Reactive species transport  ---
!
        ELSE
          IBCTX = IBCT(NSOLU+LUK+1,NB)
        ENDIF
!
!---    Zero flux boundary condition  ---
!
        IF( IBCTX.EQ.3 ) CYCLE
!
!---    Single boundary condition time  ---
!
        IF( IBCM(NB).EQ.1 ) THEN
!
!---      Solute transport  ---
!
          IF( NSL.LE.NSOLU ) THEN
            BCX(1) = BC(NSL+LBCU,1,MB)
            IF( IBCTX.EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---      Reactive species transport  ---
!
          ELSE
            BCX(1) = 0.D+0
            DO NSPX = 1,IBCSP(1,NB)
              NSP = IBCSP(NSPX+1,NB)
              MX = NSOLU+LBCU+NSPX
              BCX(NSPX+1) = BC(MX,1,MB)
!
!---          Aqueous species ---
!
              IF( NSP.LE.NSPL ) THEN
                IF( IBCT(NSOLU+LUK+1,NB).EQ.12 ) 
     &            BCX(NSPX+1) = SP_CBO(NB,NSP)
!
!---          Gas species ---
!
              ELSE
                IF( IBCT(NSOLU+LUK+2,NB).EQ.12 ) 
     &            BCX(NSPX+1) = SP_CBO(NB,NSP)
              ENDIF
            ENDDO
          ENDIF
!
!---    Multiple boundary condition times  ---
!
        ELSE
          IFIND = 0
          DO M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
              TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
              DTBC = MIN( BC(1,M,MB)-TMZ,DT )
              TFBC = (TMZ-5.D-1*DTBC-BC(1,M-1,MB))/TDBC
!
!---          Solute transport  ---
!
              IF( NSL.LE.NSOLU ) THEN
                BCX(1) = BC(NSL+LBCU,M-1,MB) +
     &            TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
                IF( IBCT(NSL+LUK,NB).EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---          Reactive species transport  ---
!
              ELSE
                BCX(1) = 0.D+0
                DO NSPX = 1,IBCSP(1,NB)
                  NSP = IBCSP(NSPX+1,NB)
                  MX = NSOLU+LBCU+NSPX
                  BCX(NSPX+1) = BC(MX,M-1,MB) +
     &              TFBC*(BC(MX,M,MB)-BC(MX,M-1,MB))
!
!---              Aqueous species ---
!
                  IF( NSP.LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.12 ) 
     &                BCX(NSPX+1) = SP_CBO(NB,NSP)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.12 ) 
     &                BCX(NSPX+1) = SP_CBO(NB,NSP)
                  ENDIF
                ENDDO
              ENDIF
              IFIND = 1
              EXIT
            ENDIF
          ENDDO
          IF( IFIND.EQ.0 ) CYCLE
        ENDIF
        N = IBCN(NB)
        N_DB = -NB
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---  Diffusion coefficients at node adjacent to boundary  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SMDLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
        DLP = TORL(2,N)*SL(2,N)*PORD(2,N)*SMDLP
        PCOR = (PG(2,N)+PATM)/PATM
        SMDGP = SMDG(NSL)*(TCOR**1.75)/PCOR
        DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SMDGP
!
!---  Phase fraction factors at node adjacent to boundary  ---
!
        XVLP = SL(2,N)*PORD(2,N)
        FCLP = 0.D+0
        IF( XVLP.GT.SMALL ) FCLP = YL(N,NSL)/XVLP
        XVGP = SG(2,N)*PORD(2,N)
        FCGP = 0.D+0
        IF( XVGP.GT.SMALL ) FCGP = YG(N,NSL)/XVGP
!
!---    Solute transport only, skip calculations for reactive
!       species transport  ---
!
        XVLB = SLB(2,NB)*PORDB(2,NB)
        XVGB = SGB(2,NB)*PORDB(2,NB)
        IF( NSL.LE.NSOLU ) THEN
!
!---      Phase fraction factors at boundary  ---
!
          IF( IPCL(NSL).EQ.2 ) THEN
            XVSB = RHOS(N)*PCSL(1,N,NSL)*(1.D+0-PORTB(2,NB))
     &        *SLB(2,NB)
          ELSE
            XVSB = RHOS(N)*PCSL(1,N,NSL)*(1.D+0-PORTB(2,NB))
          ENDIF
!
!---      Constant gas-aqueous partition coefficient  ---
!
          IF( IPCGL(NSL).EQ.0 ) THEN
            PCGLX = PCGL(1,NSL)
!
!---      Temperature dependent gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.1 ) THEN
            TK = TB(2,NB)+TABS
            PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &        + PCGL(3,NSL)*LOG(TK)
     &        + PCGL(4,NSL)*TK + PCGL(5,NSL)*TK**2 )
!
!---      Water-vapor 1 gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.2 ) THEN
            PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
          ENDIF
          PCGLX = MAX( PCGLX,1.D-20 )
          PCGLX = MIN( PCGLX,1.D+20 )
!
!---      Phase-volumetric concentration ratios  ---
!
          FCL = 1.D+0/(XVSB + XVLB + XVGB*PCGLX)
          FCG = 1.D+0/((XVSB + XVLB)/PCGLX + XVGB)
!
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = XVLB*FCL
          YGB(NB,NSL) = XVGB*FCG
!
!---      Convert boundary phase concentrations to
!         volumetric concentrations  ---
!
          IF( IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.14 .OR.
     &      IBCT(NSL+LUK,NB).EQ.23 ) THEN
            BCX(1) = BCX(1)/(FCL+SMALL)
          ELSEIF( IBCT(NSL+LUK,NB).EQ.9 .OR.
     &      IBCT(NSL+LUK,NB).EQ.15 .OR.
     &      IBCT(NSL+LUK,NB).EQ.43 ) THEN
            BCX(1) = BCX(1)/(FCG+SMALL)
          ENDIF
          CB(NB,NSL) = BCX(1)
        ELSE
!
!---      Convert species concentrations to total-component
!         concentrations  ---
!
          IF( NSL.LE.NSOLU+NEQC ) THEN
            NEQ = NSL-NSOLU
            YSPLX = 0.D+0
            YSPGX = 0.D+0
            DO NSP = 1,IEQ_C(1,NEQ)
              DO NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_C(NSP+1,NEQ) ) THEN
!
!---              Aqueous species ---
!
                  IF( IEQ_C(NSP+1,NEQ).LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+1,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPLX = YSPLX + EQ_C(NSP,NEQ)*BCX(NSPX+1)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+2,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPGX = YSPGX + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                  ENDIF                    
                  BCX(1) = BCX(1) + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
              ENDDO
            ENDDO
!
!---        Linked aqueous CO2   ---
!
            IF( ISPLK(6).EQ.NSL ) BCX(1) = 1.D+3*XLAB(2,NB)*
     &          RHOLB(2,N)*SLB(2,NB)*PORDB(2,NB)/WTMA
!
!---        Convert species concentrations to total-kinetic
!           concentrations  ---
!
          ELSEIF( NSL.LE.NSOLU+NEQC+NEQK ) THEN
            NEQ = NSL-NSOLU-NEQC
            YSPLX = 0.D+0
            DO NSP = 1,IEQ_K(1,NEQ)
              DO NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_K(NSP+1,NEQ) ) THEN
!
!---              Aqueous species ---
!
                  IF( IEQ_K(NSP+1,NEQ).LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+1,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPLX = YSPLX + EQ_K(NSP,NEQ)*BCX(NSPX+1)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+2,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                  ENDIF                    
                  BCX(1) = BCX(1) + EQ_K(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
              ENDDO
            ENDDO
!
!---        Linked aqueous CO2   ---
!
            IF( ISPLK(6).EQ.NSL ) BCX(1) = 1.D+3*XLAB(2,NB)*
     &          RHOLB(2,N)*SLB(2,NB)*PORDB(2,NB)/WTMA
          ENDIF
          IF( ABS(BCX(1))/EPSL.LT.EPSL ) THEN
            YSPLX = 0.D+0
          ELSE
            YSPLX = YSPLX/BCX(1)
          ENDIF
!
!---      Phase-volumetric concentration ratios  ---
!
          YLBX = MAX( MIN( 1.D+0,YSPLX ),0.D+0 )
          YGBX = 1.D+0-YLBX
          FCL = 0.D+0
          IF( XVLB/EPSL.GT.EPSL ) FCL = YLBX/XVLB
          FCG = 0.D+0
          IF( XVGB/EPSL.GT.EPSL ) FCG = YGBX/XVGB
!
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = YLBX
          YGB(NB,NSL) = YGBX
          CB(NB,NSL) = BCX(1)
        ENDIF
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULBX = 5.D-1*(UL(1,1,N)+UL(1,2,N))
            VLBX = 5.D-1*(VL(1,1,N)+VL(1,2,N))
            WLBX = WL(1,1,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WLBX,INDX)
              ULBX = ULBX/VMCX
              VLBX = VLBX/VMCX
              WLBX = WLBX/VMCX
            ENDIF
            ULBSQ = ULBX*ULBX
            VLBSQ = VLBX*VLBX
            WLBSQ = WLBX*WLBX
            ZLVB = SQRT(ULBSQ+VLBSQ+WLBSQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLB = (DISPL(N)*WLBSQ + DISPT(N)*(ULBSQ+VLBSQ))/
     &        (ZLVB+SMALL)
            UGBX = 5.D-1*(UG(1,1,N)+UG(1,2,N))
            VGBX = 5.D-1*(VG(1,1,N)+VG(1,2,N))
            WGBX = WG(1,1,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WGBX,INDX)
              UGBX = UGBX/VMCX
              VGBX = VGBX/VMCX
              WGBX = WGBX/VMCX
            ENDIF
            UGBSQ = UGBX*UGBX
            VGBSQ = VGBX*VGBX
            WGBSQ = WGBX*WGBX
            ZGVB = SQRT(UGBSQ+VGBSQ+WGBSQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGB = (DISPL(N)*WGBSQ + DISPT(N)*(UGBSQ+VGBSQ))/
     &        (ZGVB+SMALL)
          ELSE
            DPLB = 0.D+0
            DPGB = 0.D+0
          ENDIF
          FLB = AFZ(1,N)*WL(1,1,N)
          FGB = AFZ(1,N)*WG(1,1,N)
          CRLB = ABS( WL(1,1,N) )*DT/(DZGF(N)*XVLB+SMALL)
          CRGB = ABS( WG(1,1,N) )*DT/(DZGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLZ = DIFMN(DLB,DLP,DZGF(N),DZGF(N),WL(1,1,N),INDX)
            DLZ = AFZ(1,N)*(DLZ+DPLB)/(5.D-1*DZGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGZ = DIFMN(DGB,DGP,DZGF(N),DZGF(N),WG(1,1,N),INDX)
            DGZ = AFZ(1,N)*(DGZ+DPGB)/(5.D-1*DZGF(N))
            ALB = MAX( FLB,ZERO ) +
     &        DLZ*MAX((ONE-(TENTH*ABS(FLB)/(DLZ+SMALL)))**5,ZERO)
            AGB = MAX( FGB,ZERO ) +
     &        DGZ*MAX((ONE-(TENTH*ABS(FGB)/(DGZ+SMALL)))**5,ZERO)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLB = MIN( FLB,0.D+0 )
            FGB = MIN( FGB,0.D+0 )
            ALB = MAX( FLB,ZERO )
            AGB = MAX( FGB,ZERO )
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLB = MAX( FLB,0.D+0 )
            FGB = MAX( FGB,0.D+0 )
            ALB = MAX( FLB,ZERO )
            AGB = MAX( FGB,ZERO )
          ENDIF
          AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP
          AB = ALB*FCL + AGB*FCG
          WC(1,N,NSL) = WC(1,N,NSL) + (BCX(1)*AB - C(N,NSL)*AP)
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULSX = 5.D-1*(UL(1,1,N)+UL(1,2,N))
            VLSX = VL(1,1,N)
            WLSX = 5.D-1*(WL(1,1,N)+WL(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DYGF(N),DYGF(N),VLBX,INDX)
              ULSX = ULSX/VMCX
              VLSX = VLSX/VMCX
              WLSX = WLSX/VMCX
            ENDIF
            ULSSQ = ULSX*ULSX
            VLSSQ = VLSX*VLSX
            WLSSQ = WLSX*WLSX
            ZLVS = SQRT(ULSSQ+VLSSQ+WLSSQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLS = (DISPL(N)*VLSSQ + DISPT(N)*(ULSSQ+WLSSQ))/
     &        (ZLVS+SMALL)
            UGSX = 5.D-1*(UG(1,1,N)+UG(1,2,N))
            VGSX = VG(1,1,N)
            WGSX = 5.D-1*(WG(1,1,N)+WG(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DYGF(N),DYGF(N),VGSX,INDX)
              UGSX = UGSX/VMCX
              VGSX = VGSX/VMCX
              WGSX = WGSX/VMCX
            ENDIF
            UGSSQ = UGSX*UGSX
            VGSSQ = VGSX*VGSX
            WGSSQ = WGSX*WGSX
            ZGVS = SQRT(UGSSQ+VGSSQ+WGSSQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGS = (DISPL(N)*VGSSQ + DISPT(N)*(WGSSQ+VGSSQ))/
     &        (ZGVS+SMALL)
          ELSE
            DPLS = 0.D+0
            DPGS = 0.D+0
          ENDIF
          FLS = AFY(1,N)*VL(1,1,N)
          FGS = AFY(1,N)*VG(1,1,N)
          CRLS = ABS( VL(1,1,N) )*DT/(RP(N)*DYGF(N)*XVLB+SMALL)
          CRGS = ABS( VG(1,1,N) )*DT/(RP(N)*DYGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLY = DIFMN(DLB,DLP,DYGF(N),DYGF(N),VL(1,1,N),INDX)
            DLY = AFY(1,N)*(DLY+DPLS)/RP(N)/(5.D-1*DYGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGY = DIFMN(DGB,DGP,DYGF(N),DYGF(N),VG(1,1,N),INDX)
            DGY = AFY(1,N)*(DGY+DPGS)/RP(N)/(5.D-1*DYGF(N))
            ALS = MAX( FLS,ZERO ) +
     &        DLY*MAX((ONE-(TENTH*ABS(FLS)/(DLY+SMALL)))**5,ZERO)
            AGS = MAX( FGS,ZERO ) +
     &        DGY*MAX((ONE-(TENTH*ABS(FGS)/(DGY+SMALL)))**5,ZERO)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLS = MIN( FLS,0.D+0 )
            FGS = MIN( FGS,0.D+0 )
            ALS = MAX( FLS,ZERO )
            AGS = MAX( FGS,ZERO )
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLS = MAX( FLS,0.D+0 )
            FGS = MAX( FGS,0.D+0 )
            ALS = MAX( FLS,ZERO )
            AGS = MAX( FGS,ZERO )
          ENDIF
          AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP
          AS = ALS*FCL + AGS*FCG
          VC(1,N,NSL) = VC(1,N,NSL) + (BCX(1)*AS - C(N,NSL)*AP)
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULWX = UL(1,1,N)
            VLWX = 5.D-1*(VL(1,1,N)+VL(1,2,N))
            WLWX = 5.D-1*(WL(1,1,N)+WL(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DXGF(N),DXGF(N),ULBX,INDX)
              ULWX = ULWX/VMCX
              VLWX = VLWX/VMCX
              WLWX = WLWX/VMCX
            ENDIF
            ULWSQ = ULWX*ULWX
            VLWSQ = VLWX*VLWX
            WLWSQ = WLWX*WLWX
            ZLVW = SQRT(ULWSQ+VLWSQ+WLWSQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLW = (DISPL(N)*VLWSQ + DISPT(N)*(ULWSQ+WLWSQ))/
     &        (ZLVW+SMALL)
            UGWX = UG(1,1,N)
            VGWX = 5.D-1*(VG(1,1,N)+VG(1,2,N))
            WGWX = 5.D-1*(WG(1,1,N)+WG(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DXGF(N),DXGF(N),UGWX,INDX)
              UGWX = UGWX/VMCX
              VGWX = VGWX/VMCX
              WGWX = WGWX/VMCX
            ENDIF
            UGWSQ = UGWX*UGWX
            VGWSQ = VGWX*VGWX
            WGWSQ = WGWX*WGWX
            ZGVW = SQRT(UGWSQ+VGWSQ+WGWSQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGW = (DISPL(N)*VGWSQ + DISPT(N)*(WGWSQ+VGWSQ))/
     &        (ZGVW+SMALL)
          ELSE
            DPLW = 0.D+0
            DPGW = 0.D+0
          ENDIF
          FLW = AFX(1,N)*UL(1,1,N)
          FGW = AFX(1,N)*UG(1,1,N)
          CRLW = ABS( UL(1,1,N) )*DT/(DXGF(N)*XVLB+SMALL)
          CRGW = ABS( UG(1,1,N) )*DT/(DXGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLX = DIFMN(DLB,DLP,DXGF(N),DXGF(N),UL(1,1,N),INDX)
            DLX = AFX(1,N)*(DLX+DPLW)/(5.D-1*DXGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGX = DIFMN(DGB,DGP,DXGF(N),DXGF(N),UG(1,1,N),INDX)
            DGX = AFX(1,N)*(DGX+DPGW)/(5.D-1*DXGF(N))
            ALW = MAX(FLW,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLW)/(DLX+SMALL)))**5,ZERO)
            AGW = MAX(FGW,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGW)/(DGX+SMALL)))**5,ZERO)
            AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
            AW = ALW*FCL + AGW*FCG
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLW = MIN( FLW,0.D+0 )
            FGW = MIN( FGW,0.D+0 )
            ALW = MAX(FLW,ZERO)
            AGW = MAX(FGW,ZERO)
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLW = MAX( FLW,0.D+0 )
            FGW = MAX( FGW,0.D+0 )
              ALW = MAX(FLW,ZERO)
              AGW = MAX(FGW,ZERO)
          ENDIF
          AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
          AW = ALW*FCL + AGW*FCG
          UC(1,N,NSL) = UC(1,N,NSL) + (BCX(1)*AW - C(N,NSL)*AP)
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULEX = UL(1,2,N)
            VLEX = 5.D-1*(VL(1,1,N)+VL(1,2,N))
            WLEX = 5.D-1*(WL(1,1,N)+WL(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DXGF(N),DXGF(N),ULBX,INDX)
              ULEX = ULEX/VMCX
              VLEX = VLEX/VMCX
              WLEX = WLEX/VMCX
            ENDIF
            ULESQ = ULEX*ULEX
            VLESQ = VLEX*VLEX
            WLESQ = WLEX*WLEX
            ZLVE = SQRT(ULESQ+VLESQ+WLESQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLE = (DISPL(N)*VLESQ + DISPT(N)*(ULESQ+WLESQ))/
     &        (ZLVE+SMALL)
            UGEX = UG(1,2,N)
            VGEX = 5.D-1*(VG(1,1,N)+VG(1,2,N))
            WGEX = 5.D-1*(WG(1,1,N)+WG(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DXGF(N),DXGF(N),UGEX,INDX)
              UGEX = UGEX/VMCX
              VGEX = VGEX/VMCX
              WGEX = WGEX/VMCX
            ENDIF
            UGESQ = UGEX*UGEX
            VGESQ = VGEX*VGEX
            WGESQ = WGEX*WGEX
            ZGVE = SQRT(UGESQ+VGESQ+WGESQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGE = (DISPL(N)*VGESQ + DISPT(N)*(WGESQ+VGESQ))/
     &        (ZGVE+SMALL)
          ELSE
            DPLE = 0.D+0
            DPGE = 0.D+0
          ENDIF
          FLE = AFX(2,N)*UL(1,2,N)
          FGE = AFX(2,N)*UG(1,2,N)
          CRLE = ABS( UL(1,2,N) )*DT/(DXGF(N)*XVLB+SMALL)
          CRGE = ABS( UG(1,2,N) )*DT/(DXGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLX = DIFMN(DLP,DLB,DXGF(N),DXGF(N),UL(1,2,N),INDX)
            DLX = AFX(2,N)*(DLX+DPLE)/(5.D-1*DXGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGX = DIFMN(DGP,DGB,DXGF(N),DXGF(N),UG(1,2,N),INDX)
            DGX = AFX(2,N)*(DGX+DPGE)/(5.D-1*DXGF(N))
            ALE = MAX( -FLE,ZERO ) +
     &        DLX*MAX((ONE-(TENTH*ABS(FLE)/(DLX+SMALL)))**5,ZERO)
            AGE = MAX( -FGE,ZERO ) +
     &        DGX*MAX((ONE-(TENTH*ABS(FGE)/(DGX+SMALL)))**5,ZERO)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLE = MAX( FLE,0.D+0 )
            FGE = MAX( FGE,0.D+0 )
            ALE = MAX( -FLE,ZERO )
            AGE = MAX( -FGE,ZERO )
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLE = MIN( FLE,0.D+0 )
            FGE = MIN( FGE,0.D+0 )
            ALE = MAX( -FLE,ZERO )
            AGE = MAX( -FGE,ZERO )
          ENDIF
          AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP
          AE = ALE*FCL + AGE*FCG
          UC(2,N,NSL) = UC(2,N,NSL) + (C(N,NSL)*AP - BCX(1)*AE)
!
!---    North boundary
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULNX = 5.D-1*(UL(1,1,N)+UL(1,2,N))
            VLNX = VL(1,2,N)
            WLNX = 5.D-1*(WL(1,1,N)+WL(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DYGF(N),DYGF(N),VLBX,INDX)
              ULNX = ULNX/VMCX
              VLNX = VLNX/VMCX
              WLNX = WLNX/VMCX
            ENDIF
            ULNSQ = ULNX*ULNX
            VLNSQ = VLNX*VLNX
            WLNSQ = WLNX*WLNX
            ZLVN = SQRT(ULNSQ+VLNSQ+WLNSQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLN = (DISPL(N)*VLNSQ + DISPT(N)*(ULNSQ+WLNSQ))/
     &        (ZLVN+SMALL)
            UGNX = 5.D-1*(UG(1,1,N)+UG(1,2,N))
            VGNX = VG(1,2,N)
            WGNX = 5.D-1*(WG(1,1,N)+WG(1,2,N))
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DYGF(N),DYGF(N),VGNX,INDX)
              UGNX = UGNX/VMCX
              VGNX = VGNX/VMCX
              WGNX = WGNX/VMCX
            ENDIF
            UGNSQ = UGNX*UGNX
            VGNSQ = VGNX*VGNX
            WGNSQ = WGNX*WGNX
            ZGVN = SQRT(UGNSQ+VGNSQ+WGNSQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGN = (DISPL(N)*VGNSQ + DISPT(N)*(WGNSQ+VGNSQ))/
     &        (ZGVN+SMALL)
          ELSE
            DPLN = 0.D+0
            DPGN = 0.D+0
          ENDIF
          FLN = AFY(2,N)*VL(1,2,N)
          FGN = AFY(2,N)*VG(1,2,N)
          CRLN = ABS( VL(1,2,N) )*DT/(RP(N)*DYGF(N)*XVLB+SMALL)
          CRGN = ABS( VG(1,2,N) )*DT/(RP(N)*DYGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLY = DIFMN(DLP,DLB,DYGF(N),DYGF(N),VL(1,2,N),INDX)
            DLY = AFY(2,N)*(DLY+DPLN)/RP(N)/(5.D-1*DYGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGY = DIFMN(DGP,DGB,DYGF(N),DYGF(N),VG(1,2,N),INDX)
            DGY = AFY(2,N)*(DGY+DPGN)/RP(N)/(5.D-1*DYGF(N))
            ALN = MAX( -FLN,ZERO ) +
     &        DLY*MAX((ONE-(TENTH*ABS(FLN)/(DLY+SMALL)))**5,ZERO)
            AGN = MAX( -FGN,ZERO ) +
     &        DGY*MAX((ONE-(TENTH*ABS(FGN)/(DGY+SMALL)))**5,ZERO)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLN = MAX( FLN,0.D+0 )
            FGN = MAX( FGN,0.D+0 )
            ALN = MAX( -FLN,ZERO )
            AGN = MAX( -FGN,ZERO )
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLN = MIN( FLN,0.D+0 )
            FGN = MIN( FGN,0.D+0 )
            ALN = MAX( -FLN,ZERO )
            AGN = MAX( -FGN,ZERO )
          ENDIF
          AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP
          AN = ALN*FCL + AGN*FCG
          VC(2,N,NSL) = VC(2,N,NSL) + (C(N,NSL)*AP - BCX(1)*AN)
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            ULTX = 5.D-1*(UL(1,1,N)+UL(1,2,N))
            VLTX = 5.D-1*(VL(1,1,N)+VL(1,2,N))
            WLTX = WL(1,2,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SLB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WLTX,INDX)
              ULTX = ULTX/VMCX
              VLTX = VLTX/VMCX
              WLTX = WLTX/VMCX
            ENDIF
            ULTSQ = ULTX*ULTX
            VLTSQ = VLTX*VLTX
            WLTSQ = WLTX*WLTX
            ZLVT = SQRT(ULTSQ+VLTSQ+WLTSQ)
!
!---        Aqueous hydraulic dispersion coefficient  ---
!
            DPLT = (DISPL(N)*WLTSQ + DISPT(N)*(ULTSQ+VLTSQ))/
     &        (ZLVT+SMALL)
            UGTX = 5.D-1*(UG(1,1,N)+UG(1,2,N))
            VGTX = 5.D-1*(VG(1,1,N)+VG(1,2,N))
            WGTX = WG(1,2,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SGB(2,NB)*PORDB(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WGTX,INDX)
              UGTX = UGTX/VMCX
              VGTX = VGTX/VMCX
              WGTX = WGTX/VMCX
            ENDIF
            UGTSQ = UGTX*UGTX
            VGTSQ = VGTX*VGTX
            WGTSQ = WGTX*WGTX
            ZGVT = SQRT(UGTSQ+VGTSQ+WGTSQ)
!
!---        Gas hydraulic dispersion coefficient  ---
!
            DPGT = (DISPL(N)*WGTSQ + DISPT(N)*(UGTSQ+VGTSQ))/
     &        (ZGVT+SMALL)
          ELSE
            DPLT = 0.D+0
            DPGT = 0.D+0
          ENDIF
          FLT = AFZ(2,N)*WL(1,2,N)
          FGT = AFZ(2,N)*WG(1,2,N)
          CRLT = ABS( WL(1,2,N) )*DT/(DZGF(N)*XVLB+SMALL)
          CRGT = ABS( WG(1,2,N) )*DT/(DZGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLZ = DIFMN(DLP,DLB,DZGF(N),DZGF(N),WL(1,2,N),INDX)
            DLZ = AFZ(2,N)*(DLZ+DPLT)/(5.D-1*DZGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGZ = DIFMN(DGP,DGB,DZGF(N),DZGF(N),WG(1,2,N),INDX)
            DGZ = AFZ(2,N)*(DGZ+DPGT)/(5.D-1*DZGF(N))
            ALT = MAX( -FLT,ZERO ) +
     &        DLZ*MAX((ONE-(TENTH*ABS(FLT)/(DLZ+SMALL)))**5,ZERO)
            AGT = MAX( -FGT,ZERO ) +
     &        DGZ*MAX((ONE-(TENTH*ABS(FGT)/(DGZ+SMALL)))**5,ZERO)
!
!---      Outflow  ---
!
          ELSEIF( IBCTX.EQ.7 .OR. IBCTX.EQ.19 .OR. 
     &      IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLT = MAX( FLT,0.D+0 )
            FGT = MAX( FGT,0.D+0 )
            ALT = MAX( -FLT,ZERO )
            AGT = MAX( -FGT,ZERO )
!
!---      Inflow ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43 ) THEN
            FLT = MIN( FLT,0.D+0 )
            FGT = MIN( FGT,0.D+0 )
            ALT = MAX( -FLT,ZERO )
            AGT = MAX( -FGT,ZERO )
          ENDIF
          AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP
          AT = ALT*FCL + AGT*FCG
          WC(2,N,NSL) = WC(2,N,NSL) + (C(N,NSL)*AP - BCX(1)*AT)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFXB_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXL( NSL )
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
!     Compute solute transport flux aqueous-phase, excluding boundaries,
!     using a Patankar scheme.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 January 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
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
      REAL*8 KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFXL'
!
!---  Aqueous solute flux, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          TCOR = (T(2,N)+TABS)/TSPRF
          SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
          VMCP = SL(2,N)*PORD(2,N)
          FCLP = YL(N,NSL)/(VMCP+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFP = TORL(2,N)*VMCP*SDFLP
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFP = SDCL(1,N,NSL)*SDCL(2,N,NSL)*
     &        EXP(VMCP*SDCL(3,N,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFP = TORL(2,N)*VMCP*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFP = SDCL(1,N,NSL)*SDCL(2,N,NSL)*
     &        VMCP**SDCL(3,N,NSL)
          ENDIF
          TCOR = (T(2,NW)+TABS)/TSPRF
          SDFLW = SMDL(NSL)*TCOR*(VISRL/VISL(2,NW))
          VMCW = SL(2,NW)*PORD(2,NW)
          FCLW = YL(NW,NSL)/(VMCW+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFW = TORL(2,NW)*VMCW*SDFLW
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFW = SDCL(1,NW,NSL)*SDCL(2,NW,NSL)*
     &        EXP(VMCW*SDCL(3,NW,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFW = TORL(2,NW)*VMCW*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFW = SDCL(1,NW,NSL)*SDCL(2,NW,NSL)*
     &        VMCW**SDCL(3,NW,NSL)
          ENDIF
          IF( IDISP.EQ.1 ) THEN
            ULWX = UL(1,1,N)
            VLWX = 0.25D+0*((VL(1,1,N)+VL(1,2,N))*DXGF(NW) + 
     &        (VL(1,1,NW)+VL(1,2,NW))*DXGF(N))/DXGP(1,N)
            WLWX = 0.25D+0*((WL(1,1,N)+WL(1,2,N))*DXGF(NW) + 
     &        (WL(1,1,NW)+WL(1,2,NW))*DXGF(N))/DXGP(1,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCW = SL(2,NW)*PORD(2,NW)
              INDX = 17
              VMCX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),ULWX,INDX)
              ULWX = ULWX/VMCX
              VLWX = VLWX/VMCX
              WLWX = WLWX/VMCX
            ENDIF
            ULWSQ = ULWX*ULWX
            VLWSQ = VLWX*VLWX
            WLWSQ = WLWX*WLWX
            ZVW = SQRT(ULWSQ+VLWSQ+WLWSQ)
            INDX = 17
            DPLW = DIFMN(DISPL(NW),DISPL(N),DXGF(NW),DXGF(N),ULWX,INDX)
            DPTW = DIFMN(DISPT(NW),DISPT(N),DXGF(NW),DXGF(N),ULWX,INDX)
            DPW = (DPLW*ULWSQ + DPTW*(VLWSQ+WLWSQ))/(ZVW+SMALL)
          ELSE
            DPW = 0.D+0
          ENDIF
          INDX = 16
          DFW = DIFMN(DFW,DFP,DXGF(NW),DXGF(N),UL(1,1,N),INDX)
          DDW = (DFW+DPW)/DXGP(1,N)
          IF( ISLC(1).GE.1 )  THEN
            UC(1,N,NSL) = UC(1,N,NSL) + DDW*(C(NW,NSL)*FCLW - 
     &        C(N,NSL)*FCLP)
          ELSE
           AL = MAX( UL(1,1,N),ZERO ) +
     &       DDW*MAX( (ONE-(TENTH*ABS(UL(1,1,N))/(DDW+SMALL)))**5,ZERO )
           ALP = MAX( -UL(1,1,N),ZERO ) +
     &       DDW*MAX( (ONE-(TENTH*ABS(UL(1,1,N))/(DDW+SMALL)))**5,ZERO )
           UC(1,N,NSL) = UC(1,N,NSL) + (C(NW,NSL)*AL*FCLW - 
     &       C(N,NSL)*ALP*FCLP)
          ENDIF
          UC(2,NW,NSL) = UC(1,N,NSL)
        ELSE
          UC(1,N,NSL) = 0.D+0
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          TCOR = (T(2,N)+TABS)/TSPRF
          SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
          VMCP = SL(2,N)*PORD(2,N)
          FCLP = YL(N,NSL)/(VMCP+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFP = TORL(2,N)*VMCP*SDFLP
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFP = SDCL(1,N,NSL)*SDCL(2,N,NSL)*
     &        EXP(VMCP*SDCL(3,N,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFP = TORL(2,N)*VMCP*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFP = SDCL(1,N,NSL)*SDCL(2,N,NSL)*
     &        VMCP**SDCL(3,N,NSL)
          ENDIF
          TCOR = (T(2,NS)+TABS)/TSPRF
          SDFLS = SMDL(NSL)*TCOR*(VISRL/VISL(2,NS))
          VMCS = SL(2,NS)*PORD(2,NS)
          FCLS = YL(NS,NSL)/(VMCS+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFS = TORL(2,NS)*VMCS*SDFLS
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFS = SDCL(1,NS,NSL)*SDCL(2,NS,NSL)*
     &        EXP(VMCS*SDCL(3,NS,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFS = TORL(2,NS)*VMCS*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFS = SDCL(1,NS,NSL)*SDCL(2,NS,NSL)*
     &        VMCS**SDCL(3,NS,NSL)
          ENDIF
          IF( IDISP.EQ.1 ) THEN
            ULSX = 0.25D+0*((UL(1,1,N)+UL(1,2,N))*DYGF(NS) +
     &        (UL(1,1,NS)+UL(1,2,NS))*DYGF(N))/DYGP(1,N)
            VLSX = VL(1,1,N)
            WLSX = 0.25D+0*((WL(1,1,N)+WL(1,2,N))*DYGF(NS) + 
     &        (WL(1,1,NS)+WL(1,2,NS))*DYGF(N))/DYGP(1,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCS = SL(2,NS)*PORD(2,NS)
              INDX = 17
              VMCX = DIFMN(VMCS,VMCP,DYGF(NS),DYGF(N),VLSX,INDX)
              ULSX = ULSX/VMCX
              VLSX = VLSX/VMCX
              WLSX = WLSX/VMCX
            ENDIF
            ULSSQ = ULSX*ULSX
            VLSSQ = VLSX*VLSX
            WLSSQ = WLSX*WLSX
            ZVS = SQRT(ULSSQ+VLSSQ+WLSSQ)
            INDX = 17
            DPLS = DIFMN(DISPL(NS),DISPL(N),DYGF(NS),DYGF(N),VLSX,INDX)
            DPTS = DIFMN(DISPT(NS),DISPT(N),DYGF(NS),DYGF(N),VLSX,INDX)
            DPS = (DPLS*VLSSQ + DPTS*(ULSSQ+WLSSQ))/(ZVS+SMALL)
          ELSE
            DPS = 0.D+0
          ENDIF
          INDX = 16
          DFS = DIFMN(DFS,DFP,DYGF(NS),DYGF(N),VL(1,1,N),INDX)
          DDS = (DFS+DPS)/(DYGP(1,N)*RP(N))
          IF( ISLC(1).GE.1 )  THEN
            VC(1,N,NSL) = VC(1,N,NSL) + DDS*(C(NS,NSL)*FCLS - 
     &        C(N,NSL)*FCLP)
          ELSE
           AL = MAX( VL(1,1,N),ZERO ) +
     &       DDS*MAX( (ONE-(TENTH*ABS(VL(1,1,N))/(DDS+SMALL)))**5,ZERO )
           ALP = MAX( -VL(1,1,N),ZERO ) +
     &       DDS*MAX( (ONE-(TENTH*ABS(VL(1,1,N))/(DDS+SMALL)))**5,ZERO )
           VC(1,N,NSL) = VC(1,N,NSL) + (C(NS,NSL)*AL*FCLS - 
     &       C(N,NSL)*ALP*FCLP)
          ENDIF
          VC(2,NS,NSL) = VC(1,N,NSL)
        ELSE
          VC(1,N,NSL) = 0.D+0
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          TCOR = (T(2,N)+TABS)/TSPRF
          SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
          VMCP = SL(2,N)*PORD(2,N)
          FCLP = YL(N,NSL)/(VMCP+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFP = TORL(2,N)*VMCP*SDFLP
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFP = SDCL(1,N,NSL)*SDCL(2,N,NSL)*
     &        EXP(VMCP*SDCL(3,N,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFP = TORL(2,N)*VMCP*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFP = SDCL(1,N,NSL)*SDCL(2,N,NSL)*
     &        VMCP**SDCL(3,N,NSL)
          ENDIF
          TCOR = (T(2,NB)+TABS)/TSPRF
          SDFLB = SMDL(NSL)*TCOR*(VISRL/VISL(2,NB))
          VMCB = SL(2,NB)*PORD(2,NB)
          FCLB = YL(NB,NSL)/(VMCB+SMALL)
          IF( IEDL(NSL).EQ.1 ) THEN
            DFB = TORL(2,NB)*VMCB*SDFLB
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DFB = SDCL(1,NB,NSL)*SDCL(2,NB,NSL)*
     &        EXP(VMCB*SDCL(3,NB,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DFB = TORL(2,NB)*VMCB*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DFB = SDCL(1,NB,NSL)*SDCL(2,NB,NSL)*
     &        VMCB**SDCL(3,NB,NSL)
          ENDIF
          IF( IDISP .EQ. 1 ) THEN
            ULBX = 0.25D+0*((UL(1,1,N)+UL(1,2,N))*DZGF(NB) +
     &        (UL(1,1,NB)+UL(1,2,NB))*DZGF(N))/DZGP(1,N)
            VLBX = 0.25D+0*((VL(1,1,N)+VL(1,2,N))*DZGF(NB) + 
     &        (VL(1,1,NB)+VL(1,2,NB))*DZGF(N))/DZGP(1,N)
            WLBX = WL(1,1,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SL(2,N)*PORD(2,N)
              VMCB = SL(2,NB)*PORD(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WLBX,INDX)
              ULBX = ULBX/VMCX
              VLBX = VLBX/VMCX
              WLBX = WLBX/VMCX
            ENDIF
            ULBSQ = ULBX*ULBX
            VLBSQ = VLBX*VLBX
            WLBSQ = WLBX*WLBX
            ZVB = SQRT(ULBSQ+VLBSQ+WLBSQ)
            INDX = 17
            DPLB = DIFMN(DISPL(NB),DISPL(N),DZGF(NB),DZGF(N),WLBX,INDX)
            DPTB = DIFMN(DISPT(NB),DISPT(N),DZGF(NB),DZGF(N),WLBX,INDX)
            DPB = (DPLB*WLBSQ + DPTB*(VLBSQ+ULBSQ))/(ZVB+SMALL)
          ELSE
            DPB = 0.D+0
          ENDIF
          INDX = 16
          DFB = DIFMN(DFB,DFP,DZGF(NB),DZGF(N),WL(1,1,N),INDX)
          DDB = (DFB+DPB)/DZGP(1,N)
          IF( ISLC(1).GE.1 ) THEN
            WC(1,N,NSL) = WC(1,N,NSL) + DDB*(C(NB,NSL)*FCLB - 
     &        C(N,NSL)*FCLP)
          ELSE
           AL = MAX( WL(1,1,N),ZERO ) +
     &       DDB*MAX( (ONE-(TENTH*ABS(WL(1,1,N))/(DDB+SMALL)))**5,ZERO )
           ALP = MAX( -WL(1,1,N),ZERO ) +
     &       DDB*MAX( (ONE-(TENTH*ABS(WL(1,1,N))/(DDB+SMALL)))**5,ZERO )
           WC(1,N,NSL) = WC(1,N,NSL) + (C(NB,NSL)*AL*FCLB - 
     &       C(N,NSL)*ALP*FCLP)
          ENDIF
          WC(2,NB,NSL) = WC(1,N,NSL)
        ELSE
          WC(1,N,NSL) = 0.D+0
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFXL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXG( NSL )
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
!     Compute solute transport flux gas-phase, excluding boundaries,
!     using a Patankar scheme.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 January 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE HYST
      USE GRID
      USE FLUX
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
      REAL*8 KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFXG'
!
!---  Aqueous solute flux, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          TCOR = (T(2,N)+TABS)/TSPRF
          PCOR = (PG(2,N)+PATM)/PATM
          SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCP = SG(2,N)*PORD(2,N)
          DFP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
          FCGP = YG(N,NSL)/(VMCP+SMALL)
          TCOR = (T(2,NW)+TABS)/TSPRF
          PCOR = (PG(2,NW)+PATM)/PATM
          SDFGW = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCW = SG(2,NW)*PORD(2,NW)
          DFW = TORG(2,NW)*(SG(2,NW)-SGT(2,NW))*PORD(2,NW)*SDFGW
          FCGW = YG(NW,NSL)/(VMCW+SMALL)
          IF( IDISP.EQ.1 ) THEN
            UGWX = UG(1,1,N)
            VGWX = 0.25D+0*((VG(1,1,N)+VG(1,2,N))*DXGF(NW) + 
     &        (VG(1,1,NW)+VG(1,2,NW))*DXGF(N))/DXGP(1,N)
            WGWX = 0.25D+0*((WG(1,1,N)+WG(1,2,N))*DXGF(NW) + 
     &        (WG(1,1,NW)+WG(1,2,NW))*DXGF(N))/DXGP(1,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCW = SG(2,NW)*PORD(2,NW)
              INDX = 17
              VMCX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UGWX,INDX)
              UGWX = UGWX/VMCX
              VGWX = VGWX/VMCX
              WGWX = WGWX/VMCX
            ENDIF
            UGWSQ = UGWX*UGWX
            VGWSQ = VGWX*VGWX
            WGWSQ = WGWX*WGWX
            ZVW = SQRT(UGWSQ+VGWSQ+WGWSQ)
            INDX = 17
            DPLW = DIFMN(DISPL(NW),DISPL(N),DXGF(NW),DXGF(N),UGWX,INDX)
            DPTW = DIFMN(DISPT(NW),DISPT(N),DXGF(NW),DXGF(N),UGWX,INDX)
            DPW = (DPLW*UGWSQ + DPTW*(VGWSQ+WGWSQ))/(ZVW+SMALL)
          ELSE
            DPW = 0.D+0
          ENDIF
          INDX = 16
          DFW = DIFMN(DFW,DFP,DXGF(NW),DXGF(N),UG(1,1,N),INDX)
          DDW = (DFW+DPW)/DXGP(1,N)
          IF( ISLC(1).GE.1 )  THEN
            UC(1,N,NSL) = UC(1,N,NSL) + DDW*(C(NW,NSL)*FCGW - 
     &        C(N,NSL)*FCGP)
          ELSE
           AG = MAX( UG(1,1,N),ZERO ) +
     &       DDW*MAX( (ONE-(TENTH*ABS(UG(1,1,N))/(DDW+SMALL)))**5,ZERO )
           AGP = MAX( -UG(1,1,N),ZERO ) +
     &       DDW*MAX( (ONE-(TENTH*ABS(UG(1,1,N))/(DDW+SMALL)))**5,ZERO )
           UC(1,N,NSL) = UC(1,N,NSL) + (C(NW,NSL)*AG*FCGW - 
     &       C(N,NSL)*AGP*FCGP)
          ENDIF
          UC(2,NW,NSL) = UC(1,N,NSL)
        ELSE
          UC(1,N,NSL) = 0.D+0
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          TCOR = (T(2,N)+TABS)/TSPRF
          PCOR = (PG(2,N)+PATM)/PATM
          SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCP = SG(2,N)*PORD(2,N)
          DFP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
          FCGP = YG(N,NSL)/(VMCP+SMALL)
          TCOR = (T(2,NS)+TABS)/TSPRF
          PCOR = (PG(2,NS)+PATM)/PATM
          SDFGS = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCS = SG(2,NS)*PORD(2,NS)
          DFS = TORG(2,NS)*(SG(2,NS)-SGT(2,NS))*PORD(2,NS)*SDFGS
          FCGS = YG(NS,NSL)/(VMCS+SMALL)
          IF( IDISP.EQ.1 ) THEN
            UGSX = 0.25D+0*((UG(1,1,N)+UG(1,2,N))*DYGF(NS) +
     &        (UG(1,1,NS)+UG(1,2,NS))*DYGF(N))/DYGP(1,N)
            VGSX = VG(1,1,N)
            WGSX = 0.25D+0*((WG(1,1,N)+WG(1,2,N))*DYGF(NS) + 
     &        (WG(1,1,NS)+WG(1,2,NS))*DYGF(N))/DYGP(1,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCS = SG(2,NS)*PORD(2,NS)
              INDX = 17
              VMCX = DIFMN(VMCS,VMCP,DYGF(NS),DYGF(N),VGSX,INDX)
              UGSX = UGSX/VMCX
              VGSX = VGSX/VMCX
              WGSX = WGSX/VMCX
            ENDIF
            UGSSQ = UGSX*UGSX
            VGSSQ = VGSX*VGSX
            WGSSQ = WGSX*WGSX
            ZVS = SQRT(UGSSQ+VGSSQ+WGSSQ)
            INDX = 17
            DPLS = DIFMN(DISPL(NS),DISPL(N),DYGF(NS),DYGF(N),VGSX,INDX)
            DPTS = DIFMN(DISPT(NS),DISPT(N),DYGF(NS),DYGF(N),VGSX,INDX)
            DPS = (DPLS*VGSSQ + DPTS*(UGSSQ+WGSSQ))/(ZVS+SMALL)
          ELSE
            DPS = 0.D+0
          ENDIF
          INDX = 16
          DFS = DIFMN(DFS,DFP,DYGF(NS),DYGF(N),VG(1,1,N),INDX)
          DDS = (DFS+DPS)/(DYGP(1,N)*RP(N))
          IF( ISLC(1).GE.1 )  THEN
            VC(1,N,NSL) = VC(1,N,NSL) + DDS*(C(NS,NSL)*FCGS - 
     &        C(N,NSL)*FCGP)
          ELSE
           AG = MAX( VG(1,1,N),ZERO ) +
     &       DDS*MAX( (ONE-(TENTH*ABS(VG(1,1,N))/(DDS+SMALL)))**5,ZERO )
           AGP = MAX( -VG(1,1,N),ZERO ) +
     &       DDS*MAX( (ONE-(TENTH*ABS(VG(1,1,N))/(DDS+SMALL)))**5,ZERO )
           VC(1,N,NSL) = VC(1,N,NSL) + (C(NS,NSL)*AG*FCGS - 
     &       C(N,NSL)*AGP*FCGP)
          ENDIF
          VC(2,NS,NSL) = VC(1,N,NSL)
        ELSE
          VC(1,N,NSL) = 0.D+0
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          TCOR = (T(2,N)+TABS)/TSPRF
          PCOR = (PG(2,N)+PATM)/PATM
          SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCP = SG(2,N)*PORD(2,N)
          DFP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
          FCGP = YG(N,NSL)/(VMCP+SMALL)
          TCOR = (T(2,NB)+TABS)/TSPRF
          PCOR = (PG(2,NB)+PATM)/PATM
          SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCB = SG(2,NB)*PORD(2,NB)
          DFB = TORG(2,NB)*(SG(2,NB)-SGT(2,NB))*PORD(2,NB)*SDFGB
          FCGB = YG(NB,NSL)/(VMCB+SMALL)
          IF( IDISP .EQ. 1 ) THEN
            UGBX = 0.25D+0*((UG(1,1,N)+UG(1,2,N))*DZGF(NB) +
     &        (UG(1,1,NB)+UG(1,2,NB))*DZGF(N))/DZGP(1,N)
            VGBX = 0.25D+0*((VG(1,1,N)+VG(1,2,N))*DZGF(NB)
     &      + (VG(1,1,NB)+VG(1,2,NB))*DZGF(N))/DZGP(1,N)
            WGBX = WG(1,1,N)
            IF( ISLC(52).EQ.0 ) THEN
              VMCP = SG(2,N)*PORD(2,N)
              VMCB = SG(2,NB)*PORD(2,NB)
              INDX = 17
              VMCX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WGBX,INDX)
              UGBX = UGBX/VMCX
              VGBX = VGBX/VMCX
              WGBX = WGBX/VMCX
            ENDIF
            UGBSQ = UGBX*UGBX
            VGBSQ = VGBX*VGBX
            WGBSQ = WGBX*WGBX
            ZVB = SQRT(UGBSQ+VGBSQ+WGBSQ)
            INDX = 17
            DPLB = DIFMN(DISPL(NB),DISPL(N),DZGF(NB),DZGF(N),WGBX,INDX)
            DPTB = DIFMN(DISPT(NB),DISPT(N),DZGF(NB),DZGF(N),WGBX,INDX)
            DPB = (DPLB*WGBSQ + DPTB*(VGBSQ+UGBSQ))/(ZVB+SMALL)
          ELSE
            DPB = 0.D+0
          ENDIF
          INDX = 16
          DFB = DIFMN(DFB,DFP,DZGF(NB),DZGF(N),WG(1,1,N),INDX)
          DDB = (DFB+DPB)/DZGP(1,N)
          IF( ISLC(1).GE.1 ) THEN
            WC(1,N,NSL) = WC(1,N,NSL) + DDB*(C(NB,NSL)*FCGB - 
     &        C(N,NSL)*FCGP)
          ELSE
           AG = MAX( WG(1,1,N),ZERO ) +
     &       DDB*MAX( (ONE-(TENTH*ABS(WG(1,1,N))/(DDB+SMALL)))**5,ZERO )
           AGP = MAX( -WG(1,1,N),ZERO ) +
     &       DDB*MAX( (ONE-(TENTH*ABS(WG(1,1,N))/(DDB+SMALL)))**5,ZERO )
           WC(1,N,NSL) = WC(1,N,NSL) + (C(NB,NSL)*AG*FCGB - 
     &       C(N,NSL)*AGP*FCGP)
          ENDIF
          WC(2,NB,NSL) = WC(1,N,NSL)
        ELSE
          WC(1,N,NSL) = 0.D+0
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFXG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXZ( NSL )
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
!     STOMPX-CO2
!
!     Zero solute transport flux.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE GRID
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/SFXZ'
!
!---  Loop over local nodes  ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
        DO M = 1,2
          UC(M,N,NSL) = 0.D+0
          VC(M,N,NSL) = 0.D+0
          WC(M,N,NSL) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFXZ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SHDPG( N,DPB,DPS,DPW,DPE,DPN,DPT,NSL )
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
!     STOMPX-CO2
!
!     Calculates hydrodynamic dispersion coefficients for the gas
!     phase from phase velocities and user-specified dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/SHDPG'
!
!---  Initialize dispersion coefficients  ---
!
      M = 1
      DPB = 0.D+0
      DPS = 0.D+0
      DPW = 0.D+0
      DPE = 0.D+0
      DPN = 0.D+0
      DPT = 0.D+0
!
!---  Bottom face  ---
!
      NB = ICM(1,N)
      IF( NB.NE.0 ) THEN
        UGBX = 0.25D+0*((UG(1,1,N)+UG(1,2,N))*DZGF(NB) +
     &    (UG(1,1,NB)+UG(1,2,NB))*DZGF(N))/DZGP(1,N)
        VGBX = 0.25D+0*((VG(1,1,N)+VG(1,2,N))*DZGF(NB)
     &  + (VG(1,1,NB)+VG(1,2,NB))*DZGF(N))/DZGP(1,N)
        WGBX = WG(1,1,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SG(2,N)*PORD(2,N)
          VMCB = SG(2,NB)*PORD(2,NB)
          INDX = 17
          VMCX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WGBX,INDX)
          UGBX = UGBX/VMCX
          VGBX = VGBX/VMCX
          WGBX = WGBX/VMCX
        ENDIF
        UGBSQ = UGBX*UGBX
        VGBSQ = VGBX*VGBX
        WGBSQ = WGBX*WGBX
        ZVB = SQRT(UGBSQ+VGBSQ+WGBSQ)
        INDX = 17
        DPLB = DIFMN(DISPL(NB),DISPL(N),DZGF(NB),DZGF(N),WGBX,INDX)
        DPTB = DIFMN(DISPT(NB),DISPT(N),DZGF(NB),DZGF(N),WGBX,INDX)
        DPB = (DPLB*WGBSQ + DPTB*(VGBSQ+UGBSQ))/(ZVB+SMALL)
      ENDIF
!
!---  South face  ---
!
      NS = ICM(2,N)
      IF( NS.NE.0 ) THEN
        UGSX = 0.25D+0*((UG(1,1,N)+UG(1,2,N))*DYGF(NS) +
     &    (UG(1,1,NS)+UG(1,2,NS))*DYGF(N))/DYGP(1,N)
        VGSX = VG(1,1,N)
        WGSX = 0.25D+0*((WG(1,1,N)+WG(1,2,N))*DYGF(NS) + 
     &    (WG(1,1,NS)+WG(1,2,NS))*DYGF(N))/DYGP(1,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SG(2,N)*PORD(2,N)
          VMCS = SG(2,NS)*PORD(2,NS)
          INDX = 17
          VMCX = DIFMN(VMCS,VMCP,DYGF(NS),DYGF(N),VGSX,INDX)
          UGSX = UGSX/VMCX
          VGSX = VGSX/VMCX
          WGSX = WGSX/VMCX
        ENDIF
        UGSSQ = UGSX*UGSX
        VGSSQ = VGSX*VGSX
        WGSSQ = WGSX*WGSX
        ZVS = SQRT(UGSSQ+VGSSQ+WGSSQ)
        INDX = 17
        DPLS = DIFMN(DISPL(NS),DISPL(N),DYGF(NS),DYGF(N),VGSX,INDX)
        DPTS = DIFMN (DISPT(NS),DISPT(N),DYGF(NS),DYGF(N),VGSX,INDX)
        DPS = (DPLS*VGSSQ + DPTS*(UGSSQ+WGSSQ))/(ZVS+SMALL)
      ENDIF
!
!---  West face  ---
!
      NW = ICM(3,N)
      IF( NW.NE.0 ) THEN
        UGWX = UG(M,1,N)
        VGWX = 0.25D+0*((VG(1,1,N)+VG(1,2,N))*DXGF(NW) + 
     &    (VG(1,1,NW)+VG(1,2,NW))*DXGF(N))/DXGP(1,N)
        WGWX = 0.25D+0*((WG(1,1,N)+WG(1,2,N))*DXGF(NW) + 
     &    (WG(1,1,NW)+WG(1,2,NW))*DXGF(N))/DXGP(1,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SG(2,N)*PORD(2,N)
          VMCW = SG(2,NW)*PORD(2,NW)
          INDX = 17
          VMCX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UGWX,INDX)
          UGWX = UGWX/VMCX
          VGWX = VGWX/VMCX
          WGWX = WGWX/VMCX
        ENDIF
        UGWSQ = UGWX*UGWX
        VGWSQ = VGWX*VGWX
        WGWSQ = WGWX*WGWX
        ZVW = SQRT(UGWSQ+VGWSQ+WGWSQ)
        INDX = 17
        DPLW = DIFMN(DISPL(NW),DISPL(N),DXGF(NW),DXGF(N),UGX,INDX)
        DPTW = DIFMN(DISPT(NW),DISPT(N),DXGF(NW),DXGF(N),UGX,INDX)
        DPW = (DPLW*UGWSQ + DPTW*(VGWSQ+WGWSQ))/(ZVW+SMALL)
      ENDIF
!
!---  East face  ---
!
      NE = ICM(4,N)
      IF( NE.NE.0 ) THEN
        UGEX = UG(M,2,N)
        VGEX = 0.25D+0*((VG(1,1,N)+VG(1,2,N))*DXGF(NE) + 
     &    (VG(1,1,NE)+VG(1,2,NE))*DXGF(N))/DXGP(2,N)
        WGEX = 0.25D+0*((WG(1,1,N)+WG(1,2,N))*DXGF(NE) + 
     &    (WG(1,1,NE)+WG(1,2,NE))*DXGF(N))/DXGP(2,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SG(2,N)*PORD(2,N)
          VMCE = SG(2,NE)*PORD(2,NE)
          INDX = 17
          VMCX = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NE),UGEX,INDX)
          UGEX = UGEX/VMCX
          VGEX = VGEX/VMCX
          WGEX = WGEX/VMCX
        ENDIF
        UGESQ = UGEX*UGEX
        VGESQ = VGEX*VGEX
        WGESQ = WGEX*WGEX
        ZVE = SQRT(UGESQ+VGESQ+WGESQ)
        INDX = 17
        DPLE = DIFMN(DISPL(N),DISPL(NE),DXGF(N),DXGF(NE),UGEX,INDX)
        DPTE = DIFMN(DISPT(N),DISPT(NE),DXGF(N),DXGF(NE),UGEX,INDX)
        DPE = (DPLE*UGESQ + DPTE*(VGESQ+WGESQ))/(ZVE+SMALL)
      ENDIF
!
!---  North face  ---
!
      NN = ICM(5,N)
      IF( NN.NE.0 ) THEN
        UGNX = 0.25D+0*((UG(1,1,N)+UG(1,2,N))*DYGF(NN) +
     &    (UG(1,1,NN)+UG(1,2,NN))*DYGF(N))/DYGP(2,N)
        VGNX = VG(M,2,N)
        WGNX = 0.25D+0*((WG(1,1,N)+WG(1,2,N))*DYGF(NN) + 
     &    (WG(1,1,NN)+WG(1,2,NN))*DYGF(N))/DYGP(2,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SG(2,N)*PORD(2,N)
          VMCN = SG(2,NN)*PORD(2,NN)
          INDX = 17
          VMCX = DIFMN(VMCP,VMCN,DYGF(N),DYGF(NN),VGNX,INDX)
          UGNX = UGNX/VMCX
          VGNX = VGNX/VMCX
          WGNX = WGNX/VMCX
        ENDIF
        UGNSQ = UGNX*UGNX
        VGNSQ = VGNX*VGNX
        WGNSQ = WGNX*WGNX
        ZVN = SQRT(UGNSQ+VGNSQ+WGNSQ)
        INDX = 17
        DPLN = DIFMN(DISPL(N),DISPL(NN),DYGF(N),DYGF(NN),VGNX,INDX)
        DPTN = DIFMN(DISPT(N),DISPT(NN),DYGF(N),DYGF(NN),VGNX,INDX)
        DPN = (DPLN*VGNSQ + DPTN*(UGNSQ+WGNSQ))/(ZVN+SMALL)
      ENDIF
!
!---  Top face  ---
!
      NT = ICM(6,N)
      IF( NT.NE.0 ) THEN
        UGTX = 0.25D+0*((UG(1,1,N)+UG(1,2,N))*DZGF(NT) +
     &    (UG(1,1,NT)+UG(1,2,NT))*DZGF(N))/DZGP(2,N)
        VGTX = 0.25D+0*((VG(1,1,N)+VG(1,2,N))*DZGF(NT)
     &  + (VG(1,1,NT)+VG(1,2,NT))*DZGF(N))/DZGP(2,N)
        WGTX = WG(1,2,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SG(2,N)*PORD(2,N)
          VMCT = SG(2,NT)*PORD(2,NT)
          INDX = 17
          VMCX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NT),WGTX,INDX)
          UGTX = UGTX/VMCX
          VGTX = VGTX/VMCX
          WGTX = WGTX/VMCX
        ENDIF
        UGTSQ = UGTX*UGTX
        VGTSQ = VGTX*VGTX
        WGTSQ = WGTX*WGTX
        ZVT = SQRT(UGTSQ+VGTSQ+WGTSQ)
        INDX = 17
        DPLT = DIFMN(DISPL(N),DISPL(NT),DZGF(N),DZGF(NT),WGTX,INDX)
        DPTT = DIFMN(DISPT(N),DISPT(NT),DZGF(N),DZGF(NT),WGTX,INDX)
        DPT = (DPLT*WGTSQ + DPTT*(VGTSQ+UGTSQ))/(ZVT+SMALL)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SHDPG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SHDPL( N,DPB,DPS,DPW,DPE,DPN,DPT,NSL )
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
!     STOMPX-CO2
!
!     Calculates hydrodynamic dispersion coefficients for the aqueous
!     phase from phase velocities and user-specified dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/SHDPL'
!
!---  Initialize dispersion coefficients  ---
!
      M = 1
      DPB = 0.D+0
      DPS = 0.D+0
      DPW = 0.D+0
      DPE = 0.D+0
      DPN = 0.D+0
      DPT = 0.D+0
!
!---  Bottom face  ---
!
      NB = ICM(1,N)
      IF( NB.NE.0 ) THEN
        ULBX = 0.25D+0*((UL(1,1,N)+UL(1,2,N))*DZGF(NB) +
     &    (UL(1,1,NB)+UL(1,2,NB))*DZGF(N))/DZGP(1,N)
        VLBX = 0.25D+0*((VL(1,1,N)+VL(1,2,N))*DZGF(NB)
     &  + (VL(1,1,NB)+VL(1,2,NB))*DZGF(N))/DZGP(1,N)
        WLBX = WL(1,1,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SL(2,N)*PORD(2,N)
          VMCB = SL(2,NB)*PORD(2,NB)
          INDX = 17
          VMCX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WLBX,INDX)
          ULBX = ULBX/VMCX
          VLBX = VLBX/VMCX
          WLBX = WLBX/VMCX
        ENDIF
        ULBSQ = ULBX*ULBX
        VLBSQ = VLBX*VLBX
        WLBSQ = WLBX*WLBX
        ZVB = SQRT(ULBSQ+VLBSQ+WLBSQ)
        INDX = 17
        DPLB = DIFMN(DISPL(NB),DISPL(N),DZGF(NB),DZGF(N),WLBX,INDX)
        DPTB = DIFMN(DISPT(NB),DISPT(N),DZGF(NB),DZGF(N),WLBX,INDX)
        DPB = (DPLB*WLBSQ + DPTB*(VLBSQ+ULBSQ))/(ZVB+SMALL)
      ENDIF
!
!---  South face  ---
!
      NS = ICM(2,N)
      IF( NS.NE.0 ) THEN
        ULSX = 0.25D+0*((UL(1,1,N)+UL(1,2,N))*DYGF(NS) +
     &    (UL(1,1,NS)+UL(1,2,NS))*DYGF(N))/DYGP(1,N)
        VLSX = VL(1,1,N)
        WLSX = 0.25D+0*((WL(1,1,N)+WL(1,2,N))*DYGF(NS) + 
     &    (WL(1,1,NS)+WL(1,2,NS))*DYGF(N))/DYGP(1,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SL(2,N)*PORD(2,N)
          VMCS = SL(2,NS)*PORD(2,NS)
          INDX = 17
          VMCX = DIFMN(VMCS,VMCP,DYGF(NS),DYGF(N),VLSX,INDX)
          ULSX = ULSX/VMCX
          VLSX = VLSX/VMCX
          WLSX = WLSX/VMCX
        ENDIF
        ULSSQ = ULSX*ULSX
        VLSSQ = VLSX*VLSX
        WLSSQ = WLSX*WLSX
        ZVS = SQRT(ULSSQ+VLSSQ+WLSSQ)
        INDX = 17
        DPLS = DIFMN(DISPL(NS),DISPL(N),DYGF(NS),DYGF(N),VLSX,INDX)
        DPTS = DIFMN(DISPT(NS),DISPT(N),DYGF(NS),DYGF(N),VLSX,INDX)
        DPS = (DPLS*VLSSQ + DPTS*(ULSSQ+WLSSQ))/(ZVS+SMALL)
      ENDIF
!
!---  West face  ---
!
      NW = ICM(3,N)
      IF( NW.NE.0 ) THEN
        ULWX = UL(M,1,N)
        VLWX = 0.25D+0*((VL(1,1,N)+VL(1,2,N))*DXGF(NW) + 
     &    (VL(1,1,NW)+VL(1,2,NW))*DXGF(N))/DXGP(1,N)
        WLWX = 0.25D+0*((WL(1,1,N)+WL(1,2,N))*DXGF(NW) + 
     &    (WL(1,1,NW)+WL(1,2,NW))*DXGF(N))/DXGP(1,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SL(2,N)*PORD(2,N)
          VMCW = SL(2,NW)*PORD(2,NW)
          INDX = 17
          VMCX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),ULWX,INDX)
          ULWX = ULWX/VMCX
          VLWX = VLWX/VMCX
          WLWX = WLWX/VMCX
        ENDIF
        ULWSQ = ULWX*ULWX
        VLWSQ = VLWX*VLWX
        WLWSQ = WLWX*WLWX
        ZVW = SQRT(ULWSQ+VLWSQ+WLWSQ)
        INDX = 17
        DPLW = DIFMN(DISPL(NW),DISPL(N),DXGF(NW),DXGF(N),ULWX,INDX)
        DPTW = DIFMN(DISPT(NW),DISPT(N),DXGF(NW),DXGF(N),ULWX,INDX)
        DPW = (DPLW*ULWSQ + DPTW*(VLWSQ+WLWSQ))/(ZVW+SMALL)
      ENDIF
!
!---  East face  ---
!
      NE = ICM(4,N)
      IF( NE.NE.0 ) THEN
        ULEX = UL(M,2,N)
        VLEX = 0.25D+0*((VL(1,1,N)+VL(1,2,N))*DXGF(NE) + 
     &    (VL(1,1,NE)+VL(1,2,NE))*DXGF(N))/DXGP(2,N)
        WLEX = 0.25D+0*((WL(1,1,N)+WL(1,2,N))*DXGF(NE) + 
     &    (WL(1,1,NE)+WL(1,2,NE))*DXGF(N))/DXGP(2,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SL(2,N)*PORD(2,N)
          VMCE = SL(2,NE)*PORD(2,NE)
          INDX = 17
          VMCX = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NE),ULEX,INDX)
          ULEX = ULEX/VMCX
          VLEX = VLEX/VMCX
          WLEX = WLEX/VMCX
        ENDIF
        ULESQ = ULEX*ULEX
        VLESQ = VLEX*VLEX
        WLESQ = WLEX*WLEX
        ZVE = SQRT(ULESQ+VLESQ+WLESQ)
        INDX = 17
        DPLE = DIFMN(DISPL(N),DISPL(NE),DXGF(N),DXGF(NE),ULEX,INDX)
        DPTE = DIFMN(DISPT(N),DISPT(NE),DXGF(N),DXGF(NE),ULEX,INDX)
        DPE = (DPLE*ULESQ + DPTE*(VLESQ+WLESQ))/(ZVE+SMALL)
      ENDIF
!
!---  North face  ---
!
      NN = ICM(5,N)
      IF( NN.NE.0 ) THEN
        ULNX = 0.25D+0*((UL(1,1,N)+UL(1,2,N))*DYGF(NN) +
     &    (UL(1,1,NN)+UL(1,2,NN))*DYGF(N))/DYGP(2,N)
        VLNX = VL(M,2,N)
        WLNX = 0.25D+0*((WL(1,1,N)+WL(1,2,N))*DYGF(NN) + 
     &    (WL(1,1,NN)+WL(1,2,NN))*DYGF(N))/DYGP(2,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SL(2,N)*PORD(2,N)
          VMCN = SL(2,NN)*PORD(2,NN)
          INDX = 17
          VMCX = DIFMN(VMCP,VMCN,DYGF(N),DYGF(NN),VLNX,INDX)
          ULNX = ULNX/VMCX
          VLNX = VLNX/VMCX
          WLNX = WLNX/VMCX
        ENDIF
        ULNSQ = ULNX*ULNX
        VLNSQ = VLNX*VLNX
        WLNSQ = WLNX*WLNX
        ZVN = SQRT(ULNSQ+VLNSQ+WLNSQ)
        INDX = 17
        DPLN = DIFMN(DISPL(N),DISPL(NN),DYGF(N),DYGF(NN),VLNX,INDX)
        DPTN = DIFMN(DISPT(N),DISPT(NN),DYGF(N),DYGF(NN),VLNX,INDX)
        DPN = (DPLN*VLNSQ + DPTN*(ULNSQ+WLNSQ))/(ZVN+SMALL)
      ENDIF
!
!---  Top face  ---
!
      NT = ICM(6,N)
      IF( NT.NE.0 ) THEN
        ULTX = 0.25D+0*((UL(1,1,N)+UL(1,2,N))*DZGF(NT) +
     &    (UL(1,1,NT)+UL(1,2,NT))*DZGF(N))/DZGP(2,N)
        VLTX = 0.25D+0*((VL(1,1,N)+VL(1,2,N))*DZGF(NT)
     &  + (VL(1,1,NT)+VL(1,2,NT))*DZGF(N))/DZGP(2,N)
        WLTX = WL(1,2,N)
        IF( ISLC(52).EQ.0 ) THEN
          VMCP = SL(2,N)*PORD(2,N)
          VMCT = SL(2,NT)*PORD(2,NT)
          INDX = 17
          VMCX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NT),WLTX,INDX)
          ULTX = ULTX/VMCX
          VLTX = VLTX/VMCX
          WLTX = WLTX/VMCX
        ENDIF
        ULTSQ = ULTX*ULTX
        VLTSQ = VLTX*VLTX
        WLTSQ = WLTX*WLTX
        ZVT = SQRT(ULTSQ+VLTSQ+WLTSQ)
        INDX = 17
        DPLT = DIFMN(DISPL(N),DISPL(NT),DZGF(N),DZGF(NT),WLTX,INDX)
        DPTT = DIFMN(DISPT(N),DISPT(NT),DZGF(N),DZGF(NT),WLTX,INDX)
        DPT = (DPLT*WLTSQ + DPTT*(VLTSQ+ULTSQ))/(ZVT+SMALL)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SHDPL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBG( NSL )
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
!     STOMPX-CO2
!
!     Loads the matrix elements and solution vector for the
!     gas-phase convective-dispersive mass transport equation.
!
!     The Jacobian matrix is initially configured assuming zero-flux
!     boundary conditions.  The matrix is then updated for other
!     user-specified boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 January 2022

!----------------------LIS Modules-------------------------------------!
!
      USE LIS_STOMP
!







!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE JACOB
      USE HYST
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/SJCBG'
!
!---  Loop over local nodes, skipping inactive nodes
!     or ghost cells  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    Storage terms  ---
!
        SC = VOL(N)*DTI
        CC = 0.D+0
!
!---    Load Jacobian matrix  ---
!

        MCP = IXP(N)
        MA = 1
        MCD = KLUC(MCP-IEQC_OFFSET,MA)
        MA = MA + 1
        IF( IEQW.EQ.0 ) DLU(MCD) = DLU(MCD) + SC
!        IROW = IXP(N)
!        ICOL = IXP(N)
!        IF( IEQW.EQ.0 ) THEN
!          BUFFER = SC
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!        ENDIF








!
!---    Molecular diffusion coefficients at the nodes  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        PCOR = (PG(2,N)+PATM)/PATM
        SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
        VMCP = SG(2,N)*PORD(2,N)
        FCGP = YG(N,NSL)/(VMCP+SMALL)
        DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
!
!---    Hydrodynamic dispersion coefficients at cell faces  ---
!
        CALL SHDPG( N,DPGB,DPGS,DPGW,DPGE,DPGN,DPGT,NSL )
!
!---    Bottom face diffusion and advection terms  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          TCOR = (T(2,NB)+TABS)/TSPRF
          PCOR = (PG(2,NB)+PATM)/PATM
          SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCB = SG(2,NB)*PORD(2,NB)
          FCGB = YG(NB,NSL)/(VMCB+SMALL)
          DGB = TORG(2,NB)*(SG(2,NB)-SGT(2,NB))*PORD(2,NB)*SDFGB
          INDX = 16
          DGZ = DIFMN(DGB,DGP,DZGF(NB),DZGF(N),WG(1,1,N),INDX)
          DGZ = AFZ(1,N)*(DGZ+DPGB)/DZGP(1,N)
          FGB = AFZ(1,N)*WG(1,1,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FGB = 0.D+0
          VMCBX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WG(1,1,N),INDX)
          CRGB = ABS(WG(1,1,N))*DT/(DZGP(1,N)*VMCBX+SMALL)
!
!---      Patankar solute transport  ---
!
          AGB = MAX(FGB,ZERO)
     &      + DGZ*MAX((ONE-(TENTH*ABS(FGB)/(DGZ+SMALL)))**5,ZERO)
          AP = (AGB-FGB)*FCGP
          AB = AGB*FCGB
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          DLU(MROW) = DLU(MROW) - AB
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NB)
!          BUFFER = -AB
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )







        ENDIF
!
!---    South face diffusion and advection terms  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DYSR = RP(N)*DYGF(NS)
          TCOR = (T(2,NS)+TABS)/TSPRF
          PCOR = (PG(2,NS)+PATM)/PATM
          SDFGS = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCS = SG(2,NS)*PORD(2,NS)
          FCGS = YG(NS,NSL)/(VMCS+SMALL)
          DGS = TORG(2,NS)*(SG(2,NS)-SGT(2,NS))*PORD(2,NS)*SDFGS
          INDX = 16
          DGY = DIFMN(DGS,DGP,DYSR,DYR,VG(1,1,N),INDX)
          DGY = AFY(1,N)*(DGY+DPGS)/(DYGP(1,N)*RP(N))
          FGS = AFY(1,N)*VG(1,1,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FGS = 0.D+0
          VMCSX = DIFMN(VMCS,VMCP,DYSR,DYR,VG(1,1,N),INDX)
          CRGS = ABS(VG(1,1,N))*DT/(DYGP(1,N)*VMCSX+SMALL)/RP(N)
!
!---      Patankar solute transport  ---
!
          AGS = MAX(FGS,ZERO)
     &      + DGY*MAX((ONE-(TENTH*ABS(FGS)/(DGY+SMALL)))**5,ZERO)
          AP = (AGS-FGS)*FCGP
          AS = AGS*FCGS
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          DLU(MROW) = DLU(MROW) - AS
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NS)
!          BUFFER = -AS
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )







        ENDIF
!
!---    West face diffusion and advection terms  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          TCOR = (T(2,NW)+TABS)/TSPRF
          PCOR = (PG(2,NW)+PATM)/PATM
          SDFGW = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCW = SG(2,NW)*PORD(2,NW)
          FCGW = YG(NW,NSL)/(VMCW+SMALL)
          DGW = TORG(2,NW)*(SG(2,NW)-SGT(2,NW))*PORD(2,NW)*SDFGW
          INDX = 16
          DGX = DIFMN(DGW,DGP,DXGF(NW),DXGF(N),UG(1,1,N),INDX)
          DGX = AFX(1,N)*(DGX+DPGW)/DXGP(1,N)
          FGW = AFX(1,N)*UG(1,1,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FGW = 0.D+0
          VMCWX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UG(1,1,N),INDX)
          CRGW = ABS(UG(1,1,N))*DT/(DXGP(1,N)*VMCWX+SMALL)
!
!---      Patankar solute transport  ---
!
          AGW = MAX(FGW,ZERO)
     &      + DGX*MAX((ONE-(TENTH*ABS(FGW)/(DGX+SMALL)))**5,ZERO)
          AP = (AGW-FGW)*FCGP
          AW = AGW*FCGW
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          DLU(MROW) = DLU(MROW) - AW
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NW)
!          BUFFER = -AW
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )







        ENDIF
!
!---    East face diffusion and advection terms  ---
!
        NE = ICM(4,N)
        IF( NE.NE.0 ) THEN
          TCOR = (T(2,NE)+TABS)/TSPRF
          PCOR = (PG(2,NE)+PATM)/PATM
          SDFGE = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCE = SG(2,NE)*PORD(2,NE)
          FCGE = YG(NE,NSL)/(VMCE+SMALL)
          DGE = TORG(2,NE)*(SG(2,NE)-SGT(2,NE))*PORD(2,NE)*SDFGE
          INDX = 16
          DGX = DIFMN(DGP,DGE,DXGF(N),DXGF(NE),UG(1,2,N),INDX)
          DGX = AFX(2,N)*(DGX+DPGE)/DXGP(2,N)
          FGE = AFX(2,N)*UG(1,2,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FGE = 0.D+0
          VMCEX = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NE),UG(1,2,N),INDX)
          CRGE = ABS(UG(1,2,N))*DT/(DXGP(2,N)*VMCEX+SMALL)
!
!---      Patankar solute transport  ---
!
          AGE = MAX(-FGE,ZERO)
     &      + DGX*MAX((ONE-(TENTH*ABS(FGE)/(DGX+SMALL)))**5,ZERO)
          AP = (AGE+FGE)*FCGP
          AE = AGE*FCGE
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          DLU(MROW) = DLU(MROW) - AE
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NE)
!          BUFFER = -AE
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )







        ENDIF
!
!---    North face diffusion and advection terms  ---
!
        NN = ICM(5,N)
        IF( NN.NE.0 ) THEN
          DYNR = RP(N)*DYGF(NN)
          TCOR = (T(2,NN)+TABS)/TSPRF
          PCOR = (PG(2,NN)+PATM)/PATM
          SDFGN = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCN = SG(2,NN)*PORD(2,NN)
          FCGN = YG(NN,NSL)/(VMCN+SMALL)
          DGN = TORG(2,NN)*(SG(2,NN)-SGT(2,NN))*PORD(2,NN)*SDFGN
          INDX = 16
          DGY = DIFMN(DGP,DGN,DYNR,DYR,VG(1,2,N),INDX)
          DGY = AFY(2,N)*(DGY+DPGN)/(DYGP(2,N)*RP(N))
          FGN = AFY(2,N)*VG(1,2,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FGS = 0.D+0
          VMCNX = DIFMN(VMCP,VMCN,DYNR,DYR,VG(1,2,N),INDX)
          CRGN = ABS(VG(1,2,N))*DT/(DYGP(2,N)*VMCNX+SMALL)/RP(N)
!
!---      Patankar solute transport  ---
!
          AGN = MAX(-FGN,ZERO)
     &      + DGY*MAX((ONE-(TENTH*ABS(FGN)/(DGY+SMALL)))**5,ZERO)
          AP = (AGN+FGN)*FCGP
          AN = AGN*FCGN
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          DLU(MROW) = DLU(MROW) - AN
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NN)
!          BUFFER = -AN
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )







        ENDIF
!
!---    Top face diffusion and advection terms  ---
!
        NT = ICM(6,N)
        IF( NT.NE.0 ) THEN
          TCOR = (T(2,NT)+TABS)/TSPRF
          PCOR = (PG(2,NT)+PATM)/PATM
          SDFGT = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCT = SG(2,NT)*PORD(2,NT)
          FCGT = YG(NT,NSL)/(VMCT+SMALL)
          DGT = TORG(2,NT)*(SG(2,NT)-SGT(2,NT))*PORD(2,NT)*SDFGT
          INDX = 16
          DGZ = DIFMN(DGP,DGT,DZGF(N),DZGF(NT),WG(1,2,N),INDX)
          DGZ = AFZ(2,N)*(DGZ+DPGT)/DZGP(2,N)
          FGT = AFZ(2,N)*WG(1,2,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FGT = 0.D+0
          VMCTX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NT),WG(1,2,N),INDX)
          CRGT = ABS(WG(1,2,N))*DT/(DZGP(2,N)*VMCTX+SMALL)
!
!---      Patankar solute transport  ---
!
          AGT = MAX(-FGT,ZERO)
     &      + DGZ*MAX((ONE-(TENTH*ABS(FGT)/(DGZ+SMALL)))**5,ZERO)
          AP = (AGT+FGT)*FCGP
          AT = AGT*FCGT
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          DLU(MROW) = DLU(MROW) - AT
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NT)
!          BUFFER = -AT
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )







        ENDIF
!
!---    Solution vector  ---
!

        IF( IEQW.EQ.0 ) THEN
          BUFFER = CO(N,NSL)*SC
          CALL lis_vector_set_value_f( 1,MCP,BUFFER,
     &      T_RHS_VEC,IERR )
        ENDIF








      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBL( NSL )
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
!     STOMPX-CO2
!
!     Loads the matrix elements and solution vector for the
!     aqueous-phase convective-dispersive mass transport equation.
!
!     The Jacobian matrix is initially configured assuming zero-flux
!     boundary conditions.  The matrix is then updated for other
!     user-specified boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 January 2022

!----------------------LIS Modules------------------------------------!
!
      USE LIS_STOMP
!







!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE JACOB
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/SJCBL'
!
!---  Loop over local nodes, skipping inactive nodes
!     or ghost cells  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    Storage terms  ---
!
        SC = VOL(N)*DTI
        CC = 0.D+0
!
!---    Load Jacobian matrix  ---
!

        MCP = IXP(N)
        MA = 1
        MCD = KLUC(MCP-IEQC_OFFSET,MA)
        MA = MA + 1
        IF( MCD.LE.0 .OR. MCD.GT.NNZC ) THEN
          PRINT *,'N = ',N,'ID = ',ID
          RETURN
        ENDIF
        DLU(MCD) = DLU(MCD) + SC
!        IROW = IXP(N)
!        ICOL = IXP(N)
!        BUFFER = SC
!        IF( BUFFER.NE.BUFFER ) PRINT *,'1:ND(N) = ',ND(N),' NSL = ',NSL,
!     &    ' N = ',N,' ID = ',ID
!        CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &    T_MAT,IERR )








!
!---    Molecular diffusion coefficients at the nodes  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
        VMCP = SL(2,N)*PORD(2,N)
        FCLP = YL(N,NSL)/(VMCP+SMALL)
        IF( IEDL(NSL).EQ.2 ) THEN
          DLP = SDCL(1,N,NSL)*SDCL(2,N,NSL)*EXP(VMCP*SDCL(3,N,NSL))
        ELSEIF( IEDL(NSL).EQ.3 ) THEN
          DLP = TORL(2,N)*VMCP*SMDL(NSL)
        ELSEIF( IEDL(NSL).EQ.4 ) THEN
          DLP = SDCL(1,N,NSL)*SDCL(2,N,NSL)*VMCP**SDCL(3,N,NSL)
        ELSE
          DLP = TORL(2,N)*VMCP*SDFLP
        ENDIF
!
!---    Hydrodynamic dispersion coefficients at cell faces  ---
!
        CALL SHDPL( N,DPLB,DPLS,DPLW,DPLE,DPLN,DPLT,NSL )
!
!---    Bottom face diffusion and advection terms  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          TCOR = (T(2,NB)+TABS)/TSPRF
          SDFLB = SMDL(NSL)*TCOR*(VISRL/VISL(2,NB))
          VMCB = SL(2,NB)*PORD(2,NB)
          FCLB = YL(NB,NSL)/(VMCB+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLB = SDCL(1,NB,NSL)*SDCL(2,NB,NSL)*EXP(VMCB*SDCL(3,NB,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLB = TORL(2,NB)*VMCB*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLB = SDCL(1,NB,NSL)*SDCL(2,NB,NSL)*VMCB**SDCL(3,NB,NSL)
          ELSE
            DLB = TORL(2,NB)*VMCB*SDFLB
          ENDIF
          INDX = 16
          DLZ = DIFMN(DLB,DLP,DZGF(NB),DZGF(N),WL(1,1,N),INDX)
          DLZ = AFZ(1,N)*(DLZ+DPLB)/DZGP(1,N)
          FLB = AFZ(1,N)*WL(1,1,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FLB = 0.D+0
          VMCBX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WL(1,1,N),INDX)
          CRLB = ABS(WL(1,1,N))*DT/(DZGP(1,N)*VMCBX+SMALL)
!
!---      Patankar solute transport  ---
!
          ALB = MAX(FLB,ZERO)
     &      + DLZ*MAX((ONE-(TENTH*ABS(FLB)/(DLZ+SMALL)))**5,ZERO)
          AP = (ALB-FLB)*FCLP
          AB = ALB*FCLB
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          IF( MROW.LE.0 .OR. MROW.GT.NNZC ) THEN
            PRINT *,'NB = ',NB,'ID = ',ID
            RETURN
          ENDIF
          DLU(MROW) = DLU(MROW) - AB
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NB)
!          BUFFER = -AB
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )







        ENDIF
!
!---    South face diffusion and advection terms  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DYSR = RP(N)*DYGF(NS)
          TCOR = (T(2,NS)+TABS)/TSPRF
          SDFLS = SMDL(NSL)*TCOR*(VISRL/VISL(2,NS))
          VMCS = SL(2,NS)*PORD(2,NS)
          FCLS = YL(NS,NSL)/(VMCS+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLS = SDCL(1,NS,NSL)*SDCL(2,NS,NSL)*EXP(VMCS*SDCL(3,NS,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLS = TORL(2,NS)*VMCS*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLS = SDCL(1,NS,NSL)*SDCL(2,NS,NSL)*VMCS**SDCL(3,NS,NSL)
          ELSE
            DLS = TORL(2,NS)*VMCS*SDFLS
          ENDIF
          INDX = 16
          DLY = DIFMN(DLS,DLP,DYSR,DYR,VL(1,1,N),INDX)
          DLY = AFY(1,N)*(DLY+DPLS)/(DYGP(1,N)*RP(N))
          FLS = AFY(1,N)*VL(1,1,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FLS = 0.D+0
          VMCSX = DIFMN(VMCS,VMCP,DYSR,DYR,VL(1,1,N),INDX)
          CRLS = ABS(VL(1,1,N))*DT/(DYGP(1,N)*VMCSX+SMALL)/RP(N)
!
!---      Patankar solute transport  ---
!
          ALS = MAX(FLS,ZERO)
     &      + DLY*MAX((ONE-(TENTH*ABS(FLS)/(DLY+SMALL)))**5,ZERO)
          AP = (ALS-FLS)*FCLP
          AS = ALS*FCLS
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          IF( MROW.LE.0 .OR. MROW.GT.NNZC ) THEN
            PRINT *,'NS = ',NS,'ID = ',ID
            RETURN
          ENDIF
          DLU(MROW) = DLU(MROW) - AS
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NS)
!          BUFFER = -AS
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )







        ENDIF
!
!---    West face diffusion and advection terms  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          TCOR = (T(2,NW)+TABS)/TSPRF
          SDFLW = SMDL(NSL)*TCOR*(VISRL/VISL(2,NW))
          VMCW = SL(2,NW)*PORD(2,NW)
          FCLW = YL(NW,NSL)/(VMCW+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLW = SDCL(1,NW,NSL)*SDCL(2,NW,NSL)*EXP(VMCW*SDCL(3,NW,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLW = TORL(2,NW)*VMCW*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLW = SDCL(1,NW,NSL)*SDCL(2,NW,NSL)*VMCW**SDCL(3,NW,NSL)
          ELSE
            DLW = TORL(2,NW)*VMCW*SDFLW
          ENDIF
          INDX = 16
          DLX = DIFMN(DLW,DLP,DXGF(NW),DXGF(N),UL(1,1,N),INDX)
          DLX = AFX(1,N)*(DLX+DPLW)/DXGP(1,N)
          FLW = AFX(1,N)*UL(1,1,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FLW = 0.D+0
          VMCWX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UL(1,1,N),INDX)
          CRLW = ABS(UL(1,1,N))*DT/(DXGP(1,N)*VMCWX+SMALL)
!
!---      Patankar solute transport  ---
!
          ALW = MAX(FLW,ZERO)
     &      + DLX*MAX((ONE-(TENTH*ABS(FLW)/(DLX+SMALL)))**5,ZERO)
          AP = (ALW-FLW)*FCLP
          AW = ALW*FCLW
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          IF( MROW.LE.0 .OR. MROW.GT.NNZC ) THEN
            PRINT *,'NW = ',NW,'ID = ',ID
            RETURN
          ENDIF
          DLU(MROW) = DLU(MROW) - AW
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NW)
!          BUFFER = -AW
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )







        ENDIF
!
!---    East face diffusion and advection terms  ---
!
        NE = ICM(4,N)
        IF( NE.NE.0 ) THEN
          TCOR = (T(2,NE)+TABS)/TSPRF
          SDFLE = SMDL(NSL)*TCOR*(VISRL/VISL(2,NE))
          VMCE = SL(2,NE)*PORD(2,NE)
          FCLE = YL(NE,NSL)/(VMCE+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLE = SDCL(1,NE,NSL)*SDCL(2,NE,NSL)*EXP(VMCE*SDCL(3,NE,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLE = TORL(2,NE)*VMCE*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLE = SDCL(1,NE,NSL)*SDCL(2,NE,NSL)*VMCE**SDCL(3,NE,NSL)
          ELSE
            DLE = TORL(2,NE)*VMCE*SDFLE
          ENDIF
          INDX = 16
          DLX = DIFMN(DLP,DLE,DXGF(N),DXGF(NE),UL(1,2,N),INDX)
          DLX = AFX(2,N)*(DLX+DPLE)/DXGP(2,N)
          FLE = AFX(2,N)*UL(1,2,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FLE = 0.D+0
          VMCEX = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NE),UL(1,2,N),INDX)
          CRLE = ABS(UL(1,2,N))*DT/(DXGP(2,N)*VMCEX+SMALL)
!
!---      Patankar solute transport  ---
!
          ALE = MAX(-FLE,ZERO)
     &      + DLX*MAX((ONE-(TENTH*ABS(FLE)/(DLX+SMALL)))**5,ZERO)
          AP = (ALE+FLE)*FCLP
          AE = ALE*FCLE
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          IF( MROW.LE.0 .OR. MROW.GT.NNZC ) THEN
            PRINT *,'NE = ',NE,'ID = ',ID
            RETURN
          ENDIF
          DLU(MROW) = DLU(MROW) - AE
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NE)
!          BUFFER = -AE
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )








        ENDIF
!
!---    North face diffusion and advection terms  ---
!
        NN = ICM(5,N)
        IF( NN.NE.0 ) THEN
          DYNR = RP(N)*DYGF(NN)
          TCOR = (T(2,NN)+TABS)/TSPRF
          SDFLN = SMDL(NSL)*TCOR*(VISRL/VISL(2,NN))
          VMCN = SL(2,NN)*PORD(2,NN)
          FCLN = YL(NN,NSL)/(VMCN+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLN = SDCL(1,NN,NSL)*SDCL(2,NN,NSL)*EXP(VMCN*SDCL(3,NN,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLN = TORL(2,NN)*VMCN*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLN = SDCL(1,NN,NSL)*SDCL(2,NN,NSL)*VMCN**SDCL(3,NN,NSL)
          ELSE
            DLN = TORL(2,NN)*VMCN*SDFLN
          ENDIF
          INDX = 16
          DLY = DIFMN(DLP,DLN,DYNR,DYR,VL(1,2,N),INDX)
          DLY = AFY(2,N)*(DLY+DPLN)/(DYGP(2,N)*RP(N))
          FLN = AFY(2,N)*VL(1,2,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FLS = 0.D+0
          VMCNX = DIFMN(VMCP,VMCN,DYNR,DYR,VL(1,2,N),INDX)
          CRLN = ABS(VL(1,2,N))*DT/(DYGP(2,N)*VMCNX+SMALL)/RP(N)
!
!---      Patankar solute transport  ---
!
          ALN = MAX(-FLN,ZERO)
     &      + DLY*MAX((ONE-(TENTH*ABS(FLN)/(DLY+SMALL)))**5,ZERO)
          AP = (ALN+FLN)*FCLP
          AN = ALN*FCLN
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          IF( MROW.LE.0 .OR. MROW.GT.NNZC ) THEN
            PRINT *,'NN = ',NN,'ID = ',ID
            RETURN
          ENDIF
          DLU(MROW) = DLU(MROW) - AN
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NN)
!          BUFFER = -AN
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )







        ENDIF
!
!---    Top face diffusion and advection terms  ---
!
        NT = ICM(6,N)
        IF( NT.NE.0 ) THEN
          TCOR = (T(2,NT)+TABS)/TSPRF
          SDFLT = SMDL(NSL)*TCOR*(VISRL/VISL(2,NT))
          VMCT = SL(2,NT)*PORD(2,NT)
          FCLT = YL(NT,NSL)/(VMCT+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLT = SDCL(1,NT,NSL)*SDCL(2,NT,NSL)*EXP(VMCT*SDCL(3,NT,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLT = TORL(2,NT)*VMCT*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLT = SDCL(1,NT,NSL)*SDCL(2,NT,NSL)*VMCT**SDCL(3,NT,NSL)
          ELSE
            DLT = TORL(2,NT)*VMCT*SDFLT
          ENDIF
          INDX = 16
          DLZ = DIFMN(DLP,DLT,DZGF(N),DZGF(NT),WL(1,2,N),INDX)
          DLZ = AFZ(2,N)*(DLZ+DPLT)/DZGP(2,N)
          FLT = AFZ(2,N)*WL(1,2,N)
          IF( MOD(ISLC(23),10).EQ.1 ) FLT = 0.D+0
          VMCTX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NT),WL(1,2,N),INDX)
          CRLT = ABS(WL(1,2,N))*DT/(DZGP(2,N)*VMCTX+SMALL)
!
!---      Patankar solute transport  ---
!
          ALT = MAX(-FLT,ZERO)
     &      + DLZ*MAX((ONE-(TENTH*ABS(FLT)/(DLZ+SMALL)))**5,ZERO)
          AP = (ALT+FLT)*FCLP
          AT = ALT*FCLT
!
!---      Load Jacobian matrix  ---
!

          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          IF( MROW.LE.0 .OR. MROW.GT.NNZC ) THEN
            PRINT *,'NT = ',NT,'ID = ',ID
            RETURN
          ENDIF
          DLU(MROW) = DLU(MROW) - AT
!          ICOL = IXP(N)
!          BUFFER = AP
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )
!          ICOL = IXP(NT)
!          BUFFER = -AT
!          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &      T_MAT,IERR )


          DLU(MCD) = DLU(MCD) + AP
          MROW = KLUC(MCP-IEQC_OFFSET,MA)
          MA = MA + 1
          DLU(MROW) = DLU(MROW) - AT

        ENDIF
!
!---    Solution vector  ---
!

        BUFFER = CO(N,NSL)*SC
        CALL lis_vector_set_value_f( 1,MCP,BUFFER,
     &    T_RHS_VEC,IERR )






      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORIT_CO2( NSL )
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
!     STOMPX-CO2
!
!     Compute solute transport source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 27 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE GRID
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
      REAL*8 SRX(8+LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORIT_CO2'
!
!---  Loop over processor dependent (i.e., local) sources  ---
!
      DO NS = 1,NSR(ID+1)
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = ISRIN(NS)
        IF( TMZ.LE.SRC(1,1,MB) ) CYCLE
        IF( ISRM(NS).EQ.1 ) THEN
          DO N = 1,8
            SRX(N) = SRC(N,1,MB)
          ENDDO
        ELSE
          IFIND = 0
          DO M = 2,ISRM(NS)
            IF( TMZ.LE.SRC(1,M,MB) ) THEN
             DTSR = MIN( SRC(1,M,MB)-TMZ,DT )
             TFSR = (TMZ-0.5D+0*DTSR-SRC(1,M-1,MB))/
     &         (SRC(1,M,MB)-SRC(1,M-1,MB))
             DO N = 1,8
               SRX(N) = SRC(N,M-1,MB) + TFSR*(SRC(N,M,MB)-SRC(N,M-1,MB))
             ENDDO
             IFIND = 1
             EXIT
            ENDIF
          ENDDO
          IF( IFIND.EQ.0 ) CYCLE
        ENDIF
!
!---    Local node associated with processor-dependent (i.e., local)
!       source  ---
!
        N = ISRN(NS)
!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
!
!---    Aqueous Volumetric Sink  ---
!
        IF( MOD(ISRT(NS),100).EQ.3 .AND. SRX(4).LT.0.D+0 ) THEN
          SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &      YL(N,NSL)*DT/(PORD(2,N)*SL(2,N))
!
!---    Gas Volumetric Sink  ---
!
        ELSEIF( (MOD(ISRT(NS),100).EQ.4 .OR. 
     &    MOD(ISRT(NS),100).EQ.5) .AND. SRX(4).LT.0.D+0 ) THEN
          SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &      YG(N,NSL)*DT/(PORD(2,N)*SG(2,N))
!
!---    Aqueous Mass Sink  ---
!
        ELSEIF( MOD(ISRT(NS),100).EQ.7 .AND. 
     &    SRX(4).LT.0.D+0 ) THEN
          SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &      YL(N,NSL)*DT/(RHOL(2,N)*PORD(2,N)*SL(2,N))
!
!---    Gas Mass Sink  ---
!
        ELSEIF( (MOD(ISRT(NS),100).EQ.8 .OR. 
     &    MOD(ISRT(NS),100).EQ.9) .AND. SRX(4).LT.0.D+0 ) THEN
          SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &      YG(N,NSL)*DT/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---    Solute source  ---
!
        ELSEIF( ISRT(NS).LT.0 .AND. ISRT(NS).GE.-NSOLU ) THEN
          SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(4)*DT
!
!---   Solute source  ---
!
        ELSEIF( ISRT(NS).LT.-NSOLU .AND.
     &    ISRT(NS).GE.-2*NSOLU ) THEN
          SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(4)*DT*VOL(N)
        ENDIF
      ENDDO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORIT_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_CO2( NSL )
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
!     STOMPX-CO2
!
!     Compute solute transport source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 January 2022

!----------------------LIS Modules-------------------------------------!
!
      USE LIS_STOMP
!







!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE PROP
      USE JACOB
      USE GRID
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
      REAL*8 SRX(8+LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORT_CO2'
!
!---  Loop over processor dependent (i.e., local) sources  ---
!
      DO NS = 1,NSR(ID+1)
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = ISRIN(NS)
        IF( TMZ.LE.SRC(1,1,MB) ) CYCLE
        IF( ISRM(NS).EQ.1 ) THEN
          DO N = 1,8
            SRX(N) = SRC(N,1,MB)
          ENDDO
        ELSE
          IFIND = 0
          DO M = 2,ISRM(NS)
            IF( TMZ.LE.SRC(1,M,MB) ) THEN
             DTSR = MIN( SRC(1,M,MB)-TMZ,DT )
             TFSR = (TMZ-0.5D+0*DTSR-SRC(1,M-1,MB))/
     &         (SRC(1,M,MB)-SRC(1,M-1,MB))
             DO N = 1,8
               SRX(N) = SRC(N,M-1,MB) + TFSR*(SRC(N,M,MB)-SRC(N,M-1,MB))
             ENDDO
             IFIND = 1
             EXIT
            ENDIF
          ENDDO
          IF( IFIND.EQ.0 ) CYCLE
        ENDIF
!
!---    Local node associated with processor-dependent (i.e., local)
!       source  ---
!
        N = ISRN(NS)
        IROW = IXP(N)
        ICOL = IXP(N)
        SORTX = 0.D+0
!
!---    Aqueous Volumetric Sink  ---
!
        IF( MOD(ISRT(NS),100).EQ.3 .AND. SRX(4).LT.0.D+0 ) THEN
          SORTX = -SRX(4)*YL(N,NSL)/(PORD(2,N)*SL(2,N))
!
!---    Gas Volumetric Sink  ---
!
        ELSEIF( MOD(ISRT(NS),100).EQ.4 .AND. SRX(4).LT.0.D+0 ) THEN
          SORTX = -SRX(4)*YG(N,NSL)/(PORD(2,N)*SG(2,N))
!
!---    Gas Volumetric Sink  ---
!
        ELSEIF( MOD(ISRT(NS),100).EQ.5 .AND. SRX(4).LT.0.D+0 ) THEN
          SORTX = -SRX(4)*YG(N,NSL)/(PORD(2,N)*SG(2,N))
!
!---    Aqueous Mass Sink  ---
!
        ELSEIF( MOD(ISRT(NS),100).EQ.7 .AND. SRX(4).LT.0.D+0 ) THEN
          SORTX = -SRX(4)*YL(N,NSL)/(RHOL(2,N)*PORD(2,N)*SL(2,N))
!
!---    Gas Mass Sink  ---
!
        ELSEIF( MOD(ISRT(NS),100).EQ.8 .AND. SRX(4).LT.0.D+0 ) THEN
          SORTX = -SRX(4)*YG(N,NSL)/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---    Gas Mass Sink  ---
!
        ELSEIF( MOD(ISRT(NS),100).EQ.9 .AND. SRX(4).LT.0.D+0 ) THEN
          SORTX = -SRX(4)*YG(N,NSL)/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---    Solute source  ---
!
        ELSEIF( ISRT(NS).EQ.-NSL ) THEN

          BUFFER = SRX(4)
          CALL lis_vector_set_value_f( 1,IROW,BUFFER,
     &      T_RHS_VEC,IERR )






!
!---    Solute density source  ---
!
        ELSEIF( ISRT(NS).EQ.-(NSL+NSOLU) ) THEN

          BUFFER = SRX(4)*VOL(N)
          CALL lis_vector_set_value_f( 1,IROW,BUFFER,
     &      T_RHS_VEC,IERR )






        ENDIF
!
!---    Load Jacobian matrix  ---
!
        MP = IXP(N)
        MA = 1
        MCOL = KLUC(MP-IEQC_OFFSET,MA)
        MA = MA + 1
        DLU(MCOL) = DLU(MCOL) + SORTX
!        BUFFER = SORTX
!        CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
!     &    T_MAT,IERR )
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORT_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP_CO2( NSL )
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
!     STOMPX-CO2
!
!     Calculates the aqueous- and gas-phase solute
!     mole fractions from user-specified partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE GRID
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/SPRP_CO2'
!
!---  Loop over local nodes, skipping inactive nodes
!     but including ghost cells  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        IF( IPCL(NSL).EQ.2 ) THEN
          XVS = RHOS(N)*PCSL(1,N,NSL)*(1.D+0-PORT(2,N))*SL(2,N)
        ELSE
          XVS = RHOS(N)*PCSL(1,N,NSL)*(1.D+0-PORT(2,N))
        ENDIF
        XVL = SL(2,N)*PORD(2,N)
        XVG = SG(2,N)*PORD(2,N)
!
!---    Constant gas-aqueous partition coefficient  ---
!
        IF( IPCGL(NSL).EQ.0 ) THEN
          PCGLX = PCGL(1,NSL)
!
!---    Temperature dependent gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.1 ) THEN
          TK = T(2,N)+TABS
          PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &      + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &      + PCGL(5,NSL)*TK**2 )
!
!---    Water-vapor 1 gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.2 ) THEN
          PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
        ENDIF
        PCGLX = MAX( PCGLX,1.D-20 )
        PCGLX = MIN( PCGLX,1.D+20 )
!
!---  Phase-volumetric concentration ratios  ---
!
        YVL = 1.D+0/(XVS + XVL + XVG*PCGLX)
        YVG = 1.D+0/((XVS + XVL)/PCGLX + XVG)
!
!---  Phase mole fractions  ---
!
        YL(N,NSL) = XVL*YVL
        YG(N,NSL) = XVG*YVG
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPRP_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTC( NSL )
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
!     STOMPX-CO2
!
!     Update solute concentrations
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE JACOB
      USE GRID
      USE FDVP
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
!---  Loop over local nodes, skipping inactive nodes
!     or ghost cells  ---
!
      NMD = 0
      DO N = 1,NFCGC(ID+1)

!
!---    Skip for inactive nodes or ghost cells  ---
!
        IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE







        NMD = NMD + 1
        C(N,NSL) = BLU(NMD)
        IF( ABS(C(N,NSL)).LT.EPSL ) C(N,NSL) = 0.D+0
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDTC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTC_GC( NSL )
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
!     STOMPX-CO2
!
!     Update solute concentrations on ghost cells
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE GRID
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
      INTEGER	STATUS(MPI_STATUS_SIZE)	
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDTC_GC'
      NPVX = 1
!
!---  Load sending buffer for bottom ghost cells for processors
!     with bottom ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(1,ID+1).GT.0 ) THEN
        MCS = 0
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(1,ID+1) = ',NCGC(1,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(1,ID+1)
          SBFB(NCS+1) = C(NLSGC(M+MCS),NSL)
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(1,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(1,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Bottom Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFB,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post Bottom Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending bottom ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(6,ID+1).GT.0 ) THEN
        NRCVX = NCGC(6,ID+1)*NPVX
        IDSNDX = NPGC(6,ID+1) - 1
        IDRCVX = NPGC(1,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Bottom Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFB,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post Bottom Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,5
          MCR = MCR + NCGC(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGC(6,ID+1) = ',NCGC(6,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(6,ID+1)
          C(NLRGC(M+MCR),NSL) = RBFB(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for south ghost cells for processors
!     with south ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(2,ID+1).GT.0 ) THEN
        MCS = NCGC(1,ID+1)
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(2,ID+1) = ',NCGC(2,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(2,ID+1)
          SBFS(NCS+1) = C(NLSGC(M+MCS),NSL)
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(2,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(2,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'South Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFS,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post South Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending south ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(5,ID+1).GT.0 ) THEN
        NRCVX = NCGC(5,ID+1)*NPVX
        IDSNDX = NPGC(5,ID+1) - 1
        IDRCVX = NPGC(2,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'South Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFS,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post South Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,4
          MCR = MCR + NCGC(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGC(5,ID+1) = ',NCGC(5,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(5,ID+1)
          C(NLRGC(M+MCR),NSL) = RBFS(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for west ghost cells for processors
!     with west ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(3,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,2
          MCS = MCS + NCGC(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(3,ID+1) = ',NCGC(3,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(3,ID+1)
!          IF( ID.EQ.2 ) 
!     &      PRINT *,'S2: ND(',NLSGC(M+MCS),') = ',ND(NLSGC(M+MCS))
          SBFW(NCS+1) = C(NLSGC(M+MCS),NSL)
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(3,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(3,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'West Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFW,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post West Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending west ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(4,ID+1).GT.0 ) THEN
        NRCVX = NCGC(4,ID+1)*NPVX
        IDSNDX = NPGC(4,ID+1) - 1
        IDRCVX = NPGC(3,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'West Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFW,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post West Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,3
          MCR = MCR + NCGC(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGC(4,ID+1) = ',NCGC(4,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(4,ID+1)
!          IF( ID.EQ.1 ) 
!     &      PRINT *,'R1: ND(',NLSGC(M+MCR),') = ',ND(NLSGC(M+MCR))
          C(NLRGC(M+MCR),NSL) = RBFW(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for east ghost cells for processors
!     with east ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(4,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,3
          MCS = MCS + NCGC(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(4,ID+1) = ',NCGC(4,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(4,ID+1)
          SBFE(NCS+1) = C(NLSGC(M+MCS),NSL)
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(4,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(4,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'East Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFE,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post East Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending east ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(3,ID+1).GT.0 ) THEN
        NRCVX = NCGC(3,ID+1)*NPVX
        IDSNDX = NPGC(3,ID+1) - 1
        IDRCVX = NPGC(4,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'East Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFE,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post East Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
        DO M = 1,2
          MCR = MCR + NCGC(M,ID+1)
        ENDDO
!        PRINT *,' MCR = ',MCR,' NCGC(3,ID+1) = ',NCGC(3,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(3,ID+1)
          C(NLRGC(M+MCR),NSL) = RBFE(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for north ghost cells for processors
!     with north ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(5,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,4
          MCS = MCS + NCGC(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(5,ID+1) = ',NCGC(5,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(5,ID+1)
          SBFN(NCS+1) = C(NLSGC(M+MCS),NSL)
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(5,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(5,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'North Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFN,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post North Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending north ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(2,ID+1).GT.0 ) THEN
        NRCVX = NCGC(2,ID+1)*NPVX
        IDSNDX = NPGC(2,ID+1) - 1
        IDRCVX = NPGC(5,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'North Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFN,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post North Receive: IERR = ',IERR,' ID = ',ID
        MCR = NCGC(1,ID+1)
!        PRINT *,' MCR = ',MCR,' NCGC(2,ID+1) = ',NCGC(2,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(2,ID+1)
          C(NLRGC(M+MCR),NSL) = RBFN(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Load sending buffer for top ghost cells for processors
!     with top ghost cells, and send buffer to receiving
!     processors  ---
!
      IF( NCGC(6,ID+1).GT.0 ) THEN
        MCS = 0
        DO M = 1,5
          MCS = MCS + NCGC(M,ID+1)
        ENDDO
        NCS = 0
!        PRINT *,' MCS = ',MCS,' NCGC(6,ID+1) = ',NCGC(6,ID+1),
!     &    ' ID = ',ID
        DO M = 1,NCGC(6,ID+1)
          SBFT(NCS+1) = C(NLSGC(M+MCS),NSL)
          NCS = NCS + NPVX
        ENDDO
        NSNDX = NCGC(6,ID+1)*NPVX
        IDSNDX = ID
        IDRCVX = NPGC(6,ID+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Top Send: NSNDX = ',NSNDX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_SEND( SBFT,NSNDX,MPI_REAL8,IDRCVX,ITAGX,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'Post Top Send: IERR = ',IERR,' ID = ',ID
      ENDIF
!
!---  Receive buffer from processors sending top ghost cells,
!     and unload receiving buffer  ---
!
      IF( NCGC(1,ID+1).GT.0 ) THEN
        NRCVX = NCGC(1,ID+1)*NPVX
        IDSNDX = NPGC(1,ID+1) - 1
        IDRCVX = NPGC(6,IDSNDX+1) - 1
        ITAGX = ITAG(IDSNDX+1,IDRCVX+1)
!        PRINT *,'Top Receive: NRCVX = ',NRCVX,' IDRCVX = ',IDRCVX,
!     &    ' IDSNDX = ',IDSNDX,' ITAGX = ',ITAGX,' ID = ',ID
        CALL MPI_RECV( RBFT,NRCVX,MPI_REAL8,IDSNDX,ITAGX,
     &    MPI_COMM_WORLD,STATUS,IERR )
!        PRINT *,'Post Top Receive: IERR = ',IERR,' ID = ',ID
        MCR = 0
!        PRINT *,' MCR = ',MCR,' NCGC(1,ID+1) = ',NCGC(1,ID+1),
!     &    ' ID = ',ID
        NCR = 0
        DO M = 1,NCGC(1,ID+1)
          C(NLRGC(M+MCR),NSL) = RBFT(NCR+1)
          NCR = NCR + NPVX
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDTC_GC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTCO( NSL )
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
!     STOMPX-CO2
!
!     Load old-time-step solute concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PROP
      USE GRID
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/UPDTCO'
!
!---  Loop over local nodes, skipping inactive nodes
!     but not ghost cells  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        CO(N,NSL) = C(N,NSL)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDTCO group  ---
!
      RETURN
      END



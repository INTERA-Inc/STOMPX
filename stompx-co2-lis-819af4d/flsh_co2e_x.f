!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_11( TX,PX,PGX,RHOLSX,RHOLAX,YLSX,
     &  XLSX,XLAX,CHMSGX )
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
!     Inputs:
!       aqueous-salt concentration
!       aqueous-CO2 concentration
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 1 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 YLSY(4),XLAY(4),XLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_11'
!
!---  Guess salt brine mass fraction  ---
!
      ISRX = 1
      CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
      YLSY(2) = RHOLSX/(RHOLWX+RHOLSX)
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
      
!
!---  Guess aqueous CO2 mass fraction  ---
!
      XLAY(2) = RHOLAX/(RHOLWX+RHOLSX+RHOLAX)
!
!---  Set primary variable increments  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAMX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAMX)
      DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSY(2)))
      DXLAX = SIGN((1.D-4*XLAMX),(5.D-1*XLAMX - XLAY(2)))
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
!
!---  Two-variable Newton-Raphson loop for (YLS, XLA)  ---
!
      DO NC = 1,32
        DO M = 2,4
          YLSY(M) = YLSY(2)
          XLAY(M) = XLAY(2)
          IF( M.EQ.3 ) YLSY(M) = YLSY(2) + DYLSX
          IF( M.EQ.4 ) XLAY(M) = XLAY(2) + DXLAX
          XLSX = MIN( YLSY(M),XLSMX )
          XLSY(M) = YLSY(M)
          XLSSX = XLSY(M)
          CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLSY(M) = YLSY(M) + (XLSSX-YLSY(M))*(XLAY(M)/XLAX)
          PVAX = XLAY(M)*PGAX/XLAX
          CALL SP_B( TX,XLSX,PSBX )
          IF( PVBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
            M_ERR(1) = CHMSGX(2)
            IF( N_DB.LT.0 ) THEN
              M_ERR(2) = ' at Boundary Surface: '
            ELSE
              M_ERR(2) = ' at Node: '
            ENDIF
            CALL PATH
            R_ERR = PVBX+PVAX
            I_ERR(1) = ABS(N_DB)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
            YLSX = YLSIX
            XLSX = XLSIX
            XLAX = XLAIX
            RETURN
          ENDIF
          CALL DENS_B( TX,PX,XLSX,RHOBX )
          CALL DENS_L( TX,RHOBX,XLAY(M),RHOLX )
          GX(M,1) = XLSY(M) - RHOLSX/RHOLX
          GX(M,2) = XLAY(M) - RHOLAX/RHOLX
        ENDDO
        RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
        RX(1,2) = (GX(4,1)-GX(2,1))/DXLAX
        RX(2,1) = (GX(3,2)-GX(2,2))/DYLSX
        RX(2,2) = (GX(4,2)-GX(2,2))/DXLAX
        RPX(1) = -GX(2,1)
        RPX(2) = -GX(2,2)
        CYLSX = (RPX(2)-RPX(1)*RX(2,2)/(RX(1,2)+SMALL))/
     &    (RX(2,1)-RX(1,1)*RX(2,2)/(RX(1,2)+SMALL))
        CXLAX = (RPX(2)-RPX(1)*RX(2,1)/(RX(1,1)+SMALL))/
     &    (RX(2,2)-RX(1,2)*RX(2,1)/(RX(1,1)+SMALL))
        YLSY(2) = MAX( YLSY(2)+CYLSX,0.D+0 )
        XLAY(2) = MAX( XLAY(2)+CXLAX,0.D+0 )
        IF( ABS(CYLSX).LE.(1.D-9*XLSMX) .AND.
     &    ABS(CXLAX).LE.(1.D-9*XLAMX) ) EXIT
      ENDDO
      IF( NC.GT.32 ) THEN
        M_ERR(1) = CHMSGX(1)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = RHOLSX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
!
!---  Limit salt mass fraction  ---
!
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAX)
      PVAX = XLAY(2)*PGAX/XLAX
      CALL SP_B( TX,XLSX,PSBX )
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_11 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_12( TX,PX,PGX,RHOLSX,PHILAX,YLSX,XLSX,XLAX,
     &  CHMSGX )
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
!     Inputs:
!       aqueous-salt concentration
!       aqueous-CO2 relative saturation
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 YLSY(4),XLAY(4),XLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_12'
!
!---  Guess salt brine mass fraction  ---
!
      ISRX = 1
      CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
      YLSY(2) = RHOLSX/(RHOLWX+RHOLSX)
      YLSX = YLSY(2)
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
!
!---  Guess aqueous CO2 mass fraction  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAMX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAMX)
      CALL SOL_LS( TX,XLSMX )
      DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSY(2)))
      DXLAX = SIGN((1.D-4*XLAMX),(5.D-1*XLAMX - XLAY(2)))
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
      CALL SP_B( TX,XLSX,PSBX )
!
!---  Guess aqueous-CO2 mass fraction  ---
!
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLAY(2) = PHILAX*XLAX
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAX)
      XLAX = XLAY(2)
      XLSX = XLSY(2)
      PVAX = MIN( XLAY(2)/XLAX,1.D+0 )*PGAX
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(1)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
!
!---  Two-variable Newton-Raphson loop (YLS, XLA)  ---
!
      DO NC = 1,32
        DO M = 2,4
          YLSY(M) = YLSY(2)
          XLAY(M) = XLAY(2)
          IF( M.EQ.3 ) YLSY(M) = YLSY(2) + DYLSX
          IF( M.EQ.4 ) XLAY(M) = XLAY(2) + DXLAX
          XLSX = MIN( YLSY(M),XLSMX )
          XLSY(M) = YLSY(M)
          CALL SP_B( TX,XLSX,PSBX )
          PVBX = PSBX
          XLSSX = XLSY(M)
          CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLSY(M) = YLSY(M) + (XLSSX-YLSY(M))*(XLAY(M)/XLAX)
          PVAX = MIN( XLAY(M)/XLAX,1.D+0 )*PGAX
          IF( PVBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
            M_ERR(1) = CHMSGX(2)
            IF( N_DB.LT.0 ) THEN
              M_ERR(2) = ' at Boundary Surface: '
            ELSE
              M_ERR(2) = ' at Node: '
            ENDIF
            CALL PATH
            R_ERR = PVBX+PVAX
            I_ERR(1) = ABS(N_DB)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
            YLSX = YLSIX
            XLSX = XLSIX
            XLAX = XLAIX
            RETURN
          ENDIF
          CALL DENS_B( TX,PX,XLSX,RHOBX )
          CALL DENS_L( TX,RHOBX,XLAY(M),RHOLX )
          GX(M,1) = XLSY(M) - RHOLSX/RHOLX
          GX(M,2) = XLAY(M) - PHILAX*XLAX
        ENDDO
        RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
        RX(1,2) = (GX(4,1)-GX(2,1))/DXLAX
        RX(2,1) = (GX(3,2)-GX(2,2))/DYLSX
        RX(2,2) = (GX(4,2)-GX(2,2))/DXLAX
        RPX(1) = -GX(2,1)
        RPX(2) = -GX(2,2)
        CYLSX = (RPX(2)-RPX(1)*RX(2,2)/(RX(1,2)+SMALL))/
     &    (RX(2,1)-RX(1,1)*RX(2,2)/(RX(1,2)+SMALL))
        CXLAX = (RPX(2)-RPX(1)*RX(2,1)/(RX(1,1)+SMALL))/
     &    (RX(2,2)-RX(1,2)*RX(2,1)/(RX(1,1)+SMALL))
        YLSY(2) = MAX( YLSY(2)+CYLSX,0.D+0 )
        XLAY(2) = MAX( XLAY(2)+CXLAX,0.D+0 )
        IF( ABS(CYLSX).LE.(1.D-9*XLSMX) .AND.
     &    ABS(CXLAX).LE.(1.D-9*XLAMX) ) EXIT
      ENDDO
      IF( NC.GT.32 ) THEN
        M_ERR(1) = CHMSGX(1)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = RHOLSX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
!
!---  Limit salt mass fraction  ---
!
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAX)
      PVAX = MIN( XLAY(2)/XLAX,1.D+0 )*PGAX
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_12 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_13( TX,PX,PGX,RHOLSX,YLSX,XLSX,XLAX,CHMSGX )
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
!     Inputs:
!       aqueous-salt concentration
!       aqueous-CO2 mass fraction
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE GRID
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 YLSY(4),XLSY(4),XLAY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_13'
      XLAY(2) = XLAX
!
!---  Guess salt brine mass fraction  ---
!
      ISRX = 1
      CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
      YLSY(2) = RHOLSX/(RHOLWX+RHOLSX)
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
!
!---  Set primary variable increments  ---
!
      DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSY(2)))
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
!
!---  Single-variable Newton-Raphson loop (YLS)  ---
!
      DO NC = 1,32
        DO M = 2,3
          YLSY(M) = YLSY(2)
          IF( M.EQ.3 ) YLSY(M) = YLSY(2) + DYLSX
          XLSX = MIN( YLSY(M),XLSMX )
          XLSY(M) = YLSY(M)
          CALL SP_B( TX,XLSX,PSBX )
          PVBX = PSBX
          XLSSX = XLSY(M)
          CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLSY(M) = YLSY(M) + (XLSSX-YLSY(M))*(XLAY(2)/XLAX)
          PVAX = MIN( XLAY(2)/XLAX,1.D+0 )*PGAX
          CALL SP_B( TX,XLSX,PSBX )
          IF( PVBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
            M_ERR(1) = CHMSGX(2)
            IF( N_DB.LT.0 ) THEN
              M_ERR(2) = ' at Boundary Surface: '
            ELSE
              M_ERR(2) = ' at Node: '
            ENDIF
            CALL PATH
            R_ERR = PVBX+PVAX
            I_ERR(1) = ABS(N_DB)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
            YLSX = YLSIX
            XLSX = XLSIX
            XLAX = XLAIX
            RETURN
          ENDIF
          CALL DENS_B( TX,PX,XLSX,RHOBX )
          CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
          GX(M,1) = XLSY(M) - RHOLSX/RHOLX
        ENDDO
        RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
        RPX(1) = -GX(2,1)
        CYLSX = RPX(1)/RX(1,1)
        YLSY(2) = MAX( YLSY(2)+CYLSX,0.D+0 )
        IF( ABS(CYLSX).LE.(1.D-9*XLSMX) ) EXIT
      ENDDO
      IF( NC.GT.32 ) THEN
        M_ERR(1) = CHMSGX(1)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = RHOLSX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
!
!---  Limit salt mass fraction  ---
!
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAX)
      PVAX = MIN( XLAY(2)/XLAX,1.D+0 )*PGAX
      CALL SP_B( TX,XLSX,PSBX )
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_13 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_21( TX,PX,PGX,PHILSX,RHOLAX,YLSX,XLSX,XLAX,
     &  CHMSGX )
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
!     Inputs:
!       aqueous-salt relative saturation
!       aqueous-CO2 concentration
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE GRID
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XLSY(4),YLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_21'
!
!---  Guess salt brine mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      YLSY(2) = PHILSX*XLSMX
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
!
!---  Guess aqueous CO2 mass fraction  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      ISRX = 1
      CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
      RHOLSX = XLSX*RHOLWX
      XLAY(2) = RHOLAX/(RHOLWX+RHOLSX+RHOLAX)
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAMX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAMX)
!
!---  Set primary variable increments  ---
!
      DXLAX = SIGN((1.D-4*XLAMX),(5.D-1*XLAMX - XLAY(2)))
!
!---  Single-variable Newton-Raphson loop (XLA)  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
      DO NC = 1,32
        DO M = 2,3
          XLAY(M) = XLAY(2)
          IF( M.EQ.3 ) XLAY(M) = XLAY(2) + DXLAX
          XLSX = MIN( YLSY(2),XLSMX )
          XLSY(2) = XLSX
          XLSSX = XLSY(2)
          CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(M)/XLAX)
          PVAX = XLAY(M)*PGAX/XLAX
          CALL SP_B( TX,XLSX,PSBX )
!
!---      Check for unsaturated conditions  ---
!
          IF( PVBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
            M_ERR(1) = CHMSGX(2)
            IF( N_DB.LT.0 ) THEN
              M_ERR(2) = ' at Boundary Surface: '
            ELSE
              M_ERR(2) = ' at Node: '
            ENDIF
            CALL PATH
            R_ERR = PVBX+PVAX
            I_ERR(1) = ABS(N_DB)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
            YLSX = YLSIX
            XLSX = XLSIX
            XLAX = XLAIX
            RETURN
          ENDIF
          CALL DENS_B( TX,PX,XLSX,RHOBX )
          CALL DENS_L( TX,RHOBX,XLAY(M),RHOLX )
          GX(M,1) = XLAY(M) - RHOLAX/RHOLX
        ENDDO
        RX(1,1) = (GX(3,1)-GX(2,1))/DXLAX
        RPX(1) = -GX(2,1)
        CXLAX = RPX(1)/RX(1,1)
        XLAY(2) = MAX( XLAY(2)+CXLAX,0.D+0 )
        IF( ABS(CXLAX).LE.(1.D-9*XLAMX) ) EXIT
      ENDDO
!
!---  Convergence failure  ---
!
      IF( NC.GT.32 ) THEN
        M_ERR(1) = CHMSGX(1)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = RHOLAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        PVAX = MAX( PGX-PVBX,0.D+0 )
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_21 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_22( TX,PX,PGX,PHILSX,PHILAX,YLSX,XLSX,XLAX,
     &  CHMSGX )
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
!     Inputs:
!       aqueous-salt relative saturation
!       aqueous-CO2 relative saturation
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE GRID
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XLSY(4),YLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_22'
!
!---  Guess salt brine mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      YLSY(2) = PHILSX*XLSMX
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
!
!---  Guess aqueous CO2 mass fraction  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      ISRX = 1
      CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
      RHOLSX = XLSX*RHOLWX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAMX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLAY(2) = PHILAX*XLAMX
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAMX)
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
!
!---  Check for unsaturated conditions  ---
!
      PVAX = MIN( XLAY(2)/XLAMX,1.D+0 )*PGAX
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
!
!---  Set primary variable increments  ---
!
      DXLAX = SIGN((1.D-4*XLAMX),(5.D-1*XLAMX - XLAY(2)))
!
!---  Single-variable Newton-Raphson loop (XLA)  ---
!
      DO NC = 1,32
        DO M = 2,3
          XLAY(M) = XLAY(2)
          IF( M.EQ.3 ) XLAY(M) = XLAY(2) + DXLAX
          XLSX = MIN( YLSY(2),XLSMX )
          XLSY(2) = XLSX
          CALL SP_B( TX,XLSX,PSBX )
          PVBX = PSBX
          XLSSX = XLSY(2)
          CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(M)/XLAX)
          PVAX = MIN( XLAY(M)/XLAX,1.D+0 )*PGAX
          IF( PVBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
            M_ERR(1) = CHMSGX(2)
            IF( N_DB.LT.0 ) THEN
              M_ERR(2) = ' at Boundary Surface: '
            ELSE
              M_ERR(2) = ' at Node: '
            ENDIF
            CALL PATH
            R_ERR = PVBX+PVAX
            I_ERR(1) = ABS(N_DB)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
            YLSX = YLSIX
            XLSX = XLSIX
            XLAX = XLAIX
            RETURN
          ENDIF
          CALL DENS_B( TX,PX,XLSX,RHOBX )
          CALL DENS_L( TX,RHOBX,XLAY(M),RHOLX )
          GX(M,1) = XLAY(M) - PHILAX*XLAX
        ENDDO
        RX(1,1) = (GX(3,1)-GX(2,1))/DXLAX
        RPX(1) = -GX(2,1)
        CXLAX = RPX(1)/RX(1,1)
        XLAY(2) = MAX( XLAY(2)+CXLAX,0.D+0 )
        IF( ABS(CXLAX).LE.(1.D-9*XLAMX) ) EXIT
      ENDDO
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        M_ERR(1) = CHMSGX(1)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PHILAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_22 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_23( TX,PX,PGX,PHILSX,YLSX,XLSX,XLAX,CHMSGX )
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
!     Inputs:
!       aqueous-salt relative saturation
!       aqueous-CO2 mass fraction
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE GRID
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XLSY(4),YLSY(4)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_23'
      XLAY(2) = XLAX
      CALL SOL_LS( TX,XLSMX )
      YLSY(2) = PHILSX*XLSMX
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAX)
      PVAX = MIN( XLAY(2)/XLAX,1.D+0 )*PGAX
      CALL SP_B( TX,XLSX,PSBX )
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_23 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_31( TX,PX,PGX,RHOLAX,YLSX,XLSX,XLAX,CHMSGX )
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
!     Inputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 concentration
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XLSY(4),YLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_31'
!
!---  Guess salt brine mass fraction  ---
!
      XLSZ = YLSX
      YLSY(2) = YLSX
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
!
!---  Guess aqueous CO2 mass fraction  ---
!
      CALL DENS_B( TX,PX,XLSX,RHOBX )
      XLAY(2) = RHOLAX/(RHOBX+RHOLAX)
!
!---  Set primary variable increments  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAMX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAMX)
      DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSY(2)))
      DXLAX = SIGN((1.D-4*XLAMX),(5.D-1*XLAMX - XLAY(2)))
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
!
!---  Two-variable Newton-Raphson loop (YLS, XLA)  ---
!
      DO NC = 1,32
        DO M = 2,4
          YLSY(M) = YLSY(2)
          XLAY(M) = XLAY(2)
          IF( M.EQ.3 ) YLSY(M) = YLSY(2) + DYLSX
          IF( M.EQ.4 ) XLAY(M) = XLAY(2) + DXLAX
          XLSX = MIN( YLSY(M),XLSMX )
          XLSY(M) = YLSY(M)
          XLSSX = XLSY(M)
          CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLSY(M) = YLSY(M) + (XLSSX-YLSY(M))*(XLAY(M)/XLAX)
          PVAX = XLAY(M)*PGAX/XLAX
          CALL SP_B( TX,XLSY(M),PSBX )
          IF( PVBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
            M_ERR(1) = CHMSGX(2)
            IF( N_DB.LT.0 ) THEN
              M_ERR(2) = ' at Boundary Surface: '
            ELSE
              M_ERR(2) = ' at Node: '
            ENDIF
            CALL PATH
            R_ERR = PVBX+PVAX
            I_ERR(1) = ABS(N_DB)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
            YLSX = YLSIX
            XLSX = XLSIX
            XLAX = XLAIX
            RETURN
          ENDIF
          CALL DENS_B( TX,PX,XLSY(M),RHOBX )
          CALL DENS_L( TX,RHOBX,XLAY(M),RHOLX )
          GX(M,1) = XLSY(M) - XLSZ
          GX(M,2) = XLAY(M) - RHOLAX/RHOLX
        ENDDO
        RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
        RX(1,2) = (GX(4,1)-GX(2,1))/DXLAX
        RX(2,1) = (GX(3,2)-GX(2,2))/DYLSX
        RX(2,2) = (GX(4,2)-GX(2,2))/DXLAX
        RPX(1) = -GX(2,1)
        RPX(2) = -GX(2,2)
        CYLSX = (RPX(2)-RPX(1)*RX(2,2)/(RX(1,2)+SMALL))/
     &    (RX(2,1)-RX(1,1)*RX(2,2)/(RX(1,2)+SMALL))
        CXLAX = (RPX(2)-RPX(1)*RX(2,1)/(RX(1,1)+SMALL))/
     &    (RX(2,2)-RX(1,2)*RX(2,1)/(RX(1,1)+SMALL))
        YLSY(2) = MAX( YLSY(2)+CYLSX,0.D+0 )
        XLAY(2) = MAX( XLAY(2)+CXLAX,0.D+0 )
        IF( ABS(CYLSX).LE.(1.D-9*XLSMX) .AND.
     &    ABS(CXLAX).LE.(1.D-9*XLAMX) ) EXIT
      ENDDO
      IF( NC.GT.32 ) THEN
        M_ERR(1) = CHMSGX(1)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = YLSX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
!
!---  Limit salt mass fraction  ---
!
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAX)
      PVAX = XLAY(2)*PGAX/XLAX
      CALL SP_B( TX,XLSY(2),PSBX )
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_31 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_32( TX,PX,PGX,PHILAX,YLSX,XLSX,XLAX,CHMSGX )
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
!     Inputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 relative saturation
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XLSY(4),YLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_32'
!
!---  Guess salt brine mass fraction  ---
!
      XLSZ = YLSX
      YLSY(2) = YLSX
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
!
!---  Guess aqueous CO2 mass fraction  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAMX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLAY(2) = PHILAX*XLAMX
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAMX)
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
!
!---  Check for unsaturated conditions  ---
!
      PVAX = MIN( XLAY(2)/XLAMX,1.D+0 )*PGAX
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
!
!---  Set primary variable increments  ---
!
      DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSY(2)))
      DXLAX = SIGN((1.D-4*XLAMX),(5.D-1*XLAMX - XLAY(2)))
!
!---  Two-variable Newton-Raphson loop (YLS, XLA)  ---
!
      DO NC = 1,32
      DO 40 M = 2,4
        YLSY(M) = YLSY(2)
        XLAY(M) = XLAY(2)
        IF( M.EQ.3 ) YLSY(M) = YLSY(2) + DYLSX
        IF( M.EQ.4 ) XLAY(M) = XLAY(2) + DXLAX
        XLSX = MIN( YLSY(M),XLSMX )
        XLSY(M) = YLSY(M)
        CALL SP_B( TX,XLSX,PSBX )
        PVBX = PSBX
        XLSSX = XLSY(M)
        CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &    XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &    XMLSX,XMLWX )
        XLSY(M) = YLSY(M) + (XLSSX-YLSY(M))*(XLAY(M)/XLAX)
        PVAX = MIN( XLAY(M)/XLAX,1.D+0 )*PGAX
        IF( PVBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
          M_ERR(1) = CHMSGX(2)
          IF( N_DB.LT.0 ) THEN
            M_ERR(2) = ' at Boundary Surface: '
          ELSE
            M_ERR(2) = ' at Node: '
          ENDIF
          CALL PATH
          R_ERR = PVBX+PVAX
          I_ERR(1) = ABS(N_DB)
          I_ERR(2) = 1
          I_ERR(3) = 2
          I_ERR(4) = ID
          YLSX = YLSIX
          XLSX = XLSIX
          XLAX = XLAIX
          RETURN
        ENDIF
        CALL DENS_B( TX,PX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAY(M),RHOLX )
        GX(M,1) = XLSY(M) - XLSZ
        GX(M,2) = XLAY(M) - PHILAX*XLAX
   40 CONTINUE
      RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
      RX(1,2) = (GX(4,1)-GX(2,1))/DXLAX
      RX(2,1) = (GX(3,2)-GX(2,2))/DYLSX
      RX(2,2) = (GX(4,2)-GX(2,2))/DXLAX
      RPX(1) = -GX(2,1)
      RPX(2) = -GX(2,2)
      CYLSX = (RPX(2)-RPX(1)*RX(2,2)/(RX(1,2)+SMALL))/
     &  (RX(2,1)-RX(1,1)*RX(2,2)/(RX(1,2)+SMALL))
      CXLAX = (RPX(2)-RPX(1)*RX(2,1)/(RX(1,1)+SMALL))/
     &  (RX(2,2)-RX(1,2)*RX(2,1)/(RX(1,1)+SMALL))
      YLSY(2) = MAX( YLSY(2)+CYLSX,0.D+0 )
      XLAY(2) = MAX( XLAY(2)+CXLAX,0.D+0 )
      IF( ABS(CYLSX).LE.(1.D-9*XLSMX) .OR.
     &  ABS(CXLAX).LE.(1.D-9*XLAMX) ) EXIT
      ENDDO
      IF( NC.GT.32 ) THEN
        M_ERR(1) = CHMSGX(1)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = YLSX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
!
!---  Limit salt mass fraction  ---
!
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAX)
      PVAX = MIN( XLAY(2)/XLAX,1.D+0 )*PGAX
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  End of FLSH_32 group  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_33( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
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
!     Inputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE GRID
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XLSY(4),YLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_33'
      XLAY(2) = XLAX
!
!---  Guess salt brine mass fraction  ---
!
      XLSZ = YLSX
      YLSY(2) = YLSX
      XLAY(2) = XLAX
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
!
!---  Set primary variable increments  ---
!
      DYLSX = SIGN((1.D-4*XLSMX),(5.D-1*XLSMX - YLSY(2)))
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
!
!---  Single-variable Newton-Raphson loop (YLS)  ---
!
      DO NC = 1,32
        DO M = 2,3
          YLSY(M) = YLSY(2)
          IF( M.EQ.3 ) YLSY(M) = YLSY(2) + DYLSX
          XLSX = MIN( YLSY(M),XLSMX )
          XLSY(M) = YLSY(M)
          CALL SP_B( TX,XLSX,PSBX )
          PVBX = PSBX
          XLSSX = XLSY(M)
          CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLSY(M) = YLSY(M) + (XLSSX-YLSY(M))*(XLAY(2)/XLAX)
          PVAX = MIN( XLAY(2)/XLAX,1.D+0 )*PGAX
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( TX,PX,XLSX,RHOBX )
          CALL DENS_L( TX,RHOBX,XLAY(2),RHOLX )
          GX(M,1) = XLSY(M) - XLSZ
        ENDDO
        RX(1,1) = (GX(3,1)-GX(2,1))/DYLSX
        RPX(1) = -GX(2,1)
        CYLSX = RPX(1)/RX(1,1)
        YLSY(2) = MAX( YLSY(2)+CYLSX,0.D+0 )
        IF( ABS(CYLSX).LE.(1.D-9*XLSMX) ) EXIT
      ENDDO
      IF( NC.GT.32 ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = YLSX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        PVAX = MAX( PGX-PVBX,0.D+0 )
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_33 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_41( TX,PX,PGX,RHOLAX,YLSX,XLSX,XLAX,CHMSGX )
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
!     Inputs:
!       aqueous-salt molality (input as brine-salt mass fraction)
!       aqueous-CO2 concentration
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XLSY(4),YLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_41'
      YLSY(2) = YLSX*WTMS/(1.D+0 + YLSX*WTMS)
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
!
!---  Guess aqueous CO2 mass fraction  ---
!
      CALL DENS_B( TX,PX,XLSX,RHOBX )
      XLAY(2) = RHOLAX/(RHOBX+RHOLAX)
!
!---  Set primary variable increments  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAMX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAMX)
      DXLAX = SIGN((1.D-4*XLAMX),(5.D-1*XLAMX - XLAY(2)))
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
!
!---  Single-variable Newton-Raphson loop (XLA)  ---
!
      DO NC = 1,32
        DO M = 2,3
          XLAY(M) = XLAY(2)
          IF( M.EQ.3 ) XLAY(M) = XLAY(2) + DXLAX
          XLSX = MIN( YLSY(2),XLSMX )
          XLSY(2) = XLSX
          CALL SP_B( TX,XLSX,PSBX )
          PVBX = PSBX
          XLSSX = XLSY(2)
          CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(M)/XLAX)
          PVAX = XLAY(M)*PGAX/XLAX
          CALL SP_B( TX,XLSX,PSBX )
          IF( PVBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
            M_ERR(1) = CHMSGX(2)
            IF( N_DB.LT.0 ) THEN
              M_ERR(2) = ' at Boundary Surface: '
            ELSE
              M_ERR(2) = ' at Node: '
            ENDIF
            CALL PATH
            R_ERR = PVBX+PVAX
            I_ERR(1) = ABS(N_DB)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
            YLSX = YLSIX
            XLSX = XLSIX
            XLAX = XLAIX
            RETURN
          ENDIF
          CALL DENS_B( TX,PX,XLSX,RHOBX )
          CALL DENS_L( TX,RHOBX,XLAY(M),RHOLX )
          GX(M,1) = XLAY(M) - RHOLAX/RHOLX
        ENDDO
        RX(1,1) = (GX(3,1)-GX(2,1))/DXLAX
        RPX(1) = -GX(2,1)
        CXLAX = RPX(1)/(RX(1,1)+SMALL)
        XLAY(2) = MAX( XLAY(2)+CXLAX,0.D+0 )
        IF( ABS(CXLAX).LE.(1.D-9*XLAMX) ) EXIT
      ENDDO
      IF( NC.GT.32 ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_41 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_42( TX,PX,PGX,PHILAX,YLSX,XLSX,XLAX,CHMSGX )
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
!     Inputs:
!       aqueous-salt molality (input as brine-salt mass fraction)
!       aqueous-CO2 relative saturation
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XLSY(4),YLSY(4)
      REAL*8 GX(4,2),RX(2,2),RPX(2)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_42'
      YLSY(2) = YLSX*WTMS/(1.D+0 + YLSX*WTMS)
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
!
!---  Guess aqueous CO2 mass fraction  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAMX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLAY(2) = PHILAX*XLAMX
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAMX)
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
!
!---  Check for unsaturated conditions  ---
!
      PVAX = MIN( XLAY(2)/XLAMX,1.D+0 )*PGAX
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
!
!---  Set primary variable increments  ---
!
      DXLAX = SIGN((1.D-4*XLAMX),(5.D-1*XLAMX - XLAY(2)))
!
!---  Single-variable Newton-Raphson loop (XLA)  ---
!
      DO NC = 1,32
        DO M = 2,3
          XLAY(M) = XLAY(2)
          IF( M.EQ.3 ) XLAY(M) = XLAY(2) + DXLAX
          XLSX = MIN( YLSY(2),XLSMX )
          XLSY(2) = XLSX
          CALL SP_B( TX,XLSX,PSBX )
          PVBX = PSBX
          XLSSX = XLSY(2)
          CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &      XGAX,XGWX,XLAMX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &      XMLSX,XMLWX )
          XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(M)/XLAMX)
          PVAX = MIN( XLAY(M)/XLAMX,1.D+0 )*PGAX
          IF( PVBX+PVAX.GT.PGX .AND. M.EQ.2 ) THEN
            M_ERR(1) = CHMSGX(2)
            IF( N_DB.LT.0 ) THEN
              M_ERR(2) = ' at Boundary Surface: '
            ELSE
              M_ERR(2) = ' at Node: '
            ENDIF
            CALL PATH
            R_ERR = PVBX+PVAX
            I_ERR(1) = ABS(N_DB)
            I_ERR(2) = 1
            I_ERR(3) = 2
            I_ERR(4) = ID
            YLSX = YLSIX
            XLSX = XLSIX
            XLAX = XLAIX
            RETURN
          ENDIF
          GX(M,2) = XLAY(M) - PHILAX*XLAMX
        ENDDO
        RX(1,1) = (GX(3,1)-GX(2,1))/DXLAX
        RPX(1) = -GX(2,1)
        CXLAX = RPX(1)/(RX(1,1)+SMALL)
        XLAY(2) = MAX( XLAY(2)+CXLAX,0.D+0 )
        IF( ABS(CXLAX).LE.(1.D-9*XLAMX) ) EXIT
      ENDDO
      IF( NC.GT.32 ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF      
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_42 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLSH_43( TX,PX,PGX,YLSX,XLSX,XLAX,CHMSGX )
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
!     Inputs:
!       aqueous-salt molality (input as brine-salt mass fraction)
!       aqueous-CO2 mass fraction
!     Outputs:
!       aqueous-salt mass fraction
!       aqueous-CO2 mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XLAY(4),XLSY(4),YLSY(4)
      CHARACTER*132 CHMSGX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLSH_43'
      XLAY(2) = XLAX
      YLSY(2) = YLSX*WTMS/(1.D+0 + YLSX*WTMS)
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSY(2),XLSMX )
      XLSY(2) = XLSX
!
!---  Store initial guesses  ---
!
      YLSIX = YLSY(2)
      XLSIX = XLSY(2)
      XLAIX = XLAY(2)
!
!---  Guess aqueous CO2 mass fraction  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PVBX = PSBX
      XLSSX = XLSY(2)
      CALL EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAMX,XLSSX,XLWX,XMGAX,XMGWX,XMLAX,
     &  XMLSX,XMLWX )
      XLSY(2) = YLSY(2) + (XLSSX-YLSY(2))*(XLAY(2)/XLAMX)
!
!---  Check for unsaturated conditions  ---
!
      PVAX = MIN( XLAY(2)/XLAMX,1.D+0 )*PGAX
      IF( PVBX+PVAX.GT.PGX ) THEN
        M_ERR(1) = CHMSGX(2)
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = PVBX+PVAX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
        YLSX = YLSIX
        XLSX = XLSIX
        XLAX = XLAIX
        RETURN
      ENDIF
      XLSX = XLSY(2)
      YLSX = YLSY(2)
      XLAX = XLAY(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLSH_43 group  ---
!
      RETURN
      END



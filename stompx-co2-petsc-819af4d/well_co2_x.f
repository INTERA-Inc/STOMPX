!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE EQUIL_COUP_WELL( NCW )
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
!
!     STOMPX-CO2
!
!     Equilibrate coupled-well pressure with formation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE HYST
      USE GRID
      USE FDVP
      USE COUP_WELL
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
      SUB_LOG(ISUB_LOG) = '/EQUIL_COUP_WELL'
!
!---  Injection well, equilibrate with first well node  ---
!
      IF( IT_CW(NCW).GT.0 ) THEN
        N = IWN_CW(ID_CW(3,NCW))
        XCWX = XTP_CW(1,ID_CW(1,NCW))
        YCWX = YTP_CW(1,ID_CW(1,NCW))
        ZCWX = ZTP_CW(1,ID_CW(1,NCW))
!
!---  Withdrawl well, equilibrate with last well node  ---
!
      ELSE
        N = IWN_CW(ID_CW(4,NCW))
        XCWX = XTP_CW(2,ID_CW(2,NCW))
        YCWX = YTP_CW(2,ID_CW(2,NCW))
        ZCWX = ZTP_CW(2,ID_CW(2,NCW))
      ENDIF
      IDLX = -1
      IF( N.NE.0 ) THEN
!
!---    Aqueous unsaturated conditions  ---
!
        IF( SG(2,N)-SGT(2,N).GT.EPSL ) THEN
          IDLX = ID
          P_CW(2,NCW) = PG(2,N) - (ZCWX-ZP(N))*GRAV*RHOG(2,N)
!
!---    Aqueous saturated conditions  ---
!
        ELSE
          IDLX = ID
          P_CW(2,NCW) = PL(2,N) - (ZCWX-ZP(N))*GRAV*RHOL(2,N)
        ENDIF
      ENDIF
!
!---  Identify processor with defining well node  ---
!
      CALL MPI_ALLREDUCE( IDLX,IDX,1,MPI_INTEGER,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
!
!---  Broadcast coupled-well pressure from processor with
!     defining well node  ---
!
      CALL MPI_BCAST( P_CW(2,NCW),1,MPI_REAL8,IDX,MPI_COMM_WORLD,IERR )
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of EQUIL_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE COMST_COUP_WELL
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
!
!     STOMPX-CO2
!
!     Communicate to all processors the temperature of the 
!     field nodes containing well nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE HYST
      USE GRID
      USE FDVP
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(1:LWN_CW) :: VARX,VARGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/COMST_COUP_WELL'
!
!---  Loop over coupled wells ---
!
      DO NCW = 1,N_CW
!
!---    Loop over coupled-well nodes  ---
!
        DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          VARX(NWN) = -273.15D+0
          N = IWN_CW(NWN)
          IF( N.EQ.0 ) CYCLE
          IF( IGHC(N).EQ.1 ) CYCLE
          VARX(NWN) = T(2,N)
        ENDDO
      ENDDO
!
!---  Communicate field state to all processors  ---
!
      NVAR = LWN_CW
      CALL MPI_ALLREDUCE( VARX,VARGX,NVAR,MPI_REAL8,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
!
!---  Loop over coupled wells ---
!
      DO NCW = 1,N_CW
!
!---    Loop over coupled-well nodes  ---
!
        DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          TF_CW(NWN) = VARGX(NWN)
        ENDDO
      ENDDO
!
!---  Initialize coupled-well pressure to be in hydrostatic
!     1 with local pressure, and set the old time
!     value of the coupled-well pressure  ---
!
      DO NCW = 1,N_CW
        CALL EQUIL_COUP_WELL( NCW )
        P_CW(1,NCW) = P_CW(2,NCW)
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of COMST_COUP_WELL group  ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLUX_COUP_WELL
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
!
!     STOMPX-CO2
!
!     Mass flux of water and CO2 between coupled-well nodes and 
!     field nodes.
!
!     CO2 mass balance residuals for injection type coupled wells.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 17 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE HYST
      USE GRID
      USE FDVP
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XPX(2),YPX(2),ZPX(2)
      REAL*8 VAR_CWX(6)
      REAL*8, DIMENSION(1:(LUK_CW+1)) :: RSL_CWX
      REAL*8, DIMENSION(1:4) :: QML_CWX,QMG_CWX
      INTEGER, DIMENSION(1:(LUK+2)) :: MCW,MFD
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLUX_COUP_WELL'
      DO M = 1,ISVC+2
        IF( M.NE.ISVC+2 ) THEN
          MCW(M) = 2
        ELSE
          MCW(M) = 3
        ENDIF
        IF( M.NE.ISVC+2 ) THEN
          MFD(M) = M+1
        ELSE
          MFD(M) = 2
        ENDIF
      ENDDO       
!
!---  Initialize coupled-well parameters ---
!
      DO NCW = 1,N_CW
!
!---    Zero coupled-well fluxes ---
!
        QM_CW(1,NCW) = 0.D+0
        QM_CW(3,NCW) = 0.D+0
        QM_CW(7,NCW) = 0.D+0
!
!---    Flow controlled well ---
!
        ID_CW(8,NCW) = 0
!
!---    Loop over coupled-well nodes  ---
!
        DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---      Zero volumetric injection well fluxes
!
!         Q_CW(1,NWN) - total volumetric flux, m^3/s
!         Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!         Q_CW(3,NWN) - gas volumetric flux, m^3/s
!
          DO M = 1,3
            Q_CW(M,NWN) = 0.D+0
          ENDDO
!
!---      Loop over increment indices  ---
!
          DO M = 1,ISVC+2
            FXA_CW(M,NWN) = 0.D+0
            FXS_CW(M,NWN) = 0.D+0
            FXW_CW(M,NWN) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
!
!---  Loop over coupled wells ---
!
      CWLOOP: DO NCW = 1,N_CW
        ICHK_CWX = 0
        DQ_CWX = 1.D-6
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        NCT = 0
        DO M = 1,IM_CW(NCW)
          NCT = NCT + IMP_CW(M,NCW)
        ENDDO
        IF( ICC_CW(NCW).EQ.1 ) TMZ = MOD( TM,VAR_CW(1,NCT,NCW) )
!
!---    Coupled well is inactive set well pressure to be in 
!       1 with formation  ---
!
        IF( TMZ.LE.VAR_CW(1,1,NCW) .OR. 
     &    QM_CW(2,NCW).GE.TML_CW(NCW) ) THEN
          CALL EQUIL_COUP_WELL( NCW )
          ID_CW(8,NCW) = 1
          CYCLE CWLOOP
        ENDIF
        IF( NCT.EQ.1 ) THEN
          DO N = 2,6
            VAR_CWX(N) = VAR_CW(N,1,NCW)
          ENDDO
!
!---      Limit injection rate by total injected mass  ---
!
          VAR_CWX(2) = MIN( VAR_CWX(2),
     &      ((TML_CW(NCW)-QM_CW(2,NCW))*DTI) )
!
!---      Set well state  ---
!
          IT_CWX = IT_CW(NCW)
          IF( IT_CW(NCW).EQ.100 ) IT_CWX = ITS_CW(1,NCW)
        ELSE
          IFIND = 0
          DO M = 2,NCT
            IF( TMZ.LE.VAR_CW(1,M,NCW) ) THEN
              TD_CW = VAR_CW(1,M,NCW)-VAR_CW(1,M-1,NCW)
              DT_CW = MIN( VAR_CW(1,M,NCW)-TMZ,DT )
              TFX_CW = (TMZ-VAR_CW(1,M-1,NCW))/TD_CW
              DO N = 2,6
                VAR_CWX(N) = VAR_CW(N,M-1,NCW) + 
     &            TFX_CW*(VAR_CW(N,M,NCW)-VAR_CW(N,M-1,NCW))
              ENDDO
!
!---          Limit injection rate by total injected mass  ---
!
              VAR_CWX(2) = MIN( VAR_CWX(2),
     &          ((TML_CW(NCW)-QM_CW(2,NCW))*DTI) )
!
!---          Set well state  ---
!
              IT_CWX = IT_CW(NCW)
              IF( IT_CW(NCW).EQ.100 ) THEN
                NC = 0
                DO N = 1,IM_CW(NCW)
                  NC = NC + IMP_CW(N,NCW)
                  IF( NC.GE.M ) THEN
                    IT_CWX = ITS_CW(N,NCW)
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
              IFIND = 1
              EXIT
            ENDIF
          ENDDO
!
!---      Injection well is inactive set well pressure to be in 
!         1 with reservoir  ---
!
          IF( IFIND.EQ.0 ) THEN
            CALL EQUIL_COUP_WELL( NCW )
            ID_CW(8,NCW) = 1
            CYCLE CWLOOP
          ENDIF
        ENDIF
!
!---    Initialize local coupled-well mass residuals ---
!
        DO M = 1,(LUK_CW+1)
          RSL_CWX(M) = 0.D+0
        ENDDO
!
!---    Loop over pressure adjustments ---
!
        PALOOP: DO
!
!---    Load CO2 mass flux for use in RSDL_COUP_WELL ---
!
        FX_CW(NCW) = VAR_CWX(2)
!
!---    Load pressure limit for use in UPDT_COUP_WELL ---
!
        PL_CW(NCW) = VAR_CWX(3) - PATM
!
!---    Pressure controlled well ---
!
        IF( PL_CW(NCW)-P_CW(2,NCW).LT.EPSL ) THEN
          ID_CW(8,NCW) = 1
        ENDIF
!
!---    Excessive flow rate, pressure controlled well ---
!
        IF( VAR_CWX(2).GT.1.D+5 ) THEN
          ID_CW(8,NCW) = 1
          P_CW(2,NCW) = PL_CW(NCW)
        ENDIF
!
!---    Loop over increment indices ---
!
        DO M = 1,ISVC+2
          MW = MCW(M)
          MF = MFD(M)
!
!---      Injection well  ---
!
          IF( IT_CW(NCW).GT.0 ) THEN
            P_CWX = P_CW(MW,NCW)
!            PRINT *,'NCW = ',NCW,'VAR_CWX(2) = ',VAR_CWX(2),
!     &        'VAR_CWX(3) = ',VAR_CWX(3),'P_CWX = ',
!     &        P_CWX,'M = ',M,'ID = ',ID
            T_CWX = TF_CW(ID_CW(3,NCW))
!
!---        CO2 injection with zero solvated water  ---
!
            IF( IT_CWX.EQ.1 ) THEN
              PVW_CWX = 0.D+0
              PVA_CWX = MAX( P_CWX+PATM-PVW_CWX,0.D+0 )
              XMGA_CWX = 1.D+0
              XMGW_CWX = 0.D+0
              XGA_CWX = 1.D+0
              XGW_CWX = 0.D+0
              CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
              ISRX = 2
              CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
              RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---        CO2 injection with solvated water concentration 
!           declared as mass fraction  ---
!
            ELSEIF( IT_CWX.EQ.2 ) THEN
              XGW_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
              XGA_CWX = 1.D+0-XGW_CWX
              XMGA_CWX = (XGA_CWX/WTMA)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
              XMGW_CWX = (XGW_CWX/WTMW)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
              PVA_CWX = XMGA_CWX*(P_CWX+PATM)
              PVW_CWX = XMGW_CWX*(P_CWX+PATM)
              CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
              ISRX = 2
              CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
              RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---        CO2 injection with solvated water concentration 
!           declared as relative humidity  ---
!
            ELSEIF( IT_CWX.EQ.3 ) THEN
              CALL SP_W( T_CWX,PSW_CWX )
              PVW_CWX = PSW_CWX*VAR_CWX(4)
              PVA_CWX = MAX( P_CWX+PATM-PVW_CWX,0.D+0 )
              XMGA_CWX = PVA_CWX/(P_CWX+PATM)
              XMGW_CWX = PVW_CWX/(P_CWX+PATM)
              XGA_CWX = (XMGA_CWX*WTMA)/(XMGA_CWX*WTMA + XMGW_CWX*WTMW)
              XGW_CWX = (XMGW_CWX*WTMW)/(XMGA_CWX*WTMA + XMGW_CWX*WTMW)
              CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
              ISRX = 2
              CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
              RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---        Aqueous injection with zero dissolved CO2  ---
!
            ELSEIF( MOD(IT_CWX,10).EQ.4 ) THEN
              PL_CWX = MAX( P_CWX+PATM,0.D+0 )
              XLA_CWX = 0.D+0
!
!---          Aqueous injection with zero dissolved salt  ---
!
              IF( IT_CWX/10.EQ.1 ) THEN
                XLS_CWX = 0.D+0
!
!---          Aqueous injection with dissolved salt declared as 
!             brine mass fraction  ---
!
              ELSEIF( IT_CWX/10.EQ.2 ) THEN
                XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
!
!---          Aqueous injection with dissolved salt declared as 
!             brine relative saturation  ---
!
              ELSEIF( IT_CWX/10.EQ.3 ) THEN
                CALL SOL_LS( T_CWX,XLSMX )
                XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
              ENDIF
              CALL SP_B( T_CWX,XLS_CWX,PSBX )
              PVBX = PSBX
              CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
              PVWX = PVBX
              XLSSX = XLS_CWX
              CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &          XMGWX,XMLAX,XMLSX,XMLWX )
              XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
              XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
              CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
!
!---        Aqueous injection with dissolved CO2 declared as 
!           mass fraction  ---
!
            ELSEIF( MOD(IT_CWX,10).EQ.5 ) THEN
              PL_CWX = MAX( P_CWX+PATM,0.D+0 )
              XLA_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
!
!---          Aqueous injection with zero dissolved salt  ---
!
              IF( IT_CWX/10.EQ.1 ) THEN
                XLS_CWX = 0.D+0
!
!---          Aqueous injection with dissolved salt declared as 
!             brine mass fraction  ---
!
              ELSEIF( IT_CWX/10.EQ.2 ) THEN
                XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
!
!---          Aqueous injection with dissolved salt declared as 
!             brine relative saturation  ---
!
              ELSEIF( IT_CWX/10.EQ.3 ) THEN
                CALL SOL_LS( T_CWX,XLSMX )
                XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
              ENDIF
              CALL SP_B( T_CWX,XLS_CWX,PSBX )
              PVBX = PSBX
              CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
              PVWX = PVBX
              XLSSX = XLS_CWX
              CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &          XMGWX,XMLAX,XMLSX,XMLWX )
              XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
              XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
              CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
!
!---        Aqueous injection with dissolved CO2 declared as 
!           relative saturation  ---
!
            ELSEIF( MOD(IT_CWX,10).EQ.6 ) THEN
              PL_CWX = MAX( P_CWX+PATM,0.D+0 )
              XLA_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
!
!---          Aqueous injection with zero dissolved salt  ---
!
              IF( IT_CWX/10.EQ.1 ) THEN
                XLS_CWX = 0.D+0
!
!---          Aqueous injection with dissolved salt declared as 
!             brine mass fraction  ---
!
              ELSEIF( IT_CWX/10.EQ.2 ) THEN
                XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
 !
!---          Aqueous injection with dissolved salt declared as 
!             brine relative saturation  ---
!
              ELSEIF( IT_CWX/10.EQ.3 ) THEN
                CALL SOL_LS( T_CWX,XLSMX )
                XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
              ENDIF
              CALL SP_B( T_CWX,XLS_CWX,PSBX )
              PVBX = PSBX
              CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
              PVWX = PVBX
              XLSSX = XLS_CWX
              CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &          XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &          XMGWX,XMLAX,XMLSX,XMLWX )
              XLA_CWX = XLA_CWX*XLAX
              XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
              XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
              CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
            ENDIF
!
!---        Store top of coupled-well location in previous
!           coupled-well node location  ---
!
            XPX(1) = XTP_CW(1,ID_CW(1,NCW))
            YPX(1) = YTP_CW(1,ID_CW(1,NCW))
            ZPX(1) = ZTP_CW(1,ID_CW(1,NCW))
            ISX = ID_CW(3,NCW)
            IEX = ID_CW(4,NCW)
            IDX = 1
          ENDIF
!
!---      Loop over the nodes in the coupled well ---
!
          DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            N = IWN_CW(NWN)
            INVX = INV_CW(NWN)
            T_CWX = TF_CW(NWN)
!
!---        Coupled-well node centroids and projections ---
!
            XLX = PLX_CW(NWN)
            YLX = PLY_CW(NWN)
            ZLX = PLZ_CW(NWN)
            XPX(2) = 5.D-1*(XP_CW(2,NWN)+XP_CW(1,NWN))
            YPX(2) = 5.D-1*(YP_CW(2,NWN)+YP_CW(1,NWN))
            ZPX(2) = 5.D-1*(ZP_CW(2,NWN)+ZP_CW(1,NWN))
!
!---        Well pressure using previous coupled-well node density ---
!
            IF( IT_CWX.LE.3 ) THEN
              P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOG_CWX
            ELSE
              P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOL_CWX
            ENDIF
!
!---        Update coupled-well node density and viscosity
!           for CO2 injection, by using a fixed solvated
!           water mass fraction  ---
!
            IF( IT_CWX.LE.3 ) THEN
              XMGA_CWX = (XGA_CWX/WTMA)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
              XMGW_CWX = (XGW_CWX/WTMW)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
              PVA_CWX = XMGA_CWX*(P_CWX+PATM)
              PVW_CWX = XMGW_CWX*(P_CWX+PATM)
              CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
              ISRX = 2
              CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
              RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
              CALL VISC_A( T_CWX,RHOGA_CWX,VISGA_CWX )
              CALL VISC_W( T_CWX,PVW_CWX,RHOGW_CWX,VISGW_CWX )
              CALL VISC_G( VISGA_CWX,VISGW_CWX,XMGA_CWX,XMGW_CWX,
     &          VISG_CWX )
!
!---        Update coupled-well node density and viscosity
!           for aqueous injection, by using a fixed CO2 and salt
!           dissolved mass fraction  ---
!
            ELSE
              PL_CWX = MAX( P_CWX+PATM,0.D+0 )
              CALL SP_B( T_CWX,XLS_CWX,PSBX )
              PVAX = MAX( PL_CWX-PSBX,0.D+0 )
              CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
              WTMLX = 1.D+0/(XLA_CWX/WTMA+XLS_CWX/WTMS+XLW_CWX/WTMW)
              XMLA_CWX = XLA_CWX*WTMLX/WTMA
              CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
              ISRX = 1
              CALL DENS_W( T_CWX,PL_CWX,RHOLWX,RHOX,ISRX )
              CALL VISC_W( T_CWX,PL_CWX,RHOLWX,VISLWX )
              CALL VISC_B( T_CWX,XLS_CWX,VISLWX,VISBX )
              CALL DENS_A( T_CWX,PVAX,RHOGAX,I_VX )
              CALL VISC_A( T_CWX,RHOGAX,VISGAX )
              CALL VISC_L( XMLA_CWX,VISBX,VISGAX,VISL_CWX )
            ENDIF
!
!---        Skip for off-processor well nodes and ghost cells ---
!
            IGHOST = 0
            IF( N.NE.0 ) IGHOST = IGHC(N)
            IF( N.NE.0 .AND. IGHOST.EQ.0 ) THEN
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
                XPNX = 0.D+0
                YPNX = 0.D+0
                ZPNX = ZP(N)
!
!---          Cylindrical coordinates  ---
!
              ELSEIF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
                XPNX = XP(N)*COS(YP(N))
                YPNX = XP(N)*SIN(YP(N))
                ZPNX = ZP(N)
!
!---          Cartesian or boundary-fitted orthogonal coordinates  ---
!
              ELSE
                XPNX = XP(N)
                YPNX = YP(N)
                ZPNX = ZP(N)
              ENDIF
!
!---          Well pressure at the node centroid, used for coupled-well
!             nodal output  ---
!
              IF( M.EQ.1 ) THEN
                NWF = IWP_CW(NWN)
                IF( IT_CWX.LE.3 ) THEN
                  PF_CW(NWF) = P_CWX - (ZPNX-ZPX(1))*GRAV*RHOG_CWX
                ELSE
                  PF_CW(NWF) = P_CWX - (ZPNX-ZPX(1))*GRAV*RHOL_CWX
                ENDIF
              ENDIF
!
!---          Adjust the formation pressure to the coupled-well node
!             centroid  ---
!
              IF( (SG(MF,N)-SGT(MF,N)).GT.EPSL ) THEN
                PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(MF,N)
              ELSE
                PGFX = PG(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
              ENDIF
              PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
!
!---          Equivalent field node radius components  ---
!
              PERMX = MAX( PERM(1,N),1.D-20 )
              PERMY = MAX( PERM(2,N),1.D-20 )
              PERMZ = MAX( PERM(3,N),1.D-20 )
              RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
                ROZ = RP(N)
                RWX = MIN( RWX,9.999D-1*ROZ )
                PERMX = PERMRF(MF,N)*PERM(1,N)
                WI_CWX = 2.D+0*GPI*PERMX*ZLX/(LOG(ROZ/RWX) + 
     &            PAR_CW(1,INVX))
              ELSE
                PERMYZ = SQRT(PERMY/PERMZ)
                PERMZY = SQRT(PERMZ/PERMY)
                DXGFX = DXGF(N)/FF_CW(1,NCW)
                DYGFX = DYGF(N)*RP(N)/FF_CW(2,NCW)
                DZGFX = DZGF(N)/FF_CW(3,NCW)
                ROX = 2.8D-1*SQRT(PERMYZ*(DZGFX**2) + PERMZY*(DYGFX**2))
     &          /(SQRT(PERMYZ)+SQRT(PERMZY))
                PERMZX = SQRT(PERMZ/PERMX)
                PERMXZ = SQRT(PERMX/PERMZ)
                ROY = 2.8D-1*SQRT(PERMZX*(DXGFX**2) + PERMXZ*(DZGFX**2))
     &            /(SQRT(PERMZX)+SQRT(PERMXZ))
                PERMYX = SQRT(PERMY/PERMX)
                PERMXY = SQRT(PERMX/PERMY)
                ROZ = 2.8D-1*SQRT(PERMYX*(DXGFX**2) + PERMXY*(DYGFX**2))
     &            /(SQRT(PERMYX)+SQRT(PERMXY))
!
!---            Well index components  ---
!
                PERMX = PERMRF(MF,N)*PERM(1,N)
                PERMY = PERMRF(MF,N)*PERM(2,N)
                PERMZ = PERMRF(MF,N)*PERM(3,N)
                WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &            (LOG(ROX/RWX)+PAR_CW(1,INVX))
                WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &            (LOG(ROY/RWX)+PAR_CW(1,INVX))
                WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &            (LOG(ROZ/RWX)+PAR_CW(1,INVX))
                WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
              ENDIF
!
!---          Mass fluxes, positive into the node  ---
!
              DPGX = MAX( P_CWX-PGFX,0.D+0 )
              DPLX = MAX( P_CWX-MAX(PLFX,PGFX),0.D+0 )
!
!---          CO2 injection, flux from well to formation  ---
!
              IF( IT_CWX.LE.3 ) THEN
                FX_CWX = WI_CWX*RHOG_CWX*DPGX/VISG_CWX
                FXA_CW(M,NWN) = FX_CWX*XGA_CWX
                FXW_CW(M,NWN) = FX_CWX*XGW_CWX
!
!---            Volumetric injection well fluxes  ---
!
                Q_CW(1,NWN) = FX_CWX/RHOG_CWX
                Q_CW(3,NWN) = FX_CWX/RHOG_CWX
!
!---          Aqueous injection, flux from well to formation  ---
!
              ELSE
                FX_CWX = WI_CWX*RHOL_CWX*DPLX/VISL_CWX
                FXA_CW(M,NWN) = FX_CWX*XLA_CWX
                FXS_CW(M,NWN) = FX_CWX*XLS_CWX
                FXW_CW(M,NWN) = FX_CWX*XLW_CWX
!
!---            Volumetric injection well fluxes  ---
!
                Q_CW(1,NWN) = FX_CWX/RHOG_CWX
                Q_CW(2,NWN) = FX_CWX/RHOG_CWX
              ENDIF
!              PRINT *,'FXA_CW(',M,',',NWN,') = ',FXA_CW(M,NWN),
!     &          'NCW = ',NCW,'ND(',N,') = ',ND(N),'WI_CWX = ',WI_CWX,
!     &          'RHOG_CWX = ',RHOG_CWX,'DPGX = ',DPGX,
!     &          'VISG_CWX = ',VISG_CWX,'P_CWX = ',P_CWX,'ID = ',ID
!              PRINT *,'FXS_CW(',M,',',NWN,') = ',FXS_CW(M,NWN),
!     &          'NCW = ',NCW,'ND(',N,') = ',ND(N),'WI_CWX = ',WI_CWX,
!     &          'RHOG_CWX = ',RHOG_CWX,'DPGX = ',DPGX,
!     &          'VISG_CWX = ',VISG_CWX,'ID = ',ID
!              PRINT *,'FXW_CW(',M,',',NWN,') = ',FXW_CW(M,NWN),
!     &          'NCW = ',NCW,'ND(',N,') = ',ND(N),'WI_CWX = ',WI_CWX,
!     &          'RHOG_CWX = ',RHOG_CWX,'DPGX = ',DPGX,
!     &          'VISG_CWX = ',VISG_CWX,'ID = ',ID
            ENDIF
!
!---        Store current coupled-well node location in previous
!           coupled-well node location  ---
!
            XPX(1) = XPX(2)
            YPX(1) = YPX(2)
            ZPX(1) = ZPX(2)
!            IF( M.EQ.1 .AND. FXA_CW(1,NWN).GT.EPSL ) 
!     &        PRINT *,'FXA_CW(1,',NWN,') = ',FXA_CW(1,NWN),
!     &        ' P_CWX = ',P_CWX,' PLFX = ',PLFX,
!     &        ' PGFX = ',PGFX,' WI_CWX = ',WI_CWX,
!     &        ' RHOG_CWX = ',RHOG_CWX,' DPGX = ',DPGX,
!     &        ' VISG_CWX = ',VISG_CWX,' ID = ',ID
          ENDDO
        ENDDO
!
!---    CO2 mass balance residuals for injection type coupled well  ---
!
        NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
        NX = (NWFX*ISVC)+2
        RSL_CWX(1) = 0.D+0
        RSL_CWX(NX) = 0.D+0
        DO M = 1,4
          QML_CWX(M) = 0.D+0
          QMG_CWX(M) = 0.D+0
        ENDDO
!
!---    Loop over coupled-well nodes  ---
!
        DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          RSL_CWX(1) = RSL_CWX(1) - FXA_CW(1,NWN)
     &      - FXW_CW(1,NWN) - FXS_CW(1,NWN)
          RSL_CWX(NX) = RSL_CWX(NX) - FXA_CW(ISVC+2,NWN) 
     &      - FXW_CW(ISVC+2,NWN) - FXS_CW(ISVC+2,NWN)
          QML_CWX(1) = QML_CWX(1) + FXA_CW(1,NWN)
          QML_CWX(2) = QML_CWX(2) + FXW_CW(1,NWN)
          QML_CWX(3) = QML_CWX(3) + FXS_CW(1,NWN)
          QML_CWX(4) = QML_CWX(4) + FXA_CW(ISVC+2,NWN)
     &      + FXW_CW(ISVC+2,NWN) + FXS_CW(ISVC+2,NWN)
        ENDDO
!
!---    Global injection well mass balance residuals  ---
!
!        PRINT *,'QML_CWX(4) = ',QML_CWX(4),' ID = ',ID
        CALL MPI_ALLREDUCE( QML_CWX,QMG_CWX,4,MPI_REAL8,MPI_SUM,
     &    MPI_COMM_WORLD,IERR )
        QM_CW(1,NCW) = QMG_CWX(1)
        QM_CW(3,NCW) = QMG_CWX(2)
        QM_CW(7,NCW) = QMG_CWX(3)
        QM_CWX = QMG_CWX(4)
!        PRINT *,'QMG_CWX(4) = ',QMG_CWX(4),' ID = ',ID
!
!---    Insufficient increment in coupled-well pressure to create
!       flow from well, increase increment and well pressure  ---
!
        IF( ABS(QM_CWX/DNR_CW(NCW)).LT.1.D-7 ) THEN
          DNR_CW(NCW) = 1.25D+0*DNR_CW(NCW)
          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
          IF( P_CW(3,NCW).GE.PL_CW(NCW) ) THEN
            P_CW(2,NCW) = PL_CW(NCW)
            DNR_CW(NCW) = 1.D-1
            P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
            ID_CW(8,NCW) = 1
            CYCLE CWLOOP
          ENDIF
!
!---      Zero coupled-well fluxes ---
!
          DO M = 1,4
            QML_CWX(M) = 0.D+0
          ENDDO
!
!---      Flow controlled well ---
!
          ID_CW(8,NCW) = 0
!
!---      Loop over coupled-well nodes  ---
!
          DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---        Loop over increment indices  ---
!
            DO M = 1,ISVC+2
              FXA_CW(M,NWN) = 0.D+0
              FXS_CW(M,NWN) = 0.D+0
              FXW_CW(M,NWN) = 0.D+0
            ENDDO
          ENDDO
          CYCLE PALOOP
!
!---    Excessive well flow for pressure controlled well,
!       reduce well pressure  ---
!
        ELSEIF( ABS(ID_CW(8,NCW)).EQ.1 .AND. 
     &    (QM_CW(1,NCW)+QM_CW(3,NCW)+QM_CW(7,NCW)).GT.VAR_CWX(2) ) THEN
          DQ_CWX = 1.D+1*DQ_CWX
          DP_CWX = -DQ_CWX*PATM
          P_CW(2,NCW) = P_CW(2,NCW) + DP_CWX
          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
!
!---      Zero coupled-well fluxes ---
!
          DO M = 1,4
            QML_CWX(M) = 0.D+0
          ENDDO
!
!---      Pressure correction for pressure controlled well ---
!
          ID_CW(8,NCW) = -1
!
!---      Loop over coupled-well nodes  ---
!
          DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---        Loop over increment indices  ---
!
            DO M = 1,ISVC+2
              FXA_CW(M,NWN) = 0.D+0
              FXS_CW(M,NWN) = 0.D+0
              FXW_CW(M,NWN) = 0.D+0
            ENDDO
          ENDDO
          CYCLE PALOOP
!
!---    Acceptable pressure reduction for pressure controlled well, 
!       switch to flow controlled well  ---
!
        ELSEIF( ID_CW(8,NCW).EQ.-1 .AND. 
     &    (QM_CW(1,NCW)+QM_CW(3,NCW)+QM_CW(7,NCW)).LT.VAR_CWX(2) ) THEN
          ID_CW(8,NCW) = 0
          EXIT PALOOP
        ELSE
          EXIT PALOOP
        ENDIF
        ENDDO PALOOP
!        PRINT *,'QM_CW(:,NCW) = ',(QM_CW(M,NCW),M=1,8),' ID = ',ID
!        PRINT *,'P_CW(:,1) = ',(P_CW(M,1),M=1,3),' ID = ',ID
!
!---    Loop over field nodes that contain coupled-well nodes  ---
!
        DO NWF = ID_CW(5,NCW),ID_CW(6,NCW)
!
!---      Skip for processors without field nodes with coupled-well
!         nodes  ---
!
          IF( IWF_CW(NWF).EQ.0 ) CYCLE
          M1 = (NWF-ID_CW(5,NCW))*ISVC + 1
          DO  M2 = 1,ISVC
            M3 = M1+M2
            RSL_CWX(M3) = 0.D+0
!
!---        Loop over coupled-well nodes  ---
!
            DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---          Skip for processors without field nodes with coupled-well
!             nodes  ---
!
              IF( IWN_CW(NWN).EQ.0 ) CYCLE
!
!---          If coupled-well node is within the current field
!             node, use incremented fluxes  ---
!
              IF( IWF_CW(NWF).EQ.IWN_CW(NWN) ) THEN
                RSL_CWX(M3) = RSL_CWX(M3) - FXA_CW(M2+1,NWN)
     &            - FXW_CW(M2+1,NWN) - FXS_CW(M2+1,NWN)
!
!---          If coupled-well node is outside the current field
!             node, use un-incremented fluxes  ---
!
              ELSE
                RSL_CWX(M3) = RSL_CWX(M3) - FXA_CW(1,NWN)
     &            - FXW_CW(1,NWN) - FXS_CW(1,NWN)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!
!---    Redefine local residuals to be the change in the coupled-well
!       residual mass with respect to the change in the field-node
!       primary variables, except for the first and last residuals  ---
!
        DO NWF = ID_CW(5,NCW),ID_CW(6,NCW)
          N = IWF_CW(NWF)
          IF( N.EQ.0 ) CYCLE
          MX = (NWF-ID_CW(5,NCW))*ISVC + 1
          DO M = 1,ISVC
            RSL_CWX(MX+M) = (RSL_CWX(MX+M)-RSL_CWX(1))/DNR(M,N)
          ENDDO
        ENDDO
!
!---    Global injection well mass balance residuals  ---
!
        CALL MPI_ALLREDUCE( RSL_CWX,RS_CW(1,NCW),NX,MPI_REAL8,MPI_SUM,
     &    MPI_COMM_WORLD,IERR )
!
!---    Add mass injection rate into global injection well mass 
!       balance residuals for the first and last residuals  ---
!
        RS_CW(1,NCW) = RS_CW(1,NCW) + VAR_CWX(2)
        RS_CW(NX,NCW) = RS_CW(NX,NCW) + VAR_CWX(2)
!        DO M = 1,NX
!          PRINT *,'RS_CW(',M,',NCW) = ',RS_CW(M,NCW),' ID = ',ID
!        ENDDO
      ENDDO CWLOOP
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLUX_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_COUP_WELL
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
!     Define well nodes, determine trajectory points, and 
!     check for well trajectories within node surface planes
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_COUP_WELL'
!
!---  Loop over coupled wells ---
!
      DO NCW = 1,N_CW
!
!---    Coupled-well is an injection type well, estimate the
!       well pressure and whether the well is pressure or flow
!       controlled  ---
!
        IF( IT_CW(NCW).GT.0 .AND. ISLC(70).EQ.1 ) THEN
          CALL INJP_COUP_WELL( NCW )
        ENDIF
        DNR_CW(NCW) = 1.D-1
        P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
!        PRINT *,'P_CW(:,',NCW,') = ',(P_CW(M,NCW),M=1,3),' ID = ',ID
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INJP_COUP_WELL( NCW )
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
!
!     STOMPX-CO2
!
!     Injection coupled well model
!     
!     Rate controlled or pressure controlled
!
!     Flux of energy, water mass, and CO2 mass
!     from coupled-well nodes to field nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE HYST
      USE GRID
      USE FDVP
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XPX(2),YPX(2),ZPX(2)
      REAL*8 VAR_CWX(6)
      REAL*8 FXA_CWX(3,LWN_CW),FXW_CWX(3,LWN_CW),FXS_CWX(3,LWN_CW)
      REAL*8 QM_CWX(3),QMG_CWX(3),P_CWY(3),VAR_CW2X(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INJP_COUP_WELL'
!
!---  Loop over coupled-well nodes  ---
!
      DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!
!---    Loop over increment indices  ---
!
        DO M = 1,2
          FXA_CWX(M,NWN) = 0.D+0
          FXS_CW(M,NWN) = 0.D+0
          FXW_CWX(M,NWN) = 0.D+0
        ENDDO
      ENDDO
!
!---  Injection well time interval ---
!
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
      NCT = 0
      DO M = 1,IM_CW(NCW)
        NCT = NCT + IMP_CW(M,NCW)
      ENDDO
      IF( ICC_CW(NCW).EQ.1 ) TMZ = MOD( TM,VAR_CW(1,NCT,NCW) )
!
!---  Coupled well is inactive set well pressure to be in 
!     1 with formation  ---
!
      IF( TMZ.LE.VAR_CW(1,1,NCW) .OR. 
     &  QM_CW(2,NCW).GE.TML_CW(NCW) ) THEN
        CALL EQUIL_COUP_WELL( NCW )
        ID_CW(8,NCW) = 1
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( NCT.EQ.1 ) THEN
        DO N = 2,6
          VAR_CWX(N) = VAR_CW(N,1,NCW)
        ENDDO
!
!---    Limit injection rate by total injected mass  ---
!
        VAR_CWX(2) = MIN( VAR_CWX(2),
     &    ((TML_CW(NCW)-QM_CW(2,NCW))*DTI) )
!
!---    Set well state  ---
!
        IT_CWX = IT_CW(NCW)
        IF( IT_CW(NCW).EQ.100 ) IT_CWX = ITS_CW(1,NCW)
      ELSE
        IFIND = 0
        DO M = 2,NCT
          IF( TMZ.LE.VAR_CW(1,M,NCW) ) THEN
            TD_CW = VAR_CW(1,M,NCW)-VAR_CW(1,M-1,NCW)
            DT_CW = MIN( VAR_CW(1,M,NCW)-TMZ,DT )
            TFX_CW = (TMZ-VAR_CW(1,M-1,NCW))/TD_CW
            DO N = 2,6
              VAR_CWX(N) = VAR_CW(N,M-1,NCW) + 
     &          TFX_CW*(VAR_CW(N,M,NCW)-VAR_CW(N,M-1,NCW))
            ENDDO
!
!---        Limit injection rate by total injected mass  ---
!
            VAR_CWX(2) = MIN( VAR_CWX(2),
     &        ((TML_CW(NCW)-QM_CW(2,NCW))*DTI) )
!
!---        Set well state  ---
!
            IT_CWX = IT_CW(NCW)
            IF( IT_CW(NCW).EQ.100 ) THEN
              NC = 0
              DO N = 1,IM_CW(NCW)
                NC = NC + IMP_CW(N,NCW)
                IF( NC.GE.M ) THEN
                  IT_CWX = ITS_CW(N,NCW)
                  EXIT
                ENDIF
              ENDDO
            ENDIF
            IFIND = 1
            EXIT
          ENDIF
        ENDDO
!
!---    Injection well is inactive set well pressure to be in 
!       1 with reservoir  ---
!
        IF( IFIND.EQ.0 ) THEN
          CALL EQUIL_COUP_WELL( NCW )
          ID_CW(8,NCW) = 1
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
      ENDIF
!
!---  Load pressure limit ---
!
      PL_CW(NCW) = VAR_CWX(3) - PATM
!
!---  Upper pressure limit ---
!
      PL_CWX = VAR_CWX(3) - PATM
      P_CWY(1) = P_CW(2,NCW)
      P_CWY(2) = PL_CWX
      DP_CWX = 1.D-1
      ICHK_CWX = 0
      ML = 1
      MU = 2
      DO NC = 1,32
!
!---    Loop over increment indices ---
!
        DO M = ML,MU
          N = IWN_CW(ID_CW(3,NCW))
          P_CWX = P_CWY(M)
          T_CWX = TF_CW(ID_CW(3,NCW))
!
!---      CO2 injection with zero solvated water  ---
!
          IF( IT_CWX.EQ.1 ) THEN
            PVW_CWX = 0.D+0
            PVA_CWX = MAX( P_CWX+PATM-PVW_CWX,0.D+0 )
            XMGA_CWX = 1.D+0
            XMGW_CWX = 0.D+0
            XGA_CWX = 1.D+0
            XGW_CWX = 0.D+0
            CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
            ISRX = 2
            CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
            RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---      CO2 injection with solvated water concentration 
!         declared as mass fraction  ---
!
          ELSEIF( IT_CWX.EQ.2 ) THEN
            XGW_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
            XGA_CWX = 1.D+0-XGW_CWX
            XMGA_CWX = (XGA_CWX/WTMA)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
            XMGW_CWX = (XGW_CWX/WTMW)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
            PVA_CWX = XMGA_CWX*(P_CWX+PATM)
            PVW_CWX = XMGW_CWX*(P_CWX+PATM)
            CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
            ISRX = 2
            CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
            RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---      CO2 injection with solvated water concentration 
!         declared as relative humidity  ---
!
          ELSEIF( IT_CWX.EQ.3 ) THEN
            CALL SP_W( T_CWX,PSW_CWX )
            PVW_CWX = PSW_CWX*VAR_CWX(4)
            PVA_CWX = MAX( P_CWX+PATM-PVW_CWX,0.D+0 )
            XMGA_CWX = PVA_CWX/(P_CWX+PATM)
            XMGW_CWX = PVW_CWX/(P_CWX+PATM)
            XGA_CWX = (XMGA_CWX*WTMA)/(XMGA_CWX*WTMA + XMGW_CWX*WTMW)
            XGW_CWX = (XMGW_CWX*WTMW)/(XMGA_CWX*WTMA + XMGW_CWX*WTMW)
            CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
            ISRX = 2
            CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
            RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
!
!---      Aqueous injection with zero dissolved CO2  ---
!
          ELSEIF( MOD(IT_CWX,10).EQ.4 ) THEN
            PL_CWX = MAX( P_CWX+PATM,0.D+0 )
            XLA_CWX = 0.D+0
!
!---        Aqueous injection with zero dissolved salt  ---
!
            IF( IT_CWX/10.EQ.1 ) THEN
              XLS_CWX = 0.D+0
!
!---        Aqueous injection with dissolved salt declared as 
!           brine mass fraction  ---
!
            ELSEIF( IT_CWX/10.EQ.2 ) THEN
              XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
!
!---        Aqueous injection with dissolved salt declared as 
!           brine relative saturation  ---
!
            ELSEIF( IT_CWX/10.EQ.3 ) THEN
              CALL SOL_LS( T_CWX,XLSMX )
              XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
            ENDIF
            CALL SP_B( T_CWX,XLS_CWX,PSBX )
            PVBX = PSBX
            CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
            PVWX = PVBX
            XLSSX = XLS_CWX
            CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &        XMGWX,XMLAX,XMLSX,XMLWX )
            XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
            XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
            CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
!
!---      Aqueous injection with dissolved CO2 declared as 
!         mass fraction  ---
!
          ELSEIF( MOD(IT_CWX,10).EQ.5 ) THEN
            PL_CWX = MAX( P_CWX+PATM,0.D+0 )
            XLA_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
!
!---        Aqueous injection with zero dissolved salt  ---
!
            IF( IT_CWX/10.EQ.1 ) THEN
              XLS_CWX = 0.D+0
!
!---        Aqueous injection with dissolved salt declared as 
!           brine mass fraction  ---
!
            ELSEIF( IT_CWX/10.EQ.2 ) THEN
              XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
!
!---        Aqueous injection with dissolved salt declared as 
!           brine relative saturation  ---
!
            ELSEIF( IT_CWX/10.EQ.3 ) THEN
              CALL SOL_LS( T_CWX,XLSMX )
              XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
            ENDIF
            CALL SP_B( T_CWX,XLS_CWX,PSBX )
            PVBX = PSBX
            CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
            PVWX = PVBX
            XLSSX = XLS_CWX
            CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &        XMGWX,XMLAX,XMLSX,XMLWX )
            XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
            XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
            CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
!
!---      Aqueous injection with dissolved CO2 declared as 
!         relative saturation  ---
!
          ELSEIF( MOD(IT_CWX,10).EQ.6 ) THEN
            PL_CWX = MAX( P_CWX+PATM,0.D+0 )
            XLA_CWX = MIN( MAX( VAR_CWX(4),0.D+0 ),1.D+0 )
!
!---        Aqueous injection with zero dissolved salt  ---
!
            IF( IT_CWX/10.EQ.1 ) THEN
              XLS_CWX = 0.D+0
!
!---        Aqueous injection with dissolved salt declared as 
!           brine mass fraction  ---
!
            ELSEIF( IT_CWX/10.EQ.2 ) THEN
              XLS_CWX = MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
 !
!---        Aqueous injection with dissolved salt declared as 
!           brine relative saturation  ---
!
            ELSEIF( IT_CWX/10.EQ.3 ) THEN
              CALL SOL_LS( T_CWX,XLSMX )
              XLS_CWX = XLSMX*MIN( MAX( VAR_CWX(6),0.D+0 ),1.D+0 )
            ENDIF
            CALL SP_B( T_CWX,XLS_CWX,PSBX )
            PVBX = PSBX
            CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
            PVWX = PVBX
            XLSSX = XLS_CWX
            CALL EQUIL( T_CWX,PL_CWX,PGAX,PGWX,PSBX,PVWX,
     &        XGAX,XGWX,XLAX,XLSSX,XLWX,XMGAX,
     &        XMGWX,XMLAX,XMLSX,XMLWX )
            XLA_CWX = XLA_CWX*XLAX
            XLS_CWX = XLS_CWX + (XLSSX-XLS_CWX)*(XLA_CWX/XLAX)
            XLW_CWX = 1.D+0-XLA_CWX-XLS_CWX
            CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
          ENDIF
!
!---      Mass flow rate, kg/s  ---
!
          VAR_CW2X(M) = VAR_CWX(2)
!
!---      Store top of coupled-well location in previous
!         coupled-well node location  ---
!
          XPX(1) = XTP_CW(1,ID_CW(1,NCW))
          YPX(1) = YTP_CW(1,ID_CW(1,NCW))
          ZPX(1) = ZTP_CW(1,ID_CW(1,NCW))
!
!---      Loop over the nodes in the coupled well ---
!
          DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            N = IWN_CW(NWN)
            INVX = INV_CW(NWN)
            T_CWX = TF_CW(NWN)
!
!---        Coupled-well node centroids and projections ---
!
            XLX = PLX_CW(NWN)
            YLX = PLY_CW(NWN)
            ZLX = PLZ_CW(NWN)
            XPX(2) = 5.D-1*(XP_CW(2,NWN)+XP_CW(1,NWN))
            YPX(2) = 5.D-1*(YP_CW(2,NWN)+YP_CW(1,NWN))
            ZPX(2) = 5.D-1*(ZP_CW(2,NWN)+ZP_CW(1,NWN))
!
!---        Well pressure using previous coupled-well node density ---
!
            IF( IT_CWX.LE.3 ) THEN
              P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOG_CWX
            ELSE
              P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOL_CWX
            ENDIF
!
!---        Update coupled-well node density and viscosity
!           for CO2 injection, by using a fixed solvated
!           water mass fraction  ---
!
            IF( IT_CWX.LE.3 ) THEN
              XMGA_CWX = (XGA_CWX/WTMA)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
              XMGW_CWX = (XGW_CWX/WTMW)/((XGA_CWX/WTMA)+(XGW_CWX/WTMW))
              PVA_CWX = XMGA_CWX*(P_CWX+PATM)
              PVW_CWX = XMGW_CWX*(P_CWX+PATM)
              CALL DENS_A( T_CWX,PVA_CWX,RHOGA_CWX,I_VX )
              ISRX = 2
              CALL DENS_W( T_CWX,PVW_CWX,RHOLW_CWX,RHOGW_CWX,ISRX )
              RHOG_CWX = XGA_CWX*RHOGA_CWX + XGW_CWX*RHOGW_CWX
              CALL VISC_A( T_CWX,RHOGA_CWX,VISGA_CWX )
              CALL VISC_W( T_CWX,PVW_CWX,RHOGW_CWX,VISGW_CWX )
              CALL VISC_G( VISGA_CWX,VISGW_CWX,XMGA_CWX,XMGW_CWX,
     &          VISG_CWX )
!
!---        Update coupled-well node density and viscosity
!           for aqueous injection, by using a fixed CO2 and salt
!           dissolved mass fraction  ---
!
            ELSE
              PL_CWX = MAX( P_CWX+PATM,0.D+0 )
              CALL SP_B( T_CWX,XLS_CWX,PSBX )
              PVAX = MAX( PL_CWX-PSBX,0.D+0 )
              CALL DENS_B( T_CWX,PL_CWX,XLS_CWX,RHOBX )
              WTMLX = 1.D+0/(XLA_CWX/WTMA+XLS_CWX/WTMS+XLW_CWX/WTMW)
              XMLA_CWX = XLA_CWX*WTMLX/WTMA
              CALL DENS_L( T_CWX,RHOBX,XLA_CWX,RHOL_CWX )
              ISRX = 1
              CALL DENS_W( T_CWX,PL_CWX,RHOLWX,RHOX,ISRX )
              CALL VISC_W( T_CWX,PL_CWX,RHOLWX,VISLWX )
              CALL VISC_B( T_CWX,XLS_CWX,VISLWX,VISBX )
              CALL DENS_A( T_CWX,PVAX,RHOGAX,I_VX )
              CALL VISC_A( T_CWX,RHOGAX,VISGAX )
              CALL VISC_L( XMLA_CWX,VISBX,VISGAX,VISL_CWX )
            ENDIF
!
!---        Skip for off-processor well nodes and ghost cells ---
!
            IGHOST = 0
            IF( N.NE.0 ) IGHOST = IGHC(N)
            IF( N.NE.0 .AND. IGHOST.NE.0 ) THEN
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
                XPNX = 0.D+0
                YPNX = 0.D+0
                ZPNX = ZP(N)
!
!---          Cylindrical coordinates  ---
!
              ELSEIF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
                XPNX = XP(N)*COS(YP(N))
                YPNX = XP(N)*SIN(YP(N))
                ZPNX = ZP(N)
!
!---          Cartesian or boundary-fitted orthogonal coordinates  ---
!
              ELSE
                XPNX = XP(N)
                YPNX = YP(N)
                ZPNX = ZP(N)
              ENDIF
!
!---          Adjust the formation pressure to the coupled-well node
!             centroid  ---
!
              IF( (SG(2,N)-SGT(2,N)).GT.EPSL ) THEN
                PGFX = PG(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOG(2,N)
              ELSE
                PGFX = PG(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
              ENDIF
              PLFX = PL(2,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(2,N)
!
!---          Equivalent field node radius components  ---
!
              PERMX = MAX( PERM(1,N),1.D-20 )
              PERMY = MAX( PERM(2,N),1.D-20 )
              PERMZ = MAX( PERM(3,N),1.D-20 )
              RWX = MAX( PAR_CW(2,INVX),1.D-20 )
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
                ROZ = RP(N)
                RWX = MIN( RWX,9.999D-1*ROZ )
                PERMX = PERMRF(2,N)*PERM(1,N)
                WI_CWX = 2.D+0*GPI*PERMX*ZLX/(LOG(ROZ/RWX) + 
     &            PAR_CW(1,INVX))
              ELSE
                PERMYZ = SQRT(PERMY/PERMZ)
                PERMZY = SQRT(PERMZ/PERMY)
                DXGFX = DXGF(N)/FF_CW(1,NCW)
                DYGFX = DYGF(N)*RP(N)/FF_CW(2,NCW)
                DZGFX = DZGF(N)/FF_CW(3,NCW)
                ROX = 2.8D-1*SQRT(PERMYZ*(DZGFX**2) + PERMZY*(DYGFX**2))
     &          /(SQRT(PERMYZ)+SQRT(PERMZY))
                PERMZX = SQRT(PERMZ/PERMX)
                PERMXZ = SQRT(PERMX/PERMZ)
                ROY = 2.8D-1*SQRT(PERMZX*(DXGFX**2) + PERMXZ*(DZGFX**2))
     &            /(SQRT(PERMZX)+SQRT(PERMXZ))
                PERMYX = SQRT(PERMY/PERMX)
                PERMXY = SQRT(PERMX/PERMY)
                ROZ = 2.8D-1*SQRT(PERMYX*(DXGFX**2) + PERMXY*(DYGFX**2))
     &            /(SQRT(PERMYX)+SQRT(PERMXY))
!
!---            Well index components  ---
!
                PERMX = PERMRF(2,N)*PERM(1,N)
                PERMY = PERMRF(2,N)*PERM(2,N)
                PERMZ = PERMRF(2,N)*PERM(3,N)
                WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &            (LOG(ROX/RWX)+PAR_CW(1,INVX))
                WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &            (LOG(ROY/RWX)+PAR_CW(1,INVX))
                WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &            (LOG(ROZ/RWX)+PAR_CW(1,INVX))
                WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
              ENDIF
!
!---          Mass fluxes, positive into the node  ---
!
              DPGX = MAX( P_CWX-PGFX,0.D+0 )
              DPLX = MAX( P_CWX-MAX(PLFX,PGFX),0.D+0 )
!
!---          Zero fluxes from well to reservoir  ---
!
              FXW_CWX(M,NWN) = 0.D+0
              FXS_CWX(M,NWN) = 0.D+0
              FXA_CWX(M,NWN) = 0.D+0
!
!---          CO2 injection, flux from well to formation  ---
!
              IF( IT_CWX.LE.3 ) THEN
                FX_CWX = WI_CWX*RHOG_CWX*DPGX/VISG_CWX
                FXA_CW(M,NWN) = FX_CWX*XGA_CWX
                FXW_CW(M,NWN) = FX_CWX*XGW_CWX
!
!---          Aqueous injection, flux from well to formation  ---
!
              ELSE
                FX_CWX = WI_CWX*RHOL_CWX*DPLX/VISL_CWX
                FXA_CW(M,NWN) = FX_CWX*XLA_CWX
                FXS_CW(M,NWN) = FX_CWX*XLS_CWX
                FXW_CW(M,NWN) = FX_CWX*XLW_CWX
              ENDIF
!
!---          Store current coupled-well node location in previous
!             coupled-well node location  ---
!
              XPX(1) = XPX(2)
              YPX(1) = YPX(2)
              ZPX(1) = ZPX(2)
            ENDIF
          ENDDO
!
!---      Mass balance residuals for injection type coupled well  ---
!
          QM_CWX(M) = 0.D+0
!
!---      Loop over coupled-well nodes  ---
!
          DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
            QM_CWX(M) = QM_CWX(M) + FXW_CWX(M,NWN)
     &       + FXA_CWX(M,NWN) + FXS_CWX(M,NWN)
          ENDDO
        ENDDO
!
!---    Global injection well mass balance residuals  ---
!
!        PRINT *,'QM_CWX = ',(QM_CWX(M),M=1,3),' ID = ',ID
        CALL MPI_ALLREDUCE( QM_CWX,QMG_CWX,3,MPI_REAL8,MPI_SUM,
     &    MPI_COMM_WORLD,IERR )
!        PRINT *,'QMG_CWX = ',(QMG_CWX(M),M=1,3),' ID = ',ID
!
!---    Consider current well pressure and well pressure limit  ---
!
        IF( ICHK_CWX.EQ.0 ) THEN
          ICHK_CWX = 1
!
!---      Hold pressure controlled option  ---
!
          IF( ID_CW(8,NCW).EQ.1 .AND. NITER.GE.4 ) THEN
            P_CW(2,NCW) = PL_CW(NCW)
            DNR_CW(NCW) = 1.D-1
            P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
            ID_CW(8,NCW) = 1
            ISUB_LOG = ISUB_LOG-1
            RETURN
!
!---      Well-limit pressure insufficient to produce specified
!         rate, well is pressure controlled  ---
!
          ELSEIF( QMG_CWX(2).LT.VAR_CW2X(2) ) THEN
            P_CW(2,NCW) = PL_CW(NCW)
            DNR_CW(NCW) = 1.D-1
            P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
            ID_CW(8,NCW) = 1
            ISUB_LOG = ISUB_LOG-1
            RETURN
!
!---      Well-limit pressure and current well pressure yield
!         flow above specified rate, well is flow controlled.  Find
!         well pressure that yields positive flow below specified rate ---
!
          ELSEIF( QMG_CWX(1).GE.VAR_CW2X(1) ) THEN
            DP_CWX = 1.D+1*DP_CWX
            P_CWY(1) = P_CWY(1) - DP_CWX
            ICHK_CWX = 0
            CYCLE
!
!---      Well-limit pressure yields flow above specified rate,
!         and current well pressure yields positive flow, well is flow
!         controlled  ---
!
          ELSEIF( QMG_CWX(1).GT.EPSL ) THEN
            DNR_CW(NCW) = 1.D-1
            P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
            ID_CW(8,NCW) = 0
            ISUB_LOG = ISUB_LOG-1
            RETURN
!
!---      Well limit pressure yields flow above specified rate,
!         and current well pressure yields zero flow, use bisection to
!         determine a well pressure that yields the specified rate  ---
!
          ELSE
            P_CWY(3) = 5.D-1*(P_CWY(1)+P_CWY(2))
            ML = 3
            MU = 3
            CYCLE
          ENDIF
!
!---      Use bisection to determine a well pressure that yields
!         the specified rate  ---
!
        ELSE
          IF( (ABS(QMG_CWX(3)-VAR_CW2X(3)).LT.(1.D-3*VAR_CW2X(3))
     &      .OR. (P_CWY(2)-P_CWY(1)).LT.1.D-1) .OR. NC.GE.32 ) THEN
            P_CW(2,NCW) = P_CWY(3)
            DNR_CW(NCW) = 1.D-1
            P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
            ID_CW(8,NCW) = 0
            ISUB_LOG = ISUB_LOG-1
            RETURN
          ELSEIF( ((QMG_CWX(3)-VAR_CW2X(3)).LE.0.D+0 .AND.
     &      (QMG_CWX(1)-VAR_CW2X(1)).LE.0.D+0) .OR. 
     &      ((QMG_CWX(3)-VAR_CW2X(3)).GT.0.D+0 .AND.
     &      (QMG_CWX(1)-VAR_CW2X(1)).GT.0.D+0) ) THEN
            P_CWY(1) = P_CWY(3)
            QMG_CWX(1) = QMG_CWX(3)
            VAR_CW2X(1) = VAR_CW2X(3)
            P_CWY(3) = 5.D-1*(P_CWY(1)+P_CWY(2))
            CYCLE
          ELSE
            P_CWY(2) = P_CWY(3)
            QMG_CWX(2) = QMG_CWX(3)
            VAR_CW2X(2) = VAR_CW2X(3)
            P_CWY(3) = 5.D-1*(P_CWY(1)+P_CWY(2))
            CYCLE
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INJP_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCB_COUP_WELL
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
!
!     STOMPX-CO2
!
!     Modify Jacobian matrix for the coupled-well equations
!     and load Jacobian matrix for the coupled-well equations
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 December 2021.
!







!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
!

!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FDVP
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)





!
!
!  Include file for Fortran use of the Mat package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!




!
!
!  Include file for Fortran use of the Vec package in PETSc
!




!
!
!  Include file for Fortran use of the AO (application ordering) package in PETSc
!




!
!
!  Include file for Fortran use of the IS (index set) package in PETSc
!




!
!
!  Part of the base include file for Fortran use of PETSc.
!  Note: This file should contain only define statements and
!  not the declaration of variables.

! No spaces for #defines as some compilers (PGI) also adds
! those additional spaces during preprocessing - bad for fixed format
!
























































!
!  Include file for Fortran use of the PetscViewer package in PETSc
!













!
!  No includes needed for logging

!
!  Include file for Fortran use of the Bag package in PETSc
!







!
! The real*8,complex*16 notatiton is used so that the
! PETSc double/complex variables are not affected by
! compiler options like -r4,-r8, sometimes invoked
! by the user. NAG compiler does not like integer*4,real*8
























!
! Fortran does not support unsigned, though ISO_C_BINDING
! supports INTEGER(KIND=C_SIZE_T). We don't use that here
! only to avoid importing the module.





!

!






!

!


!



!
!     Macro for templating between real and complex
!










!
!    Allows the matrix Fortran Kernels to work with single precision
!    matrix data structures
!

!
!     PetscLogDouble variables are used to contain double precision numbers
!     that are not used in the numerical computations, but rather in logging,
!     timing etc.
!


!
!     Macros for error checking
!




















!
!  Include file for Fortran use of the type(tPetscViewer) package in PETSc
!




































































!
!  Matrix types
!


!
! MatMFFDType values
!



!
! MatSolverTypes
!


!
! GPU Storage Formats for CUSPARSE
!



!
! GPU Storage Formats for HIPSPARSE
!



!
! sparsity reducing ordering for STRUMPACK
!



!
!
!  Include file for Fortran use of the type(tVec) package in PETSc
!

!
!
!  Include file for Fortran use of the KSP package in PETSc
!




!
!
!  Include file for Fortran use of the PC (preconditioner) package in PETSc
!




!
!
!  Include file for Fortran use of the type(tMat) package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!


!
!  Include file for Fortran use of the DM package in PETSc
!




!
!
!  Include file for Fortran use of the type(tIS) (index set) package in PETSc
!

!
!
!  Include file for Fortran use of the type(tVec) package in PETSc
!

!
!
!  Include file for Fortran use of the type(tMat) package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!



















!
! GAMG types
!



!
! GAMG classical types
!



!
! Various preconditioners
!









!
!  Various Krylov subspace methods
!

!
!  Various Initial guesses for Krylov subspace methods
!


!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCB_COUP_WELL'

!
!---  PETSc solver  ---
!
!
!---  Loop over coupled wells ---
!
      DO NCW = 1,N_CW
        NC = 0
!
!---    Loop over coupled-well well nodes ---
!
        NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
        DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
          IF( N.EQ.0 ) CYCLE
          NMD = (IXP(N)-1)*ISVC
!
!---      Water mass balance equation at field node ---
!
          MP = NMD + IEQW
!          IROW = NMD + IEQW
!
!---      Change in water mass flux into field node with respect
!         to change in field node primary variables  ---
!
          MA = 0
          DO M = 1,ISVC
            MCOL = KLU(MP-IEQ_OFFSET,M+MA)
            DLU(MCOL) = DLU(MCOL) + 
     &        (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNR(M,N)
          ENDDO
!
!---      Change in water mass flux into field node with respect
!         to change in coupled-well pressure  ---
!
          MCOL = KLU1_CW(IEQW,NWN,NCW)
          MX = ISVC+2
          DLU(MCOL) = DLU(MCOL) + 
     &      (FXW_CW(1,NWN)-FXW_CW(MX,NWN))/DNR_CW(NCW)
!
!---      Water mass flux into field node, kg/s  ---
!
          BUFFER = FXW_CW(1,NWN)
          IROW_P = MP - 1
          CALL VecSetValues( F_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,IERR )
          RSDL(IEQW,N) = RSDL(IEQW,N) + BUFFER
!
!---      CO2 mass balance equation at field node ---
!
          MP = NMD + IEQA
!
!---      Change in CO2 mass flow with respect
!         to change in field node primary variables  ---
!
          MA = 0 
          DO M = 1,ISVC
            MCOL = KLU(MP-IEQ_OFFSET,M+MA)
            DLU(MCOL) = DLU(MCOL) + 
     &        (FXA_CW(1,NWN)-FXA_CW(M+1,NWN))/DNR(M,N)
          ENDDO
!
!---      Change in CO2 mass flow with respect
!         to change in coupled-well pressure  ---
!
          MCOL = KLU1_CW(IEQA,NWN,NCW)
          MX = ISVC+2
          DLU(MCOL) = DLU(MCOL) + 
     &      (FXA_CW(1,NWN)-FXA_CW(MX,NWN))/DNR_CW(NCW)
!
!---      CO2 mass flow into field node, kg/s  ---
!
          BUFFER = FXA_CW(1,NWN)
          IROW_P = MP - 1
          CALL VecSetValues( F_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,IERR )
          RSDL(IEQA,N) = RSDL(IEQA,N) + BUFFER
!
!---      Skip for isobrine simulations  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt mass balance equation at field node ---
!
            MP = NMD + IEQS
!            IROW = NMD + IEQS
!
!---        Change in salt mass flux into field node with respect
!           to change in field node primary variables  ---
!
            MA = 0
            DO M = 1,ISVC
              MCOL = KLU(MP-IEQ_OFFSET,M+MA)
              DLU(MCOL) = DLU(MCOL) + 
     &          (FXS_CW(1,NWN)-FXS_CW(M+1,NWN))/DNR(M,N)
            ENDDO
!
!---        Change in salt mass flux into field node with respect
!           to change in coupled-well pressure  ---
!
            MCOL = KLU1_CW(IEQS,NWN,NCW)
            MX = ISVC+2
            DLU(MCOL) = DLU(MCOL) + 
     &        (FXS_CW(1,NWN)-FXS_CW(MX,NWN))/DNR_CW(NCW)
!
!---        Salt mass flux into field node, kg/s  ---
!
            BUFFER = FXS_CW(1,NWN)
            IROW_P = MP - 1
            CALL VecSetValues( F_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &        IERR )
            RSDL(IEQS,N) = RSDL(IEQS,N) + BUFFER
          ENDIF
        ENDDO
!
!---    Coupled-well mass balance, coupled-well equations located
!       on last processor  ---
!
        IROW = JM_CW(NCW)
        IF( ID.EQ.(NP-1) ) THEN
          BUFFER = -RS_CW(1,NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) BUFFER = 0.D+0
          IROW_P = IROW - 1
          CALL VecSetValues( F_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,IERR )
!
!---      Change in coupled-well mass balance with respect to
!         change in coupled-well pressure  ---
!
          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
          MX = (NWFX*ISVC)+2
          NC = NC + 1
          MCOL = KLU2_CW(NC,NCW)
          DLU(MCOL) = DLU(MCOL) + 
     &      (RS_CW(MX,NCW)-RS_CW(1,NCW))/DNR_CW(NCW)
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 1.D+0
!
!---      Loop over field nodes with coupled-well nodes ---
!
          DO NWF = ID_CW(5,NCW),ID_CW(6,NCW)
            N = IWF_CW(NWF)
            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
!
!---        Change in coupled-well mass balance with respect to
!           change in field node primary variables  ---
!
            DO M = 1,ISVC
              NC = NC + 1
              MCOL = KLU2_CW(NC,NCW)
              DLU(MCOL) = DLU(MCOL) + RS_CW(MX+M,NCW)
!
!---          Pressure controlled coupled well  ---
!
              IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 0.D+0
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCB_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_COUP_WELL
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
!
!     STOMP-CO2
!
!     Coupled-well equation residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE OUTPU
      USE JACOB
      USE GRID
      USE FILES
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER, DIMENSION(1:2) :: IVARX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL_COUP_WELL'
!
!---  Zero maximum residuals  ---
!
      RSD_CW = 0.D+0
      NSD_CW = 0
!
!---  Loop over coupled wells ---
!
      DO NCW = 1,N_CW
!
!---    Injection well ---
!
        IF( IT_CW(NCW).GT.0 ) THEN
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) THEN
            RSDX = 0.D+0
!
!---      Flow controlled coupled well  ---
!
          ELSE
            RSDFX = 0.D+0
            RSDPX = ABS(DP_CW(NCW))/(ABS(P_CW(2,NCW))+PATM)
            RSDX = MAX( RSDPX,RSDFX )
          ENDIF
!
!---    Withdrawl well ---
!
        ELSEIF( IT_CW(NCW).LT.0 ) THEN
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) THEN
            RSDX = 0.D+0
!
!---      Flow controlled coupled well  ---
!
          ELSE
            RSDFX = 0.D+0
            RSDPX = ABS(DP_CW(NCW))/(ABS(P_CW(2,NCW))+PATM)
            RSDX = MAX( RSDPX,RSDFX )
          ENDIF
        ENDIF
        IF( RSDX.GT.RSD_CW ) THEN
          RSD_CW = RSDX
          NSD_CW = NCW
        ENDIF
      ENDDO
      IF( RSD_CW.GT.RSDMX ) ICNV = 2
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) THEN
        IF( RSD_CW.GE.1.D+2 ) THEN
          IF( ID.EQ.0 ) THEN
            PRINT *,'           ---  Excessive Residual  ---'
            PRINT *,'  Coupled Well Maximum Residual = ',RSD_CW,
     &        ': Coupled Well Number = ',NSD_CW
          ENDIF
          IVARX(1) = -7
        ELSE
          IF( ID.EQ.0 ) THEN
            PRINT *,'           ---  Convergence Failure  ---'
            PRINT *,'  Coupled Well Maximum Residual = ',RSD_CW,
     &        ': Coupled Well Number = ',NSD_CW
          ENDIF
          IVARX(1) = -8
        ENDIF
!
!---    Write a convergence failure index of -5 in the NSTEP location
!       of the output.bin file plus write the well number and 
!       maximum coupled-well residual  ---
!
        OFFSET = IOFFSET_REF
        NVAR = 2
        IVARX(2) = NSD_CW
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,IVARX,NVAR,
     &    MPI_INTEGER,STATUS,IERR)
        OFFSET = OFFSET + NVAR*NBYTI
        IOFFSET_REF = IOFFSET_REF + NVAR*NBYTI
        NVAR = 1
        IF( ID.EQ.0 ) CALL MPI_FILE_WRITE_AT( IWR,OFFSET,RSD_CW,NVAR,
     &    MPI_REAL8,STATUS,IERR)
        OFFSET = OFFSET + NVAR*NBYTR
        IOFFSET_REF = IOFFSET_REF + NVAR*NBYTR
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_COUP_WELL
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
!
!     STOMP-CO2
!
!     Load old-time-step values for coupled-well arrays.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LDO_COUP_WELL'
!
!---  Loop over coupled wells ---
!
      DO NCW = 1,N_CW
        P_CW(1,NCW) = P_CW(2,NCW)
        PL_CW(NCW) = P_CW(2,NCW)
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_COUP_WELL( NSL )
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
!
!     STOMP-EOR
!
!     Compute solute source transport terms for coupled wells.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 26 January 2022.
!







!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
!

!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)





!
!
!  Include file for Fortran use of the type(tMat) package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!

!
!
!  Include file for Fortran use of the type(tVec) package in PETSc
!

!
!
!  Include file for Fortran use of the type(tKSP) package in PETSc
!


!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORT_COUP_WELL'
!
!---  Loop over coupled wells  ---
!
      L1: DO NCW = 1,N_CW
!
!---    Skip for conservation or kinetic components  ---
!
        IF( NSL.GT.NSOLU ) CYCLE L1
!
!---    Check for solute in well  ---
!
        IFIND = 0
        L2: DO M = 1,NSOLU
          IF( ISOLU_CW(M,NCW).EQ.0 ) EXIT L2
          IF( ISOLU_CW(M,NCW).EQ.NSL ) THEN
            IFIND = 1
            NC = M
            EXIT L2
          ENDIF
        ENDDO L2
!
!---    Solute not found in well, cycle to next well ---
!
        IF( IFIND.EQ.0 ) CYCLE L1
!
!---    Coupled well time interval ---
!
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---    Cyclic time periods  ---
!
        IF( ICC_CW(NCW).EQ.1 ) THEN
!
!---      Loop over the coupled well time periods, 
!         to find the final well time  ---
!
          NTX = 0
          L3: DO NTP = 1,IM_CW(NCW)
            NTX = NTX + IMP_CW(NTP,NCW)
          ENDDO L3
!
!---      Determine time with the cyclic time period  ---
!
          TMZ = MOD( TM,VAR_CW(1,NTX,NCW) )
          IF( TM.GT.VAR_CW(1,NTX,NCW) ) THEN
            IF( TMZ.LT.EPSL ) TMZ = VAR_CW(1,NTX,NCW)
          ENDIF
        ENDIF
!
!---    Coupled well is inactive  ---
!
        IF( TMZ.LE.VAR_CW(1,1,NCW) ) CYCLE L1
!
!---    Loop over the coupled well time periods  ---
!
        NS = 1
        IFIND = 0
        L4: DO NTP = 1,IM_CW(NCW)
!
!---      Coupled well time period only has one time (start time)  ---
!
          IF( IMP_CW(NTP,NCW).EQ.1 ) THEN
!
!---        Time prior to start time, coupled well is inactive,
!           cycle to next well  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) ) CYCLE L1
!
!---        Time after start time, coupled well is active  ---
!
            VARC_CWX = VARC_CW(NSL,1,NCW)
            IFIND = 1
            EXIT L4
!
!---      Coupled well time period has multiple times  ---
!
          ELSE
            NE = NS + IMP_CW(NTP,NCW) - 1
!
!---        Time outside of coupled well time period, go to next 
!           coupled well time period  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) .OR. 
     &        TMZ.GT.VAR_CW(1,NE,NCW) ) THEN
              NS = NS + IMP_CW(NTP,NCW)
              CYCLE L4
            ENDIF
!
!---        Coupled well time period has multiple time points, use  
!           linear interpolation of well parameters between 
!           time points  ---
!
            L5: DO M = 2,IMP_CW(NTP,NCW)
              MX = NS + M - 1
              IF( TMZ.LE.VAR_CW(1,MX,NCW) ) THEN
                TD_CW = VAR_CW(1,MX,NCW)-VAR_CW(1,MX-1,NCW)
                DT_CW = MIN( VAR_CW(1,MX,NCW)-TMZ,DT )
                TF_CWX = (TMZ-VAR_CW(1,MX-1,NCW))/TD_CW
                VARC_CWX = VARC_CW(NSL,MX-1,NCW) + 
     &            TF_CWX*(VARC_CW(NSL,MX,NCW)-VARC_CW(NSL,MX-1,NCW))
                IFIND = 1
                EXIT L4
              ENDIF
            ENDDO L5
          ENDIF
          NS = NS + IMP_CW(NTP,NCW)
        ENDDO L4
!
!---    Coupled well is inactive, cycle to next well  ---
!
        IF( IFIND.EQ.0 ) CYCLE L1
!
!---    Loop over coupled-well nodes  ---
!
        L6: DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
!
!---      Skip for off-processor well nodes, inactive nodes
!         and ghost cells ---
!
          IF( N.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
          IROW = IXP(N)
          SORTX = 0.D+0
!
!---      Injection well (volumetric fluxes are positive from well) ---
!
          IF( IT_CW(NCW).GT.0 ) THEN






            BUFFER = Q_CW(1,NWN)*VARC_CWX
            IROW_P = IROW - 1
            CALL VecSetValues( T_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &        IERR )

!
!---      Production well (volumetric fluxes are positive into well) ---
!
          ELSEIF( IT_CW(NCW).LT.0 ) THEN
!
!---        Solute produced via aqueous phase production  ---
!
            CLX = PORD(2,N)*SL(2,N)
            IF( CLX.GT.EPSL ) SORTX = SORTX + Q_CW(2,NWN)/CLX
!
!---        Solute produced via gas phase production  ---
!
            CGX = PORD(2,N)*SG(2,N)
            IF( CGX.GT.EPSL ) SORTX = SORTX + Q_CW(3,NWN)/CGX
          ENDIF
!
!---      Load Jacobian  ---
!







          IROW_P = IROW - 1
          ICOL_P = IXP(N) - 1
          BUFFER = SORTX
          CALL MatSetValue( T_MAT,IROW_P,ICOL_P,BUFFER,ADD_VALUES,IERR )

        ENDDO L6
      ENDDO L1
      NEQ = NSL - NSOLU
!
!---  Loop over coupled wells ---
!
      L11: DO NCW = 1,N_CW
!
!---    Skip for passive solutes  ---
!
        IF( NSL.LE.NSOLU ) CYCLE L11
!
!---    Skip for non-conservation solutes  ---
!
        IF( NEQ.LT.1 .OR. NEQ.GT.NEQC ) CYCLE L11
!
!---    Check for conservation component solute in well  ---
!
        IF( NEQ.GT.0 .AND. NEQ.LE.NEQC ) THEN
          IF( ISOLC_CW(NEQ,NCW).EQ.0 ) CYCLE L11
        ENDIF
!
!---    Coupled well time interval ---
!
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---    Cyclic time periods  ---
!
        IF( ICC_CW(NCW).EQ.1 ) THEN
!
!---      Loop over the coupled well time periods, 
!         to find the final well time  ---
!
          NTX = 0
          L13: DO NTP = 1,IM_CW(NCW)
            NTX = NTX + IMP_CW(NTP,NCW)
          ENDDO L13
!
!---      Determine time with the cyclic time period  ---
!
          TMZ = MOD( TM,VAR_CW(1,NTX,NCW) )
          IF( TM.GT.VAR_CW(1,NTX,NCW) ) THEN
            IF( TMZ.LT.EPSL ) TMZ = VAR_CW(1,NTX,NCW)
          ENDIF
        ENDIF
!
!---    Coupled well is inactive  ---
!
        IF( TMZ.LE.VAR_CW(1,1,NCW) ) CYCLE L11
!
!---    Loop over the coupled well time periods  ---
!
        NS = 1
        IFIND = 0
        L14: DO NTP = 1,IM_CW(NCW)
!
!---      Coupled well time period only has one time (start time)  ---
!
          IF( IMP_CW(NTP,NCW).EQ.1 ) THEN
            ITS_CWX = MOD( ITS_CW(NTP,NCW),100 )
!
!---        Time prior to start time, coupled well is inactive,
!           cycle to next well  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) ) CYCLE L11
!
!---        Time after start time, coupled well is active  ---
!
            VARC_CWX = 0.D+0
!
!---        Aqueous species only
!
            IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
              DO NSPCW = 1,NSP_CW(NCW)
                NSP = ISPC_CW(NSPCW,NCW)
                IF( NSP.EQ.0 ) EXIT
                IF( NSP.LE.NSPL ) THEN
                  DO IX = 1,IEQ_C(1,NEQ)
                    NSP_C = IEQ_C(IX+1,NEQ)
                    IF( NSP_C.EQ.NSP ) THEN
                      VARC_CWX = VARC_CWX + 
     &                  EQ_C(IX,NEQ)*VARSP_CW(NSPCW,1,NCW)
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF   
              ENDDO             
!
!---        Gas species only
!      
            ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN
              DO NSPCW = 1,NSP_CW(NCW)
                NSP = ISPC_CW(NSPCW,NCW)
                IF( NSP.EQ.0 ) EXIT
                IF( NSP.GT.NSPL+NSPS .AND. 
     &            NSP.LE.NSPL+NSPS+NSPG ) THEN
                  DO IX = 1,IEQ_C(1,NEQ)
                    NSP_C = IEQ_C(IX+1,NEQ)
                    IF( NSP_C.EQ.NSP ) THEN
                      VARC_CWX = VARC_CWX + 
     &                  EQ_C(IX,NEQ)*VARSP_CW(NSPCW,1,NCW)
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF   
              ENDDO             
            ENDIF
            IFIND = 1
            EXIT L14
!
!---      Coupled well time period has multiple times  ---
!
          ELSE
            NE = NS + IMP_CW(NTP,NCW) - 1
            ITS_CWX = MOD( ITS_CW(NTP,NCW),100 )
!
!---        Time outside of coupled well time period, go to next 
!           coupled well time period  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) .OR. 
     &        TMZ.GT.VAR_CW(1,NE,NCW) ) THEN
              NS = NS + IMP_CW(NTP,NCW)
              CYCLE L14
            ENDIF
!
!---        Coupled well time period has multiple time points, use  
!           linear interpolation of well parameters between 
!           time points  ---
!
            L15: DO M = 2,IMP_CW(NTP,NCW)
              MX = NS + M - 1
              IF( TMZ.LE.VAR_CW(1,MX,NCW) ) THEN
                TD_CW = VAR_CW(1,MX,NCW)-VAR_CW(1,MX-1,NCW)
                DT_CW = MIN( VAR_CW(1,MX,NCW)-TMZ,DT )
                TF_CWX = (TMZ-VAR_CW(1,MX-1,NCW))/TD_CW
!
!---            Time after start time, coupled well is active  ---
!
                VARC_CWX = 0.D+0
!
!---            Aqueous species only
!
                IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
                  DO NSPCW = 1,NSP_CW(NCW)
                    NSP = ISPC_CW(NSPCW,NCW)
                    IF( NSP.EQ.0 ) EXIT
                    IF( NSP.LE.NSPL ) THEN
                      DO IX = 1,IEQ_C(1,NEQ)
                        NSP_C = IEQ_C(IX+1,NEQ)
                        IF( NSP_C.EQ.NSP ) THEN
                          VARC_CWX = VARC_CWX + EQ_C(IX,NEQ)*
     &                      (VARSP_CW(NSPCW,MX-1,NCW) + TF_CWX*
     &                      (VARSP_CW(NSPCW,MX,NCW) - 
     &                      VARSP_CW(NSPCW,MX-1,NCW)))
                          EXIT
                        ENDIF
                      ENDDO
                    ENDIF   
                  ENDDO             
!
!---            Gas species only
!
                ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN
                  DO NSPCW = 1,NSP_CW(NCW)
                    NSP = ISPC_CW(NSPCW,NCW)
                    IF( NSP.EQ.0 ) EXIT
                    IF( NSP.GT.NSPL+NSPS .AND. 
     &                NSP.LE.NSPL+NSPS+NSPG ) THEN
                      DO IX = 1,IEQ_C(1,NEQ)
                        NSP_C = IEQ_C(IX+1,NEQ)
                        IF( NSP_C.EQ.NSP ) THEN
                          VARC_CWX = VARC_CWX + EQ_C(IX,NEQ)*
     &                    (VARSP_CW(NSPCW,MX-1,NCW) + TF_CWX*
     &                    (VARSP_CW(NSPCW,MX,NCW) - 
     &                    VARSP_CW(NSPCW,MX-1,NCW)))
                          EXIT
                        ENDIF
                      ENDDO
                    ENDIF   
                  ENDDO             
                ENDIF
                IFIND = 1
                EXIT L14
              ENDIF
            ENDDO L15
          ENDIF
          NS = NS + IMP_CW(NTP,NCW)
        ENDDO L14
!
!---    Coupled well is inactive, cycle to next well  ---
!
        IF( IFIND.EQ.0 ) CYCLE L11
!
!---    Loop over coupled-well nodes  ---
!
        L16: DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
          IF( IXP(N).EQ.0 ) CYCLE
          IROW = IXP(N)
          SORTX = 0.D+0
!
!---      Injection well (volumetric fluxes are positive from well) ---
!
!         Q_CW(1,NWN) - total volumetric flux, m^3/s
!         Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!         Q_CW(3,NWN) - gas volumetric flux, m^3/s
!
          IF( IT_CW(NCW).GT.0 ) THEN
!
!---        Aqueous species only
!
            IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN






              BUFFER = Q_CW(2,NWN)*VARC_CWX
              IROW_P = IROW - 1
              CALL VecSetValues( T_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &         IERR )

!
!---        Gas species only
!
            ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN






              BUFFER = Q_CW(3,NWN)*VARC_CWX
              IROW_P = IROW - 1
              CALL VecSetValues( T_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &         IERR )

            ENDIF
!
!---      Production well (volumetric fluxes are positive into well) ---
!
          ELSEIF( IT_CW(NCW).LT.0 ) THEN
!
!---        Solute produced via aqueous phase production  ---
!
            CLX = PORD(2,N)*SL(2,N)
            IF( CLX.GT.EPSL ) SORTX = SORTX + Q_CW(2,NWN)/CLX
!
!---        Solute produced via gas phase production  ---
!
            CGX = PORD(2,N)*SG(2,N)
            IF( CGX.GT.EPSL ) SORTX = SORTX + Q_CW(3,NWN)/CGX
          ENDIF
!
!---      Load Jacobian  ---
!







          IROW_P = IROW - 1
          ICOL_P = IXP(N) - 1
          BUFFER = SORTX
          CALL MatSetValue( T_MAT,IROW_P,ICOL_P,BUFFER,ADD_VALUES,IERR )

        ENDDO L16
      ENDDO L11
      NEQ = NEQ - NEQC
!
!---  Loop over coupled wells ---
!
      L21: DO NCW = 1,N_CW
!
!---    Skip for passive solutes  ---
!
        IF( NSL.LE.NSOLU ) CYCLE L21
!
!---    Skip for non-kinetic solutes  ---
!
        IF( NEQ.LT.1 .OR. NEQ.GT.NEQK ) CYCLE L21
!
!---    Coupled well time interval ---
!
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
!
!---    Cyclic time periods  ---
!
        IF( ICC_CW(NCW).EQ.1 ) THEN
!
!---      Loop over the coupled well time periods, 
!         to find the final well time  ---
!
          NTX = 0
          L23: DO NTP = 1,IM_CW(NCW)
            NTX = NTX + IMP_CW(NTP,NCW)
          ENDDO L23
!
!---      Determine time with the cyclic time period  ---
!
          TMZ = MOD( TM,VAR_CW(1,NTX,NCW) )
          IF( TM.GT.VAR_CW(1,NTX,NCW) ) THEN
            IF( TMZ.LT.EPSL ) TMZ = VAR_CW(1,NTX,NCW)
          ENDIF
        ENDIF
!
!---    Coupled well is inactive  ---
!
        IF( TMZ.LE.VAR_CW(1,1,NCW) ) CYCLE L21
!
!---    Loop over the coupled well time periods  ---
!
        NS = 1
        IFIND = 0
        L24: DO NTP = 1,IM_CW(NCW)
!
!---      Coupled well time period only has one time (start time)  ---
!
          IF( IMP_CW(NTP,NCW).EQ.1 ) THEN
            ITS_CWX = MOD( ITS_CW(NTP,NCW),100 )
!
!---        Time prior to start time, coupled well is inactive,
!           cycle to next well  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) ) CYCLE L21
!
!---        Time after start time, coupled well is active  ---
!
            VARC_CWX = 0.D+0
!
!---        Aqueous species only
!
            IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
              DO NSPCW = 1,NSP_CW(NCW)
                NSP = ISPC_CW(NSPCW,NCW)
                IF( NSP.EQ.0 ) EXIT
                IF( NSP.LE.NSPL ) THEN
                  DO IX = 1,IEQ_K(1,NEQ)
                    NSP_K = IEQ_K(IX+1,NEQ)
                    IF( NSP_K.EQ.NSP ) THEN
                      VARC_CWX = VARC_CWX + 
     &                  EQ_K(IX,NEQ)*VARSP_CW(NSPCW,1,NCW)
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF   
              ENDDO             
!
!---        Gas species only
!
            ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN
              DO NSPCW = 1,NSP_CW(NCW)
                NSP = ISPC_CW(NSPCW,NCW)
                IF( NSP.EQ.0 ) EXIT
                IF( NSP.GT.NSPL+NSPS .AND. 
     &            NSP.LE.NSPL+NSPS+NSPG ) THEN
                  DO IX = 1,IEQ_K(1,NEQ)
                    NSP_K = IEQ_K(IX+1,NEQ)
                    IF( NSP_K.EQ.NSP ) THEN
                      VARC_CWX = VARC_CWX + 
     &                  EQ_K(IX,NEQ)*VARSP_CW(NSPCW,1,NCW)
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF   
              ENDDO             
            ENDIF
            IFIND = 1
            EXIT L24
!
!---      Coupled well time period has multiple times  ---
!
          ELSE
            NE = NS + IMP_CW(NTP,NCW) - 1
            ITS_CWX = MOD( ITS_CW(NTP,NCW),100 )
!
!---        Time outside of coupled well time period, go to next 
!           coupled well time period  ---
!
            IF( TMZ.LE.VAR_CW(1,NS,NCW) .OR. 
     &        TMZ.GT.VAR_CW(1,NE,NCW) ) THEN
              NS = NS + IMP_CW(NTP,NCW)
              CYCLE L24
            ENDIF
!
!---        Coupled well time period has multiple time points, use  
!           linear interpolation of well parameters between 
!           time points  ---
!
            L25: DO M = 2,IMP_CW(NTP,NCW)
              MX = NS + M - 1
              IF( TMZ.LE.VAR_CW(1,MX,NCW) ) THEN
                TD_CW = VAR_CW(1,MX,NCW)-VAR_CW(1,MX-1,NCW)
                DT_CW = MIN( VAR_CW(1,MX,NCW)-TMZ,DT )
                TF_CWX = (TMZ-VAR_CW(1,MX-1,NCW))/TD_CW
!
!---            Time after start time, coupled well is active  ---
!
                VARC_CWX = 0.D+0
!
!---            Aqueous species only
!
                IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
                  DO NSPCW = 1,NSP_CW(NCW)
                    NSP = ISPC_CW(NSPCW,NCW)
                    IF( NSP.EQ.0 ) EXIT
                    IF( NSP.LE.NSPL ) THEN
                      DO IX = 1,IEQ_K(1,NEQ)
                        NSP_K = IEQ_K(IX+1,NEQ)
                        IF( NSP_K.EQ.NSP ) THEN
                          VARC_CWX = VARC_CWX + EQ_K(IX,NEQ)*
     &                      (VARSP_CW(NSPCW,MX-1,NCW) + TF_CWX*
     &                      (VARSP_CW(NSPCW,MX,NCW) - 
     &                      VARSP_CW(NSPCW,MX-1,NCW)))
                          EXIT
                        ENDIF
                      ENDDO
                    ENDIF   
                  ENDDO             
!
!---            Gas species only
!
                ELSEIF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN
                  DO NSPCW = 1,NSP_CW(NCW)
                    NSP = ISPC_CW(NSPCW,NCW)
                    IF( NSP.EQ.0 ) EXIT
                    IF( NSP.GT.NSPL+NSPS .AND. 
     &                NSP.LE.NSPL+NSPS+NSPG ) THEN
                      DO IX = 1,IEQ_K(1,NEQ)
                        NSP_K = IEQ_K(IX+1,NEQ)
                        IF( NSP_K.EQ.NSP ) THEN
                          VARC_CWX = VARC_CWX + EQ_K(IX,NEQ)*
     &                      (VARSP_CW(NSPCW,MX-1,NCW) + TF_CWX*
     &                      (VARSP_CW(NSPCW,MX,NCW) - 
     &                      VARSP_CW(NSPCW,MX-1,NCW)))
                          EXIT
                        ENDIF
                      ENDDO
                    ENDIF   
                  ENDDO             
                ENDIF
                IFIND = 1
                EXIT L24
              ENDIF
            ENDDO L25
          ENDIF
          NS = NS + IMP_CW(NTP,NCW)
        ENDDO L24
!
!---    Coupled well is inactive, cycle to next well  ---
!
        IF( IFIND.EQ.0 ) CYCLE L21
!
!---    Loop over coupled-well nodes  ---
!
        L26: DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          N = IWN_CW(NWN)
          IF( IXP(N).EQ.0 ) CYCLE
          IROW = IXP(N)
          SORTX = 0.D+0
!
!---      Injection well (volumetric fluxes are positive from well) ---
!
!         Q_CW(1,NWN) - total volumetric flux, m^3/s
!         Q_CW(2,NWN) - aqueous volumetric flux, m^3/s
!         Q_CW(3,NWN) - gas volumetric flux, m^3/s
!
          IF( IT_CW(NCW).GT.0 ) THEN
!
!---        Aqueous species only
!
            IF( ITS_CWX.GE.4 .AND. ITS_CWX.LE.6 ) THEN






              BUFFER = Q_CW(2,NWN)*VARC_CWX
              IROW_P = IROW - 1
              CALL VecSetValues( T_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &          IERR )

!
!---        Gas species only
!
            ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN






              BUFFER = Q_CW(3,NWN)*VARC_CWX
              IROW_P = IROW - 1
              CALL VecSetValues( T_RHS_VEC,1,IROW_P,BUFFER,ADD_VALUES,
     &          IERR )

            ENDIF
!
!---      Production well (volumetric fluxes are positive into well) ---
!
          ELSEIF( IT_CW(NCW).LT.0 ) THEN
!
!---        Solute produced via aqueous phase production  ---
!
            CLX = PORD(2,N)*SL(2,N)
            IF( CLX.GT.EPSL ) SORTX = SORTX + Q_CW(2,NWN)/CLX
!
!---        Solute produced via gas phase production  ---
!
            CGX = PORD(2,N)*SG(2,N)
            IF( CGX.GT.EPSL ) SORTX = SORTX + Q_CW(3,NWN)/CGX
          ENDIF
!
!---      Load Jacobian  ---
!







          IROW_P = IROW - 1
          ICOL_P = IXP(N) - 1
          BUFFER = SORTX
          CALL MatSetValue( T_MAT,IROW_P,ICOL_P,BUFFER,ADD_VALUES,IERR )

        ENDDO L26
      ENDDO L21
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORT_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_COUP_WELL
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
!
!     STOMP-CO2
!
!     Update coupled-well pressure.  Injection wells are limited
!     by a high-pressure limit, and withdrawl wells are limited by a 
!     low-pressure limit.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_COUP_WELL'
!
!---  Well equations located on last processor, broadcast updates
!     to coupled-well pressures to all processors ---
!
      IDX = NP-1
      IF( ID.EQ.IDX ) THEN
        DO NCW = 1,N_CW





          M = (NUKFP(ID+1)-N_CW) + NCW
          DP_CW(NCW) = BLU(M)

        ENDDO
      ENDIF
      CALL MPI_BCAST( DP_CW,N_CW,MPI_REAL8,IDX,MPI_COMM_WORLD,IERR )
!
!---  Loop over coupled wells ---
!
      DO NCW = 1,N_CW
        DPX = 5.D+5
        DP_CWX = SIGN( MIN(ABS(DPX),ABS(DP_CW(NCW))),DP_CW(NCW) )
        P_CW(2,NCW) = P_CW(2,NCW) + DP_CWX
!
!---    Limit coupled-well pressure to upper limit for injection
!       wells or lower limit for withdrawl wells  ---
!
        IF( IT_CW(NCW).GT.0 ) THEN
          P_CW(2,NCW) = MIN( PL_CW(NCW),P_CW(2,NCW) )
        ELSEIF( IT_CW(NCW).LT.0 ) THEN
          P_CW(2,NCW) = MAX( PL_CW(NCW),P_CW(2,NCW) )
        ENDIF
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_COUP_WELL group  ---
!
      RETURN
      END



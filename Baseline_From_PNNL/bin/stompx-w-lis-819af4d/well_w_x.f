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
!     STOMPX-W (Water) Mode
!
!     Equilibrate coupled-well pressure with formation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 26 February 2023.
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
        IF( SG(2,N).GT.EPSL ) THEN
          IDLX = ID
          P_CW(2,NCW) = PG(2,N)
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
!     STOMPX-W (Water) Mode
!
!     Communicate to all processors the temperature of the 
!     field nodes containing well nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 26 February 2023.
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
!     STOMPX-W (Water) Mode
!
!     Mass flux of water between coupled-well nodes and 
!     field nodes.
!
!     Water mass balance residuals for injection type coupled wells.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 27 January 2023.
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
      REAL*8, DIMENSION(1:3) :: QML_CWX,QMG_CWX
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
!         Q_CW(1,NWN) - total (aqueous) volumetric flux, m^3/s
!
          Q_CW(1,NWN) = 0.D+0
!
!---      Loop over increment indices  ---
!
          DO M = 1,ISVC+2
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
          DO N = 2,3
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
              DO N = 2,3
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
!---    Volumetric injection or production rate, convert
!       volumetric rate to mass rate, based on bottom
!       pressure and temperature  ---
!
        IF( ABS(IT_CW(NCW)).EQ.2 ) THEN
          N = IWN_CW(ID_CW(4,NCW))
          P_CWX = P_CW(2,NCW)
          T_CWX = T(2,N)
          PL_CWX = P_CWX + PATM
          CALL WATLQD( T_CWX,PL_CWX,RHOL_CWX )
          VAR_CWX(2) = VAR_CWX(2)*RHOL_CWX
        ENDIF
!
!---    Loop over pressure adjustments ---
!
        PALOOP: DO
!
!---    Load water mass flux for use in RSDL_COUP_WELL ---
!
        FX_CW(NCW) = VAR_CWX(2)
!
!---    Load pressure limit for use in UPDT_COUP_WELL ---
!
        PL_CW(NCW) = VAR_CWX(3) - PATM
!
!---    Injection well, flop to pressure control for
!       for well pressures > limit  ---
!
        IF( IT_CW(NCW).GT.0 ) THEN
          IF( P_CW(2,NCW)-PL_CW(NCW).LT.EPSL ) THEN
            ID_CW(8,NCW) = 1
          ENDIF
!
!---    Withdrawl well, flop to pressure control for
!       for well pressures < limit  ---
!
        ELSE
          IF( PL_CW(NCW)-P_CW(2,NCW).GT.EPSL ) THEN
            ID_CW(8,NCW) = 1
          ENDIF
          VAR_CWX(2) = -VAR_CWX(2)
        ENDIF
!
!---    Loop over increment indices ---
!
        DO M = 1,ISVC+2
          MW = MCW(M)
          MF = MFD(M)
          N = IWN_CW(ID_CW(4,NCW))
          P_CWX = P_CW(MW,NCW)
          T_CWX = T(2,N)
          PL_CWX = P_CWX + PATM
          CALL WATLQD( T_CWX,PL_CWX,RHOL_CWX )
!
!---      Store bottom of coupled-well location in previous
!         coupled-well node location  ---
!
          XPX(1) = XTP_CW(2,ID_CW(1,NCW))
          YPX(1) = YTP_CW(2,ID_CW(1,NCW))
          ZPX(1) = ZTP_CW(2,ID_CW(1,NCW))
!
!---      Loop over the nodes in the coupled well ---
!
          DO NWN = ID_CW(4,NCW),ID_CW(3,NCW),-1
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
            P_CWX = P_CWX - (ZPX(2)-ZPX(1))*GRAV*RHOL_CWX
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
                PF_CW(NWF) = P_CWX - (ZPNX-ZPX(1))*GRAV*RHOL_CWX
              ENDIF
!
!---          Adjust the formation pressure to the coupled-well node
!             centroid  ---
!
              PGFX = PG(MF,N)
              PLFX = PL(MF,N) - (ZPX(2)-ZPNX)*GRAV*RHOL(MF,N)
              PL_CWX = MAX( P_CWX+PATM,0.D+0 )
              CALL WATLQD( T_CWX,PL_CWX,RHOL_CWX )
              CALL WATSP( T_CWX,PSW_CWX )
              CALL WATLQV( T_CWX,PL_CWX,PSW_CWX,VISL_CWX )
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
!---          Injection well, limit flow direction to be into
!             formation, with mass fluxes being positive into the
!             formation  ---
!
              IF( IT_CW(NCW).GT.0 ) THEN
                DPLX = MAX( P_CWX-MAX(PLFX,PGFX),0.D+0 )
!
!---          Withdrawl well, limit flow direction to be from
!             formation, with mass fluxes being positive into the
!             formation  ---
!
              ELSE
                DPLX = MIN( P_CWX-PLFX,0.D+0 )
              ENDIF
!
!---          Flux from well to formation, negative values indicate
!             flux form formation to well  ---
!
              FX_CWX = WI_CWX*RHOL_CWX*DPLX/VISL_CWX
              FXW_CW(M,NWN) = FX_CWX
!
!---          Volumetric well fluxes  ---
!
              IF( M.EQ.1 ) Q_CW(1,NWN) = FX_CWX/RHOL_CWX
!              PRINT *,'FXW_CW(',M,',',NWN,') = ',FXW_CW(M,NWN),
!     &          ' ID = ',ID
            ENDIF
!
!---        Store current coupled-well node location in previous
!           coupled-well node location  ---
!
            XPX(1) = XPX(2)
            YPX(1) = YPX(2)
            ZPX(1) = ZPX(2)
          ENDDO
        ENDDO
!
!---    CO2 mass balance residuals for injection type coupled well  ---
!
        NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
        NX = (NWFX*ISVC)+2
        RSL_CWX(1) = 0.D+0
        RSL_CWX(NX) = 0.D+0
        DO M = 1,3
          QML_CWX(M) = 0.D+0
          QMG_CWX(M) = 0.D+0
        ENDDO
!
!---    Loop over coupled-well nodes  ---
!
        DO NWN = ID_CW(3,NCW),ID_CW(4,NCW)
          RSL_CWX(1) = RSL_CWX(1) - FXW_CW(1,NWN)
          RSL_CWX(NX) = RSL_CWX(NX) - FXW_CW(ISVC+2,NWN)
          QML_CWX(1) = QML_CWX(1) + FXW_CW(1,NWN)
          QML_CWX(2) = QML_CWX(2) + Q_CW(1,NWN)
          QML_CWX(3) = QML_CWX(3) + FXW_CW(ISVC+2,NWN)
        ENDDO
!
!---    Global injection well mass balance residuals  ---
!
!        PRINT *,'QML_CWX(3) = ',QML_CWX(3),' ID = ',ID
        CALL MPI_ALLREDUCE( QML_CWX,QMG_CWX,3,MPI_REAL8,MPI_SUM,
     &    MPI_COMM_WORLD,IERR )
        QM_CW(1,NCW) = QMG_CWX(1)
        QM_CW(3,NCW) = QMG_CWX(2)
        QM_CWX = QMG_CWX(3)
!        PRINT *,'QMG_CWX(3) = ',QMG_CWX(3),' ID = ',ID
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
          DO M = 1,3
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
              FXW_CW(M,NWN) = 0.D+0
            ENDDO
          ENDDO
          CYCLE PALOOP
!
!---    Excessive well flow for pressure controlled well,
!       reduce well pressure  ---
!
        ELSEIF( ABS(ID_CW(8,NCW)).EQ.1 .AND. 
     &    ABS(QM_CW(1,NCW)).GT.ABS(VAR_CWX(2)) ) THEN
          DQ_CWX = 1.D+1*DQ_CWX
          DP_CWX = -DQ_CWX*PATM
          P_CW(2,NCW) = P_CW(2,NCW) + DP_CWX
          P_CW(3,NCW) = P_CW(2,NCW) + DNR_CW(NCW)
!
!---      Zero coupled-well fluxes ---
!
          DO M = 1,3
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
              FXW_CW(M,NWN) = 0.D+0
            ENDDO
          ENDDO
          CYCLE PALOOP
!
!---    Acceptable pressure reduction for pressure controlled well, 
!       switch to flow controlled well  ---
!
        ELSEIF( ID_CW(8,NCW).EQ.-1 .AND. 
     &    ABS(QM_CW(1,NCW)).LT.VAR_CWX(2) ) THEN
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
                RSL_CWX(M3) = RSL_CWX(M3) - FXW_CW(M2+1,NWN)
!
!---          If coupled-well node is outside the current field
!             node, use un-incremented fluxes  ---
!
              ELSE
                RSL_CWX(M3) = RSL_CWX(M3) - FXW_CW(1,NWN)
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
!     STOMPX-W (Water) Mode
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
!     STOMPX-W (Water) Mode
!
!     Modify Jacobian matrix for the coupled-well equations
!     and load Jacobian matrix for the coupled-well equations
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 December 2021.
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
      SUB_LOG(ISUB_LOG) = '/JCB_COUP_WELL'

!
!---  Lis solver  ---
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
          CALL lis_vector_set_value_f( 1,MP,BUFFER,
     &      F_RHS_VEC,IERR )
          RSDL(IEQW,N) = RSDL(IEQW,N) + BUFFER
        ENDDO
!
!---    Coupled-well mass balance, coupled-well equations located
!       on last processor  ---
!
        IROW = JM_CW(NCW)
        IF( ID.EQ.(NP-1) ) THEN
          BUFFER = -RS_CW(1,NCW)
          CALL lis_vector_set_value_f( 1,IROW,BUFFER,
     &      F_RHS_VEC,IERR )
!
!---      Pressure controlled coupled well  ---
!
          IF( ID_CW(8,NCW).EQ.1 ) THEN
            BUFFER = 0.D+0 
            CALL lis_vector_set_value_f( 0,IROW,BUFFER,
     &        F_RHS_VEC,IERR )
          ENDIF
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

!----------------------LIS Modules-------------------------------------!
!
      USE LIS_STOMP
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
            CALL lis_vector_set_value_f( 1,IROW,BUFFER,
     &        T_RHS_VEC,IERR )







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

          ICOL = IXP(N)
          BUFFER = SORTX
          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
     &      T_MAT,IERR )







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
              CALL lis_vector_set_value_f( 1,IROW,BUFFER,
     &          T_RHS_VEC,IERR )







!
!---        Gas species only
!
            ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN

              BUFFER = Q_CW(3,NWN)*VARC_CWX
              CALL lis_vector_set_value_f( 1,IROW,BUFFER,
     &          T_RHS_VEC,IERR )







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

          ICOL = IXP(N)
          BUFFER = SORTX
          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
     &      T_MAT,IERR )







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
              CALL lis_vector_set_value_f( 1,IROW,BUFFER,
     &          T_RHS_VEC,IERR )







!
!---        Gas species only
!
            ELSEIF( ITS_CWX.GE.1 .AND. ITS_CWX.LE.3 ) THEN

              BUFFER = Q_CW(3,NWN)*VARC_CWX
              CALL lis_vector_set_value_f( 1,IROW,BUFFER,
     &          T_RHS_VEC,IERR )







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

          ICOL = IXP(N)
          BUFFER = SORTX
          CALL lis_matrix_set_value_f( 1,IROW,ICOL,BUFFER,
     &      T_MAT,IERR )







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

          M = JM_CW(NCW) - NUKFO(IDX+1)
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



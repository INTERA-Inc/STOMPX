!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP_EOR_M( SGX,SLX,SNX,CPGTX,CPNLX,IZN )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Compute the gas to total-liquid capillary pressure,
!     and aqueous to nonaqueous-liquid capillary pressure, given
!     the gas saturation, aqueous saturation, nonaqueous-liquid
!     saturation and contact angle.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, July 17, 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CAP_EOR_M'
!
!---  van Genuchten saturation function
!
      IF( ISCHR(IZN).EQ.1 ) THEN
        CALL CAP_VG( SGX,SLX,SNX,CPGTX,CPNLX,IZN )
!
!---  Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 ) THEN
        CALL CAP_BC( SGX,SLX,SNX,CPGTX,CPNLX,IZN )
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP_EOR_M group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_BO_M
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 20 May 2014
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE DUAL_POR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_BO_M'
      TOLX = 2.5D-2
!
!---  Phase options, compute phase condition   ---
!
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(EPSL,GCPP,IXP,IZ,N_DB,NFLD,NPHAZ_M,PATM,PG_M,PL_M,
!$NO_OMP&    PN_M,PORD_M,POSM_M,PSO_M,PVA_M,RHOG_M,RHOL_M,RHON_M,SG_M,
!$NO_OMP&    SL_M,SN_M,T_M,TABS,TMC_M,VOL_M,XGC_M,XGW_M,XLA_M,XLS_M,
!$NO_OMP&    XLW_M,XNC_M,YLS_M) 
!$NO_OMP&  PRIVATE(CPGTX,CPNLX,IZN,N,PBPX,POSMX,PSBX,PTX,PVAX,PVBX,
!$NO_OMP&    PX,RHOBX,RHOGX,RHOOX,SGNX,SLX,TKX,TMGOX,TMLX,VGX,
!$NO_OMP&    VOX,XBSMX,XBSX,XLAX,XLSX,XLWX)
!$NO_OMP&  FIRSTPRIVATE(TOLX)
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 100
        IZN = IZ(N)
        N_DB = N
!
!---    Absolute temperature  ---
!
        TKX = T_M(2,N)+TABS
!
!---    Phase condition #1  ---
!
!       SL_M = 1.0
!       SG_M+SN_M = 0.0
!       Aqueuous saturated
!
!       water mass - aqueous (system) pressure, Pa gauge - PL_M,PSO_M
!       gas mass - total gas mass - TMC_M(1,2,N)
!       oil mass - total oil mass - TMC_M(2,2,N)
!       salt mass - total salt brine mass fraction - YLS_M
!
        IF( MOD(NPHAZ_M(2,N),10).EQ.1 ) THEN
!
!---      System pressure  ---
!
          PX = PSO_M(2,N) + PATM
!
!---      Maximum solubility of salt in pure water  ---
!
          CALL SOL_LS( T_M(2,N),XBSMX )
!
!---      Salt mass fraction in brine  ---
!
          XBSX = MIN( YLS_M(2,N),XBSMX )
          CALL SP_B( T_M(2,N),XBSX,PSBX )
!
!---      Maximum solubility of CO2 in aqueous phase  ---
!
          PVAX = 0.D+0
          CALL EQUIL( T_M(2,N),PX,PVAX,XBSX,XLAX,XLSX,XLWX )
!
!---      Solubility of CO2 in aqueous phase  ---
!
          CALL EQUIL( T_M(2,N),PX,PVA_M(2,N),XBSX,XLA_M(2,N),XLS_M(2,N),
     &      XLW_M(2,N) )
!
!---      Brine density  ---
!
          CALL DENS_B( T_M(2,N),PX,XBSX,RHOBX )
!
!---      Aqueous density  ---
!
          CALL DENS_L( T_M(2,N),RHOBX,XLA_M(2,N),RHOL_M(2,N) )
!
!---      Aqueous mass  ---
!
          TMLX = RHOL_M(2,N)*PORD_M(2,N)*VOL_M(N)
!
!---      Gas + oil mass  ---
!
          TMGOX = TMC_M(1,2,N) + TMC_M(2,2,N)
!
!---      Gas + oil mass execeeds 1.e-5 of the
!         aqueous phase mass, transition to phase condition #2 or #3 ---
!
          IF( TMGOX.GT.(1.D-5*TMLX) ) THEN
            CALL SP_B( T_M(2,N),XBSX,PSBX )
            CALL VPL( PSBX,PVBX,SLX,T_M(2,N),IZN )
!
!---        Bubble-point pressure, assuming all gas is dissolved  ---
!
            RHOGX = 1.225D+0*GCPP(37,1)
            VGX = TMC_M(1,2,N)/RHOGX
            RHOOX = 1.D+3*GCPP(37,2)
            VOX = TMC_M(2,2,N)/RHOOX
            POSMX = VGX/VOX
            CALL BO_BP( PBPX,POSMX,T_M(2,N) )
!
!---        Pressure above bubble-point pressure  ---
!
            IF( PX.GT.PBPX ) THEN
              SN_M(2,N) = 1.D-4
              SL_M(2,N) = 1.D+0 - SN_M(2,N)
              SG_M(2,N) = 0.D+0
              POSM_M(2,N) = POSMX
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 2
!
!---        Pressure below bubble-point pressure  ---
!
            ELSE
              CALL BO_GOR( PX,POSM_M(2,N),T_M(2,N) )
              CALL BO_EQ_GAS( PX,PVBX,RHOG_M(2,N),T_M(2,N),XGC_M(1,2,N),
     &          XGW_M(2,N) )
              CALL BO_EQ_OIL( RHON_M(2,N),PX,POSM_M(2,N),T_M(2,N),
     &          XNC_M(1,2,N) )
              SGNX = 1.D-4
              SG_M(2,N) = SGNX*RHON_M(2,N)*(TMC_M(1,2,N)*XNC_M(2,2,N) - 
     &          TMC_M(2,2,N) + TMC_M(2,2,N)*XNC_M(2,2,N))/
     &          (RHON_M(2,N)*TMC_M(1,2,N)*XNC_M(2,2,N) - 
     &          RHON_M(2,N)*TMC_M(2,2,N)
     &          + RHON_M(2,N)*TMC_M(2,2,N)*XNC_M(2,2,N) + 
     &          RHOG_M(2,N)*TMC_M(2,2,N)*XGC_M(1,2,N))
              SN_M(2,N) = SGNX - SG_M(2,N)
              SL_M(2,N) = 1.D+0-SGNX
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 3
            ENDIF
          ENDIF
!
!---      Phase pressure from phase capillary pressures and system
!         pressure  ---
!
          PG_M(2,N) = PSO_M(2,N)
          PTX = PG_M(2,N) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_M(2,N) = PTX
            PL_M(2,N) = PN_M(2,N) - CPNLX
          ELSE
            PL_M(2,N) = PTX
            PN_M(2,N) = PL_M(2,N) + CPNLX
          ENDIF
!
!---    Phase condition #2  ---
!
!       SG_M = 0.0
!       SL_M > 0.0
!       SN_M > 0.0
!       Aqueous and nonaqueous liquid saturated
!
!       water mass - aqueous saturation, SL_M
!       gas mass - gas to oil ratio, POSM_M
!       oil mass - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!       salt mass - total salt brine mass fraction - YLS_M
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.2 ) THEN
          PX = PSO_M(2,N) + PATM
!
!---      Maximum solubility of salt in pure water  ---
!
          CALL SOL_LS( T_M(2,N),XBSMX )
!
!---      Salt mass fraction in brine  ---
!
          XBSX = MIN( YLS_M(2,N),XBSMX )
!
!---      Nonaqueous flash calculation  ---
!
          CALL SP_B( T_M(2,N),XBSX,PSBX )
          CALL VPL( PSBX,PVBX,SL_M(2,N),T_M(2,N),IZN )
!
!---      Bubble-point pressure as a function of gas to oil ratio
!         (m^3/m^3) and temperature (C)  ---
!
          CALL BO_BP( PBPX,POSM_M(2,N),T_M(2,N) )
!
!---      Pressure below bubble-point pressure  ---
!
          IF( PX.LT.(PBPX-TOLX*GCPP(38,2)) ) THEN
            SG_M(2,N) = 1.D-4
            SN_M(2,N) = SN_M(2,N)*(1.D+0-SG_M(2,N))/
     &        (SN_M(2,N)+SL_M(2,N))
            SL_M(2,N) = 1.D+0-SG_M(2,N)-SN_M(2,N)
!
!---      Pressure above bubble-point pressure  ---
!
          ELSE
            SG_M(2,N) = 0.D+0
            SN_M(2,N) = 1.D+0-SL_M(2,N)
          ENDIF
!
!---      Phase capillary pressures from phase saturations  ---
!
          CALL CAP_EOR_M( SG_M(2,N),SL_M(2,N),SN_M(2,N),
     &      CPGTX,CPNLX,IZN )
!
!---      Nonaqueous phases disappears, transition to 
!         phase condition #1, conserving aqueous pressure  ---
!
          IF( (SG_M(2,N)+SN_M(2,N)).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 1
!
!---      Aqueous phases disappear, transition to 
!         phase condition #4  ---
!
          ELSEIF( (SG_M(2,N)+SL_M(2,N)).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 4
!
!---      Aqueous phase disappears and gas phase appears, transition to 
!         phase condition #5  ---
!
          ELSEIF( SL_M(2,N).LT.EPSL .AND. SG_M(2,N).GT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 5
!
!---      Gas phase appears, transition to phase condition #3  ---
!
          ELSEIF( SG_M(2,N).GT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 3
!
!---      No transition from phase condition #2  ---
!
          ELSE
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 2
          ENDIF
!
!---      Phase pressure from phase capillary pressures and system
!         pressure  ---
!
          PG_M(2,N) = PSO_M(2,N)
          PTX = PG_M(2,N) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_M(2,N) = PTX
            PL_M(2,N) = PN_M(2,N) - CPNLX
          ELSE
            PL_M(2,N) = PTX
            PN_M(2,N) = PL_M(2,N) + CPNLX
          ENDIF
!
!---    Phase condition #3  ---
!
!       SG_M > 0.0
!       SL_M > 0.0
!       SN_M > 0.0
!       Three-phase
!
!       water mass - aqueous saturation, SL_M
!       gas mass - gas saturation, SG_M
!       oil mass - gas (system) pressure, Pa gauge - PG_M,PSO_M
!       salt mass - total salt brine mass fraction - YLS_M
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.3 ) THEN
          PX = PSO_M(2,N) + PATM
!
!---      Maximum solubility of salt in pure water  ---
!
          CALL SOL_LS( T_M(2,N),XBSMX )
!
!---      Salt mass fraction in brine  ---
!
          XBSX = MIN( YLS_M(2,N),XBSMX )
!
!---      Nonaqueous flash calculation  ---
!
          CALL SP_B( T_M(2,N),XBSX,PSBX )
          CALL VPL( PSBX,PVBX,SL_M(2,N),T_M(2,N),IZN )
!
!---      Aqueous disappears  ---
!
          IF( SL_M(2,N).LT.EPSL ) SL_M(2,N) = 0.D+0
!
!---      Gas disappears  ---
!
          IF( SG_M(2,N).LT.EPSL ) THEN
            CALL BO_GOR( PX,POSM_M(2,N),T_M(2,N) )
            SG_M(2,N) = 0.D+0
          ENDIF
!
!---      Nonaqueous-liquid  ---
!
          SN_M(2,N) = 1.D+0-SL_M(2,N)-SG_M(2,N)
!
!---      Capillary pressure from saturations  ---
!
          CALL CAP_EOR_M( SG_M(2,N),SL_M(2,N),SN_M(2,N),
     &      CPGTX,CPNLX,IZN )
!
!---      Nonaqueous phases disappears, transition to 
!         phase condition #1  ---
!
          IF( (SG_M(2,N)+SN_M(2,N)).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 1
!
!---      Gas and aqueous phases disappear, transition to 
!         phase condition #4  ---
!
          ELSEIF( (SG_M(2,N)+SL_M(2,N)).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 4
!
!---      Aqueous phase disappears, transition to 
!         phase condition #5  ---
!
          ELSEIF( SL_M(2,N).LT.EPSL .AND. SG_M(2,N).GT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 5
!
!---      Gas phase disappears, transition to phase condition #2  ---
!
          ELSEIF( SG_M(2,N).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 2
!
!---      No transition from phase condition #3  ---
!
          ELSE
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 3
          ENDIF
!
!---      Phase pressure from phase capillary pressures and system
!         pressure  ---
!
          PG_M(2,N) = PSO_M(2,N)
          PTX = PG_M(2,N) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_M(2,N) = PTX
            PL_M(2,N) = PN_M(2,N) - CPNLX
          ELSE
            PL_M(2,N) = PTX
            PN_M(2,N) = PL_M(2,N) + CPNLX
          ENDIF
!
!---    Phase condition #4  ---
!
!       SG_M = 0.0
!       SL_M = 0.0
!       SN_M = 1.0
!       Nonaqueous liquid saturated
!
!       water mass - aqueous saturation, SL_M
!       gas mass - gas to oil ratio, POSM_M
!       oil mass - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!       salt mass - total salt mass - TMS_M(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.4 ) THEN
          PX = PSO_M(2,N) + PATM
!
!---      Maximum solubility of salt in pure water  ---
!
          CALL SOL_LS( T_M(2,N),XBSMX )
!
!---      Salt mass fraction in brine  ---
!
          XBSX = MIN( YLS_M(2,N),XBSMX )
!
!---      Nonaqueous flash calculation  ---
!
          CALL SP_B( T_M(2,N),XBSX,PSBX )
          CALL VPL( PSBX,PVBX,SL_M(2,N),T_M(2,N),IZN )
!
!---      Bubble-point pressure as a function of gas to oil ratio
!         (m^3/m^3) and temperature (C)  ---
!
          CALL BO_BP( PBPX,POSM_M(2,N),T_M(2,N) )
!
!---      Pressure below bubble-point pressure  ---
!
          IF( PX.LT.(PBPX-TOLX*GCPP(38,2)) ) THEN
            SG_M(2,N) = 1.D-4
            SN_M(2,N) = SN_M(2,N)*(1.D+0-SG_M(2,N))/
     &        (SN_M(2,N)+SL_M(2,N))
            SL_M(2,N) = 1.D+0-SG_M(2,N)-SN_M(2,N)
!
!---      Pressure above bubble-point pressure  ---
!
          ELSE
            SG_M(2,N) = 0.D+0
            SN_M(2,N) = 1.D+0-SL_M(2,N)
          ENDIF
!
!---      Aqueous phase appears, transition to 
!         phase condition #2 or #3  ---
!
          IF( SL_M(2,N).GT.1.D-6 ) THEN
!
!---        Gas phase appears, transition to phase condition #3  ---
!
            IF( SG_M(2,N).GT.EPSL ) THEN
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 3
!
!---        No gas phase, transition to phase condition #2  ---
!
            ELSE
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 2
            ENDIF
!
!---      Gas phase appears, transition to phase condition #5  ---
!
          ELSEIF( SG_M(2,N).GT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 5
!
!---      No transition from phase condition #4  ---
!
          ELSE
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 4
          ENDIF
!
!---      Phase pressure from phase capillary pressures and system
!         pressure  ---
!
          PG_M(2,N) = PSO_M(2,N)
          PTX = PG_M(2,N) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_M(2,N) = PTX
            PL_M(2,N) = PN_M(2,N) - CPNLX
          ELSE
            PL_M(2,N) = PTX
            PN_M(2,N) = PL_M(2,N) + CPNLX
          ENDIF
!
!---    Phase condition #5  ---
!
!       SG_M > 0.0
!       SL_M = 0.0
!       SN_M > 0.0
!       Nonaqueous liquid and gas saturated
!
!       water mass - aqueous saturation, SL_M
!       gas mass - gas saturation, SG_M
!       oil mass - gas (system) pressure, Pa gauge - PG_M,PSO_M
!       salt mass - total salt mass - TMS_M(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.5 ) THEN
          PX = PG_M(2,N) + PATM
!
!---      Maximum solubility of salt in pure water  ---
!
          CALL SOL_LS( T_M(2,N),XBSMX )
!
!---      Salt mass fraction in brine  ---
!
          XBSX = MIN( YLS_M(2,N),XBSMX )
!
!---      Nonaqueous flash calculation  ---
!
          CALL SP_B( T_M(2,N),XBSX,PSBX )
          CALL VPL( PSBX,PVBX,SL_M(2,N),T_M(2,N),IZN )
!
!---      Aqueous disappears  ---
!
          IF( SL_M(2,N).LT.EPSL ) SL_M(2,N) = 0.D+0
!
!---      Gas disappears  ---
!
          IF( SG_M(2,N).LT.EPSL ) THEN
            CALL BO_GOR( PX,POSM_M(2,N),T_M(2,N) )
            SG_M(2,N) = 0.D+0
          ENDIF
!
!---      Nonaqueous-liquid  ---
!
          SN_M(2,N) = 1.D+0-SL_M(2,N)-SG_M(2,N)
!
!---      Capillary pressure from saturations  ---
!
          CALL CAP_EOR_M( SG_M(2,N),SL_M(2,N),SN_M(2,N),
     &      CPGTX,CPNLX,IZN )
!
!---      Aqueous phase appears, transition to 
!         phase condition #2 or #3  ---
!
          IF( SL_M(2,N).GT.1.D-6 ) THEN
!
!---        Gas phase disappears, transition to phase condition #2,
!           conserving nonaqueous-liquid pressure  ---
!
            IF( SG_M(2,N).LT.EPSL ) THEN
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 2
!
!---        Gas phase remains, transition to phase condition #3,
!           conserving gas pressure  ---
!
            ELSE
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 3
            ENDIF
!
!---      Gas phase disappears, transition to phase condition #4, 
!         conserving nonaqueous-liquid pressure  ---
!
          ELSEIF( SG_M(2,N).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 4
!
!---      No transition from phase condition #5  ---
!
          ELSE
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 5
          ENDIF
!
!---      Phase pressure from phase capillary pressures and system
!         pressure  ---
!
          PG_M(2,N) = PSO_M(2,N)
          PTX = PG_M(2,N) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_M(2,N) = PTX
            PL_M(2,N) = PN_M(2,N) - CPNLX
          ELSE
            PL_M(2,N) = PTX
            PN_M(2,N) = PL_M(2,N) + CPNLX
          ENDIF
        ENDIF
  100 CONTINUE
!$NO_OMP END PARALLEL DO
!
!--- Compute increments  ---
!
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(DNR_M,IEQGC,IEQS,IEQW,ISLC,ISVC,IXP,IZ,N_DB,NFLD,
!$NO_OMP&    NGC,NPHAZ_M,PATM,PG_M,PL_M,PN_M,POSM_M,PSO_M,PVA_M,
!$NO_OMP&    PVW_M,SG_M,SL_M,SN_M,T_M,TMC_M,TMS_M,YLS_M) 
!$NO_OMP&  PRIVATE(IGC,IZN,M,N,XLSMX,XLSX)
      DO 300 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 300
        IZN = IZ(N)
        N_DB = N
!
!---    Phase condition #1  ---
!
!       SL_M = 1.0
!       SG_M+SN_M = 0.0
!       Aqueuous saturated
!
!       water mass - aqueous (system) pressure, Pa gauge - PL_M,PSO_M
!       gas mass - total gas mass - TMC_M(1,2,N)
!       oil mass - total oil mass - TMC_M(2,2,N)
!       salt mass - total salt brine mass fraction - YLS_M
!
        IF( MOD(NPHAZ_M(2,N),10).EQ.1 ) THEN
!
!---      Water mass (system pressure), increment  ---
!
          DNR_M(IEQW,N) = MAX( 1.D-1,1.D-7*(PL_M(2,N)+PATM) )
!
!---      Gas mass (total gas mass), increment  ---
!
          DNR_M(IEQGC(1),N) = MAX( 1.D-7,1.D-7*TMC_M(1,2,N) )
!
!---      Oil mass (total oil mass), increment  ---
!
          DNR_M(IEQGC(2),N) = MAX( 1.D-7,1.D-7*TMC_M(2,2,N) )
!
!---      Salt mass (total salt brine mass fraction), increment  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            XLSX = MIN( YLS_M(2,N),XLSMX )
            DNR_M(IEQS,N) = 1.D-7*XLSMX
          ENDIF
!
!---    Phase condition #2  ---
!
!       SG_M = 0.0
!       SL_M > 0.0
!       SN_M > 0.0
!       Aqueous and nonaqueous liquid saturated
!
!       water mass - aqueous saturation, SL_M
!       gas mass - gas to oil ratio, POSM_M
!       oil mass - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!       salt mass - total salt brine mass fraction - YLS_M
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.2 ) THEN
!
!---      Water mass (aqueous saturation), increment or decrement  ---
!
          DNR_M(IEQW,N) = SIGN( 1.D-6,5.D-1-SL_M(2,N) )
!
!---      Gas mass (gas to oil ratio) increment  ---
!
          DNR_M(IEQGC(1),N) = 1.D-6
!
!---      Oil mass (nonaqueous (system) pressure) increment   ---
!
          DNR_M(IEQGC(2),N) = MAX( 1.D-1,1.D-7*(PN_M(2,N)+PATM) )
!
!---      Salt mass (total salt brine mass fraction), increment  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            XLSX = MIN( YLS_M(2,N),XLSMX )
            DNR_M(IEQS,N) = 1.D-7*XLSMX
          ENDIF
!
!---    Phase condition #3  ---
!
!       SG_M > 0.0
!       SL_M > 0.0
!       SN_M > 0.0
!       Three-phase
!
!       water mass - aqueous saturation, SL_M
!       gas mass - gas saturation, SG_M
!       oil mass - gas (system) pressure, Pa gauge - PG_M,PSO_M
!       salt mass - total salt brine mass fraction - YLS_M
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.3 ) THEN
!
!---      Water mass (aqueous saturation), increment or decrement  ---
!
          DNR_M(IEQW,N) = SIGN( 1.D-6,5.D-1-SL_M(2,N) )
!
!---      Gas mass (gas (system) pressure) increment  ---
!
          DNR_M(IEQGC(1),N) = SIGN( 1.D-6,5.D-1-SG_M(2,N) )
!
!---      Oil mass (nonaqueous-liquid saturation) increment  ---
!
          DNR_M(IEQGC(2),N) = -MAX( 1.D-1,1.D-7*(PG_M(2,N)+PATM) )
!
!---      Salt mass (total salt brine mass fraction)  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            XLSX = MIN( YLS_M(2,N),XLSMX )
            DNR_M(IEQS,N) = 1.D-7*XLSMX
          ENDIF
!
!---    Phase condition #4  ---
!
!       SG_M = 0.0
!       SL_M = 0.0
!       SN_M = 1.0
!       Nonaqueous liquid saturated
!
!       water mass - aqueous saturation, SL_M
!       gas mass - gas to oil ratio, POSM_M
!       oil mass - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!       salt mass - total salt mass - TMS_M(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.4 ) THEN
!
!---      Water mass (aqueous saturation), increment or decrement  ---
!
          DNR_M(IEQW,N) = SIGN( 1.D-6,5.D-1-SL_M(2,N) )
!
!---      Gas mass (gas to oil ratio) increment  ---
!
          DNR_M(IEQGC(1),N) = 1.D-6
!
!---      Oil mass (nonaqueous (system) pressure) increment  ---
!
          DNR_M(IEQGC(2),N) = MAX( 1.D-1,1.D-7*(PN_M(2,N)+PATM) )
!
!---      Salt mass (total salt mass), increment  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            DNR_M(IEQS,N) = 1.D-7
          ENDIF
!
!---    Phase condition #5  ---
!
!       SG_M > 0.0
!       SL_M = 0.0
!       SN_M > 0.0
!       Nonaqueous liquid and gas saturated
!
!       water mass - aqueous saturation, SL_M
!       gas mass - gas saturation, SG_M
!       oil mass - gas (system) pressure, Pa gauge - PG_M,PSO_M
!       salt mass - total salt mass - TMS_M(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.5 ) THEN
!
!---      Water mass (aqueous saturation), increment or decrement  ---
!
          DNR_M(IEQW,N) = SIGN( 1.D-6,5.D-1-SL_M(2,N) )
!
!---      Gas mass (gas (system) pressure) increment  ---
!
          DNR_M(IEQGC(1),N) = SIGN( 1.D-6,5.D-1-SG_M(2,N) )
!
!---      Oil mass (nonaqueous-liquid saturation) increment  ---
!
          DNR_M(IEQGC(2),N) = -MAX( 1.D-1,1.D-7*(PG_M(2,N)+PATM) )
!
!---      Salt mass (total salt mass)  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            DNR_M(IEQS,N) = 1.D-7
          ENDIF
        ENDIF
!
!---    Increment the primary variables  ---
!
        DO 280 M = 3,ISVC+2
          T_M(M,N) = T_M(2,N)
          PG_M(M,N) = PG_M(2,N)
          PL_M(M,N) = PL_M(2,N)
          PN_M(M,N) = PN_M(2,N)
          POSM_M(M,N) = POSM_M(2,N)
          PSO_M(M,N) = PSO_M(2,N)
          PVA_M(M,N) = PVA_M(2,N)
          PVW_M(M,N) = PVW_M(2,N)
          SG_M(M,N) = SG_M(2,N)
          SL_M(M,N) = SL_M(2,N)
          SN_M(M,N) = SN_M(2,N)
          TMS_M(M,N) = TMS_M(2,N)
          YLS_M(M,N) = YLS_M(2,N)
          DO 200 IGC = 1,NGC+2
            TMC_M(IGC,M,N) = TMC_M(IGC,2,N)
  200     CONTINUE
!
!---      Phase condition #1  ---
!
!         SL_M = 1.0
!         SG_M+SN_M = 0.0
!         Aqueuous saturated
!
!         water mass - aqueous (system) pressure, Pa gauge - PL_M,PSO_M
!         gas mass - total gas mass - TMC_M(1,2,N)
!         oil mass - total oil mass - TMC_M(2,2,N)
!         salt mass - total salt brine mass fraction - YLS_M
!
          IF( MOD(NPHAZ_M(2,N),10).EQ.1 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
              TMC_M(1,M,N) = TMC_M(1,M,N) + DNR_M(IEQGC(1),N)
            ELSEIF( M.EQ.IEQGC(2)+2 ) THEN
              TMC_M(2,M,N) = TMC_M(2,M,N) + DNR_M(IEQGC(2),N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS_M(M,N) = YLS_M(M,N) + DNR_M(IEQS,N)
            ENDIF
!
!---      Phase condition #2  ---
!
!         SG_M = 0.0
!         SL_M > 0.0
!         SN_M > 0.0
!         Aqueous and nonaqueous liquid saturated
!
!         water mass - aqueous saturation, SL_M
!         gas mass - gas to oil ratio, POSM_M
!         oil mass - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         salt mass - total salt brine mass fraction - YLS_M
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.2 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              SL_M(M,N) = SL_M(M,N) + DNR_M(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
              POSM_M(M,N) = POSM_M(M,N) + DNR_M(IEQGC(1),N)
            ELSEIF( M.EQ.IEQGC(2)+2 ) THEN
              PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(2),N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS_M(M,N) = YLS_M(M,N) + DNR_M(IEQS,N)
            ENDIF
!
!---      Phase condition #3  ---
!
!         SG_M > 0.0
!         SL_M > 0.0
!         SN_M > 0.0
!         Three-phase
!
!         water mass - aqueous saturation, SL_M
!         gas mass - gas saturation, SG_M
!         oil mass - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         salt mass - total salt brine mass fraction - YLS_M
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.3 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              SL_M(M,N) = SL_M(M,N) + DNR_M(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
              SG_M(M,N) = SG_M(M,N) + DNR_M(IEQGC(1),N)
            ELSEIF( M.EQ.IEQGC(2)+2 ) THEN
              PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(2),N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS_M(M,N) = YLS_M(M,N) + DNR_M(IEQS,N)
            ENDIF
!
!---      Phase condition #4  ---
!
!         SG_M = 0.0
!         SL_M = 0.0
!         SN_M = 1.0
!         Nonaqueous liquid saturated
!
!         water mass - aqueous saturation, SL_M
!         gas mass - gas to oil ratio, POSM_M
!         oil mass - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         salt mass - total salt mass - TMS_M(2,N)
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.4 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              SL_M(M,N) = SL_M(M,N) + DNR_M(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
              POSM_M(M,N) = POSM_M(M,N) + DNR_M(IEQGC(1),N)
            ELSEIF( M.EQ.IEQGC(2)+2 ) THEN
              PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(2),N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              TMS_M(M,N) = TMS_M(M,N) + DNR_M(IEQS,N)
            ENDIF
!
!---    Phase condition #5  ---
!
!       SG_M > 0.0
!       SL_M = 0.0
!       SN_M > 0.0
!       Nonaqueous liquid and gas saturated
!
!       water mass - aqueous saturation, SL_M
!       gas mass - gas saturation, SG_M
!       oil mass - gas (system) pressure, Pa gauge - PG_M,PSO_M
!       salt mass - total salt mass - TMS_M(2,N)
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.5 ) THEN
            IF( M.EQ.IEQW+2 ) THEN
              SL_M(M,N) = SL_M(M,N) + DNR_M(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
              SG_M(M,N) = SG_M(M,N) + DNR_M(IEQGC(1),N)
            ELSEIF( M.EQ.IEQGC(2)+2 ) THEN
              PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(2),N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              TMS_M(M,N) = TMS_M(M,N) + DNR_M(IEQS,N)
            ENDIF
          ENDIF
  280   CONTINUE
  300 CONTINUE
!$NO_OMP END PARALLEL DO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_BO_M group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_EOR_M
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, July 11, 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FDVGC
      USE DUAL_POR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 FKX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_EOR_M'
!
!---  Black-oil option  ---
!
      IF( ISLC(68).EQ.0 ) THEN
        CALL INCRM_BO_M
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compositional option  ---
!
!
!---  K-factors computed from fugacity coefficients  ---
!
      IFK = 0
!
!---  Phase options, compute phase condition   ---
!
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(EPSL,GCPP,IXP,IZ,IZMC,N_DB,NFLD,NGC,NPHAZ_M,
!$NO_OMP&    PATM,PG_M,PL_M,PN_M,PORD_M,PSO_M,PVA_M,RHOL_M,RHOMG_M,
!$NO_OMP&    RHOMN_M,SG_M,SL_M,SN_M,T_M,TABS,TMC_M,VOL_M,XLA_M,XLS_M,
!$NO_OMP&    XLW_M,XMGC_M,XMGW_M,XMNC_M,YLS_M,ZMC_M) 
!$NO_OMP&  PRIVATE(BETAX,CPGTX,CPNLX,FKX,IGC,IZN,N,PSBX,PTX,PVAX,
!$NO_OMP&    PVBX,PVX,PX,RHOBX,SGNX,SGX,SLX,SNX,STMCX,STMMCX,TKX,
!$NO_OMP&    TMLX,XBSMX,XBSX,XLAX,XLSX,XLWX,ZGX,ZMCMX,ZNX)
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 100
        IZN = IZ(N)
        N_DB = N
!
!---    Absolute temperature  ---
!
        TKX = T_M(2,N)+TABS
!
!---    Phase condition #1  ---
!
!       SL_M = 1.0, SG_M+SN_M = 0.0
!
!       energy - temperature, C - T_M
!       water mass - aqueous (system) pressure, Pa gauge - PL_M,PSO_M
!       CO2 mass - CO2 vapor pressure - PVA_M
!       CH4 mass - total CH4 mass - TMC_M
!       petroleum component mass - total petroleum comp. mass - TMC_M
!       salt mass - total salt brine mass fraction - YLS_M
!
        IF( MOD(NPHAZ_M(2,N),10).EQ.1 ) THEN
!
!---      System pressure  ---
!
          PX = PSO_M(2,N) + PATM
!
!---      Maximum solubility of salt in pure water  ---
!
          CALL SOL_LS( T_M(2,N),XBSMX )
!
!---      Salt mass fraction in brine  ---
!
          XBSX = MIN( YLS_M(2,N),XBSMX )
          CALL SP_B( T_M(2,N),XBSX,PSBX )
!
!---      Maximum solubility of CO2 in aqueous phase  ---
!
          PVAX = PX - PSBX
          CALL EQUIL( T_M(2,N),PX,PVAX,XBSX,XLAX,XLSX,XLWX )
!
!---      Solubility of CO2 in aqueous phase  ---
!
          CALL EQUIL( T_M(2,N),PX,PVA_M(2,N),XBSX,XLA_M(2,N),XLS_M(2,N),
     &      XLW_M(2,N) )
!
!---      Brine density  ---
!
          CALL DENS_B( T_M(2,N),PX,XBSX,RHOBX )
!
!---      Aqueous density  ---
!
          CALL DENS_L( T_M(2,N),RHOBX,XLA_M(2,N),RHOL_M(2,N) )
!
!---      Total mass of CO2  ---
!
          TMC_M(1,2,N) = XLA_M(2,N)*RHOL_M(2,N)*PORD_M(2,N)*VOL_M(N)
!
!---      Non-CO2 petroleum mass and total petroleum molar mass  ---
!
          STMCX = -TMC_M(1,2,N)
          STMMCX = 0.D+0
          DO 10 IGC = 1,NGC+2
            STMCX = STMCX + TMC_M(IGC,2,N)
            STMMCX = STMMCX + TMC_M(IGC,2,N)*GCPP(1,IGC)
   10     CONTINUE
          IF( STMMCX.GT.EPSL ) THEN
            ZMCMX = 0.D+0
            DO 12 IGC = 1,NGC+2
              ZMC_M(IGC,2,N) = TMC_M(IGC,2,N)*GCPP(1,IGC)/STMMCX
              IF( ZMC_M(IGC,2,N).GT.ZMCMX ) THEN
                IZMC(N) = IGC
                ZMCMX = ZMC_M(IGC,2,N)
              ENDIF
   12       CONTINUE
          ELSE
            ZMC_M(1,2,N) = 1.D+0
            IZMC(N) = 1
            DO 14 IGC = 2,NGC+2
              ZMC_M(IGC,2,N) = 0.D+0
   14       CONTINUE
          ENDIF
!
!---      Aqueous mass  ---
!
          TMLX = RHOL_M(2,N)*PORD_M(2,N)*VOL_M(N)
!
!---      CO2 concentration in aqueous phase exceeds solubility limit
!         or total petroleum component mass execeeds 1.e-5 of the
!         aqueous phase mass, transition to phase condition #2 or #3 ---
!
          IF( XLA_M(2,N).GT.XLAX .OR. STMCX.GT.(1.D-5*TMLX) ) THEN
            CALL SP_B( T_M(2,N),XBSX,PSBX )
            CALL VPL( PSBX,PVBX,SLX,T_M(2,N),IZN )
            PVX = PX - PVBX
            BETAX = -1.D+0
            CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVBX,RHOMG_M(2,N),
     &        RHOMN_M(2,N),SGX,SNX,T_M(2,N),XMGC_M(1,2,N),XMGW_M(2,N),
     &        XMNC_M(1,2,N),ZGX,ZMC_M(1,2,N),ZNX,NPHAZ_M(2,N) )
            SGNX = 1.D-4
            SL_M(2,N) = 1.D+0-SGNX
            SG_M(2,N) = (1.D+0-SL_M(2,N))*SGX
            SN_M(2,N) = (1.D+0-SL_M(2,N))*SNX
            CALL CAP_EOR_M( SG_M(2,N),SL_M(2,N),SN_M(2,N),
     &        CPGTX,CPNLX,IZN )
            IF( SG_M(2,N).LT.EPSL ) THEN
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 2
            ELSE
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 3
            ENDIF
          ELSE
            CPGTX = 0.D+0
            CPNLX = 0.D+0
          ENDIF
!
!---      Phase pressure from phase capillary pressures and system
!         pressure  ---
!
          PG_M(2,N) = PSO_M(2,N)
          PTX = PG_M(2,N) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_M(2,N) = PTX
            PL_M(2,N) = PN_M(2,N) - CPNLX
          ELSE
            PL_M(2,N) = PTX
            PN_M(2,N) = PL_M(2,N) + CPNLX
          ENDIF
!
!---    Phase condition #2  ---
!
!       SL_M < 1.0, ST > 0.0, SG_M = 0.0
!
!       energy - temperature, C - T_M
!       water mass - aqueous saturation, SL_M
!       CO2 mass
!         if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!         endif
!       salt mass - total salt brine mass fraction - YLS_M
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.2 ) THEN
          PX = PSO_M(2,N) + PATM
!
!---      Maximum solubility of salt in pure water  ---
!
          CALL SOL_LS( T_M(2,N),XBSMX )
!
!---      Salt mass fraction in brine  ---
!
          XBSX = MIN( YLS_M(2,N),XBSMX )
!
!---      Nonaqueous flash calculation  ---
!
          CALL SP_B( T_M(2,N),XBSX,PSBX )
          CALL VPL( PSBX,PVBX,SL_M(2,N),T_M(2,N),IZN )
          PVX = PX - PVBX
          BETAX = -1.D+0
          CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVBX,RHOMG_M(2,N),
     &      RHOMN_M(2,N),SGX,SNX,T_M(2,N),XMGC_M(1,2,N),XMGW_M(2,N),
     &      XMNC_M(1,2,N),ZGX,ZMC_M(1,2,N),ZNX,NPHAZ_M(2,N) )
          SG_M(2,N) = (1.D+0-SL_M(2,N))*SGX
          SN_M(2,N) = (1.D+0-SL_M(2,N))*SNX
          CALL CAP_EOR_M( SG_M(2,N),SL_M(2,N),SN_M(2,N),
     &      CPGTX,CPNLX,IZN )
!
!---      Assign maximum concentration of petroleum component  ---
!
          ZMCMX = 0.D+0
          DO 20 IGC = 1,NGC+2
            IF( ZMC_M(IGC,2,N).GT.ZMCMX ) THEN
              IZMC(N) = IGC
              ZMCMX = ZMC_M(IGC,2,N)
            ENDIF
   20     CONTINUE
!
!---      Nonaqueous phases disappears, transition to 
!         phase condition #1, conserving aqueous pressure  ---
!
          IF( (SG_M(2,N)+SN_M(2,N)).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 1
!
!---      Aqueous phases disappear, transition to 
!         phase condition #4  ---
!
          ELSEIF( (SG_M(2,N)+SL_M(2,N)).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 4
!
!---      Aqueous phase disappears and gas phase appears, transition to 
!         phase condition #5  ---
!
          ELSEIF( SL_M(2,N).LT.EPSL .AND. SG_M(2,N).GT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 5
!
!---      Gas phase appears, transition to phase condition #3  ---
!
          ELSEIF( SG_M(2,N).GT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 3
!
!---      No transition from phase condition #2  ---
!
          ELSE
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 2
          ENDIF
!
!---      Phase pressure from phase capillary pressures and system
!         pressure  ---
!
          PG_M(2,N) = PSO_M(2,N)
          PTX = PG_M(2,N) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_M(2,N) = PTX
            PL_M(2,N) = PN_M(2,N) - CPNLX
          ELSE
            PL_M(2,N) = PTX
            PN_M(2,N) = PL_M(2,N) + CPNLX
          ENDIF
!
!---    Phase condition #3  ---
!
!       SL_M < 1.0, ST >= 0.0, SG_M > 0.0
!
!       energy - temperature, C - T_M
!       water mass - aqueous saturation, SL_M
!       CO2 mass
!         if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!         endif
!       salt mass - total salt brine mass fraction - YLS_M
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.3 ) THEN
          PX = PSO_M(2,N) + PATM
!
!---      Maximum solubility of salt in pure water  ---
!
          CALL SOL_LS( T_M(2,N),XBSMX )
!
!---      Salt mass fraction in brine  ---
!
          XBSX = MIN( YLS_M(2,N),XBSMX )
!
!---      Nonaqueous flash calculation  ---
!
          CALL SP_B( T_M(2,N),XBSX,PSBX )
          CALL VPL( PSBX,PVBX,SL_M(2,N),T_M(2,N),IZN )
          PVX = PX - PVBX
          BETAX = -1.D+0
          CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVBX,RHOMG_M(2,N),
     &      RHOMN_M(2,N),SGX,SNX,T_M(2,N),XMGC_M(1,2,N),XMGW_M(2,N),
     &      XMNC_M(1,2,N),ZGX,ZMC_M(1,2,N),ZNX,NPHAZ_M(2,N) )
          SG_M(2,N) = (1.D+0-SL_M(2,N))*SGX
          SN_M(2,N) = (1.D+0-SL_M(2,N))*SNX
          CALL CAP_EOR_M( SG_M(2,N),SL_M(2,N),SN_M(2,N),
     &      CPGTX,CPNLX,IZN )
!
!---      Assign maximum concentration of petroleum component  ---
!
          ZMCMX = 0.D+0
          DO 30 IGC = 1,NGC+2
            IF( ZMC_M(IGC,2,N).GT.ZMCMX ) THEN
              IZMC(N) = IGC
              ZMCMX = ZMC_M(IGC,2,N)
            ENDIF
   30     CONTINUE
!
!---      Nonaqueous phases disappears, transition to 
!         phase condition #1  ---
!
          IF( (SG_M(2,N)+SN_M(2,N)).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 1
!
!---      Gas and aqueous phases disappear, transition to 
!         phase condition #4  ---
!
          ELSEIF( (SG_M(2,N)+SL_M(2,N)).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 4
!
!---      Aqueous phase disappears, transition to 
!         phase condition #5  ---
!
          ELSEIF( SL_M(2,N).LT.EPSL .AND. SG_M(2,N).GT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 5
!
!---      Gas phase disappears, transition to phase condition #2  ---
!
          ELSEIF( SG_M(2,N).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 2
!
!---      No transition from phase condition #3  ---
!
          ELSE
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 3
          ENDIF
!
!---      Phase pressure from phase capillary pressures and system
!         pressure  ---
!
          PG_M(2,N) = PSO_M(2,N)
          PTX = PG_M(2,N) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_M(2,N) = PTX
            PL_M(2,N) = PN_M(2,N) - CPNLX
          ELSE
            PL_M(2,N) = PTX
            PN_M(2,N) = PL_M(2,N) + CPNLX
          ENDIF
!
!---    Phase condition #4  ---
!
!       SL_M = 0.0, ST = 1.0, SG_M = 0.0
!
!       energy - temperature, C - T_M
!       water mass - aqueous saturation, SL_M
!       CO2 mass
!         if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!         endif
!       salt mass - total salt mass - TMS_M(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.4 ) THEN
          PX = PSO_M(2,N) + PATM
!
!---      Maximum solubility of salt in pure water  ---
!
          CALL SOL_LS( T_M(2,N),XBSMX )
!
!---      Salt mass fraction in brine  ---
!
          XBSX = MIN( YLS_M(2,N),XBSMX )
!
!---      Nonaqueous flash calculation  ---
!
          CALL SP_B( T_M(2,N),XBSX,PSBX )
          CALL VPL( PSBX,PVBX,SL_M(2,N),T_M(2,N),IZN )
          PVX = PX - PVBX
          BETAX = -1.D+0
          CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVBX,RHOMG_M(2,N),
     &      RHOMN_M(2,N),SGX,SNX,T_M(2,N),XMGC_M(1,2,N),XMGW_M(2,N),
     &      XMNC_M(1,2,N),ZGX,ZMC_M(1,2,N),ZNX,NPHAZ_M(2,N) )
          SG_M(2,N) = (1.D+0-SL_M(2,N))*SGX
          SN_M(2,N) = (1.D+0-SL_M(2,N))*SNX
          CALL CAP_EOR_M( SG_M(2,N),SL_M(2,N),SN_M(2,N),
     &      CPGTX,CPNLX,IZN )
          ZMCMX = 0.D+0
          DO 40 IGC = 1,NGC+2
            IF( ZMC_M(IGC,2,N).GT.ZMCMX ) THEN
              IZMC(N) = IGC
              ZMCMX = ZMC_M(IGC,2,N)
            ENDIF
   40     CONTINUE
!
!---      Aqueous phase appears, transition to 
!         phase condition #2 or #3  ---
!
          IF( SL_M(2,N).GT.1.D-6 ) THEN
!
!---        Gas phase appears, transition to phase condition #3  ---
!
            IF( SG_M(2,N).GT.EPSL ) THEN
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 3
!
!---        No gas phase, transition to phase condition #2  ---
!
            ELSE
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 2
            ENDIF
!
!---      Gas phase appears, transition to phase condition #5  ---
!
          ELSEIF( SG_M(2,N).GT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 5
!
!---      No transition from phase condition #4  ---
!
          ELSE
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 4
          ENDIF
!
!---      Phase pressure from phase capillary pressures and system
!         pressure  ---
!
          PG_M(2,N) = PSO_M(2,N)
          PTX = PG_M(2,N) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_M(2,N) = PTX
            PL_M(2,N) = PN_M(2,N) - CPNLX
          ELSE
            PL_M(2,N) = PTX
            PN_M(2,N) = PL_M(2,N) + CPNLX
          ENDIF
!
!---    Phase condition #5  ---
!
!       SL_M = 0.0, ST >= 0.0, SG_M > 0.0
!
!       energy - temperature, C - T_M
!       water mass - aqueous saturation, SL_M
!       CO2 mass
!         if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!         endif
!       salt mass - total salt mass - TMS_M(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.5 ) THEN
          PX = PG_M(2,N) + PATM
!
!---      Maximum solubility of salt in pure water  ---
!
          CALL SOL_LS( T_M(2,N),XBSMX )
!
!---      Salt mass fraction in brine  ---
!
          XBSX = MIN( YLS_M(2,N),XBSMX )
!
!---      Nonaqueous flash calculation  ---
!
          CALL SP_B( T_M(2,N),XBSX,PSBX )
          CALL VPL( PSBX,PVBX,SL_M(2,N),T_M(2,N),IZN )
          PVX = PX - PVBX
          BETAX = -1.D+0
          CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVBX,RHOMG_M(2,N),
     &      RHOMN_M(2,N),SGX,SNX,T_M(2,N),XMGC_M(1,2,N),XMGW_M(2,N),
     &      XMNC_M(1,2,N),ZGX,ZMC_M(1,2,N),ZNX,NPHAZ_M(2,N) )
          SG_M(2,N) = (1.D+0-SL_M(2,N))*SGX
          SN_M(2,N) = (1.D+0-SL_M(2,N))*SNX
          CALL CAP_EOR_M( SG_M(2,N),SL_M(2,N),SN_M(2,N),
     &      CPGTX,CPNLX,IZN )
          ZMCMX = 0.D+0
          DO 50 IGC = 1,NGC+2
            IF( ZMC_M(IGC,2,N).GT.ZMCMX ) THEN
              IZMC(N) = IGC
              ZMCMX = ZMC_M(IGC,2,N)
            ENDIF
   50     CONTINUE
!
!---      Aqueous phase appears, transition to 
!         phase condition #2 or #3  ---
!
          IF( SL_M(2,N).GT.1.D-6 ) THEN
!
!---        Gas phase disappears, transition to phase condition #2,
!           conserving nonaqueous-liquid pressure  ---
!
            IF( SG_M(2,N).LT.EPSL ) THEN
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 2
!
!---        Gas phase remains, transition to phase condition #3,
!           conserving gas pressure  ---
!
            ELSE
              NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 3
            ENDIF
!
!---      Gas phase disappears, transition to phase condition #4, 
!         conserving nonaqueous-liquid pressure  ---
!
          ELSEIF( SG_M(2,N).LT.EPSL ) THEN
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 4
!
!---      No transition from phase condition #5  ---
!
          ELSE
            NPHAZ_M(2,N) = 10*(NPHAZ_M(2,N)/10) + 5
          ENDIF
!
!---      Phase pressure from phase capillary pressures and system
!         pressure  ---
!
          PG_M(2,N) = PSO_M(2,N)
          PTX = PG_M(2,N) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_M(2,N) = PTX
            PL_M(2,N) = PN_M(2,N) - CPNLX
          ELSE
            PL_M(2,N) = PTX
            PN_M(2,N) = PL_M(2,N) + CPNLX
          ENDIF
        ENDIF
  100 CONTINUE
!$NO_OMP END PARALLEL DO
!
!--- Compute increments  ---
!
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(DNR_M,GCPP,IEQGC,IEQS,IEQT,IEQW,ISLC,ISVC,IXP,IZ,
!$NO_OMP&    IZMC,N_DB,NFLD,NGC,NPHAZ_M,PATM,PG_M,PL_M,PN_M,PSO_M,
!$NO_OMP&    PVA_M,PVW_M,SL_M,T_M,TMC_M,TMS_M,YLS_M,ZMC_M,TABS) 
!$NO_OMP&  PRIVATE(IGC,IZN,M,N,TKX,XLSMX,XLSX)
      DO 300 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 300
        IZN = IZ(N)
        N_DB = N
!
!---    Absolute temperature  ---
!
        TKX = T_M(2,N)+TABS
!
!---    Phase condition #1  ---
!
!       SL_M = 1.0, SG_M+SN_M = 0.0
!
!       energy - temperature, C - T_M
!       water mass - aqueous (system) pressure, Pa gauge - PL_M,PSO_M
!       CO2 mass - CO2 vapor pressure - PVA_M
!       CH4 mass - total CH4 mass - TMC_M
!       petroleum component mass - total petroleum comp. mass - TMC_M
!       salt mass - total salt brine mass fraction - YLS_M
!
        IF( MOD(NPHAZ_M(2,N),10).EQ.1 ) THEN
!
!---      Energy (temperature), increment  ---
!
          IF( ISLC(30).EQ.0 ) DNR_M(IEQT,N) = -1.D-7
!
!---      Water mass (system pressure), increment  ---
!
          DNR_M(IEQW,N) = MAX( 1.D-1,1.D-7*(PL_M(2,N)+PATM) )
!
!---      CO2 mass (CO2 vapor pressure), increment  ---
!
          IF( ISLC(45).EQ.0 ) DNR_M(IEQGC(1),N) = 1.D-2
!
!---      CH4 mass (total CH4 mass), increment  ---
!
          IF( ISLC(64).EQ.0 ) 
     &      DNR_M(IEQGC(2),N) = MAX( 1.D-7,1.D-7*TMC_M(2,2,N) )
!
!---      Petroleum component mass (total petroleum component mass),
!         increment  ---
!
          DO 112 IGC = 3,NGC+2
            DNR_M(IEQGC(IGC),N) = MAX( 1.D-7,1.D-7*TMC_M(IGC,2,N) )
  112     CONTINUE
!
!---      Salt mass (total salt brine mass fraction), increment  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            XLSX = MIN( YLS_M(2,N),XLSMX )
            DNR_M(IEQS,N) = 1.D-7*XLSMX
          ENDIF
!
!---    Phase condition #2  ---
!
!       SL_M < 1.0, ST > 0.0, SG_M = 0.0
!
!       energy - temperature, C - T_M
!       water mass - aqueous saturation, SL_M
!       CO2 mass
!         if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!         endif
!       salt mass - total salt brine mass fraction - YLS_M
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.2 ) THEN
!
!---      Energy (temperature), increment  ---
!
          IF( ISLC(30).EQ.0 ) DNR_M(IEQT,N) = -1.D-7
!
!---      Water mass (aqueous saturation), increment or decrement  ---
!
          DNR_M(IEQW,N) = SIGN(1.D-6,5.D-1-SL_M(2,N))
!
!---      CO2 mass increment  ---
!
          IF( ISLC(45).EQ.0 ) THEN
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.1 ) THEN
              DNR_M(IEQGC(1),N) = MAX( 1.D-1,1.D-7*(PN_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Decrease mole fraction of volatile components,
!             and increase mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,1) ) THEN
                DNR_M(IEQGC(1),N) = MAX( 1.D-6*ZMC_M(1,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(1),N) = MAX( 1.D-6*ZMC_M(1,2,N),1.D-9 )
              ENDIF
            ENDIF
          ENDIF
!
!---      CH4 mass increment  ---
!
          IF( ISLC(64).EQ.0 ) THEN
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.2 ) THEN
              DNR_M(IEQGC(2),N) = MAX( 1.D-1,1.D-7*(PN_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Decrease mole fraction of volatile components,
!             and increase mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,2) ) THEN
                DNR_M(IEQGC(2),N) = MAX( 1.D-6*ZMC_M(2,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(2),N) = MAX( 1.D-6*ZMC_M(2,2,N),1.D-9 )
              ENDIF
            ENDIF
          ENDIF
!
!---      Petroleum component mass (total petroleum component mass),
!         increment  ---
!
          DO 122 IGC = 3,NGC+2
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.IGC ) THEN
              DNR_M(IEQGC(IGC),N) = MAX( 1.D-1,1.D-7*(PN_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Decrease mole fraction of volatile components,
!             and increase mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,IGC) ) THEN
                DNR_M(IEQGC(IGC),N) = MAX( 1.D-6*ZMC_M(IGC,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(IGC),N) = MAX( 1.D-6*ZMC_M(IGC,2,N),1.D-9 )
              ENDIF
            ENDIF
  122     CONTINUE
!
!---      Salt mass (total salt brine mass fraction), increment  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            XLSX = MIN( YLS_M(2,N),XLSMX )
            DNR_M(IEQS,N) = 1.D-7*XLSMX
          ENDIF
!
!---    Phase condition #3  ---
!
!       SL_M < 1.0, ST > 0.0, SG_M > 0.0
!
!       energy - temperature, C - T_M
!       water mass - aqueous saturation, SL_M
!       CO2 mass
!         if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!         endif
!       salt mass - total salt brine mass fraction - YLS_M
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.3 ) THEN
!
!---      Energy (temperature), increment  ---
!
          IF( ISLC(30).EQ.0 ) DNR_M(IEQT,N) = 1.D-7
!
!---      Water mass (aqueous saturation), increment or decrement  ---
!
          DNR_M(IEQW,N) = SIGN(1.D-6,5.D-1-SL_M(2,N))
!
!---      CO2 mass increment  ---
!
          IF( ISLC(45).EQ.0 ) THEN
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.1 ) THEN
              DNR_M(IEQGC(1),N) = -MAX( 1.D-1,1.D-7*(PG_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Increase mole fraction of volatile components,
!             and decrease mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,1) ) THEN
                DNR_M(IEQGC(1),N) = MAX( 1.D-6*ZMC_M(1,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(1),N) = MAX( 1.D-6*ZMC_M(1,2,N),1.D-9 )
              ENDIF
            ENDIF
          ENDIF
!
!---      CH4 mass increment  ---
!
          IF( ISLC(64).EQ.0 ) THEN
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.2 ) THEN
              DNR_M(IEQGC(2),N) = -MAX( 1.D-1,1.D-7*(PG_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Increase mole fraction of volatile components,
!             and decrease mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,2) ) THEN
                DNR_M(IEQGC(2),N) = MAX( 1.D-6*ZMC_M(2,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(2),N) = MAX( 1.D-6*ZMC_M(2,2,N),1.D-9 )
              ENDIF
            ENDIF
          ENDIF
!
!---      Petroleum component mass (total petroleum component mass),
!         increment  ---
!
          DO 132 IGC = 3,NGC+2
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.IGC ) THEN
              DNR_M(IEQGC(IGC),N) = -MAX( 1.D-1,1.D-7*(PG_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Increase mole fraction of volatile components,
!             and decrease mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,IGC) ) THEN
                DNR_M(IEQGC(IGC),N) = MAX( 1.D-6*ZMC_M(IGC,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(IGC),N) = MAX( 1.D-6*ZMC_M(IGC,2,N),1.D-9 )
              ENDIF
            ENDIF
  132     CONTINUE
!
!---      Salt mass (total salt brine mass fraction)  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            XLSX = MIN( YLS_M(2,N),XLSMX )
            DNR_M(IEQS,N) = 1.D-7*XLSMX
          ENDIF
!
!---    Phase condition #4  ---
!
!       SL_M = 0.0, ST = 1.0, SG_M = 0.0
!
!       energy - temperature, C - T_M
!       water mass - aqueous saturation, SL_M
!       CO2 mass
!         if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!           - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!         endif
!       salt mass - total salt mass - TMS_M(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.4 ) THEN
!
!---      Energy (temperature), increment  ---
!
          IF( ISLC(30).EQ.0 ) DNR_M(IEQT,N) = -1.D-7
!
!---      Water mass (aqueous saturation), increment or decrement  ---
!
          DNR_M(IEQW,N) = SIGN(1.D-6,5.D-1-SL_M(2,N))
!
!---      CO2 mass increment  ---
!
          IF( ISLC(45).EQ.0 ) THEN
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.1 ) THEN
              DNR_M(IEQGC(1),N) = MAX( 1.D-1,1.D-7*(PN_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Decrease mole fraction of volatile components,
!             and increase mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,1) ) THEN
                DNR_M(IEQGC(1),N) = MAX( 1.D-6*ZMC_M(1,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(1),N) = MAX( 1.D-6*ZMC_M(1,2,N),1.D-9 )
              ENDIF
            ENDIF
          ENDIF
!
!---      CH4 mass increment  ---
!
          IF( ISLC(64).EQ.0 ) THEN
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.2 ) THEN
              DNR_M(IEQGC(2),N) = MAX( 1.D-1,1.D-7*(PN_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Decrease mole fraction of volatile components,
!             and increase mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,2) ) THEN
                DNR_M(IEQGC(2),N) = MAX( 1.D-6*ZMC_M(2,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(2),N) = MAX( 1.D-6*ZMC_M(2,2,N),1.D-9 )
              ENDIF
            ENDIF
          ENDIF
!
!---      Petroleum component mass (total petroleum component mass),
!         increment  ---
!
          DO 142 IGC = 3,NGC+2
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.IGC ) THEN
              DNR_M(IEQGC(IGC),N) = MAX( 1.D-1,1.D-7*(PN_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Decrease mole fraction of volatile components,
!             and increase mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,IGC) ) THEN
                DNR_M(IEQGC(IGC),N) = MAX( 1.D-6*ZMC_M(IGC,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(IGC),N) = MAX( 1.D-6*ZMC_M(IGC,2,N),1.D-9 )
              ENDIF
            ENDIF
  142     CONTINUE
!
!---      Salt mass (total salt mass), increment  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            DNR_M(IEQS,N) = 1.D-7
          ENDIF
!
!---    Phase condition #5  ---
!
!       SL_M = 0.0, ST >= 0.0, SG_M > 0.0
!
!       energy - temperature, C - T_M
!       water mass - aqueous saturation, SL_M
!       CO2 mass
!         if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!           - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!         endif
!       salt mass - total salt mass - TMS_M(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.5 ) THEN
!
!---      Energy (temperature), increment  ---
!
          IF( ISLC(30).EQ.0 ) DNR_M(IEQT,N) = -1.D-7
!
!---      Water mass (aqueous saturation), increment or decrement  ---
!
          DNR_M(IEQW,N) = SIGN(1.D-6,5.D-1-SL_M(2,N))
!
!---      CO2 mass increment  ---
!
          IF( ISLC(45).EQ.0 ) THEN
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.1 ) THEN
              DNR_M(IEQGC(1),N) = -MAX( 1.D-1,1.D-7*(PG_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Increase mole fraction of volatile components,
!             and decrease mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,1) ) THEN
                DNR_M(IEQGC(1),N) = MAX( 1.D-6*ZMC_M(1,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(1),N) = MAX( 1.D-6*ZMC_M(1,2,N),1.D-9 )
              ENDIF
            ENDIF
          ENDIF
!
!---      CH4 mass increment  ---
!
          IF( ISLC(64).EQ.0 ) THEN
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.2 ) THEN
              DNR_M(IEQGC(2),N) = -MAX( 1.D-1,1.D-7*(PG_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
 !
!---          Increase mole fraction of volatile components,
!             and decrease mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,2) ) THEN
                DNR_M(IEQGC(2),N) = MAX( 1.D-6*ZMC_M(2,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(2),N) = MAX( 1.D-6*ZMC_M(2,2,N),1.D-9 )
              ENDIF
            ENDIF
          ENDIF
!
!---      Petroleum component mass (total petroleum component mass),
!         increment  ---
!
          DO 152 IGC = 3,NGC+2
!
!---        System pressure, increment  ---
!
            IF( IZMC(N).EQ.IGC ) THEN
              DNR_M(IEQGC(IGC),N) = -MAX( 1.D-1,1.D-7*(PG_M(2,N)+PATM) )
!
!---        Nonaqueous mole fraction, increment  ---
!
            ELSE
!
!---          Increase mole fraction of volatile components,
!             and decrease mole fraction of non-volatile components  ---
!
              IF( TKX.GT.GCPP(3,IGC) ) THEN
                DNR_M(IEQGC(IGC),N) = MAX( 1.D-6*ZMC_M(IGC,2,N),1.D-9 )
              ELSE
                DNR_M(IEQGC(IGC),N) = MAX( 1.D-6*ZMC_M(IGC,2,N),1.D-9 )
              ENDIF
            ENDIF
  152     CONTINUE
!
!---      Salt mass (total salt mass)  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            DNR_M(IEQS,N) = 1.D-7
          ENDIF
        ENDIF
!
!---    Increment the primary variables  ---
!
        DO 280 M = 3,ISVC+2
          T_M(M,N) = T_M(2,N)
          PG_M(M,N) = PG_M(2,N)
          PL_M(M,N) = PL_M(2,N)
          PN_M(M,N) = PN_M(2,N)
          PSO_M(M,N) = PSO_M(2,N)
          PVA_M(M,N) = PVA_M(2,N)
          PVW_M(M,N) = PVW_M(2,N)
          SL_M(M,N) = SL_M(2,N)
          TMS_M(M,N) = TMS_M(2,N)
          YLS_M(M,N) = YLS_M(2,N)
          DO 200 IGC = 1,NGC+2
            TMC_M(IGC,M,N) = TMC_M(IGC,2,N)
            ZMC_M(IGC,M,N) = ZMC_M(IGC,2,N)
  200     CONTINUE
!
!---      Phase condition #1  ---
!
!         SL_M = 1.0, SG_M+SN_M = 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous (system) pressure, Pa gauge - PL_M,PSO_M
!         CO2 mass - CO2 vapor pressure - PVA_M
!         CH4 mass - total CH4 mass - TMC_M
!         petroleum component mass - total petroleum comp. mass - TMC_M
!         salt mass - total salt brine mass fraction - YLS_M
!
          IF( MOD(NPHAZ_M(2,N),10).EQ.1 ) THEN
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T_M(M,N) = T_M(M,N) + DNR_M(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 .AND. ISLC(45).EQ.0 ) THEN
              PVA_M(M,N) = PVA_M(M,N) + DNR_M(IEQGC(1),N)
            ELSEIF( M.EQ.IEQGC(2)+2 .AND. ISLC(64).EQ.0 ) THEN
              TMC_M(2,M,N) = TMC_M(2,M,N) + DNR_M(IEQGC(2),N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS_M(M,N) = YLS_M(M,N) + DNR_M(IEQS,N)
            ELSE
              DO 210 IGC = 3,NGC+2
                IF( M.EQ.IEQGC(IGC)+2 ) THEN
                  TMC_M(IGC,M,N) = TMC_M(IGC,M,N) + DNR_M(IEQGC(IGC),N)
                ENDIF
  210         CONTINUE
            ENDIF
!
!---      Phase condition #2  ---
!
!         SL_M < 1.0, ST > 0.0, SG_M = 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous saturation, SL_M
!         CO2 mass
!           if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!           endif
!         CH4 mass
!           if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!           endif
!         petroleum component mass
!           if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!           endif
!         salt mass - total salt brine mass fraction - YLS_M
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.2 ) THEN
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T_M(M,N) = T_M(M,N) + DNR_M(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              SL_M(M,N) = SL_M(M,N) + DNR_M(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 .AND. ISLC(45).EQ.0 ) THEN
              IF( IZMC(N).EQ.1 ) THEN
                PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(1),N)
              ELSE
                ZMC_M(1,M,N) = ZMC_M(1,M,N) + DNR_M(IEQGC(1),N)
              ENDIF
            ELSEIF( M.EQ.IEQGC(2)+2 .AND. ISLC(64).EQ.0 ) THEN
              IF( IZMC(N).EQ.2 ) THEN
                PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(2),N)
              ELSE
                ZMC_M(2,M,N) = ZMC_M(2,M,N) + DNR_M(IEQGC(2),N)
              ENDIF
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS_M(M,N) = YLS_M(M,N) + DNR_M(IEQS,N)
            ELSE
              DO 220 IGC = 3,NGC+2
                IF( M.EQ.IEQGC(IGC)+2 ) THEN
                  IF( IZMC(N).EQ.IGC ) THEN
                    PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(IGC),N)
                  ELSE
                    ZMC_M(IGC,M,N) = ZMC_M(IGC,M,N)+DNR_M(IEQGC(IGC),N)
                  ENDIF
                ENDIF
  220         CONTINUE
            ENDIF
!
!---      Phase condition #3  ---
!
!         SL_M < 1.0, ST > 0.0, SG_M > 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous saturation, SL_M
!         CO2 mass
!           if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!           endif
!         CH4 mass
!           if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!           endif
!         petroleum component mass
!           if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!           endif
!         salt mass - total salt brine mass fraction - YLS_M
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.3 ) THEN
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T_M(M,N) = T_M(M,N) + DNR_M(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              SL_M(M,N) = SL_M(M,N) + DNR_M(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 .AND. ISLC(45).EQ.0 ) THEN
              IF( IZMC(N).EQ.1 ) THEN
                PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(1),N)
              ELSE
                ZMC_M(1,M,N) = ZMC_M(1,M,N) + DNR_M(IEQGC(1),N)
              ENDIF
            ELSEIF( M.EQ.IEQGC(2)+2 .AND. ISLC(64).EQ.0 ) THEN
              IF( IZMC(N).EQ.2 ) THEN
                PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(2),N)
              ELSE
                ZMC_M(2,M,N) = ZMC_M(2,M,N) + DNR_M(IEQGC(2),N)
              ENDIF
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS_M(M,N) = YLS_M(M,N) + DNR_M(IEQS,N)
            ELSE
              DO 230 IGC = 3,NGC+2
                IF( M.EQ.IEQGC(IGC)+2 ) THEN
                  IF( IZMC(N).EQ.IGC ) THEN
                    PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(IGC),N)
                  ELSE
                    ZMC_M(IGC,M,N) = ZMC_M(IGC,M,N)+DNR_M(IEQGC(IGC),N)
                  ENDIF
                ENDIF
  230         CONTINUE
            ENDIF
!
!---      Phase condition #4  ---
!
!         SL_M = 0.0, ST = 1.0, SG_M = 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous saturation, SL_M
!         CO2 mass
!           if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!           endif
!         CH4 mass
!           if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!           endif
!         petroleum component mass
!           if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!           endif
!         salt mass - total salt mass - TMS_M(2,N)
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.4 ) THEN
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T_M(M,N) = T_M(M,N) + DNR_M(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              SL_M(M,N) = SL_M(M,N) + DNR_M(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 .AND. ISLC(45).EQ.0 ) THEN
              IF( IZMC(N).EQ.1 ) THEN
                PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(1),N)
              ELSE
                ZMC_M(1,M,N) = ZMC_M(1,M,N) + DNR_M(IEQGC(1),N)
              ENDIF
            ELSEIF( M.EQ.IEQGC(2)+2 .AND. ISLC(64).EQ.0 ) THEN
              IF( IZMC(N).EQ.2 ) THEN
                PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(2),N)
              ELSE
                ZMC_M(2,M,N) = ZMC_M(2,M,N) + DNR_M(IEQGC(2),N)
              ENDIF
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              TMS_M(M,N) = TMS_M(M,N) + DNR_M(IEQS,N)
            ELSE
              DO 240 IGC = 3,NGC+2
                IF( M.EQ.IEQGC(IGC)+2 ) THEN
                  IF( IZMC(N).EQ.IGC ) THEN
                    PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(IGC),N)
                  ELSE
                    ZMC_M(IGC,M,N) = ZMC_M(IGC,M,N)+DNR_M(IEQGC(IGC),N)
                  ENDIF
                ENDIF
  240         CONTINUE
            ENDIF
!
!---      Phase condition #5  ---
!
!         SL_M = 0.0, ST >= 0.0, SG_M > 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous saturation, SL_M
!         CO2 mass
!           if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!           endif
!         CH4 mass
!           if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!           endif
!         petroleum component mass
!           if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!           endif
!         salt mass - total salt mass - TMS_M(2,N)
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.5 ) THEN
            IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
              T_M(M,N) = T_M(M,N) + DNR_M(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              SL_M(M,N) = SL_M(M,N) + DNR_M(IEQW,N)
            ELSEIF( M.EQ.IEQGC(1)+2 .AND. ISLC(45).EQ.0 ) THEN
              IF( IZMC(N).EQ.1 ) THEN
                PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(1),N)
              ELSE
                ZMC_M(1,M,N) = ZMC_M(1,M,N) + DNR_M(IEQGC(1),N)
              ENDIF
            ELSEIF( M.EQ.IEQGC(2)+2 .AND. ISLC(64).EQ.0 ) THEN
              IF( IZMC(N).EQ.2 ) THEN
                PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(2),N)
              ELSE
                ZMC_M(2,M,N) = ZMC_M(2,M,N) + DNR_M(IEQGC(2),N)
              ENDIF
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              TMS_M(M,N) = TMS_M(M,N) + DNR_M(IEQS,N)
            ELSE
              DO 250 IGC = 3,NGC+2
                IF( M.EQ.IEQGC(IGC)+2 ) THEN
                  IF( IZMC(N).EQ.IGC ) THEN
                    PSO_M(M,N) = PSO_M(M,N) + DNR_M(IEQGC(IGC),N)
                  ELSE
                    ZMC_M(IGC,M,N) = ZMC_M(IGC,M,N)+DNR_M(IEQGC(IGC),N)
                  ENDIF
                ENDIF
  250         CONTINUE
            ENDIF
          ENDIF
  280   CONTINUE
  300 CONTINUE
!$NO_OMP END PARALLEL DO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_EOR_M group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBBM_F( RSS,RSP,RST,N,MEQ )
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
!     Load the Jacobian matrix for the fracture-matrix connection 
!      of the dual porosity domain.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 December 2015.
!
#ifdef petsc
!-----------------------PETSc Module-----------------------------------!
!
      USE PETSC_STOMP
!
#endif
#ifdef lis
!------------------------Lis Module------------------------------------!
!
      USE STOMP_LIS_MODULE
!
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE DUAL_POR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#endif
#ifdef lis
#include "include/lisf.h"
#endif
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP(LUK),RST(LUK)
#ifdef petsc
      REAL*8 BUFFER(LUK)
      INTEGER IROW(1), ICOL(LUK)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBBM_F'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Fracture Node  ---
!
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        DO 100 M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSP(M)-RSS)/DNR(M,N)
  100   CONTINUE
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Adjacent Matrix Node ---
!
        DO 200 M = 1,ISVC
          MCOL = IM_M(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RST(M)-RSS)/DNR_M(M,N)
  200   CONTINUE
!
!---  SPLib or Lis solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    Fracture Node  ---
!
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        MA = 3*ISVC
        DO 2100 M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR(M,N)
 2100   CONTINUE
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Adjacent Matrix Node ---
!
        MA = 7*ISVC
        DO 2200 M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RST(M)-RSS)/DNR_M(M,N)
 2200     CONTINUE
#ifdef petsc
!
!---  PETSc solver  ---
!
      ELSEIF( ILES.EQ.5 ) THEN
!
!---    Fracture Node  ---
!
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        IROW(1) = MP-1
        DO 5100 M = 1,ISVC
          ICOL(M) = IM(M,NMD)-1
          BUFFER(M) = (RSP(M)-RSS)/DNR(M,N)
 5100   CONTINUE
        CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                    ADD_VALUES,IERR)
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Adjacent Matrix Node ---
!
        DO 5200 M = 1,ISVC
          ICOL(M) = IM_M(M,NMD)-1
          BUFFER(M) = (RST(M)-RSS)/DNR_M(M,N)
 5200   CONTINUE
        CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,ADD_VALUES,IERR)
#endif
      ELSE
        INDX = 3
        CHMSG = 'Unknown Linear Equation Solver'
        CALL WRMSGS( INDX )
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBBM_F group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBBM_M( RSS,RSP,RST,N,MEQ )
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
!     Load the Jacobian matrix for the matrix component of the 
!     dual porosity domain.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 December 2015.
!
#ifdef petsc
!-----------------------PETSc Module-----------------------------------!
!
      USE PETSC_STOMP
!
#endif
#ifdef lis
!------------------------Lis Module------------------------------------!
!
      USE STOMP_LIS_MODULE
!
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE DUAL_POR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#endif
#ifdef lis
#include "include/lisf.h"
#endif
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP(LUK),RST(LUK)
#ifdef petsc
      REAL*8 BUFFER(LUK)
      INTEGER IROW(1), ICOL(LUK)
#endif
!#ifdef lis
!      LIS_SCALAR VALUEX
!      LIS_INTEGER IROWX,ICOLX,IERR
!#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBBM_M'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Matrix Node  ---
!
        NMD = IXP(N)
        MP = IM_M(MEQ,NMD)
        DO 100 M = 1,ISVC
          MCOL = IM_M(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSP(M)-RSS)/DNR_M(M,N)
  100   CONTINUE
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Adjacent Fracture Node ---
!
        DO 200 M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RST(M)-RSS)/DNR(M,N)
  200   CONTINUE
!
!---  SPLib or Lis solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    Matrix Node  ---
!
        NMD = IXP(N)
        MP = IM_M(MEQ,NMD)
        DO 2100 M = 1,ISVC
          MCOL = KLU(MP,M)
          DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR_M(M,N)
 2100   CONTINUE
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Adjacent Fracture Node ---
!
        MA = ISVC
        DO 2200 M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RST(M)-RSS)/DNR(M,N)
 2200     CONTINUE
#ifdef petsc
!
!---  PETSc solver  ---
!
      ELSEIF( ILES.EQ.5 ) THEN
!
!---    Matrix Node  ---
!
        NMD = IXP(N)
        MP = IM_M(MEQ,NMD)
        IROW(1) = MP-1
        DO 5100 M = 1,ISVC
          ICOL(M) = IM_M(M,NMD)-1
          BUFFER(M) = (RSP(M)-RSS)/DNR_M(M,N)
 5100   CONTINUE
        CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                    ADD_VALUES,IERR)
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Adjacent Fracture Node ---
!
        DO 5200 M = 1,ISVC
          ICOL(M) = IM(M,NMD)-1
          BUFFER(M) = (RST(M)-RSS)/DNR(M,N)
 5200   CONTINUE
        CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,ADD_VALUES,IERR)
#endif
      ELSE
        INDX = 3
        CHMSG = 'Unknown Linear Equation Solver'
        CALL WRMSGS( INDX )
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBBM_M group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBC_EOR_M( IGC )
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
!     (zero flux boundary conditions) for the matrix component of the
!     dual porosity domain.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 December 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE DUAL_POR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STCX(LUK+1),RCPX(LUK),RCTX(LUK)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBC_EOR_M'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  FIRSTPRIVATE(IGC)
!$NO_OMP&  SHARED(DTI,IEQGC,ISVC,IXP,NFLD,QGC_FM,QNC_FM,TMC_M) 
!$NO_OMP&  PRIVATE(M,MM,MP,N,RCPX,RCSX,RCTX,STCX)
      DO 1000 N = 1,NFLD
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    First-order, forward-difference, time differential  ---
!
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STCX(M) = DTI*(TMC_M(IGC,MP,N)-TMC_M(IGC,1,N))
  100   CONTINUE
!
!---    Compute petroleum component equation residuals  ---
!
        RCSX = STCX(1) - QGC_FM(IGC,1,N) - QNC_FM(IGC,1,N)
        DO 920 M = 1,ISVC
          MM = 2*M
          RCPX(M) = STCX(M+1) - QGC_FM(IGC,MM,N) - QNC_FM(IGC,MM,N)
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          RCTX(M) = RCSX + QGC_FM(IGC,1,N) + QNC_FM(IGC,1,N)
     &       - QGC_FM(IGC,MM,N) - QNC_FM(IGC,MM,N)
  940   CONTINUE
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBBM_M( RCSX,RCPX,RCTX,N,IEQGC(IGC) )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
 !$NO_OMP END PARALLEL DO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBC_EOR_M group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_EOR_M
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
!     aqueous contributions (zero flux boundary conditions) 
!     for the matrix component of the dual porosity domain.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 December 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE DUAL_POR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STSX(LUK+1),RSPX(LUK),RSTX(LUK)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBS_EOR_M'
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(NFLD,IXP,ISVC,TMS_M,VOL_M,DTI,QS_FM,IEQS) 
!$NO_OMP&  PRIVATE(N,M,MP,STSX,RSSX,MM,RSPX,RSTX)
      DO 1000 N = 1,NFLD
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    First-order, forward-difference, time differential  ---
!
        DO 10 M = 1,ISVC+1
          MP = M + 1
          STSX(M) = (TMS_M(MP,N)-TMS_M(1,N))*VOL_M(N)*DTI
   10   CONTINUE
!
!---    Salt equation residuals  ---
!
        RSSX = STSX(1) - QS_FM(1,N)
        DO 920 M = 1,ISVC
          MM = 2*M
          RSPX(M) = STSX(M+1) - QS_FM(MM,N)
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          RSTX(M) = RSSX + QS_FM(1,N) - QS_FM(MM,N)
  940   CONTINUE
!
!---  Jacobian matrix loader  --
!
        CALL JCBBM_M( RSSX,RSPX,RSTX,N,IEQS )
 1000 CONTINUE
 !$NO_OMP END PARALLEL DO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_EOR_M group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_EOR_M
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
!     (zero flux boundary conditions) for the matrix component of the
!     dual porosity domain.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 December 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE JACOB
      USE GRID
      USE GLB_PAR
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
      REAL*8 STTX(LUK+1),RTPX(LUK),RTTX(LUK)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT_EOR_M'
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(NFLD,IXP,IZ,RHOS,CPS,PORD,T,SL,RHOL,UEL,SG,RHOG,UEG,SN,
!$NO_OMP&    RHON,UEN,SS,RHOSP,HSP,ISVC,PORD_M,T_M,SL_M,RHOL_M,UEL_M,SG_M,
!$NO_OMP&    RHOG_M,UEG_M,SN_M,RHON_M,UEN_M,SS_M,RHOSP_M,HSP_M,DTI,VOL_M,
!$NO_OMP&    QQ_FM,IEQT) 
!$NO_OMP&  PRIVATE(N,IZN,USTOX,USLOX,USGOX,USNOX,USSOX,M,MP,USTX,USLX,
!$NO_OMP&    USGX,USNX,USSX,DUSTX,DUSLX,DUSGX,DUSNX,DUSSX,STTX,RTSX,
!$NO_OMP&    MM,RTPX,RTTX)
      DO 1000 N = 1,NFLD
        IF( IXP(N).EQ.0 ) CYCLE
        IZN = IZ(N)
!
!---  First-order, forward-difference, time differential  ---
!
        USTOX = RHOS(IZN)*CPS(IZN)*(1.D+0-PORD(1,N))*T(1,N)
        USLOX = PORD(1,N)*SL(1,N)*RHOL(1,N)*UEL(1,N)
        USGOX = PORD(1,N)*SG(1,N)*RHOG(1,N)*UEG(1,N)
        USNOX = PORD(1,N)*SN(1,N)*RHON(1,N)*UEN(1,N)
        USSOX = PORD(1,N)*SS(1,N)*RHOSP(1,N)*HSP(1,N)
        DO 100 M = 1,ISVC+1
          MP = M + 1
          USTX = RHOS(IZN)*CPS(IZN)*(1.D+0-PORD_M(MP,N))*T_M(MP,N)
          USLX = PORD_M(MP,N)*SL_M(MP,N)*RHOL_M(MP,N)*UEL_M(MP,N)
          USGX = PORD_M(MP,N)*SG_M(MP,N)*RHOG_M(MP,N)*UEG_M(MP,N)
          USNX = PORD_M(MP,N)*SN_M(MP,N)*RHON_M(MP,N)*UEN_M(MP,N)
          USSX = PORD_M(MP,N)*SS_M(MP,N)*RHOSP_M(MP,N)*HSP_M(MP,N)
          DUSTX = (USTX-USTOX)
          DUSLX = (USLX-USLOX)
          DUSGX = (USGX-USGOX)
          DUSNX = (USNX-USNOX)
          DUSSX = (USSX-USSOX)
          STTX(M) = DUSTX + DUSLX + DUSGX + DUSNX + DUSSX
          STTX(M) = STTX(M)*DTI*VOL_M(N)
  100   CONTINUE
!
!---  Compute energy equation residuals  ---
!
        RTSX = STTX(1) - QQ_FM(1,N)
        DO 920 M = 1,ISVC
          MM = 2*M
          RTPX(M) = STTX(M+1) - QQ_FM(MM,N)
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          RTTX(M) = RTSX + QQ_FM(1,N) - QQ_FM(MM,N)
  940   CONTINUE
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBBM_M( RTSX,RTPX,RTTX,N,IEQT )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
 !$NO_OMP END PARALLEL DO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_EOR_M group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_EOR_M
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
!     (zero flux boundary conditions) for the matrix component of the
!     dual porosity domain.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 December 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE DUAL_POR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STWX(LUK+1),RWPX(LUK),RWTX(LUK)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBW_EOR_M'
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(NFLD,IXP,PORD_M,XLW_M,RHOL_M,SL_M,XGW_M,RHOG_M,SG_M,
!$NO_OMP&    ISVC,DTI,VOL_M,QLW_FM,QGW_FM,IEQW) 
!$NO_OMP&  PRIVATE(N,USLOX,USGOX,M,MP,USLX,USGX,DUSLX,DUSGX,STWX,RWSX,
!$NO_OMP&    MM,RWPX,RWTX)
      DO 1000 N = 1,NFLD
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    First-order, forward-difference, time differential  ---
!
        USLOX = PORD_M(1,N)*XLW_M(1,N)*RHOL_M(1,N)*SL_M(1,N)
        USGOX = PORD_M(1,N)*XGW_M(1,N)*RHOG_M(1,N)*SG_M(1,N)
        DO 100 M = 1,ISVC+1
          MP = M + 1
          USLX = PORD_M(MP,N)*XLW_M(MP,N)*RHOL_M(MP,N)*SL_M(MP,N)
          USGX = PORD_M(MP,N)*XGW_M(MP,N)*RHOG_M(MP,N)*SG_M(MP,N)
          DUSLX = (USLX-USLOX)
          DUSGX = (USGX-USGOX)
          STWX(M) = DUSLX + DUSGX
          STWX(M) = STWX(M)*DTI*VOL_M(N)
  100   CONTINUE
!
!---  Compute water equation residuals  ---
!
        RWSX = STWX(1) - QLW_FM(1,N) - QGW_FM(1,N)
        DO 920 M = 1,ISVC
          MM = 2*M
          RWPX(M) = STWX(M+1) - QLW_FM(MM,N) - QGW_FM(MM,N)         
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          RWTX(M) = RWSX + QLW_FM(1,N) + QGW_FM(1,N)
     &      - QLW_FM(MM,N) - QGW_FM(MM,N)
  940   CONTINUE
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBBM_M( RWSX,RWPX,RWTX,N,IEQW )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
 !$NO_OMP END PARALLEL DO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_EOR_M group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_EOR_M
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Load the current time step values into the old time step
!     variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, March 17, 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE HYST
      USE GRID
      USE FDVGC
      USE DUAL_POR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LDO_EOR_M'
!
!---  Assign old time step values  ---
!
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(ASL,ASLMIN,C,CO,DFGC,DFGW_M,DFLA_M,DFLS_M,DFNC_M,
!$NO_OMP&    HG_M,HGA_M,HGW_M,HL_M,HLW_M,HN_M,HSP_M,ISLC,NEQC,NEQK,
!$NO_OMP&    NFLD,NGC,NPHAZ_M,NSOLU,NSPR,PERMRF_M,PG_M,PL_M,PN_M,
!$NO_OMP&    PORD_M,PORT_M,POSM_M,PSO_M,PVA_M,PVW_M,RHOG_M,RHOL_M,
!$NO_OMP&    RHOMG_M,RHOML_M,RHOMN_M,RHON_M,RHOSP_M,RKG_M,RKL_M,
!$NO_OMP&    RKN_M,SG_M,SL_M,SN_M,SP_C,SP_CO,SS_M,T_M,THKG_M,
!$NO_OMP&    THKL_M,THKN_M,TMC_M,TMS_M,TORG_M,TORL_M,TORN_M,
!$NO_OMP&    UEG_M,UEL_M,UEN_M,VISG_M,VISL_M,VISN_M,XGC_M,XGW_M,
!$NO_OMP&    XLA_M,XLS_M,XLW_M,XMGC_M,XMGW_M,XMLA_M,XMLS_M,XMLW_M,
!$NO_OMP&    XMNC_M,XNC_M,YLS_M,ZG_M,ZMC_M,ZN_M) 
!$NO_OMP&  PRIVATE(IGC,N,NEQ,NSL,NSP)
      DO 100 N = 1,NFLD
        DO 10 IGC = 1,NGC+2
          DFGC(IGC,1,N) = DFGC(IGC,2,N)
          DFNC_M(IGC,1,N) = DFNC_M(IGC,2,N)
          TMC_M(IGC,1,N) = TMC_M(IGC,2,N)
          XGC_M(IGC,1,N) = XGC_M(IGC,2,N)
          XMGC_M(IGC,1,N) = XMGC_M(IGC,2,N)
          XMNC_M(IGC,1,N) = XMNC_M(IGC,2,N)
          XNC_M(IGC,1,N) = XNC_M(IGC,2,N)
          ZMC_M(IGC,1,N) = ZMC_M(IGC,2,N)
   10   CONTINUE
        ASLMIN(1,N) = MIN( ASL(N),ASLMIN(2,N) )
        ASLMIN(2,N) = ASLMIN(1,N)
        DFGW_M(1,N) = DFGW_M(2,N)
        DFLA_M(1,N) = DFLA_M(2,N)
        DFLS_M(1,N) = DFLS_M(2,N)
!
!---    Isothermal option  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          HG_M(1,N) = HG_M(2,N)
          HGA_M(1,N) = HGA_M(2,N)
          HGW_M(1,N) = HGW_M(2,N)
          HL_M(1,N) = HL_M(2,N)
          HLW_M(1,N) = HLW_M(2,N)
          HN_M(1,N) = HN_M(2,N)
          HSP_M(1,N) = HSP_M(2,N)
          UEG_M(1,N) = UEG_M(2,N)
          UEL_M(1,N) = UEL_M(2,N)
          UEN_M(1,N) = UEN_M(2,N)
          THKG_M(1,N) = THKG_M(2,N)
          THKL_M(1,N) = THKL_M(2,N)
          THKN_M(1,N) = THKN_M(2,N)
        ENDIF
        NPHAZ_M(1,N) = NPHAZ_M(2,N)
        PERMRF_M(1,N) = PERMRF_M(2,N)
        PG_M(1,N) = PG_M(2,N)
        PL_M(1,N) = PL_M(2,N)
        PN_M(1,N) = PN_M(2,N)
        POSM_M(1,N) = POSM_M(2,N)
        PSO_M(1,N) = PSO_M(2,N)
        PORD_M(1,N) = PORD_M(2,N)
        PORT_M(1,N) = PORT_M(2,N)
        PVA_M(1,N) = PVA_M(2,N)
        PVW_M(1,N) = PVW_M(2,N)
        RHOG_M(1,N) = RHOG_M(2,N)
        RHOL_M(1,N) = RHOL_M(2,N)
        RHOMG_M(1,N) = RHOMG_M(2,N)
        RHOML_M(1,N) = RHOML_M(2,N)
        RHOMN_M(1,N) = RHOMN_M(2,N)
        RHON_M(1,N) = RHON_M(2,N)
        RHOSP_M(1,N) = RHOSP_M(2,N)
        RKG_M(1,N) = RKG_M(2,N)
        RKL_M(1,1,N) = RKL_M(1,2,N)
        RKL_M(2,1,N) = RKL_M(2,2,N)
        RKL_M(3,1,N) = RKL_M(3,2,N)
        RKN_M(1,N) = RKN_M(2,N)
        SG_M(1,N) = SG_M(2,N)
        SL_M(1,N) = SL_M(2,N)
        SN_M(1,N) = SN_M(2,N)
        SS_M(1,N) = SS_M(2,N)
        T_M(1,N) = T_M(2,N)
        TMS_M(1,N) = TMS_M(2,N)
        TORG_M(1,N) = TORG_M(2,N)
        TORL_M(1,N) = TORL_M(2,N)
        TORN_M(1,N) = TORN_M(2,N)
        VISG_M(1,N) = VISG_M(2,N)
        VISL_M(1,N) = VISL_M(2,N)
        VISN_M(1,N) = VISN_M(2,N)
        XGW_M(1,N) = XGW_M(2,N)
        XLA_M(1,N) = XLA_M(2,N)
        XLS_M(1,N) = XLS_M(2,N)
        XLW_M(1,N) = XLW_M(2,N)
        XMGW_M(1,N) = XMGW_M(2,N)
        XMLA_M(1,N) = XMLA_M(2,N)
        XMLS_M(1,N) = XMLS_M(2,N)
        XMLW_M(1,N) = XMLW_M(2,N)
        YLS_M(1,N) = YLS_M(2,N)
        ZG_M(1,N) = ZG_M(2,N)
        ZN_M(1,N) = ZN_M(2,N)
        DO 90 NSL = 1,NSOLU
          CO(N,NSL) = C(N,NSL)
   90   CONTINUE
#ifdef ecke
        DO 92 NEQ = 1,NEQC+NEQK
          NSL = NEQ + NSOLU 
          CO(N,NSL) = C(N,NSL)
   92   CONTINUE
        DO 94 NSP = 1,NSPR
          SP_CO(N,NSP) = SP_C(N,NSP)
   94   CONTINUE
#endif
  100 CONTINUE
!$NO_OMP END PARALLEL DO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_EOR_M group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_BO_M
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Compute hydrologic, thermodynamic and physical properties
!     for the black-oil option
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 20 May 2014
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FDVP
      USE DUAL_POR
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
      SUB_LOG(ISUB_LOG) = '/PROP_BO_M'
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(DFLS_M,EPSL,ISLC,ISVC,IXP,IZ,N_DB,NFLD,NPHAZ_M,
!$NO_OMP&    PATM,PCMP,PERMRF_M,PG_M,PL_M,PN_M,PORD_M,PORT_M,POSM_M,
!$NO_OMP&    PSO_M,PVW_M,RHOG_M,RHOL_M,RHON_M,RHOSP_M,RKG_M,RKL_M,
!$NO_OMP&    RKN_M,SG_M,SL_M,SN_M,SS_M,T_M,TABS,TMC_M,TMS_M,TORG_M,
!$NO_OMP&    TORL_M,TORN_M,VISG_M,VISL_M,VISN_M,VOL_M,WTMS,XGC_M,
!$NO_OMP&    XGW_M,XLS_M,XLW_M,XNC_M,YLS_M) 
!$NO_OMP&  PRIVATE(CMLAX,CMLSX,CPGTX,CPNLX,IZN,M,N,PORDX,PTX,PX,
!$NO_OMP&    RHOBX,TKX,XBSMX)
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        IZN = IZ(N)
        N_DB = N
!
!---    Loop over increment indices  ---
!
        DO 100 M = 2,ISVC+2
          TKX = T_M(M,N) + TABS
!
!---      Phase condition #1  ---
!
!         SL_M = 1.0, SG_M+SN_M = 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous (system) pressure, Pa gauge - PL_M,PSO_M
!         CO2 mass - CO2 vapor pressure - PVA_M
!         CH4 mass - total CH4 mass - TMC_M
!         petroleum component mass - total petroleum comp. mass - TMC_M
!         salt mass - total salt brine mass fraction - YLS_M
!
          IF( MOD(NPHAZ_M(2,N),10).EQ.1 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_M(M,N) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY( N,PX,PCMP(N),PORD_M(M,N),PORT_M(M,N) )
            PORD_M(M,N) = MAX( PORD_M(M,N),EPSL )
            PORT_M(M,N) = MAX( PORT_M(M,N),PORD_M(M,N) )
!
!---        Flash calculations  ---
!
            CALL BO_PC1( PORD_M(M,N),POSM_M(M,N),PVW_M(M,N),PX,
     &        RHOG_M(M,N),RHOL_M(M,N),RHON_M(M,N),SG_M(M,N),SL_M(M,N),
     &        SN_M(M,N),TMS_M(M,N),T_M(M,N),VOL_M(N),
     &        XGC_M(1,M,N),XGW_M(M,N),XLS_M(M,N),XLW_M(M,N),
     &        XNC_M(1,M,N),YLS_M(M,N) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR_M( SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_M(M,N) = PSO_M(M,N)
            PTX = PG_M(M,N) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_M(M,N) = PTX
              PL_M(M,N) = PN_M(M,N) - CPNLX
            ELSE
              PL_M(M,N) = PTX
              PN_M(M,N) = PL_M(M,N) + CPNLX
            ENDIF
!
!---      Phase condition #2  ---
!
!         SG_M = 0.0
!         SL_M > 0.0
!         SN_M > 0.0
!         Aqueous and nonaqueous liquid saturated
!
!         water mass - aqueous saturation, SL_M
!         gas mass - gas to oil ratio, POSM_M
!         oil mass - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         salt mass - total salt brine mass fraction - YLS_M
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.2 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_M(M,N) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY( N,PX,PCMP(N),PORD_M(M,N),PORT_M(M,N) )
            PORD_M(M,N) = MAX( PORD_M(M,N),EPSL )
            PORT_M(M,N) = MAX( PORT_M(M,N),PORD_M(M,N) )
!
!---        Flash calculations  ---
!
            CALL BO_PC2( PORD_M(M,N),POSM_M(M,N),PVW_M(M,N),PX,
     &        RHOG_M(M,N),RHOL_M(M,N),RHON_M(M,N),SG_M(M,N),SL_M(M,N),
     &        SN_M(M,N),TMC_M(1,M,N),TMS_M(M,N),T_M(M,N),VOL_M(N),
     &        XGC_M(1,M,N),XGW_M(M,N),XLS_M(M,N),XLW_M(M,N),
     &        XNC_M(1,M,N),YLS_M(M,N),IZ(N) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR_M( SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_M(M,N) = PSO_M(M,N)
            PTX = PG_M(M,N) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_M(M,N) = PTX
              PL_M(M,N) = PN_M(M,N) - CPNLX
            ELSE
              PL_M(M,N) = PTX
              PN_M(M,N) = PL_M(M,N) + CPNLX
            ENDIF
!
!---      Phase condition #3  ---
!
!         SG_M > 0.0
!         SL_M > 0.0
!         SN_M > 0.0
!         Three-phase
!
!         water mass - aqueous saturation, SL_M
!         gas mass - gas saturation, SG_M
!         oil mass - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         salt mass - total salt brine mass fraction - YLS_M
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.3 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_M(M,N) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY( N,PX,PCMP(N),PORD_M(M,N),PORT_M(M,N) )
            PORD_M(M,N) = MAX( PORD_M(M,N),EPSL )
            PORT_M(M,N) = MAX( PORT_M(M,N),PORD_M(M,N) )
!
!---        Flash calculations  ---
!
            CALL BO_PC3( PORD_M(M,N),POSM_M(M,N),PVW_M(M,N),PX,
     &        RHOG_M(M,N),RHOL_M(M,N),RHON_M(M,N),SG_M(M,N),SL_M(M,N),
     &        SN_M(M,N),TMC_M(1,M,N),TMS_M(M,N),T_M(M,N),VOL_M(N),
     &        XGC_M(1,M,N),XGW_M(M,N),XLS_M(M,N),XLW_M(M,N),
     &        XNC_M(1,M,N),YLS_M(M,N),IZ(N) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR_M( SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_M(M,N) = PSO_M(M,N)
            PTX = PG_M(M,N) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_M(M,N) = PTX
              PL_M(M,N) = PN_M(M,N) - CPNLX
            ELSE
              PL_M(M,N) = PTX
              PN_M(M,N) = PL_M(M,N) + CPNLX
            ENDIF
!
!---      Phase condition #4  ---
!
!         SG_M = 0.0
!         SL_M = 0.0
!         SN_M = 1.0
!         Nonaqueous liquid saturated
!
!         water mass - aqueous saturation, SL_M
!         gas mass - gas to oil ratio, POSM_M
!         oil mass - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!         salt mass - total salt mass - TMS_M(2,N)
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.4 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_M(M,N) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY( N,PX,PCMP(N),PORD_M(M,N),PORT_M(M,N) )
            PORD_M(M,N) = MAX( PORD_M(M,N),EPSL )
            PORT_M(M,N) = MAX( PORT_M(M,N),PORD_M(M,N) )
!
!---        Flash calculations  ---
!
            CALL BO_PC4( PORD_M(M,N),POSM_M(M,N),PVW_M(M,N),PX,
     &        RHOG_M(M,N),RHOL_M(M,N),RHON_M(M,N),SG_M(M,N),SL_M(M,N),
     &        SN_M(M,N),TMC_M(1,M,N),TMS_M(M,N),T_M(M,N),VOL_M(N),
     &        XGC_M(1,M,N),XGW_M(M,N),XLS_M(M,N),XLW_M(M,N),
     &        XNC_M(1,M,N),YLS_M(M,N),IZ(N) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR_M( SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_M(M,N) = PSO_M(M,N)
            PTX = PG_M(M,N) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_M(M,N) = PTX
              PL_M(M,N) = PN_M(M,N) - CPNLX
            ELSE
              PL_M(M,N) = PTX
              PN_M(M,N) = PL_M(M,N) + CPNLX
            ENDIF
!
!---      Phase condition #5  ---
!
!         SG_M > 0.0
!         SL_M = 0.0
!         SN_M > 0.0
!         Nonaqueous liquid and gas saturated
!
!         water mass - aqueous saturation, SL_M
!         gas mass - gas saturation, SG_M
!         oil mass - gas (system) pressure, Pa gauge - PG_M,PSO_M
!         salt mass - total salt mass - TMS_M(2,N)
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.5 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_M(M,N) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY( N,PX,PCMP(N),PORD_M(M,N),PORT_M(M,N) )
            PORD_M(M,N) = MAX( PORD_M(M,N),EPSL )
            PORT_M(M,N) = MAX( PORT_M(M,N),PORD_M(M,N) )
!
!---        Flash calculations  ---
!
            CALL BO_PC5( PORD_M(M,N),POSM_M(M,N),PVW_M(M,N),PX,
     &        RHOG_M(M,N),RHOL_M(M,N),RHON_M(M,N),SG_M(M,N),SL_M(M,N),
     &        SN_M(M,N),TMC_M(1,M,N),TMS_M(M,N),T_M(M,N),VOL_M(N),
     &        XGC_M(1,M,N),XGW_M(M,N),XLS_M(M,N),XLW_M(M,N),
     &        XNC_M(1,M,N),YLS_M(M,N),IZ(N) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR_M( SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_M(M,N) = PSO_M(M,N)
            PTX = PG_M(M,N) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_M(M,N) = PTX
              PL_M(M,N) = PN_M(M,N) - CPNLX
            ELSE
              PL_M(M,N) = PTX
              PN_M(M,N) = PL_M(M,N) + CPNLX
            ENDIF
          ENDIF
!
!---      Relative permeability  ---
!
          CALL RKG_EOR( RKG_M(M,N),SG_M(M,N),SL_M(M,N),SN_M(M,N),IZN )
          CALL RKL_EOR( RKL_M(1,M,N),SG_M(M,N),SL_M(M,N),SN_M(M,N),IZN )
          CALL RKN_EOR( RKG_M(M,N),RKL_M(1,M,N),RKN_M(M,N),SG_M(M,N),
     &      SL_M(M,N),SN_M(M,N),IZN )
!
!---      Gas viscosity, Pa s  ---
!
          CALL BO_VISG( RHOG_M(M,N),VISG_M(M,N) )
!
!---      Aqueous viscosity, Pa s
!         concentration of CO2 in mol/kg,
!         concentration of salt in mol/kg  ---
!
          CMLAX = 0.D+0
          CMLSX = 1.D+3*XLS_M(M,N)*RHOL_M(M,N)/WTMS
          CALL VISC_L( CMLAX,CMLSX,PX,T_M(M,N),VISL_M(M,N) )
!
!---      Nonaqueous-liquid viscosity, Pa s  ---
!
          CALL BO_VISN( PX,POSM_M(M,N),T_M(M,N),VISN_M(M,N) )
!
!---      Salt aqueous diffusion coefficient  ---
!
          CALL DIFC_LS( DFLS_M(M,N),T_M(M,N),VISL_M(M,N),XLS_M(M,N) )
!
!---      Precipitated NaCl density, saturation, and enthalpy  ---
!
          CALL DENS_S( T_M(M,N),PX,RHOSP_M(M,N) )
          CALL SOL_LS( T_M(M,N),XBSMX )
          CALL DENS_B( T_M(M,N),PX,XBSMX,RHOBX )
          SS_M(M,N) = MAX(YLS_M(M,N)-XBSMX,0.D+0)*RHOBX*SL_M(M,N)/
     &      RHOSP_M(M,N)
!
!---      Kozney-Carmen permeability reduction  ---
!
          PORDX = PORD_M(M,N)*(1.D+0-SS_M(M,N))
          PERMRF_M(M,N) = ((PORDX**3)/(1.D+0+PORDX))/
     &      ((PORD_M(M,N)**3)/(1.D+0+PORD_M(M,N)))
!
!---      Change in porosity as a function of precipitated salt  ---
!
          PORD_M(M,N) = PORD_M(M,N)*(1.D+0-SS_M(M,N))
!
!---      Aqueous, gas, and nonaqueous-liquid tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL_M(M,N),SG_M(M,N),
     &      SN_M(M,N),PORD_M(M,N),TORL_M(M,N),TORG_M(M,N),TORN_M(M,N) )
  100   CONTINUE
  200 CONTINUE
!$NO_OMP END PARALLEL DO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_BO_M group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_EOR_M
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Compute hydrologic, thermodynamic and physical properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, March 17, 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FDVP
      USE FDVGC
      USE DUAL_POR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 FKX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP_EOR_M'
!
!---  Black-oil option  ---
!
      IF( ISLC(68).EQ.0 ) THEN
        CALL PROP_BO_M
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compositional option  ---
!
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(DFGC,DFGW_M,DFLA_M,DFLS_M,DFNC_M,EPSL,FK,HG_M,HGA_M,
!$NO_OMP&    HGW_M,HL_M,HLW_M,HN_M,HSP_M,IFK,IPRF,ISLC,ISVC,IXP,IZ,IZMC,
!$NO_OMP&    N_DB,NFLD,NGC,NITER,NPHAZ_M,PATM,PCMP,PERMRF_M,PG_M,PL_M,
!$NO_OMP&    PN_M,POR,POR0_M,PORD_M,PORT_M,PSO_M,PVA_M,PVW_M,RHOG_M,
!$NO_OMP&    RHOL_M,RHOMG_M,RHOML_M,RHOMN_M,RHON_M,RHOSP_M,RKG_M,
!$NO_OMP&    RKL_M,RKN_M,SG_M,SL_M,SN_M,SS_M,T_M,TABS,THKG_M,THKL_M,
!$NO_OMP&    THKN_M,TMC_M,TMS_M,TORG_M,TORL_M,TORN_M,UEG_M,UEL_M,
!$NO_OMP&    UEN_M,VISG_M,VISL_M,VISN_M,VOL_M,XGC_M,XGW_M,XLA_M,
!$NO_OMP&    XLS_M,XLW_M,XMGC_M,XMGW_M,XMLA_M,XMLS_M,XMLW_M,XMNC_M,
!$NO_OMP&    XNC_M,YLS_M,ZG_M,ZMC_M,ZN_M) 
!$NO_OMP&  PRIVATE(BETAX,CMLAX,CMLSX,CPGTX,CPNLX,FGAX,FKX,FNAX,HAX,HBX,
!$NO_OMP&    HGWX,HNAX,IGC,ISRX,IZN,M,N,PEX,PSWX,PTX,PX,RHOBX,
!$NO_OMP&    THKLWX,TKX,XBSMX)
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        IZN = IZ(N)
        N_DB = N
#ifdef ecke
        POR0_M(1,N) = POR0_M(1,N)
        POR0_M(2,N) = POR0_M(2,N)
#endif
!
!---    K-factors fixed at third iterate values  ---
!
        IF( NITER.GT.3 ) THEN
          IFK = 1
!
!---    K-factors computed from fugacity coefficients  ---
!
        ELSE
          IFK = 0
        ENDIF
!
!---    Loop over increment indices  ---
!
        DO 100 M = 2,ISVC+2
          TKX = T_M(M,N) + TABS
!
!---      K-factors fixed at third iterate values  ---
!
          IF( IFK.EQ.1 ) THEN
            DO 10 IGC = 1,NGC+2
              FKX(IGC) = FK(IGC,N)
   10       CONTINUE
          ENDIF
!
!---      Phase condition #1  ---
!
!         SL_M = 1.0, SG_M+SN_M = 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous (system) pressure, Pa gauge - PL_M,PSO_M
!         CO2 mass - CO2 vapor pressure - PVA_M
!         CH4 mass - total CH4 mass - TMC_M
!         petroleum component mass - total petroleum comp. mass - TMC_M
!         salt mass - total salt brine mass fraction - YLS_M
!
          IF( MOD(NPHAZ_M(2,N),10).EQ.1 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_M(M,N) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY( N,PX,PCMP(N),PORD_M(M,N),PORT_M(M,N) )
            PORD_M(M,N) = MAX( PORD_M(M,N),EPSL )
            PORT_M(M,N) = MAX( PORT_M(M,N),PORD_M(M,N) )
!
!---        Flash calculations  ---
!
            CALL FLH_PC1( BETAX,FKX,GNIFTX,PORD_M(M,N),PVA_M(M,N),
     &        PVW_M(M,N),PX,RHOG_M(M,N),RHOL_M(M,N),RHOMG_M(M,N),
     &        RHOML_M(M,N),RHOMN_M(M,N),RHON_M(M,N),SG_M(M,N),SL_M(M,N),
     &        SN_M(M,N),TMC_M(1,M,N),TMS_M(M,N),T_M(M,N),VOL_M(N),
     &        XGC_M(1,M,N),XGW_M(M,N),XLA_M(M,N),XLS_M(M,N),XLW_M(M,N),
     &        XNC_M(1,M,N),XMGC_M(1,M,N),XMGW_M(M,N),XMLA_M(M,N),
     &        XMLS_M(M,N),XMLW_M(M,N),XMNC_M(1,M,N),YLS_M(M,N),
     &        ZG_M(M,N),ZMC_M(1,M,N),ZN_M(M,N),IZ(N),NPHAZ_M(2,N) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR_M( SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_M(M,N) = PSO_M(M,N)
            PTX = PG_M(M,N) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_M(M,N) = PTX
              PL_M(M,N) = PN_M(M,N) - CPNLX
            ELSE
              PL_M(M,N) = PTX
              PN_M(M,N) = PL_M(M,N) + CPNLX
            ENDIF
!
!---      Phase condition #2  ---
!
!         SL_M < 1.0, ST > 0.0, SG_M = 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous saturation, SL_M
!         CO2 mass
!           if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!           endif
!         CH4 mass
!           if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!           endif
!         petroleum component mass
!           if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!           endif
!         salt mass - total salt brine mass fraction - YLS_M
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.2 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_M(M,N) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY( N,PX,PCMP(N),PORD_M(M,N),PORT_M(M,N) )
            PORD_M(M,N) = MAX( PORD_M(M,N),EPSL )
            PORT_M(M,N) = MAX( PORT_M(M,N),PORD_M(M,N) )
!
!---        Flash calculations  ---
!
            CALL FLH_PC2( BETAX,FKX,GNIFTX,PORD_M(M,N),PVA_M(M,N),
     &        PVW_M(M,N),PX,RHOG_M(M,N),RHOL_M(M,N),RHOMG_M(M,N),
     &        RHOML_M(M,N),RHOMN_M(M,N),RHON_M(M,N),SG_M(M,N),SL_M(M,N),
     &        SN_M(M,N),TMC_M(1,M,N),TMS_M(M,N),T_M(M,N),VOL_M(N),
     &        XGC_M(1,M,N),XGW_M(M,N),XLA_M(M,N),XLS_M(M,N),XLW_M(M,N),
     &        XNC_M(1,M,N),XMGC_M(1,M,N),XMGW_M(M,N),XMLA_M(M,N),
     &        XMLS_M(M,N),XMLW_M(M,N),XMNC_M(1,M,N),YLS_M(M,N),
     &        ZG_M(M,N),ZMC_M(1,M,N),ZN_M(M,N),IZMC(N),IZ(N),
     &        NPHAZ_M(2,N) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR_M( SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_M(M,N) = PSO_M(M,N)
            PTX = PG_M(M,N) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_M(M,N) = PTX
              PL_M(M,N) = PN_M(M,N) - CPNLX
            ELSE
              PL_M(M,N) = PTX
              PN_M(M,N) = PL_M(M,N) + CPNLX
            ENDIF
!
!---      Phase condition #3  ---
!
!         SL_M < 1.0, ST > 0.0, SG_M > 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous saturation, SL_M
!         CO2 mass
!           if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!           endif
!         CH4 mass
!           if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!           endif
!         petroleum component mass
!           if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!           endif
!         salt mass - total salt brine mass fraction - YLS_M
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.3 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_M(M,N) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY( N,PX,PCMP(N),PORD_M(M,N),PORT_M(M,N) )
            PORD_M(M,N) = MAX( PORD_M(M,N),EPSL )
            PORT_M(M,N) = MAX( PORT_M(M,N),PORD_M(M,N) )
!
!---        Flash calculations  ---
!
            CALL FLH_PC3( BETAX,FKX,GNIFTX,PORD_M(M,N),PVA_M(M,N),
     &        PVW_M(M,N),
     &        PX,RHOG_M(M,N),RHOL_M(M,N),RHOMG_M(M,N),RHOML_M(M,N),
     &        RHOMN_M(M,N),RHON_M(M,N),SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        TMC_M(1,M,N),TMS_M(M,N),T_M(M,N),VOL_M(N),XGC_M(1,M,N),
     &        XGW_M(M,N),XLA_M(M,N),XLS_M(M,N),XLW_M(M,N),XNC_M(1,M,N),
     &        XMGC_M(1,M,N),XMGW_M(M,N),XMLA_M(M,N),XMLS_M(M,N),
     &        XMLW_M(M,N),XMNC_M(1,M,N),YLS_M(M,N),ZG_M(M,N),
     &        ZMC_M(1,M,N),ZN_M(M,N),IZMC(N),IZ(N),NPHAZ_M(2,N) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR_M( SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_M(M,N) = PSO_M(M,N)
            PTX = PG_M(M,N) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_M(M,N) = PTX
              PL_M(M,N) = PN_M(M,N) - CPNLX
            ELSE
              PL_M(M,N) = PTX
              PN_M(M,N) = PL_M(M,N) + CPNLX
            ENDIF
!
!---      Phase condition #4  ---
!
!         SL_M = 0.0, ST = 1.0, SG_M = 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous saturation, SL_M
!         CO2 mass
!           if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!           endif
!         CH4 mass
!           if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!           endif
!         petroleum component mass
!           if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!             - nonaqueous (system) pressure, Pa gauge - PN_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!           endif
!         salt mass - total salt mass - TMS_M(2,N)
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.4 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_M(M,N) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY( N,PX,PCMP(N),PORD_M(M,N),PORT_M(M,N) )
            PORD_M(M,N) = MAX( PORD_M(M,N),EPSL )
            PORT_M(M,N) = MAX( PORT_M(M,N),PORD_M(M,N) )
!
!---        Flash calculations  ---
!
            CALL FLH_PC4( BETAX,FKX,GNIFTX,PORD_M(M,N),PVA_M(M,N),
     &        PVW_M(M,N),
     &        PX,RHOG_M(M,N),RHOL_M(M,N),RHOMG_M(M,N),RHOML_M(M,N),
     &        RHOMN_M(M,N),RHON_M(M,N),SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        TMC_M(1,M,N),TMS_M(M,N),T_M(M,N),VOL_M(N),XGC_M(1,M,N),
     &        XGW_M(M,N),XLA_M(M,N),XLS_M(M,N),XLW_M(M,N),XNC_M(1,M,N),
     &        XMGC_M(1,M,N),XMGW_M(M,N),XMLA_M(M,N),XMLS_M(M,N),
     &        XMLW_M(M,N),XMNC_M(1,M,N),YLS_M(M,N),ZG_M(M,N),
     &        ZMC_M(1,M,N),ZN_M(M,N),IZMC(N),IZ(N),NPHAZ_M(2,N) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR_M( SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_M(M,N) = PSO_M(M,N)
            PTX = PG_M(M,N) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_M(M,N) = PTX
              PL_M(M,N) = PN_M(M,N) - CPNLX
            ELSE
              PL_M(M,N) = PTX
              PN_M(M,N) = PL_M(M,N) + CPNLX
            ENDIF
!
!---      Phase condition #5  ---
!
!         SL_M = 0.0, ST >= 0.0, SG_M > 0.0
!
!         energy - temperature, C - T_M
!         water mass - aqueous saturation, SL_M
!         CO2 mass
!           if( ZMC_M(1,2,N) > max( ZMC_M(J,2,N) ) J /= 1
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(1,2,N)
!           endif
!         CH4 mass
!           if( ZMC_M(2,2,N) > max( ZMC_M(J,2,N) ) J /= 2
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(2,2,N)
!           endif
!         petroleum component mass
!           if( ZMC_M(I,2,N) > max( ZMC_M(J,2,N) ) J /= I
!             - gas (system) pressure, Pa gauge - PG_M,PSO_M
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_M(I,2,N)
!           endif
!         salt mass - total salt mass - TMS_M(2,N)
!
          ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.5 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_M(M,N) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY( N,PX,PCMP(N),PORD_M(M,N),PORT_M(M,N) )
            PORD_M(M,N) = MAX( PORD_M(M,N),EPSL )
            PORT_M(M,N) = MAX( PORT_M(M,N),PORD_M(M,N) )
!
!---        Flash calculations  ---
!
            CALL FLH_PC5( BETAX,FKX,GNIFTX,PORD_M(M,N),PVA_M(M,N),
     &        PVW_M(M,N),
     &        PX,RHOG_M(M,N),RHOL_M(M,N),RHOMG_M(M,N),RHOML_M(M,N),
     &        RHOMN_M(M,N),RHON_M(M,N),SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        TMC_M(1,M,N),TMS_M(M,N),T_M(M,N),VOL_M(N),XGC_M(1,M,N),
     &        XGW_M(M,N),XLA_M(M,N),XLS_M(M,N),XLW_M(M,N),XNC_M(1,M,N),
     &        XMGC_M(1,M,N),XMGW_M(M,N),XMLA_M(M,N),XMLS_M(M,N),
     &        XMLW_M(M,N),XMNC_M(1,M,N),YLS_M(M,N),ZG_M(M,N),
     &        ZMC_M(1,M,N),ZN_M(M,N),IZMC(N),IZ(N),NPHAZ_M(2,N) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR_M( SG_M(M,N),SL_M(M,N),SN_M(M,N),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_M(M,N) = PSO_M(M,N)
            PTX = PG_M(M,N) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_M(M,N) = PTX
              PL_M(M,N) = PN_M(M,N) - CPNLX
            ELSE
              PL_M(M,N) = PTX
              PN_M(M,N) = PL_M(M,N) + CPNLX
            ENDIF
          ENDIF
!
!---      Relative permeability  ---
!
!          RKG_M(M,N) = SG_M(M,N)**2
!          RKL_M(1,M,N) = SL_M(M,N)**2
!          RKL_M(2,M,N) = SL_M(M,N)**2
!          RKL_M(3,M,N) = SL_M(M,N)**2
!          RKN_M(M,N) = SN_M(M,N)**2
          CALL RKG_EOR( RKG_M(M,N),SG_M(M,N),SL_M(M,N),SN_M(M,N),IZN )
          CALL RKL_EOR( RKL_M(1,M,N),SG_M(M,N),SL_M(M,N),SN_M(M,N),IZN )
          CALL RKN_EOR( RKG_M(M,N),RKL_M(1,M,N),RKN_M(M,N),SG_M(M,N),
     &      SL_M(M,N),SN_M(M,N),IZN )
!
!---      Gas viscosity, Pa s  ---
!
          CALL VISC_G( RHOMG_M(M,N),T_M(M,N),VISG_M(M,N),XMGC_M(1,M,N),
     &      XMGW_M(M,N) )
!
!---      Aqueous viscosity, Pa s
!         concentration of CO2 in mol/kg,
!         concentration of salt in mol/kg  ---
!
          CMLAX = 1.D+3*XMLA_M(M,N)*RHOML_M(M,N)/RHOL_M(M,N)
          CMLSX = 1.D+3*XMLS_M(M,N)*RHOML_M(M,N)/RHOL_M(M,N)
          CALL VISC_L( CMLAX,CMLSX,PX,T_M(M,N),VISL_M(M,N) )
!
!---      Nonaqueous-liquid viscosity, Pa s  ---
!
          CALL VISC_N( RHOMN_M(M,N),T_M(M,N),VISN_M(M,N),XMNC_M(1,M,N) )
!
!---      Component gas diffusion coefficients  ---
!
          CALL DIFC_G( DFGC(1,M,N),DFGW_M(M,N),PX,T_M(M,N),
     &      XMGC_M(1,M,N),XMGW_M(M,N) )
!
!---      CO2 aqueous diffusion coefficient  ---
!
          CALL DIFC_LA( DFLA_M(M,N),T_M(M,N),VISL_M(M,N) )
!
!---      Salt aqueous diffusion coefficient  ---
!
          CALL DIFC_LS( DFLS_M(M,N),T_M(M,N),VISL_M(M,N),XLS_M(M,N) )
!
!---      Petroleum component nonaqueous-liquid diffusion 
!         coefficients  ---
!
          CALL DIFC_N( DFNC_M(1,M,N),T_M(M,N),VISN_M(M,N) )
!
!---      Gas thermal conductivity, W/m K  ---
!
          CALL THK_G( PX,T_M(M,N),THKG_M(M,N),XMGC_M(1,M,N),
     &      XMGW_M(M,N) )
!
!---      Aqueous thermal conductivity  ---
!
          CALL THK_W( PX,T_M(M,N),THKLWX )
          CALL THK_B( T_M(M,N),XLS_M(M,N),THKLWX,THKL_M(M,N) )
!
!---      Nonaqueous thermal conductivity, W/m K  ---
!
          CALL THK_N( PX,T_M(M,N),THKN_M(M,N),XMNC_M(1,M,N) )
!
!---      Gas enthalpy and internal energy, W/kg  ---
!
          CALL ENTH_G( HG_M(M,N),HGA_M(M,N),HGW_M(M,N),PX,T_M(M,N),
     &      RHOG_M(M,N),UEG_M(M,N),XGC_M(1,M,N),XGW_M(M,N),
     &      XMGC_M(1,M,N),ZG_M(M,N) )
!
!---      Nonaqueous-liquid enthalpy and internal energy  ---
!
          CALL ENTH_N( HN_M(M,N),HNAX,PX,T_M(M,N),RHON_M(M,N),
     &      UEN_M(M,N),XNC_M(1,M,N),XMNC_M(1,M,N),ZN_M(M,N) )
!
!---      Pure water liquid enthalpy  ---
!
          CALL SP_W( T_M(M,N),PSWX )
          PEX = MAX( PSWX,PX )
          ISRX = 1
          CALL ENTH_W( T_M(M,N),PEX,HLW_M(M,N),HGWX,ISRX )
!
!---      Brine enthalpy  ---
!
          CALL ENTH_B( T_M(M,N),XLS_M(M,N),HLW_M(M,N),HBX )
!!
!!---      Dissolved CO2 contribution  ---
!!
!          SGX = (BETAX/RHOMG_M(M,N))/
!     &      (BETAX/RHOMG_M(M,N) + (1.D+0-BETAX)/RHOMN_M(M,N))
!          IF( SGX.LT.1.D-12 ) SGX = 0.D+0
!          SNX = 1.D+0 - SGX
!          FGAX = SGX*RHOG_M(M,N)*XGC_M(1,M,N)
!          FNAX = SNX*RHON_M(M,N)*XNC_M(1,M,N)
          IF( SG_M(M,N).GT.0.D+0 ) THEN
            FGAX = SG_M(M,N)*RHOG_M(M,N)*XGC_M(1,M,N)
          ELSE
            FGAX = 0.D+0
          ENDIF
          IF( SN_M(M,N).GT.0.D+0 ) THEN
            FNAX = SN_M(M,N)*RHON_M(M,N)*XNC_M(1,M,N)
          ELSE
            FNAX = 0.D+0
          ENDIF
          IF( (FGAX+FNAX)/EPSL.LT.EPSL ) THEN
            HAX = 0.D+0
          ELSE
            HAX = (HGA_M(M,N)*FGAX + HNAX*FNAX)/(FGAX + FNAX)
          ENDIF
!
!---      Aqueous enthalpy and internal energy  ---
!
          CALL ENTH_L( HAX,HBX,HL_M(M,N),T_M(M,N),XLA_M(M,N),
     &      XLS_M(M,N) )
          UEL_M(M,N) = HL_M(M,N) - PX/RHOL_M(M,N)
!
!---      Precipitated NaCl density, saturation, and enthalpy  ---
!
          CALL DENS_S( T_M(M,N),PX,RHOSP_M(M,N) )
          CALL SOL_LS( T_M(M,N),XBSMX )
          CALL DENS_B( T_M(M,N),PX,XBSMX,RHOBX )
          SS_M(M,N) = MAX(YLS_M(M,N)-XBSMX,0.D+0)*RHOBX*SL_M(M,N)/
     &      RHOSP_M(M,N)
!
!---      Precipitated NaCl enthalpy  ---
!
          CALL ENTH_S( T_M(M,N),HSP_M(M,N) )
!
!---      Permeability reduction factor  ---
!
          IF( IPRF(IZN).EQ.1 ) CALL PERM_R( SS_M(M,N),PERMRF_M(M,N),
     &      PORD_M(M,N),IZN )
!
!---      Change in porosity as a function of precipitated salt  ---
!
          PORD_M(M,N) = PORD_M(M,N)*(1.D+0-SS_M(M,N))
!
!---      Aqueous, gas, and nonaqueous-liquid tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL_M(M,N),SG_M(M,N),
     &      SN_M(M,N),PORD_M(M,N),TORL_M(M,N),TORG_M(M,N),TORN_M(M,N) )
!
!---      Save unincremented K-factors  ---
!
          IF( M.EQ.2 .AND. NITER.LE.3 ) THEN
            DO 90 IGC = 1,NGC+2
              FK(IGC,N) = FKX(IGC)
   90       CONTINUE
          ENDIF
  100   CONTINUE
  200 CONTINUE
!$NO_OMP END PARALLEL DO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_EOR_M group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKG_EOR_F( RKGX,SLX,SNX,IZN )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Gas relative permeability for fractures
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 23 December 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE DUAL_POR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RKG_EOR_F'
!
!---  Stone-I 3-phase relative permeability model  ---
!
      IF( IRPN(IZN).EQ.1 ) THEN
!
!---    van Genuchten (1980), Parker et al., (1987) submodel
!       with exponent from fit of Wang and Narasimham (1985) fracture
!       relative permeability function ---
!
        STX = SLX + SNX
        STRX = 0.D+0
        SGRX = 0.D+0
        ESTX = MIN( MAX( (STX-STRX)/(1.D+0-STRX-SGRX),0.D+0 ),1.D+0 )
        RKGX = SQRT(1.D+0-ESTX)*((1.D+0-(ESTX**(1.D+0/FRAC_P(17,IZN))))
     &    **(2.D+0*FRAC_P(17,IZN)))
!
!---  Stone-II 3-phase relative permeability model  ---
!
      ELSEIF( IRPN(IZN).EQ.2 ) THEN
!
!---    van Genuchten (1980), Parker et al., (1987) submodel
!       with exponent from fit of Wang and Narasimham (1985) fracture
!       relative permeability function ---
!
        STX = SLX + SNX
        STRX = 0.D+0
        SGRX = 0.D+0
        ESTX = MIN( MAX( (STX-STRX)/(1.D+0-STRX-SGRX),0.D+0 ),1.D+0 )
        RKGX = SQRT(1.D+0-ESTX)*((1.D+0-(ESTX**(1.D+0/FRAC_P(17,IZN))))
     &    **(2.D+0*FRAC_P(17,IZN)))
!
!---  Baker 3-phase relative permeability model  ---
!
      ELSEIF( IRPN(IZN).EQ.3 ) THEN
!
!---    van Genuchten (1980), Parker et al., (1987) submodel
!       with exponent from fit of Wang and Narasimham (1985) fracture
!       relative permeability function ---
!
        STX = SLX + SNX
        STRX = 0.D+0
        SGRX = 0.D+0
        ESTX = MIN( MAX( (STX-STRX)/(1.D+0-STRX-SGRX),0.D+0 ),1.D+0 )
        RKGX = SQRT(1.D+0-ESTX)*((1.D+0-(ESTX**(1.D+0/FRAC_P(17,IZN))))
     &    **(2.D+0*FRAC_P(17,IZN)))
!
!---  Unrecognized 3-Phase relative permeability model  ---
!
      ELSE
        CHMSG = 'Unrecognized 3-Phase Relative Permeability' //
     &      ' Model: IRPN(IZN) = '
        INDX = 12
        IMSG = IRPN(IZN)
        CALL WRMSGS( INDX )
      ENDIF
      RKGX = MAX( RKGX,0.D+0 )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKG_EOR_F group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKL_EOR_F( RKLX,SLX,IZN )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Aqueous relative permeability
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, March 17, 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE DUAL_POR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RKL_EOR_F'
!
!---  Stone-I 3-phase relative permeability model  ---
!
      IF( IRPN(IZN).EQ.1 ) THEN
!
!---    van Genuchten (1980), Parker et al., (1987) submodel
!       with exponents from fit of Wang and Narasimham (1985) fracture
!       relative permeability function ---
!
        SLRX = 0.D+0
        SNRX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SNRX),0.D+0 ),1.D+0 )
        RKLX(1) = (ESLX**FRAC_P(16,IZN))*((1.D+0-((1.D+0-(ESLX**
     &    (1.D+0/FRAC_P(17,IZN))))**FRAC_P(17,IZN)))**2)
!
!---  Stone-II 3-phase relative permeability model  ---
!
      ELSEIF( IRPN(IZN).EQ.2 ) THEN
!
!---    van Genuchten (1980), Parker et al., (1987) submodel
!       with exponents from fit of Wang and Narasimham (1985) fracture
!       relative permeability function ---
!
        SLRX = 0.D+0
        SNRX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SNRX),0.D+0 ),1.D+0 )
        RKLX(1) = (ESLX**FRAC_P(16,IZN))*((1.D+0-((1.D+0-(ESLX**
     &    (1.D+0/FRAC_P(17,IZN))))**FRAC_P(17,IZN)))**2)
!
!---  Baker 3-phase relative permeability model  ---
!
      ELSEIF( IRPN(IZN).EQ.3 ) THEN
!
!---    van Genuchten (1980), Parker et al., (1987) submodel
!       with exponents from fit of Wang and Narasimham (1985) fracture
!       relative permeability function ---
!
        SLRX = 0.D+0
        SNRX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SNRX),0.D+0 ),1.D+0 )
        RKLX(1) = (ESLX**FRAC_P(16,IZN))*((1.D+0-((1.D+0-(ESLX**
     &    (1.D+0/FRAC_P(17,IZN))))**FRAC_P(17,IZN)))**2)
!
!---  Unrecognized 3-phase permeability model  ---
!
      ELSE
        CHMSG = 'Unrecognized 3-Phase Relative Permeability' //
     &      ' Model: IRPN(IZN) = '
        INDX = 12
        IMSG = IRPN(IZN)
        CALL WRMSGS( INDX )
      ENDIF
!
!---  End of RKL_EOR_F group  ---
!
      RKLX(1) = MAX( RKLX(1),0.D+0 )
      RKLX(2) = MAX( RKLX(1),0.D+0 )
      RKLX(3) = MAX( RKLX(1),0.D+0 )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKL_EOR_F group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKN_EOR_F( RKGX,RKLX,RKNX,SGX,SLX,SNX,IZN )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Nonaqueous-liquid relative permeability
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, March 17, 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GLB_PAR
      USE DUAL_POR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RKN_EOR_F'
!
!---  Stone-I 3-phase relative permeability model  ---
!
      IF( IRPN(IZN).EQ.1 ) THEN
!
!---    Nonaqueous liquid-gas van Genuchten (1980), Parker et al. (1987)
!       submodel with exponents from fit of Wang and Narasimham (1985) 
!       fracture relative permeability function ---
!
        STX = SLX + SNX
        STRX = 0.D+0
        SGRX = 0.D+0
        ESTX = MIN( MAX( (STX-STRX)/(1.D+0-STRX-SGRX),0.D+0 ),1.D+0 )
        RKNGX = (ESTX**FRAC_P(16,IZN))*((1.D+0-((1.D+0-(ESTX**
     &    (1.D+0/FRAC_P(17,IZN))))**FRAC_P(17,IZN)))**2)
!
!---    Nonaqueous liquid-aqueous van Genuchten (1980), 
!       Parker et al. (1987) submodel with exponent from fit of 
!       Wang and Narasimham (1985) fracture relative permeability 
!       function ---
!
        SLRX = 0.D+0
        SNRX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SNRX),0.D+0 ),1.D+0 )
        RKNLX = SQRT(1.D+0-ESLX)*((1.D+0-(ESLX**(1.D+0/FRAC_P(17,IZN))))
     &    **(2.D+0*FRAC_P(17,IZN)))
!
!---    Stone-I 3-phase relative permeability model  ---
!
        SNMX = MIN( SLRX/3.D+0,SLRX,STRX )
        SNSX = MAX( MIN( (SNX-SNMX)/(1.D+0-SLRX-SNMX),1.D+0 ),0.D+0 )
        SLSX = MAX( MIN( (SLX-SLRX)/(1.D+0-SLRX-SNMX),1.D+0 ),0.D+0 )
        SGSX = MAX( MIN( SGX/(1.D+0-SLRX-SNMX),1.D+0 ),0.D+0 )
        IF( SNSX.LT.EPSL ) THEN
          RKNX = 0.D+0
        ELSE
          RKNX = SNSX*RKNLX*RKNGX/((1.D+0-SLSX)*(1.D+0-SGSX)*
     &      RPLC(7,IZN))
        ENDIF
!
!---  Stone-II 3-phase relative permeability model  ---
!
      ELSEIF( IRPN(IZN).EQ.2 ) THEN
!
!---    Nonaqueous liquid-gas van Genuchten (1980), Parker et al. (1987)
!       submodel with exponents from fit of Wang and Narasimham (1985) 
!       fracture relative permeability function ---
!
        STX = SLX + SNX
        STRX = 0.D+0
        SGRX = 0.D+0
        ESTX = MIN( MAX( (STX-STRX)/(1.D+0-STRX-SGRX),0.D+0 ),1.D+0 )
        RKNGX = (ESTX**FRAC_P(16,IZN))*((1.D+0-((1.D+0-(ESTX**
     &    (1.D+0/FRAC_P(17,IZN))))**FRAC_P(17,IZN)))**2)
!
!---    Nonaqueous liquid-aqueous van Genuchten (1980), 
!       Parker et al. (1987) submodel with exponent from fit of 
!       Wang and Narasimham (1985) fracture relative permeability 
!       function ---
!
        SLRX = 0.D+0
        SNRX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SNRX),0.D+0 ),1.D+0 )
        RKNLX = SQRT(1.D+0-ESLX)*((1.D+0-(ESLX**(1.D+0/FRAC_P(17,IZN))))
     &    **(2.D+0*FRAC_P(17,IZN)))
!
!---    Stone-II 3-phase relative permeability model  ---
!
        RKNX = RPLC(7,IZN)*(((RKNLX/RPLC(7,IZN))+RKLX(1))*
     &    ((RKNGX/RPLC(7,IZN))+RKGX) - (RKLX(1)+RKGX))
!
!---  Baker 3-phase relative permeability model  ---
!
      ELSEIF( IRPN(IZN).EQ.3 ) THEN
!
!---    Nonaqueous liquid-gas van Genuchten (1980), Parker et al. (1987)
!       submodel with exponents from fit of Wang and Narasimham (1985) 
!       fracture relative permeability function ---
!
        STX = SLX + SNX
        STRX = 0.D+0
        SGRX = 0.D+0
        ESTX = MIN( MAX( (STX-STRX)/(1.D+0-STRX-SGRX),0.D+0 ),1.D+0 )
        RKNGX = (ESTX**FRAC_P(16,IZN))*((1.D+0-((1.D+0-(ESTX**
     &    (1.D+0/FRAC_P(17,IZN))))**FRAC_P(17,IZN)))**2)
!
!---    Nonaqueous liquid-aqueous van Genuchten (1980), 
!       Parker et al. (1987) submodel with exponent from fit of 
!       Wang and Narasimham (1985) fracture relative permeability 
!       function ---
!
        SLRX = 0.D+0
        SNRX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SNRX),0.D+0 ),1.D+0 )
        RKNLX = SQRT(1.D+0-ESLX)*((1.D+0-(ESLX**(1.D+0/FRAC_P(17,IZN))))
     &    **(2.D+0*FRAC_P(17,IZN)))
!
!---    Baker 3-phase relative permeability model  ---
!
        DSLX = MAX((SLX-SLRX),0.D+0)
        DSGX = MAX((SGX-SGRX),0.D+0)
        IF( (DSLX+DSGX).GT.EPSL ) THEN
          RKNX = (DSLX*RKNLX + DSGX*RKNGX)/(DSLX+DSGX)
        ELSE
          RKNX = 5.D-1*(RKNLX+RKNGX)
        ENDIF
!
!---  Unrecognized 3-Phase relative permeability model  ---
!
      ELSE
        CHMSG = 'Unrecognized 3-Phase Relative Permeability' //
     &      ' Model: IRPN(IZN) = '
        INDX = 12
        IMSG = IRPN(IZN)
        CALL WRMSGS( INDX )
      ENDIF
      RKNX = MAX( RKNX,0.D+0 )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKN_EOR_F group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRANSF_EOR
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Inter-porosity transfer fluxes
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 25 November 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FDVS
      USE FDVP
      USE FDVN
      USE FDVGC
      USE FDVG
      USE DUAL_POR
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
      SUB_LOG(ISUB_LOG) = '/TRANSF_EOR'
!
!---  Loop over active nodes  ---
!
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(FRAC_P,ISVF,IXP,IZ,MNEG,MPOS,NFLD,NGC,ONE,PERM,
!$NO_OMP&    PG,PG_M,PL,PL_M,PN,PN_M,QGC_FM,QGW_FM,QLA_FM,QLW_FM,QNC_FM,
!$NO_OMP&    QS_FM,RHOG,RHOG_M,RHOL,RHOL_M,RHON,RHON_M,RKG,RKG_M,RKL,
!$NO_OMP&    RKL_M,RKN,RKN_M,THIRD,VISG,VISG_M,VISL,VISL_M,VISN,
!$NO_OMP&    VISN_M,XGC,XGC_M,XGW,XGW_M,XLA,XLA_M,XLS,XLS_M,
!$NO_OMP&    XLW,XLW_M,XNC,XNC_M) 
!$NO_OMP&  PRIVATE(DPGRX,DPGX,DPLRX,DPLX,DPNRX,DPNX,IGC,INDX,IZN,
!$NO_OMP&    M,MF,MM,N,PERMX,QMGX,QMLX,QMNX,RHOGX,RHOLX,RHONX,
!$NO_OMP&    RKGM,RKLM,RKNM,SIGMAX,VGX,VLX,VNX,XGCX,XGWX,
!$NO_OMP&    XLAX,XLSX,XLWX,XNCX)
      DO N = 1,NFLD
        IF( IXP(N).EQ.0 ) CYCLE
        IZN = IZ(N)
!
!---    Shape factor  ---
!
        SIGMAX = 4.D+0*((1.D+0/(FRAC_P(1,IZN)**2)) + 
     &    (1.D+0/(FRAC_P(3,IZN)**2)) + (1.D+0/(FRAC_P(5,IZN)**2)))
!
!---    Aqueous, gas, and nonaqueous-liquid mass flux between fracture
!       and matrix kg/s m^3 matrix volume  ---
!
        DPGRX = PG(2,N)-PG_M(2,N)
        DPLRX = PL(2,N)-PL_M(2,N)
        DPNRX = PN(2,N)-PN_M(2,N)
        DO M = 1,ISVF
          MF = MPOS(M)
          MM = MNEG(M)
          PERMX = (PERM(4,IZN)*PERM(5,IZN)*PERM(6,IZN))**THIRD
!
!---      Aqueous mass flux from fracture to matrix  ---
!
          DPLX = PL(MF,N)-PL_M(MM,N)
          INDX = 8
          RKLM = DIFMN(RKL(1,MF,N),RKL_M(1,MM,N),ONE,ONE,DPLRX,INDX)
          INDX = 5
          VLX = DIFMN(VISL(MF,N),VISL_M(MM,N),ONE,ONE,DPLRX,INDX)
          INDX = 2
          RHOLX = DIFMN( RHOL(MF,N),RHOL_M(MM,N),ONE,ONE,DPLRX,INDX)
          QMLX = SIGMAX*RHOLX*RKLM*PERMX*DPLX/VLX
!
!---      Gas mass flux from fracture to matrix  ---
!
          DPGX = PG(MF,N)-PG_M(MM,N)
          INDX = 9
          RKGM = DIFMN(RKG(MF,N),RKG_M(MM,N),ONE,ONE,DPGRX,INDX)
          INDX = 6
          VGX = DIFMN(VISG(MF,N),VISG_M(MM,N),ONE,ONE,DPGRX,INDX)
          INDX = 3
          RHOGX = DIFMN( RHOG(MF,N),RHOG_M(MM,N),ONE,ONE,DPGRX,INDX)
          QMGX = SIGMAX*RHOGX*RKGM*PERMX*DPGX/VGX
!
!---      Nonaqueous-liquid mass flux from fracture to matrix  ---
!
          DPNX = PN(MF,N)-PN_M(MM,N)
          INDX = 10
          RKNM = DIFMN(RKN(MF,N),RKN_M(MM,N),ONE,ONE,DPNRX,INDX)
          INDX = 7
          VNX = DIFMN(VISN(MF,N),VISN_M(MM,N),ONE,ONE,DPNRX,INDX)
          INDX = 4
          RHONX = DIFMN( RHON(MF,N),RHON_M(MM,N),ONE,ONE,DPNRX,INDX)
          QMNX = SIGMAX*RHONX*RKNM*PERMX*DPNX/VNX
!
!---      CO2 aqueous mass advective flux from fracture to matrix  ---
!
          INDX = 2
          XLAX = DIFMN(XLA(MF,N),XLA_M(MM,N),ONE,ONE,DPLRX,INDX)
          QLA_FM(M,N) = XLAX*QMLX
!
!---      Water aqueous mass advective flux from fracture to matrix  ---
!
          INDX = 2
          XLWX = DIFMN(XLW(MF,N),XLW_M(MM,N),ONE,ONE,DPGRX,INDX)
          QLW_FM(M,N) = XLWX*QMLX
!
!---      Salt aqueous mass advective flux from fracture to matrix  ---
!
          INDX = 2
          XLSX = DIFMN(XLS(MF,N),XLS_M(MM,N),ONE,ONE,DPGRX,INDX)
          QS_FM(M,N) = XLSX*QMLX
!
!---      Water gas mass advective flux from fracture to matrix  ---
!
          INDX = 3
          XGWX = DIFMN(XGW(MF,N),XGW_M(MM,N),ONE,ONE,DPGRX,INDX)
          QGW_FM(M,N) = XGWX*QMGX
!
!---      Component gas mass advective flux from fracture to matrix  ---
!
          DO IGC = 1,NGC+2
            INDX = 3
            XGCX = DIFMN(XGC(IGC,MF,N),XGC_M(IGC,MM,N),ONE,ONE,
     &        DPGRX,INDX)
            QGC_FM(IGC,M,N) = XGCX*QMGX
          ENDDO
!
!---      Component nonaqueous advective mass flux from fracture 
!         to matrix  ---
!
          DO IGC = 1,NGC+2
            INDX = 4
            XNCX = DIFMN(XNC(IGC,MF,N),XNC_M(IGC,MM,N),ONE,ONE,
     &        DPNRX,INDX)
            QNC_FM(IGC,M,N) = XNCX*QMNX
          ENDDO
        ENDDO
      ENDDO
!$NO_OMP END PARALLEL DO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TRANSF_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_BO_M
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Update the primary variables for the black-oil option for the
!     matrix system.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 23 December 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE OUTPU
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FILES
      USE FDVP
      USE DUAL_POR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*128 PH_CND(5)
      REAL*8 DPC(LNGC)
      INTEGER MPC(LNGC)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Phase Condition #1, SL = 1.0, SG+SN = 0.0',
     &  'Phase Condition #2, SL > 0.0, SN > 0.0, SG = 0.0',
     &  'Phase Condition #3, SL > 0.0, SN > 0.0, SG > 0.0',
     &  'Phase Condition #4, SL = 0.0, SN > 0.0, SG = 0.0',
     &  'Phase Condition #5, SL = 0.0, SN > 0.0, SG > 0.0'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_BO_M'
      IF( ICNV.EQ.1 ) GOTO 300
      IERR = 0
!
!---  Update primary variables
!
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(BLU,EPSL,IEQGC,IEQS,IEQW,IM_M,ISLC,IXP,IZ,N_DB,NFLD,
!$NO_OMP&    NPHAZ_M,NSD,PATM,PG_M,PL_M,PN_M,POSM_M,PSO_M,SG_M,SL_M,T_M,
!$NO_OMP&    TABS,TCRW,TMC_M,TMS_M,YLS_M) 
!$NO_OMP&  PRIVATE(DPC,DPS,DPW,DPX,IERR,IZN,MPC,MPS,MPW,N,NMD,PX,TKX,XLSMX)
      DO 200 N = 1,NFLD
        N_DB = N
        IF( IXP(N).EQ.0 ) GOTO 200
        NMD = IXP(N)
        IZN = IZ(N)
!
!---    Water mass equation  ---
!
        MPW = IM_M(IEQW,NMD)
        DPW = BLU(MPW)
!
!---    Gas mass equation  ---
!
        MPC(1) = IM_M(IEQGC(1),NMD)
        DPC(1) = BLU(MPC(1))
!
!---    Oil mass equation  ---
!
        MPC(2) = IM_M(IEQGC(2),NMD)
        DPC(2) = BLU(MPC(2))
!
!---    Salt mass equation  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = IM_M(IEQS,NMD)
          DPS = BLU(MPS)
        ELSE
          DPS = 0.D+0
        ENDIF
!
!---    Phase condition #1  ---
!
!       SL = 1.0
!       SG+SN = 0.0
!       Aqueuous saturated
!
!       water mass - aqueous (system) pressure, Pa gauge - PL,PSO
!       gas mass - total gas mass - TMC(1,2,N)
!       oil mass - total oil mass - TMC(2,2,N)
!       salt mass - total salt brine mass fraction - YLS
!
        IF( MOD(NPHAZ_M(2,N),10).EQ.1 ) THEN
!
!---      Limit changes in total pressure,
!         excessive primary variable change limit at 100 MPa  ---
!
          IF( ABS(DPW).GT.1.D+8 ) IERR = 1
          DPX = 2.5D-2*(PSO_M(2,N)+PATM)
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          PSO_M(2,N) = PSO_M(2,N) + DPW
          PSO_M(2,N) = MIN( PSO_M(2,N),5.D+8 )
!
!---      Zero negative corrections for zero gas mass  ---
!
          IF( TMC_M(1,2,N).LT.EPSL .AND. DPC(1).LT.EPSL ) THEN
            BLU(MPC(1)) = 0.D+0
            DPC(1) = 0.D+0
          ENDIF
          TMC_M(1,2,N) = MAX( (TMC_M(1,2,N)+DPC(1)),0.D+0 )
!
!---      Zero negative corrections for zero oil mass  ---
!
          IF( TMC_M(2,2,N).LT.EPSL .AND. DPC(2).LT.EPSL ) THEN
            BLU(MPC(2)) = 0.D+0
            DPC(2) = 0.D+0
          ENDIF
          TMC_M(2,2,N) = MAX( (TMC_M(2,2,N)+DPC(2)),0.D+0 )
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            IF( YLS_M(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_M(2,N)/EPSL.LT.EPSL .AND. DPS/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_M(2,N) = YLS_M(2,N)+DPS
            IF( YLS_M(2,N).LT.1.D-12 ) YLS_M(2,N) = 0.D+0
          ENDIF
!
!---    Phase condition #2  ---
!
!       SG = 0.0
!       SL > 0.0
!       SN > 0.0
!       Aqueous and nonaqueous liquid saturated
!
!       water mass - aqueous saturation, SL
!       gas mass - gas to oil ratio, POSM
!       oil mass - nonaqueous (system) pressure, Pa gauge - PN,PSO
!       salt mass - total salt brine mass fraction - YLS
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.2 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_M(2,N) = SL_M(2,N) + DPW
          IF( SL_M(2,N).LT.EPSL ) SL_M(2,N) = 0.D+0
!
!---      Limit changes in gas to oil ratio, and zero negative 
!         corrections for zero gas to oil ratios  ---
!
          DPX = 5.D-2
          DPC(1) = SIGN( MIN(ABS(DPX),ABS(DPC(1))),DPC(1) )
          POSM_M(2,N) = POSM_M(2,N) + DPC(1)
!
!---      Limit changes in nonaqueous-liquid pressure (system pressure),
!         excessive primary variable change limit at 100 MPa  ---
!
          IF( ABS(DPC(2)).GT.1.D+8 ) IERR = 1
          DPX = MAX( 1.D+6,1.D-1*(PN_M(2,N)-PL_M(2,N)) )
          DPX = 2.5D-2*(PSO_M(2,N)+PATM)
          DPC(2) = SIGN( MIN(ABS(DPX),ABS(DPC(2))),DPC(2) )
!
!---      Bubble-point pressure as a function of gas to oil ratio
!         (m^3/m^3) and temperature (C)  ---
!
!          CALL BO_BP( PBPX,POSM_M(2,N),T_M(2,N) )
!
!---      Pressure below bubble-point pressure  ---
!
!          IF( PSO_M(2,N)+DPC(2).LT.PBPX-PATM ) DPC(2) = 6.D-1*DPC(2)
          PSO_M(2,N) = PSO_M(2,N) + DPC(2)
          PSO_M(2,N) = MIN( PSO_M(2,N),5.D+8 )
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            IF( YLS_M(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_M(2,N)/EPSL.LT.EPSL .AND. DPS/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_M(2,N) = YLS_M(2,N)+DPS
            IF( YLS_M(2,N).LT.1.D-12 ) YLS_M(2,N) = 0.D+0
          ENDIF
!
!---    Phase condition #3  ---
!
!       SG > 0.0
!       SL > 0.0
!       SN > 0.0
!       Three-phase
!
!       water mass - aqueous saturation, SL
!       gas mass - gas saturation, SG
!       oil mass - gas (system) pressure, Pa gauge - PG,PSO
!       salt mass - total salt brine mass fraction - YLS
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.3 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_M(2,N) = SL_M(2,N) + DPW
          IF( SL_M(2,N).LT.EPSL ) SL_M(2,N) = 0.D+0
!
!---      Limit changes in gas saturation  ---
!
          DPX = 2.5D-2
          DPC(1) = SIGN( MIN(ABS(DPX),ABS(DPC(1))),DPC(1) )
!          IF( SG_M(2,N)+DPC(1).LT.0.D+0 ) DPC(1) = 6.D-1*DPC(1)
          SG_M(2,N) = SG_M(2,N) + DPC(1)
          IF( SG_M(2,N).LT.EPSL ) SG_M(2,N) = 0.D+0
!
!---      Limit changes in gas pressure, excessive primary variable 
!         change limit at 100 MPa  ---
!
          IF( ABS(DPC(1)).GT.1.D+8 ) IERR = 1
          DPX = MAX( 1.D+6,1.D-1*(PG_M(2,N)-MAX(PL_M(2,N),PN_M(2,N))) )
          DPX = 2.5D-2*(PSO_M(2,N)+PATM)
          DPC(2) = SIGN( MIN(ABS(DPX),ABS(DPC(2))),DPC(2) )
          PSO_M(2,N) = PSO_M(2,N) + DPC(2)
          PSO_M(2,N) = MIN( PSO_M(2,N),5.D+8 )
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            IF( YLS_M(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_M(2,N)/EPSL.LT.EPSL .AND. DPS/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_M(2,N) = YLS_M(2,N)+DPS
            IF( YLS_M(2,N).LT.1.D-12 ) YLS_M(2,N) = 0.D+0
          ENDIF
!
!---    Phase condition #4  ---
!
!       SG = 0.0
!       SL = 0.0
!       SN = 1.0
!       Nonaqueous liquid saturated
!
!       water mass - aqueous saturation, SL
!       gas mass - gas to oil ratio, POSM
!       oil mass - nonaqueous (system) pressure, Pa gauge - PN,PSO
!       salt mass - total salt mass - TMS(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.4 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_M(2,N) = SL_M(2,N) + DPW
          IF( SL_M(2,N).LT.EPSL ) SL_M(2,N) = 0.D+0
!
!---      Limit changes in gas to oil ratio, and zero negative 
!         corrections for zero gas to oil ratios  ---
!
          DPX = 5.D-2
          DPC(1) = SIGN( MIN(ABS(DPX),ABS(DPC(1))),DPC(1) )
          POSM_M(2,N) = POSM_M(2,N) + DPC(1)
!
!---      Limit changes in nonaqueous-liquid pressure (system pressure),
!         excessive primary variable change limit at 100 MPa  ---
!
          IF( ABS(DPC(2)).GT.1.D+8 ) IERR = 1
          DPX = MAX( 1.D+6,1.D-1*(PN_M(2,N)-PL_M(2,N)) )
          DPX = 2.5D-2*(PSO_M(2,N)+PATM)
          DPC(2) = SIGN( MIN(ABS(DPX),ABS(DPC(2))),DPC(2) )
!
!---      Bubble-point pressure as a function of gas to oil ratio
!         (m^3/m^3) and temperature (C)  ---
!
!          CALL BO_BP( PBPX,POSM_M(2,N),T_M(2,N) )
!
!---      Pressure below bubble-point pressure  ---
!
!          IF( PSO_M(2,N)+DPC(2).LT.PBPX-PATM ) DPC(2) = 6.D-1*DPC(2)
          PSO_M(2,N) = PSO_M(2,N) + DPC(2)
          PSO_M(2,N) = MIN( PSO_M(2,N),5.D+8 )
!
!---      Isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Zero negative corrections for zero total salt mass  ---
!
            IF( TMS_M(2,N).LT.EPSL .AND. DPS.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            TMS_M(2,N) = TMS_M(2,N) + DPS
            IF( TMS_M(2,N).LT.EPSL ) TMS_M(2,N) = 0.D+0
          ENDIF
!
!---    Phase condition #5  ---
!
!       SG > 0.0
!       SL = 0.0
!       SN > 0.0
!       Nonaqueous liquid and gas saturated
!
!       water mass - aqueous saturation, SL
!       gas mass - gas saturation, SG
!       oil mass - gas (system) pressure, Pa gauge - PG,PSO
!       salt mass - total salt mass - TMS(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.5 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_M(2,N) = SL_M(2,N) + DPW
          IF( SL_M(2,N).LT.EPSL ) SL_M(2,N) = 0.D+0
!
!---      Limit changes in gas saturation  ---
!
          DPX = 2.5D-2
          DPC(1) = SIGN( MIN(ABS(DPX),ABS(DPC(1))),DPC(1) )
!          IF( SG_M(2,N)+DPC(1).LT.0.D+0 ) DPC(1) = 6.D-1*DPC(1)
          SG_M(2,N) = SG_M(2,N) + DPC(1)
          IF( SG_M(2,N).LT.EPSL ) SG_M(2,N) = 0.D+0
!
!---      Limit changes in gas pressure, excessive primary variable 
!         change limit at 100 MPa  ---
!
          IF( ABS(DPC(1)).GT.1.D+8 ) IERR = 1
          DPX = MAX( 1.D+6,1.D-1*(PG_M(2,N)-MAX(PL_M(2,N),PN_M(2,N))) )
          DPX = 2.5D-2*(PSO_M(2,N)+PATM)
          DPC(2) = SIGN( MIN(ABS(DPX),ABS(DPC(2))),DPC(2) )
          PSO_M(2,N) = PSO_M(2,N) + DPC(2)
          PSO_M(2,N) = MIN( PSO_M(2,N),5.D+8 )
!
!---      Isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Zero negative corrections for zero total salt mass  ---
!
            IF( TMS_M(2,N).LT.EPSL .AND. DPS.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            TMS_M(2,N) = TMS_M(2,N) + DPS
            IF( TMS_M(2,N).LT.EPSL ) TMS_M(2,N) = 0.D+0
          ENDIF
        ENDIF
!
!---    Check for excessive pressure or temperature   ---
!
        PX = PSO_M(2,N)+PATM
        TKX = T_M(2,N)+TABS
        IF( PX.GT.8.D+8 .OR. PX.LT.0.D+0 ) IERR = 1
        IF( TKX.GT.TCRW .OR. TKX.LT.TABS ) IERR = 1
        IF( IERR.EQ.1 ) NSD(1) = N
  200 CONTINUE
!$NO_OMP END PARALLEL DO
!
!---  Reduce time step for excessive changes in primary variables   ---
!
      IF( IERR.EQ.1 ) THEN
        ICNV = 1
        N = NSD(1)
        WRITE(ISC,'(10X,A)') '---  Excessive Matrix System Primary ' //
     &    'Variable Change  ---'
        WRITE(IWR,'(10X,A)') '---  Excessive Matrix System Primary ' //
     &    'Variable Change  ---'
        WRITE(ISC,'(4X,A,I6)') 'Node = ',N
        WRITE(IWR,'(4X,A,I6)') 'Node = ',N
        NPX = MOD(NPHAZ_M(2,N),10)
        WRITE(ISC,'(4X,2A)') 'Phase Condition = ',PH_CND(NPX)
        WRITE(IWR,'(4X,2A)') 'Phase Condition = ',PH_CND(NPX)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Temperature = ',T_M(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Temperature = ',T_M(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'System Pressure = ',PSO_M(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'System Pressure = ',PSO_M(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',PL_M(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',PL_M(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Pressure = ',PG_M(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Pressure = ',PG_M(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Nonaqueous Liquid Pressure = ',
     &    PN_M(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Nonaqueous Liquid Pressure = ',
     &    PN_M(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Saturation = ',SL_M(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Saturation = ',SL_M(2,N)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG_M(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG_M(2,N)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Nonaqueous-Liquid Saturation = ',
     &    SN_M(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Nonaqueous-Liquid Saturation = ',
     &    SN_M(2,N)
      ENDIF
!
!---  Reduce time step  ---
!
  300   CONTINUE
      IF( ICNV.EQ.1 ) THEN
        IF( NTSR.LT.4 .OR. (DTCF*DT).GT.DTMN ) THEN
          NTSR = NTSR + 1
          DTX = DT
          TM = TM - (1.D+0-DTCF)*DT
          DT = DTCF*DT
          DTO = DT
          DTI = 1.D+0/DT
          VAR = DT
          VARX = DTX
          IF( UNTM.NE.'null' ) THEN
            INDX = 1
            IUNS = 1
            CALL RDUNIT(UNTM,VAR,INDX)
            IUNS = 1
            CALL RDUNIT(UNTM,VARX,INDX)
            NCH = INDEX( UNTM,'  ')-1
          ENDIF
          WRITE(ISC,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
          WRITE(IWR,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  SHARED(NFLD,NGC,NPHAZ_M,PG_M,PL_M,PN_M,POSM_M,PSO_M,
!$NO_OMP&    SG_M,SL_M,SN_M,T_M,TMC_M,TMS_M,YLS_M) 
!$NO_OMP&  PRIVATE(IGC,N)
          DO 410 N = 1,NFLD
            T_M(2,N) = T_M(1,N)
            PL_M(2,N) = PL_M(1,N)
            PG_M(2,N) = PG_M(1,N)
            PN_M(2,N) = PN_M(1,N)
            POSM_M(2,N) = POSM_M(1,N)
            PSO_M(2,N) = PSO_M(1,N)
            SG_M(2,N) = SG_M(1,N)
            SL_M(2,N) = SL_M(1,N)
            SN_M(2,N) = SN_M(1,N)
            YLS_M(2,N) = YLS_M(1,N)
            TMS_M(2,N) = TMS_M(1,N)
            NPHAZ_M(2,N) = NPHAZ_M(1,N)
            DO 400 IGC = 1,NGC+2
              TMC_M(IGC,2,N) = TMC_M(IGC,1,N)
  400       CONTINUE
  410     CONTINUE
!$NO_OMP END PARALLEL DO
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  SHARED(NFLD,NGC,NPHAZ_M,PG_M,PL_M,PN_M,POSM_M,PSO_M,
!$NO_OMP&    SG_M,SL_M,SN_M,T_M,TMC_M,TMS_M,YLS_M) 
!$NO_OMP&  PRIVATE(IGC,N)
          DO 430 N = 1,NFLD
            T_M(2,N) = T_M(1,N)
            PL_M(2,N) = PL_M(1,N)
            PG_M(2,N) = PG_M(1,N)
            PN_M(2,N) = PN_M(1,N)
            POSM_M(2,N) = POSM_M(1,N)
            PSO_M(2,N) = PSO_M(1,N)
            SG_M(2,N) = SG_M(1,N)
            SL_M(2,N) = SL_M(1,N)
            SN_M(2,N) = SN_M(1,N)
            YLS_M(2,N) = YLS_M(1,N)
            TMS_M(2,N) = TMS_M(1,N)
            NPHAZ_M(2,N) = NPHAZ_M(1,N)
            DO 420 IGC = 1,NGC+2
              TMC_M(IGC,2,N) = TMC_M(IGC,1,N)
  420       CONTINUE
  430     CONTINUE
!$NO_OMP END PARALLEL DO
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          ICNV = 4
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_BO group
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_EOR_M
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Update the primary variables for the matrix system
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 23 December 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE OUTPU
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FILES
      USE FDVP
      USE FDVGC
      USE DUAL_POR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*128 PH_CND(5)
      REAL*8 DPC(LNGC)
      INTEGER MPC(LNGC)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Phase Condition #1, SL = 1.0, SG+SN = 0.0',
     &  'Phase Condition #2, SL < 1.0, ST > 0.0, SG = 0.0',
     &  'Phase Condition #3, SL < 1.0, ST > 0.0, SG > 0.0',
     &  'Phase Condition #4, SL = 0.0, ST = 1.0, SG = 0.0',
     &  'Phase Condition #5, SL = 0.0, ST = 1.0, SG > 0.0'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_EOR_M'
!
!---  Black-oil option  ---
!
      IF( ISLC(68).EQ.0 ) THEN
        CALL UPDT_BO
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compositional option  ---
!
      IF( ICNV.EQ.1 ) GOTO 300
      IERR = 0
!
!---  Update primary variables
!
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  COPYIN(ISUB_LOG)
!$NO_OMP&  SHARED(BLU,EPSL,IEQGC,IEQS,IEQT,IEQW,IM_M,ISLC,IXP,IZ,IZMC,
!$NO_OMP&    N_DB,NFLD,NGC,NPHAZ_M,NSD,PATM,PG_M,PL_M,PN_M,PSO_M,PVA_M,
!$NO_OMP&    SL_M,T_M,TABS,TCRW,TMC_M,TMS_M,YLS_M,ZMC_M) 
!$NO_OMP&  PRIVATE(DPC,DPS,DPT,DPW,DPX,IERR,IGC,IZN,MPC,MPS,MPT,MPW,
!$NO_OMP&    N,NMD,PX,SZMCX,TKX,XLSMX)
      DO 200 N = 1,NFLD
        N_DB = N
        IF( IXP(N).EQ.0 ) GOTO 200
        NMD = IXP(N)
        IZN = IZ(N)
!
!---    Energy equation  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          MPT = IM_M(IEQT,NMD)
          DPT = BLU(MPT)
        ELSE
          DPT = 0.D+0
        ENDIF
!
!---    Water mass equation  ---
!
        MPW = IM_M(IEQW,NMD)
        DPW = BLU(MPW)
!
!---    CO2 mass equation  ---
!
        IF( ISLC(45).EQ.0 ) THEN
          MPC(1) = IM_M(IEQGC(1),NMD)
          DPC(1) = BLU(MPC(1))
        ELSE
          DPC(1) = 0.D+0
        ENDIF
!
!---    CH4 mass equation  ---
!
        IF( ISLC(64).EQ.0 ) THEN
          MPC(2) = IM_M(IEQGC(2),NMD)
          DPC(2) = BLU(MPC(2))
        ELSE
          DPC(2) = 0.D+0
        ENDIF
!
!---    Petroleum component mass equation  ---
!
        DO 10 IGC = 3,NGC+2
          MPC(IGC) = IM_M(IEQGC(IGC),NMD)
          DPC(IGC) = BLU(MPC(IGC))
   10   CONTINUE
!
!---    Salt mass equation  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = IM_M(IEQS,NMD)
          DPS = BLU(MPS)
        ELSE
          DPS = 0.D+0
        ENDIF
!
!---    Limit changes in temperature to 0.5 C,
!       excessive primary variable change limit at 100 C  ---
!
!        IF( ABS(DPT).GT.1.D+2 ) IERR = 1
        DPT = SIGN( MIN( 5.D-1,ABS(DPT) ),DPT )
        T_M(2,N) = MAX( (T_M(2,N)+DPT),1.D-1 )
!
!---    Phase condition #1  ---
!
!       SL = 1.0, SG+SN = 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous (system) pressure, Pa gauge - PL,PSO
!       CO2 mass - CO2 vapor pressure - PVA
!       CH4 mass - total CH4 mass - TMC
!       petroleum component mass - total petroleum component mass - TMC
!       salt mass - total salt brine mass fraction - YLS
!
        IF( MOD(NPHAZ_M(2,N),10).EQ.1 ) THEN
!
!---      Limit changes in total pressure,
!         excessive primary variable change limit at 100 MPa  ---
!
!          IF( ABS(DPW).GT.1.D+8 ) IERR = 1
          DPX = MAX( 2.5D-2*(PSO_M(2,N)+PATM),1.D+6 )
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          PSO_M(2,N) = PSO_M(2,N) + DPW
          PSO_M(2,N) = MIN( PSO_M(2,N),5.D+8 )
!
!---      Zero negative corrections for zero aqueous CO2  ---
!
          IF( PVA_M(2,N).LT.EPSL .AND. BLU(MPC(1)).LT.EPSL ) THEN
            BLU(MPC(1)) = 0.D+0
            DPC(1) = 0.D+0
          ENDIF
          PVA_M(2,N) = MAX( (PVA_M(2,N)+DPC(1)),0.D+0 )
!
!---      Zero negative corrections for zero total petroleum
!         component mass  ---
!
          DO 20 IGC = 2,NGC+2
            IF( TMC_M(IGC,2,N).LT.EPSL .AND. BLU(MPC(IGC)).LT.EPSL) THEN
              BLU(MPC(IGC)) = 0.D+0
              DPC(IGC) = 0.D+0
            ENDIF
            TMC_M(IGC,2,N) = MAX( (TMC_M(IGC,2,N)+DPC(IGC)),0.D+0 )
   20     CONTINUE
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            IF( YLS_M(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_M(2,N)/EPSL.LT.EPSL .AND. DPS/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_M(2,N) = YLS_M(2,N)+DPS
            IF( YLS_M(2,N).LT.1.D-12 ) YLS_M(2,N) = 0.D+0
          ENDIF
!
!---    Phase condition #2  ---
!
!       SL < 1.0, ST > 0.0, SG = 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC(1,2,N) > max( ZMC(J,2,N) ) J /= 1
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC(2,2,N) > max( ZMC(J,2,N) ) J /= 2
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC(I,2,N) > max( ZMC(J,2,N) ) J /= I
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(I,2,N)
!         endif
!       salt mass - total salt brine mass fraction - YLS
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.2 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_M(2,N) = SL_M(2,N) + DPW
          IF( SL_M(2,N).LT.EPSL ) SL_M(2,N) = 0.D+0
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          SZMCX = 0.D+0
          DO 30 IGC = 1,NGC+2
!
!---        Maximum petroleum component, limit changes in nonaqueous
!           liquid pressure, excessive primary variable change 
!           limit at 100 MPa  ---
!
            IF( IGC.EQ.IZMC(N) ) THEN
!              IF( ABS(DPC(IGC)).GT.1.D+8 ) IERR = 1
              DPX = MAX( 1.D+6,1.D-1*(PN_M(2,N)-PL_M(2,N)) )
              DPX = 2.5D-2*(PSO_M(2,N)+PATM)
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              PSO_M(2,N) = PSO_M(2,N) + DPC(IGC)
              PSO_M(2,N) = MIN( PSO_M(2,N),5.D+8 )
!
!---        Limit changes in petroleum component mole fraction, 
!           and zero negative corrections for zero mole fractions---
!
            ELSE
              DPX = 5.D-2
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              ZMC_M(IGC,2,N) = ZMC_M(IGC,2,N) + DPC(IGC)
              IF( ZMC_M(IGC,2,N).LT.EPSL ) ZMC_M(IGC,2,N) = 0.D+0
              SZMCX = SZMCX + ZMC_M(IGC,2,N)
            ENDIF
   30     CONTINUE
          ZMC_M(IZMC(N),2,N) = MAX( 1.D+0-SZMCX,0.D+0 )
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            IF( YLS_M(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_M(2,N)/EPSL.LT.EPSL .AND. DPS/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_M(2,N) = YLS_M(2,N)+DPS
            IF( YLS_M(2,N).LT.1.D-12 ) YLS_M(2,N) = 0.D+0
          ENDIF
!
!---    Phase condition #3  ---
!
!       SL < 1.0, ST > 0.0, SG > 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC(1,2,N) > max( ZMC(J,2,N) ) J /= 1
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC(2,2,N) > max( ZMC(J,2,N) ) J /= 2
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC(I,2,N) > max( ZMC(J,2,N) ) J /= I
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(I,2,N)
!         endif
!       salt mass - total salt brine mass fraction - YLS
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.3 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_M(2,N) = SL_M(2,N) + DPW
          IF( SL_M(2,N).LT.EPSL ) SL_M(2,N) = 0.D+0
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          SZMCX = 0.D+0
          DO 40 IGC = 1,NGC+2
!
!---        Maximum petroleum component, limit changes in
!           gas pressure, excessive primary variable change 
!           limit at 100 MPa  ---
!
            IF( IGC.EQ.IZMC(N) ) THEN
!              IF( ABS(DPC(IGC)).GT.1.D+8 ) IERR = 1
              DPX = MAX( 1.D+6,1.D-1*(PG_M(2,N)-
     &          MAX(PL_M(2,N),PN_M(2,N))) )
              DPX = 2.5D-2*(PSO_M(2,N)+PATM)
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              PSO_M(2,N) = PSO_M(2,N) + DPC(IGC)
              PSO_M(2,N) = MIN( PSO_M(2,N),5.D+8 )
!
!---        Limit changes in petroleum component mole fraction, 
!           and zero negative corrections for zero mole fractions---
!
            ELSE
              DPX = 5.D-2
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              ZMC_M(IGC,2,N) = ZMC_M(IGC,2,N) + DPC(IGC)
              IF( ZMC_M(IGC,2,N).LT.EPSL ) ZMC_M(IGC,2,N) = 0.D+0
              SZMCX = SZMCX + ZMC_M(IGC,2,N)
            ENDIF
   40     CONTINUE
          ZMC_M(IZMC(N),2,N) = MAX( 1.D+0-SZMCX,0.D+0 )
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_M(2,N),XLSMX )
            IF( YLS_M(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_M(2,N)/EPSL.LT.EPSL .AND. DPS/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_M(2,N) = YLS_M(2,N)+DPS
            IF( YLS_M(2,N).LT.1.D-12 ) YLS_M(2,N) = 0.D+0
          ENDIF
!
!---    Phase condition #4  ---
!
!       SL = 0.0, ST = 1.0, SG = 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC(1,2,N) > max( ZMC(J,2,N) ) J /= 1
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC(2,2,N) > max( ZMC(J,2,N) ) J /= 2
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC(I,2,N) > max( ZMC(J,2,N) ) J /= I
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(I,2,N)
!         endif
!       salt mass - total salt mass - TMS(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.4 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_M(2,N) = SL_M(2,N) + DPW
          IF( SL_M(2,N).LT.EPSL ) SL_M(2,N) = 0.D+0
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          SZMCX = 0.D+0
          DO 50 IGC = 1,NGC+2
!
!---        Maximum petroleum component, limit changes in nonaqueous
!           liquid pressure, excessive primary variable change 
!           limit at 100 MPa  ---
!
            IF( IGC.EQ.IZMC(N) ) THEN
!              IF( ABS(DPC(IGC)).GT.1.D+8 ) IERR = 1
              DPX = MAX( 1.D+6,1.D-1*(PN_M(2,N)-PL_M(2,N)) )
              DPX = 2.5D-2*(PSO_M(2,N)+PATM)
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              PSO_M(2,N) = PSO_M(2,N) + DPC(IGC)
              PSO_M(2,N) = MIN( PSO_M(2,N),5.D+8 )
!
!---        Limit changes in petroleum component mole fraction, 
!           and zero negative corrections for zero mole fractions---
!
            ELSE
              DPX = 5.D-2
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              ZMC_M(IGC,2,N) = ZMC_M(IGC,2,N) + DPC(IGC)
              IF( ZMC_M(IGC,2,N).LT.EPSL ) ZMC_M(IGC,2,N) = 0.D+0
              SZMCX = SZMCX + ZMC_M(IGC,2,N)
            ENDIF
   50     CONTINUE
          ZMC_M(IZMC(N),2,N) = MAX( 1.D+0-SZMCX,0.D+0 )
!
!---      Isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Zero negative corrections for zero total salt mass  ---
!
            IF( TMS_M(2,N).LT.EPSL .AND. DPS.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            TMS_M(2,N) = TMS_M(2,N) + DPS
            IF( TMS_M(2,N).LT.EPSL ) TMS_M(2,N) = 0.D+0
          ENDIF
!
!---    Phase condition #5  ---
!
!       SL = 0.0, ST >= 0.0, SG > 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC(1,2,N) > max( ZMC(J,2,N) ) J /= 1
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(1,2,N)
!         endif
!       CH4 mass
!         if( ZMC(2,2,N) > max( ZMC(J,2,N) ) J /= 2
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(2,2,N)
!         endif
!       petroleum component mass
!         if( ZMC(I,2,N) > max( ZMC(J,2,N) ) J /= I
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC(I,2,N)
!         endif
!       salt mass - total salt mass - TMS(2,N)
!
        ELSEIF( MOD(NPHAZ_M(2,N),10).EQ.5 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_M(2,N) = SL_M(2,N) + DPW
          IF( SL_M(2,N).LT.EPSL ) SL_M(2,N) = 0.D+0
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          SZMCX = 0.D+0
          DO 60 IGC = 1,NGC+2
!
!---        Maximum petroleum component, limit changes in
!           gas pressure, excessive primary variable change 
!           limit at 100 MPa  ---
!
            IF( IGC.EQ.IZMC(N) ) THEN
!              IF( ABS(DPC(IGC)).GT.1.D+8 ) IERR = 1
              DPX = MAX( 1.D+6,1.D-1*(PG_M(2,N)-
     &          MAX(PL_M(2,N),PN_M(2,N))) )
              DPX = 2.5D-2*(PSO_M(2,N)+PATM)
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              PSO_M(2,N) = PSO_M(2,N) + DPC(IGC)
              PSO_M(2,N) = MIN( PSO_M(2,N),5.D+8 )
!
!---        Limit changes in petroleum component mole fraction, 
!           and zero negative corrections for zero mole fractions---
!
            ELSE
              DPX = 5.D-2
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              ZMC_M(IGC,2,N) = ZMC_M(IGC,2,N) + DPC(IGC)
              IF( ZMC_M(IGC,2,N).LT.EPSL ) ZMC_M(IGC,2,N) = 0.D+0
              SZMCX = SZMCX + ZMC_M(IGC,2,N)
            ENDIF
   60     CONTINUE
          ZMC_M(IZMC(N),2,N) = MAX( 1.D+0-SZMCX,0.D+0 )
!
!---      Isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Zero negative corrections for zero total salt mass  ---
!
            IF( TMS_M(2,N).LT.EPSL .AND. DPS.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            TMS_M(2,N) = TMS_M(2,N) + DPS
            IF( TMS_M(2,N).LT.EPSL ) TMS_M(2,N) = 0.D+0
          ENDIF
        ENDIF
!
!---    Check for excessive pressure or temperature   ---
!
        PX = PSO_M(2,N)+PATM
        TKX = T_M(2,N)+TABS
        IF( PX.GT.8.D+8 .OR. PX.LT.0.D+0 ) IERR = 1
        IF( TKX.GT.TCRW .OR. TKX.LT.TABS ) IERR = 1
        IF( IERR.EQ.1 ) NSD(1) = N
  200 CONTINUE
!$NO_OMP END PARALLEL DO
!
!---  Reduce time step for excessive changes in primary variables   ---
!
      IF( IERR.EQ.1 ) THEN
        ICNV = 1
        N = NSD(1)
        WRITE(ISC,'(10X,A)') '---  Excessive Matrix System Primary ' //
     &    'Variable Change  ---'
        WRITE(IWR,'(10X,A)') '---  Excessive Matrix System Primary ' //
     &    'Variable Change  ---'
        WRITE(ISC,'(4X,A,I6)') 'Node = ',N
        WRITE(IWR,'(4X,A,I6)') 'Node = ',N
        NPX = MOD(NPHAZ_M(2,N),10)
        WRITE(ISC,'(4X,2A)') 'Phase Condition = ',PH_CND(NPX)
        WRITE(IWR,'(4X,2A)') 'Phase Condition = ',PH_CND(NPX)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Temperature = ',T_M(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Temperature = ',T_M(2,N)
        WRITE(ISC,'(4X,A,1PE12.5)') 'System Pressure = ',PSO_M(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'System Pressure = ',PSO_M(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',PL_M(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',PL_M(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Pressure = ',PG_M(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Pressure = ',PG_M(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Nonaqueous Liquid Pressure = ',
     &    PN_M(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Nonaqueous Liquid Pressure = ',
     &    PN_M(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Saturation = ',SL_M(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Saturation = ',SL_M(2,N)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG_M(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG_M(2,N)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Nonaqueous-Liquid Saturation = ',
     &    SN_M(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Nonaqueous-Liquid Saturation = ',
     &    SN_M(2,N)
      ENDIF
!
!---  Reduce time step  ---
!
  300   CONTINUE
      IF( ICNV.EQ.1 ) THEN
        IF( NTSR.LT.4 .OR. (DTCF*DT).GT.DTMN ) THEN
          NTSR = NTSR + 1
          DTX = DT
          TM = TM - (1.D+0-DTCF)*DT
          DT = DTCF*DT
          DTO = DT
          DTI = 1.D+0/DT
          VAR = DT
          VARX = DTX
          IF( UNTM.NE.'null' ) THEN
            INDX = 1
            IUNS = 1
            CALL RDUNIT(UNTM,VAR,INDX)
            IUNS = 1
            CALL RDUNIT(UNTM,VARX,INDX)
            NCH = INDEX( UNTM,'  ')-1
          ENDIF
          WRITE(ISC,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
          WRITE(IWR,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  SHARED(NFLD,NGC,NPHAZ_M,PG_M,PL_M,PN_M,PSO_M,PVA_M,SG_M,
!$NO_OMP&    SL_M,SN_M,T_M,TMC_M,TMS_M,YLS_M,ZMC_M) 
!$NO_OMP&  PRIVATE(IGC,N)
          DO 410 N = 1,NFLD
            T_M(2,N) = T_M(1,N)
            PL_M(2,N) = PL_M(1,N)
            PG_M(2,N) = PG_M(1,N)
            PN_M(2,N) = PN_M(1,N)
            PSO_M(2,N) = PSO_M(1,N)
            PVA_M(2,N) = PVA_M(1,N)
            SG_M(2,N) = SG_M(1,N)
            SL_M(2,N) = SL_M(1,N)
            SN_M(2,N) = SN_M(1,N)
            YLS_M(2,N) = YLS_M(1,N)
            TMS_M(2,N) = TMS_M(1,N)
            NPHAZ_M(2,N) = NPHAZ_M(1,N)
            DO 400 IGC = 1,NGC+2
              TMC_M(IGC,2,N) = TMC_M(IGC,1,N)
              ZMC_M(IGC,2,N) = ZMC_M(IGC,1,N)
  400       CONTINUE
  410     CONTINUE
!$NO_OMP END PARALLEL DO
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
!$NO_OMP PARALLEL DO
!$NO_OMP&  DEFAULT(NONE)
!$NO_OMP&  SHARED(NFLD,NGC,NPHAZ_M,PG_M,PL_M,PN_M,PSO_M,PVA_M,SG_M,SL_M,
!$NO_OMP&    SN_M,T_M,TMC_M,TMS_M,YLS_M,ZMC_M) 
!$NO_OMP&  PRIVATE(IGC,N)
          DO 430 N = 1,NFLD
            T_M(2,N) = T_M(1,N)
            PL_M(2,N) = PL_M(1,N)
            PG_M(2,N) = PG_M(1,N)
            PN_M(2,N) = PN_M(1,N)
            PSO_M(2,N) = PSO_M(1,N)
            PVA_M(2,N) = PVA_M(1,N)
            SG_M(2,N) = SG_M(1,N)
            SL_M(2,N) = SL_M(1,N)
            SN_M(2,N) = SN_M(1,N)
            YLS_M(2,N) = YLS_M(1,N)
            TMS_M(2,N) = TMS_M(1,N)
            NPHAZ_M(2,N) = NPHAZ_M(1,N)
            DO 420 IGC = 1,NGC+2
              TMC_M(IGC,2,N) = TMC_M(IGC,1,N)
              ZMC_M(IGC,2,N) = ZMC_M(IGC,1,N)
  420       CONTINUE
  430     CONTINUE
!$NO_OMP END PARALLEL DO
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          ICNV = 4
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_EOR_M group
!
      RETURN
      END

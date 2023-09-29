!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ATMOS_C( PA,RH,RN,TA,UZ,RF,IERR )
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
!     Atmosphseric conditions
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 August 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
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
      SUB_LOG(ISUB_LOG) = '/ATMOS_C'
!
!---  Set local time variable  ---
!
      IERR = 0
      TMZ = TM
      TMZO = TM-DT
      IF( NSTEP-NRST.EQ.0 ) THEN
        TMZ = TMZ*(1.D+0+EPSL)+EPSL
        TMZO = TMZ
      ENDIF
!
!---  Cyclic atmospheric conditions  ---
!
      IF( IATM_C.EQ.1 ) THEN
        TMZ = MOD( TM,ATMOS(NATM_T,1) )
        TMZO =  MOD( TM-DT,ATMOS(NATM_T,1) )
      ENDIF
!
!---  Simulation time less than starting atmospheric
!     conditions time  ---
!
      IF( TMZ.LE.ATMOS(1,1) ) THEN
        IERR = 1
        GOTO 150
      ENDIF
!
!---  Single atmospheric conditions time  ---
!
      IF( NATM_T.EQ.1 ) THEN
        TA = ATMOS(1,2)
        PA = ATMOS(1,3)
        RH = ATMOS(1,4)
        RN = ATMOS(1,5)
        UZ = ATMOS(1,6)
!
!---    Assume precipitation based on 1 hour recording  ---
!
        RF = ATMOS(1,7)/3.6D+3
!
!---  Multiple atmospheric conditions times  ---
!
      ELSE
        CALL LOCATE( ATMOS(1,1),NATM_T,TMZ,NT )
!
!---    Simulation time prior to first atmospheric conditions time  ---
!
        IF( NT.EQ.0 ) THEN
          IERR = 1
          GOTO 150
!
!---    Simulation time after last atmospheric conditions time  ---
!
        ELSEIF( NT.EQ.NATM_T ) THEN
          IERR = 1
          GOTO 150
!
!---    Simulation time within atmospheric conditions time limits  ---
!
        ELSE
          TDX = ATMOS(NT+1,1)-ATMOS(NT,1)
          TFX = (TMZ-ATMOS(NT,1))/TDX
          TA = ATMOS(NT,2) + TFX*(ATMOS(NT+1,2)-ATMOS(NT,2))
          PA = ATMOS(NT,3) + TFX*(ATMOS(NT+1,3)-ATMOS(NT,3))
          RH = ATMOS(NT,4) + TFX*(ATMOS(NT+1,4)-ATMOS(NT,4))
          RN = ATMOS(NT,5) + TFX*(ATMOS(NT+1,5)-ATMOS(NT,5))
          UZ = ATMOS(NT,6) + TFX*(ATMOS(NT+1,6)-ATMOS(NT,6))
          RF = TFX*ATMOS(NT+1,7)
        ENDIF
      ENDIF
  150 CONTINUE
!
!---  Minimum wind speed, 1 cm/s  ---
!
      UZ = MAX( UZ,1.D-2 )
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ATMOS_C group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CLOUD( RN,F_CC,J_DAY,DZ_ANG )
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
!     Fractional cloud cover.
!
!     Llasat, M.C., and R.L. Snyder. 1998.  Data error effects on net
!     radiation and evapotranspiration estimation.  Agricultural
!     and Forest Meterology, 91:209-221.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 October 2003.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Data Statements---------------------------------!
!
      DATA G_SC / 1.367D+3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CLOUD'
!
!---  Incident solar radiation zero, skip calculations  ---
!
      IF( RN/EPSL.LT.EPSL ) GOTO 100
!
!---  Equation of time, hr  ---
!
!      B_EQT = 2.D+0*GPI*REAL(J_DAY-81)/364.D+0
!      T_EQT = 1.645D-1*SIN(2.D+0*B_EQT) - 1.255D-1*COS(B_EQT) -
!     &  2.5D-2*SIN(B_EQT)
      IF( J_DAY.LT.181 ) THEN
        EJ_DAY = J_DAY*1.D-2
        T_EQT = -0.04056D+0 - 0.74503D+0*EJ_DAY + 0.08823D+0*(EJ_DAY**2)
     &    + 2.0516D+0*(EJ_DAY**3) - 1.8111D+0*(EJ_DAY**4)
     &    + 0.42832D+0*(EJ_DAY**5)
      ELSE
        EJ_DAY = (J_DAY-180)*1.D-2
        T_EQT = -0.05039D+0 - 0.33954D+0*EJ_DAY + 0.04084D+0*(EJ_DAY**2)
     &    + 1.8928D+0*(EJ_DAY**3) - 1.7619D+0*(EJ_DAY**4)
     &    + 0.4224D+0*(EJ_DAY**5)
      ENDIF
!
!---  Time of solar noon, hr
!     ATMC(3) longitude, radians
!     ATMC(4) latitude, radians
!     ATMC(5) meridian, radians  ---
!
      T_SN = 1.2D+1 + (ATMC(3)-ATMC(5))*1.2D+1/GPI - T_EQT
!
!---  Time of day, hr  ---
!
      T_DAY = MOD(TM+ATMST,8.64D+4)/3.6D+3
!
!---  Solar time angle, radians  ---
!
      S_TANG = (T_DAY-T_SN)*GPI/1.2D+1
!
!---  Solar declination, radians  ---
!
      REALX = REAL(J_DAY-172)
      S_DEC = 4.101D-1*COS(2.D+0*GPI*REALX/3.65D+2)
      SJ_DAY = J_DAY*1.D-2
      S_DEC = -0.37726D+0 - 0.10564D+0*SJ_DAY + 1.2458D+0*(SJ_DAY**2)
     &  - 0.75478*(SJ_DAY**3) + 0.13627*(SJ_DAY**4)
     &  - 0.00572*(SJ_DAY**5)
      S_DEC = ASIN(S_DEC)
!
!---  Sun altitude, radians  ---
!
      S_ANG = ASIN( SIN(S_DEC)*SIN(ATMC(4)) +
     &  COS(S_DEC)*COS(ATMC(4))*COS(S_TANG) )
!
!---  Sun altitude, deg  ---
!
      DS_ANG = S_ANG*1.8D+2/GPI
!
!---  Solar zenith angle, deg  ---
!
      DZ_ANG = MIN( MAX( 9.D+1-DS_ANG,0.D+0 ),9.D+1 )
!
!---  Sun angle less than 10 deg, skip calculations ---
!
      IF( DS_ANG.LT.1.D+1 ) GOTO 100
!
!---  Incident clear sky irradiance ---
!     RN_EX - extra terrestrial radiation on a horizontal surface, W/m^2
!
      RN_EX = 1367.D+0*SIN(S_ANG)
      RN_A = (0.79D+0 - 3.75D+0/DS_ANG)*RN_EX
      F_CC = 1.088D+0*(MAX( (1.D+0 - (RN/RN_A)),0.D+0 )**(0.294D+0))
      F_CC = MIN( F_CC,1.D+0 )
      F_CC = MAX( F_CC,0.D+0 )
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CLOUD group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CONNMAP_SFC
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
!     Surface cover connection map
!
!     ICM_SFC( connection number,surface cover number) = node number
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 November 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
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
      SUB_LOG(ISUB_LOG) = '/CONNMAP_SFC'
!
!---  Find maximum plant root depth  ---
!
      PRDMX = 0.D+0
      DO IP = 1,NPLANT
        PRDMX = MAX( PRDMX,PARMS_P(1,IP) )
      ENDDO
!
!---  Loop over number of surface cover areas  ---
!
      DO NSCA = 1,NSFCA
!
!---    Loop over surface cover nodes  ---
!
        DO NSCN = 1,NSFCN
          NC = 1
          PRDX = PRDMX
          N = ICM_SFC(NC,NSCN)
!
!---      Continuous loop from surface cover node down field nodes  ---
!
          DO
!
!---        Inactive node, exit  ---
!
            IF( IXP(N).EQ.0 ) EXIT
!
!---        Roots continue through to next deeper field node  ---
!
            IF( PRDX.GT.DZGF(N) ) THEN
              NC = NC + 1
              PRDX = PRDX - DZGF(N)
              N = N-IJFLD
!
!---          Bottom of domain, exit  ---
!
              IF( N.LE.0 ) EXIT
!
!---          Inactive node, exit  ---
!
              IF( IXP(N).EQ.0 ) EXIT
              ICM_SFC(NC,NSCN) = N
            ELSE
!
!---          Maximum root depth ends within node  ---
!
              EXIT
            ENDIF
          ENDDO
          NSFCC = MAX( NC,NSFCC )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CONNMAP_SFC group ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DEW_PT( RH,TA,T_DPX )
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
!     Dew point
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 October 2003.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 PSX(2),FX(2)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DEW_PT'
!
!---  Guess dew-point temperature  ---
!
      T_DPX = TA
      DT_DP = -1.D-7
      INDX = 0
      CALL REGION_4( T_DPX,PSAX,INDX )
!
!---  Newton iteration to compute dew-point temperature  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      DO 200 M = 1,2
        T_DPX = T_DPX
        IF( M.EQ.2 ) T_DPX = T_DPX + DT_DP
        INDX = 0
        CALL REGION_4( T_DPX,PSX(M),INDX )
        FX(M) = PSAX*RH - PSX(M)
  200 CONTINUE
      DFX = (FX(2)-FX(1))/DT_DP
      DT_DPX = -FX(1)/DFX
      T_DPX = T_DPX + DT_DPX
!
!---  No convergence on phase compositions  ---
!
      IF( NC.GT.32 ) THEN
        INDX = 14
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'DEW_PT'
        CHMSGX = 'Unconverged Dew Point Temperature'
        RLMSGX = T_DPX
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
      IF( ABS(DT_DPX).GT.1.D-6 ) GOTO 100
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DEW_PT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLUX_SFC
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
!     RE_GS - Residual of the energy balance at the ground surface,
!             ignoring dissolved air mass in aqueous,
!             W/m^2 ground
!             Dimensioned (5+ISVC,NSCN)
!     RW_GS - Residual of the water balance at the ground surface,
!             kg/s m^2 ground
!             Dimensioned (5+ISVC,NSCN)
!     RE_PL - Residual of the energy balance at the plant,
!             W/m^2 ground
!             Dimensioned (5+(NSFCC*ISVC),NSCN)
!     RE_CP - Residual of the energy balance at the canopy height,
!             W/m^2 ground
!             Dimensioned (5+(NSFCC*ISVC),NSCN)
!     RW_CP - Residual of the water balance at the canopy height,
!             kg/s m^2 ground
!             Dimensioned (5+(NSFCC*ISVC),NSCN)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 November 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE PLT_ATM
      USE JACOB
      USE HYST
      USE GRID
      USE FDVT
      USE FDVP
      USE FDVG
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      EXTERNAL VRUGT
      REAL(KIND=DP) :: L_SAPX,L_SAPMX,L_NSX,L_SPX
      REAL(KIND=DP) :: RFIM_PX(2)
      REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: PAIX,LAIX
      REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: RN_SNPX,RN_LNPX
      REAL(KIND=DP) :: D
      REAL(KIND=DP) :: PA,RH,RN,TA,UZ,RF,T_DPX,DZ_ANG,F_CC,TA_KX
      REAL(KIND=DP) :: EM_CSX,SIGMAX
      REAL(KIND=DP) :: RN_LDX
      REAL(KIND=DP) :: CPG_AX,E_AX,ED_DC,HGA_AX,HGW_AX,HL_AX,
     &  K_VKX,KVIS_AX,PGA_AX,PGW_AX,PR_AX,PSW_AX,
     &  RHO_AX,RHOGA_AX,RHOGW_AX,THD_AX,THK_AX,UGA_AX,UGW_AX,VIS_AX,
     &  XGA_AX,XGW_AX,XLA_AX,XLW_AX,XMLA_AX,XMLW_AX
      REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: L_PAIX,L_SAX,
     &  RN_SDSX,RN_LDSX,
     &  RA_SCX,RA_CAX,RA_BSX,PAI_PX
      REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: C_COEFX,
     &  L_PAX,RA_PCX,RFIC_PX,RS_PCX,RWU_IX,T_PCX
      REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:) :: RWU_PX
      REAL(KIND=DP) :: NU_MCF
      INTEGER :: J_DAY
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPLANT
      INTEGER, DIMENSION(1:10) :: M_ND,M_GS
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Data Statements---------------------------------!
!
      DATA K_VKX,ED_DC / 4.1D-1,2.5D+0 /
      DATA M_ND / 2,2,2,2,2,2,3,4,5,6 /
      DATA M_GS / 2,3,4,5,6,7,2,2,2,2 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLUX_SFC'
!
!---  Allocate saved arrays ---
!
      IF( .NOT. ALLOCATED(PAIX) ) THEN
        ALLOCATE( PAIX(1:LPLANT,1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: PAIX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(LAIX) ) THEN
        ALLOCATE( LAIX(1:LPLANT,1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: LAIX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RN_SNPX) ) THEN
        ALLOCATE( RN_SNPX(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RN_SNPX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RN_LNPX) ) THEN
        ALLOCATE( RN_LNPX(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RN_LNPX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(L_PAIX) ) THEN
        ALLOCATE( L_PAIX(1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: L_PAIX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(L_SAX) ) THEN
        ALLOCATE( L_SAX(1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: L_SAX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RN_SDSX) ) THEN
        ALLOCATE( RN_SDSX(1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RN_SDSX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RN_LDSX) ) THEN
        ALLOCATE( RN_LDSX(1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RN_LDSX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RA_SCX) ) THEN
        ALLOCATE( RA_SCX(1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RA_SCX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RA_CAX) ) THEN
        ALLOCATE( RA_CAX(1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RA_CAX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RA_BSX) ) THEN
        ALLOCATE( RA_BSX(1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RA_BSX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(PAI_PX) ) THEN
        ALLOCATE( PAI_PX(1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: PAI_PX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(C_COEFX) ) THEN
        ALLOCATE( C_COEFX(1:LPLANT,1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: C_COEFX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(L_PAX) ) THEN
        ALLOCATE( L_PAX(1:LPLANT,1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: L_PAX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RA_PCX) ) THEN
        ALLOCATE( RA_PCX(1:LPLANT,1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RA_PCX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RFIC_PX) ) THEN
        ALLOCATE( RFIC_PX(1:LPLANT,1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RFIC_PX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RS_PCX) ) THEN
        ALLOCATE( RS_PCX(1:LPLANT,1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RS_PCX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RWU_IX) ) THEN
        ALLOCATE( RWU_IX(1:LPLANT,1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RWU_IX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(T_PCX) ) THEN
        ALLOCATE( T_PCX(6+LUK,1:LSFCN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: T_PCX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(RWU_PX) ) THEN
        ALLOCATE( RWU_PX(1:LSV,1:LPLANT,1:LSFCN*LSFCC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: RWU_PX'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      IF( .NOT. ALLOCATED(IPLANT) ) THEN
        ALLOCATE( IPLANT(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'FLUX_SFC'
          CHMSGX = 'Allocation Error: IPLANT'
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
!
!---  Initialize variables  ---
!
      DO IP = 1,LPLANT
        RN_LNPX(IP) = 0.D+0
        IPLANT(IP) = 0
      ENDDO
      L1: DO NSCN = 1,NSFCN
        L_PAIX(NSCN) = 0.D+0
        L_SAX(NSCN) = 0.D+0
        RN_SDSX(NSCN) = 0.D+0
        RN_LDSX(NSCN) = 0.D+0
        RA_SCX(NSCN) = 0.D+0
        RA_CAX(NSCN) = 0.D+0
        RA_BSX(NSCN) = 0.D+0
        PAI_PX(NSCN) = 0.D+0
!
!---    Ground-surface equation residuals  ---
!
        DO M = 1,6+LUK
          RE_GS(M,NSCN) = 0.D+0
          RW_GS(M,NSCN) = 0.D+0
        ENDDO
!
!---    Plant and canopy equation residuals  ---
!
        DO M = 1,LUK*LSFCC+6
          RE_PL(M,NSCN) = 0.D+0
          RE_CP(M,NSCN) = 0.D+0
          RW_CP(M,NSCN) = 0.D+0
        ENDDO
!
!---    Initialize leaf-area index and plant-area index ---
!
        DO IP = 1,LPLANT
          LAIX(IP,NSCN) = 0.D+0
          PAIX(IP,NSCN) = 0.D+0
          C_COEFX(IP,NSCN) = 0.D+0
          L_PAX(IP,NSCN) = 0.D+0
          RA_PCX(IP,NSCN) = 0.D+0
          RFIC_PX(IP,NSCN) = 0.D+0
          RS_PCX(IP,NSCN) = 0.D+0
          RWU_IX(IP,NSCN) = 0.D+0
          DO NC = 1,LSFCN*LSFCC
            DO M = 1,LSV
              RWU_PX(M,IP,NC) = 0.D+0
            ENDDO
          ENDDO
        ENDDO
!
!---    Transpiration for all plants on surface cover node  ---
!
        DO M = 1,6+LUK
          T_PCX(M,NSCN) = 0.D+0
        ENDDO
!
!---  Bottom of surface-cover node loop  ---
!
      ENDDO L1
!
!---  Time-only dependent calculations, only execute during the
!     first iteration of a time step  ---
!
      IF( NITER.LE.1 ) THEN
!
!---    Time of year  ---
!
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        D = TMZ + ATMST
        MYR = 1900
        SEC_YR = 3.1536D+7
        SEC_DAY = 8.64D+4
        L2: DO
          MYR = MYR+1
          D_X = D
          D = D - SEC_YR
          NDPYR = 365
          IF( MOD(MYR,1000).EQ.0 ) THEN
            IF( MOD(MYR,400).EQ.0 ) THEN
              D = D - SEC_DAY
              NDPYR = 366
            ENDIF
          ELSEIF( MOD(MYR,4).EQ.0 ) THEN
            D = D - SEC_DAY
            NDPYR = 366
          ENDIF
          IF( D.LE.0.D+0 ) EXIT
        ENDDO L2
        D = MOD( D_X,SEC_YR )
        J_DAY = INT(D/SEC_DAY) + 1
        J_DAY = MOD( J_DAY,NDPYR )
!
!---    Set atmospheric conditions ---
!
        CALL ATMOS_C( PA,RH,RN,TA,UZ,RF,IERR )
!
!---    Simulation time outside of atmospheric conditions time limit
!       skip Shuttleworth-Wallace calculations  ---
!
        IF( IERR.EQ.1 ) THEN
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
!
!---    Atmospheric properties  ---
!
!       E_AX - water-vapor density, kg/m^3
!       PSW_AX - saturated water-vapor pressure, Pa
!       PGW_AX - atmospheric water-vapor pressure, Pa
!       PGA_AX - atmospheric air partial pressure, Pa
!       CPG_AX - atmospheric specific heat, J/kg
!       RHOGW_AX - atmospheric water-vapor density, kg/m^3
!       RHOGA_AX - atmospheric air density, kg/m^3
!       RHOG_AX - atmospheric gas density, kg/m^3
!       HGW_AX - atmospheric water-vapor enthalpy, J/kg
!       UGW_AX - atmospheric water-vapor internal energy, J/kg
!       HGA_AX - atmospheric air enthalpy, J/kg
!       UGA_AX - atmospheric air internal energy, J/kg
!       RHO_AX - dry-air density, kg/m^3
!       VIS_AX - dry-air viscosity, Pa s
!       KVIS_AX - dry-air kinematic viscosity, m^2/s
!       CP_AX - dry-air specific heat, J/kg
!       THK_AX - dry-air thermal conductivity, W/m K
!       THD_AX - dry-air thermal diffusivity, m^2/s
!       PR_AX - dry-air Prandtl number
!
        INDX = 0
        CALL REGION_4( TA,PSW_AX,INDX )
        PGW_AX = RH*PSW_AX
        CALL P_IAPWS( TA,PGW_AX,RHOGW_AX,RHOX,HGW_AX,HX,UGW_AX,UX )
        E_AX = RHOGW_AX
        PGA_AX = MAX( PA-PGW_AX,0.D+0 )
        CALL AIRGSD( TA,PGA_AX,RHOGA_AX )
        RHOG_AX = RHOGW_AX+RHOGA_AX
        XGW_AX = RHOGW_AX/RHOG_AX
        XGA_AX = RHOGA_AX/RHOG_AX
        CALL AIRGSH( TA,PGA_AX,HGA_AX,UGA_AX )
        HG_AX = XGW_AX*HGW_AX + XGA_AX*HGA_AX
        TAIX = TA + 1.D-6
        CALL REGION_4( TAIX,PSW_AI,INDX )
        PGW_AI = RH*PSW_AI
        CALL P_IAPWS( TAIX,PGW_AX,RHOGW_AI,RHOX,HGW_AI,HX,UGW_AI,UX )
        CALL AIRGSD( TAIX,PGA_AX,RHOGA_AI )
        RHOG_AI = RHOGW_AI+RHOGA_AI
        XGW_AI = RHOGW_AI/RHOG_AI
        XGA_AI = RHOGA_AI/RHOG_AI
        CALL AIRGSH( TAIX,PGA_AX,HGA_AI,UGA_AI )
        HG_AI = XGW_AI*HGW_AI + XGA_AI*HGA_AI
        CPG_AX = (HG_AI-HG_AX)/1.D-6
        CALL AIRGSD( TA,PA,RHO_AX )
        CALL AIRGSV( TA,VIS_AX )
        KVIS_AX = VIS_AX/RHO_AX
        CALL AIRGSC( TA,CP_AX )
        CALL AIRGSK( TA,THK_AX )
        THD_AX = THK_AX/(RHO_AX*CP_AX)
        PR_AX = KVIS_AX/THD_AX
!
!---    Precipitation properties  ---
!
        XMLA_AX = PGA_AX/HCAW
        XMLW_AX = MAX( 1.D+0-XMLA_AX,0.D+0 )
        XLW_AX = XMLW_AX*WTMW/(XMLA_AX*WTMA + XMLW_AX*WTMW)
        XLA_AX = MAX( 1.D+0-XLW_AX,0.D+0 )
        TA_Y = MAX( TA,1.D-2 )
        CALL P_IAPWS( TA_Y,PA,RHOX,RHOL_AX,HX,HLW_AX,UX,UX )
        HL_AX = XLW_AX*HLW_AX + XLA_AX*HGA_AX
!
!---    Downward long-wave radiation,
!       W/m^2 ground surface  ---
!
!       RN_LDX - downward long-wave radiation, W/m^2
!       EM_CSX - clear-sky emissivity
!       F_CC - fractional cloud cover
!       SIGMAX - Stefan-Boltzmann constant, W/m^2 K^4
!       TA_KX - atmospheric temperature, K
!       T_DPX - atmospheric dew-point temperature, C
!       DZ_ANG - solar zenith angle, deg
!
        CALL DEW_PT( RH,TA,T_DPX )
        CALL CLOUD( RN,F_CC,J_DAY,DZ_ANG )
!
!---    Clear sky emissivity from dew point temperature, C
!       from the Berkeley equation (Berdahl and Fromberg)  ---
!
        EM_CSX = MIN( 0.741D+0 + 6.2D-3*T_DPX,1.D+0 )
        SIGMAX = 5.67D-8
        TA_KX = TA + TABS
        RN_LDX = EM_CSX*(1.D+0-0.84D+0*F_CC)*SIGMAX*(TA_KX**4)
     &    + 0.84D+0*F_CC*SIGMAX*(TA_KX**4)
!
!---    Loop over surface cover nodes  ---
!
        L3: DO NSCN = 1,NSFCN
!
!---      Surface cover area index and surface cover type index  ---
!
          NSCA = ID_SFC(NSCN)
          NSCT = NSFCT(NSCA)
!
!---      Leaf area index and plant area index as a function of time   ---
!
          IF( TMZ.LE.PARMS_SFC(1,1,NSCA) ) THEN
            DO IP = 1,NPLANT
              IPX = (IP-1)*2
              LAIX(IP,NSCN) = PARMS_SFC(IPX+2,1,NSCA)
              PAIX(IP,NSCN) = PARMS_SFC(IPX+3,1,NSCA)
            ENDDO
          ELSEIF( TMZ.GE.PARMS_SFC(1,NSCT,NSCA) ) THEN
            DO IP = 1,NPLANT
              IPX = (IP-1)*2
              LAIX(IP,NSCN) = PARMS_SFC(IPX+2,NSCT,NSCA)
              PAIX(IP,NSCN) = PARMS_SFC(IPX+3,NSCT,NSCA)
            ENDDO
          ELSE
            DO M = 2,NSFCT(NSCA)
              IF( TMZ.LE.PARMS_SFC(1,M,NSCA) ) THEN
                TDX = PARMS_SFC(1,M,NSCA)-PARMS_SFC(1,M-1,NSCA)
                DTX = MIN( PARMS_SFC(1,M,NSCA)-TMZ,DT )
                TFX = (TMZ-PARMS_SFC(1,M-1,NSCA))/TDX
                DO IP = 1,NPLANT
                 IPX = (IP-1)*2
                 LAIX(IP,NSCN) = PARMS_SFC(IPX+2,M-1,NSCA) + TFX*
     &             (PARMS_SFC(IPX+2,M,NSCA)-PARMS_SFC(IPX+2,M-1,NSCA))
                 PAIX(IP,NSCN) = PARMS_SFC(IPX+3,M-1,NSCA) + TFX*
     &             (PARMS_SFC(IPX+3,M,NSCA)-PARMS_SFC(IPX+3,M-1,NSCA))
                ENDDO
                EXIT
              ENDIF
            ENDDO
          ENDIF
          DO IP = 1,NPLANT
            IF( PAIX(IP,NSCN)*LAIX(IP,NSCN).LT.EPSL ) THEN
              IPLANT(IP) = 0
            ELSE
              IPLANT(IP) = 1
            ENDIF
          ENDDO
!
!---      Net short-wave radiation at ground surface, considering
!         interception by plants, W/m^2 ground  ---
!
!         RN - downward, short-wave radiation, W/m^2 ground
!         RN_SDSX(NSCN) - downward short-wave radiation
!           at ground surface, W/m^2 ground
!         RN_SUSX - upward short-wave radiation
!           at ground surface, W/m^2 ground
!         RN_SNSX - net short-wave radiation
!           at ground surface, W/m^2 ground
!         LAIX(IP,NSCN) - leaf-area index of plant specie IP,
!           m^2 leaf/m^2 plant ground
!         PAIX(IP,NSCN) - plant-area index of plant specie IP,
!           m^2 plant ground/m^2 ground
!         ALB_SX - short-wave albedo of ground surface
!
          PAI_PX(NSCN) = 0.D+0
          RN_SDSX(NSCN) = 0.D+0
          DO IP = 1,NPLANT
            IF( IPLANT(IP).EQ.0 ) CYCLE
            PAI_PX(NSCN) = PAI_PX(NSCN) + PAIX(IP,NSCN)
            RN_SDSX(NSCN) = RN_SDSX(NSCN) + PAIX(IP,NSCN)*
     &        EXP(-7.D-1*LAIX(IP,NSCN))
          ENDDO
          RN_SDSX(NSCN) = RN*((1.D+0-PAI_PX(NSCN)) + RN_SDSX(NSCN))
!
!---      Downward long-wave radiation to the ground surface,
!         excluding emission from plants, W/m^2 ground surface  ---
!
!         RN_LDX - downward long-wave radiation, W/m^2
!         RN_LDSX(NSCN) - downward long-wave radiation to 
!           ground surface, W/m^2
!         PAIX(IP,NSCN) - plant-area index of plant specie IP
!         LAIX(IP,NSCN) - leaf-area index of plant specie IP
!
          RN_LDSX(NSCN) = 0.D+0
          DO IP = 1,NPLANT
            IF( IPLANT(IP).EQ.0 ) CYCLE
            RN_LDSX(NSCN) = RN_LDSX(NSCN) + PAIX(IP,NSCN)*
     &        EXP(-7.D-1*LAIX(IP,NSCN))
          ENDDO
          RN_LDSX(NSCN) = RN_LDX*((1.D+0-PAI_PX(NSCN)) + RN_LDSX(NSCN))
!
!---      Eddy diffusion resistance above and below the average
!         plant height  ---
!
!         PH_AVG - average plant height, m
!         PAIX(IP,NSCN) - plant-area index of plant specie IP
!         PARMS_P(11,IP) - height of plant specie IP
!         ZP_D - zero-plane displacement, m
!         RL_PX - plant roughness length, m
!         RA_SCX(NSCN) - aerodynamic resistance between ground and
!           plant height, s/m
!         RA_CAX(NSCN) - aerodynamic resistance between mean canopy 
!           flow and reference height, s/m
!         ATMC(1) - wind speed reference height, m
!         K_VKX - von Karman constant, 0.4
!         ED_DC - eddy diffusion decay constant
!         UZ  Wind speed, m/s
!         ATMC(6) - momentum transport surface roughness length, m
!         ATMC(7) - heat and mass transport surface roughness length, m
!         RA_BSMX - maximum aerodynamic resistance between ground and
!           reference height for bare surface, s/m
!
          PH_AVG = 0.D+0
          VARX = 0.D+0
          DO IP = 1,NPLANT
            IF( IPLANT(IP).EQ.0 ) CYCLE
            PH_AVG = PH_AVG + PAIX(IP,NSCN)*PARMS_P(11,IP)
            VARX = VARX + PAIX(IP,NSCN)
          ENDDO
          IF( VARX.GT.EPSL ) THEN
            PH_AVG = PH_AVG/VARX
          ELSE
            PH_AVG = 0.D+0
          ENDIF
          IF( PH_AVG.GT.EPSL ) THEN
            ZP_D = 0.63D+0*PH_AVG
            RL_PX = 0.13D+0*PH_AVG
            RA_X = (LOG((ATMC(1)-ZP_D)/RL_PX)/((K_VKX**2)*UZ))*
     &        (PH_AVG/(ED_DC*(PH_AVG-ZP_D)))
            RA_SCX(NSCN) = RA_X*
     &        (EXP(ED_DC) - EXP(ED_DC*(1.D+0-((ZP_D+RL_PX)/PH_AVG))))
            RA_CAX(NSCN) = RA_X*
     &        ((EXP(ED_DC*(1.D+0-((ZP_D+RL_PX)/PH_AVG)))-1.D+0) +
     &        LOG((ATMC(1)-ZP_D)/(PH_AVG-ZP_D)))
            RA_BSMX = RHOG_AX*CPG_AX*ATMC(2)/THK_AX
            VLG_1X = LOG((ATMC(1)+ATMC(6))/ATMC(6))
            VLG_2X = LOG((ATMC(2)+ATMC(7))/ATMC(7))
            UZ_FX = K_VKX*UZ/VLG_1X
            RA_BSX(NSCN) = MIN( VLG_2X/(K_VKX*UZ_FX),RA_BSMX )
          ELSE
!
!---        Eddy diffusion resistance for bare surface  ---
!
!           ATMC(1) - wind speed reference height, m
!           K_VKX - von Karman constant, 0.4
!           UZ  Wind speed, m/s
!           UZ_FX  Friction velocity, m/s
!           RA_BSX(NSCN) - aerodynamic resistance between ground and
!             reference height for bare surface, s/m
!           RA_BSMX - maximum aerodynamic resistance between ground and
!             reference height for bare surface, s/m
!           ATMC(6) - momentum transport surface roughness length, m
!           ATMC(7) - heat and mass transport surface roughness 
!             length, m
!
            RA_BSMX = RHOG_AX*CPG_AX*ATMC(2)/THK_AX
            VLG_1X = LOG((ATMC(1)+ATMC(6))/ATMC(6))
            VLG_2X = LOG((ATMC(2)+ATMC(7))/ATMC(7))
            UZ_FX = K_VKX*UZ/VLG_1X
            RA_BSX(NSCN) = MIN( VLG_2X/(K_VKX*UZ_FX),RA_BSMX )
            RA_SCX(NSCN) = 5.D-1*RA_BSX(NSCN)
            RA_CAX(NSCN) = 5.D-1*RA_BSX(NSCN)
          ENDIF
!
!---      Aerodynamic resistance between the plant leaves and
!         mean canopy flow  ---
!
!         UZ_MCF - wind speed at the mean canopy flow height, m/s
!         RE_MCF - Reynolds number at the mean canopy flow height
!         PR_AX - dry-air Prandtl number
!         NU_MCF - Nusselt number at the mean canopy flow height
!         RHO_AX - dry-air density, kg/m^3
!         VIS_AX - dry-air viscosity, Pa s
!         KVIS_AX - dry-air kinematic viscosity, m^2/s
!         CPG_AX - dry-air specific heat, J/kg
!         THK_AX - dry-air thermal conductivity, W/m K
!         THD_AX - dry-air thermal diffusivity, m^2/s
!         RB_PX - leaf boundary layer resistance, s/m
!         RA_PCX - aerodynamic resistance between plant leaves and
!                 mean canopy flow, s/m
!         WL_PX - leaf width, m
!
          IF( PH_AVG.GT.EPSL ) THEN
            UZ_FX = UZ*K_VKX/LOG((ATMC(1)-ZP_D+RL_PX)/RL_PX)
            UZ_MCF = (UZ_FX/K_VKX)*LOG((PH_AVG-ZP_D+RL_PX)/RL_PX)
!
!---      Use dry air properties to compute Reynolds, Prandtl,
!         and Nusselt numbers ---
!
            DO IP = 1,NPLANT
              IF( IPLANT(IP).EQ.0 ) CYCLE
              WL_PX = 1.D-2*LAIX(IP,NSCN)
              RE_MCF = WL_PX*UZ_MCF/KVIS_AX
              NU_MCF = (2.D+0/3.D+0)*(RE_MCF**5.D-1)*(PR_AX**(THIRD))
              RB_PX = 7.D-1*WL_PX/(THD_AX*NU_MCF)
              RA_PCX(IP,NSCN) = RB_PX/LAIX(IP,NSCN)
            ENDDO
          ENDIF
!
!---      Stomatal resistance  ---
!
!         RS_PCX - stomatal resistance between plant leaf interior
!                 and exterior, s/m
!         RS_MINX - minimum stomatal resistance between plant leaf 
!                 interior and exterior, s/m
!         F_1X - influence of solar radiation
!
!         Hicks stomatal resistance model
!
!         PARMS_P(2,IP) - minimum stomatal resistance, s/m
!         PARMS_P(3,IP) - light response coefficient, W/m^2
!         PARMS_P(27,IP) - minimum temperature for stomatal opening, K
!         PARMS_P(28,IP) - maximum temperature for stomatal opening, K
!         PARMS_P(29,IP) - optimum temperature for stomatal opening, K
!
!---      Hicks stomatal resistance model  ---
!
          DO IP = 1,NPLANT
            IF( ISRM_P(IP).EQ.1 ) THEN
              RS_MINX = PARMS_P(2,IP)
              BETAX = PARMS_P(3,IP)
              TKAX = TA+TABS
              TKCX = PARMS_P(27,IP)
              TKHX = PARMS_P(28,IP)
              TKOX = PARMS_P(29,IP)
              F_1X = ((TKAX-TKCX)/(TKOX-TKCX))*
     &          (((TKHX-TKAX)/(TKHX-TKOX))**((TKHX-TKOX)/(TKOX-TKCX)))
              F_1X = MAX( F_1X,EPSL )
              IF( IPLANT(IP).EQ.0 ) THEN
                RS_PCX(IP,NSCN) = RS_MINX
              ELSE
                RS_PCX(IP,NSCN) = RS_MINX*(1.D+0 + BETAX/(RN+SMALL))/
     &            F_1X/LAIX(IP,NSCN)
              ENDIF
            ELSE
              RS_MINX = 5.D+1
              F_1X = (4.D+2+RN)/(1.4D+0*RN+EPSL)
              IF( IPLANT(IP).EQ.0 ) THEN
                RS_PCX(IP,NSCN) = RS_MINX
              ELSE
                RS_PCX(IP,NSCN) = RS_MINX*F_1X/LAIX(IP,NSCN)
              ENDIF
            ENDIF
          ENDDO
!
!---      Crop growth index  ---
!
          DO IP = 1,NPLANT
            IF( IPLANT(IP).EQ.0 ) THEN
              C_COEFX(IP,NSCN) = 0.D+0
            ELSE
              C_COEF1X = PARMS_P(17,IP)
              D1X = PARMS_P(18,IP)
              C_COEF2X = PARMS_P(19,IP)
              D2X = PARMS_P(20,IP)
              C_COEF3X = PARMS_P(21,IP)
              D3X = PARMS_P(22,IP)
              C_COEF4X = PARMS_P(23,IP)
              D4X = PARMS_P(24,IP)
              C_COEF5X = PARMS_P(25,IP)
              D5X = PARMS_P(26,IP)
              IF( D.GE.D1X .AND. D.LT.D2X )THEN
                C_COEFX(IP,NSCN) = ((C_COEF2X-C_COEF1X)/(D2X-D1X))*
     &            (D-D1X) + C_COEF1X
              ELSEIF( D.GE.D2X .AND. D.LT.D3X )THEN
                C_COEFX(IP,NSCN) = ((C_COEF3X-C_COEF2X)/(D3X-D2X))*
     &            (D-D2X) + C_COEF2X
              ELSEIF( D.GE.D3X .AND. D.LT.D4X )THEN
                C_COEFX(IP,NSCN) = ((C_COEF4X-C_COEF3X)/(D4X-D3X))*
     &            (D-D3X) + C_COEF3X
              ELSEIF( D.GE.D4X .AND. D.LT.D5X )THEN
                C_COEFX(IP,NSCN) = ((C_COEF5X-C_COEF4X)/(D5X-D4X))*
     &            (D-D4X) + C_COEF4X
              ELSE
                C_COEFX(IP,NSCN) = 0.D+0
              ENDIF
            ENDIF
          ENDDO
!
!---      Rainfall  ---
!
!         L_SAX(NSCN) - rainfall intensity on ground surface for
!           surface cover node NSCN, m^3/s m^2 ground
!         L_PAX(IP,NSCN) - rainfall intensity for plant IP, 
!           surface cover node NSCN m^3/s m^2 ground
!         L_PAIX(NSCN) - rainfall intensity on all plants for
!           surface cover node NSCN, m^3/s m^2 ground
!
          IF( RF/EPSL.GT.EPSL ) THEN
            L_SAX(NSCN) = -RF
          ELSE
            L_SAX(NSCN) = 0.D+0
          ENDIF
          L_PAIX(NSCN) = 0.D+0
!
!---      Rainfall interception by plants  ---
!
          DO IP = 1,NPLANT
            RFIC_PX(IP,NSCN) = 0.D+0
            L_PAX(IP,NSCN) = 0.D+0
            IF( IPLANT(IP).EQ.0 ) CYCLE
!
!---        Rainfall interception capacity, kg/m^2 ground  ---
!
!           PARMS_P(16,IP) - maximum dew depth for plant specie IP, m
!           DMMX_PX -  maximum dew mass, kg/m^2 ground
!           RFIC_PX(IP,NSCN) - rainfall interception capacity
!             for plant IP, surface cover node NSCN, kg/m^2 ground
!           RFI_CF - rainfall interception capacity coefficient, s/m
!
            IF( ISLC(26).EQ.1 ) THEN
              DMMX_PX = PARMS_P(16,IP)*LAIX(IP,NSCN)*PAIX(IP,NSCN)*
     &          RHOL_AX
              RFIC_PX(IP,NSCN) = DMMX_PX
            ENDIF
!
!---        Reduce rainfall intensity by plant area index,
!           m^3/s m^2 ground
!
            L_PAX(IP,NSCN) = L_SAX(NSCN)*PAIX(IP,NSCN)
!
!---        Sum over plants, incident rainfall m^3/m^3 ground  ---
!
            L_PAIX(NSCN) = L_PAIX(NSCN) + L_PAX(IP,NSCN)
          ENDDO
          L_SAX(NSCN) = L_SAX(NSCN) - L_PAIX(NSCN)
!
!---      Potential root-water uptake based on vertical 
!         root distribution  ---
!
          L4: DO IP = 1,NPLANT
            IF( IPLANT(IP).EQ.0 ) CYCLE
            ZO = 0.D+0
            ZM = PARMS_P(1,IP)
            CALL QROMB( VRUGT,ZO,ZM,RWU_IX(IP,NSCN),IERR,IP )
            IF( IERR.EQ.1 ) THEN
              INDX = 3
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'FLUX_SFC'
              CHMSGX = 'Unconverged Romberg Integration: ' //
     &          'Root-Water Uptake Integration'
              RLMSGX = 0.D+0
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
!
!---        Actual root-water uptake based on vertical 
!           root distribution and root stress model, using
!           current capillary pressure and saturation  ---
!
            N = ICM_SFC(1,NSCN)
            DO NC = 1,NSFCC
              NX = ICM_SFC(NC,NSCN)
              IF( NX.EQ.0 ) EXIT
              ZSX = MIN( ZM,ZE(5,N)-ZE(5,NX) )
              ZEX = MIN( ZM,ZE(5,N)-ZE(1,NX) )
              DO M = 2,ISVC+2
                CALL ROOTS( UPTAKE,ZSX,ZEX,IP,NX,M )
                NZ = (NSCN-1)*NSFCC + NC
                RWU_PX(M,IP,NZ) = UPTAKE
              ENDDO
              IF( (ZM-ZEX)/EPSL.LT.EPSL ) EXIT
            ENDDO
!
!---      Bottom of plant loop  ---
!
          ENDDO L4
!
!---    Bottom of surface cover node loop  ---
!
        ENDDO L3
!
!---  Time dependent calculations, only execute during the
!     first iteration of a time step  ---
!
      ENDIF
!
!---  Secondary surface cover variables  ---
!
      CALL PROP_SFC( PA )
      PGGX = PA - PATM
!
!---  Loop over surface cover nodes  ---
!
      DO NSCN = 1,NSFCN
        N = ICM_SFC(1,NSCN)
        NQZ = NSZ(N)+IJFLD
!
!---    Surface cover area  ---
!
        NSCA = ID_SFC(NSCN)
!
!---    Loop over ground-surface indices and node indices  ---
!
        DO M = 1,ISCVF
          MN = M_ND(M)
          MS = M_GS(M)
!
!---      Net short-wave radiation intercepted by ground surface,
!         W/m^2 ground surface  ---
!
!         RN_SNSX - net short-wave radiation transitted to the
!           ground surface, W/m^2 ground
!         RN_SDSX(NSCN) - downward short-wave radiation
!           at ground surface, W/m^2 ground
!         RN_SUSX - upward short-wave radiation
!           at ground surface, W/m^2 ground
!         ALB_SX - short-wave albedo of ground surface (0.1 - 0.4)
!
!---      Moisture ground-surface solar angle model  ---
!
          IF( IALB(NSCA).GT.0 ) THEN
            ALB_SX = ALBEDO(2,NSCA) + 
     &        (ALBEDO(1,NSCA)-ALBEDO(2,NSCA))*
     &        EXP(-ALBEDO(3,NSCA)*SL_GS(MS,NSCN)*PORD(2,N))
          ENDIF
!
!---      Plem and Xiu ground-surface solar angle model  ---
!
          IF( IALB(NSCA).EQ.1 ) THEN
            ALB_SX = ALB_SX + 1.D-2*(EXP(3.286D-3*(DZ_ANG**1.5D+0))
     &        -1.D+0)
!
!---      Briegleb ground-surface solar angle model  ---
!
          ELSEIF( IALB(NSCA).EQ.2 ) THEN
            Z_ANG = DZ_ANG*GPI/1.8D+2
            ALB_SX = ALB_SX + ALBEDO(4,NSCA)*
     &        (1.D+0 + ALBEDO(5,NSCA))/
     &        (1.D+0 + 2.D+0*ALBEDO(5,NSCA)*COS(Z_ANG))
!
!---      Wang ground-surface solar angle model  ---
!
          ELSEIF( IALB(NSCA).EQ.3 ) THEN
            Z_ANG = DZ_ANG*GPI/1.8D+2
            G_1X = (-7.574D-3) + (-7.0987D-2)*(Z_ANG**2) + 
     &        (3.07588D-1)*(Z_ANG**3)
            G_2X = (-1.284909D+0) + (-1.66314D-1)*(Z_ANG**2) +
     &        (4.184D-2)*(Z_ANG**3)
            G_1RX = 2.67808D-1
            G_2RX = -1.419244D+0
            ALB_SX = ALB_SX + ALBEDO(4,NSCA)*(1.D+0 +
     &        0.346D+0*(G_1X-G_1RX) + 0.063*(G_2X-G_2RX))
          ELSE
            ALB_SX = 0.25D+0
          ENDIF
          RN_SUSX = RN_SDSX(NSCN)*ALB_SX
          RN_SNSX = (RN_SDSX(NSCN)-RN_SUSX)
!
!---      Net short-wave radiation at plants, considering
!         single reflection from ground surface, W/m^2 ground  ---
!
!         RN_SNPX(IP) - net short-wave radiation into plant
!           specie IP, W/m^2 ground
!         ALB_PX - short-wave albedo of plant specie IP
!         LAIX(IP,NSCN) - leaf-area index of plant specie IP,
!           m^2 leaf/m^2 plant ground
!         PAIX(IP,NSCN) - plant-area index of plant specie IP,
!           m^2 plant ground/m^2 ground
!         PAI_PX(NSCN) - plant-area index of all plant species
!
          DO IP = 1,NPLANT
            RN_SNPX(IP) = 0.D+0
            IF( IPLANT(IP).EQ.0 ) CYCLE
!
!---        Temporal plant solar albedo  ---
!
            IF( IALB_P(IP).EQ.1 ) THEN
              ALB_P1X = PARMS_P(5,IP)
              D1X = PARMS_P(18,IP)
              ALB_P2X = PARMS_P(6,IP)
              D2X = PARMS_P(20,IP)
              ALB_P3X = PARMS_P(8,IP)
              D3X = PARMS_P(22,IP)
              ALB_P4X = PARMS_P(9,IP)
              D4X = PARMS_P(24,IP)
              ALB_P5X = PARMS_P(10,IP)
              D5X = PARMS_P(26,IP)
              IF( D.GE.D1X .AND. D.LT.D2X )THEN
                ALB_PX = ((ALB_P2X-ALB_P1X)/(D2X-D1X))*(D-D1X)+ALB_P1X
              ELSEIF( D.GE.D2X .AND. D.LT.D3X )THEN
                ALB_PX = ((ALB_P3X-ALB_P2X)/(D3X-D2X))*(D-D2X)+ALB_P2X
              ELSEIF( D.GE.D3X .AND. D.LT.D4X )THEN
                ALB_PX = ((ALB_P4X-ALB_P3X)/(D4X-D3X))*(D-D3X)+ALB_P3X
              ELSEIF( D.GE.D4X .AND. D.LT.D5X )THEN
                ALB_PX = ((ALB_P5X-ALB_P4X)/(D5X-D4X))*(D-D4X)+ALB_P4X
              ELSE
                ALB_PX = PARMS_P(10,IP)
              ENDIF
!
!---        Constant plant solar albedo  ---
!
            ELSE
              ALB_PX = PARMS_P(10,IP)
            ENDIF
            RN_SNPX(IP) = (RN+RN_SUSX)*PAIX(IP,NSCN)*
     &        (1.D+0-EXP(-7.D-1*LAIX(IP,NSCN)))*(1.D+0-ALB_PX)
          ENDDO
!
!---      Aqueous volumetric flux from the node to the
!         ground surface, m/s  ---
!
          D_NSX = 5.D-1*DZGF(N)
          HDL_SX = PL_GS(MS,NSCN) + 
     &      GRVZ(NQZ)*D_NSX*5.D-1*(RHOL(MN,N)+RHOL_GS(MS,NSCN))
          HDL_NX = PL(MN,N)
          HDL_NSX = HDL_NX - HDL_SX
          INDX = -1
          VISL_NSX = DIFMN(VISL(MN,N),VISL_GS(MS,NSCN),D_NSX,D_NSX,
     &      HDL_NSX,INDX)
          INDX = -4
          RKL_NSX = DIFMN(RKL(3,MN,N),RKL_GS(MS,NSCN),D_NSX,D_NSX,
     &      HDL_NSX,INDX)
          L_NSX = (RKL_NSX*PERM(3,IZ(N))/VISL_NSX)*(HDL_NSX/D_NSX)
!
!---      Gas volumetric flux from the node to the ground
!         surface, m/s  ---
!
          HDG_SX = PGGX + 
     &      GRVZ(NQZ)*D_NSX*5.D-1*(RHOG(MN,N)+RHOG_GS(MS,NSCN))
          HDG_NX = PG(MN,N)
          HDG_NSX = HDG_NX - HDG_SX
          INDX = -1
          VISG_NSX = DIFMN(VISG(MN,N),VISG_GS(MS,NSCN),D_NSX,D_NSX,
     &      HDG_NSX,INDX)
          INDX = -4
          RKG_NSX = DIFMN(RKG(MN,N),RKG_GS(MS,NSCN),D_NSX,D_NSX,
     &      HDG_NSX,INDX)
          G_NSX = (RKG_NSX*PERM(3,IZ(N))*HDG_NSX)/(VISG_NSX*D_NSX)
!
!---      Diffusive water-vapor mass flux from the node to
!         the surface, kg/m^2 s  ---
!
          DF_SX = TORG_GS(MS,NSCN)*PORD(2,N)*SG_GS(MS,NSCN)*
     &      DFGW_GS(MS,NSCN)
          DF_NX = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGW(MN,N)
          INDX = -1
          DF_NSX = DIFMN(DF_NX,DF_SX,D_NSX,D_NSX,G_NSX,INDX)
          E_SX = RHOG_GS(MS,NSCN)*XGW_GS(MS,NSCN)
          E_NX = RHOG(MN,N)*XGW(MN,N)
          E_NSX = DF_NSX*(E_NX-E_SX)/D_NSX
!
!---      Successive substitution, following Campbell [1985],
!         to determine aerodynamic resistance  ---
!
          H_SAX = (1.D+0-PAI_PX(NSCN))*(T_GS(MS,NSCN)-TA)*RHOG_AX*
     &      CPG_AX/RA_BSX(NSCN)
          VARX = K_VKX*ATMC(2)*GRAV/(RHOG_AX*CPG_AX*TA_KX)
          ASPX = -VARX*H_SAX/(UZ_FX**3)
          DO NCX = 1,4
!
!---        Surface temperature higher than atmospheric
!           temperature, unstable conditions  ---
!
            IF( (T_GS(MS,NSCN)-TA).GT.EPSL ) THEN
              SPCFHX = -2.D+0*LOG(5.D-1*(1.D+0+SQRT(
     &          MAX(1.D+0-(1.6D+1*ASPX),0.D+0))))
              SPCFMX = 6.D-1*SPCFHX
!
!---        Surface temperature lower than atmospheric
!           temperature, stable conditions  ---
!
            ELSEIF( (T_GS(MS,NSCN)-TA).LT.-EPSL ) THEN
              SPCFHX = 4.7D+0*ASPX
              SPCFMX = 4.7D+0*ASPX
!
!---        Surface temperature equal to atmospheric
!           temperature, skip calculations  ---
!
            ELSE
              SPCFHX = 0.D+0
              SPCFMX = 0.D+0
            ENDIF
!
!---        Corrections limited to 20%  ---
!
            SPCFHX = SIGN( MIN( ABS(SPCFHX),2.D-1*VLG_2X ),SPCFHX )
            SPCFMX = SIGN( MIN( ABS(SPCFMX),2.D-1*VLG_1X ),SPCFMX )
            UZ_FX = K_VKX*UZ/(VLG_1X+SPCFMX)
            RA_BSX(NSCN) = MIN( (VLG_2X+SPCFHX)/(K_VKX*UZ_FX),RA_BSMX )
!
!---        Sensible heat flux from the ground surface to the
!           atmosphere, W/m^2  ---
!
            H_SAX = (1.D+0-PAI_PX(NSCN))*(T_GS(MS,NSCN)-TA)*
     &        RHOG_AX*CPG_AX/RA_BSX(NSCN)
!
!---        Atmospheric stability parameter  ---
!
            ASPX = -VARX*H_SAX/(UZ_FX**3)
          ENDDO
!
!---      Saturated water-vapor density at ground surface
!         temperature, kg/m^3  ---
!
          INDX = 0
          CALL REGION_4( T_GS(MS,NSCN),PSWX,INDX )
          CALL P_IAPWS( T_GS(MS,NSCN),PSWX,PE_SX,RHOX,HGWX,HX,UGWX,UX )
!
!---      Diffusive water-vapor mass flux from the ground
!         surface to the atmosphere, kg/m^2 s  ---
!
          E_SAX = (1.D+0-PAI_PX(NSCN))*(E_SX-E_AX)/RA_BSX(NSCN)
          PE_SAX = (1.D+0-PAI_PX(NSCN))*(PE_SX-E_AX)/RA_BSX(NSCN)
!
!---      Enthalpy of water-vapor, air and gas at the surface, J/kg  ---
!
          HG_SX = HG_GS(MS,NSCN)
          UG_SX = UEG_GS(MS,NSCN)
!
!---      Enthalpy of liquid-water and aqueous at the
!         ground surface, J/kg  ---
!
          HL_SX = HL_GS(MS,NSCN)
!
!---      Enthalpy of gas at the node, J/kg  ---
!
          HG_NX = HG(MN,N)
!
!---      Upward long-wave radiation from the ground surface,
!         including reflected long-wave radiation from the sky
!         and emission from the ground surface, W/m^2 ground  ---
!
          TS_KX = T_GS(MS,NSCN) + TABS
          EM_SX = 9.D-1 + 1.8D-1*(SL_GS(MS,NSCN)*PORD(2,N))
          RN_LUSX = (1.D+0-EM_SX)*RN_LDSX(NSCN) + 
     &      EM_SX*SIGMAX*(TS_KX**4)
!
!---      Net long-wave radiation into the ground surface,
!         W/m^2 ground  ---
!
          RN_LNSX = RN_LDSX(NSCN) - RN_LUSX
!
!---      Net long-wave radiation into the plants and
!         net downward long-wave radiation from plants to
!         ground surface, W/m^2 ground  ---
!
!         EM_PX - Emissivity of plants [Llasat and Snyder, 1998]
!
          DO IP = 1,NPLANT
            RN_LNPX(IP) = 0.D+0
            IF( IPLANT(IP).EQ.0 ) CYCLE
            TP_KX = T_PL(MS,NSCN) + TABS
            EM_PX = 9.8D-1
            RN_LNPX(IP) = PAIX(IP,NSCN)*(RN_LDX + RN_LUSX - 
     &        2.D+0*SIGMAX*(TP_KX**4))*EM_PX*
     &        (1.D+0-EXP(-7.D-1*LAIX(IP,NSCN)))
            RN_LNSX = RN_LNSX + PAIX(IP,NSCN)*SIGMAX*
     &        (TP_KX**4)*EM_PX*(1.D+0-EXP(-7.D-1*LAIX(IP,NSCN)))*EM_SX
          ENDDO
!
!---      Net long- and short-wave radiation into the
!         ground surface, W/m^2 ground  ---
!
          RN_NSX = RN_LNSX + RN_SNSX
!
!---      Net long- and short-wave radiation into the
!         plants, W/m^2 ground  ---
!
          RN_NPX = 0.D+0
          DO IP = 1,NPLANT
            IF( IPLANT(IP).EQ.0 ) CYCLE
            RN_NPX = RN_NPX + RN_LNPX(IP) + RN_SNPX(IP)
          ENDDO
!
!---      Sensible heat flux from the field node to
!         the ground surface, W/m^2  ---
!
          INDX = -1
          TK_NSX = DIFMN(THKL(MN,N),THKL_GS(MS,NSCN),D_NSX,D_NSX,
     &      ZERO,INDX)
          H_NSX = (T(MN,N)-T_GS(MS,NSCN))*TK_NSX/D_NSX
!
!---      Plant condensate density and enthalpy  ---
!
          PG_PX = PA
          INDX = 0
          TPX = MAX( T_PL(MS,NSCN),1.D-2 )
          CALL REGION_4( TPX,PSW_PX,INDX )
          CALL P_IAPWS( TPX,PSW_PX,RHOX,RHOL_PX,HX,HLW_PL(MS,NSCN),
     &      UX,ULW_PX )
          PGW_PX = PSW_PX
          PGA_PX = MAX( PG_PX-PGW_PX,0.D+0 )
          XMLA_PX = PGA_PX/HCAW
          XMLW_PX = MAX( 1.D+0-XMLA_PX,0.D+0 )
          XLW_PX = XMLW_PX*WTMW/(XMLA_PX*WTMA + XMLW_PX*WTMW)
          XLA_PX = MAX( 1.D+0-XLW_PX,0.D+0 )
          RHOLW_PX = RHOL_PX*XLW_PX
          CALL AIRGSH( TPX,PGA_PX,HGA_PX,UGA_PX )
          HL_PX = XLW_PX*HLW_PL(MS,NSCN) + XLA_PX*HGA_PX
!
!---      Component and phase density interfacial averages  ---
!
          RHOGA_NX = RHOG(MN,N)*XGA(MN,N)
          RHOGW_NX = RHOG(MN,N)*XGW(MN,N)
          RHOLW_NX = RHOL(MN,N)*XLW(MN,N)
          RHOGW_SX = RHOG_GS(MS,NSCN)*XGW_GS(MS,NSCN)
          RHOGA_SX = RHOG_GS(MS,NSCN)*XGA_GS(MS,NSCN)
          RHOLW_SX = RHOL_GS(MS,NSCN)*XLW_GS(MS,NSCN)
          RHOLW_AX = RHOL_AX*XLW_AX
          INDX = -4
          RHOGW_SAX = DIFMN(RHOGW_SX,RHOGW_AX,D_NSX,D_NSX,G_NSX,INDX)
          RHOGW_NSX = DIFMN(RHOGW_NX,RHOGW_SX,D_NSX,D_NSX,G_NSX,INDX)
          RHOGA_SAX = DIFMN(RHOGA_SX,RHOGA_AX,D_NSX,D_NSX,G_NSX,INDX)
          RHOGA_NSX = DIFMN(RHOGA_NX,RHOGA_SX,D_NSX,D_NSX,G_NSX,INDX)
          RHOG_SAX = DIFMN(RHOG_GS(MS,NSCN),RHOG_AX,D_NSX,D_NSX,
     &      G_NSX,INDX)
          RHOG_NSX = DIFMN(RHOG(MN,N),RHOG_GS(MS,NSCN),D_NSX,D_NSX,
     &      G_NSX,INDX)
          RHOLW_NSX = DIFMN(RHOLW_NX,RHOLW_SX,D_NSX,D_NSX,L_NSX,INDX)
          RHOLW_SAX = DIFMN(RHOLW_SX,RHOLW_AX,D_NSX,D_NSX,
     &      L_SAX(NSCN),INDX)
          RHOLW_PAX = DIFMN(RHOLW_PX,RHOLW_AX,D_NSX,D_NSX,
     &      L_PAIX(NSCN),INDX)
!
!---      Gas flow from ground surface to atmosphere, m^3/s
!         using air mass balance at ground surface  ---
!
          G_SAX = G_NSX*RHOGA_NSX/RHOGA_SAX
!
!---      Component and phase interfacial averages  ---
!
          INDX = -4
          HGW_SAX = DIFMN(HGW_GS(MS,NSCN),HGW_AX,D_NSX,D_NSX,E_SAX,INDX)
          HGW_NSX = DIFMN(HGW(MN,N),HGW_GS(MS,NSCN),D_NSX,D_NSX,
     &      E_NSX,INDX)
          HGA_SAX = DIFMN(HGA_GS(MS,NSCN),HGA_AX,D_NSX,D_NSX,E_SAX,INDX)
          HGA_NSX = DIFMN(HGA(MN,N),HGA_GS(MS,NSCN),D_NSX,D_NSX,
     &      E_NSX,INDX)
          HG_SAX = DIFMN(HG_GS(MS,NSCN),HG_AX,D_NSX,D_NSX,G_SAX,INDX)
          HG_NSX = DIFMN(HG(MN,N),HG_GS(MS,NSCN),D_NSX,D_NSX,G_NSX,INDX)
          HL_NSX = DIFMN(HL(MN,N),HL_GS(MS,NSCN),D_NSX,D_NSX,L_NSX,INDX)
          HL_SAX = DIFMN(HL_GS(MS,NSCN),HL_AX,D_NSX,D_NSX,
     &      L_SAX(NSCN),INDX)
          HL_PAX = DIFMN(HL_PX,HL_AX,D_NSX,D_NSX,L_PAIX(NSCN),INDX)
!
!---      Water-vapor density at canopy height, kg/m^3  ---
!
          PGA_CX = MAX( PA-PVW_CP(MS,NSCN),0.D+0 )
          XMLA_CX = PGA_CX/HCAW
          XMLW_CX = MAX( 1.D+0-XMLA_CX,0.D+0 )
          XLW_CX = XMLW_CX*WTMW/(XMLA_CX*WTMA + XMLW_CX*WTMW)
          XLA_CX = MAX( 1.D+0-XLW_CX,0.D+0 )
          INDX = 0
          CALL REGION_4( T_CP(MS,NSCN),PSWX,INDX )
          CALL P_IAPWS( T_CP(MS,NSCN),PVW_CP(MS,NSCN),RHOGW_CX,RHOX,
     &      HGW_CX,HX,UGWX,UX )
          CALL AIRGSD( T_CP(MS,NSCN),PGA_CX,RHOGA_CX )
          RHOG_CX = RHOGW_CX + RHOGA_CX
          XGW_CX = RHOGW_CX/RHOG_CX
          XGA_CX = RHOGA_CX/RHOG_CX
          E_CX = RHOGW_CX
!
!---      Water-vapor density at the plant, assuming
!         saturated vapor conditions inside the plant
!         leaf, kg/m^3  ---
!
          INDX = 0
          CALL REGION_4( T_PL(MS,NSCN),PSW_PX,INDX )
          PGW_PX = PSW_PX
          PGA_PX = MAX( PA-PGW_PX,0.D+0 )
          XMLA_PX = PGA_PX/HCAW
          XMLW_PX = MAX( 1.D+0-XMLA_PX,0.D+0 )
          XLW_PX = XMLW_PX*WTMW/(XMLA_PX*WTMA + XMLW_PX*WTMW)
          XLA_PX = MAX( 1.D+0-XLW_PX,0.D+0 )
          CALL P_IAPWS( T_PL(MS,NSCN),PGW_PX,RHOGW_PX,RHOX,
     &      HGW_PX,HX,UGW_PX,UX )
          CALL AIRGSD( T_PL(MS,NSCN),PGA_PX,RHOGA_PX )
          RHOG_PX = RHOGW_PX + RHOGA_PX
          XGW_PX = RHOGW_PX/RHOG_PX
          XGA_PX = RHOGA_PX/RHOG_PX
          E_PX = RHOGW_PX
!
!---      Diffusive water-vapor mass flux from the ground
!         surface to the canopy height, kg/m^2 s  ---
!
          PE_SCX = PAI_PX(NSCN)*(PE_SX-E_CX)/RA_SCX(NSCN)
          E_SCX = PAI_PX(NSCN)*(E_SX-E_CX)/RA_SCX(NSCN)
!
!---      Sensible heat flux from the plant to the
!         canopy height, W/m^2  ---
!
          H_PCX = 0.D+0
          DO IP = 1,NPLANT
            IF( IPLANT(IP).EQ.0 ) CYCLE
            H_PCX = H_PCX + PAIX(IP,NSCN)/RA_PCX(IP,NSCN)
          ENDDO
          H_PCX = H_PCX*(T_PL(MS,NSCN)-T_CP(MS,NSCN))*RHOG_AX*CPG_AX
!
!---      Sensible heat flux from the ground surface to the
!         canopy height, W/m^2  ---
!
          H_SCX = PAI_PX(NSCN)*(T_GS(MS,NSCN)-T_CP(MS,NSCN))*RHOG_AX*
     &      CPG_AX/RA_SCX(NSCN)
!
!---      Sensible heat flux from the canopy height
!         to the atmosphere, W/m^2  ---
!
          H_CAX = PAI_PX(NSCN)*(T_CP(MS,NSCN)-TA)*RHOG_AX*CPG_AX/
     &      RA_CAX(NSCN)
!
!---      Enthalpy of water-vapor, air and gas
!         at the canopy height, J/kg  ---
!
          CALL P_IAPWS( T_CP(MS,NSCN),PVW_CP(MS,NSCN),RHOGW_CX,RHOX,
     &      HGW_CX,HX,UGW_CX,UX )
          CALL AIRGSH( T_CP(MS,NSCN),PGA_CX,HGA_CX,UGA_CX )
          HG_CX = XGW_CX*HGW_CX + XGA_CX*HGA_CX
!
!---      Enthalpy of water-vapor, air and gas at the plant, J/kg  ---
!
          CALL AIRGSH( T_PL(MS,NSCN),PGA_PX,HGA_PX,UGA_PX )
          HG_PX = XGW_PX*HGW_PX + XGA_PX*HGA_PX
!
!---      Enthalpy of liquid-water at the node, J/kg  ---
!
          HLW_NX = (HL(MN,N)-XLA(MN,N)*HGA(MN,N))/XLW(MN,N)
!
!---      Diffusive water-vapor mass flux from the plant
!         to the canopy height, using aerodynamic and
!         stomatal resistances; attenuated by the crop
!         growth index and normalized stressed root uptake,
!         kg/m^2 s  ---
!
          E_PCX = 0.D+0
          PE_PCX = 0.D+0
          T_PCX(M,NSCN) = 0.D+0
          DO NC = 1,NSFCC
            IF( ICM_SFC(NC,NSCN).EQ.0 ) EXIT
            NCX = (NSCN-1)*NSFCC + NC
            RWU_SFC(M,NCX) = 0.D+0
          ENDDO
          IF( M.LE.6 ) TP_SFC(M,NSCN) = 0.D+0
          PT_PCX = 0.D+0
          RFIM_PX(1) = 0.D+0
          RFIM_PX(2) = 0.D+0
          L_SPX = 0.D+0
          DO IP = 1,NPLANT
            IF( IPLANT(IP).EQ.0 ) CYCLE
            DE_PCX = E_PX - E_CX
!
!---        Rainfall interception mass on plant leaves from
!           rainfall, kg/m^2 ground
!
!           RFIM_PL(2,IP,NSCN) - rainfall interception mass on plant
!                             leaves at current time for plant IP,
!                             for surface node NSCN, kg/m^2 ground
!           RFIM_PL(1,IP,NSCN) - rainfall interception mass on plant
!                             leaves at previous time step for plant IP,
!                             for surface node NSCN, kg/m^2 ground
!           DT - time step, s
!           L_PAX(IP,NSCN) - rainfall intensity on plant IP,
!                            m^3/s m^2 ground
!           RHOL_A - density of atmospheric water, kg/m^3
!
!            RFIM_PL(2,IP,NSCN) = RFIM_PL(1,IP,NSCN) - 
!     &        L_PAX(IP,NSCN)*RHOL_AX*DT
!
!---        Evaporation and transpiration from the plant
!           to the canopy  ---
!
            IF( DE_PCX.GT.0.D+0 ) THEN
!
!---          Transpiration, including stomatal resistance and
!             root-water uptake stress ---
!
              PT_PLX = PAIX(IP,NSCN)*DE_PCX/
     &          (RA_PCX(IP,NSCN) + RS_PCX(IP,NSCN))
              PT_PCX = PT_PCX + PT_PLX
              DO NC = 1,NSFCC
                IF( ICM_SFC(NC,NSCN).EQ.0 ) EXIT
                NCX = (NSCN-1)*NSFCC + NC
!
!---            Transpiration, including stomatal resistance and
!               root-water uptake stress, including increments
!               in primary variables at connected nodes ---
!
                TRNSPX = PT_PLX*C_COEFX(IP,NSCN)*
     &            (RWU_PX(MN,IP,NCX)/RWU_IX(IP,NSCN))
!
!---            Nodal root water uptake for connected nodes, including 
!               increments in primary variables at connected nodes ---
!
                RWU_SFC(M,NCX) = RWU_SFC(M,NCX) + TRNSPX
!
!---            Transpiration, including stomatal resistance and
!               root-water uptake stress, without increments
!               in primary variables at connected nodes ---
!
                TRNSPX = PT_PLX*C_COEFX(IP,NSCN)*
     &            (RWU_PX(MN,IP,NCX)/RWU_IX(IP,NSCN))
!
!---            Total root water uptake over connected nodes
!               (transpiration), without increments in primary 
!               variables at connected nodes ---
!
                TP_SFC(M,NSCN) = TP_SFC(M,NSCN) + TRNSPX
              ENDDO
            ENDIF
!
!---        Rainfall and condenstation greater than the maximum
!           rainfall interception capacity, water is shed
!           to the ground surface, kg/s m^2 ground
!
!           L_SPZ - shed rainfall and condensation from plant IP,
!           kg/s m^2 ground
!
!            IF( RFIM_PL(2,IP,NSCN).GT.RFIC_PX(IP,NSCN) ) THEN
!              L_SPZ = -(RFIM_PL(2,IP,NSCN)-RFIC_PX(IP,NSCN))/DT
!              RFIM_PL(2,IP,NSCN) = RFIC_PX(IP,NSCN)
!            ELSE
!              L_SPZ = 0.D+0
!            ENDIF
!
!---        Sum over plants, shed rainfall and condensate
!           mass kg/m^2 ground  ---
!
!            L_SPX = L_SPX + L_SPZ
!
!---        Sum rainfall interception mass and shed water mass over
!           plant species
!
!            RFIM_PX(1) = RFIM_PX(1) + RFIM_PL(1,IP,NSCN)
!            RFIM_PX(2) = RFIM_PX(2) + RFIM_PL(2,IP,NSCN)
          ENDDO
!
!---      Diffusive water-vapor mass flux from the canopy height
!         to the atmosphere, kg/m^2 s  ---
!
          E_CAX = PAI_PX(NSCN)*(E_CX-E_AX)/RA_CAX(NSCN)
!
!---      Component and phase interfacial averages  ---
!
          INDX = -4
          HGW_SCX = DIFMN(HGW_GS(MS,NSCN),HGW_CX,ONE,ONE,E_SCX,INDX)
          HGWE_PCX = DIFMN(HGW_PX,HGW_CX,ONE,ONE,E_PCX,INDX)
          HGWT_PCX = DIFMN(HGW_PX,HGW_CX,ONE,ONE,PT_PCX,INDX)
          HGW_CAX = DIFMN(HGW_CX,HGW_AX,ONE,ONE,E_CAX,INDX)
          HLW_SPX = DIFMN(HLW_GS(MS,NSCN),HLW_PL(MS,NSCN),ONE,ONE,
     &      L_SPX,INDX)
          HLW_NPX = DIFMN(HLW_NX,HLW_PL(MS,NSCN),ONE,ONE,PT_PCX,INDX)
!
!---      Ponding height limit  ---
!
          L_SAPMX = PERM(3,IZ(N))*(RHOLW_SAX**2)*GRVZ(NQZ)*
     &      (1.D+0-PHMX(NSCA)/(0.5D+0*DZGF(N)))/VISL_GS(MS,NSCN)
          RWROX = 0.D+0
          L_SAPX = L_SAX(NSCN)*RHOLW_SAX + L_SPX - 
     &      RWROX*RHOLW_SAX/AFZ(NQZ)
          IF( ABS(L_SAPX).GT.ABS(L_SAPMX) ) THEN
            L_SAPX = L_SAPMX
          ENDIF
!
!---      Residual of the water balance at the canopy height,
!         kg/s m^2 ground, (transpiration rate from all plants at the
!         current conditions in the connected subsurface nodes) ---
!
!         E_SCX - diffusive water-vapor flux, via evaporation, 
!         from ground surface to canopy, kg/s m^2 ground
!
!         E_PCX - diffusive water-vapor flux, via evaporation, 
!         from all plants to canopy, kg/s m^2 ground
!
!         TP_SFC(M,NSCN) - transpiration rate from all plants
!         at the current conditions in the connected subsurface nodes, 
!         kg/s m^2 ground
!
!         E_CAX - diffusive water-vapor flux, via evaporation, 
!         from canopy to the atmosphere, kg/s m^2 ground
!
          RW_CP(M,NSCN) = E_SCX + E_PCX + TP_SFC(M,NSCN) - E_CAX
!
!---      Residual of the energy balance at the canopy height,
!         W/m^2 ground, (transpiration rate from all plants at the
!         current conditions in the connected subsurface nodes) ---
!
!         H_SCX - conductive-convective heat flux from ground surface
!         to canopy, W/m^2 ground
!
!         E_SCX*HGW_SCX - enthalpy transfer rate of diffusive 
!         water-vapor flux, via evaporation, from ground surface to 
!         canopy, W/m^2 ground
!
!         H_PCX - conductive-convective heat flux from all plants
!         to canopy, W/m^2 ground
!
!         E_PCX*HGWE_PCX - enthalpy transfer rate of diffusive 
!         water-vapor flux, via evaporation, from all plants to 
!         canopy, W/m^2 ground
!
!         TP_SFC(M,NSCN)*HGWT_PCX - enthalpy transpiration rate 
!         from all plants at the current conditions in the  
!         connected subsurface nodes, W/m^2 ground
!
!         H_CAX - conductive-convective heat flux from canopy
!         to atmosphere, W/m^2 ground
!
!         E_CAX*HGW_CAX - enthalpy transfer rate of diffusive 
!         water-vapor flux, via evaporation, from canopy to 
!         the atmosphere, W/m^2 ground
!
          RE_CP(M,NSCN) = H_SCX + E_SCX*HGW_SCX
     &      + H_PCX + E_PCX*HGWE_PCX + TP_SFC(M,NSCN)*HGWT_PCX
     &      - H_CAX - E_CAX*HGW_CAX
!
!---      Residual of the energy balance at the plant,
!         W/m^2 ground, (transpiration rate from all plants at the
!         current conditions in the connected subsurface nodes) ---
!
!         RN_NP - net short- and long-wave radiation into all 
!         plants, W/m^2 ground
!
!         T_PCX(M,NSCN)*(HLW_NPX-HGWT_PCX) - rate of heat of 
!         transpiration for all plants at current conditions in
!         the connected subsurface nodes), W/m^2 ground
!
!         H_PC - advective heat transport between plant and canopy,
!         W/m^2 ground
!
!         E_PCX*HGWE_PC - enthalpy rate for evaporation of accumulated 
!         water, W/m^2 ground
!
!         ACW_PL(MS,NSCN)*HLW_PL(MS,NSCN) - enthalpy of accumulated
!         water on plant surfaces at current time step, J/m^2 ground
!
!         ACW_PL(1,NSCN)*HLW_PL(1,NSCN) - enthalpy of accumulated
!         water on plant surfaces at previous time step, J/m^2 ground
!
!         WM_PL(MS,NSCN)*HLW_PL(MS,NSCN) - enthalpy of plant water mass
!         at current time step, J/m^2 ground
!
!         WM_PL(MS,NSCN)*HLW_PL(MS,NSCN) - enthalpy of plant water mass
!         at previous time step, J/m^2 ground
!
!         DM_PL(MS,NSCN)*CP_PLX*T_PL(MS,NSCN) - enthalpy of plant dry
!         mass at current time step, J/m^2 ground
!
!         DM_PL(1,NSCN)*CP_PLX*T_PL(1,NSCN) - enthalpy of plant dry
!         mass at previous time step, J/m^2 ground
!
          RE_PL(M,NSCN) = RN_NPX + T_PCX(M,NSCN)*(HLW_NPX-HGWT_PCX)
     &      - H_PCX - E_PCX*HGWE_PCX - (ACW_PL(MS,NSCN)*HLW_PL(MS,NSCN)
     &      - ACW_PL(1,NSCN)*HLW_PL(1,NSCN))/DT
     &      - (WM_PL(MS,NSCN)*HLW_PL(MS,NSCN)
     &      - WM_PL(1,NSCN)*HLW_PL(1,NSCN))/DT
!     &      - (DM_PL(MS,NSCN)*CP_PLX*T_PL(MS,NSCN)
!     &      - DM_PL(1,NSCN)*CP_PLX*T_PL(1,NSCN))/DT
!
!---      Residual of the water balance at the ground surface,
!         kg/s m^2 ground  ---
!
          RW_GS(M,NSCN) = E_NSX + G_NSX*RHOGW_NSX + L_NSX*RHOLW_NSX
     &      - E_SAX - G_SAX*RHOGW_SAX - L_SAPX - E_SCX
!
!---      Residual of the energy balance at the ground surface,
!         ignoring dissolved air mass in aqueous,
!         W/m^2 ground  ---
!
          RE_GS(M,NSCN) = RN_NSX + E_NSX*HGW_NSX + G_NSX*RHOG_NSX*HG_NSX
     &      + L_NSX*RHOLW_NSX*HL_NSX + H_NSX
     &      - E_SAX*HGW_SAX - G_SAX*RHOG_SAX*HG_SAX
     &      - L_SAPX*HL_SAX - H_SAX
     &      - H_SCX - E_SCX*HGW_SCX
!          print '(2(a,i6),21(a,1pe12.5))',
!     &      'm = ',m,' nscn = ',nscn,' e_nsx = ',e_nsx,
!     &      ' e_sax = ',e_sax,' e_scx = ',e_scx,' g_nsx = ',g_nsx,
!     &      ' g_sax = ',g_sax,' h_nsx = ',h_nsx,' h_sax = ',h_sax,
!     &      ' h_scx = ',h_scx,' hg_nsx = ',hg_nsx,' hg_sax = ',hg_sax,
!     &      ' hgw_nsx = ',hgw_nsx,' hgw_sax = ',hgw_sax,
!     &      ' hgw_scx = ',hgw_scx,' hl_nsx = ',hl_nsx,
!     &      ' hl_sax = ',hl_sax,' l_nsx = ',l_nsx,' l_sapx = ',l_sapx,
!     &      ' rhog_nsx = ',rhog_nsx,' rhog_sax = ',rhog_sax,
!     &      ' rholw_nsx = ',rholw_nsx,' rn_nsx = ',rn_nsx

!
!---    Bottom of ground-surface indices and node indices loop  ---
!
        ENDDO
!
!---  Bottom of surface cover node loop  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLUX_SFC group ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_SFC
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
!     Increment surface cover primary variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 November 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
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
      SUB_LOG(ISUB_LOG) = '/INCRM_SFC'
!
!---  Loop over surface cover nodes  ---
!
      DO NSCN = 1,NSFCN
!
!---    Node adjacent to surface cover node  ---
!
        N = ICM_SFC(1,NSCN)
!
!---    Bare-surface energy - ground-surface temperature
!
        DNR_SFC(1,NSCN) = 1.D-6
!
!---    Bare-surface water mass - ground-surface aqueous pressure
!
        PLX = PL_GS(2,NSCN) + PATM
        PGX = PG(2,N) + PATM
        DNR_SFC(2,NSCN) = MAX( 1.D-1,1.D-7*ABS(PGX-PLX) )
!
!---    Plant-leaf energy - plant-leaf temperature
!
        DNR_SFC(3,NSCN) = 1.D-6
!
!---    Canopy energy - canopy temperature
!
        DNR_SFC(4,NSCN) = 1.D-6
!
!---    Canopy water mass - canopy water vapor partial pressure
!
        DNR_SFC(5,NSCN) = 1.D-6*PATM
!
!---    Increment the primary variables  ---
!
        DO M = 3,7
          T_GS(M,NSCN) = T_GS(2,NSCN)
          PL_GS(M,NSCN) = PL_GS(2,NSCN)
          T_PL(M,NSCN) = T_PL(2,NSCN)
          T_CP(M,NSCN) = T_CP(2,NSCN)
          PVW_CP(M,NSCN) = PVW_CP(2,NSCN)
          IF( M.EQ.3 ) THEN
            T_GS(M,NSCN) = T_GS(M,NSCN) + DNR_SFC(1,NSCN)
          ELSEIF( M.EQ.4 ) THEN
            PL_GS(M,NSCN) = PL_GS(M,NSCN) + DNR_SFC(2,NSCN)
          ELSEIF( M.EQ.5 ) THEN
            T_PL(M,NSCN) = T_PL(M,NSCN) + DNR_SFC(3,NSCN)
          ELSEIF( M.EQ.6 ) THEN
            T_CP(M,NSCN) = T_CP(M,NSCN) + DNR_SFC(4,NSCN)
          ELSEIF( M.EQ.7 ) THEN
            PVW_CP(M,NSCN) = PVW_CP(M,NSCN) + DNR_SFC(5,NSCN)
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_SFC group ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCB_SFC
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
!     STOMP-CO2e
!
!     Modify Jacobian matrix for the surface cover equations
!     and load Jacobian matrix for the surface cover equations
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 March 2016.
!
#ifdef petsc
!
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#endif
!
!----------------------Type Declarations-------------------------------!
!
#ifdef petsc
      REAL*8 BUFFER(LUK+1)
      INTEGER IROW(1), ICOL(LUK+1)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCB_SFC'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!!
!!---    Loop over surface cover nodes ---
!!
!        DO 500 NSCN = 1,NSFCN
!!
!!---      Node adjacent to surface cover node ---
!!
!          N = ICM_SFC(1,NSCN)
!!
!!---      Surface cover area index and surface cover type index  ---
!!
!          NSCA = ID_SFC(NSCN)
!          NSCT = NSFCT(NSCA)
!!
!!---      Bare-surface energy - ground-surface temperature
!!
!          MP = JM_SFC(1,NSCN)
!!
!!---      Change in residual for bare-surface energy with respect
!!         to changes in surface cover node primary variables  ---
!!
!          DO M = 1,5
!            MCOL = JM_SFC(M,NSCN)
!            MROW = MP-MCOL+MDC
!            ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &        (RE_GS(1,NSCN)-RE_GS(M+1,NSCN))/DNR_SFC(1,NSCN)
!          ENDDO
!!
!!---      Change in residual for bare-surface energy with respect
!!         to changes in adjacent node primary variables  ---
!!
!          DO M = 1,ISVC
!            MCOL = IM(M,NMD)
!            MROW = MP-MCOL+MDC
!            ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &        (RE_GS(1,NSCN)-RE_GS(M+6,NSCN))/DNR(M,N)
!          ENDDO
!!
!!---      Residual for bare-surface water mass   ---
!!
!          BLU(MP) = BLU(MP) - RE_GS(1,NSCN)
!!
!!---      Bare-surface water mass - ground-surface aqueous pressure
!!
!          MP = JM_SFC(2,NSCN)
!!
!!---      Change in residual for bare-surface water mass with respect
!!         to changes in surface cover node primary variables  ---
!!
!          DO M = 1,5
!            MCOL = JM_SFC(M,NSCN)
!            MROW = MP-MCOL+MDC
!            ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &        (RW_GS(1,NSCN)-RW_GS(M+1,NSCN))/DNR_SFC(2,NSCN)
!          ENDDO
!!
!!---      Residual for bare-surface water mass   ---
!!
!          BLU(MP) = BLU(MP) - RW_GS(1,NSCN)
!!
!!---      Plant-leaf energy - plant-leaf temperature
!!
!          MP = JM_SFC(3,NSCN)
!!
!!---      Change in residual for plant-leaf energy with respect
!!         to changes in surface cover node primary variables  ---
!!
!          DO M = 1,5
!            MCOL = JM_SFC(M,NSCN)
!            MROW = MP-MCOL+MDC
!            ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &        (RE_PL(1,NSCN)-RE_PL(M+1,NSCN))/DNR_SFC(3,NSCN)
!          ENDDO
!!
!!---      Residual for bare-surface water mass   ---
!!
!          BLU(MP) = BLU(MP) - RE_PL(1,NSCN)
!!
!!---      Canopy water mass - canopy water vapor partial pressure
!!
!          MP = JM_SFC(4,NSCN)
!!
!!---      Change in residual for canopy water mass with respect
!!         to changes in surface cover node primary variables  ---
!!
!          DO M = 1,5
!            MCOL = JM_SFC(M,NSCN)
!            MROW = MP-MCOL+MDC
!            ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &        (RW_CP(1,NSCN)-RW_CP(M+1,NSCN))/DNR_SFC(4,NSCN)
!          ENDDO
!!
!!---      Residual for bare-surface water mass   ---
!!
!          BLU(MP) = BLU(MP) - RE_PL(1,NSCN)
!!
!!---      Canopy energy - canopy temperature
!!
!          MP = JM_SFC(5,NSCN)
!!
!!---      Change in residual for canopy energy with respect
!!         to changes in surface cover node primary variables  ---
!!
!          DO M = 1,5
!            MCOL = JM_SFC(M,NSCN)
!            MROW = MP-MCOL+MDC
!            ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &        (RE_CP(1,NSCN)-RE_CP(M+1,NSCN))/DNR_SFC(5,NSCN)
!          ENDDO
!!
!!---      Residual for bare-surface water mass   ---
!!
!          BLU(MP) = BLU(MP) - RE_CP(1,NSCN)
!
!---      Loop over coupled-well well nodes ---
!
!          DO 200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!            N = IWN_CW(NWN)
!            NMD = ABS(IXP(N))
!!
!!---        Nonisothermal simulations  ---
!!
!            IF( ISLC(30).EQ.0 ) THEN
!!
!!---          Energy balance equation at field node ---
!!
!              MP = IM(IEQT,NMD)
!!
!!---          Change in energy flow into field node with respect
!!             to change in field node primary variables  ---
!!
!              DO 90 M = 1,ISVC
!                DNRX = DNR(M,N)
!                MCOL = IM(M,NMD)
!                MROW = MP-MCOL+MDC
!                ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &            (FXE_CW(1,NWN)-FXE_CW(M+1,NWN))/DNRX
!   90         CONTINUE
!!
!!---          Change in energy flow into field node with respect
!!             to change in coupled-well pressure  ---
!!
!              MCOL = JM_CW(NCW)
!              MROW = MP-MCOL+MDC
!              MX = ISVC+2
!              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &          (FXE_CW(1,NWN)-FXE_CW(MX,NWN))/DNR_CW(NCW)
!!
!!---          Energy flow into field node, W  ---
!!
!              BLU(MP) = BLU(MP) + FXE_CW(1,NWN)
!              RSDL(IEQT,N) = BLU(MP)
!            ENDIF
!!
!!---        Water mass balance equation at field node ---
!!
!            MP = IM(IEQW,NMD)
!!
!!---        Change in water mass flux into field node with respect
!!           to change in field node primary variables  ---
!!
!            DO 100 M = 1,ISVC
!              DNRX = DNR(M,N)
!              MCOL = IM(M,NMD)
!              MROW = MP-MCOL+MDC
!              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &          (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNRX
!  100       CONTINUE
!!
!!---        Change in water mass flux into field node with respect
!!           to change in coupled-well pressure  ---
!!
!            MCOL = JM_CW(NCW)
!            MROW = MP-MCOL+MDC
!            MX = ISVC+2
!            ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &        (FXW_CW(1,NWN)-FXW_CW(MX,NWN))/DNR_CW(NCW)
!!
!!---        Water mass flux into field node, kg/s  ---
!!
!            BLU(MP) = BLU(MP) + FXW_CW(1,NWN)
!            RSDL(IEQW,N) = BLU(MP)
!!
!!---        CO2 mass balance equation at field node ---
!!
!            MP = IM(IEQA,NMD)
!!
!!---        Change in CO2 mass flux into field node with respect
!!           to change in field node primary variables  ---
!!
!            DO 110 M = 1,ISVC
!              DNRX = DNR(M,N)
!              MCOL = IM(M,NMD)
!              MROW = MP-MCOL+MDC
!              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &          (FXA_CW(1,NWN)-FXA_CW(M+1,NWN))/DNRX
!  110       CONTINUE
!!
!!---        Change in CO2 mass flux into field node with respect
!!           to change in coupled-well pressure  ---
!!
!            MCOL = JM_CW(NCW)
!            MROW = MP-MCOL+MDC
!            MX = ISVC+2
!            ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &        (FXA_CW(1,NWN)-FXA_CW(MX,NWN))/DNR_CW(NCW)
!!
!!---        CO2 mass flux into field node, kg/s  ---
!!
!            BLU(MP) = BLU(MP) + FXA_CW(1,NWN)
!            RSDL(IEQA,N) = BLU(MP)
!  200     CONTINUE
!!
!!---      Coupled-well mass balance  ---
!!
!          MP = JM_CW(NCW)
!          BLU(MP) = -RS_CW(1,NCW)
!!
!!---      Pressure controlled coupled well  ---
!!
!          IF( ID_CW(8,NCW).EQ.1 ) BLU(MP) = 0.D+0
!!
!!---      Change in coupled-well mass balance with respect to
!!         change in coupled-well pressure  ---
!!
!          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
!          MX = (NWFX*ISVC)+2
!          MCOL = JM_CW(NCW)
!          MROW = MP-MCOL+MDC
!          ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &      (RS_CW(MX,NCW)-RS_CW(1,NCW))/DNR_CW(NCW)
!!
!!---      Pressure controlled coupled well  ---
!!
!          IF( ID_CW(8,NCW).EQ.1 ) ALU(MROW,MCOL) = 1.D+0
!!
!!---      Loop over field nodes with coupled-well nodes ---
!!
!          DO 400 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
!            N = IWF_CW(NWF)
!            NMD = ABS(IXP(N))
!            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
!!
!!---        Change in coupled-well mass balance with respect to
!!           change in field node primary variables  ---
!!
!            DO 300 M = 1,ISVC
!              DNRX = DNR(M,N)
!              MCOL = IM(M,NMD)
!              MROW = MP-MCOL+MDC
!              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
!     &          (RS_CW(MX+M,NCW)-RS_CW(1,NCW))/DNRX
!!
!!---          Pressure controlled coupled well  ---
!!
!              IF( ID_CW(8,NCW).EQ.1 ) ALU(MROW,MCOL) = 0.D+0
!  300       CONTINUE
!  400     CONTINUE
!  500   CONTINUE
!
!---  SPLib solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!!
!!---    Loop over coupled wells ---
!!
!        DO 1500 NCW = 1,N_CW
!!
!!---      Loop over coupled-well well nodes ---
!!
!          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
!          DO 1200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!            N = IWN_CW(NWN)
!            NMD = ABS(IXP(N))
!!
!!---        Nonisothermal simulations  ---
!!
!            IF( ISLC(30).EQ.0 ) THEN
!!
!!---          Energy balance equation at field node ---
!!
!              MP = IM(IEQT,NMD)
!!
!!---          Change in energy flow into field node with respect
!!             to change in field node primary variables  ---
!!
!              MA = 3*ISVC 
!              DO 1090 M = 1,ISVC
!                DNRX = DNR(M,N)
!                MCOL = KLU(MP,M+MA)
!                DLU(MCOL) = DLU(MCOL) + 
!     &            (FXE_CW(1,NWN)-FXE_CW(M+1,NWN))/DNRX
! 1090         CONTINUE
!!
!!---          Change in energy flow with respect
!!             to change in coupled-well pressure  ---
!!
!!              MA = NWFX*ISVC + (IWP_CW(NWN)-1)*ISVC + IEQT + 1
!              MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + 
!     &          IEQT + 1
!              MCOL = KLU_CW(MA,NCW)
!              MX = ISVC+2
!              DLU(MCOL) = DLU(MCOL) + 
!     &          (FXE_CW(1,NWN)-FXE_CW(MX,NWN))/DNR_CW(NCW)
!!
!!---          Energy flow into field node, W  ---
!!
!              BLU(MP) = BLU(MP) + FXE_CW(1,NWN)
!              RSDL(IEQT,N) = BLU(MP)
!            ENDIF
!!
!!---        Water mass balance equation at field node ---
!!
!            MP = IM(IEQW,NMD)
!!
!!---        Change in water mass flux into field node with respect
!!           to change in field node primary variables  ---
!!
!            MA = 3*ISVC 
!            DO 1100 M = 1,ISVC
!              DNRX = DNR(M,N)
!              MCOL = KLU(MP,M+MA)
!              DLU(MCOL) = DLU(MCOL) + 
!     &          (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNRX
! 1100       CONTINUE
!!
!!---        Change in water fraction of coupled-well flux with respect
!!           to change in coupled-well pressure  ---
!!
!!            MA = NWFX*ISVC + (IWP_CW(NWN)-1)*ISVC + IEQW + 1
!            MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + IEQW + 1
!            MCOL = KLU_CW(MA,NCW)
!            MX = ISVC+2
!            DLU(MCOL) = DLU(MCOL) + 
!     &        (FXW_CW(1,NWN)-FXW_CW(MX,NWN))/DNR_CW(NCW)
!!
!!---        Water mass flux into field node, kg/s  ---
!!
!            BLU(MP) = BLU(MP) + FXW_CW(1,NWN)
!            RSDL(IEQW,N) = BLU(MP)
!!
!!---        CO2 mass balance equation at field node ---
!!
!            MP = IM(IEQA,NMD)
!!
!!---        Change in CO2 fraction of coupled-well flux with respect
!!           to change in field node primary variables  ---
!!
!            MA = 3*ISVC 
!            DO 1110 M = 1,ISVC
!              MCOL = KLU(MP,M+MA)
!              DNRX = DNR(M,N)
!              DLU(MCOL) = DLU(MCOL) + 
!     &          (FXA_CW(1,NWN)-FXA_CW(M+1,NWN))/DNRX
! 1110       CONTINUE
!!
!!---        Change in CO2 fraction of coupled-well flux with respect
!!           to change in coupled-well pressure  ---
!!
!!            MA = NWFX*ISVC + (IWP_CW(NWN)-1)*ISVC + IEQA + 1
!            MA = NWFX*ISVC + (IWP_CW(NWN)-ID_CW(5,NCW))*ISVC + IEQA + 1
!            MCOL = KLU_CW(MA,NCW)
!            MX = ISVC+2
!            DLU(MCOL) = DLU(MCOL) + 
!     &        (FXA_CW(1,NWN)-FXA_CW(MX,NWN))/DNR_CW(NCW)
!!
!!---        CO2 mass flux into field node, kg/s  ---
!!
!            BLU(MP) = BLU(MP) + FXA_CW(1,NWN)
!            RSDL(IEQA,N) = BLU(MP)
! 1200     CONTINUE
!!
!!---      Coupled-well mass balance  ---
!!
!          MP = JM_CW(NCW)
!          BLU(MP) = -RS_CW(1,NCW)
!!
!!---      Pressure controlled coupled well  ---
!!
!          IF( ID_CW(8,NCW).EQ.1 ) BLU(MP) = 0.D+0
!!
!!---      Change in coupled-well mass balance with respect to
!!         change in coupled-well pressure  ---
!!
!          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
!          MX = (NWFX*ISVC)+2
!          MA = 1
!          MCOL = KLU_CW(MA,NCW)
!          DLU(MCOL) = DLU(MCOL) + 
!     &      (RS_CW(MX,NCW)-RS_CW(1,NCW))/DNR_CW(NCW)
!!
!!---      Pressure controlled coupled well  ---
!!
!          IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 1.D+0
!!
!!---      Loop over field nodes with coupled-well nodes ---
!!
!          DO 1400 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
!            N = IWF_CW(NWF)
!            NMD = ABS(IXP(N))
!            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
!            MA = (NWF-ID_CW(5,NCW))*ISVC + 1
!!
!!---        Change in coupled-well mass balance with respect to
!!           change in field node primary variables  ---
!!
!            DO 1300 M = 1,ISVC
!              DNRX = DNR(M,N)
!              MCOL = KLU_CW(M+MA,NCW)
!              DLU(MCOL) = DLU(MCOL) + 
!     &          (RS_CW(MX+M,NCW)-RS_CW(1,NCW))/DNRX
!!
!!---          Pressure controlled coupled well  ---
!!
!              IF( ID_CW(8,NCW).EQ.1 ) DLU(MCOL) = 0.D+0
! 1300       CONTINUE
! 1400     CONTINUE
! 1500   CONTINUE
#ifdef petsc
!
!---  PETSc solver  ---
!
      ELSEIF( ILES.EQ.5 ) THEN
!!
!!---    Loop over coupled wells ---
!!
!        DO 2500 NCW = 1,N_CW
!!
!!---      Loop over coupled-well well nodes ---
!!
!          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
!          DO 2200 NWN = ID_CW(3,NCW),ID_CW(4,NCW)
!            N = IWN_CW(NWN)
!            NMD = ABS(IXP(N))
!!
!!---        Nonisothermal simulations  ---
!!
!            IF( ISLC(30).EQ.0 ) THEN
!!
!!---          Energy balance equation at field node ---
!!
!              MP = IM(IEQT,NMD)
!              IROW(1) = MP-1
!!
!!---          Change in energy flow into field node with respect
!!             to change in field node primary variables  ---
!!
!              DO 2090 M = 1,ISVC
!                DNRX = DNR(M,N)
!                ICOL(M) = IM(M,NMD)-1
!                BUFFER(M) = (FXE_CW(1,NWN)-FXE_CW(M+1,NWN))/DNRX
! 2090         CONTINUE
!!
!!---          Change in energy flow with respect
!!             to change in coupled-well pressure  ---
!!
!              MX = ISVC+2
!              M=ISVC+1
!              ICOL(M) = JM_CW(NCW)-1
!              BUFFER(M) = (FXE_CW(1,NWN)-FXE_CW(MX,NWN))/DNR_CW(NCW)
!              CALL MatSetValues(F_MAT,1,IROW,ISVC+1,ICOL,BUFFER,
!     &                          ADD_VALUES,IERR)
!!
!!---          Energy flow into field node, W  ---
!!
!              BLU(MP) = BLU(MP) + FXE_CW(1,NWN)
!              RSDL(IEQT,N) = BLU(MP)
!            ENDIF
!!
!!---        Water mass balance equation at field node ---
!!
!            MP = IM(IEQW,NMD)
!            IROW(1) = MP-1
!!
!!---        Change in water mass flow into field node with respect
!!           to change in field node primary variables  ---
!!
!            DO 2100 M = 1,ISVC
!              DNRX = DNR(M,N)
!              ICOL(M) = IM(M,NMD)-1
!              BUFFER(M) = (FXW_CW(1,NWN)-FXW_CW(M+1,NWN))/DNRX
! 2100       CONTINUE
!!
!!---        Change in water mass flow with respect
!!           to change in coupled-well pressure  ---
!!
!            MX = ISVC+2
!            M=ISVC+1
!            ICOL(M) = JM_CW(NCW)-1
!            BUFFER(M) = (FXW_CW(1,NWN)-FXW_CW(MX,NWN))/DNR_CW(NCW)
!            CALL MatSetValues(F_MAT,1,IROW,ISVC+1,ICOL,BUFFER,
!     &                        ADD_VALUES,IERR)
!!
!!---        Water mass flow into field node, kg/s  ---
!!
!            BLU(MP) = BLU(MP) + FXW_CW(1,NWN)
!            RSDL(IEQW,N) = BLU(MP)
!!
!!---        CO2 mass balance equation at field node ---
!!
!            MP = IM(IEQA,NMD)
!            IROW(1) = MP-1
!!
!!---        Change in CO2 mass flow with respect
!!           to change in field node primary variables  ---
!!
!            DO 2120 M = 1,ISVC
!              ICOL(M) = IM(M,NMD)-1
!              DNRX = DNR(M,N)
!              BUFFER(M) = (FXA_CW(1,NWN)-FXA_CW(M+1,NWN))/DNRX
! 2120       CONTINUE
!!
!!---        Change in CO2 mass flow with respect
!!           to change in coupled-well pressure  ---
!!
!            MX = ISVC+2
!            M=ISVC+1
!            ICOL(M) = JM_CW(NCW)-1
!            BUFFER(M) = (FXA_CW(1,NWN)-FXA_CW(MX,NWN))/DNR_CW(NCW)
!            CALL MatSetValues(F_MAT,1,IROW,ISVC+1,ICOL,BUFFER,
!     &                        ADD_VALUES,IERR)
!!
!!---        CO2 mass flow into field node, kg/s  ---
!!
!            BLU(MP) = BLU(MP) + FXA_CW(1,NWN)
!            RSDL(IEQA,N) = BLU(MP)
! 2200     CONTINUE
!!
!!---      Coupled-well mass balance  ---
!!
!          MP = JM_CW(NCW)
!          IROW(1) = MP-1
!          BLU(MP) = -RS_CW(1,NCW)
!!
!!---      Pressure controlled coupled well  ---
!!
!          IF( ID_CW(8,NCW).EQ.1 ) BLU(MP) = 0.D+0
!!
!!---      Change in coupled-well mass balance with respect to
!!         change in coupled-well pressure  ---
!!
!          NWFX = ID_CW(6,NCW)-ID_CW(5,NCW)+1
!          MX = (NWFX*ISVC)+2
!          M = 1
!          ICOL(M) = JM_CW(NCW)-1
!          BUFFER(M) = (RS_CW(MX,NCW)-RS_CW(1,NCW))/DNR_CW(NCW)
!!
!!---      Pressure controlled coupled well  ---
!!
!          IF( ID_CW(8,NCW).EQ.1 ) BUFFER(M) = 1.D+0
!          CALL MatSetValues(F_MAT,1,IROW,1,ICOL,BUFFER,
!     &                      ADD_VALUES,IERR)
!!
!!---      Loop over field nodes with coupled-well nodes ---
!!
!          DO 2400 NWF = ID_CW(5,NCW),ID_CW(6,NCW)
!            N = IWF_CW(NWF)
!            NMD = ABS(IXP(N))
!            MX = (NWF-ID_CW(5,NCW))*ISVC + 1
!!
!!---        Change in coupled-well mass balance with respect to
!!           change in field node primary variables  ---
!!
!            DO 2300 M = 1,ISVC
!              DNRX = DNR(M,N)
!              ICOL(M) = IM(M,NMD)-1
!              BUFFER(M) = (RS_CW(MX+M,NCW)-RS_CW(1,NCW))/DNRX
!!
!!---          Pressure controlled coupled well  ---
!!
!              IF( ID_CW(8,NCW).EQ.1 ) BUFFER(M) = 0.D+0
! 2300       CONTINUE
!            CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
!     &                        ADD_VALUES,IERR)
! 2400     CONTINUE
! 2500   CONTINUE
#endif
      ENDIF
!
!---  Reset subroutine character string ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCB_SFC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LOCPT( X0,Y0,X,Y,N,IINOUT,IPATH )
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
!     LOCPT locates whether a point is inside a polygon in a 2D plane
!
!     Method:
!
!     Given a polygonal line connecting the vertices 
!     (x(i),y(i)) (i = 1,...,n) taken in this order.  It is assumed that
!     the polygonal path is a loop, where (x(n),y(n)) = (x(1),y(1)) or 
!     there is an arc from (x(n),y(n)) to (x(1),y(1)). The polygon may 
!     cross itself any number of times.
!
!     Reference:
!
!     Fortran 66 version by A.H. Morris
!     Converted to ELF90 compatibility by Alan Miller, 15 February 1997
!     Available at http://jblevins.org/mirror/amiller/    
!
!     Author:
!
!     John Burkardt
!
!     Parameters:
!
!     Input, real X0, Y0, the coordinate of the point.
!     Input, real X(N), Y(N), the coordinates of the points defining 
!       the polygon.
!     Input, integer N, the number of the points defining the polygon.
! 
!     Output, integer IINOUT, assigned as follows:
!       L = -1   if (X0,Y0) is outside the polygonal path
!       L =  0   if (X0,Y0) lies on the polygonal path
!       L =  1   if (X0,Y0) is inside the polygonal path
!
!     Output, integer IPATH, where IPATH = 0 if (X0,Y0) is on or outside 
!       the path. If (X0,Y0) is inside the path then m is the winding 
!       number of the path around the point (X0,Y0).
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK White, PNNL, 15 September 2015.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 X(N),Y(N)
!
!----------------------Executable Lines--------------------------------!
!      
      EPSL = 1.D-14
      N0 = N
      IF (X(1) == X(N) .AND. Y(1) == Y(N)) N0 = N - 1
      PIX = ATAN2(0.D+0,-1.D+0)
      PI2 = 2.D+0 * PIX
      TOL = 4.D+0*EPSL*PIX
      IINOUT = -1
      IPATH = 0
!      
      U = X(1) - X0
      V = Y(1) - Y0
      IF (U .EQ. 0.D+0 .AND. V .EQ. 0.D+0) GO TO 20
      IF (N0 .LT. 2) RETURN
      THETA1 = ATAN2(V, U)
!      
      SUM = 0.0
      THETA = THETA1
      DO I = 2, N0
        U = X(I) - X0
        V = Y(I) - Y0
        IF (U .EQ. 0.D+0 .AND. V .EQ. 0.D+0) GO TO 20
        THETAI = ATAN2(V, U)
!        
        ANGLE = ABS(THETAI - THETA)
        IF (ABS(ANGLE - PIX) .LT. TOL) GO TO 20
        IF (ANGLE .GT. PIX) ANGLE = ANGLE - PI2
        IF (THETA .GT. THETAI) ANGLE = -ANGLE
        SUM = SUM + ANGLE
        THETA = THETAI
      END DO
!      
      ANGLE = ABS(THETA1 - THETA)
      IF (ABS(ANGLE - PIX) .LT. TOL) GO TO 20
      IF (ANGLE .GT. PIX) ANGLE = ANGLE - PI2
      IF (THETA .GT. THETA1) ANGLE = -ANGLE
      SUM = SUM + ANGLE
!      
!---  SUM = 2*PIX*IPATH WHERE IPATH IS THE WINDING NUMBER  ---
!      
      IPATH = INT( ABS(SUM)/PI2 + 0.2D+0 )
      IF (IPATH .EQ. 0) RETURN
      IINOUT = 1
      IF (SUM .LT. 0.D+0) M = -1 * IPATH
      RETURN
!      
!---  (X0, Y0) IS ON THE BOUNDARY OF THE PATH  ---
!      
   20 IINOUT = 0
!
!---  End of LOCPT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE POLINT( XA,YA,N,X,Y,DY,IERR )
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
!     Given arrays XA and YA, each of length N, and given a value X,
!     this routine returns a value Y, and an error estimate DY.
!     If P(X) is the polynomial of degree N-1 such that
!     P(XA(I)) = YA(I),
!     I = 1,N, then the returned value Y = P(X).
!
!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!     1986.  Numerical Recipes, The Art of Scientific Computing.
!     Cambridge University Press, Cambridge.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, February, 1999.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Parameter Statements----------------------------!
!
      PARAMETER( NMAX=10, EPSL=1.D-14 )
      REAL*8 XA(N),YA(N),C(NMAX),D(NMAX)
!
!----------------------Executable Lines--------------------------------!
!
      IERR = 0
      NS = 1
      DIF = ABS(X-XA(1))
!
!---  Find the index NS of the closest table entry  ---
!
      DO 11 I = 1,N
        DIFT = ABS(X-XA(I))
        IF( DIFT.LT.DIF ) THEN
          NS = I
          DIF = DIFT
        ENDIF
!
!---  Initialize the tableau fo C's and D's  ---
!
        C(I) = YA(I)
        D(I) = YA(I)
   11 CONTINUE
!
!---  Initial approximation to Y  ---
!
      Y = YA(NS)
      NS = NS-1
!
!---  Loop over the columns in the tableau and
!     update the C's and D's  ---
!
      DO 13 M = 1,N-1
        DO 12 I = 1,N-M
          HO = XA(I)-X
          HP = XA(I+M)-X
          W = C(I+1)-D(I)
          DEN = HO-HP
!
!---  Exit subroutine if two input XA's are
!     identical, within roundoff  ---
!
          IF( ABS( DEN )/EPSL.LT.EPSL ) THEN
            IERR = 1
            RETURN
          ENDIF
          DEN = W/DEN
!
!---  Update C's and D's  ---
!
          D(I) = HP*DEN
          C(I) = HO*DEN
   12   CONTINUE
!
!---  After each column in the tableau is completed, decide
!     which direction C or D to add the accumulating value
!     of Y, (i.e., which path to take through the tableau -
!     forking up or down).  Do this in such a way as
!     to take the most straight line route through the
!     tableau to its apex, updating NS accordingly to keep track.
!     This route keeps the partial approximations centered
!     (insofar as possible) on the target X.  The
!     last DY added is thus the error indication.  ---
!
        IF( 2*NS.LT.N-M ) THEN
          DY = C(NS+1)
        ELSE
          DY = D(NS)
          NS = NS-1
        ENDIF
        Y = Y+DY
   13 CONTINUE
!
!---  End of POLINT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_SFC( PGX )
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
!     Surface cover secondary variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 December 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE PLT_ATM
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
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP_SFC'
!
!---  Zero gas entrapment at ground surface  ---
!
      SGTX = 0.D+0
      ESGTX = 0.D+0
      ESGTMX = 0.D+0
      CALL SCF_GL( BTGLX,PGX )
!
!---  Zero salt concentration at ground surface  ---
!
      XLSX = 0.D+0
      XMLSX = 0.D+0
!
!---  Atmospheric guage pressure  ---
!
      PGGX = PGX - PATM
!
!---  Loop over surface cover nodes  ---
!
      DO NSCN = 1,NSFCN
!
!---    Node adjacent to surface cover node  ---
!
        N = ICM_SFC(1,NSCN)
        IZN = IZ(N)
!
!---    Loop over surface cover increments  ---
!
        DO M = 2,7
!
!---      Saturation and relative permeability at the 
!         ground surface  ---
!
          INDX = 1
          CALL KSP_GT( N,PGGX,PL_GS(M,NSCN),BTGLX,SG_GS(M,NSCN),
     &      SGTX,SL_GS(M,NSCN),SDPFX,SDPMX,RKLX(1),
     &      RKG_GS(M,NSCN),ASLX,ASLMINX,ESGTX,ESGTMX,SLRX,INDX )
          RKL_GS(M,NSCN) = RKLX(1)
          PLX = PL_GS(M,NSCN) + PATM
          PX = MAX( PGX,PLX )
          INDX = 0
          CALL REGION_4( T_GS(M,NSCN),PSWX,INDX )
          CALL P_IAPWS( T_GS(M,NSCN),PSWX,RHOX,RHOLWX,HX,
     &      HLW_GS(M,NSCN),UX,ULWX )
          PCX = MAX( PGGX-PL_GS(M,NSCN),0.D+0 )
          CALL VPL_B( T_GS(M,NSCN),PSWX,PCX,RHOLWX,PVWX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,T_GS(M,NSCN) )
          CALL P_IAPWS( T_GS(M,NSCN),PVWX,RHOGWX,RHOX,HGW_GS(M,NSCN),
     &      HX,UGWX,UX )
          PVA_GS(M,NSCN) = PGX-PVWX
          IF( PVA_GS(M,NSCN).LT.1.D-6 ) PVA_GS(M,NSCN) = 0.D+0
          PVW_GS(M,NSCN) = PVWX
          XMLA_GS(M,NSCN) = PVA_GS(M,NSCN)/HCAW
          XMGA_GS(M,NSCN) = PVA_GS(M,NSCN)/PGX
          CALL EQUIL( XGA_GS(M,NSCN),XGW_GS(M,NSCN),XLA_GS(M,NSCN),
     &      XLSX,XLW_GS(M,NSCN),XMGA_GS(M,NSCN),
     &      XMGW_GS(M,NSCN),XMLA_GS(M,NSCN),XMLSX,XMLW_GS(M,NSCN) )
!
!---      Gas density and component fractions  ---
!
          CALL AIRGSD( T_GS(M,NSCN),PVA_GS(M,NSCN),RHOGAX )
          RHOG_GS(M,NSCN) = XGA_GS(M,NSCN)*RHOGAX + 
     &      XGW_GS(M,NSCN)*RHOGWX
          WTMGX = XMGA_GS(M,NSCN)*WTMA + XMGW_GS(M,NSCN)*WTMW
          RHOMG_GS(M,NSCN) = RHOG_GS(M,NSCN)/WTMGX
!
!---      Gas viscosity  ---
!
          CALL AIRGSV( T_GS(M,NSCN),VISGAX )
          CALL VISC_W( T_GS(M,NSCN),PVBX,RHOGWX,VISGWX )
          CALL VISC_G( VISGAX,VISGWX,XMGA_GS(M,NSCN),XMGW_GS(M,NSCN),
     &      VISG_GS(M,NSCN) )
!
!---      Aqueous density and molar density  ---
!
          RHOL_GS(M,NSCN) = RHOBX
          WTMLX = XMLA_GS(M,NSCN)*WTMA + XMLW_GS(M,NSCN)*WTMW
          RHOML_GS(M,NSCN) = RHOL_GS(M,NSCN)/WTMLX
!
!---      Aqueous viscosity  ---
!
          CALL VISC_W( T_GS(M,NSCN),PX,RHOLWX,VISLWX )
          CALL VISC_B( T_GS(M,NSCN),XLSX,VISLWX,VISBX )
          CALL VISC_L( XMLA_GS(M,NSCN),VISBX,VISGAX,VISL_GS(M,NSCN) )
!
!---      Gas-water diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFGW_GS(M,NSCN) = DFGWC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL BNDFAW( T_GS(M,NSCN),PGX,DFGW_GS(M,NSCN) )
          ELSEIF( ISLC(2).EQ.3 ) THEN
            CALL BNDFAW( T_GS(M,NSCN),PGX,DFGW_GS(M,NSCN) )
            CMFF = 1.D+0 + 2.6D+0/(DFGWC**0.5)
            AMC = PORD(2,N)*SL_GS(M,NSCN)
            ENHF = 9.5D+0 + 6.D+0*(AMC) -
     &        8.5D+0/EXP((CMFF*AMC)**4)
            DFGW_GS(M,NSCN) = ENHF*DFGW_GS(M,NSCN)
          ELSEIF( ISLC(2).EQ.4 ) THEN
            CALL BNDFAW( T_GS(M,NSCN),PGX,DFGW_GS(M,NSCN) )
            ENHF = DFEF(1,IZ(N))+DFEF(2,IZ(N))*SL_GS(M,NSCN)-
     &        (DFEF(1,IZ(N))-DFEF(4,IZ(N)))
     &        *EXP(-((DFEF(3,IZ(N))*SL_GS(M,NSCN))**DFEF(5,IZ(N))))
            DFGW_GS(M,NSCN) = ENHF*DFGW_GS(M,NSCN)
          ENDIF
!
!---      Aqueous-air diffusion coefficients  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLA_GS(M,NSCN) = DFLAC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            CALL AIRDFL( T_GS(M,NSCN),VISL_GS(M,NSCN),DFLA_GS(M,NSCN) )
          ENDIF
!
!---      Aqueous and gas tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL_GS(M,NSCN),
     &      SG_GS(M,NSCN),ZERO,PORD(2,N),TORL_GS(M,NSCN),
     &      TORG_GS(M,NSCN),TORNX )
!
!---      Gas enthalpy and internal energy  ---
!
          CALL AIRGSH( T_GS(M,NSCN),PVA_GS(M,NSCN),HGA_GS(M,NSCN),UGAX )
          UEG_GS(M,NSCN) = XGA_GS(M,NSCN)*UGAX + XGW_GS(M,NSCN)*UGWX
          HG_GS(M,NSCN) = XGA_GS(M,NSCN)*HGA_GS(M,NSCN) + 
     &      XGW_GS(M,NSCN)*HGW_GS(M,NSCN)
!
!---      Gas thermal conductivity  ---
!
          CALL AIRGSK( T_GS(M,NSCN),THKGAX )
          CALL THK_W( T_GS(M,NSCN),PGX,RHOGWX,THKGWX )
          CALL THK_G( T_GS(M,NSCN),THKGAX,THKGWX,XMGA_GS(M,NSCN),
     &      XMGW_GS(M,NSCN),THKG_GS(M,NSCN) )
!
!---      Aqueous enthalpy and internal energy  ---
!
          HL_GS(M,NSCN) = MAX(1.D+0-XLA_GS(M,NSCN),0.D+0)*HLW_GS(M,NSCN)
     &      + XLA_GS(M,NSCN)*HGA_GS(M,NSCN)
!
!---      Aqueous thermal conductivity  ---
!
          CALL THK_W( T_GS(M,NSCN),PX,RHOLWX,THKL_GS(M,NSCN) )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_SFC group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE QROMB( FUNC,A,B,SSX,IERR,INDX )
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
!     Returns as S the integral of the function FUNC from A to B.
!     Integration is performed by Romberg's method of order 2K, where
!     e.g., K=2 is Simpson's rule.
!
!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!     1986.  Numerical Recipes, The Art of Scientific Computing.
!     Cambridge University Press, Cambridge.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, November 19, 1999.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Parameter Statements----------------------------!
!
      PARAMETER( EPS=1.D-6, EPSL=1.D-14 )
      PARAMETER( JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1 )
!
!----------------------Type Declarations-------------------------------!
!
      EXTERNAL FUNC
      REAL*8 S(JMAXP),H(JMAXP)
!
!----------------------Executable Lines--------------------------------!
!
      IERR = 0
      ZERO = 0.D+0
!
!---  S and H store the successive trapezodial approximations and their
!     relative step-sizes.  ---
!
      H(1) = 1.D+0
      DO 10 J = 1,JMAX
        CALL TRAPZD( FUNC,A,B,S(J),J,INDX )
        IF( J.GE.K ) THEN
          CALL POLINT( H(J-KM),S(J-KM),K,ZERO,SSX,DSS,IERR )
          IF( (ABS(DSS)-EPS*ABS(SSX))/EPSL.LT.EPSL ) RETURN
        ENDIF
        S(J+1) = S(J)
        H(J+1) = 2.5D-1*H(J)
   10 CONTINUE
      IERR = 1
!
!---  End of QROMB group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDATMOS
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
!     Read input file for atmospheric conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 June 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PLT_ATM
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,CDUM,FDUM,FMDUM,UNTS
      CHARACTER*512 CHDUM
      LOGICAL FCHK
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDATMOS'
!
!---  Read atmospheric data start time  ---
!
      CARD = 'Atmospheric Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Atmospheric Start Time: Month'
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      VARB = 'Atmospheric Start Time: Day'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IDAY)
      IF( INDEX(ADUM(1:),'jan').NE.0 .OR.
     &  INDEX(ADUM(1:),'01').NE.0 ) THEN
        ATMST = ATMST + 0.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'feb').NE.0 .OR.
     &  INDEX(ADUM(1:),'02').NE.0 ) THEN
        ATMST = ATMST + 31.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.28 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'mar').NE.0 .OR.
     &  INDEX(ADUM(1:),'03').NE.0 ) THEN
        ATMST = ATMST + 59.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'apr').NE.0 .OR.
     &  INDEX(ADUM(1:),'04').NE.0 ) THEN
        ATMST = ATMST + 90.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.30 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'may').NE.0 .OR.
     &  INDEX(ADUM(1:),'05').NE.0 ) THEN
        ATMST = ATMST + 120.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'jun').NE.0 .OR.
     &  INDEX(ADUM(1:),'06').NE.0 ) THEN
        ATMST = ATMST + 151.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.30 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'jul').NE.0 .OR.
     &  INDEX(ADUM(1:),'07').NE.0 ) THEN
        ATMST = ATMST + 181.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'aug').NE.0 .OR.
     &  INDEX(ADUM(1:),'08').NE.0 ) THEN
        ATMST = ATMST + 212.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'sep').NE.0 .OR.
     &  INDEX(ADUM(1:),'09').NE.0 ) THEN
        ATMST = ATMST + 243.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.30 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'oct').NE.0 .OR.
     &  INDEX(ADUM(1:),'10').NE.0 ) THEN
        ATMST = ATMST + 273.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'nov').NE.0 .OR.
     &  INDEX(ADUM(1:),'11').NE.0 ) THEN
        ATMST = ATMST + 304.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.30 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'dec').NE.0 .OR.
     &  INDEX(ADUM(1:),'12').NE.0 ) THEN
        ATMST = ATMST + 334.D+0*8.640D+4
        IF( IDAY.GE.1 .OR. IDAY.LE.31 ) THEN
          REALX = REAL(IDAY-1)
          ATMST = ATMST + REALX*8.640D+4
        ELSE
          INDX = 4
          IMSG = IDAY
          CHMSG = 'Out-of-Range Atmospheric Start Day: '
          CALL WRMSGS( INDX )
        ENDIF
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Atmospheric Start Month: '
     &    // ADUM(1:NCHA)
        CALL WRMSGS( INDX )
      ENDIF
      VARB = 'Atmospheric Start Time: Year'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IYEAR)
      IF( IYEAR.GE.1900 ) THEN
        REALX = REAL(IYEAR-1900)
        ATMST = ATMST + REALX*365.D+0*8.640D+4
      ELSE
        INDX = 4
        IMSG = IYEAR
        CHMSG = 'Out-of-Range Atmospheric Start Year: '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Correct for leap years  ---
!
      DO 5 M = 1901,IYEAR
        IF( MOD(M,1000).EQ.0 ) THEN
          IF( MOD(M,400).EQ.0 ) THEN
            ATMST = ATMST + 8.64D+4
          ENDIF
        ELSEIF( MOD(M,4).EQ.0 ) THEN
          ATMST = ATMST + 8.64D+4
        ENDIF
    5 CONTINUE
!
!---  Read atmospheric start time in military format  ---
!
      VARB = 'Atmospheric Start Time: Time'
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      READ( ADUM(1:2),'(I2)' ) IHOUR
      IF( IHOUR.GE.0 .OR. IHOUR.LT.24 ) THEN
        REALX = REAL(IHOUR)
        ATMST = ATMST + REALX*3600.D+0
      ELSE
        INDX = 4
        IMSG = IHOUR
        CHMSG = 'Out-of-Range Atmospheric Start Hour: '
        CALL WRMSGS( INDX )
      ENDIF
      READ( ADUM(4:5),'(I2)' ) IMIN
      IF( IMIN.GE.0 .OR. IMIN.LT.60 ) THEN
        REALX = REAL(IMIN)
        ATMST = ATMST + REALX*60.D+0
      ELSE
        INDX = 4
        IMSG = IMIN
        CHMSG = 'Out-of-Range Atmospheric Start Minute: '
        CALL WRMSGS( INDX )
      ENDIF
      READ( ADUM(7:8),'(I2)' ) ISEC
      IF( ISEC.GE.0 .OR. ISEC.LT.60 ) THEN
        REALX = REAL(ISEC)
        ATMST = ATMST + REALX
      ELSE
        INDX = 4
        IMSG = ISEC
        CHMSG = 'Out-of-Range Atmospheric Start Second: '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Read measurement height for wind speed  ---
!
      VARB = 'Wind Speed Measurement Height'
      UNTS = 'm'
      IUNM = 1
      IDFLT = 1
      ATMC(1) = 1.D+0
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(1))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',ATMC(1)
      INDX = 0
      CALL RDUNIT(UNTS,ATMC(1),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(1),', m)'
!
!---  Read measurement height for air temperature and
!     relative humidity  ---
!
      VARB = 'Air Temperature/Relative Humidity Measurement Height'
      UNTS = 'm'
      IUNM = 1
      IDFLT = 1
      ATMC(2) = 1.D+0
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(2))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',ATMC(2)
      INDX = 0
      CALL RDUNIT(UNTS,ATMC(2),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(2),', m)'
!
!---  Read local longitude  ---
!
      VARB = 'Local Longitude'
      UNTS = 'deg'
      IDFLT = 1
      ATMC(3) = 0.D+0
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(3))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',ATMC(3)
      INDX = 0
      CALL RDUNIT(UNTS,ATMC(3),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(3),', radians)'
!
!---  Read local latitude  ---
!
      VARB = 'Local Latitude'
      UNTS = 'deg'
      IDFLT = 1
      ATMC(4) = 0.D+0
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(4))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',ATMC(4)
      INDX = 0
      CALL RDUNIT(UNTS,ATMC(4),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(4),', radians)'
!
!---  Read local meridian, where meridians for the United States
!     are as follows:
!     Eastern Standard Time = 75 deg
!     Central Standard Time = 90 deg
!     Mountain Standard Time = 105 deg
!     Pacific Standard Time = 120 deg  ---
!
      VARB = 'Local Meridian'
      UNTS = 'deg'
      IDFLT = 1
      ATMC(5) = 120.D+0
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(5))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',ATMC(5)
      INDX = 0
      CALL RDUNIT(UNTS,ATMC(5),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(5),', radians)'
!
!---  Read roughness height for momentum transport  ---
!
      CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
      IF( INDX.EQ.1 ) THEN
        VARB = 'Roughness Height for Momentum Transport'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        ATMC(6) = 1.3D-3
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(6))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',ATMC(6)
        INDX = 0
        CALL RDUNIT(UNTS,ATMC(6),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(6),', m)'
      ENDIF
!
!---  Read roughness height for heat and mass transport  ---
!
      CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
      IF( INDX.EQ.1 ) THEN
        VARB = 'Roughness Height for Heat and Mass Transport'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        ATMC(7) = 1.3D-3
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMC(7))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',ATMC(7)
        INDX = 0
        CALL RDUNIT(UNTS,ATMC(7),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',ATMC(7),', m)'
      ENDIF
!
!---  Read number of atmospheric condition times  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Atmospheric Condition Times'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NATM_T)
      IF( NATM_T.LE.-3 ) THEN
        IATM_C = 1
        NATM_T = ABS(NATM_T)
        WRITE(IWR,'(A)') 'Cyclic Atmospheric Conditions'
      ELSEIF( NATM_T.GE.1 ) THEN
        IATM_C = 0
        WRITE(IWR,'(A)') 'Noncyclic Atmospheric Conditions'
      ELSEIF( NATM_T.EQ.0 ) THEN
        INDX = 2
        CHMSG = 'No Atmospheric Condition Times'
        CALL WRMSGS( INDX )
      ELSE
        INDX = 4
        CHMSG = 'Number of Cyclic Atmospheric Conditions Times < 3'
        CALL WRMSGS( INDX )
      ENDIF
      IF( NATM_T.GT.LATM ) THEN
        INDX = 5
        CHMSG = 'Number of Atmospheric Condition Times > LATM'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Read number of atmospheric condition times and variables  ---
!
      ACTMO = -SMALL
      WRITE(IWR,'(A)') 'Atmospheric Condition Times and Variables:'
      DO 100 NTM = 1,NATM_T
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
!
!---    Check for external atmospheric condition time file  ---
!
        IF( NTM.EQ.1 ) CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,CDUM)
        IF( INDEX(CDUM(1:),'file').NE.0 ) THEN
          CDUM = ''
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
          NCH = INDEX(FDUM,'  ')-1
          IRD = 32
!
!---      Check for external file  ---
!
          INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
          IF( .NOT.FCHK ) THEN
            INDX = 4
            CHMSG = 'Missing Atmospheric Condition File: ' // 
     &        FDUM(1:NCH)
            CALL WRMSGS( INDX )
          ELSEIF( FDUM.EQ.'unformatted' ) THEN
            INDX = 4
            CHMSG = 'Atmospheric Condition File Format: ' // 
     &        FDUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
          OPEN(UNIT=32,FILE=FDUM(1:NCH),STATUS='OLD',FORM='FORMATTED')
          WRITE(IWR,'(/,2A)') 'Atmospheric Condition Time File: ',
     &      FDUM(1:NCH)
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
        ENDIF
!
!---    Read atmospheric condition time  ---
!
        ISTART = 1
        VARB = 'Atmospheric Condition Time'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &    ', ',UNTS(1:NCH),': ',ATMOS(NTM,1)
        INDX = 0
        IUNS = 1
        CALL RDUNIT(UNTS,ATMOS(NTM,1),INDX)
        IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &    ATMOS(NTM,1),', C)'
!
!---    Read atmospheric condition temperature  ---
!
        VARB = 'Atmospheric Condition Temperature'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,2))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &    ', ',UNTS(1:NCH),': ',ATMOS(NTM,2)
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,ATMOS(NTM,2),INDX)
        IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &    ATMOS(NTM,2),', C)'
!
!---    Read atmospheric condition pressure  ---
!
        VARB = 'Atmospheric Condition Pressure'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,3))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &    ', ',UNTS(1:NCH),': ',ATMOS(NTM,3)
        INDX = 0
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        CALL RDUNIT(UNTS,ATMOS(NTM,3),INDX)
        IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &    ATMOS(NTM,3),', Pa)'
!
!---    Read atmospheric condition water-vapor relative
!       humidity  ---
!
        VARB = 'Atmospheric Condition Water-Vapor Relative Humidity'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,4))
        IF( IRD.EQ.21 ) WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),
     &    ': ',ATMOS(NTM,4)
!
!---    Read atmospheric condition net solar radiation  ---
!
        VARB = 'Atmospheric Condition Net Solar Radiation'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,5))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &    ', ',UNTS(1:NCH),': ',ATMOS(NTM,5)
        INDX = 0
        IUNKG = 1
        IUNS = -3
        CALL RDUNIT(UNTS,ATMOS(NTM,5),INDX)
        IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &    ATMOS(NTM,5),', W/m^2)'
!
!---    Read atmospheric condition wind speed  ---
!
        VARB = 'Atmospheric Condition Wind Speed'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,6))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &    ', ',UNTS(1:NCH),': ',ATMOS(NTM,6)
        INDX = 0
        IUNM = 1
        IUNS = -1
        CALL RDUNIT(UNTS,ATMOS(NTM,6),INDX)
        IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &    ATMOS(NTM,6),', m/s)'
!
!---    Read atmospheric condition precipitaton  ---
!
        VARB = 'Atmospheric Condition Precipitation'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ATMOS(NTM,7))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        IF( IRD.EQ.21 ) WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),
     &    ', ',UNTS(1:NCH),': ',ATMOS(NTM,7)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,ATMOS(NTM,7),INDX)
        IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &    ATMOS(NTM,7),', m)'
  100 CONTINUE
!
!---    Close boundary condition time file  ---
!
        IF( IRD.EQ.32 ) THEN
          IRD = 21
          CLOSE(32)
        ENDIF
!      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDATMOS group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDPLANT
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
!     Read input file for plant property information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Gene Freeman, PNNL, 9 April 2002.
!
!-------------------List of input parameters and units ----------------!
!
!     Maximum extent in Z (Zm), m
!     Maximum extent in X (Xm), m
!     Maximum extent in Y (Ym), m
!     Fit parameter in Z (z*), m
!     Fit parameter in X (x*), m
!     Fit parameter in Y (y*), m
!     Leaf area index, unitless
!     Plant canopy height (PCH), m
!     Water stress head 1 (h1), m
!     Water stress head 2 (h2), m
!     Water stress head 3 (h3), m
!     Water stress head 4 (h4), m
!     Water uptake reduced by 50% (h50), m
!     Crop coefficient - stage1, unitless
!     Crop start time jan.1=1, dec.31=365, day
!     Crop coefficient - stage2, unitless
!     Crop start time jan.1=1, dec.31=365, day
!     Crop coefficient - stage3, unitless
!     Crop start time jan.1=1, dec.31=365, day
!     Crop coefficient - stage4, unitless
!     Crop start time jan.1=1, dec.31=365, day
!     Evapotranspiration (ET), m/day
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS,ADUM
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDPLANT'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Plant Properties Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of plant varietals  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Plants'
      CALL RDINT( ISTART,ICOMMA,CHDUM,NPLANT )
      IF( NPLANT.GT.LPLANT ) THEN
        INDX = 5
        CHMSG = 'Number of Plant Varietals > Parameter LPLANT'
        CALL WRMSGS( INDX )
      ENDIF
      IF( NPLANT.GT.(LBCV-2) ) THEN
        INDX = 5
        CHMSG = 'Number of Plant Varietals > Parameter (LBCV-2)'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Check for multiple plant temperature option and
!     rainfall interception-condensation shedding option  ---
!
      CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
      IF( INDX.EQ.1 ) THEN
        VARB = 'Plant Temperature and Rainfall Interception Options'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!        IF( INDEX(ADUM(1:),'multiple').NE.0 ) THEN
!          ISLC(24) = ISLC(24)+2
!          WRITE(IWR,'(A)') '  Multiple Plant Temperature Option'
!        ELSE
          ISLC(24) = ISLC(24)+1
          WRITE(IWR,'(A)') '  Single Plant Temperature Option'
!        ENDIF
        IF( INDEX(ADUM(1:),'rainfall').NE.0 .OR.
     &    INDEX(ADUM(1:),'interception').NE.0 ) THEN
          ISLC(26) = 1
          WRITE(IWR,'(A)') '  Rainfall Interception'
        ELSE
          WRITE(IWR,'(A)') '  No Rainfall Interception'
        ENDIF
      ELSE
        ISLC(24) = ISLC(24)+1
        WRITE(IWR,'(A)') '  Single Plant Temperature Option'
        WRITE(IWR,'(A)') '  No Rainfall Interception'
      ENDIF
!
!---  Loop over the plants information lines  ---
!
      DO 500 IP = 1,NPLANT
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Plant Name: '
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,PLANT(IP))
        DO 100 M = 1,IP-1
          IF( PLANT(M).EQ.PLANT(IP) ) THEN
            INDX = 4
            CHMSG = 'Duplicate Plant Name: ' // PLANT(IP)
            CALL WRMSGS( INDX )
          ENDIF
  100   CONTINUE
        WRITE (IWR,'(/,2A)') ' Plant Name: ',PLANT(IP)
!
!---    Check for root stress or stomatal resistance options  ---
!
        CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
        IF( INDX.EQ.1 ) THEN
          VARB = 'Plant Options'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'vrugt').NE.0 ) THEN
            IRSM_P(IP) = 1
            WRITE(IWR,'(A)') '  Vrugt Root Stress Model'
          ELSEIF( INDEX(ADUM(1:),'jarvis').NE.0 ) THEN
            IRSM_P(IP) = 2
            WRITE(IWR,'(A)') '  Jarvis Root Stress Model'
          ELSEIF( INDEX(ADUM(1:),'stress').NE.0 ) THEN
            IRSM_P(IP) = 1
            WRITE(IWR,'(A)') '  Vrugt Root Stress Model'
          ELSE
            WRITE(IWR,'(A)') '  No Root Stress Model'
          ENDIF
          IF( INDEX(ADUM(1:),'hicks').NE.0 ) THEN
            ISRM_P(IP) = 1
            WRITE(IWR,'(A)') '  Hicks Stomatal Resistance Model'
          ELSE
            WRITE(IWR,'(A)') '  No Stomatal Resistance Model'
          ENDIF
        ELSE
          WRITE(IWR,'(A)') '  No Root Stress Model'
          WRITE(IWR,'(A)') '  No Stomatal Resistance Model'
        ENDIF
!
!---    Read root (Z) depth characteristics  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Max. Root Depth'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(1,IP) = 1.D-3
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(1,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(1,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(1,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(1,IP),', m)'
!
!---    Minimum maximum plant root depth of 1 mm  ---
!
        PARMS_P(1,IP) = MAX( 1.D-3,PARMS_P(1,IP) )
!
!---    Read root z* characteristics  ---
!
        VARB = 'Null Root Depth'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(4,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(4,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(4,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(4,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(4,IP),', m)'
!
!---    Read root pz characteristics  ---
!
        VARB = 'Root Depth fit parameter'
        IDFLT = 1
        PARMS_P(7,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(7,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(7,IP)
!
!---    Read plant short-wave albedo  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Plant Solar Albedo'
!
!---    Check for temporal input to be read in conjunction with crop
!       coefficient times  ---
!
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        IF( INDEX(ADUM,'temporal').NE.0 ) THEN
          IALB_P(IP) = 1
!
!---      Read plant solar albedo (initial stage start)  ---
!
          VARB = 'Plant Solar Albedo: Initial Stage Start'
          IDFLT = 1
          PARMS_P(5,IP) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(5,IP))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(5,IP)
!
!---      Read plant solar albedo (crop development start)  ---
!
          VARB = 'Plant Solar Albedo: Crop Development Start'
          IDFLT = 1
          PARMS_P(6,IP) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(6,IP))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &        ': ',PARMS_P(6,IP)
!
!---      Read plant solar albedo (mid-season start)  ---
!
          VARB = 'Plant Solar Albedo: Mid-Season Start'
          IDFLT = 1
          PARMS_P(8,IP) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(8,IP))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &        ': ',PARMS_P(8,IP)
!
!---      Read plant solar albedo (late-season start)  ---
!
          VARB = 'Plant Solar Albedo: Late-Season Start'
          IDFLT = 1
          PARMS_P(9,IP) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(9,IP))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &        ': ',PARMS_P(9,IP)
!
!---      Read plant solar albedo (late-season stop)  ---
!
          VARB = 'Plant Solar Albedo: Late-Season Stop'
          IDFLT = 1
          PARMS_P(10,IP) = PARMS_P(5,IP)
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(10,IP))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &        ': ',PARMS_P(10,IP)
        ELSE
          ISTART = 1
          IDFLT = 1
          PARMS_P(10,IP) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(10,IP))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',PARMS_P(10,IP)
        ENDIF
!
!---    Read plant height characteristics  ---
!
        VARB = 'Plant Canopy Height'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(11,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(11,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(11,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(11,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(11,IP),', m)'
!
!---    Minimum maximum plant canopy height of 1 mm  ---
!
        PARMS_P(11,IP) = MAX( 1.D-3,PARMS_P(11,IP) )
!
!---    Maximum condensate depth  ---
!
        IF( ISLC(26).EQ.1 ) THEN
          VARB = 'Maximum Condensate Depth'
          UNTS = 'm'
          IUNM = 1
          IDFLT = 1
          PARMS_P(16,IP) = 0.2D-3
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(16,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PARMS_P(16,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(16,IP),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(16,IP),', m)'
        ENDIF
!
!---  Vrugt root stress model
!
      IF( IRSM_P(IP).EQ.1 ) THEN
!
!---  Read first stress point head  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Root Stress-Point 1 (Capillary Head)'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(12,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(12,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(12,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(12,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(12,IP),', m)'
!
!---  Read second stress point head  ---
!
        VARB = 'Root Stress-Point 2 (Capillary Head)'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(13,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(13,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(13,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(13,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(13,IP),', m)'
!
!---  Read third stress point head  ---
!
        VARB = 'Root Stress-Point 3 (Capillary Head)'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(14,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(14,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(14,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(14,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(14,IP),', m)'
!
!---  Read fourth stress point head  ---
!
        VARB = 'Root Stress-Point 4 (Capillary Head)'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(15,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(15,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(15,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(15,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(15,IP),', m)'
!
!---  Jarvis root stress model
!
      ELSEIF( IRSM_P(IP).EQ.2 ) THEN
!
!---  Read first Jarvis stress point  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Wilting-Point Water Content'
        IDFLT = 1
        PARMS_P(12,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(12,IP))
        WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PARMS_P(12,IP)
!
!---  Read second Jarvis stress point  ---
!
        VARB = 'Normalized Soil Water Content: Critical Point 1'
        IDFLT = 1
        PARMS_P(13,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(13,IP))
        WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PARMS_P(13,IP)
!
!---  Read third Jarvis stress point  ---
!
        VARB = 'Normalized Soil Water Content: Critical Point 2'
        IDFLT = 1
        PARMS_P(14,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(14,IP))
        WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PARMS_P(14,IP)
!
!---  Read fourth Jarvis stress point  ---
!
        VARB = 'Saturated Water Content'
        IDFLT = 1
        PARMS_P(15,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(15,IP))
        WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PARMS_P(15,IP)
      ELSE
!
!---  Read head when uptake reduce 50%  ---
!
       CALL RDINPL( CHDUM )
       CALL LCASE( CHDUM )
       ISTART = 1
       VARB = 'root uptake reduced 50%'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(12,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(12,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(12,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(12,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(12,IP),', m)'
      ENDIF
!
!---  Read crop coefficient (initial stage start)  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Crop Coefficient: Initial Stage Start'
        IDFLT = 1
        PARMS_P(17,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(17,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(17,IP)
!
!---  Read crop coefficient day of year (initial stage start)  ---
!
        VARB = 'Initial Stage Start Day'
        UNTS = 'day'
        IUNS = 1
        IDFLT = 1
        PARMS_P(18,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(18,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(18,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(18,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(18,IP),', s)'
!
!---  Read crop coefficient (crop development start)  ---
!
        VARB = 'Crop Coefficient: Crop Development Start'
        IDFLT = 1
        PARMS_P(19,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(19,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(19,IP)
!
!---  Read crop coefficient day of year (crop development start)  ---
!
        VARB = 'Crop Development Start Day'
        UNTS = 'day'
        IUNS = 1
        IDFLT = 1
        PARMS_P(20,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(20,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(20,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(20,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(20,IP),', s)'
!
!---  Read crop coefficient (mid-season start)  ---
!
        VARB = 'Crop Coefficient: Mid-Season Start'
        IDFLT = 1
        PARMS_P(21,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(21,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(21,IP)
!
!---  Read crop coefficient day of year (mid-season start)  ---
!
        VARB = 'Mid-Season Start Day'
        UNTS = 'day'
        IUNS = 1
        IDFLT = 1
        PARMS_P(22,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(22,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(22,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(22,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(22,IP),', s)'
!
!---  Read crop coefficient (late-season start)  ---
!
        VARB = 'Crop Coefficient: Late-Season Start'
        IDFLT = 1
        PARMS_P(23,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(23,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(23,IP)
!
!---  Read crop coefficient day of year (late-season start)  ---
!
        VARB = 'Late-Season Start Day'
        UNTS = 'day'
        IUNS = 1
        IDFLT = 1
        PARMS_P(24,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(24,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(24,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(24,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(24,IP),', s)'
!
!---  Read crop coefficient (late-season stop)  ---
!
        VARB = 'Crop Coefficient: Late-Season Stop'
        IDFLT = 1
        PARMS_P(25,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(25,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(25,IP)
!
!---  Read crop coefficient day of year (late-season stop)  ---
!
        VARB = 'Late-Season Stop Day'
        UNTS = 'day'
        IUNS = 1
        IDFLT = 1
        PARMS_P(26,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(26,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(26,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(26,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(26,IP),', s)'
!
!---    Hicks stomatal resistance model  ---
!
        IF( ISRM_P(IP).EQ.1 ) THEN
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
!
!---      Read minimum stomatal resistance, s/m  ---
!
          VARB = 'Stomatal Resistance: Minimum Stomatal Resistance'
          UNTS = 's/m'
          IUNM = -1
          IUNS = 1
          IDFLT = 1
          PARMS_P(2,IP) = 100.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(2,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PARMS_P(2,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(2,IP),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(2,IP),', s/m)'
!
!---      Read light response coefficient, W/m^2  ---
!
          VARB = 'Stomatal Resistance: Light Response Coefficient'
          UNTS = 'w/m^2'
          IUNKG = 1
          IUNS = -3
          IDFLT = 1
          PARMS_P(3,IP) = 20.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(3,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PARMS_P(3,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(3,IP),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(3,IP),', W/m^2)'
!
!---      Read minimum temperature for stomatal opening, C  ---
!
          VARB = 'Stomatal Resistance: Minimum Temperature'
          UNTS = 'c'
          IUNK = 1
          IDFLT = 1
          PARMS_P(27,IP) = 5.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(27,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PARMS_P(27,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(27,IP),INDX)
          PARMS_P(27,IP) = PARMS_P(27,IP)+TABS
!
!---      Read maximum temperature for stomatal opening, C  ---
!
          VARB = 'Stomatal Resistance: Maximum Temperature'
          UNTS = 'c'
          IUNK = 1
          IDFLT = 1
          PARMS_P(28,IP) = 45.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(28,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PARMS_P(28,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(28,IP),INDX)
          PARMS_P(28,IP) = PARMS_P(28,IP)+TABS
!
!---      Read optimum temperature for stomatal opening, C  ---
!
          VARB = 'Stomatal Resistance: Optimum Temperature'
          UNTS = 'c'
          IUNK = 1
          IDFLT = 1
          PARMS_P(29,IP) = 25.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(29,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PARMS_P(29,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(29,IP),INDX)
          PARMS_P(29,IP) = PARMS_P(29,IP)+TABS
        ENDIF
!
!---  Read next plant type  ---
!
        IF( IP.LT.NPLANT ) WRITE(IWR,'(/)')
 500  CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDPLNT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSFCOV
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
!     Read input file for surface cover card.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK White, PNNL, 15 September 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PLT_ATM
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,FDUM,UNTS,UNTSPX
      CHARACTER*512 CHDUM
      REAL*8 XSPX(LSFCP),YSPX(LSFCP)
      INTEGER ITMPX(LFX,LFY),NTMPX(LFX,LFY)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSFCOV'
!
!---  Write card information to output file  ---
!
      CARD = 'Surface Cover Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Initialize arrays  ---
!
      DO J = 1,JFLD
        DO I = 1,IFLD
          ITMPX(I,J) = 0
        ENDDO
      ENDDO
!
!---  Read the number of surface cover areas ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      NSFCA = 0
      VARB = 'Number of Surface Cover Areas'
      CALL RDINT( ISTART,ICOMMA,CHDUM,NSFCA )
!
!---  Loop over number of surface cover areas  ---
!
      DO 500 NSC = 1,NSFCA
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Surface Cover Area Name: '
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,NAME_SFC(NSC))
        WRITE(IWR,'(/,2A)') 'Surface Cover Area Name: ',NAME_SFC(NSC)
!
!---    Read number of surface cover polygons corresponding to 
!       surface cover area  ---
!
        VARB = 'Number of Surface Cover Polygons'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NPLYX)
        IF( NPLYX.EQ.0 ) NPLYX = 1
!
!---    Maximum ponding height for surface cover area  ---
!
        IDFLT = 1
        PHMX(NSC) = 1.D+2
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PHMX(NSC))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',PHMX(NSC)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,PHMX(NSC),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PHMX(NSC),', m)'
!
!---    Check for ground-surface albedo model  ---
!
        CALL CHKCHR(ISTART,ICOMMA,CHDUM,INDX)
        IF( INDX.EQ.1 ) THEN
          VARB = 'Ground-Surface Albedo'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        ENDIF
!
!---    Ground-surface albedo  ---
!
        IF( INDEX(ADUM(1:),'albedo').NE.0 ) THEN
!
!---      Ground-surface solar angle model  ---
!
          VARB = 'Ground-Surface Solar Angle Model'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'plem').NE.0 .OR.
     &      INDEX(ADUM(1:),'xiu').NE.0 ) THEN
            IALB(NSC) = 1
          ELSEIF( INDEX(ADUM(1:),'wang').NE.0 ) THEN
            IALB(NSC) = 2
          ELSEIF( INDEX(ADUM(1:),'briegleb').NE.0 ) THEN
            IALB(NSC) = 3
          ELSEIF( INDEX(ADUM(1:),'moisture').NE.0 ) THEN
            IALB(NSC) = 4
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Ground-Surface Solar Angle Model: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Ground-surface dry-soil albedo  ---
!
          IDFLT = 1
          VARB = 'Dry-Soil Albedo'
          UNTS = 'null'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ALBEDO(1,NSC))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',ALBEDO(1,NSC)
!
!---      Ground-surface wet-soil albedo  ---
!
          IDFLT = 1
          VARB = 'Wet-Soil Albedo'
          UNTS = 'null'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ALBEDO(2,NSC))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',ALBEDO(2,NSC)
!
!---      Ground-surface albedo attenuation factor  ---
!
          IDFLT = 1
          VARB = 'Albedo Attenuation Factor'
          UNTS = 'null'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,ALBEDO(3,NSC))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',ALBEDO(3,NSC)
!
!---      Wang ground-surface solar angle model  ---
!
          IF( IALB(NSC).EQ.2 ) THEN
            IDFLT = 1
            VARB = 'Reference Albedo @ Solar Zenith = 60 deg'
            UNTS = 'null'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,ALBEDO(4,NSC))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        ALBEDO(4,NSC)
!
!---      Briegleb ground-surface solar angle model  ---
!
          ELSEIF( IALB(NSC).EQ.3 ) THEN
            IDFLT = 1
            VARB = 'Reference Albedo @ Solar Zenith = 60 deg'
            UNTS = 'null'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,ALBEDO(4,NSC))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        ALBEDO(4,NSC)
            IDFLT = 1
            VARB = 'Parameter C'
            UNTS = 'null'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,ALBEDO(5,NSC))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        ALBEDO(5,NSC)
          ENDIF
        ENDIF
!
!---    Loop over number of surface cover polygons  ---
!
        DO 400 NSCP = 1,NPLYX
!
!---      Read number of surface cover polygon definition x,y pairs  ---
!
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Number of Surface Cover Polygon Definition Pairs'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NSFCPX)
          WRITE(IWR,'(2X,A,I4)') 'Surface Cover Polygon: ',NSCP
          WRITE(IWR,'(2X,A,I9)') 'Number of Surface Cover Polygon' //
     &      ' Definition Pairs:',NSFCPX
!
!---      Read units for surface cover polygon coordinate units  ---
!
          VARB = 'Surface Cover Polygon Coordinate Units'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          UNTSPX = UNTS
          VARX = 1.D+0
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VARX,INDX)
!
!---      Read surface cover polygon coordinates  ---
!
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Surface Cover Polygon Points'
          WRITE(IWR,'(2X,2A)')  VARB(1:IVR),': '
          DO 300 I = 1,NSFCPX   
            CALL RDDPR(ISTART,ICOMMA,CHDUM,XSPX(I))
            CALL RDDPR(ISTART,ICOMMA,CHDUM,YSPX(I))
            XTMPX = XSPX(I)*VARX
            YTMPX = YSPX(I)*VARX
            XSPX(I) = XSPX(I)*VARX
            YSPX(I) = YSPX(I)*VARX
            WRITE(IWR,'(2X,1PE11.4,A,1PE11.4,3A,1PE11.4,A,1PE11.4,'//
     &        '2A)') XTMPX,', ',YTMPX,', ',UNTSPX(1:NCH),
     &       ' (',XSPX(I),',',YSPX(I),',',' m)'
  300     CONTINUE
!
!---      Check to make sure polygon is closed. If not, close it.  ---
!         
          IF ( XSPX(1).NE.XSPX(NSFCPX) .OR. 
     &       YSPX(1).NE.YSPX(NSFCPX) ) THEN
            NSFCPX = NSFCPX + 1
            XSPX(NSFCPX) = XSPX(1)
            YSPX(NSFCPX) = YSPX(1)
          ENDIF          
!
!---      Loop over nodes from top down to determine which nodes are  
!         active surface nodes. If the top node is inactive, find the  
!         top most active node. If node is an active surface node check 
!         to see if it is within the surfacecover area polygon.     
!   
          DO J = 1,JFLD
          DO I = 1,IFLD
            DO K = KFLD,1,-1
              N = ND(I,J,K)           
              IF( IXP(N).NE.0 ) THEN
                XPT = 2.5D-1*(XE(5,N)+XE(6,N)+XE(7,N)+XE(8,N))
                YPT = 2.5D-1*(YE(5,N)+YE(6,N)+YE(7,N)+YE(8,N))             
                CALL LOCPT( XPT,YPT,XSPX,YSPX,NSFCPX,IINOUT,IPATH )
                IF ( IINOUT.GE.0 ) THEN
                  IF( ITMPX(I,J).NE.0 ) THEN
                    INDX = 24
                    CHMSG = 'Surface Cover Polygon Overlap at Node'
                    IMSG = N
                  ENDIF
                  ITMPX(I,J) = NSC
                  NTMPX(I,J) = N
                ENDIF     
                EXIT
              ENDIF
            ENDDO
          ENDDO
          ENDDO
  400   CONTINUE
!
!---    Read number of surface cover times  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Surface Cover Times'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NSFCT(NSC))
        WRITE(IWR,'(2X,A)') VARB(1:IVR)
!
!---    Read, write, and convert surface cover time, variables,
!       and units  ---
!
        DO 450 NTM = 1,NSFCT(NSC)
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_SFC(1,NTM,NSC))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH)
     &      ,': ',PARMS_SFC(1,NTM,NSC)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,PARMS_SFC(1,NTM,NSC),INDX)
          PAI_PX = 0.D+0
          DO 425 IPLANT = 1,NPLANT
            NCH = INDEX( PLANT(IPLANT)(1:),'  ' ) - 1
            VARB = 'Leaf Area Index: ' // PLANT(IPLANT)(1:NCH)
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,FDUM )
            IPLX = (IPLANT-1)*2 + 2
            ISTART = ISX
            ICOMMA = ICX
            CALL RDDPR( ISTART,ICOMMA,CHDUM,PARMS_SFC(IPLX,NTM,NSC) )
            IF( IRD.EQ.21 ) WRITE(IWR,'(2X,2A,1PE11.4,A)') 
     &        VARB(1:IVR),': ',PARMS_SFC(IPLX,NTM,NSC),
     &        ' m^2 Green Leaf/m^2 Plant Surface'
            NCH = INDEX( PLANT(IPLANT)(1:),'  ' ) - 1
            VARB = 'Plant Area Index: ' // PLANT(IPLANT)(1:NCH)
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,FDUM )
            IPLX = (IPLANT-1)*2 + 3
            ISTART = ISX
            ICOMMA = ICX
            CALL RDDPR( ISTART,ICOMMA,CHDUM,PARMS_SFC(IPLX,NTM,NSC) )
            IF( IRD.EQ.21 ) WRITE(IWR,'(2X,2A,1PE11.4,A)') 
     &        VARB(1:IVR),': ',PARMS_SFC(IPLX,NTM,NSC),
     &        ' m^2 Plant Surface/m^2 Ground Surface'
            PAI_PX = PAI_PX + PARMS_SFC(IPLX,NTM,NSC)
 425      CONTINUE
          IF( PAI_PX.GT.1.D+0 ) THEN
            INDX = 14
            N_DB = N
            CHMSG = 'Combined Plant Area Index Greater Than One: '
            RLMSG = PAI_PX
            CALL WRMSGS( INDX )
          ENDIF
 450    CONTINUE
 500  CONTINUE
!
!---  Loop over nodes in xy plane and fill top node for surface
!     cover connection map index and surface cover area index
! 
      NC = 0        
      DO NSC = 1,NSFCA
        DO J = 1,JFLD
        DO I = 1,IFLD
          IF( ITMPX(I,J).EQ.NSC ) THEN
            NC = NC + 1
            ID_SFC(NC) = NSC
            ICM_SFC(1,NC) = NTMPX(I,J)
          ENDIF
        ENDDO
        ENDDO
      ENDDO
      NSFCN = NC
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSFCOV group ---
!
      RETURN
      END
 
 !----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ROOTS( UPTAKE,ZSX,ZEX,IP,N,M )
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
!     Stressed root-water uptake.
!
!----------------------Authors-----------------------------------------!
!
!     Written by EJ Freeman and MD White, PNNL, 20 June 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
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
      EXTERNAL VRUGT
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ROOTS'
!
!---  Set plant variables  ---
!
      ZM = PARMS_P(1,IP)
!      ZSX = PARMS_P(4,IP)
      PZ = PARMS_P(7,IP)
      H1 = PARMS_P(12,IP)
      H2 = PARMS_P(13,IP)
      H3 = PARMS_P(14,IP)
      H4 = PARMS_P(15,IP)
!
!---  Vrugt root-water uptake model  ---
!
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NK = ND(I,J,KFLD)
      CALL QROMB( VRUGT,ZSX,ZEX,RWUX,IERR,IP )
      IF( IERR.EQ.1 ) THEN
        INDX = 12
        IMSGX = ND(I,J,K)
        NMSGX = 0
        SUBLOGX = 'ROOTS'
        CHMSGX = 'Unconverged Romberg Integration: ' //
     &    'Source Root-Water Uptake Integration: Node: '
        RLMSGX = 0.D+0
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Vrugt water stress function  ---
!
      IF( IRSM_P(IP).EQ.1 ) THEN
        HDGL = (PG(M,N)-PL(M,N))/(RHORL*GRAV)
        H1 = PARMS_P(12,IP)
        H2 = PARMS_P(13,IP)
        H3 = PARMS_P(14,IP)
        H4 = PARMS_P(15,IP)
        IF( HDGL.GE.H4 .OR. HDGL.LE.H1 )THEN
          WSFX = 0.D+0
        ELSEIF( HDGL.LE.H3 .AND. HDGL.GE.H2 )THEN
          WSFX = 1.D+0
        ELSEIF( HDGL.GT.H1 .AND. HDGL.LT.H2 )THEN
          WSFX = (1.D+0/(H2-H1))*(HDGL-H1)
        ELSE
          WSFX = (-1.D+0/(H4-H3))*(HDGL-H3) + 1.D+0
        ENDIF
!
!---  Jarvis water stress function  ---
!
      ELSEIF( IRSM_P(IP).EQ.2 ) THEN
        WC1 = PARMS_P(12,IP)
        WC2 = PARMS_P(13,IP)
        WC3 = PARMS_P(14,IP)
        WC4 = PARMS_P(15,IP)
        WCX = SL(M,N)*PORD(M,N)
        WSFX = (WCX-WC1)/(WC4-WC1+SMALL)
        IF( WSFX.GT.WC2 .AND. (WSFX-EPSL).LE.1.D+0 ) THEN
          WSFX = (1.D+0-WSFX)/(1.D+0-WC2+SMALL)
        ELSEIF( WSFX.GE.WC1 .AND. WSFX.LE.WC2 ) THEN
          WSFX = 1.D+0
        ELSEIF( (WSFX+EPSL).GE.0.D+0 .AND. WSFX.LT.WC1 ) THEN
          WSFX = WSFX/(WC1+SMALL)
        ELSE
          INDX = 12
          IMSGX = ND(I,J,K)
          NMSGX = 0
          SUBLOGX = 'ROOTS'
          CHMSGX = 'Out of Range Stress Index: Node: '
          RLMSGX = 0.D+0
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF
      UPTAKE = MAX( WSFX*RWUX,0.D+0 )
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ROOTS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRAPZD( FUNC,A,B,S,N,INDX )
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
!     This routine computes the N'th stage of refinement of an
!     extended trapezoid rule.
!
!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!     1986.  Numerical Recipes, The Art of Scientific Computing.
!     Cambridge University Press, Cambridge.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, February, 1999.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      EXTERNAL FUNC
      SAVE IT
!
!----------------------Executable Lines--------------------------------!
!
      IF( N.EQ.1 ) THEN
        S = 5.D-1*(B-A)*(FUNC(A,INDX)+FUNC(B,INDX))
!
!---  IT is the number of points to be added on the next call  ---
!
        IT = 1
      ELSE
        REALX = REAL(IT)
        TNM = REALX
!
!---  Spacing of the points to be added.  ---
!
        DEL = (B-A)/TNM
        X = A + 5.D-1*DEL
        SUM = 0.D+0
        DO 100 J = 1,IT
          SUM = SUM + FUNC(X,INDX)
          X = X + DEL
  100   CONTINUE
!
!---  Replace S by its refined value  ---
!
        S = 5.D-1*(S+(B-A)*SUM/TNM)
        IT = 2*IT
      ENDIF
!
!---  End of TRAPZD group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION VRUGT( ZX,IP )
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
!     Vrugt one-dimensional root water uptake model function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 August 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VRUGT'
!
!---  Set maximum rooting depth, null root depth, and root depth
!     fitting parameter  ---
!
      ZM = PARMS_P(1,IP)
      ZSX = PARMS_P(4,IP)
      PZ = PARMS_P(7,IP)
      VRUGT = (1.D+0 - ZX/(ZM+1.D-20))*EXP(-(PZ/ZM)*ABS(ZSX-ZX))
!
!---  End of VRUGT group
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END


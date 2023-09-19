!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_BCV
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
!     Allocate global array memory for solution control variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_BCV'
!
!---  Allocate memory for BC  ---
!
      ALLOCATE( BC(1:LBCV,1:LBTM,1:LBCIN),STAT=ISTAT )
      CHMSG = 'BC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      LBCX = MAX( NBC(ID+1),1 )
      ALLOCATE( XPBC(1:LBCX),STAT=ISTAT )
      CHMSG = 'XPBC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YPBC(1:LBCX),STAT=ISTAT )
      CHMSG = 'YPBC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ZPBC(1:LBCX),STAT=ISTAT )
      CHMSG = 'ZPBC'
      ALLOCATE( TSBC(1:LBCIN),STAT=ISTAT )
      CHMSG = 'TSBC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XSBC(1:LBCIN),STAT=ISTAT )
      CHMSG = 'XSBC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YSBC(1:LBCIN),STAT=ISTAT )
      CHMSG = 'YSBC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ZSBC(1:LBCIN),STAT=ISTAT )
      CHMSG = 'ZSBC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCN(1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCD(1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCD'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCM(1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCIN(1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCIN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCC(1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      LX = LUK+LSOLU*LC
      ALLOCATE( IBCT(1:LX,1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCSP(1:LSPBC+1,1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCSP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NBHG(1:LPH,1:LBCX),STAT=ISTAT )
      CHMSG = 'NBHG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_BCV group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_BCVP
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
!     Allocate array memory for general boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
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
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_BCVP'
      LBCX = MAX( NBC(ID+1),1 )
      ALLOCATE( TB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'TB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PLB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'PLB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PGB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'PGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PSOB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'PSOB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RHOGB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'RHOGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RHOLB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'RHOLB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( BTGLB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'BTGLB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XGAB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'XGAB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XGWB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'XGWB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XMGAB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'XMGAB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XMGWB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'XMGWB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XLAB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'XLAB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XLWB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'XLWB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XMLAB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'XMLAB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XMLWB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'XMLWB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RKGB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'RKGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RKLB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'RKLB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PVAB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'PVAB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PVWB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'PVWB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DFLAB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'DFLAB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DFLSB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'DFLSB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SGB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'SGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SLB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'SLB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YLSB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'YLSB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XLSB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'XLSB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XMLSB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'XMLSB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RHOMLB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'RHOMLB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RHOMGB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'RHOMGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DFGAB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'DFGAB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DFGWB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'DFGWB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PORDB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'PORDB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PORTB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'PORTB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TORGB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'TORGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TORLB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'TORLB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( HGAB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'HGAB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( HGB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'HGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( UEGB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'UEGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( UELB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'UELB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( HGWB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'HGWB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( HLB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'HLB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( HLWB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'HLWB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( VISGB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'VISGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( VISLB(1:LSV,1:LBCX),STAT=ISTAT )
      CHMSG = 'VISLB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_BCVP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_COUP_WELL
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
!     CO2 (Carbon Sequestration in Deep Saline Aquifers) Mode
!
!     Allocate array memory for coupled-well variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 13 December 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_COUP_WELL'
      ALLOCATE( DNR_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'DNR_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DP_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'DP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( FF_CW(1:3,1:LN_CW),STAT=ISTAT )
      CHMSG = 'FF_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( FX_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'FX_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( FXW_CW(1:(LUK+2),1:LWN_CW),STAT=ISTAT )
      CHMSG = 'FXW_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( P_CW(1:3,1:LN_CW),STAT=ISTAT )
      CHMSG = 'P_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PF_CW(1:LWF_CW),STAT=ISTAT )
      CHMSG = 'PF_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PL_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'PL_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( Q_CW(1:4,1:LWN_CW),STAT=ISTAT )
      CHMSG = 'Q_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( QM_CW(1:8,1:LN_CW),STAT=ISTAT )
      CHMSG = 'QM_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RS_CW(1:(LUK_CW+1),1:LN_CW),STAT=ISTAT )
      CHMSG = 'RS_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PAR_CW(1:5,1:LWI_CW),STAT=ISTAT )
      CHMSG = 'PAR_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TML_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'TML_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( VAR_CW(1:7+LNGC,1:LWT_CW,1:LN_CW),STAT=ISTAT )
      CHMSG = 'VAR_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( VARC_CW(1:LSOLU_CW,1:LWT_CW,1:LN_CW),STAT=ISTAT )
      CHMSG = 'VARC_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( VARSP_CW(1:LSPC_CW,1:LWT_CW,1:LN_CW),STAT=ISTAT )
      CHMSG = 'VARSP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PLX_CW(1:LWN_CW),STAT=ISTAT )
      CHMSG = 'PLX_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PLY_CW(1:LWN_CW),STAT=ISTAT )
      CHMSG = 'PLY_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PLZ_CW(1:LWN_CW),STAT=ISTAT )
      CHMSG = 'PLZ_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XP_CW(1:2,1:LWN_CW),STAT=ISTAT )
      CHMSG = 'XP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YP_CW(1:2,1:LWN_CW),STAT=ISTAT )
      CHMSG = 'YP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ZP_CW(1:2,1:LWN_CW),STAT=ISTAT )
      CHMSG = 'ZP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XTP_CW(1:2,1:LWI_CW),STAT=ISTAT )
      CHMSG = 'XTP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YTP_CW(1:2,1:LWI_CW),STAT=ISTAT )
      CHMSG = 'YTP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ZTP_CW(1:2,1:LWI_CW),STAT=ISTAT )
      CHMSG = 'ZTP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( WNM_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'WNM_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TF_CW(1:LWN_CW),STAT=ISTAT )
      CHMSG = 'TF_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ICC_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'ICC_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ID_CW(1:10,1:LN_CW),STAT=ISTAT )
      CHMSG = 'ID_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IMP_CW(1:LWTP_CW,1:LN_CW),STAT=ISTAT )
      CHMSG = 'IMP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ITS_CW(1:LWTP_CW,1:LN_CW),STAT=ISTAT )
      CHMSG = 'ITS_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IM_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'IM_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( INV_CW(1:LWN_CW),STAT=ISTAT )
      CHMSG = 'INV_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IREF_CW(1:LVREF),STAT=ISTAT )
      CHMSG = 'IREF_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IS_CW(1:LWI_CW),STAT=ISTAT )
      CHMSG = 'IS_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISOLU_CW(1:LSOLU_CW,1:LN_CW),STAT=ISTAT )
      CHMSG = 'ISOLU_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISPC_CW(1:LSPC_CW,1:LN_CW),STAT=ISTAT )
      CHMSG = 'ISPC_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISOLC_CW(1:LEQC,1:LN_CW),STAT=ISTAT )
      CHMSG = 'ISOLC_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISOLK_CW(1:LEQK,1:LN_CW),STAT=ISTAT )
      CHMSG = 'ISOLK_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IT_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'IT_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NSP_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'NSP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IWFG_CW(1:LWF_CW),STAT=ISTAT )
      CHMSG = 'IWFG_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IWF_CW(1:LWF_CW),STAT=ISTAT )
      CHMSG = 'IWF_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IXPG_CW(1:LWF_CW),STAT=ISTAT )
      CHMSG = 'IXPG_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IWNG_CW(1:LWN_CW),STAT=ISTAT )
      CHMSG = 'IWNG_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IWN_CW(1:LWN_CW),STAT=ISTAT )
      CHMSG = 'IWN_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IWP_CW(1:LWN_CW),STAT=ISTAT )
      CHMSG = 'IWP_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( JM_CW(1:LN_CW),STAT=ISTAT )
      CHMSG = 'JM_CW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_FDVH
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Allocate global array memory for hydrate field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 January 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GRID
      USE GLB_PAR
      USE FDVH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ALLOCATE( SH(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SH'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  End of ALLOC_FDVH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_FDVP
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
!     Allocate global array memory for field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_FDVP'
      ALLOCATE( T(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'T'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PL(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PG(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PN(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PSO(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PSO'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DNR(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'DNR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RHOL(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'RHOL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XLA(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'XLA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XLW(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'XLW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YLS(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'YLS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XMLW(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'XMLW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RKL(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'RKL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TRPGL(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'TRPGL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PVW(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PVW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SG(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SL(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RHOML(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'RHOML'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PORD(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PORD'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PORT(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PORT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PERMRF(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PERMRF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( POR0(1:3,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'POR0'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TORL(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'TORL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( VISL(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'VISL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RSDL(1:LUK,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'RSDL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( BTGL(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'BTGL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SDPF(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SDPF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SDPM(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SDPM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NPHAZ(1:LSU,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'NPHAZ'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_FDVP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_FLUX
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
!     Allocate global array memory for flux variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE FLUX
      USE GRID
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_FLUX'
      ALLOCATE( UL(1:LSFV,1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'UL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( VL(1:LSFV,1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'VL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( WL(1:LSFV,1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'WL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ULW(1:LSFV,1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ULW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( VLW(1:LSFV,1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'VLW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( WLW(1:LSFV,1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'WLW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_FLUX group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_GMBC
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
!     CO2 (Carbon Sequestration in Deep Saline Aquifers) Mode
!
!     Allocate global array memory for geomechanical boundary
!     condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 4 October 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
      USE GEO_MECH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_GMBC'
      ALLOCATE( BC_GM(1:LBCV_GM,1:LBTM_GM,1:LBCIN_GM),STAT=ISTAT )
      CHMSG = 'BC_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      LBCX = MAX( NBC_GM(ID+1),1 )
      ALLOCATE( IBCT_GM(1:3,1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCT_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCM_GM(1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCM_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCC_GM(1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCC_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCN_GM(1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCN_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCIN_GM(1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCIN_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCD_GM(1:LBCX),STAT=ISTAT )
      CHMSG = 'IBCD_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCR_GM(1:LBCIN_GM),STAT=ISTAT )
      CHMSG = 'IBCR_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_GMBC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_GMEC
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
!     CO2 (Carbon Sequestration in Deep Saline Aquifers) Mode
!
!     Allocate global array memory for geomechanics variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 4 October 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
      USE GEO_MECH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_GMEC'
      ALLOCATE( PROP_GM(1:7,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PROP_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IPROP_GM(1:3,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IPROP_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IHCM_GM(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IHCM_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( EPS_GM(1:6,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'EPS_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SIG_GM(1:6,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SIG_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( P_GM(1:3,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'P_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( EPSV_GM(1:3,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'EPSV_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( EPSV_CMP(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'EPSV_CMP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SIGV_GM(1:3,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SIGV_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SIGV_CMP(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SIGV_CMP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RSDM_GM(1:LEPD),STAT=ISTAT )
      CHMSG = 'RSDM_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( U_GM(1:2,1:NFNGN(ID+1)),STAT=ISTAT )
      CHMSG = 'U_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( V_GM(1:2,1:NFNGN(ID+1)),STAT=ISTAT )
      CHMSG = 'V_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( W_GM(1:2,1:NFNGN(ID+1)),STAT=ISTAT )
      CHMSG = 'W_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IM_GM(1:NFNGN(ID+1)),STAT=ISTAT )
      CHMSG = 'IM_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NE_GM(1:8,1:NFNGN(ID+1)),STAT=ISTAT )
      CHMSG = 'NE_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NK_GM(1:8,1:8),STAT=ISTAT )
      CHMSG = 'NK_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NF_GM(1:NFNGN(ID+1)),STAT=ISTAT )
      CHMSG = 'NF_GM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IGHN(1:NFNGN(ID+1)),STAT=ISTAT )
      CHMSG = 'IGHN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_GMEC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_GRID
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
!     Allocate global array memory for grid variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_GRID'
      ALLOCATE( AFX(1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'AFX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( AFY(1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'AFY'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( AFZ(1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'AFZ'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DXGP(1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'DXGP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DYGP(1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'DYGP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DZGP(1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'DZGP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( GRVX(1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'GRVX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( GRVY(1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'GRVY'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( GRVZ(1:2,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'GRVZ'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XE(1:8,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'XE'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YE(1:8,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'YE'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ZE(1:8,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ZE'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RP(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'RP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( XP(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'XP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YP(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'YP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ZP(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ZP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( VOL(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'VOL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DXGF(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'DXGF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DYGF(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'DYGF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DZGF(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'DZGF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ICM(1:6,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ICM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( INBS(1:6,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'INBS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IXP(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IXP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IGHC(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IGHC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ND(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ND'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ITAG(1:NP,1:NP),STAT=ISTAT )
      CHMSG = 'ITAG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_GRID group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_HYST
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
!     Allocate global array memory for hysteretic k-s-P function 
!     variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
      USE HYST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_HYST'
      ALLOCATE( ASL(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ASL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ASLMIN(1:LSU,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ASLMIN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ASTMIN(1:LSU,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ASTMIN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ASGT(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ASGT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SGT(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SGT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IPH(1:LSU,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IPH'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_HYST group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_OUTPU
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
!     Allocate global array memory for output control variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE OUTPU
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_OUTPU'
      ALLOCATE( PRTM(1:LPTM),STAT=ISTAT )
      CHMSG = 'PRTM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SF(1:2,1:LSF),STAT=ISTAT )
      CHMSG = 'SF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CNVPLOT(1:LVPLOT),STAT=ISTAT )
      CHMSG = 'CNVPLOT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CNVREF(1:LVREF),STAT=ISTAT )
      CHMSG = 'CNVREF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CNVSF(1:2,1:LSF),STAT=ISTAT )
      CHMSG = 'CNVSF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IPLOT(1:LVPLOT),STAT=ISTAT )
      CHMSG = 'IPLOT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IREF(1:LVREF),STAT=ISTAT )
      CHMSG = 'IREF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NDREF(1:LREF),STAT=ISTAT )
      CHMSG = 'NDREF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NDREFL(1:LREF),STAT=ISTAT )
      CHMSG = 'NDREFL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NDREFI(1:LREF),STAT=ISTAT )
      CHMSG = 'NDREFI'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISFT(1:LSF),STAT=ISTAT )
      CHMSG = 'ISFT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISFF(1:LSF),STAT=ISTAT )
      CHMSG = 'ISFF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISFD(1:LSF),STAT=ISTAT )
      CHMSG = 'ISFD'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISFGP(1:LSF),STAT=ISTAT )
      CHMSG = 'ISFGP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( UNPLOT(1:LVPLOT),STAT=ISTAT )
      CHMSG = 'UNPLOT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( UNREF(1:LVREF),STAT=ISTAT )
      CHMSG = 'UNREF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CHREF(1:LVREF),STAT=ISTAT )
      CHMSG = 'CHREF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISF(1:LSF),STAT=ISTAT )
      CHMSG = 'ISF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IOFFSET_SF(1:LSF),STAT=ISTAT )
      CHMSG = 'IOFFSET_SF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( UNSF(1:2,1:LSF),STAT=ISTAT )
      CHMSG = 'UNSF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CHSF(1:2,1:LSF),STAT=ISTAT )
      CHMSG = 'CHSF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( FNSF(1:LSF),STAT=ISTAT )
      CHMSG = 'FNSF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_OUTPU group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_PROP
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
!     Allocate global array memory for property variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE PROP
      USE GRID
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_PROP'
      ALLOCATE( CPS(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'CPS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RHOS(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'RHOS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PCMP(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PCMP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TCMP(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'TCMP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( HCMWE(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'HCMWE'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CMP(1:4,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'CMP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( POR(1:6,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'POR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TOR(1:6,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'TOR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ITOR(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ITOR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PERM(1:9,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'PERM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IPRF(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IPRF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SCHR(1:LSCHR,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SCHR'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISCHR(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'GRVZ'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISM(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ISM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RPLC(1:LRPLC,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'RPLC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CPLY_RL(1:LPOLYC,1:LPOLYN,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'CPLY_RL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CPLY_SL(1:LPOLYC,1:LPOLYN,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'CPLY_SL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISKP(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ISKP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IRPL(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IRPL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IZ(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IZ'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISLTBL(1:4,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ISLTBL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IRLTBL(1:4,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'IRLTBL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NPLY_RL(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'NPLY_RL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NPLY_SL(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'NPLY_SL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_PROP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_REACT
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
!     CO2 (Carbon Sequestration in Deep Saline Aquifers) Mode
!
!     Allocate global array memory for reactive species variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 17 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE REACT
      USE GRID
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_REACT'
      ALLOCATE( ACTV16(1:6),STAT=ISTAT )
      CHMSG = 'ACTV16'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( C_PH(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'C_PH'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CHARG(1:LSPR),STAT=ISTAT )
      CHMSG = 'CHARG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( FACTV(1:LSPR),STAT=ISTAT )
      CHMSG = 'FACTV'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ACTVS(1:LSPR),STAT=ISTAT )
      CHMSG = 'ACTVS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( EQ_C(1:LSEC,1:LEQC),STAT=ISTAT )
      CHMSG = 'EQ_C'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( EQ_E(1:LSEE,1:LEQE),STAT=ISTAT )
      CHMSG = 'EQ_E'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( EQ_K(1:LSEK+LREK,1:LEQK),STAT=ISTAT )
      CHMSG = 'EQ_K'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RC_E(1:5,1:LRCE),STAT=ISTAT )
      CHMSG = 'RC_E'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RC_K(1:LSPK+11,1:LCKN,1:LRCK),STAT=ISTAT )
      CHMSG = 'RC_K'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SP_C(1:NFCGC(ID+1),1:LSPR),STAT=ISTAT )
      CHMSG = 'SP_C'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SP_CO(1:NFCGC(ID+1),1:LSPR),STAT=ISTAT )
      CHMSG = 'SP_CO'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SP_CI(1:NFCGC(ID+1),1:LSPR),STAT=ISTAT )
      CHMSG = 'SP_CI'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SP_CMN(1:NFCGC(ID+1),1:LSPS),STAT=ISTAT )
      CHMSG = 'SP_CMN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SP_RATE(1:NFCGC(ID+1),1:LSPS),STAT=ISTAT )
      CHMSG = 'SP_RATE'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SP_AREA(1:NFCGC(ID+1),1:LSPS),STAT=ISTAT )
      CHMSG = 'SP_AREA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SP_CBO(1:LBCC,1:LSPR),STAT=ISTAT )
      CHMSG = 'SP_CBO'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IC_SP(1:NFCGC(ID+1),1:LSPR),STAT=ISTAT )
      CHMSG = 'IC_SP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISP_MN(1:LSPR),STAT=ISTAT )
      CHMSG = 'ISP_MN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SP_SDCL(1:3),STAT=ISTAT )
      CHMSG = 'SP_SDCL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SP_L(1:3,1:LSPL),STAT=ISTAT )
      CHMSG = 'SP_L'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SP_S(1:2,1:LSPS),STAT=ISTAT )
      CHMSG = 'SP_S'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SPNME(1:LSPE),STAT=ISTAT )
      CHMSG = 'SPNME'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISP_E(1:LSPE),STAT=ISTAT )
      CHMSG = 'ISP_E'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RS_S(1:3,1:LSPS,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'RS_S'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YSPG(1:NFCGC(ID+1),1:(LEQC+LEQK)),STAT=ISTAT )
      CHMSG = 'YSPG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YSPL(1:NFCGC(ID+1),1:(LEQC+LEQK)),STAT=ISTAT )
      CHMSG = 'YSPL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SPNMC(1:LEQC),STAT=ISTAT )
      CHMSG = 'SPNMC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SPNMK(1:LEQK),STAT=ISTAT )
      CHMSG = 'SPNMK'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SPNMG(1:LSPG),STAT=ISTAT )
      CHMSG = 'SPNMG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SPNML(1:LSPL),STAT=ISTAT )
      CHMSG = 'SPNML'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SPNMS(1:LSPS),STAT=ISTAT )
      CHMSG = 'SPNMS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RCNME(1:LRCE),STAT=ISTAT )
      CHMSG = 'RCNME'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RCNMK(1:LRCK),STAT=ISTAT )
      CHMSG = 'RCNMK'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IEQ_C(1:LSEC+1,1:LEQC),STAT=ISTAT )
      CHMSG = 'IEQ_C'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IEQ_E(1:LSEE+2,1:LEQE),STAT=ISTAT )
      CHMSG = 'IEQ_E'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IEQ_K(1:LSEK+LREK+2,1:LEQK),STAT=ISTAT )
      CHMSG = 'IEQ_K'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IEQ_S(1:LSPR),STAT=ISTAT )
      CHMSG = 'IEQ_S'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IEL_LK(1:LSPE),STAT=ISTAT )
      CHMSG = 'IEL_LK'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISP_S(1:LEQE+LEQC+LEQK),STAT=ISTAT )
      CHMSG = 'ISP_S'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IRC_K(1:LSPK+3,1:LRCK),STAT=ISTAT )
      CHMSG = 'IRC_K'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IRCKT(1:LRCK),STAT=ISTAT )
      CHMSG = 'IRCKT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IRCKN(LSPK+11),STAT=ISTAT )
      CHMSG = 'IRCKN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISPLK(1:2*LNGC+LSPR+14),STAT=ISTAT )
      CHMSG = 'ISPLK'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISP_OW(1:LSPS,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'ISP_OW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IMMB(1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'IMMB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CFMX(1:LMC),STAT=ISTAT )
      CHMSG = 'CFMX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_REACT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_SOLTN
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
!     Allocate global array memory for solution control variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_SOLTN'
      ALLOCATE( TMPS(1:LEPD),STAT=ISTAT )
      CHMSG = 'TMPS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TMPE(1:LEPD),STAT=ISTAT )
      CHMSG = 'TMPE'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TMPD(1:LEPD),STAT=ISTAT )
      CHMSG = 'TMPD'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TMPX(1:LEPD),STAT=ISTAT )
      CHMSG = 'TMPX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TMPN(1:LEPD),STAT=ISTAT )
      CHMSG = 'TMPN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TMPA(1:LEPD),STAT=ISTAT )
      CHMSG = 'TMPA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TMPC(1:LEPD),STAT=ISTAT )
      CHMSG = 'TMPC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RSDM(1:LEPD),STAT=ISTAT )
      CHMSG = 'RSDM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( RSD(1:(LUK*(1+LWELL+LSPILL))),STAT=ISTAT )
      CHMSG = 'RSD'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( WFMN(1:20),STAT=ISTAT )
      CHMSG = 'WFMN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IEQGC(1:LNGC),STAT=ISTAT )
      CHMSG = 'IEQGC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ITOFF(1:LEPD),STAT=ISTAT )
      CHMSG = 'ITOFF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NRIM(1:LEPD),STAT=ISTAT )
      CHMSG = 'NRIM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NSD(1:(LUK*(1+LWELL+LSPILL))),STAT=ISTAT )
      CHMSG = 'NSD'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( HCNM(1:LNHC),STAT=ISTAT )
      CHMSG = 'HCNM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( MNOD(1:LSFV),STAT=ISTAT )
      CHMSG = 'MNOD'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( MADJ(1:LSFV),STAT=ISTAT )
      CHMSG = 'MADJ'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( MFLX(1:LSFV),STAT=ISTAT )
      CHMSG = 'MFLX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( MPOS(1:LSFV),STAT=ISTAT )
      CHMSG = 'MPOS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( MNEG(1:LSFV),STAT=ISTAT )
      CHMSG = 'MNEG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( MPOSB(1:LSV),STAT=ISTAT )
      CHMSG = 'MPOSB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( MNEGB(1:LSV),STAT=ISTAT )
      CHMSG = 'MNEGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_SOLTN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_SOURC
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
!     Allocate global array memory for source variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE SOURC
      USE GRID
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_SOURC'
      LX = 8+LSOLU+LSPT+LNGC
      LSRX = MAX( NSR(ID+1),1 )
      ALLOCATE( SRC(1:LX,1:LSTM,1:LSR),STAT=ISTAT )
      CHMSG = 'SRC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SRCW(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SRCW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SRCA(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SRCA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SRCS(1:LSV,1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SRCS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SRCIW(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SRCIW'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SRCIA(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SRCIA'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SRCIS(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'SRCIS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISRN(1:LSRX),STAT=ISTAT )
      CHMSG = 'ISRN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISRDM(1:13,1:LSRX),STAT=ISTAT )
      CHMSG = 'ISRDM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISRM(1:LSRX),STAT=ISTAT )
      CHMSG = 'ISRM'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISRIN(1:LSRX),STAT=ISTAT )
      CHMSG = 'ISRIN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ISRT(1:LSRX),STAT=ISTAT )
      CHMSG = 'ISRT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_SOURC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_TABL
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
!     Allocate global array memory for tabular data variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
      USE TABL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_TABL'
      ALLOCATE( TBLX(1:LTBL),STAT=ISTAT )
      CHMSG = 'TBLX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TBLY(1:LTBL),STAT=ISTAT )
      CHMSG = 'TBLY'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TBLDDX(1:LTBL),STAT=ISTAT )
      CHMSG = 'TBLDDX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( TBLDDY(1:LTBL),STAT=ISTAT )
      CHMSG = 'TBLDDY'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_TABL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ALLOC_TRNSPT
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
!     Allocate global array memory for solute transport variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
      USE TRNSPT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ALLOC_TRNSPT'
      ALLOCATE( C(1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'C'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CO(1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'CO'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CNL(1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'CNL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CB(1:LBCC,1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'CB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CBO(1:LBCC,1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'CBO'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YL(1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'YL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YG(1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'YG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YN(1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'YN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YLB(1:LBCC,1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'YLB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( YGB(1:LBCC,1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'YGB'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( PCSL(1:5,1:NFCGC(ID+1),1:LSOLU),STAT=ISTAT )
      CHMSG = 'PCSL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CRNTG(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'CRNTG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CRNTL(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'CRNTL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( UC(1:2,1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'UC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( VC(1:2,1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'VC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( WC(1:2,1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'WC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DISPL(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'DISPL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( DISPT(1:NFCGC(ID+1)),STAT=ISTAT )
      CHMSG = 'DISPT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SDCL(1:3,1:NFCGC(ID+1),1:LSOLU),STAT=ISTAT )
      CHMSG = 'SDCL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SMDG(1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'SMDG'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SMDL(1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'SMDL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SOLUT(1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'SOLUT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CCL_CRN(1:LSOLU),STAT=ISTAT )
      CHMSG = 'CCL_CRN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( HLF(1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'HLF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( SRCIC(1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'SRCIC'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( CHDF(1:LSOLU,1:LSOLU),STAT=ISTAT )
      CHMSG = 'CHDF'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IEDL(1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'IEDL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( ICT(1:NFCGC(ID+1),1:LSOLU+LSPT),STAT=ISTAT )
      CHMSG = 'ICT'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IPCL(1:LSOLU),STAT=ISTAT )
      CHMSG = 'IPCL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IPCSL(1:NFCGC(ID+1),1:LSOLU),STAT=ISTAT )
      CHMSG = 'IPCSL'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( N_CRN(1:LSOLU+1),STAT=ISTAT )
      CHMSG = 'N_CRN'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCDS(1:LCDC+LSOLU),STAT=ISTAT )
      CHMSG = 'IBCDS'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( NBCDP(1:LCDC),STAT=ISTAT )
      CHMSG = 'NBCDP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
      ALLOCATE( IBCDP(1:LCDS,1:LCDP,1:LCDC),STAT=ISTAT )
      CHMSG = 'IBCDP'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ALLOC_TRNSPT group  ---
!
      RETURN
      END




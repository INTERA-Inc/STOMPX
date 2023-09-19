!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_CO2
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
!     Initialize global variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 September 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE PROP
      USE OUTPU
      USE GRID
      USE GLB_PAR
      USE FILES
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER, DIMENSION(:), ALLOCATABLE :: NFCX
      INTEGER	STATUS(MPI_STATUS_SIZE)	
      INTEGER	(KIND=MPI_OFFSET_KIND) OFFSET
      INTEGER, DIMENSION(NP) :: ICOLOR
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_CO2'
!
!---  Allocate memory for flux variables  ---
!
      CALL ALLOC_FLUX
!      PRINT *,'Post ALLOC_FLUX: ID = ',ID
!
!---  Initialize memory for flux variables  ---
!
      CALL INTLZ_FLUX
!      PRINT *,'Post INTLZ_FLUX: ID = ',ID
!
!---  Create an output.bin file  ---
!
      CALL MPI_FILE_OPEN( MPI_COMM_WORLD,'output.bin',
     &  MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,IWR,IERR )
!
!---  Identify local grid cells for reference nodes  ---
!
      NREFL = 0
      DO M = 1,NREF
        DO N = 1,NFCGC(ID+1)
          IF( IXP(N).EQ.0 .OR. IGHC(N).EQ.1 ) CYCLE
          IF( NDREF(M).EQ.ND(N) ) THEN
            NREFL = NREFL + 1
            NDREFL(NREFL) = N
            NDREFI(NREFL) = M
            EXIT
          ENDIF
        ENDDO
      ENDDO
      IOFFSET_REF = 0
      DO L = 1,LSF
        IOFFSET_SF(L) = 0
      ENDDO
!
!---  Allocate memory for NFCX, temporary number of field cells
!     on each processor  ---
!
      ALLOCATE( NFCX(1:NP),STAT=ISTAT )
      CHMSG = 'NFCX'
      CALL ALLOC_ERROR( CHMSG,ISTAT )
!
!---  Set NFC array and NFC_G  ---
!
      DO N = 1,NP
        NFCX(N) = 0
      ENDDO
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for ghost cells only  ---
!
        IF( IGHC(N).EQ.1 ) CYCLE
        NFCX(ID+1) = NFCX(ID+1) + 1
      ENDDO
      CALL MPI_ALLREDUCE( NFCX,NFC,NP,MPI_INTEGER,MPI_MAX,
     &  MPI_COMM_WORLD,IERR )
      NFC_G = 0
      NFC_L = 0
      DO I = 1,ID
        NFC_L = NFC_L + NFC(I)
      ENDDO
      DO I = 1,NP
        NFC_G = NFC_G + NFC(I)
      ENDDO
!
!---  Deallocate memory for NFCX  ---
!
      DEALLOCATE( NFCX,STAT=ISTAT )
      CHMSG = 'NFCX'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Color processors in vertical groups  ---
!
      DO KP = 1,KPFLD
        DO JP = 1,JPFLD
          DO IP = 1,IPFLD
            NP = (KP-1)*JPFLD*IPFLD + (JP-1)*IPFLD + IP
            ICOLOR(NP) = (JP-1)*IPFLD + IP
!            PRINT *,'ICOLOR(',NP,') = ',ICOLOR(NP),'ID = ',ID
          ENDDO
        ENDDO
      ENDDO
!
!---  Split the global communicator based on color and use the 
!     original rank for ordering  ---
!
      CALL MPI_COMM_SPLIT( MPI_COMM_WORLD,ICOLOR(ID+1),ID,
     &  MPI_COMM_VERT,IERR )
      CALL MPI_COMM_RANK( MPI_COMM_VERT,IDV,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_VERT,NPV,IERR )
!      PRINT *,'ID = ',ID,'NP = ',NP,'IDV = ',IDV,'NPV = ',NPV
!
!---  Initialize variable properties  ---
!
      DO N = 1,NFCGC(ID+1)
        IF( IXP(N).EQ.0 ) CYCLE
        POR0(1,N) = POR(1,N)
        POR0(2,N) = POR(2,N)
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_BCV
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
!     Initialize global array memory for boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 January 2022
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_BCV'
      DO M = 1,LBCIN
        DO L = 1,LBTM
          DO K = 1,LBCV
            BC(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      LX = LUK+LSOLU*LC
      LBCX = MAX( NBC(ID+1),1 )
      DO L = 1,LBCX
        XPBC(L) = 0.D+0
        YPBC(L) = 0.D+0
        ZPBC(L) = 0.D+0
        IBCN(L) = 0
        IBCD(L) = 0
        IBCM(L) = 0
        IBCIN(L) = 0
        IBCC(L) = 0
        DO K = 1,LX
          IBCT(K,L) = 0
        ENDDO
        DO K = 1,LSPBC+1
          IBCSP(K,L) = 0
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_BCV group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_BCVP
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
!     Initialize global array memory for boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 September 2021
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_BCVP'
      LBCX = MAX( NBC(ID+1),1 )
      DO L = 1,LBCX
        DO K = 1,LSV
          TB(K,L) = 0.D+0
          PLB(K,L) = 0.D+0
          PGB(K,L) = 0.D+0
          RHOGB(K,L) = 0.D+0
          RHOLB(K,L) = 0.D+0
          BTGLB(K,L) = 0.D+0
          XGAB(K,L) = 0.D+0
          XGWB(K,L) = 0.D+0
          XMGAB(K,L) = 0.D+0
          XMGWB(K,L) = 0.D+0
          XLAB(K,L) = 0.D+0
          XLWB(K,L) = 0.D+0
          XMLAB(K,L) = 0.D+0
          XMLWB(K,L) = 0.D+0
          RKGB(K,L) = 0.D+0
          RKLB(K,L) = 0.D+0
          PVAB(K,L) = 0.D+0
          PVWB(K,L) = 0.D+0
          DFLAB(K,L) = 0.D+0
          DFLSB(K,L) = 0.D+0
          SGB(K,L) = 0.D+0
          SLB(K,L) = 0.D+0
          YLSB(K,L) = 0.D+0
          XLSB(K,L) = 0.D+0
          XMLSB(K,L) = 0.D+0
          RHOMLB(K,L) = 0.D+0
          RHOMGB(K,L) = 0.D+0
          DFGAB(K,L) = 0.D+0
          DFGWB(K,L) = 0.D+0
          PORDB(K,L) = 0.D+0
          PORTB(K,L) = 0.D+0
          TORGB(K,L) = 0.D+0
          TORLB(K,L) = 0.D+0
          HGAB(K,L) = 0.D+0
          HGB(K,L) = 0.D+0
          UEGB(K,L) = 0.D+0
          UELB(K,L) = 0.D+0
          HGWB(K,L) = 0.D+0
          HLB(K,L) = 0.D+0
          HLWB(K,L) = 0.D+0
          VISGB(K,L) = 0.D+0
          VISLB(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_BCVP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_COUP_WELL
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
!     Initialize coupled-well variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 December 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE SOURC
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_COUP_WELL'
      DO L = 1,LN_CW
        DNR_CW(L) = 0.D+0
        DP_CW(L) = 0.D+0
        FX_CW(L) = 0.D+0
        PL_CW(L) = 0.D+0
        RHOF_CW(L) = 0.D+0
        TML_CW(L) = 0.D+0
        WNM_CW(L) = '  '
        ICC_CW(L) = 0
        IM_CW(L) = 0
        IT_CW(L) = 0
        JM_CW(L) = 0
        NSP_CW(L) = 0
        DO K = 1,LSOLU_CW
          ISOLU_CW(K,L) = 0
        ENDDO
        DO K = 1,LSPC_CW
          ISPC_CW(K,L) = 0
        ENDDO
        DO K = 1,LEQC
          ISOLC_CW(K,L) = 0
        ENDDO
        DO K = 1,LEQK
          ISOLK_CW(K,L) = 0
        ENDDO
        DO K = 1,LWTP_CW
          IMP_CW(K,L) = 0
          ITS_CW(K,L) = 0
        ENDDO
        DO K = 1,10
          ID_CW(K,L) = 0
        ENDDO
        DO K = 1,2
          T_CW(K,L) = 0.D+0
        ENDDO
        DO K = 1,3
          P_CW(K,L) = -1.D+20
          FF_CW(K,L) = 0.D+0
        ENDDO
        DO K = 1,8
          QM_CW(K,L) = 0.D+0
        ENDDO
        DO K = 1,LUK_CW+1
          RS_CW(K,L) = 0.D+0
        ENDDO
        DO K = 1,LWT_CW
          DO J = 1,7+LNGC
            VAR_CW(J,K,L) = 0.D+0
          ENDDO
          DO J = 1,LSOLU_CW
            VARC_CW(J,K,L) = 0.D+0
          ENDDO
          DO J = 1,LSPC_CW
            VARSP_CW(J,K,L) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,LWN_CW
        PLX_CW(L) = 0.D+0
        PLY_CW(L) = 0.D+0
        PLZ_CW(L) = 0.D+0
        TF_CW(L) = 0.D+0
        INV_CW(L) = 0
        IWNG_CW(L) = 0
        IWN_CW(L) = 0
        IWP_CW(L) = 0
        DO K = 1,2
          XP_CW(K,L) = 0.D+0
          YP_CW(K,L) = 0.D+0
          ZP_CW(K,L) = 0.D+0
        ENDDO
        DO K = 1,(LUK+2)
          FXA_CW(K,L) = 0.D+0
          FXS_CW(K,L) = 0.D+0
          FXW_CW(K,L) = 0.D+0
        ENDDO
        DO K = 1,4
          Q_CW(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LWF_CW
        PF_CW(L) = 0.D+0
        IWFG_CW(L) = 0
        IWF_CW(L) = 0
        IXPG_CW(L) = 0
      ENDDO
      DO L = 1,LWI_CW
        IS_CW(L) = 0
        DO K = 1,5
          PAR_CW(K,L) = 0.D+0
        ENDDO
        DO K = 1,2
          XTP_CW(K,L) = 0.D+0
          YTP_CW(K,L) = 0.D+0
          ZTP_CW(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LVREF
        IREF_CW(L) = 0
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_FDVP
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
!     Initialize global array memory for hydrate field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 September 2021
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_FDVP'
      DO L = 1,NFCGC(ID+1)
        DO K = 1,LSV
          BTGL(K,L) = 0.D+0
          DNR(K,L) = 0.D+0
          RHOG(K,L) = 0.D+0
          RHOL(K,L) = 0.D+0
          XGA(K,L) = 0.D+0
          XGW(K,L) = 0.D+0
          XMGA(K,L) = 0.D+0
          XMGW(K,L) = 0.D+0
          XLA(K,L) = 0.D+0
          XLW(K,L) = 0.D+0
          XLO(K,L) = 0.D+0
          XMLA(K,L) = 0.D+0
          XMLW(K,L) = 0.D+0
          RKG(K,L) = 0.D+0
          RKL(K,L) = 0.D+0
          DFLA(K,L) = 0.D+0
          DFLS(K,L) = 0.D+0
          XLS(K,L) = 0.D+0
          XMLS(K,L) = 0.D+0
          RHOSP(K,L) = 0.D+0
          HSP(K,L) = 0.D+0
          SS(K,L) = 0.D+0
          TMS(K,L) = 0.D+0
          RHOML(K,L) = 0.D+0
          RHOMG(K,L) = 0.D+0
          DFGA(K,L) = 0.D+0
          DFGW(K,L) = 0.D+0
          PORD(K,L) = 0.D+0
          PORT(K,L) = 0.D+0
          PERMRF(K,L) = 1.D+0
          TORG(K,L) = 0.D+0
          TORL(K,L) = 0.D+0
          HGA(K,L) = 0.D+0
          HG(K,L) = 0.D+0
          HL(K,L) = 0.D+0
          HGW(K,L) = 0.D+0
          HLW(K,L) = 0.D+0
          UEG(K,L) = 0.D+0
          UEL(K,L) = 0.D+0
          VISG(K,L) = 0.D+0
          VISL(K,L) = 0.D+0
          PN(K,L) = 0.D+0
        ENDDO
        SDPF(L) = 0.D+0
        SDPM(L) = 0.D+0
        DO K = 1,3
          POR0(K,L) = 0.D+0
        ENDDO
        DO K = 1,LUK
          RSDL(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_FDVP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_FLUX
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
!     Initialize global array memory for flux variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 September 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE GLB_PAR
      USE FLUX
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_FLUX'
      DO L = 1,NFCGC(ID+1)
        DO K = 1,2
          DO J = 1,LSFV
            UDGA(J,K,L) = 0.D+0
            VDGA(J,K,L) = 0.D+0
            WDGA(J,K,L) = 0.D+0
            UGA(J,K,L) = 0.D+0
            VGA(J,K,L) = 0.D+0
            WGA(J,K,L) = 0.D+0
            UDGW(J,K,L) = 0.D+0
            VDGW(J,K,L) = 0.D+0
            WDGW(J,K,L) = 0.D+0
            UGW(J,K,L) = 0.D+0
            VGW(J,K,L) = 0.D+0
            WGW(J,K,L) = 0.D+0
            UDLA(J,K,L) = 0.D+0
            VDLA(J,K,L) = 0.D+0
            WDLA(J,K,L) = 0.D+0
            ULA(J,K,L) = 0.D+0
            VLA(J,K,L) = 0.D+0
            WLA(J,K,L) = 0.D+0
            UDLW(J,K,L) = 0.D+0
            VDLW(J,K,L) = 0.D+0
            WDLW(J,K,L) = 0.D+0
            ULW(J,K,L) = 0.D+0
            VLW(J,K,L) = 0.D+0
            WLW(J,K,L) = 0.D+0
            UDS(J,K,L) = 0.D+0
            VDS(J,K,L) = 0.D+0
            WDS(J,K,L) = 0.D+0
            US(J,K,L) = 0.D+0
            VS(J,K,L) = 0.D+0
            WS(J,K,L) = 0.D+0
            UL(J,K,L) = 0.D+0
            VL(J,K,L) = 0.D+0
            WL(J,K,L) = 0.D+0
            UG(J,K,L) = 0.D+0
            VG(J,K,L) = 0.D+0
            WG(J,K,L) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_FLUX group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_GMBC
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
!     Initialize memory for geomechanical boundary condition variables
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_GMBC'
      DO M = 1,LBCIN_GM
        DO L = 1,LBTM_GM
          DO K = 1,LBCV_GM
            BC_GM(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,NBC_GM(ID+1)
        DO K = 1,3
          IBCT_GM(K,L) = 0
        ENDDO
        IBCC_GM(L) = 0
        IBCN_GM(L) = 0
        IBCIN_GM(L) = 0
        IBCD_GM(L) = 0
        IBCM_GM(L) = 0
      ENDDO
      DO L = 1,LBCIN_GM
        IBCR_GM(L) = 0
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_GMBC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_GMEC
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
!     Initialize memory for geomechanical variables
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_GMEC'
      DO L = 1,NFCGC(ID+1)
        DO K = 1,6
          SIG_GM(1:6,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,NFCGC(ID+1)
        DO K = 1,9
          PROP_GM(K,L) = 0.D+0
        ENDDO
        DO K = 1,3
          IPROP_GM(K,L) = 0
          P_GM(K,L) = 0.D+0
          EPSV_GM(K,L) = 0.D+0
          SIGV_GM(K,L) = 0.D+0
        ENDDO
        IHCM_GM(L) = 0
        DO K = 1,6
          EPS_GM(1:6,L) = 0.D+0
        ENDDO
        EPSV_CMP(L) = 0.D+0
      ENDDO
      DO L = 1,LEPD
        RSDM_GM(L) = 0.D+0
      ENDDO
      DO L = 1,NFNGN(ID+1)
        DO K = 1,2
          U_GM(K,L) = 0.D+0
          V_GM(K,L) = 0.D+0
          W_GM(K,L) = 0.D+0
        ENDDO
        IM_GM(L) = 0
      ENDDO
      GRAV_GM = 9.81D+0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_GMEC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_HYST
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
!     Initialize global array memory for hysteretic k-s-P function 
!     variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 September 2021
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_HYST'
      DO L = 1,NFCGC(ID+1)
        DO K = 1,LSU
          ASLMIN(K,L) = 1.D+0
        ENDDO
        DO K = 1,LSV
          SGT(K,L) = 0.D+0
        ENDDO
        ASL(L) = 0.D+0
        ASGT(L) = 0.D+0
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_HYST group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_REACT
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
!     Initialize solute transport variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 6 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE SOURC
      USE REACT
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_REACT'
      ACTVC = 1.D+0
      DT_RST = 0.D+0
      DTI_RST = 0.D+0
      IACTEX = 0
      IACTV = 0
      ISP_IEDL = 0
      N_RST = 10
      ECKE_ER = .FALSE.
      NESITE = 0
      NRCE = 0
      NRCK = 0
      NRTSI = 0
      NSPC = 0
      NSPE = 0
      NSPG = 0
      NSPK = 0
      NSPL = 0
      NSPN = 0
      NSPLK = 0
      NEQC = 0
      NEQE = 0
      NEQK = 0
      NSPR = 0
      NSPS = 0
      CMIN = 0.D+0
      SP_MDG = 0.D+0
      SP_MDL = 0.D+0
      SP_MDN = 0.D+0
      TM_RST = 0.D+0
      DO K = 1,6
        ACTV16(K) = 0.D+0
      ENDDO
      DO K = 1,NFCGC(ID+1)
        C_PH(K) = 0.D+0
      ENDDO
      DO K = 1,LSPR
        CHARG(K) = 0.D+0
      ENDDO
      DO K = 1,LSPR
        FACTV(K) = 0.D+0
      ENDDO
      DO K = 1,LSPR
        ACTVS(K) = 0.D+0
      ENDDO
      DO L = 1,LEQC
        DO K = 1,LSEC
          EQ_C(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LEQE
        DO K = 1,LSEE
          EQ_E(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LEQK
        DO K = 1,LSEK+LREK
          EQ_K(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LRCE
        DO K = 1,5
          RC_E(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO M = 1,LRCK
        DO L = 1,LCKN
          DO K = 1,LSPK+11
            RC_K(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,LSPR
        DO K = 1,NFCGC(ID+1)
          SP_C(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSPR
        DO K = 1,NFCGC(ID+1)
          SP_CO(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSPR
        DO K = 1,NFCGC(ID+1)
          SP_CI(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSPS
        DO K = 1,NFCGC(ID+1)
          SP_CMN(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSPS
        DO K = 1,NFCGC(ID+1)
          SP_RATE(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSPS
        DO K = 1,NFCGC(ID+1)
          SP_AREA(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSPR
        DO K = 1,LBCC
          SP_CBO(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSPR
        DO K = 1,NFCGC(ID+1)
          IC_SP(K,L) = 0
        ENDDO
      ENDDO
      DO K = 1,LSPR
        ISP_MN(K) = 0
      ENDDO
      DO K = 1,3
        SP_SDCL(K) = 0.D+0
      ENDDO
      DO L = 1,LSPL
        DO K = 1,3
          SP_L(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSPS
        DO K = 1,2
          SP_S(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO K = 1,LSPE
        ISP_E(K) = 0
      ENDDO
      DO M = 1,NFCGC(ID+1)
        DO L = 1,LSPS
          DO K = 1,3
            RS_S(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,(LEQC+LEQK)
        DO K = 1,NFCGC(ID+1)
          YSPG(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,(LEQC+LEQK)
        DO K = 1,NFCGC(ID+1)
          YSPL(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO K = 1,LEQC
        SPNMC(K) = '  '
      ENDDO
      DO K = 1,LEQK
        SPNMK(K) = '  '
      ENDDO
      DO K = 1,LSPG
        SPNMG(K) = '  '
      ENDDO
      DO K = 1,LSPL
        SPNML(K) = '  '
      ENDDO
      DO K = 1,LSPS
        SPNMS(K) = '  '
      ENDDO
      DO K = 1,LSPE
        SPNME(K) = '  '
      ENDDO
      DO K = 1,LRCE
        RCNME(K) = '  '
      ENDDO
      DO K = 1,LRCK
        RCNMK(K) = '  '
      ENDDO
      DO L = 1,LEQC
        DO K = 1,LSEC+1
          IEQ_C(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LEQE
        DO K = 1,LSEE+2
          IEQ_E(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LEQK
        DO K = 1,LSEK+LREK+2
          IEQ_K(K,L) = 0
        ENDDO
      ENDDO
      DO K = 1,LSPR
        IEQ_S(K) = 0
      ENDDO
      DO K = 1,LSPE
        IEL_LK(K) = 0
      ENDDO
      DO K = 1,LEQE+LEQC+LEQK
        ISP_S(K) = 0
      ENDDO
      DO L = 1,LRCK
        DO K = 1,LSPK+3
          IRC_K(K,L) = 0
        ENDDO
      ENDDO
      DO K = 1,LRCK
        IRCKT(K) = 0
      ENDDO
      DO K = 1,LSPK+11
        IRCKN(K) = 0
      ENDDO
      DO K = 1,2*LNGC+LSPR+14
        ISPLK(K) = 0
      ENDDO
      DO L = 1,NFCGC(ID+1)
        DO K = 1,LSPS
          ISP_OW(K,L) = 0
        ENDDO
      ENDDO
      DO K = 1,LSOLU+LSPT
        IMMB(K) = 0
      ENDDO
      DO K = 1,LMC
        CFMX(K) = 0.D+0
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_REACT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_SOLTN
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
!     Initialize global array memory for hydrate field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 September 2021
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
!
!----------------------Executable Lines--------------------------------!
!
!      MNOD /2,3,2,4,2,5,2,6,2,7,2,8,2,9,2/
!      MADJ /2,2,3,2,4,2,5,2,6,2,7,2,8,2,9/
!      MPOS /2,3,2,4,2,5,2,6,2,7,2,8,2,9,2/
!      MNEG /2,2,3,2,4,2,5,2,6,2,7,2,8,2,9/
!      MFLX /1,3,2,5,4,7,6,9,8,11,10,13,12,15,14/
!      MPOSB /1,2,4,6,8,10,12,14/
!      MNEGB /1,3,5,7,9,11,13,15/
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_SOLTN'
      DO L = 1,LSFV
        IF( MOD(L,2).EQ.1 ) THEN
          MNOD(L) = 2
        ELSE
          MNOD(L) = L/2 + 2
        ENDIF
        IF( MOD(L,2).EQ.1 ) THEN
          MADJ(L) = (L-1)/2 + 2
        ELSE
          MADJ(L) = 2
        ENDIF
        IF( MOD(L,2).EQ.1 ) THEN
          MPOS(L) = 2
        ELSE
          MPOS(L) = L/2 + 2
        ENDIF
        IF( MOD(L,2).EQ.1 ) THEN
          MNEG(L) = (L-1)/2 + 2
        ELSE
          MNEG(L) = 2
        ENDIF
        IF( MOD(L,2).EQ.1 ) THEN
          MFLX(L) = MAX( 1,(L-1) )
        ELSE
          MFLX(L) = L + 1
        ENDIF
      ENDDO
      DO L = 1,LSV
        MPOSB(L) = MAX( 1,((2*L)-2) )
        MNEGB(L) = (2*L) - 1
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_SOLTN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_SOURC
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
!     Initialize source variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 September 2021
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_SOURC'
      DO L = 1,NFCGC(ID+1)
        DO K = 1,LSV
          SRCW(K,L) = 0.D+0
          SRCA(K,L) = 0.D+0
          SRCS(K,L) = 0.D+0
        ENDDO
        SRCIW(L) = 0.D+0
        SRCIA(L) = 0.D+0
        SRCIS(L) = 0.D+0
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_SOURC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ_TRNSPT
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
!     Initialize solute transport variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 6 January 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE TRNSPT
      USE SOLTN
      USE SOURC
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ_TRNSPT'
      DO L = 1,NFCGC(ID+1)
        CRNTG(L) = 0.D+0
        CRNTL(L) = 0.D+0
        DISPL(L) = 0.D+0
        DISPT(L) = 0.D+0
        DO M = 1,LSOLU+LSPT
          C(L,M) = 0.D+0
          CO(L,M) = 0.D+0
          CNL(L,M) = 0.D+0
          YL(L,M) = 0.D+0
          YG(L,M) = 0.D+0
          SRCIC(L,M) = 0.D+0
          ICT(L,M) = 0
          DO K = 1,2
            UC(K,L,M) = 0.D+0
            VC(K,L,M) = 0.D+0
            WC(K,L,M) = 0.D+0
          ENDDO
        ENDDO
        DO M = 1,LSOLU
          DO K = 1,5
            PCSL(K,L,M) = 0.D+0
          ENDDO
          DO K = 1,3
            SDCL(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,LSOLU
        IPCGL(L) = 0
        IPCL(L) = 0
        DO K = 1,5
          PCGL(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LBCC
        DO M = 1,LSOLU+LSPT
          CB(L,M) = 0.D+0
          CBO(L,M) = 0.D+0
          YLB(L,M) = 0.D+0
          YGB(L,M) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSOLU+LSPT
        SMDG(L) = 0.D+0
        SMDL(L) = 0.D+0
        SOLUT(L) = '  '
        HLF(L) = 0.D+0
        IEDL(L) = 0
      ENDDO
      DO L = 1,LSOLU
        CCL_CRN(L) = 0.D+0
        DO M = 1,LSOLU
          CHDF(L,M) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSOLU+1
        N_CRN(L) = 0
      ENDDO
      DO L = 1,LCDC+LSOLU
        IBCDS(L) = 0
      ENDDO
      DO M = 1,LCDC
        NBCDP(M) = 0
        DO L = 1,LCDP
          DO K = 1,LCDS
            IBCDP(K,L,M) = 0
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ_TRNSPT group  ---
!
      RETURN
      END



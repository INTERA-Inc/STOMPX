!---------------------Fortran 90 Module--------------------------------!
!
      MODULE BCV
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
!     Boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: BC
      REAL*8, DIMENSION(:), ALLOCATABLE :: XPBC,YPBC,ZPBC
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCC,IBCD,IBCIN,IBCM,IBCN
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IBCT,IBCSP
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE BCVP
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
!     Mechanical, hydraulic, saturation function, aqueous relative
!     permeability, gas relative permeability, and nonaqueous liquid
!     relative permeability variables for boundary surfaces
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: TB,PLB,PGB,PSOB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RHOGB,RHOLB,BTGLB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XGAB,XGWB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XMGAB,XMGWB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XLAB,XLWB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XMLAB,XMLWB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RKGB,RKLB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PVAB,PVWB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: DFLAB,DFLSB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SGB,SLB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: YLSB,XLSB,XMLSB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RHOMLB,RHOMGB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: DFGAB,DFGWB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PORDB,PORTB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: TORGB,TORLB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HGAB,HGB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: UEGB,UELB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HGWB,HLB,HLWB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: VISGB,VISLB
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE CONST
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
!     Constant variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 :: SMALL = 1.D-20
      REAL*8 :: BIG = 1.D+20
      REAL*8 :: ZERO = 0.D+0
      REAL*8 :: ONE = 1.D+0
      REAL*8 :: EPSL = 1.D-14
      REAL*8 :: GRAV = 9.81D+0
      REAL*8 :: TABS = 273.15D+0
      REAL*8 :: TMX = 374.14D+0
      REAL*8 :: TMN = 0.01D+0
      REAL*8 :: PATM = 101325.D+0
      REAL*8 :: PMX = 2.212D+8
      REAL*8 :: PMN = 6.1125D+2
      REAL*8 :: RHORL = 998.32142721500441D+0
      REAL*8 :: RHORG = 1.199D+0
      REAL*8 :: VCRA = 86.2269D+0
      REAL*8 :: VISRL = 1.0176489259594200D-3
      REAL*8 :: VISRG = 1.82D-5
      REAL*8 :: TSPRF = 293.15D+0
      REAL*8 :: HDOD = 1.D+5
      REAL*8 :: WTMA = 44.010D+0
      REAL*8 :: WTMN = 28.01348D+0
      REAL*8 :: WTMW = 18.015D+0
      REAL*8 :: WTMS = 58.4428D+0
      REAL*8 :: RCU = 8.31434D+3
      REAL*8 :: RCW = 461.52D+0
      REAL*8 :: RCA = 188.918D+0
      REAL*8 :: TCRA = 304.1282D+0
      REAL*8 :: TCRW = 647.096D+0
      REAL*8 :: PCRW = 22.064D+6
      REAL*8 :: PCRA = 73.773D+5
      REAL*8 :: PTPA = 517.95D+3
      REAL*8 :: TTPA = 216.592D+0
      REAL*8 :: ZCRW = 0.235D+0
      REAL*8 :: VCRW = 57.1075D+0
      REAL*8 :: PAFW = 0.344D+0
      REAL*8 :: DPMW = 1.8D+0
      REAL*8 :: RHOGI = 0.D+0
      REAL*8 :: RHOLI = 0.D+0
      REAL*8 :: VISGI = 0.D+0
      REAL*8 :: VISLI = 0.D+0
      REAL*8 :: VISNI = 0.D+0
      REAL*8 :: FSFLA = 1.D+0
      REAL*8 :: TBA = 83.35D+0
      REAL*8, DIMENSION(5) :: PTPS = (/ 0.0765D+0, 0.2664D+0,
     &  0.D+0, 0.00127D+0, 0.10898D+0 /)
      REAL*8 :: SUFW = 72.8D-3
      REAL*8 :: THKRW = 0.6068
      REAL*8 :: THKRA = 26.1D-3
      REAL*8 :: THIRD = 1.D+0/3.D+0
      REAL*8 :: TBW = 373.2D+0
      REAL*8 :: DFGAC = 0.D+0
      REAL*8 :: DFGOC = 0.D+0
      REAL*8 :: DFGNC = 0.D+0
      REAL*8 :: DFGWC = 0.D+0
      REAL*8 :: DFLAC = 0.D+0
      REAL*8 :: DFLOC = 0.D+0
      REAL*8 :: DFLNC = 0.D+0
      REAL*8 :: DFLSC = 0.D+0
      REAL*8 :: DFNAC = 0.D+0
      REAL*8 :: DFNOC = 0.D+0
      REAL*8 :: DFNNC = 0.D+0
      REAL*8 :: DFNWC = 0.D+0
      REAL*8 :: HCAW = 6.7419D+9
      REAL*8 :: GPI = 3.1415926536D+0
      REAL*8 :: TENTH = 1.D-1
      REAL*8 :: TOLN = LOG(1.D+1)
      INTEGER :: ISMALL = -32000
      INTEGER :: IBIG = 32000
!      INTEGER :: NBYTR = SIZEOF(ZERO)
!      INTEGER :: NBYTI = SIZEOF(ISMALL)
      INTEGER :: NBYTR = 8
      INTEGER :: NBYTI = 4
      INTEGER :: NBYTB = 4
      INTEGER, DIMENSION(4) :: IPTPS = (/ 1, 1, 1, 1 /)
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE COUP_WELL
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
!     Coupled well variables
!
!     ID_CW(1,NCW) - starting well interval index
!     ID_CW(2,NCW) - ending well interval index
!     ID_CW(3,NCW) - starting well node index, pointer to IWN_CW
!     ID_CW(4,NCW) - ending well node index, pointer to IWN_CW
!     ID_CW(5,NCW) - starting field node index, pointer to IWF_CW
!     ID_CW(6,NCW) - ending field node index, pointer to IWF_CW
!     ID_CW(7,NCW) - principal well node index
!     ID_CW(8,NCW) - well flow control index
!     ID_CW(9,NCW) - starting frac. triangle index, pointer to IBHT_FRC
!     ID_CW(10,NCW) - ending frac. triangle index, pointer to IBHT_FRC
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 13 December 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: DNR_CW,DP_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: FF_CW
      REAL*8, DIMENSION(:), ALLOCATABLE :: FX_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: FXA_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: FXS_CW,FXW_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: P_CW
      REAL*8, DIMENSION(:), ALLOCATABLE :: PF_CW
      REAL*8, DIMENSION(:), ALLOCATABLE :: PL_CW
      REAL*8, DIMENSION(:), ALLOCATABLE :: RHOF_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: Q_CW,QM_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RS_CW
      REAL*8 :: RSD_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PAR_CW
      REAL*8, DIMENSION(:), ALLOCATABLE :: TML_CW
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: VAR_CW
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: VARC_CW,VARSP_CW
      REAL*8, DIMENSION(:), ALLOCATABLE :: PLX_CW
      REAL*8, DIMENSION(:), ALLOCATABLE :: PLY_CW
      REAL*8, DIMENSION(:), ALLOCATABLE :: PLZ_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XP_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: YP_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ZP_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XTP_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: YTP_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ZTP_CW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: T_CW
      REAL*8, DIMENSION(:), ALLOCATABLE :: TF_CW
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: WNM_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICC_CW
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ID_CW,IMP_CW,ITS_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IM_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: INV_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IREF_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IS_CW
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISOLU_CW,ISPC_CW
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISOLC_CW,ISOLK_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IT_CW,NSP_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWFG_CW,IWF_CW,IXPG_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWNG_CW,IWN_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWP_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: JM_CW
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: KLU1_CW
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU2_CW
      INTEGER :: N_CW,N_GLW,NSD_CW,NWF_CW,NWN_CW
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVH
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
!     Hydrate field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 21 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SH
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVP
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
!     Mechanical, hydraulic, saturation function, aqueous relative
!     permeability, gas relative permeability, and nonaqueous liquid
!     relative permeability variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: T,PL,PG,PN,PSO,DNR
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RHOG,RHOL,BTGL
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XGA,XGW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XMGA,XMGW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XLA,XLW,XLO
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XMLA,XMLW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RKG,RKL
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PVA,PVW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: DFLA,DFLS
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SG,SL,SS,SI
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: YLS,XLS,XMLS,RHOSP,HSP,TMS
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RHOML,RHOMG
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: DFGA,DFGW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PORD,PORT,PERMRF,POR0
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: TORG,TORL
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HGA,HG,UEG,UEL
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HGW,HL,HLW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: VISG,VISL
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: THKG,THKL
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RSDL
      REAL*8, DIMENSION(:), ALLOCATABLE :: SDPF,SDPM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NPHAZ
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FILES
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
!     External file name and unit number variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 17 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: IRD = 21
      INTEGER :: IWR = 22
      INTEGER :: IPL = 23
      INTEGER :: IRS = 24
      INTEGER :: IPLX = 25
      INTEGER :: IPLY = 26
      INTEGER :: IPLZ = 27
      INTEGER :: IOFFSET
      INTEGER :: IPLSKP = 0
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FLUX
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
!     Flux variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UL,VL,WL
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UG,VG,WG
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UDGA,VDGA,WDGA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UGA,VGA,WGA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UDGW,VDGW,WDGW
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UGW,VGW,WGW
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UDLA,VDLA,WDLA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ULA,VLA,WLA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UDLW,VDLW,WDLW
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ULW,VLW,WLW
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UDS,VDS,WDS
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: US,VS,WS
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE GEO_MECH
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
!     Geomechanical variables
!
!     BC_GM(LBCV_GM,LBTM_GM,LBCIN_GM) geomechanical boundary condition
!       values
!     IBCC_GM(LBC_GM) geomechanical boundary condition cyclic index
!     IBCT_GM(3,LBC_GM) geomechanical boundary condition type index
!     IBCM_GM(LBC_GM) geomechanical boundary condition time index
!     GRAV_GM acceleration of gravity for geomechanics, m/s^2
!     EPS_GM(6,LFDM) strain at finite element centroid
!       1 xx
!       2 yy
!       3 zz
!       4 yz (tensorial strain, i.e., 1/2 (dv/dz + dw/dy) )
!       5 xz (tensorial strain, i.e., 1/2 (du/dz + dw/dx) )
!       6 xy (tensorial strain, i.e., 1/2 (du/dy + dv/dx) )
!     SIG_GM(6,LFDM) stress at finite element centroid, Pa
!       1 xx
!       2 yy
!       3 zz
!       4 yz
!       5 xz
!       6 xy
!     ND_GM(8,LFDM) finite element node number at hexahedral position
!       number and node number, initially read as global FE node
!       numbers, but converted to local FE node numbers
!       1 i,j,k
!       2 i+1,j,k
!       3 i,j+1,k
!       4 i+1,j+1,k
!       5 i,j,k+1
!       6 i+1,j,k+1
!       7 i,j+1,k+1
!       8 i+1,j+1,k+1  ---
!     NE_GM(8,LFEN) node number at hexaderal position number and finite
!       element number, initially read as global node numbers, but
!       converted to local node numbers
!       1 i,j,k
!       2 i+1,j,k
!       3 i,j+1,k
!       4 i+1,j+1,k
!       5 i,j,k+1
!       6 i+1,j,k+1
!       7 i,j+1,k+1
!       8 i+1,j+1,k+1  ---
!     U_GM(2,NFEN) x (global coordinates) displacements, m
!     V_GM(2,NFEN) y (global coordinates) displacements, m
!     W_GM(2,NFEN) z (global coordinates) displacements, m
!       1 reference displacement
!       2 current displacement
!     P_GM(3,LFDM) pore pressure, Pa
!       1 old time step value
!       2 k iterate level for seq. coupling of flow and geomechanics
!       3 k+1 iterate level for seq. coupling of flow and geomechanics
!     EPSV_GM(3,LFDM) volumetric strain
!       1 old time step value
!       2 k iterate level for seq. coupling of flow and geomechanics
!       3 k+1 iterate level for seq. coupling of flow and geomechanics
!     EPSV_CMP(LFDM) reference volumetric strain
!     SIGV_GM(3,LFDM) volumetric stress
!       1 old time step value
!       2 k iterate level for seq. coupling of flow and geomechanics
!       3 k+1 iterate level for seq. coupling of flow and geomechanics
!     SIGV_CMP(LFDM) reference volumetric stress
!     RSD_GM maximum residual for convergence of coupled flow and 
!       transport and geomechanics
!     NSD_GM node at which maximum residual occurs
!     K_GM iterate level for sequential flow and geomechanics
!     IREF_GM index for reference state calculations
!     IHCM_GM(LRC) index for hydrate composite model
!       1 C-Factor model
!     IPROP_GM(3,LRC) index for linked geomechanics functions
!       1 porosity and mean effective stress
!       2 permeability and mean effective stress
!       3 capillary pressure and mean effective stress
!     PROP_GM(9,LRC) geomechanical properties and linked geomechanics 
!       function parameters
!       1 Young's modulus, Pa
!       2 Poisson's ratio
!       3 Biot's parameter
!       4 Coefficent of thermal expansion, 1/K
!       5 Davis-Davis Residual Porosity at High Stress
!       6 Davis-Davis Porosity-Mean Stress Function Exponent
!       7 Davis-Davis Intrinsic Permeability-Porosity Function Exponent
!     NBC_GM(NP) number of geomechanical boundary condition finite-
!       element nodes on each processor
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 October 2021.
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: BC_GM
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PROP_GM
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: U_GM,V_GM,W_GM
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: EPS_GM,SIG_GM
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: P_GM,EPSV_GM,SIGV_GM
      REAL*8, DIMENSION(:), ALLOCATABLE :: EPSV_CMP,SIGV_CMP
      REAL*8, DIMENSION(:), ALLOCATABLE :: RSDM_GM
      REAL*8 :: GRAV_GM,RSD_GM
      INTEGER :: IREF_GM
      INTEGER :: NFEN_GM,MD_GM,MK_GM,ML_GM,MU_GM,NSD_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBC_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCC_GM,IBCM_GM,IBCN_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCIN_GM,IBCD_GM,IBCR_GM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IBCT_GM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IPROP_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IM_GM,IHCM_GM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IDP_GM,JDP_GM,KDP_GM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ND_GM,NE_GM,NK_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: NF_GM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: MLU_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: NLU_GM
      INTEGER, DIMENSION(1:2) :: K_GM = 0
      INTEGER, DIMENSION(:), ALLOCATABLE :: NFNGN,IGHE,IGHN
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NCGN,NPGN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NLSGN,NLRGN,NGHN
      INTEGER :: IEQ_OFFSET_GM,NNZ_GM,NFNGN_G
      INTEGER, DIMENSION(:), ALLOCATABLE :: NUKGL,NUKGO,NUKGP
      INTEGER :: NUKGG
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE GLB_PAR
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
!     Parameters
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: LNOTES=1,LEPD=1,LINC=100
      INTEGER :: LFX=1,LFY=1,LFZ=1
      INTEGER :: LPX_MPI=1,LPY_MPI=1,LPZ_MPI=1,LP_MPI=1
      INTEGER :: LFX_MPI=1,LFY_MPI=1,LFZ_MPI=1
      INTEGER :: LAN=1,LAD=1,LMNP=1
      INTEGER :: LT=0,LL=1,LG=0,LN=0,LC=0,LFW=0,LS=0,LD=0
      INTEGER :: LPC=0,LALC=0,LWELL=0,LDCO2=0,LM=0
      INTEGER :: LBD=0,LSP=1,LPT=0,LIS=0
      INTEGER :: LBC=1,LBCIN=1,LBTM=1
      INTEGER :: LBC_GM=1,LBCIN_GM=1,LBTM_GM=1,LBCV_GM=10
      INTEGER :: LSR=1,LSTM=1,LNW=1,LNWT=1,LNWS=1
      INTEGER :: L_CW=0,LN_CW=1,LWI_CW=1,LWT_CW=1,LWTP_CW=1
      INTEGER :: LWN_CW=0,LWF_CW=0,LUK_CW=0,LSOLU_CW=1,LSPC_CW=1
      INTEGER :: L_LW=0,LN_LW=1,LWI_LW=1,LWN_LW=0,LWF_LW=0
      INTEGER :: LRC=1,LSOLU=1,LCN=1,LCDC=1,LCDP=1,LCDS=1
      INTEGER :: LREF=1,LPTM=1,LSF=1,LSFDOM=1
      INTEGER :: LOBDT=1,LOBDS=1
      INTEGER :: LTBL=1,LCHEM=1,LRK=1,LUGR=0
      INTEGER :: LPTA=0,LPLANT=1,LSW=0,LATM=1,LNNGC=1
      INTEGER :: LSFCA=1,LSFCC=1,LSFCP=1,LSFCT=1,LSFCN=1,LUK_SFC=1
      INTEGER :: L_DP=0,LFD_DP=1,L_EC=0,LFD_EC=1,LBC_EC=1,L_SFC=0
      INTEGER :: LBAL=1,LREM=1,LREL=1,LHYD=0,LR=0
      INTEGER :: LANW=1,LGC=0,LN2=0,LNHC=1
      INTEGER :: LUK=1,LPH=1,LCMP=1,LSALC=1,LMPH=0
      INTEGER :: LFXY=1,LFYZ=1,LFZX=1
      INTEGER :: LFD=1,LBR=0,LRFN=0,LSTC=7
      INTEGER :: LHBW=1,LJA=1,LJB=1,LJC=1,LJC_GM=1,LJD=1
      INTEGER :: LJE=1,LJF=1,LJG=1,LJG_GM=1,LJH=1,LJH_GM=1
      INTEGER :: LJI=1,LJJ=1,LJK=1,LJL=1,LJM=1,LJN=1
      INTEGER :: LJO=1,LJO_GM=1,LSU=1
      INTEGER :: LSV=1,LSFV=1,LSFVGC=1
      INTEGER :: LFDT=1,LFDL=1,LFDG=1,LFDN=1,LFDN2=1
      INTEGER :: LFDC=1,LFDI=1,LFDS=1,LFDCR=1
      INTEGER :: LFDR=1,LFDRL=1,LFDRG=1,LFDRN=1
      INTEGER :: LFDD=1,LFDA=1,LFDH=1,LFDNH=1
      INTEGER :: LFDGC=1,LFDM=1,LFEN=1
      INTEGER :: LSX=0,LSRX=0
      INTEGER :: LSY=0,LSRY=0
      INTEGER :: LSZ=0,LSRZ=0
      INTEGER :: LSXT=1,LSXL=1,LSXG=1,LSXN=1,LSXN2=1
      INTEGER :: LSXC=1,LSXS=1,LSXGC=1,LSXLC=1,LSXNC=1
      INTEGER :: LSYT=1,LSYL=1,LSYG=1,LSYN=1,LSYN2=1
      INTEGER :: LSYC=1,LSYS=1,LSYGC=1,LSYLC=1,LSYNC=1
      INTEGER :: LSZT=1,LSZL=1,LSZG=1,LSZN=1,LSZN2=1
      INTEGER :: LSZC=1,LSZS=1,LSZGC=1,LSZLC=1,LSZNC=1
      INTEGER :: LRCT=1,LRCL=1,LRCG=1,LRCN=1
      INTEGER :: LRCS=1
      INTEGER :: LBCT=1,LBCL=1,LBCG=1,LBCN=1,LBCN2=1
      INTEGER :: LBCC=1,LBCI=1,LBCS=1
      INTEGER :: LBCA=1,LBCH=1,LBCGC=1
      INTEGER :: LBCU=1,LBCV=1
      INTEGER :: LOUPV=1,LVREF=1,LVPLOT=1,LFILES=1
      INTEGER :: LSCHR=18,LRPLC=12,LRPGC=6,LRPNC=6,LRPL=4
      INTEGER :: LNWN=1,LSZW=1,LWSI=1,LWTI=1
      INTEGER :: LNWV=1,LUKW=1
      INTEGER :: LP_TA=1,LT_TA=1,L_LV=1,LINH=15
      INTEGER :: LT_PH=155,LO_PH=11
      INTEGER :: LT_TH=100,LO_TH=11
      INTEGER :: LPOLYN=1,LPOLYC=4
      INTEGER :: LEQC=1,LEQE=1,LEQK=1
      INTEGER :: LRCE=1,LRCK=1,LREK=1
      INTEGER :: LSEC=1,LSEE=1,LSEK=1
      INTEGER :: LSPG=0,LSPK=1,LSPT=0,LCKN=1,LMC=1
      INTEGER :: LSPL=0,LSPN=0,LSPS=0
      INTEGER :: LSPE=0,LESITE=0,LSPLK=0
      INTEGER :: LSPR=0,LSPBC=0,LSOLSR=1
      INTEGER :: LXYZG=0
      INTEGER :: LNGC=1
      INTEGER :: LSPILL=0
      INTEGER :: LANI=1,LCAT=1,LNEU=1,LNAF=1,LNCF=1,LNNF=1,LMCG=1
      INTEGER :: LHF_HT=1,LCN_HT=1,LCP_HT=1,LPE_HT=1,LTP_HT=1,LPP_HT=1
      INTEGER :: LCH_HT=1,LHE_HT=1
      INTEGER :: LPF_EOR=1,LPCF=1
      INTEGER :: LBC_FRC=1,LF_FRC=1,LT_FRC=1,LNC_FRC=1,LFC_FRC=1
      INTEGER :: LTC_FRC=1,LSR_FRC=1,LSTM_FRC=1,LT_FRCC=1,LF_FRCC=1
      INTEGER :: L_FRC=0,LJG_FRC=1,LJH_FRC=1,LJK_FRC=1,LJL_FRC=1
      INTEGER :: LXP_FRC=0
      INTEGER :: LVIC_FRC=21
      INTEGER :: LJG_BH=1,LJH_BH=1,LJC_BH=1
      INTEGER :: LJK_BH=1,LJL_BH=1,LJN_BH=1
      INTEGER :: LJG_MCB=1,LJH_MCB=1,LJK_MCB=1,LJL_MCB=1
      INTEGER :: LJG_BCM=1,LJH_BCM=1,LJK_BCM=1,LJL_BCM=1
      INTEGER :: LJG_BCF=1,LJH_BCF=1,LJK_BCF=1,LJL_BCF=1
      INTEGER :: LJG_MCF=1,LJH_MCF=1,LJK_MCF=1,LJL_MCF=1
      INTEGER :: LJG_FCM=1,LJH_FCM=1,LJK_FCM=1,LJL_FCM=1
      INTEGER :: LJG_FCB=1,LJH_FCB=1,LJK_FCB=1,LJL_FCB=1
      INTEGER :: L_BH=0,LN_BH=1,LI_BH=1,LT_BH=1
      INTEGER :: LN_BHC=1,LBN_BHC=1
      INTEGER :: LBN_BH=0,LUK_BH=0,LSOLU_BH=1
      INTEGER :: LSR_BH=1,LSTM_BH=1
      INTEGER :: LFC_BH=0,LBC_BH=0,LCOAX_BH=0
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE GRID
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
!     Grid variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XE,YE,ZE
      REAL*8, DIMENSION(:), ALLOCATABLE :: RP,XP,YP,ZP
      REAL*8, DIMENSION(:), ALLOCATABLE :: DXGF,DYGF,DZGF
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: DXGP,DYGP,DZGP
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: AFX,AFY,AFZ
      REAL*8, DIMENSION(:), ALLOCATABLE :: VOL
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: GRVX,GRVY,GRVZ
      REAL*8, DIMENSION(:), ALLOCATABLE :: SBFB,SBFS,SBFW,SBFE,SBFN,SBFT
      REAL*8, DIMENSION(:), ALLOCATABLE :: RBFB,RBFS,RBFW,RBFE,RBFN,RBFT
      INTEGER :: ID,NP,IDV,NPV,NFCGC_G,NFC_G,NFC_L
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGHC,IXP,ND,NFCGC,NFC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICM,IDP,JDP,KDP,INBS
      INTEGER, DIMENSION(:), ALLOCATABLE :: NUKFL,NUKFO,NUKTL,NUKTO
      INTEGER, DIMENSION(:), ALLOCATABLE :: NUKFP,NUKTP
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NCGC,NPGC,ITAG
      INTEGER, DIMENSION(:), ALLOCATABLE :: NLSGC,NLRGC,NGHC
      INTEGER :: NUKFG,NUKTG,IFLD,JFLD,KFLD,NFLD,NXP,N_DB
      INTEGER :: ICS,IPFLD,JPFLD,KPFLD,NPFLD
      INTEGER :: MPI_COMM_VERT
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE HYST
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
!---  Hysteretic k-s-P function variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 23 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: ASL
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ASLMIN
      REAL*8, DIMENSION(:), ALLOCATABLE :: ASGT
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SGT
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE JACOB
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
!---  Jacobian matrix variables  ---
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 26 August 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: BLU
      REAL*8, DIMENSION(:), ALLOCATABLE :: DLU
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU,KLUC
      INTEGER, DIMENSION(:), ALLOCATABLE :: MLU,MLUC
      INTEGER, DIMENSION(:), ALLOCATABLE :: NLU,NLUC
      INTEGER :: IEQ_OFFSET,IEQC_OFFSET
      INTEGER :: NNZ,NNZC,NNZFR,NNZTR,NNZGR
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE LEAK_WELL
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
!     Coupled well variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 January 2011.
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: FF_LW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PAR_LW
      REAL*8, DIMENSION(:), ALLOCATABLE :: PLX_LW
      REAL*8, DIMENSION(:), ALLOCATABLE :: PLY_LW
      REAL*8, DIMENSION(:), ALLOCATABLE :: PLZ_LW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XP_LW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: YP_LW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ZP_LW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: XTP_LW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: YTP_LW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ZTP_LW
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: WNM_LW
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ID_LW
      INTEGER, DIMENSION(:), ALLOCATABLE :: INV_LW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IS_LW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IT_LW,IZ_LW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWF_LW
      INTEGER, DIMENSION(:), ALLOCATABLE :: NF_LW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWP_LW
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU_LW
      INTEGER, DIMENSION(:), ALLOCATABLE :: KLUC_LW
      INTEGER, DIMENSION(:), ALLOCATABLE :: ND_LW
      INTEGER :: N_LW,NWN_LW,NWF_LW
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE OUTPU
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
!     Output control variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 12 August 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: PRTM
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SF
      REAL*8, DIMENSION(:), ALLOCATABLE :: CNVPLOT,CNVREF
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CNVSF
      REAL*8 :: CNVTM,CNVLN
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPLOT
      INTEGER, DIMENSION(:), ALLOCATABLE :: IREF
      INTEGER, DIMENSION(:), ALLOCATABLE :: NDREF,NDREFL,NDREFI
      INTEGER :: NPRTM
      INTEGER :: NVPLOT
      INTEGER :: NREF,NREFL,IOFFSET_REF
      INTEGER :: NVREF
      INTEGER :: ICNO
      INTEGER :: ICNS
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFT
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFF
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFD
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFN,ISFS,NSFN,ISFB
      INTEGER :: NSF
      INTEGER :: NSFGP
      INTEGER :: IHSF
      INTEGER :: IFQS
      INTEGER :: IFQO
      INTEGER :: ISGNS
      INTEGER :: ISGNO
      INTEGER :: ISGNP
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFGP
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISF
      INTEGER	(KIND=MPI_OFFSET_KIND), DIMENSION(:), ALLOCATABLE :: 
     &  IOFFSET_SF
      CHARACTER(64) :: UNTM
      CHARACTER(64) :: UNLN
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: UNPLOT
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: UNREF
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: CHREF
      CHARACTER(64), DIMENSION(:,:), ALLOCATABLE :: UNSF
      CHARACTER(64), DIMENSION(:,:), ALLOCATABLE :: CHSF
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: FNSF
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE NCG_PT
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
!     Noncondensible gas property table variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 September 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: P_TA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: T_TA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: RHO_TA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: H_TA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: U_TA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: FUG_TA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: S_TA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: RHO_ST
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: H_ST
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: U_ST
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: FUG_ST
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: S_ST
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: T_LV
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: P_LV
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RHOL_LV
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HL_LV
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: UL_LV
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RHOV_LV
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HV_LV
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: UV_LV
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: FUG_LV
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SL_LV
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SV_LV
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPTP,IPCR
      INTEGER, DIMENSION(:), ALLOCATABLE :: IP_TA
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IT_TA
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IV_TA
      INTEGER :: INCG
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_LV
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE PROP
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
!     Mechanical, hydraulic, saturation function, aqueous relative
!     permeability, gas relative permeability, and nonaqueous liquid
!     relative permeability variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: CPS,RHOS,PCMP,TCMP
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CMP,POR,TOR
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PERM,SCHR
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RPGC,RPLC,THKS
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISLTBL,IRGTBL,IRLTBL
      INTEGER, DIMENSION(:), ALLOCATABLE :: ITOR,IPRF,ISCHR,ISM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IRPG,IRPL,IZ
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE PTZRCOEF
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
!     Pitzer activity coefficients
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 February 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: B0
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: B1
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: B2
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CMXX
      REAL*8, DIMENSION(:), ALLOCATABLE :: TCC
      REAL*8, DIMENSION(:), ALLOCATABLE :: TAA
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PSIC
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PSIA
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ALAMB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CLAMB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ELAMB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HOLAMB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: BPPR
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: BPHI
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: BPR
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: BMMX
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ATB0
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ATB1
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ATB2
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ATCMX
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ATNLAM
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ATCLAM
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ATALAM
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ATHLAM
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ATTC
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ATPC
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ATTA
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: ATPA
      REAL*8, DIMENSION(:), ALLOCATABLE :: CTCPH
      REAL*8, DIMENSION(:), ALLOCATABLE :: CTC
      REAL*8, DIMENSION(:), ALLOCATABLE :: CTCPR
      REAL*8, DIMENSION(:), ALLOCATABLE :: CTCPPR
      REAL*8, DIMENSION(:), ALLOCATABLE :: CTAPH
      REAL*8, DIMENSION(:), ALLOCATABLE :: CTA
      REAL*8, DIMENSION(:), ALLOCATABLE :: CTAPR
      REAL*8, DIMENSION(:), ALLOCATABLE :: CTAPPR
      REAL*8, DIMENSION(:), ALLOCATABLE :: ETH
      REAL*8, DIMENSION(:), ALLOCATABLE :: ETHP
      REAL*8, DIMENSION(:), ALLOCATABLE :: ETHP2
      INTEGER, DIMENSION(:), ALLOCATABLE :: IDD_PZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: JPA
      INTEGER, DIMENSION(:), ALLOCATABLE :: JPC
      INTEGER, DIMENSION(:), ALLOCATABLE :: JPN
      INTEGER :: NCC_PZ
      INTEGER :: NA_PZ
      INTEGER :: NNN_PZ
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE REACT
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
!     Reaction variables
!
!     ACTVC - constant activity coefficient
!     EQ_C(LSEC,LEQC) - conservation equation parameters
!       species coefficients, conservation equations
!     EQ_E(LSEE,LEQE) - 1 equation parameters
!       species coefficients, 1 equations
!     EQ_K(LSEK+LREK,LEQK) - kinetic equation parameters
!       species coefficients + reaction coefficients, kinetic equations
!     RC_E(5,LRCE) - 1 reaction parameters
!     RC_K(LSPK+11,LCKN,LRCK) - kinetic reaction parameters
!     SP_L(3,LSPL) - aqueous species parameters
!     SP_S(2,LSPS) - solid species parameters
!     RS_S(1,LSPS,LFDR) - initial mineral (solid species) area
!     RS_S(2,LSPS,LFDR) - initial mineral (solid species) volume frac.
!     RS_S(3,LSPS,LFDR) - current mineral (solid species) volume frac.
!     IACTV - activity coefficient model option
!       0 - B-Dot Equation
!       1 - Davies Equation
!       2 - Pitzer Equation
!       3 - Constant
!     NEQC - number of conservation equations
!     NEQE - number of 1 equations
!     NEQK - number of kinetic equations
!     NRCE - number of 1 reactions
!     NRCK - number of kinetic reactions
!     NSPC - number of component species
!     NSPG - number of gas species
!     NSPK - number of kinetic species
!     NSPLK - number of linked species
!     NSPL - number of aqueous species
!     NSPN - number of NAPL species
!     NSPR - number of reactive species
!     NSPS - number of solid species
!     NRTSI - number of reactive transport sequence iterations
!     SPNMC(LEQC) - conservation component species name
!     SPNMK(LEQC) - kinetic component species name
!     SPNMG(LSPG) - gas species name
!     SPNML(LSPL) - aqueous species name
!     SPNMN(LSPN) - NAPL species name
!     SPNMS(LSPS) - solid species name
!     SPNME(LSPE) - exchanged species name
!     YSPG(LFDRG,LEQC+LEQK) - gas transport fraction
!     YSPL(LFDRL,LEQC+LEQK) - aqueous transport fraction
!     YSPN(LFDRN,LEQC+LEQK) - NAPL transport fraction
!     SP_C(LFDR,LSPR) - species concentration, mol/m^3 node volume
!     SP_CO(LFDR,LSPR) - old species concentration, mol/m^3 node volume
!     SP_CI(LFDR,LSPR) - initial species conc., mol/m^3 node volume
!     SP_CMN(LFDR,LSPS) - base mineral species conc., mol/m^3 node vol.
!     SP_CBO(LBCC,LSPR) - old species bound. conc., mol/m^3 node volume
!     SP_RATE(LFDR,LSPS) - mineral reaction rate, mol/s
!     SP_AREA(LFDR,LSPS) - mineral surface area, m^2
!     IC_SP(LFDR,LSPR) - species initial condition type index
!     ISP_MN(LSPR) - species mineral index
!     RCNME(LRCE) - 1 reaction name
!     RCNMK(LRCK) - kinetic reaction name
!     IEQ_C(LSEC+1,LEQC) - integer conservation equation parameters
!       number of species + species number + ... +
!     IEQ_E(LSEE+1,LEQE) - integer 1 equation parameters
!       number of species + species number + ... +
!     IEQ_K(LSEK+LREK+2,LEQK) - integer kinetic equation parameters
!       number of species + species number + ... +
!       number of reactions + reaction number + ..., kinetic equations
!     IEQ_S(LSPR) - species --> equation sequence
!     ISP_E(LSPE) - exchange site
!     IEL_LK(LSPE) - exchanged/cation species link
!     ISP_S(LEQE+LEQC+LEQK) - equation --> species sequence
!     IRC_K(LSPK+3,LRCK) - integer kinetic reaction parameters
!     IRCKN(LSPK+11) - flag to read spatially variable kinetic 
!       reaction parameters
!       0 - No spatial variation
!       1 - Rxn parameters assigned on a node by node basis
!     IRCKT(LRCK) - kinetic reaction type
!       0 - dissolution-precipitation
!       1 - forward-backward
!     CFMX(1:LMC) - spatially variable mixing coefficient; default=1
!     ISPLK(14) - conservation species pointer
!     ISPGL(LSPG) - gas species to associate aqueous species pointer
!     ISP_OW(LSPS,LFDR) - solid species lithology card overwrite
!     LEQC - number of conservation equations 
!     LEQE - number of 1 equations
!     LEQK - number of kinetic equations 
!     LRCE - number of 1 reactions
!     LRCK - number of kinetic reactions 
!     LREK - number of kinetic equation reactions 
!     LSEC - number of conservation equation species
!     LSEE - number of 1 equation species
!     LSEK - number of kinetic equation species 
!     LSPG - number of gas species 
!     LSPK - number of kinetic reaction species 
!     LCKN - number of kinetic reaction parameters
!     LSPL - number of aqeuous species
!     LSPN - number of NAPL species
!     LSPS - number of solid species
!     LSPR - number of reactive species
!     LCAT - number of cations (pitzer)
!     LANI - number of anions (pitzer)
!     LNEU - number of neutrals (pitzer)
!     LNAF - array size for anion interaction parameters (pitzer)
!     LNCF - array size for cation interaction parameters (pitzer)
!     LNFN -array size for neutral interaction parameters (pitzer)
!     LMCG - maximum absolute value of species charge (pitzer)
!     CMIN - minimum concentration for eckechem
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 17 January 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 :: ACTVC,DT_RST,DTI_RST
      LOGICAL :: ECKE_ER
      REAL*8, DIMENSION(:), ALLOCATABLE :: ACTV16
      REAL*8, DIMENSION(:), ALLOCATABLE :: C_PH
      REAL*8, DIMENSION(:), ALLOCATABLE :: CHARG
      REAL*8, DIMENSION(:), ALLOCATABLE :: FACTV,ACTVS
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: EQ_C
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: EQ_E
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: EQ_K
      INTEGER :: IACTEX,IACTV
      INTEGER :: ISP_IEDL
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IC_SP
      INTEGER, DIMENSION(:), ALLOCATABLE :: IEL_LK
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IEQ_C
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IEQ_E
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IEQ_K
      INTEGER, DIMENSION(:), ALLOCATABLE :: IEQ_S
      INTEGER, DIMENSION(:), ALLOCATABLE :: IMMB
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IRC_K
      INTEGER, DIMENSION(:), ALLOCATABLE :: IRCKT
      INTEGER, DIMENSION(:), ALLOCATABLE :: IRCKN
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISP_E
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISP_MN
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISP_OW
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISP_S
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISPLK
      REAL*8, DIMENSION(:), ALLOCATABLE :: CFMX
      INTEGER :: N_RST,NEQC,NEQE,NEQK,NESITE
      INTEGER :: NRCE,NRCK,NRTSI,NSPC,NSPE
      INTEGER :: NSPG,NSPK,NSPL,NSPLK,NSPN
      INTEGER :: NSPR,NSPS
      REAL*8 :: CMIN
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RC_E
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: RC_K
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: RCNME
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: RCNMK
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: RS_S
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SP_C
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SP_CBO
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SP_CMN
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SP_RATE
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SP_AREA
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SP_CO
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SP_CI
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SP_L
      REAL*8 :: SP_MDG,SP_MDL,SP_MDN
      REAL*8, DIMENSION(:), ALLOCATABLE :: SP_SDCL
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SP_S
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNMC
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNME
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNMG
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNMK
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNML
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNMS
      REAL*8 :: TM_RST
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: YSPG
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: YSPL
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE SOLTN
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
!     Solution control variables
!
!     I_ERR(1) = error integer
!     I_ERR(2) = location of error message real
!       0 = no error message real
!       1 = real between error string 1 and 2
!       2 = real between error string 2 and 3
!       3 = real between error string 3 and 4
!     I_ERR(3) = location of error message integer
!       0 = no error message integer
!       1 = integer between error string 1 and 2
!       2 = integer between error string 2 and 3
!       3 = integer between error string 3 and 4
!     I_ERR(4) = generating processor rank (ID)
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 :: TM,TMMX,TMPR
      REAL*8 :: DT,DTI,DTMX,DTMN,DTAF,DTCF,DTO,DTSO
      REAL*8 :: RSDMX,RLXF,CRNTMXC
      REAL*8 :: R_ERR
      REAL*8, DIMENSION(:), ALLOCATABLE :: TMPS,TMPE,TMPD,TMPX,TMPN
      REAL*8, DIMENSION(:), ALLOCATABLE :: TMPA,TMPC
      REAL*8, DIMENSION(:), ALLOCATABLE :: RSD,RSDM
      REAL*8, DIMENSION(:), ALLOCATABLE :: WFMN
      REAL*8 :: lissv_time
      REAL*8 :: RTOL_PETSC,ATOL_PETSC,DTOL_PETSC
      INTEGER, DIMENSION(1:4) :: I_ERR
      INTEGER :: ICNV = 0
      INTEGER :: IVRSN,ISIC,IEO,ILES,IOM,ICODE
      INTEGER :: IEQT,IEQW,IEQA,IEQN,IEQO,IEQC
      INTEGER :: IEQS,IEQD,IEQDO,IEQHA,IEQHN
      INTEGER :: IEQHO,IEQDA
      INTEGER :: IAQU,IGAS,INAPL
      INTEGER :: NEPD,MEPD,IEPD,NRIMX
      INTEGER :: NSTEP,NRST,NITER,NTSR,NGC,MXSTEP
      INTEGER :: IUNM,IUNKG,IUNS,IUNK,IUNMOL
      INTEGER :: ISVC,ISVF,IEDLS
      INTEGER, DIMENSION(:), ALLOCATABLE :: IEQGC
      INTEGER, DIMENSION(1:100) :: ISLC = 0
      INTEGER, DIMENSION(1:20) :: IDMN = 1
      INTEGER, DIMENSION(:), ALLOCATABLE :: ITOFF,NRIM
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSD
      INTEGER, DIMENSION(:), ALLOCATABLE :: MNOD
      INTEGER, DIMENSION(:), ALLOCATABLE :: MADJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: MFLX
      INTEGER, DIMENSION(:), ALLOCATABLE :: MPOS
      INTEGER, DIMENSION(:), ALLOCATABLE :: MNEG
      INTEGER, DIMENSION(:), ALLOCATABLE :: MPOSB
      INTEGER, DIMENSION(:), ALLOCATABLE :: MNEGB
      INTEGER :: ISUB_LOG
      INTEGER :: MAXITS_PETSC
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: HCNM
      CHARACTER(128), DIMENSION(1:4) :: M_ERR
      CHARACTER(32), DIMENSION(1:32) :: SUB_LOG
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE SOURC
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
!     Source variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 June 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: SRC
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SRCW
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SRCA
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SRCS
      REAL*8, DIMENSION(:), ALLOCATABLE :: SRCIW
      REAL*8, DIMENSION(:), ALLOCATABLE :: SRCIA
      REAL*8, DIMENSION(:), ALLOCATABLE :: SRCIS
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSR
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISRIN,ISRM,ISRN,ISRT
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE TABL
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
!     Tabular data variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 October 2021.
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: TBLX
      REAL*8, DIMENSION(:), ALLOCATABLE :: TBLY
      REAL*8, DIMENSION(:), ALLOCATABLE :: TBLDDX
      REAL*8, DIMENSION(:), ALLOCATABLE :: TBLDDY
      INTEGER :: NTBL
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE TRNSPT
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
!     Solute transport variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 6 January 2022.
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: C,CO,CNL
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: YL,YG
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CB,CBO
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: YGB,YLB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PCGL
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: PCSL
      REAL*8, DIMENSION(:), ALLOCATABLE :: DISPL,DISPT
      REAL*8, DIMENSION(:), ALLOCATABLE :: CRNTG,CRNTL
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UC,VC,WC
      REAL*8, DIMENSION(:), ALLOCATABLE :: HLF
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SRCIC
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: SDCL
      REAL*8, DIMENSION(:), ALLOCATABLE :: SMDL,SMDG
      REAL*8 :: DT_CRN
      REAL*8 :: DTI_CRN
      REAL*8 :: TM_CRN
      REAL*8, DIMENSION(:), ALLOCATABLE :: CCL_CRN
      REAL*8 :: CRNTMXT
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CHDF
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SOLUT
      INTEGER :: NSOLU
      INTEGER :: IDISP
      INTEGER :: ICRNT
      INTEGER, DIMENSION(:), ALLOCATABLE :: IEDL
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPCL
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICT
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPCGL
      INTEGER, DIMENSION(:), ALLOCATABLE :: N_CRN
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCDS,NBCDP
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IBCDP
      INTEGER :: NBCDS
!
!---  End of module  ---
!
      END MODULE




!---------------------Fortran 90 Module--------------------------------!
!
      MODULE PETSC_STOMP
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
!     PETSc variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 1 December 2022.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE PETSCKSP
!
!----------------------Type Declarations-------------------------------!
!

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


      type(tMat), PUBLIC :: F_MAT
      type(tMat), PUBLIC :: T_MAT
      type(tMat), PUBLIC :: G_MAT
      type(tVec), PUBLIC :: F_SOL_VEC,F_SOL_VEC_S
      type(tVec), PUBLIC :: T_SOL_VEC,T_SOL_VEC_S
      type(tVec), PUBLIC :: G_SOL_VEC,G_SOL_VEC_S
      type(tVec), PUBLIC :: F_RHS_VEC
      type(tVec), PUBLIC :: T_RHS_VEC
      type(tVec), PUBLIC :: G_RHS_VEC
      type(tVecScatter) :: F_SCATTER,T_SCATTER,G_SCATTER
      type(tKSP), PUBLIC :: F_KSP
      type(tKSP), PUBLIC :: T_KSP
      type(tKSP), PUBLIC :: G_KSP
!
!---  End of module  ---
!
      END MODULE


!---------------------Fortran 90 Module--------------------------------!
!
      MODULE DB_PR
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
!     Define double precision.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!---------------------Type Declarations--------------------------!
!
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
!
!---  End of module  ---
!
      END MODULE

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
!     Global boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: BC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PHDL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PHDN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BCXYZG
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XPBC,YPBC,ZPBC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TSBC,XSBC,YSBC,ZSBC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCD
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCN
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCIN
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IBCT
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IBCSP
!      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCLL
!      INTEGER, DIMENSION(:), ALLOCATABLE :: JBCLL
!      INTEGER, DIMENSION(:), ALLOCATABLE :: KBCLL
!      INTEGER, DIMENSION(:), ALLOCATABLE :: MBCLL
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NBHG
      INTEGER :: NBC
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
!     Primary boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PLB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PGB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PSOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PSWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SLB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SGB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PORDB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PORTB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOMGB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOLB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOGB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHONB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORLB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISLB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RKLB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLSB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SIB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PIB
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SDPMB
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SDPFB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TMSB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: POSB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BTGLB,BTGNB,BTNLB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOMLB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOMNB
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE BCVT
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
!     Energy boundary conditions variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKLB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEGB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLSB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKGB
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE BCVG
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
!     Gas boundary conditions variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORGB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISGB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKGB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGAB
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE BCVN
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
!     NAPL boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XSOB
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE BCVI
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
!     Ice boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOIB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKIB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HIB
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE BCVS
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
!     Salt/surfactant boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLSB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLSB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: OECB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YLSB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SSB
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE BCVA
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
!     Alcohol boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XNAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XNNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XNOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XNWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMNAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMNNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMNOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMNWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFNAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFNNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFNOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFNWB
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE BCVGC
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
!     Gas component boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 November 2007.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: PVCB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XGCB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMGCB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XNCB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMNCB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XLCB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMLCB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFGCB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFNCB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFLCB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: HGCB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZGB
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ZMCB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZNB
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE BCVH
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
!     Hydrate boundary condition variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 September 2004.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XHWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XHOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XHAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XHNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOHB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKHB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HHB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SHB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEGAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PHB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZLAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVHAB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVHNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVHOB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVNB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGNB
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
!     Written by MD White, PNNL, 10 January 2011.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DNR_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: FF_CW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: FX_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: FXA_CW,FXO_CW
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: FXC_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: FXE_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: FXS_CW,FXW_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: P_CW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PF_CW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PL_CW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RHOF_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Q_CW,QM_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RS_CW
      REAL(KIND=DP) :: RSD_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PAR_CW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TML_CW
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: VAR_CW
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: VARC_CW,VARSP_CW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PLX_CW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PLY_CW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PLZ_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XP_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YP_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZP_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XTP_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YTP_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZTP_CW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: T_CW
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
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWF_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWN_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWP_CW
      INTEGER, DIMENSION(:), ALLOCATABLE :: JM_CW
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU_CW
      INTEGER :: N_CW,N_GLW,NSD_CW,NWF_CW,NWN_CW
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
!     Written by MD White, PNNL, 10 January 2011.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: FF_LW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PAR_LW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PLX_LW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PLY_LW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PLZ_LW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XP_LW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YP_LW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZP_LW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XTP_LW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YTP_LW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZTP_LW
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
!     Constants
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP) :: SMALL
      REAL(KIND=DP) :: BIG
      REAL(KIND=DP) :: ZERO
      REAL(KIND=DP) :: ONE
      REAL(KIND=DP) :: EPSL
      REAL(KIND=DP) :: GRAV
      REAL(KIND=DP) :: TABS
      REAL(KIND=DP) :: TMX
      REAL(KIND=DP) :: TMN
      REAL(KIND=DP) :: PATM
      REAL(KIND=DP) :: PMX
      REAL(KIND=DP) :: PMN
      REAL(KIND=DP) :: RCU
      REAL(KIND=DP) :: RCW
      REAL(KIND=DP) :: RCA
      REAL(KIND=DP) :: RHORL
      REAL(KIND=DP) :: RHORG
      REAL(KIND=DP) :: VCRA
      REAL(KIND=DP) :: VISRL
      REAL(KIND=DP) :: VISRG
      REAL(KIND=DP) :: TSPRF
      REAL(KIND=DP) :: HDOD
      REAL(KIND=DP) :: WTMA
      REAL(KIND=DP) :: WTMN
      REAL(KIND=DP) :: WTMW
      REAL(KIND=DP) :: WTMS
      REAL(KIND=DP) :: TCRW
      REAL(KIND=DP) :: PCRW
      REAL(KIND=DP) :: ZCRW
      REAL(KIND=DP) :: VCRW
      REAL(KIND=DP) :: PAFW
      REAL(KIND=DP) :: DPMW
      REAL(KIND=DP) :: RHOGI,RHOLI
      REAL(KIND=DP) :: VISGI,VISLI,VISNI,FSFLA
      REAL(KIND=DP) :: TBA
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PTPS
      REAL(KIND=DP) :: SUFW
      REAL(KIND=DP) :: THKRW
      REAL(KIND=DP) :: THKRA
      REAL(KIND=DP) :: THIRD
      REAL(KIND=DP) :: TBW
      REAL(KIND=DP) :: DFGAC,DFGOC,DFGNC,DFGWC
      REAL(KIND=DP) :: DFLAC,DFLOC,DFLNC,DFLSC
      REAL(KIND=DP) :: DFNAC,DFNOC,DFNNC,DFNWC
      REAL(KIND=DP) :: HCAW
      REAL(KIND=DP) :: GPI
      REAL(KIND=DP) :: TENTH
      REAL(KIND=DP) :: TOLN
      INTEGER :: ISMALL
      INTEGER :: IBIG
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPTPS
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE CCP
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
!     Compound critical property variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ZRA
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: VMC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PAF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WTM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: GCPP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: GWPP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BIJNA,BIJTD
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
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: IRD = 21
      INTEGER :: IWR = 22
      INTEGER :: IPL = 23
      INTEGER :: IRS = 24
      INTEGER :: IRS_BF = 34
      INTEGER :: ISC = 6
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISF
      INTEGER :: IBL,IVD,IRL,IRC
      CHARACTER(64) :: FNRD = 'input'
      CHARACTER(64) :: FNWR = 'output'
      CHARACTER(64) :: FNPL = 'plot'
      CHARACTER(64) :: FNRS = 'restart'
      CHARACTER(64) :: FNSR = 'screen'
      CHARACTER(64) :: FNRS_BF = 'restart_bf'
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: FNSF
      CHARACTER(64) :: FNBL,FNVD,FNRL,FNRC
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FLUXP
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
!     Primary flux variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDLN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDLO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ULA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ULN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ULW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDLN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDLO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VLN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDLN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDLO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WLN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WLW
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FLUXT
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
!     Energy flux variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UQ
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VQ
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WQ
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FLUXN
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
!     NAPL flux variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDNA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDNA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDNA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDNN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDNN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDNN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDNO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDNO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDNO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VDNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UNA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VNA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WNA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UNN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VNN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WNN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WNW
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FLUXS
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
!     Salt/surfactant flux variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: US
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WS
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FLUXD
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
!     Dissolved-oil flux variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ULO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VLO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WLO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UNO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VNO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WNO
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FLUXGC
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
!     Gas-component flux variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: UGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: UDGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: UNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: UDNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ULC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: UDLC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: VGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: VDGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: VNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: VDNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: VLC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: VDLC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: WGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: WDGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: WNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: WDNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: WLC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: WDLC
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
!     Parameter global variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!---------------------Type Declarations--------------------------!
!
      TYPE LIST_NODE
        CHARACTER(128) :: LIST_NAME
        TYPE(LIST_NODE), POINTER :: NEXT
      END TYPE LIST_NODE
      TYPE(LIST_NODE), POINTER :: ROCK_PTR,SOLUT_PTR,SFC_PTR,
     &  SPCL_PTR,SPCG_PTR
      TYPE LIST_SCALING
        CHARACTER(128) :: ROCK_NAME
        CHARACTER(128) :: SCALING_NAME
        INTEGER :: ROCK_NUM
        INTEGER :: SCALING_NUM
        TYPE(LIST_SCALING), POINTER :: NEXT
      END TYPE LIST_SCALING
      TYPE(LIST_SCALING), POINTER :: SCALING_PTR
      TYPE LIST_INTERVAL
        INTEGER :: ID1_CW
        INTEGER :: ID2_CW
        TYPE(LIST_INTERVAL), POINTER :: NEXT
      END TYPE LIST_INTERVAL
      TYPE(LIST_INTERVAL), POINTER :: IN_BH_PTR,IN_CW_PTR,IN_LW_PTR
      TYPE LIST_TRANSITION
        REAL*8 :: XTP1_CW
        REAL*8 :: XTP2_CW
        REAL*8 :: YTP1_CW
        REAL*8 :: YTP2_CW
        REAL*8 :: ZTP1_CW
        REAL*8 :: ZTP2_CW
        REAL*8 :: PAR3_CW
        TYPE(LIST_TRANSITION), POINTER :: NEXT
      END TYPE LIST_TRANSITION
      TYPE(LIST_TRANSITION), POINTER :: TP_BH_PTR,TP_CW_PTR,TP_LW_PTR
      TYPE LIST_WELL_NODE
        INTEGER :: IWN_CW
        REAL*8 :: XP1_CW
        REAL*8 :: XP2_CW
        REAL*8 :: YP1_CW
        REAL*8 :: YP2_CW
        REAL*8 :: ZP1_CW
        REAL*8 :: ZP2_CW
        TYPE(LIST_WELL_NODE), POINTER :: NEXT
      END TYPE LIST_WELL_NODE
      TYPE(LIST_WELL_NODE), POINTER :: WN_BH_PTR,WN_CW_PTR,WN_LW_PTR
      INTEGER :: IPF
      INTEGER :: LNOTES=1,LEPD=1,LINC=100
      INTEGER :: LFX=1,LFY=1,LFZ=1
      INTEGER :: LPX_MPI=1,LPY_MPI=1,LPZ_MPI=1,LP_MPI=1
      INTEGER :: LFX_MPI=1,LFY_MPI=1,LFZ_MPI=1
      INTEGER :: LCFGT_MPI=1,LCGT_MPI=1,LCFT_MPI=1
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
      INTEGER :: LUK=1
      INTEGER :: LPH=1,LCMP=1,LSALC=1,LMPH=0
      INTEGER :: LFXY=1,LFYZ=1,LFZX=1
      INTEGER :: LFD=1,LBR=0,LRFN=0,LSTC=7
      INTEGER :: LHBW=1
      INTEGER :: LJA=1
      INTEGER :: LJB=1
      INTEGER :: LJC=1,LJC_GM=1
      INTEGER :: LJD=1
      INTEGER :: LJE=1
      INTEGER :: LJF=1
      INTEGER :: LJG=1,LJG_GM=1
      INTEGER :: LJH=1,LJH_GM=1
      INTEGER :: LJI=1
      INTEGER :: LJJ=1
      INTEGER :: LJK=1
      INTEGER :: LJL=1
      INTEGER :: LJM=1
      INTEGER :: LJN=1
      INTEGER :: LJO=1,LJO_GM=1
      INTEGER :: LSU=1
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
      INTEGER :: LSCHR=1,LRPLC=12,LRPGC=6,LRPNC=6,LRPL=1
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
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XE,YE,ZE
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRFNX,SRFNY,SRFNZ
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XP,YP,ZP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XREF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DXGF,DYGF,DZGF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DXGP,DYGP,DZGP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: AFX,AFY,AFZ
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: VOL
      REAL(KIND=DP) :: GRAVX,GRAVY,GRAVZ,THXZ,THYZ
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: GRVX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: GRVY
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: GRVZ
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XSP     
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: YSP 
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ZSP 
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: AFZSP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TXSP     
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TYSP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PDXSP     
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PDYSP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DXSP     
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DYSP 
      REAL(KIND=DP), DIMENSION(3,3) :: ROTMAT 
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: ROCK
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: ROCK2
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: XREFU
      INTEGER :: IFLD,JFLD,KFLD
      INTEGER :: IP_MPI,JP_MPI,KP_MPI,NP_MPI
      INTEGER, DIMENSION(:), ALLOCATABLE :: NC_MPI,NCN_MPI,NTP_MPI
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ND_MPI
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ID_MPI,JD_MPI,KD_MPI
      INTEGER :: IJFLD
      INTEGER :: JKFLD
      INTEGER :: KIFLD
      INTEGER :: NFBN,NRFN,NFLD,NBRN
      INTEGER :: ICS
      INTEGER :: N_DB
      INTEGER :: NXP
      INTEGER, DIMENSION(:), ALLOCATABLE :: IXP
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KSPS
      INTEGER, DIMENSION(:), ALLOCATABLE :: IZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IZ2
      INTEGER, DIMENSION(:), ALLOCATABLE :: IXF,INP
      INTEGER :: NROCK
      INTEGER :: NROCK2
      INTEGER :: NDIM
      INTEGER, DIMENSION(:), ALLOCATABLE :: MDIM
      INTEGER, DIMENSION(:), ALLOCATABLE :: ID
      INTEGER, DIMENSION(:), ALLOCATABLE :: JD
      INTEGER, DIMENSION(:), ALLOCATABLE :: KD
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ND
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSX,NSY,NSZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSSX,NSSY,NSSZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IXREF
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IBR
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ICM
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ICMS
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INBS
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
!     Hysteretic k-s-P function variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ASL
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: AST
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ASLMIN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ASTMAX
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ASTMIN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ASNT
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ASNR
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ASGT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SGT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SNR
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SNT
      REAL(KIND=DP) :: BGN,BGL,BHL,BIL,BNL
      REAL(KIND=DP) :: CA_GN,CA_GL,CA_NL
      REAL(KIND=DP) :: SIG_GN,SIG_GL,SIG_NL,SIG_HL,SIG_IL
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NPHAZ
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ASGTL
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ASGTN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SGTL
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SGTN
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IPH
      INTEGER :: ISNR
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
!     Linear system variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ALU
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BLU
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CLU
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DLU
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RSDL
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IM
      INTEGER, DIMENSION(:), ALLOCATABLE :: ILU
      INTEGER, DIMENSION(:), ALLOCATABLE :: JLU
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: JM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU
      INTEGER, DIMENSION(:), ALLOCATABLE :: MLU
      INTEGER, DIMENSION(:), ALLOCATABLE :: NLU
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLUC
      INTEGER, DIMENSION(:), ALLOCATABLE :: MLUC
      INTEGER, DIMENSION(:), ALLOCATABLE :: NLUC
      INTEGER :: ISVC
      INTEGER :: ISVT
      INTEGER :: ISVF
      INTEGER :: MKC
      INTEGER :: MKT
      INTEGER :: MUC
      INTEGER :: MLC
      INTEGER :: MDC
      INTEGER :: MUT
      INTEGER :: MLT
      INTEGER :: MDT
      INTEGER :: NNZFR,NNZTR
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE NAPL
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
!     Oil, alcohol, and surfactant critical property and coefficient
!     variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP) :: PCRO
      REAL(KIND=DP) :: TCRO
      REAL(KIND=DP) :: ZCRO
      REAL(KIND=DP) :: VCRO
      REAL(KIND=DP) :: TFPO
      REAL(KIND=DP) :: WTMO
      REAL(KIND=DP) :: RCO
      REAL(KIND=DP) :: PAFO
      REAL(KIND=DP) :: TBO
      REAL(KIND=DP) :: TTPO
      REAL(KIND=DP) :: DPMO
      REAL(KIND=DP) :: PCCVO
      REAL(KIND=DP) :: WHBTO
      REAL(KIND=DP) :: RHORO
      REAL(KIND=DP) :: TDRO
      REAL(KIND=DP) :: ZRAO
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: VISCO
      REAL(KIND=DP) :: VISRN
      REAL(KIND=DP) :: RHORN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SATOC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CPOC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CIMTC
      REAL(KIND=DP) :: HCOW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: VISCS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TERDC
      REAL(KIND=DP) :: RCS
      REAL(KIND=DP) :: TCRS
      REAL(KIND=DP) :: PCRS
      REAL(KIND=DP) :: VCRS
      REAL(KIND=DP) :: ZCRS
      REAL(KIND=DP) :: ZCRA
      REAL(KIND=DP) :: PAFS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SFCSF
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PCSLS
      REAL(KIND=DP) :: TFPA
      REAL(KIND=DP) :: DPMA
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CPAC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SATAC
      REAL(KIND=DP) :: PCCVA
      REAL(KIND=DP) :: WHBTA
      REAL(KIND=DP) :: RHORA
      REAL(KIND=DP) :: TDRA
      REAL(KIND=DP) :: ZRAA
      REAL(KIND=DP) :: SUFA
      REAL(KIND=DP) :: SUFO
      REAL(KIND=DP) :: PCRA
      REAL(KIND=DP) :: TCRA
      REAL(KIND=DP) :: PAFA
      REAL(KIND=DP) :: TTPA
      REAL(KIND=DP) :: PTPA
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: VISCA
      INTEGER :: IVISO
      INTEGER :: IVISS
      INTEGER :: IVAPO
      INTEGER :: IRHOO
      INTEGER :: IMTC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPCSLS
      INTEGER :: IVAPS
      INTEGER :: IRHOS
      INTEGER :: NQA
      INTEGER :: NQO
      INTEGER :: IVISA
      INTEGER :: ITERDC
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
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PRTM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SF
      CHARACTER(64) :: UNTM
      CHARACTER(64) :: UNLN
      CHARACTER(64) :: UNAR
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: UNPLOT
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: UNREF
      CHARACTER(64), DIMENSION(:,:), ALLOCATABLE :: UNSF
      CHARACTER(6), DIMENSION(:), ALLOCATABLE :: CHREF
      CHARACTER(6), DIMENSION(:,:), ALLOCATABLE :: CHSF
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPLOT
      INTEGER, DIMENSION(:), ALLOCATABLE :: IREF
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPLOTGC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IREFGC
      INTEGER, DIMENSION(:), ALLOCATABLE :: NDREF
      INTEGER :: NPRTM
      INTEGER :: NVPLOT
      INTEGER :: NREF
      INTEGER :: NVREF
      INTEGER :: ICNO
      INTEGER :: ICNS
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFT
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFGC
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFF
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFD
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISFC
      INTEGER :: NSF
      INTEGER :: NSFGP
      INTEGER :: IHSF
      INTEGER :: IFQS
      INTEGER :: IFQO
      INTEGER :: ISGNS
      INTEGER :: ISGNO
      INTEGER :: ISGNP
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFGP
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISFSN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSFDOM
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ISFDOM
!
!---  End of module  ---
!
      END MODULE



!---------------------Fortran 90 Module--------------------------------!
!
      MODULE POINTE
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
!     Numerical scheme pointer variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER, DIMENSION(:), ALLOCATABLE :: MNOD
      INTEGER, DIMENSION(:), ALLOCATABLE :: MADJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: MFLX
      INTEGER, DIMENSION(:), ALLOCATABLE :: MPOS
      INTEGER, DIMENSION(:), ALLOCATABLE :: MNEG
      INTEGER, DIMENSION(:), ALLOCATABLE :: MPOSB
      INTEGER, DIMENSION(:), ALLOCATABLE :: MNEGB
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVP_FRC
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
!     Fracture field variables
!
!     PERM_FRC(M,NT) permeability of fracture triangle
!     PL_FRC(M,NT) aqueous pressure of fracture triangle
!     PG_FRC(M,NT) gas pressure of fracture triangle
!     RKG_FRC(M,NT) gas relative permeability of fracture triangle
!     RKL_FRC(M,NT) aqueous relative permeability of fracture triangle
!     SG_FRC(M,NT) gas saturation of fracture triangle
!     SL_FRC(M,NT) aqueous saturation of fracture triangle
!     T_FRC(M,NT) temperature of fracture triangle
!     XLA_FRC(M,NT) aqueous air mass fraction of fracture triangle
!     XLW_FRC(M,NT) aqueous water mass fraction of fracture triangle
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 August 2017.
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: T_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PG_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PN_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DNR_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PSO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SG_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SN_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PORD_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PORT_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: POR0_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PSW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOG_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHON_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLN_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLS_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SI_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PI_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOMG_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOML_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOMN_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TMS_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: POSM_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PERMRF_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PCMP_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TCMP_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: APH_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: APM_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PERM_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: VOL_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ASL_FRC,AST_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ASLMIN_FRC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NPHAZ_FRC
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVT_FRC
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
!     Energy field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKG_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKN_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HG_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HN_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLS_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEG_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEN_FRC
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVG_FRC
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
!     Gas field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: GNIFT_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISG_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORG_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKG_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGA_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICAIR_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICCO2_FRC
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVN_FRC
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
!     NAPL field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISN_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORN_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKN_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XNA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XNO_FRC
!
!---  End of module  ---
!
      END MODULE


!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVS_FRC
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
!     Salt/surfactant field variables
!
!     XLS_FRC(M,NT) aqueous salt mass fraction of fracture triangle
!     YLS_FRC(M,NT) aqueous total salt mass fraction of frac. triangle
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLS_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLS_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YLS_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOSP_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HSP_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SS_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICBRN_FRC
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVGC_FRC
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
!     Gas component field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 November 2007.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BETA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: FK_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XGC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMGC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XNC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMNC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XLC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMLC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFGC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFNC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFLC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: HGC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMVGC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMVLC_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMVGW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMVLW_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: TMC_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZG_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ZMC_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZN_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IZMC_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBETA_FRC
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FLUX_FRC
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
!     Fracture flux variables
!
!     TRNSA_FRC(M,NC) air mass transfer rate from fracture triangle to
!       matrix grid cell, kg/s
!     TRNSS_FRC(M,NC) salt mass transfer rate from fracture triangle to
!       matrix grid cell, kg/s
!     TRNSQ_FRC(M,NC) energy transfer rate from fracture triangle to
!       matrix grid cell, W
!     TRNSW_FRC(M,NC) water mass transfer rate from fracture triangle to
!       matrix grid cell, kg/s
!     UFFC(NT,NSOLU) solute flux between fracture triangles 1/m^2 s
!     UFMC(NC,NSOLU) solute flux from fracture triangle
!       to grid cell, 1/m^2 s
!     UFFDGW(M,NT) diffusive water vapor flux through gas, kmol/m^2 s
!     UFFDLA(M,NT) diffusive air flux through aqu., kmol/m^2 s
!     UFFDLW(M,NT) diffusive water flux through aqu., kmol/m^2 s
!     UFFDS(M,NT) diffusive salt flux through aqu., kg/m^2 s
!     UFFS(M,NT) salt flux through aqu., kg/m^2 s
!     UFFG(M,NT) gas fracture to fracture flux, m/s
!     UFFGW(M,NT) water vapor flux through gas, kg/m^2 s
!     UFFLA(M,NT) air/CO2 flux through aqueous, kg/m^2 s
!     UFFLW(M,NT) water flux through aqueous, kg/m^2 s
!     UFFL(M,NT) aqueous fracture to fracture flux, m/s
!     UFFQ(M,NT) heat fracture to fracture flux, W/m^2
!     UFML(NC) aqueous volum. flow rate from fracture to matrix, m^3/s
!     UFMG(NC) gas volum. flow rate fracture to matrix, m^3/s
!     UFMN(NC) NAPL volum. flow rate fracture to matrix, m^3/s
!     UFMT(NC) heat flow from fracture to matrix, W
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 September 2017.
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TRNSA_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TRNSS_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TRNSQ_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TRNSW_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFMC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFDGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFDLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFDLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFDS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UFFQ
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: UFMG
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: UFML
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: UFMN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: UFMT
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FLUXC_FRC
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
!     Gas-component fracture flux variables
!
!     TRNSGC_FRC component mass transfer rate from fracture/fault 
!       triangle to matrix grid cell, kg/s
!     UFFDGC diffusive component flux through gas, kmol/m^2 s
!     UFFGC component flux through gas, kg/m^2 s
!     UFFDNC diffusive component flux through NAPL, kmol/m^2 s
!     UFFNC component flux through NAPL, kg/m^2 s
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 April 2020.
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: TRNSC_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: UFFDGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: UFFDNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: UFFGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: UFFNC
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE GEOM_FRC
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
!     Fracture geometry variables
!
!     AFN_FRC(NFNC) surface area of fracture triangle to node connection
!     AFF_FRC(NFFC) surface area of frac. tri. to frac. tri. connection
!     APM_FRC(M,NT) mechanical aperture of fracture triangle
!     APH_FRC(M,NT) hydraulic aperture of fracture triangle
!     AF_FRC(NT) surface area of fracture triangle
!     DFN_FRC(NFNC) distance of fracture triangle to node connection
!     DFF_FRC(NFFC) distance of frac. tri. to frac. tri. connection
!     DFFM_FRC(NFFC) dist. of frac. tri. to frac. tri. connection surf.
!     SIGN_FRC x,y,z components of normal stress on fracture triangle
!       1 - x direction, Pa
!       2 - y direction, Pa
!       3 - z direction, Pa
!     SIGS_FRC x,y,z components of shear stress on fracture triangle
!       1 - x direction, Pa
!       2 - y direction, Pa
!       3 - z direction, Pa
!     SFNC_FRC(3,LFC_FRC) - x,y,z components of unit surface normal
!     SFNT_FRC(3,LT_FRC) - x,y,z components of unit surface normal
!     XE_FRC x dimension of fracture triangle vertex
!     YE_FRC y dimension of fracture triangle vertex
!     ZE_FRC z dimension of fracture triangle vertex
!     XP_FRC x dimension of fracture triangle centroid
!     YP_FRC y dimension of fracture triangle centroid
!     ZP_FRC z dimension of fracture triangle centroid
!     XPFN_FRC x dimension of fracture triangle to node centroid
!     YPFN_FRC y dimension of fracture triangle to node centroid
!     ZPFN_FRC z dimension of fracture triangle to node centroid
!     ND_FRC resident grid cell number of fracture triangle vertex
!     NM_FRC fracture name
!     INCM_FRC(NT) fracture triangle to grid cell connection map
!     ITCM_FRC(NT) fracture triangle to fracture triangle connection map
!     IBHT_FRC(NT) fracture triangle of borehole-fracture connection
!     IBHN_FRC(NT) borehole node of borehole-fracture connection
!     IP_FRC(2,NF) fracture triangle pointer
!       1 - low index
!       2 - high index
!     IPN_FRC(2,NT) fracture triangle to grid cell pointer
!       1 - low index
!       2 - high index
!     IPF_FRC(2,NT) fracture triangle to fracture triangle pointer
!       1 - low index
!       2 - high index
!     IXP_FRC active/inactive fracture/fault index
!     IZ_FRC fracture/fault triangle rock type
!     NDFGT_MPI(NC,NP) global of fracture/ghost triangles on processor
!     NDFT_MPI(NC,NP) global of fracture triangles on processor
!     NDGT_MPI(NC,NP) global of ghost triangles on processor
!     NCFGT_MPI(NP) count of fracture/ghost triangles on processor
!     NCFT_MPI(NP) count of fracture triangles on processor
!     NCGT_MPI(NP) count of ghost triangles on processor
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 August 2017.
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XE_FRC,YE_FRC,ZE_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XP_FRC,YP_FRC,ZP_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SIGN_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SIGS_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SFNC_FRC,SFNT_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: AF_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: AFN_FRC,AFF_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DFN_FRC,DFF_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DFNM_FRC,DFFM_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XPFN_FRC,YPFN_FRC,
     &  ZPFN_FRC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ND_FRC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NDFGT_MPI,NDFT_MPI
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IP_FRC,IPN_FRC,IPF_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: INCM_FRC,ITCM_FRC,NTP_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBHT_FRC,IBHN_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IXP_FRC,IZ_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: NCFGT_MPI,NCFT_MPI,NCGT_MPI
      INTEGER :: NF_FRC,NT_FRC,NXP_FRC
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: NM_FRC
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE PARM_FRC
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
!     Fracture parameter variables
!
!     VIC_FRC(NVIC_FRC,NF) fracture initial condition variables
!     NVIC_FRC number of fracture initial condition variables
!     SKF_FRC(NF) skin factor for fracture-fracture intersections
!     IC_FRC(NF) fracture initial condition state index
!       1 - hydrostatic
!       2 - aqueous saturation, gas pressure
!       3 - aqueous saturation, aqueous pressure
!       4 - aqueous pressure, gas pressure
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VIC_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SKF_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKSP_FRC,DSBB_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: JRC_FRC,JCS_FRC,
     &  UCS_FRC,THCP_FRC,DTHCP_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCA_FRC,SRCS_FRC,
     &  SRCT_FRC,SRCW_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: SRCGC_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SRCG_FRC,SRCL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RSDL_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RSD_FRC,RSDAVG_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: SRC_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IC_FRC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IC_OPT_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSD_FRC,NDREF_FRC,IJM_FRC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IM_FRC,KLU_FRC
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISRM_FRC,ISRT_FRC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISRDM_FRC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLUC_FRC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU_MCF,KLU_FCM,KLU_FCB
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLUC_MCF,KLUC_FCM,KLUC_FCB
      INTEGER :: NREF_FRC,NSR_FRC
      INTEGER :: NVIC_FRC = 21
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE TRNS_FRC
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
!     Transport variables
!
!     C_FRC(NT,NSOLU) volumetric solute or total species concentration 
!       of fracture triangle
!     CO_FRC(NT,NSOLU) volumetric solute or total species concentration 
!       of fracture triangle at old time step
!     CRNTG_FRC(NT) gas Courant number for fracture
!     CRNTL_FRC(NT) aqueous Courant number for fracture
!     CRNTN_FRC(NT) nonaqueous-liquid Courant number for fracture
!     ICT_FRC(NF,NSOLU) initial condition index of solute or total
!       species of fracture triangle
!     IC_SP_FRC(NF,NSP) initial condition index of reactive species
!       of fracture triangle
!     PCSL_FRC(5,NT,NSOLU) solid-aqueous partition coefficient parameters
!       of fracture triangle
!     SDCL_FRC(3,NT,NSOLU) solute aqueous diffusion model parameters
!       of fracture triangle
!     SP_C_FRC(NT,LSPR) - fracture species concentration, 
!       mol/m^3 node volume
!     SP_CO_FRC(NT,LSPR) - old fracture species concentration, 
!       mol/m^3 node volume
!     RHOS_FRC(NT) rock grain density of fracture triangle
!     YG_FRC(NT,NSOLU) gas mole fraction of solute or total species
!       of fracture triangle
!     YL_FRC(NT,NSOLU) aqueous mole fraction of solute or total species
!       of fracture triangle
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 June 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: C_FRC,CO_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_C_FRC,SP_CO_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RHOS_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: PCSL_FRC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: SDCL_FRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YG_FRC,YL_FRC,YN_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CRNTL_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CRNTG_FRC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CRNTN_FRC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICT_FRC,IC_SP_FRC
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE WELL_FD
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
!     Well field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHONW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLOW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLOW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLWW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGWW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLAW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGAW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLAW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGWW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: STW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YMLOW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLAW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGWW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLOW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PWLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PWGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PWNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOLWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOGWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHONWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLWWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGWWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLAWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGAWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLAWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGWWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLOWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLOWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISLWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISGWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISNWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLAWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGWWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLOWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCW_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCO_W
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RHONW_O
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SNW_O
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XLOW_O
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: STW_O
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RHOLW_O
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RHOGW_O
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SLW_O
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XLWW_O
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XGWW_O
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XGAW_O
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XLAW_O
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE WELL_FX
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
!     Well flux variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WL_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WG_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WN_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UL_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UG_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UN_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDLA_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDGW_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WDLO_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDLA_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDGW_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UDLO_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QL_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QN_W
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QT_W
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE WELL_CL
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
!     Well control variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WBR
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WBRS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WWD
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WHP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WRH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WIDA
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: WLVR
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WWDL
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WWDN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WIDO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DNRW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWCC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IXW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWL
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NPHAZW
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSZW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWN
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWT
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IWLDM
      INTEGER :: NWLS
      INTEGER :: NWLN
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE PORMED
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
!     Porous media property variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: POR
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TOR
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CMP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RHOS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CPS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PERM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SCHR
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CHML
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: GAMMA
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SCALNM
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RPLT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFEF
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RPGC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RPLC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RPNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: CPLY_SL
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: CPLY_RL
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DPLGA,DPLLA
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DPTGA,DPTLA
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: HCMWE
      INTEGER, DIMENSION(:), ALLOCATABLE :: ITOR
      INTEGER, DIMENSION(:), ALLOCATABLE :: IRPG
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAMMA
      INTEGER :: NSCALE
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISCALE
      INTEGER, DIMENSION(:), ALLOCATABLE :: IRPL
      INTEGER, DIMENSION(:), ALLOCATABLE :: IRPN
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISCHR
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISM
      INTEGER, DIMENSION(:), ALLOCATABLE :: ITHK
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IRPLT
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISLTBL
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IRLTBL
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IRGTBL
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IRLTBLT
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IRNTBL
      INTEGER, DIMENSION(:), ALLOCATABLE :: IDP
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISKP
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPRF
      INTEGER, DIMENSION(:), ALLOCATABLE :: INHYD
      INTEGER, DIMENSION(:), ALLOCATABLE :: NPLY_SL
      INTEGER, DIMENSION(:), ALLOCATABLE :: NPLY_RL
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
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP) :: TM
      REAL(KIND=DP) :: TMMX
      REAL(KIND=DP) :: TMPR
      REAL(KIND=DP) :: DT
      REAL(KIND=DP) :: DTI
      REAL(KIND=DP) :: DTMX
      REAL(KIND=DP) :: DTMN
      REAL(KIND=DP) :: DTAF
      REAL(KIND=DP) :: DTCF
      REAL(KIND=DP) :: DTO
      REAL(KIND=DP) :: DTSO
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TMPS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TMPE
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TMPD
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TMPX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TMPN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TMPA
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TMPC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RSD,RSDAVG
      REAL(KIND=DP) :: RSDMX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RSDM
      REAL(KIND=DP) :: RLXF
      REAL(KIND=DP) :: RLMSG
      REAL(KIND=DP) :: CPUMX
      REAL(KIND=DP) :: CLKMX
      REAL(KIND=DP) :: CPUSEC
      REAL(KIND=DP) :: CLKSEC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: WFMN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: VSLC
      REAL(KIND=DP) :: CRNTMXC
      REAL(KIND=DP) :: RTOL_PETSC,ATOL_PETSC,DTOL_PETSC
      INTEGER :: IVRSN
      INTEGER :: ISIC
      INTEGER :: ICNV
      INTEGER :: IEO
      INTEGER :: ILES
      INTEGER :: IOM
      INTEGER :: ICODE
      INTEGER :: IEQT
      INTEGER :: IEQW
      INTEGER :: IEQA
      INTEGER :: IEQN
      INTEGER :: IEQO
      INTEGER :: IEQC
      INTEGER :: IEQS
      INTEGER :: IEQD
      INTEGER :: IEQDO
      INTEGER :: IEQHA
      INTEGER :: IEQHN
      INTEGER :: IEQHO
      INTEGER :: IEQALC
      INTEGER :: IEQDA
      INTEGER :: IAQU
      INTEGER :: IGAS
      INTEGER :: INAPL
      INTEGER :: JSTOP
      INTEGER, DIMENSION(:), ALLOCATABLE :: IEQGC
      INTEGER, DIMENSION(1:100) :: ISLC = 0
      INTEGER, DIMENSION(1:20) :: IDMN = 1
      INTEGER :: NEPD
      INTEGER :: MEPD
      INTEGER :: IEPD
      INTEGER :: ICSN
      INTEGER :: IMSG
      INTEGER :: ICD
      INTEGER :: IVR
      INTEGER :: NRIMX
      INTEGER, DIMENSION(:), ALLOCATABLE :: ITOFF,NRIM
      INTEGER :: NSTEP
      INTEGER :: NRST
      INTEGER :: NITER
      INTEGER :: NTSR
      INTEGER :: NGC
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSD
      INTEGER :: MXSTEP
      INTEGER :: IDFLT
      INTEGER :: IDFLTD
      INTEGER :: IUNM
      INTEGER :: IUNKG
      INTEGER :: IUNS
      INTEGER :: IUNK
      INTEGER :: IUNMOL
      INTEGER :: ISUB_LOG
      INTEGER :: IPID
      INTEGER :: MAXITS_PETSC
      CHARACTER(132) :: USER
      CHARACTER(132) :: CMPNY
      CHARACTER(132) :: TITLE
      CHARACTER(132) :: INPDAT
      CHARACTER(132) :: INPTIM
      CHARACTER(64) :: CARD
      CHARACTER(132) :: VARB
      CHARACTER(8) :: CHDATE
      CHARACTER(64) :: CH_VRSN
      CHARACTER(10) :: CHTIME
      CHARACTER(132), DIMENSION(:), ALLOCATABLE :: NOTES
      CHARACTER(132) :: SUBNM
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: GCNM
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: HCNM
      CHARACTER(32), DIMENSION(:), ALLOCATABLE :: SUB_LOG
      CHARACTER(256) :: CHMSG
      CHARACTER(LEN=8) :: PID_CHAR=' '
      CHARACTER(132) :: PIDFNM = ' '
      LOGICAL :: IFXST
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
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: SRC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCD
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: SRCGC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SRCIT
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SRCIW
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SRCIA
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SRCIO
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SRCIS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SRCID
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCIGC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PLWB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QTW
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: GWSI
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: SWSI
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IWSI
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISRDM
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISRM
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISRT
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSOLSR
      INTEGER :: NSR
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE SPILL
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
!     Spill data variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 December 2007.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UNSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VNSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ULSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VLSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHONSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOLSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISNSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISLSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HNSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCOSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCWSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DNRSP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SPNORM
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SPHMIN
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
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TBLX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TBLY
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TBLDDX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TBLDDY
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
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: C
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CNL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CNLO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SMDEF
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: WC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: HLF
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHLF
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCIC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SMDL
      REAL(KIND=DP) :: SMDLS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SMDN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SMDG
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: SDCL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SDCLS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YN
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: PCSL
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: PCSN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PCSLD
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CMTLN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PCLN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PCGL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PCGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CBO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YLB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YGB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YNB
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DISPL
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DISPT
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RCHDF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SOLML
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RCHDFL
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RCHDFN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHLFL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHLFN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ELC_DCF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ELC_VCF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ELC_SCF
      REAL(KIND=DP) :: ELC_DUN
      REAL(KIND=DP) :: ELC_VUN
      REAL(KIND=DP) :: DT_CRN
      REAL(KIND=DP) :: DTI_CRN
      REAL(KIND=DP) :: TM_CRN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CCL_CRN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DPLGS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DPTRS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DPLD
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DPTD
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: D50
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PMDD
      CHARACTER(64) :: ELC_SOL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CHDF
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SOLUT
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CRNTL
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CRNTG
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CRNTN
      REAL(KIND=DP) :: CRNTMXT
      INTEGER :: NSOLU
      INTEGER :: IDISP
      INTEGER, DIMENSION(:), ALLOCATABLE :: IEDL
      INTEGER :: IEDLS
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPCL
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPCN
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IPCSL
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPCSLD
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICT
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICTN
      INTEGER :: IDSPS
      INTEGER :: IDSPD
      INTEGER :: ICRNT
      INTEGER, DIMENSION(:), ALLOCATABLE :: NCHEM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPCLN
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPCGL
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPCGN
      INTEGER, DIMENSION(:), ALLOCATABLE :: IMTLN
      INTEGER :: NSL_ELC
      INTEGER :: IDF_ELC
      INTEGER :: IVF_ELC
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
      MODULE DUAL_POR
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
!     Dual porosity model variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 November 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: FRAC_P
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: VOL_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DNR_M
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFGC_M
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFNC_M
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: TMC_M
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XGC_M
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMGC_M
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMNC_M
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XNC_M
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ZMC_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGW_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLA_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLS_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HG_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGA_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGW_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HL_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLW_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HN_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HSP_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEG_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEL_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEN_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKG_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKL_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKN_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PERMRF_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: POR0_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PG_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PL_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PN_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: POSM_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PSO_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PORD_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PORT_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVA_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVW_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOG_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOL_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOMG_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOML_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOMN_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHON_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOSP_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKG_M
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RKL_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKN_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SG_M,SG_F
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SL_M,SL_F
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SN_M,SN_F
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SS_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: T_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TMS_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORG_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORL_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORN_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISG_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISL_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISN_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGW_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLA_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLS_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLW_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGW_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLA_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLS_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLW_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YLS_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZG_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZN_M
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QLA_FM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QGW_FM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QLW_FM
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: QGC_FM
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: QNC_FM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QS_FM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: QQ_FM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NPHAZ_M
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IM_M
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
!     Primary field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: T
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DNR
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PSO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PORD
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PORT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PSW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHON
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORL
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RKL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SI
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PI
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOMG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOML
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOMN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SDPM
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SDPF
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TMS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: POSM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BTGL,BTGN,BTNL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PERMRF
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: POR0
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PCMP,TCMP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PERMV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SCHRV
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVT
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
!     Energy field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGD
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HLW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEN
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVG
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
!     Gas field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: GNIFT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISDG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TMBP_A
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TMBP_N
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TMBP_O
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICAIR
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICCO2
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVN
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
!     NAPL field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ANLT
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ANLF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: HKL
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: HKNT
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: HKNF
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XSO
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVI
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
!     Ice field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOI
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKI
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HI
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVS
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
!     Salt/surfactant field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: OEC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YLS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HSP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SS
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICBRN
      INTEGER :: ISALT
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVD
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
!     Dissolved-oil field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TRPNL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TRPGL
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: AVPVG
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: AVPVL
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVA
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
!     Alcohol field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XNA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XNN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XNO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMNA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMNN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMNO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMNW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SUGL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SUGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SUNL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFNA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFNN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFNO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFNW
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVGC
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
!     Gas component field variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 November 2007.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BETA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: FK
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: PVC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XLC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMLC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFNC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: DFLC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: HGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMVGC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: XMVLC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMVGW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMVLW
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: TMC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZG
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ZMC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZN
      INTEGER, DIMENSION(:), ALLOCATABLE :: IZMC,IBETA
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICCOMP
      INTEGER :: IFK
      
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE HYDT
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
!---  Ternary hydrate variables  ---
!
!     ZPC_HT - mole fraction of hydrate formers for phase equilibria
!     ZHC_HT - mole fraction of hydrate formers for hydrate equilibria
!     DZPC_HT - mole fraction spacing for phase equilibria 
!     DZHC_HT - mole fraction spacing for hydrate equilibria 
!     TCR_HT - critical temperature, K
!     PCR_HT - critical pressure, Pa
!     TCB_HT - cricondenbar temperature, K
!     PCB_HT - cricondenbar pressure, Pa
!     TCT_HT - cricondenterm temperature, K
!     PCT_HT - cricondenterm pressure, Pa
!     TLE_HT - temperature for lower phase envelope, K
!     PLE_HT - pressure for lower phase envelope, Pa
!     D2PLE_HT - second derivative of pressure for lower phase envelope
!     FKLE_HT - K-factor for lower phase envelope, Pa
!     D2FKLE_HT - second derivative of K-factor for lower phase envelope
!     TUE_HT - temperature for upper phase envelope, K
!     PUE_HT - pressure for upper phase envelope, Pa
!     D2PUE_HT - second derivative of pressure for upper phase envelope
!     FKUE_HT - K-factor for upper phase envelope, Pa
!     D2FKUE_HT - second derivative of K-factor for upper phase envelope
!     T2P_HT - temperature for two-phase region, K
!     P2P_HT - pressure for two-phase region, Pa
!     B2P_HT - beta for two-phase region
!     D2B2P_HT - second derivative of beta for two-phase region
!     FK2P_HT - K-factor for two-phase region
!     D2FK2P_HT - second derivative of K-factor for two-phase region
!     THE_HT - temperature for hydrate equilibrium table, K
!     PHE_HT - pressure for hydrate equilibrium table, Pa
!     XSCA_HT - CO2 small cage occupany for hydrate equilibrium table
!     XSCO_HT - CH4 small cage occupany for hydrate equilibrium table
!     XSCN_HT - N2 small cage occupany for hydrate equilibrium table
!     XLCA_HT - CO2 large cage occupany for hydrate equilibrium table
!     XLCO_HT - CH4 large cage occupany for hydrate equilibrium table
!     XLCN_HT - N2 large cage occupany for hydrate equilibrium table
!     THE2P_HT - second derivative of temperature for hydrate equil.
!     XSCA2T_HT - second derivative of CO2 small cage occupany for hyd.
!     XSCO2T_HT - second derivative of CH4 small cage occupany for hyd.
!     XSCN2T_HT - second derivative of N2 small cage occupany for hyd.
!     XLCA2T_HT - second derivative of CO2 large cage occupany for hyd.
!     XLCO2T_HT - second derivative of CH4 large cage occupany for hyd.
!     XLCN2T_HT - second derivative of N2 large cage occupany for hyd.
!     PHE2T_HT - second derivative of pressure for hydrate equilibrium
!     NPEP_HT - number of phase envelope table points
!     NHEP_HT - number of hydrate equilibrium table points
!     NTP_HT - number of two-phase temperature points
!     NPP_HT - number of two-phase presure points
!     NHF_HT - number of hydrate formers
!     NZP_HT - dimension of phase envelope table
!     NZH_HT - dimension of hydrate table
!     IZP_HT - bounding concentrations indices for phase equilibria
!     IZH_HT - bounding concentrations indices for hydrate equilibrium
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 July 2012.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZPC_HT,ZHC_HT
      REAL(KIND=DP) :: DZPC_HT,DZHC_HT
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TCR_HT,PCR_HT
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TCB_HT,PCB_HT
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TCT_HT,PCT_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TLE_HT,PLE_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: D2PLE_HT
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: FKLE_HT,D2FKLE_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TUE_HT,PUE_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: D2PUE_HT
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: FKUE_HT,D2FKUE_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: T2P_HT,P2P_HT
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: B2P_HT,D2B2P_HT
      REAL(KIND=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FK2P_HT
      REAL(KIND=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: D2FK2P_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THE_HT,PHE_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XSCA_HT,XSCO_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XSCN_HT,XLCA_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLCO_HT,XLCN_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THE2P_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XSCA2T_HT,XSCO2T_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XSCN2T_HT,XLCA2T_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLCO2T_HT,XLCN2T_HT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PHE2T_HT
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ZMIH_HT
      INTEGER, DIMENSION(:), ALLOCATABLE :: NPEP_HT,NHEP_HT
      INTEGER, DIMENSION(:), ALLOCATABLE :: NTP_HT,NPP_HT
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IZP_HT,IZH_HT
      INTEGER :: NHF_HT,NZH_HT,NZP_HT
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
!     Written by MD White, PNNL, 4 September 2004.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HCPP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XHW
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XHA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XHO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XHN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YMGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YMHGA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YMGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YMHGN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YMGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YMHGO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVHA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVHN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVHO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TMHA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TMHO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TMHN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZMCA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZMCO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZMCN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PPEL,PPEU,TCR,TCT
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CHKN
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IC_OPT
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE EOR
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
!     Enhanced oil recovery variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 March 2013.
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BIPC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BIPF
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PFPP
      REAL(KIND=DP), DIMENSION(15) :: PR0LKR,PR1LKR
      REAL(KIND=DP), DIMENSION(15,40) :: H0LKR,H1LKR
      REAL(KIND=DP), DIMENSION(40) :: TR0LKR,TR1LKR
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMPCF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: GC_MMP
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: PFNM
      INTEGER, DIMENSION(15) :: IL0LKR,IL1LKR
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IPCF
      INTEGER, DIMENSION(:), ALLOCATABLE :: NPCF
      INTEGER :: NPF
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE UCODE
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
!     Inverse (UCode) variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: R_OBDS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: R_OBDT
      REAL(KIND=DP) :: TMOB
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: C_OBDT
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: I_OBDT
      INTEGER :: IOBDEF
      INTEGER :: IOBDSF
      INTEGER :: IOBDUF
      INTEGER :: NOBDT
      INTEGER, DIMENSION(:), ALLOCATABLE :: NOBDS
      INTEGER :: NOBDP
      LOGICAL :: FLG_UNI
      LOGICAL :: FLG_EXT
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
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: P_TA
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: T_TA
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RHO_TA
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: H_TA
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: U_TA
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: FUG_TA
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: S_TA
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RHO_ST
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: H_ST
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: U_ST
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: FUG_ST
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: S_ST
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: P_PH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: T_TH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YMHO_PH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YMHO_TH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XSCA_PH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XSCA_TH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLCA_PH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLCA_TH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XSCO_PH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XSCO_TH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLCO_PH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLCO_TH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: T_PH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: P_TH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: YMGO_PH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: YMGO_TH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: T_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: P_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOL_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HL_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UL_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOV_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HV_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UV_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: FUG_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SL_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SV_LV
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CINH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XLIMX
      INTEGER :: IT_PH
      INTEGER :: IO_PH
      INTEGER :: IT_TH
      INTEGER :: IS_TH
      INTEGER :: IO_TH
      INTEGER :: I_INH
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPTP,IPCR
      INTEGER, DIMENSION(:), ALLOCATABLE :: IP_TA
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IT_TA
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IV_TA
      INTEGER :: INCG
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IC_NCG
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_LV
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: INHNM
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE PLT_ATM
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
!     Plant and atmospheric variables
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 5 November 2002.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ALBEDO
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ATMC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ATMOS
      REAL(KIND=DP) :: ATMST
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFGW_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DFLA_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DNR_SFC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RW_GS,RE_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RE_PL,RW_CP,RE_CP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HG_GS,HGA_GS,HGW_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HL_GS,HLW_GS,HLW_PL
      INTEGER, DIMENSION(:), ALLOCATABLE :: IALB
      INTEGER, DIMENSION(:), ALLOCATABLE :: IALB_P
      INTEGER :: IATM_C
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICM_SFC
      INTEGER, DIMENSION(:), ALLOCATABLE :: ID_SFC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IRSM_P
      INTEGER :: ISCVF
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISRM_P
      INTEGER :: ISWBCF
      INTEGER, DIMENSION(:), ALLOCATABLE :: ITMP_P
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: JM_SFC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU_SFC
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: NAME_SFC
      INTEGER :: NATM_T
      INTEGER :: NPLANT
      INTEGER :: NSFCA
      INTEGER :: NSFCC
      INTEGER :: NSFCN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSFCT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PARMS_P
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: PARMS_SFC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PL_GS
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: PLANT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVA_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVW_CP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PVW_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOG_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOL_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOMG_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RHOML_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ACW_PL,WM_PL,DM_PL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKG_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RKL_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SWM_PL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RWU_SFC,TP_SFC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SG_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SL_GS
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PHMX
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: T_CP
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: T_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: T_PL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKG_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: THKL_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORG_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TORL_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UEG_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISG_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VISL_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGA_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XGW_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLA_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XLW_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGA_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMGW_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLA_GS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XMLW_GS
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
!     EQ_E(LSEE,LEQE) - equilibrium equation parameters
!       species coefficients, equilibrium equations
!     EQ_K(LSEK+LREK,LEQK) - kinetic equation parameters
!       species coefficients + reaction coefficients, kinetic equations
!     RC_E(5,LRCE) - equilibrium reaction parameters
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
!     NEQE - number of equilibrium equations
!     NEQK - number of kinetic equations
!     NRCE - number of equilibrium reactions
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
!     RCNME(LRCE) - equilibrium reaction name
!     RCNMK(LRCK) - kinetic reaction name
!     IEQ_C(LSEC+1,LEQC) - integer conservation equation parameters
!       number of species + species number + ... +
!     IEQ_E(LSEE+1,LEQE) - integer equilibrium equation parameters
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
!     ISP_OW(LSPS,LFDR) - solid species lithology card overwrite
!     LEQC - number of conservation equations 
!     LEQE - number of equilibrium equations
!     LEQK - number of kinetic equations 
!     LRCE - number of equilibrium reactions
!     LRCK - number of kinetic reactions 
!     LREK - number of kinetic equation reactions 
!     LSEC - number of conservation equation species
!     LSEE - number of equilibrium equation species
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
!     Written by MD White, PNNL, 8 December 2004.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP) :: ACTVC,DT_RST,DTI_RST
      LOGICAL :: ECKE_ER
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ACTV16
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: C_PH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CHARG
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: FACTV,ACTVS
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: EQ_C
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: EQ_E
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: EQ_K
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
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CFMX
      INTEGER :: N_RST,NEQC,NEQE,NEQK,NESITE
      INTEGER :: NRCE,NRCK,NRTSI,NSPC,NSPE
      INTEGER :: NSPG,NSPK,NSPL,NSPLK,NSPN
      INTEGER :: NSPR,NSPS
      REAL(KIND=DP) :: CMIN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RC_E
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RC_K
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: RCNME
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: RCNMK
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: RS_S
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_C
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_CBO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_CMN
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_RATE
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_AREA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_CO
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_CI
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_L
      REAL(KIND=DP) :: SP_MDG,SP_MDL,SP_MDN
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: SP_SDCL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_S
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNMC
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNME
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNMG
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNMK
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNML
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNMN
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: SPNMS
      REAL(KIND=DP) :: TM_RST
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YSPG
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YSPL
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YSPN
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
!       number and node number
!       1 i,j,k
!       2 i+1,j,k
!       3 i,j+1,k
!       4 i+1,j+1,k
!       5 i,j,k+1
!       6 i+1,j,k+1
!       7 i,j+1,k+1
!       8 i+1,j+1,k+1  ---
!     NE_GM(8,LFEN) node number at hexaderal position number and finite
!       element number
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
!       8
!       9
!     IXP_GM(LFDM) inactive elements index 
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 17 March 2017 (Saint Patrick's Day).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: BC_GM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PROP_GM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: U_GM,V_GM,W_GM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: EPS_GM,SIG_GM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: P_GM,EPSV_GM,SIGV_GM
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: EPSV_CMP,SIGV_CMP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RSDM_GM
      REAL(KIND=DP) :: GRAV_GM,RSD_GM
      INTEGER :: IREF_GM,NNZR_GM
      INTEGER :: NFEN_GM,MD_GM,MK_GM,ML_GM,MU_GM,NBC_GM,NSD_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCC_GM,IBCM_GM,IBCN_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBCIN_GM,IBCD_GM,IBCR_GM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IBCT_GM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IPROP_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IM_GM,IHCM_GM,IXP_GM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ND_GM,NE_GM,NK_GM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: MLU_GM
      INTEGER, DIMENSION(:), ALLOCATABLE :: NLU_GM
      INTEGER, DIMENSION(1:2) :: K_GM = 0
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ID_MPI_GM,JD_MPI_GM,
     &  KD_MPI_GM
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
!     Written by MD White, PNNL, 20 June 2005.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: B0
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: B1
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: B2
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CMXX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TCC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: TAA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PSIC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PSIA
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ALAMB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: CLAMB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ELAMB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: HOLAMB
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BPPR
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BPHI
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BPR
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BMMX
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ATB0
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ATB1
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ATB2
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ATCMX
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ATNLAM
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ATCLAM
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ATALAM
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ATHLAM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ATTC
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ATPC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ATTA
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ATPA
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CTCPH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CTC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CTCPR
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CTCPPR
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CTAPH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CTA
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CTAPR
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CTAPPR
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ETH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ETHP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ETHP2
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
      MODULE SACIN
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
!     SAC number of events
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 June 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: INMBLT
      INTEGER :: INMREM
      INTEGER :: INMREL
      INTEGER :: KEXIT
      INTEGER, DIMENSION(:), ALLOCATABLE :: IRLK
      INTEGER, DIMENSION(:), ALLOCATABLE :: IRMSOL
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE SACTM
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
!     SAC times of events
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 June 2005.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP) :: TMBL
      REAL(KIND=DP) :: TMRM
      REAL(KIND=DP) :: TMRL
!
!---  End of module  ---
!
      END MODULE
      
!---------------------Fortran 90 Module--------------------------------!
!
      MODULE SACEV
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
!     SAC arrays of events
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 June 2005.
!     allo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DBLTIM
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DRMTIM
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DRMDEP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DRMEFF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DRLTIM
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CRMSIT
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FDVP_BH
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
!     Borehole field variables
!
!     DNR_BH(LUK,LBN_BH) - primary variable increments of borehole node
!     DFGW_BH(LUK,LBN_BH) - gas water diffusion coefficient of 
!       borehole node
!     DFLA_BH(LUK,LBN_BH) - aqueous air diffusion coefficient of 
!       borehole node
!     DFLS_BH(LUK,LBN_BH) - aqueous salt diffusion coefficient of 
!       borehole node
!     HGW_BH(LUK,LBN_BH) - water vapor enthalpy of borehole node
!     HGA_BH(LUK,LBN_BH) - air enthalpy of borehole node
!     HG_BH(LUK,LBN_BH) - gas enthalpy of borehole node
!     HL_BH(LUK,LBN_BH) - aqueous enthalpy of borehole node
!     HSP_BH(LUK,LBN_BH) - precipitated salt enthalpy of borehole node
!     PERMRF_BH(LUK,LBN_BH) - permeability reduction facture of 
!       borehole node
!     PL_BH(LUK,LBN_BH) - aqueous pressure of borehole node
!     PG_BH(LUK,LBN_BH) - gas pressure of borehole node
!     PORD_BH(LUK,LBN_BH) - diffusive porosity of borehole node
!     PSW_BH(LUK,LBN_BH) - saturated water pressure of borehole node
!     PVA_BH(LUK,LBN_BH) - partial pressure of air of borehole node
!     PVW_BH(LUK,LBN_BH) - partial pressure of water of borehole node
!     RHOG_BH(LUK,LBN_BH) - gas density of borehole node
!     RHOL_BH(LUK,LBN_BH) - aqueous density of borehole node
!     RHOMG_BH(LUK,LBN_BH) - gas molar density of borehole node
!     RHOML_BH(LUK,LBN_BH) - aqueous molar density of borehole node
!     RHOSP_BH(LUK,LBN_BH) - precipitated salt density of borehole node
!     RKG_BH(LUK,LBN_BH) - gas relative permeability of borehole node
!     RKL_BH(LUK,LBN_BH) - aqueous rel. permeability of borehole node
!     SG_BH(LUK,LBN_BH) - gas saturation of borehole node
!     SL_BH(LUK,LBN_BH) - aqueous saturation of borehole node
!     SS_BH(LUK,LBN_BH) - precipitated salt saturation of borehole node
!     T_BH(LUK,LBN_BH) - temperature of borehole node
!     THKG_BH(LUK,LBN_BH) - gas thermal conductivity of borehole node
!     THKL_BH(LUK,LBN_BH) - aqueous thermal conduc. of borehole node
!     T_BH(LUK,LBN_BH) - temperature of borehole node
!     TMS_BH(LUK,LBN_BH) - total salt mass of borehole node
!     UEG_BH(LUK,LBN_BH) - gas internal energy of borehole node
!     VISG_BH(LUK,LBN_BH) - gas viscosity of borehole node
!     VISL_BH(LUK,LBN_BH) - aqueous viscosity of borehole node
!     XGA_BH(LUK,LBN_BH) - gas air mass fraction of borehole node
!     XGW_BH(LUK,LBN_BH) - gas water mass fraction of borehole node
!     XLA_BH(LUK,LBN_BH) - aqueous air mass fraction of borehole node
!     XLS_BH(LUK,LBN_BH) - aqueous salt mass fraction of borehole node
!     XLW_BH(LUK,LBN_BH) - aqueous water mass fraction of borehole node
!     YLS_BH(LUK,LBN_BH) - aqueous total salt mass fraction of 
!       borehole node
!     XLW_BH(LUK,LBN_BH) - aqueous water mass fraction of borehole node
!     XMGA_BH(LUK,LBN_BH) - gas air mole fraction of borehole node
!     XMGW_BH(LUK,LBN_BH) - gas water mole fraction of borehole node
!     XMLA_BH(LUK,LBN_BH) - aqueous air mole fraction of borehole node
!     XMLS_BH(LUK,LBN_BH) - aqueous salt mole fraction of borehole node
!     XMLW_BH(LUK,LBN_BH) - aqueous water mole fraction of borehole node
!     NPHAZ_BH(LUK,LBN_BH) - phase condition of borehole node
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 April 2019.
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: DNR_BH,
     &  DFGW_BH,DFLA_BH,DFLS_BH,HGW_BH,HGA_BH,HG_BH,
     &  HL_BH,HSP_BH,PERMRF_BH,PL_BH,PG_BH,PSW_BH,PORD_BH,POR0_BH,
     &  PVA_BH,PVW_BH,RHOG_BH,RHOL_BH,RHOMG_BH,
     &  RHOML_BH,RHOSP_BH,RKG_BH,RKL_BH,SG_BH,SL_BH,SS_BH,
     &  T_BH,THKG_BH,THKL_BH,TMS_BH,UEG_BH,VISG_BH,VISL_BH,
     &  XGA_BH,XGW_BH,XLA_BH,XLS_BH,YLS_BH,XLW_BH,
     &  XMGA_BH,XMGW_BH,XMLA_BH,XMLS_BH,XMLW_BH
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NPHAZ_BH
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE FLUX_BH
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
!     Borehole node to borehole node flux variables
!
!     UBBC(LBC_BH,LSOLU+LSPT) solute flux 1/m^2 s
!     UBBDGW(LSFV,LBC_BH) diffusive water vapor flux through gas, 
!       kmol/m^2 s
!     UBBDLA(LSFV,LBC_BH) diffusive air flux through aqu., kmol/m^2 s
!     UBBDS(LSFV,LBC_BH) diffusive salt flux through aqu., kg/m^2 s
!     UBBS(LSFV,LBC_BH) salt flux through aqu., kg/m^2 s
!     UBBG(LSFV,LBC_BH) gas flux, m/s
!     UBBGW(LSFV,LBC_BH) water vapor flux through gas, kg/m^2 s
!     UBBL(LSFV,LBC_BH) aqueous flux, m/s
!     UBBL(LSFV,LBC_BH) nonaqueous flux, m/s
!     UBBQ(LSFV,LBC_BH) heat flux, W/m^2
!     UBBQC(LSFV,LBN_BH) coaxial heat flux, W/m^2
!
!     Borehole node to fracture triangle flux variables
!
!     UBFC(LFC_BH,LSOLU+LSPT) solute flux 1/m^2 s
!     UBFDGW(LSFV,LFC_BH) diffusive water vapor flux through gas, 
!       kmol/m^2 s
!     UBFDLA(LSFV,LFC_BH) diffusive air flux through aqu., kmol/m^2 s
!     UBFDS(LSFV,LFC_BH) diffusive salt flux through aqu., kg/m^2 s
!     UBFS(LSFV,LFC_BH) salt flux through aqu., kg/m^2 s
!     UBFG(LSFV,LFC_BH) gas flux, m/s
!     UBFGW(LSFV,LFC_BH) water vapor flux through gas, kg/m^2 s
!     UBFL(LSFV,LFC_BH) aqueous flux, m/s
!     UBFL(LSFV,LFC_BH) nonaqueous flux, m/s
!     UBFQ(LSFV,LFC_BH) heat flux, W/m^2
!
!     Borehole node to field node flux variables
!
!     UBML(LBN_BH) aqueous flow rate, m^3/s
!     UBMG(LBN_BH) gas flow rate, m^3/s
!     UBMT(LBN_BH) heat flow rate, W
!     UBMC(LBN_BH,LSOLU+LSPT) tracer flow rate, 1/s
!     TRNSA_BH(LSFV,LBN_BH) air mass transfer rate between borehole node
!       and grid cell, kg/s
!     TRNSS_BH(LSFV,LBN_BH) salt mass transfer rate between borehole 
!       node and grid cell, kg/s
!     TRNSQ_BH(LSFV,LBN_BH) energy transfer rate between borehole node
!       and grid cell, kg/s
!     TRNSW_BH(LSFV,LBN_BH) water mass transfer rate between borehole 
!       node and grid cell, kg/s
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 April 2019.
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UBBG,UBBL,UBBN,
     &  UBBGW,UBBDGW,UBBDLA,UBBDS,UBBS,UBBQ,UBBQC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UBBC
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UBFG,UBFL,UBFN,
     &  UBFGW,UBFDGW,UBFDLA,UBFDS,UBFS,UBFQ
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UBFC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: UBML,UBMG,UBMT
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: TRNSA_BH,TRNSS_BH,
     &  TRNSQ_BH,TRNSW_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: UBMC
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE TRNS_BH
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
!     Borehole transport variables
!
!     C_BH(LBN_BH,LSOLU+LSPT) volumetric solute or total 
!       species concentration of borehole node
!     CO_BH(LBN_BH,LSOLU+LSPT) volumetric solute or total 
!       species concentration of borehole node at old time step
!     CRNTG_BH(LBN_BH) gas Courant number for fracture
!     CRNTL_BH(LBN_BH) aqueous Courant number for fracture
!     CRNTN_BH(LBN_BH) nonaqueous-liquid Courant number for fracture
!     ICT_BH(LBN_BH,LSOLU+LSPT) initial condition index of solute or 
!       total parameters species of borehole node
!     SP_C_BH(LBN_BH,LSPR) - fracture species concentration, 
!       mol/m^3 node volume
!     SP_CO_BH(LBN_BH,LSPR) - old fracture species concentration, 
!       mol/m^3 node volume
!     YG_BH(LBN_BHC,LSOLU+LSPT) gas mole fraction of solute or total 
!       species of borehole node
!     YL_BH(LBN_BHC,LSOLU+LSPT) aqueous mole fraction of solute or total 
!       species of borehole node
!     YN_BH(LBN_BHC,LSOLU+LSPT) nonaqueous liquid mole fraction of 
!       solute or total species of borehole node
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: C_BH,CO_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SP_C_BH,SP_CO_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YG_BH,YL_BH,YN_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CRNTL_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CRNTG_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: CRNTN_BH
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICT_BH
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE GEOM_BH
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
!     Borehole geometry variables
!
!     VOL_BH(:) volume of borehole node
!     DBF_BH(:) distance of borehole node to fracture triangle 
!       connection
!     DBB_BH(:) distance of borehole node to borehole node connection
!     DBBM_BH(:) distance from borehole node to surface connection
!     DBN_BH(LBN_BH) distance from borehole node to field node
!     DIST_BH(LBN_BH) distance along borehole
!     PLX_BH(LBN_BH) - x-direction well projection on local grid
!     PLY_BH(LBN_BH) - y-direction well projection on local grid
!     PLZ_BH(LBN_BH) - z-direction well projection on local grid
!     SRFN_BH(3,LBC_BH) - x,y,z components of unit surface normal
!     XE_BH(:,:) x vertices (10x radius) of borehole node, m
!     YE_BH(:,:) y vertices (10x radius) of borehole node, m
!     ZE_BH(:,:) z vertices (10x radius) of borehole node, m
!     XP_BH(:,:) x starting and ending point of borehole node, m
!     YP_BH(:,:) y starting and ending point of borehole node, m
!     ZP_BH(:,:) z starting and ending point of borehole node, m
!     XTP_BH(2,LI_BH) - x-transition points for borehole interval, m
!     YTP_BH(2,LI_BH) - y-transition points for borehole interval, m
!     ZTP_BH(2,LI_BH) - z-transition points for borehole interval, m
!     XP_BF(LBC_FRC) x coor. borehole node to frac. tri. connection, m
!     YP_BF(LBC_FRC) y coor. borehole node to frac. tri. connection, m
!     ZP_BF(LBC_FRC) z coor. borehole node to frac. tri. connection, m
!     NM_BH(:) borehole name
!     IFCM_BH(LBN_BH) borehole node to grid cell connection map
!     IBCM_BH(LBC_BH) borehole node to borehole node connection map
!     IBN_BH(LBN_BH) - node number connected to borehole number
!     IT_BH(1,LN_BH) - type index for borehole start
!     IT_BH(2,LN_BH) - type index for borehole end
!     Borehole types
!     1 - gas mass injection w/ zero water
!     2 - gas mass injection w/ water mass frac.
!     3 - gas mass injection w/ water rel. humidity
!     4 - aqueous mass injection w/ zero air and salt
!     5 - aqueous mass injection w/ air and salt mass frac.
!     6 - aqueous mass injection w/ air and salt rel. sat.
!     7 - fixed pressure, state condition #1 (saturated)
!     8 - fixed pressure, state condition #2 (partially saturated)
!     9 - fixed pressure, state condition #3 (unsaturated)
!    11 - gas volumetric injection w/ zero water
!    12 - gas volumetric injection w/ water mass frac.
!    13 - gas volumetric injection w/ water rel. humidity
!    14 - aqueous volumetric injection w/ zero air and salt
!    15 - aqueous volumetric injection w/ air and salt mass frac.
!    16 - aqueous volumetric injection w/ air and salt rel. sat.
!    21 - Dirichlet energy
!    22 - Neumann energy
!    IS_BH(LI_BH) borehole interval type
!    Borehole interval types
!    1 - cased filled borehole interval
!    2 - uncased filled borehole interval
!    11 - unperforated pipe interval
!    12 - perforated pipe interval
!    21 - unperforated catalyst pipe interval
!    NBN_BH - number of borehole nodes
!    NXP_BH - number of inactive borehole nodes
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 April 2019.
!
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRFN_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: VOL_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PLX_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PLY_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PLZ_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XP_BH,XE_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YP_BH,YE_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZP_BH,ZE_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XP_BF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: YP_BF
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ZP_BF
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: XTP_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: YTP_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ZTP_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DBF_BH,DBN_BH,DIST_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DBB_BH,DBBM_BH
      INTEGER, DIMENSION(:), ALLOCATABLE :: INV_BH,IBCM_BH
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IPB_BH,IT_BH
      INTEGER, DIMENSION(:), ALLOCATABLE :: IS_BH,IZ_BH
      INTEGER, DIMENSION(:), ALLOCATABLE :: IBN_BH
      INTEGER :: N_BH,NBN_BH,NXP_BH
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: NM_BH
!
!---  End of module  ---
!
      END MODULE

!---------------------Fortran 90 Module--------------------------------!
!
      MODULE PARM_BH
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
!     Borehole parameter variables
!
!     HR_BH(LBN_BH) - heat of thermocatalytic reaction, J/m^3 aqu
!     PCMP_BH(LBN_BH) - reference pressure for compressibility
!     PAR_BH(1,LI_BH) - skin factor for borehole interval
!     PAR_BH(1,LI_BH) - skin factor for borehole interval
!     PAR_BH(2,LI_BH) - borehole outer radius for interval, m
!     PAR_BH(3,LI_BH) - borehole inner radius for interval, m
!     PAR_BH(4,LI_BH) - borehole wall thermal conductivity, W/m ˚K
!     PAR_BH(5,LI_BH) - borehole wall volumetric heat capacity, J/m^3 ˚K
!     PAR_BH(6,LI_BH) - borehole Colebrook roughness, m
!     PAR_BH(7,LI_BH) - borehole Brooks and Corey psi, m
!     PAR_BH(8,LI_BH) - borehole Brooks and Corey lambda
!     PAR_BH(9,LI_BH) - borehole Brooks and Corey slr
!     PAR_BH(10,LI_BH) - borehole Webb oven dried head, m
!     PAR_BH(11,LI_BH) - borehole Webb saturation matching point
!     PAR_BH(12,LI_BH) - borehole Webb head matching point, m
!     PTC_BH(1) - thermo-catalytic max. relative enthalpy
!     PTC_BH(2) - thermo-catalytic transition temperature, ˚C
!     PTC_BH(3) - thermo-catalytic rise rate
!     VAR_BH(1,LT_BH,LN_BH) - borehole time for time point, s
!     VAR_BH(2,LT_BH,LN_BH) - mass rate for time point, kg/s
!     VAR_BH(3,LT_BH,LN_BH) - production pressure, Pa
!     VAR_BH(4,LT_BH,LN_BH) - aqueous air conc. for time point
!     VAR_BH(4,LT_BH,LN_BH) - gas water conc. for time point
!     VAR_BH(5,LT_BH,LN_BH) - temperature for time point
!     VAR_BH(6,LT_BH,LN_BH) - aqu. salt conc. for time point
!     VIC_BH(12,NF) borehole initial condition variables
!     ITM_BH(LN_BH) - number of time points for borehole
!     ICC_BH(LN_BH) - cyclic borehole time index for borehole
!     IC_BH(NF) borehole initial condition state index
!       1 - aqueous saturation, gas pressure
!       2 - aqueous saturation, aqueous pressure
!       3 - aqueous pressure, gas pressure
!       4 - hydrostatic
!     ID_BH(1,NBH) - starting borehole interval index
!     ID_BH(2,NBH) - ending borehole interval index
!     ID_BH(3,NBH) - starting borehole node index
!     ID_BH(4,NBH) - ending borehole node index
!     ID_BH(5,NBH) - starting frac. triangle index, pointer to 
!       IBHT_FRC, and IBHN_FRC
!     ID_BH(6,NBH) - ending frac. triangle index, pointer to 
!       IBHT_FRC, and IBHN_FRC
!     ITC_BH - thermo-catalytic fluid equation form
!       22 - Zero-order reaction without catalyst with salt tracking
!       23 - First-order reaction without catalyst with salt tracking
!       24 - Zero-order reaction with catalyst with salt tracking
!       25 - First-order reaction with catalyst with salt tracking
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE DB_PR
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PCMP_BH,PTC_BH,HR_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: PAR_BH
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: VAR_BH
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: VARC_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: VIC_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: SRCA_BH,SRCS_BH,
     &  SRCT_BH,SRCW_BH
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RSDL_BH
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RSD_BH
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: SRC_BH
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISOLU_BH
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICC_BH,IC_BH,ITM_BH
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ID_BH,IM_BH
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISRM_BH
      INTEGER, DIMENSION(:), ALLOCATABLE :: NDREF_BH,NSD_BH
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISRT_BH
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISRDM_BH
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU_BH
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLUC_BH
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLU_MCB,KLU_BCM,KLU_BCF
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: KLUC_MCB,KLUC_BCM,KLUC_BCF
      INTEGER :: ITC_BH,ISPH_BH,NREF_BH,NSR_BH
!
!---  End of module  ---
!
      END MODULE


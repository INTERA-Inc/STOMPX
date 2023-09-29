!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ADVBB( PORDX,PORDBX,SX,SBX,U,V,W,UBX,VBX,WBX,N,M )
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
!     Computes the average Darcy velocity on a boundary bottom surface.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, August 1993.
!     Last Modified by MD White, Battelle, PNL, April 8, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 U(LSFV,*),V(LSFV,*),W(LSFV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ADVBB'
      UBX = 0.5D+0*(U(1,NSX(N))+U(1,NSX(N)+1))
      VBX = 0.5D+0*(V(1,NSY(N))+V(1,NSY(N)+IFLD))
      WBX = W(M,NSZ(N))
      IF( ISLC(52).EQ.0 ) THEN
        VMCP = SX*PORDX
        VMCB = SBX*PORDBX
        INDX = 17
        VMCX = DIFMN(VMCB,VMCP,DZGF(N),DZGF(N),WBX,INDX)
        UBX = UBX/VMCX
        VBX = VBX/VMCX
        WBX = WBX/VMCX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ADVBB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ADVSB( PORDX,PORDBX,SX,SBX,U,V,W,USX,VSX,WSX,N,M )
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
!     Computes the average Darcy velocity on a boundary south surface.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, August 1993.
!     Last Modified by MD White, Battelle, PNL, October 11, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 U(LSFV,*),V(LSFV,*),W(LSFV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ADVSB'
      USX = 0.5D+0*(U(1,NSX(N))+U(1,NSX(N)+1))
      VSX = V(M,NSY(N))
      WSX = 0.5D+0*(W(1,NSZ(N))+W(1,NSZ(N)+IJFLD))
      IF( ISLC(52).EQ.0 ) THEN
        VMCP = SX*PORDX
        VMCS = SBX*PORDBX
        INDX = 17
        VMCX = DIFMN(VMCS,VMCP,DYGF(N),DYGF(N),VSX,INDX)
        USX = USX/VMCX
        VSX = VSX/VMCX
        WSX = WSX/VMCX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ADVSB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ADVWB( PORDX,PORDBX,SX,SBX,U,V,W,UWX,VWX,WWX,N,M )
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
!     Computes the average Darcy velocity on a boundary west surface.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, August 1993.
!     Last Modified by MD White, Battelle, PNL, October 11, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 U(LSFV,*),V(LSFV,*),W(LSFV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ADVWB'
      UWX = U(M,NSX(N))
      VWX = 0.5D+0*(V(1,NSY(N))+V(1,NSY(N)+IFLD))
      WWX = 0.5D+0*(W(1,NSZ(N))+W(1,NSZ(N)+IJFLD))
      IF( ISLC(52).EQ.0 ) THEN
        VMCP = SX*PORDX
        VMCW = SBX*PORDBX
        INDX = 17
        VMCX = DIFMN(VMCW,VMCP,DXGF(N),DXGF(N),UWX,INDX)
        UWX = UWX/VMCX
        VWX = VWX/VMCX
        WWX = WWX/VMCX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ADVWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ADVEB( PORDX,PORDBX,SX,SBX,U,V,W,UEX,VEX,WEX,N,M )
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
!     Computes the average Darcy velocity on a boundary east surface.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, August 1993.
!     Last Modified by MD White, Battelle, PNL, October 11, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 U(LSFV,*),V(LSFV,*),W(LSFV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ADVEB'
      UEX = U(M,NSX(N)+1)
      VEX = 0.5D+0*(V(1,NSY(N))+V(1,NSY(N)+IFLD))
      WEX = 0.5D+0*(W(1,NSZ(N))+W(1,NSZ(N)+IJFLD))
      IF( ISLC(52).EQ.0 ) THEN
        VMCP = SX*PORDX
        VMCE = SBX*PORDBX
        INDX = 17
        VMCX = DIFMN(VMCP,VMCE,DXGF(N),DXGF(N),UEX,INDX)
        UEX = UEX/VMCX
        VEX = VEX/VMCX
        WEX = WEX/VMCX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ADVEB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ADVNB( PORDX,PORDBX,SX,SBX,U,V,W,UNX,VNX,WNX,N,M )
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
!     Computes the average Darcy velocity on a boundary north surface.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, August 1993.
!     Last Modified by MD White, Battelle, PNL, October 11, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 U(LSFV,*),V(LSFV,*),W(LSFV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ADVNB'
      UNX = 0.5D+0*(U(1,NSX(N))+U(1,NSX(N)+1))
      VNX = V(M,NSY(N)+IFLD)
      WNX = 0.5D+0*(W(1,NSZ(N))+W(1,NSZ(N)+IJFLD))
      IF( ISLC(52).EQ.0 ) THEN
        VMCP = SX*PORDX
        VMCN = SBX*PORDBX
        INDX = 17
        VMCX = DIFMN(VMCP,VMCN,DYGF(N),DYGF(N),VNX,INDX)
        UNX = UNX/VMCX
        VNX = VNX/VMCX
        WNX = WNX/VMCX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ADVNB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ADVTB( PORDX,PORDBX,SX,SBX,U,V,W,UTX,VTX,WTX,N,M )
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
!     Computes the average Darcy velocity on a boundary top surface.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, August 1993.
!     Last Modified by MD White, Battelle, PNL, August 17, 1993.
!     advbb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 U(LSFV,*),V(LSFV,*),W(LSFV,*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ADVTB'
      UTX = 0.5D+0*(U(1,NSX(N))+U(1,NSX(N)+1))
      VTX = 0.5D+0*(V(1,NSY(N))+V(1,NSY(N)+IFLD))
      WTX = W(M,NSZ(N)+IJFLD)
      IF( ISLC(52).EQ.0 ) THEN
        VMCP = SX*PORDX
        VMCT = SBX*PORDBX
        INDX = 17
        VMCX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(N),WTX,INDX)
        UTX = UTX/VMCX
        VTX = VTX/VMCX
        WTX = WTX/VMCX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ADVTB group  ---
!
      RETURN
      END

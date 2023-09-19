!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWB( N,NB )
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
!     Compute water vapor mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 22, 1998.
!     dffgwb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWB'
      NPZ = NSZ(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      K = KD(N)
      DXMGW = XMGWB(2,NB)-XMGW(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &   *DFGW(MP,N)*RHOMG(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)
     &   *DFGWB(MP,NB)*RHOMGB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGW,INDX )
       WDGW(M,NPZ) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))
     &   /DZGF(N)
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( WG(1,NPZ).LT.-EPSL ) THEN
           WGW(M,NPZ) = WG(M,NPZ)*FGW
         ELSE
           WGW(M,NPZ) = 0.D+0
         ENDIF
       ELSE
         WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      K = KD(N)
      DXMGW = XMGWB(2,NB)*RHOMGB(2,NB) - XMGW(2,N)*RHOMG(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGW(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGW,INDX )
       WDGW(M,NPZ) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB)
     &   - XMGW(MP,N)*RHOMG(MP,N))/DZGF(N)
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( WG(1,NPZ).LT.-EPSL ) THEN
           WGW(M,NPZ) = WG(M,NPZ)*FGW
         ELSE
           WGW(M,NPZ) = 0.D+0
         ENDIF
       ELSE
         WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
       ENDIF
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWE( N,NB )
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
!     Compute water vapor mole diffusion rates on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 22, 1998.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWE'
      NQX = NSX(N)+1
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      I = ID(N)
      DXMGW = XMGW(2,N)-XMGWB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &   *DFGW(MN,N)*RHOMG(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*
     &   DFGWB(MN,NB)*RHOMGB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGW,INDX )
       UDGW(M,NQX) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))
     &   /DXGF(N)
       FGWP = XGW(MN,N)*RHOG(MN,N)
       FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
       INDX = 3
       FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( UG(1,NQX).GT.EPSL ) THEN
           UGW(M,NQX) = UG(M,NQX)*FGW
         ELSE
           UGW(M,NQX) = 0.D+0
         ENDIF
       ELSE
         UGW(M,NQX) = UG(M,NQX)*FGW + WTMW*UDGW(M,NQX)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGW(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGW,INDX )
       UDGW(M,NQX) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N)
     &   - XMGWB(MN,NB)*RHOMGB(MN,NB))/DXGF(N)
       FGWP = XGW(MN,N)*RHOG(MN,N)
       FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
       INDX = 3
       FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( UG(1,NQX).GT.EPSL ) THEN
           UGW(M,NQX) = UG(M,NQX)*FGW
         ELSE
           UGW(M,NQX) = 0.D+0
         ENDIF
       ELSE
         UGW(M,NQX) = UG(M,NQX)*FGW + WTMW*UDGW(M,NQX)
       ENDIF
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWN( N,NB )
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
!     Compute water vapor mole diffusion rates on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 22, 1998.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWN'
      NQY = NSY(N)+IFLD
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      I = ID(N)
      J = JD(N)
      DXMGW = XMGW(2,N)-XMGWB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &   *DFGW(MN,N)*RHOMG(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*
     &   DFGWB(MN,NB)*RHOMGB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGW,INDX )
       VDGW(M,NQY) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))
     &   /(DYGF(N)*RP(I))
       FGWP = XGW(MN,N)*RHOG(MN,N)
       FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
       INDX = 3
       FGW = DIFMN( FGWP,FGWB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( VG(1,NQY).GT.EPSL ) THEN
           VGW(M,NQY) = VG(M,NQY)*FGW
         ELSE
           VGW(M,NQY) = 0.D+0
         ENDIF
       ELSE
         VGW(M,NQY) = VG(M,NQY)*FGW + WTMO*VDGO(M,NQY)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      J = JD(N)
      DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGW(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGW,INDX )
       VDGW(M,NQY) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N)
     &   - XMGWB(MN,NB)*RHOMGB(MN,NB))/(DYGF(N)*RP(I))
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWS( N,NB )
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
!     Compute water vapor mole diffusion rates on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 22, 1998.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWS'
      NPY = NSY(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      I = ID(N)
      J = JD(N)
      DXMGW = XMGWB(2,NB)-XMGW(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &   *DFGW(MP,N)*RHOMG(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*
     &   DFGWB(MP,NB)*RHOMGB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGW,INDX )
       VDGW(M,NPY) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))/
     &   (DYGF(N)*RP(I))
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( VG(1,NPY).LT.-EPSL ) THEN
           VGW(M,NPY) = VG(M,NPY)*FGW
         ELSE
           VGW(M,NPY) = 0.D+0
         ENDIF
       ELSE
         VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      J = JD(N)
      DXMGW = XMGWB(2,NB)*RHOMGB(2,NB) - XMGW(2,N)*RHOMG(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGW(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGW,INDX )
       VDGW(M,NPY) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB)
     &   - XMGW(MP,N)*RHOMG(MP,N))/(DYGF(N)*RP(I))
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( VG(1,NPY).LT.-EPSL ) THEN
           VGW(M,NPY) = VG(M,NPY)*FGW
         ELSE
           VGW(M,NPY) = 0.D+0
         ENDIF
       ELSE
         VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
       ENDIF
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWT( N,NB )
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
!     Compute water vapor mole diffusion rates on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 22, 1998.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWT'
      NQZ = NSZ(N)+IJFLD
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      K = KD(N)
      DXMGW = XMGW(2,N)-XMGWB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &   *DFGW(MN,N)*RHOMG(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*
     &   DFGWB(MN,NB)*RHOMGB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGW,INDX )
       WDGW(M,NQZ) = 2.D+0*DFM*(XMGW(MN,N)-
     &   XMGWB(MN,NB))/DZGF(N)
       FGWP = XGW(MN,N)*RHOG(MN,N)
       FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
       INDX = 3
       FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( WG(1,NQZ).GT.EPSL ) THEN
           WGW(M,NQZ) = WG(M,NQZ)*FGW
         ELSE
           WGW(M,NQZ) = 0.D+0
         ENDIF
       ELSE
         WGW(M,NQZ) = WG(M,NQZ)*FGW + WTMW*WDGW(M,NQZ)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      K = KD(N)
      DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGW(MN,N)
       DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGW,INDX )
       WDGW(M,NQZ) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N)
     &   - XMGWB(MN,NB)*RHOMGB(MN,NB))/DZGF(N)
       FGWP = XGW(MN,N)*RHOG(MN,N)
       FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
       INDX = 3
       FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( WG(1,NQZ).GT.EPSL ) THEN
           WGW(M,NQZ) = WG(M,NQZ)*FGW
         ELSE
           WGW(M,NQZ) = 0.D+0
         ENDIF
       ELSE
         WGW(M,NQZ) = WG(M,NQZ)*FGW + WTMW*WDGW(M,NQZ)
       ENDIF
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWW( N,NB )
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
!     Compute water vapor mole diffusion rates on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 22, 1998.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVG
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWW'
      NPX = NSX(N)
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
      I = ID(N)
      DXMGW = XMGWB(2,NB)-XMGW(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &   *DFGW(MP,N)*RHOMG(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*
     &   DFGWB(MP,NB)*RHOMGB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGW,INDX )
       UDGW(M,NPX) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))/DXGF(N)
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( UG(1,NPX).LT.-EPSL ) THEN
           UGW(M,NPX) = UG(M,NPX)*FGW
         ELSE
           UGW(M,NPX) = 0.D+0
         ENDIF
       ELSE
         UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
       ENDIF
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      DXMGW = XMGWB(2,NB)*RHOMGB(2,NB)-XMGW(2,N)*RHOMG(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGW(MP,N)
       DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGW,INDX )
       UDGW(M,NPX) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB)
     &   - XMGW(MP,N)*RHOMG(MP,N))/DXGF(N)
       FGWP = XGW(MP,N)*RHOG(MP,N)
       FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
       INDX = 3
       FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---   Dirichlet-Outflow boundary condition  ---
!
       IF( IBCT(IEQA,NB).EQ.19 ) THEN
         IF( UG(1,NPX).LT.-EPSL ) THEN
           UGW(M,NPX) = UG(M,NPX)*FGW
         ELSE
           UGW(M,NPX) = 0.D+0
         ENDIF
       ELSE
         UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
       ENDIF
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWW group  ---
!
      RETURN
      END

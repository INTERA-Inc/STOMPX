!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGOB( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 22, 1998.
!     dffgob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE FLUXD
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
      SUB_LOG(ISUB_LOG) = '/DFFGOB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMGO = XMGOB(2,NB)-XMGO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGO(MP,N)
     &    *RHOMG(MP,N)
        DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGOB(MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGO,INDX )
        WDGO(M,NPZ) = 2.D+0*DFM*(XMGOB(MP,NB)-XMGO(MP,N))
     &    /DZGF(N)
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGO(M,NPZ) = WG(M,NPZ)*FGO
          ELSE
            WGO(M,NPZ) = 0.D+0
          ENDIF
        ELSE
          WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGOE( N,NB )
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
!     Compute VOC vapor mole diffusion rates on an east boundary.
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
      USE FLUXD
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
      SUB_LOG(ISUB_LOG) = '/DFFGOE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMGO = XMGO(2,N)-XMGOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGO(MN,N)
     &    *RHOMG(MN,N)
        DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGOB(MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGO,INDX )
        UDGO(M,NQX) = 2.D+0*DFM*(XMGO(MN,N)-XMGOB(MN,NB))
     &    /DXGF(N)
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGO(M,NQX) = UG(M,NQX)*FGO
          ELSE
            UGO(M,NQX) = 0.D+0
          ENDIF
        ELSE
          UGO(M,NQX) = UG(M,NQX)*FGO + WTMO*UDGO(M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGON( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a north boundary.
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
      USE FLUXD
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
      SUB_LOG(ISUB_LOG) = '/DFFGON'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGO = XMGO(2,N)-XMGOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGO(MN,N)
     &    *RHOMG(MN,N)
        DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGOB(MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGO,INDX )
        VDGO(M,NQY) = 2.D+0*DFM*(XMGO(MN,N)-XMGOB(MN,NB))
     &    /(DYGF(N)*RP(I))
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGO(M,NQY) = VG(M,NQY)*FGO
          ELSE
            VGO(M,NQY) = 0.D+0
          ENDIF
        ELSE
          VGO(M,NQY) = VG(M,NQY)*FGO + WTMO*VDGO(M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGOS( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a south boundary.
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
      USE FLUXD
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
      SUB_LOG(ISUB_LOG) = '/DFFGOS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGO = XMGOB(2,NB)-XMGO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGO(MP,N)
     &    *RHOMG(MP,N)
        DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGOB(MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGO,INDX )
        VDGO(M,NPY) = 2.D+0*DFM*(XMGOB(MP,NB)-XMGO(MP,N))
     &    /(DYGF(N)*RP(I))
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGO(M,NPY) = VG(M,NPY)*FGO
          ELSE
            VGO(M,NPY) = 0.D+0
          ENDIF
        ELSE
          VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGOT( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a top boundary.
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
      USE FLUXD
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
      SUB_LOG(ISUB_LOG) = '/DFFGOT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGO = XMGO(2,N)-XMGOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGO(MN,N)
     &    *RHOMG(MN,N)
        DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGOB(MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGO,INDX )
        WDGO(M,NQZ) = 2.D+0*DFM*(XMGO(MN,N)-XMGOB(MN,NB))
     &   /DZGF(N)
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGO(M,NQZ) = WG(M,NQZ)*FGO
          ELSE
            WGO(M,NQZ) = 0.D+0
          ENDIF
        ELSE
          WGO(M,NQZ) = WG(M,NQZ)*FGO + WTMO*WDGO(M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGOW( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a west boundary.
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
      USE FLUXD
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
      SUB_LOG(ISUB_LOG) = '/DFFGOW'
      I = ID(N)
      NPX = NSX(N)
      DXMGO = XMGOB(2,NB)-XMGO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGO(MP,N)
     &    *RHOMG(MP,N)
        DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGOB(MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGO,INDX )
        UDGO(M,NPX) = 2.D+0*DFM*(XMGOB(MP,NB)-XMGO(MP,N))
     &    /DXGF(N)
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGO(M,NPX) = UG(M,NPX)*FGO
          ELSE
            UGO(M,NPX) = 0.D+0
          ENDIF
        ELSE
          UGO(M,NPX) = UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGOW group  ---
!
      RETURN
      END

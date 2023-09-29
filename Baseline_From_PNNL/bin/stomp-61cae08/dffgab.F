!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAB( N,NB )
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
!     Compute CO2 mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 3 March 2005.
!     Last Modified by MD White, PNNL, 3 March 2005.
!     dffgab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFGAB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMGA = XMGAB(2,NB)-XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGA(MP,N)
     &    *RHOMG(MP,N)
        DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGA,INDX )
        WDGA(M,NPZ) = 2.D+0*DFM*(XMGAB(MP,NB)-XMGA(MP,N))
     &    /DZGF(N)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGA(M,NPZ) = WG(M,NPZ)*FGA
          ELSE
            WGA(M,NPZ) = 0.D+0
          ENDIF
        ELSE
          WGA(M,NPZ) = WG(M,NPZ)*FGA + WTMA*WDGA(M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAE( N,NB )
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
!     Written by MD White, PNNL, 3 March 2005.
!     Last Modified by MD White, PNNL, 3 March 2005.
!     dffgab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFGAE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMGA = XMGA(2,N)-XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGA(MN,N)
     &    *RHOMG(MN,N)
        DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGA,INDX )
        UDGA(M,NQX) = 2.D+0*DFM*(XMGA(MN,N)-XMGAB(MN,NB))/DXGF(N)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGA(M,NQX) = UG(M,NQX)*FGA
          ELSE
            UGA(M,NQX) = 0.D+0
          ENDIF
        ELSE
          UGA(M,NQX) = UG(M,NQX)*FGA + WTMA*UDGA(M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAN( N,NB )
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
!     Written by MD White, PNNL, 3 March 2005.
!     Last Modified by MD White, PNNL, 3 March 2005.
!     dffgab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFGAN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGA = XMGA(2,N)-XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGA(MN,N)
     &    *RHOMG(MN,N)
        DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGA,INDX )
        VDGA(M,NQY) = 2.D+0*DFM*(XMGA(MN,N)-XMGAB(MN,NB))
     &    /(DYGF(N)*RP(I))
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGA(M,NQY) = VG(M,NQY)*FGA
          ELSE
            VGA(M,NQY) = 0.D+0
          ENDIF
        ELSE
          VGA(M,NQY) = VG(M,NQY)*FGA + WTMA*VDGA(M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAS( N,NB )
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
!     Written by MD White, PNNL, 3 March 2005.
!     Last Modified by MD White, PNNL, 3 March 2005.
!     dffgab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFGAS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGA = XMGAB(2,NB)-XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGA(MP,N)
     &    *RHOMG(MP,N)
        DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGA,INDX )
        VDGA(M,NPY) = 2.D+0*DFM*(XMGAB(MP,NB)-XMGA(MP,N))
     &    /(DYGF(N)*RP(I))
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGA(M,NPY) = VG(M,NPY)*FGA
          ELSE
            VGA(M,NPY) = 0.D+0
          ENDIF
        ELSE
          VGA(M,NPY) = VG(M,NPY)*FGA + WTMA*VDGA(M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAT( N,NB )
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
!     Written by MD White, PNNL, 3 March 2005.
!     Last Modified by MD White, PNNL, 3 March 2005.
!     dffgab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFGAT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGA = XMGA(2,N)-XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))*DFGA(MN,N)
     &    *RHOMG(MN,N)
        DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &    *RHOMGB(MN,NB)
        INDX = 13
        DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGA,INDX )
        WDGA(M,NQZ) = 2.D+0*DFM*(XMGA(MN,N)-XMGAB(MN,NB))
     &   /DZGF(N)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGA(M,NQZ) = WG(M,NQZ)*FGA
          ELSE
            WGA(M,NQZ) = 0.D+0
          ENDIF
        ELSE
          WGA(M,NQZ) = WG(M,NQZ)*FGA + WTMA*WDGA(M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAW( N,NB )
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
!     Written by MD White, PNNL, 3 March 2005.
!     Last Modified by MD White, PNNL, 3 March 2005.
!     dffgab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFGAW'
      I = ID(N)
      NPX = NSX(N)
      DXMGA = XMGAB(2,NB)-XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*DFGA(MP,N)
     &    *RHOMG(MP,N)
        DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &    *RHOMGB(MP,NB)
        INDX = 13
        DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGA,INDX )
        UDGA(M,NPX) = 2.D+0*DFM*(XMGAB(MP,NB)-XMGA(MP,N))
     &    /DXGF(N)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(IEQA,NB).EQ.19 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGA(M,NPX) = UG(M,NPX)*FGA
          ELSE
            UGA(M,NPX) = 0.D+0
          ENDIF
        ELSE
          UGA(M,NPX) = UG(M,NPX)*FGA + WTMA*UDGA(M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGB( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
!     drcvgb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVGB'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPZ = NSZ(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        HZ = PGB(MP,NB)-PG(MP,N)
     &    -5.D-1*GRVZ(NPZ)*DZGF(N)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
        WG(M,NPZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Diode boundary condition  ---
!
        IF( ABS(IBCT(IEQA,NB)).EQ.7 ) THEN
          WG(M,NPZ) = MIN( 0.D+0,WG(M,NPZ) )
!
!---    Dirichlet-outflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.26 ) THEN
          WG(M,NPZ) = MIN( 0.D+0,WG(M,NPZ) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.42 ) THEN
          WG(M,NPZ) = MAX( 0.D+0,WG(M,NPZ) )
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGS( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVGS'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPY = NSY(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        HY = PGB(MP,NB)-PG(MP,N)
     &    -5.D-1*GRVY(NPY)*DYGF(N)*RP(I)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
        VG(M,NPY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Diode boundary condition  ---
!
        IF( ABS(IBCT(IEQA,NB)).EQ.7 ) THEN
          VG(M,NPY) = MIN( 0.D+0,VG(M,NPY) )
!
!---    Dirchlet-outflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.26 ) THEN
          VG(M,NPY) = MIN( 0.D+0,VG(M,NPY) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.42 ) THEN
          VG(M,NPY) = MAX( 0.D+0,VG(M,NPY) )
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGW( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVGW'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPX = NSX(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        HX = PGB(MP,NB)-PG(MP,N)
     &    -5.D-1*GRVX(NPX)*DXGF(N)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
        UG(M,NPX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Diode boundary condition  ---
!
        IF( ABS(IBCT(IEQA,NB)).EQ.7 ) THEN
          UG(M,NPX) = MIN( 0.D+0,UG(M,NPX) )
!
!---    Dirichlet-outflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.26 ) THEN
          UG(M,NPX) = MIN( 0.D+0,UG(M,NPX) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.42 ) THEN
          UG(M,NPX) = MAX( 0.D+0,UG(M,NPX) )
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGE( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVGE'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQX = NSX(N)+1
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        HX = PG(MN,N)-PGB(MN,NB)
     &    -5.D-1*GRVX(NQX)*DXGF(N)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(1,IZ(N))
        UG(M,NQX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Diode boundary condition  ---
!
        IF( ABS(IBCT(IEQA,NB)).EQ.7 ) THEN
          UG(M,NQX) = MAX( 0.D+0,UG(M,NQX) )
!
!---    Dirichlet-outflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.26 ) THEN
          UG(M,NQX) = MAX( 0.D+0,UG(M,NQX) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.42 ) THEN
          UG(M,NQX) = MIN( 0.D+0,UG(M,NQX) )
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGN( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVGN'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQY = NSY(N)+IFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        HY = PG(MN,N)-PGB(MN,NB)
     &    -5.D-1*GRVY(NQY)*DYGF(N)*RP(I)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(2,IZ(N))
        VG(M,NQY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Diode boundary condition  ---
!
        IF( ABS(IBCT(IEQA,NB)).EQ.7 ) THEN
          VG(M,NQY) = MAX( 0.D+0,VG(M,NQY) )
!
!---    Dirichlet-outflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.26 ) THEN
          VG(M,NQY) = MAX( 0.D+0,VG(M,NQY))
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.42 ) THEN
          VG(M,NQY) = MIN( 0.D+0,VG(M,NQY))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVGT( N,NB )
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 13, 1995.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVGT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        HZ = PG(MN,N)-PGB(MN,NB)
     &    -5.D-1*GRVZ(NQZ)*DZGF(N)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(3,IZ(N))
        WG(M,NQZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Diode boundary condition  ---
!
        IF( ABS(IBCT(IEQA,NB)).EQ.7 ) THEN
          WG(M,NQZ) = MAX( 0.D+0,WG(M,NQZ) )
!
!---    Dirichlet-outflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.26 ) THEN
          WG(M,NQZ) = MAX( 0.D+0,WG(M,NQZ) )
!
!---    Dirichlet-inflow boundary condition  ---
!
        ELSEIF( ABS(IBCT(3,NB)).EQ.42 ) THEN
          WG(M,NQZ) = MIN( 0.D+0,WG(M,NQZ) )
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGT group  ---
!
      RETURN
      END

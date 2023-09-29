!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVG
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
!     and gravitational body forces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 KGM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVG'
!
!
!---  X-direction Darcy velocities, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
        DO 200 J = 1,JFLD
          DO 200 I = 2,IFLD
            N = ND(I,J,K)
            NW = N-1
            IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &        INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
            NPX = NSX(N)
            DO 100 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              HDGX = PG(MN,NW)-PG(MP,N) - 0.5D+0*GRVX(NPX)*
     &          (RHOG(MN,NW)*DXGF(N)+RHOG(MP,N)*DXGF(NW))
              IF( M.EQ.1 ) HDG = HDGX
!
!---          Permeability reduction factor
!
              PERM_WX = PERMRF(MN,NW)*PERM(1,IZ(NW))
              PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
              INDX = 11
              KGM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDG,INDX)
              IF( PERM_WX/EPSL.LT.EPSL ) KGM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
              INDX = 9
              RKGM = DIFMN(RKG(MN,NW),RKG(MP,N),DXGF(NW),DXGF(N),
     &          HDG,INDX)
              INDX = 6
              VGM = DIFMN(VISG(MN,NW),VISG(MP,N),DXGF(NW),DXGF(N),
     &          HDG,INDX)
              UG(M,NPX) = KGM*RKGM*HDGX/DXGP(NPX)/VGM
  100       CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      IF( IJFLD.GT.1 ) THEN
      DO 400 I = 1,IFLD
        DO 400 K = 1,KFLD
          DO 400 J = 2,JFLD
            N = ND(I,J,K)
            NS = N-IFLD
            IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &        INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
            NPY = NSY(N)
            DO 300 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              HDGY = PG(MN,NS)-PG(MP,N) - 0.5D+0*GRVY(NPY)*
     &          (RHOG(MN,NS)*DYGF(N)+RHOG(MP,N)*DYGF(NS))
              IF( M.EQ.1 ) HDG = HDGY
!
!---          Permeability reduction factor
!
              PERM_SX = PERMRF(MN,NS)*PERM(2,IZ(NS))
              PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
              INDX = 11
              KGM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDG,INDX)
              IF( PERM_SX/EPSL.LT.EPSL ) KGM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
              INDX = 9
              RKGM = DIFMN(RKG(MN,NS),RKG(MP,N),DYGF(NS),DYGF(N),
     &          HDG,INDX)
              INDX = 6
              VGM = DIFMN(VISG(MN,NS),VISG(MP,N),DYGF(NS),DYGF(N),
     &          HDG,INDX)
              VG(M,NPY) = KGM*RKGM*HDGY/DYGP(NPY)/VGM/RP(I)
  300       CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction Darcy velocities, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
      DO 600 J = 1,JFLD
        DO 600 I = 1,IFLD
          DO 600 K = 2,KFLD
            N = ND(I,J,K)
            NB = N-IJFLD
            IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &        INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
            NPZ = NSZ(N)
            DO 500 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              HDGZ = PG(MN,NB)-PG(MP,N) - 0.5D+0*GRVZ(NPZ)*
     &          (RHOG(MN,NB)*DZGF(N)+RHOG(MP,N)*DZGF(NB))
              IF( M.EQ.1 ) HDG = HDGZ
!
!---          Permeability reduction factor
!
              PERM_BX = PERMRF(MN,NB)*PERM(3,IZ(NB))
              PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
              INDX = 11
              KGM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDG,INDX)
              IF( PERM_BX/EPSL.LT.EPSL ) KGM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
              INDX = 9
              RKGM = DIFMN(RKG(MN,NB),RKG(MP,N),DZGF(NB),DZGF(N),
     &          HDG,INDX)
              INDX = 6
              VGM = DIFMN(VISG(MN,NB),VISG(MP,N),DZGF(NB),DZGF(N),
     &          HDG,INDX)
              WG(M,NPZ) = KGM*RKGM*HDGZ/DZGP(NPZ)/VGM
  500       CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVG group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVL
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FDVS
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
      REAL*8 KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVL'
!
!---  X-direction Darcy velocities, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
        DO 200 J = 1,JFLD
          DO 200 I = 2,IFLD
            N = ND(I,J,K)
            NW = N-1
            IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &        INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
            NPX = NSX(N)
            OECX = 0.D+0
            DO 100 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              INDX = 11
              IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &          OECX = DIFMN(OEC(MN,NW),OEC(MP,N),
     &          DXGF(NW),DXGF(N),ZERO,INDX)
              DPOS = OECX*(POSM(MN,NW)-POSM(MP,N))
              HDLX = PL(MN,NW)-PL(MP,N)+DPOS - 0.5D+0*GRVX(NPX)*
     &          (RHOL(MN,NW)*DXGF(N)+RHOL(MP,N)*DXGF(NW))
              IF( M.EQ.1 ) HDL = HDLX
!
!---          Permeability reduction factor
!
              PERM_WX = PERMRF(MN,NW)*PERM(1,IZ(NW))
              PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
              INDX = 11
              KLM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDL,INDX)
              INDX = 8
              RKLM = DIFMN(RKL(1,MN,NW),RKL(1,MP,N),DXGF(NW),
     &          DXGF(N),HDL,INDX)
              IF( PERM_WX/EPSL.LT.EPSL ) KLM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
              INDX = 5
              VLM = DIFMN(VISL(MN,NW),VISL(MP,N),DXGF(NW),DXGF(N),
     &          HDL,INDX)
              UL(M,NPX) = KLM*RKLM*HDLX/DXGP(NPX)/VLM
  100       CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
      DO 400 I = 1,IFLD
        DO 400 K = 1,KFLD
          DO 400 J = 2,JFLD
            N = ND(I,J,K)
            NS = N-IFLD
            IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &        INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
            NPY = NSY(N)
            OECY = 0.D+0
            DO 300 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              INDX = 11
              IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &          OECY = DIFMN(OEC(MN,NS),OEC(MP,N),
     &          DYGF(NS),DYGF(N),ZERO,INDX)
              DPOS = OECY*(POSM(MN,NS)-POSM(MP,N))
              HDLY = PL(MN,NS)-PL(MP,N)+DPOS - 0.5D+0*GRVY(NPY)*
     &          (RHOL(MN,NS)*DYGF(N)+RHOL(MP,N)*DYGF(NS))
              IF( M.EQ.1 ) HDL = HDLY
!
!---          Permeability reduction factor
!
              PERM_SX = PERMRF(MN,NS)*PERM(2,IZ(NS))
              PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
              INDX = 11
              KLM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDL,INDX)
              INDX = 8
              RKLM = DIFMN(RKL(2,MN,NS),RKL(2,MP,N),DYGF(NS),
     &          DYGF(N),HDL,INDX)
              IF( PERM_SX/EPSL.LT.EPSL ) KLM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
              INDX = 5
              VLM = DIFMN(VISL(MN,NS),VISL(MP,N),DYGF(NS),DYGF(N),
     &          HDL,INDX)
              VL(M,NPY) = KLM*RKLM*HDLY/DYGP(NPY)/VLM/RP(I)
  300       CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction Darcy velocities, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
      DO 600 J = 1,JFLD
        DO 600 I = 1,IFLD
          DO 600 K = 2,KFLD
            N = ND(I,J,K)
            NB = N-IJFLD
            IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &        INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
            NPZ = NSZ(N)
            OECZ = 0.D+0
            DO 500 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              INDX = 11
              IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &          OECZ = DIFMN(OEC(MN,NB),OEC(MP,N),
     &          DZGF(NB),DZGF(N),ZERO,INDX)
              DPOS = OECZ*(POSM(MN,NB)-POSM(MP,N))
              HDLZ = PL(MN,NB)-PL(MP,N)+DPOS - 0.5D+0*GRVZ(NPZ)*
     &          (RHOL(MN,NB)*DZGF(N)+RHOL(MP,N)*DZGF(NB))
              IF( M.EQ.1 ) HDL = HDLZ
!
!---          Permeability reduction factor
!
              PERM_BX = PERMRF(MN,NB)*PERM(3,IZ(NB))
              PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
              INDX = 11
              KLM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDL,INDX)
              INDX = 8
              RKLM = DIFMN(RKL(3,MN,NB),RKL(3,MP,N),DZGF(NB),
     &          DZGF(N),HDL,INDX)
              IF( PERM_BX/EPSL.LT.EPSL ) KLM = 0.D+0
              IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
              INDX = 5
              VLM = DIFMN(VISL(MN,NB),VISL(MP,N),DZGF(NB),DZGF(N),
     &          HDL,INDX)
              WL(M,NPZ) = KLM*RKLM*HDLZ/DZGP(NPZ)/VLM
  500       CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVL group  ---
!
      RETURN
      END
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
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          WG(M,NPZ) = MIN( 0.D+0,WG(M,NPZ) )
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
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          VG(M,NPY) = MIN( 0.D+0,VG(M,NPY) )
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
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          UG(M,NPX) = MIN( 0.D+0,UG(M,NPX) )
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
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          UG(M,NQX) = MAX( 0.D+0,UG(M,NQX) )
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
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          VG(M,NQY) = MAX( 0.D+0,VG(M,NQY) )
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
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          WG(M,NQZ) = MAX( 0.D+0,WG(M,NQZ) )
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVGT group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLB( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLB'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPZ = NSZ(N)
      OECZ = 0.D+0
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &    OECZ = DIFMN(OECB(MP,NB),OEC(MP,N),
     &    DZGF(N),DZGF(N),0.D+0,INDX)
        DPOS = OECZ*(POSB(MP,NB)-POSM(MP,N))
        HZ = PLB(MP,NB)-PL(MP,N)+DPOS
     &    -5.D-1*GRVZ(NPZ)*DZGF(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKLB(3,MP,NB),RKL(3,MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(3,IZ(N))
        WL(M,NSZ(N)) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          WL(M,NSZ(N)) = MIN( 0.D+0,WL(M,NSZ(N)) )
!
!---    Seepage or x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 
     &   .OR. ABS(IBCT(IEQW,NB)).EQ.45) THEN
           IF( PG(MP,N)-PL(MP,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MP,N)).GT.EPSL ) THEN
             IF( PLB(MP,NB)-PGB(MP,NB).GT.EPSL ) THEN
               WL(M,NPZ) = MAX( 0.D+0,WL(M,NPZ) )
             ELSE
               WL(M,NPZ) = MIN( 0.D+0,WL(M,NPZ)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MP,NB)-PGB(MP,NB).LT.EPSL ) THEN
               WL(M,NPZ) = MIN( 0.D+0,WL(M,NPZ) )
             ENDIF
           ENDIF
!          PLBX = PL(MP,N) + RHOL(MP,N)*5.D-1*GRVZ(NPZ)*DZGF(N)
!          IF( PLBX.GT.PLB(MP,NB) ) THEN
!            WL(M,NSZ(N)) = MIN( 0.D+0,WL(M,NSZ(N)) )
!          ELSE
!            WL(M,NSZ(N)) = 0.D+0
!          ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MP,N)-PG(2,N).GE.EPSL ) THEN
            WL(M,NSZ(N)) = MIN( 0.D+0,WL(M,NSZ(N)) )
          ELSE
            WL(M,NSZ(N)) = 0.D+0
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLS( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLS'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPY = NSY(N)
      OECY = 0.D+0
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &    OECY = DIFMN(OECB(MP,NB),OEC(MP,N),
     &    DYGF(N),DYGF(N),0.D+0,INDX)
        DPOS = OECY*(POSB(MP,NB)-POSM(MP,N))
        HY = PLB(MP,NB)-PL(MP,N)+DPOS
     &    -5.D-1*GRVY(NPY)*DYGF(N)*RP(I)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKLB(2,MP,NB),RKL(2,MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(2,IZ(N))
        VL(M,NPY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          VL(M,NPY) = MIN( 0.D+0,VL(M,NPY) )
!
!---    Seepage or x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 
     &   .OR. ABS(IBCT(IEQW,NB)).EQ.45 ) THEN
           IF( PG(MP,N)-PL(MP,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MP,N)).GT.EPSL ) THEN
             IF( PLB(MP,NB)-PGB(MP,NB).GT.EPSL ) THEN
               VL(M,NPY) = MAX( 0.D+0,VL(M,NPY) )
             ELSE
               VL(M,NPY) = MIN( 0.D+0,VL(M,NPY)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MP,NB)-PGB(MP,NB).LT.EPSL ) THEN
               VL(M,NPY) = MIN( 0.D+0,VL(M,NPY) )
             ENDIF
           ENDIF
!          PLBX = PL(MP,N) + RHOL(MP,N)*5.D-1*GRVY(NPY)*DYGF(N)*RP(I)
!          IF( PLBX.GT.PLB(MP,NB) ) THEN
!            VL(M,NPY) = MIN( 0.D+0,VL(M,NPY) )
!          ELSE
!            VL(M,NPY) = 0.D+0
!          ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MP,N)-PG(2,N).GE.EPSL ) THEN
            VL(M,NPY) = MIN( 0.D+0,VL(M,NPY) )
          ELSE
            VL(M,NPY) = 0.D+0
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLW( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLW'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NPX = NSX(N)
      OECX = 0.D+0
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &    OECX = DIFMN(OECB(MP,NB),OEC(MP,N),
     &    DXGF(N),DXGF(N),0.D+0,INDX)
        DPOS = OECX*(POSB(MP,NB)-POSM(MP,N))
        HX = PLB(MP,NB)-PL(MP,N)+DPOS
     &    -5.D-1*GRVX(NPX)*DXGF(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKLB(1,MP,NB),RKL(1,MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(1,IZ(N))
        UL(M,NPX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          UL(M,NPX) = MIN( 0.D+0,UL(M,NPX) )
!
!---    Seepage or x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 
     &   .OR. ABS(IBCT(IEQW,NB)).EQ.45 ) THEN
           IF( PG(MP,N)-PL(MP,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MP,N)).GT.EPSL ) THEN
             IF( PLB(MP,NB)-PGB(MP,NB).GT.EPSL ) THEN
               UL(M,NPX) = MAX( 0.D+0,UL(M,NPX) )
             ELSE
               UL(M,NPX) = MIN( 0.D+0,UL(M,NPX)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MP,NB)-PGB(MP,NB).LT.EPSL ) THEN
               UL(M,NPX) = MIN( 0.D+0,UL(M,NPX) )
             ENDIF
           ENDIF
!          PLBX = PL(MP,N) + RHOL(MP,N)*5.D-1*GRVX(NPX)*DXGF(N)
!          IF( PLBX.GT.PLB(MP,NB) ) THEN
!            UL(M,NPX) = MIN( 0.D+0,UL(M,NPX) )
!          ELSE
!            UL(M,NPX) = 0.D+0
!          ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MP,N)-PG(2,N).GE.EPSL ) THEN
            UL(M,NPX) = MIN( 0.D+0,UL(M,NPX) )
          ELSE
            UL(M,NPX) = 0.D+0
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLE( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLE'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQX = NSX(N)+1
      OECX = 0.D+0
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &    OECX = DIFMN(OEC(MN,N),OECB(MN,NB),
     &    DXGF(N),DXGF(N),0.D+0,INDX)
        DPOS = OECX*(POSM(MN,N)-POSB(MN,NB))
        HX = PL(MN,N)-PLB(MN,NB)+DPOS
     &    -5.D-1*GRVX(NQX)*DXGF(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKL(1,MN,N),RKLB(1,MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(1,IZ(N))
        UL(M,NQX) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          UL(M,NQX) = MAX( 0.D+0,UL(M,NQX) )
!
!---    Seepage or x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 
     &     .OR. ABS(IBCT(IEQW,NB)).EQ.45 ) THEN
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MN,N)).GT.EPSL ) THEN
             IF( PLB(MN,NB)-PGB(MN,NB).GT.EPSL ) THEN
               UL(M,NQX) = MIN( 0.D+0,UL(M,NQX) )
             ELSE
               UL(M,NQX) = MAX( 0.D+0,UL(M,NQX)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MN,NB)-PGB(MN,NB).LT.EPSL ) THEN
               UL(M,NQX) = MAX( 0.D+0,UL(M,NQX) )
             ENDIF
           ENDIF
!          PLBX = PL(MN,N) - RHOL(MN,N)*5.D-1*GRVX(NQX)*DXGF(N)
!          IF( PLBX.GT.PLB(MN,NB) ) THEN
!            UL(M,NQX) = MAX( 0.D+0,UL(M,NQX) )
!          ELSE
!            UL(M,NQX) = 0.D+0
!          ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MN,N)-PG(2,N).GE.EPSL ) THEN
            UL(M,NQX) = MAX( 0.D+0,UL(M,NQX) )
          ELSE
            UL(M,NQX) = 0.D+0
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLN( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLN'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQY = NSY(N)+IFLD
      OECY = 0.D+0
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &    OECY = DIFMN(OEC(MN,N),OECB(MN,NB),
     &    DYGF(N),DYGF(N),0.D+0,INDX)
        DPOS = OECY*(POSM(MN,N)-POSB(MN,NB))
        HY = PL(MN,N)-PLB(MN,NB)+DPOS
     &    -5.D-1*GRVY(NQY)*DYGF(N)*RP(I)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKL(2,MN,N),RKLB(2,MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(2,IZ(N))
        VL(M,NQY) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(I)*DYGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          VL(M,NQY) = MAX( 0.D+0,VL(M,NQY))
!
!---    Seepage or x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 
     &     .OR. ABS(IBCT(IEQW,NB)).EQ.45 ) THEN
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MN,N)).GT.EPSL ) THEN
             IF( PLB(MN,NB)-PGB(MN,NB).GT.EPSL ) THEN
               VL(M,NQY) = MIN( 0.D+0,VL(M,NQY))
             ELSE
               VL(M,NQY) = MAX( 0.D+0,VL(M,NQY)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MN,NB)-PGB(MN,NB).LT.EPSL ) THEN
               VL(M,NQY) = MAX( 0.D+0,VL(M,NQY))
             ENDIF
           ENDIF
!          PLBX = PL(MN,N) - RHOL(MN,N)*5.D-1*GRVY(NQY)*DYGF(N)*RP(I)
!          IF( PLBX.GT.PLB(MN,NB) ) THEN
!            VL(M,NQY) = MAX( 0.D+0,VL(M,NQY))
!          ELSE
!            VL(M,NQY) = 0.D+0
!          ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MN,N)-PG(2,N).GE.EPSL ) THEN
            VL(M,NQY) = MAX( 0.D+0,VL(M,NQY))
          ELSE
            VL(M,NQY) = 0.D+0
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLT( N,NB )
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by D.H. Bacon, PNNL, Oct. 31, 2012.
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
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLT'
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      OECZ = 0.D+0
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 )
     &    OECZ = DIFMN(OEC(MN,N),OECB(MN,NB),
     &    DZGF(N),DZGF(N),0.D+0,INDX)
        DPOS = OECZ*(POSM(MN,N)-POSB(MN,NB))
        HZ = PL(MN,N)-PLB(MN,NB)+DPOS
     &    -5.D-1*GRVZ(NQZ)*DZGF(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKL(3,MN,N),RKLB(3,MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(3,IZ(N))
        WL(M,NQZ) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          WL(M,NQZ) = MAX( 0.D+0,WL(M,NQZ) )
!
!---    Seepage or x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 
     &    .OR. ABS(IBCT(IEQW,NB)).EQ.45 ) THEN
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MN,N)).GT.EPSL ) THEN
             IF( PLB(MN,NB)-PGB(MN,NB).GT.EPSL ) THEN
               WL(M,NQZ) = MIN( 0.D+0,WL(M,NQZ) )
             ELSE
               WL(M,NQZ) = MAX( 0.D+0,WL(M,NQZ)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MN,NB)-PGB(MN,NB).LT.EPSL ) THEN
               WL(M,NQZ) = MAX( 0.D+0,WL(M,NQZ) )
             ENDIF
           ENDIF
!          PLBX = PL(MN,N) - RHOL(MN,N)*5.D-1*GRVZ(NQZ)*DZGF(N)
!          IF( PLBX.GT.PLB(MN,NB) ) THEN
!            WL(M,NQZ) = MAX( 0.D+0,WL(M,NQZ) )
!          ELSE
!            WL(M,NQZ) = 0.D+0
!          ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MN,N)-PG(2,N).GE.EPSL ) THEN
            WL(M,NQZ) = MAX( 0.D+0,WL(M,NQZ) )
          ELSE
            WL(M,NQZ) = 0.D+0
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLT group  ---
!
      RETURN
      END

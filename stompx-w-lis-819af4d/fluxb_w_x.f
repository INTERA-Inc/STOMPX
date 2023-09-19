!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWB_W( N,NB )
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
!     Diffusive water aqueous fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Water Mode (STOMPX-W)
!
!     Written by M.D. White, PNNL, 11 January 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUX
      USE FDVP
      USE CONST
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
      SUB_LOG(ISUB_LOG) = '/DFFLWB_W'
      DO M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,1,N),INDX )
        WLW(M,1,N) = WL(M,1,N)*FLW
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWB_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWE_W( N,NB )
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
!     Water Mode (STOMPX-W)
!
!     Diffusive water aqueous flux on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 January 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUX
      USE FDVP
      USE CONST
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
      SUB_LOG(ISUB_LOG) = '/DFFLWE_W'
      DO M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,2,N),INDX )
        ULW(M,2,N) = UL(M,2,N)*FLW
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWE_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWN_W( N,NB )
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
!     Water Mode (STOMPX-W)
!
!     Diffusive water aqueous fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 January 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUX
      USE FDVP
      USE CONST
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
      SUB_LOG(ISUB_LOG) = '/DFFLWN_W'
      DO M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,2,N),INDX )
        VLW(M,2,N) = VL(M,2,N)*FLW
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWN_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWS_W( N,NB )
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
!     Water Mode (STOMPX-W)
!
!     Diffusive water aqueous fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 January 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUX
      USE FDVP
      USE CONST
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
      SUB_LOG(ISUB_LOG) = '/DFFLWS_W'
      DO M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,1,N),INDX )
        VLW(M,1,N) = VL(M,1,N)*FLW
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWS_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWT_W( N,NB )
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
!     Water Mode (STOMPX-W)
!
!     Diffusive water aqueous fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 January 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUX
      USE FDVP
      USE CONST
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
      SUB_LOG(ISUB_LOG) = '/DFFLWT_W'
      DO M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,2,N),INDX )
        WLW(M,2,N) = WL(M,2,N)*FLW
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWT_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWW_W( N,NB )
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
!     Water Mode (STOMPX-W)
!
!     Diffusive water aqueous fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 January 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUX
      USE FDVP
      USE CONST
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
      SUB_LOG(ISUB_LOG) = '/DFFLWW_W'
      DO M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,1,N),INDX )
        ULW(M,1,N) = UL(M,1,N)*FLW
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWW_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLB_W( N,NB )
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
!     Written by M.D. White, PNNL, 3 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLB_W'
      OECZ = 0.D+0
      DO M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        HZ = PLB(MP,NB)-PL(MP,N)
     &    -5.D-1*GRVZ(1,N)*DZGF(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKLB(MP,NB),RKL(MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(3,N)
        WL(M,1,N) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          WL(M,1,N) = MIN( 0.D+0,WL(M,1,N) )
!
!---    Seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 ) THEN 
!
!---       Node is unsaturated  ---
!
           IF( PG(MP,N)-PL(MP,N).GT.EPSL ) THEN
!
!---         Boundary is saturated, flow allowed into node  ---
!
             IF( PLB(MP,NB)-PGB(MP,NB).GT.EPSL ) THEN
               WL(M,1,N) = MAX( 0.D+0,WL(M,1,N) )
!
!---         Boundary is unsaturated, zero flow  ---
!
             ELSE
               WL(M,1,N) = 0.D+0
             ENDIF
!
!---       Node saturated  ---
!
           ELSE
!
!---         Boundary is unsaturated, flow allowed from node,
!            boundary is saturated, flow allowed into and from node  ---
!
             IF( PLB(MP,NB)-PGB(MP,NB).LT.EPSL ) THEN
               WL(M,1,N) = MIN( 0.D+0,WL(M,1,N) )
             ENDIF
           ENDIF
!
!---    x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.45 ) THEN 
           IF( PG(MP,N)-PL(MP,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MP,N)).GT.EPSL ) THEN
             IF( PLB(MP,NB)-PGB(MP,NB).GT.EPSL ) THEN
               WL(M,1,N) = MAX( 0.D+0,WL(M,1,N) )
             ELSE
               WL(M,1,N) = MIN( 0.D+0,WL(M,1,N)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MP,NB)-PGB(MP,NB).LT.EPSL ) THEN
               WL(M,1,N) = MIN( 0.D+0,WL(M,1,N) )
             ENDIF
           ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MP,N)-PG(2,N).GE.EPSL ) THEN
            WL(M,1,N) = MIN( 0.D+0,WL(M,1,N) )
          ELSE
            WL(M,1,N) = 0.D+0
          ENDIF
!
!---    Dirichlet-Outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( WL(1,1,N).LT.-EPSL ) THEN
            WL(M,1,N) = WL(M,1,N)
          ELSE
            WL(M,1,N) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLB_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLS_W( N,NB )
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
!     Written by M.D. White, PNNL, 3 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLS_W'
      OECY = 0.D+0
      DO M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        HY = PLB(MP,NB)-PL(MP,N)
     &    -5.D-1*GRVY(1,N)*DYGF(N)*RP(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKLB(MP,NB),RKL(MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(2,N)
        VL(M,1,N) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(N)*DYGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          VL(M,1,N) = MIN( 0.D+0,VL(M,1,N) )
!
!---    Seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 ) THEN 
!
!---       Node is unsaturated  ---
!
           IF( PG(MP,N)-PL(MP,N).GT.EPSL ) THEN
!
!---         Boundary is saturated, flow allowed into node  ---
!
             IF( PLB(MP,NB)-PGB(MP,NB).GT.EPSL ) THEN
               VL(M,1,N) = MAX( 0.D+0,VL(M,1,N) )
!
!---         Boundary is unsaturated, zero flow  ---
!
             ELSE
               VL(M,1,N) = 0.D+0
             ENDIF
!
!---       Node saturated  ---
!
           ELSE
!
!---         Boundary is unsaturated, flow allowed from node,
!            boundary is saturated, flow allowed into and from node  ---
!
             IF( PLB(MP,NB)-PGB(MP,NB).LT.EPSL ) THEN
               VL(M,1,N) = MIN( 0.D+0,VL(M,1,N) )
             ENDIF
           ENDIF
!
!---    x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.45) THEN 
           IF( PG(MP,N)-PL(MP,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MP,N)).GT.EPSL ) THEN
             IF( PLB(MP,NB)-PGB(MP,NB).GT.EPSL ) THEN
               VL(M,1,N) = MAX( 0.D+0,VL(M,1,N) )
             ELSE
               VL(M,1,N) = MIN( 0.D+0,VL(M,1,N)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MP,NB)-PGB(MP,NB).LT.EPSL ) THEN
               VL(M,1,N) = MIN( 0.D+0,VL(M,1,N) )
             ENDIF
           ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MP,N)-PG(2,N).GE.EPSL ) THEN
            VL(M,1,N) = MIN( 0.D+0,VL(M,1,N) )
          ELSE
            VL(M,1,N) = 0.D+0
          ENDIF
!
!---    Dirichlet-Outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( VL(1,1,N).LT.-EPSL ) THEN
            VL(M,1,N) = VL(M,1,N)
          ELSE
            VL(M,1,N) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLS_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLW_W( N,NB )
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
!     Written by M.D. White, PNNL, 3 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLW_W'
      OECX = 0.D+0
      DO M = 1,ISVF
        MP = MPOS(M)
        INDX = 11
        HX = PLB(MP,NB)-PL(MP,N)
     &    -5.D-1*GRVX(1,N)*DXGF(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKLB(MP,NB),RKL(MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(1,N)
        UL(M,1,N) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          UL(M,1,N) = MIN( 0.D+0,UL(M,1,N) )
!
!---    Seepage seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 ) THEN 
!
!---       Node is unsaturated  ---
!
          IF( PG(MP,N)-PL(MP,N).GT.EPSL ) THEN
!
!---         Boundary is saturated, flow allowed into node  ---
!
             IF( PLB(MP,NB)-PGB(MP,NB).GT.EPSL ) THEN
               UL(M,1,N) = MAX( 0.D+0,UL(M,1,N) )
!
!---         Boundary is unsaturated, zero flow  ---
!
             ELSE
               UL(M,1,N) = 0.D+0
             ENDIF
!
!---       Node saturated  ---
!
           ELSE
!
!---         Boundary is unsaturated, flow allowed from node,
!            boundary is saturated, flow allowed into and from node  ---
!
             IF( PLB(MP,NB)-PGB(MP,NB).LT.EPSL ) THEN
               UL(M,1,N) = MIN( 0.D+0,UL(M,1,N) )
             ENDIF
           ENDIF
!
!---    x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.45 ) THEN 
           IF( PG(MP,N)-PL(MP,N).GT.EPSL ) THEN
             IF( PLB(MP,NB)-PGB(MP,NB).GT.EPSL ) THEN
               UL(M,1,N) = MAX( 0.D+0,UL(M,1,N) )
             ELSE
               UL(M,1,N) = MIN( 0.D+0,UL(M,1,N)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MP,NB)-PGB(MP,NB).LT.EPSL ) THEN
               UL(M,1,N) = MIN( 0.D+0,UL(M,1,N) )
             ENDIF
           ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MP,N)-PG(2,N).GE.EPSL ) THEN
            UL(M,1,N) = MIN( 0.D+0,UL(M,1,N) )
          ELSE
            UL(M,1,N) = 0.D+0
          ENDIF
!
!---    Dirichlet-Outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( UL(1,1,N).LT.-EPSL ) THEN
            UL(M,1,N) = UL(M,1,N)
          ELSE
            UL(M,1,N) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLW_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLE_W( N,NB )
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
!     Written by M.D. White, PNNL, 3 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLE_W'
      OECX = 0.D+0
      DO M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        HX = PL(MN,N)-PLB(MN,NB)
     &    -5.D-1*GRVX(2,N)*DXGF(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKL(MN,N),RKLB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(1,N)
        UL(M,2,N) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          UL(M,2,N) = MAX( 0.D+0,UL(M,2,N) )
!
!---    Seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 ) THEN 
!
!---       Node is unsaturated  ---
!
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
!
!---         Boundary is saturated, flow allowed into node  ---
!
             IF( PLB(MN,NB)-PGB(MN,NB).GT.EPSL ) THEN
               UL(M,2,N) = MIN( 0.D+0,UL(M,2,N) )
!
!---         Boundary is unsaturated, zero flow  ---
!
             ELSE
               UL(M,2,N) = 0.D+0
             ENDIF
!
!---       Node saturated  ---
!
           ELSE
!
!---         Boundary is unsaturated, flow allowed from node,
!            boundary is saturated, flow allowed into and from node  ---
!
             IF( PLB(MN,NB)-PGB(MN,NB).LT.EPSL ) THEN
               UL(M,2,N) = MAX( 0.D+0,UL(M,2,N) )
             ENDIF
           ENDIF
!
!---       Node is unsaturated, flow allowed into node  ---
!
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
             UL(M,2,N) = MIN( 0.D+0,UL(M,2,N) )
           ENDIF
!
!---    x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.45 ) THEN 
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MN,N)).GT.EPSL ) THEN
             IF( PLB(MN,NB)-PGB(MN,NB).GT.EPSL ) THEN
               UL(M,2,N) = MIN( 0.D+0,UL(M,2,N) )
             ELSE
               UL(M,2,N) = MAX( 0.D+0,UL(M,2,N)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MN,NB)-PGB(MN,NB).LT.EPSL ) THEN
               UL(M,2,N) = MAX( 0.D+0,UL(M,2,N) )
             ENDIF
           ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MN,N)-PG(2,N).GE.EPSL ) THEN
            UL(M,2,N) = MAX( 0.D+0,UL(M,2,N) )
          ELSE
            UL(M,2,N) = 0.D+0
          ENDIF
!
!---    Dirichlet-Outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( UL(1,2,N).GT.EPSL ) THEN
            UL(M,2,N) = UL(M,2,N)
          ELSE
            UL(M,2,N) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLE_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLN_W( N,NB )
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
!     Written by M.D. White, PNNL, 3 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLN_W'
      OECY = 0.D+0
      DO M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        HY = PL(MN,N)-PLB(MN,NB)
     &    -5.D-1*GRVY(2,N)*DYGF(N)*RP(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKL(MN,N),RKLB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(2,N)
        VL(M,2,N) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(N)*DYGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          VL(M,2,N) = MAX( 0.D+0,VL(M,2,N))
!
!---    Seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 ) THEN 
!
!---       Node is unsaturated  ---
!
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
!
!---         Boundary is saturated, flow allowed into node  ---
!
             IF( PLB(MN,NB)-PGB(MN,NB).GT.EPSL ) THEN
               VL(M,2,N) = MIN( 0.D+0,VL(M,2,N))
!
!---         Boundary is unsaturated, zero flow  ---
!
             ELSE
               VL(M,2,N) = 0.D+0
             ENDIF
!
!---       Node saturated  ---
!
           ELSE
!
!---         Boundary is unsaturated, flow allowed from node,
!            boundary is saturated, flow allowed into and from node  ---
!
             IF( PLB(MN,NB)-PGB(MN,NB).LT.EPSL ) THEN
               VL(M,2,N) = MAX( 0.D+0,VL(M,2,N))
             ENDIF
           ENDIF
!
!---       Node is unsaturated, flow allowed into node  ---
!
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
             VL(M,2,N) = MIN( 0.D+0,VL(M,2,N))
           ENDIF
!
!---    x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.45) THEN 
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MN,N)).GT.EPSL ) THEN
             IF( PLB(MN,NB)-PGB(MN,NB).GT.EPSL ) THEN
               VL(M,2,N) = MIN( 0.D+0,VL(M,2,N))
             ELSE
               VL(M,2,N) = MAX( 0.D+0,VL(M,2,N)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MN,NB)-PGB(MN,NB).LT.EPSL ) THEN
               VL(M,2,N) = MAX( 0.D+0,VL(M,2,N))
             ENDIF
           ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MN,N)-PG(2,N).GE.EPSL ) THEN
            VL(M,2,N) = MAX( 0.D+0,VL(M,2,N))
          ELSE
            VL(M,2,N) = 0.D+0
          ENDIF
!
!---    Dirichlet-Outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( VL(1,2,N).GT.EPSL ) THEN
            VL(M,2,N) = VL(M,2,N)
          ELSE
            VL(M,2,N) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLN_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVLT_W( N,NB )
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
!     Written by M.D. White, PNNL, 3 June 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DRCVLT_W'
      OECZ = 0.D+0
      DO M = 1,ISVF
        MN = MNEG(M)
        INDX = 11
        HZ = PL(MN,N)-PLB(MN,NB)
     &    -5.D-1*GRVZ(2,N)*DZGF(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKL(MN,N),RKLB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(3,N)
        WL(M,2,N) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow, diode or potential-evaporation
!       boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 .OR.
     &    ABS(IBCT(IEQW,NB)).EQ.23 ) THEN
          WL(M,2,N) = MAX( 0.D+0,WL(M,2,N) )
!
!---    Seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.17 ) THEN 
!
!---       Node is unsaturated  ---
!
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
!
!---         Boundary is saturated, flow allowed into node  ---
!
             IF( PLB(MN,NB)-PGB(MN,NB).GT.EPSL ) THEN
               WL(M,2,N) = MIN( 0.D+0,WL(M,2,N) )
!
!---         Boundary is unsaturated, zero flow  ---
!
             ELSE
               WL(M,2,N) = 0.D+0
             ENDIF
!
!---       Node saturated  ---
!
           ELSE
!
!---         Boundary is unsaturated, flow allowed from node,
!            boundary is saturated, flow allowed into and from node  ---
!
             IF( PLB(MN,NB)-PGB(MN,NB).LT.EPSL ) THEN
               WL(M,2,N) = MAX( 0.D+0,WL(M,2,N) )
             ENDIF
           ENDIF
!
!---    x-y-z seepage boundary condition  ---
!
        ELSEIF( ABS(IBCT(IEQW,NB)).EQ.45 ) THEN 
           IF( PG(MN,N)-PL(MN,N).GT.EPSL ) THEN
!           IF( ABS(1.D+0-SL(MN,N)).GT.EPSL ) THEN
             IF( PLB(MN,NB)-PGB(MN,NB).GT.EPSL ) THEN
               WL(M,2,N) = MIN( 0.D+0,WL(M,2,N) )
             ELSE
               WL(M,2,N) = MAX( 0.D+0,WL(M,2,N)*(RKM**4) )
             ENDIF
           ELSE
             IF( PLB(MN,NB)-PGB(MN,NB).LT.EPSL ) THEN
               WL(M,2,N) = MAX( 0.D+0,WL(M,2,N) )
             ENDIF
           ENDIF
!
!---    Shuttleworth-Wallace boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.22 ) THEN
          IF( PL(MN,N)-PG(2,N).GE.EPSL ) THEN
            WL(M,2,N) = MAX( 0.D+0,WL(M,2,N) )
          ELSE
            WL(M,2,N) = 0.D+0
          ENDIF
!
!---    Dirichlet-Outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQW,NB).EQ.26 ) THEN
          IF( WL(1,2,N).GT.EPSL ) THEN
            WL(M,2,N) = WL(M,2,N)
          ELSE
            WL(M,2,N) = 0.D+0
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLT_W group  ---
!
      RETURN
      END


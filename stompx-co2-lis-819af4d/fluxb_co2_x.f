!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAWB( N,NB )
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
!     Diffusive CO2 and water gas fluxes on bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAWB'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        DXMGW = XMGWB(2,NB) - XMGW(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
          FGWP = XGW(MP,N)*RHOG(MP,N)
          INDX = 3
          FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,1,N),INDX )
          FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
          FGAP = XGA(MP,N)*RHOG(MP,N)
          INDX = 3
          FGA = DIFMN( FGAB,FGAP,DZGF(N),DZGF(N),WG(1,1,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( WG(1,1,N).LT.-EPSL ) THEN
              WGW(M,1,N) = WG(M,1,N)*FGW
            ELSE
              WGW(M,1,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGW(MP,N)
     &        *RHOMG(MP,N)
            DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
     &        *RHOMGB(MP,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGW,INDX )
            WDGW(M,1,N) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))/DZGF(N)
            WGW(M,1,N) = WG(M,1,N)*FGW + WTMW*WDGW(M,1,N)
            WGA(M,1,N) = WG(M,1,N)*FGA - WTMA*WDGW(M,1,N)
          ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMGW = XMGWB(2,NB)*RHOMGB(2,NB) - XMGW(2,N)*RHOMG(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
          FGWP = XGW(MP,N)*RHOG(MP,N)
          INDX = 3
          FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,1,N),INDX )
          FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
          FGAP = XGA(MP,N)*RHOG(MP,N)
          INDX = 3
          FGA = DIFMN( FGAB,FGAP,DZGF(N),DZGF(N),WG(1,1,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( WG(1,1,N).LT.-EPSL ) THEN
              WGW(M,1,N) = WG(M,1,N)*FGW
            ELSE
              WGW(M,1,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGW(MP,N)
            DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGW,INDX )
            WDGW(M,1,N) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB) - 
     &        XMGW(MP,N)*RHOMG(MP,N))/DZGF(N)
            WGW(M,1,N) = WG(M,1,N)*FGW + WTMW*WDGW(M,1,N)
            WGA(M,1,N) = WG(M,1,N)*FGA - WTMA*WDGW(M,1,N)
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAWE( N,NB )
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
!     Diffusive CO2 and water gas fluxes on east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAWE'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        DXMGW = XMGW(2,N) - XMGWB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
          FGWP = XGW(MN,N)*RHOG(MN,N)
          INDX = 3
          FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,2,N),INDX )
          FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
          FGAP = XGA(MN,N)*RHOG(MN,N)
          INDX = 3
          FGA = DIFMN( FGAP,FGAB,DXGF(N),DXGF(N),UG(1,2,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( UG(1,2,N).GT.EPSL ) THEN
              UGW(M,2,N) = UG(M,2,N)*FGW
            ELSE
              UGW(M,2,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGW(MN,N)
     &        *RHOMG(MN,N)
            DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
     &        *RHOMGB(MN,NB)
            INDX = 12
            DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGW,INDX )
            UDGW(M,2,N) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))/DXGF(N)
            UGW(M,2,N) = UG(M,2,N)*FGW + WTMW*UDGW(M,2,N)
            UGA(M,2,N) = UG(M,2,N)*FGA - WTMA*UDGW(M,2,N)
          ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
          FGWP = XGW(MN,N)*RHOG(MN,N)
          INDX = 3
          FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,2,N),INDX )
          FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
          FGAP = XGA(MN,N)*RHOG(MN,N)
          INDX = 3
          FGA = DIFMN( FGAP,FGAB,DXGF(N),DXGF(N),UG(1,2,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( UG(1,2,N).GT.EPSL ) THEN
              UGW(M,2,N) = UG(M,2,N)*FGW
            ELSE
              UGW(M,2,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGW(MN,N)
            DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
            INDX = 12
            DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGW,INDX )
            UDGW(M,2,N) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N) - 
     &        XMGWB(MN,NB)*RHOMGB(MN,NB))/DXGF(N)
            UGW(M,2,N) = UG(M,2,N)*FGW + WTMW*UDGW(M,2,N)
            UGA(M,2,N) = UG(M,2,N)*FGA - WTMA*UDGW(M,2,N)
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAWN( N,NB )
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
!     Diffusive CO2 and water gas fluxes on north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAWN'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        DXMGW = XMGW(2,N) - XMGWB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
          FGWP = XGW(MN,N)*RHOG(MN,N)
          INDX = 3
          FGW = DIFMN( FGWP,FGWB,DYGF(N),DYGF(N),VG(1,2,N),INDX )
          FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
          FGAP = XGA(MN,N)*RHOG(MN,N)
          INDX = 3
          FGA = DIFMN( FGAP,FGAB,DYGF(N),DYGF(N),VG(1,2,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( VG(1,2,N).GT.EPSL ) THEN
              VGW(M,2,N) = VG(M,2,N)*FGW
            ELSE
              VGW(M,2,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGW(MN,N)
     &        *RHOMG(MN,N)
            DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
     &        *RHOMGB(MN,NB)
            INDX = 12
            DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGW,INDX )
            VDGW(M,2,N) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))
     &        /(DYGF(N)*RP(N))
            VGW(M,2,N) = VG(M,2,N)*FGW + WTMW*VDGW(M,2,N)
            VGA(M,2,N) = VG(M,2,N)*FGA - WTMA*VDGW(M,2,N)
          ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
          FGWP = XGW(MN,N)*RHOG(MN,N)
          INDX = 3
          FGW = DIFMN( FGWP,FGWB,DYGF(N),DYGF(N),VG(1,2,N),INDX )
          FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
          FGAP = XGA(MN,N)*RHOG(MN,N)
          INDX = 3
          FGA = DIFMN( FGAP,FGAB,DYGF(N),DYGF(N),VG(1,2,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( VG(1,2,N).GT.EPSL ) THEN
              VGW(M,2,N) = VG(M,2,N)*FGW
            ELSE
              VGW(M,2,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGW(MN,N)
            DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
            INDX = 12
            DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGW,INDX )
            VDGW(M,2,N) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N) - 
     &        XMGWB(MN,NB)*RHOMGB(MN,NB))/(DYGF(N)*RP(N))
            VGW(M,2,N) = VG(M,2,N)*FGW + WTMW*VDGW(M,2,N)
            VGA(M,2,N) = VG(M,2,N)*FGA - WTMA*VDGW(M,2,N)
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAWS( N,NB )
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
!     Diffusive CO2 and water gas fluxes on south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAWS'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        DXMGW = XMGWB(2,NB) - XMGW(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
          FGWP = XGW(MP,N)*RHOG(MP,N)
          INDX = 3
          FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,1,N),INDX )
          FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
          FGAP = XGA(MP,N)*RHOG(MP,N)
          INDX = 3
          FGA = DIFMN( FGAB,FGAP,DYGF(N),DYGF(N),VG(1,1,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( VG(1,1,N).LT.-EPSL ) THEN
              VGW(M,1,N) = VG(M,1,N)*FGW
            ELSE
              VGW(M,1,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGW(MP,N)
     &        *RHOMG(MP,N)
            DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
     &        *RHOMGB(MP,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGW,INDX )
            VDGW(M,1,N) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))
     &        /(DYGF(N)*RP(N))
            VGW(M,1,N) = VG(M,1,N)*FGW + WTMW*VDGW(M,1,N)
            VGA(M,1,N) = VG(M,1,N)*FGA - WTMA*VDGW(M,1,N)
          ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMGW = XMGWB(2,NB)*RHOMGB(2,NB) - XMGW(2,N)*RHOMG(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
          FGWP = XGW(MP,N)*RHOG(MP,N)
          INDX = 3
          FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,1,N),INDX )
          FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
          FGAP = XGA(MP,N)*RHOG(MP,N)
          INDX = 3
          FGA = DIFMN( FGAB,FGAP,DYGF(N),DYGF(N),VG(1,1,N),INDX )
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( VG(1,1,N).LT.-EPSL ) THEN
              VGW(M,1,N) = VG(M,1,N)*FGW
            ELSE
              VGW(M,1,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGW(MP,N)
            DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGW,INDX )
            VDGW(M,1,N) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB) - 
     &        XMGW(MP,N)*RHOMG(MP,N))/(DYGF(N)*RP(N))
            VGW(M,1,N) = VG(M,1,N)*FGW + WTMW*VDGW(M,1,N)
            VGA(M,1,N) = VG(M,1,N)*FGA - WTMA*VDGW(M,1,N)
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAWT( N,NB )
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
!     Diffusive CO2 and water gas fluxes on top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAWT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        DXMGW = XMGW(2,N)-XMGWB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
          FGWP = XGW(MN,N)*RHOG(MN,N)
          INDX = 3
          FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,2,N),INDX )
          FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
          FGAP = XGA(MN,N)*RHOG(MN,N)
          INDX = 3
          FGA = DIFMN( FGAP,FGAB,DZGF(N),DZGF(N),WG(1,2,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( WG(1,2,N).GT.EPSL ) THEN
              WGW(M,2,N) = WG(M,2,N)*FGW
            ELSE
              WGW(M,2,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGW(MN,N)
     &        *RHOMG(MN,N)
            DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
     &        *RHOMGB(MN,NB)
            INDX = 12
            DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGW,INDX )
            WDGW(M,2,N) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))/DZGF(N)
            WGW(M,2,N) = WG(M,2,N)*FGW + WTMW*WDGW(M,2,N)
            WGA(M,2,N) = WG(M,2,N)*FGA - WTMA*WDGW(M,2,N)
          ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMGW = XMGW(2,N)*RHOMG(2,N) - XMGWB(2,NB)*RHOMGB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
          FGWP = XGW(MN,N)*RHOG(MN,N)
          INDX = 3
          FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,2,N),INDX )
          FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
          FGAP = XGA(MN,N)*RHOG(MN,N)
          INDX = 3
          FGA = DIFMN( FGAP,FGAB,DZGF(N),DZGF(N),WG(1,2,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( WG(1,2,N).GT.EPSL ) THEN
              WGW(M,2,N) = WG(M,2,N)*FGW
            ELSE
              WGW(M,2,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGW(MN,N)
            DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
            INDX = 12
            DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGW,INDX )
            WDGW(M,2,N) = 2.D+0*DFM*(XMGW(MN,N)*RHOMG(MN,N) - 
     &        XMGWB(MN,NB)*RHOMGB(MN,NB))/DZGF(N)
            WGW(M,2,N) = WG(M,2,N)*FGW + WTMW*WDGW(M,2,N)
            WGA(M,2,N) = WG(M,2,N)*FGA - WTMA*WDGW(M,2,N)
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAWW( N,NB )
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
!     Diffusive CO2 and water gas fluxes on west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGAWW'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
        DXMGW = XMGWB(2,NB)-XMGW(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
          FGWP = XGW(MP,N)*RHOG(MP,N)
          INDX = 3
          FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,1,N),INDX )
          FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
          FGAP = XGA(MP,N)*RHOG(MP,N)
          INDX = 3
          FGA = DIFMN( FGAB,FGAP,DXGF(N),DXGF(N),UG(1,1,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( UG(1,1,N).LT.-EPSL ) THEN
              UGW(M,1,N) = UG(M,1,N)*FGW
            ELSE
              UGW(M,1,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGW(MP,N)
     &        *RHOMG(MP,N)
            DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
     &        *RHOMGB(MP,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGW,INDX )
            UDGW(M,1,N) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))/DXGF(N)
            UGW(M,1,N) = UG(M,1,N)*FGW + WTMW*UDGW(M,1,N)
            UGA(M,1,N) = UG(M,1,N)*FGA - WTMA*UDGW(M,1,N)
          ENDIF
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMGW = XMGWB(2,NB)*RHOMGB(2,NB)-XMGW(2,N)*RHOMG(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
          FGWP = XGW(MP,N)*RHOG(MP,N)
          INDX = 3
          FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,1,N),INDX )
          FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
          FGAP = XGA(MP,N)*RHOG(MP,N)
          INDX = 3
          FGA = DIFMN( FGAB,FGAP,DXGF(N),DXGF(N),UG(1,1,N),INDX )
!
!---      Dirichlet-Outflow boundary condition  ---
!
          IF( IBCT(IEQA,NB).EQ.19 ) THEN
            IF( UG(1,1,N).LT.-EPSL ) THEN
              UGW(M,1,N) = UG(M,1,N)*FGW
            ELSE
              UGW(M,1,N) = 0.D+0
            ENDIF
          ELSE
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGW(MP,N)
            DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGW,INDX )
            UDGW(M,1,N) = 2.D+0*DFM*(XMGWB(MP,NB)*RHOMGB(MP,NB) - 
     &        XMGW(MP,N)*RHOMG(MP,N))/DXGF(N)
            UGW(M,1,N) = UG(M,1,N)*FGW + WTMW*UDGW(M,1,N)
            UGA(M,1,N) = UG(M,1,N)*FGA - WTMA*UDGW(M,1,N)
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAWW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAWB( N,NB )
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
!     Diffusive CO2 and water aqueous fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAWB'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        DXMLA = XMLAB(2,NB) - XMLA(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,1,N),INDX )
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,1,N),INDX )
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,1,N) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/DZGF(N)
          WLA(M,1,N) = WL(M,1,N)*FLA + WTMA*WDLA(M,1,N)
          WLW(M,1,N) = WL(M,1,N)*FLW - 
     &      WTMW*(WDLA(M,1,N)+WDS(M,1,N)/WTMS)
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,1,N),INDX )
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,1,N),INDX )
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,1,N) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB) - 
     &      XMLA(MP,N)*RHOML(MP,N))/DZGF(N)
          WLA(M,1,N) = WL(M,1,N)*FLA + WTMA*WDLA(M,1,N)
          WLW(M,1,N) = WL(M,1,N)*FLW - 
     &      WTMW*(WDLA(M,1,N)+WDS(M,1,N)/WTMS)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAWE( N,NB )
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
!     Diffusive CO2 and water aqueous flux on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAWE'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        DXMLA = XMLA(2,N) - XMLAB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,2,N),INDX )
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,2,N),INDX )
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,2,N) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))/DXGF(N)
          ULA(M,2,N) = UL(M,2,N)*FLA + WTMA*UDLA(M,2,N)
          ULW(M,2,N) = UL(M,2,N)*FLW - 
     &      WTMW*(UDLA(M,2,N)+UDS(M,2,N)/WTMS)
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,2,N),INDX )
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,2,N),INDX )
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,2,N) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N) -
     &      XMLAB(MN,NB)*RHOMLB(MN,NB))/DXGF(N)
          ULA(M,2,N) = UL(M,2,N)*FLA + WTMA*UDLA(M,2,N)
          ULW(M,2,N) = UL(M,2,N)*FLW - 
     &      WTMW*(UDLA(M,2,N)+UDS(M,2,N)/WTMS)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAWN( N,NB )
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
!     Diffusive CO2 and water aqueous fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAWN'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        DXMLA = XMLA(2,N) - XMLAB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,2,N),INDX )
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,2,N),INDX )
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,2,N) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))
     &      /(DYGF(N)*RP(N))
          VLA(M,2,N) = VL(M,2,N)*FLA + WTMA*VDLA(M,2,N)
          VLW(M,2,N) = VL(M,2,N)*FLW - 
     &      WTMW*(VDLA(M,2,N)+VDS(M,2,N)/WTMS)
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,2,N),INDX )
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,2,N),INDX )
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,2,N) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N) -
     &      XMLAB(MN,NB)*RHOMLB(MN,NB))/(DYGF(N)*RP(N))
          VLA(M,2,N) = VL(M,2,N)*FLA + WTMA*VDLA(M,2,N)
          VLW(M,2,N) = VL(M,2,N)*FLW - 
     &      WTMW*(VDLA(M,2,N)+VDS(M,2,N)/WTMS)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAWS( N,NB )
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
!     Diffusive CO2 and water aqueous fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAWS'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        DXMLA = XMLAB(2,NB) - XMLA(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,1,N),INDX )
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,1,N),INDX )
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,1,N) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))
     &      /(DYGF(N)*RP(N))
          VLA(M,1,N) = VL(M,1,N)*FLA + WTMA*VDLA(M,1,N)
          VLW(M,1,N) = VL(M,1,N)*FLW - 
     &      WTMW*(VDLA(M,1,N)+VDS(M,1,N)/WTMS)
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,1,N),INDX )
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,1,N),INDX )
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,1,N) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB) - 
     &      XMLA(MP,N)*RHOML(MP,N))/(DYGF(N)*RP(N))
          VLA(M,1,N) = VL(M,1,N)*FLA + WTMA*VDLA(M,1,N)
          VLW(M,1,N) = VL(M,1,N)*FLW - 
     &      WTMW*(VDLA(M,1,N)+VDS(M,1,N)/WTMS)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAWT( N,NB )
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
!     Diffusive CO2 and water aqueous fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAWT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        DXMLA = XMLA(2,N)-XMLAB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,2,N),INDX )
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,2,N),INDX )
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,2,N) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))/DZGF(N)
          WLA(M,2,N) = WL(M,2,N)*FLA + WTMA*WDLA(M,2,N)
          WLW(M,2,N) = WL(M,2,N)*FLW - 
     &      WTMW*(WDLA(M,2,N)+WDS(M,2,N)/WTMS)
        ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
        DO M = 1,ISVF
          MN = MNEG(M)
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,2,N),INDX )
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
          INDX = 2
          FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,2,N),INDX )
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,2,N) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N) - 
     &      XMLAB(MN,NB)*RHOMLB(MN,NB))/DZGF(N)
          WLA(M,2,N) = WL(M,2,N)*FLA + WTMA*WDLA(M,2,N)
          WLW(M,2,N) = WL(M,2,N)*FLW - 
     &      WTMW*(WDLA(M,2,N)-WDS(M,2,N)/WTMS)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAWW( N,NB )
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
!     Diffusive CO2 and water aqueous fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLAWW'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
        DXMLA = XMLAB(2,NB)-XMLA(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,1,N),INDX )
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,1,N),INDX )
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,1,N) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/DXGF(N)
          ULA(M,1,N) = UL(M,1,N)*FLA + WTMA*UDLA(M,1,N)
          ULW(M,1,N) = UL(M,1,N)*FLW - 
     &      WTMW*(UDLA(M,1,N)+UDS(M,1,N)/WTMS)
      ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
        DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
        DO M = 1,ISVF
          MP = MPOS(M)
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,1,N),INDX )
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
          INDX = 2
          FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,1,N),INDX )
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,1,N) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB) - 
     &      XMLA(MP,N)*RHOML(MP,N))/DXGF(N)
          ULA(M,1,N) = UL(M,1,N)*FLA + WTMA*UDLA(M,1,N)
          ULW(M,1,N) = UL(M,1,N)*FLW - 
     &      WTMW*(UDLA(M,1,N)+UDS(M,1,N)/WTMS)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAWW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSB( N,NB )
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
!     Diffusive salt aqueous fluxes on bottom boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSB'
      DO M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MP,N)+TABS)/TSPRF
          SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
          TCOR = (TB(MP,NB)+TABS)/TSPRF
          SMDLB = DFLSB(MP,NB)*TCOR*(VISRL/VISLB(MP,NB))
          DFFLB = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*SMDLB
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
          DFFLB = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLSB(MP,NB)
        ENDIF
        INDX = 18
        DFFLB = DIFMN(DFFLB,DFFLP,DZGF(N),DZGF(N),WL(1,1,N),INDX)
!
!---    Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
          DDLB = DFFLB/(5.D-1*DZGF(N))
          AL = MAX( WL(M,1,N),ZERO ) +
     &     DDLB*MAX((ONE-(TENTH*ABS(WL(M,1,N))/(DDLB+SMALL)))**5,ZERO)
          ALP = MAX( -WL(M,1,N),ZERO ) +
     &     DDLB*MAX((ONE-(TENTH*ABS(WL(M,1,N))/(DDLB+SMALL)))**5,ZERO)
          WS(M,1,N) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL
     &     - XLS(MP,N)*RHOL(MP,N)*ALP)
          WDS(M,1,N) = DDLB*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &      XLS(MP,N)*RHOL(MP,N))
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( -WL(M,1,N),ZERO )
          WS(M,1,N) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          AL = MAX( WL(M,1,N),ZERO )
          WS(M,1,N) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSS( N,NB )
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
!     Diffusive salt aqueous fluxes on south boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSS'
      DO M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MP,N)+TABS)/TSPRF
          SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
          TCOR = (TB(MP,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MP,NB)*TCOR*(VISRL/VISLB(MP,NB))
          DFFLS = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*SMDLP
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
          DFFLS = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLSB(MP,NB)
        ENDIF
        INDX = 18
        DFFLS = DIFMN(DFFLS,DFFLP,DYGF(N),DYGF(N),VL(1,1,N),INDX)
!
!---    Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
          DDLS = DFFLS/RP(N)/(5.D-1*DYGF(N))
          AL = MAX( VL(M,1,N),ZERO ) +
     &     DDLS*MAX((ONE-(TENTH*ABS(VL(M,1,N))/(DDLS+SMALL)))**5,ZERO)
          ALP = MAX( -VL(M,1,N),ZERO ) +
     &     DDLS*MAX((ONE-(TENTH*ABS(VL(M,1,N))/(DDLS+SMALL)))**5,ZERO)
          VS(M,1,N) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL
     &      - XLS(MP,N)*RHOL(MP,N)*ALP)
          VDS(M,1,N) = DDLS*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &      XLS(MP,N)*RHOL(MP,N))
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( -VL(M,1,N),ZERO )
          VS(M,1,N) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          AL = MAX( VL(M,1,N),ZERO )
          VS(M,1,N) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSW( N,NB )
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
!     Diffusive salt aqueous fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSW'
      DO M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MP,N)+TABS)/TSPRF
          SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
          TCOR = (TB(MP,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MP,NB)*TCOR*(VISRL/VISLB(MP,NB))
          DFFLW = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*SMDLP
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
          DFFLW = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLSB(MP,NB)
        ENDIF
        INDX = 18
        DFFLW = DIFMN(DFFLW,DFFLP,DXGF(N),DXGF(N),UL(1,1,N),INDX)
!
!---    Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
          DDLW = DFFLW/(5.D-1*DXGF(N))
          AL = MAX( UL(M,1,N),ZERO ) +
     &     DDLW*MAX((ONE-(TENTH*ABS(UL(M,1,N))/(DDLW+SMALL)))**5,ZERO)
          ALP = MAX( -UL(M,1,N),ZERO ) +
     &     DDLW*MAX((ONE-(TENTH*ABS(UL(M,1,N))/(DDLW+SMALL)))**5,ZERO)
          US(M,1,N) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL
     &      - XLS(MP,N)*RHOL(MP,N)*ALP)
          UDS(M,1,N) = DDLW*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &      XLS(MP,N)*RHOL(MP,N))
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( -UL(M,1,N),ZERO )
          US(M,1,N) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          AL = MAX( UL(M,1,N),ZERO )
          US(M,1,N) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSE( N,NB )
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
!     Diffusive salt aqueous fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSE'
      DO M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MN,N)+TABS)/TSPRF
          SMDLP = DFLS(MN,N)*TCOR*(VISRL/VISL(MN,N))
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*SMDLP
          TCOR = (TB(MN,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MN,NB)*TCOR*(VISRL/VISLB(MN,NB))
          DFFLE = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*SMDLP
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLS(MN,N)
          DFFLE = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLSB(MN,NB)
        ENDIF
        INDX = 18
        DFFLE = DIFMN(DFFLP,DFFLE,DXGF(N),DXGF(N),UL(1,2,N),INDX)
!
!---    Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
          DDLE = DFFLE/(5.D-1*DXGF(N))
          AL = MAX( -UL(M,2,N),ZERO ) +
     &     DDLE*MAX((ONE-(TENTH*ABS(UL(M,2,N))/(DDLE+SMALL)))**5,ZERO)
          ALP = MAX( UL(M,2,N),ZERO ) +
     &     DDLE*MAX((ONE-(TENTH*ABS(UL(M,2,N))/(DDLE+SMALL)))**5,ZERO)
          US(M,2,N) = (XLS(MN,N)*RHOL(MN,N)*ALP
     &     - XLSB(MN,NB)*RHOLB(MN,NB)*AL)
          UDS(M,2,N) = DDLE*(XLS(MN,N)*RHOL(MN,N) -
     &     XLSB(MN,NB)*RHOLB(MN,NB))
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( UL(M,2,N),ZERO )
          US(M,2,N) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          AL = MAX( -UL(M,2,N),ZERO )
          US(M,2,N) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSN( N,NB )
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
!     Diffusive salt aqueous fluxes on north boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DFFLSN'
      DO M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MN,N)+TABS)/TSPRF
          SMDLP = DFLS(MN,N)*TCOR*(VISRL/VISL(MN,N))
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*SMDLP
          TCOR = (TB(MN,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MN,NB)*TCOR*(VISRL/VISLB(MN,NB))
          DFFLN = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*SMDLP
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLS(MN,N)
          DFFLN = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLSB(MN,NB)
        ENDIF
        INDX = 18
        DFFLN = DIFMN(DFFLP,DFFLN,DYGF(N),DYGF(N),VL(1,2,N),INDX)
!
!---    Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
          DDLN = DFFLN/RP(N)/(5.D-1*DYGF(N))
          AL = MAX( -VL(M,2,N),ZERO ) +
     &     DDLN*MAX((ONE-(TENTH*ABS(VL(M,2,N))/(DDLN+SMALL)))**5,ZERO)
          ALP = MAX( VL(M,2,N),ZERO ) +
     &     DDLN*MAX((ONE-(TENTH*ABS(VL(M,2,N))/(DDLN+SMALL)))**5,ZERO)
          VS(M,2,N) = (XLS(MN,N)*RHOL(MN,N)*ALP -
     &      XLSB(MN,NB)*RHOLB(MN,NB)*AL)
          VDS(M,2,N) = DDLN*(XLS(MN,N)*RHOL(MN,N) -
     &      XLSB(MN,NB)*RHOLB(MN,NB))
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( VL(M,2,N),ZERO )
          VS(M,2,N) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          AL = MAX( -VL(M,2,N),ZERO )
          VS(M,2,N) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLST( N,NB )
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
!     Diffusive salt aqueous fluxes on top boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DFFLST'
      DO M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MN,N)+TABS)/TSPRF
          SMDLP = DFLS(MN,N)*TCOR*(VISRL/VISL(MN,N))
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*SMDLP
          TCOR = (TB(MN,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MN,NB)*TCOR*(VISRL/VISLB(MN,NB))
          DFFLT = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*SMDLP
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLS(MN,N)
          DFFLT = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLSB(MN,NB)
        ENDIF
        INDX = 18
        DFFLT = DIFMN(DFFLP,DFFLT,DZGF(N),DZGF(N),WL(1,2,N),INDX)
!
!---    Dirichlet boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.1 .OR. IBCT(IEQS,NB).EQ.8 .OR.
     &    IBCT(IEQS,NB).EQ.12 .OR. (IBCT(IEQS,NB).GE.34 .AND.
     &    IBCT(IEQS,NB).LE.41) ) THEN
          DDLT = DFFLT/(5.D-1*DZGF(N))
          AL = MAX( -WL(M,2,N),ZERO ) +
     &     DDLT*MAX((ONE-(TENTH*ABS(WL(M,2,N))/(DDLT+SMALL)))**5,ZERO)
          ALP = MAX( WL(M,2,N),ZERO ) +
     &     DDLT*MAX((ONE-(TENTH*ABS(WL(M,2,N))/(DDLT+SMALL)))**5,ZERO)
          WS(M,2,N) = (XLS(MN,N)*RHOL(MN,N)*ALP
     &      - XLSB(MN,NB)*RHOLB(MN,NB)*AL)
          WDS(M,2,N) = DDLT*(XLS(MN,N)*RHOL(MN,N) -
     &      XLSB(MN,NB)*RHOLB(MN,NB))
!
!---   Outflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( WL(M,2,N),ZERO )
          WS(M,2,N) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.13 .OR. IBCT(IEQS,NB).EQ.14 ) THEN
          AL = MAX( -WL(M,2,N),ZERO )
          WS(M,2,N) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLST group  ---
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGB'
      DO M = 1,ISVF
        MP = MPOS(M)
        HZ = PGB(MP,NB)-PG(MP,N)-5.D-1*GRVZ(1,N)
     &    *DZGF(N)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(3,N)
        WG(M,1,N) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          WG(M,1,N) = MIN( 0.D+0,WG(M,1,N) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGS'
      DO M = 1,ISVF
        MP = MPOS(M)
        HY = PGB(MP,NB)-PG(MP,N)-5.D-1*GRVY(1,N)
     &    *DYGF(N)*RP(N)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(2,N)
        VG(M,1,N) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(N)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          VG(M,1,N) = MIN( 0.D+0,VG(M,1,N) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGW'
      DO M = 1,ISVF
        MP = MPOS(M)
        HX = PGB(MP,NB)-PG(MP,N)-5.D-1*GRVX(1,N)
     &    *DXGF(N)*RHOGB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 9
        RKM = DIFMN(RKGB(MP,NB),RKG(MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISGB(MP,NB),VISG(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(1,N)
        UG(M,1,N) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          UG(M,1,N) = MIN( 0.D+0,UG(M,1,N) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGE'
      DO M = 1,ISVF
        MN = MNEG(M)
        HX = PG(MN,N)-PGB(MN,NB)-5.D-1*GRVX(2,N)
     &    *DXGF(N)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(1,N)
        UG(M,2,N) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          UG(M,2,N) = MAX( 0.D+0,UG(M,2,N) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGN'
      DO M = 1,ISVF
        MN = MNEG(M)
        HY = PG(MN,N)-PGB(MN,NB)-5.D-1*GRVY(2,N)
     &    *DYGF(N)*RP(N)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(2,N)
        VG(M,2,N) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(N)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          VG(M,2,N) = MAX( 0.D+0,VG(M,2,N))
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVGT'
      DO M = 1,ISVF
        MN = MNEG(M)
        HZ = PG(MN,N)-PGB(MN,NB)-5.D-1*GRVZ(2,N)
     &    *DZGF(N)*RHOGB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 9
        RKM = DIFMN(RKG(MN,N),RKGB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 6
        VM = DIFMN(VISG(MN,N),VISGB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(3,N)
        WG(M,2,N) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          WG(M,2,N) = MAX( 0.D+0,WG(M,2,N) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLB'
      DO M = 1,ISVF
        MP = MPOS(M)
        HZ = PLB(MP,NB)-PL(MP,N)-5.D-1*GRVZ(1,N)
     &    *DZGF(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKLB(MP,NB),RKL(MP,N),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(3,N)
        WL(M,1,N) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          WL(M,1,N) = MIN( 0.D+0,WL(M,1,N) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLS'
      DO M = 1,ISVF
        MP = MPOS(M)
        HY = PLB(MP,NB)-PL(MP,N)-5.D-1*GRVY(1,N)
     &    *DYGF(N)*RP(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKLB(MP,NB),RKL(MP,N),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(2,N)
        VL(M,1,N) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(N)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          VL(M,1,N) = MIN( 0.D+0,VL(M,1,N) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLW'
      DO M = 1,ISVF
        MP = MPOS(M)
        HX = PLB(MP,NB)-PL(MP,N)-5.D-1*GRVX(1,N)
     &    *DXGF(N)*RHOLB(MP,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKLB(MP,NB),RKL(MP,N),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISLB(MP,NB),VISL(MP,N),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MP,N)*PERM(1,N)
        UL(M,1,N) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          UL(M,1,N) = MIN( 0.D+0,UL(M,1,N) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLE'
      DO M = 1,ISVF
        MN = MNEG(M)
        HX = PL(MN,N)-PLB(MN,NB)-5.D-1*GRVX(2,N)
     &    *DXGF(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HX
        INDX = 8
        RKM = DIFMN(RKL(MN,N),RKLB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DXGF(N),DXGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(1,N)
        UL(M,2,N) = 2.D+0*PERM_PX*RKM*HX/(VM*DXGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          UL(M,2,N) = MAX( 0.D+0,UL(M,2,N) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLN'
      DO M = 1,ISVF
        MN = MNEG(M)
        HY = PL(MN,N)-PLB(MN,NB)-5.D-1*GRVY(2,N)
     &    *DYGF(N)*RP(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HY
        INDX = 8
        RKM = DIFMN(RKL(MN,N),RKLB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DYGF(N),DYGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(2,N)
        VL(M,2,N) = 2.D+0*PERM_PX*RKM*HY/(VM*RP(N)*DYGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          VL(M,2,N) = MAX( 0.D+0,VL(M,2,N))
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
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
!     Written by M.D. White, PNNL, 2 November 2021.
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
      SUB_LOG(ISUB_LOG) = '/DRCVLT'
      DO M = 1,ISVF
        MN = MNEG(M)
        HZ = PL(MN,N)-PLB(MN,NB)-5.D-1*GRVZ(2,N)
     &    *DZGF(N)*RHOLB(MN,NB)
        IF( M.EQ.1 ) HD = HZ
        INDX = 8
        RKM = DIFMN(RKL(MN,N),RKLB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        INDX = 5
        VM = DIFMN(VISL(MN,N),VISLB(MN,NB),DZGF(N),DZGF(N),HD,INDX)
        PERM_PX = PERMRF(MN,N)*PERM(3,N)
        WL(M,2,N) = 2.D+0*PERM_PX*RKM*HZ/(VM*DZGF(N))
!
!---    Outflow boundary condition  ---
!
        IF( ABS(IBCT(IEQW,NB)).EQ.7 ) THEN
          WL(M,2,N) = MAX( 0.D+0,WL(M,2,N) )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVLT group  ---
!
      RETURN
      END


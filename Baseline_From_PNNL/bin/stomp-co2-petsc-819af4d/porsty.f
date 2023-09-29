!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY( N,PX,PREFX,PORDX,PORTX )
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
!     Compute diffusive and total porosities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, December, 1992.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
      USE GEO_MECH
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PORSTY'
      IZN = IZ(N)
      DPX = PX-PREFX
      ISLC50X = ABS(ISLC(50))
!
!---  Geomechanics simulations  --
!
      IF( ISLC(50).NE.0 .AND. MOD(ISLC50X,10).NE.2 .AND. 
     &  MOD(ISLC50X,10).NE.4 ) THEN
!
!---    Drained bulk modulus  --
!
        BLKDRNX = PROP_GM(1,IZN)/(3.0D+0*(1.0D+0-2.0D+0*PROP_GM(2,IZN)))
!
!---    1/N  ---
!
        OONMODX = (PROP_GM(3,IZN)-POR(1,IZN))*(1.0D+0-PROP_GM(3,IZN))
     &    /BLKDRNX
!
!---    Volumetric strain differential, at iterate level k  ---
!
        DEVKX = EPSV_GM(2,N)-EPSV_CMP(N)
!
!---    Pressure differential, at iterate level k  ---
!
        DPKX = P_GM(2,N)-PCMP(N)
!
!---    Diffusive and total mechanical porosity, at iterate level k  ---
!
        PORDMCHKX = POR(1,IZN) - PROP_GM(3,IZN)*DEVKX + OONMODX*DPKX
        PORTMCHKX = POR(2,IZN) - PROP_GM(3,IZN)*DEVKX + OONMODX*DPKX
!
!---    Pressure differential, at iterate level k+1  ---
!
        DPK1X = PX-P_GM(2,N)
!
!---    Diffusive and total flow porosity, at iterate level k+1  ---
!
        PORDX = PORDMCHKX + ((PROP_GM(3,IZN)**2)/BLKDRNX+OONMODX)*DPK1X
        PORTX = PORTMCHKX + ((PROP_GM(3,IZN)**2)/BLKDRNX+OONMODX)*DPK1X
!
!---  No geomechanics  --
!
      ELSE        
!
!---    Pore compressibility w/ fixed bulk volume  ---
!
        IF( ISLC(15).EQ.1 ) THEN
!
!---      Fracture properties (dual porosity model)  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
            PORTX = POR(3,IZN)*EXP(DPX*CMP(2,IZN))
            PORDX = POR(4,IZN)*EXP(DPX*CMP(2,IZN))
!
!---      Matrix properties (dual porosity model)  ---
!
          ELSEIF( ABS(IDP(IZN)).EQ.2 ) THEN
            PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
            PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
            ENDIF
          ENDIF
!
!---    Bulk compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.10 ) THEN
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(1,N)/
     &          POR0(1,N)))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(2,N)/
     &          POR0(2,N)))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)/
     &          POR(1,IZN)))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &          POR(2,IZN)))
            ENDIF
!
!---    Pore compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.11 ) THEN
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(1,N)))
            PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(2,N)))
          ELSE
            PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)))
            PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
          ENDIF
!
!---    Bulk compressibility w/ fixed bulk volume  ---
!
        ELSE
!
!---      Fracture properties (dual porosity model)  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
            PORTX = POR(3,IZN) + DPX*CMP(2,IZN)
            PORDX = POR(4,IZN) + DPX*CMP(2,IZN)
!
!---      Matrix properties (dual porosity model)  ---
!
          ELSEIF( ABS(IDP(IZN)).EQ.2 ) THEN
            PORTX = POR(1,IZN) + DPX*CMP(1,IZN)
            PORDX = POR(2,IZN) + DPX*CMP(1,IZN)
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N) + DPX*CMP(1,IZN)
              PORDX = POR0(2,N) + DPX*CMP(1,IZN)
            ELSE
              PORTX = POR(1,IZN) + DPX*CMP(1,IZN)
              PORDX = POR(2,IZN) + DPX*CMP(1,IZN)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      PORTX = MAX( MIN( PORTX,1.D+0 ),1.D-6 )
      PORDX = MAX( MIN( PORDX,1.D+0 ),1.D-6 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY_BH( NBNX,PX,PREFX,PORDX )
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
!     Compute diffusive and total porosities for borehole nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 06 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GEOM_BH
      USE FDVP_BH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PORSTY_BH'
      DPX = PX-PREFX
      IZN = IZ_BH(NBNX)
!
!---  Pore compressibility w/ fixed bulk volume  ---
!
      IF( ISLC(15).EQ.1 ) THEN
!
!---    Dual porosity model  ---
!
        IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---      Fracture properties (dual porosity model)  ---
!
          PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN))
!
!---      Matrix properties (dual porosity model)  ---
!
          PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
!
!---      Effective properties (dual porosity model)  ---
!
          PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---    Single porosity model  ---
!
        ELSE
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORDX = POR0_BH(2,NBNX)*EXP(DPX*CMP(1,IZN))
          ELSE
            PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
          ENDIF
        ENDIF
!
!---  Bulk compressibility w/ variable bulk volume  ---
!
      ELSEIF( ISLC(15).EQ.10 ) THEN
!
!---    Dual porosity model  ---
!
        IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---      Fracture properties (dual porosity model)  ---
!
          PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(4,IZN)/
     &      POR(4,IZN)))
!
!---      Matrix properties (dual porosity model)  ---
!
          PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &      POR(2,IZN)))
!
!---      Effective properties (dual porosity model)  ---
!
          PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---    Single porosity model  ---
!
        ELSE
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORDX = POR0_BH(2,NBNX)*EXP(DPX*CMP(1,IZN)*
     &        (1.D+0-POR0_BH(2,NBNX)/POR0_BH(2,NBNX)))
          ELSE
            PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &        POR(2,IZN)))
          ENDIF
        ENDIF
!
!---  Pore compressibility w/ variable bulk volume  ---
!
      ELSEIF( ISLC(15).EQ.11 ) THEN
!
!---    Dual porosity model  ---
!
        IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---      Fracture properties (dual porosity model)  ---
!
          PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(4,IZN)))
!
!---      Matrix properties (dual porosity model)  ---
!
          PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
!
!---      Effective properties (dual porosity model)  ---
!
          PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---    Single porosity model  ---
!
        ELSE
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORDX = POR0_BH(2,NBNX)*EXP(DPX*CMP(1,IZN)*
     &        (1.D+0-POR0_BH(2,NBNX)))
          ELSE
            PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
          ENDIF
        ENDIF
!
!---  Bulk compressibility w/ fixed bulk volume  ---
!
      ELSE
!
!---    Dual porosity model  ---
!
        IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---      Fracture properties (dual porosity model)  ---
!
          PORDFX = POR(4,IZN) + DPX*CMP(2,IZN)
!
!---      Matrix properties (dual porosity model)  ---
!
          PORDMX = POR(2,IZN) + DPX*CMP(1,IZN)
!
!---      Effective properties (dual porosity model)  ---
!
          PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---    Single porosity model  ---
!
        ELSE
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORDX = POR0_BH(2,NBNX) + DPX*CMP(1,IZN)
          ELSE
            PORDX = POR(2,IZN) + DPX*CMP(1,IZN)
          ENDIF
        ENDIF
      ENDIF
      PORDX = MAX( MIN( PORDX,1.D+0 ),1.D-6 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY_BH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY_CO2E( N,PX,PREFX,PORDX,PORTX )
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
!     Compute diffusive and total porosities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, December, 1992.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
      USE GEO_MECH
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PORSTY_CO2E'
      IZN = IZ(N)
      DPX = PX-PREFX
      ISLC50X = ABS(ISLC(50))
!
!---  Geomechanics simulations  --
!
      IF( ISLC(50).NE.0 .AND. MOD(ISLC50X,10).NE.2 .AND. 
     &  MOD(ISLC50X,10).NE.4 ) THEN
!
!---   Davis-Davis porosity versus mean effective stress model  ---
!
       IF( IPROP_GM(1,IZN).EQ.1 ) THEN
!
!---    Diffusive and total mechanical porosity, at iterate level k  ---
!
        PORDMCHKX = (POR(1,IZN)-PROP_GM(5,IZN))*
     &    EXP(-PROP_GM(6,IZN)*SIGV_GM(2,N)) +  PROP_GM(5,IZN)
        PORTMCHKX = (POR(1,IZN)-PROP_GM(5,IZN))*
     &    EXP(-PROP_GM(6,IZN)*SIGV_GM(2,N)) +  PROP_GM(5,IZN)
!
!---    Pressure differential, at iterate level k+1  ---
!
        DPK1X = PX-P_GM(2,N)
!
!---    Diffusive and total flow porosity, at iterate level k+1  ---
!
        DPORX = (POR(1,IZN)-PROP_GM(5,IZN))*PROP_GM(6,IZN)*
     &    PROP_GM(3,IZN)*EXP(-PROP_GM(6,IZN)*SIGV_GM(2,N))
        PORDX = PORDMCHKX + DPORX*DPK1X
        PORTX = PORTMCHKX + DPORX*DPK1X
!
!---   Classical model  ---
!
       ELSE
!
!---    Drained bulk modulus  --
!
        BLKDRNX = PROP_GM(1,IZN)/(3.0D+0*(1.0D+0-2.0D+0*PROP_GM(2,IZN)))
!
!---    1/N  ---
!
        OONMODX = (PROP_GM(3,IZN)-POR(1,IZN))*(1.0D+0-PROP_GM(3,IZN))
     &    /BLKDRNX
!
!---    Volumetric strain differential, at iterate level k  ---
!
        DEVKX = EPSV_GM(2,N)-EPSV_CMP(N)
!
!---    Pressure differential, at iterate level k  ---
!
        DPKX = P_GM(2,N)-PCMP(N)
!
!---    Diffusive and total mechanical porosity, at iterate level k  ---
!
        PORDMCHKX = POR(1,IZN) - PROP_GM(3,IZN)*DEVKX + OONMODX*DPKX
        PORTMCHKX = POR(2,IZN) - PROP_GM(3,IZN)*DEVKX + OONMODX*DPKX
!
!---    Pressure differential, at iterate level k+1  ---
!
        DPK1X = PX-P_GM(2,N)
!
!---    Diffusive and total flow porosity, at iterate level k+1  ---
!
        PORDX = PORDMCHKX + ((PROP_GM(3,IZN)**2)/BLKDRNX+OONMODX)*DPK1X
        PORTX = PORTMCHKX + ((PROP_GM(3,IZN)**2)/BLKDRNX+OONMODX)*DPK1X
       ENDIF
!
!---  No geomechanics  --
!
      ELSE        
!
!---    Pore compressibility w/ fixed bulk volume  ---
!
        IF( ISLC(15).EQ.1 ) THEN
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN)*EXP(DPX*CMP(2,IZN))
            PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN))
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
            PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
            ENDIF
          ENDIF
!
!---    Bulk compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.10 ) THEN
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(3,IZN)/
     &        POR(3,IZN)))
            PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(4,IZN)/
     &        POR(4,IZN)))
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)/
     &        POR(1,IZN)))
            PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &        POR(2,IZN)))
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(1,N)/
     &          POR0(1,N)))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(2,N)/
     &          POR0(2,N)))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)/
     &          POR(1,IZN)))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &          POR(2,IZN)))
            ENDIF
          ENDIF
!
!---    Pore compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.11 ) THEN
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(3,IZN)))
            PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(4,IZN)))
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)))
            PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(1,N)))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(2,N)))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
            ENDIF
          ENDIF
!
!---    Bulk compressibility w/ fixed bulk volume  ---
!
        ELSE
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN) + DPX*CMP(2,IZN)
            PORDFX = POR(4,IZN) + DPX*CMP(2,IZN)
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN) + DPX*CMP(1,IZN)
            PORDMX = POR(2,IZN) + DPX*CMP(1,IZN)
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N) + DPX*CMP(1,IZN)
              PORDX = POR0(2,N) + DPX*CMP(1,IZN)
            ELSE
              PORTX = POR(1,IZN) + DPX*CMP(1,IZN)
              PORDX = POR(2,IZN) + DPX*CMP(1,IZN)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      PORTX = MAX( MIN( PORTX,1.D+0 ),1.D-6 )
      PORDX = MAX( MIN( PORDX,1.D+0 ),1.D-6 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY_CO2E group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY_FRC( NTX,PX,PREFX,PORDX,PORTX )
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
!     Compute diffusive and total porosities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, December, 1992.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GEOM_FRC
      USE FDVP_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PORSTY_FRC'
      IZN = IZ_FRC(NTX)
      DPX = PX-PREFX
!
!---  Pore compressibility w/ fixed bulk volume  ---
!
      IF( ISLC(15).EQ.1 ) THEN
!
!---    Fracture properties (dual porosity model)  ---
!
        IF( ABS(IDP(IZN)).EQ.1 ) THEN
          PORTX = POR(3,IZN)*EXP(DPX*CMP(2,IZN))
          PORDX = POR(4,IZN)*EXP(DPX*CMP(2,IZN))
!
!---    Matrix properties (dual porosity model)  ---
!
        ELSEIF( ABS(IDP(IZN)).EQ.2 ) THEN
          PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
          PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
        ELSE
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0_FRC(1,NTX)*EXP(DPX*CMP(1,IZN))
            PORDX = POR0_FRC(2,NTX)*EXP(DPX*CMP(1,IZN))
          ELSE
            PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
            PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
          ENDIF
        ENDIF
!
!---  Bulk compressibility w/ variable bulk volume  ---
!
      ELSEIF( ISLC(15).EQ.10 ) THEN
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0_FRC(1,NTX)*EXP(DPX*CMP(1,IZN)*
     &        (1.D+0-POR0_FRC(1,NTX)/POR0_FRC(1,NTX)))
            PORDX = POR0_FRC(2,NTX)*EXP(DPX*CMP(1,IZN)*
     &        (1.D+0-POR0_FRC(2,NTX)/POR0_FRC(2,NTX)))
          ELSE
            PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)/
     &        POR(1,IZN)))
            PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &        POR(2,IZN)))
          ENDIF
!
!---  Pore compressibility w/ variable bulk volume  ---
!
      ELSEIF( ISLC(15).EQ.11 ) THEN
!
!---    Reactive transport porosity alteration  ---
!
        IF( ISLC(43).EQ.1 ) THEN
          PORTX = POR0_FRC(1,NTX)*EXP(DPX*CMP(1,IZN)*
     &      (1.D+0-POR0_FRC(1,NTX)))
          PORDX = POR0_FRC(2,NTX)*EXP(DPX*CMP(1,IZN)*
     &      (1.D+0-POR0_FRC(2,NTX)))
        ELSE
          PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)))
          PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
        ENDIF
!
!---  Bulk compressibility w/ fixed bulk volume  ---
!
      ELSE
!
!---    Fracture properties (dual porosity model)  ---
!
        IF( ABS(IDP(IZN)).EQ.1 ) THEN
          PORTX = POR(3,IZN) + DPX*CMP(2,IZN)
          PORDX = POR(4,IZN) + DPX*CMP(2,IZN)
!
!---    Matrix properties (dual porosity model)  ---
!
        ELSEIF( ABS(IDP(IZN)).EQ.2 ) THEN
          PORTX = POR(1,IZN) + DPX*CMP(1,IZN)
          PORDX = POR(2,IZN) + DPX*CMP(1,IZN)
        ELSE
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0_FRC(1,NTX) + DPX*CMP(1,IZN)
            PORDX = POR0_FRC(2,NTX) + DPX*CMP(1,IZN)
          ELSE
            PORTX = POR(1,IZN) + DPX*CMP(1,IZN)
            PORDX = POR(2,IZN) + DPX*CMP(1,IZN)
          ENDIF
        ENDIF
      ENDIF
      PORTX = MAX( MIN( PORTX,1.D+0 ),1.D-6 )
      PORDX = MAX( MIN( PORDX,1.D+0 ),1.D-6 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY_GT( N,PX,PREFX,PORDX,PORTX )
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
!     Compute diffusive and total porosities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, December, 1992.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
      USE GEO_MECH
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PORSTY_GT'
      IZN = IZ(N)
      DPX = PX-PREFX
      ISLC50X = ABS(ISLC(50))
!
!---  Geomechanics simulations  --
!
      IF( ISLC(50).NE.0 .AND. MOD(ISLC50X,10).NE.2 .AND. 
     &  MOD(ISLC50X,10).NE.4 ) THEN
!
!---    Drained bulk modulus  --
!
        BLKDRNX = PROP_GM(1,IZN)/(3.0D+0*(1.0D+0-2.0D+0*PROP_GM(2,IZN)))
!
!---    1/N  ---
!
        OONMODX = (PROP_GM(3,IZN)-POR(1,IZN))*(1.0D+0-PROP_GM(3,IZN))
     &    /BLKDRNX
!
!---    Volumetric strain differential, at iterate level k  ---
!
        DEVKX = EPSV_GM(2,N)-EPSV_CMP(N)
!
!---    Pressure differential, at iterate level k  ---
!
        DPKX = P_GM(2,N)-PCMP(N)
!
!---    Diffusive and total mechanical porosity, at iterate level k  ---
!
        PORDMCHKX = POR(1,IZN) - PROP_GM(3,IZN)*DEVKX + OONMODX*DPKX
        PORTMCHKX = POR(2,IZN) - PROP_GM(3,IZN)*DEVKX + OONMODX*DPKX
!
!---    Pressure differential, at iterate level k+1  ---
!
        DPK1X = PX-P_GM(2,N)
!
!---    Diffusive and total flow porosity, at iterate level k+1  ---
!
        PORDX = PORDMCHKX + ((PROP_GM(3,IZN)**2)/BLKDRNX+OONMODX)*DPK1X
        PORTX = PORTMCHKX + ((PROP_GM(3,IZN)**2)/BLKDRNX+OONMODX)*DPK1X
!
!---  No geomechanics  --
!
      ELSE        
!
!---    Pore compressibility w/ fixed bulk volume  ---
!
        IF( ISLC(15).EQ.1 ) THEN
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN)*EXP(DPX*CMP(2,IZN))
            PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN))
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
            PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
            ENDIF
          ENDIF
!
!---    Bulk compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.10 ) THEN
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(3,IZN)/
     &        POR(3,IZN)))
            PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(4,IZN)/
     &        POR(4,IZN)))
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)/
     &        POR(1,IZN)))
            PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &        POR(2,IZN)))
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(1,N)/
     &          POR0(1,N)))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(2,N)/
     &          POR0(2,N)))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)/
     &          POR(1,IZN)))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &          POR(2,IZN)))
            ENDIF
          ENDIF
!
!---    Pore compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.11 ) THEN
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(3,IZN)))
            PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(4,IZN)))
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)))
            PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(1,N)))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(2,N)))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
            ENDIF
          ENDIF
!
!---    Bulk compressibility w/ fixed bulk volume  ---
!
        ELSE
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN) + DPX*CMP(2,IZN)
            PORDFX = POR(4,IZN) + DPX*CMP(2,IZN)
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN) + DPX*CMP(1,IZN)
            PORDMX = POR(2,IZN) + DPX*CMP(1,IZN)
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N) + DPX*CMP(1,IZN)
              PORDX = POR0(2,N) + DPX*CMP(1,IZN)
            ELSE
              PORTX = POR(1,IZN) + DPX*CMP(1,IZN)
              PORDX = POR(2,IZN) + DPX*CMP(1,IZN)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      PORTX = MAX( MIN( PORTX,1.D+0 ),1.D-6 )
      PORDX = MAX( MIN( PORDX,1.D+0 ),1.D-6 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY_HYDT_KE( N,PX,PREFX,PORDX,PORTX,SHX )
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
!     Compute diffusive and total porosities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, 21 June 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
      USE GEO_MECH
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 PROP_GMX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PORSTY_HYDT_KE'
      IZN = IZ(N)
      DPX = PX-PREFX
      ISLC50X = ABS(ISLC(50))
!
!---  Geomechanics simulations  --
!
      IF( ISLC(50).NE.0 .AND. MOD(ISLC50X,10).NE.2 .AND. 
     &  MOD(ISLC50X,10).NE.4 ) THEN
!
!---    Hydrate C-Factor composite model  ---
!
        IF( IHCM_GM(IZN).EQ.1 ) THEN
          PROP_GMX(1) = PROP_GM(1,IZN) + 
     &      PROP_GM(5,IZN)*SHX*PROP_GM(6,IZN)
        ELSE
          PROP_GMX(1) = PROP_GM(1,IZN)
        ENDIF
        PROP_GMX(2) = PROP_GM(2,IZN)
!
!---    Drained bulk modulus  --
!
        BLKDRNX = PROP_GMX(1)/(3.0D+0*(1.0D+0-2.0D+0*PROP_GMX(2)))
!
!---    1/N  ---
!
        OONMODX = (PROP_GM(3,IZN)-POR(1,IZN))*(1.0D+0-PROP_GM(3,IZN))
     &    /BLKDRNX
!
!---    Volumetric strain differential, at iterate level k  ---
!
        DEVKX = EPSV_GM(2,N)-EPSV_CMP(N)
!
!---    Pressure differential, at iterate level k  ---
!
        DPKX = P_GM(2,N)-PCMP(N)
!
!---    Diffusive and total mechanical porosity, at iterate level k  ---
!
        PORDMCHKX = POR(1,IZN) - PROP_GM(3,IZN)*DEVKX + OONMODX*DPKX
        PORTMCHKX = POR(2,IZN) - PROP_GM(3,IZN)*DEVKX + OONMODX*DPKX
!
!---    Pressure differential, at iterate level k+1  ---
!
        DPK1X = PX-P_GM(2,N)
!
!---    Diffusive and total flow porosity, at iterate level k+1  ---
!
        PORDX = PORDMCHKX + ((PROP_GM(3,IZN)**2)/BLKDRNX+OONMODX)*DPK1X
        PORTX = PORTMCHKX + ((PROP_GM(3,IZN)**2)/BLKDRNX+OONMODX)*DPK1X
!
!---  No geomechanics  --
!
      ELSE        
!
!---    Pore compressibility w/ fixed bulk volume  ---
!
        IF( ISLC(15).EQ.1 ) THEN
!
!---      Fracture properties (dual porosity model)  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
            PORTX = POR(3,IZN)*EXP(DPX*CMP(2,IZN))
            PORDX = POR(4,IZN)*EXP(DPX*CMP(2,IZN))
!
!---      Matrix properties (dual porosity model)  ---
!
          ELSEIF( ABS(IDP(IZN)).EQ.2 ) THEN
            PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
            PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
            ENDIF
          ENDIF
!
!---    Bulk compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.10 ) THEN
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(1,N)/
     &          POR0(1,N)))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(2,N)/
     &          POR0(2,N)))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)/
     &          POR(1,IZN)))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &          POR(2,IZN)))
            ENDIF
!
!---    Pore compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.11 ) THEN
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(1,N)))
            PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(2,N)))
          ELSE
            PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)))
            PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
          ENDIF
!
!---    Bulk compressibility w/ fixed bulk volume  ---
!
        ELSE
!
!---      Fracture properties (dual porosity model)  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
            PORTX = POR(3,IZN) + DPX*CMP(2,IZN)
            PORDX = POR(4,IZN) + DPX*CMP(2,IZN)
!
!---      Matrix properties (dual porosity model)  ---
!
          ELSEIF( ABS(IDP(IZN)).EQ.2 ) THEN
            PORTX = POR(1,IZN) + DPX*CMP(1,IZN)
            PORDX = POR(2,IZN) + DPX*CMP(1,IZN)
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N) + DPX*CMP(1,IZN)
              PORDX = POR0(2,N) + DPX*CMP(1,IZN)
            ELSE
              PORTX = POR(1,IZN) + DPX*CMP(1,IZN)
              PORDX = POR(2,IZN) + DPX*CMP(1,IZN)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      PORTX = MAX( MIN( PORTX,1.D+0 ),1.D-6 )
      PORDX = MAX( MIN( PORDX,1.D+0 ),1.D-6 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY_M( N,PX,PREFX,PORDX,PORTX )
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
!     Compute matrix diffusive and total porosities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 December 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PORSTY_M'
      DPX = PX-PREFX
      IZN = IZ(N)
!
!---  Pore compressibility w/ fixed bulk volume  ---
!
      IF( ISLC(15).EQ.1 ) THEN
        PORTX = POR(3,IZN)*EXP(DPX*CMP(1,IZN))
        PORDX = POR(4,IZN)*EXP(DPX*CMP(1,IZN))
!
!---  Bulk compressibility w/ variable bulk volume  ---
!
      ELSEIF( ISLC(15).EQ.10 ) THEN
        PORTX = POR(3,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(3,IZN)/
     &    POR(3,IZN)))
        PORDX = POR(4,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(4,IZN)/
     &    POR(4,IZN)))
!
!---  Pore compressibility w/ variable bulk volume  ---
!
      ELSEIF( ISLC(15).EQ.11 ) THEN
        PORTX = POR(3,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(3,IZN)))
        PORDX = POR(4,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(4,IZN)))
!
!---  Bulk compressibility w/ fixed bulk volume  ---
!
      ELSE
        PORTX = POR(3,IZN) + DPX*CMP(1,IZN)
        PORDX = POR(4,IZN) + DPX*CMP(1,IZN)
      ENDIF
      PORTX = MAX( MIN( PORTX,1.D+0 ),1.D-6 )
      PORDX = MAX( MIN( PORDX,1.D+0 ),1.D-6 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY_M group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY_SEQ( N,PX,PREFX,PORDX,PORTX )
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
!     Compute diffusive and total porosities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, December, 1992.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
      USE GEO_MECH
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PORSTY_SEQ'
      IZN = IZ(N)
      DPX = PX-PREFX
      ISLC50X = ABS(ISLC(50))
!
!---  Geomechanics simulations  --
!
      IF( ISLC(50).NE.0 .AND. MOD(ISLC50X,10).NE.2 .AND. 
     &  MOD(ISLC50X,10).NE.4 ) THEN
!
!---   Davis-Davis porosity versus mean effective stress model  ---
!
       IF( IPROP_GM(1,IZN).EQ.1 ) THEN
!
!---    Diffusive and total mechanical porosity, at iterate level k  ---
!
        PORDMCHKX = (POR(1,IZN)-PROP_GM(5,IZN))*
     &    EXP(-PROP_GM(6,IZN)*SIGV_GM(2,N)) +  PROP_GM(5,IZN)
        PORTMCHKX = (POR(1,IZN)-PROP_GM(5,IZN))*
     &    EXP(-PROP_GM(6,IZN)*SIGV_GM(2,N)) +  PROP_GM(5,IZN)
!
!---    Pressure differential, at iterate level k+1  ---
!
        DPK1X = PX-P_GM(2,N)
!
!---    Diffusive and total flow porosity, at iterate level k+1  ---
!
        DPORX = (POR(1,IZN)-PROP_GM(5,IZN))*PROP_GM(6,IZN)*
     &    PROP_GM(3,IZN)*EXP(-PROP_GM(6,IZN)*SIGV_GM(2,N))
        PORDX = PORDMCHKX + DPORX*DPK1X
        PORTX = PORTMCHKX + DPORX*DPK1X
!
!---   Classical model  ---
!
       ELSE
!
!---    Drained bulk modulus  --
!
        BLKDRNX = PROP_GM(1,IZN)/(3.0D+0*(1.0D+0-2.0D+0*PROP_GM(2,IZN)))
!
!---    1/N  ---
!
        OONMODX = (PROP_GM(3,IZN)-POR(1,IZN))*(1.0D+0-PROP_GM(3,IZN))
     &    /BLKDRNX
!
!---    Volumetric strain differential, at iterate level k  ---
!
        DEVKX = EPSV_GM(2,N)-EPSV_CMP(N)
!
!---    Pressure differential, at iterate level k  ---
!
        DPKX = P_GM(2,N)-PCMP(N)
!
!---    Diffusive and total mechanical porosity, at iterate level k  ---
!
        PORDMCHKX = POR(1,IZN) - PROP_GM(3,IZN)*DEVKX + OONMODX*DPKX
        PORTMCHKX = POR(2,IZN) - PROP_GM(3,IZN)*DEVKX + OONMODX*DPKX
!
!---    Pressure differential, at iterate level k+1  ---
!
        DPK1X = PX-P_GM(2,N)
!
!---    Diffusive and total flow porosity, at iterate level k+1  ---
!
        PORDX = PORDMCHKX + ((PROP_GM(3,IZN)**2)/BLKDRNX+OONMODX)*DPK1X
        PORTX = PORTMCHKX + ((PROP_GM(3,IZN)**2)/BLKDRNX+OONMODX)*DPK1X
       ENDIF
!
!---  No geomechanics  --
!
      ELSE        
!
!---    Pore compressibility w/ fixed bulk volume  ---
!
        IF( ISLC(15).EQ.1 ) THEN
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN)*EXP(DPX*CMP(2,IZN))
            PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN))
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
            PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN))
            ENDIF
          ENDIF
!
!---    Bulk compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.10 ) THEN
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(3,IZN)/
     &        POR(3,IZN)))
            PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(4,IZN)/
     &        POR(4,IZN)))
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)/
     &        POR(1,IZN)))
            PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &        POR(2,IZN)))
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(1,N)/
     &          POR0(1,N)))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(2,N)/
     &          POR0(2,N)))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)/
     &          POR(1,IZN)))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)/
     &          POR(2,IZN)))
            ENDIF
          ENDIF
!
!---    Pore compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.11 ) THEN
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(3,IZN)))
            PORDFX = POR(4,IZN)*EXP(DPX*CMP(2,IZN)*(1.D+0-POR(4,IZN)))
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)))
            PORDMX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(1,N)))
              PORDX = POR0(2,N)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR0(2,N)))
            ELSE
              PORTX = POR(1,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(1,IZN)))
              PORDX = POR(2,IZN)*EXP(DPX*CMP(1,IZN)*(1.D+0-POR(2,IZN)))
            ENDIF
          ENDIF
!
!---    Bulk compressibility w/ fixed bulk volume  ---
!
        ELSE
!
!---      Dual porosity model  ---
!
          IF( ABS(IDP(IZN)).EQ.1 ) THEN
!
!---        Fracture properties (dual porosity model)  ---
!
            PORTFX = POR(3,IZN) + DPX*CMP(2,IZN)
            PORDFX = POR(4,IZN) + DPX*CMP(2,IZN)
!
!---        Matrix properties (dual porosity model)  ---
!
            PORTMX = POR(1,IZN) + DPX*CMP(1,IZN)
            PORDMX = POR(2,IZN) + DPX*CMP(1,IZN)
!
!---        Effective properties (dual porosity model)  ---
!
            PORTX = PORTFX + (1.D+0-PORTFX)*PORTMX
            PORDX = PORDFX + (1.D+0-PORDFX)*PORDMX
!
!---      Single porosity model  ---
!
          ELSE
!
!---        Reactive transport porosity alteration  ---
!
            IF( ISLC(43).EQ.1 ) THEN
              PORTX = POR0(1,N) + DPX*CMP(1,IZN)
              PORDX = POR0(2,N) + DPX*CMP(1,IZN)
            ELSE
              PORTX = POR(1,IZN) + DPX*CMP(1,IZN)
              PORDX = POR(2,IZN) + DPX*CMP(1,IZN)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      PORTX = MAX( MIN( PORTX,1.D+0 ),1.D-6 )
      PORDX = MAX( MIN( PORDX,1.D+0 ),1.D-6 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY_SEQ group  ---
!
      RETURN
      END




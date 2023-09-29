!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC_FRC_EOR
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Compute initial solute fault concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 30 November 2020
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE GEOM_FRC
      USE FDVP_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( IEQC.EQ.0 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CISC_FRC_EOR'
!
!---  Loop over solutes  ---
!
      DO NSL = 1,NSOLU
!
!---    Loop over faults  ---
!
        DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
          N_DB = NTX
          IZN = IZ_FRC(NTX)
          IF( IPCL(NSL).EQ.2 ) THEN
            XVS = SL_FRC(2,NTX)*RHOS_FRC(NTX)*PCSL(1,IZN,NSL)*
     &        (1.D+0-PORT_FRC(2,NTX))
          ELSE
            XVS = RHOS_FRC(NTX)*PCSL(1,IZN,NSL)*
     &        (1.D+0-PORT_FRC(2,NTX))
          ENDIF
          XVG = SG_FRC(2,NTX)*PORD_FRC(2,NTX)
          XVL = SL_FRC(2,NTX)*PORD_FRC(2,NTX)
          XVN = SN_FRC(2,NTX)*PORD_FRC(2,NTX)
!
!---      Constant gas-aqueous partition coefficient  ---
!
          IF( IPCGL(NSL).EQ.0 ) THEN
            PCGLX = PCGL(1,NSL)
!
!---      Temperature dependent gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.1 ) THEN
            TK = T_FRC(2,NTX)+TABS
            PCGLX = EXP( PCGL(1,NSL) 
     &        + PCGL(2,NSL)/TK + PCGL(3,NSL)*LOG(TK) 
     &        + PCGL(4,NSL)*TK + PCGL(5,NSL)*TK**2 )
!
!---      Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.2 ) THEN
            PCGLX = RHOG_FRC(2,NTX)*XGW_FRC(2,NTX)/
     &        (RHOL_FRC(2,NTX)*XLW_FRC(2,NTX))
          ENDIF
          PCGLX = MAX( PCGLX,1.D-20 )
          PCGLX = MIN( PCGLX,1.D+20 )
!
!---      Phase-volumetric concentration ratios  ---
!
          IF( ICT_FRC(NFX,NSL).EQ.2 ) THEN
            C_FRC(NTX,NSL) = 
     &        C_FRC(NTX,NSL)*(XVS + XVL + XVN/PCLN(1,NSL) + XVG*PCGLX)
          ELSEIF( ICT_FRC(NFX,NSL).EQ.3 ) THEN
            C_FRC(NTX,NSL) = 
     &        C_FRC(NTX,NSL)*((XVS + XVL + XVN)/PCGLX + XVG)
          ELSEIF( ICT_FRC(NFX,NSL).EQ.4 ) THEN
            C_FRC(NTX,NSL) = 
     &        C_FRC(NTX,NSL)*((XVS + XVL + XVG*PCGLX)*PCLN(1,NSL) + XVN)
          ENDIF
!
!---    End loop over fault triangles  ---
!
        ENDDO
!
!---    End loop over faults  ---
!
        ENDDO
!
!---  End loop over solutes  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CISC_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFX_FRC_EOR( NSL )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Compute solute fault transport flux for gas, aqueous,
!     and nonaqueous liquid phases.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 30 November 2020
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE NAPL
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
      USE FDVN_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFX_FRC_EOR'
!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
        DO NT1X = IP_FRC(1,NFX),IP_FRC(2,NFX)
          IF( IXP_FRC(NT1X).EQ.0 ) CYCLE
!
!---      Gas contributions, 
!         molecular diffusion coefficients at triangle  ---
!
          TCOR = (T_FRC(2,NT1X)+TABS)/TSPRF
          PCOR = (PG_FRC(2,NT1X)+PATM)/PATM
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SG_FRC(2,NT1X)*PORD_FRC(2,NT1X).LT.EPSL ) THEN
            TORGX = 0.D+0
          ELSE
            TORGX = (PORD_FRC(2,NT1X)*(SG_FRC(2,NT1X)**7))**(THIRD)
          ENDIF
          SDFG1 = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMC1 = SG_FRC(2,NT1X)*PORD_FRC(2,NT1X)
          FCG1 = YG_FRC(NT1X,NSL)/(VMC1+SMALL)
          DG1 = TORGX*SG_FRC(2,NT1X)*PORD_FRC(2,NT1X)*SDFG1
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            MT2X = NFLD - NXP + IXP_FRC(NT2X)
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Neglect hydrodynamic dispersion  ---
!
            DPGX = 0.D+0
!
!---        Molecular diffusion coefficients at adjacent triangle  ---
!
            TCOR = (T_FRC(2,NT2X)+TABS)/TSPRF
            PCOR = (PG_FRC(2,NT2X)+PATM)/PATM
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( SG_FRC(2,NT2X)*PORD_FRC(2,NT2X).LT.EPSL ) THEN
              TORGX = 0.D+0
            ELSE
              TORGX = (PORD_FRC(2,NT2X)*(SG_FRC(2,NT2X)**7))**(THIRD)
            ENDIF
            SDFG2 = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            VMC2 = SG_FRC(2,NT2X)*PORD_FRC(2,NT2X)
            FCG2 = YG_FRC(NT2X,NSL)/(VMC2+SMALL)
            DG2 = TORGX*SG_FRC(2,NT2X)*SDFG2*PORD_FRC(2,NT2X)
!
!---        Interfacial fault aperture  ---
!
            INDX = -1
            APX = DIFMN( APM_FRC(2,NT1X),APM_FRC(2,NT2X),
     &        DFF1X,DFF1X,ZERO,INDX )
!
!---        Interfacial diffusion/dispersion coefficient  ---
!
            INDX = 16
            DGX = DIFMN(DG1,DG2,DFF1X,DFF2X,UFFG(1,NCX),INDX)
            DGX = AFF_FRC(NCX)*APX*(DGX+DPGX)/DFF_FRC(NCX)
            FGX = AFF_FRC(NCX)*APX*UFFG(1,NCX)
            IF( MOD(ISLC(23),100)/10.EQ.1 ) FGX = 0.D+0
!
!---        Patankar solute transport  ---
!
            AGX = MAX(FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
            AGPX = MAX(-FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
            UFFC(NCX,NSL) = (C_FRC(NT1X,NSL)*AGX*FCG1 - 
     &        C_FRC(NT2X,NSL)*AGPX*FCG2)
          ENDDO
!
!---      Aqueous contributions,
!         molecular diffusion coefficients at triangle  ---
!
          TCOR = (T_FRC(2,NT1X)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SL_FRC(2,NT1X)*PORD_FRC(2,NT1X).LT.EPSL ) THEN
            TORLX = 0.D+0
          ELSE
            TORLX = (PORD_FRC(2,NT1X)*(SL_FRC(2,NT1X)**7))**(THIRD)
          ENDIF
          SDFL1 = SMDL(NSL)*TCOR*(VISRL/VISL_FRC(2,NT1X))
          VMC1 = SL_FRC(2,NT1X)*PORD_FRC(2,NT1X)
          FCL1 = YL_FRC(NT1X,NSL)/(VMC1+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DL1 = SDCL_FRC(1,NT1X,NSL)*SDCL_FRC(2,NT1X,NSL)*
     &        EXP(VMC1*SDCL_FRC(3,NT1X,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DL1 = TORLX*VMC1*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DL1 = SDCL_FRC(1,NT1X,NSL)*SDCL_FRC(2,NT1X,NSL)*
     &        VMC1**SDCL_FRC(3,NT1X,NSL)
          ELSE
            DL1 = TORLX*VMC1*SDFL1
          ENDIF
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            MT2X = NFLD - NXP + IXP_FRC(NT2X)
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Neglect hydrodynamic dispersion  ---
!
            DPLX = 0.D+0
!
!---        Molecular diffusion coefficients at adjacent triangle  ---
!
            TCOR = (T_FRC(2,NT2X)+TABS)/TSPRF
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( SL_FRC(2,NT2X)*PORD_FRC(2,NT2X).LT.EPSL ) THEN
              TORLX = 0.D+0
            ELSE
              TORLX = (PORD_FRC(2,NT2X)*(SL_FRC(2,NT2X)**7))**(THIRD)
            ENDIF
            SDFL2 = SMDL(NSL)*TCOR*(VISRL/VISL_FRC(2,NT2X))
            VMC2 = SL_FRC(2,NT2X)*PORD_FRC(2,NT2X)
            FCL2 = YL_FRC(NT2X,NSL)/(VMC2+SMALL)
            IF( IEDL(NSL).EQ.2 ) THEN
              DL2 = SDCL_FRC(1,NT2X,NSL)*SDCL_FRC(2,NT2X,NSL)*
     &          EXP(VMC2*SDCL_FRC(3,NT2X,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DL2 = TORLX*VMC2*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DL2 = SDCL_FRC(1,NT2X,NSL)*SDCL_FRC(2,NT2X,NSL)*
     &          VMC2**SDCL_FRC(3,NT2X,NSL)
            ELSE
              DL2 = TORLX*VMC2*SDFL2
            ENDIF
!
!---        Interfacial fault aperture  ---
!
            INDX = -1
            APX = DIFMN( APM_FRC(2,NT1X),APM_FRC(2,NT2X),
     &        DFF1X,DFF1X,ZERO,INDX )
!
!---        Interfacial diffusion/dispersion coefficient  ---
!
            INDX = 16
            DLX = DIFMN(DL1,DL2,DFF1X,DFF2X,UFFL(1,NCX),INDX)
            DLX = AFF_FRC(NCX)*APX*(DLX+DPLX)/DFF_FRC(NCX)
            FLX = AFF_FRC(NCX)*APX*UFFL(1,NCX)
            IF( MOD(ISLC(23),10).EQ.1 ) FLX = 0.D+0
!
!---        Patankar solute transport  ---
!
            ALX = MAX(FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
            ALPX = MAX(-FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
            UFFC(NCX,NSL) = UFFC(NCX,NSL)+ (C_FRC(NT1X,NSL)*ALX*FCL1 - 
     &        C_FRC(NT2X,NSL)*ALPX*FCL2)
          ENDDO
!
!---      Nonaqueous liquid contributions,
!         molecular diffusion coefficients at triangle  ---
!
          TCOR = (T_FRC(2,NT1X)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SN_FRC(2,NT1X)*PORD_FRC(2,NT1X).LT.EPSL ) THEN
            TORNX = 0.D+0
          ELSE
            TORNX = (PORD_FRC(2,NT1X)*(SN_FRC(2,NT1X)**7))**(THIRD)
          ENDIF
          SDFN1 = SMDN(NSL)*TCOR*(VISRN/VISN_FRC(2,NT1X))
          VMC1 = SN_FRC(2,NT1X)*PORD_FRC(2,NT1X)
          FCN1 = YN_FRC(NT1X,NSL)/(VMC1+SMALL)
          DN1 = TORNX*VMC1*SDFN1
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            MT2X = NFLD - NXP + IXP_FRC(NT2X)
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Neglect hydrodynamic dispersion  ---
!
            DPLX = 0.D+0
!
!---        Molecular diffusion coefficients at adjacent triangle  ---
!
            TCOR = (T_FRC(2,NT2X)+TABS)/TSPRF
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( SN_FRC(2,NT2X)*PORD_FRC(2,NT2X).LT.EPSL ) THEN
              TORNX = 0.D+0
            ELSE
              TORNX = (PORD_FRC(2,NT2X)*(SN_FRC(2,NT2X)**7))**(THIRD)
            ENDIF
            SDFN2 = SMDN(NSL)*TCOR*(VISRN/VISN_FRC(2,NT2X))
            VMC2 = SN_FRC(2,NT2X)*PORD_FRC(2,NT2X)
            FCN2 = YN_FRC(NT2X,NSL)/(VMC2+SMALL)
            DN2 = TORNX*VMC2*SDFN2
!
!---        Interfacial fault aperture  ---
!
            INDX = -1
            APX = DIFMN( APM_FRC(2,NT1X),APM_FRC(2,NT2X),
     &        DFF1X,DFF1X,ZERO,INDX )
!
!---        Interfacial diffusion/dispersion coefficient  ---
!
            INDX = 16
            DNX = DIFMN(DN1,DN2,DFF1X,DFF2X,UFFN(1,NCX),INDX)
            DNX = AFF_FRC(NCX)*APX*(DNX+DPNX)/DFF_FRC(NCX)
            FNX = AFF_FRC(NCX)*APX*UFFN(1,NCX)
            IF( MOD(ISLC(23),10).EQ.1 ) FNX = 0.D+0
!
!---        Patankar solute transport  ---
!
            ANX = MAX(FNX,ZERO)
     &        + DNX*MAX((ONE-(TENTH*ABS(FNX)/(DNX+SMALL)))**5,ZERO)
            ANPX = MAX(-FNX,ZERO)
     &        + DNX*MAX((ONE-(TENTH*ABS(FNX)/(DNX+SMALL)))**5,ZERO)
            UFFC(NCX,NSL) = UFFC(NCX,NSL)+ (C_FRC(NT1X,NSL)*ANX*FCN1 - 
     &        C_FRC(NT2X,NSL)*ANPX*FCN2)
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFX_FRC_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFX_FM_EOR( NSL )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Compute solute fault to matrix transport flux for gas, 
!     aqueous, and nonaqueous liquid phases.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 30 November 2020
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE NAPL
      USE HYST
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
      USE FDVP
      USE FDVN_FRC
      USE FDVN
      USE FDVG
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFX_FM_EOR'
!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
          MTX = NFLD - NXP + NTX
!
!---      Loop over fault triangle to grid cell connections  ---
!
          DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
            N = INCM_FRC(NCX)
            IZN = IZ(N)
!
!---        Gas contributions, 
!           fault side gas diffusion  ---
!
            TCOR = (T_FRC(2,NTX)+TABS)/TSPRF
            PCOR = (PG_FRC(2,NTX)+PATM)/PATM
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( PORD_FRC(2,NTX)*SG_FRC(2,NTX).LT.EPSL ) THEN
              TORGX = 0.D+0
            ELSE
              TORGX = (PORD_FRC(2,NTX)*(SG_FRC(2,NTX)**7))**(THIRD)
            ENDIF
            SDFGF = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            VMCF = PORD_FRC(2,NTX)*SG_FRC(2,NTX)
            FCGF = YG_FRC(NTX,NSL)/(VMCF+SMALL)
            DGF = TORGX*SG_FRC(2,NTX)*SDFGF
!
!---        Grid-cell side gas diffusion  ---
!
            TCOR = (T(2,N)+TABS)/TSPRF
            PCOR = (PG(2,N)+PATM)/PATM
            SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            VMCP = SG(2,N)*PORD(2,N)
            FCGP = YG(N,NSL)/(VMCP+SMALL)
            DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
            INDX = 16
            DGX = DIFMN(DGF,DGP,DFN_FRC(NCX),DFN_FRC(NCX),ZERO,INDX)
            DGX = AFN_FRC(NCX)*DGX/DFN_FRC(NCX)
            FGX = UFMG(NCX)
!
!---        Patankar solute transport between fault and grid cell  ---
!
            AGX = MAX(FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
            AGPX = MAX(-FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
            UFMC(NCX,NSL) = (C_FRC(NTX,NSL)*AGX*FCGF - 
     &        C(N,NSL)*AGPX*FCGP)
!
!---        Aqueous contributions,
!           fault side aqueous diffusion  ---
!
            TCOR = (T_FRC(2,NTX)+TABS)/TSPRF
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( PORD_FRC(2,NTX)*SL_FRC(2,NTX).LT.EPSL ) THEN
              TORLX = 0.D+0
            ELSE
              TORLX = (PORD_FRC(2,NTX)*(SL_FRC(2,NTX)**7))**(THIRD)
            ENDIF
            SDFLF = SMDL(NSL)*TCOR*(VISRL/VISL_FRC(2,NTX))
            VMCF = PORD_FRC(2,NTX)*SL_FRC(2,NTX)
            FCLF = YL_FRC(NTX,NSL)/(VMCF+SMALL)
            IF( IEDL(NSL).EQ.2 ) THEN
              DLF = SDCL_FRC(1,NTX,NSL)*SDCL_FRC(2,NTX,NSL)*
     &          EXP(VMCF*SDCL_FRC(3,NTX,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLF = TORLX*VMCF*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DLF = SDCL_FRC(1,NTX,NSL)*SDCL_FRC(2,NTX,NSL)*
     &          VMCF**SDCL_FRC(3,NTX,NSL)
            ELSE
              DLF = TORLX*VMCF*SDFLF
            ENDIF
!
!---        Grid-cell side aqueous diffusion  ---
!
            TCOR = (T(2,N)+TABS)/TSPRF
            SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
            VMCP = SL(2,N)*PORD(2,N)
            FCLP = YL(N,NSL)/(VMCP+SMALL)
            IF( IEDL(NSL).EQ.2 ) THEN
              DLP = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(VMCP*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLP = TORL(2,N)*VMCP*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DLP = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          VMCP**SDCL(3,IZN,NSL)
            ELSE
              DLP = TORL(2,N)*VMCP*SDFLP
            ENDIF
            INDX = 16
            DLX = DIFMN(DLF,DLP,DFN_FRC(NCX),DFN_FRC(NCX),ZERO,INDX)
            DLX = AFN_FRC(NCX)*DLX/DFN_FRC(NCX)
            FLX = UFML(NCX)
!
!---        Patankar solute transport between fault and grid cell  ---
!
            ALX = MAX(FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
            ALPX = MAX(-FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
            UFMC(NCX,NSL) = UFMC(NCX,NSL)+ (C_FRC(NTX,NSL)*ALX*FCLF - 
     &          C(N,NSL)*ALPX*FCLP)
!
!---        Nonaqueous liquid contributions,
!           fault side nonaqueous liquid diffusion  ---
!
            TCOR = (T_FRC(2,NTX)+TABS)/TSPRF
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( PORD_FRC(2,NTX)*SN_FRC(2,NTX).LT.EPSL ) THEN
              TORNX = 0.D+0
            ELSE
              TORNX = (PORD_FRC(2,NTX)*(SN_FRC(2,NTX)**7))**(THIRD)
            ENDIF
            SDFNF = SMDN(NSL)*TCOR*(VISRN/VISN_FRC(2,NTX))
            VMCF = PORD_FRC(2,NTX)*SN_FRC(2,NTX)
            FCNF = YN_FRC(NTX,NSL)/(VMCF+SMALL)
            DNF = TORNX*VMCF*SDFNF
!
!---        Grid-cell side aqueous diffusion  ---
!
            TCOR = (T(2,N)+TABS)/TSPRF
            SDFNP = SMDN(NSL)*TCOR*(VISRN/VISN(2,N))
            VMCP = SN(2,N)*PORD(2,N)
            FCNP = YN(N,NSL)/(VMCP+SMALL)
            DNP = TORN(2,N)*VMCP*SDFNP
            INDX = 16
            DNX = DIFMN(DNF,DNP,DFN_FRC(NCX),DFN_FRC(NCX),ZERO,INDX)
            DNX = AFN_FRC(NCX)*DNX/DFN_FRC(NCX)
            FNX = UFMN(NCX)
!
!---        Patankar solute transport between fault and grid cell  ---
!
            ANX = MAX(FNX,ZERO)
     &        + DNX*MAX((ONE-(TENTH*ABS(FNX)/(DNX+SMALL)))**5,ZERO)
            ANPX = MAX(-FNX,ZERO)
     &        + DNX*MAX((ONE-(TENTH*ABS(FNX)/(DNX+SMALL)))**5,ZERO)
            UFMC(NCX,NSL) = UFMC(NCX,NSL)+ (C_FRC(NTX,NSL)*ANX*FCNF - 
     &          C(N,NSL)*ANPX*FCNP)
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFX_FM_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXZ_FRC_EOR( NSL )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Zero solute transport fluxes for fault flow.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 30 November 2020
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GEOM_FRC
      USE FLUX_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFXZ_FRC_EOR'
!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            UFFC(NCX,NSL) = 0.D+0
          ENDDO
!
!---      Loop over fault triangle to grid cell connections  ---
!
          DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
            UFMC(NCX,NSL) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFXZ_FRC_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBG_FRC_EOR( NSL )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Loads the matrix elements and solution vector for the
!     gas convective-dispersive mass transport equation for faults.
!
!     The Jacobian matrix is initially configured assuming zero-flux
!     boundary conditions.  The matrix is then updated for other
!     user-specified boundary conditions.
!
!     Matrix elements are stored in the array ALU.
!     Elements for the right-hand-side are stored in the array BLU.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 30 November 2020
!
#ifdef petsc
!
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: IROW,ICOL
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SJCBG_FRC_EOR'
!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
        DO NT1X = IP_FRC(1,NFX),IP_FRC(2,NFX)
          IF( IXP_FRC(NT1X).EQ.0 ) CYCLE
          MT1X = NFLD - NXP + IXP_FRC(NT1X)
          MA = 1
!
!---      SPLib or Lis solvers  ---
!
          IF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MCD = KLUC_FRC(IXP_FRC(NT1X),MA)
            MA = MA + 1
          ENDIF
!
!---      Molecular diffusion coefficients at triangle  ---
!
          TCOR = (T_FRC(2,NT1X)+TABS)/TSPRF
          PCOR = (PG_FRC(2,NT1X)+PATM)/PATM
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SG_FRC(2,NT1X).LT.EPSL ) THEN
            TORGX = 0.D+0
          ELSE
            TORGX = (SG_FRC(2,NT1X)**7)**(THIRD)
          ENDIF
          SDFG1 = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMC1 = SG_FRC(2,NT1X)
          FCG1 = YG_FRC(NT1X,NSL)/(VMC1+SMALL)
          DG1 = TORGX*SG_FRC(2,NT1X)*SDFG1
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            MT2X = NFLD - NXP + IXP_FRC(NT2X)
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Neglect hydrodynamic dispersion  ---
!
            DPGX = 0.D+0
!
!---        Molecular diffusion coefficients at adjacent triangle  ---
!
            TCOR = (T_FRC(2,NT2X)+TABS)/TSPRF
            PCOR = (PG_FRC(2,NT2X)+PATM)/PATM
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( SG_FRC(2,NT2X).LT.EPSL ) THEN
              TORGX = 0.D+0
            ELSE
              TORGX = (SG_FRC(2,NT2X)**7)**(THIRD)
            ENDIF
            SDFG2 = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            VMC2 = SG_FRC(2,NT2X)
            FCG2 = YG_FRC(NT2X,NSL)/(VMC2+SMALL)
            DG2 = TORGX*SG_FRC(2,NT2X)*SDFG2
!
!---        Interfacial fault aperture  ---
!
            INDX = -1
            APX = DIFMN( APM_FRC(2,NT1X),APM_FRC(2,NT2X),
     &        DFF1X,DFF1X,ZERO,INDX )
!
!---        Interfacial diffusion/dispersion coefficient  ---
!
            INDX = 16
            DGX = DIFMN(DG1,DG2,DFF1X,DFF2X,UFFG(1,NCX),INDX)
            DGX = AFF_FRC(NCX)*APX*(DGX+DPGX)/DFF_FRC(NCX)
            FGX = AFF_FRC(NCX)*APX*UFFG(1,NCX)
            IF( MOD(ISLC(23),100)/10.EQ.1 ) FGX = 0.D+0
!
!---        Patankar solute transport  ---
!
            AGX = MAX(-FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
            A1 = (AGX+FGX)*FCG1
            A2 = AGX*FCG2
!
!---        Banded solver  ---
!
            IF( ILES.EQ.1 ) THEN
              MCOL = MT2X
              MROW = MT1X-MT2X+MDT
              ALU(MDT,MT1X) = ALU(MDT,MT1X) + A1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) - A2
!
!---        SPLib solver  ---
!
            ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
              DLU(MCD) = DLU(MCD) + A1
              MROW = KLUC_FRC(IXP_FRC(NT1X),MA)
              MA = MA + 1
              DLU(MROW) = DLU(MROW) - A2
#ifdef petsc
!
!---        PETSc solver  ---
!
            ELSEIF( ILES.EQ.5 ) THEN
              CALL MatSetValues(T_MAT,1,IROW,1,IROW,A1,
     &                        ADD_VALUES,IERR)
              ICOL = petsc_offset(MT2X)
              A2 = -A2
              CALL MatSetValues(T_MAT,1,IROW,1,ICOL,A2,
     &                        ADD_VALUES,IERR)
#endif
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBG_FRC_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBL_FRC_EOR( NSL )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Loads the matrix elements and solution vector for the aqueous 
!     convective-dispersive mass transport equation for faults.
!
!     The Jacobian matrix is initially configured assuming zero-flux
!     boundary conditions.  The matrix is then updated for other
!     user-specified boundary conditions.
!
!     Matrix elements are stored in the array ALU.
!     Elements for the right-hand-side are stored in the array BLU.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 30 November 2020
!
#ifdef petsc
!
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: IROW,ICOL
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SJCBL_FRC_EOR'
!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
        DO NT1X = IP_FRC(1,NFX),IP_FRC(2,NFX)
          IF( IXP_FRC(NT1X).EQ.0 ) CYCLE
          MT1X = NFLD - NXP + IXP_FRC(NT1X)
!
!---      Storage terms  ---
!
          VOLX = APM_FRC(2,NT1X)*AF_FRC(NT1X)
          SC = VOLX*DTI
          MA = 1
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
            ALU(MDT,MT1X) = ALU(MDT,MT1X) + SC
!
!---      SPLib or Lis solvers  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MCD = KLUC_FRC(IXP_FRC(NT1X),MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + SC
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            IROW = petsc_offset(MT1X)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,SC,
     &                      ADD_VALUES,IERR)
#endif  
          ENDIF
!
!---      Molecular diffusion coefficients at triangle  ---
!
          TCOR = (T_FRC(2,NT1X)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SL_FRC(2,NT1X).LT.EPSL ) THEN
            TORLX = 0.D+0
          ELSE
            TORLX = (SL_FRC(2,NT1X)**7)**(THIRD)
          ENDIF
          SDFL1 = SMDL(NSL)*TCOR*(VISRL/VISL_FRC(2,NT1X))
          VMC1 = SL_FRC(2,NT1X)
          FCL1 = YL_FRC(NT1X,NSL)/(VMC1+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DL1 = SDCL_FRC(1,NT1X,NSL)*SDCL_FRC(2,NT1X,NSL)*
     &        EXP(VMC1*SDCL_FRC(3,NT1X,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DL1 = TORLX*VMC1*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DL1 = SDCL_FRC(1,NT1X,NSL)*SDCL_FRC(2,NT1X,NSL)*
     &        VMC1**SDCL_FRC(3,NT1X,NSL)
          ELSE
            DL1 = TORLX*VMC1*SDFL1
          ENDIF
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            MT2X = NFLD - NXP + IXP_FRC(NT2X)
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Neglect hydrodynamic dispersion  ---
!
            DPLX = 0.D+0
!
!---        Molecular diffusion coefficients at adjacent triangle  ---
!
            TCOR = (T_FRC(2,NT2X)+TABS)/TSPRF
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( SL_FRC(2,NT2X).LT.EPSL ) THEN
              TORLX = 0.D+0
            ELSE
              TORLX = (SL_FRC(2,NT2X)**7)**(THIRD)
            ENDIF
            SDFL2 = SMDL(NSL)*TCOR*(VISRL/VISL_FRC(2,NT2X))
            VMC2 = SL_FRC(2,NT2X)
            FCL2 = YL_FRC(NT2X,NSL)/(VMC2+SMALL)
            IF( IEDL(NSL).EQ.2 ) THEN
              DL2 = SDCL_FRC(1,NT2X,NSL)*SDCL_FRC(2,NT2X,NSL)*
     &          EXP(VMC2*SDCL_FRC(3,NT2X,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DL2 = TORLX*VMC2*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DL2 = SDCL_FRC(1,NT2X,NSL)*SDCL_FRC(2,NT2X,NSL)*
     &          VMC2**SDCL_FRC(3,NT2X,NSL)
            ELSE
              DL2 = TORLX*VMC2*SDFL2
            ENDIF
!
!---        Interfacial fault aperture  ---
!
            INDX = -1
            APX = DIFMN( APM_FRC(2,NT1X),APM_FRC(2,NT2X),
     &        DFF1X,DFF1X,ZERO,INDX )
!
!---        Interfacial diffusion/dispersion coefficient  ---
!
            INDX = 16
            DLX = DIFMN(DL1,DL2,DFF1X,DFF2X,UFFL(1,NCX),INDX)
            DLX = AFF_FRC(NCX)*APX*(DLX+DPLX)/DFF_FRC(NCX)
            FLX = AFF_FRC(NCX)*APX*UFFL(1,NCX)
            IF( MOD(ISLC(23),10).EQ.1 ) FLX = 0.D+0
!
!---        Patankar solute transport  ---
!
            ALX = MAX(-FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
            A1 = (ALX+FLX)*FCL1
            A2 = ALX*FCL2
!
!---        Banded solver  ---
!
            IF( ILES.EQ.1 ) THEN
              MCOL = MT2X
              MROW = MT1X-MT2X+MDT
              ALU(MDT,MT1X) = ALU(MDT,MT1X) + A1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) - A2
!
!---        SPLib solver  ---
!
            ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
              DLU(MCD) = DLU(MCD) + A1
              MROW = KLUC_FRC(IXP_FRC(NT1X),MA)
              MA = MA + 1
              DLU(MROW) = DLU(MROW) - A2
#ifdef petsc
!
!---        PETSc solver  ---
!
            ELSEIF( ILES.EQ.5 ) THEN
              CALL MatSetValues(T_MAT,1,IROW,1,IROW,A1,
     &                        ADD_VALUES,IERR)
              ICOL = petsc_offset(MT2X)
              A2 = -A2
              CALL MatSetValues(T_MAT,1,IROW,1,ICOL,A2,
     &                        ADD_VALUES,IERR)
#endif
            ENDIF
          ENDDO
!
!---      Solution vector  ---
!
          BLU(MT1X) = BLU(MT1X) + CO_FRC(NT1X,NSL)*SC
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBL_FRC_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBN_FRC_EOR( NSL )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Loads the matrix elements and solution vector for the nonaqueous
!     liquid convective-dispersive mass transport equation for faults.
!
!     The Jacobian matrix is initially configured assuming zero-flux
!     boundary conditions.  The matrix is then updated for other
!     user-specified boundary conditions.
!
!     Matrix elements are stored in the array ALU.
!     Elements for the right-hand-side are stored in the array BLU.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 30 November 2020
!
#ifdef petsc
!
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE PARM_FRC
      USE NAPL
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
      USE FDVP
      USE FDVN_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: IROW,ICOL
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SJCBN_FRC_EOR'
!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
        DO NT1X = IP_FRC(1,NFX),IP_FRC(2,NFX)
          IF( IXP_FRC(NT1X).EQ.0 ) CYCLE
          MT1X = NFLD - NXP + IXP_FRC(NT1X)
          MA = 1
!
!---      SPLib or Lis solvers  ---
!
          IF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MCD = KLUC_FRC(IXP_FRC(NT1X),MA)
            MA = MA + 1
          ENDIF
!
!---      Molecular diffusion coefficients at triangle  ---
!
          TCOR = (T_FRC(2,NT1X)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( PORD(2,NT1X)*SN_FRC(2,NT1X).LT.EPSL ) THEN
            TORNX = 0.D+0
          ELSE
            TORLX = (PORD(2,NT1X)*(SN_FRC(2,NT1X)**7))**(THIRD)
          ENDIF
          SDFN1 = SMDN(NSL)*TCOR*(VISRN/VISN_FRC(2,NT1X))
          VMC1 = PORD(2,NT1X)*SN_FRC(2,NT1X)
          FCN1 = YN_FRC(NT1X,NSL)/(VMC1+SMALL)
          DN1 = TORNX*VMC1*SDFN1
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            MT2X = NFLD - NXP + IXP_FRC(NT2X)
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Neglect hydrodynamic dispersion  ---
!
            DPNX = 0.D+0
!
!---        Molecular diffusion coefficients at adjacent triangle  ---
!
            TCOR = (T_FRC(2,NT2X)+TABS)/TSPRF
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( PORD(2,NT2X)*SN_FRC(2,NT2X).LT.EPSL ) THEN
              TORNX = 0.D+0
            ELSE
              TORNX = (PORD(2,NT2X)*(SN_FRC(2,NT2X)**7))**(THIRD)
            ENDIF
            SDFN2 = SMDN(NSL)*TCOR*(VISRN/VISN_FRC(2,NT2X))
            VMC2 = PORD(2,NT2X)*SN_FRC(2,NT2X)
            FCN2 = YN_FRC(NT2X,NSL)/(VMC2+SMALL)
            DN2 = TORNX*VMC2*SDFN2
!
!---        Interfacial fault aperture  ---
!
            INDX = -1
            APX = DIFMN( APM_FRC(2,NT1X),APM_FRC(2,NT2X),
     &        DFF1X,DFF1X,ZERO,INDX )
!
!---        Interfacial diffusion/dispersion coefficient  ---
!
            INDX = 16
            DNX = DIFMN(DN1,DN2,DFF1X,DFF2X,UFFN(1,NCX),INDX)
            DNX = AFF_FRC(NCX)*APX*(DNX+DPNX)/DFF_FRC(NCX)
            FNX = AFF_FRC(NCX)*APX*UFFN(1,NCX)
            IF( MOD(ISLC(23),10).EQ.1 ) FNX = 0.D+0
!
!---        Patankar solute transport  ---
!
            ANX = MAX(-FNX,ZERO)
     &        + DNX*MAX((ONE-(TENTH*ABS(FNX)/(DNX+SMALL)))**5,ZERO)
            A1 = (ANX+FNX)*FCN1
            A2 = ANX*FCN2
!
!---        Banded solver  ---
!
            IF( ILES.EQ.1 ) THEN
              MCOL = MT2X
              MROW = MT1X-MT2X+MDT
              ALU(MDT,MT1X) = ALU(MDT,MT1X) + A1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) - A2
!
!---        SPLib solver  ---
!
            ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
              DLU(MCD) = DLU(MCD) + A1
              MROW = KLUC_FRC(IXP_FRC(NT1X),MA)
              MA = MA + 1
              DLU(MROW) = DLU(MROW) - A2
#ifdef petsc
!
!---        PETSc solver  ---
!
            ELSEIF( ILES.EQ.5 ) THEN
              CALL MatSetValues(T_MAT,1,IROW,1,IROW,A1,
     &                        ADD_VALUES,IERR)
              ICOL = petsc_offset(MT2X)
              A2 = -A2
              CALL MatSetValues(T_MAT,1,IROW,1,ICOL,A2,
     &                        ADD_VALUES,IERR)
#endif
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBL_FRC_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_FRC_EOR( NSL )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Fault solute sources.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 30 November 2020
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FDVP_FRC
      USE FDVN_FRC
      USE FDVG_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
#ifdef petsc
      INTEGER :: IROW
      REAL*8 :: BUFFER
      PetscErrorCode :: IERR
#endif
      REAL*8 SRX(10+LSOLU+LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORT_FRC_EOR'
!
!---  Loop over fault sources  ---
!
      DO NS = 1,NSR_FRC
        IF( TM.LE.SRC_FRC(1,1,NS) ) CYCLE
        SRX(1) = TM
        IF( ISRM_FRC(NS).EQ.1 ) THEN
          DO N = 2,10+NSOLU+NGC
            SRX(N) = SRC_FRC(N,1,NS)
          ENDDO
        ELSE
          DO M = 2,ISRM_FRC(NS)
            IF( TM.LE.SRC_FRC(1,M,NS) ) THEN
              DTSR = MIN( SRC_FRC(1,M,NS)-TM,DT )
              TFSR = (TM-0.5D+0*DTSR-SRC_FRC(1,M-1,NS))/
     &          (SRC_FRC(1,M,NS)-SRC_FRC(1,M-1,NS))
              DO N = 2,10+NSOLU+NGC
                SRX(N) = SRC_FRC(N,M-1,NS)+TFSR*
     &            (SRC_FRC(N,M,NS)-SRC_FRC(N,M-1,NS))
              ENDDO
              GOTO 110
            ENDIF
          ENDDO
          CYCLE
        ENDIF
  110   CONTINUE
!
!---    Loop over fault source domain  ---
!
        DO NTX = ISRDM_FRC(1,NS),ISRDM_FRC(2,NS)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
          MTX = NFLD - NXP + IXP_FRC(NTX)
          IF( ILES.EQ.1 ) THEN
            MCOL = MTX
            MROW = MDT
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC_FRC(IXP_FRC(NTX),MA)
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            IROW = petsc_offset(MTX)
#endif
          ENDIF
!
!---      Initialize source  ---
!
          SORTX = 0.D+0
!
!---      Fluid mass rate  ---
!
          IF( ISRT(NS)/100.GE.1 .AND. ISRT(NS)/100.LE.3 ) THEN
!
!---        Fluid withdrawl from grid cell  ---
!
            IF( SRX(2).LE.0.D+0 ) THEN
              FMGX = RKG_FRC(2,NTX)/VISG_FRC(2,NTX)
              FMLX = RKL_FRC(2,NTX)/VISL_FRC(2,NTX)
              FMNX = RKN_FRC(2,NTX)/VISN_FRC(2,NTX)
              SFMX = FMGX + FMLX + FMNX
              RPSLX = RHOL_FRC(2,NTX)*PORD_FRC(2,NTX)*SL_FRC(2,NTX)
              RPSGX = RHOG_FRC(2,NTX)*PORD_FRC(2,NTX)*SG_FRC(2,NTX)
              RPSNX = RHON_FRC(2,NTX)*PORD_FRC(2,NTX)*SN_FRC(2,NTX)
              IF( RPSLX.GT.EPSL ) SORTX = SORTX - SRX(2)*
     &            (YL_FRC(NTX,NSL)/RPSLX)*(FMLX/SFMX)
              IF( RPSGX.GT.EPSL ) SORTX = SORTX - SRX(2)*
     &            (YG_FRC(NTX,NSL)/RPSGX)*(FMGX/SFMX)
              IF( RPSNX.GT.EPSL ) SORTX = SORTX - SRX(2)*
     &            (YN_FRC(NTX,NSL)/RPSNX)*(FMNX/SFMX)
            ENDIF
!
!---      Fluid volumetric rate  ---
!
          ELSEIF( ISRT(NS)/100.GE.4 .AND. ISRT(NS)/100.LE.6 ) THEN
!
!---        Fluid withdrawl from grid cell  ---
!
            IF( SRX(2).LE.0.D+0 ) THEN
              FMGX = RKG_FRC(2,NTX)/VISG_FRC(2,NTX)
              FMLX = RKL_FRC(2,NTX)/VISL_FRC(2,NTX)
              FMNX = RKN_FRC(2,NTX)/VISN_FRC(2,NTX)
              SFMX = FMGX + FMLX + FMNX
              RPSLX = PORD_FRC(2,NTX)*SL_FRC(2,NTX)
              RPSGX = PORD_FRC(2,NTX)*SG_FRC(2,NTX)
              RPSNX = PORD_FRC(2,NTX)*SN_FRC(2,NTX)
              IF( RPSLX.GT.EPSL ) SORTX = SORTX - SRX(2)*
     &            (YL_FRC(NTX,NSL)/RPSLX)*(FMLX/SFMX)
              IF( RPSGX.GT.EPSL ) SORTX = SORTX - SRX(2)*
     &            (YG_FRC(NTX,NSL)/RPSGX)*(FMGX/SFMX)
              IF( RPSNX.GT.EPSL ) SORTX = SORTX - SRX(2)*
     &            (YN_FRC(NTX,NSL)/RPSNX)*(FMNX/SFMX)
            ENDIF
!
!---      Pressure sink  ---
!
          ELSEIF( ISRT_FRC(NS).EQ.35 ) THEN
            PGX = PG_FRC(2,NTX) + PATM
            PLX = PL_FRC(2,NTX) + PATM
            PNX = PN_FRC(2,NTX) + PATM
!
!---        Borehole pressuure, radius, and skin factor  ---
!
            PBHX = SRX(4) + PATM
            RBHX = SRX(3)
            SBHX = SRX(2)
!
!---        Fault triangle equivalent radius  ---
!
            RFEX = SQRT( AF_FRC(NTX)/GPI )
!
!---        Well index  ---
!
            WIX = 2.D+0*GPI*PERM_FRC(2,NTX)*APM_FRC(2,NTX)/
     &        (LOG(RFEX/RBHX)+SBHX)
!
!---        Aqueous mass flux into borehole (only negative flux 
!           allowed)  ---
!
            FLX = MIN( (PBHX-PLX),0.D+0 )*WIX*RKL_FRC(2,NTX)*
     &        RHOL_FRC(2,NTX)/VISL_FRC(2,NTX)
            IF( SL_FRC(2,NTX).GT.EPSL ) SORTX = SORTX -
     &        FLX*YL_FRC(NTX,NSL)/(RHOL_FRC(2,NTX)*SL_FRC(2,NTX))
!
!---        Gas mass flux into borehole (only negative flux 
!           allowed)  ---
!
            FGX = MIN( (PBHX-PGX),0.D+0 )*WIX*RKG_FRC(2,NTX)*
     &        RHOG_FRC(2,NTX)/VISG_FRC(2,NTX)
            IF( SG_FRC(2,NTX).GT.EPSL ) SORTX = SORTX - 
     &        FGX*YG_FRC(NTX,NSL)/(RHOG_FRC(2,NTX)*SG_FRC(2,NTX))
!
!---        Nonaqueous liquid mass flux into borehole (only negative
!           flux allowed)  ---
!
            FNX = MIN( (PBHX-PNX),0.D+0 )*WIX*RKN_FRC(2,NTX)*
     &        RHON_FRC(2,NTX)/VISN_FRC(2,NTX)
            IF( SN_FRC(2,NTX).GT.EPSL ) SORTX = SORTX - 
     &        FNX*YN_FRC(NTX,NSL)/(RHON_FRC(2,NTX)*SN_FRC(2,NTX))
!
!---      Solute source  ---
!
          ELSEIF( ISRT_FRC(NS).EQ.-NSL ) THEN
            BLU(MTX) = BLU(MTX) + SRX(4)
!
!---      Solute density source  ---
!
          ELSEIF( ISRT_FRC(NS).EQ.-(NSL+NSOLU) ) THEN
            BLU(MTX) = BLU(MTX) + SRX(4)*APM_FRC(2,NTX)*AF_FRC(NTX)
          ENDIF
        ENDDO
!
!---    Load Jacobian  ---
!
        IF( ILES.EQ.1 ) THEN
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
        ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
          DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef petsc
        ELSEIF( ILES.EQ.5 ) THEN
        CALL MatSetValues(T_MAT,1,IROW,1,IROW,SORTX,ADD_VALUES,IERR)
#endif
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORT_FRC_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP_FRC_EOR( NSL )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Calculates the aqueous- and gas-phase solute mole fractions
!     from user-specified partition coefficients for fault flow.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 30 November 2020
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE GEOM_FRC
      USE FDVP_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SPRP_FRC_EOR'
!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
          IF( IPCL(NSL).EQ.2 ) THEN
            XVS = SL_FRC(2,NTX)*RHOS_FRC(NTX)*(1.D+0-PORD_FRC(2,NTX))*
     &        PCSL_FRC(1,NTX,NSL)
          ELSE
            XVS = RHOS_FRC(NTX)*(1.D+0-PORD_FRC(2,NTX))*
     &        PCSL_FRC(1,NTX,NSL)
          ENDIF
          XVL = SL_FRC(2,NTX)*PORD_FRC(2,NTX)
          XVG = SG_FRC(2,NTX)*PORD_FRC(2,NTX)
          XVN = SN_FRC(2,NTX)*PORD_FRC(2,NTX)
!
!---      Constant gas-aqueous partition coefficient  ---
!
          IF( IPCGL(NSL).EQ.0 ) THEN
            PCGLX = PCGL(1,NSL)
!
!---      Temperature dependent gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.1 ) THEN
            TK = T_FRC(2,NTX)+TABS
            PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &        + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &        + PCGL(5,NSL)*TK**2 )
!
!---      Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.2 ) THEN
            PCGLX = RHOG_FRC(2,NTX)*XGW_FRC(2,NTX)/
     &        (RHOL_FRC(2,NTX)*XLW_FRC(2,NTX))
          ENDIF
          PCGLX = MAX( PCGLX,1.D-20 )
          PCGLX = MIN( PCGLX,1.D+20 )
!
!---      Constant aqueous-nonaqueous liquid partition coefficient  ---
!
          PCLNX = PCLN(1,NSL)
          PCLNX = MAX( PCLNX,1.D-20 )
          PCLNX = MIN( PCLNX,1.D+20 )
!
!---      Phase-volumetric concentration ratios  ---
!
          YVL = 1.D+0/(XVS + XVL + XVG*PCGLX + XVN/PCLNX)
          YVG = PCGLX/(XVS + XVL + XVG*PCGLX + XVN/PCLNX)
          YVN = 1.D+0/((XVS + XVL + XVG*PCGLX + XVN/PCLNX)*PCLNX)
!
!---      Phase mole fractions  ---
!
          YL_FRC(NTX,NSL) = XVL*YVL
          YG_FRC(NTX,NSL) = XVG*YVG
          YN_FRC(NTX,NSL) = XVN*YVN
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPRP_FRC_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRNSC_FM_EOR( NSL )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Solute/reactive species transport fault to matrix transfer
!     functions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 30 November 2020
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE HYST
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
      USE FDVP
      USE FDVN
      USE FDVG
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TRNSC_FM_EOR'
!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fault triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        IF( IXP_FRC(NTX).EQ.0 ) CYCLE
        MTX = NFLD - NXP + NTX
!
!---    Loop over fault triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
          N = INCM_FRC(NCX)
          IZN = IZ(N)
          MCP = IXP(N)
!
!---      Fault side gas diffusion  ---
!
          TCOR = (T_FRC(2,NTX)+TABS)/TSPRF
          PCOR = (PG_FRC(2,NTX)+PATM)/PATM
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SG_FRC(2,NTX).LT.EPSL ) THEN
            TORGX = 0.D+0
          ELSE
            TORGX = (SG_FRC(2,NTX)**7)**(THIRD)
          ENDIF
          SDFGF = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCF = SG_FRC(2,NTX)
          FCGF = YG_FRC(NTX,NSL)/(VMCF+SMALL)
          DGF = TORGX*SG_FRC(2,NTX)*SDFGF
!
!---      Grid-cell side gas diffusion  ---
!
          TCOR = (T(2,N)+TABS)/TSPRF
          PCOR = (PG(2,N)+PATM)/PATM
          SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCP = SG(2,N)*PORD(2,N)
          FCGP = YG(N,NSL)/(VMCP+SMALL)
          DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
          INDX = 16
          DGX = DIFMN(DGF,DGP,DFN_FRC(NCX),DFN_FRC(NCX),ZERO,INDX)
          DGX = AFN_FRC(NCX)*DGX/DFN_FRC(NCX)
          FGX = UFMG(NCX)
!
!---      Patankar solute transport between fault and grid cell  ---
!
          AGFX = MAX(-FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          AGFFX = (AGFX+FGX)*FCGF
          AGFPX = AGFX*FCGP
          AGPX = MAX(FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          AGPPX = (AGPX-FGX)*FCGP
          AGPFX = AGPX*FCGF
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fault equations  ---
!
            MCOL = MTX
            MROW = MCP-MTX+MDT
            ALU(MDT,MTX) = ALU(MDT,MTX) + AGFFX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AGFPX
!
!---        Matrix equations  ---
!
            MCOL = MCP
            MROW = MTX-MCP+MDT
            ALU(MDT,MCP) = ALU(MDT,MCP) + AGPPX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AGPFX
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fault equations  ---
!
            MA = 1
            MCD = KLUC_FRC(IXP_FRC(NTX),MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + AGFFX
            MROW = KLUC_FCM(NCX,1)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AGFPX
!
!---        Matrix equations  ---
!
            MA = 1
            MCD = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + AGPPX
            MROW = KLUC_MCF(NCX,1)
            DLU(MROW) = DLU(MROW) - AGPFX
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
!
!---        Fault equations  ---
!
            IROW = petsc_offset(MTX)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AGFFX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MCP)
            AGFPX = -AGFPX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AGFPX,
     &                      ADD_VALUES,IERR)
!
!---        Matrix equations  ---
!
            IROW = petsc_offset(MCP)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AGPPX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MTX)
            AGPFX = -AGPFX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AGPFX,
     &                      ADD_VALUES,IERR)
#endif
          ENDIF
!
!---      Fault side aqueous diffusion  ---
!
          TCOR = (T_FRC(2,NTX)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SL_FRC(2,NTX).LT.EPSL ) THEN
            TORLX = 0.D+0
          ELSE
            TORLX = (SL_FRC(2,NTX)**7)**(THIRD)
          ENDIF
          SDFLF = SMDL(NSL)*TCOR*(VISRL/VISL_FRC(2,NTX))
          VMCF = SL_FRC(2,NTX)
          FCLF = YL_FRC(NTX,NSL)/(VMCF+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLF = SDCL_FRC(1,NTX,NSL)*SDCL_FRC(2,NTX,NSL)*
     &        EXP(VMCF*SDCL_FRC(3,NTX,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLF = TORLX*VMCF*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLF = SDCL_FRC(1,NTX,NSL)*SDCL_FRC(2,NTX,NSL)*
     &        VMCF**SDCL_FRC(3,NTX,NSL)
          ELSE
            DLF = TORLX*VMCF*SDFLF
          ENDIF
!
!---      Grid-cell side gas diffusion  ---
!
          TCOR = (T(2,N)+TABS)/TSPRF
          SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
          VMCP = SL(2,N)*PORD(2,N)
          FCLP = YL(N,NSL)/(VMCP+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLP = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &        EXP(VMCP*SDCL(3,IZN,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLP = TORL(2,N)*VMCP*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLP = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &        VMCP**SDCL(3,IZN,NSL)
          ELSE
            DLP = TORL(2,N)*VMCP*SDFLP
          ENDIF
          INDX = 16
          DLX = DIFMN(DLF,DLP,DFN_FRC(NCX),DFN_FRC(NCX),ZERO,INDX)
          DLX = AFN_FRC(NCX)*DLX/DFN_FRC(NCX)
          FLX = UFML(NCX)
!
!---      Patankar solute transport between fault and grid celll  ---
!
          ALFX = MAX(-FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          ALFFX = (ALFX+FLX)*FCLF
          APLFX = ALFX*FCLP
          ALPX = MAX(FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          ALPPX = (ALPX-FLX)*FCLP
          ALPFX = ALPX*FCLF
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fault equations  ---
!
            MCOL = MTX
            MROW = MCP-MTX+MDT
            ALU(MDT,MTX) = ALU(MDT,MTX) + ALFFX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - APLFX
!
!---        Matrix equations  ---
!
            MCOL = MCP
            MROW = MTX-MCP+MDT
            ALU(MDT,MCP) = ALU(MDT,MCP) + ALPPX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - ALPFX
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fault equations  ---
!
            MA = 1
            MCD = KLUC_FRC(IXP_FRC(NTX),MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + ALFFX
            MROW = KLUC_FCM(NCX,1)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - APLFX
!
!---        Matrix equations  ---
!
            MA = 1
            MCD = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + ALPPX
            MROW = KLUC_MCF(NCX,1)
            DLU(MROW) = DLU(MROW) - ALPFX
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
!
!---        Fault equations  ---
!
            IROW = petsc_offset(MTX)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,ALFFX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MCP)
            APLFX = -APLFX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,APLFX,
     &                      ADD_VALUES,IERR)
!
!---        Matrix equations  ---
!
            IROW = petsc_offset(MCP)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,ALPPX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MTX)
            ALPFX = -ALPFX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,ALPFX,
     &                      ADD_VALUES,IERR)
#endif
          ENDIF
!
!---      Fault side nonaqueous diffusion  ---
!
          TCOR = (T_FRC(2,NTX)+TABS)/TSPRF
          PCOR = (PN_FRC(2,NTX)+PATM)/PATM
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SN_FRC(2,NTX)*PORD_FRC(2,NTX).LT.EPSL ) THEN
            TORNX = 0.D+0
          ELSE
            TORNX = (PORD_FRC(2,NTX)*(SN_FRC(2,NTX)**7))**(THIRD)
          ENDIF
          SDFNF = SMDN(NSL)*(TCOR**1.75D+0)/PCOR
          VMCF = SN_FRC(2,NTX)*PORD_FRC(2,NTX)
          FCNF = YN_FRC(NTX,NSL)/(VMCF+SMALL)
          DNF = TORNX*SN_FRC(2,NTX)*PORD_FRC(2,NTX)*SDFNF
!
!---      Grid-cell side gas diffusion  ---
!
          TCOR = (T(2,N)+TABS)/TSPRF
          PCOR = (PN(2,N)+PATM)/PATM
          SDFNP = SMDN(NSL)*(TCOR**1.75D+0)/PCOR
          VMCP = SN(2,N)*PORD(2,N)
          FCNP = YN(N,NSL)/(VMCP+SMALL)
          DNP = TORN(2,N)*(SN(2,N)-SNT(2,N))*PORD(2,N)*SDFNP
          INDX = 16
          DNX = DIFMN(DNF,DNP,DFN_FRC(NCX),DFN_FRC(NCX),ZERO,INDX)
          DNX = AFN_FRC(NCX)*DNX/DFN_FRC(NCX)
          FNX = UFMN(NCX)
!
!---      Patankar solute transport between fault and grid cell  ---
!
          ANFX = MAX(-FNX,ZERO)
     &        + DNX*MAX((ONE-(TENTH*ABS(FNX)/(DNX+SMALL)))**5,ZERO)
          ANFFX = (ANFX+FNX)*FCNF
          ANFPX = ANFX*FCNP
          ANPX = MAX(FNX,ZERO)
     &        + DNX*MAX((ONE-(TENTH*ABS(FNX)/(DNX+SMALL)))**5,ZERO)
          ANPPX = (ANPX-FNX)*FCNP
          ANPFX = ANPX*FCNF
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fault equations  ---
!
            MCOL = MTX
            MROW = MCP-MTX+MDT
            ALU(MDT,MTX) = ALU(MDT,MTX) + ANFFX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - ANFPX
!
!---        Matrix equations  ---
!
            MCOL = MCP
            MROW = MTX-MCP+MDT
            ALU(MDT,MCP) = ALU(MDT,MCP) + ANPPX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - ANPFX
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fault equations  ---
!
            MA = 1
            MCD = KLUC_FRC(IXP_FRC(NTX),MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + ANFFX
            MROW = KLUC_FCM(NCX,1)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - ANFPX
!
!---        Matrix equations  ---
!
            MA = 1
            MCD = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + ANPPX
            MROW = KLUC_MCF(NCX,1)
            DLU(MROW) = DLU(MROW) - ANPFX
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
!
!---        Fault equations  ---
!
            IROW = petsc_offset(MTX)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,ANFFX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MCP)
            ANFPX = -ANFPX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,ANFPX,
     &                      ADD_VALUES,IERR)
!
!---        Matrix equations  ---
!
            IROW = petsc_offset(MCP)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,ANPPX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MTX)
            ANPFX = -ANPFX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,ANPFX,
     &                      ADD_VALUES,IERR)
#endif
          ENDIF
!
!---    Loop over fault triangle to grid cell connections  ---
!
        ENDDO
!
!---  Loop over fault triangles  ---
!
      ENDDO
!
!---  Loop over faults  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TRNSC_FM_EOR group  ---
!
      RETURN
      END


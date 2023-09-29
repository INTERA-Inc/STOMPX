!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC_FRC_GT
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
!     Geothermal Mode (STOMP-GT)
!
!     Compute initial solute fracture concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 29 June 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE PORMED
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 PCSLX(5),SDCLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      IF( IEQC.EQ.0 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CISC_FRC_GT'
!
!---  Loop over solutes  ---
!
      DO NSL = 1,NSOLU
!
!---    Loop over fractures  ---
!
        DO NFX = 1,NF_FRC
!
!---      Loop over fracture triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---        Loop over fracture triangle to grid cell connections,
!           to determine area-weighted average rock properties  ---
!
            AFNX = 0.D+0
            RHOSX = 0.D+0
            PORDX = 0.D+0
            DO M = 1,5
              PCSLX(M) = 0.D+0
            ENDDO
            DO M = 1,3
              SDCLX(M) = 0.D+0
            ENDDO
            DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
              IZN = IZ(INCM_FRC(NCX))
              AFNX = AFNX + AFN_FRC(NCX)
              RHOSX = RHOSX + AFN_FRC(NCX)*RHOS(IZN)
              PORDX = PORDX + AFN_FRC(NCX)*POR(2,IZN)
              DO M = 1,5
                PCSLX(M) = PCSLX(M) + AFN_FRC(NCX)*PCSL(M,IZN,NSL)
              ENDDO
              DO M = 1,3
                SDCLX(M) = SDCLX(M) + AFN_FRC(NCX)*SDCL(M,IZN,NSL)
              ENDDO
            ENDDO
            RHOS_FRC(NTX) = RHOSX/(AFNX+SMALL)
            PORD_FRC(2,NTX) = PORDX/(AFNX+SMALL)
            DO M = 1,5
              PCSL_FRC(M,NTX,NSL) = PCSL_FRC(M,NTX,NSL) + 
     &          PCSLX(M)/(AFNX+SMALL)
            ENDDO
            DO M = 1,3
              SDCL_FRC(M,NTX,NSL) = SDCL_FRC(M,NTX,NSL) + 
     &          SDCLX(M)/(AFNX+SMALL)
            ENDDO
!
!---        Stationary solute  ---
!
            IF( IEDL(NSL).EQ.4  ) THEN
              XVS = 0.D+0
            ELSEIF( IPCL(NSL).EQ.2 ) THEN
              XVS = SL_FRC(2,NTX)*RHOS_FRC(NTX)*
     &          (1.D+0-PORD_FRC(2,NTX))*PCSL_FRC(1,NTX,NSL)
            ELSE
              XVS = RHOS_FRC(NTX)*PCSL_FRC(1,NTX,NSL)*
     &          (1.D+0-PORD_FRC(2,NTX))
            ENDIF
            XVL = SL_FRC(2,NTX)
            XVG = SG_FRC(2,NTX)
!
!---        Stationary solute  ---
!
            IF( IEDL(NSL).EQ.4  ) THEN
              PCGLX = 1.D+0
!
!---        Constant gas-aqueous partition coefficient  ---
!
            ELSEIF( IPCGL(NSL).EQ.0 ) THEN
              PCGLX = PCGL(1,NSL)
!
!---        Temperature dependent gas-aqueous partition coefficient  ---
!
            ELSEIF( IPCGL(NSL).EQ.1 ) THEN
              TK = T_FRC(2,NTX)+TABS
              PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &          + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &          + PCGL(5,NSL)*TK**2 )
!
!---        Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
            ELSEIF( IPCGL(NSL).EQ.2 ) THEN
              PCGLX = RHOG_FRC(2,NTX)*XGW_FRC(2,NTX)/
     &          (RHOL_FRC(2,NTX)*XLW_FRC(2,NTX))
            ENDIF
            PCGLX = MAX( PCGLX,1.D-20 )
            PCGLX = MIN( PCGLX,1.D+20 )
!
!---        Phase-volumetric concentration ratios  ---
!
            YVL = 1.D+0/(XVS + XVL + XVG*PCGLX)
            YVG = 1.D+0/((XVS + XVL)/PCGLX + XVG)
!
!---        Phase mole fractions  ---
!
            YL_FRC(NTX,NSL) = XVL*YVL
            YG_FRC(NTX,NSL) = XVG*YVG
!
!---        Phase-volumetric concentration ratios  ---
!
            IF( ICT_FRC(NFX,NSL).EQ.2 ) THEN
              C_FRC(NTX,NSL) = C_FRC(NTX,NSL)*(XVS + XVL + XVG*PCGLX)
            ELSEIF( ICT_FRC(NFX,NSL).EQ.-1 ) THEN
              C_FRC(NTX,NSL) = C_FRC(NTX,NSL)*RHOS_FRC(NTX)
            ELSEIF( ICT_FRC(NFX,NSL).EQ.3 ) THEN
              C_FRC(NTX,NSL) = C_FRC(NTX,NSL)*((XVS + XVL)/PCGLX + XVG)
            ELSEIF( ICT_FRC(NFX,NSL).EQ.4 ) THEN
              C_FRC(NTX,NSL) = C_FRC(NTX,NSL)/(LOG(2.D+0)/HLF(NSL))
     &          /VOL_FRC(NTX)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CISC_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC_BH_GT
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
!     Geothermal Mode (STOMP-GT)
!
!     Compute initial solute borehole concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 6 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE PARM_BH
      USE GEOM_BH
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/CISC_BH_GT'
!
!---  Loop over solutes  ---
!
      DO NSL = 1,NSOLU
!
!---    Loop over boreholes  ---
!
        DO NBH = 1,N_BH
!
!---      Skip for boundary condition type boreholes  ---
!
          IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---      Loop over borehole nodes  ---
!
          DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
            IZN = IZ_BH(NBN)
            INVX = INV_BH(NBN)
!
!---        Pipe flow  ---
!
            IF( IS_BH(INVX).GT.10  ) THEN
              XVS = 0.D+0
            ELSEIF( IPCL(NSL).EQ.2 ) THEN
              XVS = SL_BH(2,NBN)*RHOS(IZN)*PCSL(1,IZN,NSL)*
     &          (1.D+0-PORD_BH(2,NBN))
!
!---        Stationary solute  ---
!
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              XVS = 0.D+0
            ELSE
              XVS = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORD_BH(2,NBN))
            ENDIF
            XVL = SL_BH(2,NBN)*PORD_BH(2,NBN)
            XVG = SG_BH(2,NBN)*PORD_BH(2,NBN)
!
!---        Stationary solute  ---
!
            IF( IEDL(NSL).EQ.4 ) THEN
              PCGLX = 1.D+0
!
!---        Constant gas-aqueous partition coefficient  ---
!
            ELSEIF( IPCGL(NSL).EQ.0 ) THEN
              PCGLX = PCGL(1,NSL)
!
!---        Temperature dependent gas-aqueous partition coefficient  ---
!
            ELSEIF( IPCGL(NSL).EQ.1 ) THEN
              TK = T_BH(2,NBN)+TABS
              PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &          + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &          + PCGL(5,NSL)*TK**2 )
!
!---        Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
            ELSEIF( IPCGL(NSL).EQ.2 ) THEN
              PCGLX = RHOG_BH(2,NBN)*XGW_BH(2,NBN)/
     &          (RHOL_BH(2,NBN)*XLW_BH(2,NBN))
            ENDIF
            PCGLX = MAX( PCGLX,1.D-20 )
            PCGLX = MIN( PCGLX,1.D+20 )
!
!---        Phase-volumetric concentration ratios  ---
!
            YVL = 1.D+0/(XVS + XVL + XVG*PCGLX)
            YVG = 1.D+0/((XVS + XVL)/PCGLX + XVG)
!
!---        Phase mole fractions  ---
!
            YL_BH(NBN,NSL) = XVL*YVL
            YG_BH(NBN,NSL) = XVG*YVG
!
!---        Phase-volumetric concentration ratios  ---
!
            IF( ICT_BH(NBH,NSL).EQ.2 ) THEN
              C_BH(NBN,NSL) = C_BH(NBN,NSL)*(XVS + XVL + XVG*PCGLX)
            ELSEIF( ICT_BH(NBH,NSL).EQ.-1 ) THEN
              C_BH(NBN,NSL) = C_BH(NBN,NSL)*RHOS(IZN)*
     &          (1.D+0-PORD_BH(2,NBN))
            ELSEIF( ICT_BH(NBH,NSL).EQ.3 ) THEN
              C_BH(NBN,NSL) = C_BH(NBN,NSL)*((XVS + XVL)/PCGLX + XVG)
            ELSEIF( ICT_BH(NBH,NSL).EQ.4 ) THEN
              C_BH(NBN,NSL) = C_BH(NBN,NSL)/(LOG(2.D+0)/HLF(NSL))
     &          /VOL_BH(NBN)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CISC_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFX_FRC_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Compute solute fracture transport flux for gas and aqueous phases.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 October 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE GEOM_FRC
      USE FLUX_FRC
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFX_FRC_GT'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fracture triangles  ---
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
!---      Loop over fracture triangle connections  ---
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
!---        Interfacial fracture aperture  ---
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
!---      Loop over fracture triangle connections  ---
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
!---        Interfacial fracture aperture  ---
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
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFX_FRC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFX_BH_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Compute solute borehole transport flux for gas and aqueous phases.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 October 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/SFX_BH_GT'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN1 = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---      Gas contributions, 
!         molecular diffusion coefficients at borehole node  ---
!
          TCOR = (T_BH(2,NBN1)+TABS)/TSPRF
          PCOR = (PG_BH(2,NBN1)+PATM)/PATM
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SG_BH(2,NBN1).LT.EPSL ) THEN
            TORGX = 0.D+0
          ELSE
            TORGX = (SG_BH(2,NBN1)**7)**(THIRD)
          ENDIF
          SDFG1 = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMC1 = SG_BH(2,NBN1)*PORD_BH(2,NBN1)
          FCG1 = YG_BH(NBN1,NSL)/(VMC1+SMALL)
          DG1 = TORGX*SG_BH(2,NBN1)*SDFG1
!
!---      Loop over borehole node to node connections ---
!
          DO NCX = IPB_BH(1,NBN1),IPB_BH(2,NBN1)
            NBN2 = IBCM_BH(NCX)
            MBN2 = NFLD - NXP + NT_FRC + NBN2
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---        Neglect hydrodynamic dispersion  ---
!
            DPGX = 0.D+0
!
!---        Molecular diffusion coefficients at adjacent triangle  ---
!
            TCOR = (T_BH(2,NBN2)+TABS)/TSPRF
            PCOR = (PG_BH(2,NBN2)+PATM)/PATM
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( SG_BH(2,NBN2).LT.EPSL ) THEN
              TORGX = 0.D+0
            ELSE
              TORGX = (SG_BH(2,NBN2)**7)**(THIRD)
            ENDIF
            SDFG2 = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            VMC2 = SG_BH(2,NBN2)*PORD_BH(2,NBN2)
            FCG2 = YG_BH(NBN2,NSL)/(VMC2+SMALL)
            DG2 = TORGX*SG_BH(2,NBN2)*SDFG2
!
!---        Interfacial fracture aperture  ---
!
            INDX = -1
            INV1 = INV_BH(NBN1)
            IF( IS_BH(INV1).GT.10 ) THEN
              RBX = MAX( PAR_BH(3,INV1),1.D-20 )
            ELSE
              RBX = MAX( PAR_BH(2,INV1),1.D-20 )
            ENDIF
            AP1X = GPI*(RBX**2)
            INV2 = INV_BH(NBN2)
            IF( IS_BH(INV2).GT.10 ) THEN
              RBX = MAX( PAR_BH(3,INV2),1.D-20 )
            ELSE
              RBX = MAX( PAR_BH(2,INV2),1.D-20 )
            ENDIF
            AP2X = GPI*(RBX**2)
            APX = DIFMN( AP1X,AP2X,DBB1X,DBB1X,ZERO,INDX )
!
!---        Interfacial diffusion/dispersion coefficient  ---
!
            INDX = 16
            DGX = DIFMN(DG1,DG2,DBB1X,DBB2X,UBBG(1,NCX),INDX)
            DGX = APX*(DGX+DPGX)/DBB_BH(NCX)
            FGX = APX*UBBG(1,NCX)
            IF( MOD(ISLC(23),10).EQ.1 ) FGX = 0.D+0
!
!---        Patankar solute transport  ---
!
            AGX = MAX(FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
            AGPX = MAX(-FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
            UBBC(NCX,NSL) = (C_BH(NBN1,NSL)*AGX*FCG1 - 
     &        C_BH(NBN2,NSL)*AGPX*FCG2)
          ENDDO
!
!---      Aqueous contributions,
!         molecular diffusion coefficients at triangle  ---
!
          TCOR = (T_BH(2,NBN1)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SL_BH(2,NBN1).LT.EPSL ) THEN
            TORLX = 0.D+0
          ELSE
            TORLX = (SL_BH(2,NBN1)**7)**(THIRD)
          ENDIF
          SDFL1 = SMDL(NSL)*TCOR*(VISRL/VISL_BH(2,NBN1))
          VMC1 = SL_BH(2,NBN1)*PORD_BH(2,NBN1)
          FCL1 = YL_BH(NBN1,NSL)/(VMC1+SMALL)
          IZN1 = IZ_BH(NBN1)
!
!---      Pipe flow  ---
!
          IF( IS_BH(INV1).GT.10 ) THEN
            DL1 = TORLX*VMC1*SDFL1
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DL1 = SDCL(1,IZN1,NSL)*SDCL(2,IZN1,NSL)*
     &        EXP(VMC1*SDCL(3,IZN1,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DL1 = TORLX*VMC1*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DL1 = SDCL(1,IZN1,NSL)*SDCL(2,IZN1,NSL)*
     &        VMC1**SDCL(3,IZN1,NSL)
          ELSE
            DL1 = TORLX*VMC1*SDFL1
          ENDIF
!
!---      Loop over borehole node to node connections ---
!
          DO NCX = IPB_BH(1,NBN1),IPB_BH(2,NBN1)
            NBN2 = IBCM_BH(NCX)
            MBN2 = NFLD - NXP + NT_FRC + NBN2
            IZN2 = IZ_BH(NBN2)
            INV2 = INV_BH(NBN2)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---        Neglect hydrodynamic dispersion  ---
!
            DPLX = 0.D+0
!
!---        Molecular diffusion coefficients at adjacent triangle  ---
!
            TCOR = (T_BH(2,NBN2)+TABS)/TSPRF
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( SL_BH(2,NBN2).LT.EPSL ) THEN
              TORLX = 0.D+0
            ELSE
              TORLX = (SL_BH(2,NBN2)**7)**(THIRD)
            ENDIF
            SDFL2 = SMDL(NSL)*TCOR*(VISRL/VISL_BH(2,NBN2))
            VMC2 = SL_BH(2,NBN2)*PORD_BH(2,NBN2)
            FCL2 = YL_BH(NBN2,NSL)/(VMC2+SMALL)
!
!---        Pipe flow  ---
!
            IF( IS_BH(INV2).GT.10 ) THEN
              DL2 = TORLX*VMC2*SDFL2
            ELSEIF( IEDL(NSL).EQ.2 ) THEN
              DL2 = SDCL(1,IZN2,NSL)*SDCL(2,IZN2,NSL)*
     &          EXP(VMC2*SDCL(3,IZN2,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DL2 = TORLX*VMC2*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DL2 = SDCL(1,IZN2,NSL)*SDCL(2,IZN2,NSL)*
     &          VMC2**SDCL(3,IZN2,NSL)
            ELSE
              DL2 = TORLX*VMC2*SDFL2
            ENDIF
!
!---        Interfacial fracture aperture  ---
!
            INDX = -1
            INV1 = INV_BH(NBN1)
            IF( IS_BH(INV1).GT.10 ) THEN
              RBX = MAX( PAR_BH(3,INV1),1.D-20 )
            ELSE
              RBX = MAX( PAR_BH(2,INV1),1.D-20 )
            ENDIF
            AP1X = GPI*(RBX**2)
            INV2 = INV_BH(NBN2)
            IF( IS_BH(INV2).GT.10 ) THEN
              RBX = MAX( PAR_BH(3,INV2),1.D-20 )
            ELSE
              RBX = MAX( PAR_BH(2,INV2),1.D-20 )
            ENDIF
            AP2X = GPI*(RBX**2)
            APX = DIFMN( AP1X,AP2X,DBB1X,DBB1X,ZERO,INDX )
!
!---        Interfacial diffusion/dispersion coefficient  ---
!
            INDX = 16
            DLX = DIFMN(DL1,DL2,DBB1X,DBB2X,UBBL(1,NCX),INDX)
            DLX = APX*(DLX+DPLX)/DBB_BH(NCX)
            FLX = APX*UBBL(1,NCX)
            IF( MOD(ISLC(23),10).EQ.1 ) FLX = 0.D+0
!
!---        Patankar solute transport  ---
!
            ALX = MAX(FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
            ALPX = MAX(-FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
            UBBC(NCX,NSL) = UBBC(NCX,NSL)+ (C_BH(NBN1,NSL)*ALX*FCL1 - 
     &        C_BH(NBN2,NSL)*ALPX*FCL2)
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFX_BH_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFX_FM_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Compute solute fracture to matrix transport flux for gas 
!     and aqueous phases.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 October 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE HYST
      USE GRID
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/SFX_FM_GT'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fracture triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        IF( IXP_FRC(NTX).EQ.0 ) CYCLE
        MTX = NFLD - NXP + IXP_FRC(NTX)
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
          N = INCM_FRC(NCX)
          IZN = IZ(N)
!
!---      Gas contributions, 
!         fracture side gas diffusion  ---
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
!---      Patankar solute transport between fracture and grid cell  ---
!
          AGX = MAX(FGX,ZERO)
     &      + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          AGPX = MAX(-FGX,ZERO)
     &      + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          UFMC(NCX,NSL) = (C_FRC(NTX,NSL)*AGX*FCGF - 
     &      C(N,NSL)*AGPX*FCGP)
!
!---      Aqueous contributions,
!         fracture side aqueous diffusion  ---
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
!---      Patankar solute transport between fracture and grid cell  ---
!
          ALX = MAX(FLX,ZERO)
     &      + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          ALPX = MAX(-FLX,ZERO)
     &      + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          UFMC(NCX,NSL) = UFMC(NCX,NSL)+ (C_FRC(NTX,NSL)*ALX*FCLF - 
     &        C(N,NSL)*ALPX*FCLP)
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFX_FM_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFX_BF_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Compute solute borehole to fracture transport flux for gas 
!     and aqueous phases.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 October 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/SFX_BF_GT'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole node to fracture triangle connections  ---
!
        DO NCX  = ID_BH(5,NBH),ID_BH(6,NBH)
          NTX = IBHT_FRC(NCX)
          NBN = IBHN_FRC(NCX)
          INVX = INV_BH(NBN)
          IZN = IZ_BH(NBN)
          MTX = NFLD - NXP + NTX
          MBN = NFLD - NXP + NT_FRC + NBN 
!
!---      Borehole node centroid  ---
!
          XBHX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
          YBHX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
          ZBHX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
!
!---      Distance between borehole centroid and borehole-fracture
!         intersection  ---
!
          DBF1X = SQRT((XBHX-XP_BF(NCX))**2 + (YBHX-YP_BF(NCX))**2 +
     &      (ZBHX-ZP_BF(NCX))**2)
!
!---      Distance between fracture centroid and borehole-fracture
!         intersection  ---
!
          DBF2X = SQRT((XP_FRC(NTX)-XP_BF(NCX))**2 + 
     &      (YP_FRC(NTX)-YP_BF(NCX))**2 + (ZP_FRC(NTX)-ZP_BF(NCX))**2)
!
!---      Distance between borehole node centroid and fracture
!         triangle centroid  ---
!
          DBFX = MAX( DBF1X+DBF2X,1.D-6 )
!
!---      Gas contributions, 
!         borehole side gas diffusion  ---
!
          TCOR = (T_BH(2,NBN)+TABS)/TSPRF
          PCOR = (PG_BH(2,NBN)+PATM)/PATM
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SG_BH(2,NBN).LT.EPSL ) THEN
            TORGX = 0.D+0
          ELSE
            TORGX = (SG_BH(2,NBN)**7)**(THIRD)
          ENDIF
          SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCB = SG_BH(2,NBN)*PORD_BH(2,NBN)
          FCGB = YG_BH(NBN,NSL)/(VMCB+SMALL)
          DGB = TORGX*SG_BH(2,NBN)*PORD_BH(2,NBN)*SDFGB
!
!---      Fracture side gas diffusion  ---
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
!---      Interfacial gas diffusion between borehole and fracture  ---
!
          INDX = 16
          DGX = DIFMN(DGB,DGF,DBF1X,DBF2X,ZERO,INDX)
          APX = GPI*2.D+1*PAR_BH(2,INVX)*APM_FRC(2,NTX)
          DGX = APX*DGX/DBFX
          FGX = UBFG(1,NCX)
!
!---      Patankar solute transport between borehole node and
!         fracture triangle  ---
!
          AGX = MAX(FGX,ZERO)
     &      + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          AGPX = MAX(-FGX,ZERO)
     &      + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          UBFC(NCX,NSL) = (C_BH(NBN,NSL)*AGX*FCGB - 
     &      C_FRC(NTX,NSL)*AGPX*FCGP)
!
!---      Aqueous contributions,
!         borehole side aqueous diffusion  ---
!
          TCOR = (T_BH(2,NBN)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SL_BH(2,NBN).LT.EPSL ) THEN
            TORLX = 0.D+0
          ELSE
            TORLX = (SL_BH(2,NBN)**7)**(THIRD)
          ENDIF
          SDFLB = SMDL(NSL)*TCOR*(VISRL/VISL_BH(2,NBN))
          VMCB = SL_BH(2,NBN)*PORD_BH(2,NBN)
          FCLB = YL_BH(NBN,NSL)/(VMCB+SMALL)
!
!---      Pipe flow  ---
!
          IF( IS_BH(INVX).GT.10 ) THEN
            DLB = TORLX*VMCB*SDFLB
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &        EXP(VMCB*SDCL(3,IZN,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLB = TORLX*VMCB*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &        VMCB**SDCL(3,IZN,NSL)
          ELSE
            DLB = TORLX*VMCB*SDFLB
          ENDIF
!
!---      Fracture side aqueous diffusion  ---
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
!---      Interfacial aqueous diffusion between borehole and 
!         fracture  ---
!
          INDX = 16
          DLX = DIFMN(DLB,DLF,DBF1X,DBF2X,ZERO,INDX)
          APX = GPI*2.D+1*PAR_BH(2,INVX)*APM_FRC(2,NTX)
          DLX = APX*DLX/DBFX
          FLX = UBFL(1,NCX)
!
!---      Patankar solute transport between borehole node and 
!         fracture triangle  ---
!
          ALX = MAX(FLX,ZERO)
     &      + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          ALPX = MAX(-FLX,ZERO)
     &      + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          UBFC(NCX,NSL) = UBFC(NCX,NSL) + (C_BH(NBN,NSL)*ALX*FCLB - 
     &        C_FRC(NTX,NSL)*ALPX*FCLP)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFX_BF_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFX_BM_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Compute solute borehole to matrix transport flux for gas 
!     and aqueous phases.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 October 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE HYST
      USE GRID
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_BH
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/SFX_BM_GT'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---      Field node connected to borehole node  ---
!
          N = IBN_BH(NBN)
!
!---      Borehole node centroid and surface area  ---
!
          XP_BHX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
          YP_BHX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
          ZP_BHX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
          INVX = INV_BH(NBN)
          RBX = MAX( PAR_BH(2,INVX),1.D-20 )    
          AFB_BHX = 2.D+0*GPI*RBX*SQRT( (XP_BH(2,NBN)-XP_BH(1,NBN))**2 + 
     &      (YP_BH(2,NBN)-YP_BH(1,NBN))**2 + 
     &      (ZP_BH(2,NBN)-ZP_BH(1,NBN))**2 )
!
!---      Gas contributions, 
!         borehole side gas diffusion  ---
!
          TCOR = (T_BH(2,NBN)+TABS)/TSPRF
          PCOR = (PG_BH(2,NBN)+PATM)/PATM
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SG_BH(2,NBN).LT.EPSL ) THEN
            TORGX = 0.D+0
          ELSE
            TORGX = (SG_BH(2,NBN)**7)**(THIRD)
          ENDIF
          SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCB = SG_BH(2,NBN)*PORD_BH(2,NBN)
          FCGB = YG_BH(NBN,NSL)/(VMCB+SMALL)
          DGB = TORGX*SG_BH(2,NBN)*PORD_BH(2,NBN)*SDFGB
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
          DGX = DIFMN(DGB,DGP,DBN_BH(NBN),DBN_BH(NBN),ZERO,INDX)
          DGX = AFB_BHX*DGX/DBN_BH(NBN)
          FGX = UBMG(NBN)
!
!---      Patankar solute transport between borehole and grid cell  ---
!
          AGX = MAX(FGX,ZERO)
     &      + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          AGPX = MAX(-FGX,ZERO)
     &      + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          UBMC(NBN,NSL) = (C_BH(NBN,NSL)*AGX*FCGB - 
     &      C(N,NSL)*AGPX*FCGP)
!
!---      Aqueous contributions,
!         borehole side aqueous diffusion  ---
!
          TCOR = (T_BH(2,NBN)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SL_BH(2,NBN).LT.EPSL ) THEN
            TORLX = 0.D+0
          ELSE
            TORLX = (SL_BH(2,NBN)**7)**(THIRD)
          ENDIF
          SDFLB = SMDL(NSL)*TCOR*(VISRL/VISL_BH(2,NBN))
          VMCB = SL_BH(2,NBN)*PORD_BH(2,NBN)
          FCLB = YL_BH(NBN,NSL)/(VMCB+SMALL)
!
!---      Pipe flow  ---
!
          IZN = IZ_BH(NBN)
          IF( IS_BH(INVX).GT.10 ) THEN
            DLB = TORLX*VMCB*SDFLB
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &        EXP(VMCB*SDCL(3,IZN,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLB = TORLX*VMCB*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &        VMCB**SDCL(3,IZN,NSL)
          ELSE
            DLB = TORLX*VMCB*SDFLB
          ENDIF
!
!---      Grid-cell side aqueous diffusion  ---
!
          TCOR = (T(2,N)+TABS)/TSPRF
          SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
          VMCP = SL(2,N)*PORD(2,N)
          FCLP = YL(N,NSL)/(VMCP+SMALL)
          IZN = IZ(N)
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
          DLX = DIFMN(DLB,DLP,DBN_BH(NBN),DBN_BH(NBN),ZERO,INDX)
          DLX = AFB_BHX*DLX/DBN_BH(NBN)
          FLX = UBML(NBN)
!
!---      Patankar solute transport between borehole node and 
!         matrix grid cell  ---
!
          ALX = MAX(FLX,ZERO)
     &      + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          ALPX = MAX(-FLX,ZERO)
     &      + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          UBMC(NBN,NSL) = UBMC(NBN,NSL) + (C_BH(NBN,NSL)*ALX*FCLB - 
     &        C(N,NSL)*ALPX*FCLP)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFX_BM_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXZ_FRC_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Zero solute transport fluxes for fracture flow.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 July 2018.
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
      SUB_LOG(ISUB_LOG) = '/SFXZ_FRC_GT'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fracture triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---      Loop over fracture triangle connections  ---
!
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            UFFC(NCX,NSL) = 0.D+0
          ENDDO
!
!---      Loop over fracture triangle to grid cell connections  ---
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
!---  End of SFXZ_FRC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXZ_BH_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Zero solute transport fluxes for borehole flow.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_BH
      USE GEOM_BH
      USE FLUX_BH
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
      SUB_LOG(ISUB_LOG) = '/SFXZ_BH_GT'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---      Loop over borehole node to borehole node connections ---
!
          DO NCX = IPB_BH(1,NBN),IPB_BH(2,NBN)
            UBBC(NCX,NSL) = 0.D+0
          ENDDO
          UBMC(NBN,NSL) = 0.D+0
!
!---      Loop over borehole node to fracture triangle connections  ---
!
          DO NCX  = ID_BH(5,NBH),ID_BH(6,NBH)
            UBFC(NCX,NSL) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFXZ_BH_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBG_FRC_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Loads the matrix elements and solution vector for the
!     gas convective-dispersive mass transport equation for fractures.
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
!     Written by MD White, PNNL, 6 July 2018.
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
      SUB_LOG(ISUB_LOG) = '/SJCBG_FRC_GT'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fracture triangles  ---
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
            IF( IAQU.EQ.0 ) ALU(MDT,MT1X) = ALU(MDT,MT1X) + SC
!
!---      SPLib or Lis solvers  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MCD = KLUC_FRC(IXP_FRC(NT1X),MA)
            MA = MA + 1
            IF( IAQU.EQ.0 ) DLU(MCD) = DLU(MCD) + SC
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            IROW = petsc_offset(MT1X)
            IF( IAQU.EQ.0 ) CALL MatSetValues(T_MAT,1,IROW,1,IROW,SC,
     &        ADD_VALUES,IERR)
#endif  
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
!---      Loop over fracture triangle connections  ---
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
!---        Interfacial fracture aperture  ---
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
!
!---      Skip for active aqueous  ---
!
          IF( IAQU.EQ.0 ) THEN
!!
!!---        Solution vector and parent radioactive decay  ---
!!
            BLU(MT1X) = BLU(MT1X) + CO_FRC(NT1X,NSL)*SC
!     &        - 6.931D-1*CO_FRC(NT1X,NSL)*VOLX/HLF(NSL)
!!
!!---        Daughter radioactive chain decay  ---
!!
!            IF( NSL.LE.NSOLU ) THEN
!              DO NPSL = 1,NSOLU
!                IF( NPSL.EQ.NSL ) CYCLE
!                BLU(MT1X) = BLU(MT1X) + CHDF(NPSL,NSL)*
!     &            6.931D-1*CO_FRC(NT1X,NPSL)*VOLX/HLF(NPSL)
!              ENDDO
!            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBG_FRC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBG_BH_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Loads the matrix elements and solution vector for the
!     gas convective-dispersive mass transport equation for boreholes.
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
!     Written by MD White, PNNL, 6 May 2019.
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
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/SJCBG_BH_GT'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN1 = ID_BH(3,NBH),ID_BH(4,NBH)
          MBN1 = NFLD - NXP + NT_FRC + NBN1
!
!---      Storage terms  ---
!
          SC = VOL_BH(NBN1)*DTI
          MA = 1
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
            IF( IAQU.EQ.0 ) ALU(MDT,MBN1) = ALU(MDT,MBN1) + SC
!
!---      SPLib or Lis solvers  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MCD = KLUC_BH(NBN1,MA)
            MA = MA + 1
            IF( IAQU.EQ.0 ) DLU(MCD) = DLU(MCD) + SC
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            IROW = petsc_offset(MBN1)
            IF( IAQU.EQ.0 ) CALL MatSetValues(T_MAT,1,IROW,1,IROW,SC,
     &        ADD_VALUES,IERR)
#endif  
          ENDIF
!
!---      Molecular diffusion coefficients at borehole node  ---
!
          TCOR = (T_BH(2,NBN1)+TABS)/TSPRF
          PCOR = (PG_BH(2,NBN1)+PATM)/PATM
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SG_BH(2,NBN1).LT.EPSL ) THEN
            TORGX = 0.D+0
          ELSE
            TORGX = (SG_BH(2,NBN1)**7)**(THIRD)
          ENDIF
          SDFG1 = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMC1 = SG_BH(2,NBN1)*PORD_BH(2,NBN1)
          FCG1 = YG_BH(NBN1,NSL)/(VMC1+SMALL)
          DG1 = TORGX*SG_BH(2,NBN1)*SDFG1
!
!---      Loop over borehole node to node connections ---
!
          DO NCX = IPB_BH(1,NBN1),IPB_BH(2,NBN1)
            NBN2 = IBCM_BH(NCX)
            MBN2 = NFLD - NXP + NT_FRC + NBN2
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---        Neglect hydrodynamic dispersion  ---
!
            DPGX = 0.D+0
!
!---        Molecular diffusion coefficients at adjacent triangle  ---
!
            TCOR = (T_BH(2,NBN2)+TABS)/TSPRF
            PCOR = (PG_BH(2,NBN2)+PATM)/PATM
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( SG_BH(2,NBN2).LT.EPSL ) THEN
              TORGX = 0.D+0
            ELSE
              TORGX = (SG_BH(2,NBN2)**7)**(THIRD)
            ENDIF
            SDFG2 = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            VMC2 = SG_BH(2,NBN2)*PORD_BH(2,NBN2)
            FCG2 = YG_BH(NBN2,NSL)/(VMC2+SMALL)
            DG2 = TORGX*SG_BH(2,NBN2)*SDFG2
!
!---        Interfacial fracture aperture  ---
!
            INDX = -1
            INV1 = INV_BH(NBN1)
            IF( IS_BH(INV1).GT.10 ) THEN
              RBX = MAX( PAR_BH(3,INV1),1.D-20 )
            ELSE
              RBX = MAX( PAR_BH(2,INV1),1.D-20 )
            ENDIF    
            AP1X = GPI*(RBX**2)
            INV2 = INV_BH(NBN2)
            IF( IS_BH(INV2).GT.10 ) THEN
              RBX = MAX( PAR_BH(3,INV2),1.D-20 )
            ELSE
              RBX = MAX( PAR_BH(2,INV2),1.D-20 )
            ENDIF    
            AP2X = GPI*(RBX**2)
            APX = DIFMN( AP1X,AP2X,DBB1X,DBB1X,ZERO,INDX )
!
!---        Interfacial diffusion/dispersion coefficient  ---
!
            INDX = 16
            DGX = DIFMN(DG1,DG2,DBB1X,DBB2X,UBBG(1,NCX),INDX)
            DGX = APX*(DGX+DPGX)/DBB_BH(NCX)
            FGX = APX*UBBG(1,NCX)
            IF( MOD(ISLC(23),10).EQ.1 ) FGX = 0.D+0
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
              MCOL = MBN2
              MROW = MBN1-MBN2+MDT
              ALU(MDT,MBN1) = ALU(MDT,MBN1) + A1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) - A2
!
!---        SPLib solver  ---
!
            ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
              DLU(MCD) = DLU(MCD) + A1
              MROW = KLUC_BH(NBN1,MA)
              MA = MA + 1
              DLU(MROW) = DLU(MROW) - A2
#ifdef petsc
!
!---        PETSc solver  ---
!
            ELSEIF( ILES.EQ.5 ) THEN
              CALL MatSetValues(T_MAT,1,IROW,1,IROW,A1,
     &                        ADD_VALUES,IERR)
              ICOL = petsc_offset(MBN2)
              A2 = -A2
              CALL MatSetValues(T_MAT,1,IROW,1,ICOL,A2,
     &                        ADD_VALUES,IERR)
#endif
            ENDIF
          ENDDO
!
!---      Skip for active aqueous  ---
!
          IF( IAQU.EQ.0 ) THEN
!
!---        Solution vector and parent radioactive decay  ---
!
            BLU(MBN1) = BLU(MBN1) + CO_BH(NBN1,NSL)*SC
!     &        - 6.931D-1*CO_BH(NBN1,NSL)*VOL_BH(NBN1)/HLF(NSL)
!!
!!---        Daughter radioactive chain decay  ---
!!
!            IF( NSL.LE.NSOLU ) THEN
!              DO NPSL = 1,NSOLU
!                IF( NPSL.EQ.NSL ) CYCLE
!                BLU(MBN1) = BLU(MBN1) + CHDF(NPSL,NSL)*
!     &            6.931D-1*CO_BH(NBN1,NPSL)*VOL_BH(NBN1)/HLF(NPSL)
!              ENDDO
!            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBG_BH_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBL_FRC_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Loads the matrix elements and solution vector for the aqueous 
!     convective-dispersive mass transport equation for fractures.
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
!     Written by MD White, PNNL, 4 July 2018.
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
      SUB_LOG(ISUB_LOG) = '/SJCBL_FRC_GT'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fracture triangles  ---
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
!---      Loop over fracture triangle connections  ---
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
!---        Interfacial fracture aperture  ---
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
!---      Solution vector and parent radioactive decay  ---
!
          BLU(MT1X) = BLU(MT1X) + CO_FRC(NT1X,NSL)*SC
!     &      - 6.931D-1*CO_FRC(NT1X,NSL)*VOLX/HLF(NSL)
!!
!!---      Daughter radioactive chain decay  ---
!!
!          IF( NSL.LE.NSOLU ) THEN
!            DO NPSL = 1,NSOLU
!              IF( NPSL.EQ.NSL ) CYCLE
!              BLU(MT1X) = BLU(MT1X) +
!     &          CHDF(NPSL,NSL)*6.931D-1*CO_FRC(NT1X,NPSL)*VOLX/HLF(NPSL)
!            ENDDO
!          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBL_FRC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBL_BH_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Loads the matrix elements and solution vector for the aqueous 
!     convective-dispersive mass transport equation for boreholes.
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
!     Written by MD White, PNNL, 6 May 2019.
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
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/SJCBL_BH_GT'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN1 = ID_BH(3,NBH),ID_BH(4,NBH)
          MBN1 = NFLD - NXP + NT_FRC + NBN1
          IZN1 = IZ_BH(NBN1)
          INV1 = INV_BH(NBN1)
!
!---      Storage terms  ---
!
          SC = VOL_BH(NBN1)*DTI
          MA = 1
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
            ALU(MDT,MBN1) = ALU(MDT,MBN1) + SC
!
!---      SPLib or Lis solvers  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MCD = KLUC_BH(NBN1,MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + SC
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            IROW = petsc_offset(MBN1)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,SC,
     &                      ADD_VALUES,IERR)
#endif  
          ENDIF
!
!---      Molecular diffusion coefficients at triangle  ---
!
          TCOR = (T_BH(2,NBN1)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SL_BH(2,NBN1).LT.EPSL ) THEN
            TORLX = 0.D+0
          ELSE
            TORLX = (SL_BH(2,NBN1)**7)**(THIRD)
          ENDIF
          SDFL1 = SMDL(NSL)*TCOR*(VISRL/VISL_BH(2,NBN1))
          VMC1 = SL_BH(2,NBN1)*PORD_BH(2,NBN1)
          FCL1 = YL_BH(NBN1,NSL)/(VMC1+SMALL)
!
!---      Pipe flow  ---
!
          IF( IS_BH(INV1).GT.10 ) THEN
            DL1 = TORLX*VMC1*SDFL1
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DL1 = SDCL(1,IZN1,NSL)*SDCL(2,IZN1,NSL)*
     &        EXP(VMC1*SDCL(3,IZN1,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DL1 = TORLX*VMC1*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DL1 = SDCL(1,IZN1,NSL)*SDCL(2,IZN1,NSL)*
     &        VMC1**SDCL(3,IZN1,NSL)
          ELSE
            DL1 = TORLX*VMC1*SDFL1
          ENDIF
!
!---      Loop over borehole node to node connections ---
!
          DO NCX = IPB_BH(1,NBN1),IPB_BH(2,NBN1)
            NBN2 = IBCM_BH(NCX)
            MBN2 = NFLD - NXP + NT_FRC + NBN2
            IZN2 = IZ_BH(NBN2)
            INV2 = INV_BH(NBN2)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---        Neglect hydrodynamic dispersion  ---
!
            DPLX = 0.D+0
!
!---        Molecular diffusion coefficients at adjacent triangle  ---
!
            TCOR = (T_BH(2,NBN2)+TABS)/TSPRF
!
!---        Millington and Quirk tortuosity model  ---
!
            IF( SL_BH(2,NBN2).LT.EPSL ) THEN
              TORLX = 0.D+0
            ELSE
              TORLX = (SL_BH(2,NBN2)**7)**(THIRD)
            ENDIF
            SDFL2 = SMDL(NSL)*TCOR*(VISRL/VISL_BH(2,NBN2))
            VMC2 = SL_BH(2,NBN2)*PORD_BH(2,NBN2)
            FCL2 = YL_BH(NBN2,NSL)/(VMC2+SMALL)
!
!---        Pipe flow  ---
!
            IF( IS_BH(INV2).GT.10 ) THEN
              DL2 = TORLX*VMC2*SDFL2
            ELSEIF( IEDL(NSL).EQ.2 ) THEN
              DL2 = SDCL(1,IZN2,NSL)*SDCL(2,IZN2,NSL)*
     &          EXP(VMC2*SDCL(3,IZN2,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DL2 = TORLX*VMC2*SMDL(NSL)
            ELSEIF( IEDL(NSL).EQ.4 ) THEN
              DL2 = SDCL(1,IZN2,NSL)*SDCL(2,IZN2,NSL)*
     &          VMC2**SDCL(3,IZN2,NSL)
            ELSE
              DL2 = TORLX*VMC2*SDFL2
            ENDIF
!
!---        Interfacial fracture aperture  ---
!
            INDX = -1
            INV1 = INV_BH(NBN1)
            IF( IS_BH(INV1).GT.10 ) THEN
              RBX = MAX( PAR_BH(3,INV1),1.D-20 )
            ELSE
              RBX = MAX( PAR_BH(2,INV1),1.D-20 )
            ENDIF    
            AP1X = GPI*(RBX**2)
            INV2 = INV_BH(NBN2)
            IF( IS_BH(INV2).GT.10 ) THEN
              RBX = MAX( PAR_BH(3,INV2),1.D-20 )
            ELSE
              RBX = MAX( PAR_BH(2,INV2),1.D-20 )
            ENDIF    
            AP2X = GPI*(RBX**2)
            APX = DIFMN( AP1X,AP2X,DBB1X,DBB1X,ZERO,INDX )
!
!---        Interfacial diffusion/dispersion coefficient  ---
!
            INDX = 16
            DLX = DIFMN(DL1,DL2,DBB1X,DBB2X,UBBL(1,NCX),INDX)
            DLX = APX*(DLX+DPLX)/DBB_BH(NCX)
            FLX = APX*UBBL(1,NCX)
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
              MCOL = MBN2
              MROW = MBN1-MBN2+MDT
              ALU(MDT,MBN1) = ALU(MDT,MBN1) + A1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) - A2
!
!---        SPLib solver  ---
!
            ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
              DLU(MCD) = DLU(MCD) + A1
              MROW = KLUC_BH(NBN1,MA)
              MA = MA + 1
              DLU(MROW) = DLU(MROW) - A2
#ifdef petsc
!
!---        PETSc solver  ---
!
            ELSEIF( ILES.EQ.5 ) THEN
              CALL MatSetValues(T_MAT,1,IROW,1,IROW,A1,
     &                        ADD_VALUES,IERR)
              ICOL = petsc_offset(MBN2)
              A2 = -A2
              CALL MatSetValues(T_MAT,1,IROW,1,ICOL,A2,
     &                        ADD_VALUES,IERR)
#endif
            ENDIF
          ENDDO
!
!---      Solution vector and parent radioactive decay  ---
!
          BLU(MBN1) = BLU(MBN1) + CO_BH(NBN1,NSL)*SC
!     &      - 6.931D-1*CO_BH(NBN1,NSL)*VOL_BH(NBN1)/HLF(NSL)
!!
!!---      Daughter radioactive chain decay  ---
!!
!          IF( NSL.LE.NSOLU ) THEN
!            DO NPSL = 1,NSOLU
!              IF( NPSL.EQ.NSL ) CYCLE
!              BLU(MBN1) = BLU(MBN1) + CHDF(NPSL,NSL)*
!     &          6.931D-1*CO_BH(NBN1,NPSL)*VOL_BH(NBN1)/HLF(NPSL)
!            ENDDO
!          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBL_BH_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_BH_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Borehole solute sources.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 November 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE SOLTN
      USE PORMED
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FDVP_BH
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
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR_BHX(11)
      REAL*8 VARC_BHX(2)
      REAL*8 KGM,KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORT_BH_GT'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        NCT = ITM_BH(NBH)
        IF( ICC_BH(NBH).EQ.1 ) TMZ = MOD( TM,VAR_BH(1,NCT,NBH) )
!
!---    Borehole is inactive  ---
!
        IF( TMZ.LE.VAR_BH(1,1,NBH) ) CYCLE
!
!---    Borehole is active w/ one time stamp  ---
!
        IF( NCT.EQ.1 ) THEN
          DO N = 2,11
            VAR_BHX(N) = VAR_BH(N,1,NBH)
          ENDDO
          NSLX = (NSL-1)*2 + 1
          VARC_BHX(1) = VARC_BH(NSLX,1,NBH)
          NSLX = (NSL-1)*2 + 2
          VARC_BHX(2) = VARC_BH(NSLX,1,NBH)
!
!---    Borehole is active w/ multiple time stamps  ---
!
        ELSE
          DO M = 2,NCT
            IF( TMZ.LE.VAR_BH(1,M,NBH) ) THEN
              TD_BH = VAR_BH(1,M,NBH)-VAR_BH(1,M-1,NBH)
              DT_BH = MIN( VAR_BH(1,M,NBH)-TMZ,DT )
              TF_BH = (TMZ-VAR_BH(1,M-1,NBH))/TD_BH
              DO N = 2,11
                VAR_BHX(N) = VAR_BH(N,M-1,NBH) + 
     &            TF_BH*(VAR_BH(N,M,NBH)-VAR_BH(N,M-1,NBH))
              ENDDO
              NSLX = (NSL-1)*2 + 1
              VARC_BHX(1) = VARC_BH(NSLX,M-1,NBH) + 
     &            TF_BH*(VARC_BH(NSLX,M,NBH)-VARC_BH(NSLX,M-1,NBH))
              NSLX = (NSL-1)*2 + 2
              VARC_BHX(2) = VARC_BH(NSLX,M-1,NBH) + 
     &            TF_BH*(VARC_BH(NSLX,M,NBH)-VARC_BH(NSLX,M-1,NBH))
              EXIT
            ENDIF
          ENDDO
        ENDIF
!    
!---    Loop over start and end boreholes  ---
!    
        DO NBSE = 1,2
!    
!---      If only starting borehole is assigned a type, 
!         skip ending borehole assignment  ---
!    
          IF ( IT_BH(2,NBH).EQ.0 ) CYCLE
!    
!---      Borehole source is applied to first ID_BH(3,NBH) or last 
!         ID_BH(4,NBH) borehole node  ---
!    
          NBN = ID_BH(NBSE+2,NBH)
          IZN = IZ_BH(NBN)
          INVX = INV_BH(NBN)
          MBN = NFLD - NXP + NT_FRC + NBN
          IF( ILES.EQ.1 ) THEN
            MCOL = MBN
            MROW = MDT
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC_BH(NBN,MA)
#ifdef   petsc
          ELSEIF( ILES.EQ.5 ) THEN
            IROW = petsc_offset(MBN)
#endif  
          ENDIF
!    
!---      Define adder for VAR_BH index based on starting (0) 
!         or ending (5) borehole node ---
!               
          NBSE_ID = 0
          IF( NBSE.EQ.2 ) NBSE_ID = 5
!
!---      Initialize source  ---
!
          SORTX = 0.D+0
          PGX = PG_BH(2,NBN) + PATM
          PLX = PL_BH(2,NBN) + PATM
!
!---      Gas mass rate w/ zero water vapor   ---
!
          IF( IT_BH(NBSE,NBH).EQ.1 ) THEN
!
!---        Gas mass source ---
!
            IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
              PX = MAX( PLX,PGX )
              TX = T_BH(2,NBN)
              IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
              XGWX = 0.D+0
              XGAX = MAX( 1.D+0-XGWX,0.D+0 )
              XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
              INDX = 0
              CALL REGION_4( TX,PSWX,INDX )
              PVWX = MIN( PX*XMGWX,PSWX )
              PVAX = MAX( PX-PVWX,0.D+0 )
              INDX = 0
              CALL WATGSD( TX,PVWX,RHOW,INDX )
              CALL AIRGSD( TX,PVAX,RHOA )
              XGWX = RHOW/(RHOW+RHOA)
              XGAX = RHOA/(RHOW+RHOA)
              RHOGX = XGWX*RHOW + XGAX*RHOA
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)/
     &          RHOGX
!
!---        Gas mass sink  ---
!
            ELSE
              IF( SG_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          VAR_BHX(2+NBSE_ID)*YG_BH(NBN,NSL)/
     &          (RHOG_BH(2,NBN)*SG_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
!
!---      Gas mass rate w/ water vapor mass fraction   ---
!
          ELSEIF( IT_BH(NBSE,NBH).EQ.2 ) THEN
!
!---        Gas mass source ---
!
            IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
              PX = MAX( PLX,PGX )
              TX = T_BH(2,NBN)
              IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
              XGWX = VAR_BHX(4+NBSE_ID)
              XGAX = MAX( 1.D+0-XGWX,0.D+0 )
              XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
              INDX = 0
              CALL REGION_4( TX,PSWX,INDX )
              PVWX = MIN( PX*XMGWX,PSWX )
              PVAX = MAX( PX-PVWX,0.D+0 )
              INDX = 0
              CALL WATGSD( TX,PVWX,RHOW,INDX )
              CALL AIRGSD( TX,PVAX,RHOA )
              CALL WATGSH( TX,PSWX,RHOW,HGWX,UGWX )
              CALL AIRGSH( TX,PVAX,HGAX,UGAX )
              XGWX = RHOW/(RHOW+RHOA)
              XGAX = RHOA/(RHOW+RHOA)
              RHOGX = XGWX*RHOW + XGAX*RHOA
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)/
     &          RHOGX
!
!---        Gas mass sink  ---
!
            ELSE
              IF( SG_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          VAR_BHX(2+NBSE_ID)*YG_BH(NBN,NSL)/
     &          (RHOG_BH(2,NBN)*SG_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
!
!---      Gas mass rate w/ water vapor relative humidity   ---
!
          ELSEIF( IT_BH(NBSE,NBH).EQ.3 ) THEN
!
!---        Gas mass source ---
!
            IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
              PX = MAX( PGX,PLX )
              TX = T_BH(2,NBN)
              IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
              INDX = 0
              CALL REGION_4( TX,PSWX,INDX )
              PVWX = PSWX*VAR_BHX(4+NBSE_ID)
              PVAX = MAX( PGX-PVWX,0.D+0 )
              INDX = 0
              CALL WATGSD( TX,PVWX,RHOW,INDX )
              CALL AIRGSD( TX,PVAX,RHOA )
              CALL WATGSH( TX,PSWX,RHOW,HGWX,UGWX )
              CALL AIRGSH( TX,PVAX,HGAX,UGAX )
              XGWX = RHOW/(RHOW+RHOA)
              XGAX = RHOA/(RHOW+RHOA)
              RHOGX = XGWX*RHOW + XGAX*RHOA
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)/
     &          RHOGX
!
!---        Gas mass sink  ---
!
            ELSE
              IF( SG_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          VAR_BHX(2+NBSE_ID)*YG_BH(NBN,NSL)/
     &          (RHOG_BH(2,NBN)*SG_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
!
!---      Aqueous mass rate w/ zero air and salt mass ---
!
          ELSEIF( IT_BH(NBSE,NBH).EQ.4 ) THEN
!
!---        Aqueous Mass Source ---
!
            IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
              TX = T_BH(2,NBN)
              IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
              XLAX = 0.D+0
              XLSX = 0.D+0
              INDX = 0
              CALL REGION_4( TX,PSWX,INDX )
              PX = MAX( PGX,PLX,PSWX )
              CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
              CALL SOL_BRNS( TX,PX,XLSMX )
              XLSX = MIN( XLSMX,XLSX )
              XLWX = MAX( 1.D+0-XLAX-XLSX,0.D+0 )
              XMLAX = (XLAX/WTMA)/((XLAX/WTMA) + (XLSX/WTMS) + 
     &          (XLWX/WTMW))
              CALL SP_B( TX,XLSX,PSBX )
              PVAX = HCAW*XMLAX
              XMGAX = PVAX/(PVAX+PSBX)
              CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &          XMLAX,XMLSX,XMLWX )
              CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)/
     &          RHOLX
!
!---        Aqueous Mass Sink ---
!
            ELSE
              IF( SL_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          VAR_BHX(2+NBSE_ID)*YL_BH(NBN,NSL)/
     &          (RHOL_BH(2,NBN)*SL_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
!
!---      Aqueous mass rate w/ air and salt mass fractions  ---
!
          ELSEIF( IT_BH(NBSE,NBH).EQ.5 ) THEN
!
!---        Aqueous Mass Source ---
!
            IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
              TX = T_BH(2,NBN)
              IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
              XLAX = VAR_BHX(4+NBSE_ID)
              XLSX = VAR_BHX(6+NBSE_ID)
              INDX = 0
              CALL REGION_4( TX,PSWX,INDX )
              PX = MAX( PGX,PLX,PSWX )
              CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
              CALL SOL_BRNS( TX,PX,XLSMX )
              XLSX = MIN( XLSMX,XLSX )
              XLWX = MAX( 1.D+0-XLAX-XLSX,0.D+0 )
              XMLAX = (XLAX/WTMA)/((XLAX/WTMA) + (XLSX/WTMS) + 
     &          (XLWX/WTMW))
              CALL SP_B( TX,XLSX,PSBX )
              PVAX = HCAW*XMLAX
              XMGAX = PVAX/(PVAX+PSBX)
              CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &          XMLAX,XMLSX,XMLWX )
              CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)/
     &          RHOLX
!
!---        Aqueous Mass Sink ---
!
            ELSE
              IF( SL_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          VAR_BHX(2+NBSE_ID)*YL_BH(NBN,NSL)/
     &          (RHOL_BH(2,NBN)*SL_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
!
!---      Aqueous mass rate w/ air and salt relative saturations  ---
!
          ELSEIF( IT_BH(NBSE,NBH).EQ.6 ) THEN
!
!---        Aqueous Mass Source ---
!
            IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
              TX = T_BH(2,NBN)
              IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
              WLAX = VAR_BHX(4+NBSE_ID)
              WLSX = VAR_BHX(6+NBSE_ID)
              INDX = 0
              CALL REGION_4( TX,PSWX,INDX )
              PX = MAX( PGX,PLX,PSWX )
              CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
              CALL SOL_BRNS( TX,PX,XLSMX )
              XLSX = WLSX*XLSMX
              CALL SP_B( TX,XLSX,PSBX )
              PVAX = WLAX*MAX( PX-PSBX,0.D+0 )
              XMLAX = PVAX/HCAW
              XMGAX = PVAX/(PVAX+PSBX)
              CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &          XMLAX,XMLSX,XMLWX )
              CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)/
     &          RHOLX
!
!---        Aqueous Mass Sink ---
!
            ELSE
              IF( SL_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          VAR_BHX(2+NBSE_ID)*YL_BH(NBN,NSL)/
     &          (RHOL_BH(2,NBN)*SL_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
!
!---      Gas volumetric rate ---
!
          ELSEIF( IT_BH(NBSE,NBH).GE.11.AND.IT_BH(NBSE,NBH).LE.13 ) THEN
!
!---        Gas volumetric source ---
!
            IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)
!
!---        Gas volumetric sink  ---
!
            ELSE
              IF( SG_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          VAR_BHX(2+NBSE_ID)*YG_BH(NBN,NSL)/
     &          SG_BH(2,NBN)*PORD_BH(2,NBN)
            ENDIF
!
!---      Aqueous volumetric rate ---
!
          ELSEIF( IT_BH(NBSE,NBH).GE.14.AND.IT_BH(NBSE,NBH).LE.16 ) THEN
!
!---        Aqueous volumetric source ---
!
            IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)
!
!---        Aqueous volumetric sink  ---
!
            ELSE
              IF( SL_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          VAR_BHX(3+NBSE_ID)*YL_BH(NBN,NSL)/
     &          (SL_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
!
!---      Pressure borehole, state condition #1 (saturated)  ---
!
          ELSEIF( IT_BH(NBSE,NBH).EQ.7 ) THEN
!
!---        Distance between borehole node centroid and
!           pressure boundary  ---
!
            DBBX = 5.D-1*SQRT( (XP_BH(1,NBN)-XP_BH(2,NBN))**2 + 
     &        (YP_BH(1,NBN)-YP_BH(2,NBN))**2 + 
     &        (ZP_BH(1,NBN)-ZP_BH(2,NBN))**2 )
!
!---        Borehole pressure  ---
!
            PL_BHX = VAR_BHX(3+NBSE_ID)
            PG_BHX = PL_BHX
            ZP_BHX = ZP_BH(1,NBN)
            RKL_BHX = 1.D+0
            IF( IS_BH(INVX).GT.10 ) THEN
              AF_BHX = GPI*(PAR_BH(3,INVX)**2)
            ELSE
              AF_BHX = GPI*(PAR_BH(2,INVX)**2)
            ENDIF    
!
!---        Aqueous flux (m^3/s) into first borehole node of 
!           borehole  ---
!
            PLX = PL_BH(2,NBN) + PATM
            PGX = PG_BH(2,NBN) + PATM
            ZPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
            HDLX = PL_BHX - PLX + 
     &        GRAV*(ZP_BHX-ZPX)*RHOL_BH(2,NBN)
            IF( M.EQ.1 ) HDL = HDLX
            INDX = 8
            RKLM = DIFMN(RKL_BHX,RKL_BH(2,NBN),DBBX,DBBX,HDL,INDX)
            VLM = VISL_BH(2,NBN)
!
!---        Pipe flow  ---
!
            IF( IS_BH(INVX).GT.10 ) THEN
!
!---          Hydraulic diameter (m), twice interval radius ---
!
              DHX = 2.D+0*PAR_BH(3,INVX)
!
!---          Roughness, m  ---
!
              RFX = PAR_BH(6,INVX)
!
!---          Pipe flow velocity (m/s), first estimating from  
!             Hazen-Williamsequation, adjusting for fluid viscosity, 
!             then solving a nonlinear system of equations, combining  
!             the Colebrook equation for the friction coefficient and 
!             the Darcy-Weisbach equation for the 
!             pressure drop (i.e., head loss)  ---
!
              CALL PIPE_FLOW_VEL_GT( DBBX,DHX,HDLX,RFX,RHOL_BH(2,NBN),
     &          ULX,VLM,NBN )           
              FLX = AF_BHX*SIGN(ULX,HDLX)*RKLM
            ELSE
              KLM = PERM(1,IZN)
              IF( PERM(1,IZN)/EPSL.LT.EPSL ) KLM = 0.D+0
              FLX = AF_BHX*KLM*RKLM*HDLX/DBBX/VLM
            ENDIF
!
!---        Aqueous volumetric source ---
!
            IF( FLX.GE.0.D+0 ) THEN
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)
!
!---        Aqueous volumetric sink  ---
!
            ELSE
              IF( SL_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          FLX*YL_BH(NBN,NSL)/(SL_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
!
!---      Pressure borehole, state condition #2 
!         (partially saturated)  ---
!
          ELSEIF( IT_BH(NBSE,NBH).EQ.8 ) THEN
!
!---        Distance between borehole node centroid and
!           pressure boundary  ---
!
            DBBX = 5.D-1*SQRT( (XP_BH(1,NBN)-XP_BH(2,NBN))**2 + 
     &        (YP_BH(1,NBN)-YP_BH(2,NBN))**2 + 
     &        (ZP_BH(1,NBN)-ZP_BH(2,NBN))**2 )
!
!---        Borehole pressure  ---
!
            PG_BHX = VAR_BHX(3+NBSE_ID)
            SL_BHX = VAR_BHX(2+NBSE_ID)
            SG_BHX = 1.D+0 - SL_BHX
            CALL CAP_BH_GT( CPGLX,SL_BHX,NBN )
            PL_BHX = PG_BHX - CPGLX
            CALL RKSP_BH_GT( PG_BHX,PL_BHX,RKG_BHX,RKL_BHX,SGX,SLX,
     &        NBN )
            ZP_BHX = ZP_BH(1,NBN)
            IF( IS_BH(INVX).GT.10 ) THEN
              AF_BHX = GPI*(PAR_BH(3,INVX)**2)
            ELSE
              AF_BHX = GPI*(PAR_BH(2,INVX)**2)
            ENDIF    
!
!---        Aqueous flux (m^3/s) into first borehole node of 
!           borehole  ---
!
            PLX = PL_BH(2,NBN) + PATM
            ZPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
            HDLX = PL_BHX - PLX + 
     &        GRAV*(ZP_BHX-ZPX)*RHOL_BH(2,NBN)
            IF( M.EQ.1 ) HDL = HDLX
            INDX = 8
            RKLM = DIFMN(RKL_BHX,RKL_BH(2,NBN),DBBX,DBBX,HDL,INDX)
            VLM = VISL_BH(2,NBN)
!
!---        Pipe flow  ---
!
            IF( IS_BH(INVX).GT.10 ) THEN
!
!---          Hydraulic diameter (m), twice interval radius ---
!
              DHX = 2.D+0*PAR_BH(3,INVX)
!
!---          Roughness, m  ---
!
              RFX = PAR_BH(6,INVX)
!
!---          Pipe flow velocity (m/s), first estimating from  
!             Hazen-Williamsequation, adjusting for fluid viscosity, 
!             then solving a nonlinear system of equations, combining  
!             the Colebrook equation for the friction coefficient and 
!             the Darcy-Weisbach equation for the 
!             pressure drop (i.e., head loss)  ---
!
              CALL PIPE_FLOW_VEL_GT( DBBX,DHX,HDLX,RFX,RHOL_BH(2,NBN),
     &          ULX,VLM,NBN )           
              FLX = AF_BHX*SIGN(ULX,HDLX)*RKLM
            ELSE
              KLM = PERM(1,IZN)
              IF( PERM(1,IZN)/EPSL.LT.EPSL ) KLM = 0.D+0
              FLX = AF_BHX*KLM*RKLM*HDLX/DBBX/VLM
            ENDIF
!
!---        Aqueous volumetric source ---
!
            IF( FLX.GE.0.D+0 ) THEN
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)*
     &          SL_BHX
!
!---        Aqueous volumetric sink  ---
!
            ELSE
              IF( SL_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          FLX*YL_BH(NBN,NSL)/(SL_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
!
!---        Gas flux (m^3/s) into first borehole node of 
!           borehole  ---
!
            PGX = PG_BH(2,NBN) + PATM
            ZPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
            HDGX = PG_BHX - PGX + 
     &        GRAV*(ZP_BHX-ZPX)*RHOG_BH(2,NBN)
            IF( M.EQ.1 ) HDG = HDGX
            INDX = 8
            RKGM = DIFMN(RKG_BHX,RKG_BH(2,NBN),DBBX,DBBX,HDG,INDX)
            VGM = VISG_BH(2,NBN)
!
!---        Pipe flow  ---
!
            IF( IS_BH(INVX).GT.10 ) THEN
!
!---          Hydraulic diameter (m), twice interval radius ---
!
              DHX = 2.D+0*PAR_BH(3,INVX)
!
!---          Roughness, m  ---
!
              RFX = PAR_BH(6,INVX)
!
!---          Pipe flow velocity (m/s), first estimating from  
!             Hazen-Williamsequation, adjusting for fluid viscosity, 
!             then solving a nonlinear system of equations, combining  
!             the Colebrook equation for the friction coefficient and 
!             the Darcy-Weisbach equation for the 
!             pressure drop (i.e., head loss)  ---
!
              CALL PIPE_FLOW_VEL_GT( DBBX,DHX,HDGX,RFX,RHOG_BH(2,NBN),
     &          UGX,VGM,NBN )           
              FGX = AF_BHX*SIGN(UGX,HDGX)*RKGM
            ELSE
              KGM = PERM(1,IZN)
              IF( PERM(1,IZN)/EPSL.LT.EPSL ) KLM = 0.D+0
              FGX = AF_BHX*KGM*RKGM*HDGX/DBBX/VGM
            ENDIF
!
!---        Gas volumetric source ---
!
            IF( FGX.GE.0.D+0 ) THEN
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)*
     &          SG_BHX
!
!---        Gas volumetric sink  ---
!
            ELSE
              IF( SG_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          FGX*YG_BH(NBN,NSL)/(SG_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
!
!---      Pressure borehole, state condition #3 (unsaturated)  ---
!
          ELSEIF( IT_BH(NBSE,NBH).EQ.9 ) THEN
!
!---        Distance between borehole node centroid and
!           pressure boundary  ---
!
            DBBX = 5.D-1*SQRT( (XP_BH(1,NBN)-XP_BH(2,NBN))**2 + 
     &        (YP_BH(1,NBN)-YP_BH(2,NBN))**2 + 
     &        (ZP_BH(1,NBN)-ZP_BH(2,NBN))**2 )
!
!---        Borehole pressure  ---
!
            PG_BHX = VAR_BHX(3+NBSE_ID)
            SL_BHX = 0.D+0
            SG_BHX = 1.D+0 - SL_BHX
            RKG_BHX = 1.D+0
            ZP_BHX = ZP_BH(1,NBN)
            IF( IS_BH(INVX).GT.10 ) THEN
              AF_BHX = GPI*(PAR_BH(3,INVX)**2)
            ELSE
              AF_BHX = GPI*(PAR_BH(2,INVX)**2)
            ENDIF    
!
!---        Gas flux (m^3/s) into first borehole node of 
!           borehole  ---
!
            PGX = PG_BH(2,NBN) + PATM
            ZPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
            HDGX = PG_BHX - PGX + 
     &        GRAV*(ZP_BHX-ZPX)*RHOG_BH(2,NBN)
            IF( M.EQ.1 ) HDG = HDGX
            INDX = 8
            RKGM = DIFMN(RKG_BHX,RKG_BH(2,NBN),DBBX,DBBX,HDG,INDX)
            VGM = VISG_BH(2,NBN)
!
!---        Pipe flow  ---
!
            IF( IS_BH(INV_BH(NBN)).GT.10 ) THEN
!
!---          Hydraulic diameter (m), sum of interval radii ---
!
              DHX = 2.D+0*PAR_BH(3,INVX)
!
!---          Roughness, m  ---
!
              RFX = PAR_BH(6,INVX)
!
!---          Pipe flow velocity (m/s), first estimating from  
!             Hazen-Williamsequation, adjusting for fluid viscosity, 
!             then solving a nonlinear system of equations, combining  
!             the Colebrook equation for the friction coefficient and 
!             the Darcy-Weisbach equation for the 
!             pressure drop (i.e., head loss)  ---
!
              CALL PIPE_FLOW_VEL_GT( DBBX,DHX,HDGX,RFX,RHOG_BH(2,NBN),
     &          UGX,VGM,NBN )           
              FLX = AF_BHX*SIGN(UGX,HDGX)*RKGM
            ELSE
              KGM = PERM(1,IZN)
              IF( PERM(1,IZN)/EPSL.LT.EPSL ) KGM = 0.D+0
              FGX = AF_BHX*KGM*RKGM*HDGX/DBBX/VGM
            ENDIF
!
!---        Gas volumetric source ---
!
            IF( FGX.GE.0.D+0 ) THEN
              BLU(MBN) = BLU(MBN) + VARC_BHX(NBSE)*VAR_BHX(2+NBSE_ID)
!
!---        Gas volumetric sink  ---
!
            ELSE
              IF( SG_BH(2,NBN).GT.EPSL ) SORTX = SORTX - 
     &          FGX*YG_BH(NBN,NSL)/(SG_BH(2,NBN)*PORD_BH(2,NBN))
            ENDIF
          ENDIF
!
!---      Load Jacobian  ---
!
          IF( ILES.EQ.1 ) THEN
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef   petsc
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,SORTX,ADD_VALUES,IERR)
#endif  
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORT_BH_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_FRC_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Fracture solute sources.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 July 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE FDVP_FRC
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
      REAL*8 SRX(7+LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORT_FRC_GT'
!
!---  Loop over fracture sources  ---
!
      DO NS = 1,NSR_FRC
        IF( TM.LE.SRC_FRC(1,1,NS) ) CYCLE
        SRX(1) = TM
        IF( ISRM_FRC(NS).EQ.1 ) THEN
          DO N = 2,7+NSOLU
            SRX(N) = SRC_FRC(N,1,NS)
          ENDDO
        ELSE
          DO M = 2,ISRM_FRC(NS)
            IF( TM.LE.SRC_FRC(1,M,NS) ) THEN
              DTSR = MIN( SRC_FRC(1,M,NS)-TM,DT )
              TFSR = (TM-0.5D+0*DTSR-SRC_FRC(1,M-1,NS))/
     &          (SRC_FRC(1,M,NS)-SRC_FRC(1,M-1,NS))
              DO N = 2,7+NSOLU
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
!---    Loop over fracture source domain  ---
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
!---      Aqueous Volumetric Rate Sink ---
!
          IF( (ISRT_FRC(NS).EQ.11 .OR. ISRT_FRC(NS).EQ.12)
     &      .AND. SRX(4).LT.0.D+0 ) THEN
            IF( SL_FRC(2,NTX).GT.EPSL ) SORTX = -SRX(4)*YL_FRC(NTX,NSL)/
     &        SL_FRC(2,NTX)
!
!---      Aqueous Mass Rate Sink ---
!
          ELSEIF( (ISRT_FRC(NS).EQ.13 .OR. ISRT_FRC(NS).EQ.14)
     &      .AND. SRX(4).LT.0.D+0 ) THEN
            IF( SL_FRC(2,NTX).GT.EPSL ) SORTX = -SRX(4)*YL_FRC(NTX,NSL)/
     &        (RHOL_FRC(2,NTX)*SL_FRC(2,NTX))
!
!---      Pressure sink  ---
!
          ELSEIF( ISRT_FRC(NS).EQ.35 ) THEN
            PGX = PG_FRC(2,NTX) + PATM
            PLX = PL_FRC(2,NTX) + PATM
!
!---        Borehole pressuure, radius, and skin factor  ---
!
            PBHX = SRX(4) + PATM
            RBHX = SRX(3)
            SBHX = SRX(2)
!
!---        Fracture triangle equivalent radius  ---
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
            PLX = PL_FRC(2,NTX) + PATM
            FLX = MIN( (PBHX-PLX),0.D+0 )*WIX*RKL_FRC(2,NTX)*
     &        RHOL_FRC(2,NTX)/VISL_FRC(2,NTX)
            IF( SL_FRC(2,NTX).GT.EPSL ) SORTX = SORTX -
     &        FLX*YL_FRC(NTX,NSL)/(RHOL_FRC(2,NTX)*SL_FRC(2,NTX))
!
!---        Gas mass flux into borehole (only negative flux 
!           allowed)  ---
!
            PGX = PG_FRC(2,NTX) + PATM
            FGX = MIN( (PBHX-PGX),0.D+0 )*WIX*RKG_FRC(2,NTX)*
     &        RHOG_FRC(2,NTX)/VISG_FRC(2,NTX)
            IF( SG_FRC(2,NTX).GT.EPSL ) SORTX = SORTX - 
     &        FGX*YG_FRC(NTX,NSL)/(RHOG_FRC(2,NTX)*SG_FRC(2,NTX))
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
!---  End of SORT_FRC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP_FRC_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Calculates the aqueous- and gas-phase solute mole fractions
!     from user-specified partition coefficients for fracture flow.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 3 July 2018.
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
      SUB_LOG(ISUB_LOG) = '/SPRP_FRC_GT'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fracture triangles  ---
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
          XVL = SL_FRC(2,NTX)
          XVG = SG_FRC(2,NTX)
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
!---      Phase-volumetric concentration ratios  ---
!
          YVL = 1.D+0/(XVS + XVL + XVG*PCGLX)
          YVG = 1.D+0/((XVS + XVL)/PCGLX + XVG)
!
!---      Phase mole fractions  ---
!
          YL_FRC(NTX,NSL) = XVL*YVL
          YG_FRC(NTX,NSL) = XVG*YVG
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPRP_FRC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP_BH_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Calculates the aqueous- and gas-phase solute mole fractions
!     from user-specified partition coefficients for borehole flow.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE PARM_BH
      USE GEOM_BH
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/SPRP_BH_GT'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21 .AND. IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
          IZN = IZ_BH(NBN)
          INVX = INV_BH(NBN)
!
!---      Pipe flow  ---
!
          IF( IS_BH(INVX).GT.10 ) THEN
            XVS = 0.D+0
          ELSEIF( IPCL(NSL).EQ.2 ) THEN
            XVS = SL_BH(2,NBN)*RHOS(IZN)*(1.D+0-PORD_BH(2,NBN))*
     &        PCSL(1,IZN,NSL)
          ELSE
            XVS = RHOS(IZN)*(1.D+0-PORD_BH(2,NBN))*
     &        PCSL(1,IZN,NSL)
          ENDIF
          XVL = SL_BH(2,NBN)*PORD_BH(2,NBN)
          XVG = SG_BH(2,NBN)*PORD_BH(2,NBN)
!
!---      Constant gas-aqueous partition coefficient  ---
!
          IF( IPCGL(NSL).EQ.0 ) THEN
            PCGLX = PCGL(1,NSL)
!
!---      Temperature dependent gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.1 ) THEN
            TK = T_BH(2,NBN)+TABS
            PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &        + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &        + PCGL(5,NSL)*TK**2 )
!
!---      Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.2 ) THEN
            PCGLX = RHOG_BH(2,NBN)*XGW_BH(2,NBN)/
     &        (RHOL_BH(2,NBN)*XLW_BH(2,NBN))
          ENDIF
          PCGLX = MAX( PCGLX,1.D-20 )
          PCGLX = MIN( PCGLX,1.D+20 )
!
!---      Phase-volumetric concentration ratios  ---
!
          YVL = 1.D+0/(XVS + XVL + XVG*PCGLX)
          YVG = 1.D+0/((XVS + XVL)/PCGLX + XVG)
!
!---      Phase mole fractions  ---
!
          YL_BH(NBN,NSL) = XVL*YVL
          YG_BH(NBN,NSL) = XVG*YVG
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPRP_BH_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRNSC_BF_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Solute/reactive species transport borehole to fracture transfer
!     functions, borehole and fracture equations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 October 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/TRNSC_BF_GT'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21 .AND. IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole node to fracture triangle connections  ---
!
        DO NCX  = ID_BH(5,NBH),ID_BH(6,NBH)
          NTX = IBHT_FRC(NCX)
          NBN = IBHN_FRC(NCX)
          IZN = IZ_BH(NBN)
          INVX = INV_BH(NBN)
          MTX = NFLD - NXP + NTX
          MBN = NFLD - NXP + NT_FRC + NBN 
!
!---      Borehole node centroid  ---
!
          XBHX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
          YBHX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
          ZBHX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
!
!---      Distance between borehole centroid and borehole-fracture
!         intersection  ---
!
          DBF1X = SQRT((XBHX-XP_BF(NCX))**2 + (YBHX-YP_BF(NCX))**2 +
     &      (ZBHX-ZP_BF(NCX))**2)
!
!---      Distance between fracture centroid and borehole-fracture
!         intersection  ---
!
          DBF2X = SQRT((XP_FRC(NTX)-XP_BF(NCX))**2 + 
     &      (YP_FRC(NTX)-YP_BF(NCX))**2 + (ZP_FRC(NTX)-ZP_BF(NCX))**2)
!
!---      Distance between borehole node centroid and fracture
!         triangle centroid  ---
!
          DBFX = MAX( DBF1X+DBF2X,1.D-6 )
!
!---      Borehole side gas diffusion  ---
!
          TCOR = (T_BH(2,NBN)+TABS)/TSPRF
          PCOR = (PG_BH(2,NBN)+PATM)/PATM
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SG_BH(2,NBN).LT.EPSL ) THEN
            TORGX = 0.D+0
          ELSE
            TORGX = (SG_BH(2,NBN)**7)**(THIRD)
          ENDIF
          SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCB = SG_BH(2,NBN)*PORD_BH(2,NBN)
          FCGB = YG_BH(NBN,NSL)/(VMCB+SMALL)
          DGB = TORGX*SG_BH(2,NBN)*PORD_BH(2,NBN)*SDFGB
!
!---      Fracture side gas diffusion  ---
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
!---      Interfacial gas diffusion between borehole and fracture  ---
!
          INDX = 16
          DGX = DIFMN(DGB,DGF,DBF1X,DBF2X,ZERO,INDX)
          APX = GPI*2.D+1*PAR_BH(2,INVX)*APM_FRC(2,NTX)
          DGX = APX*DGX/DBFX
          FGX = APX*UBFG(1,NCX)
!
!---      Patankar solute transport between borehole and fracture  ---
!
          AGBX = MAX(-FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          AGBBX = (AGBX+FGX)*FCGB
          AGBFX = AGBX*FCGF
          AGFX = MAX(FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          AGFFX = (AGFX-FGX)*FCGF
          AGFBX = AGFX*FCGB
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Gas transport between borehole and fracture, borehole
!           equations  ---
!
            MCOL = MBN
            MROW = MTX-MBN+MDT
            ALU(MDT,MBN) = ALU(MDT,MBN) + AGBBX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AGBFX
!
!---        Gas transport between borehole and fracture, fracture
!           equations  ---
!
            MCOL = MTX
            MROW = MBN-MTX+MDT
            ALU(MDT,MTX) = ALU(MDT,MTX) + AGFFX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AGFBX
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Gas transport between borehole and fracture, borehole
!           equations  ---
!
            MA = 1
            MCD = KLUC_BH(NBN,MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + AGBBX
            MROW = KLUC_BCF(NCX,1)
            DLU(MROW) = DLU(MROW) - AGBFX
!
!---        Gas transport between borehole and fracture, fracture
!           equations  ---
!
            MA = 1
            MCD = KLUC_FRC(IXP_FRC(NTX),MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + AGFFX
            MROW = KLUC_FCB(NCX,1)
            DLU(MROW) = DLU(MROW) - AGFBX
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
!
!---        Gas transport between borehole and fracture, borehole
!           equations  ---
!
            IROW = petsc_offset(MBN)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AGBBX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MTX)
            AGBFX = -AGBFX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AGBFX,
     &                      ADD_VALUES,IERR)
!
!---        Gas transport between borehole and fracture, fracture
!           equations  ---
!
            IROW = petsc_offset(MTX)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AGFFX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MBN)
            AGFBX = -AGFBX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AGFBX,
     &                      ADD_VALUES,IERR)
#endif
          ENDIF
!
!---      Borehole side aqueous diffusion  ---
!
          TCOR = (T_BH(2,NBN)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SL_BH(2,NBN).LT.EPSL ) THEN
            TORLX = 0.D+0
          ELSE
            TORLX = (SL_BH(2,NBN)**7)**(THIRD)
          ENDIF
          SDFLB = SMDL(NSL)*TCOR*(VISRL/VISL_BH(2,NBN))
          VMCB = SL_BH(2,NBN)*PORD_BH(2,NBN)
          FCLB = YL_BH(NBN,NSL)/(VMCB+SMALL)
!
!---      Pipe flow  ---
!
          IF( IS_BH(INVX).GT.10 ) THEN
            DLB = TORLX*VMCB*SDFLB
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &        EXP(VMCB*SDCL(3,IZN,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLB = TORLX*VMCB*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &        VMCB**SDCL(3,IZN,NSL)
          ELSE
            DLB = TORLX*VMCB*SDFLB
          ENDIF
!
!---      Fracture side aqueous diffusion  ---
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
!---      Interfacial aqueous diffusion between borehole and 
!         fracture  ---
!
          INDX = 16
          DLX = DIFMN(DLB,DLF,DBF1X,DBF2X,ZERO,INDX)
          APX = GPI*2.D+1*PAR_BH(2,INVX)*APM_FRC(2,NTX)
          DLX = APX*DLX/DBFX
          FLX = APX*UBFL(1,NCX)
!
!---      Patankar solute transport between borehole and fracture  ---
!
          ALBX = MAX(-FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          ALBBX = (ALBX+FLX)*FCLB
          ALBFX = ALBX*FCLF
          ALFX = MAX(FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          ALFFX = (ALFX-FLX)*FCLF
          ALFBX = ALFX*FCLB
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Aqueous transport between borehole and fracture, borehole
!           equations  ---
!
            MCOL = MBN
            MROW = MTX-MBN+MDT
            ALU(MDT,MBN) = ALU(MDT,MBN) + ALBBX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - ALBFX
!
!---        Aqueous transport between borehole and fracture, fracture
!           equations  ---
!
            MCOL = MTX
            MROW = MBN-MTX+MDT
            ALU(MDT,MTX) = ALU(MDT,MTX) + ALFFX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - ALFBX
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Aqueous transport between borehole and fracture, borehole
!           equations  ---
!
            MA = 1
            MCD = KLUC_BH(NBN,MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + ALBBX
            MROW = KLUC_BCF(NCX,1)
            DLU(MROW) = DLU(MROW) - ALBFX
!
!---        Aqueous transport between borehole and fracture, fracture
!           equations  ---
!
            MA = 1
            MCD = KLUC_FRC(IXP_FRC(NTX),MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + ALFFX
            MROW = KLUC_FCB(NCX,1)
            DLU(MROW) = DLU(MROW) - ALFBX
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
!
!---        Aqueous transport between borehole and fracture, borehole
!           equations  ---
!
            IROW = petsc_offset(MBN)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,ALBBX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MTX)
            ALBFX = -ALBFX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,ALBFX,
     &                      ADD_VALUES,IERR)
!
!---        Aqueous transport between borehole and fracture, fracture
!           equations  ---
!
            IROW = petsc_offset(MTX)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,ALFFX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MBN)
            ALFBX = -ALFBX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,ALFBX,
     &                      ADD_VALUES,IERR)
#endif
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TRNSC_BF_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRNSC_BM_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Solute/reactive species transport borehole to matrix transfer
!     functions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 October 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE JACOB
      USE HYST
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_BH
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/TRNSC_BM_GT'
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21 .AND. IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
          INVX = INV_BH(NBN)
!
!---      Cased borehole  ---
!
          IF( MOD(IS_BH(INVX),10).EQ.1 ) CYCLE
!
!---      Field node connected to borehole node  ---
!
          N = IBN_BH(NBN)
          MCP = IXP(N)
          MBN = NFLD - NXP + NT_FRC + NBN 
!
!---      Borehole node centroid and surface area  ---
!
          XP_BHX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
          YP_BHX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
          ZP_BHX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
          RBX = MAX( PAR_BH(2,INVX),1.D-20 )    
          AFB_BHX = 2.D+0*GPI*RBX*SQRT( (XP_BH(2,NBN)-XP_BH(1,NBN))**2 + 
     &      (YP_BH(2,NBN)-YP_BH(1,NBN))**2 + 
     &      (ZP_BH(2,NBN)-ZP_BH(1,NBN))**2 )
!
!---      Borehole side gas diffusion  ---
!
          TCOR = (T_BH(2,NBN)+TABS)/TSPRF
          PCOR = (PG_BH(2,NBN)+PATM)/PATM
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SG_BH(2,NBN).LT.EPSL ) THEN
            TORGX = 0.D+0
          ELSE
            TORGX = (SG_BH(2,NBN)**7)**(THIRD)
          ENDIF
          SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCB = SG_BH(2,NBN)*PORD_BH(2,NBN)
          FCGB = YG_BH(NBN,NSL)/(VMCB+SMALL)
          DGB = TORGX*SG_BH(2,NBN)*PORD_BH(2,NBN)*SDFGB
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
          DGX = DIFMN(DGB,DGP,DBN_BH(NBN),DBN_BH(NBN),ZERO,INDX)
          DGX = AFB_BHX*DGX/DBN_BH(NBN)
          FGX = UBMG(NBN)
!
!---      Patankar solute transport between borehole node and 
!         matrix grid cell  ---
!
          AGBX = MAX(-FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          AGBBX = (AGBX+FGX)*FCGB
          AGBPX = AGBX*FCGP
          AGPX = MAX(FGX,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGX)/(DGX+SMALL)))**5,ZERO)
          AGPPX = (AGBX-FGX)*FCGP
          AGPBX = AGPX*FCGB
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Borehole equations  ---
!
            MCOL = MBN
            MROW = MCP-MBN+MDT
            ALU(MDT,MBN) = ALU(MDT,MBN) + AGBBX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AGBPX
!
!---        Matrix equations  ---
!
            MCOL = MCP
            MROW = MBN-MCP+MDT
            ALU(MDT,MCP) = ALU(MDT,MCP) + AGPPX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AGPBX
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Borehole equations  ---
!
            MA = 1
            MCD = KLUC_BH(NBN,MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + AGBBX
            MROW = KLUC_BCM(NBN,1)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AGBPX
!
!---        Matrix equations  ---
!
            MA = 1
            MCD = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + AGPPX
            MROW = KLUC_MCB(NBN,1)
            DLU(MROW) = DLU(MROW) - AGPBX
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
!
!---        Borehole equations  ---
!
            IROW = petsc_offset(MBN)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AGBBX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MCP)
            AGBPX = -AGBPX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AGBPX,
     &                      ADD_VALUES,IERR)
!
!---        Matrix equations  ---
!
            IROW = petsc_offset(MCP)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AGPPX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MBN)
            AGPBX = -AGPBX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AGPBX,
     &                      ADD_VALUES,IERR)
#endif
          ENDIF
!
!---      Borehole side aqueous diffusion  ---
!
          TCOR = (T_BH(2,NBN)+TABS)/TSPRF
!
!---      Millington and Quirk tortuosity model  ---
!
          IF( SL_BH(2,NBN).LT.EPSL ) THEN
            TORLX = 0.D+0
          ELSE
            TORLX = (SL_BH(2,NBN)**7)**(THIRD)
          ENDIF
          SDFLB = SMDL(NSL)*TCOR*(VISRL/VISL_BH(2,NBN))
          VMCB = SL_BH(2,NBN)*PORD_BH(2,NBN)
          FCLB = YL_BH(NBN,NSL)/(VMCB+SMALL)
          IZN = IZ_BH(NBN)
          INVX = INV_BH(NBN)
!
!---      Pipe flow  ---
!
          IF( IS_BH(INVX).GT.10 ) THEN
            DLB = TORLX*VMCB*SDFLB
          ELSEIF( IEDL(NSL).EQ.2 ) THEN
            DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &        EXP(VMCB*SDCL(3,IZN,NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLB = TORLX*VMCB*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &        VMCB**SDCL(3,IZN,NSL)
          ELSE
            DLB = TORLX*VMCB*SDFLB
          ENDIF
!
!---      Grid-cell side aqueous diffusion  ---
!
          TCOR = (T(2,N)+TABS)/TSPRF
          SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
          VMCP = SL(2,N)*PORD(2,N)
          FCLP = YL(N,NSL)/(VMCP+SMALL)
          IZN = IZ(N)
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
          DLX = DIFMN(DLB,DLP,DBN_BH(NBN),DBN_BH(NBN),ZERO,INDX)
          DLX = AFB_BHX*DLX/DBN_BH(NBN)
          FLX = UBML(NBN)
!
!---      Patankar solute transport between borehole node and 
!         matrix grid cell  ---
!
          ALBX = MAX(-FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          ALBBX = (ALBX+FLX)*FCLB
          ALBPX = ALBX*FCLP
          ALPX = MAX(FLX,ZERO)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLX)/(DLX+SMALL)))**5,ZERO)
          ALPPX = (ALPX-FLX)*FCLP
          ALPBX = ALPX*FCLB
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Borehole equations  ---
!
            MCOL = MBN
            MROW = MCP-MBN+MDT
            ALU(MDT,MBN) = ALU(MDT,MBN) + ALBBX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - ALBPX
!
!---        Matrix equations  ---
!
            MCOL = MCP
            MROW = MBN-MCP+MDT
            ALU(MDT,MCP) = ALU(MDT,MCP) + ALPPX
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - ALPBX
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Borehole equations  ---
!
            MA = 1
            MCD = KLUC_BH(NBN,MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + ALBBX
            MROW = KLUC_BCM(NBN,1)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - ALBPX
!
!---        Matrix equations  ---
!
            MA = 1
            MCD = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MCD) = DLU(MCD) + ALPPX
            MROW = KLUC_MCB(NBN,1)
            DLU(MROW) = DLU(MROW) - ALPBX
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
!
!---        Borehole equations  ---
!
            IROW = petsc_offset(MBN)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,ALBBX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MCP)
            ALBPX = -ALBPX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,ALBPX,
     &                      ADD_VALUES,IERR)
!
!---        Matrix equations  ---
!
            IROW = petsc_offset(MCP)
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,ALPPX,
     &                      ADD_VALUES,IERR)
            ICOL = petsc_offset(MBN)
            ALPBX = -ALPBX
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,ALPBX,
     &                      ADD_VALUES,IERR)
#endif
          ENDIF
!
!---    Loop over borehole nodes  ---
!
        ENDDO
!
!---  Loop over boreholes  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TRNSC_BM_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRNSC_FM_GT( NSL )
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
!     Geothermal Mode (STOMP-GT)
!
!     Solute/reactive species transport fracture to matrix transfer
!     functions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 3 August 2018.
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
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/TRNSC_FM_GT'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fracture triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        IF( IXP_FRC(NTX).EQ.0 ) CYCLE
        MTX = NFLD - NXP + NTX
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
          N = INCM_FRC(NCX)
          IZN = IZ(N)
          MCP = IXP(N)
!
!---      Fracture side gas diffusion  ---
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
!---      Patankar solute transport between fracture and grid cell  ---
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
!---        Fracture equations  ---
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
!---        Fracture equations  ---
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
!---        Fracture equations  ---
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
!---      Fracture side aqueous diffusion  ---
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
!---      Patankar solute transport between fracture and grid celll  ---
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
!---        Fracture equations  ---
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
!---        Fracture equations  ---
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
!---        Fracture equations  ---
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
!---    Loop over fracture triangle to grid cell connections  ---
!
        ENDDO
!
!---  Loop over fracture triangles  ---
!
      ENDDO
!
!---  Loop over fractures  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TRNSC_FM_GT group  ---
!
      RETURN
      END


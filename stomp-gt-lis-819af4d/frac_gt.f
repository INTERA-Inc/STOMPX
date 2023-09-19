!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP_BH_GT( CPGLX,SLX,NBNX )
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
!     Compute the gas/aqueous capillary pressure from the aqueous
!     saturation for a borehole node.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 21 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_BH
      USE GEOM_BH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CAP_BH_GT'
!
!---  Pipe flow  ---
!
      INVX = INV_BH(NBNX)
      IF( IS_BH(INVX).GT.10 ) THEN
        SLRX = PAR_BH(9,INVX)
        ESLX = (SLX-SLRX)/(1.D+0-SLRX)
        CL = MAX( PAR_BH(8,INVX),SMALL )
        SMPX = PAR_BH(11,INVX)
!
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( SLX.LT.SMPX ) THEN
          HMPX = PAR_BH(12,INVX)
          DMPX = -(LOG10(PAR_BH(10,INVX))-LOG10(HMPX))/SMPX
          HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---    Aqueous saturation at or above the matching point,
!       use Brooks and Corey function
!
        ELSE
          HDGL = PAR_BH(7,INVX)*(1.D+0/ESLX)**(1.D+0/CL)
        ENDIF
        CPGLX = HDGL*RHORL*GRAV
!
!---  Borehole filled with porous media flow  ---
!
      ELSE
        IZN = IZ_BH(NBNX)
        ASLMINX = 0.D+0
        SGTX = 0.D+0
        BTGLX = 1.D+0
        CALL CAP_GT( ASLMINX,BTGLX,CPGLX,SLX,SGTX,IZN )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP_FRC_GT( CPGLX,SLX,NFX )
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
!     Compute the gas/aqueous capillary pressure from the aqueous
!     saturation for a fracture triangle.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 7 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
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
      SUB_LOG(ISUB_LOG) = '/CAP_FRC_GT'
      SLRX = RKSP_FRC(4,NFX)
      ESLX = (SLX-SLRX)/(1.D+0-SLRX)
      CL = MAX( RKSP_FRC(3,NFX),SMALL )
      SMPX = RKSP_FRC(6,NFX)
!
!---  Aqueous saturation below the matching point,
!     use Webb extension  ---
!
      IF( SLX.LT.SMPX ) THEN
        HMPX = RKSP_FRC(7,NFX)
        DMPX = -(LOG10(RKSP_FRC(5,NFX))-LOG10(HMPX))/SMPX
        HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---  Aqueous saturation at or above the matching point,
!     use Brooks and Corey function
!
      ELSE
        HDGL = RKSP_FRC(1,NFX)*(1.D+0/ESLX)**(1.D+0/CL)
      ENDIF
      CPGLX = HDGL*RHORL*GRAV
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_BH_GT
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
!     Check the thermodynamic and hydrologic states declared through
!     user inputs for borehole nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, April 21 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE PARM_BH
      USE JACOB
      USE GEOM_BH
      USE FDVP_BH
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK_BH_GT'
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
          XLS_BH(2,NBN) = YLS_BH(2,NBN)
          XLS_BH(1,NBN) = XLS_BH(2,NBN)
!
!---      Skip for pipe flow  ---
!
          IF( IS_BH(INV_BH(NBN)).GT.10 ) CYCLE
          IZN = IZ_BH(NBN)
!
!---      Check for conflicts in saturation functions with
!         gas entrapment and gas relative permeability functions
!         with residual gas saturaitons  ---
!
          IF( ISCHR(IZN).EQ.12 .OR. ISCHR(IZN).EQ.13 .OR.
     &      ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 .OR.
     &      ISCHR(IZN).EQ.201 .OR. ISCHR(IZN).EQ.202 ) THEN
!
!---        Modified-Corey gas relative permeability function  ---
!
            IF( IRPG(IZN).EQ.3 ) THEN
              IF( RPGC(3,IZN).GT.EPSL ) THEN
                INDX = 17
                IMSGX = 0
                NMSGX = 0
                SUBLOGX = 'CHK_BH_GT'
                RLMSGX = RPGC(3,IZN)
                CHMSGX = 'Modified-Corey Gas Relative Permeability: ' //
     &            'Non-zero Residual Gas Saturation with '  //
     &            'Gas Entrapment in Saturation Function'
                CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
              ENDIF
!
!---        Free-Corey gas relative permeability function  ---
!
            ELSEIF( IRPG(IZN).EQ.7 ) THEN
              IF( RPGC(4,IZN).GT.EPSL ) THEN
                INDX = 17
                IMSGX = 0
                NMSGX = 0
                SUBLOGX = 'CHK_BH_GT'
                RLMSGX = RPGC(4,IZN)
                CHMSGX = 'Free-Corey Gas Relative Permeability: ' //
     &            'Non-zero Residual Gas Saturation with '  //
     &            'Gas Entrapment in Saturation Function'
                CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
              ENDIF
!
!---         van Genuchten gas relative permeability function  ---
!
            ELSEIF( IRPG(IZN).EQ.9 ) THEN
              IF( RPGC(3,IZN).GT.EPSL ) THEN
                INDX = 17
                IMSGX = 0
                NMSGX = 0
                SUBLOGX = 'CHK_BH_GT'
                RLMSGX = RPGC(4,IZN)
                CHMSGX = 'van Genuchten Gas Relative Permeability: ' //
     &            'Non-zero Residual Gas Saturation with '  //
     &            'Gas Entrapment in Saturation Function'
                CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
              ENDIF
!
!---        Classical-Corey gas relative permeability function  ---
!
            ELSEIF( IRPG(IZN).EQ.17 ) THEN
              IF( RPGC(4,IZN).GT.EPSL ) THEN
                INDX = 17
                IMSGX = 0
                NMSGX = 0
                SUBLOGX = 'CHK_BH_GT'
                RLMSGX = RPGC(3,IZN)
                CHMSGX = 'Classical-Corey Gas Relative ' //
     &            'Permeability: Non-zero Residual Gas Saturation '  //
     &            'with Gas Entrapment in Saturation Function'
                CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
              ENDIF
            ENDIF
          ENDIF
!
!---      Webb saturation and capillary pressure matching points  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
!
!---        van Genuchten moisture retension function  ---
!
            IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 .OR.
     &        ISCHR(IZN).EQ.201 ) THEN
              CNX = MAX( SCHR(3,IZN),SMALL )
              IF( SCHR(14,IZN).LE.0.D+0 ) THEN
                IF( IRPN(IZN).EQ.2 ) THEN
                  SCHR(14,IZN) = 1.D+0 - 2.D+0/CNX
                ELSE
                  SCHR(14,IZN) = 1.D+0 - 1.D+0/CNX
                ENDIF
              ENDIF
              SRX = SCHR(4,IZN)
              ALPHAX = SCHR(1,IZN)
              CMX = SCHR(14,IZN)
              CALL WEBB_VG( ALPHAX,CMX,CNX,HMPX,SMPX,SRX )
              SCHR(8,IZN) = SMPX
              SCHR(9,IZN) = HMPX
              IF( ISCHR(IZN).EQ.201 ) THEN
                CNX = MAX( SCHR(5,IZN),SMALL )
                IF( SCHR(13,IZN).LE.0.D+0 ) THEN
                  IF( IRPN(IZN).EQ.2 ) THEN
                    SCHR(13,IZN) = 1.D+0 - 2.D+0/CNX
                  ELSE
                    SCHR(13,IZN) = 1.D+0 - 1.D+0/CNX
                  ENDIF
                ENDIF
                ALPHAX = SCHR(2,IZN)
                CMX = SCHR(13,IZN)
                CALL WEBB_VG( ALPHAX,CMX,CNX,HMPX,SMPX,SRX )
                SCHR(10,IZN) = SMPX
                SCHR(11,IZN) = HMPX
              ENDIF
!
!---        Brooks and Corey moisture retension function  ---
!
            ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &        ISCHR(IZN).EQ.202 ) THEN
              SRX = SCHR(4,IZN)
              PSIX = SCHR(1,IZN)
              CLX = SCHR(3,IZN)
              CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
              SCHR(8,IZN) = SMPX
              SCHR(9,IZN) = HMPX
              IF( ISCHR(IZN).EQ.202 ) THEN
                PSIX = SCHR(5,IZN)
                CLX = SCHR(6,IZN)
                CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
                SCHR(10,IZN) = SMPX
                SCHR(11,IZN) = HMPX
              ENDIF
!
!---        Dual porosity van Genuchten moisture retension function  ---
!
            ELSEIF( ISCHR(IZN).EQ.3 ) THEN
!
!---          Matrix matching point  ---
!
              CNX = MAX( SCHR(3,IZN),SMALL )
              IF( SCHR(14,IZN).LE.0.D+0 ) THEN
                IF( IRPN(IZN).EQ.2 ) THEN
                  SCHR(14,IZN) = 1.D+0 - 2.D+0/CNX
                ELSE
                  SCHR(14,IZN) = 1.D+0 - 1.D+0/CNX
                ENDIF
              ENDIF
              SRX = SCHR(4,IZN)
              ALPHAX = SCHR(1,IZN)
              CMX = SCHR(14,IZN)
              CALL WEBB_VG( ALPHAX,CMX,CNX,HMPX,SMPX,SRX )
              SCHR(8,IZN) = SMPX
              SCHR(9,IZN) = HMPX
!
!---          Frature matching point  ---
!
              CNX = MAX( SCHR(6,IZN),SMALL )
              IF( SCHR(15,IZN).LE.0.D+0 ) THEN
                IF( IRPN(IZN).EQ.2 ) THEN
                  SCHR(15,IZN) = 1.D+0 - 2.D+0/CNX
                ELSE
                  SCHR(15,IZN) = 1.D+0 - 1.D+0/CNX
                ENDIF
              ENDIF
              SRX = SCHR(7,IZN)
              ALPHAX = SCHR(5,IZN)
              CMX = SCHR(15,IZN)
              CALL WEBB_VG( ALPHAX,CMX,CNX,HMPX,SMPX,SRX )
              SCHR(10,IZN) = SMPX
              SCHR(11,IZN) = HMPX
!
!---        Dual porosity Brooks and Corey moisture retension
!           function  ---
!
            ELSEIF( ISCHR(IZN).EQ.4 ) THEN
!
!---          Matrix matching point  ---
!
              SRX = SCHR(4,IZN)
              PSIX = SCHR(1,IZN)
              CLX = SCHR(3,IZN)
              CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
              SCHR(8,IZN) = SMPX
              SCHR(9,IZN) = HMPX
!
!---          Matrix matching point  ---
!
              SRX = SCHR(7,IZN)
              PSIX = SCHR(5,IZN)
              CLX = SCHR(6,IZN)
              CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
              SCHR(10,IZN) = SMPX
              SCHR(11,IZN) = HMPX
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
!---  Skip for restart simulations  ---
!
      IF( IEO.EQ.2 .AND. ISLC(87).EQ.0 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Inner Aqueous Mass Injection, Outer Aqueous Pressure  ---
!
        IF( IC_BH(NBH).EQ.11 ) THEN
!
!---      Loop over the boreholes nodes starting with the
!         last borehole node (i.e., outer coaxial node)  ---
!
          DO NBN1 = ID_BH(4,NBH),ID_BH(3,NBH),-1
            INV1 = INV_BH(NBN1)
!
!---        Borehole node midpoint  ---
!
            XP_BHX = 5.D-1*(XP_BH(1,NBN1)+XP_BH(2,NBN1))
            YP_BHX = 5.D-1*(YP_BH(1,NBN1)+YP_BH(2,NBN1))
            ZP_BHX = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
!
!---        Temperature at borehole node centroid  ---
!
            IF( VIC_BH(12,NBH).GT.-TABS )
     &        T_BH(2,NBN1) = VIC_BH(12,NBH) +
     &          (XP_BHX-VIC_BH(1,NBH))*VIC_BH(13,NBH) +
     &          (YP_BHX-VIC_BH(2,NBH))*VIC_BH(14,NBH) +
     &          (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(15,NBH)
!
!---        Last borehole node only (i.e., outer coaxial node)  ---
!
            IF( NBN1.EQ.ID_BH(4,NBH) ) THEN
              ZP_BH1X = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
!
!---          Aqueous pressure at borehole node ending point  ---
!
              P_BHX = VIC_BH(4,NBH) - PATM
              PX = VIC_BH(4,NBH)
              CALL P_IAPWS( T_BH(2,NBN1),PX,RHOG_BHX,RHOL_BHX,
     &          HGWX,HLWX,UGWX,ULWX )
!
!---          Aqueous viscosity  ---
!
              CALL VISC_W( T_BH(2,NBN1),PX,RHOL_BHX,VISL_BHX )
!
!---          Hydraulic diameter (m)  ---
!
              DHX = 2.D+0*SQRT(PAR_BH(3,INV1)**2 - PAR_BH(2,INV1)**2)
!
!---          Hydraulic area (m^2) ---
!
              AHX = GPI*((5.D-1*DHX)**2)
!
!---          Roughness, m  ---
!
              RFX = PAR_BH(16,INV1)
!
!---          Aqueous velocity, m/s  ---
!
              ULX = VIC_BH(5,NBH)/(RHOL_BHX*AHX)
!
!---          Half borehole length, m  ---
!
              DIST1X = SQRT( (XP_BH(2,NBN1)-XP_BH(1,NBN1))**2 + 
     &          (YP_BH(2,NBN1)-YP_BH(1,NBN1))**2 + 
     &          (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))**2 )
              DBB1X = 5.D-1*DIST1X
!
!---          Compute pipe pressure drop  ---
!
              CALL PIPE_FLOW_DP_GT( DBB1X,DHX,HDLX,RFX,RHOL_BHX,
     &          ULX,VISL_BHX,NBN1 )
!
!---          Pressure at borehole node equal to sum of pressure at
!             outlet plus pressure drop plus hydrostatic pressure  ---
!
              PL_BH(2,NBN1) = P_BHX - HDLX + 
     &          (ZP_BH(2,NBN1)-ZP_BH1X)*RHOL_BHX*GRAV
              PG_BH(2,NBN1) = PL_BH(2,NBN1)
              SL_BH(2,NBN1) = 1.D+0
              YLS_BH(2,NBN1) = 0.D+0
              NPHAZ_BH(2,NBN1) = 1
            ELSE
              NCX = IPB_BH(2,NBN1)
              NBN2 = IBCM_BH(NCX)
              ZP_BH1X = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
              ZP_BH2X = 5.D-1*(ZP_BH(1,NBN2)+ZP_BH(2,NBN2))
              INV2 = INV_BH(NBN2)
              DBB1X = DBBM_BH(NCX)
              DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---          Aqueous pressure at borehole node ending point  ---
!
              P_BHX = PL_BH(2,NBN2)
              PX = PL_BH(2,NBN2) + PATM
              CALL P_IAPWS( T_BH(2,NBN1),PX,RHOG_BHX,RHOL_BHX,
     &          HGWX,HLWX,UGWX,ULWX )
!
!---          Aqueous viscosity  ---
!
              CALL VISC_W( T_BH(2,NBN1),PX,RHOL_BHX,VISL_BHX )
!
!---          Outer pipe, hydraulic diameter (m) and roughness (m)  ---
!
              IF( IBN_BH(NBN1).NE.0 ) THEN
                DH1X = 2.D+0*SQRT(PAR_BH(3,INV1)**2 - PAR_BH(2,INV1)**2)
                RF1X = PAR_BH(16,INV1)
!
!---          Inner pipe, hydraulic diameter (m) and roughness (m)  ---
!
              ELSE
                DH1X = 2.D+0*PAR_BH(1,INV1)
                RF1X = PAR_BH(13,INV1)
              ENDIF
!
!---          Outer pipe, hydraulic diameter (m) and roughness (m)  ---
!
              IF( IBN_BH(NBN2).NE.0 ) THEN
                DH2X = 2.D+0*SQRT(PAR_BH(3,INV2)**2 - PAR_BH(2,INV2)**2)
                RF2X = PAR_BH(16,INV2)
!
!---          Inner pipe, hydraulic diameter (m) and roughness (m)  ---
!
              ELSE
                DH2X = 2.D+0*PAR_BH(1,INV2)
                RF2X = PAR_BH(13,INV2)
              ENDIF
              INDX = -3
              DHX = DIFMN( DH1X,DH2X,DBB1X,DBB2X,ZERO,INDX )
              RFX = DIFMN( RF1X,RF2X,DBB1X,DBB2X,ZERO,INDX )
!
!---          Hydraulic area (m^2) ---
!
              AHX = GPI*((5.D-1*DHX)**2)
!
!---          Aqueous velocity, m/s  ---
!
              ULX = VIC_BH(5,NBH)/(RHOL_BHX*AHX)
!
!---          Compute pipe pressure drop  ---
!
              CALL PIPE_FLOW_DP_GT( DBB_BH(NCX),DHX,HDLX,RFX,RHOL_BHX,
     &          ULX,VISL_BHX,NBN1 )
!
!---          Pressure at borehole node equal to sum of pressure at
!             outlet plus pressure drop plus hydrostatic pressure  ---
!
              PL_BH(2,NBN1) = PL_BH(2,NBN2) - HDLX + 
     &          (ZP_BH2X-ZP_BH1X)*RHOL_BHX*GRAV
              PG_BH(2,NBN1) = PL_BH(2,NBN1)
              SL_BH(2,NBN1) = 1.D+0
              YLS_BH(2,NBN1) = 0.D+0
              NPHAZ_BH(2,NBN1) = 1
            ENDIF
          ENDDO
          CYCLE
!
!---    Outer Aqueous Mass Injection, Inner Aqueous Pressure  ---
!
        ELSEIF( IC_BH(NBH).EQ.12 ) THEN
!
!---      Loop over the boreholes nodes starting with the
!         first borehole node (i.e., inner coaxial node)  ---
!
          DO NBN1 = ID_BH(3,NBH),ID_BH(4,NBH)
            INV1 = INV_BH(NBN1)
!
!---        Borehole node midpoint  ---
!
            XP_BHX = 5.D-1*(XP_BH(1,NBN1)+XP_BH(2,NBN1))
            YP_BHX = 5.D-1*(YP_BH(1,NBN1)+YP_BH(2,NBN1))
            ZP_BHX = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
!
!---        Temperature at borehole node centroid  ---
!
            IF( VIC_BH(12,NBH).GT.-TABS )
     &        T_BH(2,NBN1) = VIC_BH(12,NBH) +
     &          (XP_BHX-VIC_BH(1,NBH))*VIC_BH(13,NBH) +
     &          (YP_BHX-VIC_BH(2,NBH))*VIC_BH(14,NBH) +
     &          (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(15,NBH)
!
!---        First borehole node only (i.e., inner coaxial node)  ---
!
            IF( NBN1.EQ.ID_BH(3,NBH) ) THEN
              ZP_BH1X = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
!
!---          Aqueous pressure at borehole node ending point  ---
!
              P_BHX = VIC_BH(4,NBH) - PATM
              PX = VIC_BH(4,NBH)
              CALL P_IAPWS( T_BH(2,NBN1),PX,RHOG_BHX,RHOL_BHX,
     &          HGWX,HLWX,UGWX,ULWX )
!
!---          Aqueous viscosity  ---
!
              CALL VISC_W( T_BH(2,NBN1),PX,RHOL_BHX,VISL_BHX )
!
!---          Hydraulic diameter (m) ---
!
              DHX = 2.D+0*PAR_BH(1,INV1)
!
!---          Hydraulic area (m^2) ---
!
              AHX = GPI*((5.D-1*DHX)**2)
!
!---          Roughness, m  ---
!
              RFX = PAR_BH(13,INV1)
!
!---          Aqueous velocity, m/s  ---
!
              ULX = VIC_BH(5,NBH)/(RHOL_BHX*AHX)
!
!---          Half borehole length, m  ---
!
              DIST1X = SQRT( (XP_BH(2,NBN1)-XP_BH(1,NBN1))**2 + 
     &          (YP_BH(2,NBN1)-YP_BH(1,NBN1))**2 + 
     &          (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))**2 )
              DBB1X = 5.D-1*DIST1X
!
!---          Compute pipe pressure drop  ---
!
              CALL PIPE_FLOW_DP_GT( DBB1X,DHX,HDLX,RFX,RHOL_BHX,
     &          ULX,VISL_BHX,NBN1 )
!
!---          Pressure at borehole node equal to sum of pressure at
!             outlet plus pressure drop plus hydrostatic pressure  ---
!
              PL_BH(2,NBN1) = P_BHX + HDLX + 
     &          (ZP_BH(1,NBN1)-ZP_BH1X)*RHOL_BHX*GRAV
              PG_BH(2,NBN1) = PL_BH(2,NBN1)
              SL_BH(2,NBN1) = 1.D+0
              YLS_BH(2,NBN1) = 0.D+0
              NPHAZ_BH(2,NBN1) = 1
            ELSE
              NCX = IPB_BH(1,NBN1)
              NBN2 = IBCM_BH(NCX)
              ZP_BH1X = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
              ZP_BH2X = 5.D-1*(ZP_BH(1,NBN2)+ZP_BH(2,NBN2))
              INV2 = INV_BH(NBN2)
              DBB1X = DBBM_BH(NCX)
              DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---          Aqueous pressure at borehole node ending point  ---
!
              P_BHX = PL_BH(2,NBN2)
              PX = PL_BH(2,NBN2) + PATM
              CALL P_IAPWS( T_BH(2,NBN1),PX,RHOG_BHX,RHOL_BHX,
     &          HGWX,HLWX,UGWX,ULWX )
!
!---          Aqueous viscosity  ---
!
              CALL VISC_W( T_BH(2,NBN1),PX,RHOL_BHX,VISL_BHX )
!
!---          Outer pipe, hydraulic diameter (m) and roughness (m)  ---
!
              IF( IBN_BH(NBN1).NE.0 ) THEN
                DH1X = 2.D+0*SQRT(PAR_BH(3,INV1)**2 - PAR_BH(2,INV1)**2)
                RF1X = PAR_BH(16,INV1)
!
!---          Inner pipe, hydraulic diameter (m) and roughness (m)  ---
!
              ELSE
                DH1X = 2.D+0*PAR_BH(1,INV1)
                RF1X = PAR_BH(13,INV1)
              ENDIF
!
!---          Outer pipe, hydraulic diameter (m) and roughness (m)  ---
!
              IF( IBN_BH(NBN2).NE.0 ) THEN
                DH2X = 2.D+0*SQRT(PAR_BH(3,INV2)**2 - PAR_BH(2,INV2)**2)
                RF2X = PAR_BH(16,INV2)
!
!---          Inner pipe, hydraulic diameter (m) and roughness (m)  ---
!
              ELSE
                DH2X = 2.D+0*PAR_BH(1,INV2)
                RF2X = PAR_BH(13,INV2)
              ENDIF
              INDX = -3
              DHX = DIFMN( DH1X,DH2X,DBB1X,DBB2X,ZERO,INDX )
              RFX = DIFMN( RF1X,RF2X,DBB1X,DBB2X,ZERO,INDX )
!
!---          Hydraulic area (m^2) ---
!
              AHX = GPI*((5.D-1*DHX)**2)
!
!---          Aqueous velocity, m/s  ---
!
              ULX = VIC_BH(5,NBH)/(RHOL_BHX*AHX)
!
!---          Compute pipe pressure drop  ---
!
              CALL PIPE_FLOW_DP_GT( DBB_BH(NCX),DHX,HDLX,RFX,RHOL_BHX,
     &          ULX,VISL_BHX,NBN1 )
!
!---          Pressure at borehole node equal to sum of pressure at
!             outlet plus pressure drop plus hydrostatic pressure  ---
!
              PL_BH(2,NBN1) = PL_BH(2,NBN2) + HDLX + 
     &          (ZP_BH2X-ZP_BH1X)*RHOL_BHX*GRAV
              PG_BH(2,NBN1) = PL_BH(2,NBN1)
              SL_BH(2,NBN1) = 1.D+0
              YLS_BH(2,NBN1) = 0.D+0
              NPHAZ_BH(2,NBN1) = 1
            ENDIF
          ENDDO
          CYCLE
!
!---    Start Aqueous Mass Injection, End Aqueous Pressure  ---
!
        ELSEIF( IC_BH(NBH).EQ.13 ) THEN
!
!---      Loop over the boreholes nodes starting with the
!         last borehole node  ---
!
          DO NBN1 = ID_BH(4,NBH),ID_BH(3,NBH),-1
            INV1 = INV_BH(NBN1)
!
!---        Borehole node midpoint  ---
!
            XP_BHX = 5.D-1*(XP_BH(1,NBN1)+XP_BH(2,NBN1))
            YP_BHX = 5.D-1*(YP_BH(1,NBN1)+YP_BH(2,NBN1))
            ZP_BHX = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
!
!---        Temperature at borehole node centroid  ---
!
            IF( VIC_BH(12,NBH).GT.-TABS )
     &        T_BH(2,NBN1) = VIC_BH(12,NBH) +
     &          (XP_BHX-VIC_BH(1,NBH))*VIC_BH(13,NBH) +
     &          (YP_BHX-VIC_BH(2,NBH))*VIC_BH(14,NBH) +
     &          (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(15,NBH)
!
!---        Last borehole node only  ---
!
            IF( NBN1.EQ.ID_BH(4,NBH) ) THEN
              ZP_BH1X = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
!
!---          Aqueous pressure at borehole node ending point  ---
!
              P_BHX = VIC_BH(4,NBH) - PATM
              PX = VIC_BH(4,NBH)
              CALL P_IAPWS( T_BH(2,NBN1),PX,RHOG_BHX,RHOL_BHX,
     &          HGWX,HLWX,UGWX,ULWX )
!
!---          Aqueous viscosity  ---
!
              CALL VISC_W( T_BH(2,NBN1),PX,RHOL_BHX,VISL_BHX )
!
!---          Hydraulic diameter (m)  ---
!
              DHX = 2.D+0*PAR_BH(3,INV1)
!
!---          Hydraulic area (m^2) ---
!
              AHX = GPI*((5.D-1*DHX)**2)
!
!---          Roughness, m  ---
!
              RFX = PAR_BH(6,INV1)
!
!---          Aqueous velocity, m/s  ---
!
              ULX = VIC_BH(5,NBH)/(RHOL_BHX*AHX)
!
!---          Half borehole length, m  ---
!
              DIST1X = SQRT( (XP_BH(2,NBN1)-XP_BH(1,NBN1))**2 + 
     &          (YP_BH(2,NBN1)-YP_BH(1,NBN1))**2 + 
     &          (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))**2 )
              DBB1X = 5.D-1*DIST1X
!
!---          Compute pipe pressure drop  ---
!
              CALL PIPE_FLOW_DP_GT( DBB1X,DHX,HDLX,RFX,RHOL_BHX,
     &          ULX,VISL_BHX,NBN1 )
!
!---          Pressure at borehole node equal to sum of pressure at
!             outlet plus pressure drop plus hydrostatic pressure  ---
!
              PL_BH(2,NBN1) = P_BHX - HDLX + 
     &          (ZP_BH(2,NBN1)-ZP_BH1X)*RHOL_BHX*GRAV
              PG_BH(2,NBN1) = PL_BH(2,NBN1)
              SL_BH(2,NBN1) = 1.D+0
              YLS_BH(2,NBN1) = 0.D+0
              NPHAZ_BH(2,NBN1) = 1
            ELSE
              NCX = IPB_BH(2,NBN1)
              NBN2 = IBCM_BH(NCX)
              ZP_BH1X = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
              ZP_BH2X = 5.D-1*(ZP_BH(1,NBN2)+ZP_BH(2,NBN2))
              INV2 = INV_BH(NBN2)
              DBB1X = DBBM_BH(NCX)
              DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---          Aqueous pressure at borehole node ending point  ---
!
              P_BHX = PL_BH(2,NBN2)
              PX = PL_BH(2,NBN2) + PATM
              CALL P_IAPWS( T_BH(2,NBN1),PX,RHOG_BHX,RHOL_BHX,
     &          HGWX,HLWX,UGWX,ULWX )
!
!---          Aqueous viscosity  ---
!
              CALL VISC_W( T_BH(2,NBN1),PX,RHOL_BHX,VISL_BHX )
!
!---          Hydraulic diameter (m) and roughness (m)  ---
!
              DH1X = 2.D+0*PAR_BH(3,INV1)
              RF1X = PAR_BH(6,INV1)
!
!---          Hydraulic diameter (m) and roughness (m)  ---
!
              DH2X = 2.D+0*PAR_BH(3,INV2)
              RF2X = PAR_BH(6,INV2)
              INDX = -3
              DHX = DIFMN( DH1X,DH2X,DBB1X,DBB2X,ZERO,INDX )
              RFX = DIFMN( RF1X,RF2X,DBB1X,DBB2X,ZERO,INDX )
!
!---          Hydraulic area (m^2) ---
!
              AHX = GPI*((5.D-1*DHX)**2)
!
!---          Aqueous velocity, m/s  ---
!
              ULX = VIC_BH(5,NBH)/(RHOL_BHX*AHX)
!
!---          Compute pipe pressure drop  ---
!
              CALL PIPE_FLOW_DP_GT( DBB_BH(NCX),DHX,HDLX,RFX,RHOL_BHX,
     &          ULX,VISL_BHX,NBN1 )
!
!---          Pressure at borehole node equal to sum of pressure at
!             outlet plus pressure drop plus hydrostatic pressure  ---
!
              PL_BH(2,NBN1) = PL_BH(2,NBN2) - HDLX + 
     &          (ZP_BH2X-ZP_BH1X)*RHOL_BHX*GRAV
              PG_BH(2,NBN1) = PL_BH(2,NBN1)
              SL_BH(2,NBN1) = 1.D+0
              YLS_BH(2,NBN1) = 0.D+0
              NPHAZ_BH(2,NBN1) = 1
            ENDIF
          ENDDO
          CYCLE
!
!---    Start Aqueous Mass Injection, End Aqueous Pressure  ---
!
        ELSEIF( IC_BH(NBH).EQ.14 ) THEN
!
!---      Loop over the boreholes nodes starting with the
!         first borehole node  ---
!
          DO NBN1 = ID_BH(3,NBH),ID_BH(4,NBH)
            INV1 = INV_BH(NBN1)
!
!---        Borehole node midpoint  ---
!
            XP_BHX = 5.D-1*(XP_BH(1,NBN1)+XP_BH(2,NBN1))
            YP_BHX = 5.D-1*(YP_BH(1,NBN1)+YP_BH(2,NBN1))
            ZP_BHX = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
!
!---        Temperature at borehole node centroid  ---
!
            IF( VIC_BH(12,NBH).GT.-TABS )
     &        T_BH(2,NBN1) = VIC_BH(12,NBH) +
     &          (XP_BHX-VIC_BH(1,NBH))*VIC_BH(13,NBH) +
     &          (YP_BHX-VIC_BH(2,NBH))*VIC_BH(14,NBH) +
     &          (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(15,NBH)
!
!---        First borehole node only  ---
!
            IF( NBN1.EQ.ID_BH(3,NBH) ) THEN
              ZP_BH1X = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
!
!---          Aqueous pressure at borehole node ending point  ---
!
              P_BHX = VIC_BH(4,NBH) - PATM
              PX = VIC_BH(4,NBH)
              CALL P_IAPWS( T_BH(2,NBN1),PX,RHOG_BHX,RHOL_BHX,
     &          HGWX,HLWX,UGWX,ULWX )
!
!---          Aqueous viscosity  ---
!
              CALL VISC_W( T_BH(2,NBN1),PX,RHOL_BHX,VISL_BHX )
!
!---          Hydraulic diameter (m) ---
!
              DHX = 2.D+0*PAR_BH(3,INV1)
!
!---          Hydraulic area (m^2) ---
!
              AHX = GPI*((5.D-1*DHX)**2)
!
!---          Roughness, m  ---
!
              RFX = PAR_BH(6,INV1)
!
!---          Aqueous velocity, m/s  ---
!
              ULX = VIC_BH(5,NBH)/(RHOL_BHX*AHX)
!
!---          Half borehole length, m  ---
!
              DIST1X = SQRT( (XP_BH(2,NBN1)-XP_BH(1,NBN1))**2 + 
     &          (YP_BH(2,NBN1)-YP_BH(1,NBN1))**2 + 
     &          (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))**2 )
              DBB1X = 5.D-1*DIST1X
!
!---          Compute pipe pressure drop  ---
!
              CALL PIPE_FLOW_DP_GT( DBB1X,DHX,HDLX,RFX,RHOL_BHX,
     &          ULX,VISL_BHX,NBN1 )
!
!---          Pressure at borehole node equal to sum of pressure at
!             outlet plus pressure drop plus hydrostatic pressure  ---
!
              PL_BH(2,NBN1) = P_BHX + HDLX + 
     &          (ZP_BH(1,NBN1)-ZP_BH1X)*RHOL_BHX*GRAV
              PG_BH(2,NBN1) = PL_BH(2,NBN1)
              SL_BH(2,NBN1) = 1.D+0
              YLS_BH(2,NBN1) = 0.D+0
              NPHAZ_BH(2,NBN1) = 1
            ELSE
              NCX = IPB_BH(1,NBN1)
              NBN2 = IBCM_BH(NCX)
              ZP_BH1X = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
              ZP_BH2X = 5.D-1*(ZP_BH(1,NBN2)+ZP_BH(2,NBN2))
              INV2 = INV_BH(NBN2)
              DBB1X = DBBM_BH(NCX)
              DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---          Aqueous pressure at borehole node ending point  ---
!
              P_BHX = PL_BH(2,NBN2)
              PX = PL_BH(2,NBN2) + PATM
              CALL P_IAPWS( T_BH(2,NBN1),PX,RHOG_BHX,RHOL_BHX,
     &          HGWX,HLWX,UGWX,ULWX )
!
!---          Aqueous viscosity  ---
!
              CALL VISC_W( T_BH(2,NBN1),PX,RHOL_BHX,VISL_BHX )
!
!---          Hydraulic diameter (m) and roughness (m)  ---
!
              DH1X = 2.D+0*PAR_BH(3,INV1)
              RF1X = PAR_BH(6,INV1)
!
!---          Hydraulic diameter (m) and roughness (m)  ---
!
              DH2X = 2.D+0*PAR_BH(3,INV2)
              RF2X = PAR_BH(6,INV2)
              INDX = -3
              DHX = DIFMN( DH1X,DH2X,DBB1X,DBB2X,ZERO,INDX )
              RFX = DIFMN( RF1X,RF2X,DBB1X,DBB2X,ZERO,INDX )
!
!---          Hydraulic area (m^2) ---
!
              AHX = GPI*((5.D-1*DHX)**2)
!
!---          Aqueous velocity, m/s  ---
!
              ULX = VIC_BH(5,NBH)/(RHOL_BHX*AHX)
!
!---          Compute pipe pressure drop  ---
!
              CALL PIPE_FLOW_DP_GT( DBB_BH(NCX),DHX,HDLX,RFX,RHOL_BHX,
     &          ULX,VISL_BHX,NBN1 )
!
!---          Pressure at borehole node equal to sum of pressure at
!             outlet plus pressure drop plus hydrostatic pressure  ---
!
              PL_BH(2,NBN1) = PL_BH(2,NBN2) + HDLX + 
     &          (ZP_BH2X-ZP_BH1X)*RHOL_BHX*GRAV
              PG_BH(2,NBN1) = PL_BH(2,NBN1)
              SL_BH(2,NBN1) = 1.D+0
              YLS_BH(2,NBN1) = 0.D+0
              NPHAZ_BH(2,NBN1) = 1
            ENDIF
          ENDDO
          CYCLE
        ENDIF
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
          IZN = IZ_BH(NBN)
!
!---      Borehole node midpoint  ---
!
          XP_BHX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
          YP_BHX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
          ZP_BHX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
!
!---      Dissolved air relative saturation  ---
!
          PVA_BH(2,NBN) = VIC_BH(20,NBH)
!
!---      Hydrostatic initial conditions  ---
!
          IF( IC_BH(NBH).EQ.4 ) THEN
            CALL HYDST_BH_GT( VIC_BH(1,NBH),PG_BH(2,NBN),
     &        PL_BH(2,NBN),T_BH(2,NBN),XLS_BH(2,NBN),YLS_BH(2,NBN),
     &        ZP_BHX )
!
!---      Aqueous saturation + gas pressure initial conditions  ---
!
          ELSEIF( IC_BH(NBH).EQ.1 ) THEN
            PG_BH(2,NBN) = VIC_BH(4,NBH) +
     &        (XP_BHX-VIC_BH(1,NBH))*VIC_BH(5,NBH) +
     &        (YP_BHX-VIC_BH(2,NBH))*VIC_BH(6,NBH) +
     &        (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(7,NBH)
            SL_BH(2,NBN) = VIC_BH(8,NBH) +
     &        (XP_BHX-VIC_BH(1,NBH))*VIC_BH(9,NBH) +
     &        (YP_BHX-VIC_BH(2,NBH))*VIC_BH(10,NBH) +
     &        (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(11,NBH)
!
!---        Temperature at borehole node centroid  ---
!
            IF( VIC_BH(12,NBH).GT.-TABS )
     &        T_BH(2,NBN) = VIC_BH(12,NBH) +
     &          (XP_BHX-VIC_BH(1,NBH))*VIC_BH(13,NBH) +
     &          (YP_BHX-VIC_BH(2,NBH))*VIC_BH(14,NBH) +
     &          (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(15,NBH)
            YLS_BH(2,NBN) = VIC_BH(16,NBH) +
     &        (XP_BHX-VIC_BH(1,NBH))*VIC_BH(17,NBH) +
     &        (YP_BHX-VIC_BH(2,NBH))*VIC_BH(18,NBH) +
     &        (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(19,NBH)
            CALL CAP_BH_GT( CPGLX,SL_BH(2,NBN),NBN )
            PL_BH(2,NBN) = PG_BH(2,NBN) - CPGLX
!
!---      Aqueous saturation + aqueous pressure initial conditions  ---
!
          ELSEIF( IC_BH(NBH).EQ.2 ) THEN
            PL_BH(2,NBN) = VIC_BH(4,NBH) +
     &        (XP_BHX-VIC_BH(1,NBH))*VIC_BH(5,NBH) +
     &        (YP_BHX-VIC_BH(2,NBH))*VIC_BH(6,NBH) +
     &        (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(7,NBH)
            SL_BH(2,NBN) = VIC_BH(8,NBH) +
     &        (XP_BHX-VIC_BH(1,NBH))*VIC_BH(9,NBH) +
     &        (YP_BHX-VIC_BH(2,NBH))*VIC_BH(10,NBH) +
     &        (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(11,NBH)
!
!---        Temperature at borehole node centroid  ---
!
            IF( VIC_BH(12,NBH).GT.-TABS )
     &        T_BH(2,NBN) = VIC_BH(12,NBH) +
     &          (XP_BHX-VIC_BH(1,NBH))*VIC_BH(13,NBH) +
     &          (YP_BHX-VIC_BH(2,NBH))*VIC_BH(14,NBH) +
     &          (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(15,NBH)
            YLS_BH(2,NBN) = VIC_BH(16,NBH) +
     &        (XP_BHX-VIC_BH(1,NBH))*VIC_BH(17,NBH) +
     &        (YP_BHX-VIC_BH(2,NBH))*VIC_BH(18,NBH) +
     &        (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(19,NBH)
            CALL CAP_BH_GT( CPGLX,SL_BH(2,NBN),NBN )
            PG_BH(2,NBN) = PL_BH(2,NBN) + CPGLX
!
!---      Aqueous pressure + gas pressure initial conditions  ---
!
          ELSEIF( IC_BH(NBH).EQ.3 ) THEN
            PG_BH(2,NBN) = VIC_BH(4,NBH) +
     &        (XP_BHX-VIC_BH(1,NBH))*VIC_BH(5,NBH) +
     &        (YP_BHX-VIC_BH(2,NBH))*VIC_BH(6,NBH) +
     &        (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(7,NBH)
            PL_BH(2,NBN) = VIC_BH(8,NBH) +
     &        (XP_BHX-VIC_BH(1,NBH))*VIC_BH(9,NBH) +
     &        (YP_BHX-VIC_BH(2,NBH))*VIC_BH(10,NBH) +
     &        (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(11,NBH)
!
!---        Temperature at borehole node centroid  ---
!
            IF( VIC_BH(12,NBH).GT.-TABS )
     &        T_BH(2,NBN) = VIC_BH(12,NBH) +
     &          (XP_BHX-VIC_BH(1,NBH))*VIC_BH(13,NBH) +
     &          (YP_BHX-VIC_BH(2,NBH))*VIC_BH(14,NBH) +
     &          (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(15,NBH)
            YLS_BH(2,NBN) = VIC_BH(16,NBH) +
     &        (XP_BHX-VIC_BH(1,NBH))*VIC_BH(17,NBH) +
     &        (YP_BHX-VIC_BH(2,NBH))*VIC_BH(18,NBH) +
     &        (ZP_BHX-VIC_BH(3,NBH))*VIC_BH(19,NBH)
          ENDIF
!
!---      Saturation and relative permeability  ---
!
          CALL RKSP_BH_GT( PG_BH(2,NBN),PL_BH(2,NBN),
     &      RKG_BH(2,NBN),RKL_BH(2,NBN),SG_BH(2,NBN),SL_BH(2,NBN),
     &      NBN )
!
!---      Check for out of range parameters  ---
!
          IF( T_BH(2,NBN).GT.2.D+3 .OR. T_BH(2,NBN).LT.0.01D+0 ) THEN
            INDX = 17
            IMSGX = 0
            NMSGX = 0
            SUBLOGX = 'CHK_BH_GT'
            RLMSGX = T_BH(2,NBN)
            CHMSGX = 'Out of Range Initial Fracture Temperature(C) = '
            CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ENDIF
          IF( PL_BH(2,NBN).GT.100.D+6-PATM ) THEN
            INDX = 17
            IMSGX = 0
            NMSGX = 0
            SUBLOGX = 'CHK_BH_GT'
            RLMSGX = PL_BH(2,NBN)+PATM
            CHMSGX = 'Out of Range Initial Fracture Aqu. ' // 
     &        'Pressure(Pa) = '
            CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ENDIF
          IF( PG_BH(2,NBN).GT.100.D+6-PATM .OR.
     &      PG_BH(2,NBN).LT.6.1125D+2-PATM ) THEN
            INDX = 17
            IMSGX = 0
            NMSGX = 0
            SUBLOGX = 'CHK_BH_GT'
            RLMSGX = PG_BH(2,NBN)+PATM
            CHMSGX = 'Out of Range Initial Fracture Gas Pressure(Pa) = '
            CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ENDIF
          IF( SL_BH(2,NBN).GT.1.D+0 .OR. SL_BH(2,NBN).LT.0.D+0 ) THEN
            INDX = 17
            IMSGX = 0
            NMSGX = 0
            SUBLOGX = 'CHK_BH_GT'
            RLMSGX = SL_BH(2,NBN)
            CHMSGX = 'Out of Range Initial Fracture Aqu. Saturation = '
            CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ENDIF
!
!---      Assign initial phase condition  ---
!
          PLX = PL_BH(2,NBN)+PATM
          PGX = PG_BH(2,NBN)+PATM
          TKX = T_BH(2,NBN)+TABS
          TX = T_BH(2,NBN)
!
!---      Saturated conditions with no trapped gas  ---
!
          IF( (1.D+0-SL_BH(2,NBN))/EPSL.LT.EPSL ) THEN
            NPHAZ_BH(2,NBN) = 1
!
!---      Unsaturated conditions  ---
!
          ELSEIF( SL_BH(2,NBN)/EPSL.LT.EPSL .OR. TKX.GT.TCRW ) THEN
            NPHAZ_BH(2,NBN) = 4
!
!---      Partially saturated conditions  ---
!
          ELSE
            NPHAZ_BH(2,NBN) = 2
          ENDIF
!
!---      Saturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - trapped gas saturation
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ_BH(2,NBN).EQ.1 ) THEN
!
!---        Check dissolved salt mass fraction  ---
!
            CALL SOL_BRNS( TX,PLX,XLSMX )
            IF( YLS_BH(2,NBN).GT.XLSMX ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_BH_GT'
              RLMSGX = XLSMX
              CHMSGX = 'Saturated Fracture Initial Conditions ' //
     &          'Aqu. Salt Mass Frac. > Solubility Limit : ' //
     &          'Aqu. Salt Mass Frac. Solubility Limit = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
            XLS_BH(2,NBN) = YLS_BH(2,NBN)
!
!---        Check pressure exceeds water vapor pressure  ---
!
            CALL SP_B( TX,XLS_BH(2,NBN),PSBX )
            IF( PLX.LT.PSBX ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_BH_GT'
              RLMSGX = PSBX
              CHMSGX = 'Saturated Fracture Initial Conditions ' //
     &          'Aqueous Pressure < Saturated Pressure : ' //
     &          'Saturated Pressure, Pa = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
!
!---        Dissolved air mole fraction  ---
!
            XMLA_BH(2,NBN) = VIC_BH(20,NBH)*PLX/HCAW
!
!---      Unsaturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ_BH(2,NBN).EQ.2 ) THEN
!
!---        Saturated water vapor pressure  ---
!
            INDX = 0
            CALL REGION_4( TX,PSWX,INDX )
!
!---        Check dissolved salt mass fraction  ---
!
            CALL SOL_BRNS( TX,PSWX,XLSMX )
            IF( YLS_BH(2,NBN).GT.XLSMX ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_BH_GT'
              RLMSGX = XLSMX
              CHMSGX = 'Unsaturated Fracture Initial Conditions ' //
     &          'Aqu. Salt Mass Frac. > Solubility Limit : ' //
     &          'Aqu. Salt Mass Frac. Solubility Limit = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
            XLS_BH(2,NBN) = YLS_BH(2,NBN)
!
!---        Check gas pressure exceeds water vapor pressure  ---
!
            CALL P_IAPWS( TX,PSWX,RHOX,RHOLWX,HX,HLWX,UX,ULWX )
            CALL DENS_B( XLS_BH(2,NBN),RHOBX,RHOLWX,TX )
            CALL SP_B( TX,XLS_BH(2,NBN),PSBX )
            PCX = PGX - PLX
            CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX )
            IF( PGX.LT.PVBX ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_BH_GT'
              RLMSGX = PVBX
              CHMSGX = 'Unsaturated Fracture Initial Conditions ' //
     &          'Gas Pressure < Water Vapor Pressure : ' //
     &          'Water Vapor Pressure, Pa = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
!
!---        Dissolved air mole fraction  ---
!
            PVA_BH(2,NBN) = MAX( PGX-PVBX,0.D+0 )
            XMLA_BH(2,NBN) = PVA_BH(2,NBN)/HCAW
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - gas pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ_BH(2,NBN).EQ.4 ) THEN
!
!---        Check for non-zero dissolved salt concentration  ---
!
            IF( YLS_BH(2,NBN).GT.EPSL ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_BH_GT'
              RLMSGX = YLS_BH(2,NBN)
              CHMSGX = 'Fully Unsaturated Fracture Initial' //
     &          'Conditions Non-Zero Dissolved Salt Mass Fraction ' //
     &          'Declared : Salt Concentration = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
            XLS_BH(2,NBN) = 0.D+0
            YLS_BH(2,NBN) = 0.D+0
!
!---        Water vapor pressure  ---
!
            INDX = 0
            CALL REGION_4( TX,PSWX,INDX )
            SL_BH(2,NBN) = 0.D+0
            CALL CAP_BH_GT( CPGLX,SL_BH(2,NBN),NBN )
            PLX = PGX - CPGLX
            PCX = PGX - PLX
            CALL VPL_B( TX,PSWX,PCX,RHOBX,PVWX )
            PVW_BH(2,NBN) = VIC_BH(21,NBH)*PVWX
          ENDIF
          PG_BH(2,NBN) = PGX - PATM
          PL_BH(2,NBN) = PLX - PATM
        ENDDO
      ENDDO
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
!---      Small concentration limits  ---
!
          IF( YLS_BH(2,NBN).LT.1.D-15 ) YLS_BH(2,NBN) = 0.D+0
          IF( XLS_BH(2,NBN).LT.1.D-15 ) XLS_BH(2,NBN) = 0.D+0
          IF( XLA_BH(2,NBN).LT.1.D-15 ) XLA_BH(2,NBN) = 0.D+0
          XLS_BH(1,NBN) = XLS_BH(2,NBN)
!
!---      Coaxial pipe flow  ---
!
          INV = INV_BH(NBN)
          IF( IS_BH(INV).EQ.10000 ) THEN
!
!---        Initialize reference pressure for compressibility  ---
!
            PCMP_BH(NBN) = MAX( PL_BH(2,NBN),PG_BH(2,NBN) )+PATM
            IF( IBN_BH(NBN).NE.0 ) THEN
              PORD_BH(2,NBN) = MIN( 1.D+0,MAX(
     &          ((PAR_BH(3,INV)**2)-(PAR_BH(2,INV)**2))/
     &          ((PAR_BH(4,INV)**2)-(PAR_BH(2,INV)**2)),EPSL ) )
            ELSE
              PORD_BH(2,NBN) = MIN( 1.D+0,MAX( (PAR_BH(1,INV)**2)/
     &          (PAR_BH(2,INV)**2),EPSL ) )
            ENDIF
!
!---        Pipe porosity  ---
!
            DO M = 3,ISVC+2
              PORD_BH(M,NBN) = PORD_BH(2,NBN)
            ENDDO
!
!---      Pipe flow  ---
!
          ELSEIF( IS_BH(INV).GT.10 ) THEN
!
!---        Initialize reference pressure for compressibility  ---
!
            PCMP_BH(NBN) = MAX( PL_BH(2,NBN),PG_BH(2,NBN) )+PATM
            PORD_BH(2,NBN) = MIN( 1.D+0,MAX( (PAR_BH(3,INV)**2)/
     &        (PAR_BH(2,INV)**2),EPSL ) )
!
!---        Pipe porosity  ---
!
            DO M = 3,ISVC+2
              PORD_BH(M,NBN) = PORD_BH(2,NBN)
            ENDDO
!
!---      Porous media flow  ---
!
          ELSE
            IZN = IZ_BH(NBN)
!
!---        Initialize reference pressure for compressibility  ---
!
            IF( CMP(3,IZN).GT.PATM ) THEN
              PCMP_BH(NBN) = CMP(3,IZN)
            ELSEIF( ISLC(61).EQ.0 ) THEN
              PCMP_BH(NBN) = MAX( PL_BH(2,NBN),PG_BH(2,NBN) )+PATM
            ENDIF
!
!---        Porous-media porosity  ---
!
            PX = MAX( PGX,PLX )
            CALL PORSTY_BH( NBN,PX,PCMP_BH(NBN),PORD_BH(2,NBN) )
            PORD_BH(2,NBN) = MAX( PORD_BH(2,NBN),EPSL )
          ENDIF

          POR0_BH(1,NBN) = POR0_BH(1,NBN)
          POR0_BH(2,NBN) = POR0_BH(2,NBN)

        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_FRC_GT
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
!     Check the thermodynamic and hydrologic states declared through
!     user inputs for fracture triangles.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 6 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE GRID
      USE GEOM_FRC
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP
      USE FDVG_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK_FRC_GT'
!
!---  Thermal capacitance of rock surrounding fracture, loop over
!     fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fracture triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        IF( IXP_FRC(NTX).EQ.0 ) THEN
!          PRINT *,' Inactive Fracture Triangle = ',NTX
          CYCLE
        ENDIF
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        AFN_FRCX = 0.D+0
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
          N = INCM_FRC(NCX)
          IZN = IZ(N)
          THCP_FRC(NTX) = THCP_FRC(NTX)+AFN_FRC(NCX)*RHOS(IZN)*CPS(IZN)
          AFN_FRCX = AFN_FRCX + AFN_FRC(NCX)
        ENDDO
        THCP_FRC(NTX) = THCP_FRC(NTX)/AFN_FRCX
      ENDDO
      ENDDO
!
!---  Skip for restart simulations  ---
!
      IF( IEO.EQ.2 .AND. ISLC(87).EQ.0 ) THEN
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
            XLS_FRC(2,NTX) = YLS_FRC(2,NTX)
            XLS_FRC(1,NTX) = XLS_FRC(2,NTX)
          ENDDO
        ENDDO
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
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
!---      Dissolved air relative saturation  ---
!
          PVA_FRC(2,NTX) = VIC_FRC(20,NFX)
!
!---      Hydrostatic initial conditions  ---
!
          IF( IC_FRC(NFX).EQ.4 ) THEN
            CALL HYDST_FRC_GT( VIC_FRC(1,NFX),PG_FRC(2,NTX),
     &        PL_FRC(2,NTX),T_FRC(2,NTX),XLS_FRC(2,NTX),YLS_FRC(2,NTX),
     &        ZP_FRC(NTX) )
!
!---      Aqueous saturation + gas pressure initial conditions  ---
!
          ELSEIF( IC_FRC(NFX).EQ.1 ) THEN
            PG_FRC(2,NTX) = VIC_FRC(4,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(5,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(6,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(7,NFX)
            SL_FRC(2,NTX) = VIC_FRC(8,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(9,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(10,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(11,NFX)
            T_FRC(2,NTX) = VIC_FRC(12,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(13,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(14,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(15,NFX)
            YLS_FRC(2,NTX) = VIC_FRC(16,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(17,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(18,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(19,NFX)
            CALL CAP_FRC_GT( CPGLX,SL_FRC(2,NTX),NFX )
            PL_FRC(2,NTX) = PG_FRC(2,NTX) - CPGLX
!
!---      Aqueous saturation + aqueous pressure initial conditions  ---
!
          ELSEIF( IC_FRC(NFX).EQ.2 ) THEN
            PL_FRC(2,NTX) = VIC_FRC(4,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(5,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(6,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(7,NFX)
            SL_FRC(2,NTX) = VIC_FRC(8,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(9,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(10,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(11,NFX)
            T_FRC(2,NTX) = VIC_FRC(12,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(13,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(14,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(15,NFX)
            YLS_FRC(2,NTX) = VIC_FRC(16,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(17,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(18,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(19,NFX)
            CALL CAP_FRC_GT( CPGLX,SL_FRC(2,NTX),NFX )
            PG_FRC(2,NTX) = PL_FRC(2,NTX) + CPGLX
!
!---      Aqueous pressure + gas pressure initial conditions  ---
!
          ELSEIF( IC_FRC(NFX).EQ.3 ) THEN
            PG_FRC(2,NTX) = VIC_FRC(4,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(5,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(6,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(7,NFX)
            PL_FRC(2,NTX) = VIC_FRC(8,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(9,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(10,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(11,NFX)
            T_FRC(2,NTX) = VIC_FRC(12,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(13,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(14,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(15,NFX)
            YLS_FRC(2,NTX) = VIC_FRC(16,NFX) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(17,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(18,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(19,NFX)
!
!---      Equilibrium initial conditions  ---
!
          ELSEIF( IC_FRC(NFX).EQ.0 ) THEN
!
!---        Loop over connected nodes to determine the total
!           connected distance  ---
!
            DISTX = 0.D+0
            DO K = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
              DISTX = DISTX + MAX( DFN_FRC(K),1.D-9 )
            ENDDO
!
!---        Distance weighted equilibrium values for primary
!           variables  ---
!
            PG_FRC(2,NTX) = 0.D+0
            PL_FRC(2,NTX) = 0.D+0
            T_FRC(2,NTX) = 0.D+0
            YLS_FRC(2,NTX) = 0.D+0
            DO K = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
              N = INCM_FRC(K)
              PG_FRC(2,NTX) = PG_FRC(2,NTX) +
     &          PG(2,N)*MAX( DFN_FRC(K),1.D-9 )/DISTX
              PL_FRC(2,NTX) = PL_FRC(2,NTX) +
     &          PL(2,N)*MAX( DFN_FRC(K),1.D-9 )/DISTX
              T_FRC(2,NTX) = T_FRC(2,NTX) +
     &          T(2,N)*MAX( DFN_FRC(K),1.D-9 )/DISTX
              YLS_FRC(2,NTX) = YLS_FRC(2,NTX) +
     &          YLS(2,N)*MAX( DFN_FRC(K),1.D-9 )/DISTX
            ENDDO
          ENDIF
!
!---      Saturation and relative permeability  ---
!
          CALL RKSP_FRC_GT( PG_FRC(2,NTX),PL_FRC(2,NTX),
     &      RKG_FRC(2,NTX),RKL_FRC(2,NTX),SG_FRC(2,NTX),SL_FRC(2,NTX),
     &      NFX )
!
!---      Check for out of range parameters  ---
!
          IF( T_FRC(2,NTX).GT.2.D+3 .OR. T_FRC(2,NTX).LT.0.01D+0 ) THEN
            INDX = 17
            IMSGX = 0
            NMSGX = 0
            SUBLOGX = 'CHK_FRC_GT'
            RLMSGX = T_FRC(2,NTX)
            CHMSGX = 'Out of Range Initial Fracture Temperature(C) = '
            CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ENDIF
          IF( PL_FRC(2,NTX).GT.100.D+6-PATM ) THEN
            INDX = 17
            IMSGX = 0
            NMSGX = 0
            SUBLOGX = 'CHK_FRC_GT'
            RLMSGX = PL_FRC(2,NTX)+PATM
            CHMSGX = 'Out of Range Initial Fracture Aqu. ' //
     &        'Pressure(Pa) = '
            CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ENDIF
          IF( PG_FRC(2,NTX).GT.100.D+6-PATM .OR.
     &      PG_FRC(2,NTX).LT.6.1125D+2-PATM ) THEN
            IMSGX = 0
            NMSGX = 0
            SUBLOGX = 'CHK_FRC_GT'
            INDX = 17
            RLMSGX = PG_FRC(2,NTX)+PATM
            CHMSGX = 'Out of Range Initial Fracture Gas Pressure(Pa) = '
            CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ENDIF
          IF( SL_FRC(2,NTX).GT.1.D+0 .OR. SL_FRC(2,NTX).LT.0.D+0 ) THEN
            INDX = 17
            IMSGX = 0
            NMSGX = 0
            SUBLOGX = 'CHK_FRC_GT'
            RLMSGX = SL_FRC(2,NTX)
            CHMSGX = 'Out of Range Initial Fracture Aqu. Saturation = '
            CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ENDIF
!
!---      Assign initial phase condition  ---
!
          PLX = PL_FRC(2,NTX)+PATM
          PGX = PG_FRC(2,NTX)+PATM
          TKX = T_FRC(2,NTX)+TABS
          TX = T_FRC(2,NTX)
!
!---      Saturated conditions with no trapped gas  ---
!
          IF( (1.D+0-SL_FRC(2,NTX))/EPSL.LT.EPSL ) THEN
            NPHAZ_FRC(2,NTX) = 1
!
!---      Unsaturated conditions  ---
!
          ELSEIF( SL_FRC(2,NTX)/EPSL.LT.EPSL .OR. TKX.GT.TCRW ) THEN
            NPHAZ_FRC(2,NTX) = 4
!
!---      Partially saturated conditions  ---
!
          ELSE
            NPHAZ_FRC(2,NTX) = 2
          ENDIF
!
!---      Saturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - trapped gas saturation
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ_FRC(2,NTX).EQ.1 ) THEN
!
!---        Check dissolved salt mass fraction  ---
!
            CALL SOL_BRNS( TX,PLX,XLSMX )
            IF( YLS_FRC(2,NTX).GT.XLSMX ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_FRC_GT'
              RLMSGX = XLSMX
              CHMSGX = 'Saturated Fracture Initial Conditions ' //
     &          'Aqu. Salt Mass Frac. > Solubility Limit : ' //
     &          'Aqu. Salt Mass Frac. Solubility Limit = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
            XLS_FRC(2,NTX) = YLS_FRC(2,NTX)
            XLS_FRC(1,NTX) = XLS_FRC(2,NTX)
!
!---        Check pressure exceeds water vapor pressure  ---
!
            CALL SP_B( TX,XLS_FRC(2,NTX),PSBX )
            IF( PLX.LT.PSBX ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_FRC_GT'
              RLMSGX = PSBX
              CHMSG = 'Saturated Fracture Initial Conditions ' //
     &          'Aqueous Pressure < Saturated Pressure : ' //
     &          'Saturated Pressure, Pa = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
!
!---        Dissolved air mole fraction  ---
!
            XMLA_FRC(2,NTX) = VIC_FRC(20,NFX)*PLX/HCAW
!
!---      Unsaturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ_FRC(2,NTX).EQ.2 ) THEN
!
!---        Saturated water vapor pressure  ---
!
            INDX = 0
            CALL REGION_4( TX,PSWX,INDX )
!
!---        Check dissolved salt mass fraction  ---
!
            CALL SOL_BRNS( TX,PSWX,XLSMX )
            IF( YLS_FRC(2,NTX).GT.XLSMX ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_FRC_GT'
              RLMSGX = XLSMX
              CHMSGX = 'Unsaturated Fracture Initial Conditions ' //
     &          'Aqu. Salt Mass Frac. > Solubility Limit : ' //
     &          'Aqu. Salt Mass Frac. Solubility Limit = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
            XLS_FRC(2,NTX) = YLS_FRC(2,NTX)
            XLS_FRC(1,NTX) = XLS_FRC(2,NTX)
!
!---        Check gas pressure exceeds water vapor pressure  ---
!
            CALL P_IAPWS( TX,PSWX,RHOX,RHOLWX,HX,HLWX,UX,ULWX )
            CALL DENS_B( XLS_FRC(2,NTX),RHOBX,RHOLWX,TX )
            CALL SP_B( TX,XLS_FRC(2,NTX),PSBX )
            PCX = PGX - PLX
            CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX )
            IF( PGX.LT.PVBX ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_FRC_GT'
              RLMSGX = PVBX
              CHMSGX = 'Unsaturated Fracture Initial Conditions ' //
     &          'Gas Pressure < Water Vapor Pressure : ' //
     &          'Water Vapor Pressure, Pa = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
!
!---        Dissolved air mole fraction  ---
!
            PVA_FRC(2,NTX) = MAX( PGX-PVBX,0.D+0 )
            XMLA_FRC(2,NTX) = PVA_FRC(2,NTX)/HCAW
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - gas pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ_FRC(2,NTX).EQ.4 ) THEN
!
!---        Check for non-zero dissolved salt concentration  ---
!
            IF( YLS_FRC(2,NTX).GT.EPSL ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_FRC_GT'
              RLMSGX = YLS_FRC(2,NTX)
              CHMSGX = 'Fully Unsaturated Fracture Initial ' //
     &          'Conditions Non-Zero Dissolved Salt Mass Fraction ' //
     &          'Declared : Salt Concentration = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
            XLS_FRC(2,NTX) = 0.D+0
            YLS_FRC(2,NTX) = 0.D+0
            XLS_FRC(1,NTX) = 0.D+0
!
!---        Water vapor pressure  ---
!
            INDX = 0
            CALL REGION_4( TX,PSWX,INDX )
            SL_FRC(2,NTX) = 0.D+0
            CALL CAP_FRC_GT( CPGLX,SL_FRC(2,NTX),NFX )
            PLX = PGX - CPGLX
            PCX = PGX - PLX
            CALL VPL_B( TX,PSWX,PCX,RHOBX,PVWX )
            PVW_FRC(2,NTX) = VIC_FRC(21,NFX)*PVWX
          ENDIF
          PG_FRC(2,NTX) = PGX - PATM
          PL_FRC(2,NTX) = PLX - PATM

          POR0_FRC(1,NTX) = POR0_FRC(1,NTX)
          POR0_FRC(2,NTX) = POR0_FRC(2,NTX)

        ENDDO
      ENDDO
!!
!!---  Loop over fractures  ---
!!
!      DPX = 0.D+0
!      DO NFX = 1,NF_FRC
!!
!!---    Loop over fracture triangles  ---
!!
!        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!          IF( IXP_FRC(NTX).EQ.0 ) PRINT *,'Inactive Triangle: ',NTX
!!
!!---      Skip inactive triangles  ---
!!
!          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!          IF( PL_FRC(2,NTX).LT.2.D+6 .AND. 
!     &      (IPF_FRC(2,NTX)-IPF_FRC(1,NTX)).LE.1 ) THEN
!            PRINT *,'Drift Triangle: ',NTX
!            DPX = DPX + 1.D+6
!            PL_FRC(2,NTX) = 2.D+7 + DPX
!            PG_FRC(2,NTX) = 2.D+7 + DPX
!          ENDIF
!        ENDDO
!      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HYDST_BH_GT( VICX,PGX,PLX,TX,XLSX,YLSX,ZPX )
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
!     STOMP-GT
!
!     Establish hydrostatic fracture initial conditions for
!     boreholes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 21 April 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FDVS
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VICX(21)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HYDST_BH_GT'
!
!---  Hydrostatic fracture conditions, establish aqueous pressure,
!     temperature, and salt mass fraction at fracture centroid  ---
!
      ZX = VICX(2)
      KC = MAX( 1,INT(ABS(ZX-ZPX)) )
      DISTZ = (ZX-ZPX)/REAL(KC)
      PLX = VICX(1) + PATM
      PGX = PLX
      XLAX = 0.D+0
      DO 100 K = 1,KC
        ZX = ZX - 5.D-1*DISTZ
        TX = VICX(3) + (ZX-VICX(4))*VICX(5)
        YLSX = VICX(6) + (ZX-VICX(7))*VICX(8)
        ZX = ZX - 5.D-1*DISTZ
!
!---    Check for out-of-range salt mass fraction  ---
!
        CALL SOL_BRNS( TX,PLX,XLSMX )
        IF( YLSX.LT.0.D+0 ) THEN
          INDX = 14
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'HYDST_BH_GT'
          RLMSGX = YLSX
          CHMSGX = 'Hydrostatic Borehole Condition: Salt ' //
     &      'Mass Fraction < 0.0 : XLS = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ELSEIF( YLSX.GT.XLSMX ) THEN
          INDX = 17
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'HYDST_BH_GT'
          RLMSGX = YLSX
          CHMSGX = 'Hydrostatic Borehole Condition: Salt ' //
     &      'Mass Fraction > Solubility Limit : XLS = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
        CALL P_IAPWS( TX,PLX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        XLSX = YLSX
        CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
        RHOLX = RHOBX
        PLX = PLX + RHOLX*GRAVZ*DISTZ
        PGX = PLX
  100 CONTINUE
      PLX = PLX - PATM
      PGX = PLX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HYDST_BH_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HYDST_FRC_GT( VICX,PGX,PLX,TX,XLSX,YLSX,ZPX )
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
!     STOMP-GT
!
!     Establish hydrostatic fracture initial conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 6 September 2017
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FDVS
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VICX(21)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HYDST_FRC_GT'
!
!---  Hydrostatic fracture conditions, establish aqueous pressure,
!     temperature, and salt mass fraction at fracture centroid  ---
!
      ZX = VICX(2)
      KC = MAX( 1,INT(ABS(ZX-ZPX)) )
      DISTZ = (ZX-ZPX)/REAL(KC)
      PLX = VICX(1) + PATM
      PGX = PLX
      XLAX = 0.D+0
      DO 100 K = 1,KC
        ZX = ZX - 5.D-1*DISTZ
        TX = VICX(3) + (ZX-VICX(4))*VICX(5)
        YLSX = VICX(6) + (ZX-VICX(7))*VICX(8)
        ZX = ZX - 5.D-1*DISTZ
!
!---    Check for out-of-range salt mass fraction  ---
!
        CALL SOL_BRNS( TX,PLX,XLSMX )
        IF( YLSX.LT.0.D+0 ) THEN
          INDX = 14
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'HYDST_FRC_GT'
          RLMSGX = YLSX
          CHMSGX = 'Hydrostatic Fracture Condition: Salt ' //
     &      'Mass Fraction < 0.0 : XLS = '
           CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ELSEIF( YLSX.GT.XLSMX ) THEN
          INDX = 17
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'HYDST_FRC_GT'
          RLMSGX = YLSX
          CHMSGX = 'Hydrostatic Fracture Condition: Salt ' //
     &      'Mass Fraction > Solubility Limit : XLS = '
           CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
        CALL P_IAPWS( TX,PLX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        XLSX = YLSX
        CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
        RHOLX = RHOBX
        PLX = PLX + RHOLX*GRAVZ*DISTZ
        PGX = PLX
  100 CONTINUE
      PLX = PLX - PATM
      PGX = PLX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HYDST_FRC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGW_BF_GT
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
!     Compute the water-vapor molar diffusion rates for boreholes
!     and fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 13 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
      USE FDVG_FRC
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
      SUB_LOG(ISUB_LOG) = '/DFFGW_BF_GT'
!
!---  Molar-density diffusion gradient option ---
!
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
!---      Loop over borehole node to node connections ---
!
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBN2 = IBCM_BH(NCX)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
            DXMGW = XMGW_BH(2,NBN1)*RHOMG_BH(2,NBN1)-
     &        XMGW_BH(2,NBN2)*RHOMG_BH(2,NBN2)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DF1 = SG_BH(MP,NBN1)*DFGW_BH(MP,NBN1)
              DF2 = SG_BH(MN,NBN2)*DFGW_BH(MN,NBN2)
              INDX = 12
              DFM = DIFMN( DF1,DF2,DBB1X,DBB2X,DXMGW,INDX )
              UBBDGW(M,NCX) = (XMGW_BH(MP,NBN1)*RHOMG_BH(MP,NBN1)-
     &          XMGW_BH(MN,NBN2)*RHOMG_BH(MN,NBN2))*DFM/DBB_BH(NCX)
              FGW1 = XGW_BH(MP,NBN1)*RHOG_BH(MP,NBN1)
              FGW2 = XGW_BH(MN,NBN2)*RHOG_BH(MN,NBN2)
              INDX = 3
              FGW = DIFMN( FGW1,FGW2,DBB1X,DBB2X,UBBG(1,NCX),INDX )
              UBBGW(M,NCX) = UBBG(M,NCX)*FGW +
     &          WTMW*UBBDGW(M,NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Molar-density diffusion gradient option ---
!
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
!---      Loop over fracture triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
            DXMGW = XMGW_FRC(2,NT1X)*RHOMG_FRC(2,NT1X)-
     &        XMGW_FRC(2,NT2X)*RHOMG_FRC(2,NT2X)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DF1 = SG_FRC(MP,NT1X)*DFGW_FRC(MP,NT1X)
              DF2 = SG_FRC(MN,NT2X)*DFGW_FRC(MN,NT2X)
              INDX = 12
              DFM = DIFMN( DF1,DF2,DFF1X,DFF2X,DXMGW,INDX )
              UFFDGW(M,NCX) = (XMGW_FRC(MP,NT1X)*RHOMG_FRC(MP,NT1X)-
     &          XMGW_FRC(MN,NT2X)*RHOMG_FRC(MN,NT2X))*DFM/DFF_FRC(NCX)
              FGW1 = XGW_FRC(MP,NT1X)*RHOG_FRC(MP,NT1X)
              FGW2 = XGW_FRC(MN,NT2X)*RHOG_FRC(MN,NT2X)
              INDX = 3
              FGW = DIFMN( FGW1,FGW2,DFF1X,DFF2X,UFFG(1,NCX),INDX )
              UFFGW(M,NCX) = UFFG(M,NCX)*FGW +
     &          WTMW*UFFDGW(M,NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGW_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLA_BF_GT
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
!     Compute the dissolved air molar diffusion rates through the
!     aqueous phase for boreholes and fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 13 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVP_FRC
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
      SUB_LOG(ISUB_LOG) = '/DFFLA_BF_GT'
!
!---  Molar-density diffusion gradient option ---
!
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
!---      Loop over borehole node to node connections ---
!
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBN2 = IBCM_BH(NCX)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
            DXMLA = XMLA_BH(2,NBN1)*RHOML_BH(2,NBN1)-
     &        XMLA_BH(2,NBN2)*RHOML_BH(2,NBN2)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DF1 = SL_BH(MP,NBN1)*DFLA_BH(MP,NBN1)
              DF2 = SL_BH(MN,NBN2)*DFLA_BH(MN,NBN2)
              INDX = 14
              DFM = DIFMN( DF1,DF2,DBB1X,DBB2X,DXMLA,INDX )
              UBBDLA(M,NCX) = (XMLA_BH(MP,NBN1)*RHOML_BH(MP,NBN1)-
     &          XMLA_BH(MN,NBN2)*RHOML_BH(MN,NBN2))*DFM/DBB_BH(NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Molar-density diffusion gradient option ---
!
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
!---      Loop over fracture triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
            DXMLA = XMLA_FRC(2,NT1X)*RHOML_FRC(2,NT1X)-
     &        XMLA_FRC(2,NT2X)*RHOML_FRC(2,NT2X)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DF1 = SL_FRC(MP,NT1X)*DFLA_FRC(MP,NT1X)
              DF2 = SL_FRC(MN,NT2X)*DFLA_FRC(MN,NT2X)
              INDX = 14
              DFM = DIFMN( DF1,DF2,DFF1X,DFF2X,DXMLA,INDX )
              UFFDLA(M,NCX) = (XMLA_FRC(MP,NT1X)*RHOML_FRC(MP,NT1X)-
     &          XMLA_FRC(MN,NT2X)*RHOML_FRC(MN,NT2X))*DFM/DFF_FRC(NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLA_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLS_BF_GT
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
!     Compute salt aqueous-phase fluxes for boreholes and fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 13 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVS_FRC
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
      SUB_LOG(ISUB_LOG) = '/DFFLS_BF_GT'
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
!---      Loop over borehole node to node connections ---
!
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBN2 = IBCM_BH(NCX)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
!
!---          Diffusion coefficients  ---
!
              IF( IEDLS.EQ.1 ) THEN
                TCOR = (T_BH(MP,NBN1)+TABS)/TSPRF
                SMDL1 = DFLS_BH(MP,NBN1)*TCOR*(VISRL/VISL_BH(MP,NBN1))
                DFC1 = SL_BH(MP,NBN1)*SMDL1
                TCOR = (T_BH(MN,NBN2)+TABS)/TSPRF
                SMDL2 = DFLS_BH(MN,NBN2)*TCOR*(VISRL/VISL_BH(MN,NBN2))
                DFC2 = SL_BH(MN,NBN2)*SMDL2
              ELSEIF( IEDLS.EQ.3 ) THEN
                DFC1 = SL_BH(MP,NBN1)*DFLS_BH(MP,NBN1)
                DFC2 = SL_BH(MN,NBN2)*DFLS_BH(MN,NBN2)
              ENDIF
              INDX = 18
              DFC = DIFMN(DFC1,DFC2,DBB1X,DBB2X,UBBL(1,NCX),INDX)
              DDL = DFC/DBB_BH(NCX)
              AL = MAX( UBBL(M,NCX),ZERO ) +
     &          DDL*MAX((ONE-(TENTH*ABS(UBBL(M,NCX))/
     &          (DDL+SMALL)))**5,ZERO)
              ALP = MAX( -UBBL(M,NCX),ZERO ) +
     &          DDL*MAX((ONE-(TENTH*ABS(UBBL(M,NCX))/
     &          (DDL+SMALL)))**5,ZERO)
              UBBS(M,NCX) = XLS_BH(MP,NBN1)*RHOL_BH(MP,NBN1)*AL -
     &          XLS_BH(MN,NBN2)*RHOL_BH(MN,NBN2)*ALP
              UBBDS(M,NCX) = DDL*(XLS_BH(MP,NBN1)*RHOL_BH(MP,NBN1)
     &          - XLS_BH(MN,NBN2)*RHOL_BH(MN,NBN2))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
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
!---      Loop over fracture triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
!
!---          Diffusion coefficients  ---
!
              IF( IEDLS.EQ.1 ) THEN
                TCOR = (T_FRC(MP,NT1X)+TABS)/TSPRF
                SMDL1 = DFLS_FRC(MP,NT1X)*TCOR*(VISRL/VISL_FRC(MP,NT1X))
                DFC1 = SL_FRC(MP,NT1X)*SMDL1
                TCOR = (T_FRC(MN,NT2X)+TABS)/TSPRF
                SMDL2 = DFLS_FRC(MN,NT2X)*TCOR*(VISRL/VISL_FRC(MN,NT2X))
                DFC2 = SL_FRC(MN,NT2X)*SMDL2
              ELSEIF( IEDLS.EQ.3 ) THEN
                DFC1 = SL_FRC(MP,NT1X)*DFLS_FRC(MP,NT1X)
                DFC2 = SL_FRC(MN,NT2X)*DFLS_FRC(MN,NT2X)
              ENDIF
              INDX = 18
              DFC = DIFMN(DFC1,DFC2,DFF1X,DFF2X,UFFL(1,NCX),INDX)
              DDL = DFC/DFF_FRC(NCX)
              AL = MAX( UFFL(M,NCX),ZERO ) +
     &          DDL*MAX((ONE-(TENTH*ABS(UFFL(M,NCX))/
     &          (DDL+SMALL)))**5,ZERO)
              ALP = MAX( -UFFL(M,NCX),ZERO ) +
     &          DDL*MAX((ONE-(TENTH*ABS(UFFL(M,NCX))/
     &          (DDL+SMALL)))**5,ZERO)
              UFFS(M,NCX) = XLS_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)*AL -
     &          XLS_FRC(MN,NT2X)*RHOL_FRC(MN,NT2X)*ALP
              UFFDS(M,NCX) = DDL*(XLS_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)
     &          - XLS_FRC(MN,NT2X)*RHOL_FRC(MN,NT2X))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLS_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLS_BTF_GT
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
!     Compute the water-vapor molar diffusion rates, and the dissolved
!     air molar diffusion rates through the aqueous phase, from
!     borehole nodes to fracture triangles (borehole node is considered
!     the upper node for flux indexing)
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 7 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
      USE FDVS_FRC
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
      SUB_LOG(ISUB_LOG) = '/DFFLS_BTF_GT'
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
          IZN = IZ_BH(NBN)
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
!---      Loop over flux increments ---
!
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Diffusion coefficients  ---
!
            IF( IEDLS.EQ.1 ) THEN
              TCOR = (T_BH(MP,NBN)+TABS)/TSPRF
              SMDL1 = DFLS_BH(MP,NBN)*TCOR*(VISRL/VISL_BH(MP,NBN))
              DFC1 = SL_BH(MP,NBN)*SMDL1
              TCOR = (T_FRC(MN,NTX)+TABS)/TSPRF
              SMDL2 = DFLS_FRC(MN,NTX)*TCOR*(VISRL/VISL_FRC(MN,NTX))
              DFC2 = SL_FRC(MN,NTX)*SMDL2
            ELSEIF( IEDLS.EQ.3 ) THEN
              DFC1 = SL_BH(MP,NBN)*DFLS_BH(MP,NBN)
              DFC2 = SL_FRC(MN,NTX)*DFLS_FRC(MN,NTX)
            ENDIF
            INDX = 18
            DFC = DIFMN(DFC1,DFC2,DBF1X,DBF2X,UBFL(1,NCX),INDX)
            DDL = DFC/DBFX
            AL = MAX( UBFL(M,NCX),ZERO ) +
     &        DDL*MAX((ONE-(TENTH*ABS(UBFL(M,NCX))/
     &        (DDL+SMALL)))**5,ZERO)
            ALP = MAX( -UBFL(M,NCX),ZERO ) +
     &        DDL*MAX((ONE-(TENTH*ABS(UBFL(M,NCX))/
     &        (DDL+SMALL)))**5,ZERO)
            UBFS(M,NCX) = XLS_BH(MP,NBN)*RHOL_BH(MP,NBN)*AL -
     &        XLS_FRC(MN,NTX)*RHOL_FRC(MN,NTX)*ALP
            UBFDS(M,NCX) = DDL*(XLS_BH(MP,NBN)*RHOL_BH(MP,NBN)
     &        - XLS_FRC(MN,NTX)*RHOL_FRC(MN,NTX))
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLS_BTF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFX_BTF_GT
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
!     Compute the water-vapor molar diffusion rates, and the dissolved
!     air molar diffusion rates through the aqueous phase, from
!     borehole nodes to fracture triangles (borehole node is considered
!     the upper node for flux indexing)
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 7 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
      USE FDVG_FRC
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
      SUB_LOG(ISUB_LOG) = '/DFFX_BTF_GT'
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
          IZN = IZ_BH(NBN)
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
!---      Unincremented water-vapor diffusive gradient  ---
!
          DXMGW = XMGW_BH(2,NBN)*RHOMG_BH(2,NBN)-
     &      XMGW_FRC(2,NTX)*RHOMG_FRC(2,NTX)
!
!---      Unincremented dissolved-air diffusive gradient  ---
!
          DXMLA = XMLA_BH(2,NBN)*RHOML_BH(2,NBN)-
     &      XMLA_FRC(2,NTX)*RHOML_FRC(2,NTX)
!
!---      Loop over flux increments ---
!
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Water-vapor diffusive flux  ---
!
            DF1 = SG_BH(MP,NBN)*DFGW_BH(MP,NBN)
            DF2 = SG_FRC(MN,NTX)*DFGW_FRC(MN,NTX)
            INDX = 12
            DFM = DIFMN( DF1,DF2,DBF1X,DBF2X,DXMGW,INDX )
            UBFDGW(M,NCX) = (XMGW_BH(MP,NBN)*RHOMG_BH(MP,NBN)-
     &        XMGW_FRC(MN,NTX)*RHOMG_FRC(MN,NTX))*DFM/DBFX
            FGW1 = XGW_BH(MP,NBN)*RHOG_BH(MP,NBN)
            FGW2 = XGW_FRC(MN,NTX)*RHOG_FRC(MN,NTX)
            INDX = 3
            FGW = DIFMN( FGW1,FGW2,DBF1X,DBF2X,UBFG(1,NCX),INDX )
            UBFGW(M,NCX) = UBFG(M,NCX)*FGW +
     &        WTMW*UBFDGW(M,NCX)
!
!---        Dissolved-air diffusive flux  ---
!
            DF1 = SL_BH(MP,NBN)*DFLA_BH(MP,NBN)
            DF2 = SL_FRC(MN,NTX)*DFLA_FRC(MN,NTX)
            INDX = 14
            DFM = DIFMN( DF1,DF2,DBF1X,DBF2X,DXMLA,INDX )
            UBFDLA(M,NCX) = (XMLA_BH(MP,NBN)*RHOML_BH(MP,NBN)-
     &        XMLA_FRC(MN,NTX)*RHOML_FRC(MN,NTX))*DFM/DBFX
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFX_BTF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVG_BF_GT
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces for fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 13 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
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
      REAL*8 KGM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVG_BF_GT'
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
          INV1 = INV_BH(NBN1)
!
!---      Loop over borehole node to node connections ---
!
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBN2 = IBCM_BH(NCX)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
            ZP_BH1X = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
            ZP_BH2X = 5.D-1*(ZP_BH(1,NBN2)+ZP_BH(2,NBN2))
            INV2 = INV_BH(NBN2)
!
!---        Coaxial borehole  ---
!
            IF( IT_BH(1,NBH).GE.10000 ) THEN
!
!---          Outer pipe, hydraulic diameter (m) and roughness (m)  ---
!
              IF( IBN_BH(NBN1).NE.0 ) THEN
                DH1X = 2.D+0*SQRT(PAR_BH(3,INV1)**2 - PAR_BH(2,INV1)**2)
                RF1X = PAR_BH(16,INV1)
!
!---          Inner pipe, hydraulic diameter (m) and roughness (m)  ---
!
              ELSE
                DH1X = 2.D+0*PAR_BH(1,INV1)
                RF1X = PAR_BH(13,INV1)
              ENDIF
!
!---          Outer pipe, hydraulic diameter (m) and roughness (m)  ---
!
              IF( IBN_BH(NBN2).NE.0 ) THEN
                DH2X = 2.D+0*SQRT(PAR_BH(3,INV2)**2 - PAR_BH(2,INV2)**2)
                RF2X = PAR_BH(16,INV2)
!
!---          Inner pipe, hydraulic diameter (m) and roughness (m)  ---
!
              ELSE
                DH2X = 2.D+0*PAR_BH(1,INV2)
                RF2X = PAR_BH(13,INV2)
              ENDIF
              INDX = -3
              DHX = DIFMN( DH1X,DH2X,DBB1X,DBB2X,ZERO,INDX )
              RFX = DIFMN( RF1X,RF2X,DBB1X,DBB2X,ZERO,INDX )
!
!---          Loop over flux increments ---
!
              DO M = 1,ISVF
                MN = MNEG(M)
                MP = MPOS(M)
                HDGX = PG_BH(MP,NBN1) - PG_BH(MN,NBN2) +
     &            5.D-1*GRAV*(ZP_BH1X-ZP_BH2X)*
     &           (RHOG_BH(MP,NBN1)+RHOG_BH(MN,NBN2))
                IF( M.EQ.1 ) HDG = HDGX
                INDX = 5
                VGM = DIFMN(VISG_BH(MP,NBN1),VISG_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDG,INDX)
                INDX = 8
                RKGM = DIFMN(RKG_BH(MP,NBN1),RKG_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDG,INDX)
                INDX = 8
                RHOGM = DIFMN(RHOG_BH(MP,NBN1),RHOG_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDG,INDX)
                IF( RKGM.GT.1.D-9 ) THEN
!
!---              Pipe flow velocity (m/s), first estimating from
!                 Hazen-Williams equation, adjusting for fluid 
!                 viscosity, then solving a nonlinear system of  
!                 equations, combining the Colebrook equation for the  
!                 friction coefficient and the Darcy-Weisbach equation 
!                 for the pressure drop (i.e., head loss)  ---
!
                  CALL PIPE_FLOW_VEL_GT( DBB_BH(NCX),DHX,HDGX,RFX,RHOGM,
     &              UGX,VGM,NBN1 )
                  IF( VGM.LT.0.D+0 ) THEN
                    print *,'pg_bh(mp,nbn1) = ',PG_BH(MP,NBN1)
                    print *,'pg_bh(mp,nbn2) = ',PG_BH(MP,NBN2)
                    print *,'rhog_bh(mp,nbn1) = ',RHOG_BH(MP,NBN1)
                    print *,'rhog_bh(mp,nbn2) = ',RHOG_BH(MP,NBN2)
                  ENDIF
                  UBBG(M,NCX) = RKGM*SIGN(UGX,HDGX)
                ELSE
                  UBBG(M,NCX) = 0.D+0
                ENDIF
              ENDDO
!
!---        Pipe flow when both intervals are screened or closed
!           pipes ---
!
            ELSEIF( IS_BH(INV1).GT.10 .AND. IS_BH(INV2).GT.10 ) THEN
!
!---          Hydraulic diameter (m), sum of interval radii ---
!
              DHX = (PAR_BH(3,INV1)+PAR_BH(3,INV2))
!
!---          Roughness, m  ---
!
              RFX = 5.D-1*(PAR_BH(6,INV1)+PAR_BH(6,INV2))
!
!---          Loop over flux increments ---
!
              DO M = 1,ISVF
                MN = MNEG(M)
                MP = MPOS(M)
                HDGX = PG_BH(MP,NBN1) - PG_BH(MN,NBN2) +
     &            5.D-1*GRAV*(ZP_BH1X-ZP_BH2X)*
     &           (RHOG_BH(MP,NBN1)+RHOG_BH(MN,NBN2))
                IF( M.EQ.1 ) HDG = HDGX
                INDX = 5
                VGM = DIFMN(VISG_BH(MP,NBN1),VISG_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDG,INDX)
                INDX = 8
                RKGM = DIFMN(RKG_BH(MP,NBN1),RKG_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDG,INDX)
                INDX = 8
                RHOGM = DIFMN(RHOG_BH(MP,NBN1),RHOG_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDG,INDX)
                IF( RKGM.GT.1.D-9 ) THEN
!
!---              Pipe flow velocity (m/s), first estimating from
!                 Hazen-Williamsequation, adjusting for fluid viscosity,
!                 then solving a nonlinear system of equations, 
!                 combining the Colebrook equation for the friction 
!                 coefficient and the Darcy-Weisbach equation for the
!                 pressure drop (i.e., head loss)  ---
!
                  IF( RHOGM.NE.RHOGM ) THEN
                    print *,'DRCVG_BF_GT'
                    print *,'NBN1 = ',NBN1
                    print *,'NBN2 = ',NBN2
                    print *,'PG_BH(MP,NBN1) = ',PG_BH(MP,NBN1)
                    print *,'PG_BH(MN,NBN2) = ',PG_BH(MN,NBN2)
                    print *,'RHOG_BH(MP,NBN1) = ',RHOG_BH(MP,NBN1)
                    print *,'RHOG_BH(MN,NBN2) = ',RHOG_BH(MN,NBN2)
                    print *,'RKG_BH(MP,NBN1) = ',RKG_BH(MP,NBN1)
                    print *,'RKG_BH(MN,NBN2) = ',RKG_BH(MN,NBN2)
                    stop
                  ENDIF
                  CALL PIPE_FLOW_VEL_GT( DBB_BH(NCX),DHX,HDGX,RFX,RHOGM,
     &              UGX,VGM,NBN1 )
                  IF( VGM.LT.0.D+0 ) THEN
                    print *,'pg_bh(mp,nbn1) = ',PG_BH(MP,NBN1)
                    print *,'pg_bh(mp,nbn2) = ',PG_BH(MP,NBN2)
                    print *,'rhog_bh(mp,nbn1) = ',RHOG_BH(MP,NBN1)
                    print *,'rhog_bh(mp,nbn2) = ',RHOG_BH(MP,NBN2)
                  ENDIF
                  UBBG(M,NCX) = RKGM*SIGN(UGX,HDGX)
                ELSE
                  UBBG(M,NCX) = 0.D+0
                ENDIF
              ENDDO
!
!---        Darcy flow when either interval is filled ---
!
            ELSE
              IZ1X = IZ_BH(NBN1)
              IZ2X = IZ_BH(NBN2)
!
!---          Loop over flux increments ---
!
              DO M = 1,ISVF
                MN = MNEG(M)
                MP = MPOS(M)
                HDGX = PG_BH(MP,NBN1) - PG_BH(MN,NBN2) +
     &            5.D-1*GRAV*(ZP_BH1X-ZP_BH2X)*
     &           (RHOG_BH(MP,NBN1)+RHOG_BH(MN,NBN2))
                IF( M.EQ.1 ) HDG = HDGX
                INDX = 11
                KGM = DIFMN(PERM(1,IZ1X),PERM(1,IZ2X),
     &            DBB1X,DBB2X,HDG,INDX)
                IF( PERM(1,IZ1X)/EPSL.LT.EPSL ) KGM = 0.D+0
                IF( PERM(1,IZ2X)/EPSL.LT.EPSL ) KGM = 0.D+0
                INDX = 9
                RKGM = DIFMN(RKG_BH(MP,NBN1),RKG_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDG,INDX)
                INDX = 6
                VGM = DIFMN(VISG_BH(MP,NBN1),VISG_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDG,INDX)
                UBBG(M,NCX) = KGM*RKGM*HDGX/DBB_BH(NCX)/VGM
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
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
!---      Loop over fracture triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Apply fracture-fracture skin factor for inter-fracture
!           connections  ---
!
            SKFX = 0.D+0
            DO NF2X = 1,NF_FRC
              IF( NT2X.GE.IP_FRC(1,NF2X) .AND.
     &          NT2X.LE.IP_FRC(2,NF2X) ) THEN
                EXIT
              ENDIF
            ENDDO
            IF( NFX.NE.NF2X ) THEN
              SKFX = MAX( SKF_FRC(NFX),SKF_FRC(NF2X) )
            ENDIF
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              HDGX = PG_FRC(MP,NT1X) - PG_FRC(MN,NT2X) +
     &          5.D-1*GRAV*(ZP_FRC(NT1X)-ZP_FRC(NT2X))*
     &         (RHOG_FRC(MP,NT1X)+RHOG_FRC(MN,NT2X))
              IF( M.EQ.1 ) HDG = HDGX
              INDX = 11
              KGM = DIFMN(PERM_FRC(MP,NT1X),PERM_FRC(MN,NT2X),
     &          DFF1X,DFF2X,HDG,INDX)
              IF( PERM_FRC(MP,NT1X)/EPSL.LT.EPSL ) KGM = 0.D+0
              IF( PERM_FRC(MN,NT2X)/EPSL.LT.EPSL ) KGM = 0.D+0
              INDX = 9
              RKGM = DIFMN(RKG_FRC(MP,NT1X),RKG_FRC(MN,NT2X),
     &          DFF1X,DFF2X,HDG,INDX)
              INDX = 6
              VGM = DIFMN(VISG_FRC(MP,NT1X),VISG_FRC(MN,NT2X),
     &          DFF1X,DFF2X,HDG,INDX)
              UFFG(M,NCX) = KGM*RKGM*HDGX/DFF_FRC(NCX)/VGM/(1.D+0+SKFX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVG_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVL_BF_GT
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces for fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 8 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVL_BF_GT'
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
          INV1 = INV_BH(NBN1)
!
!---      Loop over borehole node to node connections ---
!
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBN2 = IBCM_BH(NCX)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
            ZP_BH1X = 5.D-1*(ZP_BH(1,NBN1)+ZP_BH(2,NBN1))
            ZP_BH2X = 5.D-1*(ZP_BH(1,NBN2)+ZP_BH(2,NBN2))
            INV2 = INV_BH(NBN2)
!
!---        Coaxial borehole  ---
!
            IF( IT_BH(1,NBH).GE.10000 ) THEN
!
!---          Outer pipe, hydraulic diameter (m) and roughness (m)  ---
!
              IF( IBN_BH(NBN1).NE.0 ) THEN
                DH1X = 2.D+0*SQRT(PAR_BH(3,INV1)**2 - PAR_BH(2,INV1)**2)
                RF1X = PAR_BH(16,INV1)
!
!---          Inner pipe, hydraulic diameter (m) and roughness (m)  ---
!
              ELSE
                DH1X = 2.D+0*PAR_BH(1,INV1)
                RF1X = PAR_BH(13,INV1)
              ENDIF
!
!---          Outer pipe, hydraulic diameter (m) and roughness (m)  ---
!
              IF( IBN_BH(NBN2).NE.0 ) THEN
                DH2X = 2.D+0*SQRT(PAR_BH(3,INV2)**2 - PAR_BH(2,INV2)**2)
                RF2X = PAR_BH(16,INV2)
!
!---          Inner pipe, hydraulic diameter (m) and roughness (m)  ---
!
              ELSE
                DH2X = 2.D+0*PAR_BH(1,INV2)
                RF2X = PAR_BH(13,INV2)
              ENDIF
              INDX = -3
              DHX = DIFMN( DH1X,DH2X,DBB1X,DBB2X,ZERO,INDX )
              RFX = DIFMN( RF1X,RF2X,DBB1X,DBB2X,ZERO,INDX )
!
!---          Loop over flux increments ---
!
              DO M = 1,ISVF
                MN = MNEG(M)
                MP = MPOS(M)
                HDLX = PL_BH(MP,NBN1) - PL_BH(MN,NBN2) +
     &            5.D-1*GRAV*(ZP_BH1X-ZP_BH2X)*
     &           (RHOL_BH(MP,NBN1)+RHOL_BH(MN,NBN2))
                IF( M.EQ.1 ) HDL = HDLX
                INDX = 5
                VLM = DIFMN(VISL_BH(MP,NBN1),VISL_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDL,INDX)
                INDX = 8
                RKLM = DIFMN(RKL_BH(MP,NBN1),RKL_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDL,INDX)
                INDX = 8
                RHOLM = DIFMN(RHOL_BH(MP,NBN1),RHOL_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDL,INDX)
!
!---            Pipe flow velocity (m/s), first estimating from
!               Hazen-Williamsequation, adjusting for fluid viscosity,
!               then solving a nonlinear system of equations, combining
!               the Colebrook equation for the friction coefficient and
!               the Darcy-Weisbach equation for the
!               pressure drop (i.e., head loss)  ---
!
                CALL PIPE_FLOW_VEL_GT( DBB_BH(NCX),DHX,HDLX,RFX,RHOLM,
     &            ULX,VLM,NBN1 )
                IF( VLM.LT.0.D+0 ) THEN
                  print *,'pl_bh(mp,nbn1) = ',PL_BH(MP,NBN1)
                  print *,'pl_bh(mp,nbn2) = ',PL_BH(MP,NBN2)
                  print *,'rhol_bh(mp,nbn1) = ',RHOL_BH(MP,NBN1)
                  print *,'rhol_bh(mp,nbn2) = ',RHOL_BH(MP,NBN2)
                ENDIF
                UBBL(M,NCX) = RKLM*SIGN(ULX,HDLX)
              ENDDO
!
!---        Pipe flow when both intervals are screened or closed
!           pipes ---
!
            ELSEIF( IS_BH(INV1).GT.10 .AND. IS_BH(INV2).GT.10 ) THEN
!
!---          Hydraulic diameter (m), sum of interval radii ---
!
              DHX = (PAR_BH(3,INV1)+PAR_BH(3,INV2))
!
!---          Roughness, m  ---
!
              RFX = 5.D-1*(PAR_BH(6,INV1)+PAR_BH(6,INV2))
!
!---          Loop over flux increments ---
!
              DO M = 1,ISVF
                MN = MNEG(M)
                MP = MPOS(M)
                HDLX = PL_BH(MP,NBN1) - PL_BH(MN,NBN2) +
     &            5.D-1*GRAV*(ZP_BH1X-ZP_BH2X)*
     &           (RHOL_BH(MP,NBN1)+RHOL_BH(MN,NBN2))
                IF( M.EQ.1 ) HDL = HDLX
                INDX = 5
                VLM = DIFMN(VISL_BH(MP,NBN1),VISL_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDL,INDX)
                INDX = 8
                RKLM = DIFMN(RKL_BH(MP,NBN1),RKL_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDL,INDX)
                INDX = 8
                RHOLM = DIFMN(RHOL_BH(MP,NBN1),RHOL_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDL,INDX)
!
!---            Pipe flow velocity (m/s), first estimating from
!               Hazen-Williamsequation, adjusting for fluid viscosity,
!               then solving a nonlinear system of equations, combining
!               the Colebrook equation for the friction coefficient and
!               the Darcy-Weisbach equation for the
!               pressure drop (i.e., head loss)  ---
!
                IF( RHOLM.NE.RHOLM ) THEN
                  print *,'DRCVL_BF_GT'
                  print *,'NBN1 = ',NBN1
                  print *,'NBN2 = ',NBN2
                  print *,'PL_BH(MP,NBN1) = ',PL_BH(MP,NBN1)
                  print *,'PL_BH(MN,NBN2) = ',PL_BH(MN,NBN2)
                  print *,'RHOL_BH(MP,NBN1) = ',RHOL_BH(MP,NBN1)
                  print *,'RHOL_BH(MN,NBN2) = ',RHOL_BH(MN,NBN2)
                  print *,'RKL_BH(MP,NBN1) = ',RKL_BH(MP,NBN1)
                  print *,'RKL_BH(MN,NBN2) = ',RKL_BH(MN,NBN2)
                  stop
                ENDIF
                CALL PIPE_FLOW_VEL_GT( DBB_BH(NCX),DHX,HDLX,RFX,RHOLM,
     &            ULX,VLM,NBN1 )
                IF( VLM.LT.0.D+0 ) THEN
                  print *,'pl_bh(mp,nbn1) = ',PL_BH(MP,NBN1)
                  print *,'pl_bh(mp,nbn2) = ',PL_BH(MP,NBN2)
                  print *,'rhol_bh(mp,nbn1) = ',RHOL_BH(MP,NBN1)
                  print *,'rhol_bh(mp,nbn2) = ',RHOL_BH(MP,NBN2)
                ENDIF
                UBBL(M,NCX) = RKLM*SIGN(ULX,HDLX)
              ENDDO
!
!---        Darcy flow when either interval is filled ---
!
            ELSE
              IZ1X = IZ_BH(NBN1)
              IZ2X = IZ_BH(NBN2)
!
!---          Loop over flux increments ---
!
              DO M = 1,ISVF
                MN = MNEG(M)
                MP = MPOS(M)
                HDLX = PL_BH(MP,NBN1) - PL_BH(MN,NBN2) +
     &            5.D-1*GRAV*(ZP_BH1X-ZP_BH2X)*
     &           (RHOL_BH(MP,NBN1)+RHOL_BH(MN,NBN2))
                IF( M.EQ.1 ) HDL = HDLX
                INDX = 11
                KLM = DIFMN(PERM(1,IZ1X),PERM(1,IZ2X),
     &            DBB1X,DBB2X,HDL,INDX)
                IF( PERM(1,IZ1X)/EPSL.LT.EPSL ) KLM = 0.D+0
                IF( PERM(1,IZ2X)/EPSL.LT.EPSL ) KLM = 0.D+0
                INDX = 8
                RKLM = DIFMN(RKL_BH(MP,NBN1),RKL_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDL,INDX)
                INDX = 5
                VLM = DIFMN(VISL_BH(MP,NBN1),VISL_BH(MN,NBN2),
     &            DBB1X,DBB2X,HDL,INDX)
                UBBL(M,NCX) = KLM*RKLM*HDLX/DBB_BH(NCX)/VLM
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
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
!---      Loop over fracture triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Apply fracture-fracture skin factor for inter-fracture
!           connections  ---
!
            SKFX = 0.D+0
            DO NF2X = 1,NF_FRC
              IF( NT2X.GE.IP_FRC(1,NF2X) .AND.
     &          NT2X.LE.IP_FRC(2,NF2X) ) THEN
                EXIT
              ENDIF
            ENDDO
            IF( NFX.NE.NF2X ) THEN
              SKFX = MAX( SKF_FRC(NFX),SKF_FRC(NF2X) )
            ENDIF
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              HDLX = PL_FRC(MP,NT1X) - PL_FRC(MN,NT2X) +
     &          5.D-1*GRAV*(ZP_FRC(NT1X)-ZP_FRC(NT2X))*
     &         (RHOL_FRC(MP,NT1X)+RHOL_FRC(MN,NT2X))
              IF( M.EQ.1 ) HDL = HDLX
              INDX = 11
              KLM = DIFMN(PERM_FRC(MP,NT1X),PERM_FRC(MN,NT2X),
     &          DFF1X,DFF2X,HDL,INDX)
              IF( PERM_FRC(MP,NT1X)/EPSL.LT.EPSL ) KLM = 0.D+0
              IF( PERM_FRC(MN,NT2X)/EPSL.LT.EPSL ) KLM = 0.D+0
              INDX = 8
              RKLM = DIFMN(RKL_FRC(MP,NT1X),RKL_FRC(MN,NT2X),
     &          DFF1X,DFF2X,HDL,INDX)
              INDX = 5
              VLM = DIFMN(VISL_FRC(MP,NT1X),VISL_FRC(MN,NT2X),
     &          DFF1X,DFF2X,HDL,INDX)
              UFFL(M,NCX) = KLM*RKLM*HDLX/DFF_FRC(NCX)/VLM/(1.D+0+SKFX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVL_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCV_BTF_GT
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
!     Compute the aqueous and gas Darcy flux from pressure gradients
!     and gravitational body forces from borehole nodes to
!     fracture triangles (borehole node is considered
!     the upper node for flux indexing)
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 7 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
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
      REAL*8 KGM,KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCV_BTF_GT'
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
          IZN = IZ_BH(NBN)
          SKFX = PAR_BH(1,INV_BH(NBN))
          INV = INV_BH(NBN)
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
!---      Loop over flux increments ---
!
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Aqueous volumetric flux  ---
!
            HDLX = PL_BH(MP,NBN) - PL_FRC(MN,NTX) +
     &        5.D-1*GRAV*(ZBHX-ZP_FRC(NTX))*
     &       (RHOL_BH(MP,NBN)+RHOL_FRC(MN,NTX))
            IF( M.EQ.1 ) HDL = HDLX
!
!---        Perforated pipe or uncased borehole  ---
!
            IF( IS_BH(INV).EQ.2 .OR. IS_BH(INV).EQ.12 ) THEN
              KLM = PERM_FRC(MN,NTX)
!
!---        Unperforated pipe or catalyst unperforated pipe  ---
!
            ELSEIF( IS_BH(INV).EQ.11 .OR. IS_BH(INV).EQ.21 ) THEN
              KLM = 0.D+0
!
!---        Filled uncased borehole  ---
!
            ELSEIF( IS_BH(INV).EQ.11 .OR. IS_BH(INV).EQ.21 ) THEN
              INDX = 11
              KLM = DIFMN(PERM(1,IZN),PERM_FRC(MN,NTX),
     &          DBF1X,DBF2X,HDL,INDX)
              IF( PERM(1,IZN)/EPSL.LT.EPSL ) KLM = 0.D+0
            ELSE
              KLM = 0.D+0
            ENDIF
            IF( PERM_FRC(MN,NTX)/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 8
            RKLM = DIFMN(RKL_BH(MP,NBN),RKL_FRC(MN,NTX),
     &        DBF1X,DBF2X,HDL,INDX)
            INDX = 5
            VLM = DIFMN(VISL_BH(MP,NBN),VISL_FRC(MN,NTX),
     &        DBF1X,DBF2X,HDL,INDX)
            UBFL(M,NCX) = KLM*RKLM*HDLX/DBFX/VLM/(1.D+0+SKFX)
!
!---        Gas volumetric flux  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              HDGX = PG_BH(MP,NBN) - PG_FRC(MN,NTX) +
     &          5.D-1*GRAV*(ZBHX-ZP_FRC(NTX))*
     &         (RHOG_BH(MP,NBN)+RHOG_FRC(MN,NTX))
              IF( M.EQ.1 ) HDG = HDGX
!
!---          Perforated pipe or uncased borehole  ---
!
              IF( IS_BH(INV).EQ.2 .OR. IS_BH(INV).EQ.12 ) THEN
                KGM = PERM_FRC(MN,NTX)
!
!---          Unperforated pipe or catalyst unperforated pipe  ---
!
              ELSEIF( IS_BH(INV).EQ.11 .OR. IS_BH(INV).EQ.21 ) THEN
                KGM = 0.D+0
!
!---          Filled uncased borehole  ---
!
              ELSEIF( IS_BH(INV).EQ.11 .OR. IS_BH(INV).EQ.21 ) THEN
                INDX = 11
                KGM = DIFMN(PERM(1,IZN),PERM_FRC(MN,NTX),
     &            DBF1X,DBF2X,HDG,INDX)
                IF( PERM(1,IZN)/EPSL.LT.EPSL ) KGM = 0.D+0
              ELSE
                KGM = 0.D+0
              ENDIF
              IF( PERM_FRC(MN,NTX)/EPSL.LT.EPSL ) KGM = 0.D+0
              INDX = 8
              RKGM = DIFMN(RKG_BH(MP,NBN),RKG_FRC(MN,NTX),
     &          DBF1X,DBF2X,HDG,INDX)
              INDX = 5
              VGM = DIFMN(VISG_BH(MP,NBN),VISG_FRC(MN,NTX),
     &          DBF1X,DBF2X,HDG,INDX)
              UBFG(M,NCX) = KGM*RKGM*HDGX/DBFX/VGM/(1.D+0+SKFX)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCV_BTF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_BH_GT
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
!     Compute borehole primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 21 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE PARM_BH
      USE JACOB
      USE GRID
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_BH_GT'
      EPSLX = 1.D-4
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
          IZN = IZ_BH(NBN)
          INV = INV_BH(NBN)
          N_DB = NBN
!
!---      Pipe flow  ---
!
          IF( IS_BH(INV).GT.10 ) THEN
            ENPR = 0.D+0
          ELSE
!
!---      Assign gas entry pressure and minimum gas saturation
!         for transition to unsaturated conditions  ---
!
            ENPR = 0.D+0
            IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
              ENPR = SCHR(1,IZN)*RHORL*GRAV
            ELSEIF( ISCHR(IZN).EQ.4 ) THEN
              ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
            ENDIF
          ENDIF
!
!---      Absolute temperature  ---
!
          TKX = T_BH(2,NBN) + TABS
!
!---      Saturated system w/o entrapped gas
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas air partial pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ_BH(2,NBN).EQ.1 ) THEN
            PG_BH(2,NBN) = PL_BH(2,NBN) + ENPR
            PLX = PL_BH(2,NBN) + PATM
            INDX = 0
            CALL REGION_4( T_BH(2,NBN),PSWX,INDX )
            PWX = MAX( PSWX,PLX )
            CALL SCF_GL( BTGLX,PWX )
            CALL SOL_BRNS( T_BH(2,NBN),PWX,XLSMX )
            XLSX = MIN( YLS_BH(2,NBN),XLSMX )
            XLS_BH(2,NBN) = XLSX
            CALL SP_B( T_BH(2,NBN),XLSX,PSBX )
            PVBX = PSBX
            PLX = PL_BH(2,NBN) + PATM
!            PGAX = XMLA_BH(2,NBN)*HCAW
            PGAX = PVA_BH(2,NBN)
!
!---        Isoair option no transition from aqueous
!           saturated conditions  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              PGX = PVBX + PGAX
            ELSE
              PGX = PLX
            ENDIF
            PX = MAX( PGX,PLX )
            CALL RKSP_BH_GT( PGX,PLX,RKGX,RKLX,SGX,SLX,NBN )
!
!---        Transition from aqueous saturated conditions to
!           aqueous unsaturated conditions: PC1 -> PC2  ---
!
            IF( SGX.GT.1.D-3 ) THEN
              SGX = 1.D-3
              SLX = 1.D+0 - SGX
              CALL CAP_BH_GT( CPGL,SLX,NBN )
              PCX = CPGL
              CALL P_IAPWS( T_BH(2,NBN),PSWX,RHOX,RHOLWX,HGWX,HLWX,
     &          UGWX,ULWX )
              CALL VPL_B( T_BH(2,NBN),PSWX,PCX,RHOLWX,PVW_BH(2,NBN) )
              CALL SCF_GL( BTGLX,PVW_BH(2,NBN) )
              CALL SOL_BRNS( T_BH(2,NBN),PSWX,XLSMX )
              XLSX = MIN( YLS_BH(2,NBN),XLSMX )
              XLS_BH(2,NBN) = XLSX
              CALL DENS_B( XLSX,RHOBX,RHOLWX,T_BH(2,NBN) )
              CALL SP_B( T_BH(2,NBN),XLSX,PSBX )
              CALL VPL_B( T_BH(2,NBN),PSBX,PCX,RHOBX,PVBX )
              PGX = PVBX + PGAX
              PG_BH(2,NBN) = PGX - PATM
              PL_BH(2,NBN) = PG_BH(2,NBN) - CPGL
              NPHAZ_BH(2,NBN) = 2
!
!---        No transition from aqueous saturated conditions  ---
!
            ELSE
              NPHAZ_BH(2,NBN) = 1
            ENDIF
!
!---      Unsaturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ_BH(2,NBN).EQ.2 ) THEN
            PLX = PL_BH(2,NBN) + PATM
            PGX = PG_BH(2,NBN) + PATM
            PX = MAX( PGX,PLX )
            INDX = 0
            CALL REGION_4( T_BH(2,NBN),PSWX,INDX )
            CALL P_IAPWS( T_BH(2,NBN),PSWX,RHOGWX,RHOLWX,HGWX,HLWX,
     &        UGWX,ULWX )
            PCX = PGX - PLX
            CALL VPL_B( T_BH(2,NBN),PSWX,PCX,RHOLWX,PVW_BH(2,NBN) )
            CALL SCF_GL( BTGLX,PVW_BH(2,NBN) )
            CALL SOL_BRNS( T_BH(2,NBN),PSWX,XLSMX )
            XLSX = MIN( YLS_BH(2,NBN),XLSMX )
            XLS_BH(2,NBN) = XLSX
            CALL DENS_B( XLSX,RHOBX,RHOLWX,T_BH(2,NBN) )
            CALL SP_B( T_BH(2,NBN),XLSX,PSBX )
            CALL VPL_B( T_BH(2,NBN),PSBX,PCX,RHOBX,PVW_BH(2,NBN) )
            PGAX = PGX - PVW_BH(2,NBN)
            CALL RKSP_BH_GT( PGX,PLX,RKGX,RKLX,SGX,SLX,NBH )
!
!---        Aqueous saturation = 1.0, transition to saturation
!           conditions  ---
!
            IF( ABS(1.D+0-SLX).LT.EPSL ) THEN
              PX = PL_BH(2,NBN)+PATM
              CALL SP_B( T_BH(2,NBN),XLSX,PSBX )
              PVBX = PSBX
              PG_BH(2,NBN) = PL_BH(2,NBN) + ENPR - EPSLX
              NPHAZ_BH(2,NBN) = 1
!
!---        No aqueous phase, transition
!           to fully unsaturated condition: PC2 -> PC4  ---
!
            ELSEIF( ( SLX.LT.EPSL .AND. (1.D+0-SL_BH(1,NBN)).GT.EPSL )
     &        .OR. TKX.GT.TCRW ) THEN
              PGX = PG_BH(2,NBN) + PATM
              INDX = 0
              CALL REGION_4( T_BH(2,NBN),PSWX,INDX )
              CALL SOL_BRNS( T_BH(2,NBN),PSWX,XLSMX )
              IF( TMS_BH(2,NBN).GT.EPSL ) THEN
                XLSX = XLSMX
              ELSE
                XLSX = 0.D+0
              ENDIF
              XLS_BH(2,NBN) = XLSX
              SL_BH(2,NBN) = 0.D+0
              CALL CAP_BH_GT( CPGL,SL_BH(2,NBN),NBN )
              PL_BH(2,NBN) = PG_BH(2,NBN) - PCX
              PLX = PGX - PCX
              CALL SP_B( T_BH(2,NBN),XLSX,PSBX )
              PX = MAX( PGX,PSBX )
              CALL P_IAPWS( T_BH(2,NBN),PX,RHOX,RHOLWX,HGWX,HLWX,
     &          UGWX,ULWX )
              CALL DENS_B( XLSX,RHOBX,RHOLWX,T_BH(2,NBN) )
              PCX = PGX - PLX
              CALL VPL_B( T_BH(2,NBN),PSBX,PCX,RHOBX,PVW_BH(2,NBN) )
              PVA_BH(2,NBN) = PGX - PVW_BH(2,NBN)
              NPHAZ_BH(2,NBN) = 4
!
!---        Gas pressure remains above gas-entry pressure, remain
!           as unsaturated condition  ---
!
            ELSE
              NPHAZ_BH(2,NBN) = 2
            ENDIF
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - gas pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ_BH(2,NBN).EQ.4 ) THEN
            PGX = PG_BH(2,NBN) + PATM
!
!---        Maximum salt solubility if aqueous phase existed  ---
!
            INDX = 0
            CALL REGION_4( T_BH(2,NBN),PSWX,INDX )
            CALL SOL_BRNS( T_BH(2,NBN),PSWX,XLSMX )
            IF( TMS_BH(2,NBN).GT.EPSL ) THEN
              XLSX = XLSMX
            ELSE
              XLSX = 0.D+0
            ENDIF
            XLS_BH(2,NBN) = XLSX
            PVA_BH(2,NBN) = PGX - PVW_BH(2,NBN)
            CALL SCF_GL( BTGLX,PVW_BH(2,NBN) )
            SL_BH(2,NBN) = 0.D+0
            CALL CAP_BH_GT( CPGL,SL_BH(2,NBN),NBN )
            PL_BH(2,NBN) = PG_BH(2,NBN) - PCX
            PLX = PGX - PCX
            CALL P_IAPWS( T_BH(2,NBN),PVW_BH(2,NBN),RHOGWX,RHOLWX,
     &        HGWX,HLWX,UGWX,ULWX )
!
!---        Saturated water vapor pressure given temperature and
!           capillary pressure  ---
!
            CALL SP_B( T_BH(2,NBN),XLSX,PSBX )
            CALL DENS_B( XLSX,RHOBX,RHOLWX,T_BH(2,NBN) )
            PCX = PGX - PLX
            CALL VPL_B( T_BH(2,NBN),PSBX,PCX,RHOBX,PVBX )
!
!---        Aqueous phase appears, transition to
!           unsaturated conditions: PC4 -> PC2  ---
!
            IF( PVW_BH(2,NBN).GT.PSBX .AND. TKX.LE.TCRW ) THEN
              SLX = 1.D-4
              CALL CAP_BH_GT( CPGL,SLX,NBN )
              PCX = CPGL
              CALL P_IAPWS( T_BH(2,NBN),PSWX,RHOX,RHOLWX,HX,HLWX,
     &          UX,ULWX )
              CALL VPL_B( T_BH(2,NBN),PSWX,PCX,RHOLWX,PVWX )
              CALL SCF_GL( BTGLX,PVWX )
              CALL SOL_BRNS( T_BH(2,NBN),PSWX,XLSMX )
              YLS_BH(2,NBN) = TMS_BH(2,NBN)/
     &          (RHOLWX*SLX+TMS_BH(2,NBN))
              XLSX = MIN( YLS_BH(2,NBN),XLSMX )
              XLS_BH(2,NBN) = XLSX
              CALL DENS_B( XLSX,RHOBX,RHOLWX,T_BH(2,NBN) )
              CALL SP_B( T_BH(2,NBN),XLSX,PSBX )
              CALL VPL_B( T_BH(2,NBN),PSBX,PCX,RHOBX,PVBX )
              PGX = PVBX + PVA_BH(2,NBN)
              PG_BH(2,NBN) = PGX - PATM
              PL_BH(2,NBN) = PG_BH(2,NBN) - CPGL
              NPHAZ_BH(2,NBN) = 2
!
!---        No aqueous phase, no transition from
!           fully unsaturated condition  ---
!
            ELSE
              NPHAZ_BH(2,NBN) = 4
            ENDIF
          ENDIF
!
!---      Compute increments  ---
!
!
!---      Saturated system w/o entrapped gas
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - aqueous-air mole fraction
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ_BH(2,NBN).EQ.1 ) THEN
!
!---        Energy - temperature
!
            DNR_BH(IEQT,NBN) = -1.D-7
!
!---        Water mass - aqueous pressure
!
            PLX = PL_BH(2,NBN) + PATM
            DNR_BH(IEQW,NBN) = MAX( 1.D-4,1.D-12*PLX )
!
!---        Air mass - gas air partial pressure
!
            IF( ISLC(37).EQ.0 ) THEN
!              XMLAX = MAX( PL_BH(2,NBN)+PATM,PATM )/HCAW
!              IF( XMLA_BH(2,NBN).GT.(1.D-2*XMLAX) ) THEN
!                DNR_BH(IEQA,NBN) = SIGN( 1.D-4*XMLAX,
!     &            5.D-1*XMLAX-XMLA_BH(2,NBN) )
!              ELSE
!                DNR_BH(IEQA,NBN) = SIGN( 1.D-3*XMLAX,
!     &            5.D-1*XMLAX-XMLA_BH(2,NBN) )
!              ENDIF
              DNR_BH(IEQA,NBN) = 1.D-1
            ENDIF
!
!---        NaCl mass - total NaCl brine mass fraction
!
            IF( ISLC(32).EQ.0 ) THEN
              CALL SOL_BRNS( T_BH(2,NBN),PLX,XLSMX )
              DNR_BH(IEQS,NBN) = 1.D-5*XLSMX
!
!---          Thermo-catalytic fluid, salt tracking  ---
!
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  )
     &          DNR_BH(IEQS,NBN) = MAX( 1.D-4*XLS_BH(2,NBN),1.D-12 )
            ENDIF
!
!---      Unsaturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ_BH(2,NBN).EQ.2 ) THEN
!
!---        Energy - temperature
!
            DNR_BH(IEQT,NBN) = -1.D-7
!
!---        Water mass - aqueous pressure
!
            DNRX = MAX( 1.D-4,1.D-12*ABS(PG_BH(2,NBN)-PL_BH(2,NBN)) )
            DNR_BH(IEQW,NBN) = -DNRX
!
!---        Air mass - gas pressure
!
            IF( ISLC(37).EQ.0 ) THEN
              DNR_BH(IEQA,NBN) = DNRX
            ENDIF
!
!---        NaCl mass - total NaCl brine mass fraction
!
            IF( ISLC(32).EQ.0 ) THEN
              INDX = 0
              CALL REGION_4( T_BH(2,NBN),PSWX,INDX )
              CALL SOL_BRNS( T_BH(2,NBN),PSWX,XLSMX )
              DNR_BH(IEQS,NBN) = 1.D-5*XLSMX
!
!---          Thermo-catalytic fluid, salt tracking  ---
!
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  )
     &          DNR_BH(IEQS,NBN) = MAX( 1.D-4*XLS_BH(2,NBN),1.D-12 )
            ENDIF
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - gas pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ_BH(2,NBN).EQ.4 ) THEN
!
!---        Energy - temperature  ---
!
            DNR_BH(IEQT,NBN) = -1.D-7
!
!---        Water mass - water vapor partial pressure  ---
!
            DPX = SIGN( MAX( 1.D-4,1.D-6*PVW_BH(2,NBN) ),
     &        1.D+0-PVW_BH(2,NBN) )
            DNR_BH(IEQW,NBN) = DPX
!
!---        Air mass - gas pressure  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              DNR_BH(IEQA,NBN) = MAX( 1.D-4,1.D-6*PVA_BH(2,NBN) )
            ENDIF
!
!---        NaCl mass - salt mass  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              DNR_BH(IEQS,NBN) = 1.D-6
            ENDIF
          ENDIF
!
!---      Increment the primary variables  ---
!
          DO M = 3,ISVC+2
            T_BH(M,NBN) = T_BH(2,NBN)
            PL_BH(M,NBN) = PL_BH(2,NBN)
            PG_BH(M,NBN) = PG_BH(2,NBN)
            PVW_BH(M,NBN) = PVW_BH(2,NBN)
            PVA_BH(M,NBN) = PVA_BH(2,NBN)
            XMLA_BH(M,NBN) = XMLA_BH(2,NBN)
            SG_BH(M,NBN) = SG_BH(2,NBN)
            SL_BH(M,NBN) = SL_BH(2,NBN)
            YLS_BH(M,NBN) = YLS_BH(2,NBN)
            TMS_BH(M,NBN) = TMS_BH(2,NBN)
!
!---        Saturated system w/o entrapped gas
!           Energy - temperature
!           Water mass - aqueous pressure
!           Air mass - aqueous-air mole fraction
!           NaCl mass - total NaCl brine mass fraction  ---
!
            IF( NPHAZ_BH(2,NBN).EQ.1 ) THEN
              IF( M.EQ.IEQT+2 ) THEN
                T_BH(M,NBN) = T_BH(M,NBN) + DNR_BH(IEQT,NBN)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                PL_BH(M,NBN) = PL_BH(M,NBN) + DNR_BH(IEQW,NBN)
                PG_BH(M,NBN) = PL_BH(M,NBN) + ENPR - EPSLX
              ELSEIF( M.EQ.IEQA+2 .AND. ISLC(37).EQ.0 ) THEN
!                XMLA_BH(M,NBN) = XMLA_BH(M,NBN) + DNR_BH(IEQA,NBN)
                PVA_BH(M,NBN) = PVA_BH(M,NBN) + DNR_BH(IEQA,NBN)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                YLS_BH(M,NBN) = YLS_BH(M,NBN) + DNR_BH(IEQS,NBN)
              ENDIF
!
!---        Unsaturated system
!
!           Energy - temperature
!           Water mass - aqueous pressure
!           Air mass - gas pressure
!           NaCl mass - total NaCl brine mass fraction  ---
!
            ELSEIF( NPHAZ_BH(2,NBN).EQ.2 ) THEN
              IF( M.EQ.IEQT+2 ) THEN
                T_BH(M,NBN) = T_BH(M,NBN) + DNR_BH(IEQT,NBN)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                PL_BH(M,NBN) = PL_BH(M,NBN) + DNR_BH(IEQW,NBN)
              ELSEIF( M.EQ.IEQA+2 .AND. ISLC(37).EQ.0 ) THEN
                PG_BH(M,NBN) = PG_BH(M,NBN) + DNR_BH(IEQA,NBN)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                YLS_BH(M,NBN) = YLS_BH(M,NBN) + DNR_BH(IEQS,NBN)
              ENDIF
!
!---        Fully unsaturated conditions
!
!           Energy - temperature
!           Water mass - water vapor partial pressure
!           Air mass - gas pressure
!           NaCl mass - salt mass  ---
!
            ELSEIF( NPHAZ_BH(2,NBN).EQ.4 ) THEN
              IF( M.EQ.IEQT+2 ) THEN
                T_BH(M,NBN) = T_BH(M,NBN) + DNR_BH(IEQT,NBN)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                PVW_BH(M,NBN) = PVW_BH(M,NBN) + DNR_BH(IEQW,NBN)
              ELSEIF( M.EQ.IEQA+2 .AND. ISLC(37).EQ.0 ) THEN
                PG_BH(M,NBN) = PG_BH(M,NBN) + DNR_BH(IEQA,NBN)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                TMS_BH(M,NBN) = TMS_BH(M,NBN) + DNR_BH(IEQS,NBN)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_FRC_GT
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
!     Compute frature primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 8 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE FDVS_FRC
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
      SUB_LOG(ISUB_LOG) = '/INCRM_FRC_GT'
      EPSLX = 1.D-4
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Assign gas entry pressure and minimum gas saturation
!       for transition to unsaturated conditions---
!
        ENPR = RKSP_FRC(1,NFX)*RHORL*GRAV
!
!---    Loop over fracture triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---      Absolute temperature  ---
!
          TKX = T_FRC(2,NTX) + TABS
!
!---      Saturated system w/o entrapped gas
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas air partial pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ_FRC(2,NTX).EQ.1 ) THEN
            PG_FRC(2,NTX) = PL_FRC(2,NTX) + ENPR
            PLX = PL_FRC(2,NTX) + PATM
            INDX = 0
            CALL REGION_4( T_FRC(2,NTX),PSWX,INDX )
            PWX = MAX( PSWX,PLX )
            CALL SCF_GL( BTGLX,PWX )
            CALL SOL_BRNS( T_FRC(2,NTX),PWX,XLSMX )
            XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
            XLS_FRC(2,NTX) = XLSX
            CALL SP_B( T_FRC(2,NTX),XLSX,PSBX )
            PVBX = PSBX
            PLX = PL_FRC(2,NTX) + PATM
!            PGAX = XMLA_FRC(2,NTX)*HCAW
            PGAX = PVA_FRC(2,NTX)
!
!---        Isoair option no transition from aqueous
!           saturated conditions  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              PGX = PVBX + PGAX
            ELSE
              PGX = PLX
            ENDIF
            PX = MAX( PGX,PLX )
            CALL RKSP_FRC_GT( PGX,PLX,RKGX,RKLX,SGX,SLX,NFX )
!
!---        Transition from aqueous saturated conditions to
!           aqueous unsaturated conditions: PC1 -> PC2  ---
!
            IF( SGX.GT.1.D-3 ) THEN
              SGX = 1.D-3
              SLX = 1.D+0 - SGX
              CALL CAP_FRC_GT( CPGL,SLX,NFX )
              PCX = CPGL
              CALL P_IAPWS( T_FRC(2,NTX),PSWX,RHOX,RHOLWX,HGWX,HLWX,
     &          UGWX,ULWX )
              CALL VPL_B( T_FRC(2,NTX),PSWX,PCX,RHOLWX,PVW_FRC(2,NTX) )
              CALL SCF_GL( BTGLX,PVW_FRC(2,NTX) )
              CALL SOL_BRNS( T_FRC(2,NTX),PSWX,XLSMX )
              XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
              XLS_FRC(2,NTX) = XLSX
              CALL DENS_B( XLSX,RHOBX,RHOLWX,T_FRC(2,NTX) )
              CALL SP_B( T_FRC(2,NTX),XLSX,PSBX )
              CALL VPL_B( T_FRC(2,NTX),PSBX,PCX,RHOBX,PVBX )
              PGX = PVBX + PGAX
              PG_FRC(2,NTX) = PGX - PATM
              PL_FRC(2,NTX) = PG_FRC(2,NTX) - CPGL
              NPHAZ_FRC(2,NTX) = 2
!
!---        No transition from aqueous saturated conditions  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 1
            ENDIF
!
!---      Unsaturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ_FRC(2,NTX).EQ.2 ) THEN
            PLX = PL_FRC(2,NTX) + PATM
            PGX = PG_FRC(2,NTX) + PATM
            PX = MAX( PGX,PLX )
            INDX = 0
            CALL REGION_4( T_FRC(2,NTX),PSWX,INDX )
            CALL P_IAPWS( T_FRC(2,NTX),PSWX,RHOGWX,RHOLWX,HGWX,HLWX,
     &        UGWX,ULWX )
            PCX = PGX - PLX
            CALL VPL_B( T_FRC(2,NTX),PSWX,PCX,RHOLWX,PVW_FRC(2,NTX) )
            CALL SCF_GL( BTGLX,PVW_FRC(2,NTX) )
            CALL SOL_BRNS( T_FRC(2,NTX),PSWX,XLSMX )
            XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
            XLS_FRC(2,NTX) = XLSX
            CALL DENS_B( XLSX,RHOBX,RHOLWX,T_FRC(2,NTX) )
            CALL SP_B( T_FRC(2,NTX),XLSX,PSBX )
            CALL VPL_B( T_FRC(2,NTX),PSBX,PCX,RHOBX,PVW_FRC(2,NTX) )
            PGAX = PGX - PVW_FRC(2,NTX)
            CALL RKSP_FRC_GT( PGX,PLX,RKGX,RKLX,SGX,SLX,NFX )
!
!---        Aqueous saturation = 1.0, transition to saturation
!           conditions  ---
!
            IF( ABS(1.D+0-SLX).LT.EPSL ) THEN
              PX = PL_FRC(2,NTX)+PATM
              CALL SP_B( T_FRC(2,NTX),XLSX,PSBX )
              PVBX = PSBX
              PG_FRC(2,NTX) = PL_FRC(2,NTX) + ENPR - EPSLX
              NPHAZ_FRC(2,NTX) = 1
!
!---        No aqueous phase, transition
!           to fully unsaturated condition: PC2 -> PC4  ---
!
            ELSEIF( ( SLX.LT.EPSL .AND. (1.D+0-SL_FRC(1,NTX)).GT.EPSL )
     &        .OR. TKX.GT.TCRW ) THEN
              PGX = PG_FRC(2,NTX) + PATM
              INDX = 0
              CALL REGION_4( T_FRC(2,NTX),PSWX,INDX )
              CALL SOL_BRNS( T_FRC(2,NTX),PSWX,XLSMX )
              IF( TMS_FRC(2,NTX).GT.EPSL ) THEN
                XLSX = XLSMX
              ELSE
                XLSX = 0.D+0
              ENDIF
              XLS_FRC(2,NTX) = XLSX
              SL_FRC(2,NTX) = 0.D+0
              CALL CAP_FRC_GT( CPGL,SL_FRC(2,NTX),NFX )
              PL_FRC(2,NTX) = PG_FRC(2,NTX) - PCX
              PLX = PGX - PCX
              CALL SP_B( T_FRC(2,NTX),XLSX,PSBX )
              PX = MAX( PGX,PSBX )
              CALL P_IAPWS( T_FRC(2,NTX),PX,RHOX,RHOLWX,HGWX,HLWX,
     &          UGWX,ULWX )
              CALL DENS_B( XLSX,RHOBX,RHOLWX,T_FRC(2,NTX) )
              PCX = PGX - PLX
              CALL VPL_B( T_FRC(2,NTX),PSBX,PCX,RHOBX,PVW_FRC(2,NTX) )
              PVA_FRC(2,NTX) = PGX - PVW_FRC(2,NTX)
              NPHAZ_FRC(2,NTX) = 4
!
!---        Gas pressure remains above gas-entry pressure, remain
!           as unsaturated condition  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 2
            ENDIF
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - gas pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ_FRC(2,NTX).EQ.4 ) THEN
            PGX = PG_FRC(2,NTX) + PATM
!
!---        Maximum salt solubility if aqueous phase existed  ---
!
            INDX = 0
            CALL REGION_4( T_FRC(2,NTX),PSWX,INDX )
            CALL SOL_BRNS( T_FRC(2,NTX),PSWX,XLSMX )
            IF( TMS_FRC(2,NTX).GT.EPSL ) THEN
              XLSX = XLSMX
            ELSE
              XLSX = 0.D+0
            ENDIF
            XLS_FRC(2,NTX) = XLSX
            PVA_FRC(2,NTX) = PGX - PVW_FRC(2,NTX)
            CALL SCF_GL( BTGLX,PVW_FRC(2,NTX) )
            SL_FRC(2,NTX) = 0.D+0
            CALL CAP_FRC_GT( CPGL,SL_FRC(2,NTX),NFX )
            PL_FRC(2,NTX) = PG_FRC(2,NTX) - PCX
            PLX = PGX - PCX
            CALL P_IAPWS( T_FRC(2,NTX),PVW_FRC(2,NTX),RHOGWX,RHOLWX,
     &        HGWX,HLWX,UGWX,ULWX )
!
!---        Saturated water vapor pressure given temperature and
!           capillary pressure  ---
!
            CALL SP_B( T_FRC(2,NTX),XLSX,PSBX )
            CALL DENS_B( XLSX,RHOBX,RHOLWX,T_FRC(2,NTX) )
            PCX = PGX - PLX
            CALL VPL_B( T_FRC(2,NTX),PSBX,PCX,RHOBX,PVBX )
!
!---        Aqueous phase appears, transition to
!           unsaturated conditions: PC4 -> PC2  ---
!
            IF( PVW_FRC(2,NTX).GT.PSBX .AND. TKX.LE.TCRW ) THEN
              SLX = 1.D-4
              CALL CAP_FRC_GT( CPGL,SLX,NFX )
              PCX = CPGL
              CALL P_IAPWS( T_FRC(2,NTX),PSWX,RHOX,RHOLWX,HX,HLWX,
     &          UX,ULWX )
              CALL VPL_B( T_FRC(2,NTX),PSWX,PCX,RHOLWX,PVWX )
              CALL SCF_GL( BTGLX,PVWX )
              CALL SOL_BRNS( T_FRC(2,NTX),PSWX,XLSMX )
              YLS_FRC(2,NTX) = TMS_FRC(2,NTX)/
     &          (RHOLWX*SLX+TMS_FRC(2,NTX))
              XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
              XLS_FRC(2,NTX) = XLSX
              CALL DENS_B( XLSX,RHOBX,RHOLWX,T_FRC(2,NTX) )
              CALL SP_B( T_FRC(2,NTX),XLSX,PSBX )
              CALL VPL_B( T_FRC(2,NTX),PSBX,PCX,RHOBX,PVBX )
              PGX = PVBX + PVA_FRC(2,NTX)
              PG_FRC(2,NTX) = PGX - PATM
              PL_FRC(2,NTX) = PG_FRC(2,NTX) - CPGL
              NPHAZ_FRC(2,NTX) = 2
!
!---        No aqueous phase, no transition from
!           fully unsaturated condition  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 4
            ENDIF
          ENDIF
!
!---      Compute increments  ---
!
!
!---      Saturated system w/o entrapped gas
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - aqueous-air mole fraction
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ_FRC(2,NTX).EQ.1 ) THEN
!
!---        Energy - temperature
!
            DNR_FRC(IEQT,NTX) = -1.D-7
!
!---        Water mass - aqueous pressure
!
            PLX = PL_FRC(2,NTX) + PATM
            DNR_FRC(IEQW,NTX) = MAX( 1.D-4,1.D-12*PLX )
!
!---        Air mass - aqueous-air mole fraction
!
            IF( ISLC(37).EQ.0 ) THEN
!              XMLAX = MAX( PL_FRC(2,NTX)+PATM,PATM )/HCAW
!              IF( XMLA_FRC(2,NTX).GT.(1.D-2*XMLAX) ) THEN
!                DNR_FRC(IEQA,NTX) = SIGN( 1.D-4*XMLAX,
!     &            5.D-1*XMLAX-XMLA_FRC(2,NTX) )
!              ELSE
!                DNR_FRC(IEQA,NTX) = SIGN( 1.D-3*XMLAX,
!     &            5.D-1*XMLAX-XMLA_FRC(2,NTX) )
!              ENDIF
              DNR_FRC(IEQA,NTX) = 1.D-1
            ENDIF
!
!---        NaCl mass - total NaCl brine mass fraction
!
            IF( ISLC(32).EQ.0 ) THEN
              CALL SOL_BRNS( T_FRC(2,NTX),PLX,XLSMX )
              DNR_FRC(IEQS,NTX) = 1.D-5*XLSMX
!
!---          Thermo-catalytic fluid, salt tracking  ---
!
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  )
     &          DNR_FRC(IEQS,NTX) = MAX( 1.D-4*XLS_FRC(2,NTX),1.D-12 )
            ENDIF
!
!---      Unsaturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ_FRC(2,NTX).EQ.2 ) THEN
!
!---        Energy - temperature
!
            DNR_FRC(IEQT,NTX) = -1.D-7
!
!---        Water mass - aqueous pressure
!
            DNRX = MAX( 1.D-4,1.D-12*ABS(PG_FRC(2,NTX)-PL_FRC(2,NTX)) )
            DNR_FRC(IEQW,NTX) = -DNRX
!
!---        Air mass - gas pressure
!
            IF( ISLC(37).EQ.0 ) THEN
              DNR_FRC(IEQA,NTX) = DNRX
            ENDIF
!
!---        NaCl mass - total NaCl brine mass fraction
!
            IF( ISLC(32).EQ.0 ) THEN
              INDX = 0
              CALL REGION_4( T_FRC(2,NTX),PSWX,INDX )
              CALL SOL_BRNS( T_FRC(2,NTX),PSWX,XLSMX )
              DNR_FRC(IEQS,NTX) = 1.D-5*XLSMX
!
!---          Thermo-catalytic fluid, salt tracking  ---
!
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  )
     &          DNR_FRC(IEQS,NTX) = MAX( 1.D-4*XLS_FRC(2,NTX),1.D-12 )
            ENDIF
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - gas pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ_FRC(2,NTX).EQ.4 ) THEN
!
!---        Energy - temperature  ---
!
            DNR_FRC(IEQT,NTX) = -1.D-7
!
!---        Water mass - water vapor partial pressure  ---
!
            DPX = SIGN( MAX( 1.D-4,1.D-6*PVW_FRC(2,NTX) ),
     &        1.D+0-PVW_FRC(2,NTX) )
            DNR_FRC(IEQW,NTX) = DPX
!
!---        Air mass - gas pressure  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              DNR_FRC(IEQA,NTX) = MAX( 1.D-4,1.D-6*PVA_FRC(2,NTX) )
            ENDIF
!
!---        NaCl mass - salt mass  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              DNR_FRC(IEQS,NTX) = 1.D-6
            ENDIF
          ENDIF
!
!---      Increment the primary variables  ---
!
          DO M = 3,ISVC+2
            T_FRC(M,NTX) = T_FRC(2,NTX)
            PL_FRC(M,NTX) = PL_FRC(2,NTX)
            PG_FRC(M,NTX) = PG_FRC(2,NTX)
            PVW_FRC(M,NTX) = PVW_FRC(2,NTX)
            PVA_FRC(M,NTX) = PVA_FRC(2,NTX)
            XMLA_FRC(M,NTX) = XMLA_FRC(2,NTX)
            SG_FRC(M,NTX) = SG_FRC(2,NTX)
            SL_FRC(M,NTX) = SL_FRC(2,NTX)
            YLS_FRC(M,NTX) = YLS_FRC(2,NTX)
            TMS_FRC(M,NTX) = TMS_FRC(2,NTX)
!
!---        Saturated system w/o entrapped gas
!           Energy - temperature
!           Water mass - aqueous pressure
!           Air mass - aqueous-air mole fraction
!           NaCl mass - total NaCl brine mass fraction  ---
!
            IF( NPHAZ_FRC(2,NTX).EQ.1 ) THEN
              IF( M.EQ.IEQT+2 ) THEN
                T_FRC(M,NTX) = T_FRC(M,NTX) + DNR_FRC(IEQT,NTX)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                PL_FRC(M,NTX) = PL_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
                PG_FRC(M,NTX) = PL_FRC(M,NTX) + ENPR - EPSLX
              ELSEIF( M.EQ.IEQA+2 .AND. ISLC(37).EQ.0 ) THEN
!                XMLA_FRC(M,NTX) = XMLA_FRC(M,NTX) + DNR_FRC(IEQA,NTX)
                PVA_FRC(M,NTX) = PVA_FRC(M,NTX) + DNR_FRC(IEQA,NTX)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                YLS_FRC(M,NTX) = YLS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ENDIF
!
!---        Unsaturated system
!
!           Energy - temperature
!           Water mass - aqueous pressure
!           Air mass - gas pressure
!           NaCl mass - total NaCl brine mass fraction  ---
!
            ELSEIF( NPHAZ_FRC(2,NTX).EQ.2 ) THEN
              IF( M.EQ.IEQT+2 ) THEN
                T_FRC(M,NTX) = T_FRC(M,NTX) + DNR_FRC(IEQT,NTX)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                PL_FRC(M,NTX) = PL_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
              ELSEIF( M.EQ.IEQA+2 .AND. ISLC(37).EQ.0 ) THEN
                PG_FRC(M,NTX) = PG_FRC(M,NTX) + DNR_FRC(IEQA,NTX)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                YLS_FRC(M,NTX) = YLS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ENDIF
!
!---        Fully unsaturated conditions
!
!           Energy - temperature
!           Water mass - water vapor partial pressure
!           Air mass - air partial pressure
!           NaCl mass - salt mass  ---
!
            ELSEIF( NPHAZ_FRC(2,NTX).EQ.4 ) THEN
              IF( M.EQ.IEQT+2 ) THEN
                T_FRC(M,NTX) = T_FRC(M,NTX) + DNR_FRC(IEQT,NTX)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                PVW_FRC(M,NTX) = PVW_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
              ELSEIF( M.EQ.IEQA+2 .AND. ISLC(37).EQ.0 ) THEN
                PVA_FRC(M,NTX) = PVA_FRC(M,NTX) + DNR_FRC(IEQA,NTX)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                TMS_FRC(M,NTX) = TMS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBA_BF_GT
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
!     Modify the Jacobian matrix for the borehole and fracture air
!     equations with aqueous-phase and gas-phase contributions, for
!     borehole to fracture connections (borehole node is considered
!     the upper node for flux indexing).
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 8 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RAP_FRC(LUK),RAA_FRC(LUK),RAP_BH(LUK),RAA_BH(LUK),FA(LSFV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBA_BF_GT'
      IEQAX = IEQA
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
!---      Borehole node to fracture triangle fluxes  ---
!
          DO M = 1,ISVF
            M1 = MNOD(M)
            M2 = MADJ(M)
            FLA1 = XLA_BH(M1,NBN)*RHOL_BH(M1,NBN)
            FLA2 = XLA_FRC(M2,NTX)*RHOL_FRC(M2,NTX)
            INDX = 2
            FLA = DIFMN( FLA1,FLA2,DBF1X,DBF2X,UBFL(1,NCX),INDX )
            FGA1 = XGA_BH(M1,NBN)*RHOG_BH(M1,NBN)
            FGA2 = XGA_FRC(M2,NTX)*RHOG_FRC(M2,NTX)
            INDX = 3
            FGA = DIFMN( FGA1,FGA2,DBF1X,DBF2X,UBFG(1,NCX),INDX )
            APX = GPI*2.D+1*PAR_BH(2,INV_BH(NBN))*APM_FRC(M2,NTX)
            FA(M) = APX*(UBFL(M,NCX)*FLA + UBFG(M,NCX)*FGA
     &          - WTMA*UBFDGW(M,NCX) + WTMA*UBFDLA(M,NCX))
          ENDDO
!
!---      Air equation residuals, borehole perspective  ---
!
          RAS_BH = FA(1)
          DO M = 1,ISVC
            MM = 2*M
            RAP_BH(M) = FA(MM)
            MM = 2*M + 1
            RAA_BH(M) = FA(MM)
          ENDDO
!
!---      Air equation residuals, fracture perspective  ---
!
          RAS_FRC = -FA(1)
          DO M = 1,ISVC
            MM = 2*M + 1
            RAP_FRC(M) = -FA(MM)
            MM = 2*M
            RAA_FRC(M) = -FA(MM)
          ENDDO
!
!---      Modify Jacobian Matrix  ---
!
          CALL JCBL_BF_GT( RAS_BH,RAP_BH,RAA_BH,
     &      RAS_FRC,RAP_FRC,RAA_FRC,NBN,NTX,NCX,IEQAX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBA_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBA_BH_GT
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
!     Load the Jacobian matrix for the borehole air equation with
!     aqueous-phase and gas-phase contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 21 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STAX(LUK+1),RAP(LUK),RAA(LUK,3),FA(LSFV,3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBA_BH_GT'
      IEQAX = IEQA
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
          INV1 = INV_BH(NBN1)
!
!---      Coaxial pipe  ---
!
          IF( IS_BH(INV1).EQ.10000 ) THEN
!
!---        Inner pipe hydraulic area  ---
!
            IF( IBN_BH(NBN1).EQ.0 ) THEN
              AP1X = GPI*(PAR_BH(1,INV1)**2)
!
!---        Outer pipe hydraulic area  ---
!
            ELSE
              AP1X = GPI*((PAR_BH(3,INV1)**2)-(PAR_BH(2,INV1)**2))
            ENDIF
!
!---      Pipe hydraulic area  ---
!
          ELSEIF( IS_BH(INV1).GT.10 ) THEN
            AP1X = GPI*(PAR_BH(3,INV1)**2)
!
!---      Porous media hydraulic area  ---
!
          ELSE
            IZN = IZ_BH(NBN1)
            AP1X = GPI*(PAR_BH(2,INV1)**2)
          ENDIF
!
!---      First-order, forward-difference, time differential  ---
!
          STA1 = PORD_BH(1,NBN1)*
     &      (XLA_BH(1,NBN1)*RHOL_BH(1,NBN1)*SL_BH(1,NBN1)
     &        + XGA_BH(1,NBN1)*RHOG_BH(1,NBN1)*SG_BH(1,NBN1))
          DO M = 1,ISVC+1
            MP = M + 1
            STA0 = PORD_BH(MP,NBN1)*
     &       (XLA_BH(MP,NBN1)*RHOL_BH(MP,NBN1)*SL_BH(MP,NBN1) +
     &        XGA_BH(MP,NBN1)*RHOG_BH(MP,NBN1)*SG_BH(MP,NBN1))
            STAX(M) = (STA0-STA1)*VOL_BH(NBN1)*DTI
          ENDDO
!
!---      Loop over borehole node to node connections ---
!
          NBCX = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBCX = NBCX + 1
            NBN2 = IBCM_BH(NCX)
            INV2 = INV_BH(NBN2)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---        Coaxial pipe  ---
!
            IF( IS_BH(INV2).EQ.10000 ) THEN
!
!---          Inner pipe hydraulic area  ---
!
              IF( IBN_BH(NBN2).EQ.0 ) THEN
                AP2X = GPI*(PAR_BH(1,INV2)**2)
!
!---          Outer pipe hydraulic area  ---
!
              ELSE
                AP2X = GPI*((PAR_BH(3,INV2)**2)-(PAR_BH(2,INV2)**2))
              ENDIF
!
!---        Pipe hydraulic area  ---
!
            ELSEIF( IS_BH(INV2).GT.10 ) THEN
              AP2X = GPI*(PAR_BH(3,INV2)**2)
!
!---        Porous media hydraulic area  ---
!
            ELSE
              IZN = IZ_BH(NBN2)
              AP2X = GPI*(PAR_BH(2,INV2)**2)
            ENDIF
            INDX = -3
            APX = DIFMN( AP1X,AP2X,DBB1X,DBB2X,ZERO,INDX )
!
!---        Borehole node surface fluxes  ---
!
            DO M = 1,ISVF
              M1 = MNOD(M)
              M2 = MADJ(M)
              FLA1 = XLA_BH(M1,NBN1)*RHOL_BH(M1,NBN1)
              FLA2 = XLA_BH(M2,NBN2)*RHOL_BH(M2,NBN2)
              INDX = 2
              FLA = DIFMN( FLA1,FLA2,DBB1X,DBB2X,UBBL(1,NCX),INDX )
              FGA1 = XGA_BH(M1,NBN1)*RHOG_BH(M1,NBN1)
              FGA2 = XGA_BH(M2,NBN2)*RHOG_BH(M2,NBN2)
              INDX = 3
              FGA = DIFMN( FGA1,FGA2,DBB1X,DBB2X,UBBG(1,NCX),INDX )
              FA(M,NBCX) = APX*
     &          (UBBL(M,NCX)*FLA + UBBG(M,NCX)*FGA
     &          - WTMA*UBBDGW(M,NCX) + WTMA*UBBDLA(M,NCX))
            ENDDO
          ENDDO
!
!---      Compute borehole node air equation residuals  ---
!
          RAS = STAX(1) - SRCA_BH(2,NBN1)
          DO MD = 1,NBCX
            RAS = RAS + FA(1,MD)
          ENDDO
          DO M = 1,ISVC
            RAP(M) = STAX(M+1) - SRCA_BH(M+2,NBN1)
            MM = 2*M
            DO MD = 1,NBCX
              RAP(M) = RAP(M) + FA(MM,MD)
            ENDDO
          ENDDO
          DO M = 1,ISVC
            MM = 2*M + 1
            DO MD = 1,NBCX
              RAA(M,MD) = RAS - FA(1,MD) + FA(MM,MD)
            ENDDO
            RAA(M,3) = RAS
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBL_BH_GT( RAS,RAP,RAA,NBN1,IEQAX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBA_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBA_FRC_GT
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
!     Load the Jacobian matrix for the fracture air equation with
!     aqueous-phase and gas-phase contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 14 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_FRC
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STAX(LUK+1),RAP(LUK),RAA(LUK,LTC_FRC),FA(LSFV,LTC_FRC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBA_FRC_GT'
      IEQAX = IEQA
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
!---      First-order, forward-difference, time differential  ---
!
          STA1 = (XLA_FRC(1,NT1X)*RHOL_FRC(1,NT1X)*SL_FRC(1,NT1X) +
     &        XGA_FRC(1,NT1X)*RHOG_FRC(1,NT1X)*SG_FRC(1,NT1X))
          DO M = 1,ISVC+1
            MP = M + 1
            STA0 = (XLA_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)*SL_FRC(MP,NT1X) +
     &        XGA_FRC(MP,NT1X)*RHOG_FRC(MP,NT1X)*SG_FRC(MP,NT1X))
            STAX(M) = (STA0*APM_FRC(MP,NT1X)-STA1*APM_FRC(1,NT1X))*
     &        DTI*AF_FRC(NT1X)
          ENDDO
!
!---      Loop over fracture triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NTCX = NCX - IPF_FRC(1,NT1X) + 1
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Fracture triangle surface fluxes  ---
!
            DO M = 1,ISVF
              M1 = MNOD(M)
              M2 = MADJ(M)
              FLA1 = XLA_FRC(M1,NT1X)*RHOL_FRC(M1,NT1X)
              FLA2 = XLA_FRC(M2,NT2X)*RHOL_FRC(M2,NT2X)
              INDX = 2
              FLA = DIFMN( FLA1,FLA2,DFF1X,DFF2X,UFFL(1,NCX),INDX )
              FGA1 = XGA_FRC(M1,NT1X)*RHOG_FRC(M1,NT1X)
              FGA2 = XGA_FRC(M2,NT2X)*RHOG_FRC(M2,NT2X)
              INDX = 3
              FGA = DIFMN( FGA1,FGA2,DFF1X,DFF2X,UFFG(1,NCX),INDX )
              INDX = -3
              APX = DIFMN( APM_FRC(M1,NT1X),APM_FRC(M2,NT2X),
     &          DFF1X,DFF2X,ZERO,INDX )
              FA(M,NTCX) = AFF_FRC(NCX)*APX*
     &          (UFFL(M,NCX)*FLA + UFFG(M,NCX)*FGA
     &          - WTMA*UFFDGW(M,NCX) + WTMA*UFFDLA(M,NCX))
            ENDDO
          ENDDO
!
!---      Compute fracture triangle air equation residuals  ---
!
          RAS = STAX(1) - SRCA_FRC(2,NT1X)
          DO MD = 1,NTCX
            RAS = RAS + FA(1,MD)
          ENDDO
          DO M = 1,ISVC
            RAP(M) = STAX(M+1) - SRCA_FRC(M+2,NT1X)
            MM = 2*M
            DO MD = 1,NTCX
              RAP(M) = RAP(M) + FA(MM,MD)
            ENDDO
          ENDDO
          DO M = 1,ISVC
            MM = 2*M + 1
            DO MD = 1,NTCX
              RAA(M,MD) = RAS - FA(1,MD) + FA(MM,MD)
            ENDDO
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBL_FRC_GT( RAS,RAP,RAA,NT1X,IEQAX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBA_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBA_MFB_GT
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
!     Modify Jacobian matrix for matrix grid cells,
!     fracture triangles, and borehole nodes for transfer of air
!     mass between matrix grid cells and fracture triangles and
!     borehole nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/JCBA_MFB_GT'
      IEQAX = IEQA
!
!---  Matrix equations, fracture connections,
!     loop over fractures  ---
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
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fracture triangle  ---
!
          N = INCM_FRC(NCX)
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSA_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQAX,NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of air mass residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSA_FRC(1,NCX)-TRNSA_FRC(MX,NCX))/DNR(M,N)
!
!---          Partial derivative of air mass residual with
!             respect to fracture triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSA_FRC(1,NCX)-TRNSA_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSA_FRC(1,NCX)
            RSDL(IEQAX,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSA_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQAX,NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of air mass residual with
!             respect to field node primary variables  ---
!
              MCOL = KLU(MP,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSA_FRC(1,NCX)-TRNSA_FRC(MX,NCX))/DNR(M,N)
!
!---          Partial derivative of air mass residual with
!             respect to fracture triangle primary variables  ---
!
              MC = (NCX-1)*ISVC + IEQAX
              MCOL = KLU_MCF(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSA_FRC(1,NCX)-TRNSA_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSA_FRC(1,NCX)
            RSDL(IEQAX,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO
!
!---  Matrix equations, borehole connections,
!     loop over boreholes  ---
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
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSA_BH is positive
!           from the borehole node to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQAX,NMD)
            DO M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
!
!---          Partial derivative of air mass residual with
!             respect to field node primary variables  ---
!
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSA_BH(1,NBN)-TRNSA_BH(MX,NBN))/DNR(M,N)
!
!---          Partial derivative of air mass residual with
!             respect to borehole node primary variables  ---
!
              MCOL = IM_BH(M,NBN)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSA_BH(1,NBN)-TRNSA_BH(MX,NBN))/DNR_BH(M,NBN)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSA_BH(1,NBN)
            RSDL(IEQAX,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSA_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQAX,NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of air mass residual with
!             respect to field node primary variables  ---
!
              MCOL = KLU(MP,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSA_BH(1,NBN)-TRNSA_BH(MX,NBN))/DNR(M,N)
!
!---          Partial derivative of air mass residual with
!             respect to borehole node primary variables  ---
!
              MC = (NBN-1)*ISVC + IEQAX
              MCOL = KLU_MCB(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSA_BH(1,NBN)-TRNSA_BH(MX,NBN))/DNR_BH(M,NBN)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSA_BH(1,NBN)
            RSDL(IEQAX,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
!
!---  Fracture equations, loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fracture triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        MTX = IXP_FRC(NTX)
        IF( MTX.EQ.0 ) CYCLE
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fracture triangle  ---
!
          N = INCM_FRC(NCX)
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fracture triangle, noting that TRNSA_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM_FRC(IEQAX,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of air mass residual with
!             respect to fracture triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSA_FRC(MX,NCX)-TRNSA_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of air mass residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSA_FRC(MX,NCX)-TRNSA_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSA_FRC(1,NCX)
            RSDL_FRC(IEQAX,NTX) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fracture triangle, noting that TRNSA_FRC is positive
!           from the fracture triangle to the field node  ---
!
            MP = IM_FRC(IEQAX,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of air mass residual with
!             respect to fracture triangle primary variables  ---
!
              NMD = (MTX-1)*ISVC + IEQAX
              MCOL = KLU_FRC(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSA_FRC(MX,NCX)-TRNSA_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of air mass residual with
!             respect to field node primary variables  ---
!
              NMD = (NCX-1)*ISVC + IEQAX
              MCOL = KLU_FCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSA_FRC(MX,NCX)-TRNSA_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSA_FRC(1,NCX)
            RSDL_FRC(IEQAX,NTX) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO
!
!---  Borehole equations, loop over boreholes  ---
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
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Borehole node, noting that TRNSA_BH is positive
!           from the borehole node to the field node  ---
!
            NMD = IXP(N)
            MP = IM_BH(IEQAX,NBN)
            DO M = 1,ISVC
!
!---          Partial derivative of air mass residual with
!             respect to borehole primary variables  ---
!
              MCOL = IM_BH(M,NBN)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSA_BH(MX,NBN)-TRNSA_BH(1,NBN))/DNR_BH(M,NBN)
!
!---          Partial derivative of air mass residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSA_BH(MX,NBN)-TRNSA_BH(1,NBN))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSA_BH(1,NBN)
            RSDL_BH(IEQAX,NBN) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Borehole node, noting that TRNSA_BH is positive
!           from the borehole node to the field node  ---
!
            MP = IM_BH(IEQAX,NBN)
            DO M = 1,ISVC
!
!---          Partial derivative of air mass residual with
!             respect to borehole primary variables  ---
!
              NMD = (NBN-1)*ISVC + IEQAX
              MCOL = KLU_BH(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSA_BH(MX,NBN)-TRNSA_BH(1,NBN))/DNR_BH(M,NBN)
!
!---          Partial derivative of air mass residual with
!             respect to field node primary variables  ---
!
              NMD = (NBN-1)*ISVC + IEQAX
              MCOL = KLU_BCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSA_BH(MX,NBN)-TRNSA_BH(1,NBN))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSA_BH(1,NBN)
            RSDL_BH(IEQAX,NBN) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBA_MFB_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBL_BF_GT( RSS_BH,RSP_BH,RSA_BH,
     &  RSS_FRC,RSP_FRC,RSA_FRC,NBN,NTX,NCX,MEQ )
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
!     Load the borehole Jacobian matrix
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 21 April 2019.
!






!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE FDVP_FRC
      USE FDVP_BH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)







!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP_BH(LUK),RSA_BH(LUK),RSP_FRC(LUK),RSA_FRC(LUK)




      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBL_BF_GT'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Borehole node  ---
!
        MP = IM_BH(MEQ,NBN)
        DO M = 1,ISVC
          MCOL = IM_BH(M,NBN)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &      (RSP_BH(M)-RSS_BH)/DNR_BH(M,NBN)
        ENDDO
        BLU(MP) = BLU(MP) - RSS_BH
        RSDL_BH(MEQ,NBN) = BLU(MP)
!
!---    Borehole node connection with fracture triangle ---
!
        DO M = 1,ISVC
          DNRX = DNR_FRC(M,NTX)
          MCOL = IM_FRC(M,NTX)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA_BH(M)-RSS_BH)/DNRX
        ENDDO
!
!---    Fracture triangle  ---
!
        MP = IM_FRC(MEQ,NTX)
        DO M = 1,ISVC
          MCOL = IM_FRC(M,NTX)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &      (RSP_FRC(M)-RSS_FRC)/DNR_FRC(M,NTX)
        ENDDO
        BLU(MP) = BLU(MP) - RSS_FRC
        RSDL_FRC(MEQ,NTX) = BLU(MP)
!
!---    Fracture triangle connection with borehole node ---
!
        DO M = 1,ISVC
          DNRX = DNR_BH(M,NBN)
          MCOL = IM_BH(M,NBN)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA_FRC(M)-RSS_FRC)/DNRX
        ENDDO
!
!---  SPLIB or Lis solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    Borehole node  ---
!
        MB = (NBN-1)*ISVC + MEQ
        DO M = 1,ISVC
          MCOL = KLU_BH(MB,M)
          DLU(MCOL) = DLU(MCOL) + (RSP_BH(M)-RSS_BH)/DNR_BH(M,NBN)
        ENDDO
        MP = IM_BH(MEQ,NBN)
        BLU(MP) = BLU(MP) - RSS_BH
        RSDL_BH(MEQ,NBN) = BLU(MP)
!
!---    Borehole node connection with fracture triangle ---
!
        MB = (NCX-1)*ISVC + MEQ
        DO M = 1,ISVC
          DNRX = DNR_FRC(M,NTX)
          MCOL = KLU_BCF(MB,M)
          DLU(MCOL) = DLU(MCOL) + (RSA_BH(M)-RSS_BH)/DNRX
        ENDDO
!
!---    Fracture triangle  ---
!
        MT1X = IXP_FRC(NTX)
        MB = (MT1X-1)*ISVC + MEQ
        DO M = 1,ISVC
          MCOL = KLU_FRC(MB,M)
          DLU(MCOL) = DLU(MCOL) + (RSP_FRC(M)-RSS_FRC)/DNR_FRC(M,NTX)
        ENDDO
        MP = IM_FRC(MEQ,NTX)
        BLU(MP) = BLU(MP) - RSS_FRC
        RSDL_FRC(MEQ,NTX) = BLU(MP)
!
!---    Fracture triangle connection with borehole node ---
!
        MB = (NCX-1)*ISVC + MEQ
        DO M = 1,ISVC
          DNRX = DNR_BH(M,NBN)
          MCOL = KLU_FCB(MB,M)
          DLU(MCOL) = DLU(MCOL) + (RSA_FRC(M)-RSS_FRC)/DNRX
        ENDDO

      ELSE
        INDX = 3
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'JCBL_BF_GT'
        RLMSGX = 0.D+0
        CHMSGX = 'Unknown Fracture Linear Equation Solver'
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBL_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBL_BH_GT( RSS,RSP,RSA,NBN1,MEQ )
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
!     Load the borehole Jacobian matrix
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 21 April 2019.
!






!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_BH
      USE JACOB
      USE GEOM_BH
      USE FDVP_BH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)







!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP(LUK),RSA(LUK,3)




      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBL_BH_GT'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Node  ---
!
        MP = IM_BH(MEQ,NBN1)
        DO M = 1,ISVC
          MCOL = IM_BH(M,NBN1)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSP(M)-RSS)/DNR_BH(M,NBN1)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL_BH(MEQ,NBN1) = BLU(MP)
!
!---    Loop over borehole node to node connections ---
!
        NBCX = 0
        DO ICX = 1,2
          NCX = IPB_BH(ICX,NBN1)
          IF( NCX.LE.0 ) CYCLE
          NBCX = NBCX + 1
          NBN2 = IBCM_BH(NCX)
          DO M = 1,ISVC
            DNRX = DNR_BH(M,NBN2)
            MCOL = IM_BH(M,NBN2)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,NBCX)-RSS)/DNRX
          ENDDO
        ENDDO
!
!---    Coaxial borehole node to node connections ---
!
        INV1 = INV_BH(NBN1)
        NCX = IPB_BH(3,NBN1) 
        IF( IS_BH(INV1).EQ.10000 .AND. NCX.GT.0 ) THEN
          NBN2 = IBCM_BH(NCX)
          DO M = 1,ISVC
            DNRX = DNR_BH(M,NBN2)
            MCOL = IM_BH(M,NBN2)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,3)-RSS)/DNRX
          ENDDO
        ENDIF
!
!---  SPLIB or Lis solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    Node  ---
!
        MB = (NBN1-1)*ISVC + MEQ
        MA = 0
        DO M = 1,ISVC
          MCOL = KLU_BH(MB,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR_BH(M,NBN1)
        ENDDO
        MP = IM_BH(MEQ,NBN1)
        BLU(MP) = BLU(MP) - RSS
        RSDL_BH(MEQ,NBN1) = BLU(MP)
!
!---    Loop over borehole node to node connections ---
!
        NBCX = 0
        DO ICX = 1,2
          NCX = IPB_BH(ICX,NBN1)
          IF( NCX.LE.0 ) CYCLE
          NBCX = NBCX + 1
          NBN2 = IBCM_BH(NCX)
          MA = MA + ISVC
          DO M = 1,ISVC
            DNRX = DNR_BH(M,NBN2)
            MCOL = KLU_BH(MB,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,NBCX)-RSS)/DNRX
          ENDDO
        ENDDO
!
!---    Coaxial borehole node to node connections ---
!
        INV1 = INV_BH(NBN1)
        NCX = IPB_BH(3,NBN1) 
        IF( IS_BH(INV1).EQ.10000 .AND. NCX.GT.0 ) THEN
          NBN2 = IBCM_BH(NCX)
          MA = MA + ISVC
          DO M = 1,ISVC
            DNRX = DNR_BH(M,NBN2)
            MCOL = KLU_BH(MB,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,3)-RSS)/DNRX
          ENDDO
        ENDIF

      ELSE
        INDX = 3
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'JCBL_BH_GT'
        RLMSGX = 0.D+0
        CHMSGX = 'Unknown Borehole Linear Equation Solver'
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBL_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBL_FRC_GT( RSS,RSP,RSA,NT1X,MEQ )
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
!     Load the fracture Jacobian matrix
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 14 September 2017.
!






!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE GEOM_FRC
      USE FDVP_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)







!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP(LUK),RSA(LUK,LTC_FRC)




      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBL_FRC_GT'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Node  ---
!
        MP = IM_FRC(MEQ,NT1X)
        DO M = 1,ISVC
          MCOL = IM_FRC(M,NT1X)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSP(M)-RSS)/DNR_FRC(M,NT1X)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL_FRC(MEQ,NT1X) = BLU(MP)
!
!---    Loop over fracture triangle connections  ---
!
        DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
          NTCX = NCX - IPF_FRC(1,NT1X) + 1
          NT2X = ITCM_FRC(NCX)
          IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
          DO M = 1,ISVC
            DNRX = DNR_FRC(M,NT2X)
            MCOL = IM_FRC(M,NT2X)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,NTCX)-RSS)/DNRX
          ENDDO
        ENDDO
!
!---  SPLIB or Lis solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    Node  ---
!
        MT1X = IXP_FRC(NT1X)
        MB = (MT1X-1)*ISVC + MEQ
        MA = 0
        DO M = 1,ISVC
          MCOL = KLU_FRC(MB,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR_FRC(M,NT1X)
        ENDDO
        MP = IM_FRC(MEQ,NT1X)
        BLU(MP) = BLU(MP) - RSS
        RSDL_FRC(MEQ,NT1X) = BLU(MP)
!
!---    Loop over fracture triangle connections  ---
!
        DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
          NTCX = NCX - IPF_FRC(1,NT1X) + 1
          NT2X = ITCM_FRC(NCX)
          IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
          MA = MA + ISVC
          DO M = 1,ISVC
            DNRX = DNR_FRC(M,NT2X)
            MCOL = KLU_FRC(MB,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,NTCX)-RSS)/DNRX
          ENDDO
        ENDDO

      ELSE
        INDX = 3
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'JCBL_FRC_GT'
        RLMSGX = 0.D+0
        CHMSGX = 'Unknown Fracture Linear Equation Solver'
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBL_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_BF_GT
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
!     Modify the Jacobian matrix for the borehole and fracture salt
!     equations with aqueous-phase and gas-phase contributions, for
!     borehole to fracture connections.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 8 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
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
      REAL*8 RSP_FRC(LUK),RSA_FRC(LUK),RSP_BH(LUK),RSA_BH(LUK),FS(LSFV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBS_BF_GT'
      IEQSX = IEQS
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
!---      Borehole node to fracture triangle fluxes  ---
!
          DO M = 1,ISVF
            M1 = MNOD(M)
            M2 = MADJ(M)
            APX = GPI*2.D+1*PAR_BH(2,INV_BH(NBN))*APM_FRC(M2,NTX)
            FS(M) = APX*UBFS(M,NCX)
          ENDDO
!
!---      Energy equation residuals, borehole perspective  ---
!
          RSS_BH = FS(1)
          DO M = 1,ISVC
            MM = 2*M
            RSP_BH(M) = FS(MM)
            MM = 2*M + 1
            RSA_BH(M) = FS(MM)
          ENDDO
!
!---      Energy equation residuals, fracture perspective  ---
!
          RSS_FRC = -FS(1)
          DO M = 1,ISVC
            MM = 2*M + 1
            RSP_FRC(M) = -FS(MM)
            MM = 2*M
            RSA_FRC(M) = -FS(MM)
          ENDDO
!
!---      Modify Jacobian Matrix  ---
!
          CALL JCBL_BF_GT( RSS_BH,RSP_BH,RSA_BH,
     &      RSS_FRC,RSP_FRC,RSA_FRC,NBN,NTX,NCX,IEQSX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_BH_GT
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
!     Load the Jacobian matrix for the borehole salt equation with
!     aqueous-phase and gas-phase contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 21 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STSX(LUK+1),RSP(LUK),RSA(LUK,3),FS(LSFV,3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBS_BH_GT'
      IEQSX = IEQS
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
          INV1 = INV_BH(NBN1)
!
!---      Coaxial pipe  ---
!
          IF( IS_BH(INV1).EQ.10000 ) THEN
!
!---        Inner pipe hydraulic area  ---
!
            IF( IBN_BH(NBN1).EQ.0 ) THEN
              AP1X = GPI*(PAR_BH(1,INV1)**2)
!
!---        Outer pipe hydraulic area  ---
!
            ELSE
              AP1X = GPI*((PAR_BH(3,INV1)**2)-(PAR_BH(2,INV1)**2))
            ENDIF
!
!---      Pipe hydraulic area  ---
!
          ELSEIF( IS_BH(INV1).GT.10 ) THEN
            AP1X = GPI*(PAR_BH(3,INV1)**2)
!
!---      Porous media hydraulic area  ---
!
          ELSE
            IZN = IZ_BH(NBN1)
            AP1X = GPI*(PAR_BH(2,INV1)**2)
          ENDIF
!
!---      First-order, forward-difference, time differential  ---
!
          DO M = 1,ISVC+1
            MP = M + 1
            STSX(M) = (TMS_BH(MP,NBN1)-TMS_BH(1,NBN1))*VOL_BH(NBN1)*DTI
          ENDDO
!
!---      Loop over borehole node to node connections ---
!
          NBCX = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBCX = NBCX + 1
            NBN2 = IBCM_BH(NCX)
            INV2 = INV_BH(NBN2)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---        Coaxial pipe  ---
!
            IF( IS_BH(INV2).EQ.10000 ) THEN
!
!---          Inner pipe hydraulic area  ---
!
              IF( IBN_BH(NBN2).EQ.0 ) THEN
                AP2X = GPI*(PAR_BH(1,INV2)**2)
!
!---          Outer pipe hydraulic area  ---
!
              ELSE
                AP2X = GPI*((PAR_BH(3,INV2)**2)-(PAR_BH(2,INV2)**2))
              ENDIF
!
!---        Pipe hydraulic area  ---
!
            ELSEIF( IS_BH(INV2).GT.10 ) THEN
              AP2X = GPI*(PAR_BH(3,INV2)**2)
!
!---        Porous media  ---
!
            ELSE
              IZN = IZ_BH(NBN2)
              AP2X = GPI*(PAR_BH(2,INV2)**2)
            ENDIF
            INDX = -3
            APX = DIFMN( AP1X,AP2X,DBB1X,DBB2X,ZERO,INDX )
!
!---        Borehole node surface fluxes  ---
!
            DO M = 1,ISVF
              M1 = MNOD(M)
              M2 = MADJ(M)
              FS(M,NBCX) = APX*UBBS(M,NCX)
            ENDDO
          ENDDO
!
!---      Compute borehole node salt equation residuals  ---
!
          RSS = STSX(1) - SRCS_BH(2,NBN1)
          DO MD = 1,NBCX
            RSS = RSS + FS(1,MD)
          ENDDO
          DO M = 1,ISVC
            RSP(M) = STSX(M+1) - SRCS_BH(M+2,NBN1)
            MM = 2*M
            DO MD = 1,NBCX
              RSP(M) = RSP(M) + FS(MM,MD)
            ENDDO
          ENDDO
          DO M = 1,ISVC
            MM = 2*M + 1
            DO MD = 1,NBCX
              RSA(M,MD) = RSS - FS(1,MD) + FS(MM,MD)
            ENDDO
            RSA(M,3) = RSS
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBL_BH_GT( RSS,RSP,RSA,NBN1,IEQSX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_FRC_GT
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
!     Load the Jacobian matrix for the fracture salt equation with
!     aqueous-phase and gas-phase contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 14 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_FRC
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STSX(LUK+1),RSP(LUK),RSA(LUK,LTC_FRC),FS(LSFV,LTC_FRC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBS_FRC_GT'
      IEQSX = IEQS
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
!---      First-order, forward-difference, time differential  ---
!
          DO M = 1,ISVC+1
            MP = M + 1
            STSX(M) = (TMS_FRC(MP,NT1X)*APM_FRC(MP,NT1X)-
     &        TMS_FRC(1,NT1X)*APM_FRC(1,NT1X))*DTI*AF_FRC(NT1X)
          ENDDO
!
!---      Loop over fracture triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NTCX = NCX - IPF_FRC(1,NT1X) + 1
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Fracture triangle surface fluxes  ---
!
            DO M = 1,ISVF
              M1 = MNOD(M)
              M2 = MADJ(M)
              INDX = -3
              APX = DIFMN( APM_FRC(M1,NT1X),APM_FRC(M2,NT2X),
     &          DFF1X,DFF2X,ZERO,INDX )
              FS(M,NTCX) = AFF_FRC(NCX)*APX*UFFS(M,NCX)
            ENDDO
          ENDDO
!
!---      Compute fracture triangle salt equation residuals  ---
!
          RSS = STSX(1) - SRCS_FRC(2,NT1X)
          DO MD = 1,NTCX
            RSS = RSS + FS(1,MD)
          ENDDO
          DO M = 1,ISVC
            RSP(M) = STSX(M+1) - SRCS_FRC(M+2,NT1X)
            MM = 2*M
            DO MD = 1,NTCX
              RSP(M) = RSP(M) + FS(MM,MD)
            ENDDO
          ENDDO
          DO M = 1,ISVC
            MM = 2*M + 1
            DO MD = 1,NTCX
              RSA(M,MD) = RSS - FS(1,MD) + FS(MM,MD)
            ENDDO
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBL_FRC_GT( RSS,RSP,RSA,NT1X,IEQSX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_MFB_GT
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
!     Modify Jacobian matrix for matrix grid cells,
!     fracture triangles, and borehole nodes for transfer of salt
!     mass between matrix grid cells and fracture triangles and
!     borehole nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/JCBS_MFB_GT'
      IEQSX = IEQS
!
!---  Matrix equations, fracture connections,
!     loop over fractures  ---
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
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fracture triangle  ---
!
          N = INCM_FRC(NCX)
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSS_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQSX,NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of salt mass residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSS_FRC(1,NCX)-TRNSS_FRC(MX,NCX))/DNR(M,N)
!
!---          Partial derivative of salt mass residual with
!             respect to fracture triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSS_FRC(1,NCX)-TRNSS_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSS_FRC(1,NCX)
            RSDL(IEQSX,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSS_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQSX,NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of salt mass residual with
!             respect to field node primary variables  ---
!
              MCOL = KLU(MP,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSS_FRC(1,NCX)-TRNSS_FRC(MX,NCX))/DNR(M,N)
!
!---          Partial derivative of salt mass residual with
!             respect to fracture triangle primary variables  ---
!
              MC = (NCX-1)*ISVC + IEQSX
              MCOL = KLU_MCF(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSS_FRC(1,NCX)-TRNSS_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSS_FRC(1,NCX)
            RSDL(IEQSX,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO
!
!---  Matrix equations, borehole connections,
!     loop over boreholes  ---
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
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSS_BH is positive
!           from the borehole node to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQSX,NMD)
            DO M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
!
!---          Partial derivative of salt mass residual with
!             respect to field node primary variables  ---
!
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSS_BH(1,NBN)-TRNSS_BH(MX,NBN))/DNR(M,N)
!
!---          Partial derivative of salt mass residual with
!             respect to borehole node primary variables  ---
!
              MCOL = IM_BH(M,NBN)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSS_BH(1,NBN)-TRNSS_BH(MX,NBN))/DNR_BH(M,NBN)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSS_BH(1,NBN)
            RSDL(IEQSX,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSS_BH is positive
!           from the borehole node to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQSX,NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of salt mass residual with
!             respect to field node primary variables  ---
!
              MCOL = KLU(MP,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSS_BH(1,NBN)-TRNSS_BH(MX,NBN))/DNR(M,N)
!
!---          Partial derivative of salt mass residual with
!             respect to borehole node primary variables  ---
!
              MC = (NBN-1)*ISVC + IEQSX
              MCOL = KLU_MCB(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSS_BH(1,NBN)-TRNSS_BH(MX,NBN))/DNR_BH(M,NBN)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSS_BH(1,NBN)
            RSDL(IEQSX,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
!
!---  Fracture equations, loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fracture triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        MTX = IXP_FRC(NTX)
        IF( MTX.EQ.0 ) CYCLE
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fracture triangle  ---
!
          N = INCM_FRC(NCX)
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fracture triangle, noting that TRNSS_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM_FRC(IEQSX,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of salt mass residual with
!             respect to fracture triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSS_FRC(MX,NCX)-TRNSS_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of salt mass residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSS_FRC(MX,NCX)-TRNSS_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSS_FRC(1,NCX)
            RSDL_FRC(IEQSX,NTX) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fracture triangle, noting that TRNSS_FRC is positive
!           from the fracture triangle to the field node  ---
!
            MP = IM_FRC(IEQSX,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of salt mass residual with
!             respect to fracture triangle primary variables  ---
!
              NMD = (MTX-1)*ISVC + IEQSX
              MCOL = KLU_FRC(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSS_FRC(MX,NCX)-TRNSS_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of salt mass residual with
!             respect to field node primary variables  ---
!
              NMD = (NCX-1)*ISVC + IEQSX
              MCOL = KLU_FCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSS_FRC(MX,NCX)-TRNSS_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSS_FRC(1,NCX)
            RSDL_FRC(IEQSX,NTX) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO
!
!---  Borehole equations, loop over boreholes  ---
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
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Borehole node, noting that TRNSS_BH is positive
!           from the borehole node to the field node  ---
!
            NMD = IXP(N)
            MP = IM_BH(IEQSX,NBN)
            DO M = 1,ISVC
!
!---          Partial derivative of salt mass residual with
!             respect to borehole primary variables  ---
!
              MCOL = IM_BH(M,NBN)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSS_BH(MX,NBN)-TRNSS_BH(1,NBN))/DNR_BH(M,NBN)
!
!---          Partial derivative of salt mass residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSS_BH(MX,NBN)-TRNSS_BH(1,NBN))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSS_BH(1,NBN)
            RSDL_BH(IEQSX,NBN) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Borehole node, noting that TRNSS_BH is positive
!           from the borehole node to the field node  ---
!
            MP = IM_BH(IEQSX,NBN)
            DO M = 1,ISVC
!
!---          Partial derivative of salt mass residual with
!             respect to borehole primary variables  ---
!
              NMD = (NBN-1)*ISVC + IEQSX
              MCOL = KLU_BH(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSS_BH(MX,NBN)-TRNSS_BH(1,NBN))/DNR_BH(M,NBN)
!
!---          Partial derivative of salt mass residual with
!             respect to field node primary variables  ---
!
              NMD = (NBN-1)*ISVC + IEQSX
              MCOL = KLU_BCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSS_BH(MX,NBN)-TRNSS_BH(1,NBN))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSS_BH(1,NBN)
            RSDL_BH(IEQSX,NBN) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_MFB_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_BF_GT
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
!     Modify the Jacobian matrix for the borehole and fracture energy
!     equations with aqueous-phase and gas-phase contributions, for
!     borehole to fracture connections.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 8 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
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
      REAL*8 RTP_FRC(LUK),RTA_FRC(LUK),RTP_BH(LUK),RTA_BH(LUK),FT(LSFV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT_BF_GT'
      IEQTX = IEQT
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
!---      Borehole node to fracture triangle fluxes  ---
!
          DO M = 1,ISVF
            M1 = MNOD(M)
            M2 = MADJ(M)
            APX = GPI*2.D+1*PAR_BH(2,INV_BH(NBN))*APM_FRC(M2,NTX)
            FT(M) = APX*UBFQ(M,NCX)
          ENDDO
!
!---      Energy equation residuals, borehole perspective  ---
!
          RTS_BH = FT(1)
          DO M = 1,ISVC
            MM = 2*M
            RTP_BH(M) = FT(MM)
            MM = 2*M + 1
            RTA_BH(M) = FT(MM)
          ENDDO
!
!---      Energy equation residuals, fracture perspective  ---
!
          RTS_FRC = -FT(1)
          DO M = 1,ISVC
            MM = 2*M + 1
            RTP_FRC(M) = -FT(MM)
            MM = 2*M
            RTA_FRC(M) = -FT(MM)
          ENDDO
!
!---      Modify Jacobian Matrix  ---
!
          CALL JCBL_BF_GT( RTS_BH,RTP_BH,RTA_BH,
     &      RTS_FRC,RTP_FRC,RTA_FRC,NBN,NTX,NCX,IEQTX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_BH_GT
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
!     Load the Jacobian matrix for the borehole energy equation with
!     aqueous-phase and gas-phase contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 22 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_BH
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
      REAL*8 STTX(LUK+1),RTP(LUK),RTA(LUK,3),FT(LSFV,3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT_BH_GT'
      IEQTX = IEQT
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
          INV1 = INV_BH(NBN1)
!
!---      Coaxial pipe  ---
!
          IF( IS_BH(INV1).EQ.10000 ) THEN
!
!---        Inner pipe hydraulic area  ---
!
            IF( IBN_BH(NBN1).EQ.0 ) THEN
              CP_BHX = PAR_BH(6,INV1)
              AP1X = GPI*(PAR_BH(1,INV1)**2)
!
!---        Outer pipe hydraulic area  ---
!
            ELSE
              CP_BHX = PAR_BH(15,INV1)
              AP1X = GPI*((PAR_BH(3,INV1)**2)-(PAR_BH(2,INV1)**2))
            ENDIF
!
!---        First-order, forward-difference, time differential  ---
!
            STT1 = (1.D+0-PORD_BH(1,NBN1))*CP_BHX*T_BH(1,NBN1)
     &        + PORD_BH(1,NBN1)*
     &        (SS_BH(1,NBN1)*RHOSP_BH(1,NBN1)*HSP_BH(1,NBN1) +
     &        SL_BH(1,NBN1)*RHOL_BH(1,NBN1)*HL_BH(1,NBN1) +
     &        SG_BH(1,NBN1)*RHOG_BH(1,NBN1)*UEG_BH(1,NBN1))
            DO M = 1,ISVC+1
              MP = M + 1
              STT0 = (1.D+0-PORD_BH(MP,NBN1))*CP_BHX*T_BH(MP,NBN1)
     &        + PORD_BH(MP,NBN1)*
     &        (SS_BH(MP,NBN1)*RHOSP_BH(MP,NBN1)*HSP_BH(MP,NBN1) +
     &          SL_BH(MP,NBN1)*RHOL_BH(MP,NBN1)*HL_BH(MP,NBN1) +
     &          SG_BH(MP,NBN1)*RHOG_BH(MP,NBN1)*UEG_BH(MP,NBN1))
              STTX(M) = (STT0-STT1)*DTI*VOL_BH(NBN1)
            ENDDO
!
!---      Pipe hydraulic area  ---
!
          ELSEIF( IS_BH(INV1).GT.10 ) THEN
            AP1X = GPI*(PAR_BH(3,INV1)**2)
!
!---        First-order, forward-difference, time differential  ---
!
            STT1 = (1.D+0-PORD_BH(1,NBN1))*PAR_BH(5,INV1)*T_BH(1,NBN1)
     &        + PORD_BH(1,NBN1)*
     &        (SS_BH(1,NBN1)*RHOSP_BH(1,NBN1)*HSP_BH(1,NBN1) +
     &        SL_BH(1,NBN1)*RHOL_BH(1,NBN1)*HL_BH(1,NBN1) +
     &        SG_BH(1,NBN1)*RHOG_BH(1,NBN1)*UEG_BH(1,NBN1))
            DO M = 1,ISVC+1
              MP = M + 1
              STT0 = (1.D+0-PORD_BH(MP,NBN1))*PAR_BH(5,INV1)*
     &        T_BH(MP,NBN1) + PORD_BH(MP,NBN1)*
     &        (SS_BH(MP,NBN1)*RHOSP_BH(MP,NBN1)*HSP_BH(MP,NBN1) +
     &          SL_BH(MP,NBN1)*RHOL_BH(MP,NBN1)*HL_BH(MP,NBN1) +
     &          SG_BH(MP,NBN1)*RHOG_BH(MP,NBN1)*UEG_BH(MP,NBN1))
              STTX(M) = (STT0-STT1)*DTI*VOL_BH(NBN1)
            ENDDO
!
!---      Porous media hydraulic area  ---
!
          ELSE
            IZN = IZ_BH(NBN1)
            AP1X = GPI*(PAR_BH(2,INV1)**2)
!
!---        First-order, forward-difference, time differential  ---
!
            STT1 = (1.D+0-PORD_BH(1,NBN1))*RHOS(IZN)*CPS(IZN)*
     &        T_BH(1,NBN1) + PORD_BH(1,NBN1)*
     &        (SS_BH(1,NBN1)*RHOSP_BH(1,NBN1)*HSP_BH(1,NBN1) +
     &        SL_BH(1,NBN1)*RHOL_BH(1,NBN1)*HL_BH(1,NBN1) +
     &        SG_BH(1,NBN1)*RHOG_BH(1,NBN1)*UEG_BH(1,NBN1))
            DO M = 1,ISVC+1
              MP = M + 1
              STT0 = (1.D+0-PORD_BH(MP,NBN1))*RHOS(IZN)*CPS(IZN)*
     &        T_BH(MP,NBN1) + PORD_BH(MP,NBN1)*
     &        (SS_BH(MP,NBN1)*RHOSP_BH(MP,NBN1)*HSP_BH(MP,NBN1) +
     &          SL_BH(MP,NBN1)*RHOL_BH(MP,NBN1)*HL_BH(MP,NBN1) +
     &          SG_BH(MP,NBN1)*RHOG_BH(MP,NBN1)*UEG_BH(MP,NBN1))
              STTX(M) = (STT0-STT1)*DTI*VOL_BH(NBN1)
            ENDDO
          ENDIF
!
!---      Loop over borehole node to node connections ---
!
          NBCX = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBCX = NBCX + 1
            NBN2 = IBCM_BH(NCX)
            INV2 = INV_BH(NBN2)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---        Coaxial pipe  ---
!
            IF( IS_BH(INV2).EQ.10000 ) THEN
!
!---          Inner pipe hydraulic area  ---
!
              IF( IBN_BH(NBN2).EQ.0 ) THEN
                AP2X = GPI*(PAR_BH(1,INV2)**2)
!
!---          Outer pipe hydraulic area  ---
!
              ELSE
                AP2X = GPI*((PAR_BH(3,INV2)**2)-(PAR_BH(2,INV2)**2))
              ENDIF
!
!---        Pipe hydraulic area  ---
!
            ELSEIF( IS_BH(INV2).GT.10 ) THEN
              AP2X = GPI*(PAR_BH(3,INV2)**2)
!
!---        Porous media hydraulic area  ---
!
            ELSE
              IZN = IZ_BH(NBN2)
              AP2X = GPI*(PAR_BH(2,INV2)**2)
            ENDIF
            INDX = -3
            APX = DIFMN( AP1X,AP2X,DBB1X,DBB2X,ZERO,INDX )
!
!---        Borehole node surface fluxes  ---
!
            DO M = 1,ISVF
              FT(M,NBCX) = APX*UBBQ(M,NCX)
            ENDDO
          ENDDO
!
!---      Special coaxial connections  ---
!
          IF( IS_BH(INV1).EQ.10000 ) THEN
            NCX = IPB_BH(3,NBN1)
            AFB_BHX = 2.D+0*GPI*PAR_BH(2,INV1)*
     &        SQRT( (XP_BH(2,NBN1)-XP_BH(1,NBN1))**2 +
     &        (YP_BH(2,NBN1)-YP_BH(1,NBN1))**2 +
     &        (ZP_BH(2,NBN1)-ZP_BH(1,NBN1))**2 )
            IF( NCX.EQ.-1 ) THEN
              DO M = 1,ISVF
                FT(M,1) = FT(M,1) + AFB_BHX*UBBQC(M,NBN1)
                FT(M,3) = 0.D+0
              ENDDO
            ELSEIF( NCX.EQ.-2 ) THEN
              DO M = 1,ISVF
                FT(M,2) = FT(M,2) + AFB_BHX*UBBQC(M,NBN1)
                FT(M,3) = 0.D+0
              ENDDO
            ELSE
              DO M = 1,ISVF
                FT(M,3) = AFB_BHX*UBBQC(M,NBN1)
              ENDDO
            ENDIF
          ELSE
            DO M = 1,ISVF
              FT(M,3) = 0.D+0
            ENDDO
          ENDIF
!
!---      Compute borehole node energy equation residuals  ---
!
          RTS = STTX(1) - SRCT_BH(2,NBN1)
          DO MD = 1,NBCX
            RTS = RTS + FT(1,MD)
          ENDDO
          IF( IS_BH(INV1).EQ.10000 ) RTS = RTS + FT(1,3)
          DO M = 1,ISVC
            RTP(M) = STTX(M+1) - SRCT_BH(M+2,NBN1)
            MM = 2*M
            DO MD = 1,NBCX
              RTP(M) = RTP(M) + FT(MM,MD)
            ENDDO
            IF( IS_BH(INV1).EQ.10000 ) RTP(M) = RTP(M) + FT(MM,3)
          ENDDO
          DO M = 1,ISVC
            MM = 2*M + 1
            DO MD = 1,NBCX
              RTA(M,MD) = RTS - FT(1,MD) + FT(MM,MD)
            ENDDO
            IF( IS_BH(INV1).EQ.10000 ) RTA(M,3) = RTS-FT(1,3)+FT(MM,3)
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBL_BH_GT( RTS,RTP,RTA,NBN1,IEQTX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_FRC_GT
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
!     Load the Jacobian matrix for the fracture energy equation with
!     aqueous-phase and gas-phase contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 14 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_FRC
      USE JACOB
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVT_FRC
      USE FDVS_FRC
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
      REAL*8 STTX(LUK+1),RTP(LUK),RTA(LUK,LTC_FRC),FT(LSFV,LTC_FRC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT_FRC_GT'
      IEQTX = IEQT
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
!---      First-order, forward-difference, time differential  ---
!
          STT1 = SS_FRC(1,NT1X)*RHOSP_FRC(1,NT1X)*HSP_FRC(1,NT1X) +
     &      SL_FRC(1,NT1X)*RHOL_FRC(1,NT1X)*HL_FRC(1,NT1X) +
     &      SG_FRC(1,NT1X)*RHOG_FRC(1,NT1X)*UEG_FRC(1,NT1X) +
     &      DTHCP_FRC(NFX)*AF_FRC(NT1X)*THCP_FRC(NT1X)*T_FRC(1,NT1X)
          DO M = 1,ISVC+1
            MP = M + 1
            STT0 = SS_FRC(MP,NT1X)*RHOSP_FRC(MP,NT1X)*HSP_FRC(MP,NT1X) +
     &        SL_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)*HL_FRC(MP,NT1X) +
     &        SG_FRC(MP,NT1X)*RHOG_FRC(MP,NT1X)*UEG_FRC(MP,NT1X) +
     &        DTHCP_FRC(NFX)*AF_FRC(NT1X)*THCP_FRC(NT1X)*T_FRC(MP,NT1X)
            STTX(M) = (STT0*APM_FRC(MP,NT1X)-STT1*APM_FRC(1,NT1X))*
     &        DTI*AF_FRC(NT1X)
          ENDDO
!
!---      Loop over fracture triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NTCX = NCX - IPF_FRC(1,NT1X) + 1
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Fracture triangle surface fluxes  ---
!
            DO M = 1,ISVF
              INDX = -3
              M1 = MNOD(M)
              M2 = MADJ(M)
              APX = DIFMN( APM_FRC(M1,NT1X),APM_FRC(M2,NT2X),
     &          DFF1X,DFF2X,ZERO,INDX )
              FT(M,NTCX) = AFF_FRC(NCX)*APX*UFFQ(M,NCX)
            ENDDO
          ENDDO
!
!---      Compute fracture triangle energy equation residuals  ---
!
          RTS = STTX(1) - SRCT_FRC(2,NT1X)
          DO MD = 1,NTCX
            RTS = RTS + FT(1,MD)
          ENDDO
          DO M = 1,ISVC
            RTP(M) = STTX(M+1) - SRCT_FRC(M+2,NT1X)
            MM = 2*M
            DO MD = 1,NTCX
              RTP(M) = RTP(M) + FT(MM,MD)
            ENDDO
          ENDDO
          DO M = 1,ISVC
            MM = 2*M + 1
            DO MD = 1,NTCX
              RTA(M,MD) = RTS - FT(1,MD) + FT(MM,MD)
            ENDDO
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBL_FRC_GT( RTS,RTP,RTA,NT1X,IEQTX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_MFB_GT
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
!     Modify Jacobian matrix for matrix grid cells,
!     fracture triangles, and borehole nodes for transfer of energy
!     between matrix grid cells and fracture triangles and
!     borehole nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/JCBT_MFB_GT'
      IEQTX = IEQT
!
!---  Matrix equations, fracture connections,
!     loop over fractures  ---
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
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fracture triangle  ---
!
          N = INCM_FRC(NCX)
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSQ_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQTX,NMD)
            DO M = 1,ISVC
              IF( ISLC(38).EQ.1 .AND. M.NE.IEQTX ) CYCLE
!
!---          Partial derivative of energy residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSQ_FRC(1,NCX)-TRNSQ_FRC(MX,NCX))/DNR(M,N)
!
!---          Partial derivative of energy residual with
!             respect to fracture triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSQ_FRC(1,NCX)-TRNSQ_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSQ_FRC(1,NCX)
            RSDL(IEQTX,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSQ_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQTX,NMD)
            DO M = 1,ISVC
              IF( ISLC(38).EQ.1 .AND. M.NE.IEQTX ) CYCLE
!
!---          Partial derivative of energy residual with
!             respect to field node primary variables  ---
!
              MCOL = KLU(MP,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSQ_FRC(1,NCX)-TRNSQ_FRC(MX,NCX))/DNR(M,N)
!
!---          Partial derivative of energy residual with
!             respect to fracture triangle primary variables  ---
!
              MC = (NCX-1)*ISVC + IEQTX
              MCOL = KLU_MCF(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSQ_FRC(1,NCX)-TRNSQ_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSQ_FRC(1,NCX)
            RSDL(IEQTX,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO
!
!---  Matrix equations, borehole connections,
!     loop over boreholes  ---
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
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSQ_BH is positive
!           from the borehole node to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQTX,NMD)
            DO M = 1,ISVC
              IF( ISLC(38).EQ.1 .AND. M.NE.IEQTX ) CYCLE
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
!
!---          Partial derivative of energy residual with
!             respect to field node primary variables  ---
!
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSQ_BH(1,NBN)-TRNSQ_BH(MX,NBN))/DNR(M,N)
!
!---          Partial derivative of energy residual with
!             respect to borehole node primary variables  ---
!
              MCOL = IM_BH(M,NBN)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSQ_BH(1,NBN)-TRNSQ_BH(MX,NBN))/DNR_BH(M,NBN)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSQ_BH(1,NBN)
            RSDL(IEQTX,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSQ_BH is positive
!           from the borehole node to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQTX,NMD)
            DO M = 1,ISVC
              IF( ISLC(38).EQ.1 .AND. M.NE.IEQTX ) CYCLE
!
!---          Partial derivative of energy residual with
!             respect to field node primary variables  ---
!
              MCOL = KLU(MP,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSQ_BH(1,NBN)-TRNSQ_BH(MX,NBN))/DNR(M,N)
!
!---          Partial derivative of energy residual with
!             respect to borehole node primary variables  ---
!
              MC = (NBN-1)*ISVC + IEQTX
              MCOL = KLU_MCB(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSQ_BH(1,NBN)-TRNSQ_BH(MX,NBN))/DNR_BH(M,NBN)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSQ_BH(1,NBN)
            RSDL(IEQTX,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
!
!---  Fracture equations, loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fracture triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        MTX = IXP_FRC(NTX)
        IF( MTX.EQ.0 ) CYCLE
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fracture triangle  ---
!
          N = INCM_FRC(NCX)
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fracture triangle, noting that TRNSQ_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM_FRC(IEQTX,NTX)
            DO M = 1,ISVC
              IF( ISLC(38).EQ.1 .AND. M.NE.IEQTX ) CYCLE
!
!---          Partial derivative of energy residual with
!             respect to fracture triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSQ_FRC(MX,NCX)-TRNSQ_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of energy residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSQ_FRC(MX,NCX)-TRNSQ_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSQ_FRC(1,NCX)
            RSDL_FRC(IEQTX,NTX) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fracture triangle, noting that TRNSQ_FRC is positive
!           from the fracture triangle to the field node  ---
!
            MP = IM_FRC(IEQTX,NTX)
            DO M = 1,ISVC
              IF( ISLC(38).EQ.1 .AND. M.NE.IEQTX ) CYCLE
!
!---          Partial derivative of energy residual with
!             respect to fracture triangle primary variables  ---
!
              NMD = (MTX-1)*ISVC + IEQTX
              MCOL = KLU_FRC(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSQ_FRC(MX,NCX)-TRNSQ_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of energy residual with
!             respect to field node primary variables  ---
!
              NMD = (NCX-1)*ISVC + IEQTX
              MCOL = KLU_FCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSQ_FRC(MX,NCX)-TRNSQ_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSQ_FRC(1,NCX)
            RSDL_FRC(IEQTX,NTX) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO
!
!---  Borehole equations, loop over boreholes  ---
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
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Borehole node, noting that TRNSQ_BH is positive
!           from the borehole node to the field node  ---
!
            NMD = IXP(N)
            MP = IM_BH(IEQTX,NBN)
            DO M = 1,ISVC
              IF( ISLC(38).EQ.1 .AND. M.NE.IEQTX ) CYCLE
!
!---          Partial derivative of energy residual with
!             respect to borehole primary variables  ---
!
              MCOL = IM_BH(M,NBN)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSQ_BH(MX,NBN)-TRNSQ_BH(1,NBN))/DNR_BH(M,NBN)
!
!---          Partial derivative of energy residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSQ_BH(MX,NBN)-TRNSQ_BH(1,NBN))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSQ_BH(1,NBN)
            RSDL_BH(IEQTX,NBN) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Borehole node, noting that TRNSQ_BH is positive
!           from the borehole node to the field node  ---
!
            MP = IM_BH(IEQTX,NBN)
            DO M = 1,ISVC
              IF( ISLC(38).EQ.1 .AND. M.NE.IEQTX ) CYCLE
!
!---          Partial derivative of energy residual with
!             respect to borehole primary variables  ---
!
              NMD = (NBN-1)*ISVC + IEQTX
              MCOL = KLU_BH(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSQ_BH(MX,NBN)-TRNSQ_BH(1,NBN))/DNR_BH(M,NBN)
!
!---          Partial derivative of energy residual with
!             respect to field node primary variables  ---
!
              NMD = (NBN-1)*ISVC + IEQTX
              MCOL = KLU_BCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSQ_BH(MX,NBN)-TRNSQ_BH(1,NBN))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSQ_BH(1,NBN)
            RSDL_BH(IEQTX,NBN) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_MFB_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_BF_GT
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
!     Modify the Jacobian matrix for the borehole and fracture water
!     equations with aqueous-phase and gas-phase contributions, for
!     borehole to fracture connections.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 8 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RWP_FRC(LUK),RWA_FRC(LUK),RWP_BH(LUK),RWA_BH(LUK),FW(LSFV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBW_BF_GT'
      IEQWX = IEQW
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
!---      Borehole node to fracture triangle fluxes  ---
!
          DO M = 1,ISVF
            M1 = MNOD(M)
            M2 = MADJ(M)
            FLW1 = XLW_BH(M1,NBN)*RHOL_BH(M1,NBN)
            FLW2 = XLW_FRC(M2,NTX)*RHOL_FRC(M2,NTX)
            INDX = 2
            FLW = DIFMN( FLW1,FLW2,DBF1X,DBF2X,UBFL(1,NCX),INDX )
            FGW1 = XGW_BH(M1,NBN)*RHOG_BH(M1,NBN)
            FGW2 = XGW_FRC(M2,NTX)*RHOG_FRC(M2,NTX)
            INDX = 3
            FGW = DIFMN( FGW1,FGW2,DBF1X,DBF2X,UBFG(1,NCX),INDX )
            APX = GPI*2.D+1*PAR_BH(2,INV_BH(NBN))*APM_FRC(M2,NTX)
            FW(M) = APX*(UBFL(M,NCX)*FLW + UBFG(M,NCX)*FGW
     &          + WTMW*UBFDGW(M,NCX) - WTMW*UBFDLA(M,NCX)
     &          - WTMW*UBFDS(M,NCX))
          ENDDO
!
!---      Energy equation residuals, borehole perspective  ---
!
          RWS_BH = FW(1)
          DO M = 1,ISVC
            MM = 2*M
            RWP_BH(M) = FW(MM)
            MM = 2*M + 1
            RWA_BH(M) = FW(MM)
          ENDDO
!
!---      Energy equation residuals, fracture perspective  ---
!
          RWS_FRC = -FW(1)
          DO M = 1,ISVC
            MM = 2*M + 1
            RWP_FRC(M) = -FW(MM)
            MM = 2*M
            RWA_FRC(M) = -FW(MM)
          ENDDO
!
!---      Modify Jacobian Matrix  ---
!
          CALL JCBL_BF_GT( RWS_BH,RWP_BH,RWA_BH,
     &      RWS_FRC,RWP_FRC,RWA_FRC,NBN,NTX,NCX,IEQWX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_BH_GT
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
!     Load the Jacobian matrix for the borehole water equation with
!     aqueous-phase and gas-phase contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 22 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STWX(LUK+1),RWP(LUK),RWA(LUK,3),FW(LSFV,3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBW_BH_GT'
      IEQWX = IEQW
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
          INV1 = INV_BH(NBN1)
!
!---      Coaxial pipe  ---
!
          IF( IS_BH(INV1).EQ.10000 ) THEN
!
!---        Inner pipe hydraulic area  ---
!
            IF( IBN_BH(NBN1).EQ.0 ) THEN
              AP1X = GPI*(PAR_BH(1,INV1)**2)
!
!---        Outer pipe hydraulic area  ---
!
            ELSE
              AP1X = GPI*((PAR_BH(3,INV1)**2)-(PAR_BH(2,INV1)**2))
            ENDIF
!
!---      Pipe hydraulic area  ---
!
          ELSEIF( IS_BH(INV1).GT.10 ) THEN
            AP1X = GPI*(PAR_BH(3,INV1)**2)
!
!---      Porous media hydraulic area  ---
!
          ELSE
            IZN = IZ_BH(NBN1)
            AP1X = GPI*(PAR_BH(2,INV1)**2)
          ENDIF
!
!---      First-order, forward-difference, time differential  ---
!
          STW1 = PORD_BH(1,NBN1)*
     &      (XLW_BH(1,NBN1)*RHOL_BH(1,NBN1)*SL_BH(1,NBN1) +
     &        XGW_BH(1,NBN1)*RHOG_BH(1,NBN1)*SG_BH(1,NBN1))
          DO M = 1,ISVC+1
            MP = M + 1
            STW0 = PORD_BH(MP,NBN1)*
     &      (XLW_BH(MP,NBN1)*RHOL_BH(MP,NBN1)*SL_BH(MP,NBN1) +
     &        XGW_BH(MP,NBN1)*RHOG_BH(MP,NBN1)*SG_BH(MP,NBN1))
            STWX(M) = (STW0-STW1)*DTI*VOL_BH(NBN1)
          ENDDO
!
!---      Loop over borehole node to node connections ---
!
          NBCX = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBCX = NBCX + 1
            NBN2 = IBCM_BH(NCX)
            INV2 = INV_BH(NBN2)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
!
!---        Coaxial pipe  ---
!
            IF( IS_BH(INV2).EQ.10000 ) THEN
!
!---          Inner pipe hydraulic area  ---
!
              IF( IBN_BH(NBN2).EQ.0 ) THEN
                AP2X = GPI*(PAR_BH(1,INV2)**2)
!
!---          Outer pipe hydraulic area  ---
!
              ELSE
                AP2X = GPI*((PAR_BH(3,INV2)**2)-(PAR_BH(2,INV2)**2))
              ENDIF
!
!---        Pipe hydraulic area  ---
!
            ELSEIF( IS_BH(INV2).GT.10 ) THEN
              AP2X = GPI*(PAR_BH(3,INV2)**2)
!
!---        Porous media hydraulic area  ---
!
            ELSE
              IZN = IZ_BH(NBN2)
              AP2X = GPI*(PAR_BH(2,INV2)**2)
            ENDIF
            INDX = -3
            APX = DIFMN( AP1X,AP2X,DBB1X,DBB2X,ZERO,INDX )
!
!---        Borehole node surface fluxes  ---
!
            DO M = 1,ISVF
              M1 = MNOD(M)
              M2 = MADJ(M)
              FLW1 = XLW_BH(M1,NBN1)*RHOL_BH(M1,NBN1)
              FLW2 = XLW_BH(M2,NBN2)*RHOL_BH(M2,NBN2)
              INDX = 2
              FLW = DIFMN( FLW1,FLW2,DBB1X,DBB2X,UBBL(1,NCX),INDX )
              FGW1 = XGW_BH(M1,NBN1)*RHOG_BH(M1,NBN1)
              FGW2 = XGW_BH(M2,NBN2)*RHOG_BH(M2,NBN2)
              INDX = 3
              FGW = DIFMN( FGW1,FGW2,DBB1X,DBB2X,UBBG(1,NCX),INDX )
              FW(M,NBCX) = APX*
     &          (UBBL(M,NCX)*FLW + UBBG(M,NCX)*FGW
     &          + WTMW*UBBDGW(M,NCX) - WTMW*UBBDLA(M,NCX)
     &          - WTMW*UBBDS(M,NCX))
            ENDDO
          ENDDO
!
!---      Compute borehole node water equation residuals  ---
!
          RWS = STWX(1) - SRCW_BH(2,NBN1)
          DO MD = 1,NBCX
            RWS = RWS + FW(1,MD)
          ENDDO
          DO M = 1,ISVC
            RWP(M) = STWX(M+1) - SRCW_BH(M+2,NBN1)
            MM = 2*M
            DO MD = 1,NBCX
              RWP(M) = RWP(M) + FW(MM,MD)
            ENDDO
          ENDDO
          DO M = 1,ISVC
            MM = 2*M + 1
            DO MD = 1,NBCX
              RWA(M,MD) = RWS - FW(1,MD) + FW(MM,MD)
            ENDDO
            RWA(M,3) = RWS
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBL_BH_GT( RWS,RWP,RWA,NBN1,IEQWX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_FRC_GT
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
!     Load the Jacobian matrix for the fracture water equation with
!     aqueous-phase and gas-phase contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 14 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_FRC
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STWX(LUK+1),RWP(LUK),RWA(LUK,LTC_FRC),FW(LSFV,LTC_FRC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBW_FRC_GT'
      IEQWX = IEQW
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
!---      First-order, forward-difference, time differential  ---
!
          STW1 = (XLW_FRC(1,NT1X)*RHOL_FRC(1,NT1X)*SL_FRC(1,NT1X) +
     &        XGW_FRC(1,NT1X)*RHOG_FRC(1,NT1X)*SG_FRC(1,NT1X))
          DO M = 1,ISVC+1
            MP = M + 1
            STW0 = (XLW_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)*SL_FRC(MP,NT1X) +
     &        XGW_FRC(MP,NT1X)*RHOG_FRC(MP,NT1X)*SG_FRC(MP,NT1X))
            STWX(M) = (STW0*APM_FRC(MP,NT1X)-STW1*APM_FRC(1,NT1X))*
     &        DTI*AF_FRC(NT1X)
          ENDDO
!
!---      Loop over fracture triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NTCX = NCX - IPF_FRC(1,NT1X) + 1
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Fracture triangle surface fluxes  ---
!
            DO M = 1,ISVF
              M1 = MNOD(M)
              M2 = MADJ(M)
              FLW1 = XLW_FRC(M1,NT1X)*RHOL_FRC(M1,NT1X)
              FLW2 = XLW_FRC(M2,NT2X)*RHOL_FRC(M2,NT2X)
              INDX = 2
              FLW = DIFMN( FLW1,FLW2,DFF1X,DFF2X,UFFL(1,NCX),INDX )
              FGW1 = XGW_FRC(M1,NT1X)*RHOG_FRC(M1,NT1X)
              FGW2 = XGW_FRC(M2,NT2X)*RHOG_FRC(M2,NT2X)
              INDX = 3
              FGW = DIFMN( FGW1,FGW2,DFF1X,DFF2X,UFFG(1,NCX),INDX )
              INDX = -3
              APX = DIFMN( APM_FRC(M1,NT1X),APM_FRC(M2,NT2X),
     &          DFF1X,DFF2X,ZERO,INDX )
              FW(M,NTCX) = AFF_FRC(NCX)*APX*
     &          (UFFL(M,NCX)*FLW + UFFG(M,NCX)*FGW
     &          + WTMW*UFFDGW(M,NCX) - WTMW*UFFDLA(M,NCX)
     &          - WTMW*UFFDS(M,NCX))
            ENDDO
          ENDDO
!
!---      Compute fracture triangle water equation residuals  ---
!
          RWS = STWX(1) - SRCW_FRC(2,NT1X)
          DO MD = 1,NTCX
            RWS = RWS + FW(1,MD)
          ENDDO
          DO M = 1,ISVC
            RWP(M) = STWX(M+1) - SRCW_FRC(M+2,NT1X)
            MM = 2*M
            DO MD = 1,NTCX
              RWP(M) = RWP(M) + FW(MM,MD)
            ENDDO
          ENDDO
          DO M = 1,ISVC
            MM = 2*M + 1
            DO MD = 1,NTCX
              RWA(M,MD) = RWS - FW(1,MD) + FW(MM,MD)
            ENDDO
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBL_FRC_GT( RWS,RWP,RWA,NT1X,IEQWX )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_MFB_GT
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
!     Modify Jacobian matrix for matrix grid cells,
!     fracture triangles, and borehole nodes for transfer of water
!     mass between matrix grid cells and fracture triangles and
!     borehole nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
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
      SUB_LOG(ISUB_LOG) = '/JCBW_MFB_GT'
      IEQWX = IEQW
!
!---  Matrix equations, fracture connections,
!     loop over fractures  ---
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
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fracture triangle  ---
!
          N = INCM_FRC(NCX)
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSW_FRC is positive
!           from the fracture triangle to the field node, and the
!           field node is considered the upper node for flux
!           indexing  ---
!
            NMD = IXP(N)
            MP = IM(IEQWX,NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSW_FRC(1,NCX)-TRNSW_FRC(MX,NCX))/DNR(M,N)
!
!---          Partial derivative of water mass residual with
!             respect to fracture triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSW_FRC(1,NCX)-TRNSW_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSW_FRC(1,NCX)
            RSDL(IEQWX,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSW_FRC is positive
!           from the fracture triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQWX,NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with
!             respect to field node primary variables  ---
!
              MCOL = KLU(MP,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSW_FRC(1,NCX)-TRNSW_FRC(MX,NCX))/DNR(M,N)
!
!---          Partial derivative of water mass residual with
!             respect to fracture triangle primary variables  ---
!
              MC = (NCX-1)*ISVC + IEQWX
              MCOL = KLU_MCF(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSW_FRC(1,NCX)-TRNSW_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSW_FRC(1,NCX)
            RSDL(IEQWX,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO
!
!---  Matrix equations, borehole connections,
!     loop over boreholes  ---
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
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSW_BH is positive
!           from the borehole node to the field node, and the
!           field node is considered the upper node for flux
!           indexing  ---
!
            NMD = IXP(N)
            MP = IM(IEQWX,NMD)
            DO M = 1,ISVC
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
!
!---          Partial derivative of water mass residual with
!             respect to field node primary variables  ---
!
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSW_BH(1,NBN)-TRNSW_BH(MX,NBN))/DNR(M,N)
!
!---          Partial derivative of water mass residual with
!             respect to borehole node primary variables  ---
!
              MCOL = IM_BH(M,NBN)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSW_BH(1,NBN)-TRNSW_BH(MX,NBN))/DNR_BH(M,NBN)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSW_BH(1,NBN)
            RSDL(IEQWX,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSW_BH is positive
!           from the borehole node to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQWX,NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with
!             respect to field node primary variables  ---
!
              MCOL = KLU(MP,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSW_BH(1,NBN)-TRNSW_BH(MX,NBN))/DNR(M,N)
!
!---          Partial derivative of water mass residual with
!             respect to borehole node primary variables  ---
!
              MC = (NBN-1)*ISVC + IEQWX
              MCOL = KLU_MCB(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSW_BH(1,NBN)-TRNSW_BH(MX,NBN))/DNR_BH(M,NBN)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSW_BH(1,NBN)
            RSDL(IEQWX,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
!
!---  Fracture equations, loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fracture triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        MTX = IXP_FRC(NTX)
        IF( MTX.EQ.0 ) CYCLE
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fracture triangle  ---
!
          N = INCM_FRC(NCX)
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fracture triangle, noting that TRNSW_FRC is positive
!           from the fracture triangle to the field node, and the
!           field node is considered the upper node for flux
!           indexing  ---
!
            NMD = IXP(N)
            MP = IM_FRC(IEQWX,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with
!             respect to fracture triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSW_FRC(MX,NCX)-TRNSW_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of water mass residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSW_FRC(MX,NCX)-TRNSW_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSW_FRC(1,NCX)
            RSDL_FRC(IEQWX,NTX) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fracture triangle, noting that TRNSW_FRC is positive
!           from the fracture triangle to the field node  ---
!
            MP = IM_FRC(IEQWX,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with
!             respect to fracture triangle primary variables  ---
!
              NMD = (MTX-1)*ISVC + IEQWX
              MCOL = KLU_FRC(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSW_FRC(MX,NCX)-TRNSW_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of water mass residual with
!             respect to field node primary variables  ---
!
              NMD = (NCX-1)*ISVC + IEQWX
              MCOL = KLU_FCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSW_FRC(MX,NCX)-TRNSW_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSW_FRC(1,NCX)
            RSDL_FRC(IEQWX,NTX) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO
!
!---  Borehole equations, loop over boreholes  ---
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
          IF( N.EQ.0 ) CYCLE
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Borehole node, noting that TRNSW_BH is positive
!           from the borehole node to the field node, and the
!           field node is considered the upper node for flux
!           indexing  ---
!
            NMD = IXP(N)
            MP = IM_BH(IEQWX,NBN)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with
!             respect to borehole primary variables  ---
!
              MCOL = IM_BH(M,NBN)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSW_BH(MX,NBN)-TRNSW_BH(1,NBN))/DNR_BH(M,NBN)
!
!---          Partial derivative of water mass residual with
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) +
     &          (TRNSW_BH(MX,NBN)-TRNSW_BH(1,NBN))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSW_BH(1,NBN)
            RSDL_BH(IEQWX,NBN) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Borehole node, noting that TRNSW_BH is positive
!           from the borehole node to the field node  ---
!
            MP = IM_BH(IEQWX,NBN)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with
!             respect to borehole primary variables  ---
!
              NMD = (NBN-1)*ISVC + IEQWX
              MCOL = KLU_BH(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSW_BH(MX,NBN)-TRNSW_BH(1,NBN))/DNR_BH(M,NBN)
!
!---          Partial derivative of water mass residual with
!             respect to field node primary variables  ---
!
              NMD = (NBN-1)*ISVC + IEQWX
              MCOL = KLU_BCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) +
     &          (TRNSW_BH(MX,NBN)-TRNSW_BH(1,NBN))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSW_BH(1,NBN)
            RSDL_BH(IEQWX,NBN) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_MFB_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_BH_GT
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
!     Load the current time step values into the old time step
!     variables for boreholes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 22 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PARM_BH
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
      SUB_LOG(ISUB_LOG) = '/LDO_BH_GT'
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
          T_BH(1,NBN) = T_BH(2,NBN)
          PL_BH(1,NBN) = PL_BH(2,NBN)
          PG_BH(1,NBN) = PG_BH(2,NBN)
          PORD_BH(1,NBN) = PORD_BH(2,NBN)
          SL_BH(1,NBN) = SL_BH(2,NBN)
          SG_BH(1,NBN) = SG_BH(2,NBN)
          PVA_BH(1,NBN) = PVA_BH(2,NBN)
          PVW_BH(1,NBN) = PVW_BH(2,NBN)
          XGA_BH(1,NBN) = XGA_BH(2,NBN)
          XGW_BH(1,NBN) = XGW_BH(2,NBN)
          XMGA_BH(1,NBN) = XMGA_BH(2,NBN)
          XMGW_BH(1,NBN) = XMGW_BH(2,NBN)
          XLA_BH(1,NBN) = XLA_BH(2,NBN)
          XLW_BH(1,NBN) = XLW_BH(2,NBN)
          XMLA_BH(1,NBN) = XMLA_BH(2,NBN)
          XMLW_BH(1,NBN) = XMLW_BH(2,NBN)
          RHOL_BH(1,NBN) = RHOL_BH(2,NBN)
          RHOG_BH(1,NBN) = RHOG_BH(2,NBN)
          RHOSP_BH(1,NBN) = RHOSP_BH(2,NBN)
          VISL_BH(1,NBN) = VISL_BH(2,NBN)
          VISG_BH(1,NBN) = VISG_BH(2,NBN)
          RKL_BH(1,NBN) = RKL_BH(2,NBN)
          RKG_BH(1,NBN) = RKG_BH(2,NBN)
          DFGW_BH(1,NBN) = DFGW_BH(2,NBN)
          DFLA_BH(1,NBN) = DFLA_BH(2,NBN)
          THKG_BH(1,NBN) = THKG_BH(2,NBN)
          THKL_BH(1,NBN) = THKL_BH(2,NBN)
          HL_BH(1,NBN) = HL_BH(2,NBN)
          HGW_BH(1,NBN) = HGW_BH(2,NBN)
          HGA_BH(1,NBN) = HGA_BH(2,NBN)
          HG_BH(1,NBN) = HG_BH(2,NBN)
          UEG_BH(1,NBN) = UEG_BH(2,NBN)
          HSP_BH(1,NBN) = HSP_BH(2,NBN)
          NPHAZ_BH(1,NBN) = NPHAZ_BH(2,NBN)
          SS_BH(1,NBN) = SS_BH(2,NBN)
          TMS_BH(1,NBN) = TMS_BH(2,NBN)
          XLS_BH(1,NBN) = XLS_BH(2,NBN)
          YLS_BH(1,NBN) = YLS_BH(2,NBN)
          DO NSL = 1,NSOLU
            CO_BH(NBN,NSL) = C_BH(NBN,NSL)
          ENDDO

          DO NEQ = 1,NEQC+NEQK
            NSL = NEQ + NSOLU
            CO_BH(NBN,NSL) = C_BH(NBN,NSL)
          ENDDO
          DO NSP = 1,NSPR
            SP_CO_BH(NBN,NSP) = SP_C_BH(NBN,NSP)
          ENDDO

        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_FRC_GT
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
!     Load the current time step values into the old time step
!     variables for fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 8 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE GEOM_FRC
      USE FDVT_FRC
      USE FDVS_FRC
      USE FDVP_FRC
      USE FDVG_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LDO_FRC_GT'
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
          APH_FRC(1,NTX) = APH_FRC(2,NTX)
          APM_FRC(1,NTX) = APM_FRC(2,NTX)
          T_FRC(1,NTX) = T_FRC(2,NTX)
          PL_FRC(1,NTX) = PL_FRC(2,NTX)
          PG_FRC(1,NTX) = PG_FRC(2,NTX)
          SL_FRC(1,NTX) = SL_FRC(2,NTX)
          SG_FRC(1,NTX) = SG_FRC(2,NTX)
          PVA_FRC(1,NTX) = PVA_FRC(2,NTX)
          PVW_FRC(1,NTX) = PVW_FRC(2,NTX)
          XGA_FRC(1,NTX) = XGA_FRC(2,NTX)
          XGW_FRC(1,NTX) = XGW_FRC(2,NTX)
          XMGA_FRC(1,NTX) = XMGA_FRC(2,NTX)
          XMGW_FRC(1,NTX) = XMGW_FRC(2,NTX)
          XLA_FRC(1,NTX) = XLA_FRC(2,NTX)
          XLW_FRC(1,NTX) = XLW_FRC(2,NTX)
          XMLA_FRC(1,NTX) = XMLA_FRC(2,NTX)
          XMLW_FRC(1,NTX) = XMLW_FRC(2,NTX)
          RHOL_FRC(1,NTX) = RHOL_FRC(2,NTX)
          RHOG_FRC(1,NTX) = RHOG_FRC(2,NTX)
          RHOSP_FRC(1,NTX) = RHOSP_FRC(2,NTX)
          VISL_FRC(1,NTX) = VISL_FRC(2,NTX)
          VISG_FRC(1,NTX) = VISG_FRC(2,NTX)
          RKL_FRC(1,NTX) = RKL_FRC(2,NTX)
          RKG_FRC(1,NTX) = RKG_FRC(2,NTX)
          DFGW_FRC(1,NTX) = DFGW_FRC(2,NTX)
          DFLA_FRC(1,NTX) = DFLA_FRC(2,NTX)
          THKG_FRC(1,NTX) = THKG_FRC(2,NTX)
          THKL_FRC(1,NTX) = THKL_FRC(2,NTX)
          HL_FRC(1,NTX) = HL_FRC(2,NTX)
          HGW_FRC(1,NTX) = HGW_FRC(2,NTX)
          HGA_FRC(1,NTX) = HGA_FRC(2,NTX)
          HG_FRC(1,NTX) = HG_FRC(2,NTX)
          UEG_FRC(1,NTX) = UEG_FRC(2,NTX)
          HSP_FRC(1,NTX) = HSP_FRC(2,NTX)
          NPHAZ_FRC(1,NTX) = NPHAZ_FRC(2,NTX)
          SS_FRC(1,NTX) = SS_FRC(2,NTX)
          TMS_FRC(1,NTX) = TMS_FRC(2,NTX)
          XLS_FRC(1,NTX) = XLS_FRC(2,NTX)
          YLS_FRC(1,NTX) = YLS_FRC(2,NTX)
          DO NSL = 1,NSOLU
            CO_FRC(NTX,NSL) = C_FRC(NTX,NSL)
          ENDDO

          DO NEQ = 1,NEQC+NEQK
            NSL = NEQ + NSOLU
            CO_FRC(NTX,NSL) = C_FRC(NTX,NSL)
          ENDDO
          DO NSP = 1,NSPR
            SP_CO_FRC(NTX,NSP) = SP_C_FRC(NTX,NSP)
          ENDDO

        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PIPE_FLOW_VEL_GT( DBBX,DHX,HDX,RFX,RHOX,UX,VISX,NBN )
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
!     Compute fluid velocity in a pipe
!
!     DBBX - pipe length, m
!     DHX - hydraulic diameter, m
!     HDX - pressure drop, Pa
!     RFX - roughness, m
!     RHOX - fluid density, kg/m^3
!     UX - fluid velocity, m/s
!     VISX - fluid viscosity, Pa s
!     NBN - borehole number
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 27 April, 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 GX(2)
      SAVE NCMAX
      DATA NCMAX / 0 /
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PIPE_FLOW_VEL_GT'
!
!---  Hydraulic area (m^2) ---
!
      AHX = GPI*((5.D-1*DHX)**2)
!!
!!---  Hydraulic diameter (inches), sum of interval radii ---
!!
!      DHIX = 39.370079D+0*DHX
!!
!!---  Hazen-Williams coefficient for steel, welded and seamless  ---
!!
!      CHWX = 1.D+2
!!
!!---  Estimate absolute flow rate from Hazen-Williams
!!     equation, adjusting for fluid viscosity (gal/min)  ---
!!
!      DHDX = HDX/RHORL/GRAV/DBBX
!      DH100X = 1.D+2*DHDX
!      QX = (ABS(DH100X)*(DHIX**4.8655D+0)*(VISRL/VISX)
!     &  /(0.2083D+0*(1.D+2/CHWX)**1.852))**5.3996D-1
!!
!!---  Convert to absolute flow velocity (m/s) ---
!!
!      UX = 6.30833D-5*QX/AHX
!
!---  Hazen-Williams coefficient for steel, welded and seamless  ---
!
      CHWX = 1.D+2
!
!---  Head drop, m ---
!
      HDDX = ABS(HDX/RHOX/GRAV)
!
!---  Volumetric flow rate from Hazen-Williams equation, m^3/s ---
!
      QX = (HDDX/DBBX)*(CHWX**1.852D+0)*(DHX**4.8704D+0)/10.67D+0
      QX = QX**(1.D+0/1.852D+0)
!
!---  Flow velocity estimate, m/s ---
!
      UX = QX/AHX
!
!---  Newton-Raphson Solution for absolute flow velocity,
!     combining the Colebrook equation for the friction
!     coefficient and the Darcy-Weisbach equation for the
!     pressure drop (i.e., head loss)  ---
!
      IF( UX.GT.1.D-9 ) THEN
        DO NC = 1,32
          DO MC = 1,2
            DUZ = MAX( UX*1.D-6,1.D-12 )
            UZ = UX
            IF( MC.EQ.2 ) UZ = UX + DUZ
!
!---        Darcy-Weisbach friction coefficient ---
!
            FCX = 2.D+0*ABS(HDX)*DHX/((UZ**2)*RHOX*DBBX)
!
!---        Reynolds number ---
!
            REYX = RHOX*UZ*DHX/VISX
!
!---        Colebrook friction coefficient equation ---
!
            GX(MC) = 1.D+0/SQRT(FCX) + 2.D+0*LOG10(
     &         (2.51D+0/(REYX*SQRT(FCX))) + RFX/(DHX*3.72D+0))
          ENDDO
          FX = GX(1)
          DFX = (GX(2)-GX(1))/DUZ
          DUX = -FX/DFX
          IF( DUX.LT.0.D+0 .AND. UX.LT.1.D-9 ) THEN
            UX = 0.D+0
            EXIT
          ENDIF
          UX = MAX( UX+DUX,1.D-10 )
          IF( ABS(DUX).LT.MAX(1.D-12,UX*1.D-6) ) EXIT
        ENDDO
        IF( NC.GT.32 ) THEN
!
!---      Hydraulic area (m^2) ---
!
          AHX = GPI*((5.D-1*DHX)**2)
!
!---      Hydraulic diameter (inches), sum of interval radii ---
!
          DHIX = 39.370079D+0*DHX
!
!---      Hazen-Williams coefficient for steel, welded and seamless  ---
!
          CHWX = 1.D+2
!
!---      Estimate absolute flow rate from Hazen-Williams
!         equation, adjusting for fluid viscosity (gal/min)  ---
!
          DHDX = HDX/RHORL/GRAV/DBBX
          DH100X = 1.D+2*DHDX
          QX = (ABS(DH100X)*(DHIX**4.8655D+0)*(VISRL/VISX)
     &      /(0.2083D+0*(1.D+2/CHWX)**1.852))**5.3996D-1
!
!---      Convert to absolute flow velocity (m/s) ---
!
          UX = 6.30833D-5*QX/AHX
!
!---      Newton-Raphson Solution for absolute flow velocity,
!         combining the Colebrook equation for the friction
!         coefficient and the Darcy-Weisbach equation for the
!         pressure drop (i.e., head loss)  ---
!
          IF( UX.GT.1.D-9 ) THEN
            print *,'DBBX = ',DBBX
            print *,'DHX = ',DHX
            print *,'HDX = ',HDX
            print *,'RFX = ',RFX
            print *,'RHOX = ',RHOX
            print *,'UX = ',UX
            print *,'VISX = ',VISX
            print *,'NBN = ',NBN
            DO NC = 1,32
              DO MC = 1,2
                DUZ = MAX( UX*1.D-6,1.D-12 )
                UZ = UX
                IF( MC.EQ.2 ) UZ = UX + DUZ
!
!---            Darcy-Weisbach friction coefficient ---
!
                FCX = 2.D+0*ABS(HDX)*DHX/((UZ**2)*RHOX*DBBX)
!
!---            Reynolds number ---
!
                REYX = RHOX*UZ*DHX/VISX
!
!---            Colebrook friction coefficient equation ---
!
                GX(MC) = 1.D+0/SQRT(FCX) + 2.D+0*LOG10(
     &             (2.51D+0/(REYX*SQRT(FCX))) + RFX/(DHX*3.72D+0))
              ENDDO
              print *,'DUZ = ',DUZ
              print *,'FCX = ',FCX
              print *,'REYX = ',REYX
              print *,'GX = ',GX
              print *,'UX = ',UX
              FX = GX(1)
              DFX = (GX(2)-GX(1))/DUZ
              DUX = -FX/DFX
              IF( DUX.LT.0.D+0 .AND. UX.LT.1.D-9 ) THEN
                UX = 0.D+0
                EXIT
              ENDIF
              UX = MAX( UX+DUX,1.D-10 )
              IF( ABS(DUX).LT.1.D-12 ) EXIT
            ENDDO
          ENDIF
          INDX = 12
          IMSGX = NBN
          NMSGX = 0
          SUBLOGX = 'PIPE_FLOW_VEL_GT'
          RLMSGX = 0.D+0
          CHMSGX = 'Pipe Flow Convergence Failure: ' //
     &      'Newton-Raphson @ Borehole Node: '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          print *,'DBBX = ',DBBX
          print *,'DHX = ',DHX
          print *,'HDX = ',HDX
          print *,'RFX = ',RFX
          print *,'RHOX = ',RHOX
          print *,'UX = ',UX
          print *,'VISX = ',VISX
          print *,'NBN = ',NBN
          VISX = -1.D+0
        ENDIF
      ELSE
        ULX = 0.D+0
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PIPE_FLOW_VEL_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PIPE_FLOW_DP_GT( DBBX,DHX,HDX,RFX,RHOX,UX,VISX,NBN )
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
!     Compute pressure drop in a pipe
!
!     DBBX - pipe length, m
!     DHX - hydraulic diameter, m
!     HDX - pressure drop, Pa
!     RFX - roughness, m
!     RHOX - fluid density, kg/m^3
!     UX - fluid velocity, m/s
!     VISX - fluid viscosity, Pa s
!     NBN - borehole number
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 26 April 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      SUB_LOG(ISUB_LOG) = '/PIPE_FLOW_DP_GT'
!
!---  Hydraulic area (m^2) ---
!
      AHX = GPI*((5.D-1*DHX)**2)
!
!---  Reynolds number ---
!
      REYX = RHOX*UX*DHX/VISX
!
!---  Haaland equation for Darcy-Weisbach friction factor ---
!
      FCX = -1.8D+0*LOG10( ((RFX/DHX/3.7D+0)**1.11D+0) + 6.9D+0/REYX )
      FCX = (1.D+0/FCX)**2
!
!---  Darcy-Weisbach pressure loss ---
!
      HDX = 5.D-1*DBBX*FCX*RHOX*(UX**2)/DHX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PIPE_FLOW_DP_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_BH_GT
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
!     Compute hydrologic, thermodynamic and physical properties for
!     boreholes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 22 April 2019.
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DFEFX(5)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP_BH_GT'
      NC = 0
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
          IF( NBN.EQ.ID_BH(3,NBH) .OR. NBN.EQ.ID_BH(4,NBH) ) THEN
            NC = NC + 1
          ELSE
            NC = NC + 2
          ENDIF

          POR0_BH(1,NBN) = POR0_BH(1,NBN)
          POR0_BH(2,NBN) = POR0_BH(2,NBN)

          IZN = IZ_BH(NBN)
          INV = INV_BH(NBN)
          N_DB = NBN
!
!---      Pipe flow  ---
!
          IF( IS_BH(INV).GT.10 ) THEN
            ENPR = 0.D+0
          ELSE
!
!---      Assign gas entry pressure and minimum gas saturation
!         for transition to unsaturated conditions  ---
!
            ENPR = 0.D+0
            IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
              ENPR = SCHR(1,IZN)*RHORL*GRAV
            ELSEIF( ISCHR(IZN).EQ.4 ) THEN
              ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
            ENDIF
          ENDIF
!
!---      Loop over increment indices  ---
!
          DO M = 2,ISVC+2
            PLX = PL_BH(M,NBN)+PATM
            PGX = PG_BH(M,NBN)+PATM
            INDX = 0
            CALL REGION_4( T_BH(M,NBN),PSW_BH(M,NBN),INDX )
!
!---        Saturated system w/o entrapped gas
!           Energy - temperature
!           Water mass - aqueous pressure
!           Air mass - aqueous-air mole fraction
!           NaCl mass - total NaCl brine mass fraction  ---
!
            IF( NPHAZ_BH(2,NBN).EQ.1 ) THEN
              PX = PLX
              PWX = PLX
              PCX = 0.D+0
              SL_BH(M,NBN) = 1.D+0
              CALL SOL_BRNS( T_BH(M,NBN),PWX,XLSMX )
              XLS_BH(M,NBN) = MIN(YLS_BH(M,NBN),XLSMX)
              CALL SP_B( T_BH(M,NBN),XLS_BH(M,NBN),PSBX )
              CALL P_IAPWS( T_BH(M,NBN),PWX,RHOGWX,RHOLWX,HGWX,HLWX,
     &          UGWX,ULWX )
              if( rholwx.lt.0.d+0 ) then
                print *,'nphaz = 1'
                print *,'t_bh(m,nbn) = ',T_BH(M,NBN)
                print *,'pwx = ',PWX
                stop
              endif
              CALL DENS_B( XLS_BH(M,NBN),RHOBX,RHOLWX,T_BH(M,NBN) )
              PVBX = PSBX
!              PVA_BH(M,NBN) = XMLA_BH(M,NBN)*HCAW
              XMLA_BH(M,NBN) = PVA_BH(M,NBN)/HCAW
              PVW_BH(M,NBN) = PVBX
              XMGA_BH(M,NBN) = PVA_BH(M,NBN)/
     &          (PVA_BH(M,NBN)+PVW_BH(M,NBN))
              CALL EQUIL( XGA_BH(M,NBN),XGW_BH(M,NBN),XLA_BH(M,NBN),
     &          XLS_BH(M,NBN),XLW_BH(M,NBN),XMGA_BH(M,NBN),
     &          XMGW_BH(M,NBN),XMLA_BH(M,NBN),XMLS_BH(M,NBN),
     &          XMLW_BH(M,NBN) )
!
!---        Unsaturated system w/ or w/o entrapped gas
!
!           Energy - temperature
!           Water mass - aqueous pressure
!           Air mass - gas pressure
!           NaCl mass - total NaCl brine mass fraction  ---
!
            ELSEIF( NPHAZ_BH(2,NBN).EQ.2 ) THEN
              PX = PG_BH(M,NBN) + PATM
              CALL P_IAPWS( T_BH(M,NBN),PSW_BH(M,NBN),RHOX,RHOLWX,
     &          HX,HLWX,UX,ULWX )
              if( rholwx.lt.0.d+0 ) then
                print *,'nphaz = 2'
                print *,'t_bh(m,nbn) = ',T_BH(M,NBN)
                print *,'pwx = ',PSW_BH(M,NBN)
                stop
              endif
              PCX = PG_BH(M,NBN)-PL_BH(M,NBN)
              CALL VPL_B( T_BH(M,NBN),PSW_BH(M,NBN),PCX,RHOLWX,PVWX )
              CALL SOL_BRNS( T_BH(M,NBN),PSW_BH(M,NBN),XLSMX )
              XLS_BH(M,NBN) = MIN(YLS_BH(M,NBN),XLSMX)
              CALL DENS_B( XLS_BH(M,NBN),RHOBX,RHOLWX,T_BH(M,NBN) )
              CALL SP_B( T_BH(M,NBN),XLS_BH(M,NBN),PSBX )
              CALL VPL_B( T_BH(M,NBN),PSBX,PCX,RHOBX,PVBX )
              CALL P_IAPWS( T_BH(M,NBN),PVBX,RHOGWX,RHOX,HGWX,HX,
     &          UGWX,UX )
              PVA_BH(M,NBN) = PGX-PVBX
              IF( PVA_BH(M,NBN).LT.1.D-6 ) PVA_BH(M,NBN) = 0.D+0
              PVW_BH(M,NBN) = PVBX
              XMLA_BH(M,NBN) = PVA_BH(M,NBN)/HCAW
              XMGA_BH(M,NBN) = PVA_BH(M,NBN)/PGX
              CALL EQUIL( XGA_BH(M,NBN),XGW_BH(M,NBN),XLA_BH(M,NBN),
     &          XLS_BH(M,NBN),XLW_BH(M,NBN),XMGA_BH(M,NBN),
     &          XMGW_BH(M,NBN),XMLA_BH(M,NBN),XMLS_BH(M,NBN),
     &          XMLW_BH(M,NBN) )
!
!---        Fully unsaturated conditions
!
!           Energy - temperature
!           Water mass - water vapor partial pressure
!           Air mass - gas pressure
!           NaCl mass - salt mass  ---
!
            ELSEIF( NPHAZ_BH(2,NBN).EQ.4 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PX = PGX
              INDX = 0
              CALL REGION_4( T_BH(M,NBN),PSW_BH(M,NBN),INDX )
              CALL SOL_BRNS( T_BH(M,NBN),PSW_BH(M,NBN),XLSMX )
              IF( TMS_BH(M,NBN).GT.EPSL ) THEN
                XLS_BH(M,NBN) = XLSMX
              ELSE
                XLS_BH(M,NBN) = 0.D+0
              ENDIF
              PVA_BH(M,NBN) = PGX - PVW_BH(M,NBN)
              SL_BH(M,NBN) = 0.D+0
              CALL CAP_BH_GT( CPGLX,SL_BH(M,NBN),NBN )
              PL_BH(M,NBN) = PG_BH(M,NBN) - PCX
              PLX = PGX - PCX
              CALL P_IAPWS( T_BH(M,NBN),PVW_BH(M,NBN),RHOGWX,RHOLWX,
     &          HGWX,HLWX,UGWX,ULWX )
              if( rholwx.lt.0.d+0 ) then
                print *,'nphaz = 4'
                print *,'t_bh(m,nbn) = ',T_BH(M,NBN)
                print *,'pwx = ',PVW_BH(M,NBN)
                stop
              endif
              CALL DENS_B( XLS_BH(M,NBN),RHOBX,RHOLWX,T_BH(M,NBN) )
              XMLA_BH(M,NBN) = PVA_BH(M,NBN)/HCAW
              XMGA_BH(M,NBN) = PVA_BH(M,NBN)/PGX
              CALL EQUIL( XGA_BH(M,NBN),XGW_BH(M,NBN),XLA_BH(M,NBN),
     &          XLS_BH(M,NBN),XLW_BH(M,NBN),XMGA_BH(M,NBN),
     &          XMGW_BH(M,NBN),XMLA_BH(M,NBN),XMLS_BH(M,NBN),
     &          XMLW_BH(M,NBN) )
            ENDIF
!
!---        Coaxial borehole, inner pipe  ---
!
            IF( IS_BH(INV).EQ.10000 .AND. IBN_BH(NBN).EQ.0 ) THEN
              PORD_BH(M,NBN) = MIN( 1.D+0,MAX( (PAR_BH(1,INV)**2)/
     &          (PAR_BH(2,INV)**2),EPSL ) )
!
!---        Coaxial borehole, outer pipe  ---
!
            ELSEIF( IS_BH(INV).EQ.10000 .AND. IBN_BH(NBN).NE.0 ) THEN
              PORD_BH(M,NBN) = MIN( 1.D+0,MAX(
     &          ((PAR_BH(3,INV)**2)-(PAR_BH(2,INV)**2))/
     &          ((PAR_BH(4,INV)**2)-(PAR_BH(2,INV)**2)),EPSL ) )
!
!---        Pipe flow  ---
!
            ELSEIF( IS_BH(INV).GT.10 ) THEN
              PORD_BH(M,NBN) = MIN( 1.D+0,MAX( (PAR_BH(3,INV)**2)/
     &          (PAR_BH(2,INV)**2),EPSL ) )
!
!---        Porous-media porosity  ---
!
            ELSE
              CALL PORSTY_BH( NBN,PX,PCMP_BH(NBN),PORD_BH(M,NBN) )
              PORD_BH(M,NBN) = MAX( PORD_BH(M,NBN),EPSL )
            ENDIF
!
!---        Saturation, relative permeability  ---
!
            CALL RKSP_BH_GT( PG_BH(M,NBN),PL_BH(M,NBN),
     &        RKG_BH(M,NBN),RKL_BH(M,NBN),SG_BH(M,NBN),
     &        SL_BH(M,NBN),NBN )
!
!---        Gas density and component fractions  ---
!
            CALL AIRGSD( T_BH(M,NBN),PVA_BH(M,NBN),RHOGAX )
            RHOG_BH(M,NBN) = XGA_BH(M,NBN)*RHOGAX +
     &        XGW_BH(M,NBN)*RHOGWX
            WTMGX = XMGA_BH(M,NBN)*WTMA + XMGW_BH(M,NBN)*WTMW
            RHOMG_BH(M,NBN) = RHOG_BH(M,NBN)/WTMGX
!
!---        Gas viscosity  ---
!
            CALL AIRGSV( T_BH(M,NBN),VISGAX )
            CALL VISC_W( T_BH(M,NBN),PVBX,RHOGWX,VISGWX )
            CALL VISC_G( VISGAX,VISGWX,XMGA_BH(M,NBN),XMGW_BH(M,NBN),
     &        VISG_BH(M,NBN) )
!
!---        Aqueous density and molar density  ---
!
            RHOL_BH(M,NBN) = RHOBX
            WTMLX = XMLA_BH(M,NBN)*WTMA + XMLS_BH(M,NBN)*WTMS +
     &        XMLW_BH(M,NBN)*WTMW
            RHOML_BH(M,NBN) = RHOL_BH(M,NBN)/WTMLX
!
!---        Aqueous viscosity  ---
!
            CALL VISC_W( T_BH(M,NBN),PX,RHOLWX,VISLWX )
            CALL VISC_B( T_BH(M,NBN),XLS_BH(M,NBN),VISLWX,VISBX )
            CALL VISC_L( XMLA_BH(M,NBN),VISBX,VISGAX,VISL_BH(M,NBN) )
!
!---        Gas-water diffusion coefficients  ---
!
            IF( ISLC(2).EQ.1 ) THEN
              DFGW_BH(M,NBN) = DFGWC
            ELSEIF( ISLC(2).EQ.2 ) THEN
              CALL BNDFAW( T_BH(M,NBN),PGX,DFGW_BH(M,NBN) )
            ELSEIF( ISLC(2).EQ.3 ) THEN
              CALL BNDFAW( T_BH(M,NBN),PGX,DFGW_BH(M,NBN) )
              CMFF = 1.D+0 + 2.6D+0/(DFGWC**0.5)
              AMC = SL_BH(M,NBN)
              ENHF = 9.5D+0 + 6.D+0*(AMC) -
     &          8.5D+0/EXP((CMFF*AMC)**4)
              DFGW_BH(M,NBN) = ENHF*DFGW_BH(M,NBN)
            ELSEIF( ISLC(2).EQ.4 ) THEN
              CALL BNDFAW( T_BH(M,NBN),PGX,DFGW_BH(M,NBN) )
              DFEFX(1) = 9.5D+0
              DFEFX(2) = 2.0D+0
              DFEFX(3) = 8.0D+0
              DFEFX(4) = 0.5D+0
              DFEFX(5) = 3.0D+0
              ENHF = DFEFX(1)+DFEFX(2)*SL_BH(M,NBN)-
     &          (DFEFX(1)-DFEFX(4))
     &          *EXP(-((DFEFX(3)*SL_BH(M,NBN))**DFEFX(5)))
              DFGW_BH(M,NBN) = ENHF*DFGW_BH(M,NBN)
            ENDIF
!
!---        Aqueous-air diffusion coefficients  ---
!
            IF( ISLC(4).EQ.1 ) THEN
             DFLA_BH(M,NBN) = DFLAC
            ELSEIF( ISLC(4).EQ.2 ) THEN
             CALL AIRDFL( T_BH(M,NBN),VISL_BH(M,NBN),DFLA_BH(M,NBN) )
            ENDIF
!
!---        Aqueous-salt diffusion coefficient  ---
!
            IF( ISLC(4).EQ.1 ) THEN
              DFLS_BH(M,NBN) = DFLSC
            ELSEIF( ISLC(4).EQ.2 ) THEN
              CALL DIFC_LS( T_BH(M,NBN),XLS_BH(M,NBN),VISL_BH(M,NBN),
     &          DFLS_BH(M,NBN) )
            ENDIF
!
!---        Precipitated NaCl density and saturation  ---
!
            CALL DENS_S( T_BH(M,NBN),PX,RHOSP_BH(M,NBN) )
!
!---        Fully unsaturated system w/ or w/o entrapped gas
!           Water mass - aqueous saturation
!           air mass - gas pressure
!           NaCl mass - total NaCl brine mass fraction  ---
!
            IF( NPHAZ_BH(2,NBN).EQ.4 ) THEN
              SS_BH(M,NBN) = TMS_BH(M,NBN)/RHOSP_BH(M,NBN)
              SLX = EPSL
              YLS_BH(M,NBN) = TMS_BH(M,NBN)*RHOBX*SLX
            ELSE
!
!---          Precipitated salt saturation  ---
!
              SS_BH(M,NBN) = MAX(YLS_BH(M,NBN)-XLSMX,0.D+0)*RHOBX*
     &          SL_BH(M,NBN)/RHOSP_BH(M,NBN)
!
!---          NaCl volumetric concentration  ---
!
              TMS_BH(M,NBN) = YLS_BH(M,NBN)*RHOBX*SL_BH(M,NBN)*
     &          PORD_BH(M,NBN)
            ENDIF
!
!---        Nonisothermal simulation  ---
!
            IF( ISLC(30).EQ.0 ) THEN
!
!---          Gas enthalpy and internal energy  ---
!
              CALL AIRGSH( T_BH(M,NBN),PVA_BH(M,NBN),HGA_BH(M,NBN),
     &          UEGAX )
              UEG_BH(M,NBN) = XGA_BH(M,NBN)*UEGAX+XGW_BH(M,NBN)*UGWX
              HGW_BH(M,NBN) = HGWX
              HG_BH(M,NBN) = XGA_BH(M,NBN)*HGA_BH(M,NBN) +
     &          XGW_BH(M,NBN)*HGW_BH(M,NBN)
!
!---          Gas thermal conductivity  ---
!
              CALL AIRGSK( T_BH(M,NBN),THKGAX )
              CALL THK_W( T_BH(M,NBN),PGX,RHOGWX,THKGWX )
              CALL THK_G( T_BH(M,NBN),THKGAX,THKGWX,XMGA_BH(M,NBN),
     &          XMGW_BH(M,NBN),THKG_BH(M,NBN) )
!
!---          Aqueous enthalpy and internal energy  ---
!
              CALL ENTH_B( T_BH(M,NBN),XLS_BH(M,NBN),HLWX,HBX )
              HL_BH(M,NBN) = MAX(1.D+0-XLA_BH(M,NBN),0.D+0)*HBX +
     &          XLA_BH(M,NBN)*HGA_BH(M,NBN)
!
!---          Thermo-catalytic fluid, salt tracking  ---
!
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  ) THEN
                VARX = MIN(1.D+0,MAX((PTC_BH(5)-XLS_BH(M,NBN))/
     &            PTC_BH(5),0.D+0))
                HL_BH(M,NBN) = HL_BH(M,NBN) +
     &            PTC_BH(1)*VARX/RHOL_BH(M,NBN)
              ENDIF
!
!---          Aqueous thermal conductivity  ---
!
              CALL THK_W( T_BH(M,NBN),PX,RHOLWX,THKLWX )
              CALL THK_B( T_BH(M,NBN),XLS_BH(M,NBN),THKLWX,
     &          THKL_BH(M,NBN) )
!
!---          Precipitated NaCl enthalpy  ---
!
              CALL ENTH_S( T_BH(M,NBN),HSP_BH(M,NBN) )
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_FRC_GT
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
!     Compute hydrologic, thermodynamic and physical properties for
!     fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 8 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEO_MECH
      USE FDVT_FRC
      USE FDVS_FRC
      USE FDVP_FRC
      USE FDVP
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
      REAL*8 DFEFX(5)
      REAL*8, DIMENSION(3,3) :: SIGX
      REAL*8, DIMENSION(1) :: SIGNX
      REAL*8, DIMENSION(3) :: VCX
      REAL*8, DIMENSION(4) :: PROP_GMX
      REAL*8 KNX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP_FRC_GT'
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
!
!---    Assign gas entry pressure and minimum gas saturation
!       for transition to unsaturated conditions---
!
        ENPR = RKSP_FRC(1,NFX)*RHORL*GRAV
!
!---    Loop over fracture triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE

          POR0_FRC(1,NTX) = POR0_FRC(1,NTX)
          POR0_FRC(2,NTX) = POR0_FRC(2,NTX)

!
!---      Loop over increment indices  ---
!
          DO M = 2,ISVC+2
            PLX = PL_FRC(M,NTX)+PATM
            PGX = PG_FRC(M,NTX)+PATM
            INDX = 0
            CALL REGION_4( T_FRC(M,NTX),PSW_FRC(M,NTX),INDX )
!
!---        Saturated system w/o entrapped gas
!           Energy - temperature
!           Water mass - aqueous pressure
!           Air mass - aqueous-air mole fraction
!           NaCl mass - total NaCl brine mass fraction  ---
!
            IF( NPHAZ_FRC(2,NTX).EQ.1 ) THEN
              PX = PLX
              PWX = PLX
              PCX = 0.D+0
              SL_FRC(M,NTX) = 1.D+0
              CALL SOL_BRNS( T_FRC(M,NTX),PWX,XLSMX )
              XLS_FRC(M,NTX) = MIN(YLS_FRC(M,NTX),XLSMX)
              CALL SP_B( T_FRC(M,NTX),XLS_FRC(M,NTX),PSBX )
              CALL P_IAPWS( T_FRC(M,NTX),PWX,RHOGWX,RHOLWX,HGWX,HLWX,
     &          UGWX,ULWX )
              CALL DENS_B( XLS_FRC(M,NTX),RHOBX,RHOLWX,T_FRC(M,NTX) )
              PVBX = PSBX
!              PVA_FRC(M,NTX) = XMLA_FRC(M,NTX)*HCAW
              XMLA_FRC(M,NTX) = PVA_FRC(M,NTX)/HCAW
              PVW_FRC(M,NTX) = PVBX
              XMGA_FRC(M,NTX) = PVA_FRC(M,NTX)/
     &          (PVA_FRC(M,NTX)+PVW_FRC(M,NTX))
              CALL EQUIL( XGA_FRC(M,NTX),XGW_FRC(M,NTX),XLA_FRC(M,NTX),
     &          XLS_FRC(M,NTX),XLW_FRC(M,NTX),XMGA_FRC(M,NTX),
     &          XMGW_FRC(M,NTX),XMLA_FRC(M,NTX),XMLS_FRC(M,NTX),
     &          XMLW_FRC(M,NTX) )
!
!---        Unsaturated system w/ or w/o entrapped gas
!
!           Energy - temperature
!           Water mass - aqueous pressure
!           Air mass - gas pressure
!           NaCl mass - total NaCl brine mass fraction  ---
!
            ELSEIF( NPHAZ_FRC(2,NTX).EQ.2 ) THEN
              PX = PG_FRC(M,NTX) + PATM
              CALL P_IAPWS( T_FRC(M,NTX),PSW_FRC(M,NTX),RHOX,RHOLWX,
     &          HX,HLWX,UX,ULWX )
              PCX = PG_FRC(M,NTX)-PL_FRC(M,NTX)
              CALL VPL_B( T_FRC(M,NTX),PSW_FRC(M,NTX),PCX,RHOLWX,PVWX )
              CALL SOL_BRNS( T_FRC(M,NTX),PSW_FRC(M,NTX),XLSMX )
              XLS_FRC(M,NTX) = MIN(YLS_FRC(M,NTX),XLSMX)
              CALL DENS_B( XLS_FRC(M,NTX),RHOBX,RHOLWX,T_FRC(M,NTX) )
              CALL SP_B( T_FRC(M,NTX),XLS_FRC(M,NTX),PSBX )
              CALL VPL_B( T_FRC(M,NTX),PSBX,PCX,RHOBX,PVBX )
              CALL P_IAPWS( T_FRC(M,NTX),PVBX,RHOGWX,RHOX,HGWX,HX,
     &          UGWX,UX )
              PVA_FRC(M,NTX) = PGX-PVBX
              IF( PVA_FRC(M,NTX).LT.1.D-6 ) PVA_FRC(M,NTX) = 0.D+0
              PVW_FRC(M,NTX) = PVBX
              XMLA_FRC(M,NTX) = PVA_FRC(M,NTX)/HCAW
              XMGA_FRC(M,NTX) = PVA_FRC(M,NTX)/PGX
              CALL EQUIL( XGA_FRC(M,NTX),XGW_FRC(M,NTX),XLA_FRC(M,NTX),
     &          XLS_FRC(M,NTX),XLW_FRC(M,NTX),XMGA_FRC(M,NTX),
     &          XMGW_FRC(M,NTX),XMLA_FRC(M,NTX),XMLS_FRC(M,NTX),
     &          XMLW_FRC(M,NTX) )
!
!---        Fully unsaturated conditions
!
!           Energy - temperature
!           Water mass - water vapor partial pressure
!           Air mass - gas pressure
!           NaCl mass - salt mass  ---
!
            ELSEIF( NPHAZ_FRC(2,NTX).EQ.4 ) THEN
              PGX = PG_FRC(M,NTX) + PATM
              PX = PGX
              INDX = 0
              CALL REGION_4( T_FRC(M,NTX),PSW_FRC(M,NTX),INDX )
              CALL SOL_BRNS( T_FRC(M,NTX),PSW_FRC(M,NTX),XLSMX )
              IF( TMS_FRC(M,NTX).GT.EPSL ) THEN
                XLS_FRC(M,NTX) = XLSMX
              ELSE
                XLS_FRC(M,NTX) = 0.D+0
              ENDIF
              PVA_FRC(M,NTX) = PGX - PVW_FRC(M,NTX)
              SL_FRC(M,NTX) = 0.D+0
              CALL CAP_FRC_GT( CPGLX,SL_FRC(M,NTX),NFX )
              PL_FRC(M,NTX) = PG_FRC(M,NTX) - PCX
              PLX = PGX - PCX
              CALL P_IAPWS( T_FRC(M,NTX),PVW_FRC(M,NTX),RHOGWX,RHOLWX,
     &          HGWX,HLWX,UGWX,ULWX )
              CALL DENS_B( XLS_FRC(M,NTX),RHOBX,RHOLWX,T_FRC(M,NTX) )
              XMLA_FRC(M,NTX) = PVA_FRC(M,NTX)/HCAW
              XMGA_FRC(M,NTX) = PVA_FRC(M,NTX)/PGX
              CALL EQUIL( XGA_FRC(M,NTX),XGW_FRC(M,NTX),XLA_FRC(M,NTX),
     &          XLS_FRC(M,NTX),XLW_FRC(M,NTX),XMGA_FRC(M,NTX),
     &          XMGW_FRC(M,NTX),XMLA_FRC(M,NTX),XMLS_FRC(M,NTX),
     &          XMLW_FRC(M,NTX) )
            ENDIF
!
!---        Constant fracture aperture model  ---
!
            IF( IJM_FRC(NFX).EQ.0 ) THEN
              APM_FRC(M,NTX) = MAX( JRC_FRC(NTX),1.D-9 )
              APH_FRC(M,NTX) = MAX( JRC_FRC(NTX),1.D-12 )
!
!---          Net normal pressure on fracture surfaces, MPa  ---
!
              PN_FRC(M,NTX) = PX - PATM
!
!---        Barton-Bandis Joint Model  ---
!
            ELSEIF( IJM_FRC(NFX).GE.1 .AND. IJM_FRC(NFX).LE.3 ) THEN
!
!---          Initial joint aperture (mm), eqn. (5),
!             Rock Mech Rock Eng. (2016) 49:837-853  ---
!
              AJX = 2.D-1*JRC_FRC(NTX)*(2.D-1*UCS_FRC(NTX)/
     &          JCS_FRC(NTX) - 1.D-1)
!
!---          Maximum joint closure (mm), eqn. (4),
!             Rock Mech Rock Eng.(2016) 49:837-853, using joint
!             compressive strength in untis of MPa  ---
!
              CSJX = 1.D-6*JCS_FRC(NTX)
              VMX = -2.96D-1 - 5.6D-3*JRC_FRC(NTX) +
     &          2.241D+0*((CSJX/AJX)**(-2.45D-1))
!
!---          Normal stifness (MPa/mm), eqn. (3), Rock Mech Rock Eng.
!             (2016) 49:837-853, using joint compressive strength
!             in untis of MPa  ---
!
              SNX = 1.78D-2*(CSJX/AJX) + 1.748D+0*JRC_FRC(NTX)
     &          - 7.155D+0
!
!---          Total normal stress (Pa) on the fracture  ---
!
              SIG_NMX = SQRT((SIGN_FRC(1,NTX)**2) + 
     &          (SIGN_FRC(2,NTX)**2) + (SIGN_FRC(3,NTX)**2))
!
!---          Effective normal stress (MPa) on the fracture  ---
!
              ESIGNX = 1.D-6*(SIG_NMX-PX)
!
!---          Mechanical aperture (micron)  ---
!
              APMX = 1.D+3*MAX(AJX-(VMX*ESIGNX/(VMX*SNX+ESIGNX)),0.D+0)
!
!---          Hydraulic aperture (micron), eqn. 7, International Journal
!             of Rock Mechanics & Mining Sciences 38 (2001) 317–329  ---
!
              APHX = MIN( (APMX**2)/(JRC_FRC(NTX)**2.5D+0),APMX )
!
!---          Mechanical and hydraulic aperture (m)  ---
!
              APM_FRC(M,NTX) = MAX( 1.D-6*APMX,1.D-9 )
              APH_FRC(M,NTX) = MAX( 1.D-6*APHX,1.D-12 )
!
!---        Dynamic Sneddon-Barton-Bandis Joint Model  ---
!
            ELSEIF( IJM_FRC(NFX).EQ.11 ) THEN
              IF( M.EQ.2 ) THEN
!
!---          Loop over fracture triangle to grid cell connections  ---
!
              AFSX = 0.D+0
              DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
                N = INCM_FRC(NCX)
                AFSX = AFSX + AFN_FRC(NCX)
              ENDDO
!
!---          Fracture triangle area weighted local stress tensor  ---
!
              DO J = 1,3
                DO I = 1,3
                  SIGX(I,J) = 0.D+0
                ENDDO
              ENDDO
              DO I = 1,4
                PROP_GMX(I) = 0.D+0
              ENDDO
              DELTAPX = 0.D+0
              DELTATX = 0.D+0
              DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
                N = INCM_FRC(NCX)
                IZN = IZ(N)
                SIGX(1,1) = SIGX(1,1) + AFN_FRC(NCX)*SIG_GM(1,N)/AFSX
                SIGX(1,2) = SIGX(1,2) + AFN_FRC(NCX)*SIG_GM(6,N)/AFSX
                SIGX(1,3) = SIGX(1,3) + AFN_FRC(NCX)*SIG_GM(5,N)/AFSX
                SIGX(2,1) = SIGX(2,1) + AFN_FRC(NCX)*SIG_GM(6,N)/AFSX
                SIGX(2,2) = SIGX(2,2) + AFN_FRC(NCX)*SIG_GM(2,N)/AFSX
                SIGX(2,3) = SIGX(2,3) + AFN_FRC(NCX)*SIG_GM(4,N)/AFSX
                SIGX(3,1) = SIGX(3,1) + AFN_FRC(NCX)*SIG_GM(5,N)/AFSX
                SIGX(3,2) = SIGX(3,2) + AFN_FRC(NCX)*SIG_GM(4,N)/AFSX
                SIGX(3,3) = SIGX(3,3) + AFN_FRC(NCX)*SIG_GM(3,N)/AFSX
                DO I = 1,4
                  PROP_GMX(I) = PROP_GMX(I) + 
     &              AFN_FRC(NCX)*PROP_GM(I,IZN)/AFSX
                ENDDO
                PEX = MAX( PG(1,N),PL(1,N) ) + PATM
                DELTAPX = DELTAPX + 
     &              AFN_FRC(NCX)*(PEX-PCMP(N))/AFSX
                DELTATX = DELTATX + 
     &              AFN_FRC(NCX)*(T(1,N)-TCMP(N))/AFSX
              ENDDO
              ENDIF
!
!---          Normal stress (Pa) on fracture triangle, where
!             fracture triangle normal vector has been aligned
!             with the local stress vector via the rotation matrix
!             ROTMAT defined through the Grid Rotation Card  ---
!
              CALL MATMUL( SIGN_FRC(1,NTX),SIGX,VCX,1,3,3 )
              CALL MATMUL( VCX,SIGN_FRC(1,NTX),SIGNX(1),1,3,1 )
              SIG_NMX = SIGNX(1)
!
!---          Poroelasticity  ---
!
              POROX = PROP_GMX(3)*DELTAPX
!
!---          Bulk modulus times 3  ---
!
              BLK3X = PROP_GMX(1)/(1.D+0-2.D+0*PROP_GMX(2))
!
!---          Thermoelasticity  ---
!
              THERMOX = PROP_GMX(4)*DELTATX/BLK3X
!
!---          Net normal pressure on fracture surfaces, MPa  ---
!
              PN_FRC(M,NTX) = PX - SIG_NMX - POROX - THERMOX - PATM
              SNX = 1.D-6*(PX - SIG_NMX - POROX - THERMOX)
!
!---          Aperture at zero net pressure, mm  ---
!
              BMAX = 1.D+3*DSBB_FRC(1,NFX)
!
!---          Young's modulus, MPa  ---
!
              YGMX = 1.D-6*PROP_GMX(1)
!
!---          Poisson's ratio  ---
!
              PRX = PROP_GMX(2)
!
!---          Fracture toughness, MPa/mm  ---
!
              KNX = 1.D-9*DSBB_FRC(2,NFX)
!
!---          Sneddon radius, mm  ---
!
              SRX = 1.D+3*JRC_FRC(NTX)
!
!---          Sneddon fractional radius, mm  ---
!
              SFRX = SRX*JCS_FRC(NTX)
!
!---          Sneddon radial scaling factor  ---
!
              SRSFX = SQRT(SRX**2-SFRX**2)
!
!---          Sigmoid function parameters  ---
!       
              PAX = 0.D+0
              PBX = 1.D+0
              STFX = 1.D+0/(1.D+0 + EXP(-(SNX-PAX)/PBX))
!
!---          Closure displacement  ---
!       
!              UNCX = MIN( MAX( SRSFX*(SNX*BMAX)/(SNX-BMAX*KNX)/SRX,
!     &          0.D+0 ),BMAX )
              UNCX = MIN( MAX( (SNX*BMAX)/(SNX-BMAX*KNX),
     &          0.D+0 ),BMAX )
!
!---          Closure aperture  ---
!       
              BCX = BMAX - UNCX
!
!---          Opening displacement  ---
!       
              UNOX = MAX( 8.D+0*(1.D+0-(PRX**2))*SNX*SRSFX/(GPI*YGMX),
     &          0.D+0 )
!
!---          Opening aperture  ---
!       
              BOX = BMAX + UNOX
!
!---          Model mechanical aperture, mm  ---
!
              APMX = (STFX*BOX + (1.D+0-STFX)*BCX)
!
!---          Model hydraulic aperture, mm  ---
!
              APHX = (STFX*BOX + (1.D+0-STFX)*BCX*SRSFX/SRX)
!
!---          Mechanical and hydraulic aperture, m  ---
!
              APM_FRC(M,NTX) = MAX( 1.D-3*APMX,1.D-12 )
              APH_FRC(M,NTX) = MAX( 1.D-3*APHX,1.D-12 )
            ENDIF
!
!---        Fracture permeability  ---
!
            PORD_FRC(M,NTX) = 1.D+0
            PERM_FRC(M,NTX) = (APH_FRC(M,NTX)**2)/12.D+0
!
!---        Saturation, relative permeability  ---
!
            CALL RKSP_FRC_GT( PG_FRC(M,NTX),PL_FRC(M,NTX),
     &        RKG_FRC(M,NTX),RKL_FRC(M,NTX),SG_FRC(M,NTX),
     &        SL_FRC(M,NTX),NFX )
!
!---        Gas density and component fractions  ---
!
            CALL AIRGSD( T_FRC(M,NTX),PVA_FRC(M,NTX),RHOGAX )
            RHOG_FRC(M,NTX) = XGA_FRC(M,NTX)*RHOGAX +
     &        XGW_FRC(M,NTX)*RHOGWX
            WTMGX = XMGA_FRC(M,NTX)*WTMA + XMGW_FRC(M,NTX)*WTMW
            RHOMG_FRC(M,NTX) = RHOG_FRC(M,NTX)/WTMGX
!
!---        Gas viscosity  ---
!
            CALL AIRGSV( T_FRC(M,NTX),VISGAX )
            CALL VISC_W( T_FRC(M,NTX),PVBX,RHOGWX,VISGWX )
            CALL VISC_G( VISGAX,VISGWX,XMGA_FRC(M,NTX),XMGW_FRC(M,NTX),
     &        VISG_FRC(M,NTX) )
!
!---        Aqueous density and molar density  ---
!
            RHOL_FRC(M,NTX) = RHOBX
            WTMLX = XMLA_FRC(M,NTX)*WTMA + XMLS_FRC(M,NTX)*WTMS +
     &        XMLW_FRC(M,NTX)*WTMW
            RHOML_FRC(M,NTX) = RHOL_FRC(M,NTX)/WTMLX
!
!---        Aqueous viscosity  ---
!
            CALL VISC_W( T_FRC(M,NTX),PX,RHOLWX,VISLWX )
            CALL VISC_B( T_FRC(M,NTX),XLS_FRC(M,NTX),VISLWX,VISBX )
            CALL VISC_L( XMLA_FRC(M,NTX),VISBX,VISGAX,VISL_FRC(M,NTX) )
!
!---        Gas-water diffusion coefficients  ---
!
            IF( ISLC(2).EQ.1 ) THEN
              DFGW_FRC(M,NTX) = DFGWC
            ELSEIF( ISLC(2).EQ.2 ) THEN
              CALL BNDFAW( T_FRC(M,NTX),PGX,DFGW_FRC(M,NTX) )
            ELSEIF( ISLC(2).EQ.3 ) THEN
              CALL BNDFAW( T_FRC(M,NTX),PGX,DFGW_FRC(M,NTX) )
              CMFF = 1.D+0 + 2.6D+0/(DFGWC**0.5)
              AMC = SL_FRC(M,NTX)
              ENHF = 9.5D+0 + 6.D+0*(AMC) -
     &          8.5D+0/EXP((CMFF*AMC)**4)
              DFGW_FRC(M,NTX) = ENHF*DFGW_FRC(M,NTX)
            ELSEIF( ISLC(2).EQ.4 ) THEN
              CALL BNDFAW( T_FRC(M,NTX),PGX,DFGW_FRC(M,NTX) )
              DFEFX(1) = 9.5D+0
              DFEFX(2) = 2.0D+0
              DFEFX(3) = 8.0D+0
              DFEFX(4) = 0.5D+0
              DFEFX(5) = 3.0D+0
              ENHF = DFEFX(1)+DFEFX(2)*SL_FRC(M,NTX)-
     &          (DFEFX(1)-DFEFX(4))
     &          *EXP(-((DFEFX(3)*SL_FRC(M,NTX))**DFEFX(5)))
              DFGW_FRC(M,NTX) = ENHF*DFGW_FRC(M,NTX)
            ENDIF
!
!---        Aqueous-air diffusion coefficients  ---
!
            IF( ISLC(4).EQ.1 ) THEN
             DFLA_FRC(M,NTX) = DFLAC
            ELSEIF( ISLC(4).EQ.2 ) THEN
             CALL AIRDFL( T_FRC(M,NTX),VISL_FRC(M,NTX),DFLA_FRC(M,NTX) )
            ENDIF
!
!---        Aqueous-salt diffusion coefficient  ---
!
            IF( ISLC(4).EQ.1 ) THEN
              DFLS_FRC(M,NTX) = DFLSC
            ELSEIF( ISLC(4).EQ.2 ) THEN
              CALL DIFC_LS( T_FRC(M,NTX),XLS_FRC(M,NTX),VISL_FRC(M,NTX),
     &          DFLS_FRC(M,NTX) )
            ENDIF
!
!---        Precipitated NaCl density and saturation  ---
!
            CALL DENS_S( T_FRC(M,NTX),PX,RHOSP_FRC(M,NTX) )
!
!---        Fully unsaturated system w/ or w/o entrapped gas
!           Water mass - aqueous saturation
!           air mass - gas pressure
!           NaCl mass - total NaCl brine mass fraction  ---
!
            IF( NPHAZ_FRC(2,NTX).EQ.4 ) THEN
              SS_FRC(M,NTX) = TMS_FRC(M,NTX)/RHOSP_FRC(M,NTX)
              SLX = EPSL
              YLS_FRC(M,NTX) = TMS_FRC(M,NTX)*RHOBX*SLX
            ELSE
!
!---          Precipitated salt saturation  ---
!
              SS_FRC(M,NTX) = MAX(YLS_FRC(M,NTX)-XLSMX,0.D+0)*RHOBX*
     &          SL_FRC(M,NTX)/RHOSP_FRC(M,NTX)
!
!---          NaCl volumetric concentration  ---
!
              TMS_FRC(M,NTX) = YLS_FRC(M,NTX)*RHOBX*SL_FRC(M,NTX)
            ENDIF
!
!---        Nonisothermal simulation  ---
!
            IF( ISLC(30).EQ.0 ) THEN
!
!---          Gas enthalpy and internal energy  ---
!
              CALL AIRGSH( T_FRC(M,NTX),PVA_FRC(M,NTX),HGA_FRC(M,NTX),
     &          UEGAX )
              UEG_FRC(M,NTX) = XGA_FRC(M,NTX)*UEGAX+XGW_FRC(M,NTX)*UGWX
              HGW_FRC(M,NTX) = HGWX
              HG_FRC(M,NTX) = XGA_FRC(M,NTX)*HGA_FRC(M,NTX) +
     &          XGW_FRC(M,NTX)*HGW_FRC(M,NTX)
!
!---          Gas thermal conductivity  ---
!
              CALL AIRGSK( T_FRC(M,NTX),THKGAX )
              CALL THK_W( T_FRC(M,NTX),PGX,RHOGWX,THKGWX )
              CALL THK_G( T_FRC(M,NTX),THKGAX,THKGWX,XMGA_FRC(M,NTX),
     &          XMGW_FRC(M,NTX),THKG_FRC(M,NTX) )
!
!---          Aqueous enthalpy and internal energy  ---
!
              CALL ENTH_B( T_FRC(M,NTX),XLS_FRC(M,NTX),HLWX,HBX )
              HL_FRC(M,NTX) = MAX(1.D+0-XLA_FRC(M,NTX),0.D+0)*HBX +
     &          XLA_FRC(M,NTX)*HGA_FRC(M,NTX)
!
!---          Thermo-catalytic fluid, salt tracking  ---
!
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  ) THEN
                VARX = MIN(1.D+0,MAX((PTC_BH(5)-XLS_FRC(M,NTX))/
     &            PTC_BH(5),0.D+0))
                HL_FRC(M,NTX) = HL_FRC(M,NTX) +
     &            PTC_BH(1)*VARX/RHOL_FRC(M,NTX)
              ENDIF
!
!---          Aqueous thermal conductivity  ---
!
              CALL THK_W( T_FRC(M,NTX),PX,RHOLWX,THKLWX )
              CALL THK_B( T_FRC(M,NTX),XLS_FRC(M,NTX),THKLWX,
     &          THKL_FRC(M,NTX) )
!
!---          Precipitated NaCl enthalpy  ---
!
              CALL ENTH_S( T_FRC(M,NTX),HSP_FRC(M,NTX) )
            ENDIF
          ENDDO
          VOL_FRC(NTX) = AF_FRC(NTX)*APM_FRC(2,NTX)
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RD_BH_GT
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
!
!     STOMP-GT
!
!     Reads the borehole card
!
!     ICC_BH(LN_BH) - cyclic borehole time index for borehole
!     ITM_BH(LN_BH) - number of time points for borehole
!     ID_BH(1,LN_BH) - starting borehole interval index for borehole
!     ID_BH(2,LN_BH) - ending borehole interval index for borehole
!     ID_BH(3,LN_BH) - starting borehole node index for borehole
!     ID_BH(4,LN_BH) - ending borehole node index for borehole
!     ID_BH(5,LN_BH) - principal borehole node index for borehole
!
!     IT_BH(1,LN_BH) - type index for starting borehole node
!     IT_BH(2,LN_BH) - type index for ending borehole node
!     Borehole types
!     1 - gas mass injection w/ zero water
!     2 - gas mass injection w/ water mass frac.
!     3 - gas mass injection w/ water rel. humidity
!     4 - aqueous mass injection w/ zero air and salt
!     5 - aqueous mass injection w/ air and salt mass frac.
!     6 - aqueous mass injection w/ air and salt rel. sat.
!     7 - fixed pressure, state condition #1 (saturated)
!     8 - fixed pressure, state condition #2 (partially saturated)
!     9 - fixed pressure, state condition #3 (unsaturated)
!    11 - gas volumetric injection w/ zero water
!    12 - gas volumetric injection w/ water mass frac.
!    13 - gas volumetric injection w/ water rel. humidity
!    14 - aqueous volumetric injection w/ zero air and salt
!    15 - aqueous volumetric injection w/ air and salt mass frac.
!    16 - aqueous volumetric injection w/ air and salt rel. sat.
!    21 - Dirichlet energy
!    22 - Neumann energy
!     N_BH - number of boreholes
!     PAR_BH(1,LI_BH) - skin factor for borehole interval
!     PAR_BH(2,LI_BH) - borehole outer radius for interval, m
!     PAR_BH(3,LI_BH) - borehole inner radius for interval, m
!     PAR_BH(4,LI_BH) - borehole wall thermal conductivity, W/m ˚K
!     PAR_BH(5,LI_BH) - borehole wall volumetric heat capacity, J/m^3 ˚K
!     PAR_BH(6,LI_BH) - borehole Colebrook roughness, m
!     PAR_BH(7,LI_BH) - borehole Brooks and Corey psi, m
!     PAR_BH(8,LI_BH) - borehole Brooks and Corey lambda
!     PAR_BH(9,LI_BH) - borehole Brooks and Corey slr
!     PAR_BH(10,LI_BH) - borehole Webb oven dried head, m
!     PAR_BH(11,LI_BH) - borehole Webb saturation matching point
!     PAR_BH(12,LI_BH) - borehole Webb head matching point, m
!     VAR_BH(2,LT_BH,LN_BH) - starting borehole node mass rate for
!                              time point, kg/s
!     VAR_BH(3,LT_BH,LN_BH) - starting borehole node production                             !                             pressure, Pa
!     VAR_BH(4,LT_BH,LN_BH) - starting borehole node aqueous air conc.                            !                             for time point
!     VAR_BH(4,LT_BH,LN_BH) - starting borehole node gas water conc.
!                             for time point
!     VAR_BH(5,LT_BH,LN_BH) - starting borehole node temperature for                            !                             time point
!     VAR_BH(6,LT_BH,LN_BH) - starting borehole node aqu. salt conc.
!                             for time point
!     VAR_BH(7,LT_BH,LN_BH) - ending borehole node mass rate for time
!                             point, kg/s
!     VAR_BH(8,LT_BH,LN_BH) - ending borehole node production pressure, !                             Pa
!     VAR_BH(9,LT_BH,LN_BH) - ending borehole node aqueous air conc.
!                             for time point
!     VAR_BH(9,LT_BH,LN_BH) - ending borehole node gas water conc
!                             for time point
!     VAR_BH(10,LT_BH,LN_BH) - ending borehole node temperature for
!                             time point
!     VAR_BH(11,LT_BH,LN_BH) - ending borehole node aqu. salt conc.
!                              for time point
!     VARC_BH((NSL-1)*2+1,LT_BH,LN_BH) - starting borehole node solute
!       conc. for time point
!     VARC_BH((NSL-1)*2+2,LT_BH,LN_BH) - ending borehole node solute
!       conc. for time point
!     XTP_BH(2,LI_BH) - x-transition points for borehole interval, m
!     YTP_BH(2,LI_BH) - y-transition points for borehole interval, m
!     ZTP_BH(2,LI_BH) - z-transition points for borehole interval, m
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE GRID
      USE GEOM_BH
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*6 FORM1
      CHARACTER*64 ADUM,BDUM,UNTS
      CHARACTER*512 CHDUM
      INTEGER NCH
      INTEGER, DIMENSION(:), ALLOCATABLE :: IZ_BHX
      LOGICAL FCHK
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM1
      DATA FORM1 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RD_BH_GT'
!
!---  Check for pre-processed borehole data  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
      IF( INDEX(ADUM(1:),'preprocess').NE.0 .AND.
     &  INDEX(ADUM(1:),'borehole').NE.0 ) THEN
        INQUIRE( FILE='ppbh.dat',EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          INDX = 4     
          CHMSG = 'Missing Preprocessed Borehole File: ppbh.dat'
          CALL WRMSGS( INDX )
        ELSE
          OPEN(UNIT=35, FILE='ppbh.dat',STATUS='OLD', FORM='FORMATTED')
        ENDIF
!
!---    Read borehole parameters  ---
!
        READ(35,*) LBC_BHX,LBC_FRCX,LBN_BHX,LFC_BHX,LI_BHX,LN_BHX,
     &    LT_BHX,LSOLU_BHX,LRCX
!
!---    Read range of borehole rock/soil type indices  ---
!
        READ(35,*) I1X,I2X
        NROCK = I2X
!
!---    Read range of borehole rock/soil names  ---
!
        J1X = I1X
        DO IX = I1X,I2X,10
          J2X = MIN( J1X+9,I2X )
          READ(35,*) (ROCK(JX),JX=J1X,J2X)
          J1X = J1X + 10
        ENDDO
!
!---    Read number of boreholes and borehole nodes  ---
!
        READ(35,*) N_BH,NBN_BH
!
!---    Loop over boreholes  ---
!
        DO NBH = 1,N_BH
!
!---      Read borehole preprocessed file  ---
!
          CALL RDPREP_BH( NBH )
        ENDDO     
!
!---    Close borehole preprocessed file  ---
!
        CLOSE(UNIT=35)
!
!---    Define borehole nodes, check borehole trajectory,
!       and write borehole.dat file  ---
!
        INDX = 1
        CALL CHK_BOREHOLE( INDX )
        CALL WR_BOREHOLE
!
!---    Reset subroutine string sequence  ---
!
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ELSE
        BACKSPACE(UNIT=IRD)
      ENDIF
!
!---  No pre-processed borehole data  ---
!
      NROCKX = NROCK+1
!
!---  Allocate temporary memory for rock/soil index  ---
!
      ALLOCATE( IZ_BHX(1:LI_BH),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IZ_BHX'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write card information to ouput file  ---
!
      CARD = 'Borehole Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of boreholes  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Boreholes'
      CALL RDINT(ISTART,ICOMMA,CHDUM,N_BH)
!
!---  Check number of boreholes parameter  ---
!
      IF( N_BH.GT.LN_BH ) THEN
        INDX = 5
        CHMSG = 'Number of Boreholes > LN_BH'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Loop over number of boreholes  ---
!
      NIT_BH = 0
      DO 400 NBH = 1,N_BH
!
!---    Read borehole type
!       (if start and end nodes and not specified then
!        assign to starting node) ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        JNDX = 1
        DO M = 1,2
        IF( JNDX.EQ.0 ) CYCLE
        VARB = 'Borehole Type'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
        IF( INDEX(ADUM(1:),'start').NE.0 ) THEN
          INDX = 1
          JNDX = 1
        ELSEIF( INDEX(ADUM(1:),'end').NE.0 ) THEN
          INDX = 2
          JNDX = 1
        ELSE
          INDX = 1
          JNDX = 0
        ENDIF
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'injection').NE.0 ) THEN
          IF( INDX.EQ.2 ) THEN
            WRITE(IWR,'(2X,A)') 'Ending Gas Mass Injection Borehole'
          ELSE
            WRITE(IWR,'(2X,A)') 'Starting Gas Mass Injection Borehole'
          ENDIF
          VARB = 'Water-Vapor Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          NM_BH(NBH) = 'GIBH'
          ICX = ICOUNT(NBH)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(NM_BH(NBH)(5:5+ICX-1),FORM1) NBH
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
            IT_BH(INDX,NBH) = 3
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
            IT_BH(INDX,NBH) = 2
          ELSE
            WRITE(IWR,'(2X,A)') 'Zero Water-Vapor'
            IT_BH(INDX,NBH) = 1
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'injection').NE.0 ) THEN
          IF( INDX.EQ.2 ) THEN
            WRITE(IWR,'(2X,A)') 'Ending Aqueous Mass Injection Borehole'
          ELSE
            WRITE(IWR,'(2X,A)') 'Starting Aqueous Mass Injection ' //
     &        'Borehole'
          ENDIF
          VARB = 'Dissolved Air and Salt Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          NM_BH(NBH) = 'LIBH'
          ICX = ICOUNT(NBH)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(NM_BH(NBH)(5:5+ICX-1),FORM1) NBH
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Relative ' //
     &        'Saturation'
            IT_BH(INDX,NBH) = 6
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Mass Fraction'
            IT_BH(INDX,NBH) = 5
          ELSE
            WRITE(IWR,'(2X,A)') 'Zero Dissolved Air and Salt '
            IT_BH(INDX,NBH) = 4
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &    INDEX(ADUM(1:),'injection').NE.0 ) THEN
          IF( INDX.EQ.2 ) THEN
            WRITE(IWR,'(2X,A)') 'Ending Gas Volumetric Injection ' //
     &        'Borehole'
          ELSE
            WRITE(IWR,'(2X,A)') 'Starting Gas Volumetric Injection ' //
     &        'Borehole'
          ENDIF
          VARB = 'Water-Vapor Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          NM_BH(NBH) = 'GIBH'
          ICX = ICOUNT(NBH)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(NM_BH(NBH)(5:5+ICX-1),FORM1) NBH
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
            IT_BH(INDX,NBH) = 13
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
            IT_BH(INDX,NBH) = 12
          ELSE
            WRITE(IWR,'(2X,A)') 'Zero Water-Vapor'
            IT_BH(INDX,NBH) = 11
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &    INDEX(ADUM(1:),'injection').NE.0 ) THEN
          IF( INDX.EQ.2 ) THEN
            WRITE(IWR,'(2X,A)') 'Ending Aqueous Volumetric ' //
     &        'Injection Borehole'
          ELSE
            WRITE(IWR,'(2X,A)') 'Starting Aqueous Volumetric ' //
     &        'Injection Borehole'
          ENDIF
          VARB = 'Dissolved Air and Salt Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          NM_BH(NBH) = 'LIBH'
          ICX = ICOUNT(NBH)
          WRITE(FORM1(3:3),'(I1)') ICX
          WRITE(NM_BH(NBH)(5:5+ICX-1),FORM1) NBH
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Relative ' //
     &        'Saturation'
            IT_BH(INDX,NBH) = 16
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Mass Fraction'
            IT_BH(INDX,NBH) = 15
          ELSE
            WRITE(IWR,'(2X,A)') 'Zero Dissolved Air and Salt '
            IT_BH(INDX,NBH) = 14
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IF( INDX.EQ.2 ) THEN
            WRITE(IWR,'(2X,A)') 'Ending Pressure Borehole'
          ELSE
            WRITE(IWR,'(2X,A)') 'Starting Pressure Borehole'
          ENDIF
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'bc1').NE.0 .OR.
     &      INDEX(BDUM(1:),'bh1').NE.0 .OR.
     &      INDEX(BDUM(1:),'ps1').NE.0 ) THEN
            NM_BH(NBH) = 'PS1'
            ICX = ICOUNT(NBH)
            WRITE(FORM1(3:3),'(I1)') ICX
            WRITE(NM_BH(NBH)(4:4+ICX-1),FORM1) NBH
            WRITE(IWR,'(2X,A)') 'Borehole State Condition #1 ' //
     &        '(Saturated)'
            IT_BH(INDX,NBH) = 7
          ELSEIF( INDEX(BDUM(1:),'bc2').NE.0 .OR.
     &      INDEX(BDUM(1:),'bh2').NE.0 .OR.
     &      INDEX(BDUM(1:),'ps2').NE.0 ) THEN
            NM_BH(NBH) = 'PS2'
            ICX = ICOUNT(NBH)
            WRITE(FORM1(3:3),'(I1)') ICX
            WRITE(NM_BH(NBH)(4:4+ICX-1),FORM1) NBH
            WRITE(IWR,'(2X,A)') 'Borehole State Condition #2 ' //
     &        '(Partially Saturated)'
            IT_BH(INDX,NBH) = 8
          ELSEIF( INDEX(BDUM(1:),'bc3').NE.0 .OR.
     &      INDEX(BDUM(1:),'bh3').NE.0 .OR.
     &      INDEX(BDUM(1:),'ps3').NE.0 ) THEN
            NM_BH(NBH) = 'PS3'
            ICX = ICOUNT(NBH)
            WRITE(FORM1(3:3),'(I1)') ICX
            WRITE(NM_BH(NBH)(4:4+ICX-1),FORM1) NBH
            WRITE(IWR,'(2X,A)') 'Borehole State Condition #3 ' //
     &        '(Unsaturated)'
            IT_BH(INDX,NBH) = 9
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Borehole State: ' // BDUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'dirichlet').NE.0 .AND.
     &    INDEX(ADUM(1:),'energy').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') ',Dirichlet Energy Borehole'
            IT_BH(1,NBH) = 21
            IT_BH(2,NBH) = 0
            JNDX = 0
            IF( ISLC(30).NE.0 ) THEN
              INDX = 4
              CHMSG = 'Dirichlet Energy Borehole w/ ' //
     &          'Isothermal Simulation: ' // ADUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
        ELSEIF( INDEX(ADUM(1:),'neumann').NE.0 .AND.
     &    INDEX(ADUM(1:),'energy').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') ',Neumann Energy Borehole'
            IT_BH(1,NBH) = 22
            IT_BH(2,NBH) = 0
            JNDX = 0
            IF( ISLC(30).NE.0 ) THEN
              INDX = 4
              CHMSG = 'Neumann Energy Borehole w/ ' //
     &          'Isothermal Simulation: ' // ADUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
        ELSEIF( INDEX(ADUM(1:),'coaxial').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') ',Coaxial Borehole'
!
!---      Inner borehole  ---
!
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &      INDEX(ADUM(1:),'mass').NE.0 .AND.
     &      INDEX(ADUM(1:),'injection').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Inner Gas Mass Injection Borehole'
            NM_BH(NBH) = 'CIGI'
            VARB = 'Water-Vapor Option'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
              IT_BH(1,NBH) = 10003
            ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(1:),'frac').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
              IT_BH(1,NBH) = 10002
            ELSE
              WRITE(IWR,'(2X,A)') 'Zero Water-Vapor'
              IT_BH(1,NBH) = 10001
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &      INDEX(ADUM(1:),'mass').NE.0 .AND.
     &      INDEX(ADUM(1:),'injection').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Inner Aqueous Mass Injection Borehole'
            NM_BH(NBH) = 'CILI'
            VARB = 'Dissolved Air and Salt Option'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Relative ' //
     &          'Saturation'
              IT_BH(1,NBH) = 10006
            ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(1:),'frac').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Mass Fraction'
              IT_BH(1,NBH) = 10005
            ELSE
              WRITE(IWR,'(2X,A)') 'Zero Dissolved Air and Salt '
              IT_BH(1,NBH) = 10004
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &      INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &      INDEX(ADUM(1:),'injection').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Inner Gas Volumetric Injection ' //
     &          'Borehole'
            NM_BH(NBH) = 'CIGI'
            VARB = 'Water-Vapor Option'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
              IT_BH(1,NBH) = 10013
            ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(1:),'frac').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
              IT_BH(1,NBH) = 10012
            ELSE
              WRITE(IWR,'(2X,A)') 'Zero Water-Vapor'
              IT_BH(1,NBH) = 10011
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &      INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &      INDEX(ADUM(1:),'injection').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Inner Aqueous Volumetric ' //
     &        'Injection Borehole'
            NM_BH(NBH) = 'CILI'
            VARB = 'Dissolved Air and Salt Option'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Relative ' //
     &          'Saturation'
              IT_BH(1,NBH) = 10016
            ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(1:),'frac').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Mass Fraction'
              IT_BH(1,NBH) = 10015
            ELSE
              WRITE(IWR,'(2X,A)') 'Zero Dissolved Air and Salt '
              IT_BH(1,NBH) = 10014
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'pressure').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Inner Pressure Borehole'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            IF( INDEX(BDUM(1:),'bc1').NE.0 .OR.
     &        INDEX(BDUM(1:),'bh1').NE.0 .OR.
     &        INDEX(BDUM(1:),'ps1').NE.0 ) THEN
              NM_BH(NBH) = 'CIPS1'
              WRITE(IWR,'(2X,A)') 'Borehole State Condition #1 ' //
     &          '(Saturated)'
              IT_BH(1,NBH) = 10007
              IT_BH(2,NBH) = 10007
            ELSEIF( INDEX(BDUM(1:),'bc2').NE.0 .OR.
     &        INDEX(BDUM(1:),'bh2').NE.0 .OR.
     &        INDEX(BDUM(1:),'ps2').NE.0 ) THEN
              NM_BH(NBH) = 'CIPS2'
              WRITE(IWR,'(2X,A)') 'Borehole State Condition #2 ' //
     &          '(Partially Saturated)'
              IT_BH(1,NBH) = 10008
              IT_BH(2,NBH) = 10008
            ELSEIF( INDEX(BDUM(1:),'bc3').NE.0 .OR.
     &        INDEX(BDUM(1:),'bh3').NE.0 .OR.
     &        INDEX(BDUM(1:),'ps3').NE.0 ) THEN
              NM_BH(NBH) = 'CIPS2'
              WRITE(IWR,'(2X,A)') 'Borehole State Condition #3 ' //
     &          '(Unsaturated)'
              IT_BH(1,NBH) = 10009
              IT_BH(2,NBH) = 10009
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Borehole State: ' // BDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
!
!---      Outer borehole  ---
!
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &      INDEX(ADUM(1:),'mass').NE.0 .AND.
     &      INDEX(ADUM(1:),'injection').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Outer Gas Mass Injection Borehole'
            NCH = INDEX(NM_BH(NBH),'  ')-1
            WRITE(NM_BH(NBH)(NCH+1:NCH+5),'(A)') 'OGIBH'
            VARB = 'Water-Vapor Option'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
              IT_BH(2,NBH) = 10003
            ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(1:),'frac').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
              IT_BH(2,NBH) = 10002
            ELSE
              WRITE(IWR,'(2X,A)') 'Zero Water-Vapor'
              IT_BH(2,NBH) = 10001
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &      INDEX(ADUM(1:),'mass').NE.0 .AND.
     &      INDEX(ADUM(1:),'injection').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Outer Aqueous Mass Injection Borehole'
            NCH = INDEX(NM_BH(NBH),'  ')-1
            WRITE(NM_BH(NBH)(NCH+1:NCH+5),'(A)') 'OLIBH'
            VARB = 'Dissolved Air and Salt Option'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Relative ' //
     &          'Saturation'
              IT_BH(2,NBH) = 10006
            ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(1:),'frac').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Mass Fraction'
              IT_BH(2,NBH) = 10005
            ELSE
              WRITE(IWR,'(2X,A)') 'Zero Dissolved Air and Salt '
              IT_BH(2,NBH) = 10007
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &      INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &      INDEX(ADUM(1:),'injection').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Outer Gas Volumetric Injection ' //
     &          'Borehole'
            NCH = INDEX(NM_BH(NBH),'  ')-1
            WRITE(NM_BH(NBH)(NCH+1:NCH+5),'(A)') 'OGIBH'
            VARB = 'Water-Vapor Option'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
              IT_BH(2,NBH) = 10013
            ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(1:),'frac').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
              IT_BH(2,NBH) = 10012
            ELSE
              WRITE(IWR,'(2X,A)') 'Zero Water-Vapor'
              IT_BH(2,NBH) = 10011
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &      INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &      INDEX(ADUM(1:),'injection').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Outer Aqueous Volumetric ' //
     &        'Injection Borehole'
            NCH = INDEX(NM_BH(NBH),'  ')-1
            WRITE(NM_BH(NBH)(NCH+1:NCH+5),'(A)') 'OLIBH'
            VARB = 'Dissolved Air and Salt Option'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Relative ' //
     &          'Saturation'
              IT_BH(2,NBH) = 10016
            ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(1:),'frac').NE.0 ) THEN
              WRITE(IWR,'(2X,A)') 'Dissolved Air and Salt Mass Fraction'
              IT_BH(2,NBH) = 10015
            ELSE
              WRITE(IWR,'(2X,A)') 'Zero Dissolved Air and Salt '
              IT_BH(2,NBH) = 10014
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'pressure').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Inner Pressure Borehole'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            IF( INDEX(BDUM(1:),'bc1').NE.0 .OR.
     &        INDEX(BDUM(1:),'bh1').NE.0 .OR.
     &        INDEX(BDUM(1:),'ps1').NE.0 ) THEN
              NCH = INDEX(NM_BH(NBH),'  ')-1
              WRITE(NM_BH(NBH)(NCH+1:NCH+6),'(A)') 'OPS1BH'
              WRITE(IWR,'(2X,A)') 'Borehole State Condition #1 ' //
     &          '(Saturated)'
              IT_BH(2,NBH) = 10007
            ELSEIF( INDEX(BDUM(1:),'bc2').NE.0 .OR.
     &        INDEX(BDUM(1:),'bh2').NE.0 .OR.
     &        INDEX(BDUM(1:),'ps2').NE.0 ) THEN
              NCH = INDEX(NM_BH(NBH),'  ')-1
              WRITE(NM_BH(NBH)(NCH+1:NCH+6),'(A)') 'OPS2BH'
              WRITE(IWR,'(2X,A)') 'Borehole State Condition #2 ' //
     &          '(Partially Saturated)'
              IT_BH(2,NBH) = 10008
            ELSEIF( INDEX(BDUM(1:),'bc3').NE.0 .OR.
     &        INDEX(BDUM(1:),'bh3').NE.0 .OR.
     &        INDEX(BDUM(1:),'ps3').NE.0 ) THEN
              NCH = INDEX(NM_BH(NBH),'  ')-1
              WRITE(NM_BH(NBH)(NCH+1:NCH+6),'(A)') 'OPS3BH'
              WRITE(IWR,'(2X,A)') 'Borehole State Condition #3 ' //
     &          '(Unsaturated)'
              IT_BH(2,NBH) = 10009
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Borehole State: ' // BDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
          ICX = ICOUNT(NBH)
          WRITE(FORM1(3:3),'(I1)') ICX
          NCH = INDEX(NM_BH(NBH),'  ')-1
          WRITE(NM_BH(NBH)(NCH+1:NCH+ICX),FORM1) NBH
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Borehole Type: ' // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for solutes  ---
!
        IF( IEQC.NE.0 ) THEN
          NC = 0
          DO
            CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
            ISX = ISTART
            ICX = ICOMMA
            IF( INDX.EQ.0 ) EXIT
!
!---        Check for solute names  ---
!
            ADUM(1:) = ' '
            VARB = 'Solute Name'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
            IFIND = 0
            DO NSL = 1,NSOLU
!
!---          Solute name found  ---
!
              IF( SOLUT(NSL).EQ.ADUM ) THEN
                NC = NC + 1
                ISOLU_BH(NC,NBH) = NSL
                IFIND = 1
                EXIT
              ENDIF
            ENDDO
            IF( IFIND.EQ.0 ) EXIT
          ENDDO
          ISTART = ISX
          ICOMMA = ICX
        ENDIF
        ENDDO
!
!---    Read borehole name  ---
!
        CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
        VARB = 'Borehole Name'
        IF( INDX.EQ.1 ) THEN
          IDFLT = 1
          CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,NM_BH(NBH) )
        ENDIF
        WRITE(IWR,'(2X,3A)') VARB(1:IVR),': ',NM_BH(NBH)(1:NCH)
!
!---    Read number of borehole intervals  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Number of Borehole Intervals'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NI_BHX)
        WRITE(IWR,'(2X,2A,I6)') VARB(1:IVR),': ',NI_BHX
        NIT_BH = NIT_BH + NI_BHX
!
!---    Check total number of borehole intervals parameter  ---
!
        IF( NIT_BH.GT.LI_BH ) THEN
          INDX = 5
          CHMSG = 'Number of Borehole Intervals > LI_BH'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for interval increment input type  ---
!
        III_BHX = 0
        CALL CHKCHR(ISTART,ICOMMA,CHDUM,INDX )
        IF( INDX.EQ.1 ) THEN
          IDFLT = 1
          CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,BDUM )
          IF( INDEX(BDUM(1:),'increment').NE.0 ) THEN
            III_BHX = 1
!
!---        Read first x-transition point  ---
!
            VARB = 'Interval First X-Transition Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,XTP_BH(1,1))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',XTP_BH(1,1)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,XTP_BH(1,1),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',XTP_BH(1,1),', m)'
!
!---        Read first y-transition point  ---
!
            VARB = 'Interval First Y-Transition Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,YTP_BH(1,1))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',YTP_BH(1,1)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,YTP_BH(1,1),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',YTP_BH(1,1),', m)'
!
!---        Cylindrical coordinates with azimuthal symmetry  ---
!
            IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
              IF( ABS(XTP_BH(1,1))/EPSL.GT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Non-Zero Interval First X-Transition Point ' //
     &            'for Radially Symmetric Domain: XTP_BH ='
                RLMSG = XTP_BH(1,1)
                CALL WRMSGS( INDX )
              ENDIF
              IF( ABS(YTP_BH(1,1))/EPSL.GT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Non-Zero Interval First Y-Transition Point ' //
     &            'for Radially Symmetric Domain: YTP_BH ='
                RLMSG = YTP_BH(1,1)
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
!
!---        Read first z-transition point  ---
!
            VARB = 'Interval First Z-Transition Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,ZTP_BH(1,1))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',ZTP_BH(1,1)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,ZTP_BH(1,1),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',ZTP_BH(1,1),', m)'
          ENDIF
        ENDIF
!
!---    Assign borehole transition pointers  ---
!
        IF( NBH.EQ.1 ) THEN
          ID_BH(1,NBH) = 1
        ELSE
          ID_BH(1,NBH) = ID_BH(2,NBH-1) + 1
        ENDIF
        ID_BH(2,NBH) = ID_BH(1,NBH) + NI_BHX - 1
!
!---    Loop over number of borehole intervals  ---
!
        DO 100 NIBH = ID_BH(1,NBH),ID_BH(2,NBH)
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          IF( III_BHX.EQ.1 ) THEN
            IF( NIBH.GT.ID_BH(1,NBH) ) THEN
              XTP_BH(1,NIBH) = XTP_BH(2,NIBH-1)
              YTP_BH(1,NIBH) = YTP_BH(2,NIBH-1)
              ZTP_BH(1,NIBH) = ZTP_BH(2,NIBH-1)
            ENDIF
!
!---        Read x-direction interval increment  ---
!
            VARB = 'X-Direction Interval Increment'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,XTP_BH(2,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',XTP_BH(2,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,XTP_BH(2,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',XTP_BH(2,NIBH),', m)'
            XTP_BH(2,NIBH) = XTP_BH(1,NIBH) + XTP_BH(2,NIBH)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,A,1PE11.4,A)') 'Interval Second ' // 
     &          'X-Transition Point (',XTP_BH(2,NIBH),', m)'
!
!---        Read y-direction interval increment  ---
!
            VARB = 'Y-Direction Interval Increment'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,YTP_BH(2,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',YTP_BH(2,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,YTP_BH(2,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',YTP_BH(2,NIBH),', m)'
            YTP_BH(2,NIBH) = YTP_BH(1,NIBH) + YTP_BH(2,NIBH)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,A,1PE11.4,A)') 'Interval Second ' // 
     &          'Y-Transition Point (',YTP_BH(2,NIBH),', m)'
!
!---        Cylindrical coordinates with azimuthal symmetry  ---
!
            IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
              IF( ABS(XTP_BH(2,NIBH))/EPSL.GT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Non-Zero Interval Second X-Transition ' //
     &            'Point for Radially Symmetric Domain: XTP_BH ='
                RLMSG = XTP_BH(2,NIBH)
                CALL WRMSGS( INDX )
              ENDIF
              IF( ABS(YTP_BH(2,NIBH))/EPSL.GT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Non-Zero Interval Second Y-Transition ' //
     &            'Point for Radially Symmetric Domain: YTP_BH ='
                RLMSG = YTP_BH(2,NIBH)
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
!
!---        Read z-direction interval increment  ---
!
            VARB = 'Z-Direction Interval Increment'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,ZTP_BH(2,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',ZTP_BH(2,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,ZTP_BH(2,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',ZTP_BH(2,NIBH),', m)'
            ZTP_BH(2,NIBH) = ZTP_BH(1,NIBH) + ZTP_BH(2,NIBH)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,A,1PE11.4,A)') 'Interval Second ' // 
     &          'Z-Transition Point (',ZTP_BH(2,NIBH),', m)'
          ELSE
!
!---        Read first x-transition point  ---
!
            VARB = 'Interval First X-Transition Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,XTP_BH(1,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',XTP_BH(1,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,XTP_BH(1,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',XTP_BH(1,NIBH),', m)'
!
!---        Read first y-transition point  ---
!
            VARB = 'Interval First Y-Transition Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,YTP_BH(1,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',YTP_BH(1,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,YTP_BH(1,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',YTP_BH(1,NIBH),', m)'
!
!---        Cylindrical coordinates with azimuthal symmetry  ---
!
            IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
              IF( ABS(XTP_BH(1,NIBH))/EPSL.GT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Non-Zero Interval First X-Transition Point ' //
     &            'for Radially Symmetric Domain: XTP_BH ='
                RLMSG = XTP_BH(1,NIBH)
                CALL WRMSGS( INDX )
              ENDIF
              IF( ABS(YTP_BH(1,NIBH))/EPSL.GT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Non-Zero Interval First Y-Transition Point ' //
     &            'for Radially Symmetric Domain: YTP_BH ='
                RLMSG = YTP_BH(1,NIBH)
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
!
!---        Read first z-transition point  ---
!
            VARB = 'Interval First Z-Transition Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,ZTP_BH(1,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',ZTP_BH(1,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,ZTP_BH(1,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',ZTP_BH(1,NIBH),', m)'
!
!---        Read second x-transition point  ---
!
            VARB = 'Interval Second X-Transition Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,XTP_BH(2,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',XTP_BH(2,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,XTP_BH(2,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',XTP_BH(2,NIBH),', m)'
!
!---        Read second y-transition point  ---
!
            VARB = 'Interval Second Y-Transition Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,YTP_BH(2,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',YTP_BH(2,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,YTP_BH(2,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',YTP_BH(2,NIBH),', m)'
!
!---        Cylindrical coordinates with azimuthal symmetry  ---
!
            IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
              IF( ABS(XTP_BH(2,NIBH))/EPSL.GT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Non-Zero Interval Second X-Transition ' //
     &            'Point for Radially Symmetric Domain: XTP_BH ='
                RLMSG = XTP_BH(2,NIBH)
                CALL WRMSGS( INDX )
              ENDIF
              IF( ABS(YTP_BH(2,NIBH))/EPSL.GT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Non-Zero Interval Second Y-Transition ' //
     &            'Point for Radially Symmetric Domain: YTP_BH ='
                RLMSG = YTP_BH(2,NIBH)
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
!
!---        Read second z-transition point  ---
!
            VARB = 'Interval Second Z-Transition Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,ZTP_BH(2,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',ZTP_BH(2,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,ZTP_BH(2,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',ZTP_BH(2,NIBH),', m)'
          ENDIF
!
!---      Coaxial borehole  ---
!
          IF( IT_BH(1,NBH).GT.10000 ) THEN
            IS_BH(NIBH) = 10000
!
!---        Read borehole inner pipe inner radius  ---
!
            VARB = 'Interval Borehole Inner Pipe Inner Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(1,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(1,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_BH(1,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(1,NIBH),', m)'
!
!---        Read borehole inner pipe outer radius  ---
!
            VARB = 'Interval Borehole Inner Pipe Outer Radius'
            UNTS = 'm'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(2,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(2,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_BH(2,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(2,NIBH),', m)'
!
!---        Read borehole Outer Pipe Inner radius  ---
!
            VARB = 'Interval Borehole Outer Pipe Inner Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(3,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(3,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_BH(3,NIBH),INDX)
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(3,NIBH),', m)'
!
!---        Read borehole Outer Pipe Outer radius  ---
!
            VARB = 'Interval Borehole Outer Pipe Outer Radius'
            UNTS = 'm'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(4,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(4,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_BH(4,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(4,NIBH),', m)'
!
!---        Read borehole inner wall thermal conductivity  ---
!
            VARB = 'Interval Borehole Inner Wall Thermal Conductivity'
            UNTS = 'w/m k'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(5,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(5,NIBH)
            INDX = 0
            IUNM = 1
            IUNKG = 1
            IUNS = -3
            IUNK = -1
            CALL RDUNIT(UNTS,PAR_BH(5,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(5,NIBH),', W/m K)'
            IF( PAR_BH(5,NIBH).LT.1.D-3 ) THEN
              PAR_BH(5,NIBH) = 1.D-3
              WRITE(IWR,'(A,1PE11.4,A)') 'Minimum Borehole Wall Thermal'
     &         // 'Conductivity Imposed: ',PAR_BH(5,NIBH),', m'
            ENDIF
!
!---        Read borehole inner wall specific heat  ---
!
            VARB = 'Interval Borehole Inner Wall Volumetric Heat ' //
     &        'Capacity'
            UNTS = 'j/m^3 k'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(6,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(6,NIBH)
            INDX = 0
            IUNM = -1
            IUNS = -2
            IUNK = -1
            IUNKG = 1
            CALL RDUNIT(UNTS,PAR_BH(6,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(6,NIBH),', J/m^3 K)'
!
!---        Read borehole outer Colebrook roughness  ---
!
            VARB = 'Interval Borehole Inner Colebrook Roughness'
            UNTS = 'm'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(13,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(13,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_BH(13,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(13,NIBH),', m)'
            IF( PAR_BH(13,NIBH).LT.1.D-6 ) THEN
              PAR_BH(13,NIBH) = 1.D-6
              WRITE(IWR,'(A,1PE11.4,A)') 'Minimum Colebrook Roughness'
     &         // 'Imposed: ',PAR_BH(13,NIBH),', m'
            ENDIF
!
!---        Read borehole outer wall thermal conductivity  ---
!
            VARB = 'Interval Borehole Outer Wall Thermal Conductivity'
            UNTS = 'w/m k'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(14,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(14,NIBH)
            INDX = 0
            IUNM = 1
            IUNKG = 1
            IUNS = -3
            IUNK = -1
            CALL RDUNIT(UNTS,PAR_BH(14,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(14,NIBH),', W/m K)'
            IF( PAR_BH(14,NIBH).LT.1.D-3 ) THEN
              PAR_BH(14,NIBH) = 1.D-3
              WRITE(IWR,'(A,1PE11.4,A)') 'Minimum Borehole Wall Thermal'
     &         // 'Conductivity Imposed: ',PAR_BH(14,NIBH),', m'
            ENDIF
!
!---        Read borehole outer wall specific heat  ---
!
            VARB = 'Interval Borehole Outer Wall Volumetric Heat ' //
     &        'Capacity'
            UNTS = 'j/m^3 k'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(15,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(15,NIBH)
            INDX = 0
            IUNM = -1
            IUNS = -2
            IUNK = -1
            IUNKG = 1
            CALL RDUNIT(UNTS,PAR_BH(15,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(15,NIBH),', J/m^3 K)'
!
!---        Read borehole outer Colebrook roughness  ---
!
            VARB = 'Interval Borehole Outer Colebrook Roughness'
            UNTS = 'm'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(16,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(16,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_BH(16,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(16,NIBH),', m)'
            IF( PAR_BH(16,NIBH).LT.1.D-6 ) THEN
              PAR_BH(16,NIBH) = 1.D-6
              WRITE(IWR,'(A,1PE11.4,A)') 'Minimum Colebrook Roughness'
     &         // 'Imposed: ',PAR_BH(16,NIBH),', m'
            ENDIF
!
!---        Read borehole Brooks and Corey (psi)  ---
!
            VARB = 'Interval Borehole Brooks and Corey (psi)'
            UNTS = 'm'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(7,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(7,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_BH(7,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(7,NIBH),', m)'
!
!---        Read borehole Brooks and Corey (lambda)
!
            VARB = 'Interval Borehole Brooks and Corey (lambda)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(8,NIBH))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',PAR_BH(8,NIBH)
!
!---        Read borehole Brooks and Corey (slr)
!
            VARB = 'Interval Borehole Brooks and Corey (slr)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(9,NIBH))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',PAR_BH(9,NIBH)
!
!---        Webb extension read optional oven dried head  ---
!
            CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
            IF( INDX.EQ.1 ) THEN
              VARB = 'Interval Borehole Oven Dried Head'
              UNTS = 'm'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(10,NIBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(10,NIBH)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,PAR_BH(10,NIBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(10,NIBH),', m)'
            ELSE
              PAR_BH(10,NIBH) = HDOD
            ENDIF
!
!---        Webb matching point parameters  ---
!
            PSIX = PAR_BH(7,NIBH)
            CLX = PAR_BH(8,NIBH)
            SRX = PAR_BH(9,NIBH)
            CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
            PAR_BH(11,NIBH) = SMPX
            VARB = 'Interval Borehole Webb Matching Point Saturation'
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',PAR_BH(11,NIBH)
            PAR_BH(12,NIBH) = HMPX
            VARB = 'Interval Borehole Webb Matching Point Head'
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ', m : ',PAR_BH(12,NIBH)
!
!---      Simple borehole  ---
!
          ELSE
!
!---      Read borehole outer radius  ---
!
          VARB = 'Interval Borehole Outer Radius'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(2,NIBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IF( ISLC(36).EQ.1 )
     &      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',PAR_BH(2,NIBH)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,PAR_BH(2,NIBH),INDX)
          IF( ISLC(36).EQ.1 )
     &      WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(2,NIBH),', m)'
!
!---      Read borehole skin factor  ---
!
          VARB = 'Interval Skin Factor'
          IDFLT = 1
          PAR_BH(1,NIBH) = 0.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(1,NIBH))
          IF( ISLC(36).EQ.1 )
     &      WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',PAR_BH(1,NIBH)
!
!---      Read borehole interval type option  ---
!
          VARB = 'Interval Type Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'uncased').NE.0 .OR.
     &      INDEX(ADUM(1:),'open').NE.0 ) THEN
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,A)') 'Uncased Filled Borehole Interval'
            IS_BH(NIBH) = 2
          ELSEIF( INDEX(ADUM(1:),'screened pipe').NE.0 .OR.
     &      INDEX(ADUM(1:),'perforated pipe').NE.0 ) THEN
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,A)') 'Perforated Pipe Borehole Interval'
            IS_BH(NIBH) = 12
          ELSEIF( INDEX(ADUM(1:),'pipe').NE.0 .AND.
     &      INDEX(ADUM(1:),'catalyst').NE.0 ) THEN
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,A)') 'Unperforated Pipe Borehole ' //
     &        'Interval with Catalyst'
            IS_BH(NIBH) = 21
          ELSEIF( INDEX(ADUM(1:),'pipe').NE.0 ) THEN
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,A)') 'Unperforated Pipe Borehole Interval'
            IS_BH(NIBH) = 11
          ELSE
            IF( ISLC(36).EQ.1 )
     &        WRITE(IWR,'(2X,A)') 'Cased Filled Borehole Interval'
            IS_BH(NIBH) = 1
          ENDIF
!
!---      Read borehole interval rock/soil/grout/fill for filled
!         borehole interval types  ---
!
          IF( IS_BH(NIBH).LE.10 ) THEN
            VARB = 'Interval Rock/Soil Name'
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,ADUM)
            IF( INDEX(ROCK(1),'indexing').NE.0 ) THEN
              DO 210 M = NROCKX,NROCK
                IF( ROCK(M).EQ.ADUM ) THEN
                  IROCK = M
                  GOTO 220
                ENDIF
  210         CONTINUE
              NROCK = NROCK+1
              IF( NROCK.GT.LRC ) THEN
                INDX = 5
                CHMSG = 'Number of Rock/Soil Types > Parameter LRC'
                CALL WRMSGS( INDX )
              ENDIF
              ROCK(NROCK) = ADUM
              IROCK = NROCK
  220         CONTINUE
              IZ_BHX(NIBH) = IROCK
            ELSE
              DO 230 M = 1,NROCK
                IF( ROCK(M).EQ.ADUM ) THEN
                  IROCK = M
                  GOTO 240
                ENDIF
  230         CONTINUE
              NROCK = NROCK+1
              IF( NROCK.GT.LRC ) THEN
                INDX = 5
                CHMSG = 'Number of Rock/Soil Types > Parameter LRC'
                CALL WRMSGS( INDX )
              ENDIF
              ROCK(NROCK) = ADUM
              IROCK = NROCK
  240         CONTINUE
              IZ_BHX(NIBH) = IROCK
            ENDIF
!
!---      Read borehole interval pipe properties  ---
!
          ELSE
!
!---        Read borehole inner radius  ---
!
            VARB = 'Interval Borehole Inner Radius'
            UNTS = 'm'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(3,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(3,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_BH(3,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(3,NIBH),', m)'
!
!---        Read borehole wall thermal conductivity  ---
!
            VARB = 'Interval Borehole Wall Thermal Conductivity'
            UNTS = 'w/m k'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(4,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(4,NIBH)
            INDX = 0
            IUNM = 1
            IUNKG = 1
            IUNS = -3
            IUNK = -1
            CALL RDUNIT(UNTS,PAR_BH(4,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(4,NIBH),', W/m K)'
            IF( PAR_BH(4,NIBH).LT.1.D-3 ) THEN
              PAR_BH(4,NIBH) = 1.D-3
              WRITE(IWR,'(A,1PE11.4,A)') 'Minimum Borehole Wall Thermal'
     &         // 'Conductivity Imposed: ',PAR_BH(4,NIBH),', m'
            ENDIF
!
!---        Read borehole wall specific heat  ---
!
            VARB = 'Interval Borehole Wall Volumetric Heat Capacity'
            UNTS = 'j/m^3 k'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(5,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(5,NIBH)
            INDX = 0
            IUNM = -1
            IUNS = -2
            IUNK = -1
            IUNKG = 1
            CALL RDUNIT(UNTS,PAR_BH(5,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(5,NIBH),', J/m^3 K)'
!
!---        Read borehole Colebrook roughness  ---
!
            VARB = 'Interval Borehole Colebrook Roughness'
            UNTS = 'm'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(6,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(6,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_BH(6,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(6,NIBH),', m)'
            IF( PAR_BH(6,NIBH).LT.1.D-6 ) THEN
              PAR_BH(6,NIBH) = 1.D-6
              WRITE(IWR,'(A,1PE11.4,A)') 'Minimum Colebrook Roughness'
     &         // 'Imposed: ',PAR_BH(6,NIBH),', m'
            ENDIF
!
!---        Read borehole Brooks and Corey (psi)  ---
!
            VARB = 'Interval Borehole Brooks and Corey (psi)'
            UNTS = 'm'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(7,NIBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(7,NIBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,PAR_BH(7,NIBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(7,NIBH),', m)'
!
!---        Read borehole Brooks and Corey (lambda)
!
            VARB = 'Interval Borehole Brooks and Corey (lambda)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(8,NIBH))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',PAR_BH(8,NIBH)
!
!---        Read borehole Brooks and Corey (slr)
!
            VARB = 'Interval Borehole Brooks and Corey (slr)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(9,NIBH))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',PAR_BH(9,NIBH)
!
!---        Webb extension read optional oven dried head  ---
!
            CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
            IF( INDX.EQ.1 ) THEN
              VARB = 'Interval Borehole Oven Dried Head'
              UNTS = 'm'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,PAR_BH(10,NIBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PAR_BH(10,NIBH)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,PAR_BH(10,NIBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',PAR_BH(10,NIBH),', m)'
            ELSE
              PAR_BH(10,NIBH) = HDOD
            ENDIF
!
!---        Webb matching point parameters  ---
!
            PSIX = PAR_BH(7,NIBH)
            CLX = PAR_BH(8,NIBH)
            SRX = PAR_BH(9,NIBH)
            CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
            PAR_BH(11,NIBH) = SMPX
            VARB = 'Interval Borehole Webb Matching Point Saturation'
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',PAR_BH(11,NIBH)
            PAR_BH(12,NIBH) = HMPX
            VARB = 'Interval Borehole Webb Matching Point Head'
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ', m : ',PAR_BH(12,NIBH)
          ENDIF
          ENDIF
  100   CONTINUE
!
!---    Read number of borehole time points  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Number of Borehole Time Points'
        CALL RDINT(ISTART,ICOMMA,CHDUM,ITM_BH(NBH))
!
!---    Check number of borehole time points parameter  ---
!
        IF( ITM_BH(NBH).GT.LT_BH ) THEN
          INDX = 5
          CHMSG = 'Number of Borehole Time Points > LT_BH'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for cyclic borehole times  ---
!
        VARB = 'Cyclic Borehole Time Option'
        CALL CHKCHR(ISTART,ICOMMA,CHDUM,INDX)
        IF( INDX.EQ.1 ) THEN
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'cyclic').NE.0 ) THEN
            ICC_BH(NBH) = 1
            WRITE(IWR,'(2X,A)') 'Cyclic Borehole Times'
          ELSE
            ICC_BH(NBH) = 0
            WRITE(IWR,'(2X,A)') 'Noncyclic Borehole Times'
          ENDIF
        ENDIF
!
!---    Loop over the number of borehole time points  ---
!
        NC = 0
        DO 300 M = 1,ITM_BH(NBH)
!
!---      Read new line  ---
!
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
!
!---      Read borehole time  ---
!
          VARB = 'Borehole Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(1,M,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',VAR_BH(1,M,NBH)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR_BH(1,M,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(1,M,NBH),', s)'
!
!---      Coaxial type borehole  ---
!
          IF( IT_BH(1,NBH).GT.10000 ) THEN
            ITI_BHX = MOD( IT_BH(1,NBH),10000 )
            ITO_BHX = MOD( IT_BH(2,NBH),10000 )
!
!---        Inner pipe  ---
!
!
!---        Read injection mass rate for inner pipe ---
!
            IF( ITI_BHX.GE.1 .AND. ITI_BHX.LE.6 ) THEN
              VARB = 'Inner Pipe Injection Mass Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(2,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(2,M,NBH)
              IUNKG = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_BH(2,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(2,M,NBH),', kg/s)'
!
!---        Read injection volumetric rate for inner pipe ---
!
            ELSEIF( ITI_BHX.GE.11 .AND. ITI_BHX.LE.16 ) THEN
              VARB = 'Inner Pipe Injection Volumetric Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(2,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(2,M,NBH)
              IUNM = 3
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_BH(2,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(2,M,NBH),', kg/s)'
!
!---        Read pressure input for state condition #1 (saturated)
!           for inner pipe  ---
!
            ELSEIF( ITI_BHX.EQ.7 ) THEN
!
!---          Read inner pipe pressure  ---
!
              VARB = 'Inner Pipe Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(3,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(3,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(3,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(3,M,NBH),', Pa)'
!
!---          Nonisoair read dissolved air relative saturation  ---
!
              IF( ISLC(37).EQ.0 ) THEN
!
!---            Read dissolved air relative saturation  ---
!
                VARB = 'Inner Pipe Dissolved Air ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(4,M,NBH)
              ENDIF
!
!---          Nonisobrine read dissolved salt relative saturation  ---
!
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Inner Pipe Dissolved Salt ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
!
!---        Read pressure input for state condition #2
!           (partially saturated) for inner pipe ---
!
            ELSEIF( ITI_BHX.EQ.8 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Inner Pipe Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(3,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(3,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(3,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(3,M,NBH),', Pa)'
!
!---          Read aqueous saturation  ---
!
              VARB = 'Inner Pipe Aqueous Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(2,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(2,M,NBH)
              IF( VAR_BH(2,M,NBH).LT.EPSL .OR.
     &          (1.D+0-VAR_BH(2,M,NBH)).LT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Out-of-Range Inner Pipe Aqueous Saturation'
                RLMSG = VAR_BH(2,M,NBH)
                CALL WRMSGS( INDX )
              ENDIF
!
!---          Nonisobrine read dissolved salt relative saturation  ---
!
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Inner Pipe Dissolved Salt ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
!
!---        Read pressure input for state condition #3 (unsaturated)
!           for inner pipe   ---
!
            ELSEIF( ITI_BHX.EQ.9 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Inner Pipe Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(3,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(3,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(3,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(3,M,NBH),', Pa)'
!
!---          Read water vapor relative humidity  ---
!
              VARB = 'Inner Pipe Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(4,M,NBH)
            ENDIF
!
!---        Read water vapor concentration in gas for inner pipe ---
!
            IF( ITI_BHX.EQ.3 .OR. ITI_BHX.EQ.13 ) THEN
              VARB = 'Inner Pipe Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(4,M,NBH)
            ELSEIF( ITI_BHX.EQ.2 .OR. ITI_BHX.EQ.12 ) THEN
              VARB = 'Inner Pipe Water Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(4,M,NBH)
            ELSEIF( ITI_BHX.EQ.1 .OR. ITI_BHX.EQ.11 ) THEN
              VAR_BH(4,M,NBH) = 0.D+0
            ENDIF
!
!---        Read air concentration in aqueous for inner pipe ---
!
            IF( ITI_BHX.EQ.6 .OR. ITI_BHX.EQ.16 ) THEN
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Inner Pipe Dissolved Air ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(4,M,NBH)
              ENDIF
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Inner Pipe Dissolved Salt ' //
     &          'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
            ELSEIF( ITI_BHX.EQ.5 .OR. ITI_BHX.EQ.15 ) THEN
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Inner Pipe Dissolved Air Mass ' //
     &            'Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(4,M,NBH)
              ENDIF
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Inner Pipe Dissolved Salt Mass ' //
     &          'Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
            ELSEIF( ITI_BHX.EQ.4 .OR. ITI_BHX.EQ.14 ) THEN
              VAR_BH(4,M,NBH) = 0.D+0
              VAR_BH(6,M,NBH) = 0.D+0
            ENDIF
!
!---        Read injection temperature for inner pipe ---
!
            VARB = 'Inner Pipe Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(5,M,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &       UNTS(1:NCH),': ',VAR_BH(5,M,NBH)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR_BH(5,M,NBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(5,M,NBH),', C)'
!
!---        Simulation with solute transport  ---
!
            IF( IEQC.NE.0 ) THEN
              DO NSL_BH = 1,LSOLU_BH
!
!---            Solute active in borehole  ---
!
                IF( ISOLU_BH(NSL_BH,NBH).EQ.0 ) EXIT
                NSL = ISOLU_BH(NSL_BH,NBH)
!
!---            Read solute concentration in borehole fluid for
!               inner pipe  ---
!
                VARB = 'Starting Solute Concentration in Inner Pipe ' //
     &            'Fluid'
                NSLX = (NSL-1)*2 + 1
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VARC_BH(NSLX,M,NBH))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VARC_BH(NSLX,M,NBH)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VARC_BH(NSLX,M,NBH),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') ' (',VARC_BH(NSLX,M,NBH),
     &            ', 1/m^3)'
              ENDDO
            ENDIF
!
!---        Read injection mass rate for outer pipe ---
!
            IF( ITO_BHX.GE.1 .AND. ITO_BHX.LE.6 ) THEN
              VARB = 'Outer Pipe Injection Mass Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(7,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(7,M,NBH)
              IUNKG = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_BH(7,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(7,M,NBH),', kg/s)'
!
!---        Read injection volumetric rate for outer pipe ---
!
            ELSEIF( ITO_BHX.GE.11 .AND. ITO_BHX.LE.16 ) THEN
              VARB = 'Outer Pipe Injection Volumetric Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(7,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(7,M,NBH)
              IUNM = 3
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_BH(7,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(7,M,NBH),', kg/s)'
!
!---        Read pressure input for state condition #1 (saturated)
!           for outer pipe  ---
!
            ELSEIF( ITO_BHX.EQ.7 ) THEN
!
!---          Read outer pipe pressure  ---
!
              VARB = 'Outer Pipe Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(8,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(8,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(8,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(8,M,NBH),', Pa)'
!
!---          Nonisoair read dissolved air relative saturation  ---
!
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Outer Pipe Dissolved Air ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(9,M,NBH)
              ENDIF
!
!---          Nonisobrine read dissolved salt relative saturation  ---
!
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Outer Pipe Dissolved Salt ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(11,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(11,M,NBH)
              ENDIF
!
!---        Read pressure input for state condition #2
!           (partially saturated) for outer pipe ---
!
            ELSEIF( ITO_BHX.EQ.8 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Outer Pipe Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(8,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(8,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(8,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(8,M,NBH),', Pa)'
!
!---          Read aqueous saturation  ---
!
              VARB = 'Outer Pipe Aqueous Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(7,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(7,M,NBH)
              IF( VAR_BH(7,M,NBH).LT.EPSL .OR.
     &          (1.D+0-VAR_BH(7,M,NBH)).LT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Out-of-Range Outer Pipe Aqueous Saturation'
                RLMSG = VAR_BH(2,M,NBH)
                CALL WRMSGS( INDX )
              ENDIF
!
!---          Nonisobrine read dissolved salt relative saturation  ---
!
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Outer Pipe Dissolved Salt ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(11,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(11,M,NBH)
              ENDIF
!
!---        Read pressure input for state condition #3 (unsaturated)
!           for outer pipe   ---
!
            ELSEIF( ITO_BHX.EQ.9 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Outer Pipe Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(8,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(8,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(8,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(8,M,NBH),', Pa)'
!
!---          Read water vapor relative humidity  ---
!
              VARB = 'Outer Pipe Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(9,M,NBH)
            ENDIF
!
!---        Read water vapor concentration in gas for outer pipe ---
!
            IF( ITO_BHX.EQ.3 .OR. ITO_BHX.EQ.13 ) THEN
              VARB = 'Outer Pipe Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(9,M,NBH)
            ELSEIF( ITO_BHX.EQ.2 .OR. ITO_BHX.EQ.12 ) THEN
              VARB = 'Outer Pipe Water Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(9,M,NBH)
            ELSEIF( ITO_BHX.EQ.1 .OR. ITO_BHX.EQ.11 ) THEN
              VAR_BH(9,M,NBH) = 0.D+0
            ENDIF
!
!---        Read air concentration in aqueous for outer pipe ---
!
            IF( ITO_BHX.EQ.6 .OR. ITO_BHX.EQ.16 ) THEN
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Outer Pipe Dissolved Air ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(9,M,NBH)
              ENDIF
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Outer Pipe Dissolved Salt ' //
     &          'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(11,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(11,M,NBH)
              ENDIF
            ELSEIF( ITO_BHX.EQ.5 .OR. ITO_BHX.EQ.15 ) THEN
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Outer Pipe Dissolved Air Mass ' //
     &            'Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(9,M,NBH)
              ENDIF
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Outer Pipe Dissolved Salt Mass ' //
     &          'Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(11,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(11,M,NBH)
              ENDIF
            ELSEIF( ITO_BHX.EQ.4 .OR. ITO_BHX.EQ.14 ) THEN
              VAR_BH(9,M,NBH) = 0.D+0
              VAR_BH(11,M,NBH) = 0.D+0
            ENDIF
!
!---        Read injection temperature for outer pipe ---
!
            VARB = 'Outer Pipe Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(10,M,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &       UNTS(1:NCH),': ',VAR_BH(10,M,NBH)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR_BH(10,M,NBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(10,M,NBH),', C)'
!
!---        Simulation with solute transport  ---
!
            IF( IEQC.NE.0 ) THEN
              DO NSL_BH = 1,LSOLU_BH
!
!---            Solute active in borehole  ---
!
                IF( ISOLU_BH(NSL_BH,NBH).EQ.0 ) EXIT
                NSL = ISOLU_BH(NSL_BH,NBH)
!
!---            Read solute concentration in borehole fluid for
!               outer pipe  ---
!
                VARB = 'Starting Solute Concentration in Outer Pipe ' //
     &            'Fluid'
                NSLX = (NSL-1)*2 + 2
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VARC_BH(NSLX,M,NBH))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VARC_BH(NSLX,M,NBH)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VARC_BH(NSLX,M,NBH),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') ' (',VARC_BH(NSLX,M,NBH),
     &            ', 1/m^3)'
              ENDDO
            ENDIF
!
!---      Start-End type borehole  ---
!
          ELSEIF( IT_BH(2,NBH).NE.0 ) THEN
!
!---        Read injection mass rate for starting borehole node ---
!
            IF( IT_BH(1,NBH).GE.1 .AND. IT_BH(1,NBH).LE.6 ) THEN
              VARB = 'Starting Borehole Node Injection Mass Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(2,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(2,M,NBH)
              IUNKG = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_BH(2,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(2,M,NBH),', kg/s)'
!
!---        Read injection volumetric rate for starting borehole node ---
!
            ELSEIF( IT_BH(1,NBH).GE.11 .AND. IT_BH(1,NBH).LE.16 ) THEN
              VARB = 'Starting Borehole Node Injection Volumetric Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(2,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(2,M,NBH)
              IUNM = 3
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_BH(2,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(2,M,NBH),', kg/s)'
!
!---        Read pressure input for state condition #1 (saturated)
!           for starting borehole node  ---
!
            ELSEIF( IT_BH(1,NBH).EQ.7 ) THEN
!
!---        Read borehole pressure  ---
!
              VARB = 'Starting Borehole Node Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(3,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(3,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(3,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(3,M,NBH),', Pa)'
!
!---          Nonisoair read dissolved air relative saturation  ---
!
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Starting Borehole Node Dissolved Air ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(4,M,NBH)
              ENDIF
!
!---          Nonisobrine read dissolved salt relative saturation  ---
!
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Starting Borehole Node Dissolved Salt ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
!
!---        Read pressure input for state condition #2
!           (partially saturated) for starting borehole node ---
!
            ELSEIF( IT_BH(1,NBH).EQ.8 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Starting Borehole Node Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(3,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(3,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(3,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(3,M,NBH),', Pa)'
!
!---          Read aqueous saturation  ---
!
              VARB = 'Starting Borehole Node Aqueous Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(2,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(2,M,NBH)
              IF( VAR_BH(2,M,NBH).LT.EPSL .OR.
     &          (1.D+0-VAR_BH(2,M,NBH)).LT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Out-of-Range Borehole Aqueous Saturation'
                RLMSG = VAR_BH(2,M,NBH)
                CALL WRMSGS( INDX )
              ENDIF
!
!---          Nonisobrine read dissolved salt relative saturation  ---
!
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Starting Borehole Node Dissolved Salt ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
!
!---        Read pressure input for state condition #3 (unsaturated)
!           for starting borehole node   ---
!
            ELSEIF( IT_BH(1,NBH).EQ.9 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Starting Borehole Node Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(3,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(3,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(3,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(3,M,NBH),', Pa)'
!
!---          Read water vapor relative humidity  ---
!
              VARB = 'Starting Borehole Node Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(4,M,NBH)
            ENDIF
!
!---        Read water vapor concentration in gas for starting
!           borehole node ---
!
            IF( IT_BH(1,NBH).EQ.3 .OR. IT_BH(1,NBH).EQ.13 ) THEN
              VARB = 'Starting Borehole Node Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(4,M,NBH)
            ELSEIF( IT_BH(1,NBH).EQ.2 .OR. IT_BH(1,NBH).EQ.12 ) THEN
              VARB = 'Starting Borehole Node Water Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(4,M,NBH)
            ELSEIF( IT_BH(1,NBH).EQ.1 .OR. IT_BH(1,NBH).EQ.11 ) THEN
              VAR_BH(4,M,NBH) = 0.D+0
            ENDIF
!
!---        Read air concentration in aqueous for starting
!           borehole node---
!
            IF( IT_BH(1,NBH).EQ.6 .OR. IT_BH(1,NBH).EQ.16 ) THEN
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Starting Borehole Node Dissolved Air ' //
     &            'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(4,M,NBH)
              ENDIF
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Starting Borehole Node Dissolved Salt ' //
     &          'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
            ELSEIF( IT_BH(1,NBH).EQ.5 .OR. IT_BH(1,NBH).EQ.15 ) THEN
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Starting Borehole Node Dissolved Air Mass ' //
     &            'Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(4,M,NBH)
              ENDIF
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Starting Borehole Node Dissolved Salt Mass ' //
     &          'Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
            ELSEIF( IT_BH(1,NBH).EQ.4 .OR. IT_BH(1,NBH).EQ.14 ) THEN
              VAR_BH(4,M,NBH) = 0.D+0
              VAR_BH(6,M,NBH) = 0.D+0
            ENDIF
!
!---        Read injection temperature for starting borehole node ---
!
            IF( IT_BH(1,NBH).NE.21 .AND. IT_BH(1,NBH).NE.22 ) THEN
              VARB = 'Starting Borehole Node Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(5,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &         UNTS(1:NCH),': ',VAR_BH(5,M,NBH)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR_BH(5,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(5,M,NBH),', C)'
            ENDIF
!
!---        Simulation with solute transport  ---
!
            IF( IEQC.NE.0 ) THEN
              DO NSL_BH = 1,LSOLU_BH
!
!---            Solute active in borehole  ---
!
                IF( ISOLU_BH(NSL_BH,NBH).EQ.0 ) EXIT
                NSL = ISOLU_BH(NSL_BH,NBH)
!
!---            Read solute concentration in borehole fluid for
!               starting borehole node  ---
!
                VARB = 'Starting Solute Concentration in Borehole Fluid'
                NSLX = (NSL-1)*2 + 1
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VARC_BH(NSLX,M,NBH))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VARC_BH(NSLX,M,NBH)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VARC_BH(NSLX,M,NBH),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') ' (',VARC_BH(NSLX,M,NBH),
     &            ', 1/m^3)'
              ENDDO
            ENDIF
!
!---        Read ending borehole node information ---
!
!---        Read injection mass rate for ending borehole node ---
!
            IF( IT_BH(2,NBH).GE.1 .AND. IT_BH(2,NBH).LE.6 ) THEN
              VARB = 'Ending Borehole Node Injection Mass Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(7,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(7,M,NBH)
              IUNKG = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_BH(7,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(7,M,NBH),', kg/s)'
!
!---      Read injection volumetric rate for ending borehole node---
!
            ELSEIF( IT_BH(2,NBH).GE.11 .AND. IT_BH(2,NBH).LE.16 ) THEN
              VARB = 'Ending Borehole Node Injection Volumetric Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(7,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(7,M,NBH)
              IUNM = 3
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_BH(7,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(7,M,NBH),', kg/s)'
!
!---        Read pressure input for state condition #1 (saturated) for
!           ending borehole node ------
!
            ELSEIF( IT_BH(2,NBH).EQ.7 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Ending Borehole Node Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(8,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(8,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(8,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(8,M,NBH),', Pa)'
!
!---          Nonisoair read dissolved air relative saturation  ---
!
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Ending Borehole Node Dissolved Air Relative ' //
     &            'Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(9,M,NBH)
              ENDIF
!
!---          Nonisobrine read dissolved salt relative saturation  ---
!
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Ending Borehole Node Dissolved Salt ' //
     &                'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(11,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(11,M,NBH)
              ENDIF
!
!---        Read pressure input for state condition #2
!           (partially saturated) for ending borehole node  ---
!
            ELSEIF( IT_BH(2,NBH).EQ.8 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Ending Borehole Node Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(8,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(8,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(8,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(8,M,NBH),', Pa)'
!
!---          Read aqueous saturation  ---
!
              VARB = 'Ending Borehole Node Aqueous Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(7,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(7,M,NBH)
              IF( VAR_BH(7,M,NBH).LT.EPSL .OR.
     &          (1.D+0-VAR_BH(7,M,NBH)).LT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Out-of-Range Borehole Aqueous Saturation'
                RLMSG = VAR_BH(7,M,NBH)
                CALL WRMSGS( INDX )
              ENDIF
!
!---          Nonisobrine read dissolved salt relative saturation  ---
!
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Ending Borehole Node Dissolved Salt ' //
     &                'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(11,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(11,M,NBH)
              ENDIF
!
!---        Read pressure input for state condition #3 (unsaturated)
!           for ending borehole node  ------
!
            ELSEIF( IT_BH(2,NBH).EQ.9 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Ending Borehole Node Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(8,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(8,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(8,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(8,M,NBH),', Pa)'
!
!---          Read water vapor relative humidity  ---
!
              VARB = 'Ending Borehole Node Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(9,M,NBH)
!
!---        Read pressure input for state condition #3 (unsaturated)
!           for ending borehole node  ------
!
            ELSEIF( IT_BH(2,NBH).EQ.9 ) THEN
!
!---          Read pressure for ending borehole node  ---
!
              VARB = 'Ending Borehole Node Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(8,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(8,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(8,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(8,M,NBH),', Pa)'
!
!---          Read water vapor relative humidity or ending
!               borehole node ---
!
              VARB = 'Ending Borehole Node Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(9,M,NBH)
            ENDIF
!
!---        Read water vapor concentration in gas for ending
!             borehole node ---
!
            IF( IT_BH(2,NBH).EQ.3 .OR. IT_BH(2,NBH).EQ.13 ) THEN
              VARB = 'Ending Borehole Node Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(9,M,NBH)
            ELSEIF( IT_BH(2,NBH).EQ.2 .OR. IT_BH(2,NBH).EQ.12 ) THEN
              VARB = 'Ending Borehole Node Water Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(9,M,NBH)
            ELSEIF( IT_BH(2,NBH).EQ.1 .OR. IT_BH(2,NBH).EQ.11 ) THEN
              VAR_BH(9,M,NBH) = 0.D+0
            ENDIF
!
!---        Read air concentration in aqueous for ending
!             borehole node---
!
            IF( IT_BH(2,NBH).EQ.6 .OR. IT_BH(2,NBH).EQ.16 ) THEN
!
!---          Nonisoair read dissolved air relative saturation  ---
!
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Ending Borehole Node Dissolved Air Relative ' //
     &            'Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(9,M,NBH)
              ENDIF
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Ending Borehole Node Dissolved Salt ' //
     &          'Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(11,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(11,M,NBH)
              ENDIF
            ELSEIF( IT_BH(2,NBH).EQ.5 .OR. IT_BH(2,NBH).EQ.15 ) THEN
!
!---          Nonisoair read dissolved air relative saturation  ---
!
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Ending Borehole Node Dissolved Air Mass ' // 
     &            'Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(9,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(9,M,NBH)
              ENDIF
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Ending Borehole Node Dissolved Salt Mass ' //
     &            'Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(11,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(11,M,NBH)
              ENDIF
            ELSEIF( IT_BH(2,NBH).EQ.4 .OR. IT_BH(2,NBH).EQ.14 ) THEN
              VAR_BH(9,M,NBH) = 0.D+0
              VAR_BH(11,M,NBH) = 0.D+0
            ENDIF
!
!---        Read injection temperature for ending borehole node ---
!
            IF( IT_BH(2,NBH).NE.21 .AND. IT_BH(2,NBH).NE.22 ) THEN
              VARB = 'Ending Borehole Node Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(10,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &         UNTS(1:NCH),': ',VAR_BH(10,M,NBH)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR_BH(10,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(10,M,NBH),', C)'
            ENDIF
!
!---        Simulation with solute transport  ---
!
            IF( IEQC.NE.0 ) THEN
              DO NSL_BH = 1,LSOLU_BH
!
!---            Solute active in borehole  ---
!
                IF( ISOLU_BH(NSL_BH,NBH).EQ.0 ) EXIT
                NSL = ISOLU_BH(NSL_BH,NBH)
!
!---            Read solute concentration in borehole fluid for
!               starting borehole node  ---
!
                VARB = 'Ending Solute Concentration in Borehole Fluid'
                NSLX = (NSL-1)*2 + 2
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VARC_BH(NSLX,M,NBH))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &            UNTS(1:NCH),': ',VARC_BH(NSLX,M,NBH)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VARC_BH(NSLX,M,NBH),INDX)
                WRITE(IWR,'(A,1PE11.4,A)') ' (',VARC_BH(NSLX,M,NBH),
     &            ', 1/m^3)'
              ENDDO
            ENDIF
          ELSE
!
!---      Read injection mass rate for borehole ---
!
            IF( IT_BH(1,NBH).GE.1 .AND. IT_BH(1,NBH).LE.6 ) THEN
              VARB = 'Injection Mass Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(2,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(2,M,NBH)
              IUNKG = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_BH(2,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(2,M,NBH),', kg/s)'
!
!---      Read injection volumetric rate for borehole ---
!
            ELSEIF( IT_BH(1,NBH).GE.11 .AND. IT_BH(1,NBH).LE.16 ) THEN
              VARB = 'Injection Volumetric Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(2,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(2,M,NBH)
              IUNM = 3
              IUNS = -1
              CALL RDUNIT(UNTS,VAR_BH(2,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(2,M,NBH),', kg/s)'
!
!---        Read pressure input for state condition #1 (saturated)
!           for borehole  ---
!
            ELSEIF( IT_BH(1,NBH).EQ.7 ) THEN
!
!---        Read borehole pressure  ---
!
              VARB = 'Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(3,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(3,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(3,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(3,M,NBH),', Pa)'
!
!---          Nonisoair read dissolved air relative saturation  ---
!
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Dissolved Air Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(4,M,NBH)
              ENDIF
!
!---          Nonisobrine read dissolved salt relative saturation  ---
!
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Dissolved Salt Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
!
!---        Read pressure input for state condition #2
!           (partially saturated) for borehole ---
!
            ELSEIF( IT_BH(1,NBH).EQ.8 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(3,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(3,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(3,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(3,M,NBH),', Pa)'
!
!---          Read aqueous saturation  ---
!
              VARB = 'Aqueous Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(2,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(2,M,NBH)
              IF( VAR_BH(2,M,NBH).LT.EPSL .OR.
     &          (1.D+0-VAR_BH(2,M,NBH)).LT.EPSL ) THEN
                INDX = 9
                CHMSG = 'Out-of-Range Borehole Aqueous Saturation'
                RLMSG = VAR_BH(2,M,NBH)
                CALL WRMSGS( INDX )
              ENDIF
!
!---          Nonisobrine read dissolved salt relative saturation  ---
!
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Dissolved Salt Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
!
!---        Read pressure input for state condition #3 (unsaturated)
!           for borehole   ---
!
            ELSEIF( IT_BH(1,NBH).EQ.9 ) THEN
!
!---          Read borehole pressure  ---
!
              VARB = 'Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(3,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(3,M,NBH)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR_BH(3,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(3,M,NBH),', Pa)'
!
!---          Read water vapor relative humidity  ---
!
              VARB = 'Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(4,M,NBH)
            ENDIF
!
!---        Read water vapor concentration in gas  ---
!
            IF( IT_BH(1,NBH).EQ.3 .OR. IT_BH(1,NBH).EQ.13 ) THEN
              VARB = 'Water Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(4,M,NBH)
            ELSEIF( IT_BH(1,NBH).EQ.2 .OR. IT_BH(1,NBH).EQ.12 ) THEN
              VARB = 'Water Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
              WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR_BH(4,M,NBH)
            ELSEIF( IT_BH(1,NBH).EQ.1 .OR. IT_BH(1,NBH).EQ.11 ) THEN
              VAR_BH(4,M,NBH) = 0.D+0
            ENDIF
!
!---        Read air concentration in aqueous  ---
!
            IF( IT_BH(1,NBH).EQ.6 .OR. IT_BH(1,NBH).EQ.16 ) THEN
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Dissolved Air Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(4,M,NBH)
              ENDIF
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Dissolved Salt Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
            ELSEIF( IT_BH(1,NBH).EQ.5 .OR. IT_BH(1,NBH).EQ.15 ) THEN
              IF( ISLC(37).EQ.0 ) THEN
                VARB = 'Dissolved Air Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(4,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(4,M,NBH)
              ENDIF
              IF( ISLC(32).EQ.0 ) THEN
                VARB = 'Dissolved Salt Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
                WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR_BH(6,M,NBH)
              ENDIF
            ELSEIF( IT_BH(1,NBH).EQ.4 .OR. IT_BH(1,NBH).EQ.14 ) THEN
              VAR_BH(4,M,NBH) = 0.D+0
              VAR_BH(6,M,NBH) = 0.D+0
            ENDIF
!
!---        Read injection temperature for borehole ---
!
            IF( IT_BH(1,NBH).NE.21 .AND. IT_BH(1,NBH).NE.22 ) THEN
              VARB = 'Borehole Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(5,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &         UNTS(1:NCH),': ',VAR_BH(5,M,NBH)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR_BH(5,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(5,M,NBH),', C)'
            ENDIF
          ENDIF
!
!---      Nonisothermal, read injection temperature  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Dirichlet energy borehole  ---
!
            IF( IT_BH(1,NBH).EQ.21 ) THEN
              VARB = 'Dirichlet Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(5,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(5,M,NBH)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR_BH(5,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(5,M,NBH),', C)'
              VARB = 'Temperature Gradient along Borehole'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(6,M,NBH)
              INDX = 0
              IUNM = -1
              CALL RDUNIT(UNTS,VAR_BH(6,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(6,M,NBH),', 1/m)'
!
!---        Neumann energy borehole  ---
!
            ELSEIF( IT_BH(1,NBH).EQ.22 ) THEN
              VARB = 'Energy Flux'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(5,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(5,M,NBH)
              INDX = 0
              IUNKG = 1
              IUNS = -3
              CALL RDUNIT(UNTS,VAR_BH(5,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(5,M,NBH),', W/m^2)'
              VARB = 'Energy Flux Gradient along Borehole'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR_BH(6,M,NBH))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR_BH(6,M,NBH)
              INDX = 0
              IUNM = -1
              CALL RDUNIT(UNTS,VAR_BH(6,M,NBH),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR_BH(6,M,NBH),', 1/m)'
            ENDIF
          ENDIF
  300   CONTINUE
  400 CONTINUE
      WRITE(IWR,'(/)')
!
!---  Define borehole nodes, check borehole trajectory,
!     and write borehole.dat file  ---
!
      INDX = 0
      CALL CHK_BOREHOLE( INDX )
      CALL WR_BOREHOLE
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---      Set boundary condition type boreholes as inactive  ---
!
          IF( ( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) .OR.
     &      ( IT_BH(2,NBH).GE.21.AND.IT_BH(2,NBH).LE.29 ) )
     &       NXP_BH = NXP_BH+1
          INVX = INV_BH(NBN)
!
!---      Pipe flow  ---
!
          IF( IS_BH(INV_BH(NBN)).GT.10 ) THEN
            IZ_BH(NBN) = 0
          ELSE
            IZ_BH(NBN) = IZ_BHX(INVX)
          ENDIF
        ENDDO
      ENDDO
!
!---  Deallocate temporary memory for rock/soil index  ---
!
      IF( ALLOCATED(IZ_BHX) ) THEN
        DEALLOCATE( IZ_BHX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: IZ_BHX'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Open borehole preprocessed file for subsequent 
!     simulations  ---
!
      OPEN(UNIT=36, FILE='ppbh.dat', STATUS='UNKNOWN', 
     &  FORM='FORMATTED')
      CLOSE(UNIT=36,STATUS='DELETE')
      OPEN(UNIT=36, FILE='ppbh.dat', STATUS='NEW', 
     &  FORM='FORMATTED')
!
!---  Write borehole parameters  ---
!
      WRITE(36,*) LBC_BH,LBC_FRC,LBN_BH,LFC_BH,LI_BH,LN_BH,LT_BH,
     &    LSOLU_BH,LRC
!
!---  Find the range of borehole rock/soil type indices  ---
!
      I1X = LRC
      I2X = 0
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
          I1X = MIN( I1X,IZ_BH(NBN) )
          I2X = MAX( I2X,IZ_BH(NBN) )
        ENDDO
      ENDDO
!
!---  Write range of borehole rock/soil type indices  ---
!
      WRITE(36,*) I1X,I2X
!
!---  Write range of borehole rock/soil names  ---
!
      J1X = I1X
      DO IX = I1X,I2X,10
        J2X = MIN( J1X+9,I2X )
        WRITE(36,*) (TRIM(ROCK(JX)),JX=J1X,J2X)
        J1X = J1X + 10
      ENDDO
!
!---  Write number of boreholes and number of borehole nodes  ---
!
      WRITE(36,*) N_BH,NBN_BH
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Write borehole preprocessed file  ---
!
        CALL WRPREP_BH( NBH )
      ENDDO     
!
!---  Close borehole preprocessed file  ---
!
      CLOSE(UNIT=36)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RD_BH_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIC_BH_GT
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
!     Read input file for borehole initial conditions information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 23 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE GEOM_BH
      USE FDVP_BH
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VIC_BHX(4)
      CHARACTER*64 ADUM,BDUM,CDUM,DDUM,FDUM,UNTS
      CHARACTER*512 CHDUM
      LOGICAL FCHK,FBIN
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIC_BH_GT'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Borehole Initial Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
        IF( IEO.EQ.2 .AND. ISLC(87).EQ.0 ) THEN
          IC_BH(NBH) = 3
        ELSE
          IC_BH(NBH) = 0
        ENDIF
      ENDDO
!
!---  Read number of borehole inputs  ---
!
      WRITE(IWR,'(/,A)') 'Borehole Initial Condition Parameters'
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Borehole Initial Condition Inputs: '
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO NL = 1, NLIN
        ISTART = 1
        ICT_BHX = 0
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Borehole Number: '
        CALL RDINT(ISTART,ICOMMA,CHDUM,NBH)
        WRITE(IWR,'(2X,A,I9)') 'Borehole Number: ',NBH
        IF( NBH.LT.1 .OR. NBH.GT.N_BH ) THEN
          INDX = 7
          IMSG = NBH
          CHMSG = 'Out of Range Borehole Number'
          CALL WRMSGS( INDX )
        ENDIF
        VARB = 'Initial Borehole Saturation Option: '
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'overwrite').EQ.0 .AND. (IEO.EQ.2 .AND.
     &    ISLC(87).EQ.0) ) THEN
          CALL RDINPL( CHDUM )
          CYCLE
        ENDIF
        IF( INDEX(ADUM(1:),'hydrostatic').EQ.0 .AND. 
     &    INDEX(ADUM(1:),'borehole flow').EQ.0 .AND. 
     &    INDEX(ADUM(1:),'coaxial flow').EQ.0 )
     &    CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
        IF( INDEX(ADUM(1:),'hydrostatic').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Hydrostatic Borehole Initial Conditions'
          IC_BH(NBH) = 4
          VARB = 'Reference Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(1,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(1,NBH)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_BH(1,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(1,NBH),', Pa)'
          VIC_BH(1,NBH) = VIC_BH(1,NBH) - PATM
          VARB = 'Reference-Pressure Z Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(2,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(2,NBH)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VIC_BH(2,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(2,NBH),', m)'
!
!---      Nonisothermal conditions, read geothermal gradient inputs  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            VARB = 'Reference Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(7,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(7,NBH)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VIC_BH(7,NBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(7,NBH),', C)'
            VARB = 'Reference-Temperature Z Location'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(8,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(8,NBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VIC_BH(8,NBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(8,NBH),', m)'
            VARB = 'Geothermal Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(9,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(9,NBH)
            INDX = 0
            IUNK = 1
            IUNM = -1
            CALL RDUNIT(UNTS,VIC_BH(9,NBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(9,NBH),', C/m)'
          ENDIF
!
!---      Nonisobrine conditions, read salt conc. gradient inputs  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            VARB = 'Reference Salt Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(10,NBH))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        VIC_BH(10,NBH)
            VARB = 'Reference-Salinity (Salt Mass Fraction) Z Location'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(11,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(11,NBH)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VIC_BH(11,NBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(11,NBH),', m)'
            VARB = 'Geosalinity (Salt Mass Fraction) Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(12,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(12,NBH)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,VIC_BH(12,NBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(12,NBH),', 1/m)'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
          IF( INDEX(BDUM(1:),'gas pressure').NE.0 ) THEN
            IC_BH(NBH) = 1
          ELSEIF( INDEX(BDUM(1:),'aqueous pressure').NE.0 ) THEN
            IC_BH(NBH) = 2
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Initial Borehole Saturation Option: '
     &        // ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas pressure').NE.0 ) THEN
          IF( INDEX(BDUM(1:),'aqueous saturation').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Pressure-Aqueous Saturation ' //
     &        'Borehole Initial Conditions'
            IC_BH(NBH) = 1
          ELSEIF( INDEX(BDUM(1:),'aqueous pressure').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Pressure-Aqueous Pressure '
     &        // 'Borehole Initial Conditions'
            IC_BH(NBH) = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &        ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
          IF( INDEX(BDUM(1:),'aqueous saturation').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous Pressure-Aqueous Saturation '
     &        // 'Borehole Initial Conditions'
            IC_BH(NBH) = 2
          ELSEIF( INDEX(BDUM(1:),'gas pressure').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Pressure-Aqueous Pressure '
     &        // 'Borehole Initial Conditions'
            IC_BH(NBH) = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &        ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'activity').NE.0 ) THEN
            ICT_BHX = 4
          ELSEIF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
            ICT_BHX = 3
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
            ICT_BHX = 2
          ELSE
            ICT_BHX = 1
          ENDIF
          DO NSL = 1,NSOLU
            IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SOLUT(NSL)(1:IDB) ) EXIT
          ENDDO
          IF( NSL.GT.NSOLU ) THEN
            INDX = 7
            IMSG = NBH
            CHMSG = 'Unrecognized Solute: ' // BDUM(1:NCHB) //
     &        ' at Borehole'
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'coaxial flow').NE.0 ) THEN
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          CALL RDCHR(ISTART,ICOMMA,NCHC,CHDUM,CDUM)
          IF( INDEX(BDUM(1:),'aqueous mass inject').NE.0 .AND. 
     &      INDEX(CDUM(1:),'aqueous pressure').NE.0 ) THEN
            IC_BH(NBH) = 11
          ELSEIF( INDEX(CDUM(1:),'aqueous mass inject').NE.0 .AND. 
     &      INDEX(BDUM(1:),'aqueous pressure').NE.0 ) THEN
            IC_BH(NBH) = 12
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Initial Condition Option: ' //
     &        ADUM(1:NCHA) // ',' // BDUM(1:NCHB) // ',' // CDUM(1:NCHC)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'borehole flow').NE.0 ) THEN
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          CALL RDCHR(ISTART,ICOMMA,NCHC,CHDUM,CDUM)
          IF( INDEX(BDUM(1:),'aqueous mass inject').NE.0 .AND. 
     &      INDEX(CDUM(1:),'aqueous pressure').NE.0 ) THEN
            IC_BH(NBH) = 13
          ELSEIF( INDEX(CDUM(1:),'aqueous mass inject').NE.0 .AND. 
     &      INDEX(BDUM(1:),'aqueous pressure').NE.0 ) THEN
            IC_BH(NBH) = 14
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Initial Condition Option: ' //
     &        ADUM(1:NCHA) // ',' // BDUM(1:NCHB) // ',' // CDUM(1:NCHC)
            CALL WRMSGS( INDX )
          ENDIF
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Option: ' //
     &      ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ENDIF
        VARB = 'Reference X Location'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(1,NBH))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(1,NBH)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VIC_BH(1,NBH),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(1,NBH),', m)'
        VARB = 'Reference Y Location'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(2,NBH))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(2,NBH)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VIC_BH(2,NBH),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(2,NBH),', m)'
        VARB = 'Reference Z Location'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(3,NBH))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(3,NBH)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VIC_BH(3,NBH),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(3,NBH),', m)'
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        IF( IC_BH(NBH).EQ.1 .AND. ICT_BHX.EQ.0 ) THEN
          VARB = 'Reference Gas Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(4,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(4,NBH)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_BH(4,NBH),INDX)
          VARX = 1.D+0
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(4,NBH),', Pa)'
          VIC_BH(4,NBH) = VIC_BH(4,NBH) - PATM
          VARB = 'X-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(5,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(5,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(5,NBH) = VARX*VIC_BH(5,NBH)
          CALL RDUNIT(UNTS,VIC_BH(5,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(5,NBH),', Pa/m)'
          VARB = 'Y-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(6,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(6,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(6,NBH) = VARX*VIC_BH(6,NBH)
          CALL RDUNIT(UNTS,VIC_BH(6,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(6,NBH),', Pa/m)'
          VARB = 'Z-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(7,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(7,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(7,NBH) = VARX*VIC_BH(7,NBH)
          CALL RDUNIT(UNTS,VIC_BH(7,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(7,NBH),', Pa/m)'
          VARB = 'Reference Aqueous Saturation'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(8,NBH))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),', ',VIC_BH(8,NBH)
          VARB = 'X-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(9,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(9,NBH)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_BH(9,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(9,NBH),', 1/m)'
          VARB = 'Y-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(10,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(10,NBH)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_BH(10,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(10,NBH),', 1/m)'
          VARB = 'Z-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(11,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(11,NBH)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_BH(11,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(11,NBH),', 1/m)'
        ELSEIF( IC_BH(NBH).EQ.2 .AND. ICT_BHX.EQ.0 ) THEN
          VARB = 'Reference Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(4,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(4,NBH)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_BH(4,NBH),INDX)
          VARX = 1.D+0
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(4,NBH),', Pa)'
          VIC_BH(4,NBH) = VIC_BH(4,NBH) - PATM
          VARB = 'X-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(5,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(5,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(5,NBH) = VARX*VIC_BH(5,NBH)
          CALL RDUNIT(UNTS,VIC_BH(5,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(5,NBH),', Pa/m)'
          VARB = 'Y-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(6,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(6,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(6,NBH) = VARX*VIC_BH(6,NBH)
          CALL RDUNIT(UNTS,VIC_BH(6,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(6,NBH),', Pa/m)'
          VARB = 'Z-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(7,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(7,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(7,NBH) = VARX*VIC_BH(7,NBH)
          CALL RDUNIT(UNTS,VIC_BH(7,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(7,NBH),', Pa/m)'
          VARB = 'Reference Aqueous Saturation'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(8,NBH))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),', ',VIC_BH(8,NBH)
          VARB = 'X-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(9,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(9,NBH)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_BH(9,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(9,NBH),', 1/m)'
          VARB = 'Y-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(10,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(10,NBH)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_BH(10,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(10,NBH),', 1/m)'
          VARB = 'Z-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(11,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(11,NBH)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_BH(11,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(11,NBH),', 1/m)'
        ELSEIF( IC_BH(NBH).EQ.3 .AND. ICT_BHX.EQ.0 ) THEN
          VARB = 'Reference Gas Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(4,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(4,NBH)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_BH(4,NBH),INDX)
          VARX = 1.D+0
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(4,NBH),', Pa)'
          VIC_BH(4,NBH) = VIC_BH(4,NBH) - PATM
          VARB = 'X-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(5,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(5,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(5,NBH) = VARX*VIC_BH(5,NBH)
          CALL RDUNIT(UNTS,VIC_BH(5,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(5,NBH),', Pa/m)'
          VARB = 'Y-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(6,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(6,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(6,NBH) = VARX*VIC_BH(6,NBH)
          CALL RDUNIT(UNTS,VIC_BH(6,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(6,NBH),', Pa/m)'
          VARB = 'Z-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(7,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(7,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(7,NBH) = VARX*VIC_BH(7,NBH)
          CALL RDUNIT(UNTS,VIC_BH(7,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(7,NBH),', Pa/m)'
          VARB = 'Reference Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(8,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(8,NBH)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_BH(8,NBH),INDX)
          VARX = 1.D+0
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(8,NBH),', Pa)'
          VIC_BH(8,NBH) = VIC_BH(8,NBH) - PATM
          VARB = 'X-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(9,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(9,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(9,NBH) = VARX*VIC_BH(9,NBH)
          CALL RDUNIT(UNTS,VIC_BH(9,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(9,NBH),', Pa/m)'
          VARB = 'Y-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(10,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(10,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(10,NBH) = VARX*VIC_BH(10,NBH)
          CALL RDUNIT(UNTS,VIC_BH(10,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(10,NBH),', Pa/m)'
          VARB = 'Z-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(11,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(11,NBH)
          INDX = 0
          IUNM = -1
          VIC_BH(11,NBH) = VARX*VIC_BH(11,NBH)
          CALL RDUNIT(UNTS,VIC_BH(11,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(11,NBH),', Pa/m)'
        ELSEIF( IC_BH(NBH).EQ.11 ) THEN
          VARB = 'Inner Aqueous Mass Injection'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(5,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(5,NBH)
          INDX = 0
          IUNKG = 1
          IUNS = -1
          CALL RDUNIT(UNTS,VIC_BH(5,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(5,NBH),', kg/s)'
          VARB = 'Outer Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(4,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(4,NBH)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_BH(4,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(4,NBH),', Pa)'
        ELSEIF( IC_BH(NBH).EQ.12 ) THEN
          VARB = 'Inner Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(4,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(4,NBH)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_BH(4,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(4,NBH),', Pa)'
          VARB = 'Outer Aqueous Mass Injection'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(5,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(5,NBH)
          INDX = 0
          IUNKG = 1
          IUNS = -1
          CALL RDUNIT(UNTS,VIC_BH(5,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(5,NBH),', kg/s)'
        ELSEIF( IC_BH(NBH).EQ.13 ) THEN
          VARB = 'Start Aqueous Mass Injection'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(5,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(5,NBH)
          INDX = 0
          IUNKG = 1
          IUNS = -1
          CALL RDUNIT(UNTS,VIC_BH(5,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(5,NBH),', kg/s)'
          VARB = 'End Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(4,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(4,NBH)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_BH(4,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(4,NBH),', Pa)'
        ELSEIF( IC_BH(NBH).EQ.14 ) THEN
          VARB = 'Start Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(4,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(4,NBH)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_BH(4,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(4,NBH),', Pa)'
          VARB = 'End Aqueous Mass Injection'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(5,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(5,NBH)
          INDX = 0
          IUNKG = 1
          IUNS = -1
          CALL RDUNIT(UNTS,VIC_BH(5,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(5,NBH),', kg/s)'
        ENDIF
!
!---    Nonisothermal conditions, read geothermal gradient inputs  ---
!
        IF( ISLC(30).EQ.0 .AND. ICT_BHX.EQ.0 ) THEN
!
!---      Check for external file  ---
!
          ISTX = ISTART
          ICMX = ICOMMA
          IDFLT = 0
          VARB = 'External Borehole Temperature File Name'
          CALL RDCHR(ISTART,ICOMMA,NCHD,CHDUM,DDUM)
          FBIN = .FALSE.
          IF( INDEX(DDUM(1:),'binary').NE.0 .OR.
     &      INDEX(DDUM(1:),'bfile').NE.0 .OR.
     &      INDEX(DDUM(1:),'b_file').NE.0 ) FBIN = .TRUE.
          IF( INDEX(DDUM(1:),'file').NE.0 ) THEN
            ICOLON = INDEX(DDUM(1:),':') + 1
            FDUM = DDUM(ICOLON:NCHD)
            NCHF = NCHD-ICOLON+1
            INQUIRE( FILE=FDUM(1:NCHF), FORM=CDUM, EXIST=FCHK )
            IF( .NOT.FCHK ) THEN
              INDX = 4
              CHMSG = 'External Borehole Temperature File does ' // 
     &          'not Exist: ' // FDUM(1:NCHF)
              CALL WRMSGS( INDX )
            ELSEIF( CDUM.EQ.'UNFORMATTED' .AND. (.NOT. FBIN) ) THEN
              INDX = 4
              CHMSG = 'External Borehole Temperature File is ' // 
     &          'Unformatted: ' // FDUM(1:NCHF)
              CALL WRMSGS( INDX )
            ELSEIF( CDUM.EQ.'FORMATTED' .AND. FBIN ) THEN
              INDX = 4
              CHMSG = 'External Borehole Temperature File is ' //
     &          'Formatted: ' // FDUM(1:NCHF)
              CALL WRMSGS( INDX )
            ENDIF
            IF( FBIN ) THEN
              OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &          FORM='UNFORMATTED')
            ELSE
              OPEN(UNIT=26, FILE=FDUM(1:NCHF), STATUS='OLD',
     &          FORM='FORMATTED')
            ENDIF
            WRITE(IWR,'(/,2A)') 'External Borehole Temperature File: ',
     &        FDUM(1:NCHF)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,3A)') 'External Borehole Temperature ' // 
     &        'File Units: ',UNTS(1:NCH)
            VARX = 1.D+0
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VARX,INDX)
            IF( FBIN ) THEN
              READ(26) (T_BH(2,NBN),NBN=ID_BH(3,NBH),ID_BH(4,NBH))
            ELSE
              READ(26,*) (T_BH(2,NBN),NBN=ID_BH(3,NBH),ID_BH(4,NBH))
            ENDIF
            DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
              T_BH(2,NBN) = VARX*T_BH(2,NBN)
            ENDDO
            VIC_BH(12,NBH) = -1.D+3
            CLOSE(26)
          ELSE
            ISTART = ISTX
            ICOMMA = ICMX         
            VARB = 'Reference Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(12,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(12,NBH)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VIC_BH(12,NBH),INDX)
            VARX = 1.D+0
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VARX,INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(12,NBH),', C)'
            VARB = 'X-Dir. Geothermal Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(13,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(13,NBH)
            INDX = 0
            IUNM = -1
            VIC_BH(13,NBH) = VARX*VIC_BH(13,NBH)
            CALL RDUNIT(UNTS,VIC_BH(13,NBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(13,NBH),', C/m)'
            VARB = 'Y-Dir. Geothermal Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(14,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(14,NBH)
            INDX = 0
            IUNM = -1
            VIC_BH(14,NBH) = VARX*VIC_BH(14,NBH)
            CALL RDUNIT(UNTS,VIC_BH(14,NBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(14,NBH),', C/m)'
            VARB = 'Z-Dir. Geothermal Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(15,NBH))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(15,NBH)
            INDX = 0
            IUNM = -1
            VIC_BH(15,NBH) = VARX*VIC_BH(15,NBH)
            CALL RDUNIT(UNTS,VIC_BH(15,NBH),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(15,NBH),', C/m)'
          ENDIF
        ENDIF
!
!---    Nonisobrine conditions, read salt conc. gradient inputs  ---
!
        IF( ISLC(32).EQ.0 .AND. ICT_BHX.EQ.0 ) THEN
          VARB = 'Reference Salt Mass Fraction'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(16,NBH))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &      VIC_BH(10,NBH)
          VARB = 'X-Dir. Geosalinity (Salt Mass Fraction) Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(17,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(17,NBH)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_BH(17,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(17,NBH),', 1/m)'
          VARB = 'Y-Dir. Geosalinity (Salt Mass Fraction) Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(18,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(18,NBH)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_BH(18,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(18,NBH),', 1/m)'
          VARB = 'Z-Dir. Geosalinity (Salt Mass Fraction) Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(19,NBH))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BH(19,NBH)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_BH(19,NBH),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BH(19,NBH),', 1/m)'
        ENDIF
        IF( ICT_BHX.EQ.0 ) THEN
          CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
          IF( INDX.EQ.1 ) THEN
            VARB = 'Dissolved Air Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(20,NBH))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        VIC_BH(20,NBH)
            VARB = 'Water Vapor Relative Humidity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BH(21,NBH))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        VIC_BH(21,NBH)
          ENDIF
        ENDIF
!
!---    Read initial borehole node solute concentrations  ---
!
        IF( ICT_BHX.NE.0 ) THEN
          IF( ICT_BHX.EQ.4 ) THEN
            VARB = 'Reference Solute Activity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BHX(1))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BHX(1)
            INDX = 0
            CALL RDUNIT(UNTS,VIC_BHX(1),INDX)
            INDX = 0
            VARX = 1.D+0
            CALL RDUNIT(UNTS,VARX,INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BHX(1),', Bq)'
            VARB = 'X-Dir. Solute Activity Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BHX(2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BHX(2)
            INDX = 0
            IUNM = -1
            VIC_BHX(2) = VARX*VIC_BHX(2)
            CALL RDUNIT(UNTS,VIC_BHX(2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BHX(2),
     &        ', Bq/m)'
            VARB = 'Y-Dir. Solute Activity Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BHX(3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BHX(3)
            INDX = 0
            IUNM = -1
            VIC_BHX(3) = VARX*VIC_BHX(3)
            CALL RDUNIT(UNTS,VIC_BHX(3),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BHX(3),
     &        ', Bq/m)'
            VARB = 'Z-Dir. Solute Activity Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BHX(4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BHX(4)
            INDX = 0
            IUNM = -1
            VIC_BHX(4) = VARX*VIC_BHX(4)
            CALL RDUNIT(UNTS,VIC_BHX(4),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BHX(4),
     &        ', Bq/m)'
          ELSE
            VARB = 'Reference Solute Concentration'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BHX(1))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BHX(1)
            INDX = 0
            IUNM = -3
            CALL RDUNIT(UNTS,VIC_BHX(1),INDX)
            INDX = 0
            IUNM = -3
            VARX = 1.D+0
            CALL RDUNIT(UNTS,VARX,INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BHX(1),', 1/m^3)'
            VARB = 'X-Dir. Solute Concentration Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BHX(2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BHX(2)
            INDX = 0
            IUNM = -1
            VIC_BHX(2) = VARX*VIC_BHX(2)
            CALL RDUNIT(UNTS,VIC_BHX(2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BHX(2),
     &        ', (1/m^3)/m)'
            VARB = 'Y-Dir. Solute Concentration Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BHX(3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BHX(3)
            INDX = 0
            IUNM = -1
            VIC_BHX(3) = VARX*VIC_BHX(3)
            CALL RDUNIT(UNTS,VIC_BHX(3),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BHX(3),
     &        ', (1/m^3)/m)'
            VARB = 'Z-Dir. Solute Concentration Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_BHX(4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_BHX(4)
            INDX = 0
            IUNM = -1
            VIC_BHX(4) = VARX*VIC_BHX(4)
            CALL RDUNIT(UNTS,VIC_BHX(4),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_BHX(4),
     &        ', (1/m^3)/m)'
          ENDIF
!
!---      Loop over borehole nodes, initializing solute
!         concentrations in borehole nodes  ---
!
          DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
            XP_BHX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
            YP_BHX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
            ZP_BHX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
            C_BH(NBN,NSL) = VIC_BHX(1) +
     &        (XP_BHX-VIC_BH(1,NBH))*VIC_BHX(2) +
     &        (YP_BHX-VIC_BH(2,NBH))*VIC_BHX(3) +
     &        (ZP_BHX-VIC_BH(3,NBH))*VIC_BHX(4)
          ENDDO
          ICT_BH(NBH,NSL) = ICT_BHX
        ENDIF
      ENDDO
!
!---  Loop over boreholes, checking for uninitialized
!     hydrologic state conditions for non-restart simulations  ---
!
      IF( IEO.NE.2 .AND. ISLC(87).EQ.0 ) THEN
        DO NBH = 1,N_BH
!
!---      Skip for boundary condition type boreholes  ---
!
          IF( ( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) .OR.
     &      ( IT_BH(2,NBH).GE.21.AND.IT_BH(2,NBH).LE.29 ) ) CYCLE
          IF( IC_BH(NBH).EQ.0 ) THEN
            INDX = 7
            IMSG = NBH
            CHMSG = 'Uninitialized Hydrologic State for Borehole'
            CALL WRMSGS( INDX )
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIC_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIC_FRC_GT
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
!     Read input file for fracture initial conditions information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 5 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE PARM_FRC
      USE GEOM_FRC
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VIC_FRCX(4)
      CHARACTER*64 ADUM,BDUM,UNTS
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIC_FRC_GT'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Fracture Initial Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
        IF( IEO.EQ.2 .AND. ISLC(87).EQ.0 ) THEN
          IC_FRC(NFX) = 3
        ELSE
          IC_FRC(NFX) = 0
        ENDIF
      ENDDO
!
!---  Read number of fracture inputs, unspecified fractures will
!     be in thermodynamic equilibrium with connecting nodes  ---
!
      WRITE(IWR,'(/,A)') 'Fracture Initial Condition Parameters'
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Fracture Initial Condition Inputs: '
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO NL = 1, NLIN
        ISTART = 1
        ICT_FRCX = 0
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Fracture Number: '
        CALL RDINT(ISTART,ICOMMA,CHDUM,NFX)
        WRITE(IWR,'(2X,A,I9)') 'Fracture Number: ',NFX
        IF( NFX.LT.1 .OR. NFX.GT.NF_FRC ) THEN
          INDX = 7
          IMSG = NFX
          CHMSG = 'Out of Range Fracture Number'
          CALL WRMSGS( INDX )
        ENDIF
        VARB = 'Initial Fracture Saturation Option: '
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'overwrite').EQ.0 .AND. (IEO.EQ.2
     &     .AND. ISLC(87).EQ.0) ) THEN
          CALL RDINPL( CHDUM )
          CYCLE
        ENDIF
        IF( INDEX(ADUM(1:),'hydrostatic').EQ.0 )
     &    CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
        IF( INDEX(ADUM(1:),'hydrostatic').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Hydrostatic Fracture Initial Conditions'
          IC_FRC(NFX) = 4
          VARB = 'Reference Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(1,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(1,NFX)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_FRC(1,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(1,NFX),', Pa)'
          VIC_FRC(1,NFX) = VIC_FRC(1,NFX) - PATM
          VARB = 'Reference-Pressure Z Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(2,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(2,NFX)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VIC_FRC(2,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(2,NFX),', m)'
!
!---      Nonisothermal conditions, read geothermal gradient inputs  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            VARB = 'Reference Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(7,NFX))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(7,NFX)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VIC_FRC(7,NFX),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(7,NFX),', C)'
            VARB = 'Reference-Temperature Z Location'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(8,NFX))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(8,NFX)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VIC_FRC(8,NFX),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(8,NFX),', m)'
            VARB = 'Geothermal Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(9,NFX))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(9,NFX)
            INDX = 0
            IUNK = 1
            IUNM = -1
            CALL RDUNIT(UNTS,VIC_FRC(9,NFX),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(9,NFX),', C/m)'
          ENDIF
!
!---      Nonisobrine conditions, read salt conc. gradient inputs  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            VARB = 'Reference Salt Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(10,NFX))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        VIC_FRC(10,NFX)
            VARB = 'Reference-Salinity (Salt Mass Fraction) Z Location'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(11,NFX))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(11,NFX)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VIC_FRC(11,NFX),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(11,NFX),', m)'
            VARB = 'Geosalinity (Salt Mass Fraction) Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(12,NFX))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(12,NFX)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,VIC_FRC(12,NFX),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(12,NFX),', 1/m)'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
          IF( INDEX(BDUM(1:),'gas pressure').NE.0 ) THEN
            IC_FRC(NFX) = 1
          ELSEIF( INDEX(BDUM(1:),'aqueous pressure').NE.0 ) THEN
            IC_FRC(NFX) = 2
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Initial Fracture Saturation Option: '
     &        // ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas pressure').NE.0 ) THEN
          IF( INDEX(BDUM(1:),'aqueous saturation').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Pressure-Aqueous Saturation ' //
     &        'Fracture Initial Conditions'
            IC_FRC(NFX) = 1
          ELSEIF( INDEX(BDUM(1:),'aqueous pressure').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Pressure-Aqueous Pressure '
     &        // 'Fracture Initial Conditions'
            IC_FRC(NFX) = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &        ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
          IF( INDEX(BDUM(1:),'aqueous saturation').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous Pressure-Aqueous Saturation '
     &        // 'Fracture Initial Conditions'
            IC_FRC(NFX) = 2
          ELSEIF( INDEX(BDUM(1:),'gas pressure').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Pressure-Aqueous Pressure '
     &        // 'Fracture Initial Conditions'
            IC_FRC(NFX) = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &        ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'activity').NE.0 ) THEN
            ICT_FRCX = 4
          ELSEIF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
            ICT_FRCX = 3
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
            ICT_FRCX = 2
          ELSE
            ICT_FRCX = 1
          ENDIF
          DO NSL = 1,NSOLU
            IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SOLUT(NSL)(1:IDB) ) EXIT
          ENDDO
          IF( NSL.GT.NSOLU ) THEN
            INDX = 7
            IMSG = NFX
            CHMSG = 'Unrecognized Solute: ' // BDUM(1:NCHB) //
     &        ' at Fracture'
            CALL WRMSGS( INDX )
          ENDIF
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Option: ' //
     &      ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ENDIF
        VARB = 'Reference X Location'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(1,NFX))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(1,NFX)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VIC_FRC(1,NFX),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(1,NFX),', m)'
        VARB = 'Reference Y Location'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(2,NFX))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(2,NFX)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VIC_FRC(2,NFX),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(2,NFX),', m)'
        VARB = 'Reference Z Location'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(3,NFX))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(3,NFX)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VIC_FRC(3,NFX),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(3,NFX),', m)'
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        IF( IC_FRC(NFX).EQ.1 .AND. ICT_FRCX.EQ.0 ) THEN
          VARB = 'Reference Gas Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(4,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(4,NFX)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_FRC(4,NFX),INDX)
          VARX = 1.D+0
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(4,NFX),', Pa)'
          VIC_FRC(4,NFX) = VIC_FRC(4,NFX) - PATM
          VARB = 'X-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(5,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(5,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(5,NFX) = VARX*VIC_FRC(5,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(5,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(5,NFX),', Pa/m)'
          VARB = 'Y-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(6,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(6,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(6,NFX) = VARX*VIC_FRC(6,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(6,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(6,NFX),', Pa/m)'
          VARB = 'Z-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(7,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(7,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(7,NFX) = VARX*VIC_FRC(7,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(7,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(7,NFX),', Pa/m)'
          VARB = 'Reference Aqueous Saturation'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(8,NFX))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),', ',VIC_FRC(8,NFX)
          VARB = 'X-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(9,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(9,NFX)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_FRC(9,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(9,NFX),', 1/m)'
          VARB = 'Y-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(10,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(10,NFX)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_FRC(10,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(10,NFX),', 1/m)'
          VARB = 'Z-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(11,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(11,NFX)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_FRC(11,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(11,NFX),', 1/m)'
        ELSEIF( IC_FRC(NFX).EQ.2 .AND. ICT_FRCX.EQ.0 ) THEN
          VARB = 'Reference Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(4,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(4,NFX)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_FRC(4,NFX),INDX)
          VARX = 1.D+0
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(4,NFX),', Pa)'
          VIC_FRC(4,NFX) = VIC_FRC(4,NFX) - PATM
          VARB = 'X-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(5,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(5,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(5,NFX) = VARX*VIC_FRC(5,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(5,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(5,NFX),', Pa/m)'
          VARB = 'Y-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(6,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(6,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(6,NFX) = VARX*VIC_FRC(6,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(6,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(6,NFX),', Pa/m)'
          VARB = 'Z-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(7,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(7,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(7,NFX) = VARX*VIC_FRC(7,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(7,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(7,NFX),', Pa/m)'
          VARB = 'Reference Aqueous Saturation'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(8,NFX))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),', ',VIC_FRC(8,NFX)
          VARB = 'X-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(9,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(9,NFX)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_FRC(9,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(9,NFX),', 1/m)'
          VARB = 'Y-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(10,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(10,NFX)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_FRC(10,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(10,NFX),', 1/m)'
          VARB = 'Z-Dir. Aqueous Saturation Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(11,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(11,NFX)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_FRC(11,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(11,NFX),', 1/m)'
        ELSEIF( IC_FRC(NFX).EQ.3 .AND. ICT_FRCX.EQ.0 ) THEN
          VARB = 'Reference Gas Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(4,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(4,NFX)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_FRC(4,NFX),INDX)
          VARX = 1.D+0
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(4,NFX),', Pa)'
          VIC_FRC(4,NFX) = VIC_FRC(4,NFX) - PATM
          VARB = 'X-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(5,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(5,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(5,NFX) = VARX*VIC_FRC(5,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(5,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(5,NFX),', Pa/m)'
          VARB = 'Y-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(6,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(6,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(6,NFX) = VARX*VIC_FRC(6,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(6,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(6,NFX),', Pa/m)'
          VARB = 'Z-Dir. Gas Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(7,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(7,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(7,NFX) = VARX*VIC_FRC(7,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(7,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(7,NFX),', Pa/m)'
          VARB = 'Reference Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(8,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(8,NFX)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VIC_FRC(8,NFX),INDX)
          VARX = 1.D+0
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(8,NFX),', Pa)'
          VIC_FRC(8,NFX) = VIC_FRC(8,NFX) - PATM
          VARB = 'X-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(9,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(9,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(9,NFX) = VARX*VIC_FRC(9,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(9,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(9,NFX),', Pa/m)'
          VARB = 'Y-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(10,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(10,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(10,NFX) = VARX*VIC_FRC(10,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(10,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(10,NFX),', Pa/m)'
          VARB = 'Z-Dir. Aqueous Pressure Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(11,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(11,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(11,NFX) = VARX*VIC_FRC(11,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(11,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(11,NFX),', Pa/m)'
        ENDIF
!
!---    Nonisothermal conditions, read geothermal gradient inputs  ---
!
        IF( ISLC(30).EQ.0 .AND. ICT_FRCX.EQ.0 ) THEN
          VARB = 'Reference Temperature'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(12,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(12,NFX)
          INDX = 0
          IUNK = 1
          CALL RDUNIT(UNTS,VIC_FRC(12,NFX),INDX)
          VARX = 1.D+0
          INDX = 0
          IUNK = 1
          CALL RDUNIT(UNTS,VARX,INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(12,NFX),', C)'
          VARB = 'X-Dir. Geothermal Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(13,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(13,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(13,NFX) = VARX*VIC_FRC(13,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(13,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(13,NFX),', C/m)'
          VARB = 'Y-Dir. Geothermal Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(14,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(14,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(14,NFX) = VARX*VIC_FRC(14,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(14,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(14,NFX),', C/m)'
          VARB = 'Z-Dir. Geothermal Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(15,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(15,NFX)
          INDX = 0
          IUNM = -1
          VIC_FRC(15,NFX) = VARX*VIC_FRC(15,NFX)
          CALL RDUNIT(UNTS,VIC_FRC(15,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(15,NFX),', C/m)'
        ENDIF
!
!---    Nonisobrine conditions, read salt conc. gradient inputs  ---
!
        IF( ISLC(32).EQ.0 .AND. ICT_FRCX.EQ.0 ) THEN
          VARB = 'Reference Salt Mass Fraction'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(16,NFX))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &      VIC_FRC(10,NFX)
          VARB = 'X-Dir. Geosalinity (Salt Mass Fraction) Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(17,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(17,NFX)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_FRC(17,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(17,NFX),', 1/m)'
          VARB = 'Y-Dir. Geosalinity (Salt Mass Fraction) Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(18,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(18,NFX)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_FRC(18,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(18,NFX),', 1/m)'
          VARB = 'Z-Dir. Geosalinity (Salt Mass Fraction) Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(19,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(19,NFX)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,VIC_FRC(19,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(19,NFX),', 1/m)'
        ENDIF
        IF( ICT_FRCX.EQ.0 ) THEN
          CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
          IF( INDX.EQ.1 ) THEN
            VARB = 'Dissolved Air Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(20,NFX))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        VIC_FRC(20,NFX)
            VARB = 'Water Vapor Relative Humidity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(21,NFX))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &        VIC_FRC(21,NFX)
          ENDIF
        ENDIF
!
!---    Read initial fracture triangle solute concentrations  ---
!
        IF( ICT_FRCX.NE.0 ) THEN
          IF( ICT_FRCX.EQ.4 ) THEN
            VARB = 'Reference Solute Activity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRCX(1))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRCX(1)
            INDX = 0
            CALL RDUNIT(UNTS,VIC_FRCX(1),INDX)
            INDX = 0
            VARX = 1.D+0
            CALL RDUNIT(UNTS,VARX,INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRCX(1),', Bq)'
            VARB = 'X-Dir. Solute Activity Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRCX(2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRCX(2)
            INDX = 0
            IUNM = -1
            VIC_FRCX(2) = VARX*VIC_FRCX(2)
            CALL RDUNIT(UNTS,VIC_FRCX(2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRCX(2),
     &        ', Bq/m)'
            VARB = 'Y-Dir. Solute Activity Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRCX(3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRCX(3)
            INDX = 0
            IUNM = -1
            VIC_FRCX(3) = VARX*VIC_FRCX(3)
            CALL RDUNIT(UNTS,VIC_FRCX(3),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRCX(3),
     &        ', Bq/m)'
            VARB = 'Z-Dir. Solute Activity Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRCX(4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRCX(4)
            INDX = 0
            IUNM = -1
            VIC_FRCX(4) = VARX*VIC_FRCX(4)
            CALL RDUNIT(UNTS,VIC_FRCX(4),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRCX(4),
     &        ', Bq/m)'
          ELSE
            VARB = 'Reference Solute Concentration'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRCX(1))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRCX(1)
            INDX = 0
            IUNM = -3
            CALL RDUNIT(UNTS,VIC_FRCX(1),INDX)
            INDX = 0
            IUNM = -3
            VARX = 1.D+0
            CALL RDUNIT(UNTS,VARX,INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRCX(1),', 1/m^3)'
            VARB = 'X-Dir. Solute Concentration Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRCX(2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRCX(2)
            INDX = 0
            IUNM = -1
            VIC_FRCX(2) = VARX*VIC_FRCX(2)
            CALL RDUNIT(UNTS,VIC_FRCX(2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRCX(2),
     &        ', (1/m^3)/m)'
            VARB = 'Y-Dir. Solute Concentration Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRCX(3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRCX(3)
            INDX = 0
            IUNM = -1
            VIC_FRCX(3) = VARX*VIC_FRCX(3)
            CALL RDUNIT(UNTS,VIC_FRCX(3),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRCX(3),
     &        ', (1/m^3)/m)'
            VARB = 'Z-Dir. Solute Concentration Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRCX(4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRCX(4)
            INDX = 0
            IUNM = -1
            VIC_FRCX(4) = VARX*VIC_FRCX(4)
            CALL RDUNIT(UNTS,VIC_FRCX(4),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRCX(4),
     &        ', (1/m^3)/m)'
          ENDIF
!
!---      Loop over fracture triangles, initializing solute
!         concentrations in fracture triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
            C_FRC(NTX,NSL) = VIC_FRCX(1) +
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRCX(2) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRCX(3) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRCX(4)
          ENDDO
          ICT_FRC(NFX,NSL) = ICT_FRCX
        ENDIF
      ENDDO
!
!---  Loop over fractures, checking for uninitialized
!     hydrologic state conditions for non-restart simulations  ---
!
      IF( IEO.NE.2 .AND. ISLC(87).EQ.0 ) THEN
        DO NFX = 1,NF_FRC
          IF( IC_FRC(NFX).EQ.0 ) THEN
            INDX = 7
            IMSG = NFX
            CHMSG = 'Uninitialized Hydrologic State for Fracture'
            CALL WRMSGS( INDX )
          ENDIF
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIC_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDPROP_FRC_GT
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
!     Read input file for fracture properties information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 7 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE GEOM_FRC
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER IFRCX(LF_FRC)
      CHARACTER*64 UNTS
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDPROP_FRC_GT'
       DO K = 1,LF_FRC
         IFRCX(K) = 0
       ENDDO
!
!---  Write card information to ouput file  ---
!
      CARD = 'Fracture Properties Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Loop over the number of fractures  ---
!
      DO K = 1,NF_FRC
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Fracture Number: '
        CALL RDINT(ISTART,ICOMMA,CHDUM,NFX)
        WRITE(IWR,'(/,2X,A,I9)') 'Fracture Number: ',NFX
        IF( NFX.LT.1 .OR. NFX.GT.NF_FRC ) THEN
          INDX = 7
          IMSG = NFX
          CHMSG = 'Out of Range Fracture Number'
          CALL WRMSGS( INDX )
        ENDIF
        IFRCX(NFX) = 1
        VARB = 'Fracture Brooks and Corey (psi)'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,RKSP_FRC(1,NFX))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',RKSP_FRC(1,NFX)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,RKSP_FRC(1,NFX),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',RKSP_FRC(1,NFX),', m)'
        VARB = 'Fracture Brooks and Corey (lambda)'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,RKSP_FRC(3,NFX))
        WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),RKSP_FRC(3,NFX)
        VARB = 'Fracture Brooks and Corey (slr)'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,RKSP_FRC(4,NFX))
        WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),RKSP_FRC(4,NFX)
        VARB = 'Fracture Thermal Capacitance Thickness'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DTHCP_FRC(NFX))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',DTHCP_FRC(NFX)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DTHCP_FRC(NFX),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',DTHCP_FRC(NFX),', m)'
!
!---    Webb extension read optional oven dried head  ---
!
        CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
        IF( INDX.EQ.1 ) THEN
          CALL RDDPR(ISTART,ICOMMA,CHDUM,RKSP_FRC(5,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',RKSP_FRC(5,NFX)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,RKSP_FRC(5,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',RKSP_FRC(5,NFX),', m)'
        ELSE
          RKSP_FRC(5,NFX) = HDOD
        ENDIF
!
!---    Webb matching point parameters  ---
!
        SRX = RKSP_FRC(4,NFX)
        PSIX = RKSP_FRC(1,NFX)
        CLX = RKSP_FRC(3,NFX)
        CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
        RKSP_FRC(6,NFX) = SMPX
        RKSP_FRC(7,NFX) = HMPX
      ENDDO
!
!---  Check for fractures with unspecified parameters  ---
!
      DO K = 1,NF_FRC
        IF( IFRCX(K).EQ.0 ) THEN
          INDX = 7
          IMSG = NFX
          CHMSG = 'Unspecified Fracture Properties: Fracture Number'
          CALL WRMSGS( INDX )
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDPROP_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSR_FRC_GT
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
!     Read input file for fracture source information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 10 October 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PARM_FRC
      USE GEOM_FRC
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,BDUM,UNTS
      CHARACTER*512 CHDUM
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE ::  VAR
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSR_FRC_GT'
!
!---  Dynamic memory allocation  ---
!
      ALLOCATE( VAR(1:LSTM_FRC,1:7),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: VAR'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write card information to ouput file  ---
!
      CARD = 'Fracture Source Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      NSR_FRC = 0
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Fracture Sources'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
!
!---  Loop over the number of fracture sources  ---
!
      DO NS = 1,NLIN
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
!
!---  Read fracture source type  ---
!
        VARB = 'Fracture Source Type'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'aqueous volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Volumetric Source'
          IF( INDEX(ADUM(1:),'relative').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous Volumetric Source w/ ' //
     &        'Relative Saturation'
            ISRTX = 11
          ELSEIF( INDEX(ADUM(1:),'mass frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous Volumetric Source w/ ' //
     &        'Mass Fraction'
            ISRTX = 12
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Source Type: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Mass Source'
          IF( INDEX(ADUM(1:),'relative').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous Mass Source w/ ' //
     &        'Relative Saturation'
            ISRTX = 13
          ELSEIF( INDEX(ADUM(1:),'mass frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous Mass Source w/ ' //
     &        'Mass Fraction'
            ISRTX = 14
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Source Type: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous press').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Pressure Source'
          IF( INDEX(ADUM(1:),'relative').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous Pressure Source w/ ' //
     &        'Relative Saturation'
            ISRTX = 15
          ELSEIF( INDEX(ADUM(1:),'mass frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Aqueous Pressure Source w/ ' //
     &        'Mass Fraction'
            ISRTX = 16
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Source Type: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Gas Volumetric Source'
          IF( INDEX(ADUM(1:),'relative').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Volumetric Source w/ ' //
     &        'Relative Saturation'
            ISRTX = 21
          ELSEIF( INDEX(ADUM(1:),'mass frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Volumetric Source w/ ' //
     &        'Mass Fraction'
            ISRTX = 22
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Source Type: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Gas Mass Source'
          IF( INDEX(ADUM(1:),'relative').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Mass Source w/ ' //
     &        'Relative Saturation'
            ISRTX = 23
          ELSEIF( INDEX(ADUM(1:),'mass frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Mass Source w/ ' //
     &        'Mass Fraction'
            ISRTX = 24
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Source Type: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas press').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Gas Pressure Source'
          IF( INDEX(ADUM(1:),'relative').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Pressure Source w/ ' //
     &        'Relative Saturation'
            ISRTX = 25
          ELSEIF( INDEX(ADUM(1:),'mass frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Pressure Source w/ ' //
     &        'Mass Fraction'
            ISRTX = 26
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Source Type: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'sink').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Pressure Sink'
          ISRTX = 35
        ELSEIF( IEQC.NE.0 .AND. INDEX(ADUM(1:),'solute').NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          DO 30 NSL = 1,NSOLU
            IDB = INDEX(SOLUT(NSL)(1:),'  ')
            IF( INDEX(BDUM(1:),SOLUT(NSL)(1:IDB)).NE.0 ) THEN
              IF( INDEX(ADUM(1:),'density').NE.0 ) THEN
                ISRTX = -(NSL+NSOLU)
                WRITE(IWR,'(2X,2A)')'Solute Source Density: ',SOLUT(NSL)
              ELSE
                ISRTX = -NSL
                WRITE(IWR,'(2X,2A)')'Solute Source: ',SOLUT(NSL)
              ENDIF
              GOTO 40
            ENDIF
   30     CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Source Solute Name: '//BDUM
            CALL WRMSGS( INDX )
   40     CONTINUE
!#ifdef 1
!        ELSEIF( IEQC.NE.0 .AND. INDEX(ADUM(1:),'specie').NE.0 ) THEN
!          VARB = 'Species Name: '
!          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
!!
!!---      Aqueous species  ---
!!
!          DO 42 NSP = 1,NSPL
!            IDB = INDEX(SPNML(NSP)(1:),'  ')
!            IF( BDUM(1:IDB).EQ.SPNML(NSP)(1:IDB) ) THEN
!              IF( INDEX(ADUM(1:),'density').NE.0 ) THEN
!                ISRTX = 100+NSPL+NSPS+NSP
!                WRITE(IWR,'(2X,2A)')'Solute Source Density: ',
!     &            SPNML(NSP)(1:IDB)
!              ELSE
!                ISRTX = 100+NSP
!                WRITE(IWR,'(2X,2A)')'Species Source: ',
!     &            SPNML(NSP)(1:IDB)
!              ENDIF
!              GOTO 44
!            ENDIF
!   42     CONTINUE
!            INDX = 4
!            CHMSG = 'Unrecognized Source Species Name: '//BDUM
!            CALL WRMSGS( INDX )
!   44     CONTINUE
!#endif
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Source Type: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read source domain indices  ---
!
        VARB = 'Fracture Source Domain Indices'
        ISX = ISTART
        CALL RDINT(ISTART,ICOMMA,CHDUM,IF1X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,IT1X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,IF2X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,IT2X)
        ICX = ISTART
        WRITE(IWR,'(2X,A)') 'Fracture Source Domain:'
        WRITE(IWR,'(4X,2(A,I6),A)') 'Fracture = ',IF1X,
     &    ' Triangle = ',IT1X,' to '
        WRITE(IWR,'(4X,2(A,I6))') 'Fracture = ',IF2X,
     &    ' Triangle = ',IT2X
!
!---    Check for ill-defined fracture source domains  ---
!
        IF( IF1X.LT.1 .OR. IF1X.GT.NF_FRC .OR. IF2X.LT.1 .OR.
     &    IF2X.GT.NF_FRC .OR. IF2X.LT.IF1X ) THEN
          INDX = 4
          CHMSG = 'Invalid Fracture Source Domain: ' // CHDUM(ISX:ICX)
          CALL WRMSGS( INDX )
        ENDIF
        IF( IT1X.LT.1 .OR. IT1X.GT.NTP_FRC(IF1X) .OR. IT2X.LT.1 .OR.
     &    IT2X.GT.NTP_FRC(IF2X) .OR.
     &    ( IF1X.EQ.IF2X .AND. IT2X.LT.IT1X ) ) THEN
          INDX = 4
          CHMSG = 'Invalid Fracture Source Domain: ' // CHDUM(ISX:ICX)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Read number of fracture source times  ---
!
        VARB = 'Number of Fracture Source Times'
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISRM_FRC(NS))
        IF( ISRM_FRC(NS).GT.LSTM ) THEN
          INDX = 5
          CHMSG = 'Number of Fracture Source Times > Parameter LSTM'
          CALL WRMSGS( INDX )
        ENDIF
        SRTMO = -SMALL
        DO NTM = 1,ISRM_FRC(NS)
          DO M = 1,7
            VAR(NTM,M) = 0.D+0
          ENDDO
!
!---     Read and write fracture source values and units  ---
!
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Fracture Source Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',VAR(NTM,1)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR(NTM,1),INDX)
!
!---      Aqueous volumetric fracture source,
!         w/ relative saturation  ---
!
          IF( ISRTX.EQ.11 ) THEN
            VARB = 'Fracture Source Aqueous Volumetric Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
!
!---      Aqueous volumetric fracture source,
!         w/ mass fraction  ---
!
          ELSEIF( ISRTX.EQ.12 ) THEN
            VARB = 'Fracture Source Aqueous Volumetric Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
!
!---      Aqueous mass fracture source,
!         w/ relative saturation  ---
!
          ELSEIF( ISRTX.EQ.13 ) THEN
            VARB = 'Fracture Source Aqueous Mass Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
!
!---      Aqueous mass fracture source,
!         w/ mass fraction  ---
!
          ELSEIF( ISRTX.EQ.14 ) THEN
            VARB = 'Fracture Source Aqueous Mass Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
!
!---      Aqueous pressure source,
!         w/ relative saturation  ---
!
          ELSEIF( ISRTX.EQ.15 ) THEN
            VARB = 'Fracture Source Aqueous Pressure: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VAR(NTM,4) = VAR(NTM,4) - PATM
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
            VARB = 'Fracture Source Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
!
!---      Aqueous pressure source,
!         w/ mass fraction  ---
!
          ELSEIF( ISRTX.EQ.16 ) THEN
            VARB = 'Fracture Source Aqueous Pressure: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VAR(NTM,4) = VAR(NTM,4) - PATM
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
            VARB = 'Fracture Source Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
!
!---      Gas volumetric fracture source,
!         w/ relative saturation  ---
!
          ELSEIF( ISRTX.EQ.21 ) THEN
            VARB = 'Fracture Source Gas Volumetric Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
!
!---      Gas volumetric fracture source,
!         w/ mass fraction  ---
!
          ELSEIF( ISRTX.EQ.22 ) THEN
            VARB = 'Fracture Source Gas Volumetric Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Source Dissolved-Salt Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
!
!---      Gas mass fracture source,
!         w/ relative saturation  ---
!
          ELSEIF( ISRTX.EQ.23 ) THEN
            VARB = 'Fracture Source Gas Mass Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
!
!---      Gas mass fracture source,
!         w/ mass fraction  ---
!
          ELSEIF( ISRTX.EQ.24 ) THEN
            VARB = 'Fracture Source Gas Mass Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
!
!---      Gas pressure source,
!         w/ relative saturation  ---
!
          ELSEIF( ISRTX.EQ.25 ) THEN
            VARB = 'Fracture Source Gas Pressure: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VAR(NTM,4) = VAR(NTM,4) - PATM
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
            VARB = 'Fracture Source Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
!
!---      Gas pressure source,
!         w/ mass fraction  ---
!
          ELSEIF( ISRTX.EQ.26 ) THEN
            VARB = 'Fracture Source Gas Pressure: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VAR(NTM,4) = VAR(NTM,4) - PATM
            VARB = 'Fracture Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Fracture Source Dissolved-Air Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Fracture Source Dissolved-Salt Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
            VARB = 'Fracture Source Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
!
!---      Pressure sink  ---
!
          ELSEIF( ISRTX.EQ.35 ) THEN
            VARB = 'Fracture Sink Pressure: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VAR(NTM,4) = VAR(NTM,4) - PATM
            VARB = 'Borehole Skin Factor'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,2)
          ELSEIF( ISRTX.LT.0 .AND. ISRTX.GE.-NSOLU ) THEN
            VARB = 'Source Solute Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ENDIF
!
!---  Check for nonascending source times  ---
!
          IF( VAR(NTM,1).LT.SRTMO ) THEN
            INDX = 4
            CHMSG = 'Invalid Fracture Source Time Sequencing'
            CALL WRMSGS( INDX )
          ENDIF
          SRTMO = VAR(NTM,1)
        ENDDO
!
!---    Assign values to source variables  ---
!
        NSR_FRC = NSR_FRC + 1
        IF( NSR_FRC.GT.LSR_FRC ) THEN
          INDX = 5
          CHMSG = 'Number of Fracture Sources > Parameter LSR_FRC'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Starting triangle  ---
!
        NC1X = 0
        L1: DO KFX = 1,IF1X
          DO KTX = 1,NTP_FRC(KFX)
            NC1X = NC1X + 1
            IF( KFX.EQ.IF1X .AND. KTX.EQ.IT1X ) EXIT L1
          ENDDO
        ENDDO L1
        ISRDM_FRC(1,NSR_FRC) = NC1X
 !
!---    Ending triangle  ---
!
        NC2X = 0
        L2: DO KFX = 1,IF2X
          DO KTX = 1,NTP_FRC(KFX)
            NC2X = NC2X + 1
            IF( KFX.EQ.IF2X .AND. KTX.EQ.IT2X ) EXIT L2
          ENDDO
        ENDDO L2
        ISRDM_FRC(2,NSR_FRC) = NC2X
        ISRT_FRC(NSR_FRC) = ISRTX
        DO NTM = 1,ISRM_FRC(NS)
          DO M = 1,7
            SRC_FRC(M,NTM,NSR_FRC) = VAR(NTM,M)
          ENDDO
        ENDDO
      ENDDO
!
!---  Deallocate memory  ---
!
      IF( ALLOCATED(VAR) ) THEN
      DEALLOCATE( VAR,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: VAR'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSR_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTHMOC_GT
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
!     Read input file for fluid thermo-catalytic properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 4 May 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,UNTS
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDTHMOC_GT'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Thermo-catalytic Borehole Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read thermo-catalytic fluid equation form  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Thermo-catalytic Fluid Equation'
      CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
!
!---  First-order rate equation, with catalyst and salt tracking
!
!     hl = hl_aqu + rhol*min( hl_rate, k_rate*dt/
!       (1 + exp((T_trans-T)/rate)) )  ---
!
      IF( INDEX(ADUM(1:),'first').NE.0 .AND.
     &  INDEX(ADUM(1:),'order').NE.0 .AND.
     &  INDEX(ADUM(1:),'catalyst').NE.0 ) THEN
        WRITE(IWR,'(/,2X,A)') 'First-Order Rate Equation, with '
     &      // 'Catalyst with Salt Tracking'
        ITC_BH = 24
!
!---    Read enthalpy of reaction density  ---
!
        VARB = 'Enthalpy of Reaction Density'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PTC_BH(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &   UNTS(1:NCH),': ',PTC_BH(1)
        INDX = 0
        IUNKG = 1
        IUNM = -1
        IUNS = -2
        CALL RDUNIT(UNTS,PTC_BH(1),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PTC_BH(1),', J/m^3 aqu)'
!
!---    Read reaction rate (i.e., k) ---
!
        VARB = 'Reaction Rate Constant'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PTC_BH(2))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &   UNTS(1:NCH),': ',PTC_BH(2)
        INDX = 0
        IUNS = -1
        CALL RDUNIT(UNTS,PTC_BH(2),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PTC_BH(2),', 1/s)'
!
!---    Read transition temperature ---
!
        VARB = 'Mid-Point Transition Temperature'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PTC_BH(3))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &   UNTS(1:NCH),': ',PTC_BH(3)
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,PTC_BH(3),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PTC_BH(3),', C)'
!
!---    Read transition rise rate ---
!
        VARB = 'Transition Rise Rate'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PTC_BH(4))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &   UNTS(1:NCH),': ',PTC_BH(4)
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,PTC_BH(4),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PTC_BH(4),', C)'
!
!---    Read salt mass fraction basis ---
!
        VARB = 'Aqueous Salt Mass Fraction Basis'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PTC_BH(5))
        WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),': ',PTC_BH(5)
!
!---  First-order rate equation, without catalyst and salt tracking
!
!     hl = hl_aqu + rhol*min( hl_rate, k_rate*dt/
!       (1 + exp((T_trans-T)/rate)) )  ---
!
      ELSEIF( INDEX(ADUM(1:),'first').NE.0 .AND.
     &  INDEX(ADUM(1:),'order').NE.0 ) THEN
        WRITE(IWR,'(/,2X,A)') 'First-Order Rate Equation, without '
     &    // 'Catalyst with Salt Tracking'
        ITC_BH = 22
!
!---    Read enthalpy of reaction density  ---
!
        VARB = 'Enthalpy of Reaction Density'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PTC_BH(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &   UNTS(1:NCH),': ',PTC_BH(1)
        INDX = 0
        IUNKG = 1
        IUNM = -1
        IUNS = -2
        CALL RDUNIT(UNTS,PTC_BH(1),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PTC_BH(1),', J/m^3 aqu)'
!
!---    Read reaction rate (i.e., k) ---
!
        VARB = 'Reaction Rate Constant'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PTC_BH(2))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &   UNTS(1:NCH),': ',PTC_BH(2)
        INDX = 0
        IUNS = -1
        CALL RDUNIT(UNTS,PTC_BH(2),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PTC_BH(2),', 1/s)'
!
!---    Read transition temperature ---
!
        VARB = 'Mid-Point Transition Temperature'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PTC_BH(3))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &   UNTS(1:NCH),': ',PTC_BH(3)
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,PTC_BH(3),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PTC_BH(3),', C)'
!
!---    Read transition rise rate ---
!
        VARB = 'Transition Rise Rate'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PTC_BH(4))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &   UNTS(1:NCH),': ',PTC_BH(4)
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,PTC_BH(4),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PTC_BH(4),', C)'
!
!---    Read salt mass fraction basis ---
!
        VARB = 'Aqueous Salt Mass Fraction Basis'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PTC_BH(5))
        WRITE(IWR,'(4X,2A,1PE11.4)') VARB(1:IVR),': ',PTC_BH(5)
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Thermo-catalytic Fluid Equation Form: '
     &    // ADUM(1:NCH)
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTHMOC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKSP_BH_GT( PGX,PLX,RKGX,RKLX,SGX,SLX,NBN )
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
!     Aqueous saturation, gas saturation, aqueous relative permeability,
!     and gas relative permeability for boreholes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, May 2, 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLZ(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RKSP_BH_GT'
!
!---  Pipe flow  ---
!
      IF( IS_BH(INV_BH(NBN)).GT.10 ) THEN
        INVX = INV_BH(NBN)
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        HMPX = PAR_BH(12,INVX)
        SLRX = PAR_BH(9,INVX)
        CL = MAX( PAR_BH(8,INVX),SMALL )
!
!---    Capillary head above the matching point head,
!       use Webb extension  ---
!
        IF( HDGL.GT.HMPX ) THEN
          SMPX = PAR_BH(11,INVX)
          HDGL = MIN( HDGL,PAR_BH(10,INVX) )
          DMPX = SMPX/(LOG10(PAR_BH(10,INVX))-LOG10(HMPX))
          SLX = -(LOG10(HDGL)-LOG10(PAR_BH(10,INVX)))*DMPX
          ASLX = SLX
!
!---    Capillary head at or below the matching point head,
!       use Brooks-Corey function
!
        ELSE
          IF( HDGL.LE.PAR_BH(7,INVX) ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (PAR_BH(7,INVX)/HDGL)**CL
          ENDIF
          SLX = ASLX*(1.D+0-SLRX) + SLRX
          ASLX = SLX
        ENDIF
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
!
!---    Aqueous relative permeability, Brooks and Corey function  ---
!
        RKLX = ESLX**((2.D+0 + 3.D+0*CL)/CL)
!
!---    Gas relative permeability, Brooks and Corey function  ---
!
        RKGX = ((1.D+0-ESLX)**2)*(1.D+0-(ESLX**((2.D+0+CL)/CL)))
!
!---  Borehole filled with porous media flow  ---
!
      ELSE
        IZN = IZ_BH(NBN)
        INDX = 1
        ASL_F = 0.D+0
        ASL_M = 0.D+0
        ASLX = 0.D+0
        ASLMINX = 0.D+0
        ESGTX = 0.D+0
        SGTX = 0.D+0
        BTGLX = 1.D+0
!
!---    Aqueous and gas saturation  ---
!
        IF( NPHAZ_BH(2,NBN).EQ.4 ) THEN
          ESLX = 0.D+0
          SLX = 0.D+0
          SGX = 1.D+0
        ELSE
          CALL SP_GT( ASL_F,ASL_M,ASLX,ASLMINX,ESLX,ESGTX,ESGTMX,
     &      PGX,PLX,BTGLX,SGX,SGTX,SLX,SLFX,SLMX,SLRX,INDX,IZN )
        ENDIF
!
!---    Aqueous relative permeability  ---
!
        CALL RKLS_GT( ASL_F,ASL_M,ESLX,PGX,PLX,RKLZ,SLX,IZN )
        RKLX = RKLZ(1)
!
!---    Gas relative permeability  ---
!
        CALL RKGS_GT( ASL_F,ASL_M,ASLX,ASLMINX,ESGTX,PGX,PLX,RKGX,
     &    SGX,SGTX,SLX,IZN )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKSP_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKSP_FRC_GT( PGX,PLX,RKGX,RKLX,SGX,SLX,NFX )
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
!     Relative permeability and saturations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 7 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
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
      SUB_LOG(ISUB_LOG) = '/RKSP_FRC_GT'
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        HMPX = RKSP_FRC(7,NFX)
        SLRX = RKSP_FRC(4,NFX)
        CL = MAX( RKSP_FRC(3,NFX),SMALL )
!
!---    Capillary head above the matching point head,
!       use Webb extension  ---
!
        IF( HDGL.GT.HMPX ) THEN
          SMPX = RKSP_FRC(6,NFX)
          HDGL = MIN( HDGL,RKSP_FRC(5,NFX) )
          DMPX = SMPX/(LOG10(RKSP_FRC(5,NFX))-LOG10(HMPX))
          SLX = -(LOG10(HDGL)-LOG10(RKSP_FRC(5,NFX)))*DMPX
          ASLX = SLX
!
!---    Capillary head at or below the matching point head,
!       use Brooks-Corey function
!
        ELSE
          IF( HDGL.LE.RKSP_FRC(1,NFX) ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (RKSP_FRC(1,NFX)/HDGL)**CL
          ENDIF
          SLRX = RKSP_FRC(4,NFX)
          SLX = ASLX*(1.D+0-SLRX) + SLRX
          ASLX = SLX
        ENDIF
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
!
!---    Aqueous relative permeability, Brooks and Corey function  ---
!
        RKLX = ESLX**((2.D+0 + 3.D+0*CL)/CL)
!
!---    Gas relative permeability, Brooks and Corey function  ---
!
        RKGX = ((1.D+0-ESLX)**2)*(1.D+0-(ESLX**((2.D+0+CL)/CL)))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKSP_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_BF_GT
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
!     Compute the maximum fracture and borehole relative residuals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 2 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE PARM_BH
      USE OUTPU
      USE JACOB
      USE HYST
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FILES
      USE FDVT_FRC
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP_BH
      USE FDVP
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 PH_CND(5)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Saturated w/o Entrapped Gas',
     &   'Unsaturated w/ or w/o Entrapped Gas',
     &   'Saturated w/ Trapped Gas',
     &   'Fully Unsaturated',
     &   'Supercritical Water'/
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL_BF_GT'
!
!---  Zero maximum residuals  ---
!
      DO M = 1,ISVC
        RSD_FRC(M) = 0.D+0
        NSD_FRC(M) = 0
        RSD_BH(M) = 0.D+0
        NSD_BH(M) = 0
      ENDDO
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
          N_DB = NTX
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) MPT = IM_FRC(IEQT,NTX)
          MPL = IM_FRC(IEQW,NTX)
          IF( ISLC(37).EQ.0 ) MPG = IM_FRC(IEQA,NTX)
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) MPS = IM_FRC(IEQS,NTX)
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Energy equation  ---
!
            ACP = (SL_FRC(2,NTX)*RHOL_FRC(2,NTX)*HL_FRC(2,NTX) +
     &        SS_FRC(2,NTX)*RHOSP_FRC(2,NTX)*HSP_FRC(2,NTX) +
     &        SG_FRC(2,NTX)*RHOG_FRC(2,NTX)*UEG_FRC(2,NTX))*DTI*
     &        AF_FRC(NTX)*APM_FRC(2,NTX)
            RSDX = MIN( ABS(BLU(MPT))/TABS,
     &        ABS(RSDL_FRC(IEQT,NTX)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD_FRC(IEQT) ) THEN
              RSD_FRC(IEQT) = RSDX
              NSD_FRC(IEQT) = NTX
            ENDIF
          ENDIF
!
!---      Saturated system w/o entrapped gas
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - aqueous-air mole fraction
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ_FRC(2,NTX).EQ.1 ) THEN
!
!---        Water mass equation  ---
!
            ACP = (RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &        XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &        XGW_FRC(2,NTX))*DTI*AF_FRC(NTX)*APM_FRC(2,NTX)
            RSDX = MIN( ABS(BLU(MPL))/(ABS(PL_FRC(2,NTX))+PATM),
     &        ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
            RSDX = 1.D-1*RSDX
            IF( RSDX.GT.RSD_FRC(IEQW) ) THEN
              RSD_FRC(IEQW) = RSDX
              NSD_FRC(IEQW) = NTX
            ENDIF
!
!---      Air mass equation, ignore residual for small aqueous-air  ---
!
            INDX = 0
            CALL REGION_4( T_FRC(2,NTX),PSWX,INDX )
            PWX = MAX( PSWX,PL_FRC(2,NTX)+PATM )
            CALL SOL_BRNS( T_FRC(2,NTX),PWX,XLSMX )
            XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
            XLS_FRC(2,NTX) = XLSX
            CALL SP_B( T_FRC(2,NTX),XLS_FRC(2,NTX),PSBX )
            IF( ISLC(37).EQ.0 ) THEN
              PVAX = MAX( PWX-PSBX,0.D+0 )
              XMLAX = PVAX/HCAW
              XMLA_FRC(2,NTX) = PVA_FRC(2,NTX)/HCAW
              IF( XMLA_FRC(2,NTX).GT.(1.D-6*XMLAX) ) THEN
                ACP = (RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &            XGA_FRC(2,NTX) + RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &            XLA_FRC(2,NTX))*DTI*AF_FRC(NTX)*APM_FRC(2,NTX)
                RSDX = MIN( ABS(BLU(MPG))/MAX( XMLAX,PATM/HCAW ),
     &            ABS(RSDL_FRC(IEQA,NTX)/(ACP+SMALL)) )
                IF( RSDX.GT.RSD_FRC(IEQA) ) THEN
                  RSD_FRC(IEQA) = RSDX
                  NSD_FRC(IEQA) = NTX
                ENDIF
              ENDIF
            ENDIF
!
!---        Salt mass equation, isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              ACP = TMS_FRC(2,NTX)*DTI*AF_FRC(NTX)*APM_FRC(2,NTX)
              RSDX = MIN( (ABS(BLU(MPS))/XLSMX),
     &          ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
              RSDX = RSDX*1.D-1
              IF( RSDX.GT.RSD_FRC(IEQS) ) THEN
                RSD_FRC(IEQS) = RSDX
                NSD_FRC(IEQS) = NTX
              ENDIF
            ENDIF
!
!---      Unsaturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ_FRC(2,NTX).EQ.2 ) THEN
!
!---        Water mass equation  ---
!
            ACP = (RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &        XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &        XGW_FRC(2,NTX))*DTI*AF_FRC(NTX)*APM_FRC(2,NTX)
            RSDX = MIN( ABS(BLU(MPL))/(ABS(PL_FRC(2,NTX))+PATM),
     &        ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD_FRC(IEQW) ) THEN
              RSD_FRC(IEQW) = RSDX
              NSD_FRC(IEQW) = NTX
            ENDIF
!
!---        Air mass equation  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              IF( SG_FRC(2,NTX).GT.1.D-3 ) THEN
                ACP = (RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &            XGA_FRC(2,NTX) + RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &            XLA_FRC(2,NTX))*DTI*AF_FRC(NTX)*APM_FRC(2,NTX)
                RSDX = MIN( ABS(BLU(MPG))/(ABS(PG_FRC(2,NTX))+PATM),
     &            ABS(RSDL_FRC(IEQA,NTX)/(ACP+SMALL)) )
                IF( RSDX.GT.RSD_FRC(IEQA) ) THEN
                  RSD_FRC(IEQA) = RSDX
                  NSD_FRC(IEQA) = NTX
                ENDIF
              ENDIF
            ENDIF
!
!---        Salt mass equation, isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              ACP = TMS_FRC(2,NTX)*DTI*AF_FRC(NTX)*APM_FRC(2,NTX)
              INDX = 0
              CALL REGION_4( T_FRC(2,NTX),PWX,INDX )
              CALL SOL_BRNS( T_FRC(2,NTX),PWX,XLSMX )
              RSDX = MIN( (ABS(BLU(MPS))/XLSMX),
     &          ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
              RSDX = RSDX*1.D-1
              IF( RSDX.GT.RSD_FRC(IEQS) ) THEN
                RSD_FRC(IEQS) = RSDX
                NSD_FRC(IEQS) = NTX
              ENDIF
            ENDIF
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - gas pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ_FRC(2,NTX).EQ.4 ) THEN
!
!---        Water mass equation  ---
!
            INDX = 0
            CALL REGION_4( T_FRC(2,NTX),PWX,INDX )
            PWX = MIN( PWX,PCRW )
            ACP = RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &        XGW_FRC(2,NTX)*DTI*AF_FRC(NTX)*APM_FRC(2,NTX)
            RSDX = MIN( ABS(BLU(MPL))/PWX,
     &        ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD_FRC(IEQW) ) THEN
              RSD_FRC(IEQW) = RSDX
              NSD_FRC(IEQW) = NTX
            ENDIF
!
!---        Air mass equation  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              PGX = PG_FRC(2,NTX) + PATM
              ACP = RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &          XGA_FRC(2,NTX)*DTI*AF_FRC(NTX)*APM_FRC(2,NTX)
              RSDX = MIN( ABS(BLU(MPG))/ABS(PGX),
     &          ABS(RSDL_FRC(IEQA,NTX)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD_FRC(IEQA) ) THEN
                RSD_FRC(IEQA) = RSDX
                NSD_FRC(IEQA) = NTX
              ENDIF
            ENDIF
!
!---        Salt mass equation, isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              CALL DENS_S( T_FRC(2,NTX),PGX,RHOSPX )
              ACP = TMS_FRC(2,NTX)*DTI*AF_FRC(NTX)*APM_FRC(2,NTX)
              RSDX = MIN( (ABS(BLU(MPS))/RHOSPX),
     &          ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD_FRC(IEQS) ) THEN
                RSD_FRC(IEQS) = RSDX
                NSD_FRC(IEQS) = NTX
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
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
          N_DB = NBN
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) MPT = IM_BH(IEQT,NBN)
          MPL = IM_BH(IEQW,NBN)
          IF( ISLC(37).EQ.0 ) MPG = IM_BH(IEQA,NBN)
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) MPS = IM_BH(IEQS,NBN)
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Energy equation  ---
!
            ACP = (SL_BH(2,NBN)*RHOL_BH(2,NBN)*HL_BH(2,NBN) +
     &        SS_BH(2,NBN)*RHOSP_BH(2,NBN)*HSP_BH(2,NBN) +
     &        SG_BH(2,NBN)*RHOG_BH(2,NBN)*UEG_BH(2,NBN))*DTI*
     &        VOL_BH(NBN)
            RSDX = MIN( ABS(BLU(MPT))/TABS,
     &        ABS(RSDL_BH(IEQT,NBN)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD_BH(IEQT) ) THEN
              RSD_BH(IEQT) = RSDX
              NSD_BH(IEQT) = NBN
            ENDIF
          ENDIF
!
!---      Saturated system w/o entrapped gas
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - aqueous-air mole fraction
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ_BH(2,NBN).EQ.1 ) THEN
!
!---        Water mass equation  ---
!
            ACP = (RHOL_BH(2,NBN)*SL_BH(2,NBN)*
     &        XLW_BH(2,NBN) + RHOG_BH(2,NBN)*SG_BH(2,NBN)*
     &        XGW_BH(2,NBN))*DTI*VOL_BH(NBN)
            RSDX = MIN( ABS(BLU(MPL))/(ABS(PL_BH(2,NBN))+PATM),
     &        ABS(RSDL_BH(IEQW,NBN)/(ACP+SMALL)) )
            RSDX = 1.D-1*RSDX
            IF( RSDX.GT.RSD_BH(IEQW) ) THEN
              RSD_BH(IEQW) = RSDX
              NSD_BH(IEQW) = NBN
            ENDIF
!
!---        Air mass equation, ignore residual for small aqueous-air  ---
!
            INDX = 0
            CALL REGION_4( T_BH(2,NBN),PSWX,INDX )
            PWX = MAX( PSWX,PL_BH(2,NBN)+PATM )
            CALL SOL_BRNS( T_BH(2,NBN),PWX,XLSMX )
            XLSX = MIN( YLS_BH(2,NBN),XLSMX )
            XLS_BH(2,NBN) = XLSX
            CALL SP_B( T_BH(2,NBN),XLS_BH(2,NBN),PSBX )
            IF( ISLC(37).EQ.0 ) THEN
              PVAX = MAX( PWX-PSBX,0.D+0 )
              XMLAX = PVAX/HCAW
              XMLA_BH(2,NBN) = PVA_BH(2,NBN)/HCAW
              IF( XMLA_BH(2,NBN).GT.(1.D-6*XMLAX) ) THEN
                ACP = (RHOG_BH(2,NBN)*SG_BH(2,NBN)*
     &            XGA_BH(2,NBN) + RHOL_BH(2,NBN)*SL_BH(2,NBN)*
     &            XLA_BH(2,NBN))*DTI*VOL_BH(NBN)
                RSDX = MIN( ABS(BLU(MPG))/MAX( XMLAX,PATM/HCAW ),
     &            ABS(RSDL_BH(IEQA,NBN)/(ACP+SMALL)) )
                IF( RSDX.GT.RSD_BH(IEQA) ) THEN
                  RSD_BH(IEQA) = RSDX
                  NSD_BH(IEQA) = NBN
                ENDIF
              ENDIF
            ENDIF
!
!---        Salt mass equation, isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              ACP = TMS_BH(2,NBN)*DTI*VOL_BH(NBN)
              RSDX = MIN( (ABS(BLU(MPS))/XLSMX),
     &          ABS(RSDL_BH(IEQS,NBN)/(ACP+SMALL)) )
              RSDX = RSDX*1.D-1
              IF( RSDX.GT.RSD_BH(IEQS) ) THEN
                RSD_BH(IEQS) = RSDX
                NSD_BH(IEQS) = NBN
              ENDIF
            ENDIF
!
!---      Unsaturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ_BH(2,NBN).EQ.2 ) THEN
!
!---        Water mass equation  ---
!
            ACP = (RHOL_BH(2,NBN)*SL_BH(2,NBN)*
     &        XLW_BH(2,NBN) + RHOG_BH(2,NBN)*SG_BH(2,NBN)*
     &        XGW_BH(2,NBN))*DTI*VOL_BH(NBN)
            RSDX = MIN( ABS(BLU(MPL))/(ABS(PL_BH(2,NBN))+PATM),
     &        ABS(RSDL_BH(IEQW,NBN)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD_BH(IEQW) ) THEN
              RSD_BH(IEQW) = RSDX
              NSD_BH(IEQW) = NBN
            ENDIF
!
!---        Air mass equation  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              IF( SG_BH(2,NBN).GT.1.D-3 ) THEN
                ACP = (RHOG_BH(2,NBN)*SG_BH(2,NBN)*
     &            XGA_BH(2,NBN) + RHOL_BH(2,NBN)*SL_BH(2,NBN)*
     &            XLA_BH(2,NBN))*DTI*VOL_BH(NBN)
                RSDX = MIN( ABS(BLU(MPG))/(ABS(PG_BH(2,NBN))+PATM),
     &            ABS(RSDL_BH(IEQA,NBN)/(ACP+SMALL)) )
                IF( RSDX.GT.RSD_BH(IEQA) ) THEN
                  RSD_BH(IEQA) = RSDX
                  NSD_BH(IEQA) = NBN
                ENDIF
              ENDIF
            ENDIF
!
!---        Salt mass equation, isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              ACP = TMS_BH(2,NBN)*DTI*VOL_BH(NBN)
              INDX = 0
              CALL REGION_4( T_BH(2,NBN),PWX,INDX )
              CALL SOL_BRNS( T_BH(2,NBN),PWX,XLSMX )
              RSDX = MIN( (ABS(BLU(MPS))/XLSMX),
     &          ABS(RSDL_BH(IEQS,NBN)/(ACP+SMALL)) )
              RSDX = RSDX*1.D-1
              IF( RSDX.GT.RSD_BH(IEQS) ) THEN
                RSD_BH(IEQS) = RSDX
                NSD_BH(IEQS) = NBN
              ENDIF
            ENDIF
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - gas pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ_BH(2,NBN).EQ.4 ) THEN
!
!---        Water mass equation  ---
!
            INDX = 0
            CALL REGION_4( T_BH(2,NBN),PWX,INDX )
            PWX = MIN( PWX,PCRW )
            ACP = RHOG_BH(2,NBN)*SG_BH(2,NBN)*
     &        XGW_BH(2,NBN)*DTI*VOL_BH(NBN)
            RSDX = MIN( ABS(BLU(MPL))/PWX,
     &        ABS(RSDL_BH(IEQW,NBN)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD_BH(IEQW) ) THEN
              RSD_BH(IEQW) = RSDX
              NSD_BH(IEQW) = NBN
            ENDIF
!
!---        Air mass equation  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              PGX = PG_BH(2,NBN) + PATM
              ACP = RHOG_BH(2,NBN)*SG_BH(2,NBN)*
     &          XGA_BH(2,NBN)*DTI*VOL_BH(NBN)
              RSDX = MIN( ABS(BLU(MPG))/ABS(PGX),
     &          ABS(RSDL_BH(IEQA,NBN)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD_BH(IEQA) ) THEN
                RSD_BH(IEQA) = RSDX
                NSD_BH(IEQA) = NBN
              ENDIF
            ENDIF
!
!---        Salt mass equation, isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              CALL DENS_S( T_BH(2,NBN),PGX,RHOSPX )
              ACP = TMS_BH(2,NBN)*DTI*VOL_BH(NBN)
              RSDX = MIN( (ABS(BLU(MPS))/RHOSPX),
     &          ABS(RSDL_BH(IEQS,NBN)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD_BH(IEQS) ) THEN
                RSD_BH(IEQS) = RSDX
                NSD_BH(IEQS) = NBN
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
!---  Assign a convergence index  ---
!
      RSDX = 1.D-20
      DO M = 1,ISVC
        IF( RSD_FRC(M).GT.RSDMX ) ICNV = 2
        IF( RSD_BH(M).GT.RSDMX ) ICNV = 2
        RSDX = MAX( RSD_FRC(M),RSD_BH(M),RSDX )
      ENDDO
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) ICNV = 1
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.1 ) THEN
        IF( RSDX.GE.1.D+2 ) THEN
          WRITE(ISC,'(10X,A)') '---  Excessive Residual in Fracture' //
     &      ' or Borehole  ---'
          WRITE(IWR,'(10X,A)') '---  Excessive Residual in Fracture' //
     &      ' or Borehole  ---'
        ELSE
          WRITE(ISC,'(10X,A)') '---  Convergence Failure in Fracture' //
     &      ' or Borehole  ---'
          WRITE(IWR,'(10X,A)') '---  Convergence Failure in Fracture' //
     &      ' or Borehole  ---'
        ENDIF
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          NT = NSD_FRC(IEQT)
          IF( NT.GT.0 ) THEN
            NPT = NPHAZ_FRC(2,NT)
            NCHT = INDEX( PH_CND(NPT),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Energy Equation Maximum Residual = ',RSD_FRC(IEQT),
     &        ': Fracture Triangle = ',NT,
     &        ': Phase Condition = ',PH_CND(NPT)(1:NCHT)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Energy Equation Maximum Residual = ',RSD_FRC(IEQT),
     &        ': Fracture Triangle = ',NT,
     &        ': Phase Condition = ',PH_CND(NPT)(1:NCHT)
          ENDIF
          NT = NSD_BH(IEQT)
          IF( NT.GT.0 ) THEN
            NPT = NPHAZ_BH(2,NT)
            NCHT = INDEX( PH_CND(NPT),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Energy Equation Maximum Residual = ',RSD_BH(IEQT),
     &        ': Borehole Node = ',NT,
     &        ': Phase Condition = ',PH_CND(NPT)(1:NCHT)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Energy Equation Maximum Residual = ',RSD_BH(IEQT),
     &        ': Borehole Node = ',NT,
     &        ': Phase Condition = ',PH_CND(NPT)(1:NCHT)
          ENDIF
        ENDIF
        NW = NSD_FRC(IEQW)
        IF( NW.GT.0 ) THEN
          NPW = NPHAZ_FRC(2,NW)
          NCHW = INDEX( PH_CND(NPW),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Equation Maximum Residual = ',RSD_FRC(IEQW),
     &      ': Fracture Triangle = ',NW,
     &      ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Equation Maximum Residual = ',RSD_FRC(IEQW),
     &      ': Fracture Triangle = ',NW,
     &      ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
!
!---      Extended output on convergence failure  ---
!
          IF( ISLC(62).EQ.1 ) THEN
            WRITE(ISC,'(4X,2A)')
     &        'Current Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Current Aqueous Saturation = ',SL_FRC(2,NW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Current Salt Saturation = ',SS_FRC(2,NW)
            WRITE(IWR,'(4X,2A)')
     &        'Current Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Current Aqueous Saturation = ',SL_FRC(2,NW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Current Salt Saturation = ',SS_FRC(2,NW)
            NPW = NPHAZ_FRC(1,NW)
            NCHW = INDEX( PH_CND(NPW),'  ') - 1
            WRITE(ISC,'(4X,2A)')
     &        'Previous Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Previous Aqueous Saturation = ',SL_FRC(1,NW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Previous Salt Saturation = ',SS_FRC(1,NW)
            WRITE(IWR,'(4X,2A)')
     &        'Previous Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Previous Aqueous Saturation = ',SL_FRC(1,NW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Previous Salt Saturation = ',SS_FRC(1,NW)
          ENDIF
        ENDIF
        NW = NSD_BH(IEQW)
        IF( NW.GT.0 ) THEN
          NPW = NPHAZ_BH(2,NW)
          NCHW = INDEX( PH_CND(NPW),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Equation Maximum Residual = ',RSD_BH(IEQW),
     &      ': Borehole Node = ',NW,
     &      ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Equation Maximum Residual = ',RSD_BH(IEQW),
     &      ': Borehole Node = ',NW,
     &      ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
!
!---      Extended output on convergence failure  ---
!
          IF( ISLC(62).EQ.1 ) THEN
            WRITE(ISC,'(4X,2A)')
     &        'Current Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Current Aqueous Saturation = ',SL_BH(2,NW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Current Salt Saturation = ',SS_BH(2,NW)
            WRITE(IWR,'(4X,2A)')
     &        'Current Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Current Aqueous Saturation = ',SL_BH(2,NW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Current Salt Saturation = ',SS_BH(2,NW)
            NPW = NPHAZ_BH(1,NW)
            NCHW = INDEX( PH_CND(NPW),'  ') - 1
            WRITE(ISC,'(4X,2A)')
     &        'Previous Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Previous Aqueous Saturation = ',SL_BH(1,NW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Previous Salt Saturation = ',SS_BH(1,NW)
            WRITE(IWR,'(4X,2A)')
     &        'Previous Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Previous Aqueous Saturation = ',SL_BH(1,NW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Previous Salt Saturation = ',SS_BH(1,NW)
          ENDIF
        ENDIF
        IF( ISLC(37).EQ.0 ) THEN
          NA = NSD_FRC(IEQA)
          IF( NA.GT.0 ) THEN
            NPA = NPHAZ_FRC(2,NA)
            NCHA = INDEX( PH_CND(NPA),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Air Equation Maximum Residual = ',RSD_FRC(IEQA),
     &        ': Fracture Triangle = ',NA,
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Air Equation Maximum Residual = ',RSD_FRC(IEQA),
     &        ': Fracture Triangle = ',NA,
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
!
!---        Extended output on convergence failure  ---
!
            IF( ISLC(62).EQ.1 ) THEN
              WRITE(ISC,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL_FRC(2,NA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS_FRC(2,NA)
              WRITE(IWR,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL_FRC(2,NA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS_FRC(2,NA)
              NPA = NPHAZ_FRC(1,NA)
              NCHA = INDEX( PH_CND(NPA),'  ') - 1
              WRITE(ISC,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL_FRC(1,NA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS_FRC(1,NA)
              WRITE(IWR,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL_FRC(1,NA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS_FRC(1,NA)
            ENDIF
          ENDIF
          NA = NSD_BH(IEQA)
          IF( NA.GT.0 ) THEN
            NPA = NPHAZ_BH(2,NA)
            NCHA = INDEX( PH_CND(NPA),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Air Equation Maximum Residual = ',RSD_BH(IEQA),
     &        ': Borehole Node = ',NA,
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Air Equation Maximum Residual = ',RSD_BH(IEQA),
     &        ': Borehole Node = ',NA,
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
!
!---        Extended output on convergence failure  ---
!
            IF( ISLC(62).EQ.1 ) THEN
              WRITE(ISC,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL_BH(2,NA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS_BH(2,NA)
              WRITE(IWR,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL_BH(2,NA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS_BH(2,NA)
              NPA = NPHAZ_BH(1,NA)
              NCHA = INDEX( PH_CND(NPA),'  ') - 1
              WRITE(ISC,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL_BH(1,NA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS_BH(1,NA)
              WRITE(IWR,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL_BH(1,NA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS_BH(1,NA)
            ENDIF
          ENDIF
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          NS = NSD_FRC(IEQS)
          IF( NS.GT.0 ) THEN
            NPS = NPHAZ_FRC(2,NS)
            NCHS = INDEX( PH_CND(NPS),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD_FRC(IEQS),
     &        ': Fracture Triangle = ',NS,
     &        ': Phase Condition = ',PH_CND(NPS)(1:NCHS)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD_FRC(IEQS),
     &        ': Fracture Triangle = ',NS,
     &        ': Phase Condition = ',PH_CND(NPS)(1:NCHS)
!
!---        Extended output on convergence failure  ---
!
            IF( ISLC(62).EQ.1 ) THEN
              WRITE(ISC,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL_FRC(2,NS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS_FRC(2,NS)
              WRITE(IWR,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL_FRC(2,NS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS_FRC(2,NS)
              NPS = NPHAZ_FRC(1,NS)
              NCHS = INDEX( PH_CND(NPS),'  ') - 1
              WRITE(ISC,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL_FRC(1,NS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS_FRC(1,NS)
              WRITE(IWR,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL_FRC(1,NS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS_FRC(1,NS)
            ENDIF
          ENDIF
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          NS = NSD_BH(IEQS)
          IF( NS.GT.0 ) THEN
            NPS = NPHAZ_BH(2,NS)
            NCHS = INDEX( PH_CND(NPS),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD_BH(IEQS),
     &        ': Borehole Node = ',NS,
     &        ': Phase Condition = ',PH_CND(NPS)(1:NCHS)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD_BH(IEQS),
     &        ': Borehole Node = ',NS,
     &        ': Phase Condition = ',PH_CND(NPS)(1:NCHS)
!
!---        Extended output on convergence failure  ---
!
            IF( ISLC(62).EQ.1 ) THEN
              WRITE(ISC,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL_BH(2,NS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS_BH(2,NS)
              WRITE(IWR,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL_BH(2,NS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS_BH(2,NS)
              NPS = NPHAZ_BH(1,NS)
              NCHS = INDEX( PH_CND(NPS),'  ') - 1
              WRITE(ISC,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL_BH(1,NS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS_BH(1,NS)
              WRITE(IWR,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL_BH(1,NS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS_BH(1,NS)
            ENDIF
          ENDIF
        ENDIF
!
!---    Reduce time step  ---
!
        IF( NTSR.LT.4 .OR. (DTCF*DT).GT.DTMN ) THEN
          NTSR = NTSR + 1
          DTX = DT
          TM = TM - (1.D+0-DTCF)*DT
          DT = DTCF*DT
          DTO = DT
          DTI = 1.D+0/DT
          VAR = DT
          VARX = DTX
          IF( UNTM.NE.'null' ) THEN
            INDX = 1
            IUNS = 1
            CALL RDUNIT(UNTM,VAR,INDX)
            IUNS = 1
            CALL RDUNIT(UNTM,VARX,INDX)
            NCH = INDEX( UNTM,'  ')-1
          ENDIF
          WRITE(ISC,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
          WRITE(IWR,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
!
!---      Loop over fractures  ---
!
          DO NFX = 1,NF_FRC
!
!---        Loop over fracture triangles  ---
!
            DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---          Skip inactive triangles  ---
!
              IF( IXP_FRC(NTX).EQ.0 ) CYCLE
              T_FRC(2,NTX) = T_FRC(1,NTX)
              PL_FRC(2,NTX) = PL_FRC(1,NTX)
              PG_FRC(2,NTX) = PG_FRC(1,NTX)
              PVA_FRC(2,NTX) = PVA_FRC(1,NTX)
              PVW_FRC(2,NTX) = PVW_FRC(1,NTX)
              XMLA_FRC(2,NTX) = XMLA_FRC(1,NTX)
              SL_FRC(2,NTX) = SL_FRC(1,NTX)
              SG_FRC(2,NTX) = SG_FRC(1,NTX)
              YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
              TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
              NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
            ENDDO
          ENDDO
!
!---      Loop over boreholes  ---
!
          DO NBH = 1,N_BH
!
!---      Skip for boundary condition type boreholes  ---
!
          IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---        Loop over borehole nodes  ---
!
            DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
              T_BH(2,NBN) = T_BH(1,NBN)
              PL_BH(2,NBN) = PL_BH(1,NBN)
              PG_BH(2,NBN) = PG_BH(1,NBN)
              PVA_BH(2,NBN) = PVA_BH(1,NBN)
              PVW_BH(2,NBN) = PVW_BH(1,NBN)
              XMLA_BH(2,NBN) = XMLA_BH(1,NBN)
              SL_BH(2,NBN) = SL_BH(1,NBN)
              SG_BH(2,NBN) = SG_BH(1,NBN)
              YLS_BH(2,NBN) = YLS_BH(1,NBN)
              TMS_BH(2,NBN) = TMS_BH(1,NBN)
              NPHAZ_BH(2,NBN) = NPHAZ_BH(1,NBN)
            ENDDO
          ENDDO
!
!---      Loop over field nodes  ---
!
          DO N = 1,NFLD
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PVA(2,N) = PVA(1,N)
            PVW(2,N) = PVW(1,N)
            XMLA(2,N) = XMLA(1,N)
            SL(2,N) = SL(1,N)
            SG(2,N) = SG(1,N)
            SGT(2,N) = SGT(1,N)
            ASLMIN(2,N) = ASLMIN(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
          ENDDO
!
!---      Coupled-well pressure  ---
!
          IF( L_CW.EQ.1 ) THEN
            DO NCW = 1,N_CW
              P_CW(2,NCW) = P_CW(1,NCW)
            ENDDO
          ENDIF
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          ICNV = 4
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORC_BH_GT
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
!     Fracture sources and sinks from intersecting boreholes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 9 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_BH
      USE FDVT
      USE FDVP_BH
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
      REAL*8 VAR_BHX(11)
      REAL*8 KGM,KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC_BH_GT'
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
!---      Zero source terms for borehole nodes  ---
!
          DO M = 2,ISVC+2
            SRCT_BH(M,NBN) = 0.D+0
            SRCA_BH(M,NBN) = 0.D+0
            SRCW_BH(M,NBN) = 0.D+0
            SRCS_BH(M,NBN) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
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
          IF ( IT_BH(NBSE,NBH).EQ.0 ) CYCLE
!
!---      Coaxial borehole  ---
!
          IF( IT_BH(NBSE,NBH).GT.10000 ) THEN
!
!---        Borehole source is applied to first inner or outer
!           borehole node  ---
!
            NBN = ID_BH(3,NBH)+(NBSE-1)*(ID_BH(4,NBH)-ID_BH(3,NBH))
!
!---      Simple borehole  ---
!
          ELSE
!
!---        Borehole source is applied to first ID_BH(3,NBH) or last
!           ID_BH(4,NBH) borehole node  ---
!
            NBN = ID_BH(NBSE+2,NBH)
          ENDIF
          IZN = IZ_BH(NBN)
!
!---      Define adder for VAR_BH index based on starting (0)
!         or ending (5) borehole node ---
!
          NBSE_ID = 0
          IF( NBSE.EQ.2 ) NBSE_ID = 5
!
!---      Loop over increment indices  ---
!
          DO M = 2,ISVC+2
!
!---        Gas mass injection w/ zero water vapor  ---
!
            IF( MOD(IT_BH(NBSE,NBH),10000).EQ.1 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                PX = MAX( PLX,PGX )
                TX = T_BH(M,NBN)
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
                CALL WATGSH( TX,PSWX,RHOW,HGWX,UGWX )
                CALL AIRGSH( TX,PVAX,HGAX,UGAX )
                XGWX = RHOW/(RHOW+RHOA)
                XGAX = RHOA/(RHOW+RHOA)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*(XGWX*HGWX+XGAX*HGAX)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGWX
              ELSE
                HGX = XGW_BH(M,NBN)*HGW_BH(M,NBN) +
     &            XGA_BH(M,NBN)*HGA_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*HGX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGW_BH(M,NBN)
              ENDIF
!
!---        Gas mass injection w/ water vapor mass fraction  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.2 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                PX = MAX( PLX,PGX )
                TX = T_BH(M,NBN)
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
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*(XGWX*HGWX+XGAX*HGAX)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGWX
              ELSE
                HGX = XGW_BH(M,NBN)*HGW_BH(M,NBN) +
     &            XGA_BH(M,NBN)*HGA_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*HGX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGW_BH(M,NBN)
              ENDIF
!
!---        Gas mass injection w/ water vapor relative humidity  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.3 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                PX = MAX( PGX,PLX )
                TX = T_BH(M,NBN)
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
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*(XGWX*HGWX+XGAX*HGAX)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGWX
              ELSE
                HGX = XGW_BH(M,NBN)*HGW_BH(M,NBN) +
     &            XGA_BH(M,NBN)*HGA_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*HGX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XGW_BH(M,NBN)
              ENDIF
!
!---        Aqueous mass injection w/ zero air and salt mass  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.4 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                TX = T_BH(M,NBN)
                IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
                XLAX = 0.D+0
                XLSX = 0.D+0
                INDX = 0
                CALL REGION_4( TX,PSWX,INDX )
                PX = MAX( PGX,PLX,PSWX )
                CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,
     &            ULWX )
                CALL SOL_BRNS( TX,PX,XLSMX )
                XLSX = MIN( XLSMX,XLSX )
                XLWX = MAX( 1.D+0-XLAX-XLSX,0.D+0 )
                XMLAX = (XLAX/WTMA)/((XLAX/WTMA) + (XLSX/WTMS) +
     &            (XLWX/WTMW))
                CALL SP_B( TX,XLSX,PSBX )
                PVAX = HCAW*XMLAX
                XMGAX = PVAX/(PVAX+PSBX)
                CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &            XMLAX,XMLSX,XMLWX )
                CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLSX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLWX
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*HLX
              ELSE
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLS_BH(M,NBN)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLW_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*HL_BH(M,NBN)
              ENDIF
!
!---        Aqueous mass injection w/ air and salt mass fractions  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.5 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                TX = T_BH(M,NBN)
                IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
                XLAX = VAR_BHX(4+NBSE_ID)
                XLSX = VAR_BHX(6+NBSE_ID)
                INDX = 0
                CALL REGION_4( TX,PSWX,INDX )
                PX = MAX( PGX,PLX,PSWX )
                CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,
     &            ULWX )
                CALL SOL_BRNS( TX,PX,XLSMX )
                XLSX = MIN( XLSMX,XLSX )
                XLWX = MAX( 1.D+0-XLAX-XLSX,0.D+0 )
                XMLAX = (XLAX/WTMA)/((XLAX/WTMA) + (XLSX/WTMS) +
     &            (XLWX/WTMW))
                CALL SP_B( TX,XLSX,PSBX )
                PVAX = HCAW*XMLAX
                XMGAX = PVAX/(PVAX+PSBX)
                CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &            XMLAX,XMLSX,XMLWX )
                CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLSX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLWX
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*HLX
              ELSE
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLS_BH(M,NBN)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLW_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*HL_BH(M,NBN)
              ENDIF
!
!---        Aqueous mass injection w/ air and salt relative
!           saturations  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.6 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                TX = T_BH(M,NBN)
                IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
                WLAX = VAR_BHX(4+NBSE_ID)
                WLSX = VAR_BHX(6+NBSE_ID)
                INDX = 0
                CALL REGION_4( TX,PSWX,INDX )
                PX = MAX( PGX,PLX,PSWX )
                CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,
     &            ULWX )
                CALL SOL_BRNS( TX,PX,XLSMX )
                XLSX = WLSX*XLSMX
                CALL SP_B( TX,XLSX,PSBX )
                PVAX = WLAX*MAX( PX-PSBX,0.D+0 )
                XMLAX = PVAX/HCAW
                XMGAX = PVAX/(PVAX+PSBX)
                CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &            XMLAX,XMLSX,XMLWX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLSX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLWX
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*HLX
              ELSE
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLS_BH(M,NBN)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*XLW_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*HL_BH(M,NBN)
              ENDIF
!
!---        Gas volumetric injection w/ zero water vapor  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.11 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                PX = MAX( PLX,PGX )
                TX = T_BH(M,NBN)
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
                CALL WATGSH( TX,PSWX,RHOW,HGWX,UGWX )
                CALL AIRGSH( TX,PVAX,HGAX,UGAX )
                XGWX = RHOW/(RHOW+RHOA)
                XGAX = RHOA/(RHOW+RHOA)
                RHOGX = XGWX*RHOW + XGAX*RHOA
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOGX*(XGWX*HGWX+XGAX*HGAX)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOGX*XGAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOGX*XGWX
              ELSE
                HGX = XGW_BH(M,NBN)*HGW_BH(M,NBN) +
     &            XGA_BH(M,NBN)*HGA_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            RHOG_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*HGX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            RHOG_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XGA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            RHOG_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XGW_BH(M,NBN)
              ENDIF
!
!---        Gas volumetric injection w/ water vapor mass fraction  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.12 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                PX = MAX( PLX,PGX )
                TX = T_BH(M,NBN)
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
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOGX*(XGWX*HGWX+XGAX*HGAX)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOGX*XGAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOGX*XGWX
              ELSE
                HGX = XGW_BH(M,NBN)*HGW_BH(M,NBN) +
     &            XGA_BH(M,NBN)*HGA_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            RHOG_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*HGX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            RHOG_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XGA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            RHOG_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XGW_BH(M,NBN)
              ENDIF
!
!---        Gas volumetric injection w/ water vapor relative
!           humidity ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.13 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                PX = MAX( PGX,PLX )
                TX = T_BH(M,NBN)
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
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOGX*(XGWX*HGWX+XGAX*HGAX)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOGX*XGAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOGX*XGWX
              ELSE
                HGX = XGW_BH(M,NBN)*HGW_BH(M,NBN) +
     &            XGA_BH(M,NBN)*HGA_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            RHOG_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*HGX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            RHOG_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XGA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            RHOG_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XGW_BH(M,NBN)
              ENDIF
!
!---        Aqueous volumetric injection w/ zero air and
!           salt mass  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.14 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                TX = T_BH(M,NBN)
                IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
                XLAX = 0.D+0
                XLSX = 0.D+0
                INDX = 0
                CALL REGION_4( TX,PSWX,INDX )
                PX = MAX( PGX,PLX,PSWX )
                CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,
     &            ULWX )
                CALL SOL_BRNS( TX,PX,XLSMX )
                XLSX = MIN( XLSMX,XLSX )
                XLWX = MAX( 1.D+0-XLAX-XLSX,0.D+0 )
                XMLAX = (XLAX/WTMA)/((XLAX/WTMA) + (XLSX/WTMS) +
     &            (XLWX/WTMW))
                CALL SP_B( TX,XLSX,PSBX )
                PVAX = HCAW*XMLAX
                XMGAX = PVAX/(PVAX+PSBX)
                CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &            XMLAX,XMLSX,XMLWX )
                CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOLX*XLSX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOLX*XLAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &             VAR_BHX(2+NBSE_ID)*RHOLX*XLWX
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &             VAR_BHX(2+NBSE_ID)*RHOLX*HLX
              ELSE
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XLS_BH(M,NBN)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XLA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XLW_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*HL_BH(M,NBN)
              ENDIF
!
!---        Aqueous volumetric injection w/ air and
!           salt mass fractions  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.15 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                TX = T_BH(M,NBN)
                IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
                XLAX = VAR_BHX(4+NBSE_ID)
                XLSX = VAR_BHX(6+NBSE_ID)
                INDX = 0
                CALL REGION_4( TX,PSWX,INDX )
                PX = MAX( PGX,PLX,PSWX )
                CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,
     &            ULWX )
                CALL SOL_BRNS( TX,PX,XLSMX )
                XLSX = MIN( XLSMX,XLSX )
                XLWX = MAX( 1.D+0-XLAX-XLSX,0.D+0 )
                XMLAX = (XLAX/WTMA)/((XLAX/WTMA) + (XLSX/WTMS) +
     &            (XLWX/WTMW))
                CALL SP_B( TX,XLSX,PSBX )
                PVAX = HCAW*XMLAX
                XMGAX = PVAX/(PVAX+PSBX)
                CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &            XMLAX,XMLSX,XMLWX )
                CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOLX*XLSX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOLX*XLAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2)*RHOLX*XLWX
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2)*RHOLX*HLX
              ELSE
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XLS_BH(M,NBN)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XLA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XLW_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*HL_BH(M,NBN)
              ENDIF
!
!---        Aqueous volumetric injection w/ air and salt
!           relative saturations  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.16 ) THEN
              PGX = PG_BH(M,NBN) + PATM
              PLX = PL_BH(M,NBN) + PATM
              IF( VAR_BHX(2+NBSE_ID).GE.0.D+0 ) THEN
                TX = T_BH(M,NBN)
                IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
                WLAX = VAR_BHX(4+NBSE_ID)
                WLSX = VAR_BHX(6+NBSE_ID)
                INDX = 0
                CALL REGION_4( TX,PSWX,INDX )
                PX = MAX( PGX,PLX,PSWX )
                CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,
     &            ULWX )
                CALL SOL_BRNS( TX,PX,XLSMX )
                XLSX = WLSX*XLSMX
                CALL SP_B( TX,XLSX,PSBX )
                PVAX = WLAX*MAX( PX-PSBX,0.D+0 )
                XMLAX = PVAX/HCAW
                XMGAX = PVAX/(PVAX+PSBX)
                CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &            XMLAX,XMLSX,XMLWX )
                CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOLX*XLSX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOLX*XLAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOLX*XLWX
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            VAR_BHX(2+NBSE_ID)*RHOLX*HLX
              ELSE
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XLS_BH(M,NBN)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XLA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*XLW_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) +
     &            RHOL_BH(M,NBN)*VAR_BHX(2+NBSE_ID)*HL_BH(M,NBN)
              ENDIF
!
!---        Pressure borehole, state condition #1 (saturated)  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.7 ) THEN
!
!---          Distance between borehole node centroid and
!             pressure boundary  ---
!
              DBBX = 5.D-1*SQRT( (XP_BH(1,NBN)-XP_BH(2,NBN))**2 +
     &          (YP_BH(1,NBN)-YP_BH(2,NBN))**2 +
     &          (ZP_BH(1,NBN)-ZP_BH(2,NBN))**2 )
!
!---          Borehole pressure  ---
!
              PL_BHX = VAR_BHX(3+NBSE_ID)
              PG_BHX = PL_BHX
              RKL_BHX = 1.D+0
              INV = INV_BH(NBN)
!
!---          Inner coaxial borehole  ---
!
              IF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.1 ) THEN
                AF_BHX = GPI*(PAR_BH(1,INV)**2)
                ZP_BHX = ZP_BH(1,NBN)
!
!---          Outer coaxial borehole  ---
!
              ELSEIF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.2 ) THEN
                AF_BHX = GPI*((PAR_BH(3,INV)**2)-(PAR_BH(2,INV)**2))
                ZP_BHX = ZP_BH(2,NBN)
!
!---          Cased or pipe borehole  ---
!
              ELSEIF( IS_BH(INV).GT.10 ) THEN
                AF_BHX = GPI*(PAR_BH(3,INV)**2)
                ZP_BHX = ZP_BH(1,NBN)
!
!---          Uncased  ---
!
              ELSE
                AF_BHX = GPI*(PAR_BH(2,INV)**2)
                ZP_BHX = ZP_BH(1,NBN)
              ENDIF
!
!---          Aqueous flux (m^3/s) into first borehole node of
!             borehole  ---
!
              PLX = PL_BH(M,NBN) + PATM
              PGX = PG_BH(M,NBN) + PATM
              ZPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
              HDLX = PL_BHX - PLX +
     &          GRAV*(ZP_BHX-ZPX)*RHOL_BH(M,NBN)
              IF( M.EQ.1 ) HDL = HDLX
!
!---          Pipe flow  ---
!
              INV = INV_BH(NBN)
              IF( IS_BH(INV).GT.10 ) THEN
!
!---            Inner coaxial borehole  ---
!
                IF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.1 ) THEN
!
!---              Hydraulic diameter (m) ---
!
                  DHX = 2.D+0*PAR_BH(1,INV)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(13,INV)
!
!---            Outer coaxial borehole  ---
!
                ELSEIF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.2 ) THEN
!
!---              Hydraulic diameter (m) ---
!
                  DHX = 2.D+0*SQRT(PAR_BH(3,INV)**2 - PAR_BH(2,INV)**2)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(16,INV)
!
!---            Cased or pipe borehole  ---
!
                ELSE
!
!---              Hydraulic diameter (m), sum of interval radii ---
!
                  DHX = 2.D+0*PAR_BH(3,INV)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(6,INV)
                ENDIF
!
!---            Pipe flow velocity (m/s), first estimating from
!               Hazen-Williamsequation, adjusting for fluid viscosity,
!               then solving a nonlinear system of equations,
!               combining the Colebrook equation for the friction
!               coefficient and the Darcy-Weisbach equation for the
!               pressure drop (i.e., head loss)  ---
!
                VISLX = VISL_BH(M,NBN)
                CALL PIPE_FLOW_VEL_GT( DBBX,DHX,HDLX,RFX,RHOL_BH(M,NBN),
     &            ULX,VISLX,NBN )
                IF( VISLX.LT.0.D+0 ) THEN
                  print *,'pl_bhx = ',PL_BHX
                  print *,'plx = ',PLX
                  print *,'zp_bhx = ',ZP_BHX
                  print *,'zpx = ',ZPX
                  print *,'rhol_bh(m,nbn) = ',RHOL_BH(M,NBN)
                  print *,'nbn = ',NBN
                ENDIF
                INDX = 8
                RKLM = DIFMN(RKL_BHX,RKL_BH(M,NBN),DBBX,DBBX,HDL,INDX)
                FLX = AF_BHX*RKLM*SIGN(ULX,HDLX)
!
!---          Darcy flow ---
!
              ELSE
                KLM = PERM(1,IZN)
                IF( PERM(1,IZN)/EPSL.LT.EPSL ) KLM = 0.D+0
                INDX = 8
                RKLM = DIFMN(RKL_BHX,RKL_BH(M,NBN),DBBX,DBBX,HDL,INDX)
                VLM = VISL_BH(M,NBN)
                FLX = AF_BHX*KLM*RKLM*HDLX/DBBX/VLM
              ENDIF
!
!---          Aqueous flux into borehole node  ---
!
              IF( FLX.GE.0.D+0 ) THEN
                TX = T_BH(M,NBN)
                IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
                WLAX = VAR_BHX(4+NBSE_ID)
                WLSX = VAR_BHX(6+NBSE_ID)
                INDX = 0
                CALL REGION_4( TX,PSWX,INDX )
                PX = MAX( PL_BHX,PSWX )
                CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,
     &            ULWX )
                CALL SOL_BRNS( TX,PX,XLSMX )
                XLSX = WLSX*XLSMX
                CALL SP_B( TX,XLSX,PSBX )
                PVAX = WLAX*MAX( PX-PSBX,0.D+0 )
                XMLAX = PVAX/HCAW
                XMGAX = PVAX/(PVAX+PSBX)
                CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &            XMLAX,XMLSX,XMLWX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
                RHOLX = RHOBX
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) + FLX*RHOLX*XLSX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) + FLX*RHOLX*XLAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) + FLX*RHOLX*XLWX
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) + FLX*RHOLX*HLX
!
!---          Aqueous flux from borehole node  ---
!
              ELSE
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) + FLX*RHOL_BH(M,NBN)*
     &            XLS_BH(M,NBN)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) + FLX*RHOL_BH(M,NBN)*
     &            XLA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) + FLX*RHOL_BH(M,NBN)*
     &            XLW_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) + FLX*RHOL_BH(M,NBN)*
     &            HL_BH(M,NBN)
              ENDIF
!
!---        Pressure borehole, state condition #2
!           (partially saturated)  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.8 ) THEN
!
!---          Distance between borehole node centroid and
!             pressure boundary  ---
!
              DBBX = 5.D-1*SQRT( (XP_BH(1,NBN)-XP_BH(2,NBN))**2 +
     &          (YP_BH(1,NBN)-YP_BH(2,NBN))**2 +
     &          (ZP_BH(1,NBN)-ZP_BH(2,NBN))**2 )
!
!---          Borehole pressure  ---
!
              PG_BHX = VAR_BHX(3+NBSE_ID)
              SL_BHX = VAR_BHX(2+NBSE_ID)
              SG_BHX = 1.D+0 - SL_BHX
              CALL CAP_BH_GT( CPGLX,SL_BHX,NBN )
              PL_BHX = PG_BHX - CPGLX
              CALL RKSP_BH_GT( PG_BHX,PL_BHX,RKG_BHX,RKL_BHX,SGX,SLX,
     &          NBN )
              INV = INV_BH(NBN)
!
!---          Inner coaxial borehole  ---
!
              IF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.1 ) THEN
                AF_BHX = GPI*(PAR_BH(1,INV)**2)
                ZP_BHX = ZP_BH(1,NBN)
!
!---          Outer coaxial borehole  ---
!
              ELSEIF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.2 ) THEN
                AF_BHX = GPI*((PAR_BH(3,INV)**2)-(PAR_BH(2,INV)**2))
                ZP_BHX = ZP_BH(2,NBN)
!
!---          Cased or pipe borehole  ---
!
              ELSEIF( IS_BH(INV).GT.10 ) THEN
                AF_BHX = GPI*(PAR_BH(3,INV)**2)
                ZP_BHX = ZP_BH(1,NBN)
!
!---          Uncased  ---
!
              ELSE
                AF_BHX = GPI*(PAR_BH(2,INV)**2)
                ZP_BHX = ZP_BH(1,NBN)
              ENDIF
!
!---          Aqueous flux (m^3/s) into first borehole node of
!             borehole  ---
!
              PLX = PL_BH(M,NBN) + PATM
              ZPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
              HDLX = PL_BHX - PLX +
     &          GRAV*(ZP_BHX-ZPX)*RHOL_BH(M,NBN)
              IF( M.EQ.1 ) HDL = HDLX
!
!---          Pipe flow  ---
!
              INV = INV_BH(NBN)
              IF( IS_BH(INV).GT.10 ) THEN
!
!---            Inner coaxial borehole  ---
!
                IF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.1 ) THEN
!
!---              Hydraulic diameter (m) ---
!
                  DHX = 2.D+0*PAR_BH(1,INV)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(13,INV)
!
!---            Outer coaxial borehole  ---
!
                ELSEIF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.2 ) THEN
!
!---              Hydraulic diameter (m) ---
!
                  DHX = 2.D+0*SQRT(PAR_BH(3,INV)**2-PAR_BH(2,INV)**2)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(16,INV)
!
!---            Cased or pipe borehole  ---
!
                ELSE
!
!---              Hydraulic diameter (m), sum of interval radii ---
!
                  DHX = 2.D+0*PAR_BH(3,INV)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(6,INV)
                ENDIF
!
!---            Pipe flow velocity (m/s), first estimating from
!               Hazen-Williamsequation, adjusting for fluid viscosity,
!               then solving a nonlinear system of equations,
!               combining the Colebrook equation for the friction
!               coefficient and the Darcy-Weisbach equation for the
!               pressure drop (i.e., head loss)  ---
!
                VISLX = VISL_BH(M,NBN)
                CALL PIPE_FLOW_VEL_GT( DBBX,DHX,HDLX,RFX,RHOL_BH(M,NBN),
     &            ULX,VISLX,NBN )
                IF( VISLX.LT.0.D+0 ) THEN
                  print *,'pl_bhx = ',PL_BHX
                  print *,'plx = ',PLX
                  print *,'zp_bhx = ',ZP_BHX
                  print *,'zpx = ',ZPX
                  print *,'rhol_bh(m,nbn) = ',RHOL_BH(M,NBN)
                  print *,'nbn = ',NBN
                ENDIF
                INDX = 8
                RKLM = DIFMN(RKL_BHX,RKL_BH(M,NBN),DBBX,DBBX,HDL,INDX)
                FLX = AF_BHX*RKLM*SIGN(ULX,HDLX)
!
!---          Darcy flow ---
!
              ELSE
                KLM = PERM(1,IZN)
                IF( PERM(1,IZN)/EPSL.LT.EPSL ) KLM = 0.D+0
                INDX = 8
                RKLM = DIFMN(RKL_BHX,RKL_BH(M,NBN),DBBX,DBBX,HDL,INDX)
                VLM = VISL_BH(M,NBN)
                FLX = AF_BHX*KLM*RKLM*HDLX/DBBX/VLM
              ENDIF
!
!---          Gas flux (m^3/s) into first borehole node of
!             borehole  ---
!
              PGX = PG_BH(M,NBN) + PATM
              ZPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
              HDGX = PG_BHX - PGX +
     &          GRAV*(ZP_BHX-ZPX)*RHOG_BH(M,NBN)
              IF( M.EQ.1 ) HDG = HDGX
!
!---          Pipe flow  ---
!
              INV = INV_BH(NBN)
              IF( IS_BH(INV).GT.10 ) THEN
!
!---            Inner coaxial borehole  ---
!
                IF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.1 ) THEN
!
!---              Hydraulic diameter (m) ---
!
                  DHX = 2.D+0*PAR_BH(1,INV)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(13,INV)
!
!---            Outer coaxial borehole  ---
!
                ELSEIF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.2 ) THEN
!
!---              Hydraulic diameter (m) ---
!
                  DHX = 2.D+0*SQRT(PAR_BH(3,INV)**2-PAR_BH(2,INV)**2)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(16,INV)
!
!---            Cased or pipe borehole  ---
!
                ELSE
!
!---              Hydraulic diameter (m), sum of interval radii ---
!
                  DHX = 2.D+0*PAR_BH(3,INV)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(6,INV)
                ENDIF
!
!---            Pipe flow velocity (m/s), first estimating from
!               Hazen-Williamsequation, adjusting for fluid viscosity,
!               then solving a nonlinear system of equations,
!               combining the Colebrook equation for the friction
!               coefficient and the Darcy-Weisbach equation for the
!               pressure drop (i.e., head loss)  ---
!
                VISGX = VISG_BH(M,NBN)
                CALL PIPE_FLOW_VEL_GT( DBBX,DHX,HDGX,RFX,RHOG_BH(M,NBN),
     &            UGX,VISGX,NBN )
                IF( VISGX.LT.0.D+0 ) THEN
                  print *,'pg_bhx = ',PG_BHX
                  print *,'pgx = ',PGX
                  print *,'zp_bhx = ',ZP_BHX
                  print *,'zpx = ',ZPX
                  print *,'rhog_bh(m,nbn) = ',RHOG_BH(M,NBN)
                  print *,'nbn = ',NBN
                ENDIF
                INDX = 8
                RKGM = DIFMN(RKG_BHX,RKG_BH(M,NBN),DBBX,DBBX,HDG,INDX)
                FGX = AF_BHX*RKGM*SIGN(UGX,HDGX)
!
!---          Darcy flow ---
!
              ELSE
                KGM = PERM(1,IZN)
                IF( PERM(1,IZN)/EPSL.LT.EPSL ) KLM = 0.D+0
                INDX = 8
                RKGM = DIFMN(RKG_BHX,RKG_BH(M,NBN),DBBX,DBBX,HDG,INDX)
                VGM = VISG_BH(M,NBN)
                FGX = AF_BHX*KGM*RKGM*HDGX/DBBX/VGM
              ENDIF
!
!---          Influx properties  ---
!
              TX = T_BH(M,NBN)
              IF( ISLC(30).EQ.0 ) TX = VAR_BHX(5+NBSE_ID)
              WLSX = VAR_BHX(6+NBSE_ID)
              INDX = 0
              CALL REGION_4( TX,PSWX,INDX )
              PX = MAX( PG_BHX,PL_BHX,PSWX )
              PCX = PG_BHX-PL_BHX
              CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
              CALL SOL_BRNS( TX,PX,XLSMX )
              XLSX = WLSX*XLSMX
              CALL SP_B( TX,XLSX,PSBX )
              CALL VPL_B( TX,PSBX,PCX,RHOLWX,PVWX )
              PVAX = WLAX*MAX( PX-PVWX,0.D+0 )
              IF( PVAX.LT.1.D-6 ) PVAX = 0.D+0
              XMLAX = PVAX/HCAW
              XMGAX = PVAX/(PVAX+PVWX)
              CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &          XMLAX,XMLSX,XMLWX )
              CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
              CALL ENTH_B( TX,XLSX,HLWX,HBX )
              HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
              HGX = XGAX*HGAX + XGWX*HGWX
              CALL AIRGSD( TX,PVAX,RHOGAX )
              RHOGX = XGAX*RHOGAX + XGWX*RHOGWX
!
!---          Aqueous flux into borehole node  ---
!
              IF( FLX.GE.0.D+0 ) THEN
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) + FLX*RHOLX*XLSX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) + FLX*RHOLX*XLAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) + FLX*RHOLX*XLWX
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) + FLX*RHOLX*HLX
!
!---          Aqueous flux from borehole node  ---
!
              ELSE
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) + FLX*RHOL_BH(M,NBN)*
     &            XLS_BH(M,NBN)
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) + FLX*RHOL_BH(M,NBN)*
     &            XLA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) + FLX*RHOL_BH(M,NBN)*
     &            XLW_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) + FLX*RHOL_BH(M,NBN)*
     &            HL_BH(M,NBN)
              ENDIF
!
!---          Gas flux into borehole node  ---
!
              IF( FGX.GE.0.D+0 ) THEN
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) + FGX*RHOGX*XGAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) + FGX*RHOGX*XGWX
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) + FGX*RHOGX*HGX
!
!---          Gas flux from borehole node  ---
!
              ELSE
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) + FGX*RHOG_BH(M,NBN)*
     &            XGA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) + FGX*RHOG_BH(M,NBN)*
     &            XGW_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) + FGX*RHOG_BH(M,NBN)*
     &            HG_BH(M,NBN)
              ENDIF
!
!---        Pressure borehole, state condition #3 (unsaturated)  ---
!
            ELSEIF( MOD(IT_BH(NBSE,NBH),10000).EQ.9 ) THEN
!
!---          Distance between borehole node centroid and
!             pressure boundary  ---
!
              DBBX = 5.D-1*SQRT( (XP_BH(1,NBN)-XP_BH(2,NBN))**2 +
     &          (YP_BH(1,NBN)-YP_BH(2,NBN))**2 +
     &          (ZP_BH(1,NBN)-ZP_BH(2,NBN))**2 )
!
!---          Borehole pressure  ---
!
              PG_BHX = VAR_BHX(3+NBSE_ID)
              SL_BHX = 0.D+0
              SG_BHX = 1.D+0 - SL_BHX
              RKG_BHX = 1.D+0
              INV = INV_BH(NBN)
!
!---          Inner coaxial borehole  ---
!
              IF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.1 ) THEN
                AF_BHX = GPI*(PAR_BH(1,INV)**2)
                ZP_BHX = ZP_BH(1,NBN)
!
!---          Outer coaxial borehole  ---
!
              ELSEIF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.2 ) THEN
                AF_BHX = GPI*((PAR_BH(3,INV)**2)-(PAR_BH(2,INV)**2))
                ZP_BHX = ZP_BH(2,NBN)
!
!---          Cased or pipe borehole  ---
!
              ELSEIF( IS_BH(INV).GT.10 ) THEN
                AF_BHX = GPI*(PAR_BH(3,INV)**2)
                ZP_BHX = ZP_BH(1,NBN)
!
!---          Uncased  ---
!
              ELSE
                AF_BHX = GPI*(PAR_BH(2,INV)**2)
                ZP_BHX = ZP_BH(1,NBN)
              ENDIF
!
!---          Gas flux (m^3/s) into first borehole node of
!             borehole  ---
!
              PGX = PG_BH(M,NBN) + PATM
              ZPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
              HDGX = PG_BHX - PGX +
     &          GRAV*(ZP_BHX-ZPX)*RHOG_BH(M,NBN)
              IF( M.EQ.1 ) HDG = HDGX
!
!---          Pipe flow  ---
!
              INV = INV_BH(NBN)
              IF( IS_BH(INV).GT.10 ) THEN
!
!---            Inner coaxial borehole  ---
!
                IF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.1 ) THEN
!
!---              Hydraulic diameter (m) ---
!
                  DHX = 2.D+0*PAR_BH(1,INV)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(13,INV)
!
!---            Outer coaxial borehole  ---
!
                ELSEIF( IS_BH(INV).EQ.10000 .AND. NBSE.EQ.2 ) THEN
!
!---              Hydraulic diameter (m) ---
!
                  DHX = 2.D+0*SQRT(PAR_BH(3,INV)**2-PAR_BH(2,INV)**2)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(16,INV)
!
!---            Cased or pipe borehole  ---
!
                ELSE
!
!---              Hydraulic diameter (m), sum of interval radii ---
!
                  DHX = 2.D+0*PAR_BH(3,INV)
!
!---              Roughness, m  ---
!
                  RFX = PAR_BH(6,INV)
                ENDIF
!
!---            Pipe flow velocity (m/s), first estimating from
!               Hazen-Williamsequation, adjusting for fluid viscosity,
!               then solving a nonlinear system of equations,
!               combining the Colebrook equation for the friction
!               coefficient and the Darcy-Weisbach equation for the
!               pressure drop (i.e., head loss)  ---
!
                VISGX = VISG_BH(M,NBN)
                CALL PIPE_FLOW_VEL_GT( DBBX,DHX,HDGX,RFX,RHOG_BH(M,NBN),
     &            UGX,VISGX,NBN )
                IF( VISGX.LT.0.D+0 ) THEN
                  print *,'pg_bhx = ',PG_BHX
                  print *,'pgx = ',PGX
                  print *,'zp_bhx = ',ZP_BHX
                  print *,'zpx = ',ZPX
                  print *,'rhog_bh(m,nbn) = ',RHOG_BH(M,NBN)
                  print *,'nbn = ',NBN
                ENDIF
                INDX = 8
                RKGM = DIFMN(RKG_BHX,RKG_BH(M,NBN),DBBX,DBBX,HDG,INDX)
                FGX = AF_BHX*RKGM*SIGN(UGX,HDGX)
!
!---          Darcy flow ---
!
              ELSE
                KGM = PERM(1,IZN)
                IF( PERM(1,IZN)/EPSL.LT.EPSL ) KLM = 0.D+0
                INDX = 8
                RKGM = DIFMN(RKG_BHX,RKG_BH(M,NBN),DBBX,DBBX,HDG,INDX)
                VGM = VISG_BH(M,NBN)
                FGX = AF_BHX*KGM*RKGM*HDGX/DBBX/VGM
              ENDIF
!
!---          Gas flux into borehole node  ---
!
              IF( FGX.GE.0.D+0 ) THEN
                PX = PG_BHX
                TX = T_BH(M,NBN)
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
                CALL AIRGSD( TX,PVAX,RHOGAX )
                RHOGX = XGAX*RHOGAX + XGWX*RHOGWX
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) + FGX*RHOGX*XGAX
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) + FGX*RHOGX*XGWX
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) + FGX*RHOGX*HGX
!
!---          Gas flux from borehole node  ---
!
              ELSE
                SRCA_BH(M,NBN) = SRCA_BH(M,NBN) + FGX*RHOG_BH(M,NBN)*
     &            XGA_BH(M,NBN)
                SRCW_BH(M,NBN) = SRCW_BH(M,NBN) + FGX*RHOG_BH(M,NBN)*
     &            XGW_BH(M,NBN)
                SRCT_BH(M,NBN) = SRCT_BH(M,NBN) + FGX*RHOG_BH(M,NBN)*
     &            HG_BH(M,NBN)
              ENDIF
!
!---        Dirichlet energy borehole type  ---
!
            ELSEIF( IT_BH(NBSE,NBH).EQ.21 ) THEN
!
!---          Borehole collar coordinates  ---
!
              DISTX = 0.D+0
              NBN = ID_BH(3,NBH)
              XP_BH0X = XP_BH(1,NBN)
              YP_BH0X = YP_BH(1,NBN)
              ZP_BH0X = ZP_BH(1,NBN)
!
!---          Loop over borehole nodes  ---
!
              DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---            Field node connected to borehole node  ---
!
                N = IBN_BH(NBN)
!
!---            Borehole node projections  ---
!
                XLX = ABS(XP_BH(1,NBN)-XP_BH(2,NBN))
                YLX = ABS(YP_BH(1,NBN)-YP_BH(2,NBN))
                ZLX = ABS(ZP_BH(1,NBN)-ZP_BH(2,NBN))
!
!---            Borehole node centroid and surface area  ---
!
                XP_BHX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
                YP_BHX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
                ZP_BHX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
                INV = INV_BH(NBN)
                IF( IS_BH(INV).GT.10 ) THEN
                  RBX = MAX( PAR_BH(3,INV_BH(NBN)),1.D-20 )
                ELSE
                  RBX = MAX( PAR_BH(2,INV_BH(NBN)),1.D-20 )
                ENDIF
                AFB_BHX = 2.D+0*GPI*RBX*SQRT(
     &            (XP_BH(2,NBN)-XP_BH(1,NBN))**2 +
     &            (YP_BH(2,NBN)-YP_BH(1,NBN))**2 +
     &            (ZP_BH(2,NBN)-ZP_BH(1,NBN))**2 )
!
!---            Distance from borehole collar to borehole node
!               centroid  ---
!
                DISTX = DISTX + SQRT( (XP_BHX-XP_BH0X)**2 +
     &            (YP_BHX-YP_BH0X)**2 + (ZP_BHX-ZP_BH0X)**2 )
                XP_BH0X = XP_BHX
                YP_BH0X = YP_BHX
                ZP_BH0X = ZP_BHX
!
!---            Borehole node temperature  ---
!
                T_BH(M,NBN) = VAR_BHX(5+NBSE_ID) +
     &            DISTX*VAR_BHX(6+NBSE_ID)
!
!---            Thermal conductivity  ---
!
                THKXY = SQRT(THKS(1,IZ(N))*THKS(2,IZ(N)))
                THKXZ = SQRT(THKS(1,IZ(N))*THKS(3,IZ(N)))
                THKYZ = SQRT(THKS(2,IZ(N))*THKS(3,IZ(N)))
                THKX = SQRT((XLX*THKYZ)**2 + (YLX*THKXZ)**2 +
     &            (ZLX*THKXY)**2)/SQRT(XLX**2 + YLX**2 + ZLX**2)
                THKX = THKX*MAX(1.D+0-PORD(M,N),0.D+0) +
     &            PORD(M,N)*(THKL(M,N)*SL(M,N) + THKG(M,N)*SG(M,N))
!
!---            Conduction heat transfer coefficient, W/m^2 ---
!
                HCX = THKX/DBN_BH(NBN)
!
!---            Heat flux into field node  ---
!
                SRCT(M,N) = SRCT(M,N) +
     &            AFB_BHX*HCX*(T_BH(M,NBN)-T(M,N))
              ENDDO
!
!---        Neumann energy borehole type  ---
!
            ELSEIF( IT_BH(NBSE,NBH).EQ.22 ) THEN
!
!---          Borehole collar coordinates  ---
!
              DISTX = 0.D+0
              NBN = ID_BH(3,NBH)
              XP_BH0X = XP_BH(1,NBN)
              YP_BH0X = YP_BH(1,NBN)
              ZP_BH0X = ZP_BH(1,NBN)
!
!---          Loop over borehole nodes  ---
!
              DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---            Field node connected to borehole node  ---
!
                N = IBN_BH(NBN)
!
!---            Borehole node centroid and surface area  ---
!
                XP_BHX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
                YP_BHX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
                ZP_BHX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
                INV = INV_BH(NBN)
                IF( IS_BH(INV).GT.10 ) THEN
                  RBX = MAX( PAR_BH(3,INV_BH(NBN)),1.D-20 )
                ELSE
                  RBX = MAX( PAR_BH(2,INV_BH(NBN)),1.D-20 )
                ENDIF
                AFB_BHX = 2.D+0*GPI*RBX*SQRT(
     &            (XP_BH(2,NBN)-XP_BH(1,NBN))**2 +
     &            (YP_BH(2,NBN)-YP_BH(1,NBN))**2 +
     &            (ZP_BH(2,NBN)-ZP_BH(1,NBN))**2 )
!
!---            Distance from borehole collar to borehole node
!               centroid  ---
!
                DISTX = DISTX + SQRT( (XP_BHX-XP_BH0X)**2 +
     &            (YP_BHX-YP_BH0X)**2 + (ZP_BHX-ZP_BH0X)**2 )
                XP_BH0X = XP_BHX
                YP_BH0X = YP_BHX
                ZP_BH0X = ZP_BHX
!
!---            Borehole node heat flux  ---
!
                Q_BHX = VAR_BHX(5+NBSE_ID) + DISTX*VAR_BHX(6+NBSE_ID)
!
!---            Heat flux into field node  ---
!
                SRCT(M,N) = SRCT(M,N) + AFB_BHX*Q_BHX
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Thermocatalytic reactions with salt tracking ---
!
      IF( ITC_BH.GT.20 .AND. ITC_BH.LE.30 ) THEN
!
!---    Loop over boreholes  ---
!
        DO NBH = 1,N_BH
!
!---      Loop over borehole nodes  ---
!
          DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
!
!---        Thermo-catalytic fluid, temperature-transition
!           equation via a first-order rate form, (J/kg aqu)/(m^3 aqu)
!           with salt tracking  ---
!
            IF( ITC_BH.EQ.22 ) THEN
              DO M = 2,ISVC+2
!
!---            Salt mass source equal to 1.e-4 times water mass in
!               borehole node  ---
!
                VARX = XLS_BH(M,NBN)*VOL_BH(NBN)*
     &            SL_BH(M,NBN)*RHOL_BH(M,NBN)*PORD_BH(M,NBN)
                VAR1X = PTC_BH(2)
                VAR2X = 1.D+0/(1.D+0 + EXP((PTC_BH(3)-T_BH(M,NBN))/
     &            PTC_BH(4)))
                SRCS_BH(M,NBN) = SRCS_BH(M,NBN) - VARX*VAR1X*VAR2X
              ENDDO
!
!---        Thermo-catalytic fluid, temperature-transition
!           equation via a first-order rate form, (J/kg aqu)/(m^3 aqu)
!           with salt tracking, under catalyst contact only  ---
!
            ELSEIF( ITC_BH.EQ.24 ) THEN
              INVX = INV_BH(NBN)
!
!---          Borehole node within catalyst interval  ---
!
              IF( IS_BH(INVX).EQ.21 ) THEN
                DO M = 2,ISVC+2
!
!---              Salt mass source equal to 1.e-4 times water mass in
!                 borehole node  ---
!
                  VARX = XLS_BH(M,NBN)*VOL_BH(NBN)*
     &              SL_BH(M,NBN)*RHOL_BH(M,NBN)*PORD_BH(M,NBN)
                  VAR1X = PTC_BH(2)
                  VAR2X = 1.D+0/(1.D+0 + EXP((PTC_BH(3)-T_BH(M,NBN))/
     &              PTC_BH(4)))
                  SRCS_BH(M,NBN) = SRCS_BH(M,NBN) - VARX*VAR1X*VAR2X
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORC_BH_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORC_FRC_GT
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
!     Compute fracture source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 10 October 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PARM_FRC
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE FDVT_FRC
      USE FDVS_FRC
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
      REAL*8 SRX(7+LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC_FRC_GT'
!
!---  Zero fracture source terms, looping over fractures  ---
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
          DO M = 2,ISVC+2
            SRCT_FRC(M,NTX) = 0.D+0
            SRCA_FRC(M,NTX) = 0.D+0
            SRCW_FRC(M,NTX) = 0.D+0
            SRCS_FRC(M,NTX) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
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
          DO M = 2,ISVC+2
            PGX = PG_FRC(M,NTX) + PATM
            PLX = PL_FRC(M,NTX) + PATM
!
!---        Aqueous Volumetric Rate w/ Relative Saturations ---
!
            IF( ISRT_FRC(NS).EQ.11 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                TX = SRX(2)
                WLAX = SRX(5)
                WLSX = SRX(6)
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
     &            XMLAX,XMLSX,XMLWX )
                CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + SRX(4)*RHOLX*XLSX
                SRCA_FRC(M,NTX) = SRCA_FRC(M,NTX) + SRX(4)*RHOLX*XLAX
                SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) + SRX(4)*RHOLX*XLWX
                SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) + SRX(4)*RHOLX*HLX
                IF( M.EQ.2 ) SRCL_FRC(NTX) = SRX(4)*RHOLX
              ELSE
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) +
     &            SRX(4)*RHOL_FRC(M,NTX)*XLS_FRC(M,NTX)
                SRCA_FRC(M,NTX) = SRCA_FRC(M,NTX) +
     &            SRX(4)*RHOL_FRC(M,NTX)*XLA_FRC(M,NTX)
                SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) +
     &            SRX(4)*RHOL_FRC(M,NTX)*XLW_FRC(M,NTX)
                SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) +
     &            SRX(4)*RHOL_FRC(M,NTX)*HL_FRC(M,NTX)
                IF( M.EQ.2 ) SRCL_FRC(NTX) = SRX(4)*RHOL_FRC(M,NTX)
              ENDIF
!
!---        Aqueous Volumetric Rate w/ Mass Fractions ---
!
            ELSEIF( ISRT_FRC(NS).EQ.12 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                TX = SRX(2)
                XLAX = SRX(5)
                XLSX = SRX(6)
                INDX = 0
                CALL REGION_4( TX,PSWX,INDX )
                PX = MAX( PGX,PLX,PSWX )
                CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
                CALL SOL_BRNS( TX,PX,XLSMX )
                XLSX = MIN( XLSMX,XLSX )
                XLWX = MAX( 1.D+0-XLAX-XLSX,0.D+0 )
                XMLAX = (XLAX/WTMA)/((XLAX/WTMA) + (XLSX/WTMS) +
     &            (XLWX/WTMW))
                CALL SP_B( TX,XLSX,PSBX )
                PVAX = HCAW*XMLAX
                XMGAX = PVAX/(PVAX+PSBX)
                CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &            XMLAX,XMLSX,XMLWX )
                CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + SRX(4)*RHOLX*XLSX
                SRCA_FRC(M,NTX) = SRCA_FRC(M,NTX) + SRX(4)*RHOLX*XLAX
                SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) + SRX(4)*RHOLX*XLWX
                SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) + SRX(4)*RHOLX*HLX
                IF( M.EQ.2 ) SRCL_FRC(NTX) = SRX(4)*RHOLX
              ELSE
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) +
     &            SRX(4)*RHOL_FRC(M,NTX)*XLS_FRC(M,NTX)
                SRCA_FRC(M,NTX) = SRCA_FRC(M,NTX) +
     &            SRX(4)*RHOL_FRC(M,NTX)*XLA_FRC(M,NTX)
                SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) +
     &            SRX(4)*RHOL_FRC(M,NTX)*XLW_FRC(M,NTX)
                SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) +
     &            SRX(4)*RHOL_FRC(M,NTX)*HL_FRC(M,NTX)
                IF( M.EQ.2 ) SRCL_FRC(NTX) = SRX(4)*RHOL_FRC(M,NTX)
              ENDIF
!
!---        Aqueous Mass Rate w/ Relative Saturations ---
!
            ELSEIF( ISRT_FRC(NS).EQ.13 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                TX = SRX(2)
                WLAX = SRX(5)
                WLSX = SRX(6)
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
     &            XMLAX,XMLSX,XMLWX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + SRX(4)*XLSX
                SRCA_FRC(M,NTX) = SRCA_FRC(M,NTX) + SRX(4)*XLAX
                SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) + SRX(4)*XLWX
                SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) + SRX(4)*HLX
              ELSE
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) +
     &            SRX(4)*XLS_FRC(M,NTX)
                SRCA_FRC(M,NTX) = SRCA_FRC(M,NTX) +
     &            SRX(4)*XLA_FRC(M,NTX)
                SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) +
     &            SRX(4)*XLW_FRC(M,NTX)
                SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) +
     &            SRX(4)*HL_FRC(M,NTX)
              ENDIF
              IF( M.EQ.2 ) SRCL_FRC(NTX) = SRX(4)
!
!---        Aqueous Mass Rate w/ Mass Fractions ---
!
            ELSEIF( ISRT_FRC(NS).EQ.14 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                TX = SRX(2)
                XLAX = SRX(5)
                XLSX = SRX(6)
                INDX = 0
                CALL REGION_4( TX,PSWX,INDX )
                PX = MAX( PGX,PLX,PSWX )
                CALL P_IAPWS( TX,PX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
                CALL SOL_BRNS( TX,PX,XLSMX )
                XLSX = MIN( XLSMX,XLSX )
                XLWX = MAX( 1.D+0-XLAX-XLSX,0.D+0 )
                XMLAX = (XLAX/WTMA)/((XLAX/WTMA) + (XLSX/WTMS) +
     &            (XLWX/WTMW))
                CALL SP_B( TX,XLSX,PSBX )
                PVAX = HCAW*XMLAX
                XMGAX = PVAX/(PVAX+PSBX)
                CALL EQUIL( XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,
     &            XMLAX,XMLSX,XMLWX )
                CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
                CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                CALL ENTH_B( TX,XLSX,HLWX,HBX )
                HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*HGAX
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + SRX(4)*XLSX
                SRCA_FRC(M,NTX) = SRCA_FRC(M,NTX) + SRX(4)*XLAX
                SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) + SRX(4)*XLWX
                SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) + SRX(4)*HLX
              ELSE
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) +
     &            SRX(4)*XLS_FRC(M,NTX)
                SRCA_FRC(M,NTX) = SRCA_FRC(M,NTX) +
     &            SRX(4)*XLA_FRC(M,NTX)
                SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) +
     &            SRX(4)*XLW_FRC(M,NTX)
                SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) +
     &            SRX(4)*HL_FRC(M,NTX)
              ENDIF
              IF( M.EQ.2 ) SRCL_FRC(NTX) = SRX(4)
!
!---        Pressure sink  ---
!
            ELSEIF( ISRT_FRC(NS).EQ.35 ) THEN
!
!---          Sink pressure and skin factor  ---
!
              PBHX = SRX(4) + PATM
              SBHX = SRX(2)
!
!---          Fracture triangle equivalent radius  ---
!
              RFEX = SQRT( AF_FRC(NTX)/GPI )
!
!---          Flow index  ---
!
              WIX = PERM_FRC(M,NTX)*APM_FRC(M,NTX)/(1.D+0+SBHX)
!
!---          Aqueous mass flux into borehole (only negative flux
!             allowed)  ---
!
              PLX = PL_FRC(M,NTX) + PATM
              FLX = MIN( (PBHX-PLX),0.D+0 )*WIX*RKL_FRC(M,NTX)*
     &          RHOL_FRC(M,NTX)/VISL_FRC(M,NTX)
              IF( M.EQ.2 ) SRCL_FRC(NTX) = FLX
!
!---          Gas mass flux into borehole (only negative flux
!             allowed)  ---
!
              PGX = PG_FRC(M,NTX) + PATM
              FGX = MIN( (PBHX-PGX),0.D+0 )*WIX*RKG_FRC(M,NTX)*
     &          RHOG_FRC(M,NTX)/VISG_FRC(M,NTX)
              IF( M.EQ.2 ) SRCG_FRC(NTX) = FGX
!
!---          Sources  ---
!
              SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + FLX*XLS_FRC(M,NTX)
              SRCA_FRC(M,NTX) = SRCA_FRC(M,NTX) + FLX*XLA_FRC(M,NTX) +
     &          FGX*XGA_FRC(M,NTX)
              SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) + FLX*XLW_FRC(M,NTX) +
     &          FGX*XGW_FRC(M,NTX)
              SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) + FLX*HL_FRC(M,NTX) +
     &          FGX*HG_FRC(M,NTX)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Thermocatalytic reactions with salt tracking ---
!
      IF( ITC_BH.GT.20 .AND. ITC_BH.LE.30 ) THEN
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
!---        Thermo-catalytic fluid, temperature-transition
!           equation via a first-order rate form, (J/kg aqu)/(m^3 aqu)
!           with salt tracking  ---
!
            IF( ITC_BH.EQ.22 ) THEN
              DO M = 2,ISVC+2
!
!---            Salt mass source equal to 1.e-4 times water mass in
!               borehole node  ---
!
                VARX = XLS_FRC(M,NTX)*APM_FRC(M,NTX)*AF_FRC(NTX)*
     &            SL_FRC(M,NTX)*RHOL_FRC(M,NTX)
                VAR1X = PTC_BH(2)
                VAR2X = 1.D+0/(1.D+0 + EXP((PTC_BH(3)-T_FRC(M,NTX))/
     &            PTC_BH(4)))
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) - VARX*VAR1X*VAR2X
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORC_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THA_BF_GT
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
!     Compute the contribution to the energy flux by aqueous and
!     gas advection for boreholes and fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 13 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVT_FRC
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
      SUB_LOG(ISUB_LOG) = '/THA_BF_GT'
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
!---      Loop over borehole node to node connections ---
!
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBN2 = IBCM_BH(NCX)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              HG1 = HG_BH(MP,NBN1)*RHOG_BH(MP,NBN1)
              HG2 = HG_BH(MN,NBN2)*RHOG_BH(MN,NBN2)
              HL1 = HL_BH(MP,NBN1)*RHOL_BH(MP,NBN1)
              HL2 = HL_BH(MN,NBN2)*RHOL_BH(MN,NBN2)
              UBBQ(M,NCX) = UBBQ(M,NCX)
     &          + HG1*MAX(UBBG(M,NCX),ZERO)
     &          - HG2*MAX(-UBBG(M,NCX),ZERO)
     &          + HL1*MAX(UBBL(M,NCX),ZERO)
     &          - HL2*MAX(-UBBL(M,NCX),ZERO)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
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
!---      Loop over fracture triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              HG1 = HG_FRC(MP,NT1X)*RHOG_FRC(MP,NT1X)
              HG2 = HG_FRC(MN,NT2X)*RHOG_FRC(MN,NT2X)
              HL1 = HL_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)
              HL2 = HL_FRC(MN,NT2X)*RHOL_FRC(MN,NT2X)
              UFFQ(M,NCX) = UFFQ(M,NCX)
     &          + HG1*MAX(UFFG(M,NCX),ZERO)
     &          - HG2*MAX(-UFFG(M,NCX),ZERO)
     &          + HL1*MAX(UFFL(M,NCX),ZERO)
     &          - HL2*MAX(-UFFL(M,NCX),ZERO)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THA_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THD_BF_GT
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
!     Compute the contribution to the energy flux by thermal conduction
!     for boreholes and fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 13 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVT_FRC
      USE FDVP_FRC
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
      REAL*8 NUL1X,NUG1X,NUL2X,NUG2X
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THD_BF_GT'
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
          INV1 = INV_BH(NBN1)
!
!---      Loop over borehole node to node connections ---
!
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBN2 = IBCM_BH(NCX)
            DBB1X = DBBM_BH(NCX)
            DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
            DTK = T_BH(2,NBN1)-T_BH(2,NBN2)
            INV2 = INV_BH(NBN2)
!
!---        Coaxial borehole  ---
!
            IF( IT_BH(1,NBH).GE.10000 ) THEN
!
!---          Both borehole nodes outer pipes  ---
!
              IF( IBN_BH(NBN1).NE.0 .AND. IBN_BH(NBN2).NE.0 ) THEN
!
!---            Annular area of pipe wall ---
!
                APW1X = GPI*(PAR_BH(4,INV1)**2-PAR_BH(3,INV1)**2)
!
!---            Fluid area of pipe s ---
!
                APF1X = GPI*((PAR_BH(3,INV1)**2)-(PAR_BH(2,INV1)**2))
                THK1X = PAR_BH(14,INV1)
!
!---            Annular area of pipe wall ---
!
                APW2X = GPI*(PAR_BH(4,INV2)**2-PAR_BH(3,INV2)**2)
!
!---            Fluid area of pipe s ---
!
                APF2X = GPI*((PAR_BH(3,INV2)**2)-(PAR_BH(2,INV2)**2))
                THK2X = PAR_BH(14,INV2)
!
!---          Borehole node 1 is outer pipe
!
              ELSEIF( IBN_BH(NBN1).NE.0 ) THEN
!
!---            Annular area of pipe wall ---
!
                APW1X = GPI*(PAR_BH(4,INV1)**2-PAR_BH(3,INV1)**2)
!
!---            Fluid area of pipe s ---
!
                APF1X = GPI*((PAR_BH(3,INV1)**2)-(PAR_BH(2,INV1)**2))
                THK1X = PAR_BH(14,INV1)
!
!---            Annular area of pipe wall ---
!
                APW2X = GPI*(PAR_BH(2,INV2)**2-PAR_BH(1,INV2)**2)
!
!---            Fluid area of pipe s ---
!
                APF2X = GPI*(PAR_BH(1,INV2)**2)
                THK2X = PAR_BH(5,INV2)
!
!---          Borehole node 2 is outer pipe
!
              ELSEIF( IBN_BH(NBN2).NE.0 ) THEN
!
!---            Annular area of pipe wall ---
!
                APW1X = GPI*(PAR_BH(2,INV1)**2-PAR_BH(1,INV1)**2)
!
!---            Fluid area of pipe s ---
!
                APF1X = GPI*(PAR_BH(1,INV1)**2)
                THK1X = PAR_BH(5,INV1)
!
!---            Annular area of pipe wall ---
!
                APW2X = GPI*(PAR_BH(4,INV2)**2-PAR_BH(3,INV2)**2)
!
!---            Fluid area of pipe s ---
!
                APF2X = GPI*((PAR_BH(3,INV2)**2)-(PAR_BH(2,INV2)**2))
                THK2X = PAR_BH(14,INV2)
!
!---          Both borehole node inner pipes  ---
!
              ELSE
!
!---            Annular area of pipe wall ---
!
                APW1X = GPI*(PAR_BH(2,INV1)**2-PAR_BH(1,INV1)**2)
!
!---            Fluid area of pipe s ---
!
                APF1X = GPI*(PAR_BH(1,INV1)**2)
                THK1X = PAR_BH(5,INV1)
!
!---            Annular area of pipe wall ---
!
                APW2X = GPI*(PAR_BH(2,INV2)**2-PAR_BH(1,INV2)**2)
!
!---            Fluid area of pipe s ---
!
                APF2X = GPI*(PAR_BH(1,INV2)**2)
                THK2X = PAR_BH(5,INV2)
              ENDIF
            ENDIF
!
!---        Pipe flow in interval of borehole node 1 ---
!
            IF( IS_BH(INV1).GT.10 .AND. IS_BH(INV1).LT.100 ) THEN
!
!---          Annular area of pipe wall ---
!
              APW1X = GPI*(PAR_BH(2,INV1)**2-PAR_BH(3,INV1)**2)
!
!---          Fluid area of pipe s ---
!
              APF1X = GPI*(PAR_BH(3,INV1)**2)
              THK1X = PAR_BH(4,INV1)
            ENDIF
!
!---        Pipe flow in interval of borehole node 2 ---
!
            IF( IS_BH(INV2).GT.10 .AND. IS_BH(INV2).LT.100 ) THEN
!
!---          Annular area of pipe wall ---
!
              APW2X = GPI*(PAR_BH(2,INV2)**2-PAR_BH(3,INV2)**2)
!
!---          Fluid area of pipe s ---
!
              APF2X = GPI*(PAR_BH(3,INV2)**2)
              THK2X = PAR_BH(4,INV2)
            ENDIF
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
!
!---          Pipe in interval of borehole node 1 ---
!
              IF( IS_BH(INV1).GT.10 ) THEN
                TK1 = APW1X*THK1X + APF1X*(THKL_BH(MP,NBN1)*
     &            SL_BH(MP,NBN1) + THKG_BH(MP,NBN1)*SG_BH(MP,NBN1))/
     &            (APW1X+APF1X)
!
!---          Porous media in interval of borehole node 1 ---
!
              ELSE
                IZ1X = IZ_BH(NBN1)
!
!---            Parallel function  ---
!
                IF( ITHK(IZ1X).EQ.2 ) THEN
                  TK1 = MAX(1.D+0-PORD_BH(MP,NBN1),0.D+0)*THKS(1,IZ1X) +
     &              PORD_BH(MP,NBN1)*(THKL_BH(MP,NBN1)*SL_BH(MP,NBN1) +
     &              THKG_BH(MP,NBN1)*SG_BH(MP,NBN1))
!
!---            Somerton function  ---
!
                ELSEIF( ITHK(IZ1X).EQ.4 ) THEN
                  TK1 = THKS(1,IZ1X) +
     &              SQRT(SL_BH(MP,NBN1))*(THKS(4,IZ1X)-THKS(1,IZ1X))
                ENDIF
              ENDIF
!
!---          Pipe in interval of borehole node 2 ---
!
              IF( IS_BH(INV1).GT.10 ) THEN
                TK2 = APW2X*THK2X + APF2X*(THKL_BH(MN,NBN2)*
     &            SL_BH(MN,NBN2) + THKG_BH(MN,NBN2)*SG_BH(MN,NBN2))/
     &            (APW2X+APF2X)
!
!---          Porous media in interval 2 ---
!
              ELSE
                IZ2X = IZ_BH(NBN2)
!
!---            Parallel function  ---
!
                IF( ITHK(IZ2X).EQ.2 ) THEN
                  TK2 = MAX(1.D+0-PORD_BH(MN,NBN2),0.D+0)*THKS(1,IZ2X) +
     &              PORD_BH(MN,NBN2)*(THKL_BH(MN,NBN2)*SL_BH(MN,NBN2) +
     &              THKG_BH(MN,NBN2)*SG_BH(MN,NBN2))
!
!---            Somerton function  ---
!
                ELSEIF( ITHK(IZ2X).EQ.4 ) THEN
                  TK2 = THKS(1,IZ2X) +
     &              SQRT(SL_BH(MN,NBN2))*(THKS(4,IZ2X)-THKS(1,IZ2X))
                ENDIF
              ENDIF
              INDX = 1
              TK = DIFMN( TK1,TK2,DBB1X,DBB2X,DTK,INDX )
              UBBQ(M,NCX) = TK*(T_BH(MP,NBN1)-T_BH(MN,NBN2))/
     &          DBB_BH(NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Coaxial boreholes connection between inner and outer
!     boreholes, loop over boreholes  ---
!
      DO NBH = 1,N_BH
!
!---    Skip for non-coaxial type boreholes  ---
!
        IF( IT_BH(1,NBH).LT.10000 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN1 = ID_BH(3,NBH),ID_BH(4,NBH)
          INV1 = INV_BH(NBN1)
!
!---      Inner borehole node to outer borehole node connections ---
!
          NCX = IPB_BH(3,NBN1)
          IF( NCX.EQ.0 ) CYCLE
!
!---      Last outer borehole node ---
!
          IF( NCX.EQ.-1 ) THEN
            NCX = IPB_BH(1,NBN1)
!
!---      Last inner borehole node ---
!
          ELSEIF( NCX.EQ.-2 ) THEN
            NCX = IPB_BH(2,NBN1)
          ENDIF
          NBN2 = IBCM_BH(NCX)
          INV2 = INV_BH(NBN2)
!
!---      Borehole node 1 is the outer pipe, hydraulic radius equals
!         square root of flow area divided by pi ---
!
          IF( IBN_BH(NBN1).NE.0 ) THEN
            RB1X = SQRT(PAR_BH(3,INV1)**2-PAR_BH(2,INV1)**2)
          ELSE
!
!---      Borehole node 1 is the inner pipe, hydraulic radius equals
!         inner radius of inner pipe ---
!
            RB1X = PAR_BH(1,INV1)
          ENDIF
!
!---      Borehole node 2 is the outer pipe, hydraulic radius equals
!         square root of flow area divided by pi ---
!
          IF( IBN_BH(NBN2).NE.0 ) THEN
            RB2X = SQRT(PAR_BH(3,INV2)**2-PAR_BH(2,INV2)**2)
!
!---      Borehole node 2 is the inner pipe, hydraulic radius equals
!         inner radius of inner pipe ---
!
          ELSE
            RB2X = PAR_BH(1,INV2)
          ENDIF
!
!---      Loop over flux increments ---
!
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Loop over borehole node to node connections ---
!
            UBL1X = 0.D+0
            UBG1X = 0.D+0
            NC = 0
            DO ICX = 1,2
              MCX = IPB_BH(ICX,NBN1)
              IF( MCX.LE.0 ) CYCLE
              UBL1X = UBL1X + ABS(UBBL(1,MCX))
              UBG1X = UBG1X + ABS(UBBG(1,MCX))
              NC = NC + 1
            ENDDO
!
!---        Aqueous and gas fluid velocity ---
!
            UBL1X = UBL1X/(REAL(NC) + SMALL)
            UBG1X = UBG1X/(REAL(NC) + SMALL)
!
!---        Hydraulic diameter  ---
!
            DB1X = 2.D+0*RB1X
!
!---        Aqueous and gas Reynolds number  ---
!
            REL1X = UBL1X*DB1X*RHOL_BH(MP,NBN1)/VISL_BH(MP,NBN1)
            REG1X = UBG1X*DB1X*RHOG_BH(MP,NBN1)/VISG_BH(MP,NBN1)
!
!---        Aqueous and gas Prandtl number ---
!
            CPL1X = (HL_BH(3,NBN1)-HL_BH(2,NBN1))/
     &        (T_BH(3,NBN1)-T_BH(2,NBN1))
            CPG1X = (HG_BH(3,NBN1)-HG_BH(2,NBN1))/
     &        (T_BH(3,NBN1)-T_BH(2,NBN1))
            PRL1X = CPL1X*VISL_BH(MP,NBN1)/THKL_BH(MP,NBN1)
            PRG1X = CPG1X*VISG_BH(MP,NBN1)/THKG_BH(MP,NBN1)
!
!---        Aqueous and gas Darcy friction factors ---
!
            FL1X = 1.D+0/((0.79D+0*LOG(REL1X)-1.64D+0)**2)
            FG1X = 1.D+0/((0.79D+0*LOG(REG1X)-1.64D+0)**2)
!
!---        Aqueous and gas Nusselt number ---
!
            NUL1X = 8.7D-5*(REL1X**0.92D+0)*(PRL1X**1.89D+0)
            NUG1X = 8.7D-5*(REG1X**0.92D+0)*(PRG1X**1.89D+0)
            NUL1X = (1.25D-1*FL1X)*(REL1X-1.D+3)*PRL1X/
     &        (1.D+0 + 1.27D+1*SQRT((1.25D-1*FL1X))*
     &        ((PRL1X**(2.D+0/3.D+0))-1.D+0))
            NUG1X = (1.25D-1*FG1X)*(REG1X-1.D+3)*PRG1X/
     &        (1.D+0 + 1.27D+1*SQRT((1.25D-1*FG1X))*
     &        ((PRG1X**(2.D+0/3.D+0))-1.D+0))
!
!---        Aqueous and gas heat transfer coefficient ---
!
            IF( DB1X.GE.EPSL ) THEN
              HCL1X = NUL1X*THKL_BH(MP,NBN1)/DB1X
              HCG1X = NUG1X*THKG_BH(MP,NBN1)/DB1X
            ELSE
              HCL1X = 0.D+0
              HCG1X = 0.D+0
            ENDIF
!
!---        Loop over borehole node to node connections ---
!
            UBL2X = 0.D+0
            UBG2X = 0.D+0
            NC = 0
            MF = MFLX(M)
            DO ICX = 1,2
              MCX = IPB_BH(ICX,NBN2)
              IF( MCX.LE.0 ) CYCLE
              UBL2X = UBL2X + ABS(UBBL(1,MCX))
              UBG2X = UBG2X + ABS(UBBG(1,MCX))
              NC = NC + 1
            ENDDO
!
!---        Aqueous and gas fluid velocity ---
!
            UBL2X = UBL2X/(REAL(NC) + SMALL)
            UBG2X = UBG2X/(REAL(NC) + SMALL)
!
!---        Hydraulic diameter  ---
!
            DB2X = 2.D+0*RB2X
!
!---        Aqueous and gas Reynolds number  ---
!
            REL2X = UBL2X*DB2X*RHOL_BH(MN,NBN2)/VISL_BH(MN,NBN2)
            REG2X = UBG2X*DB2X*RHOG_BH(MN,NBN2)/VISG_BH(MN,NBN2)
!
!---        Aqueous and gas Prandtl number ---
!
            CPL2X = (HL_BH(3,NBN2)-HL_BH(2,NBN2))/
     &        (T_BH(3,NBN2)-T_BH(2,NBN2))
            CPG2X = (HG_BH(3,NBN2)-HG_BH(2,NBN2))/
     &        (T_BH(3,NBN2)-T_BH(2,NBN2))
            PRL2X = CPL2X*VISL_BH(MN,NBN2)/THKL_BH(MN,NBN2)
            PRG2X = CPG2X*VISG_BH(MN,NBN2)/THKG_BH(MN,NBN2)
!
!---        Aqueous and gas Darcy friction factors ---
!
            FL2X = 1.D+0/((0.79D+0*LOG(REL2X)-1.64D+0)**2)
            FG2X = 1.D+0/((0.79D+0*LOG(REG2X)-1.64D+0)**2)
!
!---        Aqueous and gas Nusselt number ---
!
            NUL2X = 8.7D-5*(REL2X**0.92D+0)*(PRL2X**1.89D+0)
            NUG2X = 8.7D-5*(REG2X**0.92D+0)*(PRG2X**1.89D+0)
            NUL2X = (1.25D-1*FL2X)*(REL2X-1.D+3)*PRL2X/
     &        (1.D+0 + 1.27D+1*SQRT((1.25D-1*FL2X))*
     &        ((PRL2X**(2.D+0/3.D+0))-1.D+0))
            NUG2X = (1.25D-1*FG2X)*(REG2X-1.D+3)*PRG2X/
     &        (1.D+0 + 1.27D+1*SQRT((1.25D-1*FG2X))*
     &        ((PRG2X**(2.D+0/3.D+0))-1.D+0))
!
!---        Aqueous and gas heat transfer coefficient ---
!
            IF( DB2X.GE.EPSL ) THEN
              HCL2X = NUL2X*THKL_BH(MN,NBN2)/DB2X
              HCG2X = NUG2X*THKG_BH(MN,NBN2)/DB2X
            ELSE
              HCL2X = 0.D+0
              HCG2X = 0.D+0
            ENDIF
!
!---        Advective heat transfer coefficient, W/m^2 K, based
!           on area of outer diameter of inner pipe. Area conversion
!           from hydraulic diameter to outer diameter of inner pipe  ---
!
            HCA1X = SL_BH(MP,NBN1)*HCL1X + SG_BH(MP,NBN1)*HCG1X
            HCA1X = HCA1X*RB1X/PAR_BH(2,INV1)
            HCA2X = SL_BH(MN,NBN2)*HCL2X + SG_BH(MN,NBN2)*HCG2X
            HCA2X = HCA2X*RB2X/PAR_BH(2,INV1)
!
!---        Conduction across borehole wall for pipes, W/m^2 K, based
!           on area of outer diameter of inner pipe ---
!
            HCPX = 0.D+0
            IF( PAR_BH(2,INV1)-PAR_BH(1,INV1).GT.EPSL ) THEN
              HCPX = PAR_BH(5,INV1)/
     &          (PAR_BH(2,INV1)*LOG(PAR_BH(2,INV1)/PAR_BH(1,INV1)))
            ENDIF
!
!---        Overall heat transfer coefficient, W/m^2 K, based
!           on area of outer diameter of inner pipe ---
!
            HCOX = 0.D+0
            IF( HCA1X.GE.EPSL ) HCOX = HCOX + 1.D+0/HCA1X
            IF( HCA2X.GE.EPSL ) HCOX = HCOX + 1.D+0/HCA2X
            IF( HCPX.GE.EPSL ) HCOX = HCOX + 1.D+0/HCPX
            IF( HCOX.GE.EPSL ) HCOX = 1.D+0/HCOX
!
!---        Heat flux between inner and outer pipes, W/m^2, based
!           on area of outer diameter of inner pipe  ---
!
            UBBQC(M,NBN1) = HCOX*(T_BH(MP,NBN1)-T_BH(MN,NBN2))
          ENDDO
        ENDDO
      ENDDO
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
!---      Loop over fracture triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
            DTK = T_FRC(2,NT1X)-T_FRC(2,NT2X)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              TK1 = THKL_FRC(MP,NT1X)*SL_FRC(MP,NT1X) +
     &          THKG_FRC(MP,NT1X)*SG_FRC(MP,NT1X)
              TK2 = THKL_FRC(MN,NT2X)*SL_FRC(MN,NT2X) +
     &          THKG_FRC(MN,NT2X)*SG_FRC(MN,NT2X)
              INDX = 1
              TK = DIFMN( TK1,TK2,DFF1X,DFF2X,DTK,INDX )
              UFFQ(M,NCX) = TK*(T_FRC(MP,NT1X)-T_FRC(MN,NT2X))/
     &          DFF_FRC(NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THD_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDG_BF_GT
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
!     Compute the contribution to the energy flux energy flux by
!     water vapor and air molecular diffusion for fractures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 13 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_FRC
      USE FLUX_BH
      USE FDVT_FRC
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
      SUB_LOG(ISUB_LOG) = '/THDG_BF_GT'
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
!---      Loop over borehole node to node connections ---
!
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN1)
            IF( NCX.LE.0 ) CYCLE
            NBN2 = IBCM_BH(NCX)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              UDGAX = -UBBDGW(M,NCX)
              UBBQ(M,NCX) = UBBQ(M,NCX)
     &          + HGW_BH(MP,NBN1)*MAX(UBBDGW(M,NCX),ZERO)*WTMW
     &          - HGW_BH(MN,NBN2)*MAX(-UBBDGW(M,NCX),ZERO)*WTMW
     &          + HGA_BH(MP,NBN1)*MAX(UDGAX,ZERO)*WTMA
     &          - HGA_BH(MN,NBN2)*MAX(-UDGAX,ZERO)*WTMA
            ENDDO
          ENDDO
        ENDDO
      ENDDO
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
!---      Loop over fracture triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              UDGAX = -UFFDGW(M,NCX)
              UFFQ(M,NCX) = UFFQ(M,NCX)
     &          + HGW_FRC(MP,NT1X)*MAX(UFFDGW(M,NCX),ZERO)*WTMW
     &          - HGW_FRC(MN,NT2X)*MAX(-UFFDGW(M,NCX),ZERO)*WTMW
     &          + HGA_FRC(MP,NT1X)*MAX(UDGAX,ZERO)*WTMA
     &          - HGA_FRC(MN,NT2X)*MAX(-UDGAX,ZERO)*WTMA
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDG_BF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THFX_BTF_GT
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
!     Compute the contribution to the energy flux by thermal conduction,
!     by aqueous and  gas advection from, and water vapor and air
!     molecular diffusion for borehole nodes to fracture triangles
!     (borehole node is considered the upper node for flux indexing).
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 7 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GEOM_FRC
      USE GEOM_BH
      USE FLUX_BH
      USE FDVT_FRC
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
      SUB_LOG(ISUB_LOG) = '/THFX_BTF_GT'
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
          INV = INV_BH(NBN)
!
!---      Pipe flow in interval 1 ---
!
          IF( IS_BH(INV).GT.10 ) THEN
!
!---        Annular area of pipe wall ---
!
            APW1X = GPI*(PAR_BH(2,INV)**2-PAR_BH(3,INV)**2)
!
!---        Fluid area of pipe s ---
!
            APF1X = GPI*(PAR_BH(3,INV)**2)
          ENDIF
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
!---      Unincremented temperature difference between borehole
!         node and fracture triangle  ---
!
          DTK = T_BH(2,NBN)-T_FRC(2,NTX)
!
!---      Loop over flux increments ---
!
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Borehole thermal conduction, as pipe  ---
!
            IF( IS_BH(INV).GT.10 ) THEN
              TK1 = APW1X*PAR_BH(4,INV) + APF1X*(THKL_BH(MP,NBN)*
     &          SL_BH(MP,NBN) + THKG_BH(MP,NBN)*SG_BH(MP,NBN))/
     &          (APW1X+APF1X)
!
!---        Borehole thermal conduction, as porous media  ---
!
            ELSE
              IZ1X = IZ_BH(NBN)
!
!---          Parallel function  ---
!
              IF( ITHK(IZ1X).EQ.2 ) THEN
                TK1 = MAX(1.D+0-PORD_BH(MP,NBN),0.D+0)*THKS(1,IZ1X) +
     &            PORD_BH(MP,NBN)*(THKL_BH(MP,NBN)*SL_BH(MP,NBN) +
     &            THKG_BH(MP,NBN)*SG_BH(MP,NBN))
!
!---          Somerton function  ---
!
              ELSEIF( ITHK(IZ1X).EQ.4 ) THEN
                TK1 = THKS(1,IZ1X) +
     &            SQRT(SL_BH(MP,NBN))*(THKS(4,IZ1X)-THKS(1,IZ1X))
              ENDIF
            ENDIF
!
!---        Fracture thermal conduction  ---
!
            TK2 = THKL_FRC(MN,NTX)*SL_FRC(MN,NTX) +
     &        THKG_FRC(MN,NTX)*SG_FRC(MN,NTX)
            INDX = 1
            TK = DIFMN( TK1,TK2,DBF1X,DBF2X,DTK,INDX )
            UBFQ(M,NCX) = TK*(T_BH(MP,NBN)-T_FRC(MN,NTX))/DBFX
!
!---        Thermal advection  ---
!
            HG1 = HG_BH(MP,NBN)*RHOG_BH(MP,NBN)
            HG2 = HG_FRC(MN,NTX)*RHOG_FRC(MN,NTX)
            HL1 = HL_BH(MP,NBN)*RHOL_BH(MP,NBN)
            HL2 = HL_FRC(MN,NTX)*RHOL_FRC(MN,NTX)
            UBFQ(M,NCX) = UBFQ(M,NCX)
     &        + HG1*MAX(UBFG(M,NCX),ZERO)
     &        - HG2*MAX(-UBFG(M,NCX),ZERO)
     &        + HL1*MAX(UBFL(M,NCX),ZERO)
     &        - HL2*MAX(-UBFL(M,NCX),ZERO)
!
!---        Thermal gas diffusion  ---
!
            UDGAX = -UBFDGW(M,NCX)
            UBFQ(M,NCX) = UBFQ(M,NCX)
     &        + HGW_BH(MP,NBN)*MAX(UBFDGW(M,NCX),ZERO)*WTMW
     &        - HGW_FRC(MN,NTX)*MAX(-UBFDGW(M,NCX),ZERO)*WTMW
     &        + HGA_BH(MP,NBN)*MAX(UDGAX,ZERO)*WTMA
     &        - HGA_FRC(MN,NTX)*MAX(-UDGAX,ZERO)*WTMA
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THFX_BTF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRNS_BH_GT
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
!     Borehole to matrix transfer functions (borehole node is considered
!     the lower node for flux indexing).
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 27 August 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE PARM_BH
      USE JACOB
      USE GRID
      USE GEOM_BH
      USE FLUX_BH
      USE FDVT
      USE FDVS
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 NULX,NUGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TRNS_BH_GT'
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
!---      Coaxial borehole, only consider outer pipe  ---
!
          IF( IT_BH(1,NBH).GE.10000 .AND. IBN_BH(NBN).EQ.0 ) CYCLE
          IZB = IZ_BH(NBN)
          INV = INV_BH(NBN)
!
!---      Borehole node centroid and surface area  ---
!
          XP_BHX = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
          YP_BHX = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
          ZP_BHX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
!
!---      Outer pipe radius  ---
!
          INV = INV_BH(NBN)
          IF( IS_BH(INV).EQ.10000 ) THEN
            RBX = MAX( PAR_BH(4,INV),1.D-20 )
          ELSEIF( IS_BH(INV).GT.10 ) THEN
            RBX = MAX( PAR_BH(3,INV),1.D-20 )
          ELSE
            RBX = MAX( PAR_BH(2,INV),1.D-20 )
          ENDIF
!
!---      Outer coaxial pipe, hydraulic radius equals
!         square root of flow area divided by pi ---
!
          IF( IBN_BH(NBN).NE.0 ) THEN
            RB2X = SQRT(PAR_BH(3,INV)**2-PAR_BH(2,INV)**2)
!
!---      Inner pipe, hydraulic radius equals
!         inner radius of inner pipe ---
!
          ELSE
            RB2X = PAR_BH(1,INV)
          ENDIF
          AFB_BHX = 2.D+0*GPI*RBX*SQRT( (XP_BH(2,NBN)-XP_BH(1,NBN))**2 +
     &      (YP_BH(2,NBN)-YP_BH(1,NBN))**2 +
     &      (ZP_BH(2,NBN)-ZP_BH(1,NBN))**2 )
!
!---      Borehole node projections  ---
!
          XLX = ABS(XP_BH(1,NBN)-XP_BH(2,NBN))
          YLX = ABS(YP_BH(1,NBN)-YP_BH(2,NBN))
          ZLX = ABS(ZP_BH(1,NBN)-ZP_BH(2,NBN))
!
!---      Field node connected to borehole node  ---
!
          N = IBN_BH(NBN)
          IZN = IZ(N)
          PERMX = MAX( PERM(1,IZN),1.D-20 )
          PERMY = MAX( PERM(2,IZN),1.D-20 )
          PERMZ = MAX( PERM(3,IZN),1.D-20 )
          PERMYZ = SQRT(PERMY/PERMZ)
          PERMZY = SQRT(PERMZ/PERMY)
          DXGFX = DXGF(N)
          DYGFX = DYGF(N)
          DZGFX = DZGF(N)
          ROX = 2.8D-1*SQRT(PERMYZ*(DZGFX**2) + PERMZY*(DYGFX**2))
     &      /(SQRT(PERMYZ)+SQRT(PERMZY))
          PERMZX = SQRT(PERMZ/PERMX)
          PERMXZ = SQRT(PERMX/PERMZ)
          ROY = 2.8D-1*SQRT(PERMZX*(DXGFX**2) + PERMXZ*(DZGFX**2))
     &      /(SQRT(PERMZX)+SQRT(PERMXZ))
          PERMYX = SQRT(PERMY/PERMX)
          PERMXY = SQRT(PERMX/PERMY)
          ROZ = 2.8D-1*SQRT(PERMYX*(DXGFX**2) + PERMXY*(DYGFX**2))
     &      /(SQRT(PERMYX)+SQRT(PERMXY))
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            THKX = MAX( THKS(1,IZN),1.D-20 )
            THKY = MAX( THKS(2,IZN),1.D-20 )
            THKZ = MAX( THKS(3,IZN),1.D-20 )
            THKYZ = SQRT(THKY/THKZ)
            THKZY = SQRT(THKZ/THKY)
            ROTX = 2.8D-1*SQRT(THKYZ*(DZGFX**2) + THKZY*(DYGFX**2))
     &        /(SQRT(THKYZ)+SQRT(THKZY))
            THKZX = SQRT(THKZ/THKX)
            THKXZ = SQRT(THKX/THKZ)
            ROTY = 2.8D-1*SQRT(THKZX*(DXGFX**2) + THKXZ*(DZGFX**2))
     &        /(SQRT(THKZX)+SQRT(THKXZ))
            THKYX = SQRT(THKY/THKX)
            THKXY = SQRT(THKX/THKY)
            ROTZ = 2.8D-1*SQRT(THKYX*(DXGFX**2) + THKXY*(DYGFX**2))
     &        /(SQRT(THKYX)+SQRT(THKXY))
          ENDIF
!
!---      Loop over flux increments ---
!
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Cased borehole, unperforated pipe, or coaxial borehole  ---
!
            IF( MOD(IS_BH(INV),10).EQ.1 .OR. IS_BH(INV).EQ.10000 ) THEN
              QLX = 0.D+0
              QGX = 0.D+0
!
!---        Uncased borehole or perforated pipe  ---
!
            ELSE
!
!---          Well index components  ---
!
              PERMX = PERMRF(MP,N)*PERM(1,IZN)
              PERMY = PERMRF(MP,N)*PERM(2,IZN)
              PERMZ = PERMRF(MP,N)*PERM(3,IZN)
              WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/LOG(ROX/RBX)
              WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/LOG(ROY/RBX)
              WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/LOG(ROZ/RBX)
              WI_BHX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
!
!---          Matrix aqueous relative permeability  ---
!
              RKLX = RKL(1,MP,N)
!
!---          Aqueous head difference for borehole node and
!             grid cell, m  ---
!
              HDLX = PL_BH(MN,NBN) - PL(MP,N) +
     &          5.D-1*GRAV*(ZP_BHX-ZP(N))*(RHOL_BH(MN,NBN)+RHOL(MP,N))
!
!---          Aqueous relative permeability  ---
!
              INDX = 11
              RKLM = DIFMN( RKL_BH(MN,NBN),RKLX,DBN_BH(NBN),
     &          DBN_BH(NBN),HDLX,INDX )
!
!---          Aqueous viscosity Pa s  ---
!
              INDX = 5
              VLM = DIFMN( VISL_BH(MN,NBN),VISL(MP,N),DBN_BH(NBN),
     &          DBN_BH(NBN),HDLX,INDX )
!
!---          Aqueous volumetric rate borehole to matrix, m^3/s  ---
!
              IF( ISLC(38).EQ.0 ) THEN
                QLX = WI_BHX*RKLM*HDLX/VLM
              ELSE
                QLX = 0.D+0
              ENDIF
!
!---          Gas head difference for borehole node and field node, m  ---
!
              HDGX = PG_BH(MN,NBN) - PG(MP,N) +
     &          5.D-1*GRAV*(ZP_BHX-ZP(N))*(RHOG_BH(MN,NBN)+RHOG(MP,N))
!
!---          Gas relative permeability  ---
!
              INDX = 11
              RKGM = DIFMN( RKG_BH(MN,NBN),RKG(MP,N),DBN_BH(NBN),
     &          DBN_BH(NBN),HDGX,INDX )
!
!---          Gas viscosity Pa s  ---
!
              INDX = 5
              VGM = DIFMN( VISG_BH(MN,NBN),VISG(MP,N),DBN_BH(NBN),
     &          DBN_BH(NBN),HDGX,INDX )
!
!---          Gas volumetric rate borehole to matrix, m^3/s  ---
!
              IF( ISLC(38).EQ.0 ) THEN
                QGX = WI_BHX*RKGM*HDGX/VGM
              ELSE
                QGX = 0.D+0
              ENDIF
            ENDIF
!
!---        Nonisothermal simulations  ---
!
            IF( ISLC(30).EQ.0 ) THEN
!
!---          Thermal conductivity well index, parallel function  ---
!
              IF( ITHK(IZN).EQ.2 ) THEN
                THKX = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZN) +
     &            PORD(MP,N)*(THKL(MP,N)*SL(MP,N) +
     &            THKG(MP,N)*SG(MP,N))
                THKY = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZN) +
     &            PORD(MP,N)*(THKL(MP,N)*SL(MP,N) +
     &            THKG(MP,N)*SG(MP,N))
                THKZ = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZN) +
     &            PORD(MP,N)*(THKL(MP,N)*SL(MP,N) +
     &            THKG(MP,N)*SG(MP,N))
!
!---          Thermal conductivity well index, Somerton function  ---
!
              ELSEIF( ITHK(IZN).EQ.4 ) THEN
                THKX = THKS(1,IZN) +
     &            SQRT(SL(MP,N))*(THKS(4,IZN)-THKS(1,IZN))
                THKY = THKS(2,IZN) +
     &            SQRT(SL(MP,N))*(THKS(5,IZN)-THKS(2,IZN))
                THKZ = THKS(3,IZN) +
     &            SQRT(SL(MP,N))*(THKS(6,IZN)-THKS(3,IZN))
              ENDIF
              WIX = 2.D+0*GPI*SQRT(THKY*THKZ)*XLX/LOG(ROTX/RBX)
              WIY = 2.D+0*GPI*SQRT(THKX*THKZ)*YLX/LOG(ROTY/RBX)
              WIZ = 2.D+0*GPI*SQRT(THKX*THKY)*ZLX/LOG(ROTZ/RBX)
              WI_BHX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
!
!---          Heat transfer coefficient, Zhang et al. 2015. "The
!             analytical solution of the water-rock heat transfer
!             coefficient and sensitivity analyses of paramters."
!             Proceedings World Geothermal Congress 2015, Melbourne,
!             Australia, 19-25, April 2015  ---
!
              UBLX = 0.D+0
              UBGX = 0.D+0
!
!---          Loop over borehole node to node connections ---
!
              NC = 0
              DO ICX = 1,2
                MCX = IPB_BH(ICX,NBN)
                IF( MCX.LE.0 ) CYCLE
                UBLX = UBLX + ABS(UBBL(1,MCX))
                UBGX = UBGX + ABS(UBBG(1,MCX))
                NC = NC + 1
              ENDDO
!
!---          Aqueous and gas fluid velocity ---
!
              UBLX = UBLX/(REAL(NC) + SMALL)
              UBGX = UBGX/(REAL(NC) + SMALL)
!
!---          Aqueous and gas Reynolds number  ---
!
              DBX = 2.D+0*RB2X
              RELX = UBLX*DBX*RHOL_BH(MN,NBN)/VISL_BH(MN,NBN)
              REGX = UBGX*DBX*RHOG_BH(MN,NBN)/VISG_BH(MN,NBN)
!
!---          Aqueous and gas Prandtl number ---
!
              CPLX = (HL_BH(3,NBN)-HL_BH(2,NBN))/
     &          (T_BH(3,NBN)-T_BH(2,NBN))
              CPGX = (HG_BH(3,NBN)-HG_BH(2,NBN))/
     &          (T_BH(3,NBN)-T_BH(2,NBN))
              PRLX = CPLX*VISL_BH(MN,NBN)/THKL_BH(MN,NBN)
              PRGX = CPGX*VISG_BH(MN,NBN)/THKG_BH(MN,NBN)
!
!---          Aqueous and gas Darcy friction factors ---
!
              FLX = 1.D+0/((0.79D+0*LOG(RELX)-1.64D+0)**2)
              FGX = 1.D+0/((0.79D+0*LOG(REGX)-1.64D+0)**2)
!
!---          Aqueous and gas Nusselt number ---
!
              NULX = 8.7D-5*(RELX**0.92D+0)*(PRLX**1.89D+0)
              NUGX = 8.7D-5*(REGX**0.92D+0)*(PRGX**1.89D+0)
              NULX = (1.25D-1*FLX)*(RELX-1.D+3)*PRLX/
     &          (1.D+0 + 1.27D+1*SQRT((1.25D-1*FLX))*
     &          ((PRLX**(2.D+0/3.D+0))-1.D+0))
              NUGX = (1.25D-1*FGX)*(REGX-1.D+3)*PRGX/
     &          (1.D+0 + 1.27D+1*SQRT((1.25D-1*FGX))*
     &          ((PRGX**(2.D+0/3.D+0))-1.D+0))
!
!---          Aqueous and gas Nusselt number, assuming fully
!             developed laminar flow in a circular tube ---
!
!              NULX = 3.66D+0
!              NUGX = 3.66D+0
!
!---          Aqueous and gas heat transfer coefficient ---
!
              IF( DBX.GE.EPSL ) THEN
                HCLX = NULX*THKL_BH(MN,NBN)/DBX
                HCGX = NUGX*THKG_BH(MN,NBN)/DBX
              ELSE
                HCLX = 0.D+0
                HCGX = 0.D+0
              ENDIF
!
!---          Advective heat transfer coefficient, W/m^2 K ---
!
              HCAX = SL_BH(MN,NBN)*HCLX + SG_BH(MN,NBN)*HCGX
!
!---          Conduction across borehole wall for pipes ---
!
              IF( IS_BH(INV).EQ.10000 ) THEN
                HCAX = HCAX*RB2X/PAR_BH(4,INV)
                HCPX = 0.D+0
                IF( PAR_BH(4,INV)-PAR_BH(3,INV).GT.EPSL ) THEN
                  HCPX = PAR_BH(14,INV)/
     &              (PAR_BH(4,INV)*LOG(PAR_BH(4,INV)/PAR_BH(3,INV)))
                ENDIF
              ELSEIF( IS_BH(INV).GT.10 ) THEN
                HCAX = HCAX*RB2X/PAR_BH(3,INV)
                HCPX = 0.D+0
                IF( PAR_BH(2,INV)-PAR_BH(3,INV).GT.EPSL ) THEN
                  HCPX = PAR_BH(4,INV)/
     &              (PAR_BH(3,INV)*LOG(PAR_BH(2,INV)/PAR_BH(3,INV)))
                ENDIF
              ELSE
                HCAX = 0.D+0
                HCPX = 0.D+0
              ENDIF
!!
!!---          Conduction across borehole wall for pipes ---
!!
!              HCPX = 0.D+0
!              IF( IS_BH(INV).GT.10 ) THEN
!                AOX = GPI*(PAR_BH(2,INV)**2)
!                AIX = GPI*(PAR_BH(3,INV)**2)
!                HCAX = HCAX*AIX/AOX
!                IF( PAR_BH(2,INV)-PAR_BH(3,INV).GT.EPSL ) THEN
!                  HCPX = (PAR_BH(2,INV)-PAR_BH(3,INV))/(PAR_BH(4,INV)*
!     &              (AOX-AIX)/LOG(AOX/AIX))
!                  HCPX = 1.D+0/(HCPX+SMALL)
!                ENDIF
!              ENDIF
!
!---          Conduction heat transfer coefficient from
!             borehole wall to matrix node, W/m^2 K ---
!
              HCCX = WI_BHX/AFB_BHX
!
!---          Overall heat transfer coefficient, W/m^2 ---
!
              HCOX = 0.D+0
              IF( HCAX.GE.EPSL ) HCOX = HCOX + 1.D+0/HCAX
              IF( HCCX.GE.EPSL ) HCOX = HCOX + 1.D+0/HCCX
              IF( HCPX.GE.EPSL ) HCOX = HCOX + 1.D+0/HCPX
              IF( HCOX.GE.EPSL ) HCOX = 1.D+0/HCOX
!
!---          Overall heat transfer, W  ---
!
              QTCX = AFB_BHX*HCOX*(T_BH(MN,NBN)-T(MP,N))
            ENDIF
!
!---        Store aqueous, gas, and heat flow from fracture
!           to grid cell  ---
!
            IF( M.EQ.1 ) THEN
              UBML(NBN) = QLX
              UBMG(NBN) = QGX
              UBMT(NBN) = QTCX
            ENDIF
!
!---        Air mass flow rate from borehole node
!           to field node, kg/s  ---
!
            TRNSA_BH(M,NBN) =
     &        XLA_BH(MN,NBN)*RHOL_BH(MN,NBN)*MAX(QLX,ZERO) -
     &        XLA(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO) +
     &        XGA_BH(MN,NBN)*RHOG_BH(MN,NBN)*MAX(QGX,ZERO) -
     &        XGA(MP,N)*RHOG(MP,N)*MAX(-QGX,ZERO)
!
!---        Water mass flow rate from borehole node
!           to field node, kg/s  ---
!
            TRNSW_BH(M,NBN) =
     &        XLW_BH(MN,NBN)*RHOL_BH(MN,NBN)*MAX(QLX,ZERO) -
     &        XLW(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO) +
     &        XGW_BH(MN,NBN)*RHOG_BH(MN,NBN)*MAX(QGX,ZERO) -
     &        XGW(MP,N)*RHOG(MP,N)*MAX(-QGX,ZERO)
!
!---        Salt mass flow rate from borehole node
!           to field node, kg/s  ---
!
            TRNSS_BH(M,NBN) =
     &        XLS_BH(MN,NBN)*RHOL_BH(MN,NBN)*MAX(QLX,ZERO) -
     &        XLS(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO)
!
!---        Total heat transfer rate from borehole node
!           to field node, W  ---
!
!---        Nonisothermal simulations  ---
!
            IF( ISLC(30).EQ.0 ) THEN
              TRNSQ_BH(M,NBN) = QTCX +
     &          HL_BH(MN,NBN)*RHOL_BH(MN,NBN)*MAX(QLX,ZERO) -
     &          HL(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO) +
     &          HG_BH(MN,NBN)*RHOG_BH(MN,NBN)*MAX(QGX,ZERO) -
     &          HG(MP,N)*RHOG(MP,N)*MAX(-QGX,ZERO)
            ENDIF
!
!---      Loop over node increment indices  ---
!
          ENDDO
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
!---  End of TRNS_BH_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRNS_FRC_GT
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
!     Fracture to matrix transfer functions (fracture triangle is
!     considered the lower node for flux indexing).
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 27 August 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVT_FRC
      USE FDVT
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP
      USE FDVG_FRC
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
      REAL*8 NULX,NUGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TRNS_FRC_GT'
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
!
!---    Loop over fracture triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
          N = INCM_FRC(NCX)
          IZN = IZ(N)
!
!---      Unit surface normal magnitudes ---
!
          XLX = ABS(SFNT_FRC(1,NTX))
          YLX = ABS(SFNT_FRC(2,NTX))
          ZLX = ABS(SFNT_FRC(3,NTX))
!
!---      Loop over flux increments ---
!
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Matrix permeability, m^2  ---
!
            PERMX = SQRT((XLX*PERM(1,IZN))**2 + (YLX*PERM(2,IZN))**2 +
     &        (ZLX*PERM(3,IZN))**2)/SQRT(XLX**2 + YLX**2 + ZLX**2)
            PERMX = PERMX*PERMRF(MP,N)
!
!---        Matrix aqueous relative permeability  ---
!
            RKLX = RKL(1,MP,N)
!
!---        Aqueous head difference for fracture and grid cell, m  ---
!
            HDLX = PL_FRC(MN,NTX) - PL(MP,N) +
     &        5.D-1*GRAV*(ZP_FRC(NTX)-ZP(N))*
     &        (RHOL_FRC(MN,NTX)+RHOL(MP,N))
!
!---        Aqueous relative permeability  ---
!
            INDX = 8
            RKLM = DIFMN( RKL_FRC(MN,NTX),RKLX,DFN_FRC(NCX),
     &        DFN_FRC(NCX),HDLX,INDX )
!
!---        Aqueous viscosity Pa s  ---
!
            INDX = 5
            VLM = DIFMN( VISL_FRC(MN,NTX),VISL(MP,N),DFN_FRC(NCX),
     &        DFN_FRC(NCX),HDLX,INDX )
!
!---        Aqueous volumetric rate fracture to matrix, m^3/s  ---
!
            IF( ISLC(38).EQ.0 ) THEN
              QLX = 2.D+0*AFN_FRC(NCX)*PERMX*RKLM*HDLX/DFN_FRC(NCX)/VLM
            ELSE
              QLX = 0.D+0
            ENDIF
!
!---        Gas head difference for fracture and grid cell, m  ---
!
            HDGX = PG_FRC(MN,NTX) - PG(MP,N) +
     &        5.D-1*GRAV*(ZP_FRC(NTX)-ZP(N))*
     &        (RHOG_FRC(MN,NTX)+RHOG(MP,N))
!
!---        Gas relative permeability  ---
!
            INDX = 9
            RKGM = DIFMN( RKG_FRC(MN,NTX),RKG(MP,N),DFN_FRC(NCX),
     &        DFN_FRC(NCX),HDGX,INDX )
!
!---        Gas viscosity Pa s  ---
!
            INDX = 6
            VGM = DIFMN( VISG_FRC(MN,NTX),VISG(MP,N),DFN_FRC(NCX),
     &        DFN_FRC(NCX),HDGX,INDX )
!
!---        Gas volumetric rate fracture to matrix, m^3/s  ---
!
            IF( ISLC(38).EQ.0 ) THEN
              QGX = 2.D+0*AFN_FRC(NCX)*PERMX*RKGM*HDGX/DFN_FRC(NCX)/VGM
            ELSE
              QGX = 0.D+0
            ENDIF
!
!---        Nonisothermal simulations  ---
!
            IF( ISLC(30).EQ.0 ) THEN
!
!---          Thermal conductivity  ---
!
              THKX = SQRT((XLX*THKS(1,IZN))**2 + (YLX*THKS(2,IZN))**2 +
     &          (ZLX*THKS(3,IZN))**2)/SQRT(XLX**2 + YLX**2 + ZLX**2)
              THKX = THKX*MAX(1.D+0-PORD(MP,N),0.D+0) +
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + THKG(MP,N)*SG(MP,N))
!
!---          Heat transfer coefficient, Zhang et al. 2015. "The
!             analytical solution of the water-rock heat transfer
!             coefficient and sensitivity analyses of paramters."
!             Proceedings World Geothermal Congress 2015, Melbourne,
!             Australia, 19-25, April 2015  ---
!
              UFLX = 0.D+0
              UFGX = 0.D+0
!
!---          Loop over fracture triangle to triangle connections ---
!
              NC = 0
              DO NTTX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
                UFLX = UFLX + ABS(UFFL(1,NTTX))
                UFGX = UFGX + ABS(UFFG(1,NTTX))
                NC = NC + 1
              ENDDO
!
!---          Aqueous and gas fluid velocity ---
!
              UFLX = UFLX/(REAL(NC) + SMALL)
              UFGX = UFGX/(REAL(NC) + SMALL)
!
!---          Aqueous and gas Reynolds number ---
!
              RELX = UFLX*APH_FRC(MN,NTX)*RHOL_FRC(MN,NTX)/
     &          VISL_FRC(MN,NTX)
              REGX = UFGX*APH_FRC(MN,NTX)*RHOG_FRC(MN,NTX)/
     &          VISG_FRC(MN,NTX)
!
!---          Aqueous and gas Prandtl number ---
!
              CPLX = (HL_FRC(3,NTX)-HL_FRC(2,NTX))/
     &          (T_FRC(3,NTX)-T_FRC(2,NTX))
              CPGX = (HG_FRC(3,NTX)-HG_FRC(2,NTX))/
     &          (T_FRC(3,NTX)-T_FRC(2,NTX))
              PRLX = CPLX*VISL_FRC(MN,NTX)/THKL_FRC(MN,NTX)
              PRGX = CPGX*VISG_FRC(MN,NTX)/THKG_FRC(MN,NTX)
!
!---          Aqueous and gas Nusselt number ---
!
              NULX = 8.7D-5*(RELX**0.92D+0)*(PRLX**1.89D+0)
              NUGX = 8.7D-5*(RELX**0.92D+0)*(PRLX**1.89D+0)
!
!---          Aqueous and gas Nusselt number for laminar flow
!             in channels with infinite cross section ratios ---
!
              NULX = 7.54D+0
              NUGX = 7.54D+0
!
!---          Aqueous and gas heat transfer coefficient ---
!
              IF( APH_FRC(MN,NTX).GE.EPSL ) THEN
                HCLX = NULX*THKL_FRC(MN,NTX)/(2.D+0*APH_FRC(MN,NTX))
                HCGX = NUGX*THKG_FRC(MN,NTX)/(2.D+0*APH_FRC(MN,NTX))
              ELSE
                HCLX = 0.D+0
                HCGX = 0.D+0
              ENDIF
!
!---          Advective heat transfer coefficient, W/m^2 ---
!
              HCAX = SL_FRC(MN,NTX)*HCLX + SG_FRC(MN,NTX)*HCGX
!
!---          Conduction heat transfer coefficient, W/m^2 ---
!
              HCCX = THKX/DFN_FRC(NCX)
!
!---          Overall heat transfer coefficient, W/m^2 ---
!
              IF( HCAX.GE.EPSL .AND. HCCX.GE.EPSL ) THEN
                HCOX = 1.D+0/((1.D+0/HCAX)+(1.D+0/HCCX))
              ELSEIF( HCAX.GE.EPSL ) THEN
                HCOX = HCAX
              ELSEIF( HCCX.GE.EPSL ) THEN
                HCOX = HCCX
              ELSE
                HCOX = 0.D+0
              ENDIF
!
!---          Overall heat transfer, W  ---
!
              QTCX = 2.D+0*AFN_FRC(NCX)*HCOX*(T_FRC(MN,NTX)-T(MP,N))
            ENDIF
!
!---        Store aqueous, gas, and heat flow from fracture
!           to grid cell  ---
!
            IF( M.EQ.1 ) THEN
              UFMG(NCX) = QGX
              UFML(NCX) = QLX
              UFMT(NCX) = QTCX
            ENDIF
!
!---        Air mass flow rate from fracture to grid cell, kg/s  ---
!
            TRNSA_FRC(M,NCX) =
     &        XLA_FRC(MN,NTX)*RHOL_FRC(MN,NTX)*MAX(QLX,ZERO) -
     &        XLA(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO) +
     &        XGA_FRC(MN,NTX)*RHOG_FRC(MN,NTX)*MAX(QGX,ZERO) -
     &        XGA(MP,N)*RHOG(MP,N)*MAX(-QGX,ZERO)
!
!---        Water mass flow rate from fracture to grid cell, kg/s  ---
!
            TRNSW_FRC(M,NCX) =
     &        XLW_FRC(MN,NTX)*RHOL_FRC(MN,NTX)*MAX(QLX,ZERO) -
     &        XLW(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO) +
     &        XGW_FRC(MN,NTX)*RHOG_FRC(MN,NTX)*MAX(QGX,ZERO) -
     &        XGW(MP,N)*RHOG(MP,N)*MAX(-QGX,ZERO)
!
!---        Salt mass flow rate from fracture to grid cell, kg/s  ---
!
            TRNSS_FRC(M,NCX) =
     &        XLS_FRC(MN,NTX)*RHOL_FRC(MN,NTX)*MAX(QLX,ZERO) -
     &        XLS(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO)
!
!---        Total heat transfer rate from fracture to grid cell, W  ---
!
            IF( ISLC(30).EQ.0 ) THEN
              TRNSQ_FRC(M,NCX) = QTCX +
     &          HL_FRC(MN,NTX)*RHOL_FRC(MN,NTX)*MAX(QLX,ZERO) -
     &          HL(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO) +
     &          HG_FRC(MN,NTX)*RHOG_FRC(MN,NTX)*MAX(QGX,ZERO) -
     &          HG(MP,N)*RHOG(MP,N)*MAX(-QGX,ZERO)
            ENDIF
!
!---      Loop over node increment indices  ---
!
          ENDDO
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
!---  End of TRNS_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_FRC_GT
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
!     Update the fracture primary variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 18 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE PARM_BH
      USE PARM_FRC
      USE OUTPU
      USE JACOB
      USE HYST
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FILES
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP_BH
      USE FDVP
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 PH_CND(5)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Saturated w/o Entrapped Gas',
     &   'Unsaturated w/ or w/o Entrapped Gas',
     &   'Saturated w/ Trapped Gas',
     &   'Fully Unsaturated',
     &   'Supercritical Water'/
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_FRC_GT'
      IERR = 0
!
!---  Loop over fractures  ---
!
      DO NFX = 1,NF_FRC
        IF( IERR.NE.0 ) CYCLE
!
!---    Loop over fracture triangles  ---
!
        DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---      Skip inactive triangles  ---
!
          IF( IXP_FRC(NTX).EQ.0 ) CYCLE
          IF( IERR.NE.0 ) CYCLE
          MPL = IM_FRC(IEQW,NTX)
          IF( ISLC(37).EQ.0 ) THEN
            MPG = IM_FRC(IEQA,NTX)
            DPG = BLU(MPG)
          ENDIF
          DPL = BLU(MPL)
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            MPS = IM_FRC(IEQS,NTX)
            DPS = BLU(MPS)
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            MPT = IM_FRC(IEQT,NTX)
            DPT = BLU(MPT)
            DPT = SIGN( MIN( 1.D+0,ABS(DPT) ),DPT )
            T_FRC(2,NTX) = T_FRC(2,NTX) + DPT
          ENDIF
!
!---      Saturated system w/o entrapped gas
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - aqueous-air mole fraction
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ_FRC(2,NTX).EQ.1 ) THEN
            INDX = 0
            CALL REGION_4( T_FRC(2,NTX),PSWX,INDX )
            DPX = 1.D+6
            DPL = SIGN( MIN( ABS(DPL),DPX ),DPL )
            PL_FRC(2,NTX) = PL_FRC(2,NTX) + DPL
!
!---        Zero negative corrections for zero aqueous air  ---
!
            IF( ISLC(37).EQ.0 ) THEN
!              IF( XMLA_FRC(2,NTX)/EPSL.LT.EPSL .AND.
!     &          BLU(MPG)/EPSL.LT.EPSL ) THEN
!                BLU(MPG) = 0.D+0
!                DPG = 0.D+0
!              ENDIF
!              XMLA_FRC(2,NTX) = XMLA_FRC(2,NTX) + DPG
!              IF( XMLA_FRC(2,NTX).LT.1.D-12 ) XMLA_FRC(2,NTX) = 0.D+0
              PVA_FRC(2,NTX) = MAX( PVA_FRC(2,NTX)+DPG,0.D+0 )
              IF( PVA_FRC(2,NTX).LT.EPSL ) PVA_FRC(2,NTX) = 0.D+0
            ENDIF
!
!---        Isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Limit salt mass fraction changes to 0.25 of the
!             maximum value if salt mass fraction is less than
!             the maximum   ---
!
              PWX = MAX( PSWX,PL_FRC(2,NTX)+PATM )
              CALL SOL_BRNS( T_FRC(2,NTX),PWX,XLSMX )
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  ) THEN
                DPX = ABS(2.5D-1*PTC_BH(5))
              ELSE
                DPX = ABS(2.5D-1*XLSMX)
              ENDIF
              IF( YLS_FRC(2,NTX).LT.XLSMX ) THEN
                DPS = SIGN( MIN( DPX,ABS(DPS) ),DPS )
              ENDIF
!
!---          Zero negative corrections for zero aqueous salt  ---
!
              IF( YLS_FRC(2,NTX)/EPSL.LT.EPSL .AND.
     &          BLU(MPS)/EPSL.LT.EPSL ) THEN
                BLU(MPS) = 0.D+0
                DPS = 0.D+0
              ENDIF
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  ) THEN
                IF( YLS_FRC(2,NTX)+DPS.GT.PTC_BH(5) ) DPS = 6.D-1*DPS
                YLS_FRC(2,NTX) = YLS_FRC(2,NTX) + DPS
                IF( YLS_FRC(2,NTX).LT.EPSL ) YLS_FRC(2,NTX) = 0.D+0
                IF( YLS_FRC(2,NTX).GT.PTC_BH(5) )
     &            YLS_FRC(2,NTX) = PTC_BH(5)
              ELSE
                IF( YLS_FRC(2,NTX)+DPS.LT.0.D+0 ) DPS = 6.D-1*DPS
                YLS_FRC(2,NTX) = YLS_FRC(2,NTX) + DPS
                IF( YLS_FRC(2,NTX).LT.EPSL ) YLS_FRC(2,NTX) = 0.D+0
              ENDIF
              XLS_FRC(2,NTX) = MIN( YLS_FRC(2,NTX),XLSMX )
            ENDIF
!
!---      Unsaturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ_FRC(2,NTX).EQ.2 ) THEN
            INDX = 0
            CALL REGION_4( T_FRC(2,NTX),PSWX,INDX )
!
!---        Isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Limit salt mass fraction changes to 0.25 of the
!             maximum value if salt mass fraction is less than
!             the maximum   ---
!
              CALL SOL_BRNS( T_FRC(2,NTX),PSWX,XLSMX )
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  ) THEN
                DPX = ABS(2.5D-1*PTC_BH(5))
              ELSE
                DPX = ABS(2.5D-1*XLSMX)
              ENDIF
              IF( YLS_FRC(2,NTX).LT.XLSMX ) THEN
                DPS = SIGN( MIN( DPX,ABS(DPS) ),DPS )
              ENDIF
!
!---          Zero negative corrections for zero aqueous salt  ---
!
              IF( YLS_FRC(2,NTX)/EPSL.LT.EPSL .AND.
     &          BLU(MPS)/EPSL.LT.EPSL ) THEN
                BLU(MPS) = 0.D+0
                DPS = 0.D+0
              ENDIF
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  ) THEN
                IF( YLS_FRC(2,NTX)+DPS.GT.PTC_BH(5) ) DPS = 6.D-1*DPS
                YLS_FRC(2,NTX) = YLS_FRC(2,NTX) + DPS
                IF( YLS_FRC(2,NTX).LT.EPSL ) YLS_FRC(2,NTX) = 0.D+0
                IF( YLS_FRC(2,NTX).GT.PTC_BH(5) )
     &            YLS_FRC(2,NTX) = PTC_BH(5)
              ELSE
                IF( YLS_FRC(2,NTX)+DPS.LT.0.D+0 ) DPS = 6.D-1*DPS
                YLS_FRC(2,NTX) = YLS_FRC(2,NTX) + DPS
                IF( YLS_FRC(2,NTX).LT.EPSL ) YLS_FRC(2,NTX) = 0.D+0
              ENDIF
              XLS_FRC(2,NTX) = MIN( YLS_FRC(2,NTX),XLSMX )
            ENDIF
!
!---        Assign gas-entry pressure  ---
!
            ENPR = RKSP_FRC(1,NFX)*RHORL*GRAV
!
!---        Limit changes in pressure  ---
!
            DPX = MAX( 1.D+6,2.5D-1*(PG_FRC(2,NTX)-PL_FRC(2,NTX)) )
            DPCX = MAX( 1.D+4,2.5D-1*(PG_FRC(2,NTX)-PL_FRC(2,NTX)) )
            IF( (DPG-DPL).GT.DPCX ) THEN
              DPG = DPG - 5.D-1*DPCX
              DPL = DPL + 5.D-1*DPCX
            ELSE
              DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
              DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
            ENDIF
!
!---        Relax pressure updates when transitioning to unsaturated
!           conditions  ---
!
            IF( (PG_FRC(2,NTX)+DPG)-(PL_FRC(2,NTX)+DPL).LT.ENPR ) THEN
              DPG = 6.D-1*DPG
              DPL = 6.D-1*DPL
            ENDIF
            PG_FRC(2,NTX) = PG_FRC(2,NTX) + DPG
            PL_FRC(2,NTX) = PL_FRC(2,NTX) + DPL
            PL_FRC(2,NTX) = MAX( PL_FRC(2,NTX),
     &        (PG_FRC(2,NTX)-RKSP_FRC(5,NFX)*RHORL*GRAV) )
!
!---        Maintain the gas pressure above or at the water vapor
!           pressure  ---
!
            CALL P_IAPWS( T_FRC(2,NTX),PSWX,RHOGWX,RHOLWX,HGWX,
     &        HLWX,UGWX,ULWX )
            CALL DENS_B( XLS_FRC(2,NTX),RHOBX,RHOLWX,T_FRC(2,NTX) )
            CALL SP_B( T_FRC(2,NTX),XLS_FRC(2,NTX),PSBX )
            PGX = PG_FRC(2,NTX) + PATM
            PLX = PL_FRC(2,NTX) + PATM
            PCX = PGX - PLX
            CALL VPL_B( T_FRC(2,NTX),PSBX,PCX,RHOBX,PVBX )
            PG_FRC(2,NTX) = MAX( PG_FRC(2,NTX),(PVBX-PATM) )
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - gas pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ_FRC(2,NTX).EQ.4 ) THEN
!
!---        Limit changes in water vapor partial pressure  ---
!
            INDX = 0
            CALL REGION_4( T_FRC(2,NTX),PWX,INDX )
            PWX = MIN( PWX,PCRW )
            DPX = 2.5D-1*PWX
            DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
            PVW_FRC(2,NTX) = MAX( PVW_FRC(2,NTX)+DPL,1.D-6 )
!
!---        Limit changes in gas pressure  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              DPX = 2.5D-1*MAX(PG_FRC(2,NTX)+PATM,1.D+6)
              DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
              PG_FRC(2,NTX) = PG_FRC(2,NTX) + DPG
            ENDIF
!
!---        Isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Zero negative corrections for salt volumetric conc.  ---
!
              IF( TMS_FRC(2,NTX)/EPSL.LT.EPSL .AND.
     &          BLU(MPS)/EPSL.LT.EPSL ) THEN
                BLU(MPS) = 0.D+0
                DPS = 0.D+0
              ENDIF
              TMS_FRC(2,NTX) = TMS_FRC(2,NTX) + DPS
              IF( TMS_FRC(2,NTX).LT.1.D-9 ) TMS_FRC(2,NTX) = 0.D+0
            ENDIF
          ENDIF
!
!---      Check for excessive pressure or temperature   ---
!
          PGX = PG_FRC(2,NTX)+PATM
          PLX = PL_FRC(2,NTX)+PATM
          TKX = T_FRC(2,NTX)+TABS
          IF( PGX.GT.100.D+6 .OR. PGX.LT.0.D+0 .AND. IERR.EQ.0 ) 
     &      IERR = 1
          IF( PLX.GT.100.D+6 .AND. IERR.EQ.0 ) IERR = 1
          IF( TKX.GT.1073.15D+0 .OR. TKX.LT.TABS .AND. IERR.EQ.0 ) 
     &      IERR = 1
          IF( IERR.EQ.1 ) NSD_FRC(1) = NTX
        ENDDO
      ENDDO
!
!---  Reduce time step for excessive changes in primary variables   ---
!
      IF( IERR.EQ.1 ) THEN
        ICNV = 1
        NTX = NSD_FRC(1)
        WRITE(ISC,'(10X,A)') '---  Excessive Primary Variable Change '
     &    // 'for Fracture Flow ---'
        WRITE(IWR,'(10X,A)') '---  Excessive Primary Variable Change '
     &    // 'for Fracture Flow ---'
        WRITE(ISC,'(4X,A,I6)') 'Fracture Triangle = ',NTX
        WRITE(IWR,'(4X,A,I6)') 'Fracture Triangle = ',NTX
        WRITE(ISC,'(4X,2A)') 'Phase Condition = ',
     &    PH_CND(NPHAZ_FRC(2,NTX))
        WRITE(IWR,'(4X,2A)') 'Phase Condition = ',
     &    PH_CND(NPHAZ_FRC(2,NTX))
        WRITE(ISC,'(4X,A,1PE12.5)') 'Temperature = ',
     &    T_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Temperature = ',
     &    T_FRC(2,NTX)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',
     &    PL_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',
     &    PL_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Pressure = ',
     &    PG_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Pressure = ',
     &    PG_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG_FRC(2,NTX)
        WRITE(ISC,'(4X,A,1PE12.5)')
     &    'Aqueous-Air Mole Fraction = ',XMLA_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)')
     &    'Aqueous-Air Mole Fraction = ',XMLA_FRC(2,NTX)
        WRITE(ISC,'(4X,A,1PE12.5)')
     &    'Water Vapor Pressure = ',PVW_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)')
     &    'Water Vapor Pressure = ',PVW_FRC(2,NTX)
        WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &    'Total-Salt Aqu. Mass Fraction = ',YLS_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &    'Total-Salt Aqu. Mass Fraction = ',YLS_FRC(2,NTX)
        WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &    'Salt Volumetric Concentration = ',TMS_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &    'Salt Volumetric Concentration = ',TMS_FRC(2,NTX)
      ENDIF
!
!---  Reduce time step  ---
!
  300 CONTINUE
      IF( ICNV.EQ.1 ) THEN
        IF( NTSR.LT.4 .OR. (DTCF*DT).GT.DTMN ) THEN
          NTSR = NTSR + 1
          DTX = DT
          TM = TM - (1.D+0-DTCF)*DT
          DT = DTCF*DT
          DTO = DT
          DTI = 1.D+0/DT
          VAR = DT
          VARX = DTX
          IF( UNTM.NE.'null' ) THEN
            INDX = 1
            IUNS = 1
            CALL RDUNIT(UNTM,VAR,INDX)
            IUNS = 1
            CALL RDUNIT(UNTM,VARX,INDX)
            NCH = INDEX( UNTM,'  ')-1
          ENDIF
          WRITE(ISC,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
          WRITE(IWR,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
!
!---      Loop over fractures  ---
!
          DO NFX = 1,NF_FRC
!
!---        Loop over fracture triangles  ---
!
            DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---          Skip inactive triangles  ---
!
              IF( IXP_FRC(NTX).EQ.0 ) CYCLE
              T_FRC(2,NTX) = T_FRC(1,NTX)
              PL_FRC(2,NTX) = PL_FRC(1,NTX)
              PG_FRC(2,NTX) = PG_FRC(1,NTX)
              PVW_FRC(2,NTX) = PVW_FRC(1,NTX)
              PVA_FRC(2,NTX) = PVA_FRC(1,NTX)
              XMLA_FRC(2,NTX) = XMLA_FRC(1,NTX)
              SL_FRC(2,NTX) = SL_FRC(1,NTX)
              SG_FRC(2,NTX) = SG_FRC(1,NTX)
              YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
              TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
              NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
            ENDDO
          ENDDO
!
!---      Borehole flow  ---
!
          IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) THEN
!
!---        Loop over boreholes  ---
!
            DO NBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
                T_BH(2,NBN) = T_BH(1,NBN)
                PL_BH(2,NBN) = PL_BH(1,NBN)
                PG_BH(2,NBN) = PG_BH(1,NBN)
                PVW_BH(2,NBN) = PVW_BH(1,NBN)
                PVA_BH(2,NBN) = PVA_BH(1,NBN)
                XMLA_BH(2,NBN) = XMLA_BH(1,NBN)
                SL_BH(2,NBN) = SL_BH(1,NBN)
                SG_BH(2,NBN) = SG_BH(1,NBN)
                YLS_BH(2,NBN) = YLS_BH(1,NBN)
                TMS_BH(2,NBN) = TMS_BH(1,NBN)
                NPHAZ_BH(2,NBN) = NPHAZ_BH(1,NBN)
              ENDDO
            ENDDO
          ENDIF
!
!---      Matrix flow  ---
!
          DO N = 1,NFLD
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PVW(2,N) = PVW(1,N)
            PVA(2,N) = PVA(1,N)
            XMLA(2,N) = XMLA(1,N)
            SL(2,N) = SL(1,N)
            SG(2,N) = SG(1,N)
            SGT(2,N) = SGT(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            ASLMIN(2,N) = ASLMIN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
          ENDDO
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
!
!---      Loop over fractures  ---
!
          DO NFX = 1,NF_FRC
!
!---        Loop over fracture triangles  ---
!
            DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---          Skip inactive triangles  ---
!
              IF( IXP_FRC(NTX).EQ.0 ) CYCLE
              T_FRC(2,NTX) = T_FRC(1,NTX)
              PL_FRC(2,NTX) = PL_FRC(1,NTX)
              PG_FRC(2,NTX) = PG_FRC(1,NTX)
              PVW_FRC(2,NTX) = PVW_FRC(1,NTX)
              PVA_FRC(2,NTX) = PVA_FRC(1,NTX)
              XMLA_FRC(2,NTX) = XMLA_FRC(1,NTX)
              SL_FRC(2,NTX) = SL_FRC(1,NTX)
              SG_FRC(2,NTX) = SG_FRC(1,NTX)
              YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
              TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
              NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
            ENDDO
          ENDDO
!
!---      Borehole flow  ---
!
          IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) THEN
!
!---        Loop over boreholes  ---
!
            DO NBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
                T_BH(2,NBN) = T_BH(1,NBN)
                PL_BH(2,NBN) = PL_BH(1,NBN)
                PG_BH(2,NBN) = PG_BH(1,NBN)
                PVW_BH(2,NBN) = PVW_BH(1,NBN)
                PVA_BH(2,NBN) = PVA_BH(1,NBN)
                XMLA_BH(2,NBN) = XMLA_BH(1,NBN)
                SL_BH(2,NBN) = SL_BH(1,NBN)
                SG_BH(2,NBN) = SG_BH(1,NBN)
                YLS_BH(2,NBN) = YLS_BH(1,NBN)
                TMS_BH(2,NBN) = TMS_BH(1,NBN)
                NPHAZ_BH(2,NBN) = NPHAZ_BH(1,NBN)
              ENDDO
            ENDDO
          ENDIF
!
!---      Matrix flow  ---
!
          DO N = 1,NFLD
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PVW(2,N) = PVW(1,N)
            PVA(2,N) = PVA(1,N)
            XMLA(2,N) = XMLA(1,N)
            SL(2,N) = SL(1,N)
            SG(2,N) = SG(1,N)
            SGT(2,N) = SGT(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            ASLMIN(2,N) = ASLMIN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
          ENDDO
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     &from Fracture Flow  ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     &from Fracture Flow  ---'
          ICNV = 4
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_BH_GT
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
!     Update the borehole primary variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 2 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE PARM_BH
      USE OUTPU
      USE JACOB
      USE HYST
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FILES
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP_BH
      USE FDVP
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 PH_CND(5)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Saturated w/o Entrapped Gas',
     &   'Unsaturated w/ or w/o Entrapped Gas',
     &   'Saturated w/ Trapped Gas',
     &   'Fully Unsaturated',
     &   'Supercritical Water'/
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_BH_GT'
      IERR = 0
!
!---  Loop over boreholes  ---
!
      DO NBH = 1,N_BH
        IF( IERR.NE.0 ) CYCLE
!
!---    Skip for boundary condition type boreholes  ---
!
        IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---    Loop over borehole nodes  ---
!
        DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
          IF( IERR.NE.0 ) CYCLE
          IZN = IZ_BH(NBN)
          INV = INV_BH(NBN)
          MPL = IM_BH(IEQW,NBN)
          IF( ISLC(37).EQ.0 ) THEN
            MPG = IM_BH(IEQA,NBN)
            DPG = BLU(MPG)
          ELSE
            DPG = 0.D+0
          ENDIF
          DPL = BLU(MPL)
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            MPS = IM_BH(IEQS,NBN)
            DPS = BLU(MPS)
          ELSE
            DPS = 0.D+0
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            MPT = IM_BH(IEQT,NBN)
            DPT = BLU(MPT)
            DPT = SIGN( MIN( 1.D+0,ABS(DPT) ),DPT )
            T_BH(2,NBN) = T_BH(2,NBN) + DPT
          ELSE
            DPT = 0.D+0
          ENDIF
!
!---      Saturated system w/o entrapped gas
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - aqueous-air mole fraction
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ_BH(2,NBN).EQ.1 ) THEN
            INDX = 0
            CALL REGION_4( T_BH(2,NBN),PSWX,INDX )
!
!---        Coaxial borehole  ---
!
            IF( IT_BH(1,NBH).GE.10000 ) THEN
              DPX = 2.5D-2*(PL_BH(2,NBN)+PATM)
            ELSE
              DPX = 1.D+6
            ENDIF
            DPL = SIGN( MIN( ABS(DPL),DPX ),DPL )
            PL_BH(2,NBN) = PL_BH(2,NBN) + DPL
!
!---        Zero negative corrections for zero aqueous air  ---
!
            IF( ISLC(37).EQ.0 ) THEN
!              IF( XMLA_BH(2,NBN)/EPSL.LT.EPSL .AND.
!     &          BLU(MPG)/EPSL.LT.EPSL ) THEN
!                BLU(MPG) = 0.D+0
!                DPG = 0.D+0
!              ENDIF
!              XMLA_BH(2,NBN) = XMLA_BH(2,NBN) + DPG
!              IF( XMLA_BH(2,NBN).LT.1.D-12 ) XMLA_BH(2,NBN) = 0.D+0
              PVA_BH(2,NBN) = MAX( PVA_BH(2,NBN)+DPG,0.D+0 )
              IF( PVA_BH(2,NBN).LT.EPSL ) PVA_BH(2,NBN) = 0.D+0
            ENDIF
!
!---        Isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Limit salt mass fraction changes to 0.25 of the
!             maximum value if salt mass fraction is less than
!             the maximum   ---
!
              PWX = MAX( PSWX,PL_BH(2,NBN)+PATM )
              CALL SOL_BRNS( T_BH(2,NBN),PWX,XLSMX )
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  ) THEN
                DPX = ABS(2.5D-1*PTC_BH(5))
              ELSE
                DPX = ABS(2.5D-1*XLSMX)
              ENDIF
              IF( YLS_BH(2,NBN).LT.XLSMX ) THEN
                DPS = SIGN( MIN( DPX,ABS(DPS) ),DPS )
              ENDIF
!
!---          Zero negative corrections for zero aqueous salt  ---
!
              IF( YLS_BH(2,NBN)/EPSL.LT.EPSL .AND.
     &          BLU(MPS)/EPSL.LT.EPSL ) THEN
                BLU(MPS) = 0.D+0
                DPS = 0.D+0
              ENDIF
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  ) THEN
                IF( YLS_BH(2,NBN)+DPS.GT.PTC_BH(5) ) DPS = 6.D-1*DPS
                YLS_BH(2,NBN) = YLS_BH(2,NBN) + DPS
                IF( YLS_BH(2,NBN).LT.EPSL ) YLS_BH(2,NBN) = 0.D+0
                IF( YLS_BH(2,NBN).GT.PTC_BH(5) )
     &            YLS_BH(2,NBN) = PTC_BH(5)
              ELSE
                IF( YLS_BH(2,NBN)+DPS.LT.0.D+0 ) DPS = 6.D-1*DPS
                YLS_BH(2,NBN) = YLS_BH(2,NBN) + DPS
                IF( YLS_BH(2,NBN).LT.EPSL ) YLS_BH(2,NBN) = 0.D+0
              ENDIF
              XLS_BH(2,NBN) = MIN( YLS_BH(2,NBN),XLSMX )
            ENDIF
!
!---      Unsaturated system
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ_BH(2,NBN).EQ.2 ) THEN
            INDX = 0
            CALL REGION_4( T_BH(2,NBN),PSWX,INDX )
!
!---        Isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Limit salt mass fraction changes to 0.25 of the
!             maximum value if salt mass fraction is less than
!             the maximum   ---
!
              CALL SOL_BRNS( T_BH(2,NBN),PSWX,XLSMX )
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  ) THEN
                DPX = ABS(2.5D-1*PTC_BH(5))
              ELSE
                DPX = ABS(2.5D-1*XLSMX)
              ENDIF
              IF( YLS_BH(2,NBN).LT.XLSMX ) THEN
                DPS = SIGN( MIN( DPX,ABS(DPS) ),DPS )
              ENDIF
!
!---          Zero negative corrections for zero aqueous salt  ---
!
              IF( YLS_BH(2,NBN)/EPSL.LT.EPSL .AND.
     &          BLU(MPS)/EPSL.LT.EPSL ) THEN
                BLU(MPS) = 0.D+0
                DPS = 0.D+0
              ENDIF
              IF( ITC_BH.GE.20 .AND. ITC_BH.LT.30  ) THEN
                IF( YLS_BH(2,NBN)+DPS.GT.PTC_BH(5) ) DPS = 6.D-1*DPS
                YLS_BH(2,NBN) = YLS_BH(2,NBN) + DPS
                IF( YLS_BH(2,NBN).LT.EPSL ) YLS_BH(2,NBN) = 0.D+0
                IF( YLS_BH(2,NBN).GT.PTC_BH(5) )
     &            YLS_BH(2,NBN) = PTC_BH(5)
              ELSE
                IF( YLS_BH(2,NBN)+DPS.LT.0.D+0 ) DPS = 6.D-1*DPS
                YLS_BH(2,NBN) = YLS_BH(2,NBN) + DPS
                IF( YLS_BH(2,NBN).LT.EPSL ) YLS_BH(2,NBN) = 0.D+0
              ENDIF
              XLS_BH(2,NBN) = MIN( YLS_BH(2,NBN),XLSMX )
            ENDIF
!
!---        Pipe flow  ---
!
            IF( IS_BH(INV).GT.10 ) THEN
              ENPR = 0.D+0
            ELSE
!
!---        Assign gas entry pressure and minimum gas saturation
!           for transition to unsaturated conditions  ---
!
              ENPR = 0.D+0
              IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
                ENPR = SCHR(1,IZN)*RHORL*GRAV
              ELSEIF( ISCHR(IZN).EQ.4 ) THEN
                ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
              ENDIF
            ENDIF
!
!---        Limit changes in pressure  ---
!
            IF( ABS(DPL).GT.2.5D-1*(PL_BH(2,NBN)+PATM)
     &        .AND. IERR.EQ.0 ) IERR = 1
            IF( ABS(DPG).GT.2.5D-1*(PG_BH(2,NBN)+PATM)
     &        .AND. IERR.EQ.0 ) IERR = 1
            DPX = MAX( 1.D+6,2.5D-1*(PG_BH(2,NBN)-PL_BH(2,NBN)) )
            DPCX = MAX( 1.D+4,2.5D-1*(PG_BH(2,NBN)-PL_BH(2,NBN)) )
            IF( (DPG-DPL).GT.DPCX ) THEN
              DPG = DPG - 5.D-1*DPCX
              DPL = DPL + 5.D-1*DPCX
            ELSE
              DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
              DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
            ENDIF
!
!---        Relax pressure updates when transitioning to unsaturated
!           conditions  ---
!
            IF( (PG_BH(2,NBN)+DPG)-(PL_BH(2,NBN)+DPL).LT.ENPR ) THEN
              DPG = 6.D-1*DPG
              DPL = 6.D-1*DPL
            ENDIF
            PG_BH(2,NBN) = PG_BH(2,NBN) + DPG
            PL_BH(2,NBN) = PL_BH(2,NBN) + DPL
!
!---        Pipe flow  ---
!
            IF( IS_BH(INV).GT.10 ) THEN
              PL_BH(2,NBN) = MAX( PL_BH(2,NBN),
     &          (PG_BH(2,NBN)-HDOD*RHORL*GRAV) )
            ELSE
              PL_BH(2,NBN) = MAX( PL_BH(2,NBN),
     &          (PG_BH(2,NBN)-SCHR(12,IZN)*RHORL*GRAV) )
            ENDIF
!
!---        Maintain the gas pressure above or at the water vapor
!           pressure  ---
!
            CALL P_IAPWS( T_BH(2,NBN),PSWX,RHOGWX,RHOLWX,HGWX,
     &        HLWX,UGWX,ULWX )
            CALL DENS_B( XLS_BH(2,NBN),RHOBX,RHOLWX,T_BH(2,NBN) )
            CALL SP_B( T_BH(2,NBN),XLS_BH(2,NBN),PSBX )
            PGX = PG_BH(2,NBN) + PATM
            PLX = PL_BH(2,NBN) + PATM
            PCX = PGX - PLX
            CALL VPL_B( T_BH(2,NBN),PSBX,PCX,RHOBX,PVBX )
            PG_BH(2,NBN) = MAX( PG_BH(2,NBN),(PVBX-PATM) )
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - gas pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ_BH(2,NBN).EQ.4 ) THEN
!
!---        Limit changes in water vapor partial pressure  ---
!
            INDX = 0
            CALL REGION_4( T_BH(2,NBN),PWX,INDX )
            PWX = MIN( PWX,PCRW )
            DPX = 2.5D-1*PWX
            DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
            PVW_BH(2,NBN) = MAX( PVW_BH(2,NBN)+DPL,1.D-6 )
!
!---        Limit changes in gas pressure  ---
!
            IF( ISLC(37).EQ.0 ) THEN
              DPX = 2.5D-1*MAX(PG_BH(2,NBN)+PATM,1.D+6)
              DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
              PG_BH(2,NBN) = PG_BH(2,NBN) + DPG
            ENDIF
!
!---        Isobrine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Zero negative corrections for salt volumetric conc.  ---
!
              IF( TMS_BH(2,NBN)/EPSL.LT.EPSL .AND.
     &          BLU(MPS)/EPSL.LT.EPSL ) THEN
                BLU(MPS) = 0.D+0
                DPS = 0.D+0
              ENDIF
              TMS_BH(2,NBN) = TMS_BH(2,NBN) + DPS
              IF( TMS_BH(2,NBN).LT.1.D-9 ) TMS_BH(2,NBN) = 0.D+0
            ENDIF
          ENDIF
!
!---      Check for excessive pressure or temperature   ---
!
          PGX = PG_BH(2,NBN)+PATM
          PLX = PL_BH(2,NBN)+PATM
          TKX = T_BH(2,NBN)+TABS
          IF( PGX.GT.100.D+6 .OR. PGX.LT.0.D+0
     &      .AND. IERR.EQ.0 ) IERR = 1
          IF( PLX.GT.100.D+6 .AND. IERR.EQ.0 ) IERR = 1
          IF( TKX.GT.1073.15D+0 .OR. TKX.LT.TABS 
     &      .AND. IERR.EQ.0 ) IERR = 1
          IF( IERR.EQ.1 ) NSD_BH(1) = NBN
        ENDDO
      ENDDO
!
!---  Reduce time step for excessive changes in primary variables   ---
!
      IF( IERR.EQ.1 ) THEN
        ICNV = 1
        NBN = NSD_BH(1)
        WRITE(ISC,'(10X,A)') '---  Excessive Primary Variable Change '
     &    // 'for Borehole Flow ---'
        WRITE(IWR,'(10X,A)') '---  Excessive Primary Variable Change '
     &    // 'for Borehole Flow ---'
        WRITE(ISC,'(4X,A,I6)') 'Borehole Node = ',NBN
        WRITE(IWR,'(4X,A,I6)') 'Borehole Node = ',NBN
        WRITE(ISC,'(4X,2A)') 'Phase Condition = ',
     &    PH_CND(NPHAZ_BH(2,NBN))
        WRITE(IWR,'(4X,2A)') 'Phase Condition = ',
     &    PH_CND(NPHAZ_BH(2,NBN))
        WRITE(ISC,'(4X,A,1PE12.5)') 'Temperature = ',
     &    T_BH(2,NBN)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Temperature = ',
     &    T_BH(2,NBN)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',
     &    PL_BH(2,NBN)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',
     &    PL_BH(2,NBN)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Pressure = ',
     &    PG_BH(2,NBN)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Pressure = ',
     &    PG_BH(2,NBN)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG_BH(2,NBN)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG_BH(2,NBN)
        WRITE(ISC,'(4X,A,1PE12.5)')
     &    'Aqueous-Air Mole Fraction = ',XMLA_BH(2,NBN)
        WRITE(IWR,'(4X,A,1PE12.5)')
     &    'Aqueous-Air Mole Fraction = ',XMLA_BH(2,NBN)
        WRITE(ISC,'(4X,A,1PE12.5)')
     &    'Water Vapor Pressure = ',PVW_BH(2,NBN)
        WRITE(IWR,'(4X,A,1PE12.5)')
     &    'Water Vapor Pressure = ',PVW_BH(2,NBN)
        WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &    'Total-Salt Aqu. Mass Fraction = ',YLS_BH(2,NBN)
        WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &    'Total-Salt Aqu. Mass Fraction = ',YLS_BH(2,NBN)
        WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &    'Salt Volumetric Concentration = ',TMS_BH(2,NBN)
        WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &    'Salt Volumetric Concentration = ',TMS_BH(2,NBN)
      ENDIF
!
!---  Reduce time step  ---
!
  300 CONTINUE
      IF( ICNV.EQ.1 ) THEN
        IF( NTSR.LT.4 .OR. (DTCF*DT).GT.DTMN ) THEN
          NTSR = NTSR + 1
          DTX = DT
          TM = TM - (1.D+0-DTCF)*DT
          DT = DTCF*DT
          DTO = DT
          DTI = 1.D+0/DT
          VAR = DT
          VARX = DTX
          IF( UNTM.NE.'null' ) THEN
            INDX = 1
            IUNS = 1
            CALL RDUNIT(UNTM,VAR,INDX)
            IUNS = 1
            CALL RDUNIT(UNTM,VARX,INDX)
            NCH = INDEX( UNTM,'  ')-1
          ENDIF
          WRITE(ISC,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
          WRITE(IWR,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
!
!---      Loop over boreholes  ---
!
          DO NBH = 1,N_BH
!
!---        Skip for boundary condition type boreholes  ---
!
            IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---        Loop over borehole nodes  ---
!
            DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
              T_BH(2,NBN) = T_BH(1,NBN)
              PL_BH(2,NBN) = PL_BH(1,NBN)
              PG_BH(2,NBN) = PG_BH(1,NBN)
              PVW_BH(2,NBN) = PVW_BH(1,NBN)
              PVA_BH(2,NBN) = PVA_BH(1,NBN)
              XMLA_BH(2,NBN) = XMLA_BH(1,NBN)
              SL_BH(2,NBN) = SL_BH(1,NBN)
              SG_BH(2,NBN) = SG_BH(1,NBN)
              YLS_BH(2,NBN) = YLS_BH(1,NBN)
              TMS_BH(2,NBN) = TMS_BH(1,NBN)
              NPHAZ_BH(2,NBN) = NPHAZ_BH(1,NBN)
            ENDDO
          ENDDO
!
!---      Fracture flow and transport solution  ---
!
          IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) THEN
!
!---        Loop over fractures  ---
!
            DO NFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---            Skip inactive triangles  ---
!
                IF( IXP_FRC(NTX).EQ.0 ) CYCLE
                T_FRC(2,NTX) = T_FRC(1,NTX)
                PL_FRC(2,NTX) = PL_FRC(1,NTX)
                PG_FRC(2,NTX) = PG_FRC(1,NTX)
                PVW_FRC(2,NTX) = PVW_FRC(1,NTX)
                PVA_FRC(2,NTX) = PVA_FRC(1,NTX)
                XMLA_FRC(2,NTX) = XMLA_FRC(1,NTX)
                SL_FRC(2,NTX) = SL_FRC(1,NTX)
                SG_FRC(2,NTX) = SG_FRC(1,NTX)
                YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
                TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
                NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
              ENDDO
            ENDDO
          ENDIF
!
!---      Matrix flow  ---
!
          DO N = 1,NFLD
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PVW(2,N) = PVW(1,N)
            PVA(2,N) = PVA(1,N)
            XMLA(2,N) = XMLA(1,N)
            SL(2,N) = SL(1,N)
            SG(2,N) = SG(1,N)
            SGT(2,N) = SGT(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            ASLMIN(2,N) = ASLMIN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
          ENDDO
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
!
!---      Loop over boreholes  ---
!
          DO NBH = 1,N_BH
!
!---        Skip for boundary condition type boreholes  ---
!
           IF( IT_BH(1,NBH).GE.21.AND.IT_BH(1,NBH).LE.29 ) CYCLE
!
!---       Loop over borehole nodes  ---
!
            DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
              T_BH(2,NBN) = T_BH(1,NBN)
              PL_BH(2,NBN) = PL_BH(1,NBN)
              PG_BH(2,NBN) = PG_BH(1,NBN)
              PVW_BH(2,NBN) = PVW_BH(1,NBN)
              PVA_BH(2,NBN) = PVA_BH(1,NBN)
              XMLA_BH(2,NBN) = XMLA_BH(1,NBN)
              SL_BH(2,NBN) = SL_BH(1,NBN)
              SG_BH(2,NBN) = SG_BH(1,NBN)
              YLS_BH(2,NBN) = YLS_BH(1,NBN)
              TMS_BH(2,NBN) = TMS_BH(1,NBN)
              NPHAZ_BH(2,NBN) = NPHAZ_BH(1,NBN)
            ENDDO
          ENDDO
!
!---      Fracture flow and transport solution  ---
!
          IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) THEN
!
!---        Loop over fractures  ---
!
            DO NFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---            Skip inactive triangles  ---
!
                IF( IXP_FRC(NTX).EQ.0 ) CYCLE
                T_FRC(2,NTX) = T_FRC(1,NTX)
                PL_FRC(2,NTX) = PL_FRC(1,NTX)
                PG_FRC(2,NTX) = PG_FRC(1,NTX)
                PVW_FRC(2,NTX) = PVW_FRC(1,NTX)
                PVA_FRC(2,NTX) = PVA_FRC(1,NTX)
                XMLA_FRC(2,NTX) = XMLA_FRC(1,NTX)
                SL_FRC(2,NTX) = SL_FRC(1,NTX)
                SG_FRC(2,NTX) = SG_FRC(1,NTX)
                YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
                TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
                NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
              ENDDO
            ENDDO
          ENDIF
!
!---      Matrix flow  ---
!
          DO N = 1,NFLD
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PVW(2,N) = PVW(1,N)
            PVA(2,N) = PVA(1,N)
            XMLA(2,N) = XMLA(1,N)
            SL(2,N) = SL(1,N)
            SG(2,N) = SG(1,N)
            SGT(2,N) = SGT(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            ASLMIN(2,N) = ASLMIN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
          ENDDO
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     &from Borehole Flow  ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     &from Borehole Flow  ---'
          ICNV = 4
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_FRC_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WROBDA_EGS1
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
!     Write observational data for PEST for EGS Collab Experiment 1
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 20 January 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE UCODE
      USE TRNSPT
      USE SOLTN
      USE PARM_BH
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_FRC
      USE FDVP_BH
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/WROBDA_EGS1'
!
!---  Open observational data file for PEST  ---
!
      OPEN(UNIT=36,FILE='egs_collab_exp1_sim.dat',STATUS='UNKNOWN',
     &  FORM='FORMATTED')
      CLOSE(UNIT=36,STATUS='DELETE')
      OPEN(UNIT=36, FILE='egs_collab_exp1_sim.dat', STATUS='NEW',
     &  FORM='FORMATTED')
!
!---  E1-PB Flow Rate, ml/min  ---
!
      NCX = 1
      NBN1X = 1
      NBN2X = 2
      DBB1X = DBBM_BH(NCX)
      DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
      INV1X = INV_BH(NBN1X)
      INV2X = INV_BH(NBN2X)
      IF( IS_BH(INV1X).GT.10 ) THEN
        AP1X = GPI*(PAR_BH(3,INV1X)**2)
      ELSE
        AP1X = GPI*(PAR_BH(2,INV1X)**2)
      ENDIF
      IF( IS_BH(INV2X).GT.10 ) THEN
        AP2X = GPI*(PAR_BH(3,INV2X)**2)
      ELSE
        AP2X = GPI*(PAR_BH(2,INV2X)**2)
      ENDIF
      INDX = -3
      APX = DIFMN( AP1X,AP2X,DBB1X,DBB2X,ZERO,INDX )
      SFX = UBBL(1,NCX)*APX
      VARX = -60.D+6*SFX
      WRITE(36,'(A)') 'E1-PB Flow Rate, ml/min '
      WRITE(36,'(1PE14.7)') VARX
!
!---  E1-PI Flow Rate, ml/min  ---
!
      NCX = 21
      NBN1X = 12
      NBN2X = 13
      DBB1X = DBBM_BH(NCX)
      DBB2X = (DBB_BH(NCX)-DBBM_BH(NCX))
      INV1X = INV_BH(NBN1X)
      INV2X = INV_BH(NBN2X)
      IF( IS_BH(INV1X).GT.10 ) THEN
        AP1X = GPI*(PAR_BH(3,INV1X)**2)
      ELSE
        AP1X = GPI*(PAR_BH(2,INV1X)**2)
      ENDIF
      IF( IS_BH(INV2X).GT.10 ) THEN
        AP2X = GPI*(PAR_BH(3,INV2X)**2)
      ELSE
        AP2X = GPI*(PAR_BH(2,INV2X)**2)
      ENDIF
      INDX = -3
      APX = DIFMN( AP1X,AP2X,DBB1X,DBB2X,ZERO,INDX )
      SFX = UBBL(1,NCX)*APX
      VARX = -60.D+6*SFX
      WRITE(36,'(A)') 'E1-PI Flow Rate, ml/min '
      WRITE(36,'(1PE14.7)') VARX
!
!---  E1-I Pressure, MPa  ---
!
      VARX = 1.D-6*(PL_BH(2,29)+PATM)
      WRITE(36,'(A)') 'E1-I Pressure, MPa '
      WRITE(36,'(1PE14.7)') VARX
!
!---  Hydraulic Fracture Pressure, MPa  ---
!
      VARX = 1.D-6*(PL_FRC(2,144)+PATM)
      WRITE(36,'(A)') 'Hydraulic Fracture Pressure, MPa '
      WRITE(36,'(1PE14.7)') VARX
!
!---  E1-PB Temperature, ˚C  ---
!
      VARX = T_BH(2,6)
      WRITE(36,'(A)') 'E1-PB Temperature, ˚C '
      WRITE(36,'(1PE14.7)') VARX
!
!---  E1-PI Temperature, ˚C  ---
!
      VARX = T_BH(2,21)
      WRITE(36,'(A)') 'E1-PI Temperature, ˚C '
      WRITE(36,'(1PE14.7)') VARX
!
!---  E1-PB Peak Tracer Concentration Arrival Time, hr  ---
!
      VARX = MAX(R_OBDT(2,1)-R_OBDT(5,1),0.D+0)/3600.D+0
      WRITE(36,'(A)') 'E1-PB Peak Tracer Concentration Arrival Time, hr'
      WRITE(36,'(1PE14.7)') VARX
!
!---  E1-PI Peak Tracer Concentration Arrival Time, hr  ---
!
      VARX = MAX(R_OBDT(4,1)-R_OBDT(5,1),0.D+0)/3600.D+0
      WRITE(36,'(A)') 'E1-PI Peak Tracer Concentration Arrival Time, hr'
      WRITE(36,'(1PE14.7)') VARX
!
!---  Close observational data file for PEST  ---
!
      CLOSE(UNIT=36)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WROBDA_EGS1 group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDPREP_BH( NBH )
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
!     Read preprocessed input file for boreholes
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 23 February 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_BH
      USE OUTPU
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
!
!----------------------Data Statements---------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDPREP_BH'
!
!---  Borehole type and borehole name  ---
!
      READ(35,*) (ID_BH(M,NBH),M=1,6),(IT_BH(M,NBH),M=1,2),ITM_BH(NBH),
     &  NM_BH(NBH)
!
!---  Loop over number of borehole intervals  ---
!
      DO NIBH = ID_BH(1,NBH),ID_BH(2,NBH)
        READ(35,*) (XTP_BH(M,NIBH),M=1,2),(YTP_BH(M,NIBH),M=1,2),
     &    (ZTP_BH(M,NIBH),M=1,2),(PAR_BH(M,NIBH),M=1,16)
      ENDDO
!
!---  Borehole sources  ---
!
      READ(35,*) ((VAR_BH(M,N,NBH),M=1,11),N=1,ITM_BH(NBH))
      READ(35,*) ((VARC_BH(M,N,NBH),M=1,2*LSOLU_BH),N=1,ITM_BH(NBH))
!
!---  Loop over borehole nodes  ---
!
      DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
        READ(35,*) INV_BH(NBN)
        READ(35,*) (XP_BH(M,NBN),M=1,2),(YP_BH(M,NBN),M=1,2),
     &    (ZP_BH(M,NBN),M=1,2),IBN_BH(NBN),(IPB_BH(M,NBN),M=1,3),
     &    VOL_BH(NBN),DIST_BH(NBN),DBN_BH(NBN),PLX_BH(NBN),
     &    PLY_BH(NBN),PLZ_BH(NBN),(XE_BH(M,NBN),M=1,8),
     &    (YE_BH(M,NBN),M=1,8),(ZE_BH(M,NBN),M=1,8),
     &    IZ_BH(NBN),IS_BH(INV_BH(NBN))
        DO ICX = 1,2
          NCX = IPB_BH(ICX,NBN)
          IF( NCX.LE.0 ) CYCLE
          READ(35,*) IBCM_BH(NCX),(SRFN_BH(M,NCX),M=1,3),
     &      DBB_BH(NCX),DBBM_BH(NCX)
        ENDDO
      ENDDO
!
!---  Skip for boundary condition type boreholes  ---
!
      IF( IT_BH(1,NBH).LT.21.OR.IT_BH(1,NBH).GT.29 ) THEN
!
!---    Loop over borehole node to fracture triangle connections  ---
!
        DO NCX = ID_BH(5,NBH),ID_BH(6,NBH)
          READ(35,*) IBHN_FRC(NCX),IBHT_FRC(NCX),
     &      XP_BF(NCX),YP_BF(NCX),ZP_BF(NCX)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDPREP_BH group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WRPREP_BH( NBH )
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
!     Write preprocessed input file for boreholes
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 22 February 2023.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_BH
      USE OUTPU
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
!
!----------------------Data Statements---------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WRPREP_BH'
!
!---  Borehole type and borehole name  ---
!
      WRITE(36,*) (ID_BH(M,NBH),M=1,6),(IT_BH(M,NBH),M=1,2),ITM_BH(NBH),
     &  NM_BH(NBH)
!
!---  Loop over number of borehole intervals  ---
!
      DO NIBH = ID_BH(1,NBH),ID_BH(2,NBH)
        WRITE(36,*) (XTP_BH(M,NIBH),M=1,2),(YTP_BH(M,NIBH),M=1,2),
     &    (ZTP_BH(M,NIBH),M=1,2),(PAR_BH(M,NIBH),M=1,16)
      ENDDO
!
!---  Borehole sources  ---
!
      WRITE(36,*) ((VAR_BH(M,N,NBH),M=1,11),N=1,ITM_BH(NBH))
      WRITE(36,*) ((VARC_BH(M,N,NBH),M=1,2*LSOLU_BH),N=1,ITM_BH(NBH))
!
!---  Loop over borehole nodes  ---
!
      DO NBN = ID_BH(3,NBH),ID_BH(4,NBH)
        WRITE(36,*) INV_BH(NBN)
        WRITE(36,*) (XP_BH(M,NBN),M=1,2),(YP_BH(M,NBN),M=1,2),
     &    (ZP_BH(M,NBN),M=1,2),IBN_BH(NBN),(IPB_BH(M,NBN),M=1,3),
     &    VOL_BH(NBN),DIST_BH(NBN),DBN_BH(NBN),PLX_BH(NBN),
     &    PLY_BH(NBN),PLZ_BH(NBN),(XE_BH(M,NBN),M=1,8),
     &    (YE_BH(M,NBN),M=1,8),(ZE_BH(M,NBN),M=1,8),
     &    IZ_BH(NBN),IS_BH(INV_BH(NBN))
        DO ICX = 1,2
          NCX = IPB_BH(ICX,NBN)
          IF( NCX.LE.0 ) CYCLE
          WRITE(36,*) IBCM_BH(NCX),(SRFN_BH(M,NCX),M=1,3),
     &      DBB_BH(NCX),DBBM_BH(NCX)
        ENDDO
      ENDDO
!
!---  Skip for boundary condition type boreholes  ---
!
      IF( IT_BH(1,NBH).LT.21.OR.IT_BH(1,NBH).GT.29 ) THEN
!
!---    Loop over borehole node to fracture triangle connections  ---
!
        DO NCX = ID_BH(5,NBH),ID_BH(6,NBH)
          WRITE(36,*) IBHN_FRC(NCX),IBHT_FRC(NCX),
     &      XP_BF(NCX),YP_BF(NCX),ZP_BF(NCX)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WRPREP_BH group ---
!
      RETURN
      END



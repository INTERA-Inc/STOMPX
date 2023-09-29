!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_FRC_EOR
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
!     Check the thermodynamic and hydrologic states declared through
!     user inputs for fault triangles.
!
!     VIC_FRC(1,NFX) - x reference point, m
!     VIC_FRC(2,NFX) - y reference point, m
!     VIC_FRC(3,NFX) - z reference point, m
!     VIC_FRC(5,NFX) - temperature, C
!     VIC_FRC(6,NFX) - x gradient, C/m
!     VIC_FRC(7,NFX) - y gradient, C/m
!     VIC_FRC(8,NFX) - z gradient, C/m
!     VIC_FRC(9,NFX) - pressure, Pa
!     VIC_FRC(10,NFX) - x gradient, Pa/m
!     VIC_FRC(11,NFX) - y gradient, Pa/m
!     VIC_FRC(12,NFX) - z gradient, Pa/m
!     VIC_FRC(13,NFX) - aqueous saturation or water vapor rel. sat.
!     VIC_FRC(14,NFX) - x gradient, 1/m
!     VIC_FRC(15,NFX) - y gradient, 1/m
!     VIC_FRC(16,NFX) - z gradient, 1/m
!     VIC_FRC(17,NFX) - salt mass frac. or salt total mass, kg
!     VIC_FRC(18,NFX) - x gradient, 1/m or kg/m
!     VIC_FRC(19,NFX) - y gradient, 1/m or kg/m
!     VIC_FRC(20,NFX) - z gradient, 1/m or kg/m
!     VIC_FRC(21,NFX) - CO2 rel. sat. or mass frac.
!     VIC_FRC(22,NFX) - x gradient, 1/m
!     VIC_FRC(23,NFX) - y gradient, 1/m
!     VIC_FRC(24,NFX) - z gradient, 1/m
!     VIC_FRC(21+(IGC*4),NFX) - petroleum component total mole fraction
!     VIC_FRC(22+(IGC*4),NFX) - x gradient, 1/m
!     VIC_FRC(23+(IGC*4),NFX) - y gradient, 1/m
!     VIC_FRC(24+(IGC*4),NFX) - z gradient, 1/m
!     VIC_FRC(21+(NGC*4)+(NSL*4),NFX) - solute concentration, 1/m^3
!     VIC_FRC(21+(NGC*4)+(NSL*4),NFX) - x gradient, 1/m^4
!     VIC_FRC(21+(NGC*4)+(NSL*4),NFX) - y gradient, 1/m^4
!     VIC_FRC(21+(NGC*4)+(NSL*4),NFX) - z gradient, 1/m^4
!     VIC_FRC(21+(NGC*4)+(NSOLU*4)+(NSP*4),NFX) - species conc., 1/m^3
!     VIC_FRC(21+(NGC*4)+(NSOLU*4)+(NSP*4),NFX) - x gradient, 1/m^4
!     VIC_FRC(21+(NGC*4)+(NSOLU*4)+(NSP*4),NFX) - y gradient, 1/m^4
!     VIC_FRC(21+(NGC*4)+(NSOLU*4)+(NSP*4),NFX) - z gradient, 1/m^4
!
!----------------------Authors-----------------------------------------!
!
!     Written by Ramesh Sarathi, PNNL, 6 March 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP
      USE FDVGC_FRC
      USE FDVGC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*16 FORM1
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 / '(A,I1,A,I1,A,I1)' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK_FRC_EOR'
!
!---  Set fault width (i.e., aperture) and volume, 
!     loop over faults  ---
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
!---      Constant fault width (i.e., aperture) model  ---
!
          DO M = 1,ISVC+2
            APM_FRC(M,NTX) = MAX( JRC_FRC(NTX),1.D-9 )
            APH_FRC(M,NTX) = MAX( JRC_FRC(NTX),1.D-12 )
          ENDDO
          VOL_FRC(NTX) = AF_FRC(NTX)*APM_FRC(2,NTX)
!
!---      Loop over fault triangle to grid-cell connections, removing
!         fault volume from grid-cell volume, checking for null
!         grid-cell volume  ---
!
          DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
            N = INCM_FRC(NCX)
            VOL(N) = VOL(N) - AFN_FRC(NCX)*APM_FRC(2,NTX)
            IF( VOL(N).LT.EPSL ) THEN
              CARD = 'Fault Geometry Card'
              INDX = 4
              WRITE(FORM1(5:5),'(I1)') ICOUNT( N )
              WRITE(FORM1(10:10),'(I1)') ICOUNT( NFX )
              WRITE(FORM1(15:15),'(I1)') ICOUNT( NTX )
              WRITE(CHMSG(1:),FORM1) 'Zero Node Volume w/ Removal ' //
     &          'of Fault Volume @ Node: ',N,' Fault: ',NFX,
     &          ' Fault Triangle: ',NTX
              CALL WRMSGS( INDX )
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!---  Skip for restart simulations  ---
!
      IF( IEO.EQ.2 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compositional option  ---
!
!
!---  K-factors computed from fugacity coefficients  ---
!
      IFK = 0
!
!---  Check that the proper combination of variables is specified  ---
!

!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Check IC1  ---
!
        IF( IC_OPT_FRC(1,NFX).EQ.1 ) THEN
!
!---      Temperature unspecified  ---
!
          IF( IC_OPT_FRC(9,NFX).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 1: ' //
     &        'Temperature Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Pressure unspecified  ---
!
          IF( IC_OPT_FRC(10,NFX).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 1: ' //
     &        'Pressure Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous Saturation specified  ---
!
          IF( IC_OPT_FRC(3,NFX).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 1: ' //
     &        'Aqueous Saturation was specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Petroleum component total mole fractions specified  ---
!
          IF( IC_OPT_FRC(13,NFX).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 1: ' //
     &        'Component Total Mole Fractions were specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Default the aqueous CO2 and salt DOF to mass frac ---
!
          IF( IC_OPT_FRC(7,NFX).EQ.0) IC_OPT_FRC(7,NFX) = 3
          IF( IC_OPT_FRC(8,NFX).EQ.0) IC_OPT_FRC(8,NFX) = 3
          
!
!---    Check IC2  ---
!
        ELSEIF( IC_OPT_FRC(1,NFX).EQ.2 ) THEN
!
!---      Temperature unspecified  ---
!
          IF( IC_OPT_FRC(9,NFX).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 2: ' //
     &        'Temperature Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Pressure unspecified  ---
!
          IF( IC_OPT_FRC(10,NFX).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 2: ' //
     &        'Pressure Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous Saturation unspecified  ---
!
          IF( IC_OPT_FRC(3,NFX).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 2: ' //
     &        'Aqueous Saturation Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous CO2 Concentration specified  ---
!
          IF( IC_OPT_FRC(7,NFX).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 2: ' //
     &        'CO2 Aqueous rel sat/mass frac specified ' //
     &        'when only the total mole frac should be.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Petroleum Component Total Mole Fractions unspecified  ---
!
          IF( IC_OPT_FRC(13,NFX).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 2: ' //
     &        'Component Total Mole Fractions Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Default the aqueous salt DOF option to mass frac ---
!
          IF( IC_OPT_FRC(8,NFX).EQ.0) IC_OPT_FRC(8,NFX) = 3
          
!
!---    Check and load IC3  ---
!
        ELSEIF( IC_OPT_FRC(1,NFX).EQ.3 ) THEN
!
!---      Temperature unspecified  ---
!
          IF( IC_OPT_FRC(9,NFX).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 3: ' //
     &        'Temperature Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Pressure unspecified  ---
!
          IF( IC_OPT_FRC(10,NFX).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 3: ' //
     &        'Pressure Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous Saturation specified  ---
!
          IF( IC_OPT_FRC(3,NFX).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 3: ' //
     &        'Aqueous Saturation was specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Water Vapor Relative Saturation unspecified  ---
!
          IF( IC_OPT_FRC(11,NFX).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 3: ' //
     &        'Water Vapor Relative Saturation unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous CO2 Concentration Specified  ---
!
          IF( IC_OPT_FRC(7,NFX).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 3: ' //
     &        'CO2 Aqueous rel sat/mass frac specified ' //
     &        'when only the total mole frac should be.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Petroleum component total mole fractions unspecified  ---
!
          IF( IC_OPT_FRC(13,NFX).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 3: ' //
     &        'Component Total Mole Fractions Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Set the aqueous salt DOF default option to mass total ---
!
          IF( IC_OPT_FRC(8,NFX).EQ.0 .OR. IC_OPT_FRC(8,NFX).EQ.4) THEN
            IC_OPT_FRC(8,NFX) = 4
          ELSE
            INDX = 4
            CHMSG = 'Specification Error: ' //
     &        'Initial Condition State 3: ' //
     &        'Salt Aqueous rel sat/mass frac specified ' //
     &        'when the total mass should be.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---    End Check IC  ---
!
        ENDIF
!
!---  End loop over faults  ---
!
      ENDDO

!
!---  Load IC values into primary variable arrays  ---
!

!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Load IC1  ---
!
        IF( IC_OPT_FRC(1,NFX).EQ.1 ) THEN
!
!---      Loop over fault triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
            ! Temperature
            II = 1+1*4
            T_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
          
            ! Pressure
            II = 1+2*4
            PSO_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
            
            ! Saturation
            SG_FRC(2,NTX) = 0.D0
            SL_FRC(2,NTX) = 1.D0
            SN_FRC(2,NTX) = 0.D0
          
            ! Salt rel_sat/mass_frac/mass_total !rss_ask TMS or YLS
            II = 1+4*4
            TMS_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
          
            ! CO2 rel_sat/mass_frac
            II = 1+5*4
            PVA_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
          
!
!---      End loop over fault triangles  ---
!
          ENDDO
!
!---    Load IC2  ---
!
        ELSEIF( IC_OPT_FRC(1,NFX).EQ.2 ) THEN
!
!---      Loop over fault triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
            ! Temperature
            II = 1+1*4
            T_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
          
            ! Pressure !rss_ask should this be PL or PSO?
            II = 1+2*4
            PSO_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
            
            ! Saturation
            II = 1+3*4
            SL_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
            
            SG_FRC(2,NTX) = 0.D0
            SN_FRC(2,NTX) = 1.D0 - SL_FRC(2,NTX)
          
            ! Salt rel_sat/mass_frac/mass_total !rss_ask TMS or YLS
            II = 1+4*4
            TMS_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
            
!
!---        Set Component Mole Fractions  ---
!
            DO IGC = 1,NGC+2 ! plus two for CO2, CH4
              II = 1+(7+IGC-1)*4
              ZMC_FRC(IGC,2,NTX) = VIC_FRC(II,NFX) + 
     &          (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &          (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &          (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
            ENDDO
!
!---        Initial mole fractions from specified mass fractions and
!           verify mole fraction summation  ---
!
            IF( IC_OPT_FRC(6,NFX).NE.0 ) IOPTX = -1
            IF( IC_OPT_FRC(13,NFX).NE.0 ) IOPTX = 1
            N_DB = NTX
            CALL INIT_FRAC( ZMC_FRC(1,2,NTX),IOPTX )
            
!
!---      End loop over fault triangles  ---
!
          ENDDO
!
!---    Load IC3  ---
!
        ELSEIF( IC_OPT_FRC(1,NFX).EQ.3 ) THEN
!
!---      Loop over fault triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---        Load IC values into primary variable arrays  ---
!
            ! Temperature
            II = 1+1*4
            T_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
          
            ! Pressure
            II = 1+2*4
            PSO_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
            
            ! Water Vapor Relative Saturation
            II = 1+3*4
            PVW_FRC(2,NTX) = VIC_FRC(II,NFX) + 
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
            
            ! Saturation
            SL_FRC(2,NTX) = 0.D0
            SG_FRC(2,NTX) = 0.D0
            SN_FRC(2,NTX) = 1.D0
          
            ! Salt mass_total
            II = 1+4*4
            TMS_FRC(2,NTX) = VIC_FRC(II,NFX) + !rss_ask TMS or YLS?
     &        (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &        (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &        (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
            
!
!---        Set Component Mole Fractions  ---
!
            DO IGC = 1,NGC+2 ! plus two for CO2, CH4
              II = 1+(7+IGC-1)*4
              ZMC_FRC(IGC,2,NTX) = VIC_FRC(II,NFX) + 
     &          (XP_FRC(NTX)-VIC_FRC(1,NFX))*VIC_FRC(II+1,NFX) +
     &          (YP_FRC(NTX)-VIC_FRC(2,NFX))*VIC_FRC(II+2,NFX) +
     &          (ZP_FRC(NTX)-VIC_FRC(3,NFX))*VIC_FRC(II+3,NFX)
            ENDDO
!
!---        Initial mole fractions from specified mass fractions and
!           verify mole fraction summation  ---
!
            IF( IC_OPT_FRC(6,NFX).NE.0 ) IOPTX = -1
            IF( IC_OPT_FRC(13,NFX).NE.0 ) IOPTX = 1
            N_DB = NTX
            CALL INIT_FRAC( ZMC_FRC(1,2,NTX),IOPTX )
!
!---      End loop over fault triangles  ---
!
          ENDDO
!
!---    Equilibrium initial conditions  ---
!
        ELSEIF( IC_OPT_FRC(1,NFX).EQ.10 ) THEN
!
!---      Loop over fault triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---        Loop over connected nodes to determine the total
!           connected distance  ---
!
            DISTX = 0.D+0
            DO K = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
              DISTX = DISTX + MAX( DFN_FRC(K),1.D-9 )
            ENDDO
            DISTIX = 0.D+0
            DO K = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
              DISTIX = DISTIX + DISTX/MAX( DFN_FRC(K),1.D-9 )
            ENDDO
!
!---        Inverse distance weighted equilibrium values for primary
!           variables  ---
!
            DO K = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
              N = INCM_FRC(K)
              FRACX = (DISTX/MAX( DFN_FRC(K),1.D-9 ))/DISTIX
              PG_FRCX = PG(2,N) - GRAV*(ZP_FRC(NTX)-ZP(N))*RHOG(2,N)
              PN_FRCX = PN(2,N) - GRAV*(ZP_FRC(NTX)-ZP(N))*RHON(2,N)
              PL_FRCX = PL(2,N) - GRAV*(ZP_FRC(NTX)-ZP(N))*RHOL(2,N)
              T_FRC(2,NTX) = T_FRC(2,NTX) + T(2,N)*FRACX
              PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + 
     &          MAX(PG_FRCX,PN_FRCX,PL_FRCX)*FRACX
              SL_FRC(2,NTX) = SL_FRC(2,NTX) + SL(2,N)*FRACX
              PVW_FRC(2,NTX) = PVW_FRC(2,NTX) + PVW(2,N)*FRACX
              TMS_FRC(2,NTX) = TMS_FRC(2,NTX) + TMS(2,N)*FRACX
              YLS_FRC(2,NTX) = YLS_FRC(2,NTX) + YLS(2,N)*FRACX
              SI_FRC(1,NTX) = SI_FRC(1,NTX) + SI(1,N)*FRACX
              SI_FRC(2,NTX) = SI_FRC(2,NTX) + SI(2,N)*FRACX
              DO IGC = 1,NGC+2 ! plus two for CO2, CH4
                ZMC_FRC(IGC,2,NTX) = ZMC_FRC(IGC,2,NTX) + 
     &            ZMC(IGC,2,N)*FRACX
              ENDDO
            ENDDO
!
!---        Normalize petroleum component mole fractions  ---
!
            SZMCX = 0.D+0
            DO IGC = 1,NGC+2 ! plus two for CO2, CH4
              SZMCX = SZMCX + ZMC_FRC(IGC,2,NTX)
            ENDDO
            DO IGC = 1,NGC+2 ! plus two for CO2, CH4
              ZMC_FRC(IGC,2,NTX) = ZMC_FRC(IGC,2,NTX)/SZMCX
            ENDDO
!
!---      End loop over fault triangles  ---
!
          ENDDO
!
!---    End Check IC  ---
!
        ENDIF
!
!---  End loop over faults  ---
!
      ENDDO
      
!
!---  Check variable ranges and flash  ---
!
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
!---      Check for out of range parameters  ---
!
          IF( T_FRC(2,NTX).GT.374.14D+0 .OR. 
     &      T_FRC(2,NTX).LT.0.01D+0 ) THEN
            INDX = 17
            RLMSG = T_FRC(2,NTX)
            CHMSG = 'Out of Range Initial Fault Temperature(C) = '
            CALL WRMSGS( INDX )
          ENDIF
          IF( PSO_FRC(2,NTX).GT.8.D+8-PATM ) THEN
            INDX = 17
            RLMSG = PSO_FRC(2,NTX)+PATM
            CHMSG = 'Out of Range Initial Fault Pressure(Pa) = '
            CALL WRMSGS( INDX )
          ENDIF
          IF( SL_FRC(2,NTX).GT.1.D+0 .OR.
     &      SL_FRC(2,NTX).LT.0.D+0 ) THEN
            INDX = 17
            RLMSG = SL_FRC(2,NTX)
            CHMSG = 'Out of Range Initial Fault Aqu. Saturation = '
            CALL WRMSGS( INDX )
          ENDIF

!
!---      IC1 State  ---
!
          IF( SL_FRC(2,NTX).EQ.1.D0 ) THEN
!
!---        Flash calculation for initial condition #1  ---
!
            PX = PSO_FRC(2,NTX) + PATM
            N_DB = NTX
            
            CALL FLH_IC1( PX,PVA_FRC(2,NTX),
     &        RHOL_FRC(2,NTX),RHOML_FRC(2,NTX),TMS_FRC(2,NTX),
     &        T_FRC(2,NTX),XLA_FRC(2,NTX),XLS_FRC(2,NTX),
     &        XLW_FRC(2,NTX),XMLA_FRC(2,NTX),XMLS_FRC(2,NTX),
     &        XMLW_FRC(2,NTX),YLS_FRC(2,NTX),ZMC_FRC(1,2,NTX),
     &        IC_OPT_FRC(7,NFX),IC_OPT_FRC(8,NFX) )
!
!---        Phase condition  ---
!
            NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 1
!
!---      IC2 State  ---
!
          ELSEIF( SL_FRC(2,NTX).LT.1.D0 .AND. 
     &      SL_FRC(2,NTX).GT.0.D0 ) THEN
!
!---        Flash calculation for initial condition #2  ---
!
            PX = PSO_FRC(2,NTX) + PATM
            N_DB = NTX
            
            CALL FLH_IC2( PX,RHOG_FRC(2,NTX),
     &        RHOMG_FRC(2,NTX),RHOL_FRC(2,NTX),RHOML_FRC(2,NTX),
     &        RHON_FRC(2,NTX),RHOMN_FRC(2,NTX),
     &        SG_FRC(2,NTX),SL_FRC(2,NTX),SN_FRC(2,NTX),TMS_FRC(2,NTX),
     &        T_FRC(2,NTX),XGC_FRC(1,2,NTX),XMGC_FRC(1,2,NTX),
     &        XGW_FRC(2,NTX),XMGW_FRC(2,NTX),XLA_FRC(2,NTX),
     &        XLS_FRC(2,NTX),XLW_FRC(2,NTX),XMLA_FRC(2,NTX),
     &        XMLS_FRC(2,NTX),XMLW_FRC(2,NTX),
     &        XNC_FRC(1,2,NTX),XMNC_FRC(1,2,NTX),YLS_FRC(2,NTX),
     &        ZG_FRC(2,NTX),ZMC_FRC(1,2,NTX),ZN_FRC(2,NTX),
     &        IC_OPT_FRC(8,NFX),IZ_FRC(NTX),NPHAZ_FRC(2,NTX) )
!
!---        Phase condition  ---
!
            IF( SG_FRC(2,NTX).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
            ELSE
              PG_FRC(2,NTX) = PX - PATM
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
            ENDIF
!
!---      IC3 State  ---
!
          ELSEIF( SL_FRC(2,NTX).EQ.0.D0 ) THEN
!
!---        Flash calculation for initial condition #3  ---
!
            PX = PSO_FRC(2,NTX) + PATM
            N_DB = NTX
            
            !rss_ask - this does not make sense
            ! we read in TMS, not, YLS, so this appears wrong
            TMS_FRC(2,NTX) = YLS_FRC(2,NTX)
            IF( TMS_FRC(2,NTX).GT.0.D+0 ) THEN
              CALL SOL_LS( T_FRC(2,NTX),YLS_FRC(2,NTX) )
            ELSE
              YLS_FRC(2,NTX) = 0.D+0
            ENDIF
            
            CALL FLH_IC3( PVW_FRC(2,NTX),PX,
     &        RHOG_FRC(2,NTX),RHOMG_FRC(2,NTX),RHON_FRC(2,NTX),
     &        RHOMN_FRC(2,NTX),SG_FRC(2,NTX),SN_FRC(2,NTX),
     &        T_FRC(2,NTX),XGC_FRC(1,2,NTX),XGW_FRC(2,NTX),
     &        XMGC_FRC(1,2,NTX),XMGW_FRC(2,NTX),XNC_FRC(1,2,NTX),
     &        XMNC_FRC(1,2,NTX),YLS_FRC(2,NTX),
     &        ZG_FRC(2,NTX),ZMC_FRC(1,2,NTX),ZN_FRC(2,NTX),
     &        IZ_FRC(NTX),NPHAZ_FRC(2,NTX) )
!
!---        Phase condition  ---
!
            IF( SG_FRC(2,NTX).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 4
            ELSE
              PG_FRC(2,NTX) = PX - PATM
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 5
            ENDIF
!
!---      End IC State Check  ---
!
          ENDIF
!
!---      Capillary pressure calculation  ---
!
          CALL CAP_EOR( SG_FRC(2,NTX),SL_FRC(2,NTX),SN_FRC(2,NTX),
     &      CPGTX,CPNLX,IZ_FRC(NTX) )         
!
!---      Phase pressure from phase capillary pressures and system
!           pressure  ---
!
          PG_FRC(2,NTX) = PSO_FRC(2,NTX)
          PTX = PG_FRC(2,NTX) - CPGTX
          IF( CPNLX.GT.0.D+0 ) THEN
            PN_FRC(2,NTX) = PTX
            PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
          ELSE
            PL_FRC(2,NTX) = PTX
            PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
          ENDIF
!
!---      Initialize reference pressure for compressibility  ---
!
          IF( CMP(3,IZ_FRC(NTX)).GT.PATM ) THEN
            PCMP_FRC(NTX) = CMP(3,IZ_FRC(NTX))
          ELSEIF( ISLC(61).EQ.0 ) THEN
            PCMP_FRC(NTX) = MAX( PL_FRC(2,NTX),PG_FRC(2,NTX) )+PATM
          ENDIF
!
!---      System pressure, absolute  ---
!
          PX = PSO_FRC(M,NTX) + PATM
!
!---      Porous-media porosity  ---
!
          CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(2,NTX),
     &      PORT_FRC(2,NTX) )
          PORD_FRC(2,NTX) = MAX( PORD_FRC(2,NTX),EPSL )
          PORT_FRC(2,NTX) = MAX( PORT_FRC(2,NTX),PORD_FRC(2,NTX) )
!
!---    Petroleum component total mass non-restart simulations  ---
!
          IF( IEO.NE.2 ) THEN
            TMC_FRC(1,2,NTX) = PORD_FRC(2,NTX)*VOL_FRC(NTX)*
     &        (XLA_FRC(2,NTX)*RHOL_FRC(2,NTX)*SL_FRC(2,NTX) +
     &        XGC_FRC(1,2,NTX)*RHOG_FRC(2,NTX)*SG_FRC(2,NTX) + 
     &        XNC_FRC(1,2,NTX)*RHON_FRC(2,NTX)*SN_FRC(2,NTX))
            DO IGC = 2,NGC+2
              TMC_FRC(IGC,2,NTX) = PORD_FRC(2,NTX)*VOL_FRC(NTX)*
     &          (XGC_FRC(IGC,2,NTX)*RHOG_FRC(2,NTX)*SG_FRC(2,NTX)
     &          + XNC_FRC(IGC,2,NTX)*RHON_FRC(2,NTX)*SN_FRC(2,NTX))
            ENDDO
          ENDIF
!
!---      End loop over fault triangles  ---
!
          ENDDO
!
!---  End loop over faults  ---
!
      ENDDO







!rss_ask - these require things in PROP_FRC_EOR to be calculated first,
! especially the fault aperture
! also, the initialize-old-time-step seems to be done in 
! LDO_FRC_EOR and 
!---  Initializations for normal and restart simulations  ---
! 
!---    Porous-media porosity  ---
!---    Petroleum component total mass  ---
!---    Initialize old-time-step values  ---
!---    Assign dissolved-salt mass fractions, isobrine option
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
      SUBROUTINE DFFGC_FRC_EOR( IGC )
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
!     Diffusive petroleum component gas fluxes for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GEOM_FRC
      USE FLUX_FRC
      USE FLUXC_FRC
      USE FDVP_FRC
      USE FDVG_FRC
      USE FDVGC_FRC
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGC_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
            DXMGC = XMGC_FRC(IGC,2,NT1X) - XMGC_FRC(IGC,2,NT2X)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DF1 = TORG_FRC(MP,NT1X)*PORD_FRC(MP,NT1X)*SG_FRC(MP,NT1X)*
     &          DFGC_FRC(IGC,MP,NT1X)*RHOMG_FRC(MP,NT1X)
              DF2 = TORG_FRC(MN,NT2X)*PORD_FRC(MN,NT2X)*SG_FRC(MN,NT2X)*
     &          DFGC_FRC(IGC,MN,NT2X)*RHOMG_FRC(MN,NT2X)
              INDX = 12
              DFM = DIFMN( DF1,DF2,DFF1X,DFF2X,DXMGC,INDX )
              UFFDGC(IGC,M,NCX) = DFM*(XMGC_FRC(IGC,MP,NT1X) 
     &          - XMGC_FRC(IGC,MN,NT2X))/DFF_FRC(NCX)
              FGC1 = XGC_FRC(IGC,MP,NT1X)*RHOG_FRC(MP,NT1X)
              FGC2 = XGC_FRC(IGC,MN,NT2X)*RHOG_FRC(MN,NT2X)
              INDX = 3
              FGC = DIFMN( FGC1,FGC2,DFF1X,DFF2X,UFFG(1,NCX),INDX )
              UFFGC(IGC,M,NCX) = UFFG(M,NCX)*FGC 
     &          + GCPP(1,IGC)*UFFDGC(IGC,M,NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGC_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGW_FRC_EOR
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
!     Compute the water-vapor molar diffusion rates for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 15 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
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
      SUB_LOG(ISUB_LOG) = '/DFFGW_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
            DXMGW = XMGW_FRC(2,NT1X) - XMGW_FRC(2,NT2X)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DF1 = TORG_FRC(MP,NT1X)*PORD_FRC(MP,NT1X)*SG_FRC(MP,NT1X)*
     &          DFGW_FRC(MP,NT1X)*RHOMG_FRC(MP,NT1X)
              DF2 = TORG_FRC(MN,NT2X)*PORD_FRC(MN,NT2X)*SG_FRC(MN,NT2X)*
     &          DFGW_FRC(MN,NT2X)*RHOMG_FRC(MN,NT2X)
              INDX = 12
              DFM = DIFMN( DF1,DF2,DFF1X,DFF2X,DXMGW,INDX )
              UFFDGW(M,NCX) = DFM*(XMGW_FRC(MP,NT1X)-XMGW_FRC(MN,NT2X))
     &          /DFF_FRC(NCX)
              FGW1 = XGW_FRC(MP,NT1X)*RHOG_FRC(MP,NT1X)
              FGW2 = XGW_FRC(MN,NT2X)*RHOG_FRC(MN,NT2X)
              INDX = 3
              FGW = DIFMN( FGW1,FGW2,DFF1X,DFF2X,UFFG(1,NCX),INDX )
              UFFGW(M,NCX) = UFFG(M,NCX)*FGW + WTMW*UFFDGW(M,NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGW_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLA_FRC_EOR
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
!     Diffusive CO2 aqueous fluxes for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLA_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
            DXMLA = XMLA_FRC(2,NT1X) - XMLA_FRC(2,NT2X)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DF1 = TORL_FRC(MP,NT1X)*PORD_FRC(MP,NT1X)*SL_FRC(MP,NT1X)*
     &          DFLA_FRC(MP,NT1X)*RHOML_FRC(MP,NT1X)
              DF2 = TORL_FRC(MN,NT2X)*PORD_FRC(MN,NT2X)*SL_FRC(MN,NT2X)*
     &          DFLA_FRC(MN,NT2X)*RHOML_FRC(MN,NT2X)
              INDX = 12
              DFM = DIFMN( DF1,DF2,DFF1X,DFF2X,DXMLA,INDX )
              UFFDLA(M,NCX) = DFM*(XMLA_FRC(MP,NT1X) 
     &          - XMLA_FRC(MN,NT2X))/DFF_FRC(NCX)
              FLA1 = XLA_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)
              FLA2 = XLA_FRC(MN,NT2X)*RHOL_FRC(MN,NT2X)
              INDX = 3
              FLA = DIFMN( FLA1,FLA2,DFF1X,DFF2X,UFFL(1,NCX),INDX )
              UFFLA(M,NCX) = UFFL(M,NCX)*FLA + WTMA*UFFDLA(M,NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLA_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLS_FRC_EOR
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
!     Diffusive salt aqueous fluxes for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GEOM_FRC
      USE FLUX_FRC
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
      SUB_LOG(ISUB_LOG) = '/DFFLS_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
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
!---  End of DFFLS_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLW_FRC_EOR
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
!     Diffusive water aqueous fluxes for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLW_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
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
              UFFDLW(M,NCX) = -(UFFDLA(M,NCX) + UFFDS(M,NCX)/WTMS)
              FLW1 = XLW_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)
              FLW2 = XLW_FRC(MN,NT2X)*RHOL_FRC(MN,NT2X)
              INDX = 2
              FLW = DIFMN( FLW1,FLW2,DFF1X,DFF2X,UFFL(1,NCX),INDX)
              UFFLW(M,NCX) = UFFL(M,NCX)*FLW + WTMW*UFFDLW(M,NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLW_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNC_FRC_EOR( IGC )
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
!     Diffusive petroleum component nonaqueous-liquid fluxes for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GEOM_FRC
      USE FLUX_FRC
      USE FLUXC_FRC
      USE FDVP_FRC
      USE FDVN_FRC
      USE FDVGC_FRC
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFNC_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
            DXMGC = XMNC_FRC(IGC,2,NT1X) - XMNC_FRC(IGC,2,NT2X)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DF1 = TORN_FRC(MP,NT1X)*PORD_FRC(MP,NT1X)*SN_FRC(MP,NT1X)*
     &          DFNC_FRC(IGC,MP,NT1X)*RHOMN_FRC(MP,NT1X)
              DF2 = TORN_FRC(MN,NT2X)*PORD_FRC(MN,NT2X)*SN_FRC(MN,NT2X)*
     &          DFNC_FRC(IGC,MN,NT2X)*RHOMN_FRC(MN,NT2X)
              INDX = 12
              DFM = DIFMN( DF1,DF2,DFF1X,DFF2X,DXMNC,INDX )
              UFFDNC(IGC,M,NCX) = DFM*(XMNC_FRC(IGC,MP,NT1X) 
     &          - XMNC_FRC(IGC,MN,NT2X))/DFF_FRC(NCX)
              FNC1 = XNC_FRC(IGC,MP,NT1X)*RHON_FRC(MP,NT1X)
              FNC2 = XNC_FRC(IGC,MN,NT2X)*RHON_FRC(MN,NT2X)
              INDX = 3
              FNC = DIFMN( FNC1,FNC2,DFF1X,DFF2X,UFFN(1,NCX),INDX )
              UFFNC(IGC,M,NCX) = UFFN(M,NCX)*FNC 
     &          + GCPP(1,IGC)*UFFDNC(IGC,M,NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNC_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVG_FRC_EOR
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 15 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE PARM_FRC
      USE JACOB
      USE GLB_PAR
      USE GEOM_FRC
      USE FLUX_FRC
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
      REAL*8 KGM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVG_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Apply fault-fault skin factor for inter-fault
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
              PERM_NT1X = PERMRF_FRC(MP,NT1X)*PERM(1,IZ_FRC(NT1X))
              PERM_NT2X = PERMRF_FRC(MN,NT2X)*PERM(1,IZ_FRC(NT2X))
              INDX = 11
              KGM = DIFMN(PERM_NT1X,PERM_NT2X,DFF1X,DFF2X,HDG,INDX)
              IF( PERM_NT1X/EPSL.LT.EPSL ) KGM = 0.D+0
              IF( PERM_NT2X/EPSL.LT.EPSL ) KGM = 0.D+0
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
!---  End of DRCVG_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVL_FRC_EOR
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 15 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE PARM_FRC
      USE JACOB
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
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVL_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Apply fault-fault skin factor for inter-fault
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
              PERM_NT1X = PERMRF_FRC(MP,NT1X)*PERM(1,IZ_FRC(NT1X))
              PERM_NT2X = PERMRF_FRC(MN,NT2X)*PERM(1,IZ_FRC(NT2X))
              INDX = 11
              KLM = DIFMN(PERM_NT1X,PERM_NT2X,DFF1X,DFF2X,HDL,INDX)
              IF( PERM_NT1X/EPSL.LT.EPSL ) KLM = 0.D+0
              IF( PERM_NT2X/EPSL.LT.EPSL ) KLM = 0.D+0
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
!---  End of DRCVL_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVN_FRC_EOR
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
!     Compute the NAPL Darcy flux from pressure gradients
!     and gravitational body forces for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 15 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE PARM_FRC
      USE JACOB
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 KNM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVN_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Apply fault-fault skin factor for inter-fault
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
              HDNX = PN_FRC(MP,NT1X) - PN_FRC(MN,NT2X) + 
     &          5.D-1*GRAV*(ZP_FRC(NT1X)-ZP_FRC(NT2X))*
     &         (RHON_FRC(MP,NT1X)+RHON_FRC(MN,NT2X))
              IF( M.EQ.1 ) HDN = HDNX
              PERM_NT1X = PERMRF_FRC(MP,NT1X)*PERM(1,IZ_FRC(NT1X))
              PERM_NT2X = PERMRF_FRC(MN,NT2X)*PERM(1,IZ_FRC(NT2X))
              INDX = 11
              KNM = DIFMN(PERM_NT1X,PERM_NT2X,DFF1X,DFF2X,HDN,INDX)
              IF( PERM_NT1X/EPSL.LT.EPSL ) KNM = 0.D+0
              IF( PERM_NT2X/EPSL.LT.EPSL ) KNM = 0.D+0
              INDX = 10
              RKNM = DIFMN(RKN_FRC(MP,NT1X),RKN_FRC(MN,NT2X),
     &          DFF1X,DFF2X,HDN,INDX)
              INDX = 7
              VNM = DIFMN(VISN_FRC(MP,NT1X),VISN_FRC(MN,NT2X),
     &          DFF1X,DFF2X,HDN,INDX)
              UFFN(M,NCX) = KNM*RKNM*HDNX/DFF_FRC(NCX)/VNM/(1.D+0+SKFX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVN_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_FRC_BO
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
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 13 April 2020
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FDVS_FRC
      USE FDVP_FRC
      USE FDVGC_FRC
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_FRC_BO'
      TOLX = 2.5D-2
!
!---  Phase options, compute phase condition, loop over faults  ---
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
!
!---      Absolute temperature  ---
!
          TKX = T_FRC(2,NTX)+TABS
!
!---      Phase condition #1  ---
!
!         SL_FRC = 1.0
!         SG_FRC+SN_FRC = 0.0
!         Aqueuous saturated
!
!         water mass - aqueous (system) pressure, Pa gauge - 
!           PL_FRC,PSO_FRC
!         gas mass - total gas mass - TMC_FRC(1,2,NTX)
!         oil mass - total oil mass - TMC_FRC(2,2,NTX)
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
!
!---        System pressure  ---
!
            PX = PSO_FRC(2,NTX) + PATM
!
!---        Maximum solubility of salt in pure water  ---
!
            CALL SOL_LS( T_FRC(2,NTX),XBSMX )
!
!---        Salt mass fraction in brine  ---
!
            XBSX = MIN( YLS_FRC(2,NTX),XBSMX )
            CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
!
!---        Maximum solubility of CO2 in aqueous phase  ---
!
            PVAX = 0.D+0
            CALL EQUIL( T_FRC(2,NTX),PX,PVAX,XBSX,XLAX,XLSX,XLWX )
!
!---        Solubility of CO2 in aqueous phase  ---
!
            CALL EQUIL( T_FRC(2,NTX),PX,PVA_FRC(2,NTX),XBSX,
     &        XLA_FRC(2,NTX),XLS_FRC(2,NTX),XLW_FRC(2,NTX) )
!
!---        Brine density  ---
!
            CALL DENS_B( T_FRC(2,NTX),PX,XBSX,RHOBX )
!
!---        Aqueous density  ---
!
            CALL DENS_L( T_FRC(2,NTX),RHOBX,XLA_FRC(2,NTX),
     &        RHOL_FRC(2,NTX) )
!
!---        Aqueous mass  ---
!
            TMLX = RHOL_FRC(2,NTX)*PORD_FRC(2,NTX)*VOL_FRC(NTX)
!
!---        Gas + oil mass  ---
!
            TMGOX = TMC_FRC(1,2,NTX) + TMC_FRC(2,2,NTX)
!
!---        Gas + oil mass execeeds 1.e-5 of the aqueous phase mass, 
!           transition to phase condition #2 or #3 ---
!
            IF( TMGOX.GT.(1.D-5*TMLX) ) THEN
              CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
              CALL VPL( PSBX,PVBX,SLX,T_FRC(2,NTX),IZ_FRC(NTX) )
!
!---          Bubble-point pressure, assuming all gas is dissolved  ---
!
              RHOGX = 1.225D+0*GCPP(37,1)
              VGX = TMC_FRC(1,2,NTX)/RHOGX
              RHOOX = 1.D+3*GCPP(37,2)
              VOX = TMC_FRC(2,2,NTX)/RHOOX
              POSMX = VGX/VOX
              CALL BO_BP( PBPX,POSMX,T_FRC(2,NTX) )
!
!---          Pressure above bubble-point pressure  ---
!
              IF( PX.GT.PBPX ) THEN
                SN_FRC(2,NTX) = 1.D-4
                SL_FRC(2,NTX) = 1.D+0 - SN_FRC(2,NTX)
                SG_FRC(2,NTX) = 0.D+0
                POSM_FRC(2,NTX) = POSMX
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
!
!---          Pressure below bubble-point pressure  ---
!
              ELSE
                CALL BO_GOR( PX,POSM_FRC(2,NTX),T_FRC(2,NTX) )
                CALL BO_EQ_GAS( PX,PVBX,RHOG_FRC(2,NTX),T_FRC(2,NTX),
     &            XGC_FRC(1,2,NTX),XGW_FRC(2,NTX) )
                CALL BO_EQ_OIL( RHON_FRC(2,NTX),PX,POSM_FRC(2,NTX),
     &            T_FRC(2,NTX),XNC_FRC(1,2,NTX) )
                SGNX = 1.D-4
                SG_FRC(2,NTX) = SGNX*RHON_FRC(2,NTX)*
     &            (TMC_FRC(1,2,NTX)*XNC_FRC(2,2,NTX) - 
     &            TMC_FRC(2,2,NTX) + TMC_FRC(2,2,NTX)*XNC_FRC(2,2,NTX))/
     &            (RHON_FRC(2,NTX)*TMC_FRC(1,2,NTX)*XNC_FRC(2,2,NTX) 
     &            - RHON_FRC(2,NTX)*TMC_FRC(2,2,NTX)
     &            + RHON_FRC(2,NTX)*TMC_FRC(2,2,NTX)*XNC_FRC(2,2,NTX) + 
     &            RHOG_FRC(2,NTX)*TMC_FRC(2,2,NTX)*XGC_FRC(1,2,NTX))
                SN_FRC(2,NTX) = SGNX - SG_FRC(2,NTX)
                SL_FRC(2,NTX) = 1.D+0-SGNX
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
              ENDIF
            ENDIF
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(2,NTX) = PSO_FRC(2,NTX)
            PTX = PG_FRC(2,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(2,NTX) = PTX
              PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
            ELSE
              PL_FRC(2,NTX) = PTX
              PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #2  ---
!
!         SG_FRC = 0.0
!         SL_FRC > 0.0
!         SN_FRC > 0.0
!         Aqueous and nonaqueous liquid saturated
!
!         water mass - aqueous saturation, SL_FRC
!         gas mass - gas to oil ratio, POSM_FRC
!         oil mass - nonaqueous (system) pressure, Pa gauge - 
!           PN_FRC,PSO_FRC
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
            PX = PSO_FRC(2,NTX) + PATM
!
!---        Maximum solubility of salt in pure water  ---
!
            CALL SOL_LS( T_FRC(2,NTX),XBSMX )
!
!---        Salt mass fraction in brine  ---
!
            XBSX = MIN( YLS_FRC(2,NTX),XBSMX )
!
!---        Nonaqueous flash calculation  ---
!
            CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
            CALL VPL( PSBX,PVBX,SL_FRC(2,NTX),T_FRC(2,NTX),IZ_FRC(NTX) )
!
!---        Bubble-point pressure as a function of gas to oil ratio
!           (m^3/m^3) and temperature (C)  ---
!
            CALL BO_BP( PBPX,POSM_FRC(2,NTX),T_FRC(2,NTX) )
!
!---        Pressure below bubble-point pressure  ---
!
            IF( PX.LT.(PBPX-TOLX*GCPP(38,2)) ) THEN
              SG_FRC(2,NTX) = 1.D-4
              SN_FRC(2,NTX) = SN_FRC(2,NTX)*(1.D+0-SG_FRC(2,NTX))/
     &          (SN_FRC(2,NTX)+SL_FRC(2,NTX))
              SL_FRC(2,NTX) = 1.D+0-SG_FRC(2,NTX)-SN_FRC(2,NTX)
!
!---        Pressure above bubble-point pressure  ---
!
            ELSE
              SG_FRC(2,NTX) = 0.D+0
              SN_FRC(2,NTX) = 1.D+0-SL_FRC(2,NTX)
            ENDIF
!
!---        Phase capillary pressures from phase saturations  ---
!
            CALL CAP_EOR( SG_FRC(2,NTX),SL_FRC(2,NTX),SN_FRC(2,NTX),
     &        CPGTX,CPNLX,IZ_FRC(NTX) )
!
!---        Nonaqueous phases disappears, transition to 
!           phase condition #1, conserving aqueous pressure  ---
!
            IF( (SG_FRC(2,NTX)+SN_FRC(2,NTX)).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 1
!
!---        Aqueous phases disappear, transition to 
!           phase condition #4  ---
!
            ELSEIF( (SG_FRC(2,NTX)+SL_FRC(2,NTX)).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 4
!
!---        Aqueous phase disappears and gas phase appears, 
!           transition to phase condition #5  ---
!
            ELSEIF( SL_FRC(2,NTX).LT.EPSL .AND. 
     &        SG_FRC(2,NTX).GT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 5
!
!---        Gas phase appears, transition to phase condition #3  ---
!
            ELSEIF( SG_FRC(2,NTX).GT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
!
!---        No transition from phase condition #2  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
            ENDIF
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(2,NTX) = PSO_FRC(2,NTX)
            PTX = PG_FRC(2,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(2,NTX) = PTX
              PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
            ELSE
              PL_FRC(2,NTX) = PTX
              PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #3  ---
!
!         SG_FRC > 0.0
!         SL_FRC > 0.0
!         SN_FRC > 0.0
!         Three-phase
!
!         water mass - aqueous saturation, SL_FRC
!         gas mass - gas saturation, SG_FRC
!         oil mass - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
            PX = PSO_FRC(2,NTX) + PATM
!
!---        Maximum solubility of salt in pure water  ---
!
            CALL SOL_LS( T_FRC(2,NTX),XBSMX )
!
!---        Salt mass fraction in brine  ---
!
            XBSX = MIN( YLS_FRC(2,NTX),XBSMX )
!
!---        Nonaqueous flash calculation  ---
!
            CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
            CALL VPL( PSBX,PVBX,SL_FRC(2,NTX),T_FRC(2,NTX),IZ_FRC(NTX) )
!
!---        Aqueous disappears  ---
!
            IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
!
!---        Gas disappears  ---
!
            IF( SG_FRC(2,NTX).LT.EPSL ) THEN
              CALL BO_GOR( PX,POSM_FRC(2,NTX),T_FRC(2,NTX) )
              SG_FRC(2,NTX) = 0.D+0
            ENDIF
!
!---        Nonaqueous-liquid  ---
!
            SN_FRC(2,NTX) = 1.D+0-SL_FRC(2,NTX)-SG_FRC(2,NTX)
!
!---        Capillary pressure from saturations  ---
!
            CALL CAP_EOR( SG_FRC(2,NTX),SL_FRC(2,NTX),SN_FRC(2,NTX),
     &        CPGTX,CPNLX,IZ_FRC(NTX) )
!
!---        Nonaqueous phases disappears, transition to 
!           phase condition #1  ---
!
            IF( (SG_FRC(2,NTX)+SN_FRC(2,NTX)).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 1
!
!---        Gas and aqueous phases disappear, transition to 
!           phase condition #4  ---
!
            ELSEIF( (SG_FRC(2,NTX)+SL_FRC(2,NTX)).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 4
!
!---        Aqueous phase disappears, transition to 
!           phase condition #5  ---
!
            ELSEIF( SL_FRC(2,NTX).LT.EPSL .AND. 
     &        SG_FRC(2,NTX).GT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 5
!
!---        Gas phase disappears, transition to phase condition #2  ---
!
            ELSEIF( SG_FRC(2,NTX).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
!
!---        No transition from phase condition #3  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
            ENDIF
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(2,NTX) = PSO_FRC(2,NTX)
            PTX = PG_FRC(2,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(2,NTX) = PTX
              PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
            ELSE
              PL_FRC(2,NTX) = PTX
              PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #4  ---
!
!         SG_FRC = 0.0
!         SL_FRC = 0.0
!         SN_FRC = 1.0
!         Nonaqueous liquid saturated
!
!         water mass - aqueous saturation, SL_FRC
!         gas mass - gas to oil ratio, POSM_FRC
!         oil mass - nonaqueous (system) pressure, Pa gauge - 
!           PN_FRC,PSO_FRC
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
            PX = PSO_FRC(2,NTX) + PATM
!
!---        Maximum solubility of salt in pure water  ---
!
            CALL SOL_LS( T_FRC(2,NTX),XBSMX )
!
!---        Salt mass fraction in brine  ---
!
            XBSX = MIN( YLS_FRC(2,NTX),XBSMX )
!
!---        Nonaqueous flash calculation  ---
!
            CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
            CALL VPL( PSBX,PVBX,SL_FRC(2,NTX),T_FRC(2,NTX),IZ_FRC(NTX) )
!
!---        Bubble-point pressure as a function of gas to oil ratio
!           (m^3/m^3) and temperature (C)  ---
!
            CALL BO_BP( PBPX,POSM_FRC(2,NTX),T_FRC(2,NTX) )
!
!---        Pressure below bubble-point pressure  ---
!
            IF( PX.LT.(PBPX-TOLX*GCPP(38,2)) ) THEN
              SG_FRC(2,NTX) = 1.D-4
              SN_FRC(2,NTX) = SN_FRC(2,NTX)*(1.D+0-SG_FRC(2,NTX))/
     &          (SN_FRC(2,NTX)+SL_FRC(2,NTX))
              SL_FRC(2,NTX) = 1.D+0-SG_FRC(2,NTX)-SN_FRC(2,NTX)
!
!---        Pressure above bubble-point pressure  ---
!
            ELSE
              SG_FRC(2,NTX) = 0.D+0
              SN_FRC(2,NTX) = 1.D+0-SL_FRC(2,NTX)
            ENDIF
!
!---        Aqueous phase appears, transition to 
!           phase condition #2 or #3  ---
!
            IF( SL_FRC(2,NTX).GT.1.D-6 ) THEN
!
!---          Gas phase appears, transition to phase condition #3  ---
!
              IF( SG_FRC(2,NTX).GT.EPSL ) THEN
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
!
!---          No gas phase, transition to phase condition #2  ---
!
              ELSE
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
              ENDIF
!
!---        Gas phase appears, transition to phase condition #5  ---
!
            ELSEIF( SG_FRC(2,NTX).GT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 5
!
!---        No transition from phase condition #4  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 4
            ENDIF
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(2,NTX) = PSO_FRC(2,NTX)
            PTX = PG_FRC(2,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(2,NTX) = PTX
              PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
            ELSE
              PL_FRC(2,NTX) = PTX
              PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #5  ---
!
!         SG_FRC > 0.0
!         SL_FRC = 0.0
!         SN_FRC > 0.0
!         Nonaqueous liquid and gas saturated
!
!         water mass - aqueous saturation, SL_FRC
!         gas mass - gas saturation, SG_FRC
!         oil mass - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
            PX = PG_FRC(2,NTX) + PATM
!
!---        Maximum solubility of salt in pure water  ---
!
            CALL SOL_LS( T_FRC(2,NTX),XBSMX )
!
!---        Salt mass fraction in brine  ---
!
            XBSX = MIN( YLS_FRC(2,NTX),XBSMX )
!
!---        Nonaqueous flash calculation  ---
!
            CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
            CALL VPL( PSBX,PVBX,SL_FRC(2,NTX),T_FRC(2,NTX),IZ_FRC(NTX) )
!
!---        Aqueous disappears  ---
!
            IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
!
!---        Gas disappears  ---
!
            IF( SG_FRC(2,NTX).LT.EPSL ) THEN
              CALL BO_GOR( PX,POSM_FRC(2,NTX),T_FRC(2,NTX) )
              SG_FRC(2,NTX) = 0.D+0
            ENDIF
!
!---        Nonaqueous-liquid  ---
!
            SN_FRC(2,NTX) = 1.D+0-SL_FRC(2,NTX)-SG_FRC(2,NTX)
!
!---        Capillary pressure from saturations  ---
!
            CALL CAP_EOR( SG_FRC(2,NTX),SL_FRC(2,NTX),SN_FRC(2,NTX),
     &        CPGTX,CPNLX,IZ_FRC(NTX) )
!
!---        Aqueous phase appears, transition to 
!           phase condition #2 or #3  ---
!
            IF( SL_FRC(2,NTX).GT.1.D-6 ) THEN
!
!---          Gas phase disappears, transition to phase condition #2,
!             conserving nonaqueous-liquid pressure  ---
!
              IF( SG_FRC(2,NTX).LT.EPSL ) THEN
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
!
!---          Gas phase remains, transition to phase condition #3,
!             conserving gas pressure  ---
!
              ELSE
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
              ENDIF
!
!---        Gas phase disappears, transition to phase condition #4, 
!           conserving nonaqueous-liquid pressure  ---
!
            ELSEIF( SG_FRC(2,NTX).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 4
!
!---        No transition from phase condition #5  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 5
            ENDIF
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(2,NTX) = PSO_FRC(2,NTX)
            PTX = PG_FRC(2,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(2,NTX) = PTX
              PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
            ELSE
              PL_FRC(2,NTX) = PTX
              PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
!---  Compute increments, loop over faults  ---
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
          IZN = IZ_FRC(NTX)
          N_DB = NTX
!
!---      Phase condition #1  ---
!
!         SL_FRC = 1.0
!         SG_FRC+SN_FRC = 0.0
!         Aqueuous saturated
!
!         water mass - aqueous (system) pressure, Pa gauge - 
!           PL_FRC,PSO_FRC
!         gas mass - total gas mass - TMC_FRC(1,2,NTX)
!         oil mass - total oil mass - TMC_FRC(2,2,NTX)
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
!
!---        Water mass (system pressure), increment  ---
!
            DNR_FRC(IEQW,NTX) = MAX( 1.D-1,1.D-7*(PL_FRC(2,NTX)+PATM) )
!
!---        Gas mass (total gas mass), increment  ---
!
            DNR_FRC(IEQGC(1),NTX) = MAX( 1.D-7,1.D-7*TMC_FRC(1,2,NTX) )
!
!---        Oil mass (total oil mass), increment  ---
!
            DNR_FRC(IEQGC(2),NTX) = MAX( 1.D-7,1.D-7*TMC_FRC(2,2,NTX) )
!
!---        Salt mass (total salt brine mass fraction), increment  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              CALL SOL_LS( T_FRC(2,NTX),XLSMX )
              XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
              DNR_FRC(IEQS,NTX) = 1.D-7*XLSMX
            ENDIF
!
!---      Phase condition #2  ---
!
!         SG_FRC = 0.0
!         SL_FRC > 0.0
!         SN_FRC > 0.0
!         Aqueous and nonaqueous liquid saturated
!
!         water mass - aqueous saturation, SL_FRC
!         gas mass - gas to oil ratio, POSM_FRC
!         oil mass - nonaqueous (system) pressure, Pa gauge - 
!           PN_FRC,PSO_FRC
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
!
!---        Water mass (aqueous saturation), increment or decrement  ---
!
            DNR_FRC(IEQW,NTX) = SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
!
!---        Gas mass (gas to oil ratio) increment  ---
!
            DNR_FRC(IEQGC(1),NTX) = 1.D-6
!
!---        Oil mass (nonaqueous (system) pressure) increment   ---
!
            DNR_FRC(IEQGC(2),NTX) = MAX( 1.D-1,
     &        1.D-7*(PN_FRC(2,NTX)+PATM) )
!
!---        Salt mass (total salt brine mass fraction), increment  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              CALL SOL_LS( T_FRC(2,NTX),XLSMX )
              XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
              DNR_FRC(IEQS,NTX) = 1.D-7*XLSMX
            ENDIF
!
!---      Phase condition #3  ---
!
!         SG_FRC > 0.0
!         SL_FRC > 0.0
!         SN_FRC > 0.0
!         Three-phase
!
!         water mass - aqueous saturation, SL_FRC
!         gas mass - gas saturation, SG_FRC
!         oil mass - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
!
!---        Water mass (aqueous saturation), increment or decrement  ---
!
            DNR_FRC(IEQW,NTX) = SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
!
!---        Gas mass (gas (system) pressure) increment  ---
!
            DNR_FRC(IEQGC(1),NTX) = SIGN( 1.D-6,5.D-1-SG_FRC(2,NTX) )
!
!---        Oil mass (nonaqueous-liquid saturation) increment  ---
!
            DNR_FRC(IEQGC(2),NTX) = -MAX( 1.D-1,
     &        1.D-7*(PG_FRC(2,NTX)+PATM) )
!
!---        Salt mass (total salt brine mass fraction)  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              CALL SOL_LS( T_FRC(2,NTX),XLSMX )
              XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
              DNR_FRC(IEQS,NTX) = 1.D-7*XLSMX
            ENDIF
!
!---      Phase condition #4  ---
!
!         SG_FRC = 0.0
!         SL_FRC = 0.0
!         SN_FRC = 1.0
!         Nonaqueous liquid saturated
!
!         water mass - aqueous saturation, SL_FRC
!         gas mass - gas to oil ratio, POSM_FRC
!         oil mass - nonaqueous (system) pressure, Pa gauge - 
!           PN_FRC,PSO_FRC
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
!
!---        Water mass (aqueous saturation), increment or decrement  ---
!
            DNR_FRC(IEQW,NTX) = SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
!
!---        Gas mass (gas to oil ratio) increment  ---
!
            DNR_FRC(IEQGC(1),NTX) = 1.D-6
!
!---        Oil mass (nonaqueous (system) pressure) increment  ---
!
            DNR_FRC(IEQGC(2),NTX) = MAX( 1.D-1,
     &        1.D-7*(PN_FRC(2,NTX)+PATM) )
!
!---        Salt mass (total salt mass), increment  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              DNR_FRC(IEQS,NTX) = 1.D-7
            ENDIF
!
!---      Phase condition #5  ---
!
!         SG_FRC > 0.0
!         SL_FRC = 0.0
!         SN_FRC > 0.0
!         Nonaqueous liquid and gas saturated
!
!         water mass - aqueous saturation, SL_FRC
!         gas mass - gas saturation, SG_FRC
!         oil mass - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
!
!---        Water mass (aqueous saturation), increment or decrement  ---
!
            DNR_FRC(IEQW,NTX) = SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
!
!---        Gas mass (gas (system) pressure) increment  ---
!
            DNR_FRC(IEQGC(1),NTX) = SIGN( 1.D-6,5.D-1-SG_FRC(2,NTX) )
!
!---        Oil mass (nonaqueous-liquid saturation) increment  ---
!
            DNR_FRC(IEQGC(2),NTX) = -MAX( 1.D-1,
     &        1.D-7*(PG_FRC(2,NTX)+PATM) )
!
!---        Salt mass (total salt mass)  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              DNR_FRC(IEQS,NTX) = 1.D-7
            ENDIF
          ENDIF
!
!---      Increment the primary variables  ---
!
          DO M = 3,ISVC+2
            T_FRC(M,NTX) = T_FRC(2,NTX)
            PG_FRC(M,NTX) = PG_FRC(2,NTX)
            PL_FRC(M,NTX) = PL_FRC(2,NTX)
            PN_FRC(M,NTX) = PN_FRC(2,NTX)
            POSM_FRC(M,NTX) = POSM_FRC(2,NTX)
            PSO_FRC(M,NTX) = PSO_FRC(2,NTX)
            PVA_FRC(M,NTX) = PVA_FRC(2,NTX)
            PVW_FRC(M,NTX) = PVW_FRC(2,NTX)
            SG_FRC(M,NTX) = SG_FRC(2,NTX)
            SL_FRC(M,NTX) = SL_FRC(2,NTX)
            SN_FRC(M,NTX) = SN_FRC(2,NTX)
            TMS_FRC(M,NTX) = TMS_FRC(2,NTX)
            YLS_FRC(M,NTX) = YLS_FRC(2,NTX)
            DO IGC = 1,NGC+2
              TMC_FRC(IGC,M,NTX) = TMC_FRC(IGC,2,NTX)
            ENDDO
!
!---        Phase condition #1  ---
!
!           SL_FRC = 1.0
!           SG_FRC+SN_FRC = 0.0
!           Aqueuous saturated
!
!           water mass - aqueous (system) pressure, Pa gauge - 
!             PL_FRC,PSO_FRC
!           gas mass - total gas mass - TMC_FRC(1,2,NTX)
!           oil mass - total oil mass - TMC_FRC(2,2,NTX)
!           salt mass - total salt brine mass fraction - YLS_FRC
!
            IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
              IF( M.EQ.IEQW+2 ) THEN
                PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
              ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
                TMC_FRC(1,M,NTX) = TMC_FRC(1,M,NTX) + 
     &            DNR_FRC(IEQGC(1),NTX)
              ELSEIF( M.EQ.IEQGC(2)+2 ) THEN
                TMC_FRC(2,M,NTX) = TMC_FRC(2,M,NTX) + 
     &            DNR_FRC(IEQGC(2),NTX)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                YLS_FRC(M,NTX) = YLS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ENDIF
!
!---        Phase condition #2  ---
!
!           SG_FRC = 0.0
!           SL_FRC > 0.0
!           SN_FRC > 0.0
!           Aqueous and nonaqueous liquid saturated
!
!           water mass - aqueous saturation, SL_FRC
!           gas mass - gas to oil ratio, POSM_FRC
!           oil mass - nonaqueous (system) pressure, Pa gauge - PN_FRC,
!             PSO_FRC
!           salt mass - total salt brine mass fraction - YLS_FRC
!
            ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
              IF( M.EQ.IEQW+2 ) THEN
                SL_FRC(M,NTX) = SL_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
              ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
                POSM_FRC(M,NTX) = POSM_FRC(M,NTX) + 
     &            DNR_FRC(IEQGC(1),NTX)
              ELSEIF( M.EQ.IEQGC(2)+2 ) THEN
                PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + DNR_FRC(IEQGC(2),NTX)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                YLS_FRC(M,NTX) = YLS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ENDIF
!
!---        Phase condition #3  ---
!
!           SG_FRC > 0.0
!           SL_FRC > 0.0
!           SN_FRC > 0.0
!           Three-phase
!
!           water mass - aqueous saturation, SL_FRC
!           gas mass - gas saturation, SG_FRC
!           oil mass - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           salt mass - total salt brine mass fraction - YLS_FRC
!
            ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
              IF( M.EQ.IEQW+2 ) THEN
                SL_FRC(M,NTX) = SL_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
              ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
                SG_FRC(M,NTX) = SG_FRC(M,NTX) + DNR_FRC(IEQGC(1),NTX)
              ELSEIF( M.EQ.IEQGC(2)+2 ) THEN
                PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + DNR_FRC(IEQGC(2),NTX)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                YLS_FRC(M,NTX) = YLS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ENDIF
!
!---        Phase condition #4  ---
!
!           SG_FRC = 0.0
!           SL_FRC = 0.0
!           SN_FRC = 1.0
!           Nonaqueous liquid saturated
!
!           water mass - aqueous saturation, SL_FRC
!           gas mass - gas to oil ratio, POSM_FRC
!           oil mass - nonaqueous (system) pressure, Pa gauge - 
!             PN_FRC,PSO_FRC
!           salt mass - total salt mass - TMS_FRC(2,NTX)
!
            ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
              IF( M.EQ.IEQW+2 ) THEN
                SL_FRC(M,NTX) = SL_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
              ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
                POSM_FRC(M,NTX) = POSM_FRC(M,NTX) 
     &            + DNR_FRC(IEQGC(1),NTX)
              ELSEIF( M.EQ.IEQGC(2)+2 ) THEN
                PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + DNR_FRC(IEQGC(2),NTX)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                TMS_FRC(M,NTX) = TMS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ENDIF
!
!---      Phase condition #5  ---
!
!         SG_FRC > 0.0
!         SL_FRC = 0.0
!         SN_FRC > 0.0
!         Nonaqueous liquid and gas saturated
!
!         water mass - aqueous saturation, SL_FRC
!         gas mass - gas saturation, SG_FRC
!         oil mass - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
            ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
              IF( M.EQ.IEQW+2 ) THEN
                SL_FRC(M,NTX) = SL_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
              ELSEIF( M.EQ.IEQGC(1)+2 ) THEN
                SG_FRC(M,NTX) = SG_FRC(M,NTX) + DNR_FRC(IEQGC(1),NTX)
              ELSEIF( M.EQ.IEQGC(2)+2 ) THEN
                PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + DNR_FRC(IEQGC(2),NTX)
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
!---  End of INCRM_FRC_BO group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_FRC_EOR
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
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Ramesh Sarathi, PNNL, 6 March 2020
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FDVS_FRC
      USE FDVP_FRC
      USE FDVGC_FRC
      USE FDVGC
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 FKX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_FRC_EOR'
!
!---  Black-oil option  ---
!
      IF( ISLC(68).EQ.0 ) THEN
        CALL INCRM_FRC_BO
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compositional option, loop over faults  ---
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
!
!---      K-factors fixed at second iterate values  ---
!
          IF( NITER.GT.2 ) THEN
            IFK = 1
            DO IGC = 1,NGC+2
              FKX(IGC) = FK_FRC(IGC,NTX)
            ENDDO
!
!---      K-factors computed from fugacity coefficients  ---
!
          ELSE
            IFK = 0
          ENDIF
!
!---      Absolute temperature  ---
!
          TKX = T_FRC(2,NTX)+TABS
!
!---      Phase condition #1  ---
!
!         SL_FRC = 1.0, SG_FRC+SN_FRC = 0.0
!
!         energy - temperature, C - T_FRC
!         water mass - aqueous (system) pressure, Pa gauge - 
!           PL_FRC,PSO_FRC
!         CO2 mass - CO2 vapor pressure - PVA_FRC
!         CH4 mass - total CH4 mass - TMC_FRC
!         petroleum component mass - total petroleum component mass - 
!           TMC_FRC
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
!
!---        System pressure  ---
!
            PX = PSO_FRC(2,NTX) + PATM
!
!---        Maximum solubility of salt in pure water  ---
!
            CALL SOL_LS( T_FRC(2,NTX),XBSMX )
!
!---        Salt mass fraction in brine  ---
!
            XBSX = MIN( YLS_FRC(2,NTX),XBSMX )
            CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
!
!---        Maximum solubility of CO2 in aqueous phase  ---
!
            PVAX = PX - PSBX
            CALL EQUIL( T_FRC(2,NTX),PX,PVAX,XBSX,XLAX,XLSX,XLWX )
!
!---        Solubility of CO2 in aqueous phase  ---
!
            CALL EQUIL( T_FRC(2,NTX),PX,PVA_FRC(2,NTX),XBSX,
     &        XLA_FRC(2,NTX),XLS_FRC(2,NTX),XLW_FRC(2,NTX) )
!
!---        Brine density  ---
!
            CALL DENS_B( T_FRC(2,NTX),PX,XBSX,RHOBX )
!
!---        Aqueous density  ---
!
            CALL DENS_L( T_FRC(2,NTX),RHOBX,XLA_FRC(2,NTX),
     &        RHOL_FRC(2,NTX) )
!
!---        Total mass of CO2  ---
!
            TMC_FRC(1,2,NTX) = XLA_FRC(2,NTX)*RHOL_FRC(2,NTX)*
     &        PORD_FRC(2,NTX)*VOL_FRC(NTX)
!
!---        Non-CO2 petroleum mass and total petroleum molar mass  ---
!
            STMCX = -TMC_FRC(1,2,NTX)
            STMMCX = 0.D+0
            DO IGC = 1,NGC+2
              STMCX = STMCX + TMC_FRC(IGC,2,NTX)
              STMMCX = STMMCX + TMC_FRC(IGC,2,NTX)/GCPP(1,IGC)
            ENDDO
            IF( STMMCX.GT.EPSL ) THEN
              ZMCMX = 0.D+0
              DO IGC = 1,NGC+2
                ZMC_FRC(IGC,2,NTX) = TMC_FRC(IGC,2,NTX)/GCPP(1,IGC)
     &            /STMMCX
                IF( ZMC_FRC(IGC,2,NTX).GT.ZMCMX ) THEN
                  IZMC_FRC(NTX) = IGC
                  ZMCMX = ZMC_FRC(IGC,2,NTX)
                ENDIF
              ENDDO
            ELSE
              ZMC_FRC(1,2,NTX) = 1.D+0
              IZMC_FRC(NTX) = 1
              DO IGC = 2,NGC+2
                ZMC_FRC(IGC,2,NTX) = 0.D+0
              ENDDO
            ENDIF
!
!---        Aqueous mass  ---
!
            TMLX = RHOL_FRC(2,NTX)*PORD_FRC(2,NTX)*VOL_FRC(NTX)
!
!---        Mixture critical properties ---
!
            CALL CRIT_PROP( PAFX,PCRX,TBPX,TCRX,VMCRX,
     &        ZMC_FRC(1,2,NTX),ZCRX )
!
!---        Mixture vapor pressure (Pa), Riedel equation ---
!
            TRX = TKX/TCRX
            TBRX = TKX/TBPX
            PSI_BX = -3.5D+1 + 3.6D+1/TBRX + 4.2D+1*LOG(TBRX) - TBRX**6
            ALPHA_CX = (3.15D-1*PSI_BX + LOG(PCRX/PATM))/
     &        (8.38D-2*PSI_BX - LOG(TBRX))
            QX = 8.38D-2*(3.758D+0 - ALPHA_CX)
            AX = -3.5D+1*QX
            BX = -3.6D+1*QX
            CX = 4.2D+1*QX + ALPHA_CX
            DX = -QX
            PVPX = PATM*EXP( AX - BX/TRX + CX*LOG(TRX) + DX*(TRX**6) )
!
!---        Air concentration from vapor pressure  ---
!
            CASX = PVPX/(RCU*TKX)
!
!---        Carbon number from air concentration from Eastcott et al.  ---
!
            CNX = -(LOG10(CASX) - 4.35D+0)/5.7D-1
!
!---        Aqueous solubility (mol/cm^3) as a function of carbon 
!           number for aromatics from Eastcott et al.  ---
!
            CLSX = 1.D+1**(1.88D+0 - 2.4D-1*CNX)
!
!---        Aqueous concentration of organics (mol/cm^3)  ---
!
            CLCX = 1.D-3*STMMCX/(PORD_FRC(2,NTX)*VOL_FRC(NTX))
!
!---        CO2 concentration in aqueous phase exceeds solubility limit
!           or total petroleum component mass execeeds solubility limit,
!           transition to phase condition #2 or #3 ---
!
             IF( XLA_FRC(2,NTX).GT.XLAX .OR. CLCX.GT.CLSX ) THEN
              CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
              CALL VPL( PSBX,PVBX,SLX,T_FRC(2,NTX),IZ_FRC(NTX) )
              PVX = PX - PVBX
              IBETAX = 0
              BETAX = -1.D+0
              GNIFTX = 0.D+0
              CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVBX,RHOMG_FRC(2,NTX),
     &          RHOMN_FRC(2,NTX),SGX,SNX,T_FRC(2,NTX),XMGC_FRC(1,2,NTX),
     &          XMGW_FRC(2,NTX),XMNC_FRC(1,2,NTX),ZGX,ZMC_FRC(1,2,NTX),
     &          ZNX,IBETAX,NPHAZ_FRC(2,NTX) )
              SGNX = 1.D-4
              SL_FRC(2,NTX) = 1.D+0-SGNX
              SG_FRC(2,NTX) = (1.D+0-SL_FRC(2,NTX))*SGX
              SN_FRC(2,NTX) = (1.D+0-SL_FRC(2,NTX))*SNX
              AST_FRC(NTX) = SGX/(SGX+SNX+SMALL)
              CALL CAP_EOR( SG_FRC(2,NTX),SL_FRC(2,NTX),SN_FRC(2,NTX),
     &          CPGTX,CPNLX,IZ_FRC(NTX) )
              IF( SG_FRC(2,NTX).LT.EPSL ) THEN
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
              ELSE
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
              ENDIF
            ELSE
              CPGTX = 0.D+0
              CPNLX = 0.D+0
            ENDIF
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(2,NTX) = PSO_FRC(2,NTX)
            PTX = PG_FRC(2,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(2,NTX) = PTX
              PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
            ELSE
              PL_FRC(2,NTX) = PTX
              PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #2  ---
!
!         SL_FRC < 1.0, ST > 0.0, SG_FRC = 0.0
!
!         energy - temperature, C - T_FRC
!         water mass - aqueous saturation, SL_FRC
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
            PX = PSO_FRC(2,NTX) + PATM
!
!---        Maximum solubility of salt in pure water  ---
!
            CALL SOL_LS( T_FRC(2,NTX),XBSMX )
!
!---        Salt mass fraction in brine  ---
!
            XBSX = MIN( YLS_FRC(2,NTX),XBSMX )
!
!---        Nonaqueous flash calculation  ---
!
            CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
            CALL VPL( PSBX,PVBX,SL_FRC(2,NTX),T_FRC(2,NTX),IZ_FRC(NTX) )
            PVX = PX - PVBX
            IBETAX = 0
            BETAX = -1.D+0
            GNIFTX = 0.D+0
            CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVBX,RHOMG_FRC(2,NTX),
     &        RHOMN_FRC(2,NTX),SGX,SNX,T_FRC(2,NTX),XMGC_FRC(1,2,NTX),
     &        XMGW_FRC(2,NTX),XMNC_FRC(1,2,NTX),ZGX,ZMC_FRC(1,2,NTX),
     &        ZNX,IBETAX,NPHAZ_FRC(2,NTX) )
!            IF( NITER.GE.2 ) THEN
!              SGX = AST_FRC(NTX)
!              SNX = 1.D+0 - AST_FRC(NTX)
!            ELSE
!              AST_FRC(NTX) = SGX/(SNX+SGX+SMALL)
!            ENDIF
            SG_FRC(2,NTX) = (1.D+0-SL_FRC(2,NTX))*SGX
            SN_FRC(2,NTX) = (1.D+0-SL_FRC(2,NTX))*SNX
            CALL CAP_EOR( SG_FRC(2,NTX),SL_FRC(2,NTX),SN_FRC(2,NTX),
     &        CPGTX,CPNLX,IZ_FRC(NTX) )
!
!---        Assign maximum concentration of petroleum component  ---
!
            ZMCMX = 0.D+0
            DO IGC = 1,NGC+2
              IF( ZMC_FRC(IGC,2,NTX).GT.ZMCMX ) THEN
                IZMC_FRC(NTX) = IGC
                ZMCMX = ZMC_FRC(IGC,2,NTX)
              ENDIF
            ENDDO
!
!---        Nonaqueous phases disappears, transition to 
!           phase condition #1, conserving aqueous pressure  ---
!
            IF( (SG_FRC(2,NTX)+SN_FRC(2,NTX)).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 1
!
!---        Aqueous phases disappear, transition to 
!           phase condition #4  ---
!
            ELSEIF( (SG_FRC(2,NTX)+SL_FRC(2,NTX)).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 4
!
!---        Aqueous phase disappears and gas phase appears, transition
!           to phase condition #5  ---
!
            ELSEIF( SL_FRC(2,NTX).LT.EPSL .AND. 
     &        SG_FRC(2,NTX).GT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 5
!
!---        Gas phase appears, transition to phase condition #3  ---
!
            ELSEIF( SG_FRC(2,NTX).GT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
!
!---        No transition from phase condition #2  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
            ENDIF
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(2,NTX) = PSO_FRC(2,NTX)
            PTX = PG_FRC(2,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(2,NTX) = PTX
              PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
            ELSE
              PL_FRC(2,NTX) = PTX
              PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #3  ---
!
!         SL_FRC < 1.0, ST >= 0.0, SG_FRC > 0.0
!
!         energy - temperature, C - T_FRC
!         water mass - aqueous saturation, SL_FRC
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
            PX = PSO_FRC(2,NTX) + PATM
!
!---        Maximum solubility of salt in pure water  ---
!
            CALL SOL_LS( T_FRC(2,NTX),XBSMX )
!
!---        Salt mass fraction in brine  ---
!
            XBSX = MIN( YLS_FRC(2,NTX),XBSMX )
!
!---        Nonaqueous flash calculation  ---
!
            CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
            CALL VPL( PSBX,PVBX,SL_FRC(2,NTX),T_FRC(2,NTX),IZ_FRC(NTX) )
            PVX = PX - PVBX
            IBETAX = 0
            BETAX = -1.D+0
            GNIFTX = 0.D+0
            CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVBX,RHOMG_FRC(2,NTX),
     &        RHOMN_FRC(2,NTX),SGX,SNX,T_FRC(2,NTX),XMGC_FRC(1,2,NTX),
     &        XMGW_FRC(2,NTX),XMNC_FRC(1,2,NTX),ZGX,ZMC_FRC(1,2,NTX),
     &        ZNX,IBETAX,NPHAZ_FRC(2,NTX) )
!            IF( NITER.GE.2 ) THEN
!              SGX = AST_FRC(NTX)
!              SNX = 1.D+0 - AST_FRC(NTX)
!            ELSE
!              AST_FRC(NTX) = SGX/(SNX+SGX+SMALL)
!            ENDIF
            SG_FRC(2,NTX) = (1.D+0-SL_FRC(2,NTX))*SGX
            SN_FRC(2,NTX) = (1.D+0-SL_FRC(2,NTX))*SNX
            CALL CAP_EOR( SG_FRC(2,NTX),SL_FRC(2,NTX),SN_FRC(2,NTX),
     &        CPGTX,CPNLX,IZ_FRC(NTX) )
!
!---        Assign maximum concentration of petroleum component  ---
!
            ZMCMX = 0.D+0
            DO IGC = 1,NGC+2
              IF( ZMC_FRC(IGC,2,NTX).GT.ZMCMX ) THEN
                IZMC_FRC(NTX) = IGC
                ZMCMX = ZMC_FRC(IGC,2,NTX)
              ENDIF
            ENDDO
!
!---        Nonaqueous phases disappears, transition to 
!           phase condition #1  ---
!
            IF( (SG_FRC(2,NTX)+SN_FRC(2,NTX)).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 1
!
!---        Gas and aqueous phases disappear, transition to 
!           phase condition #4  ---
!
            ELSEIF( (SG_FRC(2,NTX)+SL_FRC(2,NTX)).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 4
!
!---        Aqueous phase disappears, transition to 
!           phase condition #5  ---
!
            ELSEIF( SL_FRC(2,NTX).LT.EPSL .AND. 
     &        SG_FRC(2,NTX).GT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 5
!
!---        Gas phase disappears, transition to phase condition #2  ---
!
            ELSEIF( SG_FRC(2,NTX).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
!
!---        No transition from phase condition #3  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
            ENDIF
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(2,NTX) = PSO_FRC(2,NTX)
            PTX = PG_FRC(2,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(2,NTX) = PTX
              PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
            ELSE
              PL_FRC(2,NTX) = PTX
              PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #4  ---
!
!         SL_FRC = 0.0, ST = 1.0, SG_FRC = 0.0
!
!         energy - temperature, C - T_FRC
!         water mass - aqueous saturation, SL_FRC
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
            PX = PSO_FRC(2,NTX) + PATM
!
!---        Maximum solubility of salt in pure water  ---
!
            CALL SOL_LS( T_FRC(2,NTX),XBSMX )
!
!---        Salt mass fraction in brine  ---
!
            XBSX = MIN( YLS_FRC(2,NTX),XBSMX )
!
!---        Nonaqueous flash calculation  ---
!
            CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
            CALL VPL( PSBX,PVBX,SL_FRC(2,NTX),T_FRC(2,NTX),IZ_FRC(NTX) )
            PVX = PX - PVBX
            IBETAX = 0
            BETAX = -1.D+0
            GNIFTX = 0.D+0
            CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVBX,RHOMG_FRC(2,NTX),
     &        RHOMN_FRC(2,NTX),SGX,SNX,T_FRC(2,NTX),XMGC_FRC(1,2,NTX),
     &        XMGW_FRC(2,NTX),XMNC_FRC(1,2,NTX),ZGX,ZMC_FRC(1,2,NTX),
     &        ZNX,IBETAX,NPHAZ_FRC(2,NTX) )
!            IF( NITER.GE.2 ) THEN
!              SGX = AST_FRC(NTX)
!              SNX = 1.D+0 - AST_FRC(NTX)
!            ELSE
!              AST_FRC(NTX) = SGX/(SNX+SGX+SMALL)
!            ENDIF
            SG_FRC(2,NTX) = (1.D+0-SL_FRC(2,NTX))*SGX
            SN_FRC(2,NTX) = (1.D+0-SL_FRC(2,NTX))*SNX
            CALL CAP_EOR( SG_FRC(2,NTX),SL_FRC(2,NTX),SN_FRC(2,NTX),
     &        CPGTX,CPNLX,IZ_FRC(NTX) )
!
!---        Assign maximum concentration of petroleum component  ---
!
            ZMCMX = 0.D+0
            DO IGC = 1,NGC+2
              IF( ZMC_FRC(IGC,2,NTX).GT.ZMCMX ) THEN
                IZMC_FRC(NTX) = IGC
                ZMCMX = ZMC_FRC(IGC,2,NTX)
              ENDIF
            ENDDO
!
!---        Aqueous phase appears, transition to 
!           phase condition #2 or #3  ---
!
            IF( SL_FRC(2,NTX).GT.1.D-6 ) THEN
!
!---          Gas phase appears, transition to phase condition #3  ---
!
              IF( SG_FRC(2,NTX).GT.EPSL ) THEN
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
!
!---          No gas phase, transition to phase condition #2  ---
!
              ELSE
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
              ENDIF
!
!---        Gas phase appears, transition to phase condition #5  ---
!
            ELSEIF( SG_FRC(2,NTX).GT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 5
!
!---        No transition from phase condition #4  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 4
            ENDIF
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(2,NTX) = PSO_FRC(2,NTX)
            PTX = PG_FRC(2,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(2,NTX) = PTX
              PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
            ELSE
              PL_FRC(2,NTX) = PTX
              PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #5  ---
!
!         SL_FRC = 0.0, ST >= 0.0, SG_FRC > 0.0
!
!         energy - temperature, C - T_FRC
!         water mass - aqueous saturation, SL_FRC
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
            PX = PG_FRC(2,NTX) + PATM
!
!---        Maximum solubility of salt in pure water  ---
!
            CALL SOL_LS( T_FRC(2,NTX),XBSMX )
!
!---        Salt mass fraction in brine  ---
!
            XBSX = MIN( YLS_FRC(2,NTX),XBSMX )
!
!---        Nonaqueous flash calculation  ---
!
            CALL SP_B( T_FRC(2,NTX),XBSX,PSBX )
            CALL VPL( PSBX,PVBX,SL_FRC(2,NTX),T_FRC(2,NTX),IZ_FRC(NTX) )
            PVX = PX - PVBX
            IBETAX = 0
            BETAX = -1.D+0
            GNIFTX = 0.D+0
            CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVBX,RHOMG_FRC(2,NTX),
     &        RHOMN_FRC(2,NTX),SGX,SNX,T_FRC(2,NTX),XMGC_FRC(1,2,NTX),
     &        XMGW_FRC(2,NTX),XMNC_FRC(1,2,NTX),ZGX,ZMC_FRC(1,2,NTX),
     &        ZNX,IBETAX,NPHAZ_FRC(2,NTX) )
!            IF( NITER.GE.2 ) THEN
!              SGX = AST_FRC(NTX)
!              SNX = 1.D+0 - AST_FRC(NTX)
!            ELSE
!              AST_FRC(NTX) = SGX/(SNX+SGX+SMALL)
!            ENDIF
            SG_FRC(2,NTX) = (1.D+0-SL_FRC(2,NTX))*SGX
            SN_FRC(2,NTX) = (1.D+0-SL_FRC(2,NTX))*SNX
            CALL CAP_EOR( SG_FRC(2,NTX),SL_FRC(2,NTX),SN_FRC(2,NTX),
     &        CPGTX,CPNLX,IZ_FRC(NTX) )
!
!---        Assign maximum concentration of petroleum component  ---
!
            ZMCMX = 0.D+0
            DO IGC = 1,NGC+2
              IF( ZMC_FRC(IGC,2,NTX).GT.ZMCMX ) THEN
                IZMC_FRC(NTX) = IGC
                ZMCMX = ZMC_FRC(IGC,2,NTX)
              ENDIF
            ENDDO
!
!---        Aqueous phase appears, transition to 
!           phase condition #2 or #3  ---
!
            IF( SL_FRC(2,NTX).GT.1.D-6 ) THEN
!
!---          Gas phase disappears, transition to phase condition #2,
!             conserving nonaqueous-liquid pressure  ---
!
              IF( SG_FRC(2,NTX).LT.EPSL ) THEN
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 2
!
!---          Gas phase remains, transition to phase condition #3,
!             conserving gas pressure  ---
!
              ELSE
                NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 3
              ENDIF
!
!---        Gas phase disappears, transition to phase condition #4, 
!           conserving nonaqueous-liquid pressure  ---
!
            ELSEIF( SG_FRC(2,NTX).LT.EPSL ) THEN
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 4
!
!---        No transition from phase condition #5  ---
!
            ELSE
              NPHAZ_FRC(2,NTX) = 10*(NPHAZ_FRC(2,NTX)/10) + 5
            ENDIF
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(2,NTX) = PSO_FRC(2,NTX)
            PTX = PG_FRC(2,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(2,NTX) = PTX
              PL_FRC(2,NTX) = PN_FRC(2,NTX) - CPNLX
            ELSE
              PL_FRC(2,NTX) = PTX
              PN_FRC(2,NTX) = PL_FRC(2,NTX) + CPNLX
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
!---  Compute increments, loop over faults  ---
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
          IZN = IZ_FRC(NTX)
          N_DB = NTX
!
!---      Absolute temperature  ---
!
          TKX = T_FRC(2,NTX)+TABS
!
!---      Phase condition #1  ---
!
!         SL_FRC = 1.0, SG_FRC+SN_FRC = 0.0
!
!         energy - temperature, C - T_FRC
!         water mass - aqueous (system) pressure, Pa gauge - 
!           PL_FRC,PSO_FRC
!         CO2 mass - CO2 vapor pressure - PVA_FRC
!         CH4 mass - total CH4 mass - TMC_FRC
!         petroleum component mass - total petroleum component mass - 
!           TMC_FRC
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
!
!---        Energy (temperature), increment  ---
!
            IF( ISLC(30).EQ.0 ) DNR_FRC(IEQT,NTX) = -1.D-7
!
!---        Water mass (system pressure), increment  ---
!
            DNR_FRC(IEQW,NTX) = MAX( 1.D-1,1.D-7*(PL_FRC(2,NTX)+PATM))
!
!---        CO2 mass (CO2 vapor pressure), increment  ---
!
            IF( ISLC(45).EQ.0 ) DNR_FRC(IEQGC(1),NTX) = 1.D-2
!
!---        Aqueous mass  ---
!
            TMLX = RHOL_FRC(2,NTX)*PORD_FRC(2,NTX)*VOL_FRC(NTX)
!
!---        CH4 mass (total CH4 mass), increment  ---
!
            IF( ISLC(64).EQ.0 ) 
     &        DNR_FRC(IEQGC(2),NTX) = MAX( 1.D-7*TMLX,1.D-6 )
!
!---        Petroleum component mass (total petroleum component mass),
!           increment  ---
!
            DO IGC = 3,NGC+2
              DNR_FRC(IEQGC(IGC),NTX) = MAX( 1.D-7*TMLX,1.D-6 )
            ENDDO
!
!---        Salt mass (total salt brine mass fraction), increment  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              CALL SOL_LS( T_FRC(2,NTX),XLSMX )
              XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
              DNR_FRC(IEQS,NTX) = 1.D-7*XLSMX
            ENDIF
!
!---      Phase condition #2  ---
!
!         SL_FRC < 1.0, ST > 0.0, SG_FRC = 0.0
!
!         energy - temperature, C - T_FRC
!         water mass - aqueous saturation, SL_FRC
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
!
!---        Energy (temperature), increment  ---
!
            IF( ISLC(30).EQ.0 ) DNR_FRC(IEQT,NTX) = -1.D-7
!
!---        Water mass (aqueous saturation), increment or decrement  ---
!
!            IF( SL_FRC(2,NTX).GT.5.D-1 ) THEN
              NPHAZ_FRC(2,NTX) = MOD( NPHAZ_FRC(2,NTX),100) + 100
              DNR_FRC(IEQW,NTX) = MAX( 1.D-1,1.D-7*PSO_FRC(2,NTX) )
!            ELSE
!              NPHAZ_FRC(2,NTX) = MOD( NPHAZ_FRC(2,NTX),100 )
!              DNR_FRC(IEQW,NTX) = SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
!            ENDIF
!
!---        CO2 mass increment  ---
!
            IF( ISLC(45).EQ.0 ) THEN
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.1 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(1),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(1),NTX) = 
     &              MAX( 1.D-1,1.D-7*(PN_FRC(2,NTX)+PATM) )
                ENDIF
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Decrease mole fraction of volatile components, and 
!               increase mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(1,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(1),NTX) = -MAX( 1.D-6*ZMC_FRC(1,2,NTX),
     &              1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(1),NTX) = MAX( 1.D-6*ZMC_FRC(1,2,NTX),
     &              1.D-9 )
                ENDIF
              ENDIF
            ENDIF
!
!---        CH4 mass increment  ---
!
            IF( ISLC(64).EQ.0 ) THEN
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.2 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(2),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(2),NTX) = 
     &              MAX( 1.D-1,1.D-7*(PN_FRC(2,NTX)+PATM) )
                ENDIF            
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Decrease mole fraction of volatile components, and
!               increase mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(2,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(2),NTX) = -MAX( 1.D-6*ZMC_FRC(2,2,NTX),
     &              1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(2),NTX) = MAX( 1.D-6*ZMC_FRC(2,2,NTX),
     &              1.D-9 )
                ENDIF
              ENDIF
            ENDIF
!
!---        Petroleum component mass (total petroleum component mass),
!           increment  ---
!
            DO IGC = 3,NGC+2
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.IGC ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              MAX( 1.D-1,1.D-7*(PN_FRC(2,NTX)+PATM) )
                ENDIF            
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Decrease mole fraction of volatile components, and
!               increase mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(IGC,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              -MAX(1.D-6*ZMC_FRC(IGC,2,NTX),1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              MAX( 1.D-6*ZMC_FRC(IGC,2,NTX),1.D-9 )
                ENDIF
              ENDIF
            ENDDO
!
!---        Salt mass (total salt brine mass fraction), increment  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              CALL SOL_LS( T_FRC(2,NTX),XLSMX )
              XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
              DNR_FRC(IEQS,NTX) = 1.D-7*XLSMX
            ENDIF
!
!---      Phase condition #3  ---
!
!         SL_FRC < 1.0, ST > 0.0, SG_FRC > 0.0
!
!         energy - temperature, C - T_FRC
!         water mass - aqueous saturation, SL_FRC
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt brine mass fraction - YLS_FRC
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
!
!---        Energy (temperature), increment  ---
!
            IF( ISLC(30).EQ.0 ) DNR_FRC(IEQT,NTX) = 1.D-7
!
!---        Water mass (aqueous saturation), increment or decrement  ---
!
!            IF( SL_FRC(2,NTX).GT.5.D-1 ) THEN
              NPHAZ_FRC(2,NTX) = MOD( NPHAZ_FRC(2,NTX),100) + 100
              DNR_FRC(IEQW,NTX) = MAX( 1.D-1,1.D-7*PSO_FRC(2,NTX) )
!            ELSE
!              NPHAZ_FRC(2,NTX) = MOD( NPHAZ_FRC(2,NTX),100 )
!              DNR_FRC(IEQW,NTX) = SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
!            ENDIF
!
!---        CO2 mass increment  ---
!
            IF( ISLC(45).EQ.0 ) THEN
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.1 ) THEN
              IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                DNR_FRC(IEQGC(1),NTX) = 
     &            SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
              ELSE
                DNR_FRC(IEQGC(1),NTX) = 
     &            -MAX( 1.D-1,1.D-5*(PG_FRC(2,NTX)+PATM) )
              ENDIF
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Increase mole fraction of volatile components, and
!               decrease mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(1,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(1),NTX) = -MAX( 1.D-6*ZMC_FRC(1,2,NTX),
     &              1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(1),NTX) = MAX( 1.D-6*ZMC_FRC(1,2,NTX),
     &              1.D-9 )
                ENDIF
              ENDIF
            ENDIF
!
!---        CH4 mass increment  ---
!
            IF( ISLC(64).EQ.0 ) THEN
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.2 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(2),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(2),NTX) = 
     &              -MAX( 1.D-1,1.D-5*(PG_FRC(2,NTX)+PATM) )
                ENDIF
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Increase mole fraction of volatile components, and
!               decrease mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(2,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(2),NTX) = -MAX( 1.D-6*ZMC_FRC(2,2,NTX),
     &              1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(2),NTX) = MAX( 1.D-6*ZMC_FRC(2,2,NTX),
     &              1.D-9 )
                ENDIF
              ENDIF
            ENDIF
!
!---        Petroleum component mass (total petroleum component mass),
!           increment  ---
!
            DO IGC = 3,NGC+2
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.IGC ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              -MAX( 1.D-1,1.D-5*(PG_FRC(2,NTX)+PATM) )
                ENDIF
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Increase mole fraction of volatile components, and
!               decrease mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(IGC,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              -MAX(1.D-6*ZMC_FRC(IGC,2,NTX),1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              MAX( 1.D-6*ZMC_FRC(IGC,2,NTX),1.D-9 )
                ENDIF
              ENDIF
            ENDDO
!
!---        Salt mass (total salt brine mass fraction)  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              CALL SOL_LS( T_FRC(2,NTX),XLSMX )
              XLSX = MIN( YLS_FRC(2,NTX),XLSMX )
              DNR_FRC(IEQS,NTX) = 1.D-7*XLSMX
            ENDIF
!
!---      Phase condition #4  ---
!
!         SL_FRC = 0.0, ST = 1.0, SG_FRC = 0.0
!
!         energy - temperature, C - T_FRC
!         water mass - aqueous saturation, SL_FRC
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - nonaqueous (system) pressure, Pa gauge - PN_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
!
!---        Energy (temperature), increment  ---
!
            IF( ISLC(30).EQ.0 ) DNR_FRC(IEQT,NTX) = -1.D-7
!
!---        Water mass (aqueous saturation), increment or decrement  ---
!
!            IF( SL_FRC(2,NTX).GT.5.D-1 ) THEN
              NPHAZ_FRC(2,NTX) = MOD( NPHAZ_FRC(2,NTX),100) + 100
              DNR_FRC(IEQW,NTX) = MAX( 1.D-1,1.D-7*PSO_FRC(2,NTX) )
!            ELSE
!              NPHAZ_FRC(2,NTX) = MOD( NPHAZ_FRC(2,NTX),100 )
!              DNR_FRC(IEQW,NTX) = SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
!            ENDIF
!
!---        CO2 mass increment  ---
!
            IF( ISLC(45).EQ.0 ) THEN
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.1 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(1),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(1),NTX) = 
     &              MAX( 1.D-1,1.D-7*(PN_FRC(2,NTX)+PATM) )
                ENDIF
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Decrease mole fraction of volatile components, and
!               increase mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(1,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(1),NTX) = -MAX( 1.D-6*ZMC_FRC(1,2,NTX),
     &              1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(1),NTX) = MAX( 1.D-6*ZMC_FRC(1,2,NTX),
     &              1.D-9 )
                ENDIF
              ENDIF
            ENDIF
!
!---        CH4 mass increment  ---
!
            IF( ISLC(64).EQ.0 ) THEN
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.2 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(2),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(2),NTX) = 
     &              MAX( 1.D-1,1.D-7*(PN_FRC(2,NTX)+PATM) )
                ENDIF
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Decrease mole fraction of volatile components, and
!               increase mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(2,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(2),NTX) = -MAX( 1.D-6*ZMC_FRC(2,2,NTX),
     &              1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(2),NTX) = MAX( 1.D-6*ZMC_FRC(2,2,NTX),
     &              1.D-9 )
                ENDIF
              ENDIF
            ENDIF
!
!---        Petroleum component mass (total petroleum component mass),
!           increment  ---
!
            DO IGC = 3,NGC+2
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.IGC ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              MAX( 1.D-1,1.D-7*(PN_FRC(2,NTX)+PATM) )
                ENDIF
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Decrease mole fraction of volatile components, and
!               increase mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(IGC,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              -MAX(1.D-6*ZMC_FRC(IGC,2,NTX),1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              MAX( 1.D-6*ZMC_FRC(IGC,2,NTX),1.D-9 )
                ENDIF
              ENDIF
            ENDDO
!
!---        Salt mass (total salt mass), increment  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              DNR_FRC(IEQS,NTX) = 1.D-7
            ENDIF
!
!---      Phase condition #5  ---
!
!         SL_FRC = 0.0, ST >= 0.0, SG_FRC > 0.0
!
!         energy - temperature, C - T_FRC
!         water mass - aqueous saturation, SL_FRC
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
!
!---        Energy (temperature), increment  ---
!
            IF( ISLC(30).EQ.0 ) DNR_FRC(IEQT,NTX) = -1.D-7
!
!---        Water mass (aqueous saturation), increment or decrement  ---
!
!            IF( SL_FRC(2,NTX).GT.5.D-1 ) THEN
              NPHAZ_FRC(2,NTX) = MOD( NPHAZ_FRC(2,NTX),100) + 100
              DNR_FRC(IEQW,NTX) = MAX( 1.D-1,1.D-7*PSO_FRC(2,NTX) )
!            ELSE
!              NPHAZ_FRC(2,NTX) = MOD( NPHAZ_FRC(2,NTX),100 )
!              DNR_FRC(IEQW,NTX) = SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
!            ENDIF
!
!---        CO2 mass increment  ---
!
            IF( ISLC(45).EQ.0 ) THEN
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.1 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(1),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(1),NTX) = 
     &              -MAX( 1.D-1,1.D-7*(PG_FRC(2,NTX)+PATM) )
                ENDIF
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Increase mole fraction of volatile components, and
!               decrease mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(1,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(1),NTX) = -MAX(1.D-6*ZMC_FRC(1,2,NTX),
     &              1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(1),NTX) = MAX( 1.D-6*ZMC_FRC(1,2,NTX),
     &              1.D-9 )
                ENDIF
              ENDIF
            ENDIF
!
!---        CH4 mass increment  ---
!
            IF( ISLC(64).EQ.0 ) THEN
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.2 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(2),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(2),NTX) = 
     &              -MAX( 1.D-1,1.D-7*(PG_FRC(2,NTX)+PATM) )
                ENDIF
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
 !
!---            Increase mole fraction of volatile components, and
!               decrease mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(2,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(2),NTX) = -MAX( 1.D-6*ZMC_FRC(2,2,NTX),
     &              1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(2),NTX) = MAX( 1.D-6*ZMC_FRC(2,2,NTX),
     &              1.D-9 )
                ENDIF
              ENDIF
            ENDIF
!
!---        Petroleum component mass (total petroleum component mass),
!           increment  ---
!
            DO IGC = 3,NGC+2
!
!---          System pressure, increment  ---
!
              IF( IZMC_FRC(NTX).EQ.IGC ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              SIGN( 1.D-6,5.D-1-SL_FRC(2,NTX) )
                ELSE
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              -MAX( 1.D-1,1.D-7*(PG_FRC(2,NTX)+PATM) )
                ENDIF
!
!---          Nonaqueous mole fraction, increment  ---
!
              ELSE
!
!---            Increase mole fraction of volatile components, and
!               decrease mole fraction of non-volatile components  ---
!
                IF( ZMC_FRC(IGC,2,NTX).GT.9.D-1 ) THEN
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              -MAX(1.D-6*ZMC_FRC(IGC,2,NTX),1.D-9 )
                ELSE
                  DNR_FRC(IEQGC(IGC),NTX) = 
     &              MAX( 1.D-6*ZMC_FRC(IGC,2,NTX),1.D-9 )
                ENDIF
              ENDIF
            ENDDO
!
!---        Salt mass (total salt mass)  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              DNR_FRC(IEQS,NTX) = 1.D-7
            ENDIF
          ENDIF
!
!---      Increment the primary variables  ---
!
          DO M = 3,ISVC+2
            T_FRC(M,NTX) = T_FRC(2,NTX)
            PG_FRC(M,NTX) = PG_FRC(2,NTX)
            PL_FRC(M,NTX) = PL_FRC(2,NTX)
            PN_FRC(M,NTX) = PN_FRC(2,NTX)
            PSO_FRC(M,NTX) = PSO_FRC(2,NTX)
            PVA_FRC(M,NTX) = PVA_FRC(2,NTX)
            PVW_FRC(M,NTX) = PVW_FRC(2,NTX)
            SL_FRC(M,NTX) = SL_FRC(2,NTX)
            TMS_FRC(M,NTX) = TMS_FRC(2,NTX)
            YLS_FRC(M,NTX) = YLS_FRC(2,NTX)
            DO IGC = 1,NGC+2
              TMC_FRC(IGC,M,NTX) = TMC_FRC(IGC,2,NTX)
              ZMC_FRC(IGC,M,NTX) = ZMC_FRC(IGC,2,NTX)
            ENDDO
!
!---        Phase condition #1  ---
!
!           SL_FRC = 1.0, SG_FRC+SN_FRC = 0.0
!
!           energy - temperature, C - T_FRC
!           water mass - aqueous (system) pressure, Pa gauge - 
!             PL_FRC,PSO_FRC
!           CO2 mass - CO2 vapor pressure - PVA_FRC
!           CH4 mass - total CH4 mass - TMC_FRC
!           petroleum component mass - total petroleum comp. mass - 
!             TMC_FRC
!           salt mass - total salt brine mass fraction - YLS_FRC
!
            IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
              IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
                T_FRC(M,NTX) = T_FRC(M,NTX) + DNR_FRC(IEQT,NTX)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
              ELSEIF( M.EQ.IEQGC(1)+2 .AND. ISLC(45).EQ.0 ) THEN
                PVA_FRC(M,NTX) = PVA_FRC(M,NTX) + DNR_FRC(IEQGC(1),NTX)
              ELSEIF( M.EQ.IEQGC(2)+2 .AND. ISLC(64).EQ.0 ) THEN
                TMC_FRC(2,M,NTX) = TMC_FRC(2,M,NTX) + 
     &            DNR_FRC(IEQGC(2),NTX)
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                YLS_FRC(M,NTX) = YLS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ELSE
                DO IGC = 3,NGC+2
                  IF( M.EQ.IEQGC(IGC)+2 ) THEN
                    TMC_FRC(IGC,M,NTX) = TMC_FRC(IGC,M,NTX) + 
     &              DNR_FRC(IEQGC(IGC),NTX)
                  ENDIF
                ENDDO
              ENDIF
!
!---        Phase condition #2  ---
!
!           SL_FRC < 1.0, ST > 0.0, SG_FRC = 0.0
!
!           energy - temperature, C - T_FRC
!           water mass - aqueous saturation, SL_FRC
!           CO2 mass
!             if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!               - nonaqueous (system) press., Pa gauge - PN_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!             endif
!           CH4 mass
!             if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!               - nonaqueous (system) press., Pa gauge - PN_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!             endif
!           petroleum component mass
!             if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!               - nonaqueous (system) press., Pa gauge - PN_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!             endif
!           salt mass - total salt brine mass fraction - YLS_FRC
!
            ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
              IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
                T_FRC(M,NTX) = T_FRC(M,NTX) + DNR_FRC(IEQT,NTX)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
                ELSE
                  SL_FRC(M,NTX) = SL_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
                ENDIF
              ELSEIF( M.EQ.IEQGC(1)+2 .AND. ISLC(45).EQ.0 ) THEN
                IF( IZMC_FRC(NTX).EQ.1 ) THEN
                  IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                    SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(1),NTX)
                  ELSE
                    PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(1),NTX)
                  ENDIF
                ELSE
                  ZMC_FRC(1,M,NTX) = ZMC_FRC(1,M,NTX) + 
     &              DNR_FRC(IEQGC(1),NTX)
                ENDIF
              ELSEIF( M.EQ.IEQGC(2)+2 .AND. ISLC(64).EQ.0 ) THEN
                IF( IZMC_FRC(NTX).EQ.2 ) THEN
                  IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                    SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(2),NTX)
                  ELSE
                    PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(2),NTX)
                  ENDIF
                ELSE
                  ZMC_FRC(2,M,NTX) = ZMC_FRC(2,M,NTX) + 
     &              DNR_FRC(IEQGC(2),NTX)
                ENDIF
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                YLS_FRC(M,NTX) = YLS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ELSE
                DO IGC = 3,NGC+2
                  IF( M.EQ.IEQGC(IGC)+2 ) THEN
                    IF( IZMC_FRC(NTX).EQ.IGC ) THEN
                      IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                        SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                    DNR_FRC(IEQGC(IGC),NTX)
                      ELSE
                        PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                    DNR_FRC(IEQGC(IGC),NTX)
                      ENDIF
                    ELSE
                      ZMC_FRC(IGC,M,NTX) = ZMC_FRC(IGC,M,NTX) + 
     &                  DNR_FRC(IEQGC(IGC),NTX)
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
!
!---        Phase condition #3  ---
!
!           SL_FRC < 1.0, ST > 0.0, SG_FRC > 0.0
!
!           energy - temperature, C - T_FRC
!           water mass - aqueous saturation, SL_FRC
!           CO2 mass
!             if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!               - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!             endif
!           CH4 mass
!             if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!               - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!             endif
!           petroleum component mass
!             if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!               - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!             endif
!           salt mass - total salt brine mass fraction - YLS_FRC
!
            ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
              IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
                T_FRC(M,NTX) = T_FRC(M,NTX) + DNR_FRC(IEQT,NTX)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
                ELSE
                  SL_FRC(M,NTX) = SL_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
                ENDIF
              ELSEIF( M.EQ.IEQGC(1)+2 .AND. ISLC(45).EQ.0 ) THEN
                IF( IZMC_FRC(NTX).EQ.1 ) THEN
                  IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                    SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(1),NTX)
                  ELSE
                    PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(1),NTX)
                  ENDIF
                ELSE
                  ZMC_FRC(1,M,NTX) = ZMC_FRC(1,M,NTX) + 
     &              DNR_FRC(IEQGC(1),NTX)
                ENDIF
              ELSEIF( M.EQ.IEQGC(2)+2 .AND. ISLC(64).EQ.0 ) THEN
                IF( IZMC_FRC(NTX).EQ.2 ) THEN
                  IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                    SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(2),NTX)
                  ELSE
                    PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(2),NTX)
                  ENDIF
                ELSE
                  ZMC_FRC(2,M,NTX) = ZMC_FRC(2,M,NTX) + 
     &              DNR_FRC(IEQGC(2),NTX)
                ENDIF
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                YLS_FRC(M,NTX) = YLS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ELSE
                DO IGC = 3,NGC+2
                  IF( M.EQ.IEQGC(IGC)+2 ) THEN
                    IF( IZMC_FRC(NTX).EQ.IGC ) THEN
                      IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                        SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                    DNR_FRC(IEQGC(IGC),NTX)
                      ELSE
                        PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                    DNR_FRC(IEQGC(IGC),NTX)
                      ENDIF
                    ELSE
                      ZMC_FRC(IGC,M,NTX) = ZMC_FRC(IGC,M,NTX) + 
     &                  DNR_FRC(IEQGC(IGC),NTX)
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
!
!---        Phase condition #4  ---
!
!           SL_FRC = 0.0, ST = 1.0, SG_FRC = 0.0
!
!           energy - temperature, C - T_FRC
!           water mass - aqueous saturation, SL_FRC
!           CO2 mass
!             if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!               - nonaqueous (system) press., Pa gauge - PN_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!             endif
!           CH4 mass
!             if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!               - nonaqueous (system) press., Pa gauge - PN_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!             endif
!           petroleum component mass
!             if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!               - nonaqueous (system) press., Pa gauge - PN_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!             endif
!           salt mass - total salt mass - TMS_FRC(2,NTX)
!
            ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
              IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
                T_FRC(M,NTX) = T_FRC(M,NTX) + DNR_FRC(IEQT,NTX)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  SL_FRC(M,NTX) = SL_FRC(M,NTX) + DNR_FRC(IEQGC(1),NTX)
                ELSE
                  PSO_FRC(M,NTX) = PSO_FRC(M,NTX)+DNR_FRC(IEQGC(1),NTX)
                ENDIF
              ELSEIF( M.EQ.IEQGC(1)+2 .AND. ISLC(45).EQ.0 ) THEN
                IF( IZMC_FRC(NTX).EQ.1 ) THEN
                  IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                    SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(1),NTX)
                  ELSE
                    PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(1),NTX)
                  ENDIF
                ELSE
                  ZMC_FRC(1,M,NTX) = ZMC_FRC(1,M,NTX) + 
     &              DNR_FRC(IEQGC(1),NTX)
                ENDIF
              ELSEIF( M.EQ.IEQGC(2)+2 .AND. ISLC(64).EQ.0 ) THEN
                IF( IZMC_FRC(NTX).EQ.2 ) THEN
                  IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                    SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(2),NTX)
                  ELSE
                    PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(2),NTX)
                  ENDIF
                ELSE
                  ZMC_FRC(2,M,NTX) = ZMC_FRC(2,M,NTX) +  
     &              DNR_FRC(IEQGC(2),NTX)
                ENDIF
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                TMS_FRC(M,NTX) = TMS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ELSE
                DO IGC = 3,NGC+2
                  IF( M.EQ.IEQGC(IGC)+2 ) THEN
                    IF( IZMC_FRC(NTX).EQ.IGC ) THEN
                      IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                        SL_FRC(M,NT) = SL_FRC(M,NTX) + 
     &                    DNR_FRC(IEQGC(IGC),NTX)
                      ELSE
                        PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                    DNR_FRC(IEQGC(IGC),NTX)
                      ENDIF
                    ELSE
                      ZMC_FRC(IGC,M,NTX) = ZMC_FRC(IGC,M,NTX) + 
     &                  DNR_FRC(IEQGC(IGC),NTX)
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
!
!---        Phase condition #5  ---
!
!           SL_FRC = 0.0, ST >= 0.0, SG_FRC > 0.0
!
!           energy - temperature, C - T_FRC
!           water mass - aqueous saturation, SL_FRC
!           CO2 mass
!             if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!               - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!             endif
!           CH4 mass
!             if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!               - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!             endif
!           petroleum component mass
!             if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!               - gas (system) pressure, Pa gauge - PG_FRC,PSO_FRC
!             else
!               - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!             endif
!           salt mass - total salt mass - TMS_FRC(2,NTX)
!
            ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
              IF( M.EQ.IEQT+2 .AND. ISLC(30).EQ.0 ) THEN
                T_FRC(M,NTX) = T_FRC(M,NTX) + DNR_FRC(IEQT,NTX)
              ELSEIF( M.EQ.IEQW+2 ) THEN
                IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                  PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
                ELSE
                  SL_FRC(M,NTX) = SL_FRC(M,NTX) + DNR_FRC(IEQW,NTX)
                ENDIF
              ELSEIF( M.EQ.IEQGC(1)+2 .AND. ISLC(45).EQ.0 ) THEN
                IF( IZMC_FRC(NTX).EQ.1 ) THEN
                  IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                    SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(1),NTX)
                  ELSE
                    PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(1),NTX)
                  ENDIF
                ELSE
                  ZMC_FRC(1,M,NTX) = ZMC_FRC(1,M,NTX) +  
     &              DNR_FRC(IEQGC(1),NTX)
                ENDIF
              ELSEIF( M.EQ.IEQGC(2)+2 .AND. ISLC(64).EQ.0 ) THEN
                IF( IZMC_FRC(NTX).EQ.2 ) THEN
                  IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                    SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(2),NTX)
                  ELSE
                    PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                DNR_FRC(IEQGC(2),NTX)
                  ENDIF
                ELSE
                  ZMC_FRC(2,M,NTX) = ZMC_FRC(2,M,NTX) +  
     &              DNR_FRC(IEQGC(2),NTX)
                ENDIF
              ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
                TMS_FRC(M,NTX) = TMS_FRC(M,NTX) + DNR_FRC(IEQS,NTX)
              ELSE
                DO IGC = 3,NGC+2
                  IF( M.EQ.IEQGC(IGC)+2 ) THEN
                    IF( IZMC_FRC(NTX).EQ.IGC ) THEN
                      IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                        SL_FRC(M,NTX) = SL_FRC(M,NTX) + 
     &                    DNR_FRC(IEQGC(IGC),NTX)
                      ELSE
                        PSO_FRC(M,NTX) = PSO_FRC(M,NTX) + 
     &                    DNR_FRC(IEQGC(IGC),NTX)
                      ENDIF
                    ELSE
                      ZMC_FRC(IGC,M,NTX) = ZMC_FRC(IGC,M,NTX) + 
     &                  DNR_FRC(IEQGC(IGC),NTX)
                    ENDIF
                  ENDIF
                ENDDO
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
!---  End of INCRM_EOR group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBC_FRC_EOR( IGC )
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
!     Load the Jacobian matrix for the fault petroleum component
!     equation with aqueous, gas, and nonaqueous-liquid contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 17 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_FRC
      USE JACOB
      USE GLB_PAR
      USE GEOM_FRC
      USE FLUX_FRC
      USE FLUXC_FRC
      USE FDVP_FRC
      USE FDVGC_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STCX(LUK+1),RCPX(LUK),RCAX(LUK,LTC_FRC),FC(LSFV,LTC_FRC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBC_FRC_GT'
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
!---      First-order, forward-difference, time differential  ---
!
          DO M = 1,ISVC+1
            MP = M + 1
            STCX(M) = DTI*(TMC_FRC(IGC,MP,NT1X)-TMC_FRC(IGC,1,NT1X))
          ENDDO
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NTCX = NCX - IPF_FRC(1,NT1X) + 1
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Fault triangle surface fluxes  ---
!
            DO M = 1,ISVF
              M1 = MNOD(M)
              M2 = MADJ(M)
              INDX = -1
              APX = DIFMN( APM_FRC(M1,NT1X),APM_FRC(M2,NT2X),
     &          DFF1X,DFF2X,ZERO,INDX )
              IF( IGC.EQ.1 ) THEN
                FC(M,NTCX) = AFF_FRC(NCX)*APX*(UFFLA(M,NCX) 
     &            + UFFGC(IGC,M,NCX) + UFFNC(IGC,M,NCX))
              ELSE
                FC(M,NTCX) = AFF_FRC(NCX)*APX*(UFFGC(IGC,M,NCX)
     &            + UFFNC(IGC,M,NCX))
              ENDIF
            ENDDO
          ENDDO
!
!---      Compute fault triangle petroleum component
!         equation residuals  ---
!
          RCSX = STCX(1) - SRCGC_FRC(IGC,2,NT1X)
          NTCX = IPF_FRC(2,NT1X) - IPF_FRC(1,NT1X) + 1
          DO MD = 1,NTCX
            RCSX = RCSX + FC(1,MD)
          ENDDO
          DO M = 1,ISVC
            RCPX(M) = STCX(M+1) - SRCGC_FRC(IGC,M+2,NT1X)
            MM = 2*M
            DO MD = 1,NTCX
              RCPX(M) = RCPX(M) + FC(MM,MD)
            ENDDO
          ENDDO
          DO M = 1,ISVC
            MM = 2*M + 1
            DO MD = 1,NTCX
              RCAX(M,MD) = RCSX - FC(1,MD) + FC(MM,MD)
            ENDDO
          ENDDO
!
!---      Load Jacobian Matrix  ---
!
          CALL JCBL_FRC_EOR( RCSX,RCPX,RCAX,NT1X,IEQGC(IGC) )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBC_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBC_MF_EOR( IGC )
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
!     Modify Jacobian matrix for matrix grid cells and 
!     fault triangles for transfer of petroleum component
!     mass between matrix grid cells and fault triangles.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 17 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE FLUXC_FRC
      USE FDVP_FRC
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
      SUB_LOG(ISUB_LOG) = '/JCBC_MF_EOR'
!
!---  Matrix equations, fault connections,
!     loop over faults  ---
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
!
!---    Loop over fault triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fault triangle  ---
!
          N = INCM_FRC(NCX)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSC_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQGC(IGC),NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with 
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + (TRNSC_FRC(IGC,1,NCX) 
     &          - TRNSC_FRC(IGC,MX,NCX))/DNR(M,N)
!
!---          Partial derivative of water mass residual with 
!             respect to fault triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + (TRNSC_FRC(IGC,1,NCX)
     &          - TRNSC_FRC(IGC,MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSC_FRC(IGC,1,NCX)
            RSDL(IEQGC(IGC),N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSC_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQGC(IGC),NMD)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with 
!             respect to field node primary variables  ---
!
              MCOL = KLU(MP,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) + (TRNSC_FRC(IGC,1,NCX)
     &          - TRNSC_FRC(IGC,MX,NCX))/DNR(M,N)
!
!---          Partial derivative of water mass residual with 
!             respect to fault triangle primary variables  ---
!
              MC = (NCX-1)*ISVC + IEQGC(IGC)
              MCOL = KLU_MCF(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) + (TRNSC_FRC(IGC,1,NCX)
     &          - TRNSC_FRC(IGC,MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSC_FRC(IGC,1,NCX)
            RSDL(IEQGC(IGC),N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO  
!
!---  Fault equations, loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fault triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        MTX = IXP_FRC(NTX)
        IF( MTX.EQ.0 ) CYCLE
!
!---    Loop over fault triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fault triangle  ---
!
          N = INCM_FRC(NCX)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fault triangle, noting that TRNSC_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM_FRC(IEQGC(IGC),NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with 
!             respect to fault triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + (TRNSC_FRC(IGC,MX,NCX)
     &          - TRNSC_FRC(IGC,1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of water mass residual with 
!             respect to field node primary variables  ---
!
              MCOL = IM(M,NMD)
              MROW = MP-MCOL+MDC
              MX = 2*M
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + (TRNSC_FRC(IGC,MX,NCX)
     &          - TRNSC_FRC(IGC,1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSC_FRC(IGC,1,NCX)
            RSDL(IEQGC(IGC),N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fault triangle, noting that TRNSC_FRC is positive
!           from the fault triangle to the field node  ---
!
            MP = IM_FRC(IEQGC(IGC),NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with 
!             respect to fault triangle primary variables  ---
!
              NMD = (MTX-1)*ISVC + IEQGC(IGC)
              MCOL = KLU_FRC(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) + (TRNSC_FRC(IGC,MX,NCX)
     &          - TRNSC_FRC(IGC,1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of water mass residual with 
!             respect to field node primary variables  ---
!
              NMD = (NCX-1)*ISVC + IEQGC(IGC)
              MCOL = KLU_FCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) + (TRNSC_FRC(IGC,MX,NCX)
     &          - TRNSC_FRC(IGC,1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSC_FRC(IGC,1,NCX)
            RSDL(IEQGC(IGC),N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO  
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBC_MF_EOR group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBL_FRC_EOR( RSS,RSP,RSA,NT1X,MEQ )
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
!     Load the fault Jacobian matrix
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 17 April 2020.
!
#ifdef petsc
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
!
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE GLB_PAR
      USE GEOM_FRC
      USE FDVP_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#endif
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP(LUK),RSA(LUK,LTC_FRC)
#ifdef petsc
      REAL*8 BUFFER(LUK)
      INTEGER IROW(1),ICOL(LUK)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBL_FRC_EOR'
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
!---    Loop over fault triangle connections  ---
!
        DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
          NTCX = NCX - IPF_FRC(1,NT1X) + 1
          NT2X = ITCM_FRC(NCX)
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
!---    Loop over fault triangle connections  ---
!
        DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
          NTCX = NCX - IPF_FRC(1,NT1X) + 1
          NT2X = ITCM_FRC(NCX)
          MA = MA + ISVC
          DO M = 1,ISVC
            DNRX = DNR_FRC(M,NT2X)
            MCOL = KLU_FRC(MB,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,NTCX)-RSS)/DNRX
          ENDDO
        ENDDO
#ifdef petsc
!
!---  PETSc solver  ---
!
      ELSEIF( ILES.EQ.5 ) THEN
!
!---    Node  ---
!
        MP = IM_FRC(MEQ,NT1X)
        IROW(1) = MP-1
        DO M = 1,ISVC
          ICOL(M) = IM_FRC(M,NT1X)-1
          BUFFER(M) = (RSP(M)-RSS)/DNR_FRC(M,NT1X)
        ENDDO
        CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                    ADD_VALUES,IERR)
        BLU(MP) = BLU(MP) - RSS
        RSDL_FRC(MEQ,NT1X) = BLU(MP)
!
!---    Loop over fault triangle connections  ---
!
        DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
          NT2X = ITCM_FRC(NCX)
          DO M = 1,ISVC
            DNRX = DNR_FRC(M,NT2X)
            ICOL(M) = IM_FRC(M,NT1X)-1
            BUFFER(M) = (RSA(M,1)-RSS)/DNRX
          ENDDO
          CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                      ADD_VALUES,IERR)
        ENDDO
#endif
      ELSE
        INDX = 3
        CHMSG = 'Unknown Fault Linear Equation Solver'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBL_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_FRC_EOR
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
!     Load the Jacobian matrix for the fault salt equation with
!     aqueous-phase contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 17 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_FRC
      USE JACOB
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
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STSX(LUK+1),RSP(LUK),RSA(LUK,LTC_FRC),FS(LSFV,LTC_FRC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBS_FRC_EOR'
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
!---      First-order, forward-difference, time differential  ---
!
          DO M = 1,ISVC+1
            MP = M + 1
            STSX(M) = (TMS_FRC(MP,NT1X)-TMS_FRC(1,NT1X))*DTI
          ENDDO
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NTCX = NCX - IPF_FRC(1,NT1X) + 1
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Fault triangle surface fluxes  ---
!
            DO M = 1,ISVF
              M1 = MNOD(M)
              M2 = MADJ(M)
              INDX = -1
              APX = DIFMN( APM_FRC(M1,NT1X),APM_FRC(M2,NT2X),
     &          DFF1X,DFF2X,ZERO,INDX )
              FS(M,NTCX) = AFF_FRC(NCX)*APX*UFFS(M,NCX)
            ENDDO
          ENDDO
!
!---      Compute fault triangle salt equation residuals  ---
!
          RSS = STSX(1) - SRCS_FRC(2,NT1X)
          NTCX = IPF_FRC(2,NT1X) - IPF_FRC(1,NT1X) + 1
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
          CALL JCBL_FRC_EOR( RSS,RSP,RSA,NT1X,IEQS )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_MF_EOR
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
!     Modify Jacobian matrix for matrix grid cells and 
!     fault triangles for transfer of salt
!     mass between matrix grid cells and fault triangles.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 17 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
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
      SUB_LOG(ISUB_LOG) = '/JCBS_MF_EOR'
!
!---  Matrix equations, fault connections,
!     loop over faults  ---
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
!
!---    Loop over fault triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fault triangle  ---
!
          N = INCM_FRC(NCX)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSS_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQS,NMD)
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
!             respect to fault triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (TRNSS_FRC(1,NCX)-TRNSS_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSS_FRC(1,NCX)
            RSDL(IEQS,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSS_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQS,NMD)
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
!             respect to fault triangle primary variables  ---
!
              MC = (NCX-1)*ISVC + IEQS
              MCOL = KLU_MCF(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) + 
     &          (TRNSS_FRC(1,NCX)-TRNSS_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSS_FRC(1,NCX)
            RSDL(IEQS,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO  
!
!---  Fault equations, loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fault triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        MTX = IXP_FRC(NTX)
        IF( MTX.EQ.0 ) CYCLE
!
!---    Loop over fault triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fault triangle  ---
!
          N = INCM_FRC(NCX)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fault triangle, noting that TRNSS_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM_FRC(IEQS,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of salt mass residual with 
!             respect to fault triangle primary variables  ---
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
            RSDL(IEQS,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fault triangle, noting that TRNSS_FRC is positive
!           from the fault triangle to the field node  ---
!
            MP = IM_FRC(IEQS,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of salt mass residual with 
!             respect to fault triangle primary variables  ---
!
              NMD = (MTX-1)*ISVC + IEQS
              MCOL = KLU_FRC(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) + 
     &          (TRNSS_FRC(MX,NCX)-TRNSS_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of salt mass residual with 
!             respect to field node primary variables  ---
!
              NMD = (NCX-1)*ISVC + IEQS
              MCOL = KLU_FCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) + 
     &          (TRNSS_FRC(MX,NCX)-TRNSS_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSS_FRC(1,NCX)
            RSDL(IEQS,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO  
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_MF_EOR group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_FRC_EOR
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
!     Load the Jacobian matrix for the fault energy equation with
!     aqueous, gas, and nonaqueous-liquid contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 17 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_FRC
      USE JACOB
      USE GLB_PAR
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
      SUB_LOG(ISUB_LOG) = '/JCBT_FRC_EOR'
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
!---      First-order, forward-difference, time differential  ---
!
          STT1 = SS_FRC(1,NT1X)*RHOSP_FRC(1,NT1X)*HSP_FRC(1,NT1X) + 
     &      SL_FRC(1,NT1X)*RHOL_FRC(1,NT1X)*HL_FRC(1,NT1X) +
     &      SG_FRC(1,NT1X)*RHOG_FRC(1,NT1X)*UEG_FRC(1,NT1X)
          DO M = 1,ISVC+1
            MP = M + 1
            STT0 = SS_FRC(MP,NT1X)*RHOSP_FRC(MP,NT1X)*HSP_FRC(MP,NT1X) +
     &        SL_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)*HL_FRC(MP,NT1X) +
     &        SG_FRC(MP,NT1X)*RHOG_FRC(MP,NT1X)*UEG_FRC(MP,NT1X)
            STTX(M) = (STT0-STT1)*VOL_FRC(NT1X)*DTI
          ENDDO
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NTCX = NCX - IPF_FRC(1,NT1X) + 1
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Fault triangle surface fluxes  ---
!
            DO M = 1,ISVF
              INDX = -1
              M1 = MNOD(M)
              M2 = MADJ(M)
              APX = DIFMN( APM_FRC(M1,NT1X),APM_FRC(M2,NT2X),
     &          DFF1X,DFF2X,ZERO,INDX )
              FT(M,NTCX) = AFF_FRC(NCX)*APX*UFFQ(M,NCX)
            ENDDO
          ENDDO
!
!---      Compute fault triangle energy equation residuals  ---
!
          RTS = STTX(1) - SRCT_FRC(2,NT1X)
          NTCX = IPF_FRC(2,NT1X) - IPF_FRC(1,NT1X) + 1
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
          CALL JCBL_FRC_EOR( RTS,RTP,RTA,NT1X,IEQT )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_MF_EOR
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
!     Modify Jacobian matrix for matrix grid cells and 
!     fault triangles for transfer of energy
!     between matrix grid cells and fault triangles.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 17 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
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
      SUB_LOG(ISUB_LOG) = '/JCBT_MF_EOR'
!
!---  Matrix equations, fault connections,
!     loop over faults  ---
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
!
!---    Loop over fault triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fault triangle  ---
!
          N = INCM_FRC(NCX)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSQ_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQT,NMD)
            DO M = 1,ISVC
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
!             respect to fault triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (TRNSQ_FRC(1,NCX)-TRNSQ_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSQ_FRC(1,NCX)
            RSDL(IEQT,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSQ_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQT,NMD)
            DO M = 1,ISVC
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
!             respect to fault triangle primary variables  ---
!
              MC = (NCX-1)*ISVC + IEQT
              MCOL = KLU_MCF(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) + 
     &          (TRNSQ_FRC(1,NCX)-TRNSQ_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSQ_FRC(1,NCX)
            RSDL(IEQT,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO  
!
!---  Fault equations, loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fault triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        MTX = IXP_FRC(NTX)
        IF( MTX.EQ.0 ) CYCLE
!
!---    Loop over fault triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fault triangle  ---
!
          N = INCM_FRC(NCX)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fault triangle, noting that TRNSQ_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM_FRC(IEQT,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of energy residual with 
!             respect to fault triangle primary variables  ---
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
            RSDL(IEQT,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fault triangle, noting that TRNSQ_FRC is positive
!           from the fault triangle to the field node  ---
!
            MP = IM_FRC(IEQT,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of energy residual with 
!             respect to fault triangle primary variables  ---
!
              NMD = (MTX-1)*ISVC + IEQT
              MCOL = KLU_FRC(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) + 
     &          (TRNSQ_FRC(MX,NCX)-TRNSQ_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of energy residual with 
!             respect to field node primary variables  ---
!
              NMD = (NCX-1)*ISVC + IEQT
              MCOL = KLU_FCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) + 
     &          (TRNSQ_FRC(MX,NCX)-TRNSQ_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSQ_FRC(1,NCX)
            RSDL(IEQT,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO  
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_MF_EOR group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_FRC_EOR
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
!     Load the Jacobian matrix for the fault water equation with
!     aqueous and gas contributions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 17 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE PARM_FRC
      USE JACOB
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
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STWX(LUK+1),RWP(LUK),RWA(LUK,LTC_FRC),FW(LSFV,LTC_FRC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBW_FRC_EOR'
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
!---      First-order, forward-difference, time differential  ---
!
          STW1 = (XLW_FRC(1,NT1X)*RHOL_FRC(1,NT1X)*SL_FRC(1,NT1X) +
     &        XGW_FRC(1,NT1X)*RHOG_FRC(1,NT1X)*SG_FRC(1,NT1X))
          DO M = 1,ISVC+1
            MP = M + 1
            STW0 = (XLW_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)*SL_FRC(MP,NT1X) +
     &        XGW_FRC(MP,NT1X)*RHOG_FRC(MP,NT1X)*SG_FRC(MP,NT1X))
            STWX(M) = (STW0-STW1)*DTI*VOL_FRC(NT1X)
          ENDDO
!
!---      Loop over fault triangle connections  ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NTCX = NCX - IPF_FRC(1,NT1X) + 1
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
!
!---        Fault triangle surface fluxes  ---
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
              INDX = -1
              APX = DIFMN( APM_FRC(M1,NT1X),APM_FRC(M2,NT2X),
     &          DFF1X,DFF2X,ZERO,INDX )
              FW(M,NTCX) = AFF_FRC(NCX)*APX*
     &          (UFFL(M,NCX)*FLW + UFFG(M,NCX)*FGW
     &          + WTMW*UFFDGW(M,NCX) + WTMW*UFFDLA(M,NCX)
     &          - WTMW*UFFDS(M,NCX))
            ENDDO
          ENDDO
!
!---      Compute fault triangle water equation residuals  ---
!
          RWS = STWX(1) - SRCW_FRC(2,NT1X)
          NTCX = IPF_FRC(2,NT1X) - IPF_FRC(1,NT1X) + 1
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
          CALL JCBL_FRC_EOR( RWS,RWP,RWA,NT1X,IEQW )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_MF_EOR
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
!     Modify Jacobian matrix for matrix grid cells and 
!     fault triangles for transfer of water
!     mass between matrix grid cells and fault triangles.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 17 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVP_FRC
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
      SUB_LOG(ISUB_LOG) = '/JCBW_MF_EOR'
!
!---  Matrix equations, fault connections,
!     loop over faults  ---
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
!
!---    Loop over fault triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fault triangle  ---
!
          N = INCM_FRC(NCX)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Field node, noting that TRNSW_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQW,NMD)
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
!             respect to fault triangle primary variables  ---
!
              MCOL = IM_FRC(M,NTX)
              MROW = MP-MCOL+MDC
              MX = 2*M + 1
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + 
     &          (TRNSW_FRC(1,NCX)-TRNSW_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSW_FRC(1,NCX)
            RSDL(IEQW,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Field node, noting that TRNSW_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM(IEQW,NMD)
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
!             respect to fault triangle primary variables  ---
!
              MC = (NCX-1)*ISVC + IEQW
              MCOL = KLU_MCF(MC,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) + 
     &          (TRNSW_FRC(1,NCX)-TRNSW_FRC(MX,NCX))/DNR_FRC(M,NTX)
            ENDDO
            BLU(MP) = BLU(MP) + TRNSW_FRC(1,NCX)
            RSDL(IEQW,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO  
!
!---  Fault equations, loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---  Loop over fault triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        MTX = IXP_FRC(NTX)
        IF( MTX.EQ.0 ) CYCLE
!
!---    Loop over fault triangle to grid cell connections  ---
!
        DO NCX = IPN_FRC(1,NTX),IPN_FRC(2,NTX)
!
!---      Field node connected to fault triangle  ---
!
          N = INCM_FRC(NCX)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
!
!---        Fault triangle, noting that TRNSW_FRC is positive
!           from the fault triangle to the field node  ---
!
            NMD = IXP(N)
            MP = IM_FRC(IEQW,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with 
!             respect to fault triangle primary variables  ---
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
            RSDL(IEQW,N) = BLU(MP)
!
!---      SPLIB or Lis solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Fault triangle, noting that TRNSW_FRC is positive
!           from the fault triangle to the field node  ---
!
            MP = IM_FRC(IEQW,NTX)
            DO M = 1,ISVC
!
!---          Partial derivative of water mass residual with 
!             respect to fault triangle primary variables  ---
!
              NMD = (MTX-1)*ISVC + IEQW
              MCOL = KLU_FRC(NMD,M)
              MX = 2*M + 1
              DLU(MCOL) = DLU(MCOL) + 
     &          (TRNSW_FRC(MX,NCX)-TRNSW_FRC(1,NCX))/DNR_FRC(M,NTX)
!
!---          Partial derivative of water mass residual with 
!             respect to field node primary variables  ---
!
              NMD = (NCX-1)*ISVC + IEQW
              MCOL = KLU_FCM(NMD,M)
              MX = 2*M
              DLU(MCOL) = DLU(MCOL) + 
     &          (TRNSW_FRC(MX,NCX)-TRNSW_FRC(1,NCX))/DNR(M,N)
            ENDDO
            BLU(MP) = BLU(MP) - TRNSW_FRC(1,NCX)
            RSDL(IEQW,N) = BLU(MP)
          ENDIF
        ENDDO
      ENDDO
      ENDDO  
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_MF_EOR group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_FRC_EOR
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
!     Load the current time step values into the old time step
!     variables for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Ramesh Sarathi, PNNL, 6 March 2020.
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
      USE FDVN_FRC
      USE FDVG_FRC
      USE FDVGC_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LDO_FRC_EOR'
!
!---  Assign old time step values, loop over faults  ---
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
        DO IGC = 1,NGC+2
          DFGC_FRC(IGC,1,NTX) = DFGC_FRC(IGC,2,NTX)
          DFNC_FRC(IGC,1,NTX) = DFNC_FRC(IGC,2,NTX)
          TMC_FRC(IGC,1,NTX) = TMC_FRC(IGC,2,NTX)
          XGC_FRC(IGC,1,NTX) = XGC_FRC(IGC,2,NTX)
          XMGC_FRC(IGC,1,NTX) = XMGC_FRC(IGC,2,NTX)
          XLC_FRC(IGC,1,NTX) = XLC_FRC(IGC,2,NTX)
          XMLC_FRC(IGC,1,NTX) = XMLC_FRC(IGC,2,NTX)
          XMNC_FRC(IGC,1,NTX) = XMNC_FRC(IGC,2,NTX)
          XNC_FRC(IGC,1,NTX) = XNC_FRC(IGC,2,NTX)
          ZMC_FRC(IGC,1,NTX) = ZMC_FRC(IGC,2,NTX)
        ENDDO
        ASLMIN_FRC(1,NTX) = MIN( ASL_FRC(NTX),ASLMIN_FRC(2,NTX) )
        ASLMIN_FRC(2,NTX) = ASLMIN_FRC(1,NTX)
        DO M = 1,6
          BETA_FRC(M,NTX) = 0.D+0
        ENDDO
        IBETA_FRC(NTX) = 0
        DFGW_FRC(1,NTX) = DFGW_FRC(2,NTX)
        DFLA_FRC(1,NTX) = DFLA_FRC(2,NTX)
        DFLS_FRC(1,NTX) = DFLS_FRC(2,NTX)
        GNIFT_FRC(1,NTX) = GNIFT_FRC(2,NTX)
!
!---    Isothermal option  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          HG_FRC(1,NTX) = HG_FRC(2,NTX)
          HGA_FRC(1,NTX) = HGA_FRC(2,NTX)
          HGW_FRC(1,NTX) = HGW_FRC(2,NTX)
          HL_FRC(1,NTX) = HL_FRC(2,NTX)
          HLW_FRC(1,NTX) = HLW_FRC(2,NTX)
          HN_FRC(1,NTX) = HN_FRC(2,NTX)
          HSP_FRC(1,NTX) = HSP_FRC(2,NTX)
          UEG_FRC(1,NTX) = UEG_FRC(2,NTX)
          UEL_FRC(1,NTX) = UEL_FRC(2,NTX)
          UEN_FRC(1,NTX) = UEN_FRC(2,NTX)
          THKG_FRC(1,NTX) = THKG_FRC(2,NTX)
          THKL_FRC(1,NTX) = THKL_FRC(2,NTX)
          THKN_FRC(1,NTX) = THKN_FRC(2,NTX)
        ENDIF
        NPHAZ_FRC(1,NTX) = NPHAZ_FRC(2,NTX)
        PERMRF_FRC(1,NTX) = PERMRF_FRC(2,NTX)
        PG_FRC(1,NTX) = PG_FRC(2,NTX)
        PL_FRC(1,NTX) = PL_FRC(2,NTX)
        PN_FRC(1,NTX) = PN_FRC(2,NTX)
        POSM_FRC(1,NTX) = POSM_FRC(2,NTX)
        PSO_FRC(1,NTX) = PSO_FRC(2,NTX)
        PORD_FRC(1,NTX) = PORD_FRC(2,NTX)
        PORT_FRC(1,NTX) = PORT_FRC(2,NTX)
        PVA_FRC(1,NTX) = PVA_FRC(2,NTX)
        PVW_FRC(1,NTX) = PVW_FRC(2,NTX)
        RHOG_FRC(1,NTX) = RHOG_FRC(2,NTX)
        RHOL_FRC(1,NTX) = RHOL_FRC(2,NTX)
        RHOMG_FRC(1,NTX) = RHOMG_FRC(2,NTX)
        RHOML_FRC(1,NTX) = RHOML_FRC(2,NTX)
        RHOMN_FRC(1,NTX) = RHOMN_FRC(2,NTX)
        RHON_FRC(1,NTX) = RHON_FRC(2,NTX)
        RHOSP_FRC(1,NTX) = RHOSP_FRC(2,NTX)
        RKG_FRC(1,NTX) = RKG_FRC(2,NTX)
        RKL_FRC(1,NTX) = RKL_FRC(2,NTX)
        RKN_FRC(1,NTX) = RKN_FRC(2,NTX)
        SG_FRC(1,NTX) = SG_FRC(2,NTX)
        SL_FRC(1,NTX) = SL_FRC(2,NTX)
        SN_FRC(1,NTX) = SN_FRC(2,NTX)
        SS_FRC(1,NTX) = SS_FRC(2,NTX)
        T_FRC(1,NTX) = T_FRC(2,NTX)
        TMS_FRC(1,NTX) = TMS_FRC(2,NTX)
        TORG_FRC(1,NTX) = TORG_FRC(2,NTX)
        TORL_FRC(1,NTX) = TORL_FRC(2,NTX)
        TORN_FRC(1,NTX) = TORN_FRC(2,NTX)
        VISG_FRC(1,NTX) = VISG_FRC(2,NTX)
        VISL_FRC(1,NTX) = VISL_FRC(2,NTX)
        VISN_FRC(1,NTX) = VISN_FRC(2,NTX)
        XGW_FRC(1,NTX) = XGW_FRC(2,NTX)
        XLA_FRC(1,NTX) = XLA_FRC(2,NTX)
        XLS_FRC(1,NTX) = XLS_FRC(2,NTX)
        XLW_FRC(1,NTX) = XLW_FRC(2,NTX)
        XMGW_FRC(1,NTX) = XMGW_FRC(2,NTX)
        XMLA_FRC(1,NTX) = XMLA_FRC(2,NTX)
        XMLS_FRC(1,NTX) = XMLS_FRC(2,NTX)
        XMLW_FRC(1,NTX) = XMLW_FRC(2,NTX)
        YLS_FRC(1,NTX) = YLS_FRC(2,NTX)
        ZG_FRC(1,NTX) = ZG_FRC(2,NTX)
        ZN_FRC(1,NTX) = ZN_FRC(2,NTX)
        IBETA_FRC(NTX) = 0
        DO NSL = 1,NSOLU
          CO_FRC(NTX,NSL) = C_FRC(NTX,NSL)
        ENDDO
#ifdef ecke
        IF( ISLC(40).EQ.1 ) THEN
          DO NEQ = 1,NEQC+NEQK
            NSL = NEQ + NSOLU
            CO_FRC(NTX,NSL) = C_FRC(NTX,NSL)
          ENDDO
          DO NSP = 1,NSPR
            SP_CO_FRC(NTX,NSP) = SP_C_FRC(NTX,NSP)
          ENDDO
        ENDIF
#endif
!
!---  End loop over fault triangles  ---
!
      ENDDO
!
!---  End loop over faults  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_FRC_EOR group
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_FRC_BO
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
!     Compute hydrologic, thermodynamic and physical properties
!     for the black-oil option
!
!----------------------Authors-----------------------------------------!
!
!     Written by Ramesh Sarathi, PNNL, 6 March 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE FDVS_FRC
      USE FDVP_FRC
      USE FDVN_FRC
      USE FDVG_FRC
      USE FDVGC_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP_FRC_BO'
!
!---  Black-oil option, loop over faults  ---
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
        N_DB = NTX
        IZN = IZ_FRC(NTX)
!
!---    Loop over increment indices  ---
!
        DO M = 2,ISVC+2
          TKX = T_FRC(M,NTX) + TABS
!
!---      Phase condition #1  ---
!
!         SL = 1.0, SG+SN = 0.0
!
!         energy - temperature, C - T
!         water mass - aqueous (system) pressure, Pa gauge - PL,PSO
!         CO2 mass - CO2 vapor pressure - PVA
!         CH4 mass - total CH4 mass - TMC
!         petroleum component mass - total petroleum comp. mass - TMC
!         salt mass - total salt brine mass fraction - YLS
!
          IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(M,NTX),
     &        PORT_FRC(M,NTX) )
            PORD_FRC(M,NTX) = MAX( PORD_FRC(M,NTX),EPSL )
            PORT_FRC(M,NTX) = MAX( PORT_FRC(M,NTX),PORD_FRC(M,NTX) )
!
!---        Flash calculations  ---
!
            CALL BO_PC1( PORD_FRC(M,NTX),POSM_FRC(M,NTX),PVW_FRC(M,NTX),
     &        PX,RHOG_FRC(M,NTX),RHOL_FRC(M,NTX),RHON_FRC(M,NTX),
     &        SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),TMS_FRC(M,NTX),
     &        T_FRC(M,NTX),VOL_FRC(NTX),XGC_FRC(1,M,NTX),XGW_FRC(M,NTX),
     &        XLS_FRC(M,NTX),XLW_FRC(M,NTX),XNC_FRC(1,M,NTX),
     &        YLS_FRC(M,NTX) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR( SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(M,NTX) = PSO_FRC(M,NTX)
            PTX = PG_FRC(M,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(M,NTX) = PTX
              PL_FRC(M,NTX) = PN_FRC(M,NTX) - CPNLX
            ELSE
              PL_FRC(M,NTX) = PTX
              PN_FRC(M,NTX) = PL_FRC(M,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #2  ---
!
!         SG = 0.0
!         SL > 0.0
!         SN > 0.0
!         Aqueous and nonaqueous liquid saturated
!
!         water mass - aqueous saturation, SL
!         gas mass - gas to oil ratio, POSM
!         oil mass - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         salt mass - total salt brine mass fraction - YLS
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(M,NTX),
     &        PORT_FRC(M,NTX) )
            PORD_FRC(M,NTX) = MAX( PORD_FRC(M,NTX),EPSL )
            PORT_FRC(M,NTX) = MAX( PORT_FRC(M,NTX),PORD_FRC(M,NTX) )
!
!---        Flash calculations  ---
!
            CALL BO_PC2( PORD_FRC(M,NTX),POSM_FRC(M,NTX),PVW_FRC(M,NTX),
     &        PX,RHOG_FRC(M,NTX),RHOL_FRC(M,NTX),RHON_FRC(M,NTX),
     &        SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        TMC_FRC(1,M,NTX),TMS_FRC(M,NTX),T_FRC(M,NTX),VOL_FRC(NTX),
     &        XGC_FRC(1,M,NTX),XGW_FRC(M,NTX),XLS_FRC(M,NTX),
     &        XLW_FRC(M,NTX),XNC_FRC(1,M,NTX),YLS_FRC(M,NTX),IZN )
            CALL CAP_EOR( SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(M,NTX) = PSO_FRC(M,NTX)
            PTX = PG_FRC(M,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(M,NTX) = PTX
              PL_FRC(M,NTX) = PN_FRC(M,NTX) - CPNLX
            ELSE
              PL_FRC(M,NTX) = PTX
              PN_FRC(M,NTX) = PL_FRC(M,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #3  ---
!
!         SG > 0.0
!         SL > 0.0
!         SN > 0.0
!         Three-phase
!
!         water mass - aqueous saturation, SL
!         gas mass - gas saturation, SG
!         oil mass - gas (system) pressure, Pa gauge - PG,PSO
!         salt mass - total salt brine mass fraction - YLS
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(M,NTX),
     &        PORT_FRC(M,NTX) )
            PORD_FRC(M,NTX) = MAX( PORD_FRC(M,NTX),EPSL )
            PORT_FRC(M,NTX) = MAX( PORT_FRC(M,NTX),PORD_FRC(M,NTX) )
!
!---        Flash calculations  ---
!
            CALL BO_PC3( PORD_FRC(M,NTX),POSM_FRC(M,NTX),PVW_FRC(M,NTX),
     &        PX,RHOG_FRC(M,NTX),RHOL_FRC(M,NTX),RHON_FRC(M,NTX),
     &        SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        TMC_FRC(1,M,NTX),TMS_FRC(M,NTX),T_FRC(M,NTX),
     &        VOL_FRC(NTX),XGC_FRC(1,M,NTX),XGW_FRC(M,NTX),
     &        XLS_FRC(M,NTX),XLW_FRC(M,NTX),XNC_FRC(1,M,NTX),
     &        YLS_FRC(M,NTX),IZN )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR( SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(M,NTX) = PSO_FRC(M,NTX)
            PTX = PG_FRC(M,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(M,NTX) = PTX
              PL_FRC(M,NTX) = PN_FRC(M,NTX) - CPNLX
            ELSE
              PL_FRC(M,NTX) = PTX
              PN_FRC(M,NTX) = PL_FRC(M,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #4  ---
!
!         SG = 0.0
!         SL = 0.0
!         SN = 1.0
!         Nonaqueous liquid saturated
!
!         water mass - aqueous saturation, SL
!         gas mass - gas to oil ratio, POSM
!         oil mass - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         salt mass - total salt mass - TMS(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(M,NTX),
     &        PORT_FRC(M,NTX) )
            PORD_FRC(M,NTX) = MAX( PORD_FRC(M,NTX),EPSL )
            PORT_FRC(M,NTX) = MAX( PORT_FRC(M,NTX),PORD_FRC(M,NTX) )
!
!---        Flash calculations  ---
!
            CALL BO_PC4( PORD_FRC(M,NTX),POSM_FRC(M,NTX),PVW_FRC(M,NTX),
     &        PX,RHOG_FRC(M,NTX),RHOL_FRC(M,NTX),RHON_FRC(M,NTX),
     &        SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        TMC_FRC(1,M,NTX),TMS_FRC(M,NTX),T_FRC(M,NTX),VOL_FRC(NTX),
     &        XGC_FRC(1,M,NTX),XGW_FRC(M,NTX),XLS_FRC(M,NTX),
     &        XLW_FRC(M,NTX),XNC_FRC(1,M,NTX),YLS_FRC(M,NTX),IZN )
            CALL CAP_EOR( SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(M,NTX) = PSO_FRC(M,NTX)
            PTX = PG_FRC(M,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(M,NTX) = PTX
              PL_FRC(M,NTX) = PN_FRC(M,NTX) - CPNLX
            ELSE
              PL_FRC(M,NTX) = PTX
              PN_FRC(M,NTX) = PL_FRC(M,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #5  ---
!
!         SG > 0.0
!         SL = 0.0
!         SN > 0.0
!         Nonaqueous liquid and gas saturated
!
!         water mass - aqueous saturation, SL
!         gas mass - gas saturation, SG
!         oil mass - gas (system) pressure, Pa gauge - PG,PSO
!         salt mass - total salt mass - TMS(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(M,NTX),
     &        PORT_FRC(M,NTX) )
            PORD_FRC(M,NTX) = MAX( PORD_FRC(M,NTX),EPSL )
            PORT_FRC(M,NTX) = MAX( PORT_FRC(M,NTX),PORD_FRC(M,NTX) )
!
!---        Flash calculations  ---
!
            CALL BO_PC5( PORD_FRC(M,NTX),POSM_FRC(M,NTX),PVW_FRC(M,NTX),
     &        PX,RHOG_FRC(M,NTX),RHOL_FRC(M,NTX),RHON_FRC(M,NTX),
     &        SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        TMC_FRC(1,M,NTX),TMS_FRC(M,NTX),T_FRC(M,NTX),VOL_FRC(NTX),
     &        XGC_FRC(1,M,NTX),XGW_FRC(M,NTX),XLS_FRC(M,NTX),
     &        XLW_FRC(M,NTX),XNC_FRC(1,M,NTX),YLS_FRC(M,NTX),IZN )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR( SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(M,NTX) = PSO_FRC(M,NTX)
            PTX = PG_FRC(M,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(M,NTX) = PTX
              PL_FRC(M,NTX) = PN_FRC(M,NTX) - CPNLX
            ELSE
              PL_FRC(M,NTX) = PTX
              PN_FRC(M,NTX) = PL_FRC(M,NTX) + CPNLX
            ENDIF
          ENDIF
!
!---      Relative permeability  ---
!
#ifdef dual_porosity
          CALL RKG_EOR_F( RKG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &      IZN )
          CALL RKL_EOR_F( RKLX,SL_FRC(M,NTX),IZN )
          RKL_FRC(M,NTX) = RKLX(1)
          CALL RKN_EOR_F( RKG_FRC(M,NTX),RKLX,RKN_FRC(M,NTX),
     &      SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),IZN )
#else
          CALL RKG_EOR( RKG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),IZN )
          CALL RKL_EOR( RKLX,SL_FRC(M,NTX),IZN )
          RKL_FRC(M,NTX) = RKLX(1)
          CALL RKN_EOR( RKG_FRC(M,NTX),RKLX,RKN_FRC(M,NTX),
     &      SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),IZN )
#endif
!
!---      Gas viscosity, Pa s  ---
!
          CALL BO_VISG( RHOG_FRC(M,NTX),VISG_FRC(M,NTX) )
!
!---      Aqueous viscosity, Pa s
!         concentration of CO2 in mol/kg,
!         concentration of salt in mol/kg  ---
!
          CMLAX = 0.D+0
          CMLSX = 1.D+3*XLS_FRC(M,NTX)*RHOL_FRC(M,NTX)/WTMS
          CALL VISC_L( CMLAX,CMLSX,PX,T_FRC(M,NTX),VISL_FRC(M,NTX) )
!
!---      Nonaqueous-liquid viscosity, Pa s  ---
!
          CALL BO_VISN( PX,POSM_FRC(M,NTX),T_FRC(M,NTX),VISN_FRC(M,NTX))
!
!---      Salt aqueous diffusion coefficient  ---
!
          CALL DIFC_LS( DFLS_FRC(M,NTX),T_FRC(M,NTX),VISL_FRC(M,NTX),
     &      XLS_FRC(M,NTX) )
!
!---      Precipitated NaCl density, saturation, and enthalpy  ---
!
          CALL DENS_S( T_FRC(M,NTX),PX,RHOSP_FRC(M,NTX) )
          CALL SOL_LS( T_FRC(M,NTX),XBSMX )
          CALL DENS_B( T_FRC(M,NTX),PX,XBSMX,RHOBX )
          SS_FRC(M,NTX) = MAX(YLS_FRC(M,NTX)-XBSMX,0.D+0)*RHOBX*
     &      SL_FRC(M,NTX)/RHOSP_FRC(M,NTX)
!
!---      Kozney-Carmen permeability reduction  ---
!
          PORDX = PORD_FRC(M,NTX)*(1.D+0-SS_FRC(M,NTX))
          PERMRF_FRC(M,NTX) = ((PORDX**3)/(1.D+0+PORDX))/
     &      ((PORD_FRC(M,NTX)**3)/(1.D+0+PORD_FRC(M,NTX)))
!
!---      Change in porosity as a function of precipitated salt  ---
!
          PORD_FRC(M,NTX) = PORD_FRC(M,NTX)*(1.D+0-SS_FRC(M,NTX))
!
!---      Aqueous, gas, and nonaqueous-liquid tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL_FRC(M,NTX),
     &      SG_FRC(M,NTX),SN_FRC(M,NTX),PORD_FRC(M,NTX),TORL_FRC(M,NTX),
     &      TORG_FRC(M,NTX),TORN_FRC(M,NTX) )
        ENDDO
      ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_FRC_BO group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_FRC_EOR
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
!     Compute hydrologic, thermodynamic and physical properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, March 17, 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FDVT_FRC
      USE FDVS_FRC
      USE FDVP_FRC
      USE FDVN_FRC
      USE FDVG_FRC
      USE FDVGC_FRC
      USE FDVGC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 FKX(LNGC)
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP_FRC_EOR'
!
!---  Black-oil option  ---
!
      IF( ISLC(68).EQ.0 ) THEN
        CALL PROP_BO
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compositional option, loop over faults  ---
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
        N_DB = NTX
        IZN = IZ_FRC(NTX)
!
!---    K-factors fixed at second iterate values  ---
!
        IF( NITER.GT.2 ) THEN
          IFK = 1
!
!---    K-factors computed from fugacity coefficients  ---
!
        ELSE
          IFK = 0
        ENDIF
!
!---    Schema for fixing beta  ---
!
        IBETAX = 0
        IF( NITER.GT.7 ) THEN
          DBETA1X = BETA_FRC(2,NTX)-BETA_FRC(1,NTX)
          DBETA2X = BETA_FRC(3,NTX)-BETA_FRC(2,NTX)
          ABETAX = 0.D+0
          DO M = 1,3
            ABETAX = ABETAX + BETA_FRC(M,NTX)
          ENDDO
          ABETAX = ABETAX/3.D+0
          SBETAX = 0.D+0
          DO M = 1,3
            SBETAX = SBETAX + (BETA_FRC(M,NTX)-ABETAX)**2
          ENDDO
          SBETAX = SQRT( SBETAX/3.D+0 )
!
!---      Beta previously fixed  ---
!
          IF( IBETA_FRC(NTX).EQ.1 ) THEN
            IBETAX = 1
            BETAX = BETA_FRC(3,NTX)
!
!---      Small variance in beta  ---
!
          ELSEIF( SBETAX.LT.1.D-4 ) THEN
            BETAX = ABETAX
            IBETA_FRC(NTX) = 1
            IBETAX = 1
!
!---      Beta oscillating  ---
!
          ELSEIF( (DBETA1X*DBETA2X).LT.0.D+0 .AND. 
     &      (DBETA1X+DBETA2X).LT.1.D-4 ) THEN
            BETAX = 5.D-1*(BETA_FRC(3,NTX)+BETA_FRC(2,NTX))
            IBETA_FRC(NTX) = 1
            IBETAX = 1
!
!---      Gas-nonaqueous liquid interfacial tension oscillating  ---
!
          ELSEIF( (BETA_FRC(1,NTX)*BETA_FRC(5,NTX)*BETA_FRC(6,NTX))
     &      .LT.EPSL .AND. (BETA_FRC(4,NTX)+BETA_FRC(5,NTX)+
     &      BETA_FRC(6,NTX)).GT.EPSL ) THEN
            BETAX = 0.D+0
            IBETAX = 0
            DO M = 1,3
              IF( BETA_FRC(M+3,NTX).GT.0.D+0 ) THEN
                BETAX = BETAX + BETA_FRC(M,NTX)
                IBETAX = IBETAX + 1
              ENDIF
            ENDDO
            BETAX = BETAX/REAL(IBETAX)
            IBETA_FRC(NTX) = 1
            IBETAX = 1
          ENDIF
        ENDIF
!
!---    Loop over increment indices  ---
!
        DO M = 2,ISVC+2
          TKX = T_FRC(M,NTX) + TABS
!
!---      K-factors fixed at third iterate values  ---
!
          IF( IFK.EQ.1 ) THEN
            DO IGC = 1,NGC+2
              FKX(IGC) = FK_FRC(IGC,NTX)
            ENDDO
          ENDIF
!
!---      Phase condition #1  ---
!
!         SL = 1.0, SG+SN = 0.0
!
!         energy - temperature, C - T
!         water mass - aqueous (system) pressure, Pa gauge - PL,PSO
!         CO2 mass - CO2 vapor pressure - PVA
!         CH4 mass - total CH4 mass - TMC
!         petroleum component mass - total petroleum comp. mass - TMC
!         salt mass - total salt brine mass fraction - YLS
!
          IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(M,NTX),
     &        PORT_FRC(M,NTX) )
            PORD_FRC(M,NTX) = MAX( PORD_FRC(M,NTX),EPSL )
            PORT_FRC(M,NTX) = MAX( PORT_FRC(M,NTX),PORD_FRC(M,NTX) )
!
!---        Flash calculations  ---
!
            CALL FLH_PC1( BETAX,FKX,GNIFT_FRC(M,NTX),PORD_FRC(M,NTX),
     &        PVA_FRC(M,NTX),PVW_FRC(M,NTX),PX,RHOG_FRC(M,NTX),
     &        RHOL_FRC(M,NTX),RHOMG_FRC(M,NTX),RHOML_FRC(M,NTX),
     &        RHOMN_FRC(M,NTX),RHON_FRC(M,NTX),SG_FRC(M,NTX),
     &        SL_FRC(M,NTX),SN_FRC(M,NTX),TMC_FRC(1,M,NTX),
     &        TMS_FRC(M,NTX),T_FRC(M,NTX),VOL_FRC(NTX),XGC_FRC(1,M,NTX),
     &        XGW_FRC(M,NTX),XLA_FRC(M,NTX),XLS_FRC(M,NTX),
     &        XLW_FRC(M,NTX),XNC_FRC(1,M,NTX),XMGC_FRC(1,M,NTX),
     &        XMGW_FRC(M,NTX),XMLA_FRC(M,NTX),XMLS_FRC(M,NTX),
     &        XMLW_FRC(M,NTX),XMNC_FRC(1,M,NTX),YLS_FRC(M,NTX),
     &        ZG_FRC(M,NTX),ZMC_FRC(1,M,NTX),ZN_FRC(M,NTX),IBETAX,
     &        IZN,NPHAZ_FRC(2,NTX) )
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR( SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(M,NTX) = PSO_FRC(M,NTX)
            PTX = PG_FRC(M,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(M,NTX) = PTX
              PL_FRC(M,NTX) = PN_FRC(M,NTX) - CPNLX
            ELSE
              PL_FRC(M,NTX) = PTX
              PN_FRC(M,NTX) = PL_FRC(M,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #2  ---
!
!         SL < 1.0, ST > 0.0, SG = 0.0
!
!         energy - temperature, C - T
!         water mass - aqueous saturation, SL
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - nonaqueous (system) pressure, Pa gauge - PN,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - nonaqueous (system) pressure, Pa gauge - PN,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - nonaqueous (system) pressure, Pa gauge - PN,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt brine mass fraction - YLS
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(M,NTX),
     &        PORT_FRC(M,NTX) )
            PORD_FRC(M,NTX) = MAX( PORD_FRC(M,NTX),EPSL )
            PORT_FRC(M,NTX) = MAX( PORT_FRC(M,NTX),PORD_FRC(M,NTX) )
!
!---        Flash calculations  ---
!
            CALL FLH_PC2( BETAX,FKX,GNIFT_FRC(M,NTX),PORD_FRC(M,NTX),
     &        PVA_FRC(M,NTX),PVW_FRC(M,NTX),PX,RHOG_FRC(M,NTX),
     &        RHOL_FRC(M,NTX),RHOMG_FRC(M,NTX),RHOML_FRC(M,NTX),
     &        RHOMN_FRC(M,NTX),RHON_FRC(M,NTX),SG_FRC(M,NTX),
     &        SL_FRC(M,NTX),SN_FRC(M,NTX),TMC_FRC(1,M,NTX),
     &        TMS_FRC(M,NTX),T_FRC(M,NTX),VOL_FRC(NTX),XGC_FRC(1,M,NTX),
     &        XGW_FRC(M,NTX),XLA_FRC(M,NTX),XLS_FRC(M,NTX),
     &        XLW_FRC(M,NTX),XNC_FRC(1,M,NTX),XMGC_FRC(1,M,NTX),
     &        XMGW_FRC(M,NTX),XMLA_FRC(M,NTX),XMLS_FRC(M,NTX),
     &        XMLW_FRC(M,NTX),XMNC_FRC(1,M,NTX),YLS_FRC(M,NTX),
     &        ZG_FRC(M,NTX),ZMC_FRC(1,M,NTX),ZN_FRC(M,NTX),
     &        IBETAX,IZMC_FRC(NTX),IZN,NPHAZ_FRC(2,NTX) )
!            IF( NITER.GE.2 ) THEN
!              SG_FRC(M,NTX) = AST_FRC(NTX)*(1.D+0-SL_FRC(M,NTX))
!              SN_FRC(M,NTX) = (1.D+0-AST_FRC(NTX))*(1.D+0-SL_FRC(M,NTX))
!            ENDIF
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR( SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(M,NTX) = PSO_FRC(M,NTX)
            PTX = PG_FRC(M,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(M,NTX) = PTX
              PL_FRC(M,NTX) = PN_FRC(M,NTX) - CPNLX
            ELSE
              PL_FRC(M,NTX) = PTX
              PN_FRC(M,NTX) = PL_FRC(M,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #3  ---
!
!         SL < 1.0, ST > 0.0, SG > 0.0
!
!         energy - temperature, C - T
!         water mass - aqueous saturation, SL
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - gas (system) pressure, Pa gauge - PG,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - gas (system) pressure, Pa gauge - PG,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - gas (system) pressure, Pa gauge - PG,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt brine mass fraction - YLS
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(M,NTX),
     &        PORT_FRC(M,NTX) )
            PORD_FRC(M,NTX) = MAX( PORD_FRC(M,NTX),EPSL )
            PORT_FRC(M,NTX) = MAX( PORT_FRC(M,NTX),PORD_FRC(M,NTX) )
!
!---        Flash calculations  ---
!
            CALL FLH_PC3( BETAX,FKX,GNIFT_FRC(M,NTX),PORD_FRC(M,NTX),
     &        PVA_FRC(M,NTX),PVW_FRC(M,NTX),PX,RHOG_FRC(M,NTX),
     &        RHOL_FRC(M,NTX),RHOMG_FRC(M,NTX),RHOML_FRC(M,NTX),
     &        RHOMN_FRC(M,NTX),RHON_FRC(M,NTX),SG_FRC(M,NTX),
     &        SL_FRC(M,NTX),SN_FRC(M,NTX),TMC_FRC(1,M,NTX),
     &        TMS_FRC(M,NTX),T_FRC(M,NTX),VOL_FRC(NTX),XGC_FRC(1,M,NTX),
     &        XGW_FRC(M,NTX),XLA_FRC(M,NTX),XLS_FRC(M,NTX),
     &        XLW_FRC(M,NTX),XNC_FRC(1,M,NTX),XMGC_FRC(1,M,NTX),
     &        XMGW_FRC(M,NTX),XMLA_FRC(M,NTX),XMLS_FRC(M,NTX),
     &        XMLW_FRC(M,NTX),XMNC_FRC(1,M,NTX),YLS_FRC(M,NTX),
     &        ZG_FRC(M,NTX),ZMC_FRC(1,M,NTX),ZN_FRC(M,NTX),
     &        IBETAX,IZMC_FRC(NTX),IZN,NPHAZ_FRC(2,NTX) )
!            IF( NITER.GE.2 ) THEN
!              SG_FRC(M,NTX) = AST_FRC(NTX)*(1.D+0-SL_FRC(M,NTX))
!              SN_FRC(M,NTX) = (1.D+0-AST_FRC(NTX))*(1.D+0-SL_FRC(M,NTX))
!            ENDIF
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR( SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(M,NTX) = PSO_FRC(M,NTX)
            PTX = PG_FRC(M,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(M,NTX) = PTX
              PL_FRC(M,NTX) = PN_FRC(M,NTX) - CPNLX
            ELSE
              PL_FRC(M,NTX) = PTX
              PN_FRC(M,NTX) = PL_FRC(M,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #4  ---
!
!         SL = 0.0, ST = 1.0, SG = 0.0
!
!         energy - temperature, C - T
!         water mass - aqueous saturation, SL
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - nonaqueous (system) pressure, Pa gauge - PN,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - nonaqueous (system) pressure, Pa gauge - PN,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - nonaqueous (system) pressure, Pa gauge - PN,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(M,NTX),
     &        PORT_FRC(M,NTX) )
            PORD_FRC(M,NTX) = MAX( PORD_FRC(M,NTX),EPSL )
            PORT_FRC(M,NTX) = MAX( PORT_FRC(M,NTX),PORD_FRC(M,NTX) )
!
!---        Flash calculations  ---
!
            CALL FLH_PC4( BETAX,FKX,GNIFT_FRC(M,NTX),PORD_FRC(M,NTX),
     &        PVA_FRC(M,NTX),PVW_FRC(M,NTX),PX,RHOG_FRC(M,NTX),
     &        RHOL_FRC(M,NTX),RHOMG_FRC(M,NTX),RHOML_FRC(M,NTX),
     &        RHOMN_FRC(M,NTX),RHON_FRC(M,NTX),SG_FRC(M,NTX),
     &        SL_FRC(M,NTX),SN_FRC(M,NTX),TMC_FRC(1,M,NTX),
     &        TMS_FRC(M,NTX),T_FRC(M,NTX),VOL_FRC(NTX),XGC_FRC(1,M,NTX),
     &        XGW_FRC(M,NTX),XLA_FRC(M,NTX),XLS_FRC(M,NTX),
     &        XLW_FRC(M,NTX),XNC_FRC(1,M,NTX),XMGC_FRC(1,M,NTX),
     &        XMGW_FRC(M,NTX),XMLA_FRC(M,NTX),XMLS_FRC(M,NTX),
     &        XMLW_FRC(M,NTX),XMNC_FRC(1,M,NTX),YLS_FRC(M,NTX),
     &        ZG_FRC(M,NTX),ZMC_FRC(1,M,NTX),ZN_FRC(M,NTX),
     &        IBETAX,IZMC_FRC(NTX),IZN,NPHAZ_FRC(2,NTX) )
!            IF( NITER.GE.2 ) THEN
!              SG_FRC(M,NTX) = AST_FRC(NTX)*(1.D+0-SL_FRC(M,NTX))
!              SN_FRC(M,NTX) = (1.D+0-AST_FRC(NTX))*(1.D+0-SL_FRC(M,NTX))
!            ENDIF
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR( SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(M,NTX) = PSO_FRC(M,NTX)
            PTX = PG_FRC(M,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(M,NTX) = PTX
              PL_FRC(M,NTX) = PN_FRC(M,NTX) - CPNLX
            ELSE
              PL_FRC(M,NTX) = PTX
              PN_FRC(M,NTX) = PL_FRC(M,NTX) + CPNLX
            ENDIF
!
!---      Phase condition #5  ---
!
!         SL = 0.0, ST >= 0.0, SG > 0.0
!
!         energy - temperature, C - T
!         water mass - aqueous saturation, SL
!         CO2 mass
!           if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!             - gas (system) pressure, Pa gauge - PG,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!           endif
!         CH4 mass
!           if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!             - gas (system) pressure, Pa gauge - PG,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!           endif
!         petroleum component mass
!           if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!             - gas (system) pressure, Pa gauge - PG,PSO
!           else
!             - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!           endif
!         salt mass - total salt mass - TMS_FRC(2,NTX)
!
          ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
!
!---        System pressure, absolute  ---
!
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Porous-media porosity  ---
!
            CALL PORSTY_FRC( NTX,PX,PCMP_FRC(NTX),PORD_FRC(M,NTX),
     &        PORT_FRC(M,NTX) )
            PORD_FRC(M,NTX) = MAX( PORD_FRC(M,NTX),EPSL )
            PORT_FRC(M,NTX) = MAX( PORT_FRC(M,NTX),PORD_FRC(M,NTX) )
!
!---        Flash calculations  ---
!
            CALL FLH_PC5( BETAX,FKX,GNIFT_FRC(M,NTX),PORD_FRC(M,NTX),
     &        PVA_FRC(M,NTX),PVW_FRC(M,NTX),PX,RHOG_FRC(M,NTX),
     &        RHOL_FRC(M,NTX),RHOMG_FRC(M,NTX),RHOML_FRC(M,NTX),
     &        RHOMN_FRC(M,NTX),RHON_FRC(M,NTX),SG_FRC(M,NTX),
     &        SL_FRC(M,NTX),SN_FRC(M,NTX),TMC_FRC(1,M,NTX),
     &        TMS_FRC(M,NTX),T_FRC(M,NTX),VOL_FRC(NTX),
     &        XGC_FRC(1,M,NTX),XGW_FRC(M,NTX),XLA_FRC(M,NTX),
     &        XLS_FRC(M,NTX),XLW_FRC(M,NTX),XNC_FRC(1,M,NTX),
     &        XMGC_FRC(1,M,NTX),XMGW_FRC(M,NTX),XMLA_FRC(M,NTX),
     &        XMLS_FRC(M,NTX),XMLW_FRC(M,NTX),XMNC_FRC(1,M,NTX),
     &        YLS_FRC(M,NTX),ZG_FRC(M,NTX),ZMC_FRC(1,M,NTX),
     &        ZN_FRC(M,NTX),IBETAX,IZMC_FRC(NTX),IZN,NPHAZ_FRC(2,NTX) )
!            IF( NITER.GE.2 ) THEN
!              SG_FRC(M,NTX) = AST_FRC(NTX)*(1.D+0-SL_FRC(M,NTX))
!              SN_FRC(M,NTX) = (1.D+0-AST_FRC(NTX))*(1.D+0-SL_FRC(M,NTX))
!            ENDIF
!
!---        Phase capillary pressure from phase saturations  ---
!
            CALL CAP_EOR( SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &        CPGTX,CPNLX,IZN )
!
!---        Phase pressure from phase capillary pressures and system
!           pressure  ---
!
            PG_FRC(M,NTX) = PSO_FRC(M,NTX)
            PTX = PG_FRC(M,NTX) - CPGTX
            IF( CPNLX.GT.0.D+0 ) THEN
              PN_FRC(M,NTX) = PTX
              PL_FRC(M,NTX) = PN_FRC(M,NTX) - CPNLX
            ELSE
              PL_FRC(M,NTX) = PTX
              PN_FRC(M,NTX) = PL_FRC(M,NTX) + CPNLX
            ENDIF
          ENDIF
!
!---      Save unincremented K-factors  ---
!
          IF( M.EQ.2 .AND. NITER.LE.3 ) THEN
            DO IGC = 1,NGC+2
              FK_FRC(IGC,NTX) = FKX(IGC)
            ENDDO
          ENDIF
!
!---      Load component aqueous mass fraction for CO2  ---
!
          XLC_FRC(1,M,NTX) = XLA_FRC(M,NTX)
          IF( M.EQ.2 ) THEN
!
!---        Store recent iterative values of beta  ---
!
            BETA_FRC(1,NTX) = BETA_FRC(2,NTX)
            BETA_FRC(2,NTX) = BETA_FRC(3,NTX)
            BETA_FRC(3,NTX) = BETAX
!
!---        Store recent iterative values of gas-nonaqueous liquid
!           interfacial tension ---
!
            BETA_FRC(4,NTX) = BETA_FRC(5,NTX)
            BETA_FRC(5,NTX) = BETA_FRC(6,NTX)
            BETA_FRC(6,NTX) = GNIFT_FRC(2,NTX)
          ENDIF
        ENDDO
      ENDDO
      ENDDO
!
!---  Subroutine split, loop over faults  ---
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
        N_DB = NTX
        IZN = IZ_FRC(NTX)
!
!---    Loop over increment indices  ---
!
        DO M = 2,ISVC+2
          TKX = T_FRC(M,NTX) + TABS
!
!---      System pressure, absolute  ---
!
          PX = PSO_FRC(M,NTX) + PATM
!
!---      Relative permeability  ---
!
#ifdef dual_porosity
          CALL RKG_EOR_F( RKG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),
     &      IZN )
          CALL RKL_EOR_F( RKLX,SL_FRC(M,NTX),IZN )
          RKL_FRC(M,NTX) = RKLX(1)
          CALL RKN_EOR_F( RKG_FRC(M,NTX),RKLX,RKN_FRC(M,NTX),
     &      SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),IZN )
#else
          CALL RKG_EOR( RKG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),IZN )
          CALL RKL_EOR( RKLX,SL_FRC(M,NTX),IZN )
          RKL_FRC(M,NTX) = RKLX(1)
          CALL RKN_EOR( RKG_FRC(M,NTX),RKLX,RKN_FRC(M,NTX),
     &      SG_FRC(M,NTX),SL_FRC(M,NTX),SN_FRC(M,NTX),IZN )
#endif
!
!---      Gas viscosity, Pa s  ---
!
          CALL VISC_G( RHOMG_FRC(M,NTX),T_FRC(M,NTX),VISG_FRC(M,NTX),
     &      XMGC_FRC(1,M,NTX),XMGW_FRC(M,NTX) )
!
!---      Aqueous viscosity, Pa s
!         concentration of CO2 in mol/kg,
!         concentration of salt in mol/kg  ---
!
          CMLAX = 1.D+3*XMLA_FRC(M,NTX)*RHOML_FRC(M,NTX)/RHOL_FRC(M,NTX)
          CMLSX = 1.D+3*XMLS_FRC(M,NTX)*RHOML_FRC(M,NTX)/RHOL_FRC(M,NTX)
          CALL VISC_L( CMLAX,CMLSX,PX,T_FRC(M,NTX),VISL_FRC(M,NTX) )
!
!---      Nonaqueous-liquid viscosity, Pa s  ---
!
          CALL VISC_N( RHOMN_FRC(M,NTX),T_FRC(M,NTX),VISN_FRC(M,NTX),
     &      XMNC_FRC(1,M,NTX) )
!
!---      Component gas diffusion coefficients  ---
!
          CALL DIFC_G( DFGC_FRC(1,M,NTX),DFGW_FRC(M,NTX),PX,
     &      T_FRC(M,NTX),XMGC_FRC(1,M,NTX),XMGW_FRC(M,NTX) )
!
!---      CO2 aqueous diffusion coefficient  ---
!
          CALL DIFC_LA( DFLA_FRC(M,NTX),T_FRC(M,NTX),VISL_FRC(M,NTX) )
!
!---      Salt aqueous diffusion coefficient  ---
!
          CALL DIFC_LS( DFLS_FRC(M,NTX),T_FRC(M,NTX),VISL_FRC(M,NTX),
     &      XLS_FRC(M,NTX) )
!
!---      Petroleum component nonaqueous-liquid diffusion 
!         coefficients  ---
!
          CALL DIFC_N( DFNC_FRC(1,M,NTX),T_FRC(M,NTX),VISN_FRC(M,NTX) )
!
!---      Gas thermal conductivity, W/m K  ---
!
          CALL THK_G( PX,T_FRC(M,NTX),THKG_FRC(M,NTX),XMGC_FRC(1,M,NTX),
     &      XMGW_FRC(M,NTX) )
!
!---      Aqueous thermal conductivity  ---
!
          CALL THK_W( PX,T_FRC(M,NTX),THKLWX )
          CALL THK_B( T_FRC(M,NTX),XLS_FRC(M,NTX),THKLWX,
     &      THKL_FRC(M,NTX) )
!
!---      Nonaqueous thermal conductivity, W/m K  ---
!
          CALL THK_N( PX,T_FRC(M,NTX),THKN_FRC(M,NTX),XMNC_FRC(1,M,NTX))
!
!---      Gas enthalpy and internal energy, W/kg  ---
!
          CALL ENTH_G( HG_FRC(M,NTX),HGA_FRC(M,NTX),HGW_FRC(M,NTX),PX,
     &      T_FRC(M,NTX),RHOG_FRC(M,NTX),UEG_FRC(M,NTX),
     &      XGC_FRC(1,M,NTX),XGW_FRC(M,NTX),XMGC_FRC(1,M,NTX),
     &      ZG_FRC(M,NTX) )
!
!---      Nonaqueous-liquid enthalpy and internal energy  ---
!
          CALL ENTH_N( HN_FRC(M,NTX),HNAX,PX,T_FRC(M,NTX),
     &     RHON_FRC(M,NTX),UEN_FRC(M,NTX),XNC_FRC(1,M,NTX),
     &     XMNC_FRC(1,M,NTX),ZN_FRC(M,NTX) )
!
!---      Pure water liquid enthalpy  ---
!
          CALL SP_W( T_FRC(M,NTX),PSWX )
          PEX = MAX( PSWX,PX )
          ISRX = 1
          CALL ENTH_W( T_FRC(M,NTX),PEX,HLW_FRC(M,NTX),HGWX,ISRX )
!
!---      Brine enthalpy  ---
!
          CALL ENTH_B( T_FRC(M,NTX),XLS_FRC(M,NTX),HLW_FRC(M,NTX),HBX )
!!
!!---      Dissolved CO2 contribution  ---
!!
!          SGX = (BETAX/RHOMG_FRC(M,NTX))/
!     &      (BETAX/RHOMG_FRC(M,NTX) + (1.D+0-BETAX)/RHOMN_FRC(M,NTX))
!          IF( SGX.LT.1.D-12 ) SGX = 0.D+0
!          SNX = 1.D+0 - SGX
!          FGAX = SGX*RHOG_FRC(M,NTX)*XGC_FRC(1,M,NTX)
!          FNAX = SNX*RHON_FRC(M,NTX)*XNC_FRC(1,M,NTX)
          IF( SG_FRC(M,NTX).GT.0.D+0 ) THEN
            FGAX = SG_FRC(M,NTX)*RHOG_FRC(M,NTX)*XGC_FRC(1,M,NTX)
          ELSE
            FGAX = 0.D+0
          ENDIF
          IF( SN_FRC(M,NTX).GT.0.D+0 ) THEN
            FNAX = SN_FRC(M,NTX)*RHON_FRC(M,NTX)*XNC_FRC(1,M,NTX)
          ELSE
            FNAX = 0.D+0
          ENDIF
          IF( (FGAX+FNAX)/EPSL.LT.EPSL ) THEN
            HAX = 0.D+0
          ELSE
            HAX = (HGA_FRC(M,NTX)*FGAX + HNAX*FNAX)/(FGAX + FNAX)
          ENDIF
!
!---      Aqueous enthalpy and internal energy  ---
!
          CALL ENTH_L( HAX,HBX,HL_FRC(M,NTX),T_FRC(M,NTX),
     &      XLA_FRC(M,NTX),XLS_FRC(M,NTX) )
          UEL_FRC(M,NTX) = HL_FRC(M,NTX) - PX/RHOL_FRC(M,NTX)
!
!---      Precipitated NaCl density, saturation, and enthalpy  ---
!
          CALL DENS_S( T_FRC(M,NTX),PX,RHOSP_FRC(M,NTX) )
          CALL SOL_LS( T_FRC(M,NTX),XBSMX )
          CALL DENS_B( T_FRC(M,NTX),PX,XBSMX,RHOBX )
          SS_FRC(M,NTX) = MAX(YLS_FRC(M,NTX)-XBSMX,0.D+0)*RHOBX*
     &      SL_FRC(M,NTX)/RHOSP_FRC(M,NTX)
!
!---      Precipitated NaCl enthalpy  ---
!
          CALL ENTH_S( T_FRC(M,NTX),HSP_FRC(M,NTX) )
!
!---      Permeability reduction factor  ---
!
          IF( IPRF(IZN).EQ.1 ) CALL PERM_R( SS_FRC(M,NTX),
     &      PERMRF_FRC(M,NTX),PORD_FRC(M,NTX),IZN )
!
!---      Change in porosity as a function of precipitated salt  ---
!
          PORD_FRC(M,NTX) = PORD_FRC(M,NTX)*(1.D+0-SS_FRC(M,NTX))
!
!---      Aqueous, gas, and nonaqueous-liquid tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL_FRC(M,NTX),
     &      SG_FRC(M,NTX),SN_FRC(M,NTX),PORD_FRC(M,NTX),TORL_FRC(M,NTX),
     &      TORG_FRC(M,NTX),TORN_FRC(M,NTX) )
        ENDDO
        CPUMX = CPUMX + VOL_FRC(NTX)*PORD_FRC(2,NTX)*
     &    (SL_FRC(2,NTX)*RHOL_FRC(2,NTX) + 
     &    SG_FRC(2,NTX)*RHOG_FRC(2,NTX) + 
     &    SN_FRC(2,NTX)*RHON_FRC(2,NTX))
      ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_FRC_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIC_FRC_EOR
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
!     Read input file for fault initial conditions information.
!
!     VIC_FRC(1,NFX) - x reference point, m
!     VIC_FRC(2,NFX) - y reference point, m
!     VIC_FRC(3,NFX) - z reference point, m
!     VIC_FRC(5,NFX) - temperature, C
!     VIC_FRC(6,NFX) - x gradient, C/m
!     VIC_FRC(7,NFX) - y gradient, C/m
!     VIC_FRC(8,NFX) - z gradient, C/m
!     VIC_FRC(9,NFX) - pressure, Pa
!     VIC_FRC(10,NFX) - x gradient, Pa/m
!     VIC_FRC(11,NFX) - y gradient, Pa/m
!     VIC_FRC(12,NFX) - z gradient, Pa/m
!     VIC_FRC(13,NFX) - aqueous saturation or water vapor rel. sat.
!     VIC_FRC(14,NFX) - x gradient, 1/m
!     VIC_FRC(15,NFX) - y gradient, 1/m
!     VIC_FRC(16,NFX) - z gradient, 1/m
!     VIC_FRC(17,NFX) - salt mass frac. or salt total mass, kg
!     VIC_FRC(18,NFX) - x gradient, 1/m or kg/m
!     VIC_FRC(19,NFX) - y gradient, 1/m or kg/m
!     VIC_FRC(20,NFX) - z gradient, 1/m or kg/m
!     VIC_FRC(21,NFX) - CO2 rel. sat. or mass frac.
!     VIC_FRC(22,NFX) - x gradient, 1/m
!     VIC_FRC(23,NFX) - y gradient, 1/m
!     VIC_FRC(24,NFX) - z gradient, 1/m
!     VIC_FRC(21+(IGC*4),NFX) - petroleum component total mole fraction
!     VIC_FRC(22+(IGC*4),NFX) - x gradient, 1/m
!     VIC_FRC(23+(IGC*4),NFX) - y gradient, 1/m
!     VIC_FRC(24+(IGC*4),NFX) - z gradient, 1/m
!     VIC_FRC(21+(NGC*4)+(NSL*4),NFX) - solute concentration, 1/m^3
!     VIC_FRC(21+(NGC*4)+(NSL*4),NFX) - x gradient, 1/m^4
!     VIC_FRC(21+(NGC*4)+(NSL*4),NFX) - y gradient, 1/m^4
!     VIC_FRC(21+(NGC*4)+(NSL*4),NFX) - z gradient, 1/m^4
!     VIC_FRC(21+(NGC*4)+(NSOLU*4)+(NSP*4),NFX) - species conc., 1/m^3
!     VIC_FRC(21+(NGC*4)+(NSOLU*4)+(NSP*4),NFX) - x gradient, 1/m^4
!     VIC_FRC(21+(NGC*4)+(NSOLU*4)+(NSP*4),NFX) - y gradient, 1/m^4
!     VIC_FRC(21+(NGC*4)+(NSOLU*4)+(NSP*4),NFX) - z gradient, 1/m^4
!
!----------------------Authors-----------------------------------------!
!
!     Written by Ramesh Sarathi, PNNL, 6 March 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PARM_FRC
      USE GLB_PAR
      USE GEOM_FRC
      USE FILES
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

      CHARACTER*24 CHLB(3)
      INTEGER II ! IC variable position counter
!
!----------------------Data Statements---------------------------------!
!
      SAVE CHLB
      DATA CHLB /'X-Direction Gradient, ','Y-Direction Gradient, ',
     &           'Z-Direction Gradient, '/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIC_FRC_EOR'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Fault Initial Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '

!---  Read number of fault inputs, unspecified faults
!     will be in thermodynamic equilibrium with connecting nodes  ---
!
      WRITE(IWR,'(/,A)') 'Fault Initial Condition Parameters'
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Fault Initial Condition Inputs: '
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
!
!---  Loop over IC definition blocks (up to one block per 
!     fault/fault ) ---
!
      DO NL = 1, NLIN
!
!---  First line of block lists fault number, 
!     initial condition state, reference location ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Fault Number: '
        CALL RDINT(ISTART,ICOMMA,CHDUM,NFX)
        WRITE(IWR,'(2X,A,I9)') 'Fault Number: ',NFX
        IF( NFX.LT.1 .OR. NFX.GT.LF_FRC ) THEN
          INDX = 7
          IMSG = NFX
          CHMSG = 'Out of Range Fault Number'
          CALL WRMSGS( INDX )
        ENDIF

        VARB = 'Fault Initial Condition State: '
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
!
!---    Skip for restarts(IEO==2) unless overwrite specified  ---
!
        IF( INDEX( ADUM(1:),'overwrite').EQ.0 .AND.
     &    ( IEO.EQ.2 ) ) THEN
          CYCLE
!
!---    Read initial condition state  ---
!
        ELSEIF( INDEX( ADUM(1:),'ic1' ).NE.0 ) THEN
          VARB = 'Initial Condition State #1'
          IC_OPT_FRC(1,NFX) = 1
        ELSEIF( INDEX( ADUM(1:),'ic2' ).NE.0 ) THEN
          VARB = 'Initial Condition State #2'
          IC_OPT_FRC(1,NFX) = 2
        ELSEIF( INDEX( ADUM(1:),'ic3' ).NE.0 ) THEN
          VARB = 'Initial Condition State #3'
          IC_OPT_FRC(1,NFX) = 3
        ELSEIF( INDEX( ADUM(1:),'equilibrium' ).NE.0 ) THEN
          VARB = 'Equilibrium Initial Condition State'
          IC_OPT_FRC(1,NFX) = 10
          IC_OPT_FRC(2,NFX) = 1
          IC_OPT_FRC(3,NFX) = 1
          IC_OPT_FRC(4,NFX) = 1
          IC_OPT_FRC(5,NFX) = 1
          IC_OPT_FRC(6,NFX) = 1
          IC_OPT_FRC(7,NFX) = 3
          IC_OPT_FRC(8,NFX) = 3
          IC_OPT_FRC(9,NFX) = 1
          IC_OPT_FRC(10,NFX) = 1
          IC_OPT_FRC(11,NFX) = 1
          CYCLE
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Fault Initial ' // 
     &      'Condition Variable: ' // ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Read reference location  ---
!
        VARB = 'Reference X Location'
        II = 1 ! put the ref loc in 1:3, for consistency 4 is empty
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(II,NFX))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(II,NFX)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VIC_FRC(II,NFX),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(II,NFX),', m)'

        VARB = 'Reference Y Location'
        II = II+1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(II,NFX))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(II,NFX)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VIC_FRC(II,NFX),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(II,NFX),', m)'

        VARB = 'Reference Z Location'
        II = II+1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(II,NFX))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(II,NFX)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VIC_FRC(II,NFX),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(II,NFX),', m)'
!
!---    Read Next line, number of IC lines for the block ---
!
        WRITE(IWR,'(/,A,I4,A)') 'Fault ', NFX,
     &    'Initial Condition Inputs'
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Number of Fault Initial Condition Inputs: '
        CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN_NIC)
!
!---    Loop over number of initial condition lines for the 
!       fault  ---
!
        DO ILIN_NIC=1,NLIN_NIC
!
!---      Read input line  ---
!
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Fault Initial Condition Variable: '
          CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
          ICX = 0
          II = 0
          VARX = 1.D0
!
!---      Match IC statements label1, [label2], and set options 
!         and UCFs ---
!
          IF( INDEX( ADUM(1:),'temperature' ).NE.0 ) THEN
            VARB = 'Fault Initial Temperature'
            IC_OPT_FRC(9,NFX) = 1 ! mark temperature as specified
            II = 1+1*4 ! IC variable array index
            ! Set unit conversion factors
            INDX = 0
            IUNK = 1
          ELSEIF( INDEX( ADUM(1:),'pressure' ).NE.0 ) THEN
            VARB = 'Fault Initial Pressure'
            IC_OPT_FRC(10,NFX) = 1 ! mark pressure as specified
            II = 1+2*4
            ! Set unit conversion factors
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
          ELSEIF( INDEX( ADUM(1:),'aqueous' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'saturation' ).NE.0 ) THEN
            VARB = 'Fault Initial Aqueous Saturation'
            IC_OPT_FRC(3,NFX) = 1 ! mark aqueous saturation as specified
            II = 1+3*4
            ! no unit conversion
          ELSEIF( INDEX( ADUM(1:),'water' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'vapor' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
            VARB = 'Fault Initial Water-Vapor Relative ' // 
     &        'Saturation'
            IC_OPT_FRC(11,NFX) = 1
            II = 1+3*4
            ! no unit conversion
          ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
            VARB = 'Fault Initial Aqueous-Salt Relative ' // 
     &        'Saturation'
            IC_OPT_FRC(8,NFX) = 2 ! mark salt input dof as rel sat
            II = 1+4*4
            ! no unit conversion
          ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &      ( INDEX( ADUM(1:),'mass' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'frac' ).NE.0 ) ) THEN
            VARB = 'Fault Initial Aqueous-Salt Mass ' // 
     &        'Fraction'
            IC_OPT_FRC(8,NFX) = 3 ! mark salt input dof as mass frac
            II = 1+4*4
            ! no unit conversion
          ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &      ( INDEX( ADUM(1:),'mass' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'total' ).NE.0 ) ) THEN
            VARB = 'Fault Initial Salt Mass Total'
            IC_OPT_FRC(8,NFX) = 4 ! mark salt input dof as mass total
            II = 1+4*4
            !rss_ask should there be a unit conversion
            ! Set unit conversion factors
            INDX = 0
            IUNKG = 1
          ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
            VARB = 'Fault Initial Aqueous-CO2 Relative ' // 
     &        'Saturation'
            IC_OPT_FRC(7,NFX) = 2 ! mark co2 input dof as rel sat
            II = 1+5*4
            ! no unit conversion
          ELSEIF( INDEX( ADUM(1:),'co2' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'mass' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'frac' ).NE.0 ) THEN
            VARB = 'Fault Initial Aqueous-CO2 Mass Fraction'
            IC_OPT_FRC(7,NFX) = 3 ! mark co2 input dof as mass frac
            II = 1+5*4
            ! no unit conversion
          ELSEIF( INDEX( ADUM(1:),'total' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'mole' ).NE.0 .AND.
     &      INDEX( ADUM(1:),'fraction' ).NE.0 ) THEN
            VARB = 'Fault Petroleum Component Total Mole ' // 
     &        'Fraction'
            IC_OPT_FRC(13,NFX) = 1 ! mark total component mole frac as specified
            
            ! Read [label2], the component name
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            
            ! Find the index corresponding to the component name
            IGCX = 0
            DO IGC = 1,NGC+2 ! plus two for CO2, CH4
              IDB = INDEX(GCNM(IGC)(1:),'  ') - 1
              IF( BDUM(1:NCHB).EQ.GCNM(IGC)(1:IDB) ) IGCX = IGC
            END DO
            IF( IGCX.EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Unrecognized Fault Petroleum ' // 
     &          'Component Name: ' // BDUM(1:NCHB)
              CALL WRMSGS( INDX )
              ! BAIL
              CYCLE
            ENDIF
            ! Set VARB to include the component name
            VARB = BDUM(1:NCHB) // ' Fault Total Mole Fraction'
            ! Set the IC variable index
            II = 1+(6+IGCX)*4 ! one for each NGC+2 components
          ELSEIF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
            VARB = 'Solute Name: '
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            NSLX = 0
            DO NSL = 1,NSOLU
              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
              IF( BDUM(1:NCHB).EQ.SOLUT(NSL)(1:IDB) ) THEN
                NSLX = NSL
                EXIT
              ENDIF
            ENDDO
            IF( NSLX.EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Unrecognized Solute: ' // BDUM(1:NCHB)
              CALL WRMSGS( INDX )
            ENDIF
            IF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
              ICT_FRC(NFX,NSLX) = 3
              VARB = BDUM(1:NCHB) // 
     &        'Fault Reference Solute Gas Concentration'
            ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
              ICT_FRC(NFX,NSLX) = 2
              VARB = BDUM(1:NCHB) // 
     &        'Fault Reference Solute Aqueous Concentration'
            ELSE
              ICT_FRC(NFX,NSLX) = 1
              VARB = BDUM(1:NCHB) // 
     &        'Fault Reference Solute Volumetric Concentration'
            ENDIF
            II = 1+(6+NGC+NSLX)*4
            ! Set unit conversion factors
            INDX = 0
            IUNM = -3
#ifdef ecke
          ELSEIF( INDEX( ADUM(1:),'specie' ).NE.0 ) THEN
            VARB = 'Reactive Species Name: '
            CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
            NSL = 0
            NSP = 0
!
!---        Conservation- or kinetic-component species  ---
!
            IF( INDEX( BDUM(1:),'total_' ).NE.0 ) THEN
              DO NSLX = NSOLU+1,NSOLU+NEQC+NEQK
                IDB = INDEX(SOLUT(NSLX)(1:),'  ') - 1
                IF( BDUM(1:NCHB).EQ.SOLUT(NSLX) ) THEN
                  NSL = NSLX
                  EXIT
                ENDIF
              ENDDO
            ELSE
!
!---          Aqueous reactive species  ---
!
              DO NSPX = 1,NSPL
                IDB = INDEX(SPNML(NSPX)(1:),'  ') - 1
                IF( BDUM(1:NCHB).EQ.SPNML(NSPX)(1:IDB) ) THEN
                  NSP = NSPX
                  EXIT
                ENDIF
              ENDDO
!
!---          Solid reactive species  ---
!
              IF( NSP.EQ.0 ) THEN
                DO NSPX = 1,NSPS
                  IDB = INDEX(SPNMS(NSPX)(1:),'  ') - 1
                  IF( BDUM(1:NCHB).EQ.SPNMS(NSPX)(1:IDB) ) THEN
                    NSP = NSPX + NSPL
!
!---                Verify that solid-species is not a mineral  ---
!
                    IF( ISP_MN(NSP).EQ.1 ) THEN
                      INDX = 4
                      CHMSG = 'Solid-Species Mineral ' // 
     &                '(see Lithology Card): ' // BDUM(1:NCHB)
                      CALL WRMSGS( INDX )
                    ENDIF
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
!
!---          Gas reactive species  ---
!
              IF( NSP.EQ.0 ) THEN
                DO NSPX = 1,NSPG
                  IDB = INDEX(SPNMG(NSPX)(1:),'  ') - 1
                  IF( BDUM(1:NCHB).EQ.SPNMG(NSPX)(1:IDB) ) THEN
                    NSP = NSPX + NSPL + NSPS
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
!
!---      pH  ---
!
              IF( NSP.EQ.0 ) THEN
                IF( BDUM(1:NCHB).EQ.'ph' .AND. ISPLK(1).NE.0 ) THEN
                  NSP = MOD(ISPLK(1),1000)
                  ISPLK(1) = ISPLK(1) + 1000
                  IVAR = 2
                  IF( IEO.EQ.2 ) IVAR = IVAR+10
                  ADDER = 7.D+0
!
!---              Verify that species linked to pH is a conservation
!                 component species  ---
!
                  IPHX = 0
                  DO NEQ = 1,NEQC
                    IF( NSP.EQ.IEQ_C(2,NEQ) ) THEN
                      IPHX = 1
                      EXIT
                    ENDIF
                  ENDDO
                  IF( IPHX.EQ.0 ) THEN
                    INDX = 4
                    CHMSG = 'pH Species not a Conservation ' //
     &                'Component Species: ' // BDUM(1:NCHB)
                    CALL WRMSGS( INDX )
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
!
!---        Conservation- or kinetic-component species  ---
!
            IF( INDEX( BDUM(1:),'total_' ).NE.0 ) THEN
              IF( NSL.EQ.0 ) THEN
                INDX = 4
                CHMSG = 'Unrecognized Reactive Species: ' //BDUM(1:NCHB)
                CALL WRMSGS( INDX )
              ENDIF
!
!---          Set species units  ---
!
              IUNMOL = 1
              IF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
                VARB = BDUM(1:NCHB) // 'Fault ' //
     &          'Reference Total-Species Gas Concentration'
                ICT_FRC(NFX,NSL) = 4
                IUNM = -3
              ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
                IF( INDEX(ADUM(1:),'molal').NE.0 ) THEN
                  VARB = BDUM(1:NCHB) // 'Fault ' //
     &            'Reference Total-Species Aqueous Molal Concentration'
                  ICT_FRC(NFX,NSL) = 3
                  IUNKG = -1
                ELSE
                  VARB = BDUM(1:NCHB) // 'Fault ' //
     &            'Reference Total-Species Aqueous Concentration'
                  ICT_FRC(NFX,NSL) = 2
                  IVAR = 2
                  IUNM = -3
                ENDIF
              ELSE
                VARB = BDUM(1:NCHB) // 'Fault ' //
     &          'Reference Total-Species Volumetric Concentration'
                ICT_FRC(NFX,NSL) = 1
                IVAR = 1
                IUNM = -3
              ENDIF
              IF( IEO.EQ.2 ) IVAR = IVAR+10
              II = 1+(6+NGC+NSL)*4
            ELSE
              IF( NSP.EQ.0 ) THEN
                INDX = 4
                CHMSG = 'Unrecognized Reactive Species: ' //BDUM(1:NCHB)
                CALL WRMSGS( INDX )
              ENDIF
!
!---          Set species units  ---
!
              IUNMOL = 1
              IF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
                IC_SP_FRC(NFX,NSP) = 4
                IUNM = -3
                VARB = BDUM(1:NCHB) // 'Fault ' //
     &          'Reference Species Gas Concentration'
              ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
                IF( INDEX(ADUM(1:),'molal').NE.0 ) THEN
                  VARB = BDUM(1:NCHB) // 'Fault ' //
     &            'Reference Species Aqueous Molal Concentration'
                  IC_SP_FRC(NFX,NSP) = 3
                  IUNKG = -1
                ELSE
                  VARB = BDUM(1:NCHB) // 'Fault ' //
     &            'Reference Species Aqueous Concentration'
                  IC_SP_FRC(NFX,NSP) = 2
                  IUNM = -3
                ENDIF
              ELSE
                VARB = BDUM(1:NCHB) // 'Fault ' //
     &          'Reference Species Volumetric Concentration'
                IC_SP_FRC(NFX,NSP) = 1
                IUNM = -3
              ENDIF
              IF( IEO.EQ.2 ) IVAR = IVAR+10
              II = 1+(6+NGC+NSOLU+NSP)*4
            ENDIF
#endif
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Fault Initial Condition ' // 
     &        'Variable: ' // ADUM(1:NCHA)
            CALL WRMSGS( INDX )
!
!---      End IF parsing IC definition lines  ---
!
          ENDIF
!
!---      Parse the rest of the IC statement line, convert units  ---
!
          ! Read value, unit
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(II,NFX))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(II,NFX)
          
          ! Convert units
          CALL RDUNIT(UNTS,VIC_FRC(II,NFX),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(II,NFX),', 1/m)'
          
          ! Create unit conversion factor for the gradients
          INDX = 2 ! 2 means do not check unit dimensions
          VARX = 1.D+0
          ! Special logic for temperature scaling (but no translation)
          NCH = INDEX( UNTS,'  ' ) - 1
          IF( UNTS(1:NCH).EQ.'f' .OR. UNTS(1:NCH).EQ.'r' ) THEN
            VARX = 1.D+0/1.8D+0
          ELSEIF( UNTS(1:NCH).EQ.'c' .OR. UNTS(1:NCH).EQ.'k' ) THEN
            VARX = 1.D+0
          ELSE
            CALL RDUNIT( UNTS,VARX,INDX )
          ENDIF
          
          ! Loop over gradient directions
          DO I=1,3
            II = II+1
            ! Read value, unit
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VIC_FRC(II,NFX))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,3A,$)') VARB(1:IVR), ' ', CHLB(I)
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VIC_FRC(II,NFX)
            ! apply base unit convertion factor
            VIC_FRC(II,NFX) = VIC_FRC(II,NFX) * VARX
            ! now convert the gradient-length-unit
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,VIC_FRC(II,NFX),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VIC_FRC(II,NFX),', 1/m)'
          ENDDO ! end loop over gradients

        ! Black oil options
        ! Conservation- or kinetic-component species  ---
        ! Aqueous reactive species  ---
        ! Solid reactive species  ---
        ! Gas reactive species  ---
        ! pH  ---
!
!---    End Loop over IC definition lines  ---
!
        ENDDO
!
!---  End Loop over IC blocks  ---
!
      ENDDO
!
!---  Loop over faults, checking for uninitialized
!     hydrologic state conditions for non-restart simulations  ---
!
      IF( IEO.NE.2 ) THEN
        DO NFX = 1,NF_FRC
          IF( IC_OPT_FRC(1,NFX).EQ.0 ) THEN
            INDX = 7
            IMSG = NFX
            CHMSG = 'Uninitialized Initial Condition for Fault'
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
      SUBROUTINE RDSR_FRC_EOR
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
!     Read input file for fault source information.
!
!     Source types
!
!     1 - Power, W
!     2 - Power density, W/m^3
!     35 - Pressure sink
!     122 - Fluid mass, kg/s, state #1, salt rel. sat., CO2 rel. sat.
!     123 - Fluid mass, kg/s, state #1, salt rel. sat., CO2 mass frac.
!     132 - Fluid mass, kg/s, state #1, salt mass frac., CO2 rel. sat.
!     133 - Fluid mass, kg/s, state #1, salt mass frac., CO2 mass frac.
!     220 - Fluid mass, kg/s, state #2, salt rel. sat.
!     230 - Fluid mass, kg/s, state #2, salt mass frac.
!     300 - Fluid mass, kg/s, state #3
!     422 - Fluid vol., m^3/s, state #1, salt rel. sat., CO2 rel. sat.
!     423 - Fluid vol., m^3/s, state #1, salt rel. sat., CO2 mass frac.
!     432 - Fluid vol., m^3/s, state #1, salt mass frac., CO2 rel. sat.
!     433 - Fluid vol., m^3/s, state #1, salt mass frac., CO2 mass frac.
!     520 - Fluid vol., m^3/s, state #2, salt rel. sat.
!     530 - Fluid vol., m^3/s, state #2, salt mass frac.
!     600 - Fluid vol., m^3/s, state #3
!
!     Source variables
!
!     SRC_FRC(1,NTM,NSR) - time, s
!     SRC_FRC(2,NTM,NSR) - power, W
!     SRC_FRC(2,NTM,NSR) - power density, W/m^3
!     SRC_FRC(2,NTM,NSR) - fluid mass rate, kg/s
!     SRC_FRC(2,NTM,NSR) - fluid volumetric rate, m^3/s
!     SRC_FRC(2,NTM,NSR) - borehole skin factor,,
!     SRC_FRC(3,NTM,NSR) - fluid temperature, C
!     SRC_FRC(3,NTM,NSR) - borehole radius, m
!     SRC_FRC(4,NTM,NSR) - aqueous CO2 relative saturation
!     SRC_FRC(4,NTM,NSR) - aqueous CO2 mass fraction
!     SRC_FRC(4,NTM,NSR) - aqueous saturation
!     SRC_FRC(4,NTM,NSR) - solute rate, 1/s
!     SRC_FRC(4,NTM,NSR) - solute density rate, 1/s m^3
!     SRC_FRC(4,NTM,NSR) - sink pressure, Pa
!     SRC_FRC(5,NTM,NSR) - aqueous salt relative saturation
!     SRC_FRC(5,NTM,NSR) - aqueous salt mass fraction
!     SRC_FRC(6,NTM,NSR) - CO2 nonaqueous mole fraction
!     SRC_FRC(7,NTM,NSR) - CH4 nonaqueous mole fraction
!     SRC_FRC(7+I,NTM,NSR) - petroleum component I nonaqueous 
!       mole fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by Ramesh Sarathi, PNNL, 6 March 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PARM_FRC
      USE GLB_PAR
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
      CHARACTER*64 ADUM,BDUM(3),UNTS
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
      CARD = 'Fault Source Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      NSR_FRC = 0
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Fault Sources'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
!
!---  Loop over the number of fault sources  ---
!
      DO NS = 1,NLIN
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
!
!---    Read fault source type  ---
!
        VARB = 'Fault Source Type'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'power density').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Power Density Source'
          ISRTX = 2
        ELSEIF( INDEX(ADUM(1:),'power').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Power Source'
          ISRTX = 1
        ELSEIF( INDEX(ADUM(1:),'fluid').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Fluid Mass Source'
          VARB = 'Source State'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(1))
          IF( INDEX(BDUM(1)(1:),'sr1').NE.0 ) THEN
            ISRTX = 100
            WRITE(IWR,'(2X,A)') 'Source State: #1'
!
!---        Read aqueous CO2 option  ---
!
            VARB = 'Aqueous CO2 Option'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(2))
            IF( INDEX(BDUM(2)(1:),'co2').NE.0 .AND.
     &        INDEX(BDUM(2)(1:),'rel').NE.0 .AND.
     &        INDEX(BDUM(2)(1:),'sat').NE.0 ) THEN
              ISRTX = ISRTX + 2
              WRITE(IWR,'(2X,2A)') 'Aqueous CO2: Relative Saturation'
            ELSEIF( INDEX(BDUM(2)(1:),'co2').NE.0 .AND.
     &        INDEX(BDUM(2)(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(2)(1:),'frac').NE.0 ) THEN
              ISRTX = ISRTX + 3
              WRITE(IWR,'(2X,2A)') 'Aqueous CO2: Mass Fraction'
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Aqueous CO2 Option: '
     &          // BDUM(2)
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Read aqueous salt option  ---
!
            VARB = 'Aqueous Salt Option'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
            IF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'rel').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'sat').NE.0 ) THEN
              ISRTX = ISRTX + 20
              WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Relative Saturation'
            ELSEIF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'frac').NE.0 ) THEN
              ISRTX = ISRTX + 30
              WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Mass Fraction'
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Aqueous Salt Option: '
     &          // BDUM(3)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(1)(1:),'sr2').NE.0 ) THEN
            ISRTX = 200
            WRITE(IWR,'(2X,2A)') 'Source State: #2'
!
!---        Read aqueous salt option  ---
!
            VARB = 'Aqueous Salt Option'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
            IF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'rel').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'sat').NE.0 ) THEN
              ISRTX = ISRTX + 20
              WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Relative Saturation'
            ELSEIF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'frac').NE.0 ) THEN
              ISRTX = ISRTX + 30
              WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Mass Fraction'
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Aqueous Salt Option: '
     &          // BDUM(3)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(1)(1:),'sr3').NE.0 ) THEN
             ISRTX = 300
            WRITE(IWR,'(2X,2A)') 'Source State: #3'
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Source State: ' // BDUM(1)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'fluid').NE.0 .AND.
     &    INDEX(ADUM(1:),'volum').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Fluid Volumetric Source'
          VARB = 'Source State'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(1))
          IF( INDEX(BDUM(1)(1:),'sr1').NE.0 ) THEN
            ISRTX = 400
            WRITE(IWR,'(2X,A)') 'Source State: #1'
!
!---        Read aqueous CO2 option  ---
!
            VARB = 'Aqueous CO2 Option'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(2))
            IF( INDEX(BDUM(2)(1:),'co2').NE.0 .AND.
     &        INDEX(BDUM(2)(1:),'rel').NE.0 .AND.
     &        INDEX(BDUM(2)(1:),'sat').NE.0 ) THEN
              ISRTX = ISRTX + 2
              WRITE(IWR,'(2X,2A)') 'Aqueous CO2: Relative Saturation'
            ELSEIF( INDEX(BDUM(2)(1:),'co2').NE.0 .AND.
     &        INDEX(BDUM(2)(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(2)(1:),'frac').NE.0 ) THEN
              ISRTX = ISRTX + 3
              WRITE(IWR,'(2X,2A)') 'Aqueous CO2: Mass Fraction'
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Aqueous CO2 Option: '
     &          // BDUM(2)
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Read aqueous salt option  ---
!
            VARB = 'Aqueous Salt Option'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
            IF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'rel').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'sat').NE.0 ) THEN
              ISRTX = ISRTX + 20
              WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Relative Saturation'
            ELSEIF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'frac').NE.0 ) THEN
              ISRTX = ISRTX + 20
              WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Mass Fraction'
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Aqueous Salt Option: '
     &          // BDUM(3)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(1)(1:),'sr2').NE.0 ) THEN
            ISRTX = 500
            WRITE(IWR,'(2X,2A)') 'Source State: #2'
!
!---        Read aqueous salt option  ---
!
            VARB = 'Aqueous Salt Option'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
            IF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'rel').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'sat').NE.0 ) THEN
              ISRTX = ISRTX + 20
              WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Relative Saturation'
            ELSEIF( INDEX(BDUM(3)(1:),'salt').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'frac').NE.0 ) THEN
              ISRTX = ISRTX + 30
              WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Mass Fraction'
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Aqueous Salt Option: '
     &          // BDUM(3)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(BDUM(1)(1:),'sr3').NE.0 ) THEN
             ISRTX = 600
            WRITE(IWR,'(2X,2A)') 'Source State: #3'
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Source State: ' // BDUM(1)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'sink').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Pressure Sink'
          ISRTX = 35
        ELSEIF( IEQC.NE.0 .AND. INDEX(ADUM(1:),'solute').NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(1))
          DO NSL = 1,NSOLU
            IDB = INDEX(SOLUT(NSL)(1:),'  ')
            IF( INDEX(BDUM(1)(1:),SOLUT(NSL)(1:IDB)).NE.0 ) THEN
              IF( INDEX(ADUM(1:),'density').NE.0 ) THEN
                ISRTX = -(NSL+NSOLU)
                WRITE(IWR,'(2X,2A)')'Solute Source Density: ',SOLUT(NSL)
              ELSE
                ISRTX = -NSL
                WRITE(IWR,'(2X,2A)')'Solute Source: ',SOLUT(NSL)
              ENDIF
              GOTO 32
            ENDIF
          ENDDO
          INDX = 4
          CHMSG = 'Unrecognized Solute Name: '//ADUM
          CALL WRMSGS( INDX )
   32     CONTINUE
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Fault Source Type: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read source domain indices  ---
!
        VARB = 'Fault Source Domain Indices'
        ISX = ISTART
        CALL RDINT(ISTART,ICOMMA,CHDUM,IF1X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,IT1X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,IF2X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,IT2X)
        ICX = ISTART
        WRITE(IWR,'(2X,A)') 'Fault Source Domain:'
        WRITE(IWR,'(4X,2(A,I6),A)') 'Fault = ',IF1X,
     &    ' Triangle = ',IT1X,' to '
        WRITE(IWR,'(4X,2(A,I6))') 'Fault = ',IF2X,
     &    ' Triangle = ',IT2X
!
!---    Check for ill-defined fault source domains  ---
!
        IF( IF1X.LT.1 .OR. IF1X.GT.NF_FRC .OR. IF2X.LT.1 .OR.
     &    IF2X.GT.NF_FRC .OR. IF2X.LT.IF1X ) THEN
          INDX = 4
          CHMSG = 'Invalid Fault Source Domain: ' // CHDUM(ISX:ICX)
          CALL WRMSGS( INDX )
        ENDIF
        IF( IT1X.LT.1 .OR. IT1X.GT.NTP_FRC(IF1X) .OR. IT2X.LT.1 .OR.
     &    IT2X.GT.NTP_FRC(IF2X) .OR.
     &    ( IF1X.EQ.IF2X .AND. IT2X.LT.IT1X ) ) THEN
          INDX = 4
          CHMSG = 'Invalid Fault Source Domain: ' // CHDUM(ISX:ICX)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for sources applied to inactive fault triangles  ---
!
        DO NFX = IF1X,IF2X
          DO NTX = IP_FRC(1,NFX)-1+IT1X,IP_FRC(1,NFX)-1+IT2X
            IF( IXP_FRC(NTX).EQ.0 ) THEN
              INDX = 7
              IMSG = NTX
              CHMSG = 'Source Applied to an Inactive Fault Triangle: '
              CALL WRMSGS( INDX )
            ENDIF
          ENDDO
        ENDDO
!
!---    Read number of fault source times  ---
!
        VARB = 'Number of Fault Source Times'
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISRM_FRC(NS))
        IF( ISRM_FRC(NS).GT.LSTM ) THEN
          INDX = 5
          CHMSG = 'Number of Fault Source Times > Parameter LSTM'
          CALL WRMSGS( INDX )
        ENDIF
        SRTMO = -SMALL
        DO NTM = 1,ISRM_FRC(NS)
          DO M = 1,8+NGC+NSOLU
            VAR(NTM,M) = 0.D+0
          ENDDO
!
!---     Read and write fault source values and units  ---
!
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Fault Source Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',VAR(NTM,1)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR(NTM,1),INDX)
!
!---      Power Source  ---
!
          IF( ISRTX.EQ.1 ) THEN
            VARB = 'Source: Power'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNM = 2
            IUNKG = 1
            IUNS = -3
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
!
!---      Power-Density Source  ---
!
          ELSEIF( ISRTX.EQ.2 ) THEN
            VARB = 'Source: Power Density'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -3
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
!
!---      Fluid Mass Source  ---
!
          ELSEIF( ISRTX/100.GE.1 .AND. ISRTX/100.LE.3 ) THEN
!
!---        Check for isothermal option, otherwise read fluid
!           temperature  ---
!
            IF( ISLC(30).EQ.0 ) THEN
              VARB = 'Source: Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,3)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,3),', C)'
            ENDIF
!
!---        Fluid mass rate  ---
!
            VARB = 'Source: Fluid Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            WRITE(IWR,'(2X,A,1PE11.4,A)') ' (',VAR(NTM,2),', kg/s)'
            IF( VAR(NTM,2).LT.0.D+0 ) GOTO 190
!
!---        Source state #1  ---
!
!           SL = 1.0
!           SG+SN = 0.0
!
!           Declared variables:
!
!           SRC(4,NTM,NSR) - aqueous CO2 relative saturation, or
!                            aqueous CO2 mass fraction
!           SRC(5,NTM,NSR) - aqueous salt relative saturation, or
!                            aqueous salt mass fraction
!
            IF( ISRTX/100.EQ.1 ) THEN
!
!---          Read aqueous CO2  ---
!
              IF( MOD(ISRTX,10).EQ.2 ) THEN
                VARB = 'Aqueous CO2 Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,4)
              ELSEIF( MOD(ISRTX,10).EQ.3 ) THEN
                VARB = 'Aqueous CO2 Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,4)
              ENDIF
!
!---          Read aqueous salt  ---
!
              IF( MOD(ISRTX/10,10).EQ.2 ) THEN
                VARB = 'Aqueous Salt Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,5)
              ELSEIF( MOD(ISRTX/10,10).EQ.3 ) THEN
                VARB = 'Aqueous Salt Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,5)
              ENDIF
!
!---        Source state #2  ---
!
!           SL < 1.0
!           SG+SN > 0.0
!
!           Declared variables:
!
!           SRC(4,NTM,NSR) - aqueous saturation
!           SRC(5,NTM,NSR) - aqueous salt relative saturation, or
!                            aqueous salt mass fraction
!           SRC(6,NTM,NSR) - total mole fraction of petroleum components
!           SRC(7+NGC,NTM,NSR) - total mole fraction of petroleum comp.
!
            ELSEIF( ISRTX/100.EQ.2 ) THEN
!
!---          Read aqueous saturation  ---
!
              VARB = 'Aqueous Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,4)
!
!---          Read aqueous salt  ---
!
              IF( MOD(ISRTX/10,10).EQ.2 ) THEN
                VARB = 'Aqueous Salt Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,5)
              ELSEIF( MOD(ISRTX/10,10).EQ.3 ) THEN
                VARB = 'Aqueous Salt Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,5)
              ENDIF
!
!---          Read total mole fraction of petroleum components  ---
!
              SZMCX = 0.D+0
              DO IGC = 1,NGC+2
                JGC = 5 + IGC
                VARB = 'Total Mole Fraction of ' // GCNM(IGC)
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,JGC))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,JGC)
                SZMCX = SZMCX + VAR(NTM,JGC)
              ENDDO
              IF( ABS(1.D+0-SZMCX).GT.EPSL ) THEN
                INDX = 9
                RLMSG = SZMCX
                CHMSG = 'Summation of Total Mole Fraction ' // 
     &            'of Petroleum Components'
                CALL WRMSGS( INDX )
              ELSE
                DO IGC = 1,NGC+2
                  VAR(NTM,JGC) = VAR(NTM,JGC)/SZMCX
                ENDDO
              ENDIF
!
!---        Source state #3  ---
!
!           SL = 0.0
!           SG+SN = 1.0
!
!           Declared variables:
!
!           SRC(5,NTM,NSR) - water vapor relative saturation
!           SRC(6+,NTM,NSR) - total mole fraction of petroleum comp.
!          
            ELSEIF( ISRTX/100.EQ.3 ) THEN
!
!---          Read water vapor relative saturation  ---
!
              VARB = 'Water Vapor Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
!
!---          Read total mole fraction of petroleum components  ---
!
              SZMCX = 0.D+0
              DO IGC = 1,NGC+2
                JGC = 5 + IGC
                VARB = 'Total Mole Fraction of ' // GCNM(IGC)
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,JGC))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,JGC)
                SZMCX = SZMCX + VAR(NTM,JGC)
              ENDDO
              IF( ABS(1.D+0-SZMCX).GT.EPSL ) THEN
                INDX = 9
                RLMSG = SZMCX
                CHMSG = 'Summation of Total Mole Fraction ' // 
     &            'of Petroleum Components'
                CALL WRMSGS( INDX )
              ELSE
                DO IGC = 1,NGC+2
                  VAR(NTM,JGC) = VAR(NTM,JGC)/SZMCX
                ENDDO
              ENDIF
            ENDIF
!
!---      Fluid Volumetric Source  ---
!
          ELSEIF( ISRTX/100.GE.4 .AND. ISRTX/100.LE.6 ) THEN
            VARB = 'Source: Fluid Volumetric Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            WRITE(IWR,'(2X,A,1PE11.4,A)') ' (',VAR(NTM,2),', m^3/s)'
            IF( VAR(NTM,2).LT.0.D+0 ) GOTO 190
!
!---        Isothermal option  ---
!
            IF( ISLC(30).EQ.0 ) THEN
              VARB = 'Source: Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,3)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,3),', C)'
            ENDIF
!
!---        Source state #1  ---
!
!           SL = 1.0
!           SG+SN = 0.0
!
!           Declared variables:
!
!           SRC(4,NTM,NSR) - aqueous CO2 relative saturation, or
!                            aqueous CO2 mass fraction
!           SRC(5,NTM,NSR) - aqueous salt relative saturation, or
!                            aqueous salt mass fraction
!
            IF( ISRTX/100.EQ.4 ) THEN
!
!---          Read aqueous CO2  ---
!
              IF( MOD(ISRTX,10).EQ.2 ) THEN
                VARB = 'Aqueous CO2 Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,4)
              ELSEIF( MOD(ISRTX,10).EQ.3 ) THEN
                VARB = 'Aqueous CO2 Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,4)
              ENDIF
!
!---          Read aqueous salt  ---
!
              IF( MOD(ISRTX/10,10).EQ.2 ) THEN
                VARB = 'Aqueous Salt Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,5)
              ELSEIF( MOD(ISRTX/10,10).EQ.3 ) THEN
                VARB = 'Aqueous Salt Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,5)
              ENDIF
!
!---        Source state #2  ---
!
!           SL < 1.0
!           SG+SN > 0.0
!
!           Declared variables:
!
!           SRC(4,NTM,NSR) - aqueous saturation
!           SRC(5,NTM,NSR) - aqueous salt relative saturation, or
!                            aqueous salt mass fraction
!           SRC(6,NTM,NSR) - total mole fraction of petroleum components
!           SRC(7+NGC,NTM,NSR) - total mole fraction of petroleum comp.
!
            ELSEIF( ISRTX/100.EQ.5 ) THEN
!
!---          Read aqueous saturation  ---
!
              VARB = 'Aqueous Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,4)
!
!---          Read aqueous salt  ---
!
              IF( MOD(ISRTX/10,10).EQ.2 ) THEN
                VARB = 'Aqueous Salt Relative Saturation'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,5)
              ELSEIF( MOD(ISRTX/10,10).EQ.3 ) THEN
                VARB = 'Aqueous Salt Mass Fraction'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,5)
              ENDIF
!
!---          Read total mole fraction of petroleum components  ---
!
              SZMCX = 0.D+0
              DO IGC = 1,NGC+2
                JGC = 5 + IGC
                VARB = 'Total Mole Fraction of ' // GCNM(IGC)
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,JGC))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,JGC)
                SZMCX = SZMCX + VAR(NTM,JGC)
              ENDDO
              IF( ABS(1.D+0-SZMCX).GT.EPSL ) THEN
                INDX = 9
                RLMSG = SZMCX
                CHMSG = 'Summation of Total Mole Fraction ' // 
     &            'of Petroleum Components'
                CALL WRMSGS( INDX )
              ELSE
                DO IGC = 1,NGC+2
                  VAR(NTM,JGC) = VAR(NTM,JGC)/SZMCX
                ENDDO
              ENDIF
!
!---        Source state #3  ---
!
!           SL = 0.0
!           SG+SN = 1.0
!
!           Declared variables:
!
!           SRC(5,NTM,NSR) - water vapor relative saturation
!           SRC(6+,NTM,NSR) - total mole fraction of petroleum comp.
!          
            ELSEIF( ISRTX/100.EQ.6 ) THEN
!
!---          Read water vapor relative saturation  ---
!
              VARB = 'Water Vapor Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
!
!---          Read total mole fraction of petroleum components  ---
!
              SZMCX = 0.D+0
              DO IGC = 1,NGC+2
                JGC = 5 + IGC
                VARB = 'Total Mole Fraction of ' // GCNM(IGC)
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,JGC))
                WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &            ': ',VAR(NTM,JGC)
                SZMCX = SZMCX + VAR(NTM,JGC)
              ENDDO
              IF( ABS(1.D+0-SZMCX).GT.EPSL ) THEN
                INDX = 9
                RLMSG = SZMCX
                CHMSG = 'Summation of Total Mole Fraction ' // 
     &            'of Petroleum Components'
                CALL WRMSGS( INDX )
              ELSE
                DO IGC = 1,NGC+2
                  VAR(NTM,JGC) = VAR(NTM,JGC)/SZMCX
                ENDDO
              ENDIF
            ENDIF
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
            VARB = 'Borehole Radius'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Borehole Skin Factor'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,2)
!
!---      Solute Source  ---
!
          ELSEIF( ISRTX.LT.0 .AND. ISRTX.GE.-NSOLU ) THEN
            VARB = 'Source Solute Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ELSEIF( ISRTX.LT.-NSOLU .AND. ISRTX.GE.-2*NSOLU ) THEN
            VARB = 'Source Solute Density Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNS = -1
            IUNM = -3
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ENDIF
  190     CONTINUE
!
!---      Check for nonascending source times  ---
!
          IF( VAR(NTM,1).LT.SRTMO ) THEN
            INDX = 4
            CHMSG = 'Source Time Sequencing'
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
          CHMSG = 'Number of Fault Sources > Parameter LSR_FRC'
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
          DO M = 1,8+NGC+NSOLU
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
!---  End of RDSR_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_FRC_BO
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
!     Compute the maximum relative residuals for the black-oil option
!     for faults
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 20 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE JACOB
      USE HYST
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FILES
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP
      USE FDVGC_FRC
      USE FDVGC
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
      CHARACTER*128 PH_CND(5)
      REAL*8 DPC(LNGC)
      INTEGER MPC(LNGC)
!
!----------------------Data Statements---------------------------------!
!
      DATA PH_CND /'Phase Condition #1, SL = 1.0, SG+SN = 0.0',
     &  'Phase Condition #2, SL < 1.0, ST > 0.0, SG = 0.0',
     &  'Phase Condition #3, SL < 1.0, ST > 0.0, SG > 0.0',
     &  'Phase Condition #4, SL = 0.0, ST = 1.0, SG = 0.0',
     &  'Phase Condition #5, SL = 0.0, ST = 1.0, SG > 0.0'/
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL_BO_EOR'
!
!---  Black-oil option  ---
!
!
!---  Zero maximum residuals  ---
!
      DO M = 1,ISVC
        RSD_FRC(M) = 0.D+0
        NSD_FRC(M) = 0
      ENDDO
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
        N_DB = NTX
        IZN = IZ_FRC(NTX)
!
!---    Water mass equation  ---
!
        MPW = IM_FRC(IEQW,NMD)
        DPW = BLU(MPW)
!
!---    Gas mass equation  ---
!
        MPC(1) = IM_FRC(IEQGC(1),NMD)
        DPC(1) = BLU(MPC(1))
!
!---    Oil mass equation  ---
!
        MPC(2) = IM_FRC(IEQGC(2),NMD)
        DPC(2) = BLU(MPC(2))
!
!---    Salt mass equation  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = IM_FRC(IEQS,NMD)
          DPS = BLU(MPS)
        ELSE
          DPS = 0.D+0
        ENDIF
!
!---    Skip selected nodes in the residual calculation  ---
!
        IF( ISKP(IZN).EQ.1 ) CYCLE
!
!---    Absolute system pressure  ---
!
        PX = PSO_FRC(2,NTX) + PATM
!
!---    Phase condition #1  ---
!
!       SL = 1.0
!       SG+SN = 0.0
!       Aqueuous saturated
!
!       water mass - aqueous (system) pressure, Pa gauge - PL,PSO
!       gas mass - total gas mass - TMC_FRC(1,2,NTX)
!       oil mass - total oil mass - TMC_FRC(2,2,NTX)
!       salt mass - total salt brine mass fraction - YLS
!
        IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
!
!---      Water mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &      XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGW_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPW)/PX,
     &      ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD_FRC(IEQW) = RSDX
            NSD_FRC(IEQW) = NTX
          ENDIF
!
!---      Gas mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGC_FRC(1,2,NTX) +RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &      XNC_FRC(1,2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPC(1))/PX,
     &      ABS(RSDL_FRC(IEQGC(1),NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQGC(1)) ) THEN
            RSD_FRC(IEQGC(1)) = RSDX
            NSD_FRC(IEQGC(1)) = NTX
          ENDIF
!
!---      Oil mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &      XNC_FRC(2,2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPC(2))/PX,
     &      ABS(RSDL_FRC(IEQGC(2),NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQGC(2)) ) THEN
            RSD_FRC(IEQGC(2)) = RSDX
            NSD_FRC(IEQGC(2)) = NTX
          ENDIF
!
!---      Salt mass equation residual  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS_FRC(2,NTX)*DTI*VOL_FRC(NTX)
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            RSDX = MIN( (ABS(DPS)/XLSMX),
     &        ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            IF( RSDX.GT.RSD(IEQS) ) THEN
              RSD_FRC(IEQS) = RSDX
              NSD_FRC(IEQS) = NTX
            ENDIF
          ENDIF
!
!---    Phase condition #2  ---
!
!       SG = 0.0
!       SL > 0.0
!       SN > 0.0
!       Aqueous and nonaqueous liquid saturated
!
!       water mass - aqueous saturation, SL
!       gas mass - gas to oil ratio, POSM
!       oil mass - nonaqueous (system) pressure, Pa gauge - PN,PSO
!       salt mass - total salt brine mass fraction - YLS
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
!
!---      Water mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &      XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGW_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPW),
     &      ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD_FRC(IEQW) = RSDX
            NSD_FRC(IEQW) = NTX
          ENDIF
!
!---      Gas mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGC_FRC(1,2,NTX) + RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &      XNC_FRC(1,2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPC(1)),
     &      ABS(RSDL_FRC(IEQGC(1),NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQGC(1)) ) THEN
            RSD_FRC(IEQGC(1)) = RSDX
            NSD_FRC(IEQGC(1)) = NTX
          ENDIF
!
!---      Oil mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &      XNC_FRC(2,2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPC(2))/PX,
     &      ABS(RSDL_FRC(IEQGC(2),NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQGC(2)) ) THEN
            RSD_FRC(IEQGC(2)) = RSDX
            NSD_FRC(IEQGC(2)) = NTX
          ENDIF
!
!---      Salt mass equation residual  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS_FRC(2,NTX)*DTI*VOL_FRC(NTX)
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            RSDX = MIN( (ABS(DPS)/XLSMX),
     &        ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            IF( RSDX.GT.RSD(IEQS) ) THEN
              RSD_FRC(IEQS) = RSDX
              NSD_FRC(IEQS) = NTX
            ENDIF
          ENDIF
!
!---    Phase condition #3  ---
!
!       SG > 0.0
!       SL > 0.0
!       SN > 0.0
!       Three-phase
!
!       water mass - aqueous saturation, SL
!       gas mass - gas saturation, SG
!       oil mass - gas (system) pressure, Pa gauge - PG,PSO
!       salt mass - total salt brine mass fraction - YLS
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
!
!---      Water mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &      XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGW_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPW),
     &      ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD_FRC(IEQW) = RSDX
            NSD_FRC(IEQW) = NTX
          ENDIF
!
!---      Gas mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGC_FRC(1,2,NTX) + RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &      XNC_FRC(1,2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPC(1))/PX,
     &      ABS(RSDL_FRC(IEQGC(1),NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQGC(1)) ) THEN
            RSD_FRC(IEQGC(1)) = RSDX
            NSD_FRC(IEQGC(1)) = NTX
          ENDIF
!
!---      Oil mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &      XNC_FRC(2,2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPC(2)),
     &      ABS(RSDL_FRC(IEQGC(2),NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQGC(2)) ) THEN
            RSD_FRC(IEQGC(2)) = RSDX
            NSD_FRC(IEQGC(2)) = NTX
          ENDIF
!
!---      Salt mass equation residual  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS_FRC(2,NTX)*DTI*VOL_FRC(NTX)
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            RSDX = MIN( (ABS(DPS)/XLSMX),
     &        ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            IF( RSDX.GT.RSD(IEQS) ) THEN
              RSD_FRC(IEQS) = RSDX
              NSD_FRC(IEQS) = NTX
            ENDIF
          ENDIF
!
!---    Phase condition #4  ---
!
!       SG = 0.0
!       SL = 0.0
!       SN = 1.0
!       Nonaqueous liquid saturated
!
!       water mass - aqueous saturation, SL
!       gas mass - gas to oil ratio, POSM
!       oil mass - nonaqueous (system) pressure, Pa gauge - PN,PSO
!       salt mass - total salt mass - TMS_FRC(2,NTX)
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
!
!---      Water mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &      XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGW_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPW),
     &      ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD_FRC(IEQW) = RSDX
            NSD_FRC(IEQW) = NTX
          ENDIF
!
!---      Gas mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGC_FRC(1,2,NTX) + RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &      XNC_FRC(1,2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPC(1)),
     &      ABS(RSDL_FRC(IEQGC(1),NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQGC(1)) ) THEN
            RSD_FRC(IEQGC(1)) = RSDX
            NSD_FRC(IEQGC(1)) = NTX
          ENDIF
!
!---      Oil mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &      XNC_FRC(2,2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPC(2))/PX,
     &      ABS(RSDL_FRC(IEQGC(2),NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQGC(2)) ) THEN
            RSD_FRC(IEQGC(2)) = RSDX
            NSD_FRC(IEQGC(2)) = NTX
          ENDIF
!
!---      Salt mass equation residual  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS_FRC(2,NTX)*DTI*VOL_FRC(NTX)
            RSDX = MIN( ABS(DPS),
     &        ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            IF( RSDX.GT.RSD(IEQS) ) THEN
              RSD_FRC(IEQS) = RSDX
              NSD_FRC(IEQS) = NTX
            ENDIF
          ENDIF
!
!---    Phase condition #5  ---
!
!       SG > 0.0
!       SL = 0.0
!       SN > 0.0
!       Nonaqueous liquid and gas saturated
!
!       water mass - aqueous saturation, SL
!       gas mass - gas saturation, SG
!       oil mass - gas (system) pressure, Pa gauge - PG,PSO
!       salt mass - total salt mass - TMS_FRC(2,NTX)
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
!
!---      Water mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &      XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGW_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPW),
     &      ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD_FRC(IEQW) = RSDX
            NSD_FRC(IEQW) = NTX
          ENDIF
!
!---      Gas mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGC_FRC(1,2,NTX) + RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &      XNC_FRC(1,2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPC(1))/PX,
     &      ABS(RSDL_FRC(IEQGC(1),NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQGC(1)) ) THEN
            RSD_FRC(IEQGC(1)) = RSDX
            NSD_FRC(IEQGC(1)) = NTX
          ENDIF
!
!---      Oil mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &      XNC_FRC(2,2,NTX))*DTI*VOL_FRC(NTX)
          RSDX = MIN( ABS(DPC(2)),
     &      ABS(RSDL_FRC(IEQGC(2),NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQGC(2)) ) THEN
            RSD_FRC(IEQGC(2)) = RSDX
            NSD_FRC(IEQGC(2)) = NTX
          ENDIF
!
!---      Salt mass equation residual  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS_FRC(2,NTX)*DTI*VOL_FRC(NTX)
            RSDX = MIN( ABS(DPS),
     &        ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            IF( RSDX.GT.RSD(IEQS) ) THEN
              RSD_FRC(IEQS) = RSDX
              NSD_FRC(IEQS) = NTX
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      ENDDO
!
!---  Assign a convergence index  ---
!
      RSDX = 0.D+0
      DO M = 1,ISVC
        RSDX = MAX( RSD_FRC(M),RSDX )
        IF( RSD_FRC(M).GT.RSDMX ) ICNV = 2
      ENDDO
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) ICNV = 1
      IF( RSDX.GE.1.D+2 ) ICNV = 1
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.1 ) THEN
        IF( RSDX.GE.1.D+2 ) THEN
          WRITE(ISC,'(10X,A)') '---  Excessive Residual in Fault  ---'
          WRITE(IWR,'(10X,A)') '---  Excessive Residual in Fault  ---'
        ELSE
          WRITE(ISC,'(10X,A)') '---  Convergence Failure in Fault  ---'
          WRITE(IWR,'(10X,A)') '---  Convergence Failure in Fault  ---'
        ENDIF
!
!---    Water mass equation  ---
!
        NTX = NSD_FRC(IEQW)
        IF( NTX.GT.0 ) THEN
          NPX = NPHAZ_FRC(2,NTX)
          NCHX = INDEX( PH_CND(NPX),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Mass Equation Maximum Residual = ',RSD_FRC(IEQW),
     &      ': Fault Triangle = ',NTX,
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Mass Equation Maximum Residual = ',RSD_FRC(IEQW),
     &      ': Fault Triangle = ',NTX,
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
        ENDIF
!
!---    Gas mass equation  ---
!
        NTX = NSD_FRC(IEQGC(1))
        IF( NTX.GT.0 ) THEN
          NPX = NPHAZ_FRC(2,NTX)
          NCAX = INDEX( PH_CND(NPX),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Gas Mass Equation Maximum Residual = ',
     &      RSD_FRC(IEQGC(1)),': Fault Triangle = ',NTX,
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCAX)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Gas Mass Equation Maximum Residual = ',
     &      RSD_FRC(IEQGC(1)),': Fault Triangle = ',NTX,
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCAX)
        ENDIF
!
!---    Oil mass equation  ---
!
        NTX = NSD_FRC(IEQGC(2))
        IF( NTX.GT.0 ) THEN
          NPX = NPHAZ_FRC(2,NTX)
          NCAX = INDEX( PH_CND(NPX),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Oil Mass Equation Maximum Residual = ',
     &      RSD_FRC(IEQGC(2)),': Fault Triangle = ',NTX,
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCAX)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Oil Mass Equation Maximum Residual = ',
     &      RSD_FRC(IEQGC(2)),': Fault Triangle = ',NTX,
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCAX)
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          NTX = NSD_FRC(IEQS)
          IF( NTX.GT.0 ) THEN
            NPX = MOD(NPHAZ_FRC(2,NTX),100)+1
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD_FRC(IEQS),
     &        ': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD_FRC(IEQS),
     &         ': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---  Reduce time step  ---
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
!---      Loop over faults  ---
!
          DO NFX = 1,NF_FRC
!
!---      Loop over fault triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
            DO M = 1,6
              BETA_FRC(M,NTX) = 0.D+0
            ENDDO
            IBETA_FRC(NTX) = 0
            T_FRC(2,NTX) = T_FRC(1,NTX)
            PL_FRC(2,NTX) = PL_FRC(1,NTX)
            PG_FRC(2,NTX) = PG_FRC(1,NTX)
            PN_FRC(2,NTX) = PN_FRC(1,NTX)
            PSO_FRC(2,NTX) = PSO_FRC(1,NTX)
            SG_FRC(2,NTX) = SG_FRC(1,NTX)
            SL_FRC(2,NTX) = SL_FRC(1,NTX)
            SN_FRC(2,NTX) = SN_FRC(1,NTX)
            YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
            TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
            NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
            DO IGC = 1,NGC+2
              TMC_FRC(IGC,2,NTX) = TMC_FRC(IGC,1,NTX)
            ENDDO
          ENDDO
          ENDDO
!
!---      Loop over field nodes  ---
!
          DO N = 1,NFBN
            IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
            DO M = 1,6
              BETA(M,N) = 0.D+0
            ENDDO
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            POSM(2,N) = POSM(1,N)
            PSO(2,N) = PSO(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
            DO IGC = 1,NGC+2
              TMC(IGC,2,N) = TMC(IGC,1,N)
            ENDDO
          ENDDO
!
!---      Coupled-well pressure  ---
!
          DO NCW = 1,N_CW
            P_CW(2,NCW) = P_CW(1,NCW)
          ENDDO
          NITER = 1
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
!---  End of RSDL_BO_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_FRC_EOR
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
!     Compute the maximum fault relative residuals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 20 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE PARM_FRC
      USE OUTPU
      USE JACOB
      USE HYST
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FILES
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP
      USE FDVGC_FRC
      USE FDVGC
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
      CHARACTER*128 PH_CND(5)
      REAL*8 DPC(LNGC)
      INTEGER MPC(LNGC)
!
!----------------------Data Statements---------------------------------!
!
      DATA PH_CND /'Phase Condition #1, SL = 1.0, SG+SN = 0.0',
     &  'Phase Condition #2, SL < 1.0, ST > 0.0, SG = 0.0',
     &  'Phase Condition #3, SL < 1.0, ST > 0.0, SG > 0.0',
     &  'Phase Condition #4, SL = 0.0, ST = 1.0, SG = 0.0',
     &  'Phase Condition #5, SL = 0.0, ST = 1.0, SG > 0.0'/
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL_FRC_EOR'
!
!---  Black-oil option  ---
!
      IF( ISLC(68).EQ.0 ) THEN
        CALL RSDL_FRC_BO
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compositional option  ---
!
 !
!---  Zero maximum residuals  ---
!
      DO M = 1,ISVC
        RSD_FRC(M) = 0.D+0
        NSD_FRC(M) = 0
        RSDAVG_FRC(M) = 0.D+0
      ENDDO
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
        N_DB = NTX
        IZN = IZ_FRC(NTX)
!
!---    Energy equation  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          MPT = IM_FRC(IEQT,NTX)
          DPT = BLU(MPT)
        ELSE
          DPT = 0.D+0
        ENDIF
!
!---    Water mass equation  ---
!
        MPW = IM_FRC(IEQW,NTX)
        DPW = BLU(MPW)
!
!---    CO2 mass equation  ---
!
        IF( ISLC(45).EQ.0 ) THEN
          MPC(1) = IM_FRC(IEQGC(1),NTX)
          DPC(1) = BLU(MPC(1))
        ELSE
          DPC(1) = 0.D+0
        ENDIF
!
!---    CH4 mass equation  ---
!
        IF( ISLC(64).EQ.0 ) THEN
          MPC(2) = IM_FRC(IEQGC(2),NTX)
          DPC(2) = BLU(MPC(2))
        ELSE
          DPC(2) = 0.D+0
        ENDIF
!
!---    Petroleum component mass equation  ---
!
        DO IGC = 3,NGC+2
          MPC(IGC) = IM_FRC(IEQGC(IGC),NTX)
          DPC(IGC) = BLU(MPC(IGC))
        ENDDO
!
!---    Salt mass equation  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = IM_FRC(IEQS,NTX)
          DPS = BLU(MPS)
        ELSE
          DPS = 0.D+0
        ENDIF
!
!---    Total fluid and salt mass density, kg/m^3 rock volume  ---
!
        TMX = PORD_FRC(2,NTX)*(SL_FRC(2,NTX)*RHOL_FRC(2,NTX) + 
     &    SG_FRC(2,NTX)*RHOG_FRC(2,NTX) + 
     &    SN_FRC(2,NTX)*RHON_FRC(2,NTX) + 
     &    SS_FRC(2,NTX)*RHOSP_FRC(2,NTX))
!
!---    Skip selected nodes in the residual calculation  ---
!
        IF( ISKP(IZN).EQ.1 ) CYCLE
        IF( ZMC_FRC(1,2,NTX).GT.9.D-1 ) CYCLE
!
!---    Absolute system pressure  ---
!
        PX = PSO_FRC(2,NTX) + PATM
!
!---    Phase condition #1  ---
!
!       SL = 1.0, SG+SN = 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous (system) pressure, Pa gauge - PL,PSO
!       CO2 mass - CO2 vapor pressure - PVA
!       CH4 mass - total CH4 mass - TMC
!       petroleum component mass - total petroleum component mass - TMC
!       salt mass - total salt brine mass fraction - YLS
!
        IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
!
!---      Energy equation residual  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            RSDX = 1.D-2*ABS(DPT)/TABS
            RSDAVG_FRC(IEQT) = RSDAVG_FRC(IEQT) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQT) ) THEN
              RSD_FRC(IEQT) = RSDX
              NSD_FRC(IEQT) = NTX
            ENDIF
          ENDIF
!
!---      Water mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &      XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGW_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          ACP = TMX*DTI*VOL_FRC(NTX)
          RSDAVG_FRC(IEQW) = RSDAVG_FRC(IEQW) + RSDX
          RSDX = 1.D-2*MIN( ABS(DPW)/PX,
     &      ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD_FRC(IEQW) ) THEN
            RSD_FRC(IEQW) = RSDX
            NSD_FRC(IEQW) = NTX
          ENDIF
!
!---      CO2 mass equation residual  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( PVA_FRC(2,NTX).GT.(1.D-6*PX) ) THEN
              ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &          XGA_FRC(2,NTX) + RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &          XLA_FRC(2,NTX) + RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &          XNC_FRC(1,2,NTX))*DTI*VOL_FRC(NTX)
              ACP = TMX*DTI*VOL_FRC(NTX)
              RSDAVG_FRC(IEQGC(1)) = RSDAVG_FRC(IEQGC(1)) + RSDX
              RSDX = 1.D-2*MIN( ABS(DPC(1))/PX,
     &          1.D-1*ABS(RSDL_FRC(IEQGC(1),NTX)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD_FRC(IEQGC(1)) ) THEN
                RSD_FRC(IEQGC(1)) = RSDX
                NSD_FRC(IEQGC(1)) = NTX
              ENDIF
            ENDIF
          ENDIF
!
!---      CH4 mass equation residual  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &        XGO_FRC(2,NTX) + RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &        XNC_FRC(2,2,NTX))*DTI*VOL_FRC(NTX)
            ACP = TMX*DTI*VOL_FRC(NTX)
            RSDAVG_FRC(IEQGC(2)) = RSDAVG_FRC(IEQGC(2)) + RSDX
            RSDX = 1.D-2*MIN( ABS(DPC(2))/PX,
     &        1.D-1*ABS(RSDL_FRC(IEQGC(2),NTX)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD_FRC(IEQGC(2)) ) THEN
              RSD_FRC(IEQGC(2)) = RSDX
              NSD_FRC(IEQGC(2)) = NTX
            ENDIF
          ENDIF
!
!---      Petroleum component mass equation residual  ---
!
          DO IGC = 3,NGC+2
            ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &        XGC_FRC(IGC,2,NTX) + RHON_FRC(2,NTX)*SN_FRC(2,NTX)*
     &        XNC_FRC(IGC,2,NTX))*DTI*VOL_FRC(NTX)
            ACP = TMX*DTI*VOL_FRC(NTX)
            RSDX = 1.D-2*MIN( ABS(DPC(IGC))/PX,
     &        1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
            RSDAVG_FRC(IEQGC(IGC)) = RSDAVG_FRC(IEQGC(IGC)) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQGC(IGC)) ) THEN
              RSD_FRC(IEQGC(IGC)) = RSDX
              NSD_FRC(IEQGC(IGC)) = NTX
            ENDIF
          ENDDO
!
!---      Salt mass equation residual  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS_FRC(2,NTX)*DTI*VOL_FRC(NTX)
            ACP = TMX*DTI*VOL_FRC(NTX)
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            RSDX = 1.D-2*MIN( (ABS(DPS)/XLSMX),
     &        1.D-1*ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            RSDAVG_FRC(IEQS) = RSDAVG_FRC(IEQS) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQS) ) THEN
              RSD_FRC(IEQS) = RSDX
              NSD_FRC(IEQS) = NTX
            ENDIF
          ENDIF
!
!---    Phase condition #2  ---
!
!       SL < 1.0, ST > 0.0, SG = 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!         endif
!       CH4 mass
!         if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!         endif
!       petroleum component mass
!         if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!         endif
!       salt mass - total salt brine mass fraction - YLS
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
!
!---      Energy equation residual  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            RSDX = ABS(DPT)/TABS
            RSDAVG_FRC(IEQT) = RSDAVG_FRC(IEQT) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQT) ) THEN
              RSD_FRC(IEQT) = RSDX
              NSD_FRC(IEQT) = NTX
            ENDIF
          ENDIF
!
!---      Water mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &      XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGW_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          ACP = TMX*DTI*VOL_FRC(NTX)
          IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
            RSDX = MIN( 1.D-1*ABS(DPW)/PX,
     &        1.D-1*ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          ELSE
            RSDX = MIN( ABS(DPW),
     &        1.D-1*ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          ENDIF
          RSDAVG_FRC(IEQW) = RSDAVG_FRC(IEQW) + RSDX
          IF( RSDX.GT.RSD_FRC(IEQW) ) THEN
            RSD_FRC(IEQW) = RSDX
            NSD_FRC(IEQW) = NTX
          ENDIF
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX) +
     &      RHON_FRC(2,NTX)*SN_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          ACP = TMX*DTI*VOL_FRC(NTX)
          DO IGC = 1,NGC+2
            IF( ISLC(45).NE.0 .AND. IGC.EQ.1 ) CYCLE
            IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
!
!---        Nonaqueous-liquid pressure  ---
!
            IF( IGC.EQ.IZMC_FRC(NTX) ) THEN
              IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                RSDX = MIN( ABS(DPC(IGC)),
     &            1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
              ELSE
                RSDX = MIN( 1.D-1*ABS(DPC(IGC))/PX,
     &            1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
              ENDIF
!
!---        Petroleum component mole fraction  ---
!
            ELSE
              RSDX = MIN( 1.D-2*ABS(DPC(IGC)),
     &          1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
            ENDIF
            RSDAVG_FRC(IEQGC(IGC)) = RSDAVG_FRC(IEQGC(IGC)) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQGC(IGC)) ) THEN
              RSD_FRC(IEQGC(IGC)) = RSDX
              NSD_FRC(IEQGC(IGC)) = NTX
            ENDIF
          ENDDO
!
!---      Salt mass equation residual  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS_FRC(2,NTX)*DTI*VOL_FRC(NTX)
            ACP = TMX*DTI*VOL_FRC(NTX)
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            RSDX = MIN( (ABS(DPS)/XLSMX),
     &        1.D-1*ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            RSDAVG_FRC(IEQS) = RSDAVG_FRC(IEQS) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQS) ) THEN
              RSD_FRC(IEQS) = RSDX
              NSD_FRC(IEQS) = NTX
            ENDIF
          ENDIF
!
!---    Phase condition #3  ---
!
!       SL < 1.0, ST > 0.0, SG > 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!         endif
!       CH4 mass
!         if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!         endif
!       petroleum component mass
!         if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!         endif
!       salt mass - total salt brine mass fraction - YLS
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
!
!---      Energy equation residual  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            RSDX = ABS(DPT)/TABS
            RSDAVG_FRC(IEQT) = RSDAVG_FRC(IEQT) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQT) ) THEN
              RSD_FRC(IEQT) = RSDX
              NSD_FRC(IEQT) = NTX
            ENDIF
          ENDIF
!
!---      Water mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &      XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGW_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          ACP = TMX*DTI*VOL_FRC(NTX)
          IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
            RSDX = MIN( 1.D-1*ABS(DPW)/PX,
     &        1.D-1*ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          ELSE
            RSDX = MIN( ABS(DPW),
     &        1.D-1*ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          ENDIF
          RSDAVG_FRC(IEQW) = RSDAVG_FRC(IEQW) + RSDX
          IF( RSDX.GT.RSD_FRC(IEQW) ) THEN
            RSD_FRC(IEQW) = RSDX
            NSD_FRC(IEQW) = NTX
          ENDIF
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX) +
     &      RHON_FRC(2,NTX)*SN_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          ACP = TMX*DTI*VOL_FRC(NTX)
          DO IGC = 1,NGC+2
            IF( ISLC(45).NE.0 .AND. IGC.EQ.1 ) CYCLE
            IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
!
!---        Nonaqueous-liquid pressure  ---
!
            IF( IGC.EQ.IZMC_FRC(NTX) ) THEN
              IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                RSDX = MIN( 1.D-1*ABS(DPC(IGC))/PX,
     &            1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
              ELSE
                RSDX = MIN( 1.D-1*ABS(DPC(IGC))/PX,
     &            1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
              ENDIF
!
!---        Petroleum component mole fraction  ---
!
            ELSE
              RSDX = MIN( 1.D-2*ABS(DPC(IGC)),
     &          1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
            ENDIF
            RSDAVG_FRC(IEQGC(IGC)) = RSDAVG_FRC(IEQGC(IGC)) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQGC(IGC)) ) THEN
              RSD_FRC(IEQGC(IGC)) = RSDX
              NSD_FRC(IEQGC(IGC)) = NTX
            ENDIF
          ENDDO
!
!---      Salt mass equation residual  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS_FRC(2,NTX)*DTI*VOL_FRC(NTX)
            ACP = TMX*DTI*VOL_FRC(NTX)
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            RSDX = MIN( (ABS(DPS)/XLSMX),
     &        1.D-1*ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            RSDAVG_FRC(IEQS) = RSDAVG_FRC(IEQS) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQS) ) THEN
              RSD_FRC(IEQS) = RSDX
              NSD_FRC(IEQS) = NTX
            ENDIF
          ENDIF
!
!---    Phase condition #4  ---
!
!       SL = 0.0, ST = 1.0, SG = 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!         endif
!       CH4 mass
!         if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!         endif
!       petroleum component mass
!         if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!         endif
!       salt mass - total salt mass - TMS_FRC(2,NTX)
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
!
!---      Energy equation residual  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            RSDX = ABS(DPT)/TABS
            RSDAVG_FRC(IEQT) = RSDAVG_FRC(IEQT) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQT) ) THEN
              RSD_FRC(IEQT) = RSDX
              NSD_FRC(IEQT) = NTX
            ENDIF
          ENDIF
!
!---      Water mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &      XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGW_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          ACP = TMX*DTI*VOL_FRC(NTX)
          IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
            RSDX = MIN( 1.D-1*ABS(DPW)/PX,
     &        1.D-1*ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          ELSE
            RSDX = MIN( ABS(DPW),
     &        1.D-1*ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          ENDIF
          RSDAVG_FRC(IEQW) = RSDAVG_FRC(IEQW) + RSDX
          IF( RSDX.GT.RSD_FRC(IEQW) ) THEN
            RSD_FRC(IEQW) = RSDX
            NSD_FRC(IEQW) = NTX
          ENDIF
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX) +
     &      RHON_FRC(2,NTX)*SN_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          ACP = TMX*DTI*VOL_FRC(NTX)
          DO IGC = 1,NGC+2
            IF( ISLC(45).NE.0 .AND. IGC.EQ.1 ) CYCLE
            IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
!
!---        Nonaqueous-liquid pressure  ---
!
            IF( IGC.EQ.IZMC_FRC(NTX) ) THEN
              IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                RSDX = MIN( ABS(DPC(IGC)),
     &            1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
              ELSE
                RSDX = MIN( 1.D-1*ABS(DPC(IGC))/PX,
     &            1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
              ENDIF
!
!---        Petroleum component mole fraction  ---
!
            ELSE
              RSDX = MIN( 1.D-2*ABS(DPC(IGC)),
     &          1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
            ENDIF
            RSDAVG_FRC(IEQGC(IGC)) = RSDAVG_FRC(IEQGC(IGC)) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQGC(IGC)) ) THEN
              RSD_FRC(IEQGC(IGC)) = RSDX
              NSD_FRC(IEQGC(IGC)) = NTX
            ENDIF
          ENDDO
!
!---      Salt mass equation residual  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS_FRC(2,NTX)*DTI*VOL_FRC(NTX)
            ACP = TMX*DTI*VOL_FRC(NTX)
            RSDX = MIN( ABS(DPS),
     &        1.D-1*ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            RSDAVG_FRC(IEQS) = RSDAVG_FRC(IEQS) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQS) ) THEN
              RSD_FRC(IEQS) = RSDX
              NSD_FRC(IEQS) = NTX
            ENDIF
          ENDIF
!
!---    Phase condition #5  ---
!
!       SL = 0.0, ST >= 0.0, SG > 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!         endif
!       CH4 mass
!         if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!         endif
!       petroleum component mass
!         if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!         endif
!       salt mass - total salt mass - TMS_FRC(2,NTX)
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
!
!---      Energy equation residual  ---
!
          IF( ISLC(30).EQ.0 ) THEN
            RSDX = ABS(DPT)/TABS
            RSDAVG_FRC(IEQT) = RSDAVG_FRC(IEQT) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQT) ) THEN
              RSD_FRC(IEQT) = RSDX
              NSD_FRC(IEQT) = NTX
            ENDIF
          ENDIF
!
!---      Water mass equation residual  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOL_FRC(2,NTX)*SL_FRC(2,NTX)*
     &      XLW_FRC(2,NTX) + RHOG_FRC(2,NTX)*SG_FRC(2,NTX)*
     &      XGW_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          ACP = TMX*DTI*VOL_FRC(NTX)
          IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
            RSDX = MIN( 1.D-1*ABS(DPW)/PX,
     &        1.D-1*ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          ELSE
            RSDX = MIN( ABS(DPW),
     &        1.D-1*ABS(RSDL_FRC(IEQW,NTX)/(ACP+SMALL)) )
          ENDIF
          RSDAVG_FRC(IEQW) = RSDAVG_FRC(IEQW) + RSDX
          IF( RSDX.GT.RSD_FRC(IEQW) ) THEN
            RSD_FRC(IEQW) = RSDX
            NSD_FRC(IEQW) = NTX
          ENDIF
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          ACP = PORD_FRC(2,NTX)*(RHOG_FRC(2,NTX)*SG_FRC(2,NTX) +
     &      RHON_FRC(2,NTX)*SN_FRC(2,NTX))*DTI*VOL_FRC(NTX)
          ACP = TMX*DTI*VOL_FRC(NTX)
          DO IGC = 1,NGC+2
            IF( ISLC(45).NE.0 .AND. IGC.EQ.1 ) CYCLE
            IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
!
!---        Nonaqueous-liquid pressure  ---
!
            IF( IGC.EQ.IZMC_FRC(NTX) ) THEN
              IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                RSDX = MIN( ABS(DPC(IGC)),
     &            1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
              ELSE
                RSDX = MIN( 1.D-1*ABS(DPC(IGC))/PX,
     &            1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
              ENDIF
!
!---        Petroleum component mole fraction  ---
!
            ELSE
              RSDX = MIN( 1.D-2*ABS(DPC(IGC)),
     &          1.D-1*ABS(RSDL_FRC(IEQGC(IGC),NTX)/(ACP+SMALL)) )
            ENDIF
            RSDAVG_FRC(IEQGC(IGC)) = RSDAVG_FRC(IEQGC(IGC)) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQGC(IGC)) ) THEN
              RSD_FRC(IEQGC(IGC)) = RSDX
              NSD_FRC(IEQGC(IGC)) = NTX
            ENDIF
          ENDDO
!
!---      Salt mass equation residual  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS_FRC(2,NTX)*DTI*VOL_FRC(NTX)
            ACP = TMX*DTI*VOL_FRC(NTX)
            RSDX = MIN( ABS(DPS),
     &        1.D-1*ABS(RSDL_FRC(IEQS,NTX)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            RSDAVG_FRC(IEQS) = RSDAVG_FRC(IEQS) + RSDX
            IF( RSDX.GT.RSD_FRC(IEQS) ) THEN
              RSD_FRC(IEQS) = RSDX
              NSD_FRC(IEQS) = NTX
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      ENDDO
!
!---  Assign a convergence index  ---
!
      DO M = 1,ISVC
        RSDAVG_FRC(M) = 1.D+1*RSDAVG_FRC(M)/REAL(NT_FRC-NXP_FRC)
      ENDDO
      RSDX = 0.D+0
      DO M = 1,ISVC
        RSDX = MAX( RSD_FRC(M),RSDX )
        IF( RSD_FRC(M).GT.RSDMX ) ICNV = 2
      ENDDO
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) ICNV = 1
      IF( RSDX.GE.1.D+2 ) ICNV = 1
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.1 ) THEN
        IF( RSDX.GE.1.D+2 ) THEN
          WRITE(ISC,'(10X,A)') '---  Excessive Residual in Fault  ---'
          WRITE(IWR,'(10X,A)') '---  Excessive Residual in Fault  ---'
        ELSE
          WRITE(ISC,'(10X,A)') '---  Convergence Failure in Fault  ---'
          WRITE(IWR,'(10X,A)') '---  Convergence Failure in Fault  ---'
        ENDIF
!
!---    Energy equation  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          NTX = NSD_FRC(IEQT)
          IF( NTX.GT.0 ) THEN
            NPX = MOD(NPHAZ_FRC(2,NTX),10)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Energy Equation Maximum Residual = ',RSD_FRC(IEQT),
     &        ': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Energy Equation Maximum Residual = ',RSD_FRC(IEQT),
     &        ': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---    Water mass equation  ---
!
        NTX = NSD_FRC(IEQW)
        IF( NTX.GT.0 ) THEN
          NPX = MOD(NPHAZ_FRC(2,NTX),10)
          NCHX = INDEX( PH_CND(NPX),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Mass Equation Maximum Residual = ',RSD_FRC(IEQW),
     &      ': Fault Triangle = ',NTX,
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &      'Water Mass Equation Maximum Residual = ',RSD_FRC(IEQW),
     &      ': Fault Triangle = ',NTX,
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
        ENDIF
!
!---    CO2 mass equation  ---
!
        IF( ISLC(45).EQ.0 ) THEN
          NTX = NSD_FRC(IEQGC(1))
          IF( NTX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NTX),10)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'CO2 Mass Equation Maximum Residual = ',
     &        RSD_FRC(IEQGC(1)),': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'CO2 Mass Equation Maximum Residual = ',
     &        RSD_FRC(IEQGC(1)),': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---    CH4 mass equation  ---
!
        IF( ISLC(64).EQ.0 ) THEN
          NTX = NSD_FRC(IEQGC(2))
          IF( NTX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NTX),10)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'CH4 Mass Equation Maximum Residual = ',
     &        RSD_FRC(IEQGC(2)),': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'CH4 Mass Equation Maximum Residual = ',
     &        RSD_FRC(IEQGC(2)),': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---    Petroleum component mass equation  ---
!
        DO IGC = 3,NGC+2
          NTX = NSD_FRC(IEQGC(IGC))
          IF( NTX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NTX),10)
            NCAX = INDEX( PH_CND(NPX),'  ') - 1
            NCBX = INDEX( GCNM(IGC),'  ') - 1
            WRITE(ISC,'(4X,2A,1PE11.4,A,I6,2A)')
     &        GCNM(IGC)(1:NCBX),' Mass Equation Maximum Residual = ',
     &        RSD_FRC(IEQGC(IGC)),': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCAX)
            WRITE(IWR,'(4X,2A,1PE11.4,A,I6,2A)')
     &        GCNM(IGC)(1:NCBX),' Mass Equation Maximum Residual = ',
     &        RSD_FRC(IEQGC(IGC)),': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCAX)
          ENDIF
        ENDDO
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          NTX = NSD_FRC(IEQS)
          IF( NTX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NTX),10)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD_FRC(IEQS),
     &        ': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD_FRC(IEQS),
     &         ': Fault Triangle = ',NTX,
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---  Reduce time step  ---
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
!---      Loop over faults  ---
!
          DO NFX = 1,NF_FRC
!
!---      Loop over fault triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
            DO M = 1,6
              BETA_FRC(M,NTX) = 0.D+0
            ENDDO
            IBETA_FRC(NTX) = 0
            T_FRC(2,NTX) = T_FRC(1,NTX)
            PL_FRC(2,NTX) = PL_FRC(1,NTX)
            PG_FRC(2,NTX) = PG_FRC(1,NTX)
            PN_FRC(2,NTX) = PN_FRC(1,NTX)
            PSO_FRC(2,NTX) = PSO_FRC(1,NTX)
            PVA_FRC(2,NTX) = PVA_FRC(1,NTX)
            SG_FRC(2,NTX) = SG_FRC(1,NTX)
            SL_FRC(2,NTX) = SL_FRC(1,NTX)
            SN_FRC(2,NTX) = SN_FRC(1,NTX)
            YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
            TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
            NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
            DO IGC = 1,NGC+2
              TMC_FRC(IGC,2,NTX) = TMC_FRC(IGC,1,NTX)
              ZMC_FRC(IGC,2,NTX) = ZMC_FRC(IGC,1,NTX)
            ENDDO
          ENDDO
          ENDDO
!
!---      Loop over field nodes  ---
!
          DO N = 1,NFBN
            IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
            DO M = 1,6
              BETA(M,N) = 0.D+0
            ENDDO
            IBETA(N) = 0
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            PSO(2,N) = PSO(1,N)
            PVA(2,N) = PVA(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
            DO IGC = 1,NGC+2
              TMC(IGC,2,N) = TMC(IGC,1,N)
              ZMC(IGC,2,N) = ZMC(IGC,1,N)
            ENDDO
          ENDDO
!
!---      Coupled-well pressure  ---
!
          DO NCW = 1,N_CW
            P_CW(2,NCW) = P_CW(1,NCW)
          ENDDO
          NITER = 1
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
!---  End of RSDL_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORC_FRC_EOR
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
!     Compute fault source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PARM_FRC
      USE JACOB
      USE GLB_PAR
      USE GEOM_FRC
      USE FDVS_FRC
      USE FDVP_FRC
      USE FDVN_FRC
      USE FDVG_FRC
      USE FDVGC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 SRX(10+LSOLU+LNGC)
      REAL*8 ZMCX(LNGC)
      REAL*8 XGCX(LNGC),XNCX(LNGC),XMGCX(LNGC),XMNCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC_FRC_EOR'
!
!---  K-factors computed from fugacity coefficients  ---
!
      IFK = 0
#ifdef ecke
      DO IGC = 1,NGC+2
!
!---  Load component sources for linked aqueous species   ---
!                             or linked gaseous species   ---
!
        IF( ( ISPLK(14+NSPLK+IGC).NE.0 ) .OR.
     &      ( ISPLK(14+NSPLK+NGC+2+IGC).NE.0 ) ) THEN
!
!---      Load NGC sources associated with chemical
!         reactions, looping over faults  ---
!
          DO NFX = 1,NF_FRC
!
!---        Loop over fault triangles  ---
!
            DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---          Skip inactive triangles  ---
!
              IF( IXP_FRC(NTX).EQ.0 ) CYCLE
              SRCGCX = SRCGC_FRC(IGC,1,NTX)*DTI
              DO M = 2,ISVC+2
                SRCGC_FRC(IGC,M,NTX) = SRCGCX
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!
!---  Zero fault source terms, looping over faults  ---
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
          DO M = 2,ISVC+2
            SRCT_FRC(M,NTX) = 0.D+0
            SRCW_FRC(M,NTX) = 0.D+0
            SRCS_FRC(M,NTX) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
#else
!
!---  Zero fault source terms, looping over faults  ---
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
          DO M = 2,ISVC+2
            DO IGC = 1,NGC+2
              SRCGC_FRC(IGC,M,NTX) = 0.D+0
            ENDDO
            SRCT_FRC(M,NTX) = 0.D+0
            SRCW_FRC(M,NTX) = 0.D+0
            SRCS_FRC(M,NTX) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
#endif
!
!---  Loop over fault sources  ---
!
      DO NS = 1,NSR_FRC
        IF( TM.LE.SRC_FRC(1,1,NS) ) CYCLE
        SRX(1) = TM
        IF( ISRM_FRC(NS).EQ.1 ) THEN
          DO N = 1,10+NSOLU+NGC
            SRX(N) = SRC_FRC(N,1,NS)
          ENDDO
        ELSE
          DO M = 2,ISRM_FRC(NS)
            IF( TM.LE.SRC_FRC(1,M,NS) ) THEN
              DTSR = MIN( SRC_FRC(1,M,NS)-TM,DT )
              TFSR = (TM-0.5D+0*DTSR-SRC_FRC(1,M-1,NS))/
     &          (SRC_FRC(1,M,NS)-SRC_FRC(1,M-1,NS))
              DO N = 1,10+NSOLU+NGC
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
          DO M = 2,ISVC+2
            TX = T_FRC(M,NTX)
            PX = PSO_FRC(M,NTX) + PATM
!
!---        Power  ---
!
            IF( ISRT_FRC(NS).EQ.1 ) THEN
              SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) + SRX(2)
!
!---        Power Density  ---
!
            ELSEIF( ISRT_FRC(NS).EQ.2 ) THEN
              SRCT_FRC(M,NTX) = SRCT_FRC(M,NTX) + SRX(2)*VOL_FRC(NTX)
!
!---        Fluid mass rate  ---
!
            ELSEIF( ISRT_FRC(NS)/100.GE.1 .AND. 
     &        ISRT_FRC(NS)/100.LE.3 ) THEN
!
!---          Fluid withdrawl from grid cell  ---
!
              IF( SRX(2).LE.0.D+0 ) THEN
                FMGX = RKG_FRC(2,NTX)/VISG_FRC(2,NTX)
                FMLX = RKL_FRC(2,NTX)/VISL_FRC(2,NTX)
                FMNX = RKN_FRC(2,NTX)/VISN_FRC(2,NTX)
                SFMX = FMGX + FMLX + FMNX 
!
!---            Petroleum components from gas and 
!               nonaqueous liquid  ---
!
                DO IGC = 1,NGC+2
                  SRCGC_FRC(IGC,M,NTX) = SRCGC_FRC(IGC,M,NTX) + SRX(2)*
     &              (XGC(IGC,M,NTX)*FMGX + XNC(IGC,M,NTX)*FMNX)/SFMX
                ENDDO
!
!---            CO2 from aqueous  ---
!
                SRCGC_FRC(1,M,NTX) = SRCGC_FRC(1,M,NTX) + SRX(2)*
     &            XLA_FRC(M,NTX)*FMLX/SFMX
!
!---            Water from aqueous and gas  ---
!
                SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) + SRX(2)*
     &            (XGW_FRC(M,NTX)*FMGX + XLW_FRC(M,NTX)*FMLX)/SFMX
!
!---            Salt from aqueous  ---
!
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + SRX(2)*
     &            XLS_FRC(M,NTX)*FMLX/SFMX
!
!---          Fluid injection into grid cell  ---
!
              ELSE
                IF( ISLC(30).EQ.0 ) TX = SRX(3)
!
!---            Source state #1  ---
!
!               SL = 1.0
!               SG+SN = 0.0
!
!               Declared variables:
!
!               SRX(4) - aqueous CO2 relative saturation, or
!                        aqueous CO2 mass fraction
!               SRX(5) - aqueous salt relative saturation, or
!                        aqueous salt mass fraction
!
                IF( ISRT_FRC(NS)/100.EQ.1 ) THEN
                  PVAX = SRX(4)
                  TMSX = SRX(5)
                  ISR4X = MOD(ISRT_FRC(NS),10)
                  ISR5X = MOD(ISRT_FRC(NS)/10,10)
                  CALL FLH_IC1( PX,PVAX,RHOLX,RHOMLX,TMSX,TX,XLAX,XLSX,
     &              XLWX,XMLAX,XMLSX,XMLWX,YLSX,ZMCX,ISR4X,ISR5X )
!
!---              CO2 from aqueous  ---
!
                  SRCGC_FRC(1,M,NTX) = SRCGC_FRC(1,M,NTX) + SRX(2)*XLAX
!
!---              Salt from aqueous  ---
!
                  SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + SRX(2)*XLSX
!
!---              Water from aqueous  ---
!
                  SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) + SRX(2)*XLWX
!
!---            Source state #2  ---
!
!               SL < 1.0
!               SG+SN > 0.0
!
!               Declared variables:
!
!               SRX(4) - aqueous saturation
!               SRX(5) - aqueous salt relative saturation, or
!                         aqueous salt mass fraction
!               SRX(6+) - total mole fraction of petroleum components
!
                ELSEIF( ISRT_FRC(NS)/100.EQ.2 ) THEN
                  SLX = SRX(4)
                  TMSX = SRX(5)
                  ISR5X = MOD(ISRT_FRC(NS)/10,10)
                  DO IGC = 1,NGC+2
                    JGC = 5 + IGC
                    ZMCX(IGC) = SRX(JGC)
                  ENDDO
                  IZNX = 0
                  CALL FLH_IC2( PX,RHOGX,RHOMGX,RHOLX,RHOMLX,RHONX,
     &              RHOMNX,SGX,SLX,SNX,TMSX,TX,XGCX,XMGCX,XGWX,XMGWX,
     &              XLAX,XLSX,XLWX,XMLAX,XMLSX,XMLWX,XNCX,XMNCX,YLSX,
     &              ZGX,ZMCX,ZNX,ISR5X,IZNX,NPHAZX )
!
!---              Petroleum components from gas and 
!                 nonaqueous liquid  ---
!
                  DO IGC = 1,NGC+2
                    SRCGC_FRC(IGC,M,NTX) = SRCGC_FRC(IGC,M,NTX) + 
     &                SRX(2)*(SGX*XGCX(IGC) + SNX*XNCX(IGC))
                  ENDDO
!
!---              CO2 from aqueous  ---
!
                  SRCGC_FRC(1,M,NTX) = SRCGC_FRC(1,M,NTX) + 
     &              SRX(2)*SLX*XLAX
!
!---              Salt from aqueous  ---
!
                  SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + SRX(2)*SLX*XLSX
!
!---              Water from aqueous and gas  ---
!
                  SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) + 
     &              SRX(2)*(SGX*XGWX + SLX*XLWX)
!
!---            Source state #3  ---
!
!               SL = 0.0
!               SG+SN = 1.0
!
!               Declared variables:
!
!               SRX(5) - water vapor relative saturation
!               SRX(6+) - total mole fraction of petroleum comp.
!          
                ELSEIF( ISRT_FRC(NS)/100.EQ.3 ) THEN
                  DO IGC = 1,NGC+2
                    JGC = 6 + IGC
                    ZMCX(IGC) = SRX(JGC)
                  ENDDO
                  PVWX = SRX(5)
                  YLSX = 0.D+0
                  IZNX = 0
                  CALL FLH_IC3( PVWX,PX,RHOGX,RHOMGX,RHONX,RHOMNX,SGX,
     &              SNX,TX,XGCX,XGWX,XMGCX,XMGWX,XNCX,XMNCX,
     &              YLSX,ZGX,ZMCX,ZNX,IZNX,NPHAZX )
!
!---              Petroleum components from gas and 
!                 nonaqueous liquid  ---
!
                  DO IGC = 1,NGC+2
                    SRCGC_FRC(IGC,M,NTX) = SRCGC_FRC(IGC,M,NTX) + 
     &                SRX(2)*(SGX*XGCX(IGC) + SNX*XNCX(IGC))
                  ENDDO
                ENDIF
              ENDIF
!
!---        Fluid volumetric rate  ---
!
            ELSEIF( ISRT_FRC(NS)/100.GE.4 .AND. 
     &        ISRT_FRC(NS)/100.LE.6 ) THEN
!
!---          Fluid withdrawl from grid cell  ---
!
              IF( SRX(2).LE.0.D+0 ) THEN
                FMGX = RKG_FRC(2,NTX)/VISG_FRC(2,NTX)
                FMLX = RKL_FRC(2,NTX)/VISL_FRC(2,NTX)
                FMNX = RKN_FRC(2,NTX)/VISN_FRC(2,NTX)
                SFMX = FMGX + FMLX + FMNX 
!
!---            Petroleum components from gas and 
!               nonaqueous liquid  ---
!
                DO IGC = 1,NGC+2
                  SRCGC_FRC(IGC,M,NTX) = SRCGC_FRC(IGC,M,NTX) + SRX(2)*
     &              (XGC(IGC,M,NTX)*RHOG_FRC(M,NTX)*FMGX + 
     &              XNC(IGC,M,NTX)*RHON_FRC(M,NTX)*FMNX)/SFMX
                ENDDO
!
!---            CO2 from aqueous  ---
!
                SRCGC_FRC(1,M,NTX) = SRCGC_FRC(1,M,NTX) + SRX(2)*
     &            XLA_FRC(M,NTX)*RHOL_FRC(M,NTX)*FMLX/SFMX
!
!---            Water from aqueous and gas  ---
!
                SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) 
     &            + SRX(2)*(XGW_FRC(M,NTX)*RHOG_FRC(M,NTX)*FMGX
     &            + XLW_FRC(M,NTX)*RHOL_FRC(M,NTX)*FMLX)/SFMX
!
!---            Salt from aqueous  ---
!
                SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + SRX(2)*
     &            XLS_FRC(M,NTX)*RHOL_FRC(M,NTX)*FMLX/SFMX
!
!---          Fluid injection into grid cell  ---
!
              ELSE
                IF( ISLC(30).EQ.0 ) TX = SRX(3)
!
!---            Source state #1  ---
!
!               SL = 1.0
!               SG+SN = 0.0
!
!               Declared variables:
!
!               SRX(4) - aqueous CO2 relative saturation, or
!                        aqueous CO2 mass fraction
!               SRX(5) - aqueous salt relative saturation, or
!                        aqueous salt mass fraction
!
                IF( ISRT_FRC(NS)/100.EQ.4 ) THEN
                  PVAX = SRX(4)
                  TMSX = SRX(5)
                  ISR4X = MOD(ISRT_FRC(NS),10)
                  ISR5X = MOD(ISRT_FRC(NS)/10,10)
                  CALL FLH_IC1( PX,PVAX,RHOLX,RHOMLX,TMSX,TX,XLAX,XLSX,
     &              XLWX,XMLAX,XMLSX,XMLWX,YLSX,ZMCX,ISR4X,ISR5X )
!
!---              CO2 from aqueous  ---
!
                  SRCGC_FRC(1,M,NTX) = SRCGC_FRC(1,M,NTX) + 
     &              SRX(2)*XLAX*RHOLX
!
!---              Salt from aqueous  ---
!
                  SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + SRX(2)*XLSX*RHOLX
!
!---              Water from aqueous  ---
!
                  SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) + SRX(2)*XLWX*RHOLX
!
!---            Source state #2  ---
!
!               SL < 1.0
!               SG+SN > 0.0
!
!               Declared variables:
!
!               SRX(4) - aqueous saturation
!               SRX(5) - aqueous salt relative saturation, or
!                         aqueous salt mass fraction
!               SRX(6+) - total mole fraction of petroleum components
!
                ELSEIF( ISRT_FRC(NS)/100.EQ.5 ) THEN
                  SLX = SRX(4)
                  TMSX = SRX(5)
                  ISR5X = MOD(ISRT_FRC(NS)/10,10)
                  DO IGC = 1,NGC+2
                    JGC = 5 + IGC
                    ZMCX(IGC) = SRX(JGC)
                  ENDDO
                  IZNX = 0
                  CALL FLH_IC2( PX,RHOGX,RHOMGX,RHOLX,RHOMLX,RHONX,
     &              RHOMNX,SGX,SLX,SNX,TMSX,TX,XGCX,XMGCX,XGWX,XMGWX,
     &              XLAX,XLSX,XLWX,XMLAX,XMLSX,XMLWX,XNCX,XMNCX,YLSX,
     &              ZGX,ZMCX,ZNX,ISR5X,IZNX,NPHAZX )
!
!---              Petroleum components from gas and 
!                 nonaqueous liquid  ---
!
                  DO IGC = 1,NGC+2
                    SRCGC_FRC(IGC,M,NTX) = SRCGC_FRC(IGC,M,NTX) + 
     &                SRX(2)*(SGX*XGCX(IGC)*RHOGX + SNX*XNCX(IGC)*RHONX)
                  ENDDO
!
!---              CO2 from aqueous  ---
!
                  SRCGC_FRC(1,M,NTX) = SRCGC_FRC(1,M,NTX) + 
     &              SRX(2)*SLX*XLAX*RHOLX
!
!---              Salt from aqueous  ---
!
                  SRCS_FRC(M,NTX) = SRCS_FRC(M,NTX) + 
     &              SRX(2)*SLX*XLSX*RHOLX
!
!---              Water from aqueous and gas  ---
!
                  SRCW_FRC(M,NTX) = SRCW_FRC(M,NTX) + 
     &              SRX(2)*(SGX*XGWX*RHOGX + SLX*XLWX*RHOLX)
!
!---            Source state #3  ---
!
!               SL = 0.0
!               SG+SN = 1.0
!
!               Declared variables:
!
!               SRX(5) - water vapor relative saturation
!               SRX(6+) - total mole fraction of petroleum comp.
!          
                ELSEIF( ISRT_FRC(NS)/100.EQ.6 ) THEN
                  DO IGC = 1,NGC+2
                    JGC = 6 + IGC
                    ZMCX(IGC) = SRX(JGC)
                  ENDDO
                  PVWX = SRX(5)
                  YLSX = 0.D+0
                  IZNX = 0
                  CALL FLH_IC3( PVWX,PX,RHOGX,RHOMGX,RHONX,RHOMNX,SGX,
     &              SNX,TX,XGCX,XGWX,XMGCX,XMGWX,XNCX,XMNCX,
     &              YLSX,ZGX,ZMCX,ZNX,IZNX,NPHAZX )
!
!---              Petroleum components from gas and 
!                 nonaqueous liquid  ---
!
                  DO IGC = 1,NGC+2
                    SRCGC_FRC(IGC,M,NTX) = SRCGC_FRC(IGC,M,NTX) + 
     &                SRX(2)*(SGX*XGCX(IGC)*RHOGX + SNX*XNCX(IGC)*RHONX)
                  ENDDO
                ENDIF
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
!---  End of SORC_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THD_FRC_EOR( IGC )
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
!     Compute the contribution to the energy flux by thermal conduction
!     for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVT_FRC
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
      SUB_LOG(ISUB_LOG) = '/THD_FRC_EOR'
!
!---  Loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
        DO NT1X = IP_FRC(1,NFX),IP_FRC(2,NFX)
          IF( IXP_FRC(NT1X).EQ.0 ) CYCLE
          IZN1 = IZ_FRC(NT1X)
!
!---      Loop over fault triangle to triangle connections ---
!
          DO NCX = IPF_FRC(1,NT1X),IPF_FRC(2,NT1X)
            NT2X = ITCM_FRC(NCX)
            IF( IXP_FRC(NT2X).EQ.0 ) CYCLE
            IZN2 = IZ_FRC(NT2X)
            DFF1X = DFFM_FRC(NCX)
            DFF2X = (DFF_FRC(NCX)-DFFM_FRC(NCX))
            DTK = T_FRC(2,NT1X)-T_FRC(2,NT2X)
!
!---        Loop over flux increments ---
!
            DO M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
!
!---          Parallel function  ---
!
              IF( ITHK(IZN1).EQ.2 ) THEN
                TK1 = MAX(1.D+0-PORD_FRC(MP,NT1X),0.D+0)*THKS(1,IZN1) + 
     &            PORD_FRC(MP,NT1X)*(THKL_FRC(MP,NT1X)*SL_FRC(MP,NT1X) +
     &            THKN_FRC(MP,NT1X)*SN_FRC(MP,NT1X) + 
     &            THKG_FRC(MP,NT1X)*SG_FRC(MP,NT1X))
!
!---          Somerton function  ---
!
              ELSEIF( ITHK(IZN1).EQ.4 ) THEN
                TK1 = THKS(1,IZN1) + 
     &            SQRT(SL_FRC(MP,NT1X))*(THKS(4,IZN1)-THKS(1,IZN1))
              ENDIF
!
!---          Parallel function  ---
!
              IF( ITHK(IZN2).EQ.2 ) THEN
                TK2 = MAX(1.D+0-PORD_FRC(MN,NT2X),0.D+0)*THKS(1,IZN2) +
     &            PORD_FRC(MN,NT2X)*(THKL_FRC(MN,NT2X)*SL_FRC(MN,NT2X) +
     &            THKN_FRC(MN,NT2X)*SN_FRC(MN,NT2X) + 
     &            THKG_FRC(MN,NT2X)*SG_FRC(MN,NT2X))
!
!---          Somerton function  ---
!
              ELSEIF( ITHK(IZN2).EQ.4 ) THEN
                TK2 = THKS(1,IZN2) + 
     &            SQRT(SL_FRC(MN,NT2X))*(THKS(4,IZN2)-THKS(1,IZN2))
              ENDIF
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
!---  End of THD_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THG_FRC_EOR( IGC )
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
!     Compute the contribution to the energy flux by gas advection
!     for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVT_FRC
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
      SUB_LOG(ISUB_LOG) = '/THG_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
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
              H1 = HG_FRC(MP,NT1X)*RHOG_FRC(MP,NT1X)
              H2 = HG_FRC(MN,NT2X)*RHOG_FRC(MN,NT2X)
!
!---          Upwind mean  ---
!
              INDX = -4
              HM = DIFMN( H1,H2,DFF1X,DFF2X,UFFG(1,NCX),INDX )
              UFFQ(M,NCX) = UFFQ(M,NCX) + HM*UFFG(M,NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THG_FRC_EOR group ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THL_FRC_EOR( IGC )
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
!     Compute the contribution to the energy flux by aqueous advection
!     for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVT_FRC
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
      SUB_LOG(ISUB_LOG) = '/THL_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
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
              H1 = HL_FRC(MP,NT1X)*RHOL_FRC(MP,NT1X)
              H2 = HL_FRC(MN,NT2X)*RHOL_FRC(MN,NT2X)
!
!---          Upwind mean  ---
!
              INDX = -4
              HM = DIFMN( H1,H2,DFF1X,DFF2X,UFFL(1,NCX),INDX )
              UFFQ(M,NCX) = UFFQ(M,NCX) + HM*UFFL(M,NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THL_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THN_FRC_EOR( IGC )
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
!     Compute the contribution to the energy flux by aqueous advection
!     for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVT_FRC
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
      SUB_LOG(ISUB_LOG) = '/THN_FRC_EOR'
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
!---      Loop over fault triangle to triangle connections ---
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
              H1 = HN_FRC(MP,NT1X)*RHON_FRC(MP,NT1X)
              H2 = HN_FRC(MN,NT2X)*RHON_FRC(MN,NT2X)
!
!---          Upwind mean  ---
!
              INDX = -4
              HM = DIFMN( H1,H2,DFF1X,DFF2X,UFFN(1,NCX),INDX )
              UFFQ(M,NCX) = UFFQ(M,NCX) + HM*UFFN(M,NCX)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THN_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRNS_FRC_EOR
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
!     Fault to matrix transfer functions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 20 May 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FLUX_FRC
      USE FLUXC_FRC
      USE FDVT_FRC
      USE FDVT
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP
      USE FDVN_FRC
      USE FDVN
      USE FDVG_FRC
      USE FDVGC_FRC
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 UFFDGCX(LNGC),UFFDNCX(LNGC)
      REAL*8 KM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TRNS_FRC_EOR'
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
!
!---    Loop over fault triangle to grid cell connections  ---
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
!---        Fault permeability, m^2  ---
!
            PERM_NTX = PERMRF_FRC(MN,NTX)*PERM(1,IZ_FRC(NTX))
!
!---        Permeability  ---
!
            INDX = -1
            KM = DIFMN(PERM_NTX,PERMX,DFN_FRC(NCX),
     &        DFN_FRC(NCX),ZERO,INDX )
            IF( PERMX/EPSL.LT.EPSL ) KM = 0.D+0
            IF( PERM_NTX/EPSL.LT.EPSL ) KM = 0.D+0
!
!---        Aqueous head difference for fault and grid cell, m  ---
!
            HDLX = PL_FRC(MN,NTX) - PL(MP,N) + 
     &        5.D-1*GRAV*(ZP_FRC(NTX)-ZP(N))*
     &        (RHOL_FRC(MN,NTX)+RHOL(MP,N))
!
!---        Matrix aqueous relative permeability  ---
!
            RKLX = RKL(1,MP,N)
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
!---        Aqueous volumetric rate fault to matrix, m^3/s  ---
!
            QLX = AFN_FRC(NCX)*KM*RKLM*HDLX/DFN_FRC(NCX)/VLM
!
!---        Gas head difference for fault and grid cell, m  ---
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
!---        Gas volumetric rate fault to matrix, m^3/s  ---
!
            QGX = AFN_FRC(NCX)*KM*RKGM*HDGX/DFN_FRC(NCX)/VGM         
!
!---        NAPL head difference for fault and grid cell, m  ---
!
            HDNX = PN_FRC(MN,NTX) - PN(MP,N) + 
     &        5.D-1*GRAV*(ZP_FRC(NTX)-ZP(N))*
     &        (RHON_FRC(MN,NTX)+RHON(MP,N))
!
!---        NAPL relative permeability  ---
!
            INDX = 10
            RKNM = DIFMN( RKN_FRC(MN,NTX),RKN(MP,N),DFN_FRC(NCX),
     &        DFN_FRC(NCX),HDNX,INDX )
!
!---        NAPL viscosity Pa s  ---
!
            INDX = 7
            VNM = DIFMN( VISN_FRC(MN,NTX),VISN(MP,N),DFN_FRC(NCX),
     &        DFN_FRC(NCX),HDNX,INDX )
!
!---        NAPL volumetric rate fault to matrix, m^3/s  ---
!
            QNX = AFN_FRC(NCX)*KM*RKNM*HDNX/DFN_FRC(NCX)/VNM         
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
!---          Overall heat transfer, W  ---
!
              QTCX = AFN_FRC(NCX)*HCOX*(T_FRC(MN,NTX)-T(MP,N))
            ENDIF
!
!---        Store aqueous, gas, and heat flow from fault
!           to grid cell  ---
!
            IF( M.EQ.1 ) THEN
              UFML(NCX) = QLX
              UFMG(NCX) = QGX
              UFMN(NCX) = QNX
              UFMT(NCX) = QTCX
            ENDIF
!
!---        Petroleum component molar gas diffusive flux rate from 
!           fault to grid cell, kmol/s  ---
!
            DO IGC = 1,NGC+2
              UFFDGCX(IGC) = 0.D+0
!
!---          Iso-CO2 option  ---
!
              IF( IGC.EQ.1 .AND. ISLC(45).NE.0 ) CYCLE
!
!---          Iso-CH4 option  ---
!
              IF( IGC.EQ.2 .AND. ISLC(64).NE.0 ) CYCLE
              DXMGC = XMGC_FRC(IGC,MN,NTX) - XMGC(IGC,MP,N)
              DF_NTX = TORG_FRC(MN,NTX)*PORD_FRC(MN,NTX)*SG_FRC(MN,NTX)*
     &          DFGC_FRC(IGC,MN,NTX)*RHOMG_FRC(MN,NTX)
              DFX = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &          DFGC(IGC,MP,N)*RHOMG(MP,N)
              INDX = 12
              DFM = DIFMN( DF_NTX,DFX,DFN_FRC(NCX),DFN_FRC(NCX),
     &          DXMGC,INDX )
              UFFDGCX(IGC) = DFM*(XMGC_FRC(IGC,MN,NTX) 
     &          - XMGC(IGC,MP,N))/DFN_FRC(NCX)
            ENDDO
!
!---        CO2 molar aqueous diffusive flux rate from 
!           fault to grid cell, kmol/s, inactive for iso-CO2
!           option  ---
!
            UFFDLAX = 0.D+0
            IF( ISLC(45).EQ.0 .AND. ISLC(68).NE.0 ) THEN
              DXMLA = XMLA_FRC(2,NTX) - XMLA(2,N)
              DFX = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)*
     &          RHOML(MP,N)
              DF_NTX = TORL_FRC(MN,NTX)*PORD_FRC(MN,NTX)*SL_FRC(MN,NTX)*
     &          DFLA_FRC(MN,NTX)*RHOML_FRC(MN,NTX)
              INDX = 14
              DFM = DIFMN( DF_NTX,DFX,DFN_FRC(NCX),DFN_FRC(NCX),
     &          DXMLA,INDX )
              UFFDLAX = DFM*(XMLA_FRC(MN,NTX) - XMLA(MP,N))/DFN_FRC(NCX)
            ENDIF
!
!---        Salt diffusive aqueous flux rate from fault to 
!           grid cell, kg/s, inactive for iso-brine option  ---
!
            UFFDSX = 0.D+0
            IF( ISLC(32).EQ.0 ) THEN
!
!---          Diffusion coefficients  ---
!
              IF( IEDLS.EQ.1 ) THEN
                TCOR = (T(MP,N)+TABS)/TSPRF
                SMDL1 = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
                DFCX = SL(MP,N)*SMDL1
                TCOR = (T_FRC(MN,NTX)+TABS)/TSPRF
                SMDL2 = DFLS_FRC(MN,NTX)*TCOR*(VISRL/VISL_FRC(MN,NTX))
                DFC_NTX = SL_FRC(MN,NTX)*SMDL2
              ELSEIF( IEDLS.EQ.3 ) THEN
                DFCX = SL(MP,N)*DFLS(MP,N)
                DFC_NTX = SL_FRC(MN,NTX)*DFLS_FRC(MN,NTX)
              ENDIF
              INDX = 18
              DFC = DIFMN(DFC_NTX,DFCX,DFN_FRC(NCX),DFN_FRC(NCX),
     &          UFML(NCX),INDX)
              DDL = DFC/DFN_FRC(NCX)
              UFFDSX = DDL*(XLS_FRC(MN,NTX)*RHOL_FRC(MN,NTX)
     &          - XLS(MP,N)*RHOL(MP,N))
            ENDIF
!
!---        Water diffusive aqueous molar flux rate from 
!           fault to grid cell, kmol/s  ---
!
            UFFDLWX = -(UFFDLAX + UFFDSX/WTMS)
!
!---        Water diffusive gas molar flux rate from 
!           fault to grid cell, kmol/s  ---
!
            DXMGW = XMGW_FRC(2,NTX) - XMGW(2,N)
            DFX = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DF_NTX = TORG_FRC(MN,NTX)*PORD_FRC(MN,NTX)*SG_FRC(MN,NTX)*
     &        DFGW_FRC(MN,NTX)*RHOMG_FRC(MN,NTX)
            INDX = 12
            DFM = DIFMN( DF_NTX,DFX,DFN_FRC(NCX),DFN_FRC(NCX),
     &        DXMGW,INDX )
            UFFDGWX = DFM*(XMGW_FRC(MN,NTX)-XMGW(MP,N))
     &        /DFN_FRC(NCX)
!
!---        Petroleum component molar NAPL diffusive flux rate from 
!           fault to grid cell, kmol/s  ---
!
            DO IGC = 1,NGC+2
              UFFDNCX(IGC) = 0.D+0
!
!---          Iso-CO2 option  ---
!
              IF( IGC.EQ.1 .AND. ISLC(45).NE.0 ) CYCLE
!
!---          Iso-CH4 option  ---
!
              IF( IGC.EQ.2 .AND. ISLC(64).NE.0 ) CYCLE
              DXMNC = XMNC_FRC(IGC,MN,NTX) - XMNC(IGC,MP,N)
              DF_NTX = TORN_FRC(MN,NTX)*PORD_FRC(MN,NTX)*SN_FRC(MN,NTX)*
     &          DFNC_FRC(IGC,MN,NTX)*RHOMN_FRC(MN,NTX)
              DFX = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &          DFNC(IGC,MP,N)*RHOMN(MP,N)
              INDX = 20
              DFM = DIFMN( DF_NTX,DFX,DFN_FRC(NCX),DFN_FRC(NCX),
     &          DXMNC,INDX )
              UFFDNCX(IGC) = DFM*(XMNC_FRC(IGC,MN,NTX) 
     &          - XMNC(IGC,MP,N))/DFN_FRC(NCX)
            ENDDO
!
!---        Water mass flow rate from fault to grid cell, kg/s  ---
!
            TRNSW_FRC(M,NCX) = 
     &        XLW_FRC(MN,NTX)*RHOL_FRC(MN,NTX)*MAX(QLX,ZERO) -
     &        XLW(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO) + UFFDLWX*WTMW + 
     &        XGW_FRC(MN,NTX)*RHOG_FRC(MN,NTX)*MAX(QGX,ZERO) -
     &        XGW(MP,N)*RHOG(MP,N)*MAX(-QGX,ZERO) + UFFDGWX*WTMW
!
!---        Petroleum component mass flow rate from fault to 
!           grid cell, kg/s  ---
!
            DO IGC = 1,NGC+2
              TRNSC_FRC(IGC,M,NCX) = 
     &          XGC_FRC(IGC,MN,NTX)*RHOG_FRC(MN,NTX)*MAX(QGX,ZERO) -
     &          XGC(IGC,MP,N)*RHOG(MP,N)*MAX(-QGX,ZERO) +
     &          UFFDGCX(IGC)*GCPP(1,IGC) + 
     &          XNC_FRC(IGC,MN,NTX)*RHON_FRC(MN,NTX)*MAX(QNX,ZERO) -
     &          XNC(IGC,MP,N)*RHON(MP,N)*MAX(-QNX,ZERO) +
     &          UFFDNCX(IGC)*GCPP(1,IGC)
              IF( IGC.EQ.1 ) THEN
                TRNSC_FRC(IGC,M,NCX) = TRNSC_FRC(IGC,M,NCX) + 
     &            XLA_FRC(MN,NTX)*RHOL_FRC(MN,NTX)*MAX(QLX,ZERO) -
     &            XLA(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO) + UFFDLAX*WTMA
              ENDIF
            ENDDO
!
!---        Salt mass flow rate from fault to grid cell, kg/s
!           inactive for iso-brine option  ---
!
            IF( ISLC(32).EQ.0 ) THEN
              TRNSS_FRC(M,NCX) = 
     &          XLS_FRC(MN,NTX)*RHOL_FRC(MN,NTX)*MAX(QLX,ZERO) -
     &          XLS(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO) + UFFDSX
            ENDIF
!
!---        Total heat transfer rate from fault to grid cell, W  ---
!
            IF( ISLC(30).EQ.0 ) THEN
              TRNSQ_FRC(M,NCX) = QTCX + 
     &          HL_FRC(MN,NTX)*RHOL_FRC(MN,NTX)*MAX(QLX,ZERO) -
     &          HL(MP,N)*RHOL(MP,N)*MAX(-QLX,ZERO) +
     &          HN_FRC(MN,NTX)*RHON_FRC(MN,NTX)*MAX(QNX,ZERO) -
     &          HN(MP,N)*RHON(MP,N)*MAX(-QNX,ZERO) +
     &          HG_FRC(MN,NTX)*RHOG_FRC(MN,NTX)*MAX(QGX,ZERO) -
     &          HG(MP,N)*RHOG(MP,N)*MAX(-QGX,ZERO)
            ENDIF
!
!---      Loop over node increment indices  ---
!
          ENDDO
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
!---  End of TRNS_FRC_EOR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_FRC_BO
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
!     Update the primary variables for the black-oil option for
!     faults .
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE OUTPU
      USE JACOB
      USE HYST
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FILES
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP
      USE FDVGC_FRC
      USE FDVGC
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
      CHARACTER*128 PH_CND(5)
      REAL*8 DPC(LNGC)
      INTEGER MPC(LNGC)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Phase Condition #1, SL = 1.0, SG+SN = 0.0', 
     &  'Phase Condition #2, SL < 1.0, ST > 0.0, SG = 0.0', 
     &  'Phase Condition #3, SL < 1.0, ST > 0.0, SG > 0.0', 
     &  'Phase Condition #4, SL = 0.0, ST = 1.0, SG = 0.0', 
     &  'Phase Condition #5, SL = 0.0, ST = 1.0, SG > 0.0'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_FRC_BO'
!
!---  Black-oil option  ---
!
      IF( ICNV.EQ.1 ) GOTO 300
      IERR = 0
!
!---  Update primary variables, loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        IF( IXP_FRC(NTX).EQ.0 ) CYCLE
        IZN = IZ_FRC(NTX)
!
!---    Water mass equation  ---
!
        MPW = IM_FRC(IEQW,NMD)
        DPW = BLU(MPW)
!
!---    Gas mass equation  ---
!
        MPC(1) = IM_FRC(IEQGC(1),NMD)
        DPC(1) = BLU(MPC(1))
!
!---    Oil mass equation  ---
!
        MPC(2) = IM_FRC(IEQGC(2),NMD)
        DPC(2) = BLU(MPC(2))
!
!---    Salt mass equation  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = IM_FRC(IEQS,NMD)
          DPS = BLU(MPS)
        ELSE
          DPS = 0.D+0
        ENDIF
!
!---    Phase condition #1  ---
!
!       SL_FRC = 1.0
!       SG_FRC+SN_FRC = 0.0
!       Aqueuous saturated
!
!       water mass - aqueous (system) pressure - PL_FRC,PSO_FRC
!       gas mass - total gas mass - TMC_FRC(1,2,NTX)
!       oil mass - total oil mass - TMC_FRC(2,2,NTX)
!       salt mass - total salt brine mass fraction - YLS_FRC
!
        IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
!
!---      Limit changes in total pressure,
!         excessive primary variable change limit at 100 MPa  ---
!
          IF( ABS(DPW).GT.1.D+8 ) IERR = 1
          DPX = 2.5D-2*(PSO_FRC(2,NTX)+PATM)
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPW
          PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
!
!---      Zero negative corrections for zero gas mass  ---
!
          IF( TMC_FRC(1,2,NTX).LT.EPSL .AND. DPC(1).LT.EPSL ) THEN
            BLU(MPC(1)) = 0.D+0
            DPC(1) = 0.D+0
          ENDIF
          TMC_FRC(1,2,NTX) = MAX( (TMC_FRC(1,2,NTX)+DPC(1)),0.D+0 )
!
!---      Zero negative corrections for zero oil mass  ---
!
          IF( TMC_FRC(2,2,NTX).LT.EPSL .AND. DPC(2).LT.EPSL ) THEN
            BLU(MPC(2)) = 0.D+0
            DPC(2) = 0.D+0
          ENDIF
          TMC_FRC(2,2,NTX) = MAX( (TMC_FRC(2,2,NTX)+DPC(2)),0.D+0 )
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            IF( YLS_FRC(2,NTX).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_FRC(2,NTX)/EPSL.LT.EPSL .AND. 
     &        DPS/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_FRC(2,NTX) = YLS_FRC(2,NTX)+DPS
            IF( YLS_FRC(2,NTX).LT.1.D-12 ) YLS_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---    Phase condition #2  ---
!
!       SG_FRC = 0.0
!       SL_FRC > 0.0
!       SN_FRC > 0.0
!       Aqueous and nonaqueous liquid saturated
!
!       water mass - aqueous saturation - SL_FRC
!       gas mass - gas to oil ratio - POSM_FRC
!       oil mass - nonaqueous (system) pressure - PN_FRC,PSO_FRC
!       salt mass - total salt brine mass fraction - YLS_FRC
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPW
          IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
!
!---      Limit changes in gas to oil ratio, and zero negative 
!         corrections for zero gas to oil ratios  ---
!
          DPX = 5.D-2/REAL(NITER)
          DPC(1) = SIGN( MIN(ABS(DPX),ABS(DPC(1))),DPC(1) )
          POSM_FRC(2,NTX) = POSM_FRC(2,NTX) + DPC(1)
!
!---      Limit changes in nonaqueous-liquid pressure (system pressure),
!         excessive primary variable change limit at 100 MPa  ---
!
          IF( ABS(DPC(2)).GT.1.D+8 ) IERR = 1
          DPX = MAX( 1.D+6,1.D-1*(PN_FRC(2,NTX)-PL_FRC(2,NTX)) )
          DPX = 2.5D-2*(PSO_FRC(2,NTX)+PATM)
          DPC(2) = SIGN( MIN(ABS(DPX),ABS(DPC(2))),DPC(2) )
!
!---      Bubble-point pressure as a function of gas to oil ratio
!         (m^3/m^3) and temperature (C)  ---
!
!          CALL BO_BP( PBPX,POSM_FRC(2,NTX),T_FRC(2,NTX) )
!
!---      Pressure below bubble-point pressure  ---
!
!          IF( PSO_FRC(2,NTX)+DPC(2).LT.PBPX-PATM ) DPC(2) = 6.D-1*DPC(2)
          PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPC(2)
          PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            IF( YLS_FRC(2,NTX).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_FRC(2,NTX)/EPSL.LT.EPSL .AND. 
     &        DPS/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_FRC(2,NTX) = YLS_FRC(2,NTX)+DPS
            IF( YLS_FRC(2,NTX).LT.1.D-12 ) YLS_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---    Phase condition #3  ---
!
!       SG_FRC > 0.0
!       SL_FRC > 0.0
!       SN_FRC > 0.0
!       Three-phase
!
!       water mass - aqueous saturation - SL_FRC
!       gas mass - gas saturation - SG_FRC
!       oil mass - gas (system) pressure - PG_FRC,PSO_FRC
!       salt mass - total salt brine mass fraction - YLS_FRC
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPW
          IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
!
!---      Limit changes in gas saturation  ---
!
          DPX = 2.5D-2
          DPC(1) = SIGN( MIN(ABS(DPX),ABS(DPC(1))),DPC(1) )
!          IF( SG_FRC(2,NTX)+DPC(1).LT.0.D+0 ) DPC(1) = 6.D-1*DPC(1)
          SG_FRC(2,NTX) = SG_FRC(2,NTX) + DPC(1)
          IF( SG_FRC(2,NTX).LT.EPSL ) SG_FRC(2,NTX) = 0.D+0
!
!---      Limit changes in gas pressure, excessive primary variable 
!         change limit at 100 MPa  ---
!
          IF( ABS(DPC(1)).GT.1.D+8 ) IERR = 1
          DPX = MAX( 1.D+6,1.D-1*(PG_FRC(2,NTX)-MAX(PL_FRC(2,NTX),
     &      PN_FRC(2,NTX))) )
          DPX = 2.5D-2*(PSO_FRC(2,NTX)+PATM)
          DPC(2) = SIGN( MIN(ABS(DPX),ABS(DPC(2))),DPC(2) )
          PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPC(2)
          PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            IF( YLS_FRC(2,NTX).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_FRC(2,NTX)/EPSL.LT.EPSL .AND. 
     &        DPS/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_FRC(2,NTX) = YLS_FRC(2,NTX)+DPS
            IF( YLS_FRC(2,NTX).LT.1.D-12 ) YLS_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---    Phase condition #4  ---
!
!       SG_FRC = 0.0
!       SL_FRC = 0.0
!       SN_FRC = 1.0
!       Nonaqueous liquid saturated
!
!       water mass - aqueous saturation - SL_FRC
!       gas mass - gas to oil ratio - POSM_FRC
!       oil mass - nonaqueous (system) pressure - PN_FRC,PSO_FRC
!       salt mass - total salt mass - TMS_FRC(2,NTX)
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPW
          IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
!
!---      Limit changes in gas to oil ratio, and zero negative 
!         corrections for zero gas to oil ratios  ---
!
          DPX = 5.D-2/REAL(NITER)
          DPC(1) = SIGN( MIN(ABS(DPX),ABS(DPC(1))),DPC(1) )
          POSM_FRC(2,NTX) = POSM_FRC(2,NTX) + DPC(1)
!
!---      Limit changes in nonaqueous-liquid pressure (system pressure),
!         excessive primary variable change limit at 100 MPa  ---
!
          IF( ABS(DPC(2)).GT.1.D+8 ) IERR = 1
          DPX = MAX( 1.D+6,1.D-1*(PN_FRC(2,NTX)-PL_FRC(2,NTX)) )
          DPX = 2.5D-2*(PSO_FRC(2,NTX)+PATM)
          DPC(2) = SIGN( MIN(ABS(DPX),ABS(DPC(2))),DPC(2) )
!
!---      Bubble-point pressure as a function of gas to oil ratio
!         (m^3/m^3) and temperature (C)  ---
!
!          CALL BO_BP( PBPX,POSM_FRC(2,NTX),T_FRC(2,NTX) )
!
!---      Pressure below bubble-point pressure  ---
!
!          IF( PSO_FRC(2,NTX)+DPC(2).LT.PBPX-PATM ) DPC(2) = 6.D-1*DPC(2)
          PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPC(2)
          PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
!
!---      Isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Zero negative corrections for zero total salt mass  ---
!
            IF( TMS_FRC(2,NTX).LT.EPSL .AND. DPS.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            TMS_FRC(2,NTX) = TMS_FRC(2,NTX) + DPS
            IF( TMS_FRC(2,NTX).LT.EPSL ) TMS_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---    Phase condition #5  ---
!
!       SG_FRC > 0.0
!       SL_FRC = 0.0
!       SN_FRC > 0.0
!       Nonaqueous liquid and gas saturated
!
!       water mass - aqueous saturation, SL_FRC
!       gas mass - gas saturation, SG_FRC
!       oil mass - gas (system) pressure - PG_FRC,PSO_FRC
!       salt mass - total salt mass - TMS_FRC(2,NTX)
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPW
          IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
!
!---      Limit changes in gas saturation  ---
!
          DPX = 2.5D-2
          DPC(1) = SIGN( MIN(ABS(DPX),ABS(DPC(1))),DPC(1) )
!          IF( SG_FRC(2,NTX)+DPC(1).LT.0.D+0 ) DPC(1) = 6.D-1*DPC(1)
          SG_FRC(2,NTX) = SG_FRC(2,NTX) + DPC(1)
          IF( SG_FRC(2,NTX).LT.EPSL ) SG_FRC(2,NTX) = 0.D+0
!
!---      Limit changes in gas pressure, excessive primary variable 
!         change limit at 100 MPa  ---
!
          IF( ABS(DPC(1)).GT.1.D+8 ) IERR = 1
          DPX = MAX( 1.D+6,1.D-1*(PG_FRC(2,NTX)-MAX(PL_FRC(2,NTX),
     &      PN_FRC(2,NTX))) )
          DPX = 2.5D-2*(PSO_FRC(2,NTX)+PATM)
          DPC(2) = SIGN( MIN(ABS(DPX),ABS(DPC(2))),DPC(2) )
          PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPC(2)
          PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
!
!---      Isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Zero negative corrections for zero total salt mass  ---
!
            IF( TMS_FRC(2,NTX).LT.EPSL .AND. DPS.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            TMS_FRC(2,NTX) = TMS_FRC(2,NTX) + DPS
            IF( TMS_FRC(2,NTX).LT.EPSL ) TMS_FRC(2,NTX) = 0.D+0
          ENDIF
        ENDIF
!
!---    Check for excessive pressure or temperature   ---
!
        PX = PSO_FRC(2,NTX)+PATM
        TKX = T_FRC(2,NTX)+TABS
        IF( PX.GT.8.D+8 .OR. PX.LT.0.D+0 ) IERR = 1
        IF( TKX.GT.TCRW .OR. TKX.LT.TABS ) IERR = 1
        IF( IERR.EQ.1 ) NSD(1) = NTX
!
!---  End loop over fault triangles  ---
!
      ENDDO
!
!---  End loop over faults  ---
!
      ENDDO
!
!---  Reduce time step for excessive changes in primary variables   ---
!
      IF( IERR.EQ.1 ) THEN
        ICNV = 1
        NTX = NSD_FRC(1)
        WRITE(ISC,'(10X,A)') '---  Excessive Primary Variable ' // 
     &    'Change  ---'
        WRITE(IWR,'(10X,A)') '---  Excessive Primary Variable ' // 
     &    'Change  ---'
        WRITE(ISC,'(4X,A,I6)') 'Fault Triangle = ',NTX
        WRITE(IWR,'(4X,A,I6)') 'Fault Triangle = ',NTX
        NPX = MOD(NPHAZ_FRC(2,NTX),10)
        WRITE(ISC,'(4X,2A)') 'Phase Condition = ',PH_CND(NPX)
        WRITE(IWR,'(4X,2A)') 'Phase Condition = ',PH_CND(NPX)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Temperature = ',T_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Temperature = ',T_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'System Pressure = ',
     &    PSO_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'System Pressure = ',
     &    PSO_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',
     &    PL_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Pressure = ',
     &    PL_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Pressure = ',
     &    PG_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Pressure = ',
     &    PG_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Nonaqueous Liquid Pressure = ',
     &    PN_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Nonaqueous Liquid Pressure = ',
     &    PN_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Saturation = ',
     &    SL_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Saturation = ',
     &    SL_FRC(2,NTX)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Saturation = ',
     &    SG_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Saturation = ',
     &    SG_FRC(2,NTX)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Nonaqueous-Liquid Saturation = ',
     &    SN_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Nonaqueous-Liquid Saturation = ',
     &    SN_FRC(2,NTX)
      ENDIF
!
!---  Reduce time step  ---
!
  300   CONTINUE
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
!---      Update primary variables, loop over faults  ---
!
          DO NFX = 1,NF_FRC
!
!---      Loop over fault triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
            T_FRC(2,NTX) = T_FRC(1,NTX)
            PL_FRC(2,NTX) = PL_FRC(1,NTX)
            PG_FRC(2,NTX) = PG_FRC(1,NTX)
            PN_FRC(2,NTX) = PN_FRC(1,NTX)
            POSM_FRC(2,NTX) = POSM_FRC(1,NTX)
            PSO_FRC(2,NTX) = PSO_FRC(1,NTX)
            SG_FRC(2,NTX) = SG_FRC(1,NTX)
            SL_FRC(2,NTX) = SL_FRC(1,NTX)
            SN_FRC(2,NTX) = SN_FRC(1,NTX)
            YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
            TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
            NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
            DO IGC = 1,NGC+2
              TMC_FRC(IGC,2,NTX) = TMC_FRC(IGC,1,NTX)
            ENDDO
          ENDDO
          ENDDO
!
!---      Loop over field nodes  ---
!
          DO N = 1,NFBN
            IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            PSO(2,N) = PSO(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
            DO IGC = 1,NGC+2
              TMC(IGC,2,N) = TMC(IGC,1,N)
            ENDDO
          ENDDO
!
!---      Coupled-well pressure  ---
!
          DO NCW = 1,N_CW
            P_CW(2,NCW) = P_CW(1,NCW)
          ENDDO
          NITER = 1
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
!
!---      Update primary variables, loop over faults  ---
!
          DO NFX = 1,NF_FRC
!
!---      Loop over fault triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
            T_FRC(2,NTX) = T_FRC(1,NTX)
            PL_FRC(2,NTX) = PL_FRC(1,NTX)
            PG_FRC(2,NTX) = PG_FRC(1,NTX)
            PN_FRC(2,NTX) = PN_FRC(1,NTX)
            POSM_FRC(2,NTX) = POSM_FRC(1,NTX)
            PSO_FRC(2,NTX) = PSO_FRC(1,NTX)
            SG_FRC(2,NTX) = SG_FRC(1,NTX)
            SL_FRC(2,NTX) = SL_FRC(1,NTX)
            SN_FRC(2,NTX) = SN_FRC(1,NTX)
            YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
            TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
            NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
            DO IGC = 1,NGC+2
              TMC_FRC(IGC,2,NTX) = TMC_FRC(IGC,1,NTX)
            ENDDO
          ENDDO
          ENDDO
!
!---      Loop over field nodes  ---
!
          DO N = 1,NFBN
            IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            PSO(2,N) = PSO(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
            DO IGC = 1,NGC+2
              TMC(IGC,2,N) = TMC(IGC,1,N)
            ENDDO
          ENDDO
!
!---      Coupled-well pressure  ---
!
          DO NCW = 1,N_CW
            P_CW(2,NCW) = P_CW(1,NCW)
          ENDDO
          NITER = 1
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
!---  End of UPDT_FRC_BO group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_FRC_EOR
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
!     Update the primary variables for faults.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Ramesh Sarathi, PNNL, 6 March 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE OUTPU
      USE JACOB
      USE HYST
      USE GRID
      USE GLB_PAR
      USE GEOM_FRC
      USE FILES
      USE FDVS_FRC
      USE FDVS
      USE FDVP_FRC
      USE FDVP
      USE FDVGC_FRC
      USE FDVGC
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
      CHARACTER*128 PH_CND(5)
      REAL*8 DPC(LNGC)
      INTEGER MPC(LNGC)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PH_CND
      DATA PH_CND /'Phase Condition #1, SL = 1.0, SG+SN = 0.0', 
     &  'Phase Condition #2, SL < 1.0, ST > 0.0, SG = 0.0', 
     &  'Phase Condition #3, SL < 1.0, ST > 0.0, SG > 0.0', 
     &  'Phase Condition #4, SL = 0.0, ST = 1.0, SG = 0.0', 
     &  'Phase Condition #5, SL = 0.0, ST = 1.0, SG > 0.0'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_FRC_EOR'
!
!---  Black-oil option  ---
!
      IF( ISLC(68).EQ.0 ) THEN
        CALL UPDT_FRC_BO
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Compositional option  ---
!
      !rss_ask
      IF( ICNV.EQ.1 ) GOTO 300
      IERR = 0
!
!---  Update primary variables, loop over faults  ---
!
      DO NFX = 1,NF_FRC
!
!---    Loop over fault triangles  ---
!
      DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---    Skip inactive triangles  ---
!
        IF( IXP_FRC(NTX).EQ.0 ) CYCLE
!
!---    Skip inactive nodes  ---
!
        N_DB = NTX
        IZN = IZ_FRC(NTX)
!
!---    Energy equation  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          MPT = IM_FRC(IEQT,NTX)
          DPT = BLU(MPT)
        ELSE
          DPT = 0.D+0
        ENDIF
!
!---    Water mass equation  ---
!
        MPW = IM_FRC(IEQW,NTX)
        DPW = BLU(MPW)
!
!---    CO2 mass equation  ---
!
        IF( ISLC(45).EQ.0 ) THEN
          MPC(1) = IM_FRC(IEQGC(1),NTX)
          DPC(1) = BLU(MPC(1))
        ELSE
          MPC(1) = 0
          DPC(1) = 0.D+0
        ENDIF
!
!---    CH4 mass equation  ---
!
        IF( ISLC(64).EQ.0 ) THEN
          MPC(2) = IM_FRC(IEQGC(2),NTX)
          DPC(2) = BLU(MPC(2))
        ELSE
          MPC(2) = 0
          DPC(2) = 0.D+0
        ENDIF
!
!---    Petroleum component mass equation  ---
!
        DO IGC = 3,NGC+2
          MPC(IGC) = IM_FRC(IEQGC(IGC),NTX)
          DPC(IGC) = BLU(MPC(IGC))
        ENDDO
!
!---    Salt mass equation  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = IM_FRC(IEQS,NTX)
          DPS = BLU(MPS)
        ELSE
          DPS = 0.D+0
        ENDIF
!
!---    Limit changes in temperature to 0.5 C,
!       excessive primary variable change limit at 100 C  ---
!
!        IF( ABS(DPT).GT.1.D+2 ) IERR = 1
        DPT = SIGN( MIN( 5.D-1,ABS(DPT) ),DPT )
        T_FRC(2,NTX) = MAX( (T_FRC(2,NTX)+DPT),1.D-1 )
!
!---    Phase condition #1  ---
!
!       SL = 1.0, SG+SN = 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous (system) pressure, Pa gauge - PL,PSO
!       CO2 mass - CO2 vapor pressure - PVA
!       CH4 mass - total CH4 mass - TMC
!       petroleum component mass - total petroleum component mass - TMC
!       salt mass - total salt brine mass fraction - YLS
!
        IF( MOD(NPHAZ_FRC(2,NTX),10).EQ.1 ) THEN
!
!---      Limit changes in total pressure,
!         excessive primary variable change limit at 100 MPa  ---
!
!          IF( ABS(DPW).GT.1.D+8 ) IERR = 1
          DPX = MAX( 2.5D-2*(PSO_FRC(2,NTX)+PATM),1.D+6 )
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPW
          PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
!
!---      Zero negative corrections for zero aqueous CO2  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( PVA_FRC(2,NTX).LT.EPSL .AND. BLU(MPC(1)).LT.EPSL ) THEN
              BLU(MPC(1)) = 0.D+0
              DPC(1) = 0.D+0
            ENDIF
            DPX = 1.D+6
            DPC(1) = SIGN( MIN(ABS(DPX),ABS(DPC(1))),DPC(1) )
            IF( (PVA_FRC(2,NTX)+DPC(1)).LT.0.D+0 ) DPC(1) = 6.D-1*DPC(1)
            PVA_FRC(2,NTX) = MAX( (PVA_FRC(2,NTX)+DPC(1)),0.D+0 )
          ENDIF
!
!---      Zero negative corrections for zero total petroleum
!         component mass  ---
!
          TMWX = PORD_FRC(2,NTX)*RHOL_FRC(2,NTX)*VOL(NTX)
          DPX = 1.D-4*TMWX
          DO IGC = 2,NGC+2
            IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
            IF( TMC_FRC(IGC,2,NTX).LT.EPSL .AND. 
     &        BLU(MPC(IGC)).LT.EPSL ) THEN
              BLU(MPC(IGC)) = 0.D+0
              DPC(IGC) = 0.D+0
            ENDIF
            DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
            TMC_FRC(IGC,2,NTX) 
     &        = MAX( (TMC_FRC(IGC,2,NTX)+DPC(IGC)),0.D+0 )
           IF( TMC_FRC(IGC,2,NTX).LT.1.D-12 ) TMC_FRC(IGC,2,NTX) = 0.D+0
          ENDDO
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            IF( YLS_FRC(2,NTX).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_FRC(2,NTX)/EPSL.LT.EPSL .AND. DPS/EPSL.LT.EPSL ) 
     &        THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_FRC(2,NTX) = YLS_FRC(2,NTX)+DPS
            IF( YLS_FRC(2,NTX).LT.1.D-12 ) YLS_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---    Phase condition #2  ---
!
!       SL < 1.0, ST > 0.0, SG = 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!         endif
!       CH4 mass
!         if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!         endif
!       petroleum component mass
!         if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!         endif
!       salt mass - total salt brine mass fraction - YLS
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.2 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
            DPX = MAX( 1.D+6,1.D-1*(PN_FRC(2,NTX)-PL_FRC(2,NTX)) )
            DPX = 5.D-2*(PSO_FRC(2,NTX)+PATM)/SQRT(REAL(NITER))
            DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
            PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPW
            PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
          ELSE
            DPX = 5.D-2/SQRT(REAL(NITER))
            DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
            SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPW
            IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          SZMCX = 0.D+0
          DO IGC = 1,NGC+2
            IF( ISLC(45).NE.0 .AND. IGC.EQ.1 ) CYCLE
            IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
!
!---        Maximum petroleum component, limit changes in nonaqueous
!           liquid pressure, excessive primary variable change
!           limit at 100 MPa  ---
!
            IF( IGC.EQ.IZMC_FRC(NTX) ) THEN
              IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                DPX = 5.D-2/SQRT(REAL(NITER))
                DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
                SL_FRC(2,NTX) = SL_FRC(2,NT) + DPC(IGC)
                IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
              ELSE
!                IF( ABS(DPC(IGC)).GT.1.D+8 ) IERR = 1
                DPX = MAX( 1.D+6,1.D-1*(PN_FRC(2,NTX)-PL_FRC(2,NTX)) )
                DPX = 5.D-2*(PSO_FRC(2,NTX)+PATM)/SQRT(REAL(NITER))
                DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
                PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPC(IGC)
                PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
              ENDIF
!
!---        Limit changes in petroleum component mole fraction,
!           and zero negative corrections for zero mole fractions---
!
            ELSE
              DPX = 5.D-2/REAL(NITER)
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              ZMC_FRC(IGC,2,NTX) = ZMC_FRC(IGC,2,NTX) + DPC(IGC)
              IF( ZMC_FRC(IGC,2,NTX).LT.EPSL ) 
     &          ZMC_FRC(IGC,2,NTX) = 0.D+0
              SZMCX = SZMCX + ZMC_FRC(IGC,2,NTX)
            ENDIF
          ENDDO
          ZMC_FRC(IZMC_FRC(NTX),2,NTX) = MAX( 1.D+0-SZMCX,0.D+0 )
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            IF( YLS_FRC(2,NTX).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_FRC(2,NTX)/EPSL.LT.EPSL .AND. DPS/EPSL.LT.EPSL ) 
     &        THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_FRC(2,NTX) = YLS_FRC(2,NTX)+DPS
            IF( YLS_FRC(2,NTX).LT.1.D-12 ) YLS_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---    Phase condition #3  ---
!
!       SL < 1.0, ST > 0.0, SG > 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!         endif
!       CH4 mass
!         if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!         endif
!       petroleum component mass
!         if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!         endif
!       salt mass - total salt brine mass fraction - YLS
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.3 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
            DPX = 5.D-2*(PSO_FRC(2,NTX)+PATM)/SQRT(REAL(NITER))
            DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
            PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPW
            PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
          ELSE
            DPX = 5.D-2/SQRT(REAL(NITER))
            DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
            SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPW
            IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          SZMCX = 0.D+0
          DO IGC = 1,NGC+2
            IF( ISLC(45).NE.0 .AND. IGC.EQ.1 ) CYCLE
            IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
!
!---        Maximum petroleum component, limit changes in
!           gas pressure, excessive primary variable change
!           limit at 100 MPa  ---
!
            IF( IGC.EQ.IZMC_FRC(NTX) ) THEN
              IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                DPX = 5.D-2/SQRT(REAL(NITER))
                DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
                SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPC(IGC)
                IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
              ELSE
                DPX = 5.D-2*(PSO_FRC(2,NTX)+PATM)/SQRT(REAL(NITER))
                DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
                PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPC(IGC)
                PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
              ENDIF
!
!---        Limit changes in petroleum component mole fraction,
!           and zero negative corrections for zero mole fractions---
!
            ELSE
              DPX = 5.D-2/REAL(NITER)
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              ZMC_FRC(IGC,2,NTX) = ZMC_FRC(IGC,2,NTX) + DPC(IGC)
              IF( ZMC_FRC(IGC,2,NTX).LT.EPSL ) 
     &          ZMC_FRC(IGC,2,NTX) = 0.D+0
              SZMCX = SZMCX + ZMC_FRC(IGC,2,NTX)
            ENDIF
          ENDDO
          ZMC_FRC(IZMC_FRC(NTX),2,NTX) = MAX( 1.D+0-SZMCX,0.D+0 )
!
!---      Limit salt mass fraction changes to 0.25 of the
!         maximum value if salt mass fraction is less than
!         the maximum, isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_LS( T_FRC(2,NTX),XLSMX )
            IF( YLS_FRC(2,NTX).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero dissolved salt  ---
!
            IF( YLS_FRC(2,NTX)/EPSL.LT.EPSL .AND. DPS/EPSL.LT.EPSL ) 
     &        THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS_FRC(2,NTX) = YLS_FRC(2,NTX)+DPS
            IF( YLS_FRC(2,NTX).LT.1.D-12 ) YLS_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---    Phase condition #4  ---
!
!       SL = 0.0, ST = 1.0, SG = 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!         endif
!       CH4 mass
!         if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!         endif
!       petroleum component mass
!         if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!           - nonaqueous (system) pressure, Pa gauge - PN,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!         endif
!       salt mass - total salt mass - TMS(2,NTX)
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.4 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
            DPX = MAX( 1.D+6,1.D-1*(PN_FRC(2,NTX)-PL_FRC(2,NTX)) )
            DPX = 5.D-2*(PSO_FRC(2,NTX)+PATM)/SQRT(REAL(NITER))
            DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
            PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPW
            PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
          ELSE
            DPX = 5.D-2/SQRT(REAL(NITER))
            DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
            SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPW
            IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          SZMCX = 0.D+0
          DO IGC = 1,NGC+2
            IF( ISLC(45).NE.0 .AND. IGC.EQ.1 ) CYCLE
            IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
!
!---        Maximum petroleum component, limit changes in nonaqueous
!           liquid pressure, excessive primary variable change
!           limit at 100 MPa  ---
!
            IF( IGC.EQ.IZMC_FRC(NTX) ) THEN
              IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                DPX = 5.D-2/SQRT(REAL(NITER))
                DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
                SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPC(IGC)
                IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
              ELSE
                DPX = MAX( 1.D+6,1.D-1*(PN_FRC(2,NTX)-PL_FRC(2,NTX)) )
                DPX = 5.D-2*(PSO(2,N)+PATM)/SQRT(REAL(NITER))
                DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
                PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPC(IGC)
                PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
              ENDIF
!
!---        Limit changes in petroleum component mole fraction,
!           and zero negative corrections for zero mole fractions---
!
            ELSE
              DPX = 5.D-2/REAL(NITER)
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              ZMC_FRC(IGC,2,NTX) = ZMC_FRC(IGC,2,NTX) + DPC(IGC)
              IF( ZMC_FRC(IGC,2,NTX).LT.EPSL ) 
     &          ZMC_FRC(IGC,2,NTX) = 0.D+0
              SZMCX = SZMCX + ZMC_FRC(IGC,2,NTX)
            ENDIF
          ENDDO
          ZMC_FRC(IZMC_FRC(NTX),2,NTX) = MAX( 1.D+0-SZMCX,0.D+0 )
!
!---      Isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Zero negative corrections for zero total salt mass  ---
!
            IF( TMS_FRC(2,NTX).LT.EPSL .AND. DPS.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            TMS_FRC(2,NTX) = TMS_FRC(2,NTX) + DPS
            IF( TMS_FRC(2,NTX).LT.EPSL ) TMS_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---    Phase condition #5  ---
!
!       SL = 0.0, ST >= 0.0, SG > 0.0
!
!       energy - temperature, C - T
!       water mass - aqueous saturation, SL
!       CO2 mass
!         if( ZMC_FRC(1,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 1
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(1,2,NTX)
!         endif
!       CH4 mass
!         if( ZMC_FRC(2,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= 2
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(2,2,NTX)
!         endif
!       petroleum component mass
!         if( ZMC_FRC(I,2,NTX) > max( ZMC_FRC(J,2,NTX) ) J /= I
!           - gas (system) pressure, Pa gauge - PG,PSO
!         else
!           - nonaqueous mole fraction of CO2 - ZMC_FRC(I,2,NTX)
!         endif
!       salt mass - total salt mass - TMS(2,NTX)
!
        ELSEIF( MOD(NPHAZ_FRC(2,NTX),10).EQ.5 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
            DPX = MAX( 1.D+6,1.D-1*(PG_FRC(2,NTX) - 
     &        MAX(PL_FRC(2,NTX),PN_FRC(2,NTX))) )
            DPX = 5.D-2*(PSO_FRC(2,NTX)+PATM)/SQRT(REAL(NITER))
            DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
            PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPW
            PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
          ELSE
            DPX = 5.D-2/SQRT(REAL(NITER))
            DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
            SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPW
            IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
          ENDIF
!
!---      Loop over all petroleum components, including CO2 and CH4  ---
!
          SZMCX = 0.D+0
          DO IGC = 1,NGC+2
            IF( ISLC(45).NE.0 .AND. IGC.EQ.1 ) CYCLE
            IF( ISLC(64).NE.0 .AND. IGC.EQ.2 ) CYCLE
!
!---        Maximum petroleum component, limit changes in
!           gas pressure, excessive primary variable change
!           limit at 100 MPa  ---
!
            IF( IGC.EQ.IZMC_FRC(NTX) ) THEN
              IF( NPHAZ_FRC(2,NTX).GT.100 ) THEN
                DPX = 5.D-2/SQRT(REAL(NITER))
                DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
                SL_FRC(2,NTX) = SL_FRC(2,NTX) + DPC(IGC)
                IF( SL_FRC(2,NTX).LT.EPSL ) SL_FRC(2,NTX) = 0.D+0
              ELSE
                DPX = MAX( 1.D+6,1.D-1*(PG_FRC(2,NTX) - 
     &            MAX(PL_FRC(2,NTX),PN_FRC(2,NTX))) )
                DPX = 5.D-2*(PSO_FRC(2,NTX)+PATM)/SQRT(REAL(NITER))
                DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
                PSO_FRC(2,NTX) = PSO_FRC(2,NTX) + DPC(IGC)
                PSO_FRC(2,NTX) = MIN( PSO_FRC(2,NTX),5.D+8 )
              ENDIF
!
!---        Limit changes in petroleum component mole fraction,
!           and zero negative corrections for zero mole fractions---
!
            ELSE
              DPX = 5.D-2/REAL(NITER)
              DPC(IGC) = SIGN( MIN(ABS(DPX),ABS(DPC(IGC))),DPC(IGC) )
              ZMC_FRC(IGC,2,NTX) = ZMC_FRC(IGC,2,NTX) + DPC(IGC)
              IF( ZMC_FRC(IGC,2,NTX).LT.EPSL ) 
     &          ZMC_FRC(IGC,2,NTX) = 0.D+0
              SZMCX = SZMCX + ZMC_FRC(IGC,2,NTX)
            ENDIF
          ENDDO
          ZMC_FRC(IZMC_FRC(NTX),2,NTX) = MAX( 1.D+0-SZMCX,0.D+0 )
!
!---      Isobrine option   ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Zero negative corrections for zero total salt mass  ---
!
            IF( TMS_FRC(2,NTX).LT.EPSL .AND. DPS.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            TMS_FRC(2,NTX) = TMS_FRC(2,NTX) + DPS
            IF( TMS_FRC(2,NTX).LT.EPSL ) TMS_FRC(2,NTX) = 0.D+0
          ENDIF
        ENDIF
!
!---    Check for excessive pressure or temperature   ---
!
        PX = PSO_FRC(2,NTX)+PATM
        TKX = T_FRC(2,NTX)+TABS
        IF( PX.GT.8.D+8 .OR. PX.LT.0.D+0 ) IERR = 1
        IF( TKX.GT.TCRW .OR. TKX.LT.TABS ) IERR = 1
        IF( IERR.EQ.1 ) NSD_FRC(1) = NTX
!
!---  End loop over fault triangles  ---
!
      ENDDO
!
!---  End loop over faults  ---
!
      ENDDO
!
!---  Reduce time step for excessive changes in primary variables   ---
!
      IF( IERR.EQ.1 ) THEN
        ICNV = 1
        NTX = NSD_FRC(1)
        WRITE(ISC,'(10X,A)') '---  Excessive Fault Primary Variable ' //
     &    'Change  ---'
        WRITE(IWR,'(10X,A)') '---  Excessive Fault Primary Variable ' //
     &    'Change  ---'
        WRITE(ISC,'(4X,A,I6)') 'Fault Triangle = ',NTX
        WRITE(IWR,'(4X,A,I6)') 'Fault Triangle = ',NTX
        NPX = MOD(NPHAZ_FRC(2,NTX),10)
        WRITE(ISC,'(4X,2A)') 'Fault Phase Condition = ',PH_CND(NPX)
        WRITE(IWR,'(4X,2A)') 'Fault Phase Condition = ',PH_CND(NPX)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Fault Temperature = ',T_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Fault Temperature = ',T_FRC(2,NTX)
        WRITE(ISC,'(4X,A,1PE12.5)') 
     &    'Fault System Pressure = ',PSO_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 
     &    'Fault System Pressure = ',PSO_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 
     &    'Fault Aqueous Pressure = ',PL_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 
     &    'Fault Aqueous Pressure = ',PL_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Fault Gas Pressure = ',
     &    PG_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Fault Gas Pressure = ',
     &    PG_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 
     &    'Fault Nonaqueous Liquid Pressure = ',PN_FRC(2,NTX)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 
     &    'Fault Nonaqueous Liquid Pressure = ',PN_FRC(2,NTX)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 
     &    'Fault Aqueous Saturation = ',SL_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)') 
     &    'Fault Aqueous Saturation = ',SL_FRC(2,NTX)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Fault Gas Saturation = ',
     &    SG_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Fault Gas Saturation = ',
     &    SG_FRC(2,NTX)
        WRITE(ISC,'(4X,A,1PE12.5)') 
     &    'Fault Nonaqueous-Liquid Saturation = ',SN_FRC(2,NTX)
        WRITE(IWR,'(4X,A,1PE12.5)') 
     &    'Fault Nonaqueous-Liquid Saturation = ',SN_FRC(2,NTX)
      ENDIF
!
!---  Reduce time step  ---
!
  300   CONTINUE
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
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),
     &      ' to ',VAR,UNTM(1:NCH)
          WRITE(IWR,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)') 
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),
     &      ' to ',VAR,UNTM(1:NCH)
!
!---      Loop over faults  ---
!
          DO NFX = 1,NF_FRC
!
!---      Loop over fault triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
            DO M = 1,6
              BETA_FRC(M,NTX) = 0.D+0
            ENDDO
            IBETA_FRC(NTX) = 0
            T_FRC(2,NTX) = T_FRC(1,NTX)
            PL_FRC(2,NTX) = PL_FRC(1,NTX)
            PG_FRC(2,NTX) = PG_FRC(1,NTX)
            PN_FRC(2,NTX) = PN_FRC(1,NTX)
            PSO_FRC(2,NTX) = PSO_FRC(1,NTX)
            PVA_FRC(2,NTX) = PVA_FRC(1,NTX)
            SG_FRC(2,NTX) = SG_FRC(1,NTX)
            SL_FRC(2,NTX) = SL_FRC(1,NTX)
            SN_FRC(2,NTX) = SN_FRC(1,NTX)
            YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
            TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
            NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
            DO IGC = 1,NGC+2
              TMC_FRC(IGC,2,NTX) = TMC_FRC(IGC,1,NTX)
              ZMC_FRC(IGC,2,NTX) = ZMC_FRC(IGC,1,NTX)
            ENDDO
!
!---      End loop over fault triangles  ---
!
          ENDDO
!
!---      End loop over faults  ---
!
          ENDDO
!
!---      Loop over field nodes  ---
!
          DO N = 1,NFBN
            IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
            DO M = 1,6
              BETA(M,N) = 0.D+0
            ENDDO
            IBETA(N) = 0
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            PSO(2,N) = PSO(1,N)
            PVA(2,N) = PVA(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
            DO IGC = 1,NGC+2
              TMC(IGC,2,N) = TMC(IGC,1,N)
              ZMC(IGC,2,N) = ZMC(IGC,1,N)
            ENDDO
          ENDDO
!
!---      Coupled-well pressure  ---
!
          DO NCW = 1,N_CW
            P_CW(2,NCW) = P_CW(1,NCW)
          ENDDO
          NITER = 1
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
!
!---      Loop over faults  ---
!
          DO NFX = 1,NF_FRC
!
!---      Loop over fault triangles  ---
!
          DO NTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
!
!---        Skip inactive triangles  ---
!
            IF( IXP_FRC(NTX).EQ.0 ) CYCLE
            DO M = 1,6
              BETA_FRC(M,NTX) = 0.D+0
            ENDDO
            IBETA_FRC(NTX) = 0
            T_FRC(2,NTX) = T_FRC(1,NTX)
            PL_FRC(2,NTX) = PL_FRC(1,NTX)
            PG_FRC(2,NTX) = PG_FRC(1,NTX)
            PN_FRC(2,NTX) = PN_FRC(1,NTX)
            PSO_FRC(2,NTX) = PSO_FRC(1,NTX)
            PVA_FRC(2,NTX) = PVA_FRC(1,NTX)
            SG_FRC(2,NTX) = SG_FRC(1,NTX)
            SL_FRC(2,NTX) = SL_FRC(1,NTX)
            SN_FRC(2,NTX) = SN_FRC(1,NTX)
            YLS_FRC(2,NTX) = YLS_FRC(1,NTX)
            TMS_FRC(2,NTX) = TMS_FRC(1,NTX)
            NPHAZ_FRC(2,NTX) = NPHAZ_FRC(1,NTX)
            DO IGC = 1,NGC+2
              TMC_FRC(IGC,2,NTX) = TMC_FRC(IGC,1,NTX)
              ZMC_FRC(IGC,2,NTX) = ZMC_FRC(IGC,1,NTX)
            ENDDO
            
!
!---      End loop over fault triangles  ---
!
          ENDDO
!
!---      End loop over faults  ---
!
          ENDDO
!
!---      Loop over field nodes  ---
!
          DO N = 1,NFBN
            IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
            DO M = 1,6
              BETA(M,N) = 0.D+0
            ENDDO
            IBETA(N) = 0
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            PSO(2,N) = PSO(1,N)
            PVA(2,N) = PVA(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            TMS(2,N) = TMS(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
            DO IGC = 1,NGC+2
              TMC(IGC,2,N) = TMC(IGC,1,N)
              ZMC(IGC,2,N) = ZMC(IGC,1,N)
            ENDDO
          ENDDO
!
!---      Coupled-well pressure  ---
!
          DO NCW = 1,N_CW
            P_CW(2,NCW) = P_CW(1,NCW)
          ENDDO
          NITER = 1
          WRITE(ISC,'(10X,A)') 
     &      '---  Time Step Reduction Limit Exceeded---'
          WRITE(IWR,'(10X,A)') 
     &      '---  Time Step Reduction Limit Exceeded---'
          ICNV = 4
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_FRC_EOR group
!
      RETURN
      END

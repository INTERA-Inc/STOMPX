!----------------------Program-----------------------------------------!
!
      PROGRAM STOMP
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
!     STOMP: Subsurface Transport Over Multiple Phases
!
!     Geothermal Mode (STOMP-GT)
!
!     This engineering program numerically simulates thermal
!     and hydrologic transport phenomena in variably saturated
!     subsurface environments, contaminated with a water immiscible
!     volatile organic compound.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!

!
!----------------------Lis Modules-----------------------------------!
!
      USE STOMP_LIS_MODULE







!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE UCODE
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PLT_ATM
      USE PARM_FRC
      USE PARM_BH
      USE OUTPU
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE GEO_MECH
      USE FILES
      USE FDVT
      USE FDVP
      USE FDVG
      USE COUP_WELL
      USE CONST
      USE BCVT
      USE BCVP
      USE BCVG
      USE BCV
      
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)



















































	









































	



















































































































































































































































































































































































































      































































































!
!----------------------Type Declarations-------------------------------!
!
      LOGICAL HALT,PLOT,RESTART




      integer :: IERR

!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = 1

!
!---  Initialize Lis ---
!
      CALL lis_init_f(IERR)

!
!---  Read input file to determine memory requirements  ---
!
      CALL STEP
!
!---  Allocate memory  ---
!
      CALL ALLOC

!
!---  Number of grid cells exceeds binary file size limit  ---
!
      IF( ISLC(67).EQ.1 .AND. LFLD.GT.12000000 ) THEN
        INDX = 3
        CHMSG = 'Number of Grid Cells Yields Binary Files > 2.4 GB'
        CALL WRMSGS( INDX )
      ENDIF
      SUB_LOG(1) = 'STOMP-GT'
      ICODE = 3
!
!---  Intialize variables in common blocks and open files  ---
!
      CALL INTLZ
!
!---  Print banner on screen and output file  ---
!
      CALL BANNER
!
!---  Pressurized crack problem  ---
!
      CALL CROUCH
!
!---  Read user input and restart files  ---
!
      CALL RDINPT_GT
!!
!!---  Write a zonation file for GTO-CCS Problem #4  --
!!
!      OPEN(UNIT=26, FILE='zon_prb_4.dat', STATUS='UNKNOWN', 
!     &  FORM='FORMATTED')
!      CLOSE(UNIT=26,STATUS='DELETE')
!      OPEN(UNIT=26, FILE='zon_prb_4.dat', STATUS='NEW', 
!     &  FORM='FORMATTED')
!      DO 10 N = 1,NFLD
!        IF( JD(N).GT.1 ) THEN
!          WRITE(26,'(A)') "2"
!        ELSE
!          RADX = SQRT( (XP(N)**2) + ((ZP(N)-(-2715.D+0))**2) )
!          IF( RADX.GT.6.D+1 ) THEN
!            WRITE(26,'(A)') "2"
!          ELSE
!            WRITE(26,'(A)') "1"
!          ENDIF
!        ENDIF
!   10 CONTINUE
!      CLOSE(26)
!      STOP
!!
!!---  Write a zonation file for EGS-Collab test bed 2 thermal
!!     breakthrough  --
!!
!      OPEN(UNIT=26, FILE='zonation.dat', STATUS='UNKNOWN', 
!     &  FORM='FORMATTED')
!      CLOSE(UNIT=26,STATUS='DELETE')
!      OPEN(UNIT=26, FILE='zonation.dat', STATUS='NEW', 
!     &  FORM='FORMATTED')
!      DO N = 1,NFLD
!        IF( KD(N).GT.1 ) THEN
!          WRITE(26,'(A)') "1"
!        ELSE
!          RADX = SQRT( (XP(N)**2) + (YP(N)**2) )
!          IF( RADX.GT.10.1D+0 ) THEN
!            WRITE(26,'(A)') "1"
!          ELSE
!            CALL RANDOM_NUMBER( VARX )
!            INDICATOR = 2+INT( VARX*5.D+0 )
!            WRITE(26,'(I1)') INDICATOR
!          ENDIF
!        ENDIF
!      ENDDO
!      CLOSE(26)
!      STOP
!
!---  Assign local stress states for EGS Collab Stanford 2021 paper  ---
!
      DO N = 1,NFLD
        IF( IXP(N).EQ.0 ) CYCLE
        SIG_GM(1,N) = 35.5D+06
        SIG_GM(2,N) = 21.7D+06
        SIG_GM(3,N) = 41.8D+06
        SIG_GM(4,N) = 0.D+0
        SIG_GM(5,N) = 0.D+0
        SIG_GM(6,N) = 0.D+0
      ENDDO
!
!---  Check for internal boundary surfaces and write connectivity
!     list file  --
!
      CALL CONNLST
!
!---  Surface cover connection map  --
!
      IF( NSFCN.GT.0 ) CALL CONNMAP_SFC
!
!---  For geomechanics simulations create a finite-element node map  --
!
      IF( ISLC(50).NE.0 ) CALL CONNFEN
!
!---  Define ground-loop well nodes, check ground-loop well trajectory,
!     and write well.dat file  ---
!
      IF( L_CW.EQ.1 ) THEN
        CALL CHK_GRLP_WELL
        CALL WR_GRLP_WELL
      ENDIF
!
!---  For geomechanics simulations check and preprocess boundary
!     conditions, and set the reference volumetric stress from
!     the initial displacements stored in the restart file  ---
!
      IF( ISLC(50).NE.0 ) CALL CHK_GM
!
!---  Check initial conditions type boundary conditions  ---
!
      CALL CHK_BC_GT
!
!---  Check thermodynamic and hydrologic initial states  ---
!
      CALL CHK_GT
!
!---  Check thermodynamic and hydrologic fracture and borehole
!     initial states  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL CHK_FRC_GT
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL CHK_BH_GT
!
!---  Write bin input files for parallel processing  ---
!
      IF( ISLC(67).EQ.1 ) THEN
        CALL INCRM_GT
        CALL PROP_GT
!
!---    For geomechanics simulations compute Jacobian 
!       matrix pointers  --
!
        IF( ISLC(50).NE.0 ) CALL JCBP_GM_PPC
        CALL WRITE_BIN_GT
        WRITE(IWR,'(A)') 'NOTE: Preprocessing for STOMPX-GT'
        WRITE(ISC,'(A)') 'NOTE: Preprocessing for STOMPX-GT'
!
!---    End parallel pre-processing unless output of the flow or
!       transport linear system is requested  ---
!
        IF( ISLC(34).EQ.0 ) THEN
          WRITE(IWR,'(/,A,/)') '---  End of STOMPX-GT Preprocessing ---'
          WRITE(ISC,'(/,A,/)') '---  End of STOMPX-GT Preprocessing ---'
          STOP
        ENDIF
      ENDIF
!
!---  For geomechanics set k iterate value of pore pressure  ---
!
      IF( ISLC(50).NE.0 ) THEN
        INDX = 2
        CALL PRESS_GM( INDX )
      ENDIF
!
!---  Verify the algorithms for the IAPWS Industrial Formulation
!     1997 for the Thermodynamic Properties of Water and Steam. Lucerne, 
!     Switzerland, August 2007, available at http://www.iapws.org---
!
!      CALL V_IAPWS

!
!---  Sequence reaction equations  ---
!
      IF( ISLC(40).EQ.1 ) CALL SEQEQ

!
!---  Load old time step arrays  ---
!
      CALL LDO_GT
!
!---  Compute Jacobian matrix pointers  ---
!
      CALL JCBP
!
!---  For geomechanics simulations compute Jacobian matrix pointers  --
!
      IF( ISLC(50).NE.0 .AND. ISLC(67).EQ.0 ) CALL JCBP_GM
!
!---  Compute primary variable increments  ---
!
      CALL INCRM_GT
!
!---  Compute fracture and borehole primary variable increments  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL INCRM_FRC_GT 
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL INCRM_BH_GT
!
!---  Compute surface cover primary variable increments  --
!
      IF( NSFCN.GT.0 ) CALL INCRM_SFC
!
!---  Initial hydrologic and thermodynamic properties for active
!     field nodes  ---
!
      CALL PROP_GT
!
!---  Initial hydrologic and thermodynamic properties for active
!     fracture triangles and borehole nodes  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL PROP_FRC_GT
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL PROP_BH_GT
!
!---  Compute initial solute concentrations  ---
!
      CALL CISC_GT
!
!---  Compute initial solute concentrations for active
!     fracture triangles and borehole nodes  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL CISC_FRC_GT
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL CISC_BH_GT

!
!---  Reactive transport  ---
!
      IF( ISLC(40).EQ.1 ) THEN
!
!---    Convert initial reactive species concentrations to
!       node volume basis, mol/m^3  ---
!
        CALL FLHSP
!
!---    Temporarily store time stepping  ---
!
        DT_RST = DT
        DTI_RST = DTI
        TM_RST = TM
!
!---    Loop over number of conservation component species  ---
!
        DO 20 NEQ = 1,NEQC
          NSL = NEQ + NSOLU
!
!---      Mobile conservation component fractions   ---
!
          CALL MOBCF( NEQ )
!
!---      Add immobile conservation component fractions   ---
!
          CALL IMOBCF( NEQ )
!
!---    End of conservation component species transport  ---
!
   20   CONTINUE
!
!---    Loop over number of kinetic component species  ---
!
        DO 40 NEQ = 1,NEQK
          NSL = NEQ + NEQC + NSOLU
! 
!---      Mobile kinetic component fractions   ---
!
          CALL MOBKF( NEQ )
! 
!---      Add immobile kinetic component fractions   ---
!
          CALL IMOBKF( NEQ )
!
!---    End of conservation component species transport  ---
!
   40   CONTINUE
!
!---    Equilibrium-conservation-kinetic reaction chemistry   ---
!
        CALL ECKECHEM
!
!---    Reconstitute mineral species concentrations for initial
!       output  ---
!
        CALL RMNSP
!
!---    Reset time stepping  ---
!
        DT = DT_RST
        DTI = DTI_RST
        TM = TM_RST
      ENDIF

!
!---  Initialize SPLIB routines  ---
!
      IF( ILES.EQ.3 ) THEN
        INDX = -1
        CALL PSPLIB( 0,INDX )
      ENDIF

!
!---  Create Lis matrix, solver, and solution and problem vectors  ---
!
      IF( ILES.EQ.4 ) THEN
!
!---    Solver for coupled flow  ---
!
        INDX = 0
        CALL STOMP_LIS_CREATE(ISVC,F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,INDX)
!
!---    Solver for solute/species transport  ---
!

        NSL = NEQ + NSOLU



        IF( NSL.GT.0 ) THEN
          INDX = 1
          CALL STOMP_LIS_CREATE(0,T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,INDX)
        ENDIF
!
!---    Solver for geomechanics  ---
!
        IF( ISLC(50).NE.0 ) THEN
          INDX = 2
          CALL STOMP_LIS_CREATE(0,G_KSP,G_MAT,G_RHS_VEC,G_SOL_VEC,INDX)
        ENDIF
      ENDIF

!
!---  Initialize geomechanics  ---
!
      IF( ISLC(50).NE.0 ) THEN
!
!---    Reference state porothermoelastic geomechanics; first call
!       to STATIC_GM eliminates reference boundary conditions  ---
!
        IREF_GM = 1
        CALL STATIC_GM
!
!---    Load reference displacements at finite elment nodes  ---
!
        CALL LDDISP_GM
!
!---    Reference volumetric stresses at finite element centroids  ---
!
        IF( ISLC(50).LT.0 ) THEN
          INDX = 0
          CALL VOLSS_GM( INDX )
!
!---      Remove restart check for geomechanics options  ---
!
          ISLC(50) = ABS(ISLC(50))
        ENDIF
!
!---    Static porothermoelastic geomechanics  ---
!
        IREF_GM = 0
        CALL STATIC_GM
!
!---    Set k iterate value of pore pressure and volumetric stress
!
        INDX = 2
        CALL PRESS_GM( INDX )
        CALL VOLSS_GM( INDX )
      ENDIF
!
!---  Boundary saturation, relative permeability, and
!     thermodynamic properties  ---
!
      CALL BCP_GT
!
!---  Initial matrix, fracture, surface-cover and borehole fluxes  ---
!
      ISVF = 1
      ISCVF = 1
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL LDO_BH_GT
      CALL FLUX_GT
      ISVF = 2*ISVC+1
      ISCVF = ISVC+6
!
!---  Surface flux integrator for zero time step  ---
!
      DTX = DT
      DT = 0.D+0
      CALL SFIN
      DT = DTX
!
!---  New Time Step ---
!
  100 CONTINUE
!
!---  Reference node(s) output  ---
!
      IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 .OR.
     &  MOD( (NSTEP-NRST),IFQO ).EQ.0 ) CALL REFNOD
!
!---  Load old time step arrays  ---
!
      CALL LDO_GT
!
!---  Load old time step arrays for active fracture triangles
!     and borehole nodes  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL LDO_FRC_GT
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL LDO_BH_GT
!
!---  Update porosity and permeability in response to geomechanical
!     stress  ---
!
      IF( ISLC(50).NE.0 .AND. NSTEP.EQ.0 ) THEN
        CALL PORSTY_GM
        CALL PERMRF_GM
      ENDIF
!
!---  Load old time step arrays for the volumetric stress
!     and pore pressure  ---
!
      IF( ISLC(50).NE.0 ) THEN
        INDX = 1
        CALL LD_GM( INDX )
      ENDIF

!
!---  Normalize mineral species concentrations after initial
!     output for normal simulations  ---
!
      IF( (NSTEP-NRST).EQ.0 ) CALL NMNSP

!
!---  End of initial conditions simulations  ---
!
      IF( IEO.EQ.3 ) THEN
        INDX = 1
        CHMSG = 'Simulation Stopped:  Initial Condition'
        CALL WRMSGS( INDX )
        GOTO 900
      ENDIF
!
!---  Stop simulation if simulation time exceeds limit  ---
!
      IF( ABS(TMMX-TM).LE.1.D-6 ) THEN
        INDX = 1
        CHMSG = 'Simulation Stopped:  Simulation Time Limit'
        CALL WRMSGS( INDX )
        GOTO 900
      ENDIF
!
!---  Stop simulation if file "stop_stomp" exists  ---
!
      INQUIRE( FILE="stop_stomp", EXIST=HALT )
      IF( HALT ) THEN
        OPEN( UNIT=19, FILE="stop_stomp" )
        CLOSE( UNIT=19, STATUS='DELETE' )
        INDX = 1
        CHMSG = 'Simulation Stopped:  User Interrupt'
        CALL WRMSGS( INDX )
        ISLC(18) = 0
        GOTO 900
      ENDIF

!
!---  Generate plot file if file "plot_stomp" exists  ---
!
      INQUIRE( FILE="plot_stomp", EXIST=PLOT )
      IF( PLOT ) THEN
        OPEN( UNIT=19, FILE="plot_stomp" )
        CLOSE( UNIT=19, STATUS='DELETE' )
        CALL WRPLOT
        IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL WRPLOT_FRC
        IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL WRPLOT_BH
        IF( ISLC(18).LT.1 ) CALL WRRST 
      ENDIF
!
!---  Generate restart file if file "restart_stomp" exists  ---
!
      INQUIRE( FILE="restart_stomp", EXIST=RESTART )
      IF( RESTART ) THEN
        OPEN( UNIT=19, FILE="restart_stomp" )
        CLOSE( UNIT=19, STATUS='DELETE' )
        CALL WRRST
      ENDIF
!
!---  Restart and plot file outputs  ---
!
      IF( ABS(TMPR-TM).LE.1.D-6 ) THEN
        CALL WRPLOT
        IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL WRPLOT_FRC
        IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL WRPLOT_BH
        IF( ISLC(18).LT.1 ) CALL WRRST
      ENDIF
!
!---  Inverse output  ---
!
      IF( ISLC(20).EQ.1 .AND. ABS(TMOB-TM).LT.(1.D+1*EPSL) ) CALL WROBDA
!
!---  Compute the next time step and increment time step counter  ---
!
      DTSO = DT
      CALL TMSTEP
      IF( NSTEP.EQ.0 ) DTSO = DT
      NSTEP = NSTEP + 1
      IF( NSTEP-NRST.GT.MXSTEP ) THEN
        INDX = 1
        CHMSG = 'Simulation Stopped:  Time Step Limit'
        CALL WRMSGS( INDX )
        NSTEP = NSTEP - 1
        GOTO 900
      ENDIF
!
!---  No flow solution  ---
!
      IF( ISLC(47).EQ.1 ) THEN
        CALL BCP_GT
        GOTO 600
      ENDIF
      NTSR = 0
!
!---  Top of sequential flow and transport and geomechanics  ---
!
      K_GM(1) = 0
      K_GM(2) = 0
  110 CONTINUE
      K_GM(1) = K_GM(1) + 1
!
!---  Newton-Raphson iteration restart  ---
!
  200 CONTINUE
      NITER = 0
!
!---  Newton-Raphson iteration start  ---
!
  300 CONTINUE
      NITER = NITER + 1
      K_GM(2) = K_GM(2) + 1
!
!---  Compute source contributions  ---
!
      CALL SORC_GT
!
!---  Compute fracture source contributions  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL SORC_FRC_GT
!
!---  Compute borehole source contributions  ---
!
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL SORC_BH_GT
!
!---  Compute ground-loop-well energy contributions  ---
!
      DO NCW = 1,N_GLW
        CALL GRLP_WELL( NCW )
      ENDDO
!
!---  Compute boundary saturation, relative permeability, and
!     thermodynamic properties  ---
!
      CALL BCP_GT
!
!---  Matrix, fracture, surface-cover and borehole fluxes  ---
!
!      time_flux_start = OMP_GET_WTIME()
      CALL FLUX_GT
!      time_flux_end = OMP_GET_WTIME()
!
!---  Zero Jacobian matrix  ---
!



      INDX = 5
      CALL JCBZ( ISVC,MUC,MLC,MKC,INDX )
!
!---  Load Jacobian matrix for matrix, fractures, and boreholes  ---
!
!      time_jcb_start = OMP_GET_WTIME()
      CALL JCB_GT
!      time_jcb_end = OMP_GET_WTIME()
!
!---  Linear equation solver  ---
!
      IF( ILES.EQ.1 ) THEN
        INDX = 0
        CALL BAND( ISVC,MUC,MLC,INDX )
      ELSEIF( ILES.EQ.3 ) THEN
        INDX = 0
        CALL PSPLIB( ISVC,INDX )

      ELSEIF( ILES.EQ.4 ) THEN
!        time_lis_start = OMP_GET_WTIME()
        INDX = 0
        CALL STOMP_LIS_SOLVE(ISVC,F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,INDX)
!        time_lis_end = OMP_GET_WTIME()







      ENDIF
!
!---  Update primary variables  ---
!
      CALL UPDT_GT
!
!---  Compute convergence from maximum relative residuals  ---
!
      CALL RSDL_GT
!
!---  Increment primary variables  ---
!
!      time_incrm_start = OMP_GET_WTIME()
      CALL INCRM_GT
!      time_incrm_end = OMP_GET_WTIME()
!
!---  Secondary field variables  ---
!
!      time_prop_start = OMP_GET_WTIME()
      CALL PROP_GT
!      time_prop_end = OMP_GET_WTIME()
!      print *,'flux = ',time_flux_end-time_flux_start,
!     &  'jcb = ',time_jcb_end-time_jcb_start,
!     &  'lis = ',time_lis_end-time_lis_start,
!     &  'incrm = ',time_incrm_end-time_incrm_start,
!     &  'prop = ',time_prop_end-time_prop_start
!
!---  Fracture flow and transport solution  ---
!
      IF( ISLC(74).NE.0 ) THEN
!
!---    Update fracture primary variables  ---
!
        IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL UPDT_FRC_GT
!
!---    Update borehole primary variables  ---
!
        IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL UPDT_BH_GT
!
!---    Compute borehole and fracture flow convergence from 
!       maximum relative residuals  ---
!
        CALL RSDL_BF_GT
!
!---    Increment fracture primary variables  ---
!
        IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL INCRM_FRC_GT
!
!---    Increment borehole primary variables  ---
!
        IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL INCRM_BH_GT
!
!---    Fracture secondary field variables  ---
!
        IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL PROP_FRC_GT
!
!---    Borehole secondary field variables  ---
!
        IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL PROP_BH_GT
      ENDIF
!
!---  For geomechanics simulations alter permeability with
!     porosity  --
!
      IF( ISLC(50).NE.0 ) CALL PERMRF_GM
!      PRINT *,'NSTEP = ',NSTEP,' NITER = ',NITER
!      DO M = 1,ISVC
!        PRINT *,'RSD(',M,') = ',RSD(M),
!     &    'NSD(',M,') = ',NSD(M),' ICNV = ',ICNV
!        PRINT *,'RSD_FRC(',M,') = ',RSD_FRC(M),
!     &    'NSD_FRC(',M,') = ',NSD_FRC(M),' ICNV = ',ICNV
!        PRINT *,'RSD_BH(',M,') = ',RSD_BH(M),
!     &    'NSD_BH(',M,') = ',NSD_BH(M),' ICNV = ',ICNV
!      ENDDO
      GOTO( 200,300,600,900 ) ICNV
  600 CONTINUE
!
!---  Solve geomechanics  ---
!
      IF( ISLC(50).NE.0 ) THEN
!
!---    Set k+1 iterate value of pore pressure  ---
!
        INDX = 3
        CALL PRESS_GM( INDX )
!
!---    Static porothermoelastic geomechanics  ---
!
        CALL STATIC_GM
!
!---    Convergence check for sequential coupled flow and transport
!       and geomechanics  ---
!
        CALL RSDL_GM
        IF( RSD_GM.GT.RSDM_GM(IEPD) ) THEN
!
!---      Load k level arrays for the volumetric stress
!         and pore pressure  ---
!
          INDX = 2
          CALL LD_GM( INDX )
!
!---      Update porosity and permeability for geomechical stress  ---
!
          CALL PORSTY_GM
          CALL PERMRF_GM
          GOTO 110
        ENDIF
!
!---    Update porosity and permeability for geomechical stress  ---
!
        CALL PORSTY_GM
        CALL PERMRF_GM
      ENDIF
!
!---  Compute current fluxes for transport solutions or flux
!     integrations  ---
!
      ISVF = 1
!
!---  Matrix, fracture, surface-cover and borehole fluxes  ---
!
      CALL FLUX_GT
!
!---  Compute Local Courant Numbers  ---
!
      IF( ICRNT.EQ.1 ) CALL CRNTNB
!
!---  Compute local Courant numbers for boreholes and fractures  ---
!
      IF( ICRNT.EQ.1 .AND. ISLC(74).NE.0 ) CALL CRNTNB_BF
      ISVF = 2*ISVC+1
!
!---  Beginning of transport equation solution, skip for
!     no transport simulations or transport solution switched off  ---
!
      IF( IEQC.EQ.0 .OR. ISLC(81).EQ.1 ) GOTO 800
!
!---  PEST output for EGS Collab Experiment 1, set time for
!     start of transport calculations  ---
!
      IF( ISLC(80).EQ.1 ) THEN
        IF( R_OBDT(5,1)/EPSL.LT.EPSL ) R_OBDT(5,1) = TM
      ENDIF
!
!---  Loop over number of solutes  ---
!
      DO 700 NSL = 1,NSOLU
!
!---    Skip transport for stationary solutes  ---
!
        IF( IEDL(NSL).EQ.4 ) CYCLE
!
!---    Courant number limiting  ---
!
        N_CRN(NSL) = 1
        IF( ISLC(17).NE.0 ) CALL CRN_LIM( NSL )
!
!---    Fracture flow and transport solution  ---
!
        IF( ISLC(17).NE.0 .AND. (ISLC(74).EQ.1 .OR. ISLC(74).EQ.3) ) 
     &    CALL CRN_LIM_FRC( NSL )
        IF( ISLC(17).NE.0 .AND. (ISLC(74).EQ.2 .OR. ISLC(74).EQ.3) ) 
     &    CALL CRN_LIM_BH( NSL )
!
!---    Sub-time step loop  ---
!
        DO NC = 1,N_CRN(NSL)
          IF( ISLC(17).NE.0 ) TM = MIN( TM+DT,TM_CRN )
!
!---      Compute solute mole fractions ---
!
          CALL SPRP_GT( NSL )
!
!---      Compute solute mole fractions for fracture flow ---
!
          IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) 
     &      CALL SPRP_FRC_GT( NSL )
!
!---      Compute solute mole fractions for borehole flow ---
!
          IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) 
     &      CALL SPRP_BH_GT( NSL )
!
!---      Solute transport for matrix, fracture, and boreholes ---
!
          CALL TPORT_GT( NSL )
!
!---      Load old sub-time-step concentrations  ---
!
          IF( ISLC(17).NE.0 ) CALL UPDTCO( NSL )
!
!---      Load old sub-time-step concentrations for fracture flow  ---
!
          IF( ISLC(17).NE.0 .AND. (ISLC(74).EQ.1 .OR. ISLC(74).EQ.3) )
     &      CALL UPDTCO_FRC( NSL )
!
!---      Load old sub-time-step concentrations for borehole flow  ---
!
          IF( ISLC(17).NE.0 .AND. (ISLC(74).EQ.2 .OR. ISLC(74).EQ.3) )
     &      CALL UPDTCO_BH( NSL )
!
!---    Bottom of sub-time step loop  ---
!
        ENDDO
!
!---    Courant number limiting, reset time stepping  ---
!
        IF( ISLC(17).NE.0 ) THEN
          DT = DT_CRN
          DTI = DTI_CRN
          TM = TM_CRN
        ENDIF
!
!---  End of transport equation solution  ---
!
  700 CONTINUE
!
!---  Decay matrix, fracture, and borehole solutes via Bateman
!     chain decay solution  ---
!
      CALL CHAIN_DECAY

!
!---  Reactive transport  ---
!
      IF( ISLC(40).EQ.1 ) THEN
        N_CRN(NSOLU+1) = 1
        IF( ISLC(17).NE.0 ) CALL CRN_LIM( NSOLU+1 )
!
!---    Courant-limiting sub-time step loop  ---
!
        DO 792 NCR = 1,N_CRN(NSOLU+1)
          IF( ISLC(17).NE.0 ) TM = MIN( TM+DT,TM_CRN )
          DT_RST = DT
          DTI_RST = DTI
          TM_RST = TM
          TM = TM - DT
          N_RST = 1
  710     CONTINUE
          IF( N_RST.GT.16 ) THEN
            WRITE(ISC,'(A)') '          ---  ECKEChem ' // 
     &        'Sub-Time Step Reduction Limit Exceeded  ---'
            WRITE(IWR,'(A)') '          ---  ECKEChem ' // 
     &        'Sub-Time Step Reduction Limit Exceeded  ---'
            DT = DT_RST
            DTI = DTI_RST
            TM = TM_RST
            NSTEP = NSTEP-1
            TM = TM-DT
            DT = DTO
            CALL BCK_STP
            GOTO 900
          ENDIF
!
!---      ECKEChem sub-time step loop  ---
!
          DO 790 NC = 1,N_RST
            TM = TM + DT
!
!---        Loop over number of conservation component species  ---
!
            DO 730 NEQ = 1,NEQC
              NSL = NEQ + NSOLU
!
!---          Mobile conservation component fractions   ---
!
              CALL MOBCF( NEQ )
!
!---          Mobile conservation component fractions for 
!               fractures   ---
!
              IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) 
     &          CALL MOBCF_FRC( NEQ )
!
!---          Mobile conservation component fractions for 
!             boreholes   ---
!
              IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) 
     &          CALL MOBCF_BH( NEQ )
!
!---          Skip transport for immobile conservation component 
!             species   ---
!
              IF( IMMB(NEQ).EQ.1 ) GOTO 720
!
!---          Solute transport ---
!
              CALL TPORT_GT( NSL )
!
!---          Add immobile conservation component fractions   ---
!
  720         CONTINUE
              CALL IMOBCF( NEQ )
!
!---          Add immobile conservation component fractions for
!             fractures   ---
!
              IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) 
     &          CALL IMOBCF_FRC( NEQ )
!
!---          Add immobile conservation component fractions for
!             boreholes   ---
!
              IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) 
     &          CALL IMOBCF_BH( NEQ )
!
!---        End of conservation component species transport  ---
!
  730       CONTINUE
!
!---        Loop over number of kinetic component species  ---
!
            DO 750 NEQ = 1,NEQK
              NSL = NEQ + NEQC + NSOLU
! 
!---          Mobile kinetic component fractions   ---
!
              CALL MOBKF( NEQ )
!
!---          Mobile kinetic component fractions for fractures   ---
!
              IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) 
     &          CALL MOBKF_FRC( NEQ )
!
!---          Mobile kinetic component fractions for boreholes   ---
!
              IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) 
     &          CALL MOBKF_BH( NEQ )
!
!---          Skip transport for immobile kinetic component 
!             species   ---
!
              IF( IMMB(NEQ+NEQC).EQ.1 ) GOTO 740
!
!---          Solute transport ---
!
              CALL TPORT_GT( NSL )
! 
!---          Add immobile kinetic component fractions   ---
!
  740         CONTINUE
              CALL IMOBKF( NEQ )
!
!---          Add immobile kinetic component fractions for
!             fractures   ---
!
              IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) 
     &          CALL IMOBKF_FRC( NEQ )
!
!---          Add immobile kinetic component fractions for
!             boreholes   ---
!
              IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) 
     &          CALL IMOBKF_BH( NEQ )
!
!---        End of kinetic component species transport  ---
!
  750       CONTINUE
!
!---        Equilibrium-conservation-kinetic reaction chemistry   ---
!
            CALL ECKECHEM
            IF( ECKE_ER ) THEN
              CALL RESET_SP
              GOTO 710
            ENDIF
!
!---        Load old sub-time-step reactive species
!           concentrations and component species concentrations  ---
!
            IF( ISLC(17).NE.0 ) CALL UPDTCHEM
!
!---        Load old sub-time-step reactive species
!           concentrations and component species concentrations
!           for fractures  ---
!
            IF( ISLC(17).NE.0 .AND. (ISLC(74).EQ.1.OR.ISLC(74).EQ.3) )
     &        CALL UPDTCHEM_FRC
!
!---        Load old sub-time-step reactive species
!           concentrations and component species concentrations
!           for boreholes  ---
!
            IF( ISLC(17).NE.0 .AND. (ISLC(74).EQ.2.OR.ISLC(74).EQ.3) )
     &        CALL UPDTCHEM_BH
!
!---      Bottom of sub-time step loop  ---
!
  790     CONTINUE
!
!---      Reset time stepping  ---
!
          IF( N_RST.GT.1 ) THEN
            DT = DT_RST
            DTI = DTI_RST
            TM = TM_RST
          ENDIF
  792   CONTINUE
!
!---    Courant number limiting, reset time stepping  ---
!
        IF( ISLC(17).NE.0 ) THEN
          DT = DT_CRN
          DTI = DTI_CRN
          TM = TM_CRN
        ENDIF
      ENDIF







  800 CONTINUE
!
!---  Surface flux integrator  ---
!
      CALL SFIN
!
!---  PEST output for EGS Collab Experiment 1  ---
!
      IF( ISLC(80).EQ.1 ) THEN
!
!---    E1-PB Peak Tracer Concentration Arrival Time  ---
!
        IF( C_BH(6,1).GT.R_OBDT(1,1) ) THEN
          R_OBDT(1,1) = C_BH(6,1)
          R_OBDT(2,1) = TM
        ENDIF
!
!---    E1-PI Peak Tracer Concentration Arrival Time  ---
!
        IF( C_BH(21,1).GT.R_OBDT(3,1) ) THEN
          R_OBDT(3,1) = C_BH(21,1)
          R_OBDT(4,1) = TM
        ENDIF
      ENDIF
!
!---  Proceed to new time step  ---
!
      GOTO 100
!
!---  Write plot file, restart file, close files, and
!     terminate simulation  ---
!
  900 CONTINUE
      CALL WRPLOT
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL WRPLOT_FRC
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL WRPLOT_BH
      IF( ISLC(18).LT.2 ) CALL WRRST
!
!---  Inverse output  ---
!
      IF( ISLC(20).EQ.1 ) THEN
        IF( ABS(TMOB-TM).LT.(1.D+1*EPSL) ) CALL WROBDA
        IF( FLG_EXT ) WRITE(IOBDEF,'(A,/)') 'END'
        IF( FLG_UNI ) WRITE(IOBDUF,'(A,/)') 'END'
      ENDIF
!
!---  PEST output for EGS Collab Experiment 1  ---
!
      IF( ISLC(80).EQ.1 ) THEN
        CALL WROBDA_EGS1
      ENDIF
      WRITE(IWR,'(/,A)') '---  End of STOMP Simulation ---'
      WRITE(ISC,'(/,A)') '---  End of STOMP Simulation ---'

!
!---  Finalize Lis execution  ---
!
      CALL lis_finalize_f(IERR)

      STOP
!
!---  End of STOMP program  ---
!
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCF_GT
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
!     Compute boundary surface fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXS
      USE FLUXP
      USE CONST
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 BCX(LBCV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCF_GT'
!
!---  Zero boundary fluxes  ---
!
      DO NB = 1,NBC
        N = IBCN(NB)
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        IF( IBCD(NB).EQ.-3 ) THEN
          DO M = 1,ISVF
            WL(M,NPZ) = 0.D+0
            WG(M,NPZ) = 0.D+0
            WDGW(M,NPZ) = 0.D+0
            WDLA(M,NPZ) = 0.D+0
            WQ(M,NPZ) = 0.D+0
            WS(M,NPZ) = 0.D+0
            WDS(M,NPZ) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DO M = 1,ISVF
            VL(M,NPY) = 0.D+0
            VG(M,NPY) = 0.D+0
            VDGW(M,NPY) = 0.D+0
            VDLA(M,NPY) = 0.D+0
            VQ(M,NPY) = 0.D+0
            VS(M,NPY) = 0.D+0
            VDS(M,NPY) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DO M = 1,ISVF
            UL(M,NPX) = 0.D+0
            UG(M,NPX) = 0.D+0
            UDGW(M,NPX) = 0.D+0
            UDLA(M,NPX) = 0.D+0
            UQ(M,NPX) = 0.D+0
            US(M,NPX) = 0.D+0
            UDS(M,NPX) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DO M = 1,ISVF
            UL(M,NQX) = 0.D+0
            UG(M,NQX) = 0.D+0
            UDGW(M,NQX) = 0.D+0
            UDLA(M,NQX) = 0.D+0
            UQ(M,NQX) = 0.D+0
            US(M,NQX) = 0.D+0
            UDS(M,NQX) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DO M = 1,ISVF
            VL(M,NQY) = 0.D+0
            VG(M,NQY) = 0.D+0
            VDGW(M,NQY) = 0.D+0
            VDLA(M,NQY) = 0.D+0
            VQ(M,NQY) = 0.D+0
            VS(M,NQY) = 0.D+0
            VDS(M,NQY) = 0.D+0
          ENDDO
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DO M = 1,ISVF
            WL(M,NQZ) = 0.D+0
            WG(M,NQZ) = 0.D+0
            WDGW(M,NQZ) = 0.D+0
            WDLA(M,NQZ) = 0.D+0
            WQ(M,NQZ) = 0.D+0
            WS(M,NQZ) = 0.D+0
            WDS(M,NQZ) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
!
!---  Loop over boundary conditions  ---
!
      DO NB = 1,NBC
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) CYCLE
        IF( IBCM(NB).EQ.1 ) THEN
          DO N = 1,LBCV
            BCX(N) = BC(N,1,MB)
          ENDDO
        ELSE
          DO M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
             TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
             DTBC = MIN( BC(1,M,MB)-TMZ,DT )
             TFBC = (TMZ-5.D-1*DTBC-BC(1,M-1,MB))/TDBC
             DO N = 1,LBCV
               BCX(N) = BC(N,M-1,MB) + TFBC*(BC(N,M,MB)-BC(N,M-1,MB))
             ENDDO
             GOTO 105
            ENDIF
          ENDDO
          CYCLE
        ENDIF
  105   CONTINUE
        N = IBCN(NB)
        I = ID(N)
        J = JD(N)
        K = KD(N)
        NPZ = NSZ(N)
        NPY = NSY(N)
        NPX = NSX(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        ITFX = MOD(IBCT(2,NB),100)
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---      Fluid flow Neumann  ---
!
          IF( ITFX.EQ.2 ) THEN
!
!---        Aqueous saturated (BC1)  ---
!
            IF( IBCT(2,NB)/100.EQ.1 ) THEN
              WLX = BCX(3)
              WGX = 0.D+0
!
!---        Aqueous partially saturated (BC2)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.2 ) THEN
              WLX = BCX(3)*BCX(4)
              WGX = BCX(3)*(1.D+0-BCX(4))
!
!---        Aqueous unsaturated (BC3)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.3 ) THEN
              WLX = 0.D+0
              WGX = BCX(3)
            ENDIF
            DO M = 1,ISVF
              WL(M,NPZ) = WLX
              WG(M,NPZ) = WGX
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAB_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWB_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSB_GT( N,NB )
!
!---      Fluid flow aqueous Neumann  ---
!
          ELSEIF( ITFX.EQ.12 ) THEN
            DO M = 1,ISVF
              WL(M,NPZ) = BCX(3)
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAB_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGB_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWB_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSB_GT( N,NB )
!
!---      Fluid flow Seepage Face  ---
!
          ELSEIF( ITFX.EQ.13 ) THEN
            CALL DRCVLB_GT( N,NB )
!
!---      Fluid flow evaporative  ---
!
          ELSEIF( ITFX.EQ.25 ) THEN
            CALL DRCVGB_GT( N,NB )
            CALL DFFGWB_GT( N,NB )
!
!---      Fluid flow Dirichlet, Saturated, or Unit Gradient  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL DRCVLB_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFLAB_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGB_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWB_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSB_GT( N,NB )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO M = 1,ISVF
                WQ(M,NPZ) = BCX(2)
              ENDDO
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.6 ) THEN
              CALL THALB_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGB_GT( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDB_GT( N,NB )
              CALL THALB_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGB_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THDGB_GT( N,NB )
            ENDIF
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Flow flow Neumann  ---
!
          IF( ITFX.EQ.2 ) THEN
!
!---        Aqueous saturated (BC1)  ---
!
            IF( IBCT(2,NB)/100.EQ.1 ) THEN
              VLX = BCX(3)
              VGX = 0.D+0
!
!---        Aqueous partially saturated (BC2)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.2 ) THEN
              VLX = BCX(3)*BCX(4)
              VGX = BCX(3)*(1.D+0-BCX(4))
!
!---        Aqueous unsaturated (BC3)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.3 ) THEN
              VLX = 0.D+0
              VGX = BCX(3)
            ENDIF
            DO M = 1,ISVF
              VL(M,NPY) = VLX
              VG(M,NPY) = VGX
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAS_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWS_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSS_GT( N,NB )
!
!---      Fluid flow Aqueous Neumann  ---
!
          ELSEIF( ITFX.EQ.12 ) THEN
            DO M = 1,ISVF
              VL(M,NPY) = BCX(3)
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAS_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGS_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWS_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSS_GT( N,NB )
!
!---      Fluid flow Seepage Face  ---
!
          ELSEIF( ITFX.EQ.13 ) THEN
            CALL DRCVLS_GT( N,NB )
!
!---      Fluid flow evaporative  ---
!
          ELSEIF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL DRCVGS_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWS_GT( N,NB )
!
!---      Fluid flow Dirichlet, Saturated, or Unit Gradient  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL DRCVLS_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFLAS_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGS_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWS_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSS_GT( N,NB )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO M = 1,ISVF
                VQ(M,NPY) = BCX(2)
              ENDDO
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.6 ) THEN
              CALL THALS_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGS_GT( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDS_GT( N,NB )
              CALL THALS_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGS_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THDGS_GT( N,NB )
            ENDIF
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Fluid flow Neumann  ---
!
          IF( ITFX.EQ.2 ) THEN
!
!---        Aqueous saturated (BC1)  ---
!
            IF( IBCT(2,NB)/100.EQ.1 ) THEN
              ULX = BCX(3)
              UGX = 0.D+0
!
!---        Aqueous partially saturated (BC2)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.2 ) THEN
              ULX = BCX(3)*BCX(4)
              UGX = BCX(3)*(1.D+0-BCX(4))
!
!---        Aqueous unsaturated (BC3)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.3 ) THEN
              ULX = 0.D+0
              UGX = BCX(3)
            ENDIF
            DO M = 1,ISVF
              UL(M,NPX) = ULX
              UG(M,NPX) = UGX
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAW_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWW_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSW_GT( N,NB )
!
!---      Fluid flow Aqueous Neumann  ---
!
          ELSEIF( ITFX.EQ.12 ) THEN
            DO M = 1,ISVF
              UL(M,NPX) = BCX(3)
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAW_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGW_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWW_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSW_GT( N,NB )
!
!---      Fluid flow Seepage Face  ---
!
          ELSEIF( ITFX.EQ.13 ) THEN
            CALL DRCVLW_GT( N,NB )
!
!---      Fluid flow evaporative  ---
!
          ELSEIF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL DRCVGW_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWW_GT( N,NB )
!
!---      Fluid flow Dirichlet, Saturated, or Unit Gradient  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL DRCVLW_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFLAW_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGW_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWW_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSW_GT( N,NB )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO M = 1,ISVF
                UQ(M,NPX) = BCX(2)
              ENDDO
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.6 ) THEN
              CALL THALW_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGW_GT( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDW_GT( N,NB )
              CALL THALW_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGW_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THDGW_GT( N,NB )
            ENDIF
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Fluid flow Neumann  ---
!
          IF( ITFX.EQ.2 ) THEN
!
!---        Aqueous saturated (BC1)  ---
!
            IF( IBCT(2,NB)/100.EQ.1 ) THEN
              ULX = BCX(3)
              UGX = 0.D+0
!
!---        Aqueous partially saturated (BC2)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.2 ) THEN
              ULX = BCX(3)*BCX(4)
              UGX = BCX(3)*(1.D+0-BCX(4))
!
!---        Aqueous unsaturated (BC3)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.3 ) THEN
              ULX = 0.D+0
              UGX = BCX(3)
            ENDIF
            DO M = 1,ISVF
              UL(M,NQX) = ULX
              UG(M,NQX) = UGX
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAE_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWE_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSE_GT( N,NB )
!
!---      Fluid flow Aqueous Neumann  ---
!
          ELSEIF( ITFX.EQ.12 ) THEN
            DO M = 1,ISVF
              UL(M,NQX) = BCX(3)
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAE_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGE_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWE_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSE_GT( N,NB )
!
!---      Fluid flow Seepage Face  ---
!
          ELSEIF( ITFX.EQ.13 ) THEN
            CALL DRCVLE_GT( N,NB )
!
!---      Fluid flow evaporative  ---
!
          ELSEIF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL DRCVGE_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWE_GT( N,NB )
!
!---      Fluid flow Dirichlet, Saturated, or Unit Gradient  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL DRCVLE_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFLAE_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGE_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWE_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSE_GT( N,NB )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO M = 1,ISVF
                UQ(M,NQX) = BCX(2)
              ENDDO
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.6 ) THEN
              CALL THALE_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGE_GT( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDE_GT( N,NB )
              CALL THALE_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGE_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THDGE_GT( N,NB )
            ENDIF
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Fluid flow Neumann  ---
!
          IF( ITFX.EQ.2 ) THEN
!
!---        Aqueous saturated (BC1)  ---
!
            IF( IBCT(2,NB)/100.EQ.1 ) THEN
              VLX = BCX(3)
              VGX = 0.D+0
!
!---        Aqueous partially saturated (BC2)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.2 ) THEN
              VLX = BCX(3)*BCX(4)
              VGX = BCX(3)*(1.D+0-BCX(4))
!
!---        Aqueous unsaturated (BC3)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.3 ) THEN
              VLX = 0.D+0
              VGX = BCX(3)
            ENDIF
            DO M = 1,ISVF
              VL(M,NQY) = VLX
              VG(M,NQY) = VGX
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAN_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWN_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSN_GT( N,NB )
!
!---      Fluid flow Aqueous Neumann  ---
!
          ELSEIF( ITFX.EQ.12 ) THEN
            DO M = 1,ISVF
              VL(M,NQY) = BCX(3)
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAN_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGN_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWN_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSN_GT( N,NB )
!
!---      Fluid flow Seepage Face  ---
!
          ELSEIF( ITFX.EQ.13 ) THEN
            CALL DRCVLN_GT( N,NB )
!
!---      Fluid flow evaporative  ---
!
          ELSEIF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL DRCVGN_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWN_GT( N,NB )
!
!---      Fluid flow Dirichlet, Saturated, or Unit Gradient  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL DRCVLN_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFLAN_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGN_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWN_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLSN_GT( N,NB )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO M = 1,ISVF
                VQ(M,NQY) = BCX(2)
              ENDDO
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.6 ) THEN
              CALL THALN_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGN_GT( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDN_GT( N,NB )
              CALL THALN_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGN_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THDGN_GT( N,NB )
            ENDIF
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Fluid flow Neumann  ---
!
          IF( ITFX.EQ.2 ) THEN
!
!---        Aqueous saturated (BC1)  ---
!
            IF( IBCT(2,NB)/100.EQ.1 ) THEN
              WLX = BCX(3)
              WGX = 0.D+0
!
!---        Aqueous partially saturated (BC2)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.2 ) THEN
              WLX = BCX(3)*BCX(4)
              WGX = BCX(3)*(1.D+0-BCX(4))
!
!---        Aqueous unsaturated (BC3)  ---
!
            ELSEIF( IBCT(2,NB)/100.EQ.3 ) THEN
              WLX = 0.D+0
              WGX = BCX(3)
            ENDIF
            DO M = 1,ISVF
              WL(M,NQZ) = WLX
              WG(M,NQZ) = WGX
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAT_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWT_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLST_GT( N,NB )
!
!---      Fluid flow Aqueous Neumann  ---
!
          ELSEIF( ITFX.EQ.12 ) THEN
            DO M = 1,ISVF
              WL(M,NQZ) = BCX(3)
            ENDDO
            IF( ISLC(37).EQ.0 ) CALL DFFLAT_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGT_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWT_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLST_GT( N,NB )
!
!---      Fluid flow Seepage Face  ---
!
          ELSEIF( ITFX.EQ.13 ) THEN
            CALL DRCVLT_GT( N,NB )
!
!---      Fluid flow evaporative  ---
!
          ELSEIF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL DRCVGT_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWT_GT( N,NB )
!
!---      Fluid flow Dirichlet, Saturated, or Unit Gradient  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL DRCVLT_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFLAT_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DRCVGT_GT( N,NB )
            IF( ISLC(37).EQ.0 ) CALL DFFGWT_GT( N,NB )
            IF( ISLC(32).EQ.0 ) CALL DFFLST_GT( N,NB )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy Neumann  ---
!
            IF( IBCT(1,NB).EQ.2 ) THEN
              DO M = 1,ISVF
                WQ(M,NQZ) = BCX(2)
              ENDDO
!
!---        Energy outflow  ---
!
            ELSEIF( IBCT(1,NB).EQ.6 ) THEN
              CALL THALT_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGT_GT( N,NB )
!
!---        Energy not zero flux  ---
!
            ELSEIF( IBCT(1,NB).NE.3 ) THEN
              CALL THDT_GT( N,NB )
              CALL THALT_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THAGT_GT( N,NB )
              IF( ISLC(37).EQ.0 ) CALL THDGT_GT( N,NB )
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCF_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCJ_GT
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
!     Modify the Jacobian matrix for boundary conditions
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
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
      SUB_LOG(ISUB_LOG) = '/BCJ_GT'
!
!---  Loop over boundary conditions  ---
!
        DO NB = 1,NBC
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) CYCLE
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) CYCLE
        N = IBCN(NB)
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        ITFX = MOD(IBCT(2,NB),100)
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---      Fluid flow evaporative  ---
!
          IF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL JCBGWB_GT( N,NB,NPZ )
            IF( ISLC(37).EQ.0 ) CALL JCBGAB_GT( N,NB,NPZ )
!
!---      Fluid flow and salt transport  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL JCBLWB_GT( N,NB,NPZ )
            IF( ISLC(37).EQ.0 ) CALL JCBLAB_GT( N,NB,NPZ )
            IF( ISLC(37).EQ.0 ) CALL JCBGWB_GT( N,NB,NPZ )
            IF( ISLC(37).EQ.0 ) CALL JCBGAB_GT( N,NB,NPZ )
            IF( ISLC(32).EQ.0 ) CALL JCBSB_GT( N,NB,NPZ )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBTB_GT( N,NB,NPZ )
            ENDIF
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Fluid flow evaporative  ---
!
          IF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL JCBGWS_GT( N,NB,NPY )
            IF( ISLC(37).EQ.0 ) CALL JCBGAS_GT( N,NB,NPY )
!
!---      Fluid flow and salt transport  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL JCBLWS_GT( N,NB,NPY )
            IF( ISLC(37).EQ.0 ) CALL JCBLAS_GT( N,NB,NPY )
            IF( ISLC(37).EQ.0 ) CALL JCBGWS_GT( N,NB,NPY )
            IF( ISLC(37).EQ.0 ) CALL JCBGAS_GT( N,NB,NPY )
            IF( ISLC(32).EQ.0 ) CALL JCBSS_GT( N,NB,NPY )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBTS_GT( N,NB,NPY )
            ENDIF
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Fluid flow evaporative  ---
!
          IF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL JCBGWW_GT( N,NB,NPX )
            IF( ISLC(37).EQ.0 ) CALL JCBGAW_GT( N,NB,NPX )
!
!---      Fluid flow and salt transport  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL JCBLWW_GT( N,NB,NPX )
            IF( ISLC(37).EQ.0 ) CALL JCBLAW_GT( N,NB,NPX )
            IF( ISLC(37).EQ.0 ) CALL JCBGWW_GT( N,NB,NPX )
            IF( ISLC(37).EQ.0 ) CALL JCBGAW_GT( N,NB,NPX )
            IF( ISLC(32).EQ.0 ) CALL JCBSW_GT( N,NB,NPX )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBTW_GT( N,NB,NPX )
            ENDIF
          ENDIF
!
!---    East boundary
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Fluid flow evaporative  ---
!
          IF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL JCBGWE_GT( N,NB,NQX )
            IF( ISLC(37).EQ.0 ) CALL JCBGAE_GT( N,NB,NQX )
!
!---      Fluid flow and salt transport  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL JCBLWE_GT( N,NB,NQX )
            IF( ISLC(37).EQ.0 ) CALL JCBLAE_GT( N,NB,NQX )
            IF( ISLC(37).EQ.0 ) CALL JCBGWE_GT( N,NB,NQX )
            IF( ISLC(37).EQ.0 ) CALL JCBGAE_GT( N,NB,NQX )
            IF( ISLC(32).EQ.0 ) CALL JCBSE_GT( N,NB,NQX )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBTE_GT( N,NB,NQX )
            ENDIF
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Fluid flow evaporative  ---
!
          IF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL JCBGWN_GT( N,NB,NQY )
            IF( ISLC(37).EQ.0 ) CALL JCBGAN_GT( N,NB,NQY )
!
!---      Fluid flow and salt transport  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL JCBLWN_GT( N,NB,NQY )
            IF( ISLC(37).EQ.0 ) CALL JCBLAN_GT( N,NB,NQY )
            IF( ISLC(37).EQ.0 ) CALL JCBGWN_GT( N,NB,NQY )
            IF( ISLC(37).EQ.0 ) CALL JCBGAN_GT( N,NB,NQY )
            IF( ISLC(32).EQ.0 ) CALL JCBSN_GT( N,NB,NQY )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBTN_GT( N,NB,NQY )
            ENDIF
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Fluid flow evaporative  ---
!
          IF( ITFX.EQ.25 ) THEN
            IF( ISLC(37).EQ.0 ) CALL JCBGWT_GT( N,NB,NQZ )
            IF( ISLC(37).EQ.0 ) CALL JCBGAT_GT( N,NB,NQZ )
!
!---      Fluid flow and salt transport  ---
!
          ELSEIF( ITFX.NE.3 ) THEN
            CALL JCBLWT_GT( N,NB,NQZ )
            IF( ISLC(37).EQ.0 ) CALL JCBLAT_GT( N,NB,NQZ )
            IF( ISLC(37).EQ.0 ) CALL JCBGWT_GT( N,NB,NQZ )
            IF( ISLC(37).EQ.0 ) CALL JCBGAT_GT( N,NB,NQZ )
            IF( ISLC(32).EQ.0 ) CALL JCBST_GT( N,NB,NQZ )
          ENDIF
!
!---      Nonisothermal simulations  ---
!
          IF( ISLC(30).EQ.0 )  THEN
!
!---        Energy  ---
!
            IF( IBCT(1,NB).NE.3 ) THEN
              CALL JCBTT_GT( N,NB,NQZ )
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJ_GT group  ---
!      
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCP_GT
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
!     Compute saturation, relative permeability and thermodynamic
!     properties for boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
      USE GRID
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVT
      USE BCVS
      USE BCVP
      USE BCVG
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 BCX(LBCV),RKLX(3)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCP_GT'
!
!---  Loop over boundary conditions  ---
!
      DO 400 NB = 1,NBC
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 400
        ITFX = MOD(IBCT(2,NB),100)
!
!---  Assign local boundary condition variables  ---
!
        IF( IBCM(NB).EQ.1 ) THEN
          DO 80 L = 1,LBCV
            BCX(L) = BC(L,1,MB)
   80     CONTINUE
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
             TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
             DTBC = MIN( BC(1,M,MB)-TMZ,DT )
             TFBC = (TMZ-BC(1,M-1,MB))/TDBC
             DO 90 L = 1,LBCV
               BCX(L) = BC(L,M-1,MB) + TFBC*(BC(L,M,MB)-BC(L,M-1,MB))
   90        CONTINUE
!
!---         Energy boundary condition  ---
!
             IF( IBCT(1,NB).EQ.2 ) THEN
              BCX(2) = BCX(2)-5.D-1*DTBC*(BC(2,M,MB)-BC(2,M-1,MB))/TDBC
             ENDIF
!
!---         Fluid flow boundary condition  ---
!
             IF( ITFX.EQ.2 ) THEN
               BCX(3) = BCX(3)-5.D-1*DTBC*(BC(3,M,MB)-BC(3,M-1,MB))/TDBC
               BCX(4) = BCX(4)-5.D-1*DTBC*(BC(4,M,MB)-BC(4,M-1,MB))/TDBC
               BCX(5) = BCX(5)-5.D-1*DTBC*(BC(6,M,MB)-BC(6,M-1,MB))/TDBC
             ENDIF
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 400
        ENDIF
  110   CONTINUE
!
!---    Initial condition boundary condition check  ---
!
        IF( IBCT(1,NB).EQ.4 ) THEN
          IF( TB(1,NB).LT.-1.D+6 ) THEN
            INDX = 7
            IMSGX = IBCN(NB)
            NMSGX = 0
            SUBLOGX = 'BCP_GT'
            RLMSGX = 0.D+0
            CHMSGX = 'Missing Input for the Energy Initial '
     &        //'Condition Type Boundary: Check Initial Boundary '
     &        //'Condition Card: Node Number'
            CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ELSE
            BCX(2) = TB(1,NB)
          ENDIF
        ENDIF
        IF( ITFX.EQ.4 ) THEN
          IF( PGB(1,NB).LT.-1.D+6 ) THEN
            INDX = 7
            IMSGX = IBCN(NB)
            NMSGX = 0
            SUBLOGX = 'BCP_GT'
            RLMSGX = 0.D+0
            CHMSG = 'Missing Input for the Flow Initial '
     &        //'Condition Type Boundary: Check Initial Boundary '
     &        //'Condition Card: Node Number'
            CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ELSE
            IF( IBCT(2,NB).EQ.104 ) THEN
              BCX(3) = PLB(1,NB)
              BCX(4) = XLAB(1,NB)
              BCX(5) = YLSB(1,NB)
            ELSEIF( IBCT(2,NB).EQ.204 ) THEN
              BCX(3) = PGB(1,NB)
              BCX(4) = SLB(1,NB)
            ELSEIF( IBCT(2,NB).EQ.304 ) THEN
              BCX(3) = PGB(1,NB)
              BCX(4) = XGWB(1,NB)
            ENDIF
          ENDIF
        ENDIF
        N = IBCN(NB)
        N_DB = -NB
        IBD = ABS(IBCD(NB))
        IZN = IZ(N)

        POR0(1,N) = POR0(1,N)
        POR0(2,N) = POR0(2,N)

!
!---    Assign gas-entry pressure for non Brooks-Corey;
!       Brooks-Corey; Brooks-Corey, Dual Porosity; and
!       Brooks-Corey, Entrapment  ---
!
        IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &      ISCHR(IZN).EQ.202 ) THEN
          ENPR = SCHRV(1,N)*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          ENPR = MIN( SCHRV(1,N),SCHR(5,IZN) )*RHORL*GRAV
        ELSE
          ENPR = 0.D+0
        ENDIF
!
!---    Initial trapped gas saturation for the van Genuchten or
!       Brooks/Corey entrapment model  ---
!
        IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 .OR.
     &      ISCHR(IZN).EQ.201 .OR. ISCHR(IZN).EQ.202 ) THEN
!
!---      w/ Webb extension  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            ESGTMX = SCHR(15,IZN)
!
!---      w/o Webb extension  ---
!
          ELSE
            ESGTMX = SCHR(15,IZN)/(1.D+0-SCHR(4,IZN))
          ENDIF
        ELSE
          ESGTMX = 0.D+0
        ENDIF
!
!---  Boundary Direction  ---
!
        I = ID(N)
        J = JD(N)
        K = KD(N)
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        IF( IBCD(NB).EQ.-3 ) THEN
           DB = DZGP(NPZ)
           GB = GRVZ(NPZ)*DB
           GBX = GRVZ(NPZ)
           AFBX = AFZ(NPZ)
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
           DB = DYGP(NPY)*RP(ID(N))
           GB = GRVY(NPY)*DB
           GBX = GRVY(NPY)
           AFBX = AFY(NPY)
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
           DB = DXGP(NPX)
           GB = GRVX(NPX)*DB
           GBX = GRVX(NPX)
           AFBX = AFX(NPX)
        ELSEIF( IBCD(NB).EQ.1 ) THEN
           DB = -DXGP(NQX)
           GB = GRVX(NQX)*DB
           GBX = GRVX(NQX)
           AFBX = AFX(NQX)
        ELSEIF( IBCD(NB).EQ.2 ) THEN
           DB = -DYGP(NQY)*RP(ID(N))
           GB = GRVY(NQY)*DB
           GBX = GRVY(NQY)
           AFBX = AFY(NQY)
        ELSEIF( IBCD(NB).EQ.3 ) THEN
           DB = -DZGP(NQZ)
           GB = GRVZ(NQZ)*DB
           GBX = GRVZ(NQZ)
           AFBX = AFZ(NQZ)
        ENDIF
!
!---  Loop over secondary variable indices  ---
!
        DO 300 M = 2,ISVC+2
          TX = T(M,N)
          PLX = PL(M,N) + PATM
          PGX = PG(M,N) + PATM
!
!---      Energy Dirichlet or initial conditions  ---
!
          IF( IBCT(1,NB).EQ.1 .OR. IBCT(1,NB).EQ.4 ) THEN
            TX = BCX(2)
!
!---      Energy Neumann  ---
!
          ELSEIF( IBCT(1,NB).EQ.2 ) THEN
            INDX = ABS( IBCD(NB) )
            TKP = THKE_L( IZ(N),SL(M,N),THKL(M,N),PORD(M,N),
     &        PORT(M,N),INDX )
            TX = TX + BCX(2)*DB/TKP
!
!---      Geothermal gradient  ---
!
          ELSEIF( IBCT(1,NB).EQ.8 ) THEN
            TX = BCX(2) + BCX(7)*(ZPBC(NB)-BCX(6))
!
!---      Barometric formula  ---
!
          ELSEIF( IBCT(1,NB).EQ.12 ) THEN
!
!---        Temperature lapse rate, K/m  ---
!
            TLRX = BCX(7)
            TBX = BCX(2)
            TKBX = BCX(2) + TABS
            PGBX = BCX(4) + PATM
            TX = TBX + TLRX*(ZPBC(NB)-BCX(6))
!
!---        Non-zero lapse rate  ---
!
            IF( ABS(TLRX).GT.EPSL ) THEN
              PGX = PGBX*(TKBX/(TKBX + TLRX*(ZPBC(NB)-BCX(6))))**
     &          (GBX*WTMA/(RCU*TLRX))
!
!---        Zero lapse rate  ---
!
            ELSE
              PGX = PGBX*EXP(-GBX*WTMA*(ZPBC(NB)-BCX(6))/(RCU*TBX))
            ENDIF
!
!---      Energy Convective  ---
!
          ELSEIF( IBCT(1,NB).EQ.18 ) THEN
            INDX = ABS( IBCD(NB) )
            TKP = THKE_L( IZ(N),SL(M,N),THKL(M,N),PORD(M,N),
     &        PORT(M,N),INDX )
            TX = (TX*TKP/ABS(DB) + BCX(2)*BCX(6))
     &        /(TKP/ABS(DB) + BCX(6))
!
!---      Energy Convective-Radiative  ---
!
          ELSEIF( IBCT(1,NB).EQ.28 ) THEN
            INDX = ABS( IBCD(NB) )
            TKP = THKE_L( IZ(N),SL(M,N),THKL(M,N),PORD(M,N),
     &        PORT(M,N),INDX )
            TSX = TX
  150       CONTINUE
            TKS = TSX+TABS
            TKR = BCX(LBCU)+TABS
            TKC = BCX(2)+TABS
            TKX = TX+TABS
            HCX = 1.259D+0*(MAX( TKS-TKC,0.D+0 ))**(3.3333D-1)
            HRX = 5.6699D-8*(TKS**3+(TKS**2)*TKR+TKS*(TKR**2)+TKR**3)
            UX = (TKX*TKP/ABS(DB) + TKC*HCX + TKR*HRX)
            VX = (TKP/ABS(DB) + HCX + HRX)
            UXX = HCX*3.3333D-1/(TKS-TKC)
            VXX = 5.6699D-8*(3.D+0*(TKS**2)+2.D+0*TKS*TKR+(TKR**2))
            DUX = TKC*UXX + TKR*VXX
            DVX = UXX + VXX
            DF = 1.D+0 - DUX/VX + UX*DVX/(VX**2)
            F = TKS - UX/VX
            DTSX = -F/DF
            TSX = TSX + DTSX
            IF( ABS(DTSX)/TKS.GT.1.D-8 ) GOTO 150
            TX = TSX
!
!---      Energy Convective  ---
!
          ELSEIF( IBCT(1,NB).EQ.29 ) THEN
            INDX = ABS( IBCD(NB) )
            TKP = THKE_L( IZ(N),SL(M,N),THKL(M,N),PORD(M,N),
     &        PORT(M,N),INDX )
            TX = (TX*TKP/ABS(DB) + BCX(2)*BCX(LBCU))
     &        /(TKP/ABS(DB) + BCX(LBCU))
          ENDIF
!
!---      Boundary condition state #1  ---
!
!         SL = 1.0
!         SG = 0.0
!
!         Declared variables:
!
!         BCX(4) - aqueous air relative saturation, or
!                  aqueous air mass fraction
!         BCX(5) - aqueous salt relative saturation, or
!                  aqueous salt mass fraction
!
          IF( IBCT(2,NB)/100.EQ.1 ) THEN
            SGB(M,NB) = 0.D+0
            SLB(M,NB) = 1.D+0
            SNB(M,NB) = 0.D+0
            SGTX = 0.D+0
!
!---        Fluid Flow: Dirichlet, initial condition, Dirichlet-inflow,
!           Dirichlet-outflow, or hydraulic gradient  ---
!
            IF( ITFX.EQ.1 .OR. ITFX.EQ.4 .OR. ITFX.EQ.5 .OR.   
     &        ITFX.EQ.6 .OR. ITFX.EQ.8 ) THEN
              PLX = BCX(3) + PATM
              PGX = PLX
            ELSE
              PLX = MAX( PLX,PGX )
              PGX = MAX( PLX,PGX )
            ENDIF
            PX = MAX( PLX,PGX )
            PVAB(M,NB) = BCX(4)
            TMSX = BCX(5)
!
!---        Input aqueous pressure and aqueous saturation  ---
!
            ICSX = 2
            CALL FLH_IC1( BTGLB(M,NB),PGX,PLX,PVAB(M,NB),
     &        SGTX,SLB(M,NB),TX,TMSX,XMLAB(M,NB),YLSB(M,NB),
     &        IBCT(4,NB),IBCT(3,NB),ICSX,N )
            PGB(M,NB) = PGX - PATM
            PLB(M,NB) = PLX - PATM
            CALL SOL_BRNS( TX,PX,XLSMX )
            XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
            CALL SP_B( TX,XLSB(M,NB),PSBX )
            PWX = PLX
            CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL DENS_B( XLSB(M,NB),RHOBX,RHOLWX,TX )
            PVBX = PSBX
            PVAB(M,NB) = XMLAB(M,NB)*HCAW
            PVWB(M,NB) = PVBX
            XMGAB(M,NB) = PVAB(M,NB)/(PVAB(M,NB)+PVWB(M,NB))
            CALL EQUIL( XGAB(M,NB),XGWB(M,NB),XLAB(M,NB),XLSB(M,NB),
     &        XLWB(M,NB),XMGAB(M,NB),XMGWB(M,NB),XMLAB(M,NB),
     &        XMLSB(M,NB),XMLWB(M,NB) )
!
!---      Boundary condition state #2  ---
!
!         SL < 1.0
!         SG > 0.0
!
!         Declared variables:
!
!         BCX(4) - aqueous saturation
!         BCX(5) - aqueous salt relative saturation, or
!                  aqueous salt mass fraction
!
          ELSEIF( IBCT(2,NB)/100.EQ.2 ) THEN
            SLB(M,NB) = BCX(4)
!
!---        Fluid Flow: Dirichlet, initial condition, Dirichlet-inflow,
!           Dirichlet-outflow, or hydraulic gradient  ---
!
            IF( ITFX.EQ.1 .OR. ITFX.EQ.4 .OR. ITFX.EQ.5 .OR.   
     &        ITFX.EQ.6 .OR. ITFX.EQ.8 ) THEN
              PGX = BCX(3) + PATM
!
!---        Fluid Flow: Aqueous Neumann  ---
!
            ELSEIF( ITFX.EQ.12 ) THEN
              RKLX(1) = ABS(BCX(3))*1.0391029D-7/PERM(IBD,IZ(N))
              IF( RKLX(1).GT.1.D+0 ) THEN
                INDX = 14
                IMSGX = 0
                NMSGX = 0
                SUBLOGX = 'BCP_GT'
                RLMSGX = RKLX(1)
                CHMSGX = 'Aqueous Flux > Hydraulic Conductivity ' //
     &           '(Aqueous Flux/Hydraulic Conductivity) = '
                CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
              ENDIF
              ESLX = SQRT(RKLX(1))
              SLRX = SCHR(4,IZ(N))
              SLB(2,NB) = ESLX*(1.D+0-SLRX) + SLRX
!
!---        Fluid Flow: Aqueous Seepage Face  ---
!
            ELSEIF( ITFX.EQ.13 ) THEN
              ZX = BCX(8)
              KC = MAX( 1,INT(ABS(ZX-ZPBC(NB))) )
              DISTZ = (ZX-ZPBC(NB))/REAL(KC)
              PLX = BCX(3) + PATM
              DO K = 1,KC
                ZX = ZX - 5.D-1*DISTZ
                TLRX = BCX(7)
                TBX = BCX(2)
                TKBX = BCX(2) + TABS
                PGBX = BCX(4) + PATM
                TX = TBX + TLRX*(ZPBC(NB)-BCX(6))
                PGX = PGBX*(TKBX/(TKBX + TLRX*(ZPBC(NB)-BCX(6))))**
     &            (GBX*WTMA/(RCU*TLRX))
                CALL REGION_4( TX,PSWX,INDX )
                PWX = MAX( PSWX,PLX )
                CALL P_IAPWS( TX,PWX,RHOX,RHOLWX,HX,HLWX,UX,ULWX )
                CALL SOL_BRNS( TX,PSWX,XLSMX )
                IF( IBCT(4,NB).EQ.2 ) THEN
                  XLSX = BCX(5)*XLSMX
                ELSE
                  XLSX = MIN( XLSMX,BCX(5) )
                ENDIF
                CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
                PLX = PLX + RHOBX*GRAVZ*DISTZ
                PGX = MAX( PGX,PLX )
              ENDDO
              INDX = 0
              CALL REGION_4( TX,PSWX,INDX )
              CALL SCF_GL( BTGLX,PSWX )
              SGTX = 0.D+0
              INDX = 1
              CALL SP_GT( ASLFX,ASLMX,ASLX,ASLMINX,ESLX,ESGTX,
     &          ESGTMX,PGX,PLX,BTGLX,SGX,SGTX,SLB(2,NB),SLFX,
     &          SLMX,SLRX,INDX,IZN )
            ENDIF
            PX = PGX
            TMSX = BCX(5)
!
!---        Input gas pressure and aqueous saturation  ---
!
            ICSX = 1
            CALL FLH_IC2( BTGLB(M,NB),PGX,PLX,SGTX,SLB(2,NB),
     &        TX,TMSX,XMLAB(M,NB),YLSB(M,NB),IBCT(4,NB),ICSX,N )
            PGB(M,NB) = PGX - PATM
            PLB(M,NB) = PLX - PATM
            INDX = 0
            CALL REGION_4( TX,PSWX,INDX )
            CALL P_IAPWS( TX,PSWX,RHOX,RHOLWX,HX,HLWX,UX,ULWX )
            PCX = PGB(M,NB)-PLB(M,NB)
            CALL VPL_B( TX,PSWX,PCX,RHOLWX,PVWX )
            CALL SCF_GL( BTGLB(M,NB),PVWX )
            CALL SOL_BRNS( TX,PSWX,XLSMX )
            XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
            CALL DENS_B( XLSB(M,NB),RHOBX,RHOLWX,TX )
            CALL SP_B( TX,XLSB(M,NB),PSBX )
            CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX )
            CALL P_IAPWS( TX,PVBX,RHOGWX,RHOX,HGWX,HX,UGWX,UX )
            PVAB(M,NB) = PGX-PVBX
            IF( PVAB(M,NB).LT.1.D-6 ) PVAB(M,NB) = 0.D+0
            PVWB(M,NB) = PVBX
            XMLAB(M,NB) = PVAB(M,NB)/HCAW
            XMGAB(M,NB) = PVAB(M,NB)/PGX
            CALL EQUIL( XGAB(M,NB),XGWB(M,NB),XLAB(M,NB),XLSB(M,NB),
     &        XLWB(M,NB),XMGAB(M,NB),XMGWB(M,NB),XMLAB(M,NB),
     &        XMLSB(M,NB),XMLWB(M,NB) )
!
!---      Boundary condition state #3  ---
!
!         SL = 0.0
!         SG = 1.0
!
!         Declared variables:
!
!         BCX(4) - water-vapor relative saturation
!        
          ELSEIF( IBCT(2,NB)/100.EQ.3 ) THEN
!
!---        Fluid Flow: Dirichlet, initial condition, Dirichlet-inflow,
!           Dirichlet-outflow, or hydraulic gradient  ---
!
            IF( ITFX.EQ.1 .OR. ITFX.EQ.4 .OR. ITFX.EQ.5 .OR.   
     &        ITFX.EQ.6 .OR. ITFX.EQ.8 ) THEN
              PGX = BCX(3) + PATM
            ENDIF
            PVWB(M,NB) = MIN( MAX(BCX(4),0.D+0),1.D+0 )
            PX = PGX
            SLB(M,NB) = 0.D+0
            YLSB(M,NB) = 0.D+0
            TMSX = 0.D+0
!
!---        Input gas pressure and aqueous saturation  ---
!
            ICSX = 1
            CALL FLH_IC4( BTGLB(M,NB),PGX,PLX,PVAB(M,NB),PVWB(M,NB),
     &        SGTX,SLB(M,NB),TX,TMSX,XMGAB(M,NB),YLSB(M,NB),
     &        IBCT(3,NB),ICSX,N )
            PGB(M,NB) = PGX - PATM
            PLB(M,NB) = PLX - PATM
            INDX = 0
            CALL REGION_4( TX,PSWX,INDX )
            CALL SOL_BRNS( TX,PSWX,XLSMX )
            XLSB(M,NB) = 0.D+0
            CALL SCF_GL( BTGLB(M,NB),PVWB(M,NB) )
            SLB(M,NB) = 0.D+0
            SGTX = 0.D+0
            ASLMINX = -1.D+0
            CALL CAP_GT( ASLMINX,BTGLB(M,NB),PCX,SLB(M,NB),SGTX,IZN )
            PLB(M,NB) = PGB(M,NB) - PCX
            CALL P_IAPWS( TX,PVWB(M,NB),RHOGWX,RHOLWX,HGWX,HLWX,
     &        UGWX,ULWX )
            CALL DENS_B( XLSB(M,NB),RHOBX,RHOLWX,TX )
            XMLAB(M,NB) = PVAB(M,NB)/HCAW
            XMGAB(M,NB) = PVAB(M,NB)/PGX
            CALL EQUIL( XGAB(M,NB),XGWB(M,NB),XLAB(M,NB),XLSB(M,NB),
     &        XLWB(M,NB),XMGAB(M,NB),XMGWB(M,NB),XMLAB(M,NB),
     &        XMLSB(M,NB),XMLWB(M,NB) )
!
!---      Hydrostatic boundary condition  ---
!
          ELSEIF( IBCT(2,NB).EQ.11 ) THEN
            IF( IBCM(NB).EQ.1 .AND. (NSTEP-NRST).GT.1 ) GOTO 400
            IF( M.EQ.2 ) THEN
              N_DB = -NB
              CALL HYDST_BC_GT( BCX,PGX,PLX,TX,XLSB(M,NB),
     &          YLSB(M,NB),ZP(N) )
            ELSE
              TX = TB(2,NB)
              PLX = PLB(2,NB)
              PGX = PGB(2,NB)
              XLSB(M,NB) = XLSB(2,NB)
              YLSB(M,NB) = YLSB(2,NB)
            ENDIF
            PGB(M,NB) = PGX - PATM
            PLB(M,NB) = PLX - PATM
            CALL SP_B( TX,XLSB(M,NB),PSBX )
            PWX = PLX
            CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL DENS_B( XLSB(M,NB),RHOBX,RHOLWX,TX )
            PVBX = PSBX
            XMLAB(M,NB) = 0.D+0
            PVAB(M,NB) = XMLAB(M,NB)*HCAW
            PVWB(M,NB) = PVBX
            CALL SCF_GL( BTGLB(M,NB),PVWB(M,NB) )
            XMGAB(M,NB) = PVAB(M,NB)/(PVAB(M,NB)+PVWB(M,NB))
            CALL EQUIL( XGAB(M,NB),XGWB(M,NB),XLAB(M,NB),XLSB(M,NB),
     &        XLWB(M,NB),XMGAB(M,NB),XMGWB(M,NB),XMLAB(M,NB),
     &        XMLSB(M,NB),XMLWB(M,NB) )
!
!---      Evaporative boundary condition  ---
!
!         Set aqueous and gas saturation equal to the nodal value
!         compute evporation rate and apply rate as a water sink, and
!         compute associated cooling rate and apply rate as an 
!         energy sink.
!
          ELSEIF( IBCT(2,NB).EQ.25 ) THEN
            PGX = BCX(3) + PATM
            PX = PGX
            PGB(M,NB) = PGX - PATM
            PORDB(M,NB) = PORD(M,N)
            PORTB(M,NB) = PORT(M,N)
            RKLB(1,M,NB) = RKL(1,M,N)
            RKLB(2,M,NB) = RKL(2,M,N)
            RKLB(3,M,NB) = RKL(3,M,N)
            RKGB(M,NB) = RKG(M,N)
            BTGLB(M,NB) = BTGL(M,N)
            SGB(M,NB) = SG(M,N)
            SLB(M,NB) = SL(M,N)
            XLAB(M,NB) = XLA(M,N)
            XLSB(M,NB) = XLS(M,N)
            XLWB(M,NB) = XLW(M,N)
            YLSB(M,NB) = YLS(M,N)
            XMLAB(M,NB) = XMLA(M,N)
            XMLSB(M,NB) = XMLS(M,N)
            XMLWB(M,NB) = XMLW(M,N)
            SGTX = SGT(M,N)
            ASLMINX = ASLMIN(2,N)
            CALL CAP_GT( ASLMINX,BTGLB(M,NB),PCX,SLB(M,NB),SGTX,IZN )
            PLB(M,NB) = PGB(M,NB) - PCX
!
!---        Air relative humidity  ---
!
            RHX = BCX(5)
!
!---        Air specific humidity at air temperature and relative
!           humidity, kg water/kg moist air  ---
!
            INDX = 0
            TGX = BCX(2)
            CALL REGION_4( TGX,PSWX,INDX )
            PVWB(M,NB) = RHX*PSWX
            PVAB(M,NB) = PGX - PVWB(M,NB)
            IF( PVWB(M,NB).LT.1.D-6 ) PVWB(M,NB) = 0.D+0
            XMGAB(M,NB) = PVAB(M,NB)/PGX
            XMGWB(M,NB) = MAX( 1.D+0-XMGAB(M,NB),0.D+0 )
            WTMGX = XMGWB(M,NB)*WTMW + XMGAB(M,NB)*WTMA
            XGWB(M,NB) = XMGWB(M,NB)*WTMW/WTMGX
            XGAB(M,NB) = XMGAB(M,NB)*WTMA/WTMGX
            CALL P_IAPWS( TX,PVWB(M,NB),RHOGWX,RHOLWX,HGWX,HLWX,
     &        UGWX,ULWX )
            CALL DENS_B( XLSB(M,NB),RHOBX,RHOLWX,TX )
            GOTO 210
!
!---      Boundary variables equal those at the node  ---
!
          ELSE
            SGB(M,NB) = SG(M,N)
            SLB(M,NB) = SL(M,N)
            PVAB(M,NB) = PVA(M,N)
            BTGLB(M,NB) = BTGL(M,N)
            YLSB(M,NB) = YLS(M,N)
            PGB(M,NB) = PG(M,N)
            PLB(M,NB) = PL(M,N)
            XLSB(M,NB) = XLS(M,N)
            PVWB(M,NB) = PVW(M,N)
            XMGAB(M,NB) = XMGA(M,N)
            XGAB(M,NB) = XGA(M,N)
            XGWB(M,NB) = XGW(M,N)
            XLAB(M,NB) = XLA(M,N)
            XLSB(M,NB) = XLS(M,N)
            XLWB(M,NB) = XLW(M,N)
            XMGAB(M,NB) = XMGA(M,N)
            XMGWB(M,NB) = XMGW(M,N)
            XMLAB(M,NB) = XMLA(M,N)
            XMLSB(M,NB) = XMLS(M,N)
            XMLWB(M,NB) = XMLW(M,N)
            PORDB(M,NB) = PORD(M,N)
            PORTB(M,NB) = PORT(M,N)
            RKLB(1,M,NB) = RKL(1,M,N)
            RKLB(2,M,NB) = RKL(2,M,N)
            RKLB(3,M,NB) = RKL(3,M,N)
            RKGB(M,NB) = RKG(M,N)
            RHOGB(M,NB) = RHOG(M,N)
            RHOMGB(M,NB) = RHOMG(M,N)
            VISGB(M,NB) = VISG(M,N)
            RHOLB(M,NB) = RHOL(M,N)
            RHOMLB(M,NB) = RHOML(M,N)
            VISLB(M,NB) = VISL(M,N)
            DFGWB(M,NB) = DFGW(M,N)
            DFLAB(M,NB) = DFLA(M,N)
            DFLSB(M,NB) = DFLS(M,N)
            TORLB(M,NB) = TORL(M,N)
            TORGB(M,NB) = TORG(M,N)
            HGAB(M,NB) = HGA(M,N)
            HGWB(M,NB) = HGW(M,N)
            HGB(M,NB) = HG(M,N)
            THKGB(M,NB) = THKG(M,N)
            HLB(M,NB) = HL(M,N)
            THKLB(M,NB) = THKL(M,N)
            GOTO 290
          ENDIF
!
!---      Porous-media porosity  ---
!
          CALL PORSTY_GT( N,PX,PCMP(N),PORDB(M,NB),PORTB(M,NB) )
          PORDB(M,NB) = MAX( PORDB(M,NB),EPSL )
          PORTB(M,NB) = MAX( PORTB(M,NB),PORDB(M,NB) )
!
!---      Saturation, relative permeability  ---
!
          INDX = 0
          CALL KSP_GT( N,PGB(M,NB),PLB(M,NB),BTGLB(M,NB),SGB(M,NB),SGTX,
     &      SLB(M,NB),SDPF(N),SDPM(N),RKLB(1,M,NB),RKGB(M,NB),ASLX,
     &      ASLMINX,ESGTX,ESGTMX,SLRX,INDX )
  210     CONTINUE
!
!---      Gas density and component fractions  ---
!
          CALL AIRGSD( TX,PVAB(M,NB),RHOGAX )
          RHOGB(M,NB) = XGAB(M,NB)*RHOGAX + XGWB(M,NB)*RHOGWX
          WTMGX = XMGAB(M,NB)*WTMA + XMGWB(M,NB)*WTMW
          RHOMGB(M,NB) = RHOGB(M,NB)/WTMGX
!
!---      Gas viscosity  ---
!
          CALL AIRGSV( TX,VISGAX )
          CALL VISC_W( TX,PVBX,RHOGWX,VISGWX )
          CALL VISC_G( VISGAX,VISGWX,XMGAB(M,NB),XMGWB(M,NB),
     &      VISGB(M,NB) )
!
!---      Aqueous density and molar density  ---
!
          RHOLB(M,NB) = RHOBX
          WTMLX = XMLAB(M,NB)*WTMA + XMLSB(M,NB)*WTMS + XMLWB(M,NB)*WTMW
          RHOMLB(M,NB) = RHOLB(M,NB)/WTMLX
!
!---      Aqueous viscosity  ---
!
          CALL VISC_W( TX,PX,RHOLWX,VISLWX )
          CALL VISC_B( TX,XLSB(M,NB),VISLWX,VISBX )
          CALL VISC_L( XMLAB(M,NB),VISBX,VISGAX,VISLB(M,NB) )
!
!---      Gas-water diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFGWB(M,NB) = DFGWC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL BNDFAW( TX,PGX,DFGWB(M,NB) )
          ELSEIF( ISLC(2).EQ.3 ) THEN
            CALL BNDFAW( TX,PGX,DFGWB(M,NB) )
            CMFF = 1.D+0 + 2.6D+0/(DFGWC**0.5)
            AMC = PORDB(M,NB)*SLB(M,NB)
            ENHF = 9.5D+0 + 6.D+0*(AMC) -
     &        8.5D+0/EXP((CMFF*AMC)**4)
            DFGWB(M,NB) = ENHF*DFGWB(M,NB)
          ELSEIF( ISLC(2).EQ.4 ) THEN
            CALL BNDFAW( TX,PGX,DFGWB(M,NB) )
            ENHF = DFEF(1,IZ(N))+DFEF(2,IZ(N))*SLB(M,NB)-
     &        (DFEF(1,IZ(N))-DFEF(4,IZ(N)))
     &        *EXP(-((DFEF(3,IZ(N))*SLB(M,NB))**DFEF(5,IZ(N))))
            DFGWB(M,NB) = ENHF*DFGWB(M,NB)
          ENDIF
!
!---      Aqueous-air diffusion coefficients  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLAB(M,NB) = DFLAC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            CALL AIRDFL( TX,VISLB(M,NB),DFLAB(M,NB) )
          ENDIF
!
!---      Aqueous-salt diffusion coefficient  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLSB(M,NB) = DFLSC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            CALL DIFC_LS( TX,XLSB(M,NB),VISLB(M,NB),DFLSB(M,NB) )
          ENDIF
!
!---      Aqueous and gas tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SLB(M,NB),SGB(M,NB),ZERO,
     &      PORDB(M,NB),TORLB(M,NB),TORGB(M,NB),TORNX )
!
!---      Nonisothermal simulation  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Gas enthalpy and internal energy  ---
!
            CALL AIRGSH( TX,PVAB(M,NB),HGAB(M,NB),UEGAX )
            HGWB(M,NB) = HGWX
            HGB(M,NB) = XGAB(M,NB)*HGAB(M,NB) + XGWB(M,NB)*HGWB(M,NB)
!
!---        Gas thermal conductivity  ---
!
            CALL AIRGSK( TX,THKGAX )
            CALL THK_W( TX,PGX,RHOGWX,THKGWX )
            CALL THK_G( TX,THKGAX,THKGWX,XMGAB(M,NB),
     &        XMGWB(M,NB),THKGB(M,NB) )
!
!---        Aqueous enthalpy and internal energy  ---
!
            CALL ENTH_B( TX,XLSB(M,NB),HLWX,HBX )
            HLB(M,NB) = MAX(1.D+0-XLAB(M,NB),0.D+0)*HBX + 
     &        XLAB(M,NB)*HGAB(M,NB)
!
!---        Aqueous thermal conductivity  ---
!
            CALL THK_W( TX,PX,RHOLWX,THKLWX )
            CALL THK_B( TX,XLSB(M,NB),THKLWX,THKLB(M,NB) )
          ENDIF
!
!---      Evaporative boundary condition  ---
!
!         Set aqueous and gas saturation equal to the nodal value
!         compute evporation rate and apply rate as a water sink, and
!         compute associated cooling rate and apply rate as an 
!         energy sink.  ---
!
          IF( IBCT(2,NB).EQ.25 ) THEN
!
!---        Surface wet function from T.J. Lee and R.A. Pielke
!           1991. "Estimating the soil surface specific heat."
!           Journal of Applied Meteorology, 31:480-484.  ---
!
            IZN = IZ(N)
            SLRX = SCHR(4,IZN)
            IF( SLB(M,NB).LT.SLRX ) THEN
              BETAX = 2.5D-1*((1.D+0-COS(SLB(M,NB)*GPI/SLRX))**2)
            ELSE
              BETAX = 1.D+0
            ENDIF
!
!---        Exchange coefficient for mositure from J. Kondo, N. Saigusa,
!           and T. Sato. 1990. "A parameterization of evaporation
!           from bare soil surfaces." 
!           Journal of Applied Meteorology, 29:385-389.  ---
!
            CEX = 2.5D-2
!
!---        Air speed (m/s)  ---
!
            UGX = BCX(4)
!
!---        Air relative humidity  ---
!
            RHX = BCX(5)
!
!---        Air specific humidity at air temperature and relative
!           humidity, kg water/kg moist air  ---
!
            INDX = 0
            TGX = BCX(2)
            CALL REGION_4( TGX,PSWX,INDX )
            PVWX = RHX*PSWX
            PVAX = PGX - PVWX
            IF( PVAX.LT.1.D-6 ) PVAX = 0.D+0
            XMGAX = PVAX/PGX
            XMGWX = MAX( 1.D+0-XMGAX,0.D+0 )
            WTMGX = XMGWX*WTMW + XMGAX*WTMA
            XGWX = XMGWX*WTMW/WTMGX
            XGAX = XMGAX*WTMA/WTMGX
!
!---        Air density at air temperature and relative
!           humidity, kg/m^3  ---
!
            CALL P_IAPWS( TGX,PVWX,RHOGWX,RHOX,HGWX,HX,UGWX,UX )
            CALL AIRGSD( TGX,PVAX,RHOGAX )
            RHOGX = XGAX*RHOGAX + XGWX*RHOGWX
!
!---        Saturated specific humidity at surface temperature, 
!           kg water/kg moist air  ---
!
            CALL REGION_4( TX,PSWX,INDX )
            PSAX = PGX - PSWX
            IF( PSAX.LT.1.D-6 ) PSAX = 0.D+0
            XMSAX = PSAX/PGX
            XMSWX = MAX( 1.D+0-XMSAX,0.D+0 )
            WTMSX = XMSWX*WTMW + XMSAX*WTMA
            XSWX = XMSWX*WTMW/WTMSX
!
!---        Evaporation rate, kg/s m^2  ---
!
            EVAPX = RHOGX*CEX*UGX*BETAX*(XSWX-XGWX)
            SRCW(M,N) = SRCW(M,N) - EVAPX*AFBX
!
!---        Cooling rate, W  ---
!
            SRCT(M,N) = SRCT(M,N) + EVAPX*(HLWX-HGWX)*AFBX
          ENDIF
  290     CONTINUE
!
!---      Assign boundary primary variables  ---
!
          TB(M,NB) = TX
          PLB(M,NB) = PLX - PATM
          PGB(M,NB) = PGX - PATM
  300   CONTINUE
  400 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCP_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BD_GT( SLX,IZN )
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
!     Bone-dry saturation
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 12 June 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/BD_GT'
!
!---  van Genuchten saturation functions  ---
!
      IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101
     &   .OR. ISCHR(IZN).EQ.201 ) THEN
!
!---    van Genuchten saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          SLX = 0.D+0
!
!---    van Genuchten saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ------
!
        ELSE
          SLRX = SCHR(4,IZN)
          SLX = SLRX + 1.D-6
        ENDIF
!
!---  Brooks and Corey saturation functions  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102
     &   .OR. ISCHR(IZN).EQ.202 ) THEN
!
!---    Brooks and Corey saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          SLX = 0.D+0
!
!---    Brooks and Corey saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ---
!
        ELSE
          SLRX = SCHR(4,IZN)
          SLX = SLRX + 1.D-6
        ENDIF
!
!---  Dual porosity van Genuchten saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.3 ) THEN
        PORD_MX = (1.D+0-POR(4,IZN))*POR(2,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
        PORD_FX = POR(4,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
!
!---    van Genuchten saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          SLX = 0.D+0
!
!---    van Genuchten saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ------
!
        ELSE
          SLR_M = SCHR(4,IZN)
          SLR_F = SCHR(7,IZN)
          SLX = SLR_F*PORD_FX + SLR_M*PORD_MX + 1.D-6
        ENDIF
!
!---  Dual porosity Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.4 ) THEN
        PORD_MX = (1.D+0-POR(4,IZN))*POR(2,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
        PORD_FX = POR(4,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
!
!---    Brooks and Corey saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          SLX = 0.D+0
!
!---    Brooks and Corey saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ---
!
        ELSE
          SLR_M = SCHR(4,IZN)
          SLR_F = SCHR(7,IZN)
          SLX = SLR_F*PORD_FX + SLR_M*PORD_MX + 1.D-6
        ENDIF
!
!---  Haverkamp saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.5 ) THEN
        SLRX = SCHR(4,IZN)
        SLX = SLRX + 1.D-6
!
!---  Linear interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.10 ) THEN
        SLX = 1.D-6
!
!---  Log-linear interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.11 ) THEN
        SLX = 1.D-6
!
!---  Hysteretic linear interpolation function, using the drainage
!     curve as an initial guess  ---
!
      ELSEIF( ISCHR(IZN).EQ.12 ) THEN
        SLX = 1.D-6
!
!---  Hysteretic log-linear interpolation function, using the drainage
!     curve  ---
!
      ELSEIF( ISCHR(IZN).EQ.13 ) THEN
        SLX = 1.D-6
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BD_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP_GT( ASLMINX,BTGLX,CPGLX,SLX,SGTX,IZN )
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
!     saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 February 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CAP_GT'
!
!---  Set capillary pressure to oven dried conditions for
!     near-zero aqueous saturations  ---
!
      IF( SLX.LT.EPSL ) THEN
        CPGLX = SCHR(12,IZN)*RHORL*GRAV
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  van Genuchten saturation function w/o gas entrapment,
!     w/ or w/o Webb extension  ---
!
      IF( ISCHR(IZN).EQ.1 ) THEN
!
!---    van Genuchten saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        SLRX = SCHR(4,IZN)
        IF( ISM(IZN).EQ.2 ) THEN
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          CN = MAX( SCHR(3,IZN),SMALL )
          IF( SCHR(14,IZN).LE.ZERO ) THEN
            IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(14,IZN)
          ENDIF
          SMPX = SCHR(8,IZN)
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( SLX.LT.SMPX ) THEN
            HMPX = SCHR(9,IZN)
            DMPX = -(LOG10(SCHR(12,IZN))-LOG10(HMPX))/SMPX
            HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use van Genuchten function
!
          ELSE
            HDGL = (((1.D+0/ESLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
     &        SCHR(1,IZN)
          ENDIF
!
!---    van Genuchten saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ------
!
        ELSE
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          CN = MAX( SCHR(3,IZN),SMALL )
          IF( SCHR(14,IZN).LE.ZERO ) THEN
            IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(14,IZN)
          ENDIF
          IF( SLX.GT.SCHR(4,IZN) ) THEN
            HDGL = (((1.D+0/ESLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
     &        SCHR(1,IZN)
          ELSE
            HDGL = SCHR(12,IZN)
            SLX = SLRX + 1.D-6
          ENDIF
        ENDIF
        CPGLX = HDGL*RHORL*GRAV
!
!---  Brooks and Corey saturation function w/o gas entrapment,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 ) THEN
!
!---    Brooks and Corey saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        SLRX = SCHR(4,IZN)
        IF( ISM(IZN).EQ.2 ) THEN
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          CL = MAX( SCHR(3,IZN),SMALL )
          SMPX = SCHR(8,IZN)
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( SLX.LT.SMPX ) THEN
            HMPX = SCHR(9,IZN)
            DMPX = -(LOG10(SCHR(12,IZN))-LOG10(HMPX))/SMPX
            HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use Brooks and Corey function
!
          ELSE
            HDGL = SCHR(1,IZN)*(1.D+0/ESLX)**(1.D+0/CL)
          ENDIF
!
!---    Brooks and Corey saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ---
!
        ELSE
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          IF( (1.D+0-SLX)/EPSL.LT.EPSL ) THEN
            HDGL = SCHR(1,IZN)
            GOTO 222
          ENDIF
          CL = MAX( SCHR(3,IZN),SMALL )
          IF( SLX.GT.SCHR(4,IZN) ) THEN
            HDGL = SCHR(1,IZN)*(1.D+0/ESLX)**(1.D+0/CL)
          ELSE
            HDGL = SCHR(12,IZN)
            SLX = SLRX + 1.D-9
          ENDIF
  222     CONTINUE
        ENDIF
        CPGLX = HDGL*RHORL*GRAV
!
!---  Dual porosity van Genuchten saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.3 ) THEN
        CNM = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
            CMM = 1.D+0 - 2.D+0/CNM
          ELSE
            CMM = 1.D+0 - 1.D+0/CNM
          ENDIF
        ELSE
          CMM = SCHR(14,IZN)
        ENDIF
        CNF = MAX( SCHR(6,IZN),SMALL )
        IF( SCHR(15,IZN).LE.EPSL ) THEN
          IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
            CMF = 1.D+0 - 2.D+0/CNF
          ELSE
            CMF = 1.D+0 - 1.D+0/CNF
          ENDIF
        ELSE
          CMF = SCHR(15,IZN)
        ENDIF
        PORD_MX = (1.D+0-POR(4,IZN))*POR(2,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
        PORD_FX = POR(4,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
!
!---    Use matrix properties to generate a guess for 
!       capillary head  ---
!
!---    van Genuchten saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          SMPX = SCHR(8,IZN)
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( SLX.LT.SMPX ) THEN
            HMPX = SCHR(9,IZN)
            DMPX = -(LOG10(SCHR(12,IZN))-LOG10(HMPX))/SMPX
            HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use van Genuchten function
!
          ELSE
            SLRX = SCHR(4,IZN)
            ESLX = (SLX-SLRX)/(1.D+0-SLRX)
            HDGL = (((1.D+0/ESLX)**(1.D+0/CMM)-1.D+0)**(1.D+0/CNM))/
     &        SCHR(1,IZN)
          ENDIF
!
!---    van Genuchten saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ------
!
        ELSE
          IF( SLX.GT.SCHR(4,IZN) ) THEN
            SLRX = SCHR(4,IZN)
            ESLX = (SLX-SLRX)/(1.D+0-SLRX)
            HDGL = (((1.D+0/ESLX)**(1.D+0/CMM)-1.D+0)**(1.D+0/CNM))/
     &        SCHR(1,IZN)
          ELSE
            HDGL = SCHR(12,IZN)
          ENDIF
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  300   CONTINUE
        NC = NC + 1
        DO M = 1,2
          HDGLX = HDGL
          DHDGLX = MAX( 1.D-6*HDGL,1.D-6 )
          IF( M.EQ.2 ) HDGLX = HDGL + DHDGLX
!
!---      Matrix van Genuchten saturation function  
!         w/o gas entrapment, w/ Webb extension ASL = SL  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            HMPX = SCHR(9,IZN)
            SLR_M = SCHR(4,IZN)
!
!---        Capillary head above the matching point head,
!           use Webb extension  ---
!
            IF( HDGLX.GT.HMPX ) THEN
              SMPX = SCHR(8,IZN)
              HDGLZ = MIN( HDGLX,SCHR(12,IZN) )
              DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
              SLMX = -(LOG10(HDGLZ)-LOG10(SCHR(12,IZN)))*DMPX
!
!---        Capillary head at or below the matching point head,
!           use van Genuchten function
!
            ELSE
              ASL_M = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGLX)**CNM))**CMM
              SLMX = ASL_M*(1.D+0-SLR_M) + SLR_M
            ENDIF
!
!---      Matrix van Genuchten saturation function 
!         w/o gas entrapment, w/o Webb extension ASL = ESL  ------
!
          ELSE
            ASL_M = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGL)**CNM))**CMM
            SLR_M = SCHR(4,IZN)
            SLMX = ASL_M*(1.D+0-SLR_M) + SLR_M
          ENDIF
!
!---      Fracture van Genuchten saturation function  
!         w/o gas entrapment, w/ Webb extension ASL = SL  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            HMPX = SCHR(11,IZN)
            SLR_F = SCHR(7,IZN)
!
!---        Capillary head above the matching point head,
!           use Webb extension  ---
!
            IF( HDGLX.GT.HMPX ) THEN
              SMPX = SCHR(10,IZN)
              HDGLZ = MIN( HDGLX,SCHR(12,IZN) )
              DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
              SLFX = -(LOG10(HDGLZ)-LOG10(SCHR(12,IZN)))*DMPX
              ASL_F = SLFX
!
!---        Capillary head at or below the matching point head,
!           use van Genuchten function
!
            ELSE
              ASL_F = (1.D+0/(1.D+0 + (SCHR(5,IZN)*HDGLX)**CNF))**CMF
              SLFX = ASL_F*(1.D+0-SLR_F) + SLR_F
            ENDIF
!
!---      Fracture van Genuchten saturation function 
!         w/o gas entrapment, w/o Webb extension ASL = ESL  ------
!
          ELSE
            ASL_F = (1.D+0/(1.D+0 + (SCHR(5,IZN)*HDGLX)**CNF))**CMF
            SLR_F = SCHR(7,IZN)
            SLFX = ASL_F*(1.D+0-SLR_F) + SLR_F
          ENDIF
          GX(M) = SLX - SLFX*PORD_FX - SLMX*PORD_MX
        ENDDO
        F = GX(1)
        DF = (GX(2)-GX(1))/DHDGLX
        DH = -F/DF
        IF( HDGL+DH.LE.0.D+0 .AND. SLX.LT.1.D+0 ) DH = -6.D-1*HDGL
        IF( HDGL+DH.GE.SCHR(12,IZN) ) DH = 6.D-1*DH
        HDGL = HDGL + DH
        HDGL = MAX( HDGL,1.D-14 )
        HDGL = MIN( HDGL,SCHR(12,IZN) )
!
!---      No convergence on dual porosity van Genuchten 
!         capillary pressure  ---
!
        IF( NC.GT.16 ) THEN
          IF( HDGL.EQ.1.D-14 .OR. HDGL.EQ.1.D+6 ) GOTO 310
        ENDIF
        IF( NC.GT.32 ) THEN
          INDX = 14
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'CAP_GT'
          RLMSGX = SLX
          CHMSGX = 'Dual Porosity van Genuchten: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &         'Aqueous Saturation = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 300
  310   CONTINUE
        CPGLX = HDGL*RHORL*GRAV
!
!---  Dual porosity Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.4 ) THEN
        IF( (1.D+0-SLX)/EPSL.LT.EPSL ) THEN
          HDGL = SCHR(1,IZN)
          GOTO 410
        ENDIF
        CLM = MAX( SCHR(3,IZN),SMALL )
        CLF = MAX( SCHR(6,IZN),SMALL )
        PORD_MX = (1.D+0-POR(4,IZN))*POR(2,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
        PORD_FX = POR(4,IZN)/
     &    ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
!
!---    Use matrix properties to generate a guess for 
!       capillary head  ---
!
!---    Brooks and Corey saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          SMPX = SCHR(8,IZN)
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( SLX.LT.SMPX ) THEN
            HMPX = SCHR(9,IZN)
            DMPX = -(LOG10(SCHR(12,IZN))-LOG10(HMPX))/SMPX
            HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use Brooks and Corey function
!
          ELSE
            SLRX = SCHR(4,IZN)
            ESLX = (SLX-SLRX)/(1.D+0-SLRX)
            HDGL = SCHR(1,IZN)*(1.D+0/ESLX)**(1.D+0/CLM)
          ENDIF
!
!---    Brooks and Corey saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ---
!
        ELSE
          SLRX = SCHR(4,IZN)
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          IF( (1.D+0-SLX)/EPSL.LT.EPSL ) THEN
            HDGL = SCHR(1,IZN)
            GOTO 322
          ENDIF
          CL = MAX( SCHR(3,IZN),SMALL )
          IF( SLX.GT.SCHR(4,IZN) ) THEN
            HDGL = SCHR(1,IZN)*(1.D+0/ESLX)**(1.D+0/CL)
          ELSE
            HDGL = SCHR(12,IZN)
          ENDIF
  322     CONTINUE
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  400   CONTINUE
        NC = NC + 1
        DO M = 1,2
          HDGLX = HDGL
          DHDGLX = MAX( 1.D-6*HDGL,1.D-6 )
          IF( M.EQ.2 ) HDGLX = HDGL + DHDGLX
!
!---      Matrix Brooks and Corey saturation function  
!         w/o gas entrapment, w/ Webb extension ASL = SL  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            HMPX = SCHR(9,IZN)
            SLR_M = SCHR(4,IZN)
!
!---        Capillary head above the matching point head,
!           use Webb extension  ---
!
            IF( HDGLX.GT.HMPX ) THEN
              SMPX = SCHR(8,IZN)
              HDGLZ = MIN( HDGLX,SCHR(12,IZN) )
              DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
              SLMX = -(LOG10(HDGLZ)-LOG10(SCHR(12,IZN)))*DMPX
!
!---        Capillary head at or below the matching point head,
!           use Brooks and Corey function
!
            ELSE
              IF( HDGLX.LE.SCHR(1,IZN) ) THEN
                ASL_M = 1.D+0
              ELSE
                ASL_M = (SCHR(1,IZN)/HDGLX)**CLM
              ENDIF
              SLR_M = SCHR(4,IZN)
              SLMX = ASL_M*(1.D+0-SLR_M) + SLR_M
            ENDIF
!
!---      Matrix Brooks and Corey saturation function 
!         w/o gas entrapment, w/o Webb extension ASL = ESL  ------
!
          ELSE
            IF( HDGLX.LE.SCHR(1,IZN) ) THEN
              ASL_M = 1.D+0
            ELSE
              ASL_M = (SCHR(1,IZN)/HDGLX)**CLM
            ENDIF
            SLR_M = SCHR(4,IZN)
            SGR_M = SCHR(14,IZN)
            SLMX = ASL_M*(1.D+0-SLR_M-SGR_M) + SLR_M
          ENDIF
!
!---      Fracture Brooks and Corey saturation function  
!         w/o gas entrapment, w/ Webb extension ASL = SL  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            HMPX = SCHR(11,IZN)
            SLR_F = SCHR(7,IZN)
!
!---        Capillary head above the matching point head,
!           use Webb extension  ---
!
            IF( HDGLX.GT.HMPX ) THEN
              SMPX = SCHR(10,IZN)
              HDGLZ = MIN( HDGLX,SCHR(12,IZN) )
              DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
              SLFX = -(LOG10(HDGLZ)-LOG10(SCHR(12,IZN)))*DMPX
!
!---        Capillary head at or below the matching point head,
!           use Brooks and Corey function
!
            ELSE
              IF( HDGLX.LE.SCHR(5,IZN) ) THEN
                ASL_F = 1.D+0
              ELSE
                ASL_F = (SCHR(5,IZN)/HDGLX)**CLF
              ENDIF
              SLR_F = SCHR(7,IZN)
              SLFX = ASL_F*(1.D+0-SLR_F) + SLR_F
            ENDIF
!
!---      Fracture Brooks and Corey saturation function 
!         w/o gas entrapment, w/o Webb extension ASL = ESL  ------
!
          ELSE
            IF( HDGLX.LE.SCHR(5,IZN) ) THEN
              ASL_F = 1.D+0
            ELSE
              ASL_F = (SCHR(5,IZN)/HDGLX)**CLF
            ENDIF
            SLR_F = SCHR(7,IZN)
            SGR_F = SCHR(15,IZN)
            SLFX = ASL_F*(1.D+0-SLR_F-SGR_F) + SLR_F
          ENDIF
          GX(M) = SLX - SLFX*PORD_FX - SLMX*PORD_MX
        ENDDO
        F = GX(1)
        DF = (GX(2)-GX(1))/DHDGLX
        DH = -F/DF
        IF( HDGL+DH.LE.0.D+0 .AND. SLX.LT.1.D+0 ) DH = -6.D-1*HDGL
        IF( HDGL+DH.GE.SCHR(12,IZN) ) DH = 6.D-1*DH
        HDGL = HDGL + DH
        HDGL = MAX( HDGL,1.D-14 )
        HDGL = MIN( HDGL,SCHR(12,IZN) )
!
!---    No convergence Brooks-Corey capillary pressure  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'CAP_GT'
          RLMSGX = SLX
          CHMSGX = 'Dual Porosity Brooks and Corey: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &         'Aqueous Saturation = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 400
  410   CONTINUE
        CPGLX = HDGL*RHORL*GRAV
!
!---  Haverkamp saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.5 ) THEN
        IF( SLX.GT.SCHR(4,IZN) ) THEN
          ASLX = (SLX-SCHR(4,IZN))/(1.D+0-SCHR(4,IZN))
          HDGL = SCHR(1,IZN) + SCHR(5,IZN)*
     &      ((SCHR(2,IZN)/ASLX)-SCHR(2,IZN))**(1.D+0/SCHR(3,IZN))
        ELSE
          HDGL = SCHR(12,IZN)
          SLX = SCHR(4,IZN) + 1.D-9
        ENDIF
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
  500   CONTINUE
        NC = NC + 1
        SLRX = SCHR(4,IZN)
        ASLX = SCHR(2,IZN)/(SCHR(2,IZN)+((HDGL-SCHR(1,IZN))/
     &    SCHR(5,IZN))**SCHR(3,IZN))
        DASLX = -(SCHR(2,IZN)*SCHR(3,IZN)*
     &    (((HDGL-SCHR(1,IZN))/SCHR(5,IZN))**(SCHR(3,IZN)-1.D+0))
     &    /SCHR(5,IZN))/((SCHR(2,IZN)+((HDGL-SCHR(1,IZN))/SCHR(5,IZN))
     &    **SCHR(3,IZN))**2)
        SLZ = ASLX*(1.D+0-SLRX) + SLRX
        DSLZ = DASLX*(1.D+0-SLRX)
        F = SLX - SLZ
        DF = -DSLZ
        DH = -F/(DF+SMALL)
        HDGL = HDGL + DH
!
!---    No convergence Haverkamp capillary pressure  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'CAP_GT'
          RLMSGX = SLX
          CHMSGX = 'Haverkamp: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &         'Aqueous Saturation = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
        IF( ABS(DH).GT.1.D-7 ) GOTO 500
        CPGLX = HDGL*RHORL*GRAV
!
!---  Linear interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.10 ) THEN
        ITBX = 0
        HDGL = FNTBLX( SLX,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        CPGLX = HDGL*RHORL*GRAV
!
!---  Log-linear interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.11 ) THEN
        ITBX = 0
        HDGL = FNTBLX( SLX,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        HDGL = EXP(HDGL)
        CPGLX = HDGL*RHORL*GRAV
!
!---  Hysteretic linear interpolation function, using the drainage
!     curve as an initial guess  ---
!
      ELSEIF( ISCHR(IZN).EQ.12 ) THEN
        SLDX = MAX( SLX-SGTX,0.D+0 )
        ITBX = 0
        HDGL = FNTBLX( SLDX,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        CPGLX = HDGL*RHORL*GRAV
        
!
!---  Hysteretic log-linear interpolation function, using the drainage
!     curve  ---
!
      ELSEIF( ISCHR(IZN).EQ.13 ) THEN
        SLDX = MAX( SLX-SGTX,0.D+0 )
        ITBX = 0
        HDGL = FNTBLX( SLDX,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        HDGL = EXP(HDGL)
        CPGLX = HDGL*RHORL*GRAV
!
!---  van Genuchten saturation function w/ gas entrapment,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.101 ) THEN
!
!---    van Genuchten saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        SLRX = SCHR(4,IZN)
        IF( ISM(IZN).EQ.2 ) THEN
          ESGTMX = SCHR(15,IZN)
          SGFX = MAX( 1.D+0-SLX-SGTX,0.D+0 )
          SGFX = SGFX + MAX( SGTX-ESGTMX,0.D+0 )
          SGTX = MIN( SGTX,ESGTMX )
          SLX = 1.D+0 - SGFX - SGTX 
          ESGTX = SGTX
          ASLX = SLX + SGTX
          ESLX = (ASLX-SLRX)/(1.D+0-SLRX)
!
!---    van Genuchten saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          ESGTMX = SCHR(15,IZN)/(1.D+0-SLRX)
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          ESGTX = SGTX/(1.D+0-SLRX)
          ESGFX = MAX( 1.D+0-ESLX-ESGTX,0.D+0 )
          ESGFX = ESGFX + MAX( ESGTX-ESGTMX,0.D+0 )
          ESGTX = MIN( ESGTX,ESGTMX )
          ESLX = 1.D+0 - ESGFX - ESGTX 
          ASLX = ESLX + ESGTX
          ESLX = ASLX
        ENDIF
!
!---    Minimum apparent aqueous saturation  ---
!
        R = 1.D+0/ESGTMX - 1.D+0
        IF( ASLMINX.LT.0.D+0 ) THEN
          ASLMINX = (ESGTX*R*ASLX + ESGTX*(R**2)*ASLX + ASLX - ESGTX -
     &      2.D+0*ESGTX*R - ESGTX*(R**2))/
     &      (1.D+0 + ESGTX*(R**2)*ASLX - ESGTX*R - ESGTX*(R**2) )
          ASLMINX = MIN( MAX( ASLMINX,0.D+0 ),1.D+0 )
        ENDIF
!
!---    Capillary head  ---
!
        CN = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(14,IZN)
        ENDIF
        SMPX = SCHR(8,IZN)
!
!---    van Genuchten saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( SLX.LT.SMPX ) THEN
            HMPX = SCHR(9,IZN)
            DMPX = -(LOG10(SCHR(12,IZN))-LOG10(HMPX))/SMPX
            HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use van Genuchten function
!
          ELSE
            HDGL = (((1.D+0/ESLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
     &        SCHR(1,IZN)
          ENDIF
!
!---    van Genuchten saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          HDGL = (((1.D+0/ASLX)**(1.D+0/CM)-1.D+0)**(1.D+0/CN))/
     &      SCHR(1,IZN)
        ENDIF
        CPGLX = HDGL*RHORL*GRAV
!
!---  Brooks and Corey saturation function w/ gas entrapment,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.102 ) THEN
!
!---    Brooks and Corey saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        SLRX = SCHR(4,IZN)
        IF( ISM(IZN).EQ.2 ) THEN
          ESGTMX = SCHR(15,IZN)
          SGFX = MAX( 1.D+0-SLX-SGTX,0.D+0 )
          SGFX = SGFX + MAX( SGTX-ESGTMX,0.D+0 )
          SGTX = MIN( SGTX,ESGTMX )
          SLX = 1.D+0 - SGFX - SGTX 
          ESGTX = SGTX
          ASLX = SLX + SGTX
          ESLX = (ASLX-SLRX)/(1.D+0-SLRX)
!
!---    Brooks and Corey saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          ESGTMX = SCHR(15,IZN)/(1.D+0-SLRX)
          ESLX = (SLX-SLRX)/(1.D+0-SLRX)
          ESGTX = SGTX/(1.D+0-SLRX)
          ESGFX = MAX( 1.D+0-ESLX-ESGTX,0.D+0 )
          ESGFX = ESGFX + MAX( ESGTX-ESGTMX,0.D+0 )
          ESGTX = MIN( ESGTX,ESGTMX )
          ESLX = 1.D+0 - ESGFX - ESGTX 
          ASLX = ESLX + ESGTX
          ESLX = ASLX
        ENDIF
!
!---    Minimum apparent aqueous saturation  ---
!
        R = 1.D+0/ESGTMX - 1.D+0
        IF( ASLMINX.LT.0.D+0 ) THEN
          ASLMINX = (ESGTX*R*ASLX + ESGTX*(R**2)*ASLX + ASLX - ESGTX -
     &      2.D+0*ESGTX*R - ESGTX*(R**2))/
     &      (1.D+0 + ESGTX*(R**2)*ASLX - ESGTX*R - ESGTX*(R**2) )
          ASLMINX = MIN( MAX( ASLMINX,0.D+0 ),1.D+0 )
        ENDIF
!
!---    Capillary head  ---
!
        CL = MAX( SCHR(3,IZN),SMALL )
        SMPX = SCHR(8,IZN)
!
!---    Brooks and Corey saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( SLX.LT.SMPX ) THEN
            HMPX = SCHR(9,IZN)
            DMPX = -(LOG10(SCHR(12,IZN))-LOG10(HMPX))/SMPX
            HDGL = 1.D+1**(DMPX*(SLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use Brooks and Corey function
!
          ELSE
            HDGL = SCHR(1,IZN)*(1.D+0/ESLX)**(1.D+0/CL)
          ENDIF
!
!---    Brooks and Corey saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          HDGL = SCHR(1,IZN)*(1.D+0/ASLX)**(1.D+0/CL)
        ENDIF
        CPGLX = HDGL*RHORL*GRAV
!
!---  van Genuchten drainage-imbibition saturation function,
!     w/ or w/o extensions---
!
      ELSEIF( ISCHR(IZN).EQ.201 ) THEN
!
!---    van Genuchten drainage-imbibition saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        SLRX = SCHR(4,IZN)
        IF( ISM(IZN).EQ.2 ) THEN
          ESGTX = SGTX
          ASLX = SLX + SGTX
          ESLX = (ASLX-SLRX)/(1.D+0-SLRX)
          ESGTMX = SCHR(15,IZN)
!
!---    van Genuchten drainage-imbibition saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          ESGTX = SGTX/(1.D+0-SLRX)
          ASLX = (SLX+SGTX-SLRX)/(1.D+0-SLRX)
          ESLX = ASLX
          ESGTMX = SCHR(15,IZN)/(1.D+0-SLRX)
        ENDIF
!
!---    Minimum apparent aqueous saturation  ---
!
        R = 1.D+0/ESGTMX - 1.D+0
        IF( ASLMINX.LT.0.D+0 ) THEN
          ASLMINX = (ESGTX*R*ASLX + ESGTX*(R**2)*ASLX + ASLX - ESGTX -
     &      2.D+0*ESGTX*R - ESGTX*(R**2))/
     &      (1.D+0 + ESGTX*(R**2)*ASLX - ESGTX*R - ESGTX*(R**2) )
          ASLMINX = MIN( MAX( ASLMINX,0.D+0 ),1.D+0 )
        ENDIF
!
!---    Drainage capillary head  ---
!
        CND = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
            CMD = 1.D+0 - 2.D+0/CND
          ELSE
            CMD = 1.D+0 - 1.D+0/CND
          ENDIF
        ELSE
          CMD = SCHR(14,IZN)
        ENDIF
        SMPDX = SCHR(8,IZN)
!
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( SLX.LT.SMPDX ) THEN
          HMPDX = SCHR(9,IZN)
          DMPDX = -(LOG10(SCHR(12,IZN))-LOG10(HMPDX))/SMPDX
          HDGLD = 1.D+1**(DMPDX*(ASLX-SMPDX) + LOG10(HMPDX))
!
!---    Aqueous saturation at or above the matching point,
!       use van Genuchten function
!
        ELSE
          HDGLD = (((1.D+0/ESLX)**(1.D+0/CMD)-1.D+0)**(1.D+0/CND))/
     &      SCHR(1,IZN)
        ENDIF
        CNI = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
            CMI = 1.D+0 - 2.D+0/CNI
          ELSE
            CMI = 1.D+0 - 1.D+0/CNI
          ENDIF
        ELSE
          CMI = SCHR(14,IZN)
        ENDIF
        SMPIX = SCHR(8,IZN)
!
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( ASLX.LT.SMPIX ) THEN
          HMPIX = SCHR(11,IZN)
          DMPIX = -(LOG10(SCHR(12,IZN))-LOG10(HMPIX))/SMPIX
          HDGLI = 1.D+1**(DMPIX*(ASLX-SMPIX) + LOG10(HMPIX))
!
!---    Aqueous saturation at or above the matching point,
!       use van Genuchten function
!
        ELSE
          HDGLI = (((1.D+0/ESLX)**(1.D+0/CMI)-1.D+0)**(1.D+0/CNI))/
     &      SCHR(2,IZN)
        ENDIF
!
!---    Guess capillary head  ---
!
        HDGL = (ASLMINX/ASLX)*HDGLD + (1.D+0-(ASLMINX/ASLX))*HDGLI
        HDGLMN = MIN( HDGLD,HDGLI )
        HDGLMX = MAX( HDGLD,HDGLI )
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
 1100   CONTINUE
        NC = NC + 1
!
!---    No convergence on van Genuchten drainage-imbibition 
!       saturation function  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'CAP_GT'
          RLMSGX = SLX
          CHMSGX = 'van Genuchten Drainage-Imbibition: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &         'Aqueous Saturation = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
        DO 1110 M = 1,2
          HDGLZ = HDGL
          DHDGLZ = SIGN( MAX( 1.D-6*HDGL,1.D-6 ),
     &      (5.D-1*SCHR(12,IZN) - HDGL) )
          IF( M.EQ.2 ) HDGLZ = HDGL + DHDGLZ
          HMPDZ = SCHR(9,IZN)
          HMPIZ = SCHR(11,IZN)
          SLRZ = SCHR(4,IZN)
!
!---      van Genuchten drainage-imbibition saturation function,
!         w/ Webb extension ASL = SL + SGT  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
!
!---        Capillary head above the drainage matching point head,
!           use Webb extension  ---
!
            IF( HDGLZ.GT.HMPDZ ) THEN
              SMPDZ = SCHR(8,IZN)
              HDGLZ = MIN( HDGLZ,SCHR(12,IZN) )
              DMPDZ = SMPDZ/(LOG10(SCHR(12,IZN))-LOG10(HMPDZ))
              SLDZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,IZN)))*DMPDZ
              ASLDZ = SLDZ
              ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---          Capillary head above the imbibition matching point head,
!             use Webb extension  ---
!
              IF( HDGLZ.GT.HMPIZ ) THEN
                SMPIZ = SCHR(10,IZN)
                DMPIZ = SMPIZ/(LOG10(SCHR(12,IZN))-LOG10(HMPIZ))
                SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,IZN)))*DMPIZ
                ASLIZ = SLIZ
                IF( ASLDZ.GE.ASLIZ ) THEN
                  ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ELSE
                  ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ENDIF
!
!---          Capillary head at or below the imbibition matching 
!             point head, use van Genuchten function
!
              ELSE
                CNI = MAX( SCHR(5,IZN),SMALL )
                IF( SCHR(13,IZN).LE.ZERO ) THEN
                  IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                    CMI = 1.D+0 - 2.D+0/CNI
                  ELSE
                     CMI = 1.D+0 - 1.D+0/CNI
                  ENDIF
                ELSE
                  CMI = SCHR(13,IZN)
                ENDIF
                ASLIZ = (1.D+0/(1.D+0 + (SCHR(2,IZN)*HDGLZ)**CNI))**CMI
                ASLIZ = ASLIZ*(1.D+0-SLRZ) + SLRZ
                IF( ASLDZ.GE.ASLIZ ) THEN
                  ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ELSE
                  ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ENDIF
              ENDIF
!
!---        Capillary head at or below the drainage matching point head,
!           use van Genuchten function
!
            ELSE
              CND = MAX( SCHR(3,IZN),SMALL )
              IF( SCHR(14,IZN).LE.ZERO ) THEN
                IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                  CMD = 1.D+0 - 2.D+0/CND
                ELSE
                  CMD = 1.D+0 - 1.D+0/CND
                ENDIF
              ELSE
                CMD = SCHR(14,IZN)
              ENDIF
              ASLDZ = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGLZ)**CND))**CMD
              ASLDZ = ASLDZ*(1.D+0-SLRZ) + SLRZ
              ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---          Capillary head above the imbibition matching point head,
!             use Webb extension  ---
!
              IF( HDGLZ.GT.HMPIZ ) THEN
                SMPIZ = SCHR(10,IZN)
                DMPIZ = SMPIZ/(LOG10(SCHR(12,IZN))-LOG10(HMPIZ))
                SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,IZN)))*DMPIZ
                ASLIZ = SLIZ
                IF( ASLDZ.GE.ASLIZ ) THEN
                  ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ELSE
                  ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ENDIF
!
!---          Capillary head at or below the imbibition matching 
!             point head, use van Genuchten function
!
              ELSE
                CNI = MAX( SCHR(5,IZN),SMALL )
                IF( SCHR(13,IZN).LE.ZERO ) THEN
                  IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                    CMI = 1.D+0 - 2.D+0/CNI
                  ELSE
                    CMI = 1.D+0 - 1.D+0/CNI
                  ENDIF
                ELSE
                  CMI = SCHR(13,IZN)
                ENDIF
                ASLIZ = (1.D+0/(1.D+0 + (SCHR(2,IZN)*HDGLZ)**CNI))**CMI
                ASLIZ = ASLIZ*(1.D+0-SLRZ) + SLRZ
                IF( ASLDZ.GE.ASLIZ ) THEN
                  ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ELSE
                  ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ENDIF
              ENDIF
            ENDIF
!
!---        Compute trapped gas saturation, using the minimum
!           apparent aqueous saturation  ---
!
            ASLM = MAX( MIN( ASLZ,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLZ.GT.ASLM ) THEN
              SGTMZ = ESGTMX
              R = 1.D+0/SGTMZ - 1.D+0
              ESGTZ = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLZ)/(1.D+0 + R*(1.D+0-ASLZ))
              IF( ESGTZ.LT.EPSL ) ESGTZ = 0.D+0
            ELSE
              ESGTZ = 0.D+0
            ENDIF
            SLZ = ASLZ - ESGTZ
            SGTZ = ESGTZ
            SGZ = 1.D+0-SLZ
            IF( SGZ.LT.EPSL ) SGZ = 0.D+0
!
!---      van Genuchten drainage-imbibition saturation function,
!         w/o Webb extension ASL = ESL + ESGT  ---
!
          ELSE
            CND = MAX( SCHR(3,IZN),SMALL )
            IF( SCHR(14,IZN).LE.ZERO ) THEN
              IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                CMD = 1.D+0 - 2.D+0/CND
              ELSE
                CMD = 1.D+0 - 1.D+0/CND
              ENDIF
            ELSE
              CMD = SCHR(14,IZN)
            ENDIF
            ASLDZ = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGLZ)**CND))**CMD
            CNI = MAX( SCHR(5,IZN),SMALL )
            IF( SCHR(13,IZN).LE.ZERO ) THEN
              IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                CMI = 1.D+0 - 2.D+0/CNI
              ELSE
                CMI = 1.D+0 - 1.D+0/CNI
              ENDIF
            ELSE
              CMI = SCHR(13,IZN)
            ENDIF
            ASLIZ = (1.D+0/(1.D+0 + (SCHR(2,IZN)*HDGLZ)**CNI))**CMI
            ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
            IF( ASLDZ.GE.ASLIZ ) THEN
              ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &           + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
            ELSE
              ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
            ENDIF
!
!---        Compute trapped gas saturation, using the minimum
!           apparent aqueous saturation  ---
!
            ASLM = MAX( MIN( ASLZ,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLZ.GT.ASLM ) THEN
              R = 1.D+0/ESGTMX - 1.D+0
              ESGTZ = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLZ)/(1.D+0 + R*(1.D+0-ASLZ))
              IF( ESGTZ.LT.EPSL ) ESGTZ = 0.D+0
            ELSE
              ESGTZ = 0.D+0
            ENDIF
            ESLZ = ASLZ - ESGTZ
            SLZ = ESLZ*(1.D+0-SLRZ) + SLRZ
            SGTZ = ESGTZ*(1.D+0-SLRZ)
            SGZ = 1.D+0-SLZ
            IF( SGZ.LT.EPSL ) SGZ = 0.D+0
          ENDIF
          GX(M) = SLX - SLZ
 1110   CONTINUE
!
!---    Convergence check  ---
!
        IF( ABS(GX(1)).LT.1.D-9 ) GOTO 1120
!
!---    Solve linear system  ---
!
        FX = GX(1)
        DFX = (GX(2)-GX(1))/DHDGLZ
        DHDGL = -FX/DFX
!
!---    Update primary unknowns  ---
!
        DHX = MAX( 1.D-3,2.D-1*HDGL )
        DHDGL = SIGN( MIN( DHX,ABS(DHDGL) ),DHDGL )
        HDGL = HDGL + DHDGL
        IF( (HDGL-HDGLMN).LT.EPSL .AND. NC.GE.3 ) GOTO 1120
        IF( (HDGLMX-HDGL).LT.EPSL .AND. NC.GE.3 ) GOTO 1120
        HDGL = MIN( HDGL,HDGLMX )
        HDGL = MAX( HDGL,HDGLMN )
        GOTO 1100
!
!---    Converged solution return capillary pressure  ---
!
 1120   CONTINUE
        CPGLX = HDGL*RHORL*GRAV
!
!---  Brooks and Corey drainage-imbibition saturation function,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.202 ) THEN
!
!---    Brooks and Corey drainage-imbibition saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        SLRX = SCHR(4,IZN)
        IF( ISM(IZN).EQ.2 ) THEN
          ESGTX = SGTX
          ASLX = SLX + SGTX
          ESGTMX = SCHR(15,IZN)
!
!---    Brooks and Corey drainage-imbibition saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          ESGTX = SGTX/(1.D+0-SLRX)
          ASLX = (SLX+SGTX-SLRX)/(1.D+0-SLRX)
          ESGTMX = SCHR(15,IZN)/(1.D+0-SLRX)
        ENDIF
!
!---    Minimum apparent aqueous saturation  ---
!
        R = 1.D+0/ESGTMX - 1.D+0
        IF( ASLMINX.LT.0.D+0 ) THEN
          ASLMINX = (ESGTX*R*ASLX + ESGTX*(R**2)*ASLX + ASLX - ESGTX -
     &      2.D+0*ESGTX*R - ESGTX*(R**2))/
     &      (1.D+0 + ESGTX*(R**2)*ASLX - ESGTX*R - ESGTX*(R**2) )
          ASLMINX = MIN( MAX( ASLMINX,0.D+0 ),1.D+0 )
        ENDIF
!
!---    Drainage capillary head  ---
!
        CLD = MAX( SCHR(3,IZN),SMALL )
        SMPDX = SCHR(8,IZN)
!
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( ASLX.LT.SMPDX ) THEN
          HMPDX = SCHR(9,IZN)
          DMPDX = -(LOG10(SCHR(12,IZN))-LOG10(HMPDX))/SMPDX
          HDGLD = 1.D+1**(DMPDX*(ASLX-SMPDX) + LOG10(HMPDX))
!
!---    Aqueous saturation at or above the matching point,
!       use Brooks and Corey function
!
        ELSE
!
!---      w/ Webb extension  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            ASLZ = (ASLX-SLRX)/(1.D+0-SLRX)
            HDGLD = SCHR(1,IZN)*(1.D+0/ASLZ)**(1.D+0/CLD)
!
!---      w/o Webb extension  ---
!
          ELSE
            HDGLD = SCHR(1,IZN)*(1.D+0/ASLX)**(1.D+0/CLD)
          ENDIF
        ENDIF
!
!---    Imbibition capillary head  ---
!
        CLI = MAX( SCHR(6,IZN),SMALL )
        SMPIX = SCHR(10,IZN)
!
!---    Aqueous saturation below the matching point,
!       use Webb extension  ---
!
        IF( ASLX.LT.SMPIX ) THEN
          HMPIX = SCHR(11,IZN)
          DMPIX = -(LOG10(SCHR(12,IZN))-LOG10(HMPIX))/SMPIX
          HDGLI = 1.D+1**(DMPIX*(ASLX-SMPIX) + LOG10(HMPIX))
!
!---    Aqueous saturation at or above the matching point,
!       use Brooks and Corey function
!
        ELSE
!
!---      w/ Webb extension  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            ASLZ = (ASLX-SLRX)/(1.D+0-SLRX)
            HDGLI = SCHR(5,IZN)*(1.D+0/ASLZ)**(1.D+0/CLI)
!
!---      w/o Webb extension  ---
!
          ELSE
            HDGLI = SCHR(5,IZN)*(1.D+0/ASLX)**(1.D+0/CLI)
          ENDIF         
        ENDIF
!
!---    Guess capillary head  ---
!
        HDGL = (ASLMINX/ASLX)*HDGLD + (1.D+0-(ASLMINX/ASLX))*HDGLI
        HDGLMN = MIN( HDGLD,HDGLI )
        HDGLMX = MAX( HDGLD,HDGLI )
!
!---    Start Newton-Raphson solution  ---
!
        NC = 0
 1200   CONTINUE
        NC = NC + 1
!
!---    No convergence on Brooks and Corey drainage-imbibition 
!       saturation function  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 14
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'CAP_GT'
          RLMSGX = SLX
          CHMSGX = 'Brooks and Corey Drainage-Imbibition: '
     &      // 'No Convergence on Capillary Pressure: ' //
     &         'Aqueous Saturation = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
        DO 1210 M = 1,2
          HDGLZ = HDGL
          DHDGLZ = SIGN( MAX( 1.D-6*HDGL,1.D-6 ),
     &      (5.D-1*SCHR(12,IZN) - HDGL) )
          IF( M.EQ.2 ) HDGLZ = HDGL + DHDGLZ
          HMPDZ = SCHR(9,IZN)
          HMPIZ = SCHR(11,IZN)
          SLRZ = SCHR(4,IZN)
!
!---      Brooks and Corey drainage-imbibition saturation function,
!         w/ Webb extension ASL = SL + SGT  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
!
!---        Capillary head above the drainage matching point head,
!           use Webb extension  ---
!
            IF( HDGLZ.GT.HMPDZ ) THEN
              SMPDZ = SCHR(8,IZN)
              HDGLZ = MIN( HDGLZ,SCHR(12,IZN) )
              DMPDZ = SMPDZ/(LOG10(SCHR(12,IZN))-LOG10(HMPDZ))
              SLDZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,IZN)))*DMPDZ
              ASLDZ = SLDZ
              ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---          Capillary head above the imbibition matching point head,
!             use Webb extension  ---
!
              IF( HDGLZ.GT.HMPIZ ) THEN
                SMPIZ = SCHR(10,IZN)
                DMPIZ = SMPIZ/(LOG10(SCHR(12,IZN))-LOG10(HMPIZ))
                SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,IZN)))*DMPIZ
                ASLIZ = SLIZ
                IF( ASLDZ.GE.ASLIZ ) THEN
                  ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ELSE
                  ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ENDIF
!
!---          Capillary head at or below the imbibition matching 
!             point head, use Brooks and Corey function
!
              ELSE
                CLI = MAX( SCHR(6,IZN),SMALL )
                IF( HDGLZ.LE.SCHR(5,IZN) ) THEN
                  ASLIZ = 1.D+0
                ELSE
                  ASLIZ = (SCHR(5,IZN)/HDGLZ)**CLI
                ENDIF
                ASLIZ = ASLIZ*(1.D+0-SLRZ) + SLRZ
                IF( ASLDZ.GE.ASLIZ ) THEN
                  ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ELSE
                  ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ENDIF
              ENDIF
!
!---        Capillary head at or below the drainage matching point head,
!           use Brooks and Corey function
!
            ELSE
              CLD = MAX( SCHR(3,IZN),SMALL )
              IF( HDGLZ.LE.SCHR(1,IZN) ) THEN
                ASLDZ = 1.D+0
              ELSE
                ASLDZ = (SCHR(1,IZN)/HDGLZ)**CLD
              ENDIF
              ASLDZ = ASLDZ*(1.D+0-SLRZ) + SLRZ
              ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
!
!---          Capillary head above the imbibition matching point head,
!             use Webb extension  ---
!
              IF( HDGLZ.GT.HMPIZ ) THEN
                SMPIZ = SCHR(10,IZN)
                DMPIZ = SMPIZ/(LOG10(SCHR(12,IZN))-LOG10(HMPIZ))
                SLIZ = -(LOG10(HDGLZ)-LOG10(SCHR(12,IZN)))*DMPIZ
                ASLIZ = SLIZ
                IF( ASLDZ.GE.ASLIZ ) THEN
                  ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ELSE
                  ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ENDIF
!
!---          Capillary head at or below the imbibition matching 
!             point head, use Brooks and Corey function
!
              ELSE
                CLI = MAX( SCHR(6,IZN),SMALL )
                IF( HDGLZ.LE.SCHR(5,IZN) ) THEN
                  ASLIZ = 1.D+0
                ELSE
                  ASLIZ = (SCHR(5,IZN)/HDGLZ)**CLI
                ENDIF
                ASLIZ = ASLIZ*(1.D+0-SLRZ) + SLRZ
                IF( ASLDZ.GE.ASLIZ ) THEN
                  ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ELSE
                  ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
                ENDIF
              ENDIF
            ENDIF
!
!---        Compute trapped gas saturation, using the minimum
!           apparent aqueous saturation  ---
!
            ASLM = MAX( MIN( ASLZ,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLZ.GT.ASLM ) THEN
              SGTMZ = ESGTMX
              R = 1.D+0/SGTMZ - 1.D+0
              ESGTZ = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLZ)/(1.D+0 + R*(1.D+0-ASLZ))
              IF( ESGTZ.LT.EPSL ) ESGTZ = 0.D+0
            ELSE
              ESGTZ = 0.D+0
            ENDIF
            SLZ = ASLZ - ESGTZ
            SGTZ = ESGTZ
            SGZ = 1.D+0-SLZ
            IF( SGZ.LT.EPSL ) SGZ = 0.D+0
!
!---      Brooks and Corey drainage-imbibition saturation function,
!         w/o Webb extension ASL = ESL + ESGT  ---
!
          ELSE
            CLD = MAX( SCHR(3,IZN),SMALL )
            IF( HDGLZ.LE.SCHR(1,IZN) ) THEN
              ASLDZ = 1.D+0
            ELSE
              ASLDZ = (SCHR(1,IZN)/HDGLZ)**CLD
            ENDIF
            CLI = MAX( SCHR(6,IZN),SMALL )
            IF( HDGLZ.LE.SCHR(5,IZN) ) THEN
              ASLIZ = 1.D+0
            ELSE
              ASLIZ = (SCHR(5,IZN)/HDGLZ)**CLI
            ENDIF
            ASLMZ = MAX( MIN( ASLDZ,ASLMINX ),0.D+0 )
            IF( ASLDZ.GE.ASLIZ ) THEN
              ASLZ = 5.D-1*(ASLIZ + SQRT((ASLIZ**2)
     &           + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
            ELSE
              ASLZ = 5.D-1*(ASLIZ - SQRT((ASLIZ**2)
     &               + 4.D+0*ASLMZ*ASLDZ - 4.D+0*ASLIZ*ASLMZ))
            ENDIF
!
!---        Compute trapped gas saturation, using the minimum
!           apparent aqueous saturation  ---
!
            ASLM = MAX( MIN( ASLZ,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLZ.GT.ASLM ) THEN
              R = 1.D+0/ESGTMX - 1.D+0
              ESGTZ = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLZ)/(1.D+0 + R*(1.D+0-ASLZ))
              IF( ESGTZ.LT.EPSL ) ESGTZ = 0.D+0
            ELSE
              ESGTZ = 0.D+0
            ENDIF
            ESLZ = ASLZ - ESGTZ
            SLZ = ESLZ*(1.D+0-SLRZ) + SLRZ
            SGTZ = ESGTZ*(1.D+0-SLRZ)
            SGZ = 1.D+0-SLZ
            IF( SGZ.LT.EPSL ) SGZ = 0.D+0
          ENDIF
          GX(M) = SLX - SLZ
 1210   CONTINUE
!
!---    Convergence check  ---
!
        IF( ABS(GX(1)).LT.1.D-9 ) GOTO 1220
!
!---    Solve linear system  ---
!
        FX = GX(1)
        DFX = (GX(2)-GX(1))/DHDGLZ
        DHDGL = -FX/DFX
!
!---    Update primary unknowns  ---
!
        DHX = MAX( 1.D-3,2.D-1*HDGL )
        DHDGL = SIGN( MIN( DHX,ABS(DHDGL) ),DHDGL )
        IF( (HDGL-HDGLMN).LT.EPSL .AND. NC.GE.3 ) GOTO 1220
        IF( (HDGLMX-HDGL).LT.EPSL .AND. NC.GE.3 ) GOTO 1220
        HDGL = HDGL + DHDGL
        HDGL = MIN( HDGL,HDGLMX )
        HDGL = MAX( HDGL,HDGLMN )
        GOTO 1200
!
!---    Converged solution return capillary pressure  ---
!
 1220   CONTINUE
        CPGLX = HDGL*RHORL*GRAV
      ENDIF
      CPGLX = CPGLX*BTGLX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_BC_GT
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
!     Boundary condition check.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, November 14, 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      SUB_LOG(ISUB_LOG) = '/CHK_BC_GT'
!
!---  Loop over boundary conditions  ---
!
      DO NB = 1,NBC
        N = IBCN(NB)
!
!---    Initial condition boundary condition check  ---
!
        IF( IBCT(1,NB).EQ.4 ) THEN
          IF( TB(1,NB).LT.-1.D+6 ) THEN
            INDX = 30
            IMSG = N
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        IF( ITFX.EQ.4 ) THEN
          IF( PGB(1,NB).LT.-1.D+6 ) THEN
            INDX = 30
            IMSG = N
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_BC_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_GT
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
!     user inputs.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE PLT_ATM
      USE JACOB
      USE HYST
      USE GRID
      USE FDVS
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
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHK_GT'
      EPSLX = 1.D-4
!
!---  Scaling factors  ---
!
      IF( ISLC(19).EQ.1 ) THEN
!
!---    Simple scaling  ---
!
        DO IZN = 1,NROCK
          PERM(1,IZN) = SCALING( GAMMA(1,IZN),PERM(1,IZN),IGAMMA(1) )
          PERM(2,IZN) = SCALING( GAMMA(1,IZN),PERM(2,IZN),IGAMMA(1) )
          PERM(3,IZN) = SCALING( GAMMA(1,IZN),PERM(3,IZN),IGAMMA(1) )
          PERM(6,IZN) = SCALING( GAMMA(1,IZN),PERM(6,IZN),IGAMMA(1) )
          PERM(7,IZN) = SCALING( GAMMA(6,IZN),PERM(7,IZN),IGAMMA(1) )
          PERM(8,IZN) = SCALING( GAMMA(6,IZN),PERM(8,IZN),IGAMMA(1) )
          PERM(9,IZN) = SCALING( GAMMA(6,IZN),PERM(9,IZN),IGAMMA(1) )
          POR(1,IZN) = SCALING( GAMMA(2,IZN),POR(1,IZN),IGAMMA(2) )
          POR(2,IZN) = SCALING( GAMMA(2,IZN),POR(2,IZN),IGAMMA(2) )
          POR(3,IZN) = SCALING( GAMMA(7,IZN),POR(3,IZN),IGAMMA(2) )
          POR(4,IZN) = SCALING( GAMMA(7,IZN),POR(4,IZN),IGAMMA(2) )
          SCHR(1,IZN) = SCALING( GAMMA(3,IZN),SCHR(1,IZN),IGAMMA(3) )
          SCHR(3,IZN) = SCALING( GAMMA(4,IZN),SCHR(3,IZN),IGAMMA(4) )
!
!---      Scaling for defaulted Mualem/Burdine parameters  ---
!
          IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
            IF( SCHR(14,IZN)/EPSL.LE.EPSL ) THEN
              IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                SCHR(11,IZN) = 1.D+0 - 2.D+0/(SCHR(3,IZN)+SMALL)
              ELSE
                SCHR(11,IZN) = 1.D+0 - 1.D+0/(SCHR(3,IZN)+SMALL)
              ENDIF
            ENDIF
          ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
            IF( SCHR(14,IZN)/EPSL.LE.EPSL ) THEN
              IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                SCHR(11,IZN) = SCHR(3,IZN)
              ELSE
                SCHR(11,IZN) = SCHR(3,IZN)
              ENDIF
            ENDIF
          ENDIF
!
!---      Scaling for Mualem-Anisotropy parameters  ---
!
          IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
            IF( IRPL(IZN).EQ.301 ) THEN
              SCHR(11,IZN) = SCALING( GAMMA(11,IZN),SCHR(11,IZN),
     &          IGAMMA(4) )
              SCHR(12,IZN) = SCALING( GAMMA(12,IZN),SCHR(12,IZN),
     &          IGAMMA(4) )
              SCHR(13,IZN) = SCALING( GAMMA(13,IZN),SCHR(13,IZN),
     &          IGAMMA(4) )
            ENDIF
          ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
            IF( IRPL(IZN).EQ.301 ) THEN
              SCHR(11,IZN) = SCALING( GAMMA(11,IZN),SCHR(11,IZN),
     &          IGAMMA(4) )
              SCHR(12,IZN) = SCALING( GAMMA(12,IZN),SCHR(12,IZN),
     &          IGAMMA(4) )
              SCHR(13,IZN) = SCALING( GAMMA(13,IZN),SCHR(13,IZN),
     &          IGAMMA(4) )
            ENDIF
          ENDIF
          SCHR(4,IZN) = SCALING( GAMMA(5,IZN),SCHR(4,IZN),IGAMMA(5) )
          SCHR(5,IZN) = SCALING( GAMMA(8,IZN),SCHR(5,IZN),IGAMMA(3) )
          SCHR(6,IZN) = SCALING( GAMMA(9,IZN),SCHR(6,IZN),IGAMMA(4) )
          SCHR(7,IZN) = SCALING( GAMMA(10,IZN),SCHR(7,IZN),IGAMMA(5) )
          PERM(4,IZN) = PERM(1,IZN)
          PERM(1,IZN) = PERM(4,IZN)*(1.D+0-POR(4,IZN)) +
     &      PERM(7,IZN)*POR(4,IZN)
          PERM(5,IZN) = PERM(2,IZN)
          PERM(2,IZN) = PERM(5,IZN)*(1.D+0-POR(4,IZN)) +
     &      PERM(8,IZN)*POR(4,IZN)
          PERM(6,IZN) = PERM(3,IZN)
          PERM(3,IZN) = PERM(6,IZN)*(1.D+0-POR(4,IZN)) +
     &      PERM(9,IZN)*POR(4,IZN)
        ENDDO
      ENDIF
!
!---  Check initial temperature, aqueous pressure, gas pressure,
!     and aqueous saturation  ---
!
      INDX = 0
      DO N = 1,NFLD
        DO M = 1,ISVC+2
          T(M,N) = T(2,N)
          PG(M,N) = PG(2,N)
          PL(M,N) = PL(2,N)
          SG(M,N) = SG(2,N)
          SL(M,N) = SL(2,N)
        ENDDO
        IZN = IZ(N)
        N_DB = NB
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
        IF( T(2,N).GT.2000.D+0 .OR. T(2,N).LT.0.01D+0 ) THEN
          INDX = 17
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'CHK_GT'
          RLMSGX = T(2,N)
          CHMSGX = 'Out of Range Initial Temperature(C) = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
        IF( PL(2,N).GT.100.D+6-PATM ) THEN
          INDX = 17
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'CHK_GT'
          RLMSGX = PL(2,N)+PATM
          CHMSGX = 'Out of Range Initial Aqueous Pressure(Pa) = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
        IF( PG(2,N).GT.100.D+6-PATM .OR. 
     &    PG(2,N).LT.6.1125D+2-PATM ) THEN
          INDX = 17
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'CHK_GT'
          RLMSGX = PG(2,N)+PATM
          CHMSGX = 'Out of Range Initial Gas Pressure(Pa) = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
        IF( SL(2,N).GT.1.D+0 .OR. SL(2,N).LT.0.D+0 ) THEN
          INDX = 17
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'CHK_GT'
          RLMSGX = SL(2,N)
          CHMSGX = 'Out of Range Initial Aqueous Saturation = '
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
        ENDIF
!
!---    Check for conflicts in saturation functions with
!       gas entrapment and gas relative permeability functions
!       with residual gas saturaitons  ---
!
        IF( ISCHR(IZN).EQ.12 .OR. ISCHR(IZN).EQ.13 .OR.
     &    ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.201 .OR. ISCHR(IZN).EQ.202 ) THEN
!
!---      Modified-Corey gas relative permeability function  ---
!
          IF( IRPG(IZN).EQ.3 ) THEN
            IF( RPGC(3,IZN).GT.EPSL ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_GT'
              RLMSGX = RPGC(3,IZN)
              CHMSGX = 'Modified-Corey Gas Relative Permeability: ' // 
     &          'Non-zero Residual Gas Saturation with '  //
     &          'Gas Entrapment in Saturation Function'
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
!
!---      Free-Corey gas relative permeability function  ---
!
          ELSEIF( IRPG(IZN).EQ.7 ) THEN
            IF( RPGC(4,IZN).GT.EPSL ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_GT'
              RLMSGX = RPGC(4,IZN)
              CHMSGX = 'Free-Corey Gas Relative Permeability: ' // 
     &          'Non-zero Residual Gas Saturation with '  //
     &          'Gas Entrapment in Saturation Function'
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF         
!
!---       van Genuchten gas relative permeability function  ---
!
          ELSEIF( IRPG(IZN).EQ.9 ) THEN
            IF( RPGC(3,IZN).GT.EPSL ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_GT'
              RLMSGX = RPGC(4,IZN)
              CHMSGX = 'van Genuchten Gas Relative Permeability: ' // 
     &          'Non-zero Residual Gas Saturation with '  //
     &          'Gas Entrapment in Saturation Function'
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
!
!---      Classical-Corey gas relative permeability function  ---
!
          ELSEIF( IRPG(IZN).EQ.17 ) THEN
            IF( RPGC(4,IZN).GT.EPSL ) THEN
              INDX = 17
              IMSGX = 0
              NMSGX = 0
              SUBLOGX = 'CHK_GT'
              RLMSGX = RPGC(3,IZN)
              CHMSGX = 'Classical-Corey Gas Relative Permeability: ' // 
     &          'Non-zero Residual Gas Saturation with '  //
     &          'Gas Entrapment in Saturation Function'
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
          ENDIF
        ENDIF
!
!---    Webb saturation and capillary pressure matching points  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
!
!---      van Genuchten moisture retension function  ---
!
          IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 .OR. 
     &      ISCHR(IZN).EQ.201 ) THEN
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
              HDOD = SCHR(12,IZN)
              CALL WEBB_VG( ALPHAX,CMX,CNX,HMPX,SMPX,SRX )
              SCHR(10,IZN) = SMPX
              SCHR(11,IZN) = HMPX            
            ENDIF
!
!---      Brooks and Corey moisture retension function  ---
!
          ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR. 
     &      ISCHR(IZN).EQ.202 ) THEN
            SRX = SCHR(4,IZN)
            PSIX = SCHR(1,IZN)
            CLX = SCHR(3,IZN)
            HDOD = SCHR(12,IZN)         
            CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
            SCHR(8,IZN) = SMPX
            SCHR(9,IZN) = HMPX
            IF( ISCHR(IZN).EQ.202 ) THEN
              PSIX = SCHR(5,IZN)
              CLX = SCHR(6,IZN)
              HDOD = SCHR(12,IZN)         
              CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
              SCHR(10,IZN) = SMPX
              SCHR(11,IZN) = HMPX
            ENDIF
!
!---      Dual porosity van Genuchten moisture retension function  ---
!
          ELSEIF( ISCHR(IZN).EQ.3 ) THEN
!
!---        Matrix matching point  ---
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
            HDOD = SCHR(12,IZN)
            CALL WEBB_VG( ALPHAX,CMX,CNX,HMPX,SMPX,SRX )
            SCHR(8,IZN) = SMPX
            SCHR(9,IZN) = HMPX            
!
!---        Frature matching point  ---
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
            HDOD = SCHR(12,IZN)
            CALL WEBB_VG( ALPHAX,CMX,CNX,HMPX,SMPX,SRX )
            SCHR(10,IZN) = SMPX
            SCHR(11,IZN) = HMPX            
!
!---      Dual porosity Brooks and Corey moisture retension 
!         function  ---
!
          ELSEIF( ISCHR(IZN).EQ.4 ) THEN
!
!---        Matrix matching point  ---
!
            SRX = SCHR(4,IZN)
            PSIX = SCHR(1,IZN)
            CLX = SCHR(3,IZN)
            HDOD = SCHR(12,IZN)        
            CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
            SCHR(8,IZN) = SMPX
            SCHR(9,IZN) = HMPX
!
!---        Matrix matching point  ---
!
            SRX = SCHR(7,IZN)
            PSIX = SCHR(5,IZN)
            CLX = SCHR(6,IZN)
            HDOD = SCHR(12,IZN)         
            CALL WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
            SCHR(10,IZN) = SMPX
            SCHR(11,IZN) = HMPX
          ENDIF
        ENDIF

!
!---    Load reactive transport total and diffusive porosity  ---
!
        POR0(1,N) = POR(1,IZN)
        POR0(2,N) = POR(2,IZN)

      ENDDO
      IF( INDX.GT.0 ) STOP
!
!---    Establish reference pressure for soil compressibility  ---
!
      DO N = 1,NFLD
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
        IF( IEO.EQ.2 .AND. NPHAZ(2,N).NE.0 ) CYCLE
        TCMP(N) = T(2,N)
        IZN = IZ(N)
        IF( CMP(3,IZN).GT.PATM ) THEN
          PCMP(N) = CMP(3,IZN)
        ELSEIF( ISLC(61).EQ.0 ) THEN
          PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
        ENDIF
      ENDDO
!
!---  Hydrostatic with geothermal gradient initial conditions, either
!     zonal or entire domain, with salt mass gradient  ---
!
      IF( MOD(ISIC,10).EQ.4 ) THEN
        NHORZ = (ISIC-10)/10 + 1
!
!---    Loop over all nodes  ---
!
        DO N = 1,NFLD
!
!---      Skip inactive nodes  ---
!
          IF( IXP(N).EQ.0 ) CYCLE
          IZN = IZ(N)
!
!---      Loop over number of hydrostatic horizons  ---
!
          IF( ISIC.NE.4 ) THEN
            IFIND = 0
            DO NH = 1,NHORZ
              NHX = (NH-1)*10
!
!---          Identify the hydrostatic horizon  ---
!
              IF( ZP(N).GE.VSLC(NHX+9).AND.ZP(N).LE.VSLC(NHX+10) ) THEN
                IFIND = NH
                EXIT
              ENDIF
            ENDDO
            IF( IFIND.EQ.0 ) THEN
              INDX = 12
              IMSGX = N
              NMSGX = 0
              RLMSGX = 0.D+0
              SUBLOGX = 'CHK_GT'
              CHMSGX = 'Zonal Hydrostatic Initial Condition: ' // 
     &          'Node not within Zone: Node = '
               CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
          ELSE
            NHX = 0
          ENDIF
          DISTZ = (ZP(N)-VSLC(NHX+2))/25.D+0
!
!---      Loop over twenty five increments to compute initial
!         conditions at node  ---
!
          PLX = VSLC(NHX+1) + PATM
          ZX = VSLC(NHX+2)
          DO M = 1,25
            ZX = ZX + DISTZ
            TX = VSLC(NHX+3) + (ZX-VSLC(NHX+4))*VSLC(NHX+5)
            YLSX = VSLC(NHX+6) + (ZX-VSLC(NHX+7))*VSLC(NHX+8)
!
!---        Check for out-of-range salt mass fraction  ---
!
            INDX = 0
            CALL REGION_4( TX,PSWX,INDX )
            PSWX = MIN( PSWX,PCRW )
            PX = MAX( PLX,PSWX )
            CALL SOL_BRNS( TX,PX,XLSMX )
            IF( YLSX.LT.0.D+0 ) THEN
              INDX = 17
              IMSGX = N
              NMSGX = n
              SUBLOGX = 'CHK_GT'
              RLMSGX = YLSX
              CHMSGX = 'Hydrostatic Initial Condition: Salt ' //
     &          'Mass Fraction < 0.0 : XLS = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ELSEIF( YLSX.GT.XLSMX ) THEN
              INDX = 17
              IMSGX = N
              NMSGX = N
              SUBLOGX = 'CHK_GT'
              RLMSGX = YLSX
              CHMSGX = 'Hydrostatic Initial Condition: Salt ' //
     &          'Mass Fraction > Solubility Limit : XLS = '
              CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
            CALL P_IAPWS( TX,PLX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL DENS_B( YLSX,RHOBX,RHOLWX,TX )
            RHOLX = RHOBX
            PLX = PLX - RHOLX*GRAVZ*DISTZ
          ENDDO
!
!---      Assign initial phase condition  ---
!
          NPHAZ(2,N) = 1
!
!---      Assign gas-entry pressure for non Brooks-Corey;
!         Brooks-Corey; Brooks-Corey, Dual Porosity; and
!         Brooks-Corey, Entrapment  ---
!
          IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &      ISCHR(IZN).EQ.202 ) THEN
            ENPR = SCHR(1,IZN)*RHORL*GRAV
          ELSEIF( ISCHR(IZN).EQ.4 ) THEN
            ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
          ELSE
            ENPR = 0.D+0
          ENDIF
          PL(2,N) = PLX - PATM
          PG(2,N) = PL(2,N) + ENPR - EPSLX
          T(2,N) = TX
          YLS(2,N) = YLSX
          SL(2,N) = 1.D+0
          SG(2,N) = 0.D+0
          ICBRN(N) = 3
          ICAIR(N) = 3

          POR0(1,N) = POR0(1,N)
          POR0(2,N) = POR0(2,N)

        ENDDO
!
!---  Hydrostatic with geothermal gradient initial conditions, either
!     zonal or entire domain, with salt mass fraction
!     read from an external file  ---
!
      ELSEIF( MOD(ISIC,10).EQ.5 ) THEN
!
!---    Loop over all nodes  ---
!
        DO N = 1,NFLD
!
!---      Skip inactive nodes  ---
!
          IF( IXP(N).EQ.0 ) CYCLE
          IZN = IZ(N)
!
!---      Loop over number of hydrostatic horizons  ---
!
          IFIND = 0
          DO NH = 1,NHORZ
            NHX = (NH-1)*10
!
!---        Identify the hydrostatic horizon  ---
!
            IF( ZP(N).GE.VSLC(NHX+9) .AND. ZP(N).LE.VSLC(NHX+10) ) THEN
              IFIND = NH
              EXIT
            ENDIF
          ENDDO
          IF( IFIND.EQ.0 ) THEN
            INDX = 12
            IMSGX = N
            NMSGX = 0
            RLMSGX = 0.D+0
            SUBLOGX = 'CHK_GT'
            CHMSGX = 'Zonal Hydrostatic Initial Condition: ' // 
     &        'Node not within Zone: Node = '
             CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
          ENDIF
          DISTZ = (ZP(N)-VSLC(NHX+2))/25.D+0
!
!---      Loop over twenty five increments to compute initial
!         conditions at node  ---
!
          PLX = VSLC(NHX+1) + PATM
          ZX = VSLC(NHX+2)
          DO M = 1,25
            ZX = ZX + DISTZ
            TX = VSLC(NHX+3) + (ZX-VSLC(NHX+4))*VSLC(NHX+5)
            YLSX = YLS(2,N)
!
!---        Check for out-of-range salt mass fraction  ---
!
            INDX = 0
            CALL REGION_4( TX,PSWX,INDX )
            PSWX = MIN( PSWX,PCRW )
            PX = MAX( PLX,PSWX )
            CALL SOL_BRNS( TX,PX,XLSMX )
            IF( YLSX.LT.0.D+0 ) THEN
              INDX = 17
              NMSGX = 0
              SUBLOGX = 'CHK_GT'
              RLMSGX = YLSX
              CHMSGX = 'Hydrostatic Initial Condition: Salt ' //
     &          'Mass Fraction < 0.0 : XLS = '
               CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ELSEIF( YLSX.GT.XLSMX ) THEN
              INDX = 17
              NMSGX = 0
              SUBLOGX = 'CHK_GT'
              RLMSGX = YLSX
              CHMSGX = 'Hydrostatic Initial Condition: Salt ' //
     &          'Mass Fraction > Solubility Limit : XLS = '
               CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
            ENDIF
            CALL P_IAPWS( TX,PLX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL DENS_B( YLSX,RHOBX,RHOLWX,TX )
            RHOLX = RHOBX
            PLX = PLX - RHOLX*GRAVZ*DISTZ
          ENDDO
!
!---      Assign initial phase condition  ---
!
          NPHAZ(2,N) = 1
!
!---      Assign gas-entry pressure for non Brooks-Corey;
!         Brooks-Corey; Brooks-Corey, Dual Porosity; and
!         Brooks-Corey, Entrapment  ---
!
          IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &      ISCHR(IZN).EQ.202 ) THEN
            ENPR = SCHR(1,IZN)*RHORL*GRAV
          ELSEIF( ISCHR(IZN).EQ.4 ) THEN
            ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
          ELSE
            ENPR = 0.D+0
          ENDIF
          PL(2,N) = PLX - PATM
          PG(2,N) = PL(2,N) + ENPR - EPSLX
          T(2,N) = TX
          SL(2,N) = 1.D+0
          SG(2,N) = 0.D+0
          ICBRN(N) = 3
          ICAIR(N) = 3

          POR0(1,N) = POR0(1,N)
          POR0(2,N) = POR0(2,N)

        ENDDO
      ENDIF
!
!---  Convert initial aqueous-salt inputs into
!     aqueous mass fractions, considering initial saturation,
!     skip for hydrostatic initial conditions  ---
!
      IF( ISIC.LE.3 ) THEN
          DO N = 1,NFLD
!
!---      Skip inactive nodes  ---
!
          IF( IXP(N).EQ.0 .AND. ISIC.NE.4 ) CYCLE
          IF( IEO.EQ.2 .AND. NPHAZ(2,N).NE.0 ) CYCLE
          IZN = IZ(N)
          N_DB = N
!
!---      Assign gas-entry pressure for non Brooks-Corey;
!         Brooks-Corey; Brooks-Corey, Dual Porosity; and
!         Brooks-Corey, Entrapment  ---
!
          IF( IXP(N).EQ.0 ) THEN
            ENPR = 0.D+0
            ESGTMX = 0.D+0
          ELSE
            IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &        ISCHR(IZN).EQ.202 ) THEN
              ENPR = SCHR(1,IZN)*RHORL*GRAV
            ELSEIF( ISCHR(IZN).EQ.4 ) THEN
              ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
            ELSE
              ENPR = 0.D+0
            ENDIF
!
!---        Initial trapped gas saturation for the van Genuchten or
!           Brooks/Corey entrapment model  ---
!
            IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 .OR.
     &        ISCHR(IZN).EQ.201 .OR. ISCHR(IZN).EQ.202 ) THEN
!
!---          w/ Webb extension  ---
!
              IF( ISM(IZN).EQ.2 ) THEN
                ESGTMX = SCHR(15,IZN)
!
!---          w/o Webb extension  ---
!
              ELSE
                ESGTMX = SCHR(15,IZN)/(1.D+0-SCHR(4,IZN))
              ENDIF
            ELSE
              ESGTMX = 0.D+0
            ENDIF
          ENDIF
!
!---      Non-hydrostatic initial conditions  ---
!
          IF( ISIC.LE.3 ) THEN
!
!---        Gas and aqueous pressure specified, determine aqueous
!           saturation without scaling  ---
!
            IF( ISIC.EQ.3 ) THEN
              BTGLX = 1.D+0
              INDX = 1
              CALL SP_GT( ASLFX,ASLMX,ASLX,ASLMINX,ESLX,ESGTX,ESGTMX,
     &          PG(2,N),PL(2,N),BTGLX,SGX,SGT(2,N),SL(2,N),SLFX,
     &          SLMX,SLRX,INDX,IZN )
            ENDIF
!
!---        Assign initial phase condition  ---
!
            PLX = PL(2,N)+PATM
            PGX = PG(2,N)+PATM
            TKX = T(2,N)+TABS
!
!---        Saturated conditions with no trapped gas  ---
!
            IF( (1.D+0-SL(2,N))/EPSL.LT.EPSL ) THEN
              SL(2,N) = 1.D+0
              SGT(2,N) = 0.D+0
              NPHAZ(2,N) = 1
!
!---        Saturated conditions with trapped gas  ---
!
            ELSEIF( (1.D+0-SL(2,N)-SGT(2,N))/EPSL.LT.EPSL ) THEN
              SL(2,N) = 1.D+0-SGT(2,N)
              NPHAZ(2,N) = 3
!
!---        Unsaturated conditions  ---
!
            ELSEIF( SL(2,N)/EPSL.LT.EPSL .OR. TKX.GT.TCRW ) THEN
              SL(2,N) = 0.D+0
              SGT(2,N) = 0.D+0
              NPHAZ(2,N) = 4
!
!---        Partially saturated conditions  ---
!
            ELSE
              NPHAZ(2,N) = 2
            ENDIF
!
!---        Saturated system w/ entrapped gas
!
!           Energy - temperature
!           Water mass - aqueous pressure
!           Air mass - trapped gas saturation
!           NaCl mass - total NaCl brine mass fraction  ---
!
            IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.3 ) THEN
              CALL FLH_IC1( BTGL(2,N),PGX,PLX,PVA(2,N),SGT(2,N),SL(2,N),
     &         T(2,N),TMS(2,N),XMLA(2,N),YLS(2,N),ICBRN(N),ICAIR(N),
     &         ISIC,N )
!
!---        Unsaturated system w/ or w/o entrapped gas
!         
!           Energy - temperature
!           Water mass - aqueous pressure
!           Air mass - gas pressure
!           NaCl mass - total NaCl brine mass fraction  ---
!
            ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
              CALL FLH_IC2( BTGL(2,N),PGX,PLX,SGT(2,N),SL(2,N),
     &          T(2,N),TMS(2,N),XMLA(2,N),YLS(2,N),ICBRN(N),ISIC,N )
!
!---        Fully unsaturated conditions
!
!           Energy - temperature
!           Water mass - water vapor partial pressure
!           Air mass - gas pressure
!           NaCl mass - salt mass  ---
!
            ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
              CALL FLH_IC4( BTGL(2,N),PGX,PLX,PVA(2,N),PVW(2,N),
     &          SGT(2,N),SL(2,N),T(2,N),TMS(2,N),XMGA(2,N),YLS(2,N),
     &          ICAIR(N),ISIC,N )
            ENDIF
            PG(2,N) = PGX - PATM
            PL(2,N) = PLX - PATM
          ENDIF
        ENDDO
      ENDIF
!
!---  Compressibility, porosity, and small limits loop  ---
!
      DO N = 1,NFLD
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
        IF( IEO.EQ.2 .AND. NPHAZ(2,N).NE.0 ) CYCLE
        IZN = IZ(N)
!
!---    Small concentration limits  ---
!
        IF( YLS(2,N).LT.1.D-15 ) YLS(2,N) = 0.D+0
        IF( XLS(2,N).LT.1.D-15 ) XLS(2,N) = 0.D+0
        IF( XLA(2,N).LT.1.D-15 ) XLA(2,N) = 0.D+0
!
!---    Initialize reference pressure for compressibility  ---
!
        IF( CMP(3,IZN).GT.PATM ) THEN
          PCMP(N) = CMP(3,IZN)
        ELSEIF( ISLC(61).EQ.0 ) THEN
          PCMP(N) = MAX( PL(2,N),PG(2,N) )+PATM
        ENDIF
!
!---    Porous-media porosity  ---
!
        PX = MAX( PGX,PLX )
        CALL PORSTY_GT( N,PX,PCMP(N),PORD(2,N),PORT(2,N) )
        PORD(2,N) = MAX( PORD(2,N),EPSL )
        PORT(2,N) = MAX( PORT(2,N),PORD(2,N) )

        POR0(1,N) = POR0(1,N)
        POR0(2,N) = POR0(2,N)

      ENDDO
!
!---  Initialize surface cover primary variables  ---
!
      IF( NSFCN.GT.0 ) THEN
        DO NSCN = 1,NSFCN
          N = ICM_SFC(1,NSCN)
          NQZ = NSZ(N)+IJFLD
          DO M = 1,2
            T_GS(M,NSCN) = T(2,N)
            T_PL(M,NSCN) = T(2,N)
            T_CP(M,NSCN) = T(2,N)
            PL_GS(M,NSCN) = PL(2,N) - DZGP(NQZ)*RHOL(2,N)*GRVZ(NQZ)
            PCX = PG(2,N)-PL(2,N)
            INDX = 0
            CALL REGION_4( TX,PSWX,INDX )
            CALL P_IAPWS( T(2,N),PSWX,RHOX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL VPL_B( T(2,N),PSWX,PCX,RHOLWX,PVWX )     
            PVW_CP(M,NSCN) = PVWX
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC_GT
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
!     Compute initial solute concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE GRID
      USE FDVP
      USE CONST
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/CISC_GT'
      DO NSL = 1,NSOLU
        DO N = 1,NFLD
          IF( IXP(N).EQ.0 ) CYCLE
          IZN = IZ(N)
!
!---      Stationary solute  ---
!
          IF( IEDL(NSL).EQ.4  ) THEN
            XVS = 0.D+0
          ELSEIF( IPCL(NSL).EQ.2 ) THEN
            XVS = SL(2,N)*RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
          ELSE
            XVS = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
          ENDIF
          XVL = SL(2,N)*PORD(2,N)
          XVG = SG(2,N)*PORD(2,N)
!
!---      Stationary solute  ---
!
          IF( IEDL(NSL).EQ.4  ) THEN
            PCGLX = 1.D+0
!
!---      Constant gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.0 ) THEN
            PCGLX = PCGL(1,NSL)
!
!---      Temperature dependent gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.1 ) THEN
            TK = T(2,N)+TABS
            PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &        + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &        + PCGL(5,NSL)*TK**2 )
!
!---      Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.2 ) THEN
            PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
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
          YL(N,NSL) = XVL*YVL
          YG(N,NSL) = XVG*YVG
!
!---      Phase-volumetric concentration ratios  ---
!
          IF( ICT(N,NSL).EQ.2 ) THEN
            C(N,NSL) = C(N,NSL)*(XVS + XVL + XVG*PCGLX)
          ELSEIF( ICT(N,NSL).EQ.-1 ) THEN
            C(N,NSL) = C(N,NSL)*RHOS(IZN)*(1.D+0-PORT(2,N))
          ELSEIF( ICT(N,NSL).EQ.3 ) THEN
            C(N,NSL) = C(N,NSL)*((XVS + XVL)/PCGLX + XVG)
!
!---      Convert from activity to molar mass  ---
!
          ELSEIF( ICT(N,NSL).EQ.4 ) THEN
            C(N,NSL) = C(N,NSL)/(LOG(2.D+0)/HLF(NSL))/VOL(N)
          ENDIF
        ENDDO
!
!---    Assign boundary solute concentrations for initial condition
!       type boundary conditions  ---
!
        DO NB = 1,NBC
          IF( IBCT(NSL+LUK,NB).EQ.12 ) THEN
            N = IBCN(NB)
            CBO(NB,NSL) = C(N,NSL)
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CISC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLUX_GT
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
!     Matrix, fracture and borehole fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 17 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLUX_GT'
!      time_flux_start = OMP_GET_WTIME()
!
!---  Compute aqueous-phase volumetric flux (interior surfaces)  ---
!
      CALL DRCVL_GT
!      time_drcvl_end = OMP_GET_WTIME()
!
!---  Compute gas-phase volumetric flux (interior surfaces)  ---
!
      IF( ISLC(37).EQ.0 ) CALL DRCVG_GT
!      time_drcvg_end = OMP_GET_WTIME()
!
!---  Compute water vapor diffusion flux through the gas phase
!     (interior surfaces)  ---
!
      IF( ISLC(37).EQ.0 ) CALL DFFGW_GT
!      time_dffgw_end = OMP_GET_WTIME()
!
!---  Compute dissolved air diffusion flux through the aqueous phase
!     (interior surfaces)  ---
!
      IF( ISLC(37).EQ.0 ) CALL DFFLA_GT
!      time_dffla_end = OMP_GET_WTIME()
!
!---  Compute aqueous-salt diffusion flux through the aqueous phase
!     (interior surfaces), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS_GT
!      time_dffls_end = OMP_GET_WTIME()
!
!---  Nonisothermal simulations  ---
!
      IF( ISLC(30).EQ.0 )  THEN
!
!---    Compute diffusive thermal flux (interior surfaces)  ---
!
        CALL THD_GT
!        time_thd_end = OMP_GET_WTIME()
!
!---    Compute aqueous-phase advective thermal flux
!       (interior surfaces)  ---
!
        CALL THAL_GT
!        time_thal_end = OMP_GET_WTIME()
!
!---    Compute gas-phase advective thermal flux
!       (interior surfaces)  ---
!
        IF( ISLC(37).EQ.0 ) CALL THAG_GT
!        time_thag_end = OMP_GET_WTIME()
!
!---    Compute vapor diffusive thermal flux
!       (interior surfaces)  ---
!
        IF( ISLC(37).EQ.0 ) CALL THDG_GT
!        time_thdg_end = OMP_GET_WTIME()
      ENDIF
!
!---  Compute aqueous-phase volumetric flux, gas-phase volumetric flux,
!     water vapor mass flux, diffusive thermal flux, aqueous-phase
!     advective thermal flux, gas-phase advective thermal flux,
!     vapor diffusive thermal flux, aqueous-phase salt flux
!     (boundary surfaces)  ---
!
      CALL BCF_GT
!      time_bcf_end = OMP_GET_WTIME()
!      print *,'drcvl = ',time_drcvl_end-time_flux_start,
!     &  ' drcvg = ',time_drcvg_end-time_drcvl_end,
!     &  ' dffgw = ',time_dffgw_end-time_drcvg_end,
!     &  ' dffla = ',time_dffla_end-time_dffgw_end,
!     &  ' dffls = ',time_dffls_end-time_dffla_end,
!     &  ' thd = ',time_thd_end-time_dffls_end,
!     &  ' thal = ',time_thal_end-time_thd_end,
!     &  ' thag = ',time_thag_end-time_thal_end,
!     &  ' thdg = ',time_thdg_end-time_thag_end,
!     &  ' bcf = ',time_bcf_end-time_thdg_end
!
!---  Compute surface cover fluxes  --
!
      IF( NSFCN.GT.0 ) CALL FLUX_SFC
!
!---  Fracture flow and transport solution  ---
!
      IF( ISLC(74).NE.0 ) THEN
!
!---    Compute advective fluxes from borehole node to 
!       fracture triangle  ---
!
        IF( ISLC(74).EQ.3 ) CALL DRCV_BTF_GT
!
!---    Compute diffusive fluxes from borehole node to 
!       fracture triangle  ---
!
        IF( ISLC(74).EQ.3 .AND. ISLC(37).EQ.0 ) CALL DFFX_BTF_GT
!
!---    Compute borehole and fracture aqueous-phase volumetric flux  ---
!
        CALL DRCVL_BF_GT
!
!---    Compute borehole and frature gas-phase volumetric flux  ---
!
        IF( ISLC(37).EQ.0 ) CALL DRCVG_BF_GT
!
!---    Compute borehole and fracture water-vapor flux  ---
!
        IF( ISLC(37).EQ.0 ) CALL DFFGW_BF_GT
!
!---    Compute borehole and fracture dissolved-air flux  ---
!
        IF( ISLC(37).EQ.0 ) CALL DFFLA_BF_GT
!
!---    Isobrine option ---
!
        IF( ISLC(32).EQ.0 ) THEN
!
!---      Compute borehole and fracture advective and diffusive
!         salt flux  ---
!
          CALL DFFLS_BF_GT
!
!---      Compute borehole to fracture salt flux  ---
!
          IF( ISLC(74).EQ.3 ) CALL DFFLS_BTF_GT
        ENDIF
!
!---    Isothermal option  ---
!
        IF( ISLC(30).EQ.0 )  THEN
!
!---      Compute thermal fluxes from borehole node to 
!         fracture triangle  ---
!
          IF( ISLC(74).EQ.3 ) CALL THFX_BTF_GT
!
!---      Compute borehole and fracture thermal diffusive flux  ---
!
          CALL THD_BF_GT
!
!---      Compute borehole and fracture thermal advective flux  ---
!
          CALL THA_BF_GT
!
!---      Compute borehole and fracture thermal gas component 
!         diffusive flux  ---
!
          IF( ISLC(37).EQ.0 ) CALL THDG_BF_GT
        ENDIF
!
!---    Fracture to matrix transfer functions  ---
!
        IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL TRNS_FRC_GT
!
!---    Borehole to matrix transfer functions  ---
!
        IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL TRNS_BH_GT
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLUX_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HYDST_BC_GT( BCX,PGX,PLX,TX,XLSX,YLSX,ZPX )
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
!     Establish hydrostatic boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 26 July 2017
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
      REAL*8 BCX(*)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HYDST_BC_GT'
!
!---  Hydrostatic initial conditions, establish aqueous pressure,
!     temperature, and salt mass fraction at boundary surface  ---
!
      ZX = BCX(3)
      KC = MAX( 1,INT(ABS(ZX-ZPX)) )
      DISTZ = (ZX-ZPX)/REAL(KC)
      PLX = BCX(2) + PATM
      PGX = PLX
      XLAX = 0.D+0
      DO 100 K = 1,KC
        ZX = ZX - 5.D-1*DISTZ
        TX = BCX(4) + (ZX-BCX(5))*BCX(6)
        YLSX = BCX(7) + (ZX-BCX(8))*BCX(9)
        ZX = ZX - 5.D-1*DISTZ
!
!---    Check for out-of-range salt mass fraction  ---
!
        CALL SOL_BRNS( TX,PLX,XLSMX )
        IF( YLSX.LT.0.D+0 ) THEN
          INDX = 14
          RLMSG = YLSX
          CHMSG = 'Hydrostatic Boundary Condition: Salt ' //
     &      'Mass Fraction < 0.0 : XLS = '
           CALL WRMSGS( INDX )
        ELSEIF( YLSX.GT.XLSMX ) THEN
          INDX = 17
          RLMSG = YLSX
          CHMSG = 'Hydrostatic Boundary Condition: Salt ' //
     &      'Mass Fraction > Solubility Limit : XLS = '
           CALL WRMSGS( INDX )
        ENDIF
        CALL P_IAPWS( TX,PLX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        XLSX = YLSX
        CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
        RHOLX = RHOBX
        PLX = PLX + RHOLX*GRAVZ*DISTZ
        PGX = PLX
  100 CONTINUE
      PGX = PLX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HYDST_BC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_GT
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
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
      USE GRID
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
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM_GT'
      EPSLX = 1.D-4
!
!---  Phase options, compute phase condition   ---
!
      DO N = 1,NFLD
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
        IZN = IZ(N)
        N_DB = N
!
!---    Minimum apparent aqueous saturation  ---
!
        ASLMINX = ASLMIN(2,N)
!
!---    Assign gas entry pressure and minimum gas saturation
!       for transition to unsaturated conditions---
!
        SGMNX = 1.D-3
        ENPR = 0.D+0
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
          SGMNX = 1.D+1**(-3.D+0+LOG10(1.D+0/SCHR(1,IZN)))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
          SGMNX = 1.D+1**(-3.D+0+LOG10(SCHR(1,IZN)))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ELSEIF( ISCHR(IZN).EQ.3 ) THEN
          SCHRX = MAX( SCHR(1,IZN),SCHR(5,IZN) )
          SGMNX = 1.D+1**(-3.D+0+LOG10(1.D+0/SCHRX))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
          SCHRX = MIN( SCHR(1,IZN),SCHR(5,IZN) )
          SGMNX = 1.D+1**(-3.D+0+LOG10(SCHRX))
          SGMNX = MIN( MAX( SGMNX,1.D-4 ),1.D-3 )
        ENDIF
!
!---    Maximum effective trapped gas saturation  ---
!
        IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 ) THEN
!
!---      w/ Webb extension  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            ESGTMX = SCHR(15,IZN)
!
!---      w/o Webb extension  ---
!
          ELSE
            ESGTMX = SCHR(15,IZN)/(1.D+0-SCHR(4,IZN))
          ENDIF
!
!---    Zero maximum trapped gas saturation  ---
!
        ELSE
          ESGTMX = 0.D+0
        ENDIF
!
!---    Absolute temperature  ---
!
        TKX = T(2,N) + TABS
!
!---    Saturated system w/o entrapped gas
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - aqueous-air mole fraction
!       Air mass - gas air partial pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          PG(2,N) = PL(2,N) + ENPR
          PLX = PL(2,N) + PATM
          INDX = 0
          CALL REGION_4( T(2,N),PSWX,INDX )
          PWX = MAX( PSWX,PLX )
          CALL SCF_GL( BTGLX,PWX )
          CALL SOL_BRNS( T(2,N),PWX,XLSMX )
          XLSX = MIN( YLS(2,N),XLSMX )
          XLS(2,N) = XLSX
          CALL SP_B( T(2,N),XLSX,PSBX )
          PVBX = PSBX
          PLX = PL(2,N) + PATM
!          PGAX = XMLA(2,N)*HCAW
          PGAX = PVA(2,N)
!
!---      Isoair option no transition from aqueous 
!         saturated conditions  ---
!
          IF( ISLC(37).EQ.0 ) THEN
            PGX = PVBX + PGAX
          ELSE
            PGX = PLX
          ENDIF
          PX = MAX( PGX,PLX )
          INDX = 0
          CALL KSP_GT( N,PGX,PLX,BTGLX,SGX,SGTX,SLX,SLFX,SLMX,
     &      RKLX,RKGX,ASLX,ASLMINX,ESGTX,ESGTMX,
     &      SLRX,INDX )
!
!---      Transition from aqueous saturated conditions to
!         aqueous unsaturated conditions: PC1 -> PC2  ---
!
          IF( SGX.GT.1.D-3 ) THEN
            SGX = 1.D-3
            SLX = 1.D+0 - SGX
            CALL CAP_GT( ASLMINX,BTGLX,CPGL,SLX,SGTX,IZN )
            PCX = CPGL
            CALL P_IAPWS( T(2,N),PSWX,RHOX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL VPL_B( T(2,N),PSWX,PCX,RHOLWX,PVW(2,N) )
            CALL SCF_GL( BTGLX,PVW(2,N) )
            CALL SOL_BRNS( T(2,N),PSWX,XLSMX )
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL DENS_B( XLSX,RHOBX,RHOLWX,T(2,N) )
            CALL SP_B( T(2,N),XLSX,PSBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX )
            PGX = PVBX + PGAX
            PG(2,N) = PGX - PATM
            PL(2,N) = PG(2,N) - CPGL
            NPHAZ(2,N) = 2
!
!---      No transition from aqueous saturated conditions  ---
!
          ELSE
            NPHAZ(2,N) = 1
          ENDIF
!
!---    Unsaturated system w/ or w/o entrapped gas
!       
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          PLX = PL(2,N) + PATM
          PGX = PG(2,N) + PATM
          PX = MAX( PGX,PLX )
          INDX = 0
          CALL REGION_4( T(2,N),PSWX,INDX )
          CALL P_IAPWS( T(2,N),PSWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          PCX = PGX - PLX
          CALL VPL_B( T(2,N),PSWX,PCX,RHOLWX,PVW(2,N) )
          CALL SCF_GL( BTGLX,PVW(2,N) )
          CALL SOL_BRNS( T(2,N),PSWX,XLSMX )
          XLSX = MIN( YLS(2,N),XLSMX )
          XLS(2,N) = XLSX
          CALL DENS_B( XLSX,RHOBX,RHOLWX,T(2,N) )
          CALL SP_B( T(2,N),XLSX,PSBX )
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVW(2,N) )
          PGAX = PGX - PVW(2,N)
          INDX = 0
          CALL KSP_GT( N,PGX,PLX,BTGLX,SGX,SGTX,SLX,SLFX,SLMX,
     &      RKLX,RKGX,ASLX,ASLMINX,ESGTX,ESGTMX,
     &      SLRX,INDX )
!         
!---      Apparent aqueous saturation (i.e., aqueous saturation +
!         trapped gas saturation = 1.0), transition to saturation
!         conditions  ---
!
          IF( ABS(1.D+0-ASLX).LT.EPSL ) THEN
!
!---        Trapped gas exists transition to saturated condition
!           w/ entrapped gas: PC2 -> PC3  ---
!
            IF( ESGTMX.GT.EPSL .AND. SGTX.GT.EPSL ) THEN
              CALL CAP_GT( ASLMIN(2,N),BTGLX,PCX,SLX,SGTX,IZN )
              PG(2,N) = PL(2,N) + PCX
              SG(2,N) = SGX
              SGT(2,N) = SGX
              NPHAZ(2,N) = 3
!
!---        No trapped gas transition to saturated condition
!           w/o entrapped gas: PC2 -> PC1  ---
!
            ELSE
              PX = PL(2,N)+PATM
              CALL SP_B( T(2,N),XLSX,PSBX )
              PVBX = PSBX
              PG(2,N) = PL(2,N) + ENPR - EPSLX
              NPHAZ(2,N) = 1
            ENDIF
!
!---      No aqueous phase, transition
!         to fully unsaturated condition: PC2 -> PC4  ---
!
          ELSEIF( ( SLX.LT.EPSL .AND. (1.D+0-SL(1,N)).GT.EPSL ) .OR.
     &      TKX.GT.(TCRW+2.5D+0) ) THEN
            PGX = PG(2,N) + PATM
            INDX = 0
            CALL REGION_4( T(2,N),PSWX,INDX )
            CALL SOL_BRNS( T(2,N),PSWX,XLSMX )
            IF( TMS(2,N).GT.EPSL ) THEN
              XLSX = XLSMX
            ELSE
              XLSX = 0.D+0
            ENDIF
            XLS(2,N) = XLSX
            SL(2,N) = 0.D+0
            SGT(2,N) = 0.D+0
            CALL CAP_GT( ASLMIN(2,N),BTGLX,PCX,SL(2,N),SGT(2,N),IZN )
            PL(2,N) = PG(2,N) - PCX
            PLX = PGX - PCX
            CALL SP_B( T(2,N),XLSX,PSBX )
            PX = MAX( PGX,PSBX )
            CALL P_IAPWS( T(2,N),PX,RHOX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL DENS_B( XLSX,RHOBX,RHOLWX,T(2,N) )
            PCX = PGX - PLX
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVW(2,N) )
            PVA(2,N) = PGX - PVW(2,N)
            NPHAZ(2,N) = 4
!
!---      Gas pressure remains above gas-entry pressure, remain
!         as unsaturated condition  ---
!
          ELSE
            NPHAZ(2,N) = 2
          ENDIF
!
!---    Saturated system w/ entrapped gas
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
          SLRX = SCHR(4,IZN)
!
!---      w/ Webb extension  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            ESGTX = SGT(2,N)
            ESGTOX = SGT(1,N)
!
!---      w/o Webb extension  ---
!
          ELSE
            ESGTX = SGT(2,N)/(1.D+0-SLRX)
            ESGTOX = SGT(1,N)/(1.D+0-SLRX)
          ENDIF
!
!---      Test against previous effective trapped gas saturation when
!         previous phase condition is saturated w/ entrapped gas  ---
!
          IF( NPHAZ(1,N).EQ.3 ) THEN
            ESGTPX = ESGTOX
          ELSE
            ESGTPX = ESGTX
          ENDIF
!
!---      Reversal point floats with trapped gas saturation  ---
!
          IF( (ESGTMX-ESGTX).GT.EPSL ) THEN
            ASLMIN(2,N) = (ESGTMX-ESGTX)/(ESGTMX + ESGTMX*ESGTX - ESGTX)
          ELSE
            ASLMIN(2,N) = 0.D+0
          ENDIF
          PX = PL(2,N)+PATM
          XLS(2,N) = XLSX
          CALL SP_B( T(2,N),XLSX,PSBX )
          PVBX = PSBX
!
!---      Trapped gas saturation disappears, transition to
!         saturated condition w/o entrapped gas: PC3 -> PC1  ---
!
          IF( SG(2,N).LT.EPSL ) THEN
            SGT(2,N) = 0.D+0
            SG(2,N) = 0.D+0
            PG(2,N) = PL(2,N) + ENPR - EPSLX
            ASLMIN(2,N) = 1.D+0
            NPHAZ(2,N) = 1
!
!---      Trapped gas saturation exceeds previous trapped gas
!         saturation or maximum trapped gas saturation, transition 
!         to free gas w/ trapped gas phase condition: PC3 -> PC2  ---
!
          ELSEIF( ((ESGTX-ESGTPX).GT.EPSL) .OR.
     &      ((ESGTX-ESGTMX).GT.EPSL) ) THEN
!
!---        w/ Webb extension  ---
!
            IF( ISM(IZN).EQ.2 ) THEN
              ASLX = 1.D+0 - ESGTX + ESGTOX
              ASLMIN(2,N) = MIN( ASLMINX,ASLX )
              SGT(2,N) = ESGTOX
              SLX = ASLX - ESGTOX 
!
!---        w/o Webb extension  ---
!
            ELSE
              ASLX = 1.D+0 - ESGTX + ESGTOX
              ASLMIN(2,N) = MIN( ASLMINX,ASLX )
              SGT(2,N) = ESGTOX*(1.D+0-SLRX)
              SLX = (ASLX-ESGTOX)*(1.D+0-SLRX) + SLRX
            ENDIF
            CALL CAP_GT( ASLMIN(2,N),BTGLX,PCX,SLX,SGT(2,N),IZN )
            PCX = CPGL
            CALL P_IAPWS( T(2,N),PSWX,RHOX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL VPL_B( T(2,N),PSWX,PCX,RHOLWX,PVW(2,N) )
            CALL SCF_GL( BTGLX,PVW(2,N) )
            CALL SOL_BRNS( T(2,N),PSWX,XLSMX )
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL DENS_B( XLSX,RHOBX,RHOLWX,T(2,N) )
            CALL SP_B( T(2,N),XLSX,PSBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX )
            PGX = PVBX + PGAX
            PG(2,N) = PGX - PATM
            PL(2,N) = PG(2,N) - CPGL
            NPHAZ(2,N) = 2
!
!---      Trapped gas saturation neither disappears nor
!         exceeds previous trapped gas saturation,
!         remain as saturated w/ entrapped gas condition  ---
!
          ELSE
            IF( ISCHR(IZN).EQ.202 ) THEN
              PG(2,N) = PL(2,N) + SCHR(5,IZN) - EPSLX
            ELSE
              PG(2,N) = PL(2,N) + ENPR - EPSLX
            ENDIF
            NPHAZ(2,N) = 3
          ENDIF
!
!---    Fully unsaturated conditions
!
!       Energy - temperature
!       Water mass - water vapor partial pressure
!       Air mass - air gas partial pressure
!       NaCl mass - salt mass  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
          PGX = PVA(2,N) + PVW(2,N)
          PG(2,N) = PGX - PATM
!
!---      Maximum salt solubility if aqueous phase existed  ---
!
          INDX = 0
          CALL REGION_4( T(2,N),PSWX,INDX )
          CALL SOL_BRNS( T(2,N),PSWX,XLSMX )
          IF( TMS(2,N).GT.EPSL ) THEN
            XLSX = XLSMX
          ELSE
            XLSX = 0.D+0
          ENDIF
          XLS(2,N) = XLSX
          CALL SCF_GL( BTGLX,PVW(2,N) )
          SL(2,N) = 0.D+0
          SGT(2,N) = 0.D+0
          CALL CAP_GT( ASLMIN(2,N),BTGLX,PCX,SL(2,N),SGT(2,N),IZN )
          PL(2,N) = PG(2,N) - PCX
          PLX = PGX - PCX
          CALL P_IAPWS( T(2,N),PVW(2,N),RHOGWX,RHOLWX,HGWX,HLWX,
     &      UGWX,ULWX )
!
!---      Saturated water vapor pressure given temperature and
!         capillary pressure  ---
!
          CALL SP_B( T(2,N),XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,T(2,N) )
          PCX = PGX - PLX
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX )
!
!---      Aqueous phase appears, transition to
!         unsaturated conditions: PC4 -> PC2  ---
!
          IF( PVW(2,N).GT.PSBX .AND. TKX.LE.(TCRW-2.5D+0) ) THEN
            SLX = MIN( MAX( 1.D-4,(1.D+1*(PVW(2,N)-PSBX)/PSBX) ),1.D-1 )
            SGTX = 0.D+0
            CALL CAP_GT( ASLMIN(2,N),BTGLX,CPGL,SLX,SGTX,IZN )
            PCX = CPGL
            CALL P_IAPWS( T(2,N),PSWX,RHOX,RHOLWX,HX,HLWX,UX,ULWX )
            CALL VPL_B( T(2,N),PSWX,PCX,RHOLWX,PVWX )
            CALL SCF_GL( BTGLX,PVWX )
            CALL SOL_BRNS( T(2,N),PSWX,XLSMX )
            YLS(2,N) = TMS(2,N)/(RHOLWX*SLX*PORD(2,N)+TMS(2,N))
            XLSX = MIN( YLS(2,N),XLSMX )
            XLS(2,N) = XLSX
            CALL DENS_B( XLSX,RHOBX,RHOLWX,T(2,N) )
            CALL SP_B( T(2,N),XLSX,PSBX )
            CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX )
            PGX = PVBX + PVA(2,N)
            PG(2,N) = PGX - PATM
            PL(2,N) = PG(2,N) - CPGL
            NPHAZ(2,N) = 2
!
!---      No aqueous phase, no transition from
!         fully unsaturated condition  ---
!
          ELSE
            NPHAZ(2,N) = 4
          ENDIF
        ENDIF
!
!---    Compute increments  ---
!
!
!---    Saturated system w/o entrapped gas
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - aqueous-air mole fraction
!       NaCl mass - total NaCl brine mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
!
!---      Energy - temperature
!
          DNR(IEQT,N) = -1.D-7
!
!---      Water mass - aqueous pressure
!
          PLX = PL(2,N) + PATM
          DNR(IEQW,N) = MAX( 1.D-4,1.D-7*PLX )
!
!---      Air mass - aqueous-air mole fraction
!---      Air mass - gas air partial pressure
!
          IF( ISLC(37).EQ.0 ) THEN
!            XMLAX = MAX( PL(2,N)+PATM,PATM )/HCAW
!            IF( XMLA(2,N).GT.(1.D-2*XMLAX) ) THEN
!              DNR(IEQA,N) = SIGN( 1.D-4*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
!            ELSE
!              DNR(IEQA,N) = SIGN( 1.D-3*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
!            ENDIF
            DNR(IEQA,N) = 1.D-1
          ENDIF
!
!---      NaCl mass - total NaCl brine mass fraction
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_BRNS( T(2,N),PLX,XLSMX )
            DNR(IEQS,N) = 1.D-5*XLSMX
          ENDIF
!
!---    Unsaturated system w/ or w/o entrapped gas
!       
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
!
!---      Energy - temperature
!
          DNR(IEQT,N) = -1.D-7
!
!---      Water mass - aqueous pressure
!
          DNRX = MAX( 1.D-1,1.D-7*ABS(PG(2,N)-PL(2,N)) )
          DNR(IEQW,N) = -DNRX
!
!---      Air mass - gas pressure
!
          IF( ISLC(37).EQ.0 ) THEN
            DNRX = MAX( 1.D-3,1.D-7*(PG(2,N)+PATM) )
            DNR(IEQA,N) = DNRX
          ENDIF
!
!---      NaCl mass - total NaCl brine mass fraction
!
          IF( ISLC(32).EQ.0 ) THEN
            INDX = 0
            CALL REGION_4( T(2,N),PSWX,INDX )
            CALL SOL_BRNS( T(2,N),PSWX,XLSMX )
            DNR(IEQS,N) = 1.D-5*XLSMX
          ENDIF
!
!---    Saturated system w/ entrapped gas
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
!
!---      Energy - temperature
!
          DNR(IEQT,N) = -1.D-7
!
!---      Water mass - aqueous pressure
!
          PLX = PL(2,N)+PATM
          DNR(IEQW,N) = MAX( 1.D-1,1.D-6*PLX )
!
!---      Initial trapped gas saturation for the van Genuchten or
!         Brooks/Corey entrapment model  ---
!
          IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 .OR.
     &      ISCHR(IZN).EQ.201 .OR. ISCHR(IZN).EQ.202 ) THEN
!
!---        w/ Webb extension  ---
!
            IF( ISM(IZN).EQ.2 ) THEN
              ESGTMX = SCHR(15,IZN)
!
!---        w/o Webb extension  ---
!
            ELSE
              ESGTMX = SCHR(15,IZN)/(1.D+0-SCHR(4,IZN))
            ENDIF
          ELSE
            ESGTMX = 0.D+0
          ENDIF
!
!---      w/ Webb extension  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            SGTMX = ESGTMX
!
!---      w/o Webb extension  ---
!
          ELSE
            SGTMX = ESGTMX*(1.D+0-SCHR(4,IZN))
          ENDIF
!
!---      Air mass - trapped gas saturation
!
          IF( ISLC(37).EQ.0 ) THEN
            DNR(IEQA,N) = SIGN( 1.D-7,5.D-1*SGTMX-SG(2,N) )
          ENDIF
!
!---      NaCl mass - total NaCl brine mass fraction
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL SOL_BRNS( T(2,N),PLX,XLSMX )
            DNR(IEQS,N) = 1.D-5*XLSMX
          ENDIF
!
!---    Fully unsaturated conditions
!
!       Energy - temperature
!       Water mass - water vapor partial pressure
!       Air mass - gas pressure
!       NaCl mass - salt mass  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
!
!---      Energy - temperature  ---
!
          DNR(IEQT,N) = -1.D-7
!
!---      Water mass - water vapor partial pressure  ---
!
          DPX = SIGN( MAX( 1.D-4,1.D-6*PVW(2,N) ),1.D+0-PVW(2,N) )
          DNR(IEQW,N) = DPX
!
!---      Air mass - air gas partial pressure  ---
!
          IF( ISLC(37).EQ.0 ) THEN
            DPX = MAX( 1.D-4,1.D-6*PVA(2,N) )
            DNR(IEQA,N) = DPX
          ENDIF
!
!---      NaCl mass - salt mass  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            DNR(IEQS,N) = 1.D-6
          ENDIF
        ENDIF
!
!---    Assign gas-entry pressure for non Brooks-Corey;
!       Brooks-Corey; Brooks-Corey, Dual Porosity; and
!       Brooks-Corey, Entrapment  ---
!
        IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &      ISCHR(IZN).EQ.202 ) THEN
          ENPR = SCHR(1,IZN)*RHORL*GRAV
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          ENPR = MIN( SCHR(1,IZN),SCHR(5,IZN) )*RHORL*GRAV
        ELSE
          ENPR = 0.D+0
        ENDIF
!
!--- Increment the primary variables  ---
!
        DO M = 3,ISVC+2
          T(M,N) = T(2,N)
          PL(M,N) = PL(2,N)
          PG(M,N) = PG(2,N)
          PVW(M,N) = PVW(2,N)
          PVA(M,N) = PVA(2,N)
          XMLA(M,N) = XMLA(2,N)
          SG(M,N) = SG(2,N)
          SL(M,N) = SL(2,N)
          YLS(M,N) = YLS(2,N)
          TMS(M,N) = TMS(2,N)
!
!---      Saturated system w/o entrapped gas
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - aqueous-air mole fraction
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
            IF( M.EQ.IEQT+2 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
              PG(M,N) = PL(M,N) + ENPR - EPSLX
            ELSEIF( M.EQ.IEQA+2 .AND. ISLC(37).EQ.0 ) THEN
!              XMLA(M,N) = XMLA(M,N) + DNR(IEQA,N)
              PVA(M,N) = PVA(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Unsaturated system w/ or w/o entrapped gas
!       
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
            IF( M.EQ.IEQT+2 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQA+2 .AND. ISLC(37).EQ.0 ) THEN
              PG(M,N) = PG(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Saturated system w/ entrapped gas
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - trapped gas saturation
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
            IF( M.EQ.IEQT+2 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PL(M,N) = PL(M,N) + DNR(IEQW,N)
              IF( ISCHR(IZN).EQ.202 ) THEN
                PG(M,N) = PL(M,N) + SCHR(5,IZN) - EPSLX
              ELSE
                PG(M,N) = PL(M,N) + ENPR - EPSLX
              ENDIF
            ELSEIF( M.EQ.IEQA+2 .AND. ISLC(37).EQ.0 ) THEN
              SG(M,N) = SG(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - air partial pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
            IF( M.EQ.IEQT+2 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PVW(M,N) = PVW(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQA+2 .AND. ISLC(37).EQ.0 ) THEN
              PVA(M,N) = PVA(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
              TMS(M,N) = TMS(M,N) + DNR(IEQS,N)
            ENDIF
         ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCB_GT
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
!     Load Jacobian matrix for matrix, fractures, and boreholes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCB_GT'
!
!---  Load Jacobian matrix for the energy equation
!     (zero flux boundary), isothermal option  ---
!
      IF( ISLC(30).EQ.0 ) CALL JCBT_GT
!
!---  Load Jacobian matrix for the water equation
!     (zero flux boundary)  ---
!
      CALL JCBW_GT
!
!---  Load Jacobian matrix for the air equation
!     (zero flux boundary)  ---
!
      IF( ISLC(37).EQ.0 ) CALL JCBA_GT
!
!---  Load Jacobian matrix for the salt equation for field nodes
!     (zero flux boundary)  ---
!
      IF( ISLC(32).EQ.0 ) CALL JCBS_GT
!
!---  Modify the Jacobian matrix for boundary conditions  ---
!
      CALL BCJ_GT
!
!---  Modify the Jacobian matrix for surface cover nodes  --
!
      IF( NSFCN.GT.0 ) CALL JCB_SFC
!
!---  Fracture and borehole flow and transport solution  ---
!
      IF( ISLC(74).NE.0 ) THEN
!
!---    Load Jacobian matrix for the fracture energy equation, 
!       isothermal option  ---
!
        IF( ISLC(30).EQ.0 .AND. (ISLC(74).EQ.1 .OR. ISLC(74).EQ.3) )
     &    CALL JCBT_FRC_GT
!
!---    Load Jacobian matrix for the borehole energy equation, 
!       isothermal option  ---
!
        IF( ISLC(30).EQ.0 .AND. (ISLC(74).EQ.2 .OR. ISLC(74).EQ.3) )
     &    CALL JCBT_BH_GT
!
!---    Modify Jacobian matrix for the borehole and fracture energy 
!       equations for borehole to fracture connections, 
!       isothermal option  ---
!
        IF( ISLC(30).EQ.0 ) CALL JCBT_BF_GT
!
!---    Modify Jacobian matrix for matrix grid cells, 
!       fracture triangles, and borehole nodes for transfer of energy
!       between matrix grid cells and fracture triangles and borehole
!       nodes  ---
!
        IF( ISLC(30).EQ.0 ) CALL JCBT_MFB_GT
!
!---    Load Jacobian matrix for the fracture water equation  ---
!
        IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL JCBW_FRC_GT
!
!---    Load Jacobian matrix for the borehole water equation  ---
!
        IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL JCBW_BH_GT
!
!---    Modify Jacobian matrix for the borehole and fracture water 
!       equations for borehole to fracture connections  ---
!
        CALL JCBW_BF_GT
!
!---    Modify Jacobian matrix for matrix grid cells, 
!       fracture triangles, and borehole nodes for transfer of water
!       mass between matrix grid cells and fracture triangles and 
!       borehole nodes  ---
!
        IF( ISLC(38).EQ.0 ) CALL JCBW_MFB_GT
!
!---    Load Jacobian matrix for the fracture air equation  ---
!
        IF( ISLC(37).EQ.0 ) THEN
          IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL JCBA_FRC_GT
!
!---      Load Jacobian matrix for the borehole air equation  ---
!
          IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL JCBA_BH_GT
!
!---      Modify Jacobian matrix for the borehole and fracture air 
!         equations for borehole to fracture connections  ---
!
          CALL JCBA_BF_GT
!
!---      Modify Jacobian matrix for matrix grid cells, 
!         fracture triangles, and borehole nodes for transfer of air
!         mass between matrix grid cells and fracture triangles and 
!         borehole nodes  ---
!
          IF( ISLC(38).EQ.0 ) CALL JCBA_MFB_GT
        ENDIF
!
!---    Load Jacobian matrix for the fracture salt equation  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL JCBS_FRC_GT
!
!---      Load Jacobian matrix for the borehole salt equation  ---
!
          IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL JCBS_BH_GT
!
!---      Modify Jacobian matrix for the borehole and fracture salt 
!         equations for borehole to fracture connections, 
!         isothermal option  ---
!
          CALL JCBS_BF_GT
!
!---      Modify Jacobian matrix for matrix grid cells, 
!         fracture triangles, and borehole nodes for transfer of salt
!         mass between matrix grid cells and fracture triangles and 
!         borehole nodes  ---
!
          IF( ISLC(38).EQ.0 ) CALL JCBS_MFB_GT
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KSP_GT( N,PGX,PLX,BTGLX,SGX,SGTX,SLX,SLFX,SLMX,
     &  RKLX,RKGX,ASLX,ASLMINX,ESGTX,ESGTMX,SLRX,INDX )
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
!     and gas relative permeability.
!
!     INDX = 0 : Trapped-gas saturation computed.
!     INDX = 1 : Trapped-gas saturation given.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 February 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE GRID
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
      SUB_LOG(ISUB_LOG) = '/KSP_GT'
      IZN = IZ(N)
!
!---  Aqueous and gas saturation  ---
!
      IF( NPHAZ(2,N).EQ.4 ) THEN
        ASL_F = 0.D+0
        ASL_M = 0.D+0
        ASLX = 0.D+0
        ASLMINX = 0.D+0
        ESLX = 0.D+0
        ESGTX = 0.D+0
        SGTX = 0.D+0
        SLX = 0.D+0
        SGX = 1.D+0
      ELSE
        CALL SP_GT( ASL_F,ASL_M,ASLX,ASLMINX,ESLX,ESGTX,ESGTMX,
     &    PGX,PLX,BTGLX,SGX,SGTX,SLX,SLFX,SLMX,SLRX,INDX,IZN )
      ENDIF
!
!---  Aqueous relative permeability  ---
!
      CALL RKLS_GT( ASL_F,ASL_M,ESLX,PGX,PLX,RKLX,SLX,IZN )
!
!---  Aqueous relative permeability tensor  ---
!
      DO 100 ITX = 1,3
        IF( IRPLT(ITX,IZN).NE.0 ) CALL RKLST_GT( ASL_F,ASL_M,
     &    ESLX,PGX,PLX,RKLX(ITX),SLX,IZN,ITX )
  100 CONTINUE
!
!---  Gas relative permeability  ---
!
      RKGX = (RKLX(1)*RKLX(2)*RKLX(3))**(THIRD)
      CALL RKGS_GT( ASL_F,ASL_M,ASLX,ASLMINX,ESGTX,PGX,PLX,RKGX,
     &  SGX,SGTX,SLX,IZN )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KSP_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_GT
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
!     variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE HYST
      USE GRID
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVI
      USE FDVG
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LDO_GT'
!
!---  Assign old time step values  ---
!
      DO N = 1,NFLD
!
!---    Saturated system w/o entrapped gas
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - aqueous-air mole fraction
!       NaCl mass - total NaCl brine mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          ASLMIN(1,N) = 1.D+0
          ASLMIN(2,N) = 1.D+0
!
!---    Unsaturated system w/ or w/o entrapped gas
!       
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          ASLMIN(1,N) = MIN( ASL(N),ASLMIN(1,N),ASLMIN(2,N) )
          ASLMIN(1,N) = MAX( ASLMIN(1,N),0.D+0 )
          ASLMIN(2,N) = ASLMIN(1,N)
!
!---    Saturated system w/ entrapped gas
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
          ASLMIN(1,N) = ASLMIN(2,N)
!
!---    Fully unsaturated conditions
!
!       Energy - temperature
!       Water mass - water vapor partial pressure
!       Air mass - gas pressure
!       NaCl mass - salt mass  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
          ASLMIN(1,N) = 0.D+0
          ASLMIN(2,N) = 0.D+0
        ENDIF
        BTGL(1,N) = BTGL(2,N)
        T(1,N) = T(2,N)
        PL(1,N) = PL(2,N)
        PG(1,N) = PG(2,N)
        PORD(1,N) = PORD(2,N)
        PORT(1,N) = PORT(2,N)
        BTGL(1,N) = BTGL(2,N)
        SL(1,N) = SL(2,N)
        SG(1,N) = SG(2,N)
        SI(1,N) = SI(2,N)
        PVA(1,N) = PVA(2,N)
        PVW(1,N) = PVW(2,N)
        XGA(1,N) = XGA(2,N)
        XGW(1,N) = XGW(2,N)
        XMGA(1,N) = XMGA(2,N)
        XMGW(1,N) = XMGW(2,N)
        XLA(1,N) = XLA(2,N)
        XLW(1,N) = XLW(2,N)
        XMLA(1,N) = XMLA(2,N)
        XMLW(1,N) = XMLW(2,N)
        RHOL(1,N) = RHOL(2,N)
        RHOG(1,N) = RHOG(2,N)
        RHOSP(1,N) = RHOSP(2,N)
        RHOI(1,N) = RHOI(2,N)
        VISL(1,N) = VISL(2,N)
        VISG(1,N) = VISG(2,N)
        TORL(1,N) = TORL(2,N)
        TORG(1,N) = TORG(2,N)
        RKL(1,1,N) = RKL(1,2,N)
        RKL(2,1,N) = RKL(2,2,N)
        RKL(3,1,N) = RKL(3,2,N)
        RKG(1,N) = RKG(2,N)
        DFGW(1,N) = DFGW(2,N)
        DFLA(1,N) = DFLA(2,N)
        THKG(1,N) = THKG(2,N)
        THKL(1,N) = THKL(2,N)
        HI(1,N) = HI(2,N)
        HL(1,N) = HL(2,N)
        HGW(1,N) = HGW(2,N)
        HGA(1,N) = HGA(2,N)
        HG(1,N) = HG(2,N)
        UEG(1,N) = UEG(2,N)
        HSP(1,N) = HSP(2,N)
        NPHAZ(1,N) = NPHAZ(2,N)
        SS(1,N) = SS(2,N)
        TMS(1,N) = TMS(2,N)
        XLS(1,N) = XLS(2,N)
        YLS(1,N) = YLS(2,N)
        DO NSL = 1,NSOLU
          CO(N,NSL) = C(N,NSL)
        ENDDO

        DO NEQ = 1,NEQC+NEQK
          NSL = NEQ + NSOLU 
          CO(N,NSL) = C(N,NSL)
        ENDDO
        DO NSP = 1,NSPR
          SP_CO(N,NSP) = SP_C(N,NSP)
        ENDDO

      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_GT group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_GT
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
!     Compute hydrologic, thermodynamic and physical properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 24 February 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
      USE GRID
      USE FDVT
      USE FDVS
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
      SUB_LOG(ISUB_LOG) = '/PROP_GT'
!
!---  Loop over nodes computing secondary variables   ---
!
        DO N = 1,NFLD
        IZN = IZ(N)

        POR0(1,N) = POR0(1,N)
        POR0(2,N) = POR0(2,N)

!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
        N_DB = N
!
!---    Initial trapped gas saturation for the van Genuchten or
!       Brooks/Corey entrapment model  ---
!
        IF( ISCHR(IZN).EQ.101 .OR. ISCHR(IZN).EQ.102 .OR.
     &      ISCHR(IZN).EQ.201 .OR. ISCHR(IZN).EQ.202 ) THEN
!
!---      w/ Webb extension  ---
!
          IF( ISM(IZN).EQ.2 ) THEN
            ESGTMX = SCHR(15,IZN)
!
!---      w/o Webb extension  ---
!
          ELSE
            ESGTMX = SCHR(15,IZN)/(1.D+0-SCHR(4,IZN))
          ENDIF
        ELSE
          ESGTMX = 0.D+0
        ENDIF
!
!---    Minimum apparent aqueous saturation  ---
!
        ASLMINX = ASLMIN(2,N)
!
!---    Loop over increment indices  ---
!
        DO M = 2,ISVC+2
          PLX = PL(M,N)+PATM
          PGX = PG(M,N)+PATM
          INDX = 0
          CALL REGION_4( T(M,N),PSW(M,N),INDX )
!
!---      Saturated system w/o entrapped gas
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas air partial pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
            PX = PLX
            PWX = PLX
            PCX = 0.D+0
            SL(M,N) = 1.D+0
            CALL SCF_GL( BTGL(M,N),PWX )
            CALL SOL_BRNS( T(M,N),PWX,XLSMX )
            XLS(M,N) = MIN(YLS(M,N),XLSMX)
            CALL SP_B( T(M,N),XLS(M,N),PSBX )
            CALL P_IAPWS( T(M,N),PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL DENS_B( XLS(M,N),RHOBX,RHOLWX,T(M,N) )
            PVBX = PSBX
!            PVA(M,N) = XMLA(M,N)*HCAW
            XMLA(M,N) = PVA(M,N)/HCAW
            PVW(M,N) = PVBX
            XMGA(M,N) = PVA(M,N)/(PVA(M,N)+PVW(M,N))
            CALL EQUIL( XGA(M,N),XGW(M,N),XLA(M,N),XLS(M,N),XLW(M,N),
     &        XMGA(M,N),XMGW(M,N),XMLA(M,N),XMLS(M,N),XMLW(M,N) )
!
!---      Unsaturated system w/ or w/o entrapped gas
!       
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
            PX = PG(2,N) + PATM
            CALL P_IAPWS( T(M,N),PSW(M,N),RHOX,RHOLWX,HX,HLWX,UX,ULWX )
            PCX = PG(M,N)-PL(M,N)
            CALL VPL_B( T(M,N),PSW(M,N),PCX,RHOLWX,PVWX )
            CALL SCF_GL( BTGL(M,N),PVWX )
            CALL SOL_BRNS( T(M,N),PSW(M,N),XLSMX )
            XLS(M,N) = MIN(YLS(M,N),XLSMX)
            CALL DENS_B( XLS(M,N),RHOBX,RHOLWX,T(M,N) )
            CALL SP_B( T(M,N),XLS(M,N),PSBX )
            CALL VPL_B( T(M,N),PSBX,PCX,RHOBX,PVBX )
            CALL P_IAPWS( T(M,N),PVBX,RHOGWX,RHOX,HGWX,HX,UGWX,UX )
            PVA(M,N) = PGX-PVBX
            IF( PVA(M,N).LT.1.D-6 ) PVA(M,N) = 0.D+0
            PVW(M,N) = PVBX
            XMLA(M,N) = PVA(M,N)/HCAW
            XMGA(M,N) = PVA(M,N)/PGX
            CALL EQUIL( XGA(M,N),XGW(M,N),XLA(M,N),XLS(M,N),XLW(M,N),
     &        XMGA(M,N),XMGW(M,N),XMLA(M,N),XMLS(M,N),XMLW(M,N) )
!
!---      Saturated system w/ entrapped gas
!
!         Energy - temperature
!         Water mass - aqueous pressure
!         Air mass - trapped gas saturation
!         NaCl mass - total NaCl brine mass fraction  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
            PCX = 0.D+0
            PWX = PSW(M,N)
            CALL SCF_GL( BTGL(M,N),PWX )
            CALL P_IAPWS( T(M,N),PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL SOL_BRNS( T(M,N),PWX,XLSMX )
            XLS(M,N) = MIN(YLS(M,N),XLSMX)
            CALL SP_B( T(M,N),XLS(M,N),PSBX )
            PX = MAX( PLX,PSBX )
            CALL DENS_B( XLS(M,N),RHOBX,RHOLWX,T(M,N) )
            CALL VPL_B( T(M,N),PSBX,PCX,RHOBX,PVBX )
            PVA(M,N) = PX-PVBX
            IF( PVA(M,N).LT.1.D-6 ) PVA(M,N) = 0.D+0
            PVW(M,N) = PVBX
            XMLA(M,N) = PVA(M,N)/HCAW
            XMGA(M,N) = PVA(M,N)/(PVA(M,N)+PVW(M,N))
            CALL EQUIL( XGA(M,N),XGW(M,N),XLA(M,N),XLS(M,N),XLW(M,N),
     &        XMGA(M,N),XMGW(M,N),XMLA(M,N),XMLS(M,N),XMLW(M,N) )
!
!---      Fully unsaturated conditions
!
!         Energy - temperature
!         Water mass - water vapor partial pressure
!         Air mass - air gas partial pressure
!         NaCl mass - salt mass  ---
!
          ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
            PGX = PVA(M,N) + PVW(M,N)
            PG(M,N) = PGX - PATM
            PX = PGX
            INDX = 0
            CALL REGION_4( T(M,N),PSW(M,N),INDX )
            CALL SOL_BRNS( T(M,N),PSW(M,N),XLSMX )
            IF( TMS(M,N).GT.EPSL ) THEN
              XLS(M,N) = XLSMX
            ELSE
              XLS(M,N) = 0.D+0
            ENDIF
            CALL SCF_GL( BTGL(M,N),PVW(M,N) )
            SL(M,N) = 0.D+0
            SGT(M,N) = 0.D+0
            CALL CAP_GT( ASLMIN(2,N),BTGL(M,N),PCX,SL(M,N),SGT(M,N),IZN)
            PL(M,N) = PG(M,N) - PCX
            PLX = PGX - PCX
            CALL P_IAPWS( T(M,N),PVW(M,N),RHOGWX,RHOLWX,HGWX,HLWX,
     &        UGWX,ULWX )
            CALL DENS_B( XLS(M,N),RHOBX,RHOLWX,T(M,N) )
            XMLA(M,N) = PVA(M,N)/HCAW
            XMGA(M,N) = PVA(M,N)/PGX
            CALL EQUIL( XGA(M,N),XGW(M,N),XLA(M,N),XLS(M,N),XLW(M,N),
     &        XMGA(M,N),XMGW(M,N),XMLA(M,N),XMLS(M,N),XMLW(M,N) )
          ENDIF
!
!---      Porous-media porosity  ---
!
          CALL PORSTY_GT( N,PX,PCMP(N),PORD(M,N),PORT(M,N) )
          PORD(M,N) = MAX( PORD(M,N),EPSL )
          PORT(M,N) = MAX( PORT(M,N),PORD(M,N) )
!
!---      Saturation, relative permeability  ---
!
          INDX = 0
          IF( NPHAZ(2,N).EQ.3 ) THEN
            SGT(M,N) = SG(M,N)
            INDX = 1
          ENDIF
          CALL KSP_GT( N,PG(M,N),PL(M,N),BTGL(M,N),SG(M,N),SGT(M,N),
     &      SL(M,N),SDPF(N),SDPM(N),RKL(1,M,N),RKG(M,N),ASLX,ASLMINX,
     &      ESGTX,ESGTMX,SLRX,INDX )
          IF( M.EQ.2 ) THEN
            ASL(N) = ASLX
            ASGT(N) = ESGTX
          ENDIF
!
!---      Gas density and component fractions  ---
!
          CALL AIRGSD( T(M,N),PVA(M,N),RHOGAX )
          RHOG(M,N) = XGA(M,N)*RHOGAX + XGW(M,N)*RHOGWX
          WTMGX = XMGA(M,N)*WTMA + XMGW(M,N)*WTMW
          RHOMG(M,N) = RHOG(M,N)/WTMGX
!
!---      Gas viscosity  ---
!
          CALL AIRGSV( T(M,N),VISGAX )
          CALL VISC_W( T(M,N),PVBX,RHOGWX,VISGWX )
          CALL VISC_G( VISGAX,VISGWX,XMGA(M,N),XMGW(M,N),VISG(M,N) )
!
!---      Aqueous density and molar density  ---
!
          RHOL(M,N) = RHOBX
          WTMLX = XMLA(M,N)*WTMA + XMLS(M,N)*WTMS + XMLW(M,N)*WTMW
          RHOML(M,N) = RHOL(M,N)/WTMLX
!
!---      Aqueous viscosity  ---
!
          CALL VISC_W( T(M,N),PX,RHOLWX,VISLWX )
          CALL VISC_B( T(M,N),XLS(M,N),VISLWX,VISBX )
          CALL VISC_L( XMLA(M,N),VISBX,VISGAX,VISL(M,N) )
!
!---      Gas-water diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFGW(M,N) = DFGWC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL BNDFAW( T(M,N),PGX,DFGW(M,N) )
          ELSEIF( ISLC(2).EQ.3 ) THEN
            CALL BNDFAW( T(M,N),PGX,DFGW(M,N) )
            CMFF = 1.D+0 + 2.6D+0/(DFGWC**0.5)
            AMC = PORD(M,N)*SL(M,N)
            ENHF = 9.5D+0 + 6.D+0*(AMC) -
     &        8.5D+0/EXP(MIN(((CMFF*AMC)**4),1.D+2))
            DFGW(M,N) = ENHF*DFGW(M,N)
          ELSEIF( ISLC(2).EQ.4 ) THEN
            CALL BNDFAW( T(M,N),PGX,DFGW(M,N) )
            ENHF = DFEF(1,IZ(N))+DFEF(2,IZ(N))*SL(M,N)-
     &        (DFEF(1,IZ(N))-DFEF(4,IZ(N)))
     &        *EXP(-((DFEF(3,IZ(N))*SL(M,N))**DFEF(5,IZ(N))))
            DFGW(M,N) = ENHF*DFGW(M,N)
          ENDIF
!
!---      Aqueous-air diffusion coefficients  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLA(M,N) = DFLAC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            CALL AIRDFL( T(M,N),VISL(M,N),DFLA(M,N) )
          ENDIF
!
!---      Aqueous-salt diffusion coefficient  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLS(M,N) = DFLSC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            CALL DIFC_LS( T(M,N),XLS(M,N),VISL(M,N),DFLS(M,N) )
          ENDIF
!
!---      Precipitated NaCl density and saturation  ---
!
          CALL DENS_S( T(M,N),PX,RHOSP(M,N) )
!
!---      Fully unsaturated system w/ or w/o entrapped gas
!         Water mass - aqueous saturation
!         air mass - gas pressure
!         NaCl mass - total NaCl brine mass fraction  ---
!
          IF( NPHAZ(2,N).EQ.4 ) THEN
            SS(M,N) = TMS(M,N)/(RHOSP(M,N)*PORD(M,N))
            SLX = EPSL
            YLS(M,N) = TMS(M,N)*RHOBX*SLX*PORD(M,N)
          ELSE
!
!---        Precipitated salt saturation  ---
!
            SS(M,N) = MAX(YLS(M,N)-XLSMX,0.D+0)*RHOBX*SL(M,N)/
     &        RHOSP(M,N)
!
!---        NaCl volumetric concentration  ---
!
            TMS(M,N) = YLS(M,N)*RHOBX*SL(M,N)*PORD(M,N)
          ENDIF

!
!---      Permeability reduction factor with simplified Verma and
!         Preuss model  ---
!
          IF( IPRF(IZN)/10.EQ.4 ) THEN
!
!---        Porosity with salt precipitation  ---
!
            PORD(M,N) = MAX( PORD(M,N)*(1.D+0-SS(M,N)),
     &        PORD(M,N)*PERM(5,IZN),1.D-12 )        
!
!---        Permeability reduction factor with salt precipitation 
!           and mineral precipitation or dissolution  ---
!
            PERMRF(M,N) = ((PORD(M,N)-PERM(5,IZN))/
     &        (PERM(4,IZN)-PERM(5,IZN)))**PERM(6,IZN)
!
!---      Permeability reduction factor with Verma and Preuss model  ---
!
          ELSEIF( IPRF(IZN).EQ.1 ) THEN
            CALL PERM_R( SS(M,N),PERMRF(M,N),PORD(M,N),IZN )
          ENDIF

!
!---      Aqueous and gas tortuosity  ---
!
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL(M,N),SG(M,N),ZERO,
     &      PORD(M,N),TORL(M,N),TORG(M,N),TORNX )
!
!---      Nonisothermal simulation  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Gas enthalpy and internal energy  ---
!
            CALL AIRGSH( T(M,N),PVA(M,N),HGA(M,N),UEGAX )
            UEG(M,N) = XGA(M,N)*UEGAX + XGW(M,N)*UGWX
            HGW(M,N) = HGWX
            HG(M,N) = XGA(M,N)*HGA(M,N) + XGW(M,N)*HGW(M,N)
!
!---        Gas thermal conductivity  ---
!
            CALL AIRGSK( T(M,N),THKGAX )
            CALL THK_W( T(M,N),PGX,RHOGWX,THKGWX )
            CALL THK_G( T(M,N),THKGAX,THKGWX,XMGA(M,N),
     &        XMGW(M,N),THKG(M,N) )
!
!---        Aqueous enthalpy and internal energy  ---
!
            CALL ENTH_B( T(M,N),XLS(M,N),HLWX,HBX )
            HL(M,N) = MAX(1.D+0-XLA(M,N),0.D+0)*HBX + XLA(M,N)*HGA(M,N)
!
!---        Aqueous thermal conductivity  ---
!
            CALL THK_W( T(M,N),PX,RHOLWX,THKLWX )
            CALL THK_B( T(M,N),XLS(M,N),THKLWX,THKL(M,N) )
!
!---        Precipitated NaCl enthalpy  ---
!
            CALL ENTH_S( T(M,N),HSP(M,N) )
          ENDIF
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDBC_GT
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
!     Read input file for boundary condition information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE GRID
      USE FILES
      USE CONST
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!

      CHARACTER*64 SDUM
      INTEGER IBCSPX(LSPBC+1)

      CHARACTER*64 ADUM,FDUM,GDUM
      CHARACTER*64 UNTS
      CHARACTER*512 CHDUM
      CHARACTER*64 BDUM(LUK+LSOLU+2)
      REAL*8 VAR(LBTM,LBCV)
      REAL*8 XVX(4),YVX(4),ZVX(4)
      INTEGER ITYP(LUK+LSOLU+2)
      LOGICAL FLG_CHK
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDBC_GT'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Boundary Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      NBC = 0
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Boundary Condition Cards'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO 400 NB = 1, NLIN
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---  Read boundary orientation  ---
!
        VARB = 'Boundary Condition Orientation'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'west').NE.0) THEN
          IBCDX = -1
          WRITE(IWR,'(A)') 'X-Direction: West Surface'
        ELSEIF( INDEX(ADUM(1:),'east').NE.0) THEN
          IBCDX = 1
          WRITE(IWR,'(A)') 'X-Direction: East Surface'
        ELSEIF( INDEX(ADUM(1:),'south').NE.0) THEN
          IBCDX = -2
          WRITE(IWR,'(A)') 'Y-Direction: South Surface'
        ELSEIF( INDEX(ADUM(1:),'north').NE.0) THEN
          IBCDX = 2
          WRITE(IWR,'(A)') 'Y-Direction: North Surface'
        ELSEIF( INDEX(ADUM(1:),'bottom').NE.0) THEN
          IBCDX = -3
          WRITE(IWR,'(A)') 'Z-Direction: Bottom Surface'
        ELSEIF( INDEX(ADUM(1:),'top').NE.0) THEN
          IBCDX = 3
          WRITE(IWR,'(A)') 'Z-Direction: Top Surface'
        ELSEIF( INDEX(ADUM(1:),'file').NE.0 ) THEN
          CALL RDCHR(ISTART,ICOMMA,NCHF,CHDUM,FDUM)
          NCHF = INDEX(FDUM,'  ')-1
!
!---      Check that boundary condition domain file exists  ---
!
          INQUIRE( FILE=FDUM(1:NCHF), FORM=GDUM, EXIST=FLG_CHK )
          IF( .NOT.FLG_CHK ) THEN
            INDX = 4
            CHMSG = 'Boundary Condition Domain File: '
     &        // FDUM(1:NCHF)
            CALL WRMSGS( INDX )
          ELSEIF( GDUM.EQ.'UNFORMATTED' ) THEN
            INDX = 4
            CHMSG = 'Unformatted Boundary Condition Domain File: '
     &        // FDUM(1:NCHF)
            CALL WRMSGS( INDX )
          ENDIF
          OPEN(UNIT=26,FILE=FDUM(1:NCHF),STATUS='OLD',FORM='FORMATTED')
          WRITE(IWR,'(/,2A)') 'Boundary Condition Domain File: ',
     &      FDUM(1:NCHF)
          I1X = 1
          I2X = 1
          J1X = 1
          J2X = 1
          K1X = 1
          K2X = 0
    5     CONTINUE
          READ(26,*,END=10) IX,JX,KX,IBCDX
          K2X = K2X+1
          GOTO 5
   10     CONTINUE
          REWIND(26)
        ENDIF
!
!---    Write boundary condition type header ---
!
        WRITE(IWR,'(A)') 'Boundary Condition Type(s): '
        DO M = 1,4
          ITYP(M) = 0
        ENDDO
!
!---    Nonisothermal simulations option  ---
!
        IF( ISLC(30).EQ.0 ) THEN
!
!---      Read energy boundary condition type  ---
!
          VARB = 'Energy Boundary Condition Type'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(1))
          IF( INDEX(BDUM(1)(1:),'dirichlet').NE.0 ) THEN
            ITYP(1) = 1
            WRITE(IWR,'(2X,2A)') 'Energy: Dirichlet'
          ELSEIF( INDEX(BDUM(1)(1:),'neumann').NE.0 ) THEN
            ITYP(1) = 2
            ITYP(2) = 3
            ITFX = 3
            WRITE(IWR,'(2X,2A)') 'Energy: Neumann'
          ELSEIF( INDEX(BDUM(1)(1:),'zero flux').NE.0 ) THEN
            ITYP(1) = 3
            ITYP(2) = 3
            ITFX = 3
            WRITE(IWR,'(2X,2A)') 'Energy: Zero Flux'
          ELSEIF( INDEX(BDUM(1)(1:),'initial cond').NE.0 ) THEN
            ITYP(1) = 4
            WRITE(IWR,'(2X,2A)') 'Energy: Initial Condition'
          ELSEIF( INDEX(BDUM(1)(1:),'inflow').NE.0 ) THEN
            ITYP(1) = 5
            WRITE(IWR,'(2X,2A)') 'Energy: Inflow'
          ELSEIF( INDEX(BDUM(1)(1:),'outflow').NE.0 ) THEN
            ITYP(1) = 6
            WRITE(IWR,'(2X,2A)') 'Energy: Outflow'
          ELSEIF( INDEX(BDUM(1)(1:),'advective').NE.0 ) THEN
            ITYP(1) = 18
            WRITE(IWR,'(2X,2A)') 'Energy: Advective Only'
          ELSEIF( INDEX(BDUM(1)(1:),'barometric').NE.0  .AND.
     &      INDEX(BDUM(1)(1:),'formula').NE.0 ) THEN
            ITYP(1) = 12
            WRITE(IWR,'(2X,2A)') 'Energy: Barometric Formula'
          ELSEIF( INDEX(BDUM(1)(1:),'geothermal').NE.0 ) THEN
            ITYP(1) = 8
            WRITE(IWR,'(2X,2A)') 'Energy: Geothermal Gradient'
          ELSEIF( INDEX(BDUM(1)(1:),'ground').NE.0 ) THEN
            ITYP(1) = 29
            WRITE(IWR,'(2X,2A)') 'Energy: Convective Ground Surface'
          ELSEIF( INDEX(BDUM(1)(1:),'conv').NE.0 .AND.
     &      INDEX(BDUM(1)(1:),'rad').NE.0 ) THEN
            ITYP(1) = 28
            WRITE(IWR,'(2X,2A)') 'Energy: Convective-Radiative'
          ELSEIF( INDEX(BDUM(1)(1:),'hydrostatic').NE.0 ) THEN
            ITYP(1) = 11
            ITYP(2) = 11
            ITYP(3) = 11
            ITYP(4) = 11
            ITFX = ITYP(2)
            WRITE(IWR,'(2X,2A)') 'Energy: Convective-Radiative'
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Energy Boundary Condition Type: ' 
     &        // BDUM(1)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( ISLC(30).NE.0 ) THEN
          ITYP(1) = 0
          WRITE(IWR,'(2X,2A)') 'Energy: Isothermal Option'
        ENDIF
!
!---    Skip fluid boundary condition reads for energy Neumann type
!       and hydrostatic boundaries  ---
!
        IF( ITYP(1).NE.2 .AND. ITYP(1).NE.3 .AND. ITYP(1).NE.11 ) THEN
!
!---      Read fluid flow boundary condition type  ---
!
          VARB = 'Fluid Flow Boundary Condition Type'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(2))
          IF( INDEX(BDUM(2)(1:),'dirichlet').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'inflow').NE.0 ) THEN
            ITYP(2) = 5
            WRITE(IWR,'(2X,2A)') 'Fluid Flow: Dirichlet Inflow'
          ELSEIF( INDEX(BDUM(2)(1:),'dirichlet').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'outflow').NE.0 ) THEN
            ITYP(2) = 6
            WRITE(IWR,'(2X,2A)') 'Fluid Flow: Dirichlet Outflow'
          ELSEIF( INDEX(BDUM(2)(1:),'dirichlet').NE.0 ) THEN
            ITYP(2) = 1
            WRITE(IWR,'(2X,2A)') 'Fluid Flow: Dirichlet'
          ELSEIF( INDEX(BDUM(2)(1:),'neumann').NE.0 .AND.
     &     INDEX(BDUM(2)(1:),'inflow').NE.0 ) THEN
            ITYP(2) = 7
            WRITE(IWR,'(2X,2A)') 'Fluid Flow: Neumann Inflow'
          ELSEIF( INDEX(BDUM(2)(1:),'neumann').NE.0 .AND.
     &     INDEX(BDUM(2)(1:),'outflow').NE.0 ) THEN
            ITYP(2) = 9
            WRITE(IWR,'(2X,2A)') 'Fluid Flow: Neumann Outflow'
          ELSEIF( INDEX(BDUM(2)(1:),'aqueous').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'neumann').NE.0  ) THEN
            ITYP(2) = 12
            WRITE(IWR,'(2X,2A)') 'Fluid Flow: Aqueous Neumann'
          ELSEIF( INDEX(BDUM(2)(1:),'neumann').NE.0 ) THEN
            ITYP(2) = 2
            WRITE(IWR,'(2X,2A)') 'Fluid Flow: Neumann'
          ELSEIF( INDEX(BDUM(2)(1:),'zero flux').NE.0 ) THEN
            ITYP(2) = 3
            WRITE(IWR,'(2X,2A)') 'Fluid Flow: Zero Flux'
          ELSEIF( INDEX(BDUM(2)(1:),'initial cond').NE.0 ) THEN
            ITYP(2) = 4
            WRITE(IWR,'(2X,2A)') 'Fluid Flow: Initial Condition'
          ELSEIF( INDEX(BDUM(2)(1:),'aqueous').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'seepage').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'face').NE.0 ) THEN
            IF( ITYP(1).NE.6 ) THEN
              INDX = 4
              CHMSG = 'Incompatible Fluid Flow: Aqueous Seepage Face'
     &          // ' and Energy Boundary: ' // BDUM(1)
              CALL WRMSGS( INDX )
            ENDIF
            ITYP(2) = 13
            WRITE(IWR,'(2X,2A)') 'Fluid Flow: Aqueous Seepage Face'
          ELSEIF( INDEX(BDUM(2)(1:),'poten').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'evap').NE.0 ) THEN
            ITYP(2) = 24
          ELSEIF( INDEX(BDUM(2)(1:),'evaporative').NE.0 ) THEN
            ITYP(2) = 25
            ITYP(3) = 0
            ITYP(4) = 0
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Fluid Flow Boundary Condition Type: '
     &        // BDUM(2)
            CALL WRMSGS( INDX )
          ENDIF
          ITFX = ITYP(2)
!
!---      Skip boundary condition state read for zero flux, or 
!         initial condition, or
!         evaporative fluid flow boundary condition type  ---
!
          IF( ITYP(2).NE.3 .AND. ITYP(2).NE.4 .AND. ITYP(2).NE.25 ) THEN
!
!---        Read boundary condition state  ---
!
            VARB = 'Boundary Condition State'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
            IF( INDEX(BDUM(3)(1:),'bc1').NE.0 ) THEN
              IF( ITYP(2).EQ.12 .OR. ITYP(2).EQ.13 ) THEN
                INDX = 4
                CHMSG = 'Incompatible Fluid Flow Boundary and ' // 
     &            'Boundary Condition State: ' // BDUM(3)
                CALL WRMSGS( INDX )
              ENDIF
              ITYP(2) = ITYP(2) + 100
              WRITE(IWR,'(2X,2A)') 'Boundary Condition State: #1'
!
!---          Read aqueous air option  ---
!
              VARB = 'Aqueous Air Option'
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(4))
              IF( INDEX(BDUM(4)(1:),'air').NE.0 .AND.
     &          INDEX(BDUM(4)(1:),'rel').NE.0 .AND.
     &          INDEX(BDUM(4)(1:),'sat').NE.0 ) THEN
                ITYP(3) = 2
                WRITE(IWR,'(2X,2A)') 'Aqueous Air: Relative Saturation'
              ELSEIF( INDEX(BDUM(4)(1:),'air').NE.0 .AND.
     &          INDEX(BDUM(4)(1:),'mass').NE.0 .AND.
     &          INDEX(BDUM(4)(1:),'frac').NE.0 ) THEN
                ITYP(3) = 3
                WRITE(IWR,'(2X,2A)') 'Aqueous Air: Mass Fraction'
              ELSE
                INDX = 4
                CHMSG = 'Unrecognized Aqueous Air Option: '
     &            // BDUM(4)
                CALL WRMSGS( INDX )
              ENDIF
!
!---          Read aqueous salt option  ---
!
              VARB = 'Aqueous Salt Option'
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(5))
              IF( INDEX(BDUM(5)(1:),'salt').NE.0 .AND.
     &          INDEX(BDUM(5)(1:),'rel').NE.0 .AND.
     &          INDEX(BDUM(5)(1:),'sat').NE.0 ) THEN
                ITYP(4) = 2
                WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Relative Saturation'
              ELSEIF( INDEX(BDUM(5)(1:),'salt').NE.0 .AND.
     &          INDEX(BDUM(5)(1:),'mass').NE.0 .AND.
     &          INDEX(BDUM(5)(1:),'frac').NE.0 ) THEN
                ITYP(4) = 3
                WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Mass Fraction'
              ELSE
                INDX = 4
                CHMSG = 'Unrecognized Aqueous Salt Option: '
     &            // BDUM(5)
                CALL WRMSGS( INDX )
              ENDIF
            ELSEIF( INDEX(BDUM(3)(1:),'bc2').NE.0 ) THEN
              ITYP(2) = ITYP(2) + 200
              WRITE(IWR,'(2X,2A)') 'Boundary Condition State: #2'
!
!---          Read aqueous salt option  ---
!
              VARB = 'Aqueous Salt Option'
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(5))
              IF( INDEX(BDUM(5)(1:),'salt').NE.0 .AND.
     &          INDEX(BDUM(5)(1:),'rel').NE.0 .AND.
     &          INDEX(BDUM(5)(1:),'sat').NE.0 ) THEN
                ITYP(4) = 2
                WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Relative Saturation'
              ELSEIF( INDEX(BDUM(5)(1:),'salt').NE.0 .AND.
     &          INDEX(BDUM(5)(1:),'mass').NE.0 .AND.
     &          INDEX(BDUM(5)(1:),'frac').NE.0 ) THEN
                ITYP(4) = 3
                WRITE(IWR,'(2X,2A)') 'Aqueous Salt: Mass Fraction'
              ELSE
                INDX = 4
                CHMSG = 'Unrecognized Aqueous Salt Option: '
     &            // BDUM(5)
                CALL WRMSGS( INDX )
              ENDIF
            ELSEIF( INDEX(BDUM(3)(1:),'bc3').NE.0 ) THEN
              IF( ITYP(2).EQ.12 .OR. ITYP(2).EQ.13 ) THEN
                INDX = 4
                CHMSG = 'Incompatible Fluid Flow Boundary and ' // 
     &            'Boundary Condition State: ' // BDUM(3)
                CALL WRMSGS( INDX )
              ENDIF
              ITYP(2) = ITYP(2) + 300
              WRITE(IWR,'(2X,2A)') 'Boundary Condition State: #3'
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Boundary Condition State: ' // 
     &          BDUM(3)
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
        ENDIF
!
!---    Read solute boundary condition type(s) ---
!
        IF( IEQC.GT.0 ) THEN
          DO 12 NSL = 1,NSOLU
!
!---        Skip reads for stationary solutes  ---
!
            IF( IEDL(NSL).EQ.4 ) CYCLE
            VARB = 'Solute Boundary Condition Type'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(NSL+LUK))
            IF( INDEX(BDUM(NSL+LUK)(1:),'volum').NE.0 .AND.
     &        INDEX(BDUM(NSL+LUK)(1:),'conc').NE.0 ) THEN
              ITYP(NSL+LUK) = 1
              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
              WRITE(IWR,'(2X,A)') 'Solute (' // SOLUT(NSL)(1:IDB) // 
     &          ': Volumetric Concentration'
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'outflow').NE.0 ) THEN
              ITYP(NSL+LUK) = 7
              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
              WRITE(IWR,'(2X,A)') 'Solute (' // SOLUT(NSL)(1:IDB) // 
     &          ': Outflow'
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'zero').NE.0 .AND.
     &        INDEX(BDUM(NSL+LUK)(1:),'flux').NE.0 ) THEN
              ITYP(NSL+LUK) = 3
              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
              WRITE(IWR,'(2X,A)') 'Solute (' // SOLUT(NSL)(1:IDB) // 
     &          ': Zero Flux'
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'init').NE.0 .AND.
     &        INDEX(BDUM(NSL+LUK)(1:),'cond').NE.0 ) THEN
              ITYP(NSL+LUK) = 12
              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
              WRITE(IWR,'(2X,A)') 'Solute (' // SOLUT(NSL)(1:IDB) // 
     &          ': Initial Condition'
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'aqu').NE.0 .AND.
     &        INDEX(BDUM(NSL+LUK)(1:),'conc').NE.0 ) THEN
              ITYP(NSL+LUK) = 8
              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
              WRITE(IWR,'(2X,A)') 'Solute (' // SOLUT(NSL)(1:IDB) // 
     &          ': Aqueous Concentration'
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'gas').NE.0 .AND.
     &        INDEX(BDUM(NSL+LUK)(1:),'conc').NE.0 ) THEN
              ITYP(NSL+LUK) = 9
              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
              WRITE(IWR,'(2X,A)') 'Solute (' // SOLUT(NSL)(1:IDB) // 
     &          ': Gas Concentration'
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'aqu').NE.0 .AND.
     &        INDEX(BDUM(NSL+LUK)(1:),'inflow').NE.0 ) THEN
              ITYP(NSL+LUK) = 14
              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
              WRITE(IWR,'(2X,A)') 'Solute (' // SOLUT(NSL)(1:IDB) // 
     &          ': Aqueous Inflow'
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'gas').NE.0 .AND.
     &        INDEX(BDUM(NSL+LUK)(1:),'inflow').NE.0 ) THEN
              ITYP(NSL+LUK) = 15
              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
              WRITE(IWR,'(2X,A)') 'Solute (' // SOLUT(NSL)(1:IDB) // 
     &          ': Gas Inflow'
!            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'aqu').NE.0 .AND.
!     &        INDEX(BDUM(NSL+LUK)(1:),'advec').NE.0 ) THEN
!              ITYP(NSL+LUK) = 11
!              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
!              WRITE(IWR,'(2X,A)') 'Solute (' // SOLUT(NSL)(1:IDB) // 
!     &          ': Aqueous Advection'
!            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'gas').NE.0 .AND.
!     &        INDEX(BDUM(NSL+LUK)(1:),'advec').NE.0 ) THEN
!              ITYP(NSL+LUK) = 12
!              IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
!              WRITE(IWR,'(2X,A)') 'Solute (' // SOLUT(NSL)(1:IDB) // 
!     &          ': Gas Advection'
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Solute Boundary Condition: ' //
     &          BDUM(NSL+LUK)
              CALL WRMSGS( INDX )
            ENDIF
   12     CONTINUE
        ENDIF

        IF( ISLC(40).EQ.1 ) THEN
!
!---    Aqueous species boundary condition types, 
!       allowing for returns in input lines  ---
!
        CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
        IF( INDX.EQ.0 ) THEN
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
        ENDIF
        BDUM(NSOLU+LUK+1) = 'zero flux'
        IDFLT = 1
        VARB = 'Aqueous Reactive Species Boundary Condition Type: '
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(NSOLU+LUK+1))
!
!---    Gas species boundary condition types, 
!       allowing for returns in input lines  ---
!
        CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
        IF( INDX.EQ.0 ) THEN
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
        ENDIF
        BDUM(NSOLU+LUK+2) = 'zero flux'
        IDFLT = 1
        VARB = 'Gas Reactive Species Boundary Condition Type: '
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(NSOLU+LUK+2))
!
!---    Number of reactive species  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Number of Reactive Species'
        CALL RDINT(ISTART,ICOMMA,CHDUM,IBCSPX(1))
        DO 16 NSPX = 2,IBCSPX(1)+1
          IBCSPX(NSPX) = 0
   16   CONTINUE
!
!---    Loop over number of reactive species  ---
!
        DO 21 NSPX = 1,IBCSPX(1)
!
!---      Allow for returns in input lines  ---
!
          CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
          IF( INDX.EQ.0 ) THEN
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
          ENDIF
          VARB = 'Species Name'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,SDUM)
!
!---      Aqueous species  ---
!
          DO 17 M = 1,NSPL
            IF( SPNML(M).EQ.SDUM ) THEN
              IBCSPX(NSPX+1) = M
              GOTO 19
            ENDIF
   17     CONTINUE
!
!---      Gas species  ---
!
          DO 18 M = 1,NSPG
            IF( SPNMG(M).EQ.SDUM ) THEN
              MX = M + NSPL + NSPS
              GOTO 19
            ENDIF
   18     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Aqueous or Gas Species Name: '
     &       // SDUM(1:NCH)
          CALL WRMSGS( INDX )
   19     CONTINUE
   21   CONTINUE
        ENDIF

!
!---    Read and write boundary domain indices  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        IF( INDEX(ADUM(1:),'file').EQ.0 ) THEN
          VARB = 'Boundary Condition Domain: '
          CALL RDINT(ISTART,ICOMMA,CHDUM,I1X)
          CALL RDINT(ISTART,ICOMMA,CHDUM,I2X)
          CALL RDINT(ISTART,ICOMMA,CHDUM,J1X)
          CALL RDINT(ISTART,ICOMMA,CHDUM,J2X)
          CALL RDINT(ISTART,ICOMMA,CHDUM,K1X)
          CALL RDINT(ISTART,ICOMMA,CHDUM,K2X)
          WRITE(IWR,'(A)') VARB(1:IVR)
          WRITE(IWR, '(2X,A,I6,A,I6)') 'I = ',I1X,' to ',I2X
          WRITE(IWR, '(2X,A,I6,A,I6)') 'J = ',J1X,' to ',J2X
          WRITE(IWR, '(2X,A,I6,A,I6)') 'K = ',K1X,' to ',K2X
!
!---      Check boundary domain  ---
!
          IF( I1X.GT.I2X .OR. J1X.GT.J2X .OR. K1X.GT.K2X ) THEN
            INDX = 4
            CHMSG = 'Nonascending Boundary Condition Domain Indices'
            CALL WRMSGS( INDX )
          ENDIF
          IF( I1X.LT.1 .OR. I2X.GT.IFLD. OR. J1X.LT.1 .OR.
     &      J2X.GT.JFLD .OR. K1X.LT.1 .OR. K2X.GT.KFLD ) THEN
            INDX = 4
            CHMSG = 'Illegal Boundary Condition Domain'
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---    Read number of boundary times  ---
!
        VARB = 'Number of Boundary Condition Times'
        CALL RDINT(ISTART,ICOMMA,CHDUM,IBCMX)
        IF( IBCMX.LE.-3 ) THEN
          IBCCX = 1
          IBCMX = -IBCMX
          WRITE(IWR,'(A)') 'Cyclic Boundary Conditions'
        ELSEIF( IBCMX.GE.1 ) THEN
          IBCCX = 0
          WRITE(IWR,'(A)') 'Noncyclic Boundary Conditions'
        ELSEIF( IBCMX.EQ.0 ) THEN
          INDX = 4
          CHMSG = 'No Boundary Condition Times'
          CALL WRMSGS( INDX )
        ELSE
          INDX = 4
          CHMSG = 'Number of Cyclic Boundary Conditions Times < 3'
          CALL WRMSGS( INDX )
        ENDIF
        IF( IBCMX.GT.LBTM ) THEN
          INDX = 5
          CHMSG = 'Number of Boundary Condition Times > LBTM'
          CALL WRMSGS( INDX )
        ENDIF
        BCTMO = -SMALL
        WRITE(IWR,'(A)') 'Boundary Condition Times and Variables:'
        DO 100 NTM = 1,IBCMX
          DO 20 M = 1,LBCV
            VAR(NTM,M) = 0.D+0
   20     CONTINUE
!
!---      Read, write, and convert boundary condition time, variables,
!         and units  ---
!
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
          VARB = 'Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH)
     &      ,': ',VAR(NTM,1)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR(NTM,1),INDX)
!
!---      Nonisothermal simulations option  ---
!
          IF( ISLC(30).EQ.0 ) THEN
!
!---        Read temperature, C for Dirichlet and inflow type 
!           energy boundary conditions  ---
!
            IF( ITYP(1).EQ.1 .OR. ITYP(1).EQ.5 .OR. ITYP(1).EQ.7 ) THEN
              VARB = 'Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
!
!---        Read energy flux, W/m^2 for Neumann type 
!           energy boundary conditions  ---
!
            ELSEIF( ITYP(1).EQ.2 ) THEN
              VARB = 'Energy Flux'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNKG = 1
              IUNS = -3
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
!
!---        Read reference temperature, C, reference Z point, m, and
!           geothermal gradient, C/m for geothermal gradient type 
!           energy boundary conditions  ---
!
            ELSEIF( ITYP(1).EQ.8 ) THEN
              VARB = 'Reference Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
              VARB = 'Reference Z Point'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,6)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
              VARB = 'Geothermal Gradient'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,7)
              INDX = 0
              IUNK = 1
              IUNM = -1
              CALL RDUNIT(UNTS,VAR(NTM,7),INDX)
!
!---        Hydrostatic boundary conditions  ---
!
            ELSEIF( ITYP(1).EQ.11 ) THEN
              VARB = 'Reference Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
              VAR(NTM,2) = VAR(NTM,2) - PATM
              VARB = 'Reference Pressure Z Point'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
              VARB = 'Reference Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,4)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
              VARB = 'Reference Temperature Z Point'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,5)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,VAR(NTM,5),INDX)
              VARB = 'Geothermal Gradient'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,6)
              INDX = 0
              IUNK = 1
              IUNM = -1
              CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
              VARB = 'Reference Salt Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
              WRITE(IWR,'(2X,2A)') VARB(1:IVR),': ',VAR(NTM,7)
              VARB = 'Reference Salinity Z Point'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,8)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,VAR(NTM,8),INDX)
              VARB = 'Geosalinity Gradient'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,9))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,9)
              INDX = 0
              IUNM = -1
              CALL RDUNIT(UNTS,VAR(NTM,9),INDX)
!
!---        Read reference temperature, C, and reference Z point,
!           for barometric formula energy boundary conditions  ---
!
            ELSEIF( ITYP(1).EQ.12 ) THEN
              VARB = 'Reference Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
              VARB = 'Reference Z Point'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,6)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
              VARB = 'Thermal Lapse Rate'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,7)
              INDX = 0
              IUNK = 1
              IUNM = -1
              CALL RDUNIT(UNTS,VAR(NTM,7),INDX)
              VARB = 'Reference Gas Pressure'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,4)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
              VAR(NTM,4) = VAR(NTM,4) - PATM
            ELSEIF( ITYP(1).EQ.18 ) THEN
              VARB = 'Convective Temperature'
              ISX = ISTART
              ICX = ICOMMA
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),
     &          ': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
              VARB = 'Convective Heat Transfer Coefficient'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),
     &          ': ',VAR(NTM,6)
              INDX = 0
              IUNK = -1
              IUNS = -3
              IUNKG = 1
              CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
            ELSEIF( ITYP(1).EQ.28 ) THEN
              VARB = 'Convective Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),
     &          ': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
              VARB = 'Radiative Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),
     &          ': ',VAR(NTM,6)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
            ELSEIF( ITYP(1).EQ.29 ) THEN
              VARB = 'Air Temperature'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),
     &          ': ',VAR(NTM,2)
              INDX = 0
              IUNK = 1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            ENDIF
          ENDIF
!
!---      Read pressure, Pa for Dirichlet, Dirichlet Inflow, 
!         Dirichlet Outflow, or Hydraulic Gradient
!         type fluid flow boundary conditions  ---
!
          IF( ITFX.EQ.1 .OR. ITFX.EQ.5 .OR. ITFX.EQ.6 .OR.
     &      ITFX.EQ.8 ) THEN
            VARB = 'Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VAR(NTM,3) = VAR(NTM,3) - PATM
!
!---      Read volumetric flux, m/s for Neumann, Neumann Inflow,
!         and Neumann Outflow type fluid flow boundary conditions  ---
!
          ELSEIF( ITFX.EQ.2 .OR. ITFX.EQ.7 .OR. ITFX.EQ.9 ) THEN
            VARB = 'Volumetric Flux'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
!
!---      Read potential evaporation rate and maximum capillary head
!         for the potential evaporation type fluid flow boundary 
!         conditions  ---
!
          ELSEIF( ITFX.EQ.2 .OR. ITFX.EQ.7 .OR. ITFX.EQ.9 ) THEN
            VARB = 'Potential Evaporation Rate'
            IF( IRD.EQ.21 ) WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( IRD.EQ.21 ) WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),
     &        ': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &        VAR(NTM,3),', m/s)'
            VARB = 'Maximum Capillary Head'
            IF( IRD.EQ.21 ) WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( IRD.EQ.21 ) WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),
     &        ': ',VAR(NTM,7)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,7),INDX)
            IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &        VAR(NTM,7),', m)'
!
!---      Read pressure, wind speed and water vapor relative humidity
!         for the evaporative type fluid flow boundary condition  ---
!
          ELSEIF( ITFX.EQ.25 ) THEN
            VARB = 'Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VAR(NTM,3) = VAR(NTM,3) - PATM
            VARB = 'Wind Speed'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IF( IRD.EQ.21 ) WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            IF( IRD.EQ.21 ) WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &        VAR(NTM,4),', m/s)'
            VARB = 'Water Vapor Relative Humidity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,5)
!
!---      Read aqueous volumetric flux, m/s for aqueous Neumann
!         boundary conditions  ---
!
          ELSEIF( ITFX.EQ.12 ) THEN
            VARB = 'Aqueous Volumetric Flux'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
!
!---      Read aqueous pressure, Pa and reference elevation, m for 
!         seepage face boundary conditions  ---
!
          ELSEIF( ITFX.EQ.13 ) THEN
            VARB = 'Reference Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Reference Z Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,6)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
            VARB = 'Thermal Lapse Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,7)
            INDX = 0
            IUNK = 1
            IUNM = -1
            CALL RDUNIT(UNTS,VAR(NTM,7),INDX)
            VARB = 'Reference Gas Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VAR(NTM,4) = VAR(NTM,4) - PATM
            VARB = 'Reference Aqueous Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VAR(NTM,3) = VAR(NTM,3) - PATM
            VARB = 'Reference Z Point'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,8)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,8),INDX)
          ENDIF
!
!---      Boundary condition state #1  ---
!
!         SL = 1.0
!         SG = 0.0
!
!         Declared variables:
!
!         BC(4,NTM,NB) - aqueous air relative saturation, or
!                        aqueous air mass fraction
!         BC(5,NTM,NB) - aqueous salt relative saturation, or
!                        aqueous salt mass fraction
!
          IF( ITYP(2)/100.EQ.1 ) THEN
!
!---        Read aqueous air  ---
!
            IF( ITYP(3).EQ.2 ) THEN
              VARB = 'Aqueous Air Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,4)
            ELSEIF( ITYP(3).EQ.3 ) THEN
              VARB = 'Aqueous Air Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,4)
            ENDIF
!
!---        Read aqueous salt  ---
!
            IF( ITYP(4).EQ.2 ) THEN
              VARB = 'Aqueous Salt Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,5)
            ELSEIF( ITYP(4).EQ.3 ) THEN
              VARB = 'Aqueous Salt Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,5)
            ENDIF
!
!---      Boundary condition state #2  ---
!
!         SL < 1.0
!         SG > 0.0
!
!         Declared variables:
!
!         BC(4,NTM,NB) - aqueous fraction
!         BC(5,NTM,NB) - aqueous salt relative saturation, or
!                        aqueous salt mass fraction
!
          ELSEIF( ITYP(2)/100.EQ.2 ) THEN
!
!---        Read aqueous fraction  ---
!
            IF( ITFX.NE.12 .AND. ITFX.NE.13 ) THEN
              VARB = 'Aqueous Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,4)
            ENDIF
!
!---        Read aqueous salt  ---
!
            IF( ITYP(4).EQ.2 ) THEN
              VARB = 'Aqueous Salt Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,5)
            ELSEIF( ITYP(4).EQ.3 ) THEN
              VARB = 'Aqueous Salt Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,5)
            ENDIF
!
!---      Boundary condition state #3  ---
!
!         SL = 0.0
!         SG = 1.0
!
!         Declared variables:
!
!         BC(4,NTM,NB) - water vapor relative saturation, or
!                        relative humidity
!         BC(8,NTM,NB) - reference elevation
!         BC(9,NTM,NB) - relative humidity gradient
!        
          ELSEIF( ITYP(2)/100.EQ.3 ) THEN
            VARB = 'Water Vapor Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,4)
            ITYP(3) = 2           
          ENDIF
!
!---      Read solute transport boundary condition variables ---
!
          IF( IEQC.GT.0 ) THEN
           JGC = LBCU
           DO 60 NSL = 1,NSOLU
              IF( ITYP(NSL+LUK).EQ.1 ) THEN
                VARB = 'Volumetric Concentration'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+JGC))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
                WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL)(1:IDB),
     &            VARB(1:IVR),', ',UNTS(1:NCH),': ',VAR(NTM,NSL+JGC)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VAR(NTM,NSL+JGC),INDX)
              ELSEIF( ITYP(NSL+LUK).EQ.8 .OR. ITYP(NSL+LUK).EQ.14 ) THEN
                VARB = 'Aqueous Volumetric Concentration'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+JGC))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
                WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL)(1:IDB),
     &            VARB(1:IVR),', ',UNTS(1:NCH),': ',VAR(NTM,NSL+JGC)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VAR(NTM,NSL+JGC),INDX)
              ELSEIF( ITYP(NSL+LUK).EQ.9 .OR. ITYP(NSL+LUK).EQ.15 ) THEN
                VARB = 'Gas Volumetric Concentration'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+JGC))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
                WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL)(1:IDB),
     &            VARB(1:IVR),', ',UNTS(1:NCH),': ',VAR(NTM,NSL+JGC)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VAR(NTM,NSL+JGC),INDX)
              ENDIF
   60       CONTINUE
          ENDIF

          IF( ISLC(40).EQ.1 ) THEN
!
!---        Loop over reactive species inputs  ---
!
            DO 62 NSPX = 1,IBCSPX(1)
              NSP = IBCSPX(NSPX+1)
              M = NSOLU+LBCU+NSPX
!
!---          Initial input line  ---
!
              IF( NSPX.EQ.1 ) THEN
                CALL RDINPL( CHDUM )
                CALL LCASE( CHDUM )
                ISTART = 1
              ENDIF
!
!---          Allow for returns in input lines  ---
!
              CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
              IF( INDX.EQ.0 ) THEN
                CALL RDINPL( CHDUM )
                CALL LCASE( CHDUM )
                ISTART = 1
              ENDIF
!
!---          Aqueous species  ---
!
              IF( NSP.LE.NSPL ) THEN
                IF( ITYP(NSOLU+LUK+1).EQ.8
     &            .OR. ITYP(NSOLU+LUK+1).EQ.14
     &            .OR. ITYP(NSOLU+LUK+1).EQ.23 ) THEN
                  VARB = 'Aqueous Concentration, '
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,M))
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  IDB = INDEX( SPNML(NSP)(1:),'  ') - 1
                  WRITE(IWR,'(2X,A,2X,3A,1PE11.4,$)') 
     &              SPNML(NSP)(1:IDB),VARB(1:IVR),UNTS(1:NCH),
     &              ': ',VAR(NTM,M)
                  INDX = 0
                  IUNM = -3
                  IUNMOL = 1
                  CALL RDUNIT(UNTS,VAR(NTM,M),INDX)
!
!---              Convert aqueous concentration from kmol/m^3 to
!                 mol/m^3  ---
!
                  VAR(NTM,M) = VAR(NTM,M)*1.D+3
                  WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &              VAR(NTM,M),', mol/m^3)'
                ELSEIF( ITYP(NSOLU+LUK+1).EQ.12 ) THEN
                  VARB = 'Dummy Variable, '
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                ELSE
                  VARB = 'Dummy Variable, '
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                ENDIF
!
!---          Gas species  ---
!
              ELSEIF( NSP.GT.(NSPL+NSPS) .AND. 
     &          NSP.LE.(NSPL+NSPS+NSPG) ) THEN
                IF( ITYP(NSOLU+LUK+2).EQ.8
     &            .OR. ITYP(NSOLU+LUK+2).EQ.9
     &            .OR. ITYP(NSOLU+LUK+2).EQ.15
     &            .OR. ITYP(NSOLU+LUK+2).EQ.43 ) THEN
                  VARB = 'Gas Concentration, '
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,M))
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  IDB = INDEX( SPNML(NSP)(1:),'  ') - 1
                  WRITE(IWR,'(2X,A,2X,3A,1PE11.4,$)') 
     &              SPNML(NSP)(1:IDB),VARB(1:IVR),UNTS(1:NCH),
     &              ': ',VAR(NTM,M)
                  INDX = 0
                  IUNM = -3
                  IUNMOL = 1
                  CALL RDUNIT(UNTS,VAR(NTM,M),INDX)
!
!---              Convert gas concentration from kmol/m^3 to
!                 mol/m^3  ---
!
                  VAR(NTM,M) = VAR(NTM,M)*1.D+3
                  WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &              VAR(NTM,M),', mol/m^3)'
                ELSEIF( ITYP(NSOLU+LUK+1).EQ.12 ) THEN
                  VARB = 'Dummy Variable, '
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                ELSE
                  VARB = 'Dummy Variable, '
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                ENDIF
              ENDIF
   62       CONTINUE
          ENDIF

!
!---      Check for nonascending boundary condition times  ---
!
          IF( VAR(NTM,1).LT.BCTMO ) THEN
            INDX = 4
            CHMSG = 'Boundary Condition Time Sequencing'
            CALL WRMSGS( INDX )
          ENDIF
          BCTMO = VAR(NTM,1)
  100   CONTINUE
!
!---    Assign values to boundary variables  ---
!
        JGC = LBCU
        DO 108 NTM = 1,IBCMX
          DO 102 M = 1,JGC
            BC(M,NTM,NB) = VAR(NTM,M)
  102     CONTINUE
          IF( IEQC.GT.0 ) THEN
            DO 104 NSL = 1,NSOLU
              BC(NSL+JGC,NTM,NB) = VAR(NTM,NSL+JGC)
  104       CONTINUE
          ENDIF

          IF( ISLC(40).EQ.1 ) THEN
            DO 106 NSPX = 1,IBCSPX(1)
              M = NSOLU+LBCU+NSPX
              BC(M,NTM,NB) = VAR(NTM,M)
  106       CONTINUE
          ENDIF

  108   CONTINUE
!
!---  Assign values to boundary variables  ---
!
        NBCL = 0
        DO 320 K = K1X, K2X
          DO 310 J = J1X, J2X
            DO 300 I = I1X, I2X
              IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
                READ(26,*,END=320) IX,JX,KX,IBCDX
                N = ND(IX,JX,KX)
              ELSE
                N = ND(I,J,K)
                IX = I
                JX = J
                KX = K
              ENDIF
!
!---         Check for boundary applied to inactive nodes  ---
!
              IF( IXP(N).EQ.0 ) THEN
                WRITE(IWR,'(A,I9)') 'Boundary Condition Applied ' //
     &            'to an Inactive Node: ',N
                CYCLE
              ENDIF
!
!---          Check for boundary applied to interior surfaces  ---
!
              IERR = 0
              IF( IBCDX.EQ.-3 .AND. KX.NE.1) THEN
                IF( IXP(N-IJFLD).NE.0 .AND. INBS(1,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'Bottom Boundary'
                  WRITE(IWR,'(A)') 'Bottom Boundary'
                ENDIF
              ELSEIF( IBCDX.EQ.-2 .AND. JX.NE.1) THEN
                IF( IXP(N-IFLD).NE.0 .AND. INBS(2,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'South Boundary'
                  WRITE(IWR,'(A)') 'South Boundary'
                ENDIF
              ELSEIF( IBCDX.EQ.-1 .AND. IX.NE.1) THEN
                IF( IXP(N-1).NE.0 .AND. INBS(3,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'West Boundary'
                  WRITE(IWR,'(A)') 'West Boundary'
                ENDIF
              ELSEIF( IBCDX.EQ.1 .AND. IX.NE.IFLD) THEN
                IF( IXP(N+1).NE.0 .AND. INBS(4,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'East Boundary'
                  WRITE(IWR,'(A)') 'East Boundary'
                ENDIF
              ELSEIF( IBCDX.EQ.2 .AND. JX.NE.JFLD) THEN
                IF( IXP(N+IFLD).NE.0 .AND. INBS(5,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'North Boundary'
                  WRITE(IWR,'(A)') 'North Boundary'
                ENDIF
              ELSEIF( IBCDX.EQ.3 .AND. KX.NE.KFLD) THEN
                IF( IXP(N+IJFLD).NE.0 .AND. INBS(6,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'Top Boundary'
                  WRITE(IWR,'(A)') 'Top Boundary'
                ENDIF
              ENDIF
!
!---          Report boundary error  ---
!
              IF( IERR.EQ.1 ) THEN
                WRITE(ISC,'(A,I9)') 'Node = ',N
                WRITE(IWR,'(A,I9)') 'Node = ',N
                WRITE(ISC,'(3(A,I9))') 'I = ',I,' J = ',J,' K = ',K
                WRITE(IWR,'(3(A,I9))') 'I = ',I,' J = ',J,' K = ',K
                INDX = 7
                IMSG = NBC
                CHMSG = 'Boundary Cond. Applied to an Interior Surface'
     &            //': Boundary Number'
                CALL WRMSGS( INDX )
              ENDIF
              NBCL = NBCL + 1
              NBC = NBC + 1
              IF( NBC.GT.LBC ) THEN
                INDX = 5
                CHMSG = 'Number of Boundary Condition Surfaces > '
     &            //'Parameter LBC'
                CALL WRMSGS( INDX )
              ENDIF
              IBCN(NBC) = N
              IBCC(NBC) = IBCCX
              IBCD(NBC) = IBCDX
              IBCT(1,NBC) = ITYP(1)
              IBCT(2,NBC) = ITYP(2)
              IBCT(3,NBC) = ITYP(3)
              IBCT(4,NBC) = ITYP(4)
              IF( IEQC.GT.0 ) THEN
                DO 110 NSL = 1,NSOLU
                  IBCT(NSL+LUK,NBC) = ITYP(NSL+LUK)
  110           CONTINUE
              ENDIF

              IF( ISLC(40).EQ.1 ) THEN
                IBCT(NSOLU+LUK+1,NBC) = ITYP(NSOLU+LUK+1)
                IBCT(NSOLU+LUK+2,NBC) = ITYP(NSOLU+LUK+2)
                DO 120 NSP = 1,LSPBC+1
                  IBCSP(NSP,NBC) = IBCSPX(NSP)
  120           CONTINUE
              ENDIF

              IBCM(NBC) = IBCMX
              IBCIN(NBC) = NB
!
!---          Check for double boundary conditions  ---
!
              DO 220 M = 1,NBC-1
                MB = IBCIN(M)
                IF( IBCN(M).EQ.N .AND. IBCD(M).EQ.IBCDX ) THEN
                  IF( (VAR(1,1).GT.BC(1,1,MB) .AND.
     &              VAR(1,1).LT.BC(1,IBCM(M),MB)) .OR.
     &              (VAR(IBCMX,1).GT.BC(1,1,MB) .AND.
     &              VAR(IBCMX,1).LT.BC(1,IBCM(M),MB)) ) THEN
                      INDX = 4
                      CHMSG = 'Multiple Boundary Conditions'
                      CALL WRMSGS( INDX )
                  ENDIF
                ENDIF
  220         CONTINUE
!
!---          Boundary surface centroids  ---
!
              IERR = 0
!
!---          Bottom surface centroid  ---
!
              IF( IBCDX.EQ.-3 ) THEN
                XVX(1) = XE(1,N)
                XVX(2) = XE(2,N)
                XVX(3) = XE(4,N)
                XVX(4) = XE(3,N)
                YVX(1) = YE(1,N)
                YVX(2) = YE(2,N)
                YVX(3) = YE(4,N)
                YVX(4) = YE(3,N)
                ZVX(1) = ZE(1,N)
                ZVX(2) = ZE(2,N)
                ZVX(3) = ZE(4,N)
                ZVX(4) = ZE(3,N)
                NP = 4
                CALL PGCNTRD( NP,XVX,YVX,ZVX,
     &            XPBC(NBC),YPBC(NBC),ZPBC(NBC) )
!
!---          South surface centroid  ---
!
              ELSEIF( IBCDX.EQ.-2 ) THEN
                XVX(1) = XE(1,N)
                XVX(2) = XE(2,N)
                XVX(3) = XE(6,N)
                XVX(4) = XE(5,N)
                YVX(1) = YE(1,N)
                YVX(2) = YE(2,N)
                YVX(3) = YE(6,N)
                YVX(4) = YE(5,N)
                ZVX(1) = ZE(1,N)
                ZVX(2) = ZE(2,N)
                ZVX(3) = ZE(6,N)
                ZVX(4) = ZE(5,N)
                NP = 4
                CALL PGCNTRD( NP,XVX,YVX,ZVX,
     &            XPBC(NBC),YPBC(NBC),ZPBC(NBC) )
!
!---          West surface centroid  ---
!
              ELSEIF( IBCDX.EQ.-1 ) THEN
                XVX(1) = XE(1,N)
                XVX(2) = XE(3,N)
                XVX(3) = XE(7,N)
                XVX(4) = XE(5,N)
                YVX(1) = YE(1,N)
                YVX(2) = YE(3,N)
                YVX(3) = YE(7,N)
                YVX(4) = YE(5,N)
                ZVX(1) = ZE(1,N)
                ZVX(2) = ZE(3,N)
                ZVX(3) = ZE(7,N)
                ZVX(4) = ZE(5,N)
                NP = 4
                CALL PGCNTRD( NP,XVX,YVX,ZVX,
     &            XPBC(NBC),YPBC(NBC),ZPBC(NBC) )
!
!---          East surface centroid  ---
!
              ELSEIF( IBCDX.EQ.1 ) THEN
                XVX(1) = XE(2,N)
                XVX(2) = XE(4,N)
                XVX(3) = XE(8,N)
                XVX(4) = XE(6,N)
                YVX(1) = YE(2,N)
                YVX(2) = YE(4,N)
                YVX(3) = YE(8,N)
                YVX(4) = YE(6,N)
                ZVX(1) = ZE(2,N)
                ZVX(2) = ZE(4,N)
                ZVX(3) = ZE(8,N)
                ZVX(4) = ZE(6,N)
                NP = 4
                CALL PGCNTRD( NP,XVX,YVX,ZVX,
     &            XPBC(NBC),YPBC(NBC),ZPBC(NBC) )
!
!---          North surface centroid  ---
!
              ELSEIF( IBCDX.EQ.2 ) THEN
                XVX(1) = XE(3,N)
                XVX(2) = XE(4,N)
                XVX(3) = XE(8,N)
                XVX(4) = XE(7,N)
                YVX(1) = YE(3,N)
                YVX(2) = YE(4,N)
                YVX(3) = YE(8,N)
                YVX(4) = YE(7,N)
                ZVX(1) = ZE(3,N)
                ZVX(2) = ZE(4,N)
                ZVX(3) = ZE(8,N)
                ZVX(4) = ZE(7,N)
                NP = 4
                CALL PGCNTRD( NP,XVX,YVX,ZVX,
     &            XPBC(NBC),YPBC(NBC),ZPBC(NBC) )
!
!---          Top surface centroid  ---
!
              ELSEIF( IBCDX.EQ.3 ) THEN
                XVX(1) = XE(5,N)
                XVX(2) = XE(6,N)
                XVX(3) = XE(8,N)
                XVX(4) = XE(7,N)
                YVX(1) = YE(5,N)
                YVX(2) = YE(6,N)
                YVX(3) = YE(8,N)
                YVX(4) = YE(7,N)
                ZVX(1) = ZE(5,N)
                ZVX(2) = ZE(6,N)
                ZVX(3) = ZE(8,N)
                ZVX(4) = ZE(7,N)
                NP = 4
                CALL PGCNTRD( NP,XVX,YVX,ZVX,
     &            XPBC(NBC),YPBC(NBC),ZPBC(NBC) )
              ENDIF
!
!---          Set base surface centroids  ---
!
              IF( NBCL.EQ.1 ) THEN
                XSBC(NB) = XPBC(NBC)
                YSBC(NB) = YPBC(NBC)
                ZSBC(NB) = ZPBC(NBC)
              ENDIF
  300       CONTINUE
  310     CONTINUE
  320   CONTINUE
        IF( INDEX(ADUM(1:),'file').NE.0 ) CLOSE(UNIT=26)
  400 CONTINUE
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDBC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIBC_GT
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
!     Read input file for initial boundary conditions information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 November 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVG
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
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,BDUM,FDUM,FMDUM,UNTS
      CHARACTER*24 CHLB(3)
      CHARACTER*512 CHDUM
      CHARACTER*20 FORM5
      CHARACTER*38 FORM6
      CHARACTER*39 FORM7
      INTEGER IDOM(6)
      REAL*8 VAR(5)
      LOGICAL FCHK
!
!----------------------Data Statements---------------------------------!
!
      DATA CHLB /'X-Direction Gradient, ','Y-Direction Gradient, ',
     &           'Z-Direction Gradient, '/
      DATA FORM5 /'(10(1PE22.15,1X),I3)'/
      DATA FORM6 /'(10(1PE22.15,1X),I3,1X,1(1PE22.15,1X))'/
      DATA FORM7 /'(10(1PE22.15,1X),I3,1X,10(1PE22.15,1X))'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIBC_GT'
      IERR = 0
!
!---  Write card information to ouput file  ---
!
      CARD = 'Initial Boundary Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Initialize initial condition states for air and salt  ---
!
      DO 10 N = 1,NFLD
        ICAIR(N) = 0
        ICBRN(N) = 0
   10 CONTINUE
!
!---  Read saturation initial condition option ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Initial Saturation Option: '
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      IF( INDEX(ADUM(1:),'hydrostatic').NE.0 .AND.
     &  INDEX(ADUM(1:),'geothermal').NE.0 ) THEN
        ISIC = 4
        VARB = 'Reference Aqueous Pressure'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(1)
        INDX = 0
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        CALL RDUNIT(UNTS,VSLC(1),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(1),', Pa)'
        VSLC(1) = VSLC(1) - PATM
        IF( INDEX( CHDUM(1:),'file' ).NE.0 ) THEN
!
!---      Read salt mass fraction from an external file, read
!         x-,y-,and z-coordinates for reference pressure  ---
!
          VARB = 'Reference-Pressure X Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(6))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(6)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(6),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(6),', m)'
          VARB = 'Reference-Pressure Y Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(7))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(7)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(7),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(7),', m)'
          VARB = 'Reference-Pressure Z Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(8))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(8)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(8),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(8),', m)'
!
!---    Read salt mass fraction from reference value, z-location, 
!       and gradient, read z-coordinate for reference pressure  ---
!
        ELSE
          VARB = 'Reference-Pressure Z Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(2))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(2)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(2),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(2),', m)'
        ENDIF
        VARB = 'Reference Temperature'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(3))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(3)
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,VSLC(3),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(3),', C)'
        VARB = 'Reference-Temperature Z Location'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(4))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(4)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VSLC(4),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(4),', m)'
        VARB = 'Geothermal Gradient'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(5))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(5)
        INDX = 0
        IUNK = 1
        IUNM = -1
        CALL RDUNIT(UNTS,VSLC(5),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(5),', C/m)'
!
!---    Check for salt mass fraction read from an external file  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          ISTARTX = ISTART
          ICOMMAX = ICOMMA
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---      Read salt mass fraction from an external file  ---
!
          IF( INDEX( ADUM(1:),'file' ).NE.0 ) THEN
            ISIC = 5
!
!---        Check for external file  ---
!
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            NCH = INDEX(FDUM,'  ')-1
            INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
            IF( .NOT.FCHK ) THEN
              INDX = 4
              CHMSG = 'Missing Salt Mass Fraction File: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ELSEIF( FDUM.EQ.'unformatted' ) THEN
              INDX = 4
              CHMSG = 'Salt Mass Fraction File Format: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
            OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',FORM='FORMATTED')
            WRITE(IWR,'(/,2A)') 'Salt Mass Fraction File: ',FDUM(1:NCH)
            UNTS = 'null'
            ADDER = 0.D+0
            INDX = 2
            CALL RDINAS( YLS,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
!
!---      Read salt mass fraction from reference value, z-location, 
!         and gradient  ---
!
          ELSE
            ISTART = ISTARTX
            ICOMMA = ICOMMAX
            VARB = 'Reference Salt Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(6))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VSLC(6)
            VARB = 'Reference-Salinity (Salt Mass Fraction) Z Location'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(7))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(7)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VSLC(7),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(7),', m)'
            VARB = 'Geosalinity (Salt Mass Fraction) Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(8))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(8)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,VSLC(8),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(8),', 1/m)'
          ENDIF
        ENDIF
        IF( NSOLU+NEQC.GT.0 ) THEN
          GOTO 22
        ELSE
          GOTO 1002
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'hydrostatic').NE.0 ) THEN
        ISIC = 6
        VARB = 'Reference-Pressure I Index'
        CALL RDINT(ISTART,ICOMMA,CHDUM,IVARX)
        VSLC(6) = REAL(IVARX)
        VARB = 'Reference-Pressure J Index'
        CALL RDINT(ISTART,ICOMMA,CHDUM,JVARX)
        VSLC(7) = REAL(JVARX)
        VARB = 'Reference-Pressure K Index'
        CALL RDINT(ISTART,ICOMMA,CHDUM,KVARX)
        VSLC(8) = REAL(KVARX)
        VARB = 'Reference-Pressure'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(1)
        INDX = 0
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        CALL RDUNIT(UNTS,VSLC(1),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(1),', Pa)'
        VSLC(1) = VSLC(1) - PATM
        GOTO 22
      ENDIF      
      CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
      WRITE(IWR,'(/,A)') VARB(1:IVR)
      WRITE(IWR,'(2X,A)') ADUM(1:NCHA)
      WRITE(IWR,'(2X,A)') BDUM(1:NCHB)
      IF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
        IF( INDEX(BDUM(1:),'gas pressure').NE.0 ) THEN
          ISIC = 1
        ELSEIF( INDEX(BDUM(1:),'aqueous pressure').NE.0 ) THEN
          ISIC = 2
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &      ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'gas pressure').NE.0 ) THEN
        IF( INDEX(BDUM(1:),'aqueous saturation').NE.0 ) THEN
          ISIC = 1
        ELSEIF( INDEX(BDUM(1:),'aqueous pressure').NE.0 ) THEN
          ISIC = 3
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &      ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
        IF( INDEX(BDUM(1:),'aqueous saturation').NE.0 ) THEN
          ISIC = 2
        ELSEIF( INDEX(BDUM(1:),'gas pressure').NE.0 ) THEN
          ISIC = 3
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &      ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ENDIF
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &    ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Read initial conditions  ---
!
   22 WRITE(IWR,'(/,A)') 'Initial Condition Variable(s) and Domain(s)'
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Initial Condition Cards: '
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO 1000 NL = 1, NLIN
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Initial Condition Variable: '
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        NPZX = 1
        IF( INDEX( ADUM(1:),'aqueous pres' ).NE.0 ) THEN
          VARB = 'Initial Aqueous Pressure'
          IUNM = -1
          IUNKG = 1
          IUNS = -2
        ELSEIF( INDEX( ADUM(1:),'gas pres' ).NE.0 ) THEN
          VARB = 'Initial Gas Pressure'
          IUNM = -1
          IUNKG = 1
          IUNS = -2
        ELSEIF( INDEX( ADUM(1:),'temperature' ).NE.0 ) THEN
          VARB = 'Initial Temperature'
          IUNK = 1
        ELSEIF( INDEX( ADUM(1:),'aqueous sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous Saturation'
        ELSEIF( INDEX( ADUM(1:),'relative' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'trapped gas' ).NE.0 ) THEN
          VARB = 'Initial Relative Trapped Gas Saturation'
        ELSEIF( INDEX( ADUM(1:),'trapped gas' ).NE.0 ) THEN
          VARB = 'Initial Trapped Gas Saturation'
        ELSEIF( INDEX( ADUM(1:),'air' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-Air Relative Saturation'
        ELSEIF( INDEX( ADUM(1:),'air' ).NE.0 .AND.
     &    ( INDEX( ADUM(1:),'mass' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'frac' ).NE.0 ) ) THEN
          VARB = 'Initial Aqueous-Air Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'air' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'conc' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-Air Concentration'
          IUNKG = 1
          IUNM = -3
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-Salt Relative Saturation'
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    ( INDEX( ADUM(1:),'mass' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'frac' ).NE.0 ) ) THEN
          VARB = 'Initial Aqueous-Salt Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'conc' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-Salt Concentration'
          IUNKG = 1
          IUNM = -3
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'molality' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-Salt Molality'
          IUNMOL = 1
          IUNKG = -1
        ELSEIF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          VARB = 'Initial Solute Concentration'
          IVAR = 1
          IUNM = -3

        ELSEIF( INDEX( ADUM(1:),'specie' ).NE.0 ) THEN
          VARB = 'Reactive Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          VARB = 'Initial Reactive Species Concentration'
!
!---      Set species units  ---
!
          IUNMOL = 1
          IF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
            IVAR = 4
            IUNM = -3
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'molal').NE.0 ) THEN
              IVAR = 3
              IUNKG = -1
            ELSE
              IVAR = 2
              IUNM = -3
            ENDIF
          ELSE
            IVAR = 1
            IUNM = -3
          ENDIF

        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(1))
        CALL RDCHR(ISTART,ICOMMA,NCHU,CHDUM,UNTS)
!
!---  Read initial conditions input from an external file  ---
!
        IF( INDEX( ADUM(1:),'file' ).NE.0 ) THEN
          IF( INDEX( ADUM(1:),'binary' ).NE.0 ) THEN
            WRITE(IWR,'(2X,3A)') ADUM(1:NCHA),',',UNTS(1:NCHU)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            NCH = INDEX(FDUM,'  ')-1
!
!---        Check for external file  ---
!
            INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
            IF( .NOT.FCHK ) THEN
              INDX = 4
              CHMSG = 'Missing Initial Conditions File: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ELSEIF( FDUM.EQ.'formatted' ) THEN
              INDX = 4
              CHMSG = 'Initial Conditions File Format: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
            OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',
     &        FORM='UNFORMATTED')
            WRITE(IWR,'(/,2A)') 'Initial Conditions File: ',FDUM(1:NCH)
          ELSEIF( INDEX( ADUM(1:),'ascii' ).NE.0 ) THEN
            WRITE(IWR,'(2X,3A)') ADUM(1:NCHA),',',UNTS(1:NCHU)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            NCH = INDEX(FDUM,'  ')-1
!
!---        Check for external file  ---
!
            INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
            IF( .NOT.FCHK ) THEN
              INDX = 4
              CHMSG = 'Missing Initial Conditions File: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ELSEIF( FDUM.EQ.'unformatted' ) THEN
              INDX = 4
              CHMSG = 'Initial Conditions File Format: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
            OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',
     &        FORM='FORMATTED')
            WRITE(IWR,'(/,2A)') 'Initial Conditions File: ',FDUM(1:NCH)
          ELSE
            WRITE(IWR,'(2X,4A,1PE11.4)') ADUM(1:NCHA),
     &        ' (Default Value), ',UNTS(1:NCHU),': ',VAR(1)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            NCH = INDEX(FDUM,'  ')-1
!
!---        Check for external file  ---
!
            INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
            IF( .NOT.FCHK ) THEN
              INDX = 4
              CHMSG = 'Missing Initial Conditions File: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ELSEIF( FDUM.EQ.'unformatted' ) THEN
              INDX = 4
              CHMSG = 'Initial Conditions File Format: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
            OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',FORM='FORMATTED')
            WRITE(IWR,'(/,2A)') 'Initial Conditions File: ',FDUM(1:NCH)
            INDX = 0
            CALL RDUNIT( UNTS,VAR(1),INDX )
          ENDIF
!
!---  Read initial conditions according to rock/soil zonations  ---
!
        ELSEIF( INDEX( ADUM(1:),'rock' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'zonation' ).NE.0 ) THEN
          VARB = 'Rock/Soil Name'
          CALL RDCHR(ISTART,ICOMMA,NCHF,CHDUM,FDUM)
!
!---  Search known rock types for a matching type ---
!
          DO 30 M = 1, NROCK
            IF( FDUM .EQ. ROCK(M)) THEN
            IROCK = M
            GOTO 40
          ENDIF
   30     CONTINUE
          INDX = 2
          CHMSG = 'Unrecognized Rock/Soil Type: '//FDUM
          CALL WRMSGS( INDX )
          GOTO 1000
   40     CONTINUE
          WRITE(IWR,'(2X,3A,1PE11.4,2A)') ADUM(1:NCHA),UNTS(1:NCHU),
     &      ': ',VAR(1),' Rock/Soil Type: ',FDUM(1:NCHF)
          INDX = 0
          CALL RDUNIT( UNTS,VAR(1),INDX )
!
!---    Read initial conditions input from the input file  ---
!
        ELSE
          WRITE(IWR,'(2X,4A,1PE11.4)') ADUM(1:NCHA),', ',
     &      UNTS(1:NCHU),': ',VAR(1)
          INDX = 0
          CALL RDUNIT( UNTS,VAR(1),INDX )
          INDX = 2
          VAR(5) = 1.D+0
          NCH = INDEX( UNTS,'  ' ) - 1
          IF( UNTS(1:NCH).EQ.'f' .OR. UNTS(1:NCH).EQ.'r' ) THEN
            VAR(5) = VAR(5)/1.8D+0
          ELSEIF( UNTS(1:NCH).EQ.'c' .OR. UNTS(1:NCH).EQ.'k' ) THEN
            VAR(5) = 1.D+0
          ELSE
            CALL RDUNIT( UNTS,VAR(5),INDX )
          ENDIF
          VARB = 'Initial Condition Variable Gradient: '
          DO 100 I = 2,4
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(I))
            VAR(I) = VAR(I)*VAR(5)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') CHLB(I-1),', ',UNTS(1:NCH),
     &        ': ',VAR(I)
            INDX = 0
            IUNM = -1
            CALL RDUNIT( UNTS,VAR(I),INDX )
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(I),', 1/m)'
  100     CONTINUE
!
!---      Read domain indices  ---
!
          VARB = 'Initial Condition Domain Index: '
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(1))
          IF( IDOM(1).LT.1 .OR. IDOM(1).GT.IFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Lower I-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(1)
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(2))
          IF( IDOM(2).LT.1 .OR. IDOM(2).GT.IFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Upper I-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(2)
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(3))
          IF( IDOM(3).LT.1 .OR. IDOM(3).GT.JFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Lower J-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(3)
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(4))
          IF( IDOM(4).LT.1 .OR. IDOM(4).GT.JFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Upper J-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(4)
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(5))
          IF( IDOM(5).LT.1 .OR. IDOM(5).GT.KFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Lower K-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(5)
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(6))
          IF( IDOM(6).LT.1 .OR. IDOM(6).GT.KFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Upper K-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(6)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Read variables  ---
!
        IF( INDEX(ADUM(1:),'temp').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNK = 1
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINBR( T,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINBS( T,ADDER,UNTS,INDX )
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINAR( T,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINAS( T,ADDER,UNTS,INDX )
              ENDIF
            ELSE
              IF( NPZX.EQ.0 ) THEN
                CALL RDINFR( T,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINFS( T,VAR,ADDER,UNTS,INDX )
              ENDIF    
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( T,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( T,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( T,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( T,VAR,ADDER,IDOM,INDX )
            ENDIF    
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous press').NE.0 ) THEN
          ADDER = -PATM
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINBR( PL,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINBS( PL,ADDER,UNTS,INDX )
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINAR( PL,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINAS( PL,ADDER,UNTS,INDX )
              ENDIF
            ELSE
              IF( NPZX.EQ.0 ) THEN
                CALL RDINFR( PL,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINFS( PL,VAR,ADDER,UNTS,INDX )
              ENDIF    
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( PL,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( PL,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( PL,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( PL,VAR,ADDER,IDOM,INDX )
            ENDIF    
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas pres').NE.0 ) THEN
          ADDER = -PATM
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINBR( PG,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINBS( PG,ADDER,UNTS,INDX )
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINAR( PG,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINAS( PG,ADDER,UNTS,INDX )
              ENDIF
            ELSE
              IF( NPZX.EQ.0 ) THEN
                CALL RDINFR( PG,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINFS( PG,VAR,ADDER,UNTS,INDX )
              ENDIF    
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( PG,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( PG,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( PG,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( PG,VAR,ADDER,IDOM,INDX )
            ENDIF    
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous sat').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( SL,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( SL,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( SL,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( SL,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( SL,VAR,ADDER,IDOM,INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'rel').NE.0 .AND.
     &    INDEX(ADUM(1:),'trapped gas').NE.0 ) THEN
          ADDER = 1.D+2
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINBR( SGT,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINBS( SGT,ADDER,UNTS,INDX )
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINAR( SGT,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINAS( SGT,ADDER,UNTS,INDX )
              ENDIF
            ELSE
              IF( NPZX.EQ.0 ) THEN
                CALL RDINFR( SGT,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINFS( SGT,VAR,ADDER,UNTS,INDX )
              ENDIF    
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( SGT,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( SGT,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( SGT,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( SGT,VAR,ADDER,IDOM,INDX )
            ENDIF    
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'trapped').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINBR( SGT,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINBS( SGT,ADDER,UNTS,INDX )
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINAR( SGT,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINAS( SGT,ADDER,UNTS,INDX )
              ENDIF
            ELSE
              IF( NPZX.EQ.0 ) THEN
                CALL RDINFR( SGT,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINFS( SGT,VAR,ADDER,UNTS,INDX )
              ENDIF    
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( SGT,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( SGT,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( SGT,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( SGT,VAR,ADDER,IDOM,INDX )
            ENDIF    
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0 ) THEN
          IVAR = 3
          IF( INDEX(ADUM(1:),'molality').NE.0 ) THEN
            IVAR = 4
          ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .OR.
     &      INDEX(ADUM(1:),'frac').NE.0 ) THEN
            IVAR = 3
          ELSEIF( INDEX(ADUM(1:),'rel').NE.0 .OR.
     &      INDEX(ADUM(1:),'sat').NE.0 ) THEN
            IVAR = 2
          ELSEIF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &      INDEX(ADUM(1:),'conc').NE.0 ) THEN
            IVAR = 1
          ENDIF
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( IVAR.EQ.0 .OR. IVAR.EQ.1 ) THEN
              IUNM = -3
              IUKG = 1
            ELSEIF( IVAR.EQ.4 ) THEN
              IUNMOL = 1
              IUNKG = -1
            ENDIF
            IF( NPZX.EQ.0 ) THEN
              CALL RDINFR( TMS,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
            ELSE
              CALL RDINFS( TMS,VAR,ADDER,UNTS,INDX )
            ENDIF    
            CLOSE(UNIT=26)
            DO 206 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 206
              ICBRN(N) = IVAR
  206       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( TMS,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( TMS,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
            DO 208 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 208
              IF( IZ(N).EQ.IROCK ) ICBRN(N) = IVAR
  208       CONTINUE
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( TMS,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( TMS,VAR,ADDER,IDOM,INDX )
            ENDIF    
            DO K = IDOM(5),IDOM(6)
            DO J = IDOM(3),IDOM(4)
            DO I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) CYCLE
              ICBRN(N) = IVAR
            ENDDO
            ENDDO
            ENDDO
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'air').NE.0 ) THEN
          IVAR = 3
          IF( INDEX(ADUM(1:),'mass').NE.0 .OR.
     &      INDEX(ADUM(1:),'frac').NE.0 ) THEN
            IVAR = 3
            ADDER = 0.D+0
          ELSEIF( INDEX(ADUM(1:),'rel').NE.0 .OR.
     &      INDEX(ADUM(1:),'sat').NE.0 ) THEN
            IVAR = 2
            ADDER = 0.D+0
          ELSEIF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &      INDEX(ADUM(1:),'conc').NE.0 ) THEN
            IVAR = 1
            ADDER = 0.D+0
          ENDIF
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( IVAR.EQ.0 .OR. IVAR.EQ.1 ) THEN
              IUNM = -3
              IUKG = 1
            ELSEIF( IVAR.EQ.4 ) THEN
              IUNM = -1
              IUNKG = 1
              IUNS = -2
            ENDIF
            IF( NPZX.EQ.0 ) THEN
              CALL RDINFR( PVA,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
            ELSE
              CALL RDINFS( PVA,VAR,ADDER,UNTS,INDX )
            ENDIF    
            CLOSE(UNIT=26)
            DO 306 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 306
              ICAIR(N) = IVAR
  306       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( PVA,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( PVA,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
            DO 308 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 308
              IF( IZ(N).EQ.IROCK ) ICAIR(N) = IVAR
  308       CONTINUE
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( PVA,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( PVA,VAR,ADDER,IDOM,INDX )
            ENDIF    
            DO K = IDOM(5),IDOM(6)
            DO J = IDOM(3),IDOM(4)
            DO I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) CYCLE
              ICAIR(N) = IVAR
            ENDDO
            ENDDO
            ENDDO
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
            IVAR = 3
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
            IVAR = 2
          ELSE
            IVAR = 1
          ENDIF
          IF( INDEX( UNTS(1:),'bd' ).NE.0 ) IVAR = -IVAR
          DO 420 NSL = 1,NSOLU
            IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SOLUT(NSL)(1:IDB) ) THEN
              ADDER = 0.D+0
              IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
                IUNM = -3
                IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
                  CALL RDINBP( C(1,NSL),ADDER,ICT(1,NSL),IVAR,UNTS )
                ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
                  CALL RDINAP( C(1,NSL),ADDER,ICT(1,NSL),IVAR,UNTS )
                ELSE
                  CALL RDINFP( C(1,NSL),VAR,ADDER,ICT(1,NSL),IVAR,UNTS )
                ENDIF
                CLOSE(UNIT=26)
              ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &          INDEX(ADUM(1:),'zonation').NE.0 )  THEN
                CALL RDINZP( C(1,NSL),VAR(1),ADDER,ICT(1,NSL),
     &            IVAR,IROCK )
              ELSE
                CALL RDINIP( C(1,NSL),VAR,ADDER,ICT(1,NSL),IVAR,IDOM )
              ENDIF
              GOTO 430
            ENDIF
  420     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute: ' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
  430     CONTINUE

        ELSEIF( INDEX(ADUM(1:),'specie').NE.0 ) THEN
          ADDER = 0.D+0
!
!---      Conservation- or kinetic-component species  ---
!
          IF( INDEX( BDUM(1:),'total_' ).NE.0 ) THEN
            DO 500 NSLX = NSOLU+1,NSOLU+NEQC+NEQK
              IDB = INDEX(SOLUT(NSLX)(1:),'  ') - 1
              IF( BDUM(1:NCHB).EQ.SOLUT(NSLX) ) THEN
                NSL = NSLX
                GOTO 540
              ENDIF
  500       CONTINUE
          ENDIF
!
!---      Aqueous reactive species  ---
!
          DO 510 NSPX = 1,NSPL
            IDB = INDEX(SPNML(NSPX)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SPNML(NSPX)(1:IDB) ) THEN
              NSP = NSPX
              GOTO 540
            ENDIF
  510     CONTINUE
!
!---      Solid reactive species  ---
!
          DO 520 NSPX = 1,NSPS
            IDB = INDEX(SPNMS(NSPX)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SPNMS(NSPX)(1:IDB) ) THEN
              NSP = NSPX + NSPL
!
!---          Verify that solid-species is not a mineral  ---
!
              IF( ISP_MN(NSP).EQ.1 ) THEN
                INDX = 4
                CHMSG = 'Solid-Species Mineral ' // 
     &            '(see Lithology Card): ' // BDUM(1:NCHB)
                CALL WRMSGS( INDX )
              ENDIF
              GOTO 540
            ENDIF
  520     CONTINUE
!
!---      Gas reactive species  ---
!
          DO 530 NSPX = 1,NSPG
            IDB = INDEX(SPNMG(NSPX)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SPNMG(NSPX)(1:IDB) ) THEN
              NSP = NSPX + NSPL + NSPS
              GOTO 540
            ENDIF
  530     CONTINUE
!
!---      pH  ---
!
          IF( BDUM(1:NCHB).EQ.'ph' .AND. ISPLK(1).NE.0 ) THEN
            NSP = MOD(ISPLK(1),1000)
            ISPLK(1) = ISPLK(1) + 1000
            IVAR = 2
            ADDER = 7.D+0
!
!---        Verify that species linked to pH is a conservation
!           component species  ---
!
            DO 532 NEQ = 1,NEQC
              IF( NSP.EQ.IEQ_C(2,NEQ) ) GOTO 540
  532       CONTINUE
            INDX = 4
            CHMSG = 'pH Species not a Conservation ' //
     &        'Component Species: ' // BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
          INDX = 4
          CHMSG = 'Unrecognized Reactive Species: ' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
  540     CONTINUE
!
!---      Conservation- or kinetic-component species  ---
!
          IF( INDEX( BDUM(1:),'total_' ).NE.0 ) THEN
            IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
              IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
                CALL RDINBP( C(1,NSL),ADDER,ICT(1,NSL),IVAR,UNTS )
              ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
                CALL RDINAP( C(1,NSL),ADDER,ICT(1,NSL),IVAR,UNTS )
              ELSE
                CALL RDINFP( C(1,NSL),VAR,ADDER,ICT(1,NSL),
     &            IVAR,UNTS )
              ENDIF
              CLOSE(UNIT=26)
            ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &        INDEX(ADUM(1:),'zonation').NE.0 )  THEN
              CALL RDINZP( C(1,NSL),VAR(1),ADDER,ICT(1,NSL),
     &          IVAR,IROCK )
            ELSE
              CALL RDINIP( C(1,NSL),VAR,ADDER,ICT(1,NSL),IVAR,
     &          IDOM )
            ENDIF
          ELSE
            IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
              IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
                CALL RDINBP( SP_C(1,NSP),ADDER,IC_SP(1,NSP),IVAR,UNTS )
              ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
                CALL RDINAP( SP_C(1,NSP),ADDER,IC_SP(1,NSP),IVAR,UNTS )
              ELSE
                CALL RDINFP( SP_C(1,NSP),VAR,ADDER,IC_SP(1,NSP),
     &            IVAR,UNTS )
              ENDIF
              CLOSE(UNIT=26)
            ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &        INDEX(ADUM(1:),'zonation').NE.0 )  THEN
              CALL RDINZP( SP_C(1,NSP),VAR(1),ADDER,IC_SP(1,NSP),
     &          IVAR,IROCK )
            ELSE
              CALL RDINIP( SP_C(1,NSP),VAR,ADDER,IC_SP(1,NSP),IVAR,
     &          IDOM )
            ENDIF
          ENDIF

        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Boundary Condition ' // 
     &      'Variable: ' // ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
 1000 CONTINUE
 1002 CONTINUE
      WRITE(IWR,'(/)')
!
!---  Check thermodynamic and hydrologic initial states  ---
!
      IEOX = IEO
      IEO = 1
      CALL CHK_GT
      ISVCX = ISVC
      ISVC = 1
      CALL PROP_GT
      ISVC = ISVCX
      IEO = IEOX
      
!
!---  Assign values for initial condition type boundary conditions  ---
!
      DO NB = 1,NBC
        ITFX = MOD(IBCT(2,NB),100)
        IF( IBCT(1,NB).EQ.4 .OR. ITFX.EQ.4 ) THEN
          N = IBCN(NB)
          IF( IBCD(NB).EQ.-3 ) THEN
            NPZ = NSZ(N)
            DB = DZGP(NPZ)
            GB = GRVZ(NPZ)*DB
          ELSEIF( IBCD(NB).EQ.-2 ) THEN
            NPY = NSY(N)
            DB = DYGP(NPY)*RP(ID(N))
            GB = GRVY(NPY)*DB
          ELSEIF( IBCD(NB).EQ.-1 ) THEN
            NPX = NSX(N)
            DB = DXGP(NPX)
            GB = GRVX(NPX)*DB
          ELSEIF( IBCD(NB).EQ.1 ) THEN
            NQX = NSX(N)+1
            IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
            DB = -DXGP(NQX)
            GB = GRVX(NQX)*DB
          ELSEIF( IBCD(NB).EQ.2 ) THEN
            NQY = NSY(N)+IFLD
            IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
            DB = -DYGP(NQY)*RP(ID(N))
            GB = GRVY(NQY)*DB
          ELSEIF( IBCD(NB).EQ.3 ) THEN
            NQZ = NSZ(N)+IJFLD
            IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
            DB = -DZGP(NQZ)
            GB = GRVZ(NQZ)*DB
          ENDIF
        ENDIF
!
!---    Energy initial conditions  ---
!
        IF( IBCT(1,NB).EQ.4 ) TB(1,NB) = T(2,N)
!
!---    Fluid flow initial conditions  ---
!
        IF( ITFX.EQ.4 ) THEN
          PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
          PGB(1,NB) = PG(2,N) + RHOG(2,N)*GB
          XLAB(1,NB) = XLA(2,N)
          XGWB(1,NB) = XGW(2,N)
          YLSB(1,NB) = YLS(2,N)
          SLB(1,NB) = SL(2,N)
          IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.3 ) THEN
            IBCT(2,NB) = 104
            IBCT(3,NB) = 3
            IBCT(4,NB) = 3
          ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
            IBCT(2,NB) = 204
            IBCT(3,NB) = 3
            IBCT(4,NB) = 3
          ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
            IBCT(2,NB) = 304
            IBCT(3,NB) = 3
            IBCT(4,NB) = 3
          ENDIF
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIBC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIC_GT
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
!     Read input file for initial conditions information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE HYST
      USE GRID
      USE FILES
      USE FDVS
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
      CHARACTER*64 ADUM,BDUM,FDUM,FMDUM,UNTS
      CHARACTER*24 CHLB(3)
      CHARACTER*512 CHDUM
      CHARACTER*20 FORM5
      CHARACTER*38 FORM6
      CHARACTER*39 FORM7
      INTEGER IDOM(6)
      REAL*8 VAR(5)
      LOGICAL FCHK
!
!----------------------Data Statements---------------------------------!
!
      DATA CHLB /'X-Direction Gradient, ','Y-Direction Gradient, ',
     &           'Z-Direction Gradient, '/
      DATA FORM5 /'(10(1PE22.15,1X),I3)'/
      DATA FORM6 /'(10(1PE22.15,1X),I3,1X,1(1PE22.15,1X))'/
      DATA FORM7 /'(10(1PE22.15,1X),I3,1X,10(1PE22.15,1X))'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIC_GT'
      IERR = 0
!
!---  Write card information to ouput file  ---
!
      CARD = 'Initial Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Initialize initial condition states for air and salt  ---
!
      DO 10 N = 1,NFLD
        ICAIR(N) = 0
        ICBRN(N) = 0
   10 CONTINUE
!
!---  Restart file will be read for initial conditions  ---
!
      IF( IEO.EQ.2 ) THEN
        INDX = 2
        CALL RDRST(INDX)
        ISIC = 3
      ENDIF
!
!---  Read saturation initial condition option ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Initial Saturation Option: '
      CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
      IF( INDEX(ADUM(1:),'yz').NE.0 .AND.
     &  INDEX(ADUM(1:),'restart').NE.0 ) THEN
        ISIC = 3
!
!---    Check for external restart file  ---
!
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
        NCH = INDEX(FDUM,'  ')-1
        INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          INDX = 4
          CHMSG = 'Missing YZ Restart File: ' // FDUM(1:NCH)
          CALL WRMSGS( INDX )
        ELSEIF( FDUM.EQ.'unformatted' ) THEN
          INDX = 4
          CHMSG = 'YZ File Format: ' // FDUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',FORM='FORMATTED')
        WRITE(IWR,'(/,2A)') 'YZ Restart File: ',FDUM(1:NCH)
        DO N = 1,23
          READ(26,'(A)') CHDUM
        ENDDO
!
!---    Read timing data  ---
!

        READ(26,'(7(1PE22.15),8(I9))') TMPSX,TMPDX,TMPXX,
     &    TMPNX,TMPAX,RSDMX,TMPCX,NRIMX,NRSTX,NFLDX,NSOLUX,
     &    IROMX,IRFCX,NSPRX,NFBNX






        IF( JFLD*KFLD.NE.NFLDX ) THEN
          INDX = 4
          CHMSG = 'Incompatible YZ Restart File: ' // FDUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        NRSV = 11
        IF( NSOLUX.LE.0 ) THEN
          WRITE( FORM5(2:3),'(I2)' ) NRSV
          DO K = 1,KFLD
            DO J = 1,JFLD
              N = ND(1,J,K)
              READ(26,FORM5) T(2,N),PL(2,N),PG(2,N),PVW(2,N),SL(2,N),
     &          SGT(2,N),ASLMIN(2,N),SG(2,N),PVA(2,N),PCMP(N),TCMP(N),
     &          NPHAZ(2,N)
              DO I = 2,IFLD
                NX = ND(I,J,K)
                T(2,NX) = T(2,N)
                PL(2,NX) = PL(2,N)
                PG(2,NX) = PG(2,N)
                SL(2,NX) = SL(2,N)
                SGT(2,NX) = SGT(2,N)
                ASLMIN(2,NX) = ASLMIN(2,N)
                SG(2,NX) = SG(2,N)
                PVA(2,NX) = PVA(2,N)
                SI(2,NX) = SI(2,N)
                PI(2,NX) = PI(2,N)
                PCMP(NX) = PCMP(N)
                TCMP(NX) = TCMP(N)
                NPHAZ(2,NX) = NPHAZ(2,N)
              ENDDO
            ENDDO
          ENDDO
        ELSEIF( NSOLUX.LT.10 ) THEN
          WRITE( FORM6(2:3),'(I2)' ) NRSV
          WRITE( FORM6(23:23),'(I1)' ) NSOLUX
          DO K = 1,KFLD
            DO J = 1,JFLD
              N = ND(1,J,K)
              READ(26,FORM6) T(2,N),PL(2,N),PG(2,N),PVW(2,N),SL(2,N),
     &          SGT(2,N),ASLMIN(2,N),SG(2,N),PVA(2,N),PCMP(N),TCMP(N),
     &          NPHAZ(2,N),(C(N,M),M=1,NSOLUX)
              DO I = 2,IFLD
                NX = ND(I,J,K)
                T(2,NX) = T(2,N)
                PL(2,NX) = PL(2,N)
                PG(2,NX) = PG(2,N)
                SL(2,NX) = SL(2,N)
                SGT(2,NX) = SGT(2,N)
                ASLMIN(2,NX) = ASLMIN(2,N)
                SG(2,NX) = SG(2,N)
                PVA(2,NX) = PVA(2,N)
                SI(2,NX) = SI(2,N)
                PI(2,NX) = PI(2,N)
                PCMP(NX) = PCMP(N)
                TCMP(NX) = TCMP(N)
                NPHAZ(2,NX) = NPHAZ(2,N)
                DO M = 1,NSOLUX
                  C(NX,M) = C(N,M)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          WRITE( FORM7(2:3),'(I2)' ) NRSV
          WRITE( FORM7(23:24),'(I2)' ) NSOLUX
          DO K = 1,KFLD
            DO J = 1,JFLD
              N = ND(1,J,K)
              READ(26,FORM6) T(2,N),PL(2,N),PG(2,N),PVW(2,N),SL(2,N),
     &          SGT(2,N),ASLMIN(2,N),SG(2,N),PVA(2,N),PCMP(N),TCMP(N),
     &          NPHAZ(2,N),(C(N,M),M=1,NSOLUX)
              DO I = 2,IFLD
                NX = ND(I,J,K)
                T(2,NX) = T(2,N)
                PL(2,NX) = PL(2,N)
                PG(2,NX) = PG(2,N)
                SL(2,NX) = SL(2,N)
                SGT(2,NX) = SGT(2,N)
                ASLMIN(2,NX) = ASLMIN(2,N)
                SG(2,NX) = SG(2,N)
                PVA(2,NX) = PVA(2,N)
                SI(2,NX) = SI(2,N)
                PI(2,NX) = PI(2,N)
                PCMP(NX) = PCMP(N)
                TCMP(NX) = TCMP(N)
                NPHAZ(2,NX) = NPHAZ(2,N)
                DO M = 1,NSOLUX
                  C(NX,M) = C(N,M)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        CLOSE(UNIT=26)
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ELSEIF( INDEX(ADUM(1:),'hydrostatic').NE.0 .AND.
     &  INDEX(ADUM(1:),'geothermal').NE.0 .AND.
     &  INDEX(ADUM(1:),'zon').NE.0 ) THEN
        VARB = 'Number of Hydrostatic Geothermal Zonal Lines: '
        CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
        DO NL = 1,NLIN
          NLX = (NL-1)*10
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISIC = NL*10 + 4
          VARB = 'Interval Lower Z Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+9))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+9)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(NLX+9),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+9),', m)'
          VARB = 'Interval Upper Z Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+10))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+10)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(NLX+10),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+10),', m)'
          IF( VSLC(NLX+9).GT.VSLC(NLX+10) ) THEN
            VSLCX = VSLC(NLX+10)
            VSLC(NLX+10) = VSLC(NLX+9)
            VSLC(NLX+9) = VSLCX
          ENDIF
          VARB = 'Reference Aqueous Pressure'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+1)
          INDX = 0
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,VSLC(NLX+1),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+1),', Pa)'
          VSLC(NLX+1) = VSLC(NLX+1) - PATM
          IF( INDEX( CHDUM(1:),'file' ).NE.0 ) THEN
!
!---        Read salt mass fraction from an external file, read
!           x-,y-,and z-coordinates for reference pressure  ---
!
            VARB = 'Reference-Pressure X Location'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+6))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+6)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VSLC(NLX+6),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+6),', m)'
            VARB = 'Reference-Pressure Y Location'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+7))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+7)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VSLC(NLX+7),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+7),', m)'
            VARB = 'Reference-Pressure Z Location'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+8))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+8)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VSLC(NLX+8),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+8),', m)'
!
!---      Read salt mass fraction from reference value, z-location, 
!         and gradient, read z-coordinate for reference pressure  ---
!
          ELSE
            VARB = 'Reference-Pressure Z Location'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+2)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VSLC(NLX+2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+2),', m)'
          ENDIF
          VARB = 'Reference Temperature'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+3))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+3)
          INDX = 0
          IUNK = 1
          CALL RDUNIT(UNTS,VSLC(NLX+3),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+3),', C)'
          VARB = 'Reference-Temperature Z Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+4))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+4)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(NLX+4),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+4),', m)'
          VARB = 'Geothermal Gradient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+5))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+5)
          INDX = 0
          IUNK = 1
          IUNM = -1
          CALL RDUNIT(UNTS,VSLC(NLX+5),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+5),', C/m)'
!
!---      Check for salt mass fraction read from an external file  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ISTARTX = ISTART
            ICOMMAX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---        Read salt mass fraction from an external file  ---
!
            IF( INDEX( ADUM(1:),'file' ).NE.0 ) THEN
              ISIC = NL*10 + 5
!
!---          Check for external file  ---
!
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              NCH = INDEX(FDUM,'  ')-1
              INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
              IF( .NOT.FCHK ) THEN
                INDX = 4
                CHMSG = 'Missing Salt Mass Fraction File: ' // 
     &            FDUM(1:NCH)
                CALL WRMSGS( INDX )
              ELSEIF( FDUM.EQ.'unformatted' ) THEN
                INDX = 4
                CHMSG = 'Salt Mass Fraction File Format: ' //  
     &            FDUM(1:NCH)
                CALL WRMSGS( INDX )
              ENDIF
              OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',
     &          FORM='FORMATTED')
              WRITE(IWR,'(/,2A)') 'Salt Mass Fraction File: ',
     &          FDUM(1:NCH)
              UNTS = 'null'
              ADDER = 0.D+0
              INDX = 2
              CALL RDINAS( YLS,ADDER,UNTS,INDX )
              CLOSE(UNIT=26)
!
!---        Read salt mass fraction from reference value, z-location, 
!           and gradient  ---
!
            ELSE
              ISTART = ISTARTX
              ICOMMA = ICOMMAX
              VARB = 'Reference Salt Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+6))
              WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',
     &          VSLC(NLX+6)
              VARB = 'Reference-Salinity (Salt Mass Fraction) ' // 
     &          'Z Location'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+7))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+7)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,VSLC(NLX+7),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+7),', m)'
              VARB = 'Geosalinity (Salt Mass Fraction) Gradient'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(NLX+8))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
              WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(NLX+8)
              INDX = 0
              IUNM = -1
              CALL RDUNIT(UNTS,VSLC(NLX+8),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(NLX+8),', 1/m)'
            ENDIF
          ENDIF
        ENDDO
        IF( NSOLU+NEQC.GT.0 ) THEN
          GOTO 22
        ELSE
          GOTO 1002
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'hydrostatic').NE.0 .AND.
     &  INDEX(ADUM(1:),'geothermal').NE.0 ) THEN
        ISIC = 4
        VARB = 'Reference Aqueous Pressure'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(1)
        INDX = 0
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        CALL RDUNIT(UNTS,VSLC(1),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(1),', Pa)'
        VSLC(1) = VSLC(1) - PATM
        IF( INDEX( CHDUM(1:),'file' ).NE.0 ) THEN
!
!---      Read salt mass fraction from an external file, read
!         x-,y-,and z-coordinates for reference pressure  ---
!
          VARB = 'Reference-Pressure X Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(6))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(6)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(6),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(6),', m)'
          VARB = 'Reference-Pressure Y Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(7))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(7)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(7),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(7),', m)'
          VARB = 'Reference-Pressure Z Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(8))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(8)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(8),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(8),', m)'
!
!---    Read salt mass fraction from reference value, z-location, 
!       and gradient, read z-coordinate for reference pressure  ---
!
        ELSE
          VARB = 'Reference-Pressure Z Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(2))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
          WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(2)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,VSLC(2),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(2),', m)'
        ENDIF
        VARB = 'Reference Temperature'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(3))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(3)
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,VSLC(3),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(3),', C)'
        VARB = 'Reference-Temperature Z Location'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(4))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(4)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,VSLC(4),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(4),', m)'
        VARB = 'Geothermal Gradient'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(5))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
        WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(5)
        INDX = 0
        IUNK = 1
        IUNM = -1
        CALL RDUNIT(UNTS,VSLC(5),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(5),', C/m)'
!
!---    Check for salt mass fraction read from an external file  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          ISTARTX = ISTART
          ICOMMAX = ICOMMA
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---      Read salt mass fraction from an external file  ---
!
          IF( INDEX( ADUM(1:),'file' ).NE.0 ) THEN
            ISIC = 5
!
!---        Check for external file  ---
!
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            NCH = INDEX(FDUM,'  ')-1
            INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
            IF( .NOT.FCHK ) THEN
              INDX = 4
              CHMSG = 'Missing Salt Mass Fraction File: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ELSEIF( FDUM.EQ.'unformatted' ) THEN
              INDX = 4
              CHMSG = 'Salt Mass Fraction File Format: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
            OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',FORM='FORMATTED')
            WRITE(IWR,'(/,2A)') 'Salt Mass Fraction File: ',FDUM(1:NCH)
            UNTS = 'null'
            ADDER = 0.D+0
            INDX = 2
            CALL RDINAS( YLS,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
!
!---      Read salt mass fraction from reference value, z-location, 
!         and gradient  ---
!
          ELSE
            ISTART = ISTARTX
            ICOMMA = ICOMMAX
            VARB = 'Reference Salt Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(6))
            WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),': ',VSLC(6)
            VARB = 'Reference-Salinity (Salt Mass Fraction) Z Location'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(7))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(7)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VSLC(7),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(7),', m)'
            VARB = 'Geosalinity (Salt Mass Fraction) Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VSLC(8))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4,$)') UNTS(1:NCH),': ',VSLC(8)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,VSLC(8),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VSLC(8),', 1/m)'
          ENDIF
        ENDIF
        IF( NSOLU+NEQC.GT.0 ) THEN
          GOTO 22
        ELSE
          GOTO 1002
        ENDIF
      ENDIF      
      CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
      IF( IEO.EQ.2 ) GOTO 20
      WRITE(IWR,'(/,A)') VARB(1:IVR)
      WRITE(IWR,'(2X,A)') ADUM(1:NCHA)
      WRITE(IWR,'(2X,A)') BDUM(1:NCHB)
      IF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
        IF( INDEX(BDUM(1:),'gas pressure').NE.0 ) THEN
          ISIC = 1
        ELSEIF( INDEX(BDUM(1:),'aqueous pressure').NE.0 ) THEN
          ISIC = 2
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &      ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'gas pressure').NE.0 ) THEN
        IF( INDEX(BDUM(1:),'aqueous saturation').NE.0 ) THEN
          ISIC = 1
        ELSEIF( INDEX(BDUM(1:),'aqueous pressure').NE.0 ) THEN
          ISIC = 3
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &      ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
        IF( INDEX(BDUM(1:),'aqueous saturation').NE.0 ) THEN
          ISIC = 2
        ELSEIF( INDEX(BDUM(1:),'gas pressure').NE.0 ) THEN
          ISIC = 3
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &      ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ENDIF
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Initial Saturation Option: ' //
     &    ADUM(1:NCHA) // ',' // BDUM(1:NCHB)
        CALL WRMSGS( INDX )
      ENDIF
   20 CONTINUE
!
!---  Read initial conditions  ---
!
   22 WRITE(IWR,'(/,A)') 'Initial Condition Variable(s) and Domain(s)'
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Initial Condition Cards: '
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO 1000 NL = 1, NLIN
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Initial Condition Variable: '
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'overwrite').EQ.0 .AND.
     &    ( IEO.EQ.2 ) ) GOTO 1000
        NPZX = 1
        IF( INDEX( ADUM(1:),'overwrite').NE.0 .AND.
     &    ( IEO.EQ.2 ) ) NPZX = 0
        IF( INDEX( ADUM(1:),'aqueous pres' ).NE.0 ) THEN
          VARB = 'Initial Aqueous Pressure'
          IUNM = -1
          IUNKG = 1
          IUNS = -2
        ELSEIF( INDEX( ADUM(1:),'gas pres' ).NE.0 ) THEN
          VARB = 'Initial Gas Pressure'
          IUNM = -1
          IUNKG = 1
          IUNS = -2
        ELSEIF( INDEX( ADUM(1:),'temperature' ).NE.0 ) THEN
          VARB = 'Initial Temperature'
          IUNK = 1
        ELSEIF( INDEX( ADUM(1:),'aqueous sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous Saturation'
        ELSEIF( INDEX( ADUM(1:),'relative' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'trapped gas' ).NE.0 ) THEN
          VARB = 'Initial Relative Trapped Gas Saturation'
        ELSEIF( INDEX( ADUM(1:),'trapped gas' ).NE.0 ) THEN
          VARB = 'Initial Trapped Gas Saturation'
        ELSEIF( INDEX( ADUM(1:),'air' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-Air Relative Saturation'
        ELSEIF( INDEX( ADUM(1:),'air' ).NE.0 .AND.
     &    ( INDEX( ADUM(1:),'mass' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'frac' ).NE.0 ) ) THEN
          VARB = 'Initial Aqueous-Air Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'air' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'conc' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-Air Concentration'
          IUNKG = 1
          IUNM = -3
        ELSEIF( INDEX( ADUM(1:),'water' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &    (INDEX( ADUM(1:),'sat' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'humid' ).NE.0) ) THEN
          VARB = 'Initial Gas-Water Relative Saturation/Humidity'
        ELSEIF( INDEX( ADUM(1:),'water' ).NE.0 .AND.
     &    ( INDEX( ADUM(1:),'mass' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'frac' ).NE.0 ) ) THEN
          VARB = 'Initial Gas-Water Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'air' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'conc' ).NE.0 ) THEN
          VARB = 'Initial Gas-Water Concentration'
          IUNKG = 1
          IUNM = -3
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-Salt Relative Saturation'
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    ( INDEX( ADUM(1:),'mass' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'frac' ).NE.0 ) ) THEN
          VARB = 'Initial Aqueous-Salt Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'conc' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-Salt Concentration'
          IUNKG = 1
          IUNM = -3
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'molality' ).NE.0 ) THEN
          VARB = 'Initial Aqueous-Salt Molality'
          IUNMOL = 1
          IUNKG = -1
        ELSEIF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          VARB = 'Initial Solute Concentration'
          IF( INDEX( ADUM(1:),'activity' ).NE.0 ) THEN
            IVAR = 4
          ELSE
            IVAR = 1
            IUNM = -3
          ENDIF

        ELSEIF( INDEX( ADUM(1:),'specie' ).NE.0 ) THEN
          VARB = 'Reactive Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          VARB = 'Initial Reactive Species Concentration'
!
!---      Set species units  ---
!
          IUNMOL = 1
          IF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
            IVAR = 4
            IUNM = -3
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'molal').NE.0 ) THEN
              IVAR = 3
              IUNKG = -1
            ELSE
              IVAR = 2
              IUNM = -3
            ENDIF
          ELSE
            IVAR = 1
            IUNM = -3
          ENDIF
          IF( IEO.EQ.2 ) IVAR = IVAR+10

        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(1))
        CALL RDCHR(ISTART,ICOMMA,NCHU,CHDUM,UNTS)
!
!---  Read initial conditions input from an external file  ---
!
        IF( INDEX( ADUM(1:),'file' ).NE.0 ) THEN
          IF( INDEX( ADUM(1:),'binary' ).NE.0 ) THEN
            WRITE(IWR,'(2X,3A)') ADUM(1:NCHA),',',UNTS(1:NCHU)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            NCH = INDEX(FDUM,'  ')-1
!
!---        Check for external file  ---
!
            INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
            IF( .NOT.FCHK ) THEN
              INDX = 4
              CHMSG = 'Missing Initial Conditions File: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ELSEIF( FDUM.EQ.'formatted' ) THEN
              INDX = 4
              CHMSG = 'Initial Conditions File Format: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
            OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',
     &        FORM='UNFORMATTED')
            WRITE(IWR,'(/,2A)') 'Initial Conditions File: ',FDUM(1:NCH)
          ELSEIF( INDEX( ADUM(1:),'ascii' ).NE.0 ) THEN
            WRITE(IWR,'(2X,3A)') ADUM(1:NCHA),',',UNTS(1:NCHU)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            NCH = INDEX(FDUM,'  ')-1
!
!---        Check for external file  ---
!
            INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
            IF( .NOT.FCHK ) THEN
              INDX = 4
              CHMSG = 'Missing Initial Conditions File: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ELSEIF( FDUM.EQ.'unformatted' ) THEN
              INDX = 4
              CHMSG = 'Initial Conditions File Format: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
            OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',
     &        FORM='FORMATTED')
            WRITE(IWR,'(/,2A)') 'Initial Conditions File: ',FDUM(1:NCH)
          ELSE
            WRITE(IWR,'(2X,4A,1PE11.4)') ADUM(1:NCHA),
     &        ' (Default Value), ',UNTS(1:NCHU),': ',VAR(1)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            NCH = INDEX(FDUM,'  ')-1
!
!---        Check for external file  ---
!
            INQUIRE( FILE=FDUM(1:NCH), FORM=FMDUM, EXIST=FCHK )
            IF( .NOT.FCHK ) THEN
              INDX = 4
              CHMSG = 'Missing Initial Conditions File: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ELSEIF( FDUM.EQ.'unformatted' ) THEN
              INDX = 4
              CHMSG = 'Initial Conditions File Format: ' // FDUM(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
            OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='OLD',FORM='FORMATTED')
            WRITE(IWR,'(/,2A)') 'Initial Conditions File: ',FDUM(1:NCH)
            INDX = 0
            CALL RDUNIT( UNTS,VAR(1),INDX )
          ENDIF
!
!---  Read initial conditions according to rock/soil zonations  ---
!
        ELSEIF( INDEX( ADUM(1:),'rock' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'zonation' ).NE.0 ) THEN
          VARB = 'Rock/Soil Name'
          CALL RDCHR(ISTART,ICOMMA,NCHF,CHDUM,FDUM)
!
!---  Search known rock types for a matching type ---
!
          DO 30 M = 1, NROCK
            IF( FDUM .EQ. ROCK(M)) THEN
            IROCK = M
            GOTO 40
          ENDIF
   30     CONTINUE
          INDX = 2
          CHMSG = 'Unrecognized Rock/Soil Type: '//FDUM
          CALL WRMSGS( INDX )
          GOTO 1000
   40     CONTINUE
          WRITE(IWR,'(2X,3A,1PE11.4,2A)') ADUM(1:NCHA),UNTS(1:NCHU),
     &      ': ',VAR(1),' Rock/Soil Type: ',FDUM(1:NCHF)
          INDX = 0
          CALL RDUNIT( UNTS,VAR(1),INDX )
!
!---    Read initial conditions input from the input file  ---
!
        ELSE
          WRITE(IWR,'(2X,4A,1PE11.4)') ADUM(1:NCHA),', ',
     &      UNTS(1:NCHU),': ',VAR(1)
          INDX = 0
          CALL RDUNIT( UNTS,VAR(1),INDX )
          INDX = 2
          VAR(5) = 1.D+0
          NCH = INDEX( UNTS,'  ' ) - 1
          IF( UNTS(1:NCH).EQ.'f' .OR. UNTS(1:NCH).EQ.'r' ) THEN
            VAR(5) = VAR(5)/1.8D+0
          ELSEIF( UNTS(1:NCH).EQ.'c' .OR. UNTS(1:NCH).EQ.'k' ) THEN
            VAR(5) = 1.D+0
          ELSE
            CALL RDUNIT( UNTS,VAR(5),INDX )
          ENDIF
          VARB = 'Initial Condition Variable Gradient: '
          DO 100 I = 2,4
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(I))
            VAR(I) = VAR(I)*VAR(5)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') CHLB(I-1),', ',UNTS(1:NCH),
     &        ': ',VAR(I)
            INDX = 0
            IUNM = -1
            CALL RDUNIT( UNTS,VAR(I),INDX )
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(I),', 1/m)'
  100     CONTINUE
!
!---      Read domain indices  ---
!
          VARB = 'Initial Condition Domain Index: '
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(1))
          IF( IDOM(1).LT.1 .OR. IDOM(1).GT.IFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Lower I-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(1)
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(2))
          IF( IDOM(2).LT.1 .OR. IDOM(2).GT.IFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Upper I-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(2)
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(3))
          IF( IDOM(3).LT.1 .OR. IDOM(3).GT.JFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Lower J-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(3)
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(4))
          IF( IDOM(4).LT.1 .OR. IDOM(4).GT.JFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Upper J-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(4)
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(5))
          IF( IDOM(5).LT.1 .OR. IDOM(5).GT.KFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Lower K-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(5)
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDINT(ISTART,ICOMMA,CHDUM,IDOM(6))
          IF( IDOM(6).LT.1 .OR. IDOM(6).GT.KFLD ) THEN
            INDX = 7
            CHMSG = 'Out-of-Range Upper K-Index: ' // ADUM(1:NCHA)
            IMSG = IDOM(6)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Read variables  ---
!
        IF( INDEX(ADUM(1:),'temp').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNK = 1
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINBR( T,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINBS( T,ADDER,UNTS,INDX )
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINAR( T,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINAS( T,ADDER,UNTS,INDX )
              ENDIF
            ELSE
              IF( NPZX.EQ.0 ) THEN
                CALL RDINFR( T,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINFS( T,VAR,ADDER,UNTS,INDX )
              ENDIF    
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( T,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( T,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( T,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( T,VAR,ADDER,IDOM,INDX )
            ENDIF    
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous press').NE.0 ) THEN
          ADDER = -PATM
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINBR( PL,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINBS( PL,ADDER,UNTS,INDX )
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINAR( PL,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINAS( PL,ADDER,UNTS,INDX )
              ENDIF
            ELSE
              IF( NPZX.EQ.0 ) THEN
                CALL RDINFR( PL,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINFS( PL,VAR,ADDER,UNTS,INDX )
              ENDIF    
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( PL,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( PL,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( PL,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( PL,VAR,ADDER,IDOM,INDX )
            ENDIF    
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas pres').NE.0 ) THEN
          ADDER = -PATM
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINBR( PG,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINBS( PG,ADDER,UNTS,INDX )
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINAR( PG,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINAS( PG,ADDER,UNTS,INDX )
              ENDIF
            ELSE
              IF( NPZX.EQ.0 ) THEN
                CALL RDINFR( PG,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINFS( PG,VAR,ADDER,UNTS,INDX )
              ENDIF    
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( PG,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( PG,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( PG,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( PG,VAR,ADDER,IDOM,INDX )
            ENDIF    
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous sat').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( SL,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( SL,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( SL,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( SL,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( SL,VAR,ADDER,IDOM,INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'rel').NE.0 .AND.
     &    INDEX(ADUM(1:),'trapped gas').NE.0 ) THEN
          ADDER = 1.D+2
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINBR( SGT,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINBS( SGT,ADDER,UNTS,INDX )
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINAR( SGT,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINAS( SGT,ADDER,UNTS,INDX )
              ENDIF
            ELSE
              IF( NPZX.EQ.0 ) THEN
                CALL RDINFR( SGT,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINFS( SGT,VAR,ADDER,UNTS,INDX )
              ENDIF    
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( SGT,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( SGT,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( SGT,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( SGT,VAR,ADDER,IDOM,INDX )
            ENDIF    
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'trapped').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINBR( SGT,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINBS( SGT,ADDER,UNTS,INDX )
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              IF( NPZX.EQ.0 ) THEN
                CALL RDINAR( SGT,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINAS( SGT,ADDER,UNTS,INDX )
              ENDIF
            ELSE
              IF( NPZX.EQ.0 ) THEN
                CALL RDINFR( SGT,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
              ELSE
                CALL RDINFS( SGT,VAR,ADDER,UNTS,INDX )
              ENDIF    
            ENDIF
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( SGT,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( SGT,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( SGT,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( SGT,VAR,ADDER,IDOM,INDX )
            ENDIF    
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0 ) THEN
          IVAR = 3
          IF( INDEX(ADUM(1:),'molality').NE.0 ) THEN
            IVAR = 4
          ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .OR.
     &      INDEX(ADUM(1:),'frac').NE.0 ) THEN
            IVAR = 3
          ELSEIF( INDEX(ADUM(1:),'rel').NE.0 .OR.
     &      INDEX(ADUM(1:),'sat').NE.0 ) THEN
            IVAR = 2
          ELSEIF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &      INDEX(ADUM(1:),'conc').NE.0 ) THEN
            IVAR = 1
          ENDIF
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( IVAR.EQ.0 .OR. IVAR.EQ.1 ) THEN
              IUNM = -3
              IUKG = 1
            ELSEIF( IVAR.EQ.4 ) THEN
              IUNMOL = 1
              IUNKG = -1
            ENDIF
            IF( NPZX.EQ.0 ) THEN
              CALL RDINFR( TMS,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
            ELSE
              CALL RDINFS( TMS,VAR,ADDER,UNTS,INDX )
            ENDIF    
            CLOSE(UNIT=26)
            DO 206 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 206
              ICBRN(N) = IVAR
  206       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( TMS,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( TMS,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
            DO 208 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 208
              IF( IZ(N).EQ.IROCK ) ICBRN(N) = IVAR
  208       CONTINUE
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( TMS,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( TMS,VAR,ADDER,IDOM,INDX )
            ENDIF    
            DO K = IDOM(5),IDOM(6)
            DO J = IDOM(3),IDOM(4)
            DO I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) CYCLE
              ICBRN(N) = IVAR
            ENDDO
            ENDDO
            ENDDO
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'air').NE.0 ) THEN
          IVAR = 3
          IF( INDEX(ADUM(1:),'mass').NE.0 .OR.
     &      INDEX(ADUM(1:),'frac').NE.0 ) THEN
            IVAR = 3
            ADDER = 0.D+0
          ELSEIF( INDEX(ADUM(1:),'rel').NE.0 .OR.
     &      INDEX(ADUM(1:),'sat').NE.0 ) THEN
            IVAR = 2
            ADDER = 0.D+0
          ELSEIF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &      INDEX(ADUM(1:),'conc').NE.0 ) THEN
            IVAR = 1
            ADDER = 0.D+0
          ENDIF
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( IVAR.EQ.0 .OR. IVAR.EQ.1 ) THEN
              IUNM = -3
              IUKG = 1
            ELSEIF( IVAR.EQ.4 ) THEN
              IUNM = -1
              IUNKG = 1
              IUNS = -2
            ENDIF
            IF( NPZX.EQ.0 ) THEN
              CALL RDINFR( PVA,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
            ELSE
              CALL RDINFS( PVA,VAR,ADDER,UNTS,INDX )
            ENDIF    
            CLOSE(UNIT=26)
            DO N = 1,NFLD
              IF( IXP(N).EQ.0 ) CYCLE
              ICAIR(N) = IVAR
            ENDDO
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( PVA,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( PVA,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
            DO N = 1,NFLD
              IF( IXP(N).EQ.0 ) CYCLE
              IF( IZ(N).EQ.IROCK ) ICAIR(N) = IVAR
            ENDDO
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( PVA,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( PVA,VAR,ADDER,IDOM,INDX )
            ENDIF    
            DO K = IDOM(5),IDOM(6)
            DO J = IDOM(3),IDOM(4)
            DO I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) CYCLE
              ICAIR(N) = IVAR
            ENDDO
            ENDDO
            ENDDO
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'water').NE.0 ) THEN
          IVAR = 3
          IF( INDEX(ADUM(1:),'mass').NE.0 .OR.
     &      INDEX(ADUM(1:),'frac').NE.0 ) THEN
            IVAR = 3
            ADDER = 0.D+0
          ELSEIF( INDEX(ADUM(1:),'rel').NE.0 .OR.
     &      INDEX(ADUM(1:),'sat').NE.0 .OR.
     &      INDEX(ADUM(1:),'humid').NE.0 ) THEN
            IVAR = 2
            ADDER = 0.D+0
          ELSEIF( INDEX(ADUM(1:),'aqu').NE.0 .AND.
     &      INDEX(ADUM(1:),'conc').NE.0 ) THEN
            IVAR = 1
            ADDER = 0.D+0
          ENDIF
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( IVAR.EQ.0 .OR. IVAR.EQ.1 ) THEN
              IUNM = -3
              IUKG = 1
            ELSEIF( IVAR.EQ.4 ) THEN
              IUNM = -1
              IUNKG = 1
              IUNS = -2
            ENDIF
            IF( NPZX.EQ.0 ) THEN
              CALL RDINFR( PVW,VAR,ADDER,NPHAZ,NPZX,UNTS,INDX )
            ELSE
              CALL RDINFS( PVW,VAR,ADDER,UNTS,INDX )
            ENDIF    
            CLOSE(UNIT=26)
            DO N = 1,NFLD
              IF( IXP(N).EQ.0 ) CYCLE
              ICAIR(N) = IVAR
            ENDDO
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            IF( NPZX.EQ.0 ) THEN
              CALL RDINZR( PVW,VAR(1),ADDER,NPHAZ,NPZX,IROCK,INDX )
            ELSE
              CALL RDINZS( PVW,VAR(1),ADDER,IROCK,INDX )
            ENDIF    
            DO N = 1,NFLD
              IF( IXP(N).EQ.0 ) CYCLE
              IF( IZ(N).EQ.IROCK ) ICAIR(N) = IVAR
            ENDDO
          ELSE
            IF( NPZX.EQ.0 ) THEN
              CALL RDINIR( PVW,VAR,ADDER,NPHAZ,NPZX,IDOM,INDX )
            ELSE
              CALL RDINIS( PVW,VAR,ADDER,IDOM,INDX )
            ENDIF    
            DO K = IDOM(5),IDOM(6)
            DO J = IDOM(3),IDOM(4)
            DO I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) CYCLE
              ICAIR(N) = IVAR
            ENDDO
            ENDDO
            ENDDO
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'activity').NE.0 ) THEN
            IVAR = 4
          ELSEIF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
            IVAR = 3
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
            IVAR = 2
          ELSE
            IVAR = 1
          ENDIF
          IF( INDEX( UNTS(1:),'bd' ).NE.0 ) IVAR = -IVAR
          DO 420 NSL = 1,NSOLU
            IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SOLUT(NSL)(1:IDB) ) THEN
              ADDER = 0.D+0
              IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
                IUNM = -3
                IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
                  CALL RDINBP( C(1,NSL),ADDER,ICT(1,NSL),IVAR,UNTS )
                ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
                  CALL RDINAP( C(1,NSL),ADDER,ICT(1,NSL),IVAR,UNTS )
                ELSE
                  CALL RDINFP( C(1,NSL),VAR,ADDER,ICT(1,NSL),IVAR,UNTS )
                ENDIF
                CLOSE(UNIT=26)
              ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &          INDEX(ADUM(1:),'zonation').NE.0 )  THEN
                CALL RDINZP( C(1,NSL),VAR(1),ADDER,ICT(1,NSL),
     &            IVAR,IROCK )
              ELSE
                CALL RDINIP( C(1,NSL),VAR,ADDER,ICT(1,NSL),IVAR,IDOM )
              ENDIF
              GOTO 430
            ENDIF
  420     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute: ' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
  430     CONTINUE

        ELSEIF( INDEX(ADUM(1:),'specie').NE.0 ) THEN
          ADDER = 0.D+0
!
!---      Conservation- or kinetic-component species  ---
!
          IF( INDEX( BDUM(1:),'total_' ).NE.0 ) THEN
            DO 500 NSLX = NSOLU+1,NSOLU+NEQC+NEQK
              IDB = INDEX(SOLUT(NSLX)(1:),'  ') - 1
              IF( BDUM(1:NCHB).EQ.SOLUT(NSLX) ) THEN
                NSL = NSLX
                GOTO 540
              ENDIF
  500       CONTINUE
          ENDIF
!
!---      Aqueous reactive species  ---
!
          DO 510 NSPX = 1,NSPL
            IDB = INDEX(SPNML(NSPX)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SPNML(NSPX)(1:IDB) ) THEN
              NSP = NSPX
              GOTO 540
            ENDIF
  510     CONTINUE
!
!---      Solid reactive species  ---
!
          DO 520 NSPX = 1,NSPS
            IDB = INDEX(SPNMS(NSPX)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SPNMS(NSPX)(1:IDB) ) THEN
              NSP = NSPX + NSPL
!
!---          Verify that solid-species is not a mineral  ---
!
              IF( ISP_MN(NSP).EQ.1 ) THEN
                INDX = 4
                CHMSG = 'Solid-Species Mineral ' // 
     &            '(see Lithology Card): ' // BDUM(1:NCHB)
                CALL WRMSGS( INDX )
              ENDIF
              GOTO 540
            ENDIF
  520     CONTINUE
!
!---      Gas reactive species  ---
!
          DO 530 NSPX = 1,NSPG
            IDB = INDEX(SPNMG(NSPX)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SPNMG(NSPX)(1:IDB) ) THEN
              NSP = NSPX + NSPL + NSPS
              GOTO 540
            ENDIF
  530     CONTINUE
!
!---      pH  ---
!
          IF( BDUM(1:NCHB).EQ.'ph' .AND. ISPLK(1).NE.0 ) THEN
            NSP = MOD(ISPLK(1),1000)
            ISPLK(1) = ISPLK(1) + 1000
            IVAR = 2
            IF( IEO.EQ.2 ) IVAR = IVAR+10
            ADDER = 7.D+0
!
!---        Verify that species linked to pH is a conservation
!           component species  ---
!
            DO 532 NEQ = 1,NEQC
              IF( NSP.EQ.IEQ_C(2,NEQ) ) GOTO 540
  532       CONTINUE
            INDX = 4
            CHMSG = 'pH Species not a Conservation ' //
     &        'Component Species: ' // BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
          INDX = 4
          CHMSG = 'Unrecognized Reactive Species: ' // BDUM(1:NCHB)
          CALL WRMSGS( INDX )
  540     CONTINUE
!
!---      Conservation- or kinetic-component species  ---
!
          IF( INDEX( BDUM(1:),'total_' ).NE.0 ) THEN
            IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
              IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
                CALL RDINBP( C(1,NSL),ADDER,ICT(1,NSL),IVAR,UNTS )
              ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
                CALL RDINAP( C(1,NSL),ADDER,ICT(1,NSL),IVAR,UNTS )
              ELSE
                CALL RDINFP( C(1,NSL),VAR,ADDER,ICT(1,NSL),
     &            IVAR,UNTS )
              ENDIF
              CLOSE(UNIT=26)
            ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &        INDEX(ADUM(1:),'zonation').NE.0 )  THEN
              CALL RDINZP( C(1,NSL),VAR(1),ADDER,ICT(1,NSL),
     &          IVAR,IROCK )
            ELSE
              CALL RDINIP( C(1,NSL),VAR,ADDER,ICT(1,NSL),IVAR,
     &          IDOM )
            ENDIF
          ELSE
            IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
              IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
                CALL RDINBP( SP_C(1,NSP),ADDER,IC_SP(1,NSP),IVAR,UNTS )
              ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
                CALL RDINAP( SP_C(1,NSP),ADDER,IC_SP(1,NSP),IVAR,UNTS )
              ELSE
                CALL RDINFP( SP_C(1,NSP),VAR,ADDER,IC_SP(1,NSP),
     &            IVAR,UNTS )
              ENDIF
              CLOSE(UNIT=26)
            ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &        INDEX(ADUM(1:),'zonation').NE.0 )  THEN
              CALL RDINZP( SP_C(1,NSP),VAR(1),ADDER,IC_SP(1,NSP),
     &          IVAR,IROCK )
            ELSE
              CALL RDINIP( SP_C(1,NSP),VAR,ADDER,IC_SP(1,NSP),IVAR,
     &          IDOM )
            ENDIF
          ENDIF

        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
 1000 CONTINUE
 1002 CONTINUE
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINPT_GT
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
!     Read input file cards.
!     Direct control to card reader subroutines.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINPT_GT'
!
!---  Write header line to output file  ---
!
      WRITE(IWR,'(/,A)') ' --- Input File Record ---'
!
!---  Search input file for simulation title card  ---
!
  100 CONTINUE
  109 READ(IRD,'(A)', END=110) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 109
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'simulation').NE.0 ) THEN
        CALL RDSIMU
        REWIND(IRD)
        GOTO 200
      ELSE
        GOTO 100
      ENDIF
  110 CONTINUE
      INDX = 18
      CHMSG = 'Missing Simulation Title Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for solution control card  ---
!
  200 CONTINUE
  209 READ(IRD,'(A)', END=210) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 209
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'solution').NE.0 ) THEN
        CALL RDSOLU
        IF( IOM.NE.ICODE ) THEN
          INDX = 18
          CHMSG = 'Incompatible Operational Mode'
          CALL WRMSGS( INDX )
        ENDIF
        REWIND(IRD)
        GOTO 300
      ELSE
        GOTO 200
      ENDIF
  210 CONTINUE
      INDX = 18
      CHMSG = 'Missing Solution Control Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for grid card  ---
!
  300 CONTINUE
  302 READ(IRD,'(A)', END=304) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 302
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'grid').NE.0 .AND.
     &  INDEX(CHDUM(2:),'rotation').EQ.0 ) THEN
        CALL RDGRID
        REWIND(IRD)
        GOTO 400
      ELSE
        GOTO 302
      ENDIF
  304 CONTINUE
      INDX = 18
      CHMSG = 'Missing Grid Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for grid rotation card --
!
  400 CONTINUE
  402 READ(IRD,'(A)', END=404) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 402
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'grid').NE.0 .AND.
     &  INDEX(CHDUM(2:),'rotation').NE.0 ) THEN
        CALL RDROTATION
        REWIND(IRD)
        GOTO 410
      ELSE
        GOTO 402
      ENDIF
  404 CONTINUE
      INDX = 1
      CHMSG = 'Missing Grid Rotation Card'
      CARD = 'Grid Rotation Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for rock/soil zonation card  ---
!
  410 CONTINUE
  412 READ(IRD,'(A)', END=414) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 412
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'rock/soil').NE.0 ) THEN
        CALL RDROCK
        REWIND(IRD)
        GOTO 500
      ELSE
        GOTO 412
      ENDIF
  414 CONTINUE
      INDX = 18
      CHMSG = 'Missing Rock/Soil Zonation Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for inactive nodes card  ---
!
  500 CONTINUE
  502 READ(IRD,'(A)', END=504) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 502
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'inactive').NE.0 .AND.
     &    INDEX(CHDUM(2:),'triangle').EQ.0 ) THEN
        CALL RDINAC
        REWIND(IRD)
        GOTO 520
      ELSE
        GOTO 502
      ENDIF
  504 CONTINUE
      INDX = 1
      CHMSG = 'Missing Inactive Nodes Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for fracture geometry card --
!
  520 CONTINUE
  529 READ(IRD,'(A)', END=530) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 529
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'fracture').NE.0 .AND.
     &  INDEX(CHDUM(2:),'geometry').NE.0 ) THEN
        INDX = 1
        CALL RDGEOM_FRC( INDX )
        REWIND(IRD)
        GOTO 540
       ELSE
        GOTO 520
      ENDIF
  530 CONTINUE
!
!---  Fracture flow simulations  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) THEN
        INDX = 4
        CHMSG = 'Missing Fracture Geometry Card'
        CALL WRMSGS( INDX )
      ELSE
        INDX = 1
        CHMSG = 'Missing Fracture Geometry Card'
        CALL WRMSGS( INDX )
        REWIND(IRD)
      ENDIF
!
!---  Search input file for inactive triangles card  ---
!
  540 CONTINUE
  542 READ(IRD,'(A)', END=544) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 542
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'inactive').NE.0 .AND.
     &    INDEX(CHDUM(2:),'triangle').NE.0 ) THEN
        INDX = 1
        CALL RDINAC_FRC( INDX )
        REWIND(IRD)
        GOTO 550
      ELSE
        GOTO 542
      ENDIF
  544 CONTINUE
      INDX = 1
      CHMSG = 'Missing Inactive Fracture Triangles Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for solute/fluid interaction card --
!
  550 CONTINUE
  552 READ(IRD,'(A)', END=558) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 552
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'solute/fluid').NE.0 ) THEN
        CALL RDTF_GT
        REWIND(IRD)
        GOTO 560
      ELSE
        GOTO 550
      ENDIF
  558 CONTINUE
      IF( IEQC.EQ.0 ) THEN
        REWIND(IRD)
      ELSE
        INDX = 18
        CHMSG = 'Missing Solute/Fluid Interaction Card'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Search input file for borehole card --
!
  560 CONTINUE
  562 READ(IRD,'(A)', END= 568) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 562
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'borehole').NE.0 .AND.
     &  INDEX(CHDUM(2:),'initial').EQ.0 ) THEN
        CALL RD_BH_GT
        REWIND(IRD)
        GOTO 600
      ELSE
        GOTO 560
      ENDIF
  568 CONTINUE
      INDX = 1
      CHMSG = 'Missing Borehole Card'
      CARD = 'Borehole Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for mechanical properties card  ---
!
  600 CONTINUE
  609 READ(IRD,'(A)', END=610) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 609
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'mechanical').NE.0 ) THEN
        CALL RDMECH
        REWIND(IRD)
        GOTO 700
      ELSE
        GOTO 600
      ENDIF
  610 CONTINUE
      INDX = 18
      CHMSG = 'Missing Mechanical Properties Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for hydraulic properties card  ---
!
  700 CONTINUE
  709 READ(IRD,'(A)', END=710) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 709
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'hydraulic').NE.0 ) THEN
        CALL RDHYDR
        REWIND(IRD)
        GOTO 800
      ELSE
        GOTO 700
      ENDIF
  710 CONTINUE
      INDX = 18
      CHMSG = 'Missing Hydraulic Properties Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for saturation function card  ---
!
  800 CONTINUE
  809 READ(IRD,'(A)', END=810) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 809
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'saturation').NE.0 ) THEN
        CALL RDSP_GT
        REWIND(IRD)
        GOTO 900
      ELSE
        GOTO 800
      ENDIF
  810 CONTINUE
      INDX = 18
      CHMSG = 'Missing Saturation Function Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for aqueous relative permeability cards  ---
!
  900 CONTINUE
  909 READ(IRD,'(A)', END=910) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 909
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'aqueous rel').NE.0 .AND.
     &    INDEX(CHDUM(2:),'x-aqueous rel').EQ.0 .AND.
     &    INDEX(CHDUM(2:),'y-aqueous rel').EQ.0 .AND.
     &    INDEX(CHDUM(2:),'z-aqueous rel').EQ.0 ) THEN
        CALL RDLRP
        IRPLX = 1
        IRPLY = 1
        IRPLZ = 1
        REWIND(IRD)
        GOTO 920
      ELSE
        GOTO 900
      ENDIF
  910 CONTINUE
      REWIND( IRD )
  920 CONTINUE
  929 READ(IRD,'(A)', END=930) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 929
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'x-aqueous rel').NE.0 ) THEN
        ITX = 1
        CALL RDLRPT( ITX )
        IRPLX = 1
        REWIND(IRD)
        GOTO 940
      ELSE
        GOTO 920
      ENDIF
  930 CONTINUE
      REWIND( IRD )
  940 CONTINUE
  949 READ(IRD,'(A)', END=950) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 949
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'y-aqueous rel').NE.0 ) THEN
        ITX = 2
        CALL RDLRPT( ITX )
        IRPLY = 1
        REWIND(IRD)
        GOTO 960
      ELSE
        GOTO 940
      ENDIF
  950 CONTINUE
      REWIND( IRD )
  960 CONTINUE
  969 READ(IRD,'(A)', END=970) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 969
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'z-aqueous rel').NE.0 ) THEN
        ITX = 3
        CALL RDLRPT( ITX )
        IRPLZ = 1
        REWIND(IRD)
        GOTO 970
      ELSE
        GOTO 960
      ENDIF
  970 CONTINUE
      REWIND( IRD )
      IF( IRPLX.EQ.0 .AND. IRPLY.EQ.0 .AND. IRPLZ.EQ.0 ) THEN
        INDX = 18
        CHMSG = 'Missing Aqueous Relative Permeability Card'
        CALL WRMSGS( INDX )
      ELSE
        IF( IRPLX.EQ.0 .AND. IFLD.GT.1 ) THEN
          INDX = 18
          CHMSG = 'Missing X-Aqueous Relative Permeability Card'
          CALL WRMSGS( INDX )
        ELSEIF( IRPLY.EQ.0 .AND. JFLD.GT.1 ) THEN
          INDX = 18
          CHMSG = 'Missing Y-Aqueous Relative Permeability Card'
          CALL WRMSGS( INDX )
        ELSEIF( IRPLZ.EQ.0 .AND. KFLD.GT.1 ) THEN
          INDX = 18
          CHMSG = 'Missing Z-Aqueous Relative Permeability Card'
          CALL WRMSGS( INDX )
        ENDIF
        IF( IRPLX.EQ.0 .AND. IFLD.EQ.1 ) THEN
          INDX = 1
          CHMSG = 'Missing X-Aqueous Relative Permeability Card'
          CALL WRMSGS( INDX )
        ENDIF
        IF( IRPLY.EQ.0 .AND. JFLD.EQ.1 ) THEN
          INDX = 1
          CHMSG = 'Missing Y-Aqueous Relative Permeability Card'
          CALL WRMSGS( INDX )
        ENDIF
        IF( IRPLZ.EQ.0 .AND. KFLD.EQ.1 ) THEN
          INDX = 1
          CHMSG = 'Missing Z-Aqueous Relative Permeability Card'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Search input file for scaling card  ---
!
 1000 CONTINUE
      IF( ISLC(19).EQ.1 ) THEN
 1009   READ(IRD,'(A)', END=1010) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 1009
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &      INDEX(CHDUM(2:),'scaling').NE.0 ) THEN
          CALL RDSCLF_GT
          REWIND(IRD)
          GOTO 1100
        ELSE
          GOTO 1000
        ENDIF
 1010   CONTINUE
        INDX = 18
        CHMSG = 'Missing Scaling Card'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Search input file for gas relative permeability card  ---
!
 1100 CONTINUE
 1109 READ(IRD,'(A)', END=1110) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 1109
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'gas rel').NE.0 ) THEN
        CALL RDGRP
        REWIND(IRD)
        GOTO 1200
      ELSE
        GOTO 1100
      ENDIF
 1110 CONTINUE
      INDX = 18
      CHMSG = 'Missing Gas Relative Permeability Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for thermal properties card  ---
!
 1200 CONTINUE
 1191 READ(IRD,'(A)', END=1250) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 1191
      CALL LCASE( CHDUM )
      IF( CHDUM(1:8).EQ.'~thermal' ) THEN
        CALL RDTHER
        REWIND(IRD)
        GOTO 1300
      ELSE
        GOTO 1200
      ENDIF
 1250 CONTINUE
      INDX = 18
      CHMSG = 'Missing Thermal Properties Card'
      CALL WRMSGS( INDX )
 1300 CONTINUE
!
!---  Search input file for solute/rock interaction card --
!
 2100 CONTINUE
 2109 READ(IRD,'(A)', END=2110) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2109
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'solute/porous').NE.0 ) THEN
        CALL RDTP_GT
        REWIND(IRD)
        GOTO 2200
      ELSE
        GOTO 2100
      ENDIF
 2110 CONTINUE
      IF( IEQC.EQ.0 ) THEN
       REWIND(IRD)
      ELSE
        INDX = 18
        CHMSG = 'Missing Solute/Porous Media Interaction Card'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Skip reaction and equation cards, no
!     reactive transport  ---
!
 2200 CONTINUE

!      IF( ISLC(40).EQ.0 ) GOTO 3000
!
!---  Search input file for aqueous species card  ---
!
 2209 READ(IRD,'(A)', END=2210) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2209
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'aqueous').NE.0 .AND.
     &  INDEX(CHDUM(2:),'specie').NE.0) THEN
        CALL RDAQSP
        REWIND(IRD)
        GOTO 2300
      ELSE
        GOTO 2209
      ENDIF
 2210 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Aqueous Species Card'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Search input file for solid species card  ---
!
 2300 CONTINUE
 2309 READ(IRD,'(A)', END=2310) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2309
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'solid').NE.0 .AND.
     &  INDEX(CHDUM(2:),'specie').NE.0) THEN
        CALL RDSDSP
        REWIND(IRD)
        GOTO 2400
      ELSE
        GOTO 2300
      ENDIF
 2310 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Solid Species Card'
        CALL WRMSGS( INDX )
      ENDIF
      REWIND(IRD)
 2400 CONTINUE
!
!---  Search input file for exchanged species card  ---
!
 2409 READ(IRD,'(A)', END=2410) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2409
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'excha').NE.0 .AND.
     &  INDEX(CHDUM(2:),'specie').NE.0) THEN
        CALL RDEXSP
        REWIND(IRD)
        GOTO 2420
      ELSE
        GOTO 2400
      ENDIF
 2410 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Exchanged Species Card'
        CALL WRMSGS( INDX )
      ENDIF
      REWIND(IRD)
 2420 CONTINUE
!
!---  Search input file for equilibrium reactions card  ---
!
 2429 READ(IRD,'(A)', END=2430) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2429
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'equil').NE.0 .AND.
     &  INDEX(CHDUM(2:),'react').NE.0) THEN
        CALL RDEQRC
        REWIND(IRD)
        GOTO 2500
      ELSE
        GOTO 2420
      ENDIF
 2430 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Equilibrium Reactions Card'
        CALL WRMSGS( INDX )
      ENDIF
      REWIND(IRD)
!
!---  Search input file for kinetic reactions card  ---
!
 2500 CONTINUE
 2509 READ(IRD,'(A)', END=2510) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2509
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'kinetic').NE.0 .AND.
     &  INDEX(CHDUM(2:),'react').NE.0) THEN
        CALL RDKNRC
        REWIND(IRD)
        GOTO 2700
      ELSE
        GOTO 2500
      ENDIF
 2510 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Kinetic Reactions Card'
        CALL WRMSGS( INDX )
      ENDIF
      REWIND(IRD)
!
!---  Search input file for equilibrium equation card  ---
!
 2700 CONTINUE
 2709 READ(IRD,'(A)', END=2710) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2709
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'equil').NE.0 .AND.
     &  INDEX(CHDUM(2:),'equat').NE.0) THEN
        CALL RDEQEQ
        REWIND(IRD)
        GOTO 2800
      ELSE
        GOTO 2700
      ENDIF
 2710 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Equilibrium Equations Card'
        CALL WRMSGS( INDX )
      ENDIF
      REWIND(IRD)
!
!---  Search input file for conservation equations card  ---
!
 2800 CONTINUE
 2809 READ(IRD,'(A)', END=2810) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2809
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'conservation').NE.0 .AND.
     &  INDEX(CHDUM(2:),'equat').NE.0) THEN
        CALL RDCNEQ
        REWIND(IRD)
        GOTO 2900
      ELSE
        GOTO 2800
      ENDIF
 2810 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Conservation Equations Card'
        CALL WRMSGS( INDX )
      ENDIF
      REWIND(IRD)
!
!---  Search input file for kinetic equations card  ---
!
 2900 CONTINUE
 2909 READ(IRD,'(A)', END=2910) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2909
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'kinetic').NE.0 .AND.
     &  INDEX(CHDUM(2:),'equat').NE.0) THEN
        CALL RDKNEQ
        REWIND(IRD)
        GOTO 3000
      ELSE
        GOTO 2900
      ENDIF
 2910 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Kinetic Equations Card'
        CALL WRMSGS( INDX )
      ENDIF
      REWIND(IRD)
!
!---  Search input file for lithology card  ---
!
 3000 CONTINUE
 3009 READ(IRD,'(A)', END=3010) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3009
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'lithol').NE.0) THEN
        CALL RDLITH
        REWIND(IRD)
        GOTO 3100
      ELSE
        GOTO 3000
      ENDIF
 3010 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Lithology Card'
        CALL WRMSGS( INDX )
      ENDIF
      REWIND(IRD)
!
!---  Search input file for reactive species link card  ---
!
 3100 CONTINUE
 3109 READ(IRD,'(A)', END=3110) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3109
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'link').NE.0 .AND.
     &  INDEX(CHDUM(2:),'specie').NE.0) THEN
        CALL RDSPLK
        REWIND(IRD)
        GOTO 3200
      ELSE
        GOTO 3100
      ENDIF
 3110 CONTINUE
      REWIND(IRD)
      IF( ISLC(40).GT.0 ) THEN
        INDX = 1
        CHMSG = 'Missing Reactive Species Link Card'
        CALL WRMSGS( INDX )
      ENDIF
      REWIND(IRD)

!
!---  Search input file for salt transport card --
!
 3200 CONTINUE
 3209 READ(IRD,'(A)', END=3210) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3209
      CALL LCASE( CHDUM )
      IF( CHDUM(1:15).EQ.'~salt transport' ) THEN
        CALL RDST_GT
        REWIND(IRD)
        GOTO 3250
      ELSE
        GOTO 3200
      ENDIF
 3210 CONTINUE
      INDX = 1
      CHMSG = 'Missing Salt Transport Card'
      CARD = 'Salt Transport Card Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for ground-loop well card --
!
 3250 CONTINUE
 3259 READ(IRD,'(A)', END= 3260) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3259
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'ground').NE.0 .AND.
     &  INDEX(CHDUM(2:),'loop').NE.0 .AND.
     &  INDEX(CHDUM(2:),'well').NE.0 ) THEN
        CALL RD_GRLP_WELL
        REWIND(IRD)
        GOTO 3300
      ELSE
        GOTO 3250
      ENDIF
 3260 CONTINUE
      INDX = 1
      CHMSG = 'Missing Ground-Loop Well Card'
      CARD = 'Ground-Loop Well Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for thermo-catalytic card --
!
 3300 CONTINUE
 3309 READ(IRD,'(A)', END= 3310) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 3309
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'thermo').NE.0 .AND.
     &  INDEX(CHDUM(2:),'catalytic').NE.0 ) THEN
        CALL RDTHMOC_GT
        REWIND(IRD)
        GOTO 4000
      ELSE
        GOTO 3300
      ENDIF
 3310 CONTINUE
      INDX = 1
      CHMSG = 'Missing Thermo-catalytic Card'
      CARD = 'Thermo-catalytic Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for plant card --
!
 4000 CONTINUE
 4009 READ(IRD,'(A)', END= 4010) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4009
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'plant').NE.0 ) THEN
        CALL RDPLANT
        REWIND(IRD)
        GOTO 4100
      ELSE
        GOTO 4000
      ENDIF
 4010 CONTINUE
      INDX = 1
      CHMSG = 'Missing Plant Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for surface cover card --
!
 4100 CONTINUE
 4109 READ(IRD,'(A)', END= 4110) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4109
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'cover').NE.0 ) THEN
        CALL RDSFCOV
        REWIND(IRD)
        GOTO 4200
      ELSE
        GOTO 4100
      ENDIF
 4110 CONTINUE
      INDX = 1
      CHMSG = 'Missing Surface Cover Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for atmospheric conditions card --
!
 4200 CONTINUE
 4209 READ(IRD,'(A)', END= 4210) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4209
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'atmospheric').NE.0 ) THEN
        CALL RDATMOS
        REWIND(IRD)
        GOTO 4300
      ELSE
        GOTO 4200
      ENDIF
 4210 CONTINUE
      IF( NSFCN.GT.0 ) THEN
        INDX = 18
        CHMSG = 'Missing Atmospheric Conditions Card'
        CALL WRMSGS( INDX )
      ELSE
        INDX = 1
        CHMSG = 'Missing Atmospheric Conditions Card'
        CALL WRMSGS( INDX )
        REWIND(IRD)
      ENDIF
!
!---  Search input file for source card --
!
 4300 CONTINUE
 4309 READ(IRD,'(A)', END= 4310) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4309
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'source').NE.0 .AND.
     &    INDEX(CHDUM(2:),'fracture').EQ.0 .AND.
     &    INDEX(CHDUM(2:),'borehole').EQ.0 ) THEN
        CALL RDSR_GT
        REWIND(IRD)
        GOTO 4320
      ELSE
        GOTO 4300
      ENDIF
 4310 CONTINUE
      INDX = 1
      CHMSG = 'Missing Source Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for geomechanics property card --
!
 4320 CONTINUE
 4329 READ(IRD,'(A)', END=4330) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4329
      CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'geomech').NE.0 .AND.
     &    INDEX(CHDUM(2:),'prop').NE.0 ) THEN
        CALL RDGMP
        REWIND(IRD)
        GOTO 4400
      ELSE
        GOTO 4320
      ENDIF
 4330 CONTINUE
!
!---  Geomechanical simulations  ---
!
      IF( ISLC(50).NE.0 ) THEN
        INDX = 4
        CHMSG = 'Missing Geomechanical Properties Card'
        CALL WRMSGS( INDX )
      ELSE
        REWIND(IRD)
      ENDIF
!
!---  Search input file for boundary conditions card --
!
 4400 CONTINUE
 4409 READ(IRD,'(A)', END=4410) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4409
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'boundary').NE.0 .AND.
     &  INDEX(CHDUM(2:),'initial').EQ.0 .AND.
     &  INDEX(CHDUM(2:),'geomech').EQ.0 ) THEN
        CALL RDBC_GT
        REWIND(IRD)
        GOTO 4500
      ELSE
        GOTO 4400
      ENDIF
 4410 CONTINUE
      INDX = 1
      CHMSG = 'Missing Boundary Conditions Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for initial boundary conditions card --
!
 4500 CONTINUE
 4509 READ(IRD,'(A)', END=4510) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4509
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'initial').NE.0 .AND.
     &    INDEX(CHDUM(2:),'boundary').NE.0 .AND.
     &    INDEX(CHDUM(2:),'fracture').EQ.0 .AND.
     &    INDEX(CHDUM(2:),'borehole').EQ.0  ) THEN
        CALL RDIBC_GT
        REWIND(IRD)
        GOTO 4550
      ELSE
        GOTO 4500
      ENDIF
 4510 CONTINUE
      INDX = 1
      CHMSG = 'Missing Initial Boundary Conditions Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for initial conditions card --
!
 4550 CONTINUE
 4559 READ(IRD,'(A)', END=4560) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4559
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'initial').NE.0 .AND.
     &    INDEX(CHDUM(2:),'boundary').EQ.0 .AND.
     &    INDEX(CHDUM(2:),'fracture').EQ.0 .AND.
     &    INDEX(CHDUM(2:),'borehole').EQ.0  ) THEN
        CALL RDIC_GT
        REWIND(IRD)
        GOTO 4600
      ELSE
        GOTO 4550
      ENDIF
 4560 CONTINUE
      IF( IEO.EQ.2 ) THEN
        INDX = 1
        CHMSG = 'Missing Initial Conditions Card'
        CALL WRMSGS( INDX )
        INDX = 2
        CALL RDRST(INDX)
        ISIC = 3
        REWIND(IRD)
      ELSE
        INDX = 18
        CHMSG = 'Missing Initial Conditions Card'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Search input file for surface flux card --
!
 4600 CONTINUE
 4609 READ(IRD,'(A)', END=4610) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4609
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'surface').NE.0 .AND.
     &  INDEX(CHDUM(2:),'flux').NE.0 ) THEN
        CALL RDSF_GT
        REWIND(IRD)
        GOTO 4700
      ELSE
        GOTO 4600
      ENDIF
 4610 CONTINUE
      INDX = 1
      CHMSG = 'Missing Surface Flux Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for observed data card --
!
 4700 CONTINUE
      IF( ISLC(20).EQ.1 ) THEN
 4709   READ(IRD,'(A)', END=4710) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 4709
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &      INDEX(CHDUM(2:),'observed').NE.0 ) THEN
          CALL RDOBDA
          REWIND(IRD)
          GOTO 4800
        ELSE
          GOTO 4700
        ENDIF
 4710   CONTINUE
        IF( ISLC(20).EQ.0 ) THEN
          REWIND(IRD)
        ELSE
          INDX = 18
          CHMSG = 'Missing Observed Data Card'
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Search input file for UCode control Card --
!
 4800   CONTINUE
 4809   READ(IRD,'(A)', END=4810) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 4809
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &      INDEX(CHDUM(2:),'ucode').NE.0 ) THEN
          CALL RDUCODE
          REWIND(IRD)
          GOTO 4900
        ELSE
          GOTO 4800
        ENDIF
 4810   CONTINUE
        IF( ISLC(20).EQ.0 ) THEN
          REWIND(IRD)
        ELSE
          INDX = 18
          CHMSG = 'Missing UCode Control Card'
          CALL WRMSGS( INDX )
        ENDIF
 4900   CONTINUE
      ENDIF
!
!---  Geomechanics  ---
!
      IF( ISLC(50).NE.0 ) THEN
 7100   CONTINUE
!
!---    Search input file for inactive nodes card  ---
!
 7109   READ(IRD,'(A)', END=7110) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 7109
        CALL LCASE( CHDUM )
        IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'inactive').NE.0 .AND.
     &    INDEX(CHDUM(2:),'elements').NE.0 ) THEN
          CALL RDINAC_GM
          REWIND(IRD)
          GOTO 7200
        ELSE
          GOTO 7100
        ENDIF
 7110   CONTINUE
        INDX = 1
        CHMSG = 'Missing Inactive Elements Card'
        CARD = 'Inactive Elements Card'
        CALL WRMSGS( INDX )
        REWIND(IRD)
 7200   CONTINUE
!
!---    Search input file for geomechanics link card --
!
 7309   READ(IRD,'(A)', END=7310) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 7309
        CALL LCASE( CHDUM )
          IF( CHDUM(1:1).EQ.'~' .AND.
     &      INDEX(CHDUM(2:),'geomech').NE.0 .AND.
     &      INDEX(CHDUM(2:),'link').NE.0 ) THEN
          CALL RDGMLK
          REWIND(IRD)
          GOTO 7400
        ELSE
          GOTO 7200
        ENDIF
 7310   CONTINUE
!
!---    Geomechanical simulations  ---
!
        IF( ISLC(50).NE.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Geomechanics Link Card'
          CALL WRMSGS( INDX )
          REWIND(IRD)
        ELSE
          REWIND(IRD)
        ENDIF
 7400   CONTINUE
!
!---    Search input file for geomechanics boundary condition card --
!
 7409   READ(IRD,'(A)', END=7410) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 7409
        CALL LCASE( CHDUM )
          IF( CHDUM(1:1).EQ.'~' .AND.
     &      INDEX(CHDUM(2:),'geomech').NE.0 .AND.
     &      INDEX(CHDUM(2:),'bound').NE.0 ) THEN
          CALL RDGMBC
          REWIND(IRD)
          GOTO 7500
        ELSE
          GOTO 7400
        ENDIF
 7410   CONTINUE
!
!---    Geomechanical simulations  ---
!
        IF( ISLC(50).NE.0 ) THEN
          INDX = 1
          CHMSG = 'Missing Geomechanics Boundary Condition Card'
          CALL WRMSGS( INDX )
        ELSE
          REWIND(IRD)
        ENDIF
 7500   CONTINUE
      ENDIF
!
!---  Create a node connection map  ---
!
      CALL CONNMAP
!
!---  Search input file for fracture initial conditions card --
!
 5500 CONTINUE
 5509 READ(IRD,'(A)', END=5510) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 5509
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'fracture').NE.0 .AND.
     &  INDEX(CHDUM(2:),'initial').NE.0 ) THEN
        CALL RDIC_FRC_GT
        REWIND(IRD)
        GOTO 5600
       ELSE
        GOTO 5500
      ENDIF
 5510 CONTINUE
      INDX = 1
      CHMSG = 'Missing Fracture Initial Conditions Card'
      CARD = 'Fracture Initial Conditions Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for fracture properties card --
!
 5600 CONTINUE
 5609 READ(IRD,'(A)', END=5610) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 5609
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'fracture').NE.0 .AND.
     &  INDEX(CHDUM(2:),'prop').NE.0 ) THEN
        CALL RDPROP_FRC_GT
        REWIND(IRD)
        GOTO 5650
       ELSE
        GOTO 5600
      ENDIF
 5610 CONTINUE
!
!---  Fracture flow simulations  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) THEN
        INDX = 4
        CARD = 'Fracture Properties Card'
        CHMSG = 'Missing Fracture Properties Card'
        CALL WRMSGS( INDX )
      ELSE
        INDX = 1
        CARD = 'Fracture Properties Card'
        CHMSG = 'Missing Fracture Properties Card'
        CALL WRMSGS( INDX )
        REWIND(IRD)
      ENDIF
!
!---  Search input file for borehole initial conditions card --
!
 5650 CONTINUE
 5659 READ(IRD,'(A)', END=5660) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 5659
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'borehole').NE.0 .AND.
     &  INDEX(CHDUM(2:),'initial').NE.0 ) THEN
        CALL RDIC_BH_GT
        REWIND(IRD)
        GOTO 5700
       ELSE
        GOTO 5650
      ENDIF
 5660 CONTINUE
      INDX = 1
      CHMSG = 'Missing Borehole Initial Conditions Card'
      CARD = 'Borehole Initial Conditions Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for output control card --
!
 5700 CONTINUE
 5709 READ(IRD,'(A)', END=5710) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 5709
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'output').NE.0 ) THEN
        CALL RDOU_GT
        REWIND(IRD)
        GOTO 5800
      ELSE
        GOTO 5700
      ENDIF
 5710 CONTINUE
      INDX = 1
      CHMSG = 'Missing Output Control Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for fracture source card --
!
 5800 CONTINUE
 5809 READ(IRD,'(A)', END=5810) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 5809
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'fracture').NE.0 .AND.
     &  INDEX(CHDUM(2:),'source').NE.0 ) THEN
        CALL RDSR_FRC_GT
        REWIND(IRD)
        GOTO 5900
       ELSE
        GOTO 5800
      ENDIF
 5810 CONTINUE
      INDX = 1
      CHMSG = 'Missing Fracture Source Card'
      CARD = 'Fracture Source Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for 4850 Drift stress card --
!
 5900 CONTINUE
 5909 READ(IRD,'(A)', END=5910) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 5909
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'4850').NE.0 .AND.
     &  INDEX(CHDUM(2:),'drift').NE.0 .AND.
     &  INDEX(CHDUM(2:),'stress').NE.0 ) THEN
        CALL RD4850D_SIG
        REWIND(IRD)
        GOTO 6000
       ELSE
        GOTO 5900
      ENDIF
 5910 CONTINUE
      INDX = 1
      CHMSG = 'Missing 4850 Drift Stress Card'
      CARD = '4850 Drift Stress Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
 6000 CONTINUE
!
!---  End of input record --
!
      CARD = 'End of Input Record'
      ICD = INDEX( CARD,'  ' )-1
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINPT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDOU_GT
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
!     Read input file for output information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PLT_ATM
      USE PARM_FRC
      USE PARM_BH
      USE OUTPU
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FILES
      USE FDVS
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      EXTERNAL ICOUNT
!
!----------------------Type Declarations-------------------------------!
!

      CHARACTER*64 SPNM

      CHARACTER*64 ADUM,UNTS,SOLNM,PLNTNM
      CHARACTER*512 CHDUM
      CHARACTER*6 FORM
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM
      DATA FORM / '(I6,$)' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDOU_GT'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Output Control Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read reference node information  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Reference Nodes'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NREF)
      IF( NREF.GT.LREF ) THEN
        INDX = 5
        CHMSG = 'Number of Reference Nodes > LREF'
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(IWR,'(/,A,I6)') 'Reference Node No. and Indices: ',NREF
      DO 100 N = 1,NREF
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Reference Node Index'
        CALL RDINT(ISTART,ICOMMA,CHDUM,IRF)
        CALL RDINT(ISTART,ICOMMA,CHDUM,JRF)
        CALL RDINT(ISTART,ICOMMA,CHDUM,KRF)
        IF( IRF.LT.1 .OR. IRF.GT.IFLD .OR. JRF.LT.1 .OR.
     &    JRF.GT.JFLD. OR. KRF.LT.1 .OR. KRF.GT.KFLD) THEN
          INDX = 4
          CHMSG = 'Out-of-Range Reference Node Index'
          CALL WRMSGS( INDX )
        ENDIF
        NDREF(N) = ND(IRF,JRF,KRF)
        WRITE(FORM(3:3),'(I1)') ICOUNT(NDREF(N))
        WRITE(IWR,'(2X,A,$)') 'Reference Node No. '
        WRITE(IWR,FORM) NDREF(N)
        WRITE(FORM(3:3),'(I1)') ICOUNT(IRF)
        WRITE(IWR,'(2X,A,$)') 'I = '
        WRITE(IWR,FORM) IRF
        WRITE(FORM(3:3),'(I1)') ICOUNT(JRF)
        WRITE(IWR,'(2X,A,$)') 'J = '
        WRITE(IWR,FORM) JRF
        WRITE(FORM(3:3),'(I1)') ICOUNT(KRF)
        WRITE(IWR,'(2X,A,$)') 'K = '
        WRITE(IWR,FORM) KRF
        WRITE(IWR,'(2X,A)' ) 'Indices'
  100 CONTINUE
!
!---  Read fracture reference node information  ---
!
      IF( (ISLC(74).EQ.1 .OR. ISLC(74).EQ.3) .AND. NF_FRC.GT.0 ) THEN
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Number of Fracture Reference Nodes'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NREF_FRC)
        IF( NREF_FRC.GT.LREF ) THEN
          INDX = 5
          CHMSG = 'Number of Fracture Reference Nodes > LREF'
          CALL WRMSGS( INDX )
        ENDIF
        WRITE(IWR,'(/,A,I6)') 'Fracture Nos., Local Triangle Nos. ' //
     &    'and Global Triangle Nos.: ',NREF_FRC
        DO N = 1,NREF_FRC
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Fracture Reference Node: Fracture No.'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NFX)
          IF( NFX.LT.1 .OR. NFX.GT.NF_FRC ) THEN
            INDX = 4
            CHMSG = 'Out-of-Range Fracture Reference Node Fracture'
            CALL WRMSGS( INDX )
          ENDIF
          VARB = 'Fracture Reference Node: Local Triangle No.'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NTX)
          IF( NTX.LT.1 .OR. NTX.GT.NTP_FRC(NFX) ) THEN
            INDX = 4
            CHMSG = 'Out-of-Range Fracture Reference Node ' // 
     &        'Triangle No.'
            CALL WRMSGS( INDX )
          ENDIF
          NC = 0
          L1: DO KFX = 1,NFX
            DO KTX = 1,NTP_FRC(KFX)
                NC = NC + 1
              IF( KFX.EQ.NFX .AND. KTX.EQ.NTX ) EXIT L1
            ENDDO
          ENDDO L1
          NDREF_FRC(N) = NC
          WRITE(FORM(3:3),'(I1)') ICOUNT(NFX)
          WRITE(IWR,'(2X,A,$)') 'Fracture No. = '
          WRITE(IWR,FORM) NFX
          WRITE(FORM(3:3),'(I1)') ICOUNT(NTX)
          WRITE(IWR,'(2X,A,$)') 'Local Triangle No. = '
          WRITE(IWR,FORM) NTX
          WRITE(FORM(3:3),'(I1)') ICOUNT(NDREF_FRC(N))
          WRITE(IWR,'(2X,A,$)') 'Global Triangle No. = '
          WRITE(IWR,FORM) NDREF_FRC(N)
          WRITE(IWR,'(A)' ) ''
        ENDDO
      ENDIF
!
!---  Read borehole reference node information  ---
!
      IF( (ISLC(74).EQ.2 .OR. ISLC(74).EQ.3) .AND. N_BH.GT.0 ) THEN
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Number of Borehole Reference Nodes'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NREF_BH)
        IF( NREF_BH.GT.LREF ) THEN
          INDX = 5
          CHMSG = 'Number of Borehole Reference Nodes > LREF'
          CALL WRMSGS( INDX )
        ENDIF
        WRITE(IWR,'(/,A,I6)') 'Borehole Nos., Local Borehole Nos. ' //
     &    'and Global Borehole Nos.: ',NREF_BH
        DO N = 1,NREF_BH
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Borehole Reference Node: Borehole No.'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NBH)
          IF( NBH.LT.1 .OR. NBH.GT.N_BH ) THEN
            INDX = 4
            CHMSG = 'Out-of-Range Borehole No.'
            CALL WRMSGS( INDX )
          ENDIF
          VARB = 'Borehole Reference Node: Local Borehole Node No.'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NBN)
          IF( NBN.LT.1 .OR. NBN.GT.(ID_BH(4,NBH)-ID_BH(3,NBH)+1) ) THEN
            INDX = 4
            CHMSG = 'Out-of-Range Local Borehole Node No.'
            CALL WRMSGS( INDX )
          ENDIF
          L2: DO KHX = 1,N_BH
            DO KNX = ID_BH(3,NBH),ID_BH(4,NBH)
              IF( KHX.EQ.NBH .AND. KNX-ID_BH(3,NBH)+1.EQ.NBN ) EXIT L2
            ENDDO
          ENDDO L2
          NDREF_BH(N) = KNX
          WRITE(FORM(3:3),'(I1)') ICOUNT(NBH)
          WRITE(IWR,'(2X,A,$)') 'Borehole No. = '
          WRITE(IWR,FORM) NBH
          WRITE(FORM(3:3),'(I1)') ICOUNT(NBN)
          WRITE(IWR,'(2X,A,$)') 'Local Borehole Node No. = '
          WRITE(IWR,FORM) NBN
          WRITE(FORM(3:3),'(I1)') ICOUNT(NDREF_BH(N))
          WRITE(IWR,'(2X,A,$)') 'Global Borehole Node No. = '
          WRITE(IWR,FORM) NDREF_BH(N)
          WRITE(IWR,'(A)' ) ''
        ENDDO
      ENDIF
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      IDFLT = 1
      IFQS = IBIG
      VARB = 'Reference Node Screen Output Frequency'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IFQS)
      WRITE(IWR,'(/,2A,I6,A)') VARB(1:IVR),': Every ',IFQS,
     &' Time Step(s)'
      IF( IFQS.LE.0 ) IFQS = IBIG
      IDFLT = 1
      IFQO = IBIG
      VARB = 'Reference Node Output File Frequency'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IFQO)
      WRITE(IWR,'(2A,I6,A)') VARB(1:IVR),': Every ',IFQO,' Time Step(s)'
      IF( IFQO.LE.0 ) IFQO = IBIG
      IDFLT = 1
      VARB = 'Time Output Units'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTM)
      WRITE(IWR,'(3A)') VARB(1:IVR),': ',UNTM(1:NCH)
      IDFLT = 1
      VARB = 'Length Output Units'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNLN)
      WRITE(IWR,'(3A)') VARB(1:IVR),': ',UNLN(1:NCH)
      IF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
        IDFLT = 1
        VARB = 'Arc Output Units'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNAR)
        WRITE(IWR,'(3A)') VARB(1:IVR),': ',UNAR(1:NCH)
      ENDIF
      IDFLT = 1
      VARB = 'Screen Significant Digits'
      CALL RDINT(ISTART,ICOMMA,CHDUM,ISGNS)
      WRITE(IWR,'(2A,I2)') VARB(1:IVR),': ',ISGNS
      IDFLT = 1
      VARB = 'Output File Significant Digits'
      CALL RDINT(ISTART,ICOMMA,CHDUM,ISGNO)
      WRITE(IWR,'(2A,I2)') VARB(1:IVR),': ',ISGNO
      IDFLT = 1
      VARB = 'Plot File Significant Digits'
      CALL RDINT(ISTART,ICOMMA,CHDUM,ISGNP)
      WRITE(IWR,'(2A,I2)') VARB(1:IVR),': ',ISGNP
!
!---  Read reference node variables  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Reference Node Variables'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NVREF)
      WRITE( IWR,'(/,A,I6)') 'Reference Node Variables: ',NVREF
      NVC = 0
      NPT = 0
      DO 200 NV = 1,NVREF
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Reference Node Variable'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Reference Node Variable: Solute Name'
          CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SOLNM)
          DO 110 NSL = 1,NSOLU
            IF( SOLNM.EQ.SOLUT(NSL) ) GOTO 120
  110     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Reference-Node Solute Name: '//SOLNM
          CALL WRMSGS( INDX )
          NVC = NVC -1
          GOTO 200
  120     CONTINUE
        ENDIF

        IF( INDEX( ADUM(1:),'species' ).NE.0 ) THEN
          IF( ISLC(40).EQ.0 ) THEN
            NVC = NVC -1
            GOTO 200
          ENDIF
          VARB = 'Reference Node Variable: Reactive Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SPNM)
!
!---      Aqueous species  ---
!
          DO 122 M = 1,NSPL
            NSP = M
            IF( SPNM.EQ.SPNML(M) ) GOTO 128
  122     CONTINUE
!
!---      Solid species  ---
!
          DO 124 M = 1,NSPS
            NSP = M+NSPL
            IF( SPNM.EQ.SPNMS(M) ) GOTO 128
  124     CONTINUE
!
!---      Exchanged species  ---
!
          DO 125 M = 1,NSPE
            NSP = M+NSPL+NSPS
            IF( SPNM.EQ.SPNME(M) ) GOTO 128
  125     CONTINUE
!
!---      Conservation- or kinetic-component species  ---
!
          IF( INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
            DO 126 NSL = NSOLU+1,NSOLU+NEQC+NEQK
              IF( SPNM.EQ.SOLUT(NSL) ) GOTO 128
  126       CONTINUE
          ENDIF
!
!---      Unrecognized species name  ---
!
          INDX = 4
          CHMSG = 'Unrecognized Reference-Node Reactive Species Name: '
     &      // SPNM
          CALL WRMSGS( INDX )
          NVC = NVC -1
          GOTO 200
  128     CONTINUE
        ENDIF

!
!---    Multiple plant temperature option, output up to five plant
!       temperatures
!
        IF( INDEX( ADUM(1:),'plant temp' ).NE.0 .AND.
     &    ISLC(31).EQ.2 ) THEN
          VARB = 'Reference Node Variable: Plant Name'
          CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,PLNTNM)
          DO 130 IP = 1,NPLANT
            IF( PLNTNM.EQ.PLANT(IP) ) THEN
              IF( NPT.GE.5 ) THEN
                INDX = 2
                CHMSG = 'Number of Reference-Node Plant Temperatures' //
     &            ' > 5: Plant Name: ' // PLNTNM
                CALL WRMSGS( INDX )
                GOTO 200
              ENDIF
              ITMP_P(IP) = 1
              NPT = NPT+1
              GOTO 140
            ENDIF
  130     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Reference-Node Plant Name: '//PLNTNM
          CALL WRMSGS( INDX )
  140     CONTINUE
        ENDIF
        IF( INDEX(ADUM(1:),'atmospheric temperature').NE.0 ) THEN
          IREF(NV) = 206
        ELSEIF( INDEX(ADUM(1:),'ground').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 .AND. 
     &    INDEX(ADUM(1:),'temp').NE.0 ) THEN
          IREF(NV) = 206
          VARB = 'Reference Node Variable: Ground-Loop-Well Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,IREF_CW(NV))
          IF( IREF_CW(NV).LT.1 .OR. IREF_CW(NV).GT.N_GLW ) THEN
            INDX = 7
            CHMSG = 'Unrecognized Well Number: '
            IMSG = IREF_CW(NV)
            CALL WRMSGS( INDX )
          ENDIF
          CHREF(206) = 'TGLW'
          UNREF(206) = 'c'
          WRITE(ADUM(NCH+1:NCH+11),'(A,I3)') ': well #',IREF_CW(NV)
          NCH = NCH+11
        ELSEIF( INDEX(ADUM(1:),'atmospheric relative').NE.0 ) THEN
          IREF(NV) = 207
        ELSEIF( INDEX(ADUM(1:),'atmospheric solar').NE.0 ) THEN
          IREF(NV) = 208
        ELSEIF( INDEX(ADUM(1:),'atmospheric wind').NE.0 ) THEN
          IREF(NV) = 209
        ELSEIF( INDEX(ADUM(1:),'surface temperature').NE.0 ) THEN
          IREF(NV) = 213
        ELSEIF( INDEX(ADUM(1:),'surface vapor press').NE.0 ) THEN
          IREF(NV) = 214
        ELSEIF( INDEX(ADUM(1:),'actual evaporation rate').NE.0 ) THEN
          IREF(NV) = 215
        ELSEIF( INDEX(ADUM(1:),'potential evaporation rate').NE.0 ) THEN
          IREF(NV) = 216
        ELSEIF( INDEX(ADUM(1:),'actual transpiration').NE.0 ) THEN
          IREF(NV) = 217
        ELSEIF( INDEX(ADUM(1:),'potential transpiration').NE.0 ) THEN
          IREF(NV) = 218
        ELSEIF( INDEX(ADUM(1:),'atmospheric pressure').NE.0 ) THEN
          IREF(NV) = 224
        ELSEIF( INDEX(ADUM(1:),'surface aqueous press').NE.0 ) THEN
          IREF(NV) = 225
        ELSEIF( INDEX(ADUM(1:),'surface gas pressure').NE.0 ) THEN
          IREF(NV) = 226
        ELSEIF( INDEX(ADUM(1:),'surface aqueous saturation').NE.0 ) THEN
          IREF(NV) = 227
        ELSEIF( INDEX(ADUM(1:),'surface latent heat flux').NE.0 ) THEN
          IREF(NV) = 228
        ELSEIF( INDEX(ADUM(1:),'surface sensible heat flux').NE.0 ) THEN
          IREF(NV) = 229
        ELSEIF( INDEX(ADUM(1:),'surface net long-wave rad').NE.0 ) THEN
          IREF(NV) = 230
        ELSEIF( INDEX(ADUM(1:),'surface net short-wave rad').NE.0 ) THEN
          IREF(NV) = 231
        ELSEIF( INDEX(ADUM(1:),'surface ground heat flux').NE.0 ) THEN
          IREF(NV) = 232
        ELSEIF( INDEX(ADUM(1:),'surface water mass bal').NE.0 ) THEN
          IREF(NV) = 233
        ELSEIF( INDEX(ADUM(1:),'plant temperature').NE.0 ) THEN
          IF( ISLC(31).EQ.2 ) THEN
            IREF(NV) = 233+NPT
            CHREF(233+NPT) = 'TP  '
            IC = ICOUNT(IP)
            IF( IC.EQ.1 ) THEN
              WRITE(CHREF(233+NPT)(3:3),'(I1)') IP
            ELSEIF( IC.EQ.2 ) THEN
              WRITE(CHREF(233+NPT)(3:4),'(I2)') IP
            ELSE
              CHREF(233+NPT) = 'TPXX'
            ENDIF
          ELSE
            IREF(NV) = 234
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'rainfall interception mass').NE.0 ) THEN
          IREF(NV) = 239
        ELSEIF( INDEX(ADUM(1:),'bare-surface aero').NE.0 .OR.
     &    INDEX(ADUM(1:),'bare-soil aero').NE.0 ) THEN
          IREF(NV) = 243
        ELSEIF( INDEX(ADUM(1:),'surface volumetric precip').NE.0 ) THEN
          IREF(NV) = 244
        ELSEIF( INDEX(ADUM(1:),'surface mass precip').NE.0 ) THEN
          IF( ISLC(74).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Reference Node Variable for ' // 
     &        'Fracture or Borehole Flow: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
          IREF(NV) = 245
        ELSEIF( INDEX(ADUM(1:),'atmospheric mass precip').NE.0 ) THEN
          IF( ISLC(74).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Reference Node Variable for ' // 
     &        'Fracture or Borehole Flow: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
          IREF(NV) = 246
        ELSEIF( INDEX(ADUM(1:),'x intrinsic perm').NE.0 ) THEN
          IREF(NV) = 247
        ELSEIF( INDEX(ADUM(1:),'y intrinsic perm').NE.0 ) THEN
          IREF(NV) = 248
        ELSEIF( INDEX(ADUM(1:),'z intrinsic perm').NE.0 ) THEN
          IREF(NV) = 249
        ELSEIF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
          IREF(NV) = 1
        ELSEIF( INDEX(ADUM(1:),'gas pressure').NE.0 ) THEN
          IREF(NV) = 2
        ELSEIF( INDEX(ADUM(1:),'net').NE.0 .AND.
     &    INDEX(ADUM(1:),'normal').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          CHREF(3) = 'PNN'
          IREF(NV) = 3
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          IREF(NV) = 4
        ELSEIF( INDEX(ADUM(1:),'phase condition').NE.0 ) THEN
          IREF(NV) = 5
        ELSEIF( INDEX(ADUM(1:),'aqueous gauge pressure').NE.0 ) THEN
          IREF(NV) = 6
        ELSEIF( INDEX(ADUM(1:),'gas gauge pressure').NE.0 ) THEN
          IREF(NV) = 7
        ELSEIF( INDEX(ADUM(1:),'trapped gas sat').NE.0 ) THEN
          IREF(NV) = 105
        ELSEIF( INDEX(ADUM(1:),'minimum').NE.0 .AND.
     &    INDEX(ADUM(1:),'apparent').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IREF(NV) = 127
        ELSEIF( INDEX(ADUM(1:),'apparent aqueous sat').NE.0 ) THEN
          IREF(NV) = 9
        ELSEIF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
          IREF(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'gas saturation').NE.0 ) THEN
          IREF(NV) = 12
        ELSEIF( INDEX(ADUM(1:),'salt saturation').NE.0 ) THEN
          IREF(NV) = 264
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 ) THEN
          CHREF(205) = 'XMLS'
          IREF(NV) = 205
        ELSEIF( INDEX(ADUM(1:),'aqueous moisture cont').NE.0 ) THEN
          IREF(NV) = 15
        ELSEIF( INDEX(ADUM(1:),'diffusive porosity').NE.0 ) THEN
          IREF(NV) = 20
        ELSEIF( INDEX(ADUM(1:),'water gas mass frac').NE.0 ) THEN
          IREF(NV) = 21
        ELSEIF( INDEX(ADUM(1:),'air gas mass frac').NE.0 ) THEN
          IREF(NV) = 22
        ELSEIF( INDEX(ADUM(1:),'water aqueous mass frac').NE.0 ) THEN
          IREF(NV) = 24
        ELSEIF( INDEX(ADUM(1:),'air aqueous mass frac').NE.0 ) THEN
          IREF(NV) = 25
        ELSEIF( INDEX(ADUM(1:),'aqueous hydraulic head').NE.0 ) THEN
          IREF(NV) = 27
        ELSEIF( INDEX(ADUM(1:),'gas hydraulic head').NE.0 ) THEN
          IREF(NV) = 28
        ELSEIF( INDEX(ADUM(1:),'dist').NE.0 .AND.
     &    INDEX(ADUM(1:),'borehole').NE.0 .AND.
     &    (ISLC(74).EQ.2 .OR. ISLC(74).EQ.3) .AND. N_BH.GT.0 ) THEN
          IREF(NV) = 399
          CHREF(399) = 'DABH'
        ELSEIF ( INDEX(ADUM(1:),'rock/soil type').NE.0 ) THEN
          IREF(NV) = 30
        ELSEIF( INDEX(ADUM(1:),'aqueous relative perm').NE.0 ) THEN
          IREF(NV) = 31
        ELSEIF( INDEX(ADUM(1:),'gas relative perm').NE.0 ) THEN
          IREF(NV) = 32
        ELSEIF( INDEX(ADUM(1:),'aqueous density').NE.0 ) THEN
          IREF(NV) = 34
        ELSEIF( INDEX(ADUM(1:),'gas density').NE.0 ) THEN
          IREF(NV) = 35
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'water').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          CHREF(140) = 'TSRW'
          IREF(NV) = 140
          IREF_CW(NV) = -1
        ELSEIF( INDEX(ADUM(1:),'total water mass').NE.0 ) THEN
          IREF(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total air mass').NE.0 ) THEN
          IREF(NV) = 38
        ELSEIF( INDEX(ADUM(1:),'water mass source int').NE.0 ) THEN
          IREF(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'air mass source int').NE.0 ) THEN
          IREF(NV) = 41
        ELSEIF( INDEX(ADUM(1:),'energy source int').NE.0 ) THEN
          IREF(NV) = 43
        ELSEIF( INDEX(ADUM(1:),'x thermal cond').NE.0 ) THEN
          IREF(NV) = 44
        ELSEIF( INDEX(ADUM(1:),'y thermal cond').NE.0 ) THEN
          IREF(NV) = 45
        ELSEIF( INDEX(ADUM(1:),'z thermal cond').NE.0 ) THEN
          IREF(NV) = 46
        ELSEIF ( INDEX(ADUM(1:),'salt volumetric conc').NE.0 ) THEN
          IREF(NV) = 47
        ELSEIF ( INDEX(ADUM(1:),'salt aqueous conc').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt conc').NE.0 ) THEN
          IREF(NV) = 48
        ELSEIF( INDEX(ADUM(1:),'aqueous courant').NE.0 ) THEN
          ICRNT = 1
          IREF(NV) = 49
        ELSEIF ( INDEX(ADUM(1:),'total salt mass').NE.0 ) THEN
          IREF(NV) = 50
        ELSEIF( INDEX(ADUM(1:),'x aqueous vol').NE.0 ) THEN
          IREF(NV) = 51
        ELSEIF( INDEX(ADUM(1:),'y aqueous vol').NE.0 ) THEN
          IREF(NV) = 52
        ELSEIF( INDEX(ADUM(1:),'z aqueous vol').NE.0 ) THEN
          IREF(NV) = 53
        ELSEIF( INDEX(ADUM(1:),'x gas vol').NE.0 ) THEN
          IREF(NV) = 54
        ELSEIF( INDEX(ADUM(1:),'y gas vol').NE.0 ) THEN
          IREF(NV) = 55
        ELSEIF( INDEX(ADUM(1:),'z gas vol').NE.0 ) THEN
          IREF(NV) = 56
        ELSEIF( INDEX(ADUM(1:),'x heat flux').NE.0 ) THEN
          IREF(NV) = 60
        ELSEIF( INDEX(ADUM(1:),'y heat flux').NE.0 ) THEN
          IREF(NV) = 61
        ELSEIF( INDEX(ADUM(1:),'z heat flux').NE.0 ) THEN
          IREF(NV) = 62
        ELSEIF( INDEX(ADUM(1:),'matric potential').NE.0 ) THEN
          IREF(NV) = 63
        ELSEIF ( INDEX(ADUM(1:),'x salt flux').NE.0 .OR.
     &    INDEX(ADUM(1:),'x-salt flux').NE.0 ) THEN
          IREF(NV) = 64
        ELSEIF ( INDEX(ADUM(1:),'y salt flux').NE.0 .OR.
     &    INDEX(ADUM(1:),'y-salt flux').NE.0 ) THEN
          IREF(NV) = 65
        ELSEIF ( INDEX(ADUM(1:),'z salt flux').NE.0 .OR.
     &    INDEX(ADUM(1:),'z-salt flux').NE.0 ) THEN
          IREF(NV) = 66
        ELSEIF ( INDEX(ADUM(1:),'xnc salt flux').NE.0 ) THEN
          IREF(NV) = 67
        ELSEIF ( INDEX(ADUM(1:),'ync salt flux').NE.0 ) THEN
          IREF(NV) = 68
        ELSEIF ( INDEX(ADUM(1:),'znc salt flux').NE.0 ) THEN
          IREF(NV) = 69
        ELSEIF( INDEX(ADUM(1:),'water gas mole').NE.0 ) THEN
          IREF(NV) = 70
        ELSEIF( INDEX(ADUM(1:),'air gas mole').NE.0 ) THEN
          IREF(NV) = 71
        ELSEIF( INDEX(ADUM(1:),'water gas conc').NE.0 ) THEN
          IREF(NV) = 73
        ELSEIF( INDEX(ADUM(1:),'air gas conc').NE.0 ) THEN
          IREF(NV) = 74
        ELSEIF( INDEX(ADUM(1:),'water aqueous conc').NE.0 ) THEN
          IREF(NV) = 76
        ELSEIF( INDEX(ADUM(1:),'air aqueous conc').NE.0 ) THEN
          IREF(NV) = 77
        ELSEIF( INDEX(ADUM(1:),'gas courant').NE.0 ) THEN
          ICRNT = 1
          IREF(NV) = 79
        ELSEIF( INDEX(ADUM(1:),'aqueous matrix').NE.0 ) THEN
          IREF(NV) = 83
        ELSEIF( INDEX(ADUM(1:),'aqueous fracture').NE.0 ) THEN
          IREF(NV) = 84
        ELSEIF( INDEX(ADUM(1:),'gas matrix').NE.0 ) THEN
          IREF(NV) = 85
        ELSEIF( INDEX(ADUM(1:),'gas fracture').NE.0 ) THEN
          IREF(NV) = 86
        ELSEIF( INDEX(ADUM(1:),'xnc aqueous vol').NE.0 ) THEN
          IREF(NV) = 87
        ELSEIF( INDEX(ADUM(1:),'ync aqueous vol').NE.0 ) THEN
          IREF(NV) = 88
        ELSEIF( INDEX(ADUM(1:),'znc aqueous vol').NE.0 ) THEN
          IREF(NV) = 89
        ELSEIF( INDEX(ADUM(1:),'xnc gas vol').NE.0 ) THEN
          IREF(NV) = 90
        ELSEIF( INDEX(ADUM(1:),'ync gas vol').NE.0 ) THEN
          IREF(NV) = 91
        ELSEIF( INDEX(ADUM(1:),'znc gas vol').NE.0 ) THEN
          IREF(NV) = 92
        ELSEIF( INDEX(ADUM(1:),'xnc heat flux').NE.0 ) THEN
          IREF(NV) = 96
        ELSEIF( INDEX(ADUM(1:),'ync heat flux').NE.0 ) THEN
          IREF(NV) = 97
        ELSEIF( INDEX(ADUM(1:),'znc heat flux').NE.0 ) THEN
          IREF(NV) = 98
        ELSEIF( INDEX(ADUM(1:),'node').NE.0 .AND.
     &    INDEX(ADUM(1:),'number').NE.0 ) THEN
          IREF(NV) = 100
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mass').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mass').NE.0 ) THEN
          IREF(NV) = 110
        ELSEIF( INDEX(ADUM(1:),'water vapor part').NE.0 ) THEN
          IREF(NV) = 128
        ELSEIF( INDEX(ADUM(1:),'water relative humid').NE.0 ) THEN
          IREF(NV) = 137
        ELSEIF( INDEX(ADUM(1:),'salt mass source rate').NE.0 ) THEN
          IREF(NV) = 147
        ELSEIF( INDEX(ADUM(1:),'salt mass source int').NE.0 ) THEN
          IREF(NV) = 148
        ELSEIF( INDEX(ADUM(1:),'actual evaporation').NE.0 ) THEN
          IREF(NV) = 215
        ELSEIF( INDEX(ADUM(1:),'potential evaporation').NE.0 ) THEN
          IREF(NV) = 216
        ELSEIF( INDEX(ADUM(1:),'actual transpiration').NE.0 ) THEN
          IREF(NV) = 217
        ELSEIF( INDEX(ADUM(1:),'potential transpiration').NE.0 ) THEN
          IREF(NV) = 218
        ELSEIF( INDEX(ADUM(1:),'water mass source rate').NE.0 ) THEN
          IREF(NV) = 140
          IREF_CW(NV) = 0
        ELSEIF( INDEX(ADUM(1:),'air mass source rate').NE.0 ) THEN
          IREF(NV) = 141
        ELSEIF( INDEX(ADUM(1:),'energy source rate').NE.0 ) THEN
          IREF(NV) = 143
        ELSEIF( INDEX(ADUM(1:),'aqueous viscosity').NE.0 ) THEN
          IREF(NV) = 176
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'trapped').NE.0 .AND.
     &    INDEX(ADUM(1:),'air').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IREF(NV) = 190
          CHREF(190) = 'IMTGA'
        ELSEIF( INDEX(ADUM(1:),'integrated water mass').NE.0 ) THEN
          IREF(NV) = 191
        ELSEIF( INDEX(ADUM(1:),'integrated air mass').NE.0 ) THEN
          IREF(NV) = 192
        ELSEIF( INDEX(ADUM(1:),'integrated aqueous water').NE.0 ) THEN
          IREF(NV) = 194
        ELSEIF( INDEX(ADUM(1:),'integrated aqueous air').NE.0 ) THEN
          IREF(NV) = 195
        ELSEIF( INDEX(ADUM(1:),'integrated gas water').NE.0 ) THEN
          IREF(NV) = 197
        ELSEIF( INDEX(ADUM(1:),'integrated gas air').NE.0 ) THEN
          IREF(NV) = 198
        ELSEIF( INDEX(ADUM(1:),'fracture').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 ) THEN
          IREF(NV) = 245
          CHREF(245) = 'SRCG'
        ELSEIF( INDEX(ADUM(1:),'fracture').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 ) THEN
          IREF(NV) = 246
          CHREF(246) = 'SRCL'
        ELSEIF( INDEX(ADUM(1:),'surface volumetric runoff').NE.0 ) THEN
          IREF(NV) = 274
        ELSEIF( INDEX(ADUM(1:),'mcstress').NE.0 ) THEN
          IREF(NV) = 275
          CHREF(275) = 'MCSTR'
        ELSEIF( INDEX(ADUM(1:),'gas viscosity').NE.0 ) THEN
          IREF(NV) = 289
        ELSEIF( INDEX(ADUM(1:),'x node position').NE.0 ) THEN
          IREF(NV) = 291
        ELSEIF( INDEX(ADUM(1:),'y node position').NE.0 ) THEN
          IREF(NV) = 292
        ELSEIF( INDEX(ADUM(1:),'z node position').NE.0 ) THEN
          IREF(NV) = 293
        ELSEIF( INDEX(ADUM(1:),'aqueous enthalpy').NE.0 ) THEN
          IREF(NV) = 296
        ELSEIF( INDEX(ADUM(1:),'gas enthalpy').NE.0 ) THEN
          IREF(NV) = 297
        ELSEIF( ( INDEX(ADUM(1:),'similarity').NE.0 .OR.
     &    INDEX(ADUM(1:),'similitude').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'variable').NE.0) THEN
          IREF(NV) = 299
        ELSEIF ( INDEX(ADUM(1:),'fracture').NE.0 .AND.
     &    INDEX(ADUM(1:),'mech').NE.0 .AND.
     &    INDEX(ADUM(1:),'aperture').NE.0 ) THEN
          IREF(NV) = 397
          CHREF(397) = 'APM'
        ELSEIF ( INDEX(ADUM(1:),'fracture').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydr').NE.0 .AND.
     &    INDEX(ADUM(1:),'aperture').NE.0 ) THEN
          IREF(NV) = 398
          CHREF(398) = 'APH'
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'xx').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 80
          CHREF(80) = 'SIG-XX'
          IREFGC(NV) = 0
          IF( INDEX(ADUM(1:),'effective').NE.0 ) THEN
            CHREF(80) = 'ESIGXX'
            IREFGC(NV) = -1
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'yy').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 101
          CHREF(101) = 'SIG-YY'
          IREFGC(NV) = 0
          IF( INDEX(ADUM(1:),'effective').NE.0 ) THEN
            CHREF(101) = 'ESIGYY'
            IREFGC(NV) = -1
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'zz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 130
          CHREF(130) = 'SIG-ZZ'
          IREFGC(NV) = 0
          IF( INDEX(ADUM(1:),'effective').NE.0 ) THEN
            CHREF(130) = 'ESIGZZ'
            IREFGC(NV) = -1
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'yz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 214
          CHREF(214) = 'SIG-YZ'
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'xz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 224
          CHREF(224) = 'SIG-XZ'
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'xy').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 226
          CHREF(226) = 'SIG-XY'
        ELSEIF( INDEX(ADUM(1:),'mean').NE.0 .AND.
     &    INDEX(ADUM(1:),'eff').NE.0 .AND.
     &    INDEX(ADUM(1:),'stress').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 368
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'xx').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 369
          CHREF(369) = 'EPS-XX'
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'yy').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 370
          CHREF(370) = 'EPS-YY'
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'zz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 371
          CHREF(371) = 'EPS-ZZ'
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'yz').NE.0 ) THEN
          IREF(NV) = 372
          CHREF(372) = 'EPS-YZ'
          IF( ISLC(50).EQ.0 ) IREF(NV) = 0
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'xz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 373
          CHREF(373) = 'EPS-XZ'
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'xy').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 374
          CHREF(374) = 'EPS-XY'
        ELSEIF( INDEX(ADUM(1:),'displacement').NE.0 .AND.
     &    INDEX(ADUM(1:),'x').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 375
          IREFGC(NV) = 0
          CHREF(375) = 'DISPX'
          IF( INDEX(ADUM(1:),'fe-node').NE.0 .OR.
     &      INDEX(ADUM(1:),'vertice').NE.0 ) THEN
            VARB = 'Finite Element Node Number'
            CALL RDINT(ISTART,ICOMMA,CHDUM,IREFGC(NV))
            IF( IREFGC(NV).LT.1 .OR. IREFGC(NV).GT.8 ) THEN
              INDX = 7
              CHMSG = 'Unrecognized Finite Element Node Number: '
              IMSG = IREFGC(NV)
              CALL WRMSGS( INDX )
            ENDIF
            WRITE(CHREF(375)(6:6),'(I1)') IREFGC(NV)
            WRITE(ADUM(NCH+1:NCH+3),'(A,I1)') ': ',IREFGC(NV)
            NCH = NCH+3
            IREFGC(NV) = -IREFGC(NV)
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'displacement').NE.0 .AND.
     &    INDEX(ADUM(1:),'y').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 376
          IREFGC(NV) = 0
          CHREF(376) = 'DISPY'
          IF( INDEX(ADUM(1:),'fe-node').NE.0 .OR.
     &      INDEX(ADUM(1:),'vertice').NE.0 ) THEN
            VARB = 'Finite Element Node Number'
            CALL RDINT(ISTART,ICOMMA,CHDUM,IREFGC(NV))
            IF( IREFGC(NV).LT.1 .OR. IREFGC(NV).GT.8 ) THEN
              INDX = 7
              CHMSG = 'Unrecognized Finite Element Node Number: '
              IMSG = IREFGC(NV)
              CALL WRMSGS( INDX )
            ENDIF
            WRITE(CHREF(376)(6:6),'(I1)') IREFGC(NV)
            WRITE(ADUM(NCH+1:NCH+3),'(A,I1)') ': ',IREFGC(NV)
            NCH = NCH+3
            IREFGC(NV) = -IREFGC(NV)
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'displacement').NE.0 .AND.
     &    INDEX(ADUM(1:),'z').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 377
          IREFGC(NV) = 0
          CHREF(377) = 'DISPZ'
          IF( INDEX(ADUM(1:),'fe-node').NE.0 .OR.
     &      INDEX(ADUM(1:),'vertice').NE.0 ) THEN
            VARB = 'Finite Element Node Number'
            CALL RDINT(ISTART,ICOMMA,CHDUM,IREFGC(NV))
            IF( IREFGC(NV).LT.1 .OR. IREFGC(NV).GT.8 ) THEN
              INDX = 7
              CHMSG = 'Unrecognized Finite Element Node Number: '
              IMSG = IREFGC(NV)
              CALL WRMSGS( INDX )
            ENDIF
            WRITE(CHREF(377)(6:6),'(I1)') IREFGC(NV)
            WRITE(ADUM(NCH+1:NCH+3),'(A,I1)') ': ',IREFGC(NV)
            NCH = NCH+3
            IREFGC(NV) = -IREFGC(NV)
          ENDIF

        ELSEIF( (INDEX(ADUM(1:),'solute volumetric conc').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species volumetric conc').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 1
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 1
            CHREF(INDX) = 'SP'
            UNREF(INDX) = 'mol/m^3'
          ENDIF





        ELSEIF( (INDEX(ADUM(1:),'solute aqueous conc').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species aqueous conc').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 2
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 2
            CHREF(INDX) = 'SPL'
            UNREF(INDX) = 'mol/m^3'
          ENDIF




        ELSEIF( INDEX(ADUM(1:),'solute gas conc').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 3
        ELSEIF( INDEX(ADUM(1:),'solute aqueous mol').NE.0 )THEN
          IREF(NV) = 400 + (NSL-1)*33 + 5
        ELSEIF( INDEX(ADUM(1:),'solute gas mol').NE.0 )THEN
          IREF(NV) = 400 + (NSL-1)*33 + 6
        ELSEIF( INDEX(ADUM(1:),'x solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 8
        ELSEIF( INDEX(ADUM(1:),'y solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 9
        ELSEIF( INDEX(ADUM(1:),'z solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 10

        ELSEIF( (INDEX(ADUM(1:),'solute source').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species source').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 11
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 11
            CHREF(INDX) = 'SPSR'
            UNREF(INDX) = 'mol/s'
          ENDIF





        ELSEIF( (INDEX(ADUM(1:),'solute integrated mass').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species integrated mass').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 23
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 23
            CHREF(INDX) = 'SPIM'
            UNREF(INDX) = 'mol'
          ENDIF





        ELSEIF( INDEX(ADUM(1:),'species volumetric conc').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 1
          IREF(NV) = INDX
          CHREF(INDX) = 'SP'
          UNREF(INDX) = 'mol/m^3'
        ELSEIF( INDEX(ADUM(1:),'species aqueous conc').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 2
          IREF(NV) = INDX
          CHREF(INDX) = 'SPL'
          UNREF(INDX) = 'mol/m^3'
        ELSEIF( INDEX(ADUM(1:),'species source').NE.0 ) THEN
          IREF(NV) = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 11
        ELSEIF( INDEX(ADUM(1:),'species integrated mass').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 23
          IREF(NV) = INDX
          CHREF(INDX) = 'SPIM'
          UNREF(INDX) = 'mol'
        ELSEIF( INDEX(ADUM(1:),'mineral area').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 24
          IREF(NV) = INDX
          CHREF(INDX) = 'SPMA'
          UNREF(INDX) = 'm^2'
        ELSEIF( INDEX(ADUM(1:),'mineral rate').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 25
          IREF(NV) = INDX
          CHREF(INDX) = 'SPMR'
          UNREF(INDX) = 'mol/s'
        ELSEIF( INDEX(ADUM(1:),'volume fraction').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 26
          IREF(NV) = INDX
          CHREF(INDX) = 'SPVF'
        ELSEIF( INDEX(ADUM(1:),'ph').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 27
          IREF(NV) = INDX
          CHREF(INDX) = 'pH'

        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'activity').NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 32
          CHREF(INDX) = 'SIA'
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'activity').NE.0 ) THEN
          IREF(NV) = 400+(NSL-1)*33 + 31
          CHREF(INDX) = 'SA'
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Reference Node Variable: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for duplicate reference node variables  ---
!
        DO 190 NX = 1,NV-1
          IF( IREF(NV).EQ.IREF(NX) ) THEN
            IF( IREF_CW(NV).EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Duplicate Reference Node Variable: '//ADUM
              CALL WRMSGS( INDX )
            ELSEIF( IREF_CW(NV).EQ.IREF_CW(NX) ) THEN
              INDX = 4
              CHMSG = 'Duplicate Reference Node Variable: '//ADUM
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
  190   CONTINUE
!
!---    Reference node variable units  ---
!
        IDFLT = 1
        VARB = 'Reference Node Variable Unit'
        CALL RDCHR(ISTART,ICOMMA,NCU,CHDUM,UNREF(IREF(NV)))
        IF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          WRITE( IWR,'(2X,3A,2X,2A,I2,A)' ) ADUM(1:NCH),', ',
     &      UNREF(IREF(NV))(1:NCU),SOLNM(1:NCS),' Solute(',NSL,')'

        ELSEIF( INDEX( ADUM(1:),'species' ).NE.0 ) THEN
          WRITE( IWR,'(2X,3A,2X,2A,I2,A)' ) ADUM(1:NCH),', ',
     &      UNREF(IREF(NV))(1:NCU),SPNM(1:NCS),' Species(',NSP,')'

        ELSE
          WRITE( IWR,'(2X,3A)' ) ADUM(1:NCH),', ',UNREF(IREF(NV))(1:NCU)
        ENDIF
        CALL RDOUUN( IREF(NV) )
        VAR = 0.D+0
        INDX = 0
        CALL RDUNIT( UNREF(IREF(NV)),VAR,INDX )
  200 CONTINUE
      NVREF = NVREF + NVC
!
!---  Plot file output times  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Plot File Output Times'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NPRTM)
      WRITE(IWR,'(/,A)') ' Plot File Output Times:'
      PRTMX = 0.D+0
      IC = 0
      DO 300 N = 1,NPRTM
        IF( IC.GT.1 ) PRTMX = PRTM(IC-1)
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        ICMX = INDEX( CHDUM(ISTART:), ',' )
        IATX = INDEX( CHDUM(ISTART:), '@' )
!
!---    Sequence of plot file output times  ---
!
        IF( IATX.GT.1 .AND. IATX.LT.ICMX ) THEN
          CHDUM(IATX:IATX) = ','
          VARB = 'Count Integer'
          CALL RDINT(ISTART,ICOMMA,CHDUM,IATX )
          VARB = 'Delta Plot File Output Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DTX )
          VARB = 'Plot File Output Time Units'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,DTX,INDX)
          DO 210 II = 1,IATX
            IC = IC + 1
            IF( IC.GT.LPTM ) THEN
        INDX = 5
              CHMSG = 'Number of Output Times > Parameter LPTM'
              CALL WRMSGS( INDX )
            ENDIF
            IF( IC.EQ.1 ) THEN
              PRTM(IC) = DTX
            ELSE
              PRTM(IC) = PRTM(IC-1) + DTX
            ENDIF
            PRTMX = PRTM(IC)
            INDX = 1
            IUNS = 1
            CALL RDUNIT(UNTS,PRTMX,INDX)
            WRITE(IWR,'(2X,1PE11.4,1X,A)') PRTMX,UNTS(1:NCH)
            TMPR = MIN( TMPR,PRTM(IC) )
  210     CONTINUE
!
!---    Single plot file output time  ---
!
        ELSE
          IC = IC + 1
          IF( IC.GT.LPTM ) THEN
        INDX = 5
            CHMSG = 'Number of Output Times > Parameter LPTM'
        CALL WRMSGS( INDX )
      ENDIF
        VARB = 'Plot File Output Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PRTM(IC))
        VARB = 'Plot File Output Time Units'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,1PE11.4,1X,A)') PRTM(IC),UNTS(1:NCH)
        INDX = 0
        IUNS = 1
          CALL RDUNIT(UNTS,PRTM(IC),INDX)
          TMPR = MIN( TMPR,PRTM(IC) )
        ENDIF
 300  CONTINUE
      NPRTM = IC
      WRITE(IWR,'(2X,A)') 'After the Final Time Step'
!
!---  Read Plot File Variables  ---
!
      WRITE( IWR,'(/,A)') 'Plot File Variables:'
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Plot File Variables'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NVPLOT)
      NVC = 0
      DO 400 NV = 1,NVPLOT
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Plot File Variable'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Plot File Variable: Solute Name'
          CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SOLNM)
          DO 310 NSL = 1,NSOLU
            IF( SOLNM.EQ.SOLUT(NSL) ) GOTO 320
  310     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Plot File Solute Name: '//SOLNM
          CALL WRMSGS( INDX )
          NVC = NVC -1
          GOTO 400
  320     CONTINUE
        ENDIF

        IF( INDEX( ADUM(1:),'species' ).NE.0 ) THEN
          IF( ISLC(40).EQ.0 ) THEN
            NVC = NVC -1
            GOTO 400
          ENDIF
          VARB = 'Plot File Variable: Reactive Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SPNM)
!
!---      Conservation- or kinetic-component species  ---
!
          IF( INDEX( SPNM(1:),'total_' ).NE.0 ) THEN
            DO 330 NSL = NSOLU+1,NSOLU+NEQC+NEQK
              IF( SPNM.EQ.SOLUT(NSL) ) GOTO 350
  330       CONTINUE
          ENDIF
!
!---      Aqueous species  ---
!
          DO 332 M = 1,NSPL
            NSP = M
            IF( SPNM.EQ.SPNML(M) ) GOTO 350
  332     CONTINUE
!
!---      Solid species  ---
!
          DO 334 M = 1,NSPS
            NSP = M+NSPL
            IF( SPNM.EQ.SPNMS(M) ) GOTO 350
  334     CONTINUE
!
!---      Exchanged species  ---
!
          DO 336 M = 1,NSPE
            NSP = M + NSPL + NSPS
            IF( SPNM.EQ.SPNME(M) ) GOTO 350
  336     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Plot File Reactive Species Name: '
     &      // SPNM
          CALL WRMSGS( INDX )
          NVC = NVC -1
          GOTO 400
  350     CONTINUE
        ENDIF

        IF( INDEX(ADUM(1:),'final restart').NE.0 ) THEN
          ISLC(18) = 1
          IPLOT(NV) = 200
        ELSEIF( INDEX(ADUM(1:),'no restart').NE.0 ) THEN
          ISLC(18) = 2
          IPLOT(NV) = 200
        ELSEIF( INDEX(ADUM(1:),'3d grid').NE.0 ) THEN
          ISLC(63) = 1
          IPLOT(NV) = 200
        ELSEIF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
          IPLOT(NV) = 1
        ELSEIF( INDEX(ADUM(1:),'gas pressure').NE.0 ) THEN
          IPLOT(NV) = 2
        ELSEIF( INDEX(ADUM(1:),'net').NE.0 .AND.
     &    INDEX(ADUM(1:),'normal').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 3
        ELSEIF( INDEX(ADUM(1:),'surface vapor press').NE.0 ) THEN
          IPLOT(NV) = 214
        ELSEIF( INDEX(ADUM(1:),'surface temperature').NE.0 ) THEN
          IPLOT(NV) = 213
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          IPLOT(NV) = 4
        ELSEIF( INDEX(ADUM(1:),'phase condition').NE.0 ) THEN
          IPLOT(NV) = 5
        ELSEIF( INDEX(ADUM(1:),'aqueous gauge pressure').NE.0 ) THEN
          IPLOT(NV) = 6
        ELSEIF( INDEX(ADUM(1:),'gas gauge pressure').NE.0 ) THEN
          IPLOT(NV) = 7
        ELSEIF( INDEX(ADUM(1:),'trapped gas sat').NE.0 ) THEN
          IPLOT(NV) = 105
        ELSEIF( INDEX(ADUM(1:),'minimum').NE.0 .AND.
     &    INDEX(ADUM(1:),'apparent').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IPLOT(NV) = 127
        ELSEIF( INDEX(ADUM(1:),'apparent aqueous sat').NE.0 ) THEN
          IPLOT(NV) = 9
        ELSEIF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
          IPLOT(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'gas saturation').NE.0 ) THEN
          IPLOT(NV) = 12
        ELSEIF( INDEX(ADUM(1:),'salt saturation').NE.0 ) THEN
          IPLOT(NV) = 264
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mole frac').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mole frac').NE.0 ) THEN
          IPLOT(NV) = 205
        ELSEIF( INDEX(ADUM(1:),'aqueous moisture cont').NE.0 ) THEN
          IPLOT(NV) = 15
        ELSEIF( INDEX(ADUM(1:),'diffusive porosity').NE.0 ) THEN
          IPLOT(NV) = 20
        ELSEIF( INDEX(ADUM(1:),'water gas mass frac').NE.0 ) THEN
          IPLOT(NV) = 21
        ELSEIF( INDEX(ADUM(1:),'air gas mass frac').NE.0 ) THEN
          IPLOT(NV) = 22
        ELSEIF( INDEX(ADUM(1:),'water aqueous mass frac').NE.0 ) THEN
          IPLOT(NV) = 24
        ELSEIF( INDEX(ADUM(1:),'air aqueous mass frac').NE.0 ) THEN
          IPLOT(NV) = 25
        ELSEIF( INDEX(ADUM(1:),'aqueous hydraulic head').NE.0 ) THEN
          IPLOT(NV) = 27
        ELSEIF( INDEX(ADUM(1:),'gas hydraulic head').NE.0 ) THEN
          IPLOT(NV) = 28
        ELSEIF( INDEX(ADUM(1:),'dist').NE.0 .AND.
     &    INDEX(ADUM(1:),'borehole').NE.0 .AND.
     &    (ISLC(74).EQ.2 .OR. ISLC(74).EQ.3) .AND. N_BH.GT.0 ) THEN
          IPLOT(NV) = 399
        ELSEIF ( INDEX(ADUM(1:),'fracture').NE.0 .AND.
     &    INDEX(ADUM(1:),'mech').NE.0 .AND.
     &    INDEX(ADUM(1:),'aperture').NE.0 ) THEN
          IPLOT(NV) = 397
        ELSEIF ( INDEX(ADUM(1:),'rock/soil type').NE.0 ) THEN
          IPLOT(NV) = 30
        ELSEIF( INDEX(ADUM(1:),'aqueous relative perm').NE.0 ) THEN
          IPLOT(NV) = 31
        ELSEIF( INDEX(ADUM(1:),'gas relative perm').NE.0 ) THEN
          IPLOT(NV) = 32
        ELSEIF( INDEX(ADUM(1:),'aqueous density').NE.0 ) THEN
          IPLOT(NV) = 34
        ELSEIF( INDEX(ADUM(1:),'gas density').NE.0 ) THEN
          IPLOT(NV) = 35
        ELSEIF( INDEX(ADUM(1:),'total water mass').NE.0 ) THEN
          IPLOT(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total air mass').NE.0 ) THEN
          IPLOT(NV) = 38
        ELSEIF( INDEX(ADUM(1:),'water mass source int').NE.0 ) THEN
          IPLOT(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'air source').NE.0 ) THEN
          IPLOT(NV) = 41
        ELSEIF( INDEX(ADUM(1:),'energy source').NE.0 ) THEN
          IPLOT(NV) = 43
        ELSEIF( INDEX(ADUM(1:),'x thermal cond').NE.0 ) THEN
          IPLOT(NV) = 44
        ELSEIF( INDEX(ADUM(1:),'y thermal cond').NE.0 ) THEN
          IPLOT(NV) = 45
        ELSEIF( INDEX(ADUM(1:),'z thermal cond').NE.0 ) THEN
          IPLOT(NV) = 46
        ELSEIF ( INDEX(ADUM(1:),'salt conc').NE.0 ) THEN
          IPLOT(NV) = 47
        ELSEIF ( INDEX(ADUM(1:),'salt aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 48
        ELSEIF( INDEX(ADUM(1:),'aqueous courant').NE.0 ) THEN
          ICRNT = 1
          IPLOT(NV) = 49
        ELSEIF ( INDEX(ADUM(1:),'total salt').NE.0 ) THEN
          IPLOT(NV) = 50
        ELSEIF( INDEX(ADUM(1:),'x aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 51
        ELSEIF( INDEX(ADUM(1:),'y aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 52
        ELSEIF( INDEX(ADUM(1:),'z aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 53
        ELSEIF( INDEX(ADUM(1:),'x gas vol').NE.0 ) THEN
          IPLOT(NV) = 54
        ELSEIF( INDEX(ADUM(1:),'y gas vol').NE.0 ) THEN
          IPLOT(NV) = 55
        ELSEIF( INDEX(ADUM(1:),'z gas vol').NE.0 ) THEN
          IPLOT(NV) = 56
        ELSEIF( INDEX(ADUM(1:),'x heat flux').NE.0 ) THEN
          IPLOT(NV) = 60
        ELSEIF( INDEX(ADUM(1:),'y heat flux').NE.0 ) THEN
          IPLOT(NV) = 61
        ELSEIF( INDEX(ADUM(1:),'z heat flux').NE.0 ) THEN
          IPLOT(NV) = 62
        ELSEIF( INDEX(ADUM(1:),'matric potential').NE.0 ) THEN
          IPLOT(NV) = 63
        ELSEIF ( INDEX(ADUM(1:),'x salt flux').NE.0 .OR.
     &    INDEX(ADUM(1:),'x-salt flux').NE.0 ) THEN
          IPLOT(NV) = 64
        ELSEIF ( INDEX(ADUM(1:),'y salt flux').NE.0 .OR.
     &    INDEX(ADUM(1:),'y-salt flux').NE.0 ) THEN
          IPLOT(NV) = 65
        ELSEIF ( INDEX(ADUM(1:),'z salt flux').NE.0 .OR.
     &    INDEX(ADUM(1:),'z-salt flux').NE.0 ) THEN
          IPLOT(NV) = 66
        ELSEIF ( INDEX(ADUM(1:),'xnc salt flux').NE.0 ) THEN
          IPLOT(NV) = 67
        ELSEIF ( INDEX(ADUM(1:),'ync salt flux').NE.0 ) THEN
          IPLOT(NV) = 68
        ELSEIF ( INDEX(ADUM(1:),'znc salt flux').NE.0 ) THEN
          IPLOT(NV) = 69
        ELSEIF( INDEX(ADUM(1:),'water gas mole').NE.0 ) THEN
          IPLOT(NV) = 70
        ELSEIF( INDEX(ADUM(1:),'air gas mole').NE.0 ) THEN
          IPLOT(NV) = 71
        ELSEIF( INDEX(ADUM(1:),'water gas conc').NE.0 ) THEN
          IPLOT(NV) = 73
        ELSEIF( INDEX(ADUM(1:),'air gas conc').NE.0 ) THEN
          IPLOT(NV) = 74
        ELSEIF( INDEX(ADUM(1:),'water aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 76
        ELSEIF( INDEX(ADUM(1:),'air aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 77
        ELSEIF( INDEX(ADUM(1:),'gas courant').NE.0 ) THEN
          ICRNT = 1
          IPLOT(NV) = 79
        ELSEIF( INDEX(ADUM(1:),'aqueous matrix').NE.0 ) THEN
          IPLOT(NV) = 83
        ELSEIF( INDEX(ADUM(1:),'aqueous fracture').NE.0 ) THEN
          IPLOT(NV) = 84
        ELSEIF( INDEX(ADUM(1:),'gas matrix').NE.0 ) THEN
          IPLOT(NV) = 85
        ELSEIF( INDEX(ADUM(1:),'gas fracture').NE.0 ) THEN
          IPLOT(NV) = 86
        ELSEIF( INDEX(ADUM(1:),'xnc aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 87
        ELSEIF( INDEX(ADUM(1:),'ync aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 88
        ELSEIF( INDEX(ADUM(1:),'znc aqueous vol').NE.0 ) THEN
          IPLOT(NV) = 89
        ELSEIF( INDEX(ADUM(1:),'xnc gas vol').NE.0 ) THEN
          IPLOT(NV) = 90
        ELSEIF( INDEX(ADUM(1:),'ync gas vol').NE.0 ) THEN
          IPLOT(NV) = 91
        ELSEIF( INDEX(ADUM(1:),'znc gas vol').NE.0 ) THEN
          IPLOT(NV) = 92
        ELSEIF( INDEX(ADUM(1:),'xnc heat flux').NE.0 ) THEN
          IPLOT(NV) = 96
        ELSEIF( INDEX(ADUM(1:),'ync heat flux').NE.0 ) THEN
          IPLOT(NV) = 97
        ELSEIF( INDEX(ADUM(1:),'znc heat flux').NE.0 ) THEN
          IPLOT(NV) = 98
        ELSEIF( INDEX(ADUM(1:),'node').NE.0 .AND.
     &    INDEX(ADUM(1:),'number').NE.0 ) THEN
          IPLOT(NV) = 100
        ELSEIF( INDEX(ADUM(1:),'salt aqueous mass').NE.0 .OR.
     &    INDEX(ADUM(1:),'aqueous salt mass').NE.0 ) THEN
          IPLOT(NV) = 110
        ELSEIF( INDEX(ADUM(1:),'water vapor part').NE.0 ) THEN
          IPLOT(NV) = 128
        ELSEIF( INDEX(ADUM(1:),'water relative humid').NE.0 ) THEN
          IPLOT(NV) = 137
        ELSEIF( INDEX(ADUM(1:),'water mass source rate').NE.0 ) THEN
          IPLOT(NV) = 140
        ELSEIF ( INDEX(ADUM(1:),'fracture').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydr').NE.0 .AND.
     &    INDEX(ADUM(1:),'aperture').NE.0 ) THEN
          IPLOT(NV) = 398
        ELSEIF( INDEX(ADUM(1:),'salt mass source rate').NE.0 ) THEN
          IPLOT(NV) = 147
        ELSEIF( INDEX(ADUM(1:),'salt mass source int').NE.0 ) THEN
          IPLOT(NV) = 148
        ELSEIF( INDEX(ADUM(1:),'aqueous viscosity').NE.0 ) THEN
          IPLOT(NV) = 176
        ELSEIF( INDEX(ADUM(1:),'x intrinsic perm').NE.0 ) THEN
          IPLOT(NV) = 247
        ELSEIF( INDEX(ADUM(1:),'y intrinsic perm').NE.0 ) THEN
          IPLOT(NV) = 248
        ELSEIF( INDEX(ADUM(1:),'z intrinsic perm').NE.0 ) THEN
          IPLOT(NV) = 249
        ELSEIF( INDEX(ADUM(1:),'mcstress').NE.0 ) THEN
          IPLOT(NV) = 275
        ELSEIF( INDEX(ADUM(1:),'gas viscosity').NE.0 ) THEN
          IPLOT(NV) = 289
        ELSEIF( INDEX(ADUM(1:),'x node position').NE.0 ) THEN
          IPLOT(NV) = 291
        ELSEIF( INDEX(ADUM(1:),'y node position').NE.0 ) THEN
          IPLOT(NV) = 292
        ELSEIF( INDEX(ADUM(1:),'z node position').NE.0 ) THEN
          IPLOT(NV) = 293
        ELSEIF( INDEX(ADUM(1:),'aqueous enthalpy').NE.0 ) THEN
          IPLOT(NV) = 296
        ELSEIF( INDEX(ADUM(1:),'gas enthalpy').NE.0 ) THEN
          IPLOT(NV) = 297
        ELSEIF( ( INDEX(ADUM(1:),'similarity').NE.0 .OR.
     &    INDEX(ADUM(1:),'similitude').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'variable').NE.0) THEN
          IPLOT(NV) = 299
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'xx').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 80
          IPLOTGC(NV) = 0
          IF( INDEX(ADUM(1:),'effective').NE.0 ) THEN
            IPLOTGC(NV) = 1
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'yy').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 101
          IPLOTGC(NV) = 0
          IF( INDEX(ADUM(1:),'effective').NE.0 ) THEN
            IPLOTGC(NV) = 1
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'zz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 130
          IPLOTGC(NV) = 0
          IF( INDEX(ADUM(1:),'effective').NE.0 ) THEN
            IPLOTGC(NV) = 1
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'yz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 214
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'xz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 224
        ELSEIF( INDEX(ADUM(1:),'stress').NE.0 .AND.
     &    INDEX(ADUM(1:),'xy').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 226
        ELSEIF( INDEX(ADUM(1:),'mean').NE.0 .AND.
     &    INDEX(ADUM(1:),'eff').NE.0 .AND.
     &    INDEX(ADUM(1:),'stress').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 368
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'xx').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 369
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'yy').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 370
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'zz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 371
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'yz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 372
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'xz').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 373
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'xy').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 374
        ELSEIF( INDEX(ADUM(1:),'displacement').NE.0 .AND.
     &    INDEX(ADUM(1:),'x').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 375
          IPLOTGC(NV) = 0
          IF( INDEX(ADUM(1:),'fe-node').NE.0 .OR.
     &      INDEX(ADUM(1:),'vertice').NE.0 ) THEN
            VARB = 'Finite Element Node Number'
            CALL RDINT(ISTART,ICOMMA,CHDUM,IPLOTGC(NV))
            IF( IPLOTGC(NV).LT.1 .OR. IPLOTGC(NV).GT.8 ) THEN
              INDX = 7
              CHMSG = 'Unrecognized Finite Element Node Number: '
              IMSG = IPLOTGC(NV)
              CALL WRMSGS( INDX )
            ENDIF
            WRITE(ADUM(NCH+1:NCH+3),'(A,I1)') ': ',IPLOTGC(NV)
            NCH = NCH+3
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'displacement').NE.0 .AND.
     &    INDEX(ADUM(1:),'y').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 376
          IPLOTGC(NV) = 0
          IF( INDEX(ADUM(1:),'fe-node').NE.0 .OR.
     &      INDEX(ADUM(1:),'vertice').NE.0 ) THEN
            VARB = 'Finite Element Node Number'
            CALL RDINT(ISTART,ICOMMA,CHDUM,IPLOTGC(NV))
            IF( IPLOTGC(NV).LT.1 .OR. IPLOTGC(NV).GT.8 ) THEN
              INDX = 7
              CHMSG = 'Unrecognized Finite Element Node Number: '
              IMSG = IPLOTGC(NV)
              CALL WRMSGS( INDX )
            ENDIF
            WRITE(ADUM(NCH+1:NCH+3),'(A,I1)') ': ',IPLOTGC(NV)
            NCH = NCH+3
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'displacement').NE.0 .AND.
     &    INDEX(ADUM(1:),'z').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 377
          IPLOTGC(NV) = 0
          IF( INDEX(ADUM(1:),'fe-node').NE.0 .OR.
     &      INDEX(ADUM(1:),'vertice').NE.0 ) THEN
            VARB = 'Finite Element Node Number'
            CALL RDINT(ISTART,ICOMMA,CHDUM,IPLOTGC(NV))
            IF( IPLOTGC(NV).LT.1 .OR. IPLOTGC(NV).GT.8 ) THEN
              INDX = 7
              CHMSG = 'Unrecognized Finite Element Node Number: '
              IMSG = IPLOTGC(NV)
              CALL WRMSGS( INDX )
            ENDIF
            WRITE(ADUM(NCH+1:NCH+3),'(A,I1)') ': ',IPLOTGC(NV)
            NCH = NCH+3
          ENDIF

        ELSEIF( (INDEX(ADUM(1:),'solute volumetric conc').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species volumetric conc').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 1
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 1
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF





        ELSEIF( (INDEX(ADUM(1:),'solute aqueous conc').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species aqueous conc').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 2
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 2
            UNPLOT(INDX) = 'mol/m^3'
          ENDIF




        ELSEIF( INDEX(ADUM(1:),'solute gas conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 3
        ELSEIF( INDEX(ADUM(1:),'solute aqueous mol').NE.0 )THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 5
        ELSEIF( INDEX(ADUM(1:),'solute gas mol').NE.0 )THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 6
        ELSEIF( INDEX(ADUM(1:),'x solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 8
        ELSEIF( INDEX(ADUM(1:),'y solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 9
        ELSEIF( INDEX(ADUM(1:),'z solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 10

        ELSEIF( (INDEX(ADUM(1:),'solute source').NE.0) .OR.
     &    ((INDEX(ADUM(1:),'species source').NE.0) .AND.
     &    (INDEX( SPNM(1:),'total_' ).NE.0)) ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 11
          IF( NSL.GT.NSOLU ) THEN
            INDX = 400+(NSL-1)*33 + 11
            UNPLOT(INDX) = 'mol/s'
          ENDIF





        ELSEIF( INDEX(ADUM(1:),'species volumetric conc').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 1
          IPLOT(NV) = INDX
          UNPLOT(INDX) = 'mol/m^3'
        ELSEIF( INDEX(ADUM(1:),'species aqueous conc').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 2
          IPLOT(NV) = INDX
          UNPLOT(INDX) = 'mol/m^3'
        ELSEIF( INDEX(ADUM(1:),'species source').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 11
          IPLOT(NV) = INDX
          UNPLOT(INDX) = 'mol'
        ELSEIF( INDEX(ADUM(1:),'mineral area').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 24
          IPLOT(NV) = INDX
          UNPLOT(INDX) = 'm^2'
        ELSEIF( INDEX(ADUM(1:),'mineral rate').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 25
          IPLOT(NV) = INDX
          UNPLOT(INDX) = 'mol/s'
        ELSEIF( INDEX(ADUM(1:),'volume fraction').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 26
          IPLOT(NV) = INDX
        ELSEIF( INDEX(ADUM(1:),'ph').NE.0 ) THEN
          INDX = 400+(NSOLU*33)+((NEQC+NEQK)*33)+(NSP-1)*33 + 27
          IPLOT(NV) = INDX

        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'activity').NE.0 ) THEN
          IPLOT(NV) = 400+(NSL-1)*33 + 31
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Plot File Variable: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for duplicate plot file variables  ---
!
        DO 390 NX = 1,NV-1
          IF( IPLOT(NV).EQ.IPLOT(NX) ) THEN
            INDX = 4
            CHMSG = 'Duplicate Plot File Variable: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
  390   CONTINUE
!
!---    Reference node variable units  ---
!
        IDFLT = 1
        VARB = 'Plot File Variable Units'
        CALL RDCHR(ISTART,ICOMMA,NCU,CHDUM,UNPLOT(IPLOT(NV)))
        IF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          WRITE( IWR,'(2X,3A,2X,2A,I2,A)' ) ADUM(1:NCH),', ',
     &      UNPLOT(IPLOT(NV))(1:NCU),SOLNM(1:NCS),' Solute(',NSL,')'
        ELSE
          WRITE( IWR,'(2X,3A)' ) ADUM(1:NCH),', ',
     &      UNPLOT(IPLOT(NV))(1:NCU)
        ENDIF
        CALL RDOUUN( IPLOT(NV) )
        VAR = 0.D+0
        INDX = 0
        CALL RDUNIT( UNPLOT(IPLOT(NV)),VAR,INDX )
  400 CONTINUE
      NVPLOT = NVPLOT + NVC
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDOU_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSCLF_GT
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
!     Water Mode
!
!     Read input file for rock/soil scaling factor information.
!
!     1 -- Saturated Hydraulic Conductivity
!     2 -- Diffusive Porosity
!     3 -- van Genuchten "alpha" parameter
!     3 -- Brooks and Corey "psi" parameter
!     4 -- van Genuchten "n" parameter
!     4 -- Brooks and Corey "lambda" parameter
!     5 -- Residual saturation
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 24 April 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSCLF_GT'
!
!---  Write card information to output file  ---
!
      CARD = 'Scaling Factor Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read scaling factor equation types  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Saturated Hydraulic Conductivity Scaling Function'
      CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
      IF( INDEX(ADUM(1:),'log').NE.0 ) THEN
        IGAMMA(1) = 2
        WRITE(IWR,'(2A)') VARB(1:NCH),': Logarithmic Scaling'
      ELSE
        IGAMMA(1) = 1
        WRITE(IWR,'(2A)') VARB(1:NCH),': Linear Scaling'
      ENDIF
      VARB = 'Diffusive Porosity Scaling Function'
      CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
      IF( INDEX(ADUM(1:),'log').NE.0 ) THEN
        IGAMMA(2) = 2
        WRITE(IWR,'(2A)') VARB(1:NCH),': Logarithmic Scaling'
      ELSE
        IGAMMA(2) = 1
        WRITE(IWR,'(2A)') VARB(1:NCH),': Linear Scaling'
      ENDIF
      VARB = 'van Genuchten "alpha" or Brooks/Corey "psi" ' //
     &  'Scaling Function'
      CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
      IF( INDEX(ADUM(1:),'log').NE.0 ) THEN
        IGAMMA(3) = 2
        WRITE(IWR,'(2A)') VARB(1:NCH),': Logarithmic Scaling'
      ELSE
        IGAMMA(3) = 1
        WRITE(IWR,'(2A)') VARB(1:NCH),': Linear Scaling'
      ENDIF
      VARB = 'van Genuchten "n" or Brooks/Corey "lambda" ' //
     &  'Scaling Function'
      CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
      IF( INDEX(ADUM(1:),'log').NE.0 ) THEN
        IGAMMA(4) = 2
        WRITE(IWR,'(2A)') VARB(1:NCH),': Logarithmic Scaling'
      ELSE
        IGAMMA(4) = 1
        WRITE(IWR,'(2A)') VARB(1:NCH),': Linear Scaling'
      ENDIF
      VARB = 'Residual Saturation Scaling Function'
      CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
      IF( INDEX(ADUM(1:),'log').NE.0 ) THEN
        IGAMMA(5) = 2
        WRITE(IWR,'(2A)') VARB(1:NCH),': Logarithmic Scaling'
      ELSE
        IGAMMA(5) = 1
        WRITE(IWR,'(2A)') VARB(1:NCH),': Linear Scaling'
      ENDIF
!
!---  Read input lines until all rock/soil types are found  ---
!
      N = 0
      IJK = 0
   10 CONTINUE
      IF( N.GE.NROCK .OR. IJK.GT.0 ) GOTO 500
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Rock/Soil Name'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---  Search known rock types for a matching type ---
!
      DO 100 M = 1, NROCK
        IF( ADUM .EQ. ROCK(M)) THEN
           IROCK = M
           GOTO 200
        ENDIF
  100 CONTINUE
      INDX = 2
      CHMSG = 'Unrecognized Rock/Soil Type: '//ADUM(1:NCH)
      CALL WRMSGS( INDX )
      GOTO 10
  200 CONTINUE
      WRITE (IWR,'(/,2A)') 'Rock/Soil Name: ',ROCK(IROCK)
      N = N+1
!
!---  Read saturated hydraulic conductivity (intrinsic permeability)
!     scaling factor ---
!
      IDFLT = 1
      VARB = 'Saturated Hydraulic Conductivity Scaling Factor'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(1,IROCK))
!
!---  Read fracture saturated hydraulic conductivity
!     (intrinsic permeability) scaling factor ---
!
      IF( INDEX(ADUM(1:),'fractured').NE.0 .OR.
     &  INDEX(ADUM(1:),'dp').NE.0 .OR.
     &  INDEX(ADUM(1:),'dual').NE.0 ) THEN
        IDFLT = 1
        VARB = 'Fracture Saturated Hydraulic Conductivity' //
     &    ' Scaling Factor'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(6,IROCK))
      ENDIF
!
!---  Read diffusive porosity scaling factor  ---
!
      IDFLT = 1
      VARB = 'Diffusive Porosity Scaling Factor'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(2,IROCK))
!
!---  Read fracture diffusive porosity scaling factor  ---
!
      IF( INDEX(ADUM(1:),'fractured').NE.0 .OR.
     &  INDEX(ADUM(1:),'dp').NE.0 .OR.
     &  INDEX(ADUM(1:),'dual').NE.0 ) THEN
        IDFLT = 1
        VARB = 'Fracture Diffusive Porosity' //
     &    ' Scaling Factor'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(7,IROCK))
      ENDIF
!
!---  Read van Genuchten "alpha" or Brooks/Corey "psi"
!     scaling factor ---
!
      IDFLT = 1
      VARB = 'van Genuchten "alpha" Brooks/Corey "psi" Scaling Factor'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(3,IROCK))
!
!---  Read fracture van Genuchten "alpha" or Brooks/Corey "psi"
!     scaling factor ---
!
      IF( INDEX(ADUM(1:),'fractured').NE.0 .OR.
     &  INDEX(ADUM(1:),'dp').NE.0 .OR.
     &  INDEX(ADUM(1:),'dual').NE.0 ) THEN
        IDFLT = 1
        VARB = 'Fracture van Genuchten "alpha" Brooks/Corey "psi"' //
     &    ' Scaling Factor'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(8,IROCK))
      ENDIF
!
!---  Read van Genuchten "n" or Brooks/Corey "lambda"
!     scaling factor ---
!
      IDFLT = 1
      VARB = 'van Genuchten "n" Brooks/Corey "lambda" Scaling Factor'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(4,IROCK))
!
!---  Read fracture van Genuchten "n" or Brooks/Corey "lambda"
!     scaling factor ---
!
      IF( INDEX(ADUM(1:),'fractured').NE.0 .OR.
     &  INDEX(ADUM(1:),'dp').NE.0 .OR.
     &  INDEX(ADUM(1:),'dual').NE.0 ) THEN
        IDFLT = 1
        VARB = 'Fracture van Genuchten "n" Brooks/Corey "lambda"' //
     &    ' Scaling Factor'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(9,IROCK))
      ENDIF
!
!---  Read residual saturation scaling factor ---
!
      IDFLT = 1
      VARB = 'Residual Saturation Scaling Factor'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(5,IROCK))
!
!---  Read fracture residual saturation scaling factor ---
!
      IF( INDEX(ADUM(1:),'fractured').NE.0 .OR.
     &  INDEX(ADUM(1:),'dp').NE.0 .OR.
     &  INDEX(ADUM(1:),'dual').NE.0 ) THEN
        IDFLT = 1
        VARB = 'Fracture Residual Saturation' //
     &    ' Scaling Factor'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(10,IROCK))
      ENDIF
!
!---  Read scaling parameters for Mualem-Anisotropy
!     Relative Permeability ---
!
      IF( IRPL(IROCK).EQ.301 ) THEN
        IDFLT = 1
        VARB = 'van Genuchten "m" or Brooks/Corey "lambda"' //
     &    ' Scaling Factor'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(11,IROCK))
        IDFLT = 1
        VARB = 'Horizontal Pore-Scale Scaling Factor'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(12,IROCK))
        IDFLT = 1
        VARB = 'Vertical Pore-Scale Scaling Factor'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,GAMMA(13,IROCK))
      ENDIF
      IF( N .LT. NROCK ) WRITE(IWR,'(/)')
      GOTO 10
 500  CONTINUE
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSCLF_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSF_GT
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
!     Read input file surface flux information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PARM_BH
      USE OUTPU
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,BDUM,FDUM,GDUM
      CHARACTER*512 CHDUM,CHDUMX
      LOGICAL FLG_CHK
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSF_GT'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Surface Flux Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read surface flux card information  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Surface Flux Inputs'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NSF)
      IF( NSF.GT.LSF ) THEN
        INDX = 5
        CHMSG = 'Number of Surface Flux Domains > Parameter LSF'
        CALL WRMSGS( INDX )
      ENDIF
      NSFF = 0
      DO 100 NS = 1, NSF
        IF( NS.NE.1 ) WRITE(IWR, '(/)')
        CALL RDINPL( CHDUM )
        CHDUMX = CHDUM
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---  Check for specified surface flux filename  ---
!
        CALL CHKINT(ISTART,ICOMMA,CHDUM,INDX)
        IF( INDX .EQ. 1 ) THEN
          VARB = 'Number of Surface Flux Inputs for the Specified File'
          CALL RDINT(ISTART,ICOMMA,CHDUMX,NSFF)
          IF( NSFF.LT.1 ) THEN
            INDX = 4
            CHMSG = 'Number of Surface Flux Inputs < 1'
            CALL WRMSGS( INDX )
          ENDIF
          VARB = 'Surface Output Filename: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUMX,ADUM)
          NSFGP = NSFGP + 1
          ISFGP(NSFGP) = NSFF
          IF( NSFGP.GT.20 ) THEN
            INDX = 4
            CHMSG = 'Number of Surface Flux Files > 20'
            CALL WRMSGS( INDX )
          ENDIF
          FNSF(NSFGP) = ADUM
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
        ENDIF
        IF( NSFF.GT.0 ) THEN
          ISFF(NS) = NSFGP
          NSFF = NSFF-1
        ELSE
          NSFF = 0
          ISFF(NS) = 1
          ISFGP(1) = ISFGP(1) + 1
        ENDIF
!
!---  Read surface flux type  ---
!
        VARB = 'Surface Flux Type'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,A,$)') VARB(1:IVR),': '
        ISFT_FMX = 0
        ISFT_BHX = 0
        IF( INDEX(ADUM(1:),'aqueous').NE.0) THEN
          IF( INDEX(ADUM(1:),'volum').NE.0) THEN
            IF( INDEX(ADUM(1:),'fracture-matrix').NE.0) THEN
              ISFT_FMX = 1
              ISFT(NS) = 66
              WRITE(IWR,'(A)') 'Aqueous Volumetric Fracture-Matrix ' // 
     &          'Flux Surface'
              UNSF(1,NS) = 'm^3/s'
              UNSF(2,NS) = 'm^3'
            ELSEIF( INDEX(ADUM(1:),'borehole').NE.0) THEN
              ISFT_BHX = 1
              ISFT(NS) = 67
              WRITE(IWR,'(A)') 'Aqueous Volumetric Borehole ' // 
     &          'Flux Surface'
              UNSF(1,NS) = 'm^3/s'
              UNSF(2,NS) = 'm^3'
            ELSE
              ISFT(NS) = 2
              WRITE(IWR,'(A)') 'Aqueous Volumetric Flux Surface'
              UNSF(1,NS) = 'm^3/s'
              UNSF(2,NS) = 'm^3'
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'water').NE.0 .OR.
     &      INDEX(ADUM(1:),'h2o').NE.0 ) THEN
            ISFT(NS) = 43
            WRITE(IWR,'(A)') 'Aqueous Water Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 5
            WRITE(IWR,'(A)') 'Aqueous-Phase Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'air').NE.0) THEN
          ISFT(NS) = 30
          WRITE(IWR,'(A)') 'Air-Mass Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'volum').NE.0) THEN
            ISFT(NS) = 3
            WRITE(IWR,'(A)') 'Gas-Phase Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ELSEIF( INDEX(ADUM(1:),'advective').NE.0) THEN
            IF( INDEX(ADUM(1:),'heat').NE.0) THEN
              ISFT(NS) = 20
              WRITE(IWR,'(A)') 'Gas-Advective Heat Flux Surface'
              UNSF(1,NS) = 'w'
              UNSF(2,NS) = 'j'
            ELSEIF( INDEX(ADUM(1:),'water').NE.0) THEN
              ISFT(NS) = 21
              WRITE(IWR,'(A)') 'Gas-Advective Water-Mass Flux Surface'
              UNSF(1,NS) = 'kg/s'
              UNSF(2,NS) = 'kg'
            ELSEIF( INDEX(ADUM(1:),'air').NE.0) THEN
              ISFT(NS) = 22
              WRITE(IWR,'(A)') 'Gas-Advective Air-Mass Flux Surface'
              UNSF(1,NS) = 'kg/s'
              UNSF(2,NS) = 'kg'
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'diffusive').NE.0) THEN
            IF( INDEX(ADUM(1:),'heat').NE.0) THEN
              ISFT(NS) = 25
              WRITE(IWR,'(A)') 'Gas-Diffusive Heat Flux Surface'
              UNSF(1,NS) = 'w'
              UNSF(2,NS) = 'j'
            ELSEIF( INDEX(ADUM(1:),'water').NE.0) THEN
              ISFT(NS) = 26
              WRITE(IWR,'(A)') 'Gas-Diffusive Water-Mass Flux Surface'
              UNSF(1,NS) = 'kg/s'
              UNSF(2,NS) = 'kg'
            ELSEIF( INDEX(ADUM(1:),'air').NE.0) THEN
              ISFT(NS) = 27
              WRITE(IWR,'(A)') 'Gas-Diffusive Air-Mass Flux Surface'
              UNSF(1,NS) = 'kg/s'
              UNSF(2,NS) = 'kg'
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'water').NE.0 .OR.
     &      INDEX(ADUM(1:),'h2o').NE.0 ) THEN
            ISFT(NS) = 44
            WRITE(IWR,'(A)') 'Gas Water Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 6
            WRITE(IWR,'(A)') 'Gas-Phase Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'condensate').NE.0 ) THEN
          ISFT(NS) = 10
          WRITE(IWR,'(A)') 'Condensate Water Mass Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'heat').NE.0) THEN
          ISFT(NS) = 1
          WRITE(IWR,'(A)') 'Heat Flux Surface'
          UNSF(1,NS) = 'w'
          UNSF(2,NS) = 'j'
        ELSEIF( INDEX(ADUM(1:),'evaporation').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'actual').NE.0 ) THEN
            ISFT(NS) = 34
            WRITE(IWR,'(A)') 'Actual Evaporation Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'potential').NE.0 ) THEN
            ISFT(NS) = 35
            WRITE(IWR,'(A)') 'Potential Evaporation Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'transpiration').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'actual').NE.0 ) THEN
            ISFT(NS) = 36
            WRITE(IWR,'(A)') 'Actual Transpiration Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'potential').NE.0 ) THEN
            ISFT(NS) = 37
            WRITE(IWR,'(A)') 'Potential Transpiration Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'net').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'total').NE.0 ) THEN
            ISFT(NS) = 38
            WRITE(IWR,'(A)') 'Net Total Radiation Surface'
            UNSF(1,NS) = 'w'
            UNSF(2,NS) = 'j'
          ELSEIF( INDEX(ADUM(1:),'short-wave').NE.0 ) THEN
            ISFT(NS) = 39
            WRITE(IWR,'(A)') 'Net Short-Wave Radiation Surface'
            UNSF(1,NS) = 'w'
            UNSF(2,NS) = 'j'
          ELSEIF( INDEX(ADUM(1:),'long-wave').NE.0 ) THEN
            ISFT(NS) = 40
            WRITE(IWR,'(A)') 'Net Long-Wave Radiation Surface'
            UNSF(1,NS) = 'w'
            UNSF(2,NS) = 'j'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'water').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'balance').NE.0 ) THEN
          ISFT(NS) = 41
          WRITE(IWR,'(A)') 'Water-Mass Balance Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'rain').NE.0 .AND.
     &    INDEX(ADUM(1:),'water').NE.0 .AND.
     &    INDEX(ADUM(1:),'runoff').NE.0 ) THEN
          ISFT(NS) = 42
          WRITE(IWR,'(A)') 'Rain-Water Runoff Surface'
          UNSF(1,NS) = 'm^3/s'
          UNSF(2,NS) = 'm^3'
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          DO 10 NSL = 1,NSOLU
            IDB = INDEX(SOLUT(NSL)(1:),'  ')
            IF( INDEX(BDUM(1:),SOLUT(NSL)(1:IDB)).NE.0 ) THEN
              ISFT(NS) = NSL+100
              WRITE(IWR,'(2X,2A)') SOLUT(NSL),' Flux Surface'
              UNSF(1,NS) = 'sol/s'
              UNSF(2,NS) = 'sol'
              GOTO 12
            ENDIF
   10     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Surface Flux Solute Name: '//BDUM
          CALL WRMSGS( INDX )
   12     CONTINUE

!
!---    Conservation-component species surface flux input  ---
!
        ELSEIF( INDEX(ADUM(1:),'conservation').NE.0 .AND.
     &    INDEX(ADUM(1:),'component').NE.0 ) THEN
          VARB = 'Conservation-Component Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          DO 14 NSL = NSOLU+1,NSOLU+NEQC
            IDB = INDEX(SOLUT(NSL)(1:),'  ')
            IF( BDUM(1:IDB).EQ.SOLUT(NSL)(1:IDB) ) THEN
              ISFT(NS) = NSL+100
              WRITE(IWR,'(2X,2A)') SOLUT(NSL),' Flux Surface'
              UNSF(1,NS) = 'mol/s'
              UNSF(2,NS) = 'mol'
              GOTO 16
            ENDIF
   14     CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Conservation-Component ' // 
     &       'Species Name: '//BDUM
            CALL WRMSGS( INDX )
   16     CONTINUE
!
!---    Kinetic-component species surface flux input  ---
!
        ELSEIF( INDEX(ADUM(1:),'kinetic').NE.0 .AND.
     &    INDEX(ADUM(1:),'component').NE.0 ) THEN
          VARB = 'Kinetic-Component Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          DO 18 NSL = NSOLU+NEQC+1,NSOLU+NEQC+NEQK
            IDB = INDEX(SOLUT(NSL)(1:),'  ')
            IF( BDUM(1:IDB).EQ.SOLUT(NSL)(1:IDB) ) THEN
              ISFT(NS) = NSL+100
              WRITE(IWR,'(2X,2A)') SOLUT(NSL),' Flux Surface'
              UNSF(1,NS) = 'mol/s'
              UNSF(2,NS) = 'mol'
              GOTO 20
            ENDIF
   18     CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Kinetic-Component ' // 
     &       'Species Name: '//BDUM
            CALL WRMSGS( INDX )
   20     CONTINUE

        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Surface Flux Type: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Read surface flux variable units  ---
!
        IDFLT = 1
        VARB = 'Surface Flux Rate Variable Unit'
        CALL RDCHR(ISTART,ICOMMA,NCU,CHDUM,UNSF(1,NS))
        CALL RDSFUN( ISFT(NS) )
        VAR = 0.D+0
        INDX = 0
        CALL RDUNIT(UNSF(1,NS),VAR,INDX)
        IDFLT = 1
        VARB = 'Surface Flux Integral Variable Unit'
        CALL RDCHR(ISTART,ICOMMA,NCU,CHDUM,UNSF(2,NS))
        INDX = -ISFT(NS)
        CALL RDSFUN( INDX )
        VAR = 0.D+0
        INDX = 0
        CALL RDUNIT(UNSF(2,NS),VAR,INDX)
!
!---  Read surface flux orientation  ---
!
        VARB = 'Surface Flux Orientation'
        ISFSN(NS) = 0
!
!---    Non-fracture and non-borehole surface  ---
!
        IF( ISFT_FMX.EQ.0 .AND. ISFT_BHX.EQ.0 ) THEN
          IF( INDEX(ADUM(1:),'surface normal').NE.0 )  ISFSN(NS) = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
          WRITE(IWR,'(2A,$)') VARB(1:IVR),': '
          IF( INDEX(ADUM(1:),'west').NE.0) THEN
            ISFD(NS) = -1
            WRITE(IWR,'(A)') 'X-Direction: West Surface'
          ELSEIF( INDEX(ADUM(1:),'east').NE.0) THEN
            ISFD(NS) = 1
            WRITE(IWR,'(A)') 'X-Direction: East Surface'
          ELSEIF( INDEX(ADUM(1:),'south').NE.0) THEN
            ISFD(NS) = -2
            WRITE(IWR,'(A)') 'Y-Direction: South Surface'
          ELSEIF( INDEX(ADUM(1:),'north').NE.0) THEN
            ISFD(NS) = 2
            WRITE(IWR,'(A)') 'Y-Direction: North Surface'
          ELSEIF( INDEX(ADUM(1:),'bottom').NE.0) THEN
            ISFD(NS) = -3
            WRITE(IWR,'(A)') 'Z-Direction: Bottom Surface'
          ELSEIF( INDEX(ADUM(1:),'top').NE.0) THEN
            ISFD(NS) = 3
            WRITE(IWR,'(A)') 'Z-Direction: Top Surface'
          ELSEIF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            CALL RDCHR(ISTART,ICOMMA,NCHF,CHDUM,FDUM)
            NCHF = INDEX(FDUM,'  ')-1
!
!---        Check that surface flux domain file exists  ---
!
            INQUIRE( FILE=FDUM(1:NCHF), FORM=GDUM, EXIST=FLG_CHK )
            IF( .NOT.FLG_CHK ) THEN
              INDX = 4
              CHMSG = 'Surface-Flux-Domain File: '
     &          // FDUM(1:NCHF)
              CALL WRMSGS( INDX )
            ELSEIF( GDUM.EQ.'UNFORMATTED' ) THEN
              INDX = 4
              CHMSG = 'Unformatted Surface-Flux-Domain File: '
     &          // FDUM(1:NCHF)
              CALL WRMSGS( INDX )
            ENDIF
            OPEN(UNIT=26,FILE=FDUM(1:NCHF),STATUS='OLD',
     &        FORM='FORMATTED')
            WRITE(IWR,'(/,2A)') 'Surface-Flux-Domain File: ',
     &        FDUM(1:NCHF)
            ISFD(NS) = 4
            ISFC(1,NS) = 1
            ISFC(2,NS) = 1
            ISFC(3,NS) = 1
            ISFC(4,NS) = 1
            ISFC(5,NS) = 1
            ISFC(6,NS) = 1
            NC = 0
   30       CONTINUE
            READ(26,*,END=40) IX,JX,KX,ISFDX
            NC = NC + 1
            IF( NC.GT.LSFDOM ) THEN
              INDX = 5
              CHMSG = 'Number of Surface-Flux-Domain Surfaces ' //
     &          '> Parameter LSFDOM'
              CALL WRMSGS( INDX )
            ENDIF
            ISFDOM(1,NC,NS) = IX
            ISFDOM(2,NC,NS) = JX
            ISFDOM(3,NC,NS) = KX
            ISFDOM(4,NC,NS) = ISFDX
            GOTO 30
   40       CONTINUE
            NSFDOM(NS) = NC
            CLOSE(26)
          ENDIF
!
!---      Check surface flux domain  ---
!
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            DO 50 NC = 1,NSFDOM(NS)
              IX = ISFDOM(1,NC,NS)
              JX = ISFDOM(2,NC,NS)
              KX = ISFDOM(3,NC,NS)
              ISFDX = ISFDOM(4,NC,NS)
              IF( IX.LT.1 .OR. IX.GT.IFLD ) THEN
                INDX = 4
                CHMSG = 'Illegal Surface Flux Domain: I Index'
                CALL WRMSGS( INDX )
              ENDIF
              IF( JX.LT.1 .OR. JX.GT.JFLD ) THEN
                INDX = 4
                CHMSG = 'Illegal Surface Flux Domain: J Index'
                CALL WRMSGS( INDX )
              ENDIF
              IF( KX.LT.1 .OR. KX.GT.KFLD ) THEN
                INDX = 4
                CHMSG = 'Illegal Surface Flux Domain: K Index'
                CALL WRMSGS( INDX )
              ENDIF
   50       CONTINUE
!
!---      Read and check surface flux domain  ---
!
          ELSE
            VARB = 'Surface Flux Domain'
            CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(1,NS))
            CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(2,NS))
            CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(3,NS))
            CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(4,NS))
            CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(5,NS))
            CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(6,NS))
            IF( ISFC(1,NS).LT.1 .OR. ISFC(1,NS).GT.IFLD .OR.
     &        ISFC(2,NS).LT.1 .OR. ISFC(2,NS).GT.IFLD .OR.
     &        ISFC(1,NS).GT.ISFC(2,NS) ) THEN
              INDX = 4
              CHMSG = 'Illegal Surface Flux Domain: I Indices'
              CALL WRMSGS( INDX )
            ENDIF
            IF( ISFC(3,NS).LT.1 .OR. ISFC(3,NS).GT.JFLD .OR.
     &        ISFC(4,NS).LT.1 .OR. ISFC(4,NS).GT.JFLD .OR.
     &        ISFC(3,NS).GT.ISFC(4,NS) ) THEN
              INDX = 4
              CHMSG = 'Illegal Surface Flux Domain: J Indices'
              CALL WRMSGS( INDX )
            ENDIF
            IF( ISFC(5,NS).LT.1 .OR. ISFC(5,NS).GT.KFLD .OR.
     &        ISFC(6,NS).LT.1 .OR. ISFC(6,NS).GT.KFLD .OR.
     &        ISFC(5,NS).GT.ISFC(6,NS) ) THEN
              INDX = 4
              CHMSG = 'Illegal Surface Flux Domain: K Indices'
              CALL WRMSGS( INDX )
            ENDIF
            ISFC(1,NS) = MAX( 1,ISFC(1,NS) )
            ISFC(1,NS) = MIN( IFLD,ISFC(1,NS),ISFC(2,NS) )
            ISFC(2,NS) = MAX( 1,ISFC(1,NS),ISFC(2,NS) )
            ISFC(2,NS) = MIN( IFLD,ISFC(2,NS) )
            ISFC(3,NS) = MAX( 1,ISFC(3,NS) )
            ISFC(3,NS) = MIN( JFLD,ISFC(3,NS),ISFC(4,NS) )
            ISFC(4,NS) = MAX( 1,ISFC(3,NS),ISFC(4,NS) )
            ISFC(4,NS) = MIN( JFLD,ISFC(4,NS) )
            ISFC(5,NS) = MAX( 1,ISFC(5,NS) )
            ISFC(5,NS) = MIN( KFLD,ISFC(5,NS),ISFC(6,NS) )
            ISFC(6,NS) = MAX( 1,ISFC(5,NS),ISFC(6,NS) )
            ISFC(6,NS) = MIN( KFLD,ISFC(6,NS) )
            WRITE(IWR,'(/,2A)') VARB(1:IVR),': '
            WRITE(IWR,'(2X,2(A,I6))') 'I = ',ISFC(1,NS),' to ',
     &        ISFC(2,NS)
            WRITE(IWR,'(2X,2(A,I6))') 'J = ',ISFC(3,NS),' to ',
     &        ISFC(4,NS)
            WRITE(IWR,'(2X,2(A,I6))') 'K = ',ISFC(5,NS),' to ',
     &        ISFC(6,NS)
          ENDIF
!
!---    Fracture-matrix surface  ---
!
        ELSEIF( ISFT_FMX.EQ.1 ) THEN
          ISFD(NS) = 5
          VARB = 'Fracture-Matrix Surface Flux Fracture Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(1,NS))
          VARB = 'Fracture-Matrix Surface Flux Starting Triangle'
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(2,NS))
          VARB = 'Fracture-Matrix Surface Flux Ending Triangle'
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(3,NS))
          NFX = ISFC(1,NS)
          NT1X = ISFC(2,NS)
          NT2X = ISFC(3,NS)
          IF( NFX.LT.1 .OR. NFX.GT.NF_FRC ) THEN
            INDX = 7
            IMSG = ISFC(1,NS)
            CHMSG = 'Illegal Fracture-Matrix Surface Flux Domain: ' // 
     &        'Fracture Number: '
            CALL WRMSGS( INDX )
          ENDIF
          NTX = IP_FRC(2,NFX)-IP_FRC(1,NFX)+1
          IF( NT1X.LT.1 .OR. NT1X.GT.NTX ) THEN
            INDX = 7
            IMSG = ISFC(2,NS)
            CHMSG = 'Illegal Fracture-Matrix Surface Flux Domain: ' // 
     &        'Local Fracture Triangle: '
            CALL WRMSGS( INDX )
          ENDIF
          IF( NT1X.GT.NT2X .OR. NT2X.GT.NTX ) THEN
            INDX = 7
            IMSG = ISFC(3,NS)
            CHMSG = 'Illegal Fracture-Matrix Surface Flux Domain: ' // 
     &        'Local Fracture Triangle: '
            CALL WRMSGS( INDX )
          ENDIF
          WRITE(IWR,'(/,2A)') VARB(1:IVR),': '
          WRITE(IWR,'(2X,A,I6)') 'Fracture Number = ',ISFC(1,NS)
          WRITE(IWR,'(2X,2(A,I6))') 'Local Triangles = ',ISFC(2,NS),
     &      ' to ',ISFC(3,NS)
          WRITE(IWR,'(2X,2(A,I6))') 'Global Triangles = ',
     &      (IP_FRC(1,NFX)+ISFC(2,NS)-1),
     &      ' to ',(IP_FRC(1,NFX)+ISFC(3,NS)-1)
!
!---    Borehole surface  ---
!
        ELSEIF( ISFT_BHX.EQ.1 ) THEN
          ISFD(NS) = 6
          VARB = 'Borehole Surface Flux Borehole Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(1,NS))
          VARB = 'Borehole Surface Flux Borehole Node Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(2,NS))
          VARB = 'Borehole Surface Flux Connecting Borehole Node Number'
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(3,NS))
          NBHX = ISFC(1,NS)
          NBN1X = ISFC(2,NS)
          NBN2X = ISFC(3,NS)
          IF( NBHX.LT.1 .OR. NBHX.GT.N_BH ) THEN
            INDX = 7
            IMSG = ISFC(1,NS)
            CHMSG = 'Illegal Borehole Surface Flux Domain: ' // 
     &        'Borehole Number: '
            CALL WRMSGS( INDX )
          ENDIF
          NBNG1X = ID_BH(3,NBHX)-1+NBN1X
          NBNG2X = ID_BH(3,NBHX)-1+NBN2X
          IF( NBNG1X.LT.ID_BH(3,NBHX).OR.NBNG1X.GT.ID_BH(4,NBHX) ) THEN
            INDX = 7
            IMSG = ISFC(2,NS)
            CHMSG = 'Illegal Borehole Surface Flux Domain: ' // 
     &        'Local Borehole Number: '
            CALL WRMSGS( INDX )
          ENDIF
          DO NCX = IPB_BH(1,NBNG1X),IPB_BH(2,NBNG1X)
            IF( NCX.EQ.0 ) CYCLE
            IF( NBNG2X.EQ.IBCM_BH(NCX) ) EXIT
          ENDDO
          IF( NCX.GT.IPB_BH(2,NBNG1X) ) THEN
            INDX = 7
            IMSG = ISFC(3,NS)
            CHMSG = 'Illegal Borehole Surface Flux Domain: ' // 
     &        'Local Connecting Borehole Number: '
            CALL WRMSGS( INDX )
          ENDIF
          WRITE(IWR,'(/,2A)') VARB(1:IVR),': '
          WRITE(IWR,'(2X,A,I6)') 'Borehole Number = ',ISFC(1,NS)
          WRITE(IWR,'(2X,2(A,I6))') 'Local Borehole Numbers = ',
     &      ISFC(2,NS),' and ',ISFC(3,NS)
          WRITE(IWR,'(2X,2(A,I6))') 'Global Borehole Numbers = ',
     &      (ID_BH(3,NBHX)+ISFC(2,NS)-1),' and ',
     &      (ID_BH(3,NBHX)+ISFC(3,NS)-1)
        ENDIF
  100 CONTINUE
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSF_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSP_GT
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
!     Read input file for rock/soil saturation function information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 February 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TABL
      USE SOLTN
      USE PORMED
      USE GRID
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
      REAL*8 SCHRX(LSCHR)
      CHARACTER*4 FORM
      CHARACTER*64 ADUM,RDUM,UNTS
      CHARACTER*512 CHDUM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IVAR,JVAR
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM
      DATA FORM /'(I9)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSP_GT'
      ALLOCATE( IVAR(1:LRC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: IVAR'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( JVAR(1:LRC),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: JVAR'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Set oven-dried head to Webb value (1.e+7 cm)  ---
!
      HDOD = 1.D+5
!
!---  Write card information to ouput file  ---
!
      CARD = 'Rock/Soil Saturation Function Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Loop over the rock/soil saturation information lines  ---
!
      NR = 0
      IJK = 0
   10 CONTINUE
      IF( NR.GE.NROCK .OR. IJK.GT.0 ) GOTO 500
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Rock/Soil Name'
      DO L = 1,LSCHR
        SCHRX(L) = 0.D+0
      ENDDO
   12 CONTINUE
!
!---  Rock/Soil option for IJK Indexing  ---
!
      IF( IJK.LT.0 ) THEN
        RDUM = 'Rock/Soil #'
        NR = NR + 1
        ICX = ICOUNT(NR)
        WRITE( FORM(3:3),'(I1)' ) ICX
        NCH = 12 + ICX - 1
        WRITE( RDUM(12:NCH),FORM ) NR
        WRITE (IWR,'(/,A)') RDUM(1:NCH)
        GOTO 220
!
!---  Read rock/soil name  ---
!
      ELSE
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,RDUM)
      ENDIF
!
!---  Check for a pair of delimiting slashes in the rock/soil name,
!     indicating a pattern of rock/soil types  ---
!
      KBS = 0
      IBS = INDEX( RDUM(1:),'/' )
      IF( IBS.GT.0 ) THEN
        IBS = IBS + 1
        JBS = INDEX( RDUM(IBS:),'/')
        IF( JBS.GT.0 ) THEN
          JBS = IBS + JBS - 2
          KBS = 1
          ISBS = ISTART
        ENDIF
      ENDIF
      IROCK = 1
   20 CONTINUE
!
!---  IJK, KIJ, or JKI indexing  ---
!
      IF( INDEX(RDUM(1:),'indexing').NE.0 ) THEN
        IF( INDEX(ROCK(1)(1:),'indexing').EQ.0 ) THEN
          INDX = 4
          CHMSG = 'Indexing Option Not Declared ' // 
     &      'in Rock/Soil Zonation Card'
          CALL WRMSGS( INDX )
        ENDIF
        IF( INDEX(RDUM,'ijk').NE.0 ) THEN
          IJK = 1
        ELSEIF( INDEX(RDUM,'jki').NE.0 ) THEN
          IJK = 2
        ELSEIF( INDEX(RDUM,'kij').NE.0 ) THEN
          IJK = 3
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Indexing Option' // RDUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        IROCK = 1
        GOTO 220
      ENDIF
!
!---  Search known rock types for a matching type ---
!
      DO 100 M = IROCK,NROCK
        IF( KBS.EQ.1 ) THEN
          IF( INDEX( ROCK(M)(1:),RDUM(IBS:JBS) ).GT.0 ) THEN
            IROCK = M
            GOTO 200
          ENDIF
        ELSE
          IF( RDUM.EQ.ROCK(M) ) THEN
            IROCK = M
            GOTO 200
          ENDIF
        ENDIF
  100 CONTINUE
!
!---  Search known scaling groups for a matching type ---
!
      IF( ISLC(19).EQ.1 ) THEN
        DO 110 M = 1,NSCALE
           IF( RDUM.EQ.SCALNM(M) ) THEN
              ISGRP = M
              IROCK = 1
              GOTO 200
           ENDIF
  110   CONTINUE
        INDX = 2
        CHMSG = 'Unrecognized Rock/Soil Type or Scaling Group: '
     &    // RDUM(1:NCH)
        CALL WRMSGS( INDX )
        GOTO 10
      ENDIF
      INDX = 2
      CHMSG = 'Unrecognized Rock/Soil Type: ' // RDUM(1:NCH)
      CALL WRMSGS( INDX )
      GOTO 10
  200 CONTINUE
!
!---  Loop over rock/soils within scaling group  ---
!
      IF( ISLC(19).EQ.1 .AND. ISGRP.NE.0 ) THEN
        DO 202 M = IROCK,NROCK
          IF( ISCALE(M).EQ.ISGRP ) THEN
            IROCK = M
            GOTO 204
          ENDIF
  202   CONTINUE
      ENDIF
  204 CONTINUE
!
!---  Write rock/soil name  ---
!
      WRITE (IWR,'(/,2A)') 'Rock/Soil Name: ',ROCK(IROCK)
!
!---  Read aqueous relative permeability pressure function  ---
!
      NR = NR + 1
  220 CONTINUE
!
!---  Rock/Soil option for IJK Indexing, dissociate from
!     saturation function type  ---
!
      IF( IJK.LT.0 ) THEN
        ISCHRX = 0
      ELSE
        ISCHRX = MOD( ISCHR(IROCK),1000 )
      ENDIF
!
!---  Read saturation/capillary pressure function  ---
!
      VARB = 'Saturation Function Type: '
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---  Rock/Soil Zonation Option for IJK Indexing  ---
!
      IF( IJK.GT.0 .AND. INDEX(ADUM(1:),'rock/soil').NE.0 ) THEN
        VARB = 'Number of Rock/Soil Entries'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NROCK)
        CALL RDIJKI( ISTART,IJK,CHDUM,IZ2 )
        IJK = -IJK
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        GOTO 12
      ENDIF
!
!---  Extension model  ---
!
      ISMX = 0
      IF( INDEX(ADUM(1:),'webb').NE.0 .OR.
     &  INDEX(ADUM(1:),'extension').NE.0 .OR.
     &  INDEX(ADUM(1:),'extended').NE.0 ) THEN
        ISMX = 2
      ENDIF
!
!---  Gas entrapment options  ---
!
      IF( INDEX(ADUM(1:),'entrap').NE.0 ) THEN
        IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
          ISCHRX = 101
        ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &    INDEX(ADUM(1:),'corey').NE.0 ) THEN
          ISCHRX = 102
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Saturation Function: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Drainage-imbibition options  ---
!
      ELSEIF( INDEX(ADUM(1:),'imbibition').NE.0 .AND. 
     &  INDEX(ADUM(1:),'drainage').NE.0 ) THEN
        IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
          ISCHRX = 201
        ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &    INDEX(ADUM(1:),'corey').NE.0 ) THEN
          ISCHRX = 202
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Saturation Function: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Webb extension w/o gas entrapment  ---
!
      ELSEIF( ISMX.EQ.2 ) THEN
        IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
          ISCHRX = 1
          IF( INDEX(RDUM(1:),'fractured').NE.0 .OR.
     &      INDEX(RDUM(1:),'dp').NE.0 .OR.
     &      INDEX(RDUM(1:),'dual').NE.0 ) ISCHRX = 3
        ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &    INDEX(ADUM(1:),'corey').NE.0 ) THEN
          ISCHRX = 2
          IF( INDEX(RDUM(1:),'fractured').NE.0 .OR.
     &      INDEX(RDUM(1:),'dp').NE.0 .OR.
     &      INDEX(RDUM(1:),'dual').NE.0 ) ISCHRX = 4
        ELSE
          INDX = 4
          CHMSG = 'Saturation Function Conflict w/ Webb Extension: '
     &      // ADUM
          CALL WRMSGS( INDX )
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
        ISCHRX = 1
        IF( INDEX(RDUM(1:),'fractured').NE.0 .OR.
     &    INDEX(RDUM(1:),'dp').NE.0 .OR.
     &    INDEX(RDUM(1:),'dual').NE.0 ) ISCHRX = 3
      ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &    INDEX(ADUM(1:),'corey').NE.0 ) THEN
        ISCHRX = 2
        IF( INDEX(RDUM(1:),'fractured').NE.0 .OR.
     &    INDEX(RDUM(1:),'dp').NE.0 .OR.
     &    INDEX(RDUM(1:),'dual').NE.0 ) ISCHRX = 4
      ELSEIF( INDEX(ADUM(1:),'haverkamp').NE.0 ) THEN
        ISCHRX = 5
      ELSEIF( INDEX(ADUM(1:),'russo').NE.0 ) THEN
        ISCHRX = 9
      ELSEIF( INDEX(ADUM(1:),'tabular').NE.0 ) THEN
        IF( INDEX( ADUM(1:),'log').NE.0 ) THEN
          IF( INDEX( ADUM(1:),'hysteretic').NE.0 ) THEN
            ISCHRX = 13
          ELSE
            ISCHRX = 11
          ENDIF
        ELSE
          IF( INDEX( ADUM(1:),'hysteretic').NE.0 ) THEN
            ISCHRX = 12
          ELSE
            ISCHRX = 10
          ENDIF
        ENDIF
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Saturation Function: '//ADUM
        CALL WRMSGS( INDX )
      ENDIF
!
!---  van Genuchten Function  ---
!
      IF( ISCHRX.EQ.1 .OR. ISCHRX.EQ.101 .OR. ISCHRX.EQ.201 ) THEN
        IF( ISCHRX.EQ.1 ) THEN
          WRITE(IWR,'(A)') 'van Genuchten Function'
        ELSEIF( ISCHRX.EQ.101 ) THEN
          WRITE(IWR,'(A)') 'van Genuchten Function w/ Gas Entrapment'
        ELSEIF( ISCHRX.EQ.201 ) THEN
          WRITE(IWR,'(A)') 'van Genuchten Drainage-Imbibition'
        ENDIF
        IF( ISMX.EQ.2 ) WRITE(IWR,'(A)') 'Webb Extension Model'
        IF( ISCHRX.EQ.201 ) THEN
          VARB = 'van Genuchten (alpha) Drainage'
        ELSE
          VARB = 'van Genuchten (alpha)'
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 1
          LNDX = LSCHR
          UNTS = '1/m'
          IUNM = -1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHRX(1)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,SCHRX(1),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(1),', 1/m)'
        ENDIF
        IF( ISCHRX.EQ.201 ) THEN
          VARB = 'van Genuchten (n) Drainage'
        ELSE
          VARB = 'van Genuchten (n)'
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 3
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(3))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHRX(3)
        ENDIF
        IF( ISCHRX.EQ.201 ) THEN
          VARB = 'van Genuchten (residual saturation) Drainage'
        ELSE
          VARB = 'van Genuchten (residual saturation)'
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 4
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(4))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHRX(4)
        ENDIF
        IF( IJK.GT.0 ) THEN
          DO 240 N = 1,NFLD
            IF( (1.D+0-SCHR(4,IZ(N))).LT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Excessive ' // VARB(1:IVR)
              RLMSG = SCHR(4,IZ(N))
              CALL WRMSGS( INDX )
            ENDIF
  240     CONTINUE
        ELSE
          IF( (1.D+0-SCHRX(4)).LT.EPSL ) THEN
            INDX = 9
            CHMSG = 'Excessive ' // VARB(1:IVR)
            RLMSG = SCHRX(4)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        IF( ISCHRX.EQ.201 ) THEN
          VARB = 'van Genuchten (m) Drainage'
        ELSE
          VARB = 'van Genuchten (m)'
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 14
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(14))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(14)
        ENDIF
        IF( ISCHRX.EQ.201 ) THEN
          VARB = 'van Genuchten (alpha) Imbibition'
          IF( IJK.GT.0 ) THEN
            INDX = 2
            LNDX = LSCHR
            UNTS = '1/m'
            IUNM = -1
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHRX(2)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,SCHRX(2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(2),', 1/m)'
          ENDIF
          VARB = 'van Genuchten (n) Imbibition'
          IF( IJK.GT.0 ) THEN
            INDX = 5
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(5))
            WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHRX(5)
          ENDIF
          VARB = 'van Genuchten (m) Imbibition'
          IF( IJK.GT.0 ) THEN
            INDX = 13
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(13))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(13)
          ENDIF
        ENDIF
        IF( ISCHRX.EQ.201 ) THEN
          VARB = 'van Genuchten (Actual Maximum Trapped Gas Saturation)'
        ELSE
          VARB = 'van Genuchten (Actual Gas Residual Saturation)'
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 15
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(15))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &      SCHRX(15)
        ENDIF
!
!---    Webb extension read optional oven dried head  ---
!
        CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
        IF( ISMX.EQ.2 .AND. INDX.EQ.1 ) THEN
          VARB = 'Webb Extension Oven Dried Head'
          IF( IJK.GT.0 ) THEN
            INDX = 12
            LNDX = LSCHR
            UNTS = 'm'
            IUNM = 1
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(12))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHRX(12)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,SCHRX(12),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(12),', m)'
          ENDIF
        ELSE
          IF( IJK.GT.0 ) THEN
            DO N = 1,NFLD
              SCHR(12,N) = HDOD
            ENDDO
          ELSE
            SCHRX(12) = HDOD
          ENDIF        
        ENDIF
!
!---  Brooks and Corey Function  ---
!
      ELSEIF( ISCHRX.EQ.2 .OR. ISCHRX.EQ.102 .OR. ISCHRX.EQ.202 ) THEN
        IF( ISCHRX.EQ.2 ) THEN
          WRITE(IWR,'(A)') 'Brooks and Corey Function'
        ELSEIF( ISCHRX.EQ.102 ) THEN
          WRITE(IWR,'(A)') 'Brooks and Corey Function w/ Gas Entrapment'
        ELSEIF( ISCHRX.EQ.202 ) THEN
          WRITE(IWR,'(A)') 'Brooks and Corey Drainage-Imbibition'
        ENDIF
        IF( ISMX.EQ.2 ) WRITE(IWR,'(A)') 'Webb Extension Model'
        IF( ISCHRX.EQ.202 ) THEN
          VARB = 'Brooks and Corey (psi) Drainage'
        ELSE
          VARB = 'Brooks and Corey (psi)'
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 1
          LNDX = LSCHR
          UNTS = 'm'
          IUNM = 1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHRX(1)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SCHRX(1),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(1),', m)'
        ENDIF
        IF( ISCHRX.EQ.202 ) THEN
          VARB = 'Brooks and Corey (lambda) Drainage'
        ELSE
          VARB = 'Brooks and Corey (lambda)'
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 3
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(3))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHRX(3)
        ENDIF
        IF( ISCHRX.EQ.202 ) THEN
          VARB = 'Brooks and Corey (Aqueous Residual Sat.) Drainage'
        ELSE
          VARB = 'Brooks and Corey (Aqueous Residual Saturation)'
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 4
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(4))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHRX(4)
        ENDIF
        IF( IJK.GT.0 ) THEN
          DO 244 N = 1,NFLD
            IF( (1.D+0-SCHR(4,IZ(N))).LT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Excessive ' // VARB(1:IVR)
              RLMSG = SCHR(4,IZ(N))
              CALL WRMSGS( INDX )
            ENDIF
  244     CONTINUE
        ELSE
          IF( (1.D+0-SCHRX(4)).LT.EPSL ) THEN
            INDX = 9
              CHMSG = 'Excessive ' // VARB(1:IVR)
            RLMSG = SCHRX(4)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        IF( ISCHRX.EQ.202 ) THEN
          VARB = 'Brooks and Corey (psi) Imbibition'
          IF( IJK.GT.0 ) THEN
            INDX = 5
            LNDX = LSCHR
            UNTS = 'm'
            IUNM = 1
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(5))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHRX(5)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,SCHRX(5),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(5),', m)'
          ENDIF
          VARB = 'Brooks and Corey (lambda) Imbibition'
          IF( IJK.GT.0 ) THEN
            INDX = 6
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(6))
            WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHRX(6)
          ENDIF
        ENDIF
        IF( ISCHRX.EQ.202 ) THEN
          VARB = 'Brooks and Corey (Actual Maximum Trapped Gas Sat.)'
        ELSE
          VARB = 'Brooks and Corey (Actual Gas Residual Saturation)'
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 15
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(15))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &      SCHRX(15)
        ENDIF
!
!---    Webb extension read optional oven dried head  ---
!
        CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
        IF( ISMX.EQ.2 .AND. INDX.EQ.1 ) THEN
          VARB = 'Webb Extension Oven Dried Head'
          IF( IJK.GT.0 ) THEN
            INDX = 12
            LNDX = LSCHR
            UNTS = 'm'
            IUNM = 1
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(12))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHRX(12)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,SCHRX(12),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(12),', m)'
          ENDIF
        ELSE
          IF( IJK.GT.0 ) THEN
            DO N = 1,NFLD
              SCHR(12,N) = HDOD
            ENDDO
          ELSE
            SCHRX(12) = HDOD
          ENDIF        
        ENDIF
!
!---  Dual Porosity van Genuchten Function  ---
!
      ELSEIF( ISCHRX.EQ.3 ) THEN
        WRITE(IWR,'(A)') 'Dual Porosity van Genuchten Function'
        VARB = 'Matrix van Genuchten (alpha)'
        IF( IJK.GT.0 ) THEN
          INDX = 1
          LNDX = LSCHR
          UNTS = '1/m'
          IUNM = -1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHRX(1)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,SCHRX(1),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(1),', 1/m)'
        ENDIF
        VARB = 'Matrix van Genuchten (n): '
        IF( IJK.GT.0 ) THEN
          INDX = 3
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(3))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHRX(3)
        ENDIF
        VARB = 'Matrix van Genuchten (residual saturation): '
        IF( IJK.GT.0 ) THEN
          INDX = 4
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(4))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHRX(4)
        ENDIF
        IF( IJK.GT.0 ) THEN
          DO 248 N = 1,NFLD
            IF( (1.D+0-SCHR(4,IZ(N))).LT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Excessive ' // VARB(1:IVR)
              RLMSG = SCHR(4,IZ(N))
              CALL WRMSGS( INDX )
            ENDIF
  248     CONTINUE
        ELSE
          IF( (1.D+0-SCHRX(4)).LT.EPSL ) THEN
            INDX = 9
            CHMSG = 'Excessive ' // VARB(1:IVR)
            RLMSG = SCHRX(4)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        VARB = 'Fracture van Genuchten (alpha), '
        IF( IJK.GT.0 ) THEN
          INDX = 5
          LNDX = LSCHR
          UNTS = '1/m'
          IUNM = -1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(5))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHRX(5)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,SCHRX(5),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(5),', 1/m)'
        ENDIF
        VARB = 'Fracture van Genuchten (n): '
        IF( IJK.GT.0 ) THEN
          INDX = 6
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(6))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHRX(6)
        ENDIF
        VARB = 'Fracture van Genuchten (residual saturation): '
        IF( IJK.GT.0 ) THEN
          INDX = 7
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(7))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(7)
        ENDIF
        VARB = 'Matrix van Genuchten (m)'
        IF( IJK.GT.0 ) THEN
          INDX = 14
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(14))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(14)
        ENDIF
        VARB = 'Fracture van Genuchten (m)'
        IF( IJK.GT.0 ) THEN
          INDX = 15
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(15))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(15)
        ENDIF
!
!---    Webb extension read optional oven dried head  ---
!
        CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
        IF( ISMX.EQ.2 .AND. INDX.EQ.1 ) THEN
          VARB = 'Webb Extension Oven Dried Head'
          IF( IJK.GT.0 ) THEN
            INDX = 12
            LNDX = LSCHR
            UNTS = 'm'
            IUNM = 1
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(12))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHRX(12)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,SCHRX(12),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(12),', m)'
          ENDIF
        ELSE
          IF( IJK.GT.0 ) THEN
            DO N = 1,NFLD
              SCHR(12,N) = HDOD
            ENDDO
          ELSE
            SCHRX(12) = HDOD
          ENDIF        
        ENDIF
!
!---  Dual Porosity Brooks and Corey Function  ---
!
      ELSEIF( ISCHRX.EQ.4 ) THEN
        WRITE(IWR,'(A)') 'Dual Porosity Brooks and Corey Function'
        VARB = 'Matrix Brooks and Corey (psi)'
        IF( IJK.GT.0 ) THEN
          INDX = 1
          LNDX = LSCHR
          UNTS = 'm'
          IUNM = 1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHRX(1)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SCHRX(1),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(1),', m)'
        ENDIF
        VARB = 'Matrix Brooks and Corey (lambda): '
        IF( IJK.GT.0 ) THEN
          INDX = 3
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(3))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(3)
        ENDIF
        VARB = 'Matrix Brooks and Corey (residual saturation): '
        IF( IJK.GT.0 ) THEN
          INDX = 4
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(4))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(4)
        ENDIF
        IF( IJK.GT.0 ) THEN
          DO 250 N = 1,NFLD
            IF( (1.D+0-SCHR(4,IZ(N))).LT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Excessive ' // VARB(1:IVR)
              RLMSG = SCHR(4,IZ(N))
              CALL WRMSGS( INDX )
            ENDIF
  250     CONTINUE
        ELSE
          IF( (1.D+0-SCHRX(4)).LT.EPSL ) THEN
            INDX = 9
            CHMSG = 'Excessive ' // VARB(1:IVR)
            RLMSG = SCHRX(4)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        VARB = 'Fracture Brooks and Corey (psi)'
        IF( IJK.GT.0 ) THEN
          INDX = 5
          LNDX = LSCHR
          UNTS = 'm'
          IUNM = 1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(5))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHRX(5)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SCHRX(5),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(5),', m)'
        ENDIF
        VARB = 'Fracture Brooks and Corey (lambda): '
        IF( IJK.GT.0 ) THEN
          INDX = 6
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(6))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(6)
        ENDIF
        VARB = 'Fracture Brooks and Corey (residual saturation): '
        IF( IJK.GT.0 ) THEN
          INDX = 7
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(7))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(7)
        ENDIF
!
!---    Webb extension read optional oven dried head  ---
!
        CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
        IF( ISMX.EQ.2 .AND. INDX.EQ.1 ) THEN
          VARB = 'Webb Extension Oven Dried Head'
          IF( IJK.GT.0 ) THEN
            INDX = 12
            LNDX = LSCHR
            UNTS = 'm'
            IUNM = 1
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(12))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SCHRX(12)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,SCHRX(12),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(12),', m)'
          ENDIF
        ELSE
          IF( IJK.GT.0 ) THEN
            DO N = 1,NFLD
              SCHR(12,N) = HDOD
            ENDDO
          ELSE
            SCHRX(12) = HDOD
          ENDIF        
        ENDIF
!
!---  Haverkamp Function  ---
!
      ELSEIF( ISCHRX.EQ.5 ) THEN
        WRITE(IWR,'(A)') 'Haverkamp Function'
        VARB = 'Haverkamp Gas Entry Head (psi): '
        IF( IJK.GT.0 ) THEN
          INDX = 1
          LNDX = LSCHR
          UNTS = 'm'
          IUNM = 1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHRX(1)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SCHRX(1),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(1),', m)'
        ENDIF
        INDX = 0
        IUNM = 1
        VAR = 1.D+0
        CALL RDUNIT(UNTS,VAR,INDX)
        IF( IJK.GT.0 ) THEN
          DO 260 N = 1,NFLD
            SCHR(5,IZ(N)) = VAR
  260     CONTINUE
        ELSE
          SCHRX(5) = VAR
        ENDIF
        VARB = 'Haverkamp (alpha): '
        IF( IJK.GT.0 ) THEN
          INDX = 2
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(2))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(2)
        ENDIF
        VARB = 'Haverkamp (beta): '
        IF( IJK.GT.0 ) THEN
          INDX = 3
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(3))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(3)
        ENDIF
        VARB = 'Haverkamp (residual saturation): '
        IF( IJK.GT.0 ) THEN
          INDX = 4
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(4))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(4)
        ENDIF
        IF( IJK.GT.0 ) THEN
          DO 262 N = 1,NFLD
            IF( (1.D+0-SCHR(4,IZ(N))).LT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Excessive ' // VARB(1:IVR)
              RLMSG = SCHR(4,IZ(N))
              CALL WRMSGS( INDX )
            ENDIF
  262     CONTINUE
        ELSE
          IF( (1.D+0-SCHRX(4)).LT.EPSL ) THEN
            INDX = 9
            CHMSG = 'Excessive ' // VARB(1:IVR)
            RLMSG = SCHRX(4)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Russo Function  ---
!
      ELSEIF( ISCHRX.EQ.9 ) THEN
        WRITE(IWR,'(A)') 'Russo Function'
        VARB = 'Russo (alpha)'
        IF( IJK.GT.0 ) THEN
          INDX = 1
          LNDX = LSCHR
          UNTS = '1/m'
          IUNM = -1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHRX(1)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,SCHRX(1),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHRX(1),', 1/m)'
        ENDIF
        VARB = 'Russo (n): '
        IF( IJK.GT.0 ) THEN
          INDX = 3
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(3))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(3)
        ENDIF
        VARB = 'Russo (residual saturation): '
        IF( IJK.GT.0 ) THEN
          INDX = 4
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHRX(4))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHRX(4)
        ENDIF
        IF( IJK.GT.0 ) THEN
          DO 264 N = 1,NFLD
            IF( (1.D+0-SCHR(4,IZ(N))).LT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Excessive ' // VARB(1:IVR)
              RLMSG = SCHR(4,IZ(N))
              CALL WRMSGS( INDX )
            ENDIF
  264     CONTINUE
        ELSE
          IF( (1.D+0-SCHRX(4)).LT.EPSL ) THEN
            INDX = 9
            CHMSG = 'Excessive ' // VARB(1:IVR)
            RLMSG = SCHRX(4)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Tabular non-hysteretic  ---
!
      ELSEIF( ISCHRX.EQ.10 .OR. ISCHRX.EQ.11 ) THEN
        WRITE(IWR,'(A)') 'Tabular Aqueous Saturation versus ' //
     &    'Scaled Gas-Water Capillary Head'
        IF( ISCHRX.EQ.11 ) THEN
          WRITE(IWR,'(A)') 'Log-Linear Interpolation'
        ELSE
          WRITE(IWR,'(A)') 'Linear Interpolation'
        ENDIF
!
!---    IJK Indexing  ---
!
        IF( IJK.GT.0 ) THEN
          VARB = 'Number of Tables'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NTABLX)
          IF( NTABLX.LT.1 ) THEN
            INDX = 4
            CHMSG = 'Invalid Number of Saturation Tables'
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDIJKI( ISTART,IJK,CHDUM,IVAR )
!
!---      Loop over saturation function tables  ---
!
          DO 274 NTX = 1,NTABLX
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
            VARB = 'Number of Table Entries'
            CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
            WRITE(IWR,'(2A,I6)') VARB(1:IVR),': ',NLIN
!
!---        Loop over lines in saturation function tables  ---
!
            NTBLX = NTBL+1
            DO 270 NL = 1,NLIN
              NTBL = NTBL + 1
              IF( NTBL.GT.LTBL ) THEN
                INDX = 5
                CHMSG = 'Number of Table Values > Parameter LTBL'
                CALL WRMSGS( INDX )
              ENDIF
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              ISTART = 1
              VARB = 'Capillary Head'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLX(NTBL))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',TBLX(NTBL)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,TBLX(NTBL),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',TBLX(NTBL),', m)'
!
!---          Log-linear interpolation  ---
!
              IF( ISCHRX.EQ.11 ) THEN
                TBLX(NTBL) = LOG( MAX( TBLX(NTBL),1.D-14 ) )
              ENDIF
              VARB = 'Aqueous Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',TBLY(NTBL)
              IF( NL.EQ.2 ) THEN
                IF( TBLX(NTBL-1).LT.TBLX(NTBL) ) THEN
                  ITDX = 1
                ELSEIF( TBLX(NTBL-1).GT.TBLX(NTBL) ) THEN
                  ITDX = -1
                ELSE
                  INDX = 4
                  CHMSG = 'Invalid Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
                IF( TBLY(NTBL-1).LT.TBLY(NTBL) ) THEN
                  ITDY = 1
                ELSEIF( TBLY(NTBL-1).GT.TBLY(NTBL) ) THEN
                  ITDY = -1
                ELSE
                  INDX = 4
                  CHMSG = 'Invalid Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
              ELSEIF( NL.GT.2 ) THEN
                IF( (ITDX.EQ.1 .AND. TBLX(NTBL).LE.TBLX(NTBL-1)) .OR.
     &            (ITDX.EQ.-1 .AND. TBLX(NTBL).GE.TBLX(NTBL-1)) ) THEN
                  INDX = 4
                  CHMSG = 'Invalid Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
                IF( (ITDY.EQ.1 .AND. TBLY(NTBL).LE.TBLY(NTBL-1)) .OR.
     &            (ITDY.EQ.-1 .AND. TBLY(NTBL).GE.TBLY(NTBL-1)) ) THEN
                  INDX = 4
                  CHMSG = 'Invalid Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
  270       CONTINUE
!
!---        Correlate table numbers with nodes  ---
!
            DO 272 N = 1,NFLD
              IF( IVAR(N).EQ.NTX ) THEN
                ISLTBL(1,N) = NTBLX
                ISLTBL(2,N) = NTBL
              ELSEIF( IVAR(N).LT.1 .OR. IVAR(N).GT.NTABLX ) THEN
                INDX = 4
                CHMSG = 'Invalid Saturation Table Number'
                CALL WRMSGS( INDX )
              ENDIF
 272        CONTINUE
 274      CONTINUE
!
!---    Rock/soil zonation  ---
!
        ELSE
          VARB = 'Number of Table Entries'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
          WRITE(IWR,'(2A,I6)') VARB(1:IVR),': ',NLIN
          IF( NLIN.LT.2 ) THEN
            INDX = 4
            CHMSG = 'Saturation Invalid Table'
            CALL WRMSGS( INDX )
          ENDIF
          ISLTBL(1,IROCK) = NTBL + 1
          DO 276 NL = 1,NLIN
            NTBL = NTBL + 1
            IF( NTBL.GT.LTBL ) THEN
              INDX = 5
              CHMSG = 'Number of Table Values > Parameter LTBL'
              CALL WRMSGS( INDX )
            ENDIF
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
            VARB = 'Capillary Head'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLX(NTBL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',TBLX(NTBL)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,TBLX(NTBL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',TBLX(NTBL),', m)'
!
!---        Log-linear interpolation  ---
!
            IF( ISCHRX.EQ.11 ) THEN
              TBLX(NTBL) = LOG( MAX( TBLX(NTBL),1.D-14 ) )
            ENDIF
            VARB = 'Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',TBLY(NTBL)
            IF( NL.EQ.2 ) THEN
              IF( TBLX(NTBL-1).LT.TBLX(NTBL) ) THEN
                ITDX = 1
              ELSEIF( TBLX(NTBL-1).GT.TBLX(NTBL) ) THEN
                ITDX = -1
              ELSE
                INDX = 4
                CHMSG = 'Invalid Saturation Table'
                CALL WRMSGS( INDX )
              ENDIF
              IF( TBLY(NTBL-1).LT.TBLY(NTBL) ) THEN
                ITDY = 1
              ELSEIF( TBLY(NTBL-1).GT.TBLY(NTBL) ) THEN
                ITDY = -1
              ELSE
                INDX = 4
                CHMSG = 'Invalid Saturation Table'
                CALL WRMSGS( INDX )
              ENDIF
            ELSEIF( NL.GT.2 ) THEN
              IF( (ITDX.EQ.1 .AND. TBLX(NTBL).LE.TBLX(NTBL-1)) .OR.
     &          (ITDX.EQ.-1 .AND. TBLX(NTBL).GE.TBLX(NTBL-1)) ) THEN
                INDX = 4
                CHMSG = 'Invalid Saturation Table'
                CALL WRMSGS( INDX )
              ENDIF
              IF( (ITDY.EQ.1 .AND. TBLY(NTBL).LE.TBLY(NTBL-1)) .OR.
     &          (ITDY.EQ.-1 .AND. TBLY(NTBL).GE.TBLY(NTBL-1)) ) THEN
                INDX = 4
                CHMSG = 'Invalid Saturation Table'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
  276     CONTINUE
          ISLTBL(2,IROCK) = NTBL
        ENDIF
!
!---  Tabular hysteretic  ---
!
      ELSEIF( ISCHRX.EQ.12 .OR. ISCHRX.EQ.13 ) THEN
        WRITE(IWR,'(A)') 'Tabular Hysteretic Aqueous Saturation ' //
     &    ' versus Scaled Gas-Water Capillary Head'
        IF( ISCHRX.EQ.13 ) THEN
          WRITE(IWR,'(A)') 'Log-Linear Interpolation'
        ELSE
          WRITE(IWR,'(A)') 'Linear Interpolation'
        ENDIF
!
!---    IJK Indexing  ---
!
        IF( IJK.GT.0 ) THEN
          VARB = 'Number of Drainage Saturation Tables'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NDTABLX)
          IF( NDTABLX.LT.1 ) THEN
            INDX = 4
            CHMSG = 'Invalid Number of Drainage Saturation Tables'
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDIJKI( ISTART,IJK,CHDUM,IVAR )
          VARB = 'Number of Imbibition Saturation Tables'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NITABLX)
          IF( NITABLX.LT.1 ) THEN
            INDX = 4
            CHMSG = 'Invalid Number of Imbibition Saturation Tables'
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDIJKI( ISTART,IJK,CHDUM,JVAR )
!
!---      Loop over drainage saturation function tables  ---
!
          DO 284 NTX = 1,NDTABLX
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
            VARB = 'Number of Drainage Table Entries'
            CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
            WRITE(IWR,'(2A,I6)') VARB(1:IVR),': ',NLIN
!
!---        Loop over lines in saturation function tables  ---
!
            NTBLX = NTBL+1
            DO 280 NL = 1,NLIN
              NTBL = NTBL + 1
              IF( NTBL.GT.LTBL ) THEN
                INDX = 5
                CHMSG = 'Number of Table Entries > Parameter LTBL'
                CALL WRMSGS( INDX )
              ENDIF
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              ISTART = 1
              VARB = 'Drainage Capillary Head'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLX(NTBL))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',TBLX(NTBL)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,TBLX(NTBL),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',TBLX(NTBL),', m)'
!
!---          Log-linear interpolation  ---
!
              IF( ISCHRX.EQ.13 ) THEN
                TBLX(NTBL) = LOG( MAX( TBLX(NTBL),1.D-14 ) )
              ENDIF
              VARB = 'Drainage Aqueous Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',TBLY(NTBL)
              IF( NL.EQ.2 ) THEN
                IF( TBLX(NTBL-1).LT.TBLX(NTBL) ) THEN
                  ITDX = 1
                ELSEIF( TBLX(NTBL-1).GT.TBLX(NTBL) ) THEN
                  ITDX = -1
                ELSE
                  INDX = 4
                  CHMSG = 'Invalid Drainage Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
                IF( TBLY(NTBL-1).LT.TBLY(NTBL) ) THEN
                  ITDY = 1
                ELSEIF( TBLY(NTBL-1).GT.TBLY(NTBL) ) THEN
                  ITDY = -1
                ELSE
                  INDX = 4
                  CHMSG = 'Invalid Drainage Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
              ELSEIF( NL.GT.2 ) THEN
                IF( (ITDX.EQ.1 .AND. TBLX(NTBL).LE.TBLX(NTBL-1)) .OR.
     &            (ITDX.EQ.-1 .AND. TBLX(NTBL).GE.TBLX(NTBL-1)) ) THEN
                  INDX = 4
                  CHMSG = 'Invalid Drainage Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
                IF( (ITDY.EQ.1 .AND. TBLY(NTBL).LE.TBLY(NTBL-1)) .OR.
     &            (ITDY.EQ.-1 .AND. TBLY(NTBL).GE.TBLY(NTBL-1)) ) THEN
                  INDX = 4
                  CHMSG = 'Invalid Drainage Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
  280       CONTINUE
!!
!!---        Build cubic splines  ---
!!
!            IF( ISCHRX.EQ.13 ) THEN
!              CALL SPLINY( NTBLX,NTBL )
!              CALL SPLINX( NTBLX,NTBL )
!            ENDIF
!
!---        Correlate table numbers with nodes  ---
!
            DO 282 N = 1,NFLD
              IF( IVAR(N).EQ.NTX ) THEN
                ISLTBL(1,N) = NTBLX
                ISLTBL(2,N) = NTBL
              ELSEIF( IVAR(N).LT.1 .OR. IVAR(N).GT.NDTABLX ) THEN
                INDX = 4
                CHMSG = 'Invalid Drainage Saturation Table Number'
                CALL WRMSGS( INDX )
              ENDIF
 282        CONTINUE
 284      CONTINUE
!
!---      Loop over imbibition saturation function tables  ---
!
          DO 294 NTX = 1,NITABLX
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
            VARB = 'Number of Imbibition Table Entries'
            CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
            WRITE(IWR,'(2A,I6)') VARB(1:IVR),': ',NLIN
!
!---        Loop over lines in saturation function tables  ---
!
            NTBLX = NTBL+1
            DO 290 NL = 1,NLIN
              NTBL = NTBL + 1
              IF( NTBL.GT.LTBL ) THEN
                INDX = 5
                CHMSG = 'Number of Table Entries > Parameter LTBL'
                CALL WRMSGS( INDX )
              ENDIF
              CALL RDINPL( CHDUM )
              CALL LCASE( CHDUM )
              ISTART = 1
              VARB = 'Imbibition Capillary Head'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLX(NTBL))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',TBLX(NTBL)
              INDX = 0
              IUNM = 1
              CALL RDUNIT(UNTS,TBLX(NTBL),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',TBLX(NTBL),', m)'
!
!---          Log-linear interpolation  ---
!
              IF( ISCHRX.EQ.13 ) THEN
                TBLX(NTBL) = LOG( MAX( TBLX(NTBL),1.D-14 ) )
              ENDIF
              VARB = 'Imbibition Aqueous Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',TBLY(NTBL)
              IF( NL.EQ.2 ) THEN
                IF( TBLX(NTBL-1).LT.TBLX(NTBL) ) THEN
                  ITDX = 1
                ELSEIF( TBLX(NTBL-1).GT.TBLX(NTBL) ) THEN
                  ITDX = -1
                ELSE
                  INDX = 4
                  CHMSG = 'Invalid Imbibition Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
                IF( TBLY(NTBL-1).LT.TBLY(NTBL) ) THEN
                  ITDY = 1
                ELSEIF( TBLY(NTBL-1).GT.TBLY(NTBL) ) THEN
                  ITDY = -1
                ELSE
                  INDX = 4
                  CHMSG = 'Invalid Imbibition Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
              ELSEIF( NL.GT.2 ) THEN
                IF( (ITDX.EQ.1 .AND. TBLX(NTBL).LE.TBLX(NTBL-1)) .OR.
     &            (ITDX.EQ.-1 .AND. TBLX(NTBL).GE.TBLX(NTBL-1)) ) THEN
                  INDX = 4
                  CHMSG = 'Invalid Imbibition Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
                IF( (ITDY.EQ.1 .AND. TBLY(NTBL).LE.TBLY(NTBL-1)) .OR.
     &            (ITDY.EQ.-1 .AND. TBLY(NTBL).GE.TBLY(NTBL-1)) ) THEN
                  INDX = 4
                  CHMSG = 'Invalid Imbibition Saturation Table'
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
  290       CONTINUE
!
!---        Correlate table numbers with nodes  ---
!
            DO 292 N = 1,NFLD
              IF( JVAR(N).EQ.NTX ) THEN
                ISLTBL(3,N) = NTBLX
                ISLTBL(4,N) = NTBL
              ELSEIF( IVAR(N).LT.1 .OR. IVAR(N).GT.NITABLX ) THEN
                INDX = 4
                CHMSG = 'Invalid Imbibition Saturation Table Number'
                CALL WRMSGS( INDX )
              ENDIF
  292       CONTINUE
  294     CONTINUE
!
!---    Rock/soil zonation  ---
!
        ELSE
!
!---      Drainage table  ---
!
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
          VARB = 'Number of Drainage Table Entries'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
          WRITE(IWR,'(2A,I6)') VARB(1:IVR),': ',NLIN
          IF( NLIN.LT.2 ) THEN
            INDX = 4
            CHMSG = 'Saturation Invalid Drainage Table'
            CALL WRMSGS( INDX )
          ENDIF
          ISLTBL(1,IROCK) = NTBL + 1
          DO 300 NL = 1,NLIN
            NTBL = NTBL + 1
            IF( NTBL.GT.LTBL ) THEN
              INDX = 5
              CHMSG = 'Number of Table Entries > Parameter LTBL'
              CALL WRMSGS( INDX )
            ENDIF
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
            VARB = 'Drainage Capillary Head'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLX(NTBL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',TBLX(NTBL)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,TBLX(NTBL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',TBLX(NTBL),', m)'
!
!---        Log-linear interpolation  ---
!
            IF( ISCHRX.EQ.13 ) THEN
              TBLX(NTBL) = LOG( MAX( TBLX(NTBL),1.D-14 ) )
            ENDIF
            VARB = 'Drainage Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',TBLY(NTBL)
            IF( NL.EQ.2 ) THEN
              IF( TBLX(NTBL-1).LT.TBLX(NTBL) ) THEN
                ITDX = 1
              ELSEIF( TBLX(NTBL-1).GT.TBLX(NTBL) ) THEN
                ITDX = -1
              ELSE
                INDX = 4
                CHMSG = 'Invalid Saturation Drainage Table'
                CALL WRMSGS( INDX )
              ENDIF
              IF( TBLY(NTBL-1).LT.TBLY(NTBL) ) THEN
                ITDY = 1
              ELSEIF( TBLY(NTBL-1).GT.TBLY(NTBL) ) THEN
                ITDY = -1
              ELSE
                INDX = 4
                CHMSG = 'Invalid Saturation Drainage Table'
                CALL WRMSGS( INDX )
              ENDIF
            ELSEIF( NL.GT.2 ) THEN
              IF( (ITDX.EQ.1 .AND. TBLX(NTBL).LE.TBLX(NTBL-1)) .OR.
     &          (ITDX.EQ.-1 .AND. TBLX(NTBL).GE.TBLX(NTBL-1)) ) THEN
                INDX = 4
                CHMSG = 'Invalid Saturation Drainage Table'
                CALL WRMSGS( INDX )
              ENDIF
              IF( (ITDY.EQ.1 .AND. TBLY(NTBL).LE.TBLY(NTBL-1)) .OR.
     &          (ITDY.EQ.-1 .AND. TBLY(NTBL).GE.TBLY(NTBL-1)) ) THEN
                INDX = 4
                CHMSG = 'Invalid Saturation Drainage Table'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
  300     CONTINUE
          ISLTBL(2,IROCK) = NTBL
!
!---      Imbibition table  ---
!
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
          VARB = 'Number of Imbibition Table Entries'
          CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
          WRITE(IWR,'(2A,I6)') VARB(1:IVR),': ',NLIN
          IF( NLIN.LT.2 ) THEN
            INDX = 4
            CHMSG = 'Saturation Invalid Imbibition Table'
            CALL WRMSGS( INDX )
          ENDIF
          ISLTBL(3,IROCK) = NTBL + 1
          DO 310 NL = 1,NLIN
            NTBL = NTBL + 1
            IF( NTBL.GT.LTBL ) THEN
              INDX = 5
              CHMSG = 'Number of Table Entries > Parameter LTBL'
              CALL WRMSGS( INDX )
            ENDIF
            CALL RDINPL( CHDUM )
            CALL LCASE( CHDUM )
            ISTART = 1
            VARB = 'Imbibition Capillary Head'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLX(NTBL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',TBLX(NTBL)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,TBLX(NTBL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',TBLX(NTBL),', m)'
!
!---        Log-linear interpolation  ---
!
            IF( ISCHRX.EQ.13 ) THEN
              TBLX(NTBL) = LOG( MAX( TBLX(NTBL),1.D-14 ) )
            ENDIF
            VARB = 'Imbibition Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,TBLY(NTBL))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',TBLY(NTBL)
            IF( NL.EQ.2 ) THEN
              IF( TBLX(NTBL-1).LT.TBLX(NTBL) ) THEN
                ITDX = 1
              ELSEIF( TBLX(NTBL-1).GT.TBLX(NTBL) ) THEN
                ITDX = -1
              ELSE
                INDX = 4
                CHMSG = 'Invalid Saturation Imbibition Table'
                CALL WRMSGS( INDX )
              ENDIF
              IF( TBLY(NTBL-1).LT.TBLY(NTBL) ) THEN
                ITDY = 1
              ELSEIF( TBLY(NTBL-1).GT.TBLY(NTBL) ) THEN
                ITDY = -1
              ELSE
                INDX = 4
                CHMSG = 'Invalid Saturation Imbibition Table'
                CALL WRMSGS( INDX )
              ENDIF
            ELSEIF( NL.GT.2 ) THEN
              IF( (ITDX.EQ.1 .AND. TBLX(NTBL).LE.TBLX(NTBL-1)) .OR.
     &          (ITDX.EQ.-1 .AND. TBLX(NTBL).GE.TBLX(NTBL-1)) ) THEN
                INDX = 4
                CHMSG = 'Invalid Saturation Imbibition Table'
                CALL WRMSGS( INDX )
              ENDIF
              IF( (ITDY.EQ.1 .AND. TBLY(NTBL).LE.TBLY(NTBL-1)) .OR.
     &          (ITDY.EQ.-1 .AND. TBLY(NTBL).GE.TBLY(NTBL-1)) ) THEN
                INDX = 4
                CHMSG = 'Invalid Saturation Imbibition Table'
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
  310     CONTINUE
          ISLTBL(4,IROCK) = NTBL
        ENDIF
      ENDIF
!
!---  Translate saturation function type for IJK indexing  ---
!
      IF( IJK.GT.0 ) THEN
        DO 462 N = 1,NFLD
          ISCHR(IZ(N)) = ISCHRX
          IF( ISCHRX.EQ.7 ) IRPL(IZ(N)) = 7
          ISM(IZ(N)) = ISMX
  462   CONTINUE
!
!---  For IJK indexing with the rock/soil option, correlate rock/soil 
!     numbers with nodes for saturation function type and
!     parameters  ---
!
      ELSEIF( IJK.LT.0 ) THEN
        DO 472 N = 1,NFLD
          IF( IZ2(N).EQ.NR ) THEN
            DO 470 L = 1,LSCHR
              SCHR(L,N) = SCHRX(L)
  470       CONTINUE
            ISCHR(N) = ISCHRX
            IF( ISCHRX.EQ.7 ) IRPL(N) = 7
            ISM(N) = ISMX
          ENDIF
  472   CONTINUE
      ELSE
        ISCHR(IROCK) = ISCHRX
        IF( ISCHRX.EQ.7 ) IRPL(IROCK) = 7
        ISM(IROCK) = ISMX
        DO 480 L = 1,LSCHR
          SCHR(L,IROCK) = SCHRX(L)
  480   CONTINUE
      ENDIF
!
!---  Loop over remaining rock/soils within scaling group  ---
!
      IF( ISLC(19).EQ.1 .AND. IROCK.LT.NROCK ) THEN
        DO 494 M = IROCK+1,NROCK
          IF( ISCALE(M).EQ.ISGRP ) THEN
            NR = NR + 1
            ISM(M) = ISM(IROCK)
            ISCHR(M) = ISCHR(IROCK)
            IF( ISCHR(M).EQ.7 ) IRPL(M) = 7
            DO 490 L = 1,LSCHR
              SCHR(L,M) = SCHR(L,IROCK)
  490       CONTINUE
            DO 492 L = 1,4
              ISLTBL(L,M) = ISLTBL(L,IROCK)
  492       CONTINUE
          ENDIF
  494   CONTINUE
      ENDIF
!
!---  Read next rock/soil type or scaling group  ---
!
      IF( NR.LT.NROCK ) WRITE(IWR,'(/)')
!
!---  Continue reading rock/soil type names for a pattern match  ---
!
      IF( KBS.EQ.1 .AND. IROCK.LT.NROCK ) THEN
        IROCK = IROCK + 1
        ISTART = ISBS
        GOTO 20
      ENDIF
      GOTO 10
 500  CONTINUE
      IF( ALLOCATED(IVAR) ) THEN
      DEALLOCATE( IVAR,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: IVAR'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      IF( ALLOCATED(JVAR) ) THEN      
      DEALLOCATE( JVAR,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: JVAR'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSP_GT group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSR_GT
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
!     Read input file for source information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE REACT
      USE GRID
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
      SUB_LOG(ISUB_LOG) = '/RDSR_GT'
!
!---  Dynamic memory allocation  ---
!
      ALLOCATE( VAR(1:LSTM,1:8),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: VAR'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write card information to ouput file  ---
!
      CARD = 'Source Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      NSR = 0
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Sources'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO 140 NS = 1, NLIN
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
!
!---  Read source type  ---
!
        VARB = 'Source Type'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'power density').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Power Density Source'
          ISRTX = 2
        ELSEIF( INDEX(ADUM(1:),'power').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Power Source'
          ISRTX = 1
        ELSEIF( INDEX(ADUM(1:),'portland').NE.0 .AND. 
     &    INDEX(ADUM(1:),'cement').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Setting of Portland Cement'
          ISRTX = 11
        ELSEIF( INDEX(ADUM(1:),'aqueous volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Volumetric Source'
          ISRTX = 3
        ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .AND. 
     &    INDEX(ADUM(1:),'sink').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Mass Sink'
          ISRTX = 10
        ELSEIF( INDEX(ADUM(1:),'gas volu').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'relative').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Volumetric Source w/ Rel. Humidity'
            ISRTX = 4
          ELSEIF( INDEX(ADUM(1:),'mass frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Volumetric Source w/ Mass Frac.'
            ISRTX = 5
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'condensate' ).NE. 0 ) THEN
          WRITE(IWR,'(2X,A)') 'Undocumented Condensate Source'
          ISRTX = 6
        ELSEIF( INDEX(ADUM(1:),'aqueous mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Mass Source'
          ISRTX = 7
        ELSEIF( INDEX(ADUM(1:),'relative').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'mass frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Mass Source w/ Rel. Humidity'
            ISRTX = 8
          ELSEIF( INDEX(ADUM(1:),'gas mass').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Mass Source w/ Mass Frac.'
            ISRTX = 9
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'z-dir').NE.0 .AND.
     &    INDEX(ADUM(1:),'injec').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Z-Direction Vertical ' //
     &      'Injection Well Source'
          VARB = 'Water-Vapor Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
            ISRTX = 13
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
            ISRTX = 14
          ELSE
            ISRTX = 15
            WRITE(IWR,'(2X,A)') 'No Water-Vapor'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'x-dir').NE.0 .AND.
     &    INDEX(ADUM(1:),'injec').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'X-Direction Horizontal ' //
     &      'Injection Well Source'
          VARB = 'Water-Vapor Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
            ISRTX = 16
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
            ISRTX = 17
          ELSE
            ISRTX = 18
            WRITE(IWR,'(2X,A)') 'No Water-Vapor'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'y-dir').NE.0 .AND.
     &    INDEX(ADUM(1:),'injec').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Y-Direction Horizontal ' //
     &      'Injection Well Source'
          VARB = 'Water-Vapor Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
            ISRTX = 19
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
            ISRTX = 20
          ELSE
            ISRTX = 21
            WRITE(IWR,'(2X,A)') 'No Water-Vapor'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volum').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Gas Volumetric' //
     &      'Injection/Extraction Well Source'
          VARB = 'Water-Vapor Source Option'
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          IF( INDEX(BDUM(1:),'rel').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Relative Humidity'
            ISRTX = 23
          ELSEIF( INDEX(BDUM(1:),'mass').NE.0 .AND.
     &      INDEX(BDUM(1:),'frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Water-Vapor Gas Mass Fraction'
            ISRTX = 24
          ELSE
            ISRTX = 25
            WRITE(IWR,'(2X,A)') 'No Water-Vapor'
          ENDIF
        ELSEIF( IEQC.NE.0 .AND. INDEX(ADUM(1:),'solute').NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          DO 30 NSL = 1,NSOLU
            IDB = INDEX(SOLUT(NSL)(1:),'  ')
            IF( INDEX(BDUM(1:),SOLUT(NSL)(1:IDB)).NE.0 ) THEN
              IF( INDEX(ADUM(1:),'density').NE.0 ) THEN
                ISRTX = -(NSL+NSOLU)
                WRITE(IWR,'(2X,2A)')'Solute Source Density: ',SOLUT(NSL)
              ELSEIF( INDEX(ADUM(1:),'diffusion').NE.0 ) THEN
                IF( INDEX(ADUM(1:),'variable').NE.0 ) THEN
                  ISRTX = -(NSL+8*NSOLU)
                  WRITE(IWR,'(2X,2A)')'Variable Diffusion-Dominated ' //
     &              'Solute Release Model: ',SOLUT(NSL)
                ELSE
                  ISRTX = -(NSL+5*NSOLU)
                  WRITE(IWR,'(2X,2A)')'Diffusion-Dominated Solute ' //
     &              'Release Model: ',SOLUT(NSL)
                ENDIF
              ELSEIF( INDEX(ADUM(1:),'solubility').NE.0 ) THEN
                IF( INDEX(ADUM(1:),'salt').NE.0 .OR.
     &            INDEX(ADUM(1:),'cake').NE.0 ) THEN
                  ISRTX = -(NSL+7*NSOLU)
                  WRITE(IWR,'(2X,2A)')'Solubility-Controlled Salt ' //
     &              'Cake Release Model: ',SOLUT(NSL)
                ELSE
                  ISRTX = -(NSL+6*NSOLU)
                  WRITE(IWR,'(2X,2A)')'Solubility-Controlled Solute ' //
     &              'Release Model: ',SOLUT(NSL)
                ENDIF
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

        ELSEIF( IEQC.NE.0 .AND. INDEX(ADUM(1:),'specie').NE.0 ) THEN
          VARB = 'Species Name: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
!
!---      Aqueous species  ---
!          
          DO 42 NSP = 1,NSPL
            IDB = INDEX(SPNML(NSP)(1:),'  ')
            IF( BDUM(1:IDB).EQ.SPNML(NSP)(1:IDB) ) THEN
              IF( INDEX(ADUM(1:),'density').NE.0 ) THEN
                ISRTX = 100+NSPL+NSPS+NSP
                WRITE(IWR,'(2X,2A)')'Solute Source Density: ',
     &            SPNML(NSP)(1:IDB)
              ELSE
                ISRTX = 100+NSP
                WRITE(IWR,'(2X,2A)')'Species Source: ',
     &            SPNML(NSP)(1:IDB)
              ENDIF
              GOTO 44
            ENDIF
   42     CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Source Species Name: '//BDUM
            CALL WRMSGS( INDX )
   44     CONTINUE

        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Source Type: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Read source domain indices  ---
!
        VARB = 'Source Domain Index'
        ISX = ISTART
        CALL RDINT(ISTART,ICOMMA,CHDUM,I1X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,I2X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,J1X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,J2X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,K1X)
        CALL RDINT(ISTART,ICOMMA,CHDUM,K2X)
        ICX = ISTART
        WRITE(IWR,'(/,2X,A)') 'Source Domain:'
        WRITE(IWR,'(4X,A,I6,A,I6)') 'I = ',I1X,' to ',I2X
        WRITE(IWR,'(4X,A,I6,A,I6)') 'J = ',J1X,' to ',J2X
        WRITE(IWR,'(4X,A,I6,A,I6)') 'K = ',K1X,' to ',K2X
!
!---  Check for ill-defined source domains  ---
!
        IF( I1X.LT.1 .OR. I1X.GT.IFLD .OR. I2X.LT.1 .OR.
     &    I2X.GT.IFLD .OR. I2X.LT.I1X ) THEN
          INDX = 4
          CHMSG = 'Invalid Source Domain: ' // CHDUM(ISX:ICX)
          CALL WRMSGS( INDX )
        ENDIF
        IF( J1X.LT.1 .OR. J1X.GT.JFLD .OR. J2X.LT.1 .OR.
     &    J2X.GT.JFLD .OR. J2X.LT.J1X ) THEN
          INDX = 4
          CHMSG = 'Invalid Source Domain: ' // CHDUM(ISX:ICX)
          CALL WRMSGS( INDX )
        ENDIF
        IF( K1X.LT.1 .OR. K1X.GT.KFLD .OR. K2X.LT.1 .OR.
     &    K2X.GT.KFLD .OR. K2X.LT.K1X ) THEN
          INDX = 4
          CHMSG = 'Invalid Source Domain: ' // CHDUM(ISX:ICX)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Check for sources applied to inactive nodes  ---
!
        DO K = K1X,K2X
          DO J = J1X,J2X
            DO I = I1X,I2X
              IF( IXP(ND(I,J,K)).EQ.0 ) THEN
                INDX = 4
                CHMSG = 'Source Applied to an Inactive Node'
                CALL WRMSGS( INDX )
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!
!---  Check for z-direction well sources  ---
!
        IF( (ISRTX.GE.13.AND.ISRTX.LE.15) .OR. 
     &       (ISRTX.GE.23.AND.ISRTX.LE.25)  ) THEN
          IF( ((I2X-I1X).GT.1) .OR. ((J2X-J1X).GT.1) ) THEN
            INDX = 4
            CHMSG = 'Invalid Well Source Domain'
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Read number of source times  ---
!
        IF( NS.GT.LSR ) THEN
          INDX = 5
          CHMSG = 'Number of Sources > Parameter LSR'
          CALL WRMSGS( INDX )
        ENDIF
        VARB = 'Number of Source Times'
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISRM(NS))
        IF( ISRM(NS).GT.LSTM ) THEN
          INDX = 5
          CHMSG = 'Number of Source Times > Parameter LSTM'
          CALL WRMSGS( INDX )
        ENDIF
        SRTMO = -SMALL
        DO 100 NTM = 1,ISRM(NS)
          DO 60 M = 1,6
            VAR(NTM,M) = 0.D+0
   60     CONTINUE
!
!---  Read and write source values and units  ---
!
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Source Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',VAR(NTM,1)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR(NTM,1),INDX)
          IF( ISRTX.EQ.1 ) THEN
            VARB = 'Source Power'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 2
            IUNKG = 1
            IUNS = -3
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ELSEIF( ISRTX.EQ.2 ) THEN
            VARB = 'Source Power Density'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -3
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ELSEIF( ISRTX.EQ.11 ) THEN
            VARB = 'Source Total Heat Evolution, J/kg'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 2
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'XHALF Parameter yO'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,2)
            VARB = 'XHALF Parameter A'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,3)
            VARB = 'XHALF Parameter invTau'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'RATE Parameter yO'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
            VARB = 'RATE Parameter A'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,7)
            VARB = 'RATE Parameter invTau'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,8)
          ELSEIF( ISRTX.EQ.3 ) THEN
            VARB = 'Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Source Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Aqueous Volumetric Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Dissolved-Air Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VARB = 'Source Dissolved-Salt Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,6)
          ELSEIF( ISRTX.EQ.4 ) THEN
            VARB = 'Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Source Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Gas Volumetric Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Water Vapor Relative Humidity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
          ELSEIF( ISRTX.EQ.5 ) THEN
            VARB = 'Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Source Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Gas Volumetric Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Water Vapor Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
          ELSEIF( ISRTX.EQ.6 ) THEN
            VARB = 'Source Condensate Domain Index'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            IVAR = INT(VAR(NTM,4))
            WRITE(IWR,'(/,2A,I4)') VARB(1:IVR),': ',IVAR
            VARB = 'Source Condensate Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Source System Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Percent Saturated Dissolved Air in Condensate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            VAR(NTM,6) = 0.D+0
            DO K = K1X,K2X
              DO J = J1X,J2X
                DO I = I1X,I2X
                  VAR(NTM,6) = VAR(NTM,6) + VOL(ND(I,J,K))
                ENDDO
              ENDDO
            ENDDO
            VARB = 'Source Condensate Fractional Mass Loss Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,7)
          ELSEIF( ISRTX.EQ.7 ) THEN
            VARB = 'Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Source Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Aqueous Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Dissolved-Air Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
          ELSEIF( ISRTX.EQ.8 ) THEN
             VARB = 'Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Source Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Gas Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Water Vapor Relative Humidity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
          ELSEIF( ISRTX.EQ.9 ) THEN
             VARB = 'Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Source Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Gas Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Water Vapor Mass Fraction'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
!
!---      Mass sink  ---
!
          ELSEIF( ISRTX.EQ.10 ) THEN
            VARB = 'Mass Sink Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            IF( VAR(NTM,4).LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Negative Mass Sink'
              RLMSG = VAR(NTM,4)
              CALL WRMSGS( INDX )
            ENDIF
            VAR(NTM,4) = -VAR(NTM,4)
!
!---      Injection Well Source  ---
!
          ELSEIF( ISRTX.GE.13 .AND. ISRTX.LE.21 ) THEN
            VARB = 'Well Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,2),', Pa)'
            VARB = 'Well Diameter'
            IDFLT = 1
            VAR(NTM,3) = 1.7D-1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,3),', m)'
            VARB = 'Symmetry Factor'
            IDFLT = 1
            VAR(NTM,4) = 1.D+0
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,4)
!
!---        Convert well pressure to guage and well diameter
!           to well radius  ---
!
            VAR(NTM,2) = VAR(NTM,2)-PATM
            VAR(NTM,3) = 5.D-1*VAR(NTM,3)
!
!---        Well water vapor  ---
!
            IF( ISRTX.EQ.13.OR.ISRTX.EQ.16.OR.ISRTX.EQ.19 ) THEN
              VARB = 'Source Water Vapor Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( ISRTX.EQ.14.OR.ISRTX.EQ.17.OR.ISRTX.EQ.20 ) THEN
              VARB = 'Source Water Vapor Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( ISRTX.EQ.15.OR.ISRTX.EQ.18.OR.ISRTX.EQ.21 ) THEN
              ISRTX = ISRTX-1
              VAR(NTM,5) = 0.D+0
            ENDIF
!
!---        Minimum permeability for injection/withdrawal  ---
!
            VARB = 'Minimum Permeability'
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            IHC = INDEX( UNTS(1:),'hc' )
            IF( IHC .EQ. 0 ) THEN
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,6)
              IUNM = 2
            ELSE
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',VAR(NTM,6)
              IUNM = 1
              IUNS = -1
            ENDIF
            INDX = 0
            CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,6),', m^2)'
            VARB = 'Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,7)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,7),INDX)
!
!---      Gas Volumetric Well Source  ---
!
          ELSEIF( ISRTX.GE.23 .AND. ISRTX.LE.25 ) THEN
            VARB = 'Source Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            VARB = 'Minimum Extraction Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            VARB = 'Source Gas Volumetric Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
!
!---        Well water vapor  ---
!
            IF( ISRTX.EQ.23 ) THEN
              VARB = 'Source Water Vapor Relative Humidity'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( ISRTX.EQ.24 ) THEN
              VARB = 'Source Water Vapor Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
              WRITE(IWR,'(/,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,5)
            ELSEIF( ISRTX.EQ.25 ) THEN
              ISRTX = ISRTX-1
              VAR(NTM,5) = 0.D+0
            ENDIF
            VARB = 'Well Diameter'
            IDFLT = 1
            VAR(NTM,6) = 1.7D-1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,6)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,6),', m)'
            VARB = 'Symmetry Factor'
            IDFLT = 1
            VAR(NTM,7) = 1.D+0
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',VAR(NTM,7)
!
!---        Convert well pressure to guage and well diameter
!           to well radius  ---
!
            VAR(NTM,3) = VAR(NTM,3)-PATM
            VAR(NTM,6) = 5.D-1*VAR(NTM,6)

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
!
!---      Diffusion-dominated solute release model  ---
!
          ELSEIF( ISRTX.LT.-5*NSOLU .AND. ISRTX.GE.-6*NSOLU ) THEN
            VARB = 'Nodal Solute Inventory'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,3)
            VARB = 'Vertical Depth of Residual Waste'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,4),', m)'
            VARB = 'Diffusion Coefficient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,5)
            INDX = 0
            IUNM = 2
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,5),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,5),', m^2/s)'
!
!---      Solubility-controlled solute release model  ---
!
          ELSEIF( ISRTX.LT.-6*NSOLU .AND. ISRTX.GE.-7*NSOLU ) THEN
!
!---        SRX(2): nodal solute inventory
!           SRX(3): aqueous solubility  ---
!
            VARB = 'Nodal Solute Inventory'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,2)
            VARB = 'Aqueous Solubility of Solute'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,3)
            INDX = 0
            IUNM = -3
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,3),', 1/m^3)'
!
!---      Solubility-controlled salt cake release model  ---
!
          ELSEIF( ISRTX.LT.-7*NSOLU .AND. ISRTX.GE.-8*NSOLU ) THEN
!
!---        SRX(2): nodal solute inventory
!           SRX(3): nodal salt cake inventory
!           SRX(4): salt cake solubility  ---
!
            VARB = 'Nodal Solute Inventory'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,2)
            VARB = 'Nodal Salt Cake Inventory'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,3)
            VARB = 'Aqueous Solubility of Salt Cake'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = -3
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,4),', 1/m^3)'
!
!---      Diffusion-dominated solute release model (w/variable diffusion) ---
!
          ELSEIF( ISRTX.LT.-8*NSOLU .AND. ISRTX.GE.-9*NSOLU ) THEN
            VARB = 'Nodal Solute Inventory'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,3)
            VARB = 'Vertical Depth of Residual Waste'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,4),', m)'
            VARB = 'Diffusion Coefficient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,5)
            INDX = 0
            IUNM = 2
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,5),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,5),', m^2/s)'
            VARB = 'Constrictivity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2A,1PE11.4,$)') VARB(1:IVR),': ',VAR(NTM,6)
          ENDIF
!
!---  Check for nonascending source times  ---
!
          IF( VAR(NTM,1).LT.SRTMO ) THEN
            INDX = 4
            CHMSG = 'Invalid Source Time Sequencing'
            CALL WRMSGS( INDX )
          ENDIF
          SRTMO = VAR(NTM,1)
  100   CONTINUE
!
!---  Assign values to source variables  ---
!
        NSR = NSR + 1
        IF( NSR.GT.LSR ) THEN
          INDX = 5
          CHMSG = 'Number of Sources > Parameter LSR'
          CALL WRMSGS( INDX )
        ENDIF
        ISRDM(1,NSR) = I1X
        ISRDM(2,NSR) = I2X
        ISRDM(3,NSR) = J1X
        ISRDM(4,NSR) = J2X
        ISRDM(5,NSR) = K1X
        ISRDM(6,NSR) = K2X
        ISRT(NSR) = ISRTX
        DO 130 NTM = 1,ISRM(NS)
          DO 120 M = 1,8
            SRC(M,NTM,NSR) = VAR(NTM,M)
  120     CONTINUE
  130   CONTINUE
  140 CONTINUE
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
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSR_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDST_GT
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
!     Reads the salt transport card dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 27 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE GRID
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
      CHARACTER*64 RDUM,UNTS
      CHARACTER*512 CHDUM
      INTEGER NCH
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDST_GT'
!
!---  Write card information to output file  ---
!
      CARD = 'Salt Transport Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      IDSPS = 0
!
!---  Salt name  ---
!
      WRITE(IWR,'(A)') 'Sodium Chloride Salt'
      WTMS = 58.4428
      WRITE(IWR,'(2X,A,1PE11.4)')  'Molecular Weight, kg/kg_mol: ',WTMS
      IEDLS = 1
      IF( ISLC(4).EQ.1 ) IEDLS = 3
      WRITE( IWR,'(A)' ) 'Conventional Diffusion Model'
!
!---  Loop over the rock/soil saturation information lines  ---
!
      N = 0
      IJK = 0
   10 CONTINUE
        IF( N.GE.NROCK .OR. IJK.GT.0 ) GOTO 600
        ISTART = 1
        VARB = 'Rock Name: '
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,RDUM)
!
!---  IJK, KIJ, or JKI indexing  ---
!
        IF( INDEX(RDUM(1:),'indexing').NE.0 ) THEN
          IF( INDEX(ROCK(1)(1:),'indexing').EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Indexing Option Not Declared ' // 
     &        'in Rock/Soil Zonation Card'
            CALL WRMSGS( INDX )
          ENDIF
          IF( INDEX(RDUM,'ijk').NE.0 ) THEN
            IJK = 1
          ELSEIF( INDEX(RDUM,'jki').NE.0 ) THEN
            IJK = 2
          ELSEIF( INDEX(RDUM,'kij').NE.0 ) THEN
            IJK = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Indexing Option' // RDUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
          GOTO 220
        ENDIF
!
!---  Search known rock types for a matching type ---
!
        DO 100 M = 1, NROCK
          IF( RDUM.EQ.ROCK(M)) THEN
            IROCK = M
            GOTO 200
          ENDIF
  100   CONTINUE
!
!---  Search known scaling groups for a matching type ---
!
        IF( ISLC(19).EQ.1 ) THEN
          DO 110 M = 1,NSCALE
             IF( RDUM.EQ.SCALNM(M) ) THEN
                ISGRP = M
                IROCK = 1
                GOTO 200
             ENDIF
  110     CONTINUE
          INDX = 2
          CHMSG = 'Unrecognized Rock/Soil Type or Scaling Group: '
     &      // RDUM(1:NCH)
          CALL WRMSGS( INDX )
          GOTO 10
        ENDIF
        INDX = 2
        CHMSG = 'Unrecognized Rock/Soil Type: ' // RDUM(1:NCH)
        CALL WRMSGS( INDX )
        GOTO 10
  200   CONTINUE
!
!---  Loop over rock/soils within scaling group  ---
!
        IF( ISLC(19).EQ.1 .AND. ISGRP.NE.0 ) THEN
          DO 202 M = IROCK,NROCK
            IF( ISCALE(M).EQ.ISGRP ) THEN
              IROCK = M
              GOTO 204
            ENDIF
  202     CONTINUE
        ENDIF
  204   CONTINUE
!
!---    Write rock/soil name  ---
!
        WRITE (IWR,'(/,2A)') 'Rock/Soil Name: ',ROCK(IROCK)
        N = N + 1
  220   CONTINUE
!
!---  Longitudinal dispersivity  ---
!
        VARB = 'Longitudinal Dispersivity: '
        IF( IJK.GT.0 ) THEN
          UNTS = 'm'
          IUNM = 1
          CALL RDIJK( ISTART,IJK,CHDUM,UNTS,DPLGS )
          IDSPS = 1
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DPLGS(IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &      DPLGS(IROCK)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,DPLGS(IROCK),INDX)
          IF( DPLGS(IROCK).GE.SMALL ) IDSPS = 1
        ENDIF
!
!---  Transverse dispersivity  ---
!
        VARB = 'Transverse Dispersivity: '
        IF( IJK.GT.0 ) THEN
          UNTS = 'm'
          IUNM = 1
          CALL RDIJK( ISTART,IJK,CHDUM,UNTS,DPTRS )
          IDSPS = 1
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DPTRS(IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &      DPTRS(IROCK)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,DPTRS(IROCK),INDX)
          IF( DPTRS(IROCK).GE.SMALL ) IDSPS = 1
        ENDIF
!
!---    Read next rock/soil type or scaling group  ---
!
        IF( N.LT.NROCK ) WRITE(IWR,'(/)')
        GOTO 10
  600 CONTINUE
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDST_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTF_GT
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
!     Reads solute/fluid interaction card for diffusion and partition
!     coefficients, and internodal diffusion term averaging scheme for
!     single phase (aqueous) solute transport equation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE FILES
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
      CHARACTER*64 ADUM,BDUM,UNTS
      CHARACTER*512 CHDUM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICDSX,JCDSX
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISCX,IPSPX,IBCDPX
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICLX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IADJM,JADJM,KADJM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDTF_GT'
!
!---  Write card information to output file  ---
!
      CARD = 'Solute/Fluid Interaction Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of different solutes  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Solutes'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      NSOLU = 0
      DO 200 NL = 1, NLIN
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        ADUM(1:) = ' '
        VARB = 'Solute Name'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        DO 100 NSL = 1,NSOLU
          IF( SOLUT(NSL).EQ.ADUM ) GOTO 110
  100   CONTINUE
        NSOLU = NSOLU + 1
        IF( NSOLU.GT.LSOLU ) THEN
          INDX = 5
          CHMSG = 'Number of Solutes > Parameter LSOLU'
          CALL WRMSGS( INDX )
        ENDIF
        SOLUT(NSOLU) = ADUM
        NSL = NSOLU
  110   CONTINUE
        WRITE(IWR,'(/,3A)') VARB(1:IVR),': ',ADUM
!
!---    Aqueous effective diffusion option, first check for 
!       real input, indicating aqueous-phase molecular diffusion
!       coefficient was entered  ---
!
        VARB = 'Aqueous Effective Diffusion Option: '
        ISX = ISTART
        ICX = ICOMMA
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX(ADUM(1:),'empirical').NE.0 .OR.
     &    INDEX(ADUM(1:),'conventional').NE.0 .OR.
     &    INDEX(ADUM(1:),'constant').NE.0 .OR.
     &    INDEX(ADUM(1:),'immobile').NE.0 .OR.
     &    INDEX(ADUM(1:),'stationary').NE.0 ) THEN
          WRITE( IWR,'(/,A,$)' ) VARB(1:IVR)
          IF( INDEX(ADUM(1:),'empirical').NE.0 )  THEN
            IEDL(NSL) = 2
            WRITE( IWR,'(A)' ) 'Kemper and van Schaik ' // 
     &        'Empirical Diffusion Model'
            WRITE( IWR,'(A)' ) '  Model Parameters Entered ' // 
     &        'on the Solute/Porous Media Interaction Card'
          ELSEIF( INDEX(ADUM(1:),'conventional').NE.0 )  THEN
            IEDL(NSL) = 1
            WRITE( IWR,'(A)' ) 'Conventional Diffusion Model'
            VARB = 'Aqueous Molecular Diffusion Coefficient @ 20 C'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDL(NSL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',SMDL(NSL)
            INDX = 0
            IUNM = 2
            IUNS = -1
            CALL RDUNIT(UNTS,SMDL(NSL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SMDL(NSL),', m^2/s)'
          ELSEIF( INDEX(ADUM(1:),'constant').NE.0 )  THEN
            IEDL(NSL) = 3
            WRITE( IWR,'(A)' ) 'Constant Diffusion Model'
            VARB = 'Aqueous Molecular Diffusion Coefficient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDL(NSL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',SMDL(NSL)
            INDX = 0
            IUNM = 2
            IUNS = -1
            CALL RDUNIT(UNTS,SMDL(NSL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SMDL(NSL),', m^2/s)'
          ELSEIF( INDEX(ADUM(1:),'immobile').NE.0 .OR.
     &      INDEX(ADUM(1:),'stationary').NE.0 )  THEN
            IEDL(NSL) = 4
            WRITE( IWR,'(A)' ) 'Immobile/Stationary Solute'
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Aqueous Diffusion Option: '//ADUM
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Aqueous-phase molecular diffusion coefficient  ---
!
        ELSE
          ISTART = ISX
          ICOMMA = ICX
          IEDL(NSL) = 1
          VARB = 'Aqueous Molecular Diffusion Coefficient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDL(NSL))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',SMDL(NSL)
          INDX = 0
          IUNM = 2
          IUNS = -1
          CALL RDUNIT(UNTS,SMDL(NSL),INDX)
        ENDIF
!
!---    Gas-phase molecular diffusion coefficient  ---
!
        IF( IEDL(NSL).NE.4 ) THEN 
          VARB = 'Gas-Phase Molecular Diffusion Coefficient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDG(NSL))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',SMDG(NSL)
          INDX = 0
          IUNM = 2
          IUNS = -1
          CALL RDUNIT(UNTS,SMDG(NSL),INDX)
!
!---      Gas-aqueous partition coefficient option  ---
!
          VARB = 'Gas-Aqueous Partition Function: '
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
          WRITE( IWR,'(/,A,$)' ) VARB(1:IVR),': '
          IF( INDEX(ADUM(1:),'constant').NE.0 )  THEN
            IPCGL(NSL) = 0
            WRITE( IWR,'(A)' ) ': Constant'
            VARB = 'Gas-Aqueous Partition Coefficient'
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCGL(1,NSL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PCGL(1,NSL)
            INDX = 0
            CALL RDUNIT(UNTS,PCGL(1,NSL),INDX)
          ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 )  THEN
            IPCGL(NSL) = 1
            WRITE( IWR,'(A)' ) ': Temperature Dependent'
            WRITE( IWR,'(A)' ) 'ln( Kgl ) = a + b/T + c ln(T) ' // 
     &        '+ dT + eT^2'
            VARB = 'Gas-Aqueous Partition Function Coefficients: '
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCGL(1,NSL))
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCGL(2,NSL))
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCGL(3,NSL))
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCGL(4,NSL))
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCGL(5,NSL))
            WRITE(IWR,'(2X,A,1PE11.4)') 'Constant a: ',PCGL(1,NSL)
            WRITE(IWR,'(2X,A,1PE11.4)') 'Constant b: ',PCGL(2,NSL)
            WRITE(IWR,'(2X,A,1PE11.4)') 'Constant c: ',PCGL(3,NSL)
            WRITE(IWR,'(2X,A,1PE11.4)') 'Constant d: ',PCGL(4,NSL)
            WRITE(IWR,'(2X,A,1PE11.4)') 'Constant e: ',PCGL(5,NSL)
          ELSEIF( INDEX(ADUM(1:),'water').NE.0 .AND. 
     &      INDEX(ADUM(1:),'vapor').NE.0 )  THEN
            IPCGL(NSL) = 2
            WRITE( IWR,'(A)' ) ': Water Vapor Equilibrium'
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Gas-Aqueous Partition Option: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Solid-Aqueous Partition option  ---
!
          VARB = 'Solid-Aqueous Partition Option: '
          ADUM = 'continuous'
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
          WRITE( IWR,'(/,A,$)' ) VARB(1:IVR)
          IF( INDEX(ADUM(1:),'noncontinuous').NE.0 )  THEN
            IPCL(NSL) = 2
            WRITE( IWR,'(A)' ) 'Noncontinuous Solid Wetting'
          ELSE
            IPCL(NSL) = 1
            WRITE( IWR,'(A)' ) 'Continuous Solid Wetting'
          ENDIF
        ENDIF

!
!---    Half-life  ---
!
        IDFLT = 1
        VARB = 'Radioactive Half-Life'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,HLF(NSL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),': '
     &    ,HLF(NSL)
        INDX = 0
        IUNS = 1
        CALL RDUNIT(UNTS,HLF(NSL),INDX)
        HLF(NSL) = MAX( HLF(NSL),SMALL )
!
!---    Half-life modification  ---
!
        CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
        IF( INDX.EQ.1 ) THEN
          CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'schroth').NE.0 ) THEN
            WRITE( IWR,'(A)' ) 'Schroth Half-Life Modification Model'
            PCLN(3,NSL) = -1.D+0
            VARB = 'Minimum Saturation for Biological Activity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCLN(1,NSL))
            WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCLN(1,NSL)
            VARB = 'Half-Velocity Factor of Biological Activity'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCLN(2,NSL))
            WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCLN(2,NSL)
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Half-Life Modification Option: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  Read number of lines of chain decay information  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = ''
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      IF( NLIN.GT.0 ) THEN
        WRITE(IWR,'(/,A)') 'Chain Decay Fractions:'
        ALLOCATE( ICDSX(1:2,1:NLIN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: ICDSX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( JCDSX(1:2,1:NLIN),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: JCDSX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( ICLX(1:NSOLU),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: ICLX'
          CALL WRMSGS( INDX )
        ENDIF
        DO NL = 1, NLIN
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ADUM(1:) = ' '
          VARB = 'Parent Solute Name'
          NPSL = 0
          CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
          DO NSL = 1,NSOLU
            IF( SOLUT(NSL).EQ.ADUM ) NPSL = NSL
          ENDDO
          BDUM(1:) = ' '
          VARB = 'Daughter Solute Name'
          NDSL = 0
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          DO NSL = 1,NSOLU
            IF( SOLUT(NSL).EQ.BDUM ) NDSL = NSL
          ENDDO
          IF( NPSL.EQ.0 .OR. NDSL.EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Invalid Chain Decay: '//
     &        ADUM(1:NCHA)//': '//BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ELSEIF( NPSL.EQ.NDSL ) THEN
            INDX = 4
            CHMSG = 'Invalid Chain Decay (Parent = Progeny): '//
     &        ADUM(1:NCHA)//': '//BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ELSEIF( NPSL.GT.NDSL ) THEN
            INDX = 4
            CHMSG = 'Invalid Chain Decay (Parent # > Progeny #): '//
     &        ADUM(1:NCHA)//': '//BDUM(1:NCHB)
            CALL WRMSGS( INDX )
          ENDIF
          VARB = 'Chain Decay Fraction'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,CHDF(NPSL,NDSL))
          WRITE(IWR,'(2X,5A,1PE11.4)') 'From ',
     &      ADUM(1:NCHA),' to ',BDUM(1:NCHB),': ',CHDF(NPSL,NDSL)
          ICDSX(1,NL) = NPSL
          ICDSX(2,NL) = NDSL
          JCDSX(1,NL) = NPSL
          JCDSX(2,NL) = NDSL
        ENDDO
        DO NPSL = 1,NSOLU
          CHDFX = 0.D+0
          DO NDSL = NPSL+1,NSOLU
            CHDFX = CHDFX + CHDF(NPSL,NDSL)        
          ENDDO
          IF( ABS(CHDFX-1.D+0)/EPSL.GT.EPSL .AND. 
     &      ABS(CHDFX)/EPSL.GT.EPSL ) THEN
            INDX = 4
            CHMSG = 'Chain Decay Fraction Summation \= 1.0 and \= 0.0'
            CALL WRMSGS( INDX )
          ENDIF
        ENDDO
!
!---    Find chain decay series  ---
!
        NC = 0
        MC = 0
        DO
          MC = MC + 1
          MC0 = MC
!
!---      Zero indices of chain decay list  ---
!
          DO NSL = 1,NSOLU
            ICLX(NSL) = 0
          ENDDO
!
!---      Find the lowest numbered parent  ---
!
          NP = NSOLU+1
          DO NL = 1,NLIN
            IF( ICDSX(1,NL).LE.NP ) THEN
              NP = ICDSX(1,NL)
            ENDIF
          ENDDO
          IF( NP.EQ.NSOLU+1 ) EXIT
          NC = NC + 1
          ICLX(NP) = 1
!
!---      Follow the chain decay from the parent  ---
!
          DO NP = 1,NSOLU
            IF( ICLX(NP).EQ.0 ) CYCLE
            DO NL = 1,NLIN
              IF( NP.EQ.ICDSX(1,NL) ) THEN
                ICLX(ICDSX(1,NL)) = 1
                ICDSX(1,NL) = NSOLU+2
                ICLX(ICDSX(2,NL)) = 1
                ICDSX(2,NL) = NSOLU+2
              ENDIF
            ENDDO
          ENDDO
          DO NP = 1,NSOLU
            IF( ICLX(NP).EQ.0 ) CYCLE
            MC = MC + 1
            IBCDS(MC) = NP
          ENDDO
          IBCDS(MC0) = MC - MC0
        ENDDO
        NBCDS = NC+1
!
!---    Put all solutes not part of a chain-decay series in the
!       the last chain decay series, to be treated as individual
!       decays  ---
!
        MC0 = MC
!
!---    Loop over all solutes checking for those not in a chain-decay
!       series  ---
!
        DO NSL = 1,NSOLU
          IFIND = 0
          KC = 0
!
!---      Loop over the number of active chain-decay series  ---
!
          DO NC = 1,NBCDS-1
            KC = KC + 1
            KC0 = KC
            DO M = 1,IBCDS(KC0)
              KC = KC + 1
              IF( NSL.EQ.IBCDS(KC) ) IFIND = 1
            ENDDO
          ENDDO
          IF( IFIND.EQ.0 ) THEN
            MC = MC + 1
            IBCDS(MC) = NSL
          ENDIF
        ENDDO
        IBCDS(MC0) = MC - MC0
        IF( ALLOCATED(ICDSX) ) THEN
          DEALLOCATE( ICDSX,STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            INDX = 3
            CHMSG = 'Deallocation Error: ICDSX'
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  No chain decay series  ---
!
      ELSE
        NBCDS = 0
        IBCDS(1) = NSOLU
        DO NSL = 1,NSOLU
          IBCDS(NSL+1) = NSL
        ENDDO
      ENDIF
!
!---  Loop over chains  ---
!
      MC = 0
      DO NC = 1,NBCDS
        MC = MC + 1
        NSC = IBCDS(MC)
        IF( NSC.LE.1 ) CYCLE
        ALLOCATE( ISCX(1:NSC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: ISCX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IPSPX(1:NSC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IPSPX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IBCDPX(1:NSC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IBCDPX'
          CALL WRMSGS( INDX )
        ENDIF
        DO NS = 1,NSC
          MC = MC + 1
          ISCX(NS) = IBCDS(MC)
        ENDDO
!
!---    Adjacency matrix  ---
!
        ALLOCATE( IADJM(1:NSC,1:NSC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IADJM'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( JADJM(1:NSC,1:NSC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: JADJM'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( KADJM(1:NSC,1:NSC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: KADJM'
          CALL WRMSGS( INDX )
        ENDIF
        DO I = 1,NSC
          DO J = 1,NSC
            IADJM(I,J) = 0
            JADJM(I,J) = 0
            KADJM(I,J) = 0
          ENDDO
        ENDDO
        DO NL = 1,NLIN
          NC1 = 0
          NC2 = 0
          DO IC = 1,NSC
            IF( JCDSX(1,NL).EQ.ISCX(IC) ) NC1 = IC
            IF( JCDSX(2,NL).EQ.ISCX(IC) ) NC2 = IC
          ENDDO
          IF( NC1.NE.0 .AND. NC2.NE.0 ) THEN
            IADJM(NC1,NC2) = 1
            JADJM(NC1,NC2) = 1
          ENDIF
        ENDDO
        NBCDP(NC) = JADJM(1,NSC)
        DO NSL = 2,NSC
          CALL IMAT_MUL(IADJM,JADJM,KADJM,NSC,NSC,NSC)
          DO I = 1,NSC
            DO J = 1,NSC
              JADJM(I,J) = KADJM(I,J)
              KADJM(I,J) = 0
            ENDDO
          ENDDO
          NBCDP(NC) = NBCDP(NC) + JADJM(1,NSC)
        ENDDO
        IF( NBCDP(NC).EQ.0 ) THEN
          INDX = 7
          IMSG = NC
          CHMSG = 'No Full Decay-Chain Path: Chain'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Use the adjacency matrix to determine the decay paths  ---
!
        DO I = 1,NSC
          IBCDPX(I) = 0
          IPSPX(I) = 0
          DO J = 1,NSC
            IF( IADJM(I,J).EQ.1 ) THEN
              IPSPX(I) = J
              EXIT
            ENDIF
          ENDDO
        ENDDO
        IC = 1
        JC = 1
        IBCDPX(1) = 1
        DO NP = IC,NBCDP(NC)
          DO
            IC = IC + 1
            IBCDPX(IC) = IPSPX(JC)
            JC = IPSPX(JC)
            IF( JC.EQ.NSC ) EXIT
          ENDDO
          DO I = 1,IC
!            IBCDP(I,NP,NC) = ISCX(IBCDPX(I))
            IBCDP(I,NP,NC) = IBCDPX(I)
          ENDDO
          IF( NP.EQ.NBCDP(NC) ) EXIT
!
!---      Recursion, stepping down the adjacency matrix to find
!         additional paths  ---
!
          IFIND = 0
          DO
            IBCDPX(IC) = 0
            IC = IC - 1
            JC = IBCDPX(IC)
            IF( IC.EQ.0 ) THEN
              INDX = 7
              IMSG = NC
              CHMSG = 'Path Finding Error: Chain'
              CALL WRMSGS( INDX )
            ENDIF
            DO J = IPSPX(JC)+1,NSC
              IF( IADJM(JC,J).EQ.1 ) THEN
                IFIND = 1
                EXIT
              ENDIF
            ENDDO
            IF( IFIND.EQ.1 ) THEN
              IPSPX(JC) = J
              EXIT
            ELSE
              DO J = 1,IPSPX(JC)
                IF( IADJM(JC,J).EQ.1 ) THEN
                  IPSPX(JC) = J
                  EXIT
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        IF( ALLOCATED(IADJM) ) THEN
          DEALLOCATE( IADJM,STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            INDX = 3
            CHMSG = 'Deallocation Error: IADJM'
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        IF( ALLOCATED(JADJM) ) THEN
          DEALLOCATE( JADJM,STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            INDX = 3
            CHMSG = 'Deallocation Error: JADJM'
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        IF( ALLOCATED(KADJM) ) THEN
          DEALLOCATE( KADJM,STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            INDX = 3
            CHMSG = 'Deallocation Error: KADJM'
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        IF( ALLOCATED(ISCX) ) THEN
          DEALLOCATE( ISCX,STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            INDX = 3
            CHMSG = 'Deallocation Error: ISCX'
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        IF( ALLOCATED(IPSPX) ) THEN
          DEALLOCATE( IPSPX,STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            INDX = 3
            CHMSG = 'Deallocation Error: IPSPX'
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        IF( ALLOCATED(IBCDPX) ) THEN
          DEALLOCATE( IBCDPX,STAT=ISTAT )
          IF( ISTAT.NE.0 ) THEN
            INDX = 3
            CHMSG = 'Deallocation Error: IBCDPX'
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
      ENDDO
      IF( ALLOCATED(JCDSX) ) THEN
        DEALLOCATE( JCDSX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: JCDSX'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IF( ALLOCATED(ICLX) ) THEN
        DEALLOCATE( ICLX,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: ICLX'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF

!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTF_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTP_GT
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
!     Reads the solute/porous media interaction card for the
!     dispersivities, half-lives, and partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE GRID
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
      CHARACTER*64 ADUM,RDUM,UNTS
      CHARACTER*512 CHDUM
      INTEGER NCH
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDTP_GT'
!
!---  Write card information to output file  ---
!
      CARD = 'Solute/Porous Media Interaction Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      IDISP = 0
!
!---  Loop over the rock/soil saturation information lines  ---
!
      N = 0
      IJK = 0
      ISGRP = 0
   10 CONTINUE
        IF( N.GE.NROCK .OR. IJK.GT.0 ) GOTO 600
        ISTART = 1
        VARB = 'Rock Name: '
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,RDUM)
!
!---  IJK, KIJ, or JKI indexing  ---
!
        IF( INDEX(RDUM(1:),'indexing').NE.0 ) THEN
          IF( INDEX(ROCK(1)(1:),'indexing').EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Indexing Option Not Declared ' // 
     &        'in Rock/Soil Zonation Card'
            CALL WRMSGS( INDX )
          ENDIF
          IF( INDEX(RDUM,'ijk').NE.0 ) THEN
            IJK = 1
          ELSEIF( INDEX(RDUM,'jki').NE.0 ) THEN
            IJK = 2
          ELSEIF( INDEX(RDUM,'kij').NE.0 ) THEN
            IJK = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Indexing Option' // RDUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
          GOTO 220
        ENDIF
!
!---  Search known rock types for a matching type ---
!
        DO 100 M = 1, NROCK
          IF( RDUM.EQ.ROCK(M)) THEN
            IROCK = M
            GOTO 200
          ENDIF
  100   CONTINUE
!
!---  Search known scaling groups for a matching type ---
!
        IF( ISLC(19).EQ.1 ) THEN
          DO 110 M = 1,NSCALE
             IF( RDUM.EQ.SCALNM(M) ) THEN
                ISGRP = M
                IROCK = 1
                GOTO 200
             ENDIF
  110     CONTINUE
          INDX = 2
          CHMSG = 'Unrecognized Rock/Soil Type or Scaling Group: '
     &      // RDUM(1:NCH)
          CALL WRMSGS( INDX )
          GOTO 10
        ENDIF
        INDX = 2
        CHMSG = 'Unrecognized Rock/Soil Type: ' // RDUM(1:NCH)
        CALL WRMSGS( INDX )
        GOTO 10
  200   CONTINUE
!
!---  Loop over rock/soils within scaling group  ---
!
        IF( ISLC(19).EQ.1 .AND. ISGRP.NE.0 ) THEN
          DO 202 M = IROCK,NROCK
            IF( ISCALE(M).EQ.ISGRP ) THEN
              IROCK = M
              GOTO 204
            ENDIF
  202     CONTINUE
        ENDIF
  204   CONTINUE
!
!---    Write rock/soil name  ---
!
        WRITE (IWR,'(/,2A)') 'Rock/Soil Name: ',ROCK(IROCK)
        N = N + 1
  220   CONTINUE
!
!---  Longitudinal dispersivity  ---
!
      VARB = 'Longitudinal Dispersivity: '
      IF( IJK.GT.0 ) THEN
        UNTS = 'm'
        IUNM = 1
        CALL RDIJK( ISTART,IJK,CHDUM,UNTS,DISPL )
        IDISP = 1
      ELSE
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DISPL(IROCK))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',DISPL(IROCK)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DISPL(IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',DISPL(IROCK),', m)'
        IF( DISPL(IROCK).GE.SMALL ) IDISP = 1
      ENDIF
!
!---  Transverse dispersivity  ---
!
      VARB = 'Transverse Dispersivity: '
      IF( IJK.GT.0 ) THEN
        UNTS = 'm'
        IUNM = 1
        CALL RDIJK( ISTART,IJK,CHDUM,UNTS,DISPT )
        IDISP = 1
      ELSE
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DISPT(IROCK))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',DISPT(IROCK)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DISPT(IROCK),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',DISPT(IROCK),', m)'
        IF( DISPT(IROCK).GE.SMALL ) IDISP = 1
      ENDIF
!
!---  Loop over number of solutes or radionuclides  ---
!
        DO 500 NS = 1,NSOLU
!
!---      Skip reads for stationary solutes  ---
!
          IF( IEDL(NS).EQ.4 ) CYCLE
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          VARB = 'Solute Name'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---      Search known solutes for matching name  ---
!
          IFIND = 0
          DO NSL = 1,NSOLU
            IF( ADUM.EQ.SOLUT(NSL)) THEN
              IFIND = 1
              EXIT
            ENDIF
          ENDDO
          IF( IFIND.EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Unrecognized Solute Name: '
     &        // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
          WRITE(IWR,'(/,2A)') 'Solute Name:',SOLUT(NSL)
!
!---     Solid-aqueous partition coefficient  ---
!
          UNTS = 'm^3/kg'
          IUNKG = -1
          IUNM = 3
          IF( IJK.GT.0 ) THEN
            INDX = 1
            LNDX = 5
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,PCSL(1,1,NSL),INDX,LNDX )
            DO IROCK = 1,NFLD
              PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-20 )
            ENDDO
          ELSE
            IDFLT = 1
            VARB = 'Solid-Aqueous Partition Coefficient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PCSL(1,IROCK,NSL)
            INDX = 0
            CALL RDUNIT(UNTS,PCSL(1,IROCK,NSL),INDX)
            PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-12 )
          ENDIF
!
!---      van Schaik and Kemper Empirical Aqueous Diffusion Model  ---
!
          IF( IEDL(NSL).EQ.2 ) THEN
            IF( IJK.GT.0 ) THEN
              UNTS = 'm^2/s'
              IUNM = 2
              IUNS = -1
              INDX = 1
              LNDX = 3
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SDCL(1,1,NSL),
     &          INDX,LNDX )
              DO IROCK = 1,NFLD
                SDCL(1,IROCK,NSL) = MAX( SDCL(1,IROCK,NSL),1.D-20 )
              ENDDO
              UNTS = 'null'
              INDX = 2
              LNDX = 3
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SDCL(1,1,NSL),
     &          INDX,LNDX )
              DO IROCK = 1,NFLD
                SDCL(1,IROCK,NSL) = MAX( SDCL(1,IROCK,NSL),1.D-20 )
              ENDDO
              UNTS = 'null'
              INDX = 3
              LNDX = 3
              CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SDCL(1,1,NSL),
     &          INDX,LNDX )
              DO IROCK = 1,NFLD
                SDCL(1,IROCK,NSL) = MAX( SDCL(1,IROCK,NSL),1.D-20 )
              ENDDO
            ELSE
              UNTS = 'm^2/s'
              IUNM = 2
              IUNS = -1
              VARB = 'Aqueous Molecular Diffusion Coefficient'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SDCL(1,IROCK,NSL))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &          UNTS(1:NCH),': ',SDCL(1,IROCK,NSL)
              INDX = 0
              CALL RDUNIT(UNTS,SDCL(1,IROCK,NSL),INDX)
              WRITE(IWR,'(A,1PE11.4,A)') ' (',
     &          SDCL(1,IROCK,NSL),', m^2/s)'
              VARB = 'a Constant'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SDCL(2,IROCK,NSL))
              WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',
     &          SDCL(2,IROCK,NSL)
              VARB = 'b Constant'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,SDCL(3,IROCK,NSL))
              WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',
     &          SDCL(3,IROCK,NSL)
            ENDIF
          ENDIF
  500   CONTINUE
!
!---  Read next rock/soil type or scaling group  ---
!
      IF( N.LT.NROCK ) WRITE(IWR,'(/)')
      GOTO 10
  600 CONTINUE
      WRITE(IWR,'(/)')
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTP_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKGS_GT( ASL_F,ASL_M,ASLX,ASLMINX,ESGTX,PGX,PLX,
     &  RKGX,SGX,SGTX,SLX,IZN )
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
!     Gas relative permeability.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/RKGS_GT'
      ASLZ = MAX( ASLX,0.D+0 )
!
!---  Gas relative permeability function  ---
!
      IF( IRPG(IZN).EQ.0 ) THEN
        RKGX = RPGC(3,IZN)
!
!---  Mualem porosity distribution function  ---
!
      ELSEIF( IRPG(IZN).EQ.1 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 .OR.
     &    ISCHR(IZN).EQ.102 ) THEN
          ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
          SGP = 1.D+0 - ASLZ
          RKGX = SQRT(SGP)*((1.D+0-ASLZ**(1.D+0/RPGC(3,IZN)))
     &      **RPGC(3,IZN))**2
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 ) THEN
          ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
          SGP = 1.D+0 - ASLZ
          RKGX = SQRT(SGP)*(1.D+0-ASLZ**(1.D+0+1.D+0/RPGC(3,IZN)))**2
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
        ELSEIF( ISCHR(IZN).EQ.3 ) THEN
          SGPM = 1.D+0-ASL_M
          RKGM = SQRT(SGPM)*((1.D+0-ASL_M**(1.D+0/RPGC(3,IZN)))
     &      **RPGC(3,IZN))**2
          SGPF = 1.D+0-ASL_F
          RKGF = SQRT(SGPF)*((1.D+0-ASL_F**(1.D+0/RPGC(1,IZN)))
     &      **RPGC(1,IZN))**2
          RKGX = ( PERM(4,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(7,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(4,IZN)*(1.D+0-POR(4,IZN)) + PERM(7,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGY = ( PERM(5,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(8,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(5,IZN)*(1.D+0-POR(4,IZN)) + PERM(8,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGZ = ( PERM(6,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(9,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(6,IZN)*(1.D+0-POR(4,IZN)) + PERM(9,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGX = MAX( RKGX,RKGY,RKGZ )
!
!---    Single-pressure dual-porosity Brooks and Corey  ---
!
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          SGPM = 1.D+0-ASL_M
          RKGM = SGPM**(2.5D+0 + 2.0D+0/RPGC(3,IZN))
          SGPF = 1.D+0-ASL_F
          RKGF = SGPF**(2.5D+0 + 2.0D+0/RPGC(1,IZN))
          RKGX = ( PERM(4,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(7,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(4,IZN)*(1.D+0-POR(4,IZN)) + PERM(7,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGY = ( PERM(5,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(8,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(5,IZN)*(1.D+0-POR(4,IZN)) + PERM(8,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGZ = ( PERM(6,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(9,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(6,IZN)*(1.D+0-POR(4,IZN)) + PERM(9,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGX = MAX( RKGX,RKGY,RKGZ )
        ENDIF
!
!---  Burdine porosity distribution function  ---
!
      ELSEIF( IRPG(IZN).EQ.2 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 .OR.
     &    ISCHR(IZN).EQ.102 ) THEN
          ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
          SGP = 1.D+0 - ASLZ
          RKGX = (SGP**2)*((1.D+0-ASLZ**(1.D+0/RPGC(3,IZN)))
     &      **RPGC(3,IZN))
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 ) THEN
          ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
          SGP = 1.D+0 - ASLZ
          RKGX = (SGP**2)*(1.D+0-ASLZ**(1.D+0 + 2.D+0/RPGC(3,IZN)))
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
        ELSEIF( ISCHR(IZN).EQ.3 ) THEN
          SGPM = 1.D+0-ASL_M
          RKGM = (SGPM**2)*((1.D+0-ASL_M**(1.D+0/RPGC(3,IZN)))
     &      **RPGC(3,IZN))
          SGPF = 1.D+0-ASL_F
          RKGF = (SGPF**2)*((1.D+0-ASL_F**(1.D+0/RPGC(1,IZN)))
     &      **RPGC(1,IZN))
          RKGX = ( PERM(4,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(7,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(4,IZN)*(1.D+0-POR(4,IZN)) + PERM(7,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGY = ( PERM(5,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(8,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(5,IZN)*(1.D+0-POR(4,IZN)) + PERM(8,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGZ = ( PERM(6,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(9,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(6,IZN)*(1.D+0-POR(4,IZN)) + PERM(9,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGX = MAX( RKGX,RKGY,RKGZ )
!
!---    Single-pressure dual-porosity Brooks and Corey  ---
!
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          SGPM = 1.D+0-ASL_M
          RKGM = SGPM**(3.0D+0 + 2.0D+0/RPGC(3,IZN))
          SGPF = 1.D+0-ASL_F
          RKGF = SGPF**(3.0D+0 + 2.0D+0/RPGC(1,IZN))
          RKGX = ( PERM(4,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(7,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(4,IZN)*(1.D+0-POR(4,IZN)) + PERM(7,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGY = ( PERM(5,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(8,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(5,IZN)*(1.D+0-POR(4,IZN)) + PERM(8,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGZ = ( PERM(6,IZN)*RKGM*(1.D+0-POR(4,IZN)) +
     &      PERM(9,IZN)*RKGF*POR(4,IZN) )/
     &      ( PERM(6,IZN)*(1.D+0-POR(4,IZN)) + PERM(9,IZN)*POR(4,IZN)
     &      + SMALL )
          RKGX = MAX( RKGX,RKGY,RKGZ )
        ENDIF
!
!---  Modified-Corey relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.3 ) THEN
        SLRX = RPGC(1,IZN)
        SGRX = RPGC(3,IZN)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        ASLPX = SLPX + ESGTX
        RKGX = RPGC(2,IZN)*((1.D+0-ASLPX)**2)*(1.D+0-ASLPX**2)
!
!---  Fatt and Klikoff relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.4 ) THEN
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        RKGX = (1.D+0-ASLZ)**3
!
!---  Free-Corey relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.7 ) THEN
        SLRX = RPGC(3,IZN)
        SGRX = RPGC(4,IZN)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        ASLPX = SLPX + ESGTX
        SGPX = 1.D+0 - ASLPX
        RKGX = RPGC(1,IZN)*(SGPX**(RPGC(2,IZN)))
!
!---  Sandia relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.8 ) THEN
        RKGX = MAX( 1.D+0-RKGX,0.D+0 )
!
!---  Classical-Corey relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.17 ) THEN
        SLRX = RPGC(3,IZN)
        SGRX = RPGC(4,IZN)
        SLPX = (SLX-SLRX)/(1.D+0-SLRX-SGRX)
        SLPX = MIN( MAX( SLPX,0.D+0 ),1.D+0 )
        ASLPX = SLPX + ESGTX
        RKGX = RPGC(1,IZN)*(1.D+0-(ASLPX**(RPGC(2,IZN)/2.D+0)))*
     &    ((1.D+0-ASLPX)**(RPGC(2,IZN)/2.D+0))
!
!---  van Genuchten gas relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.9 ) THEN
        SLRX = RPGC(1,IZN)
        SGRX = RPGC(3,IZN)
        ASLZ = MIN( ASLZ+ESGTX,1.D+0 )
        SGX = (1.D+0-ASLZ)*(1.D+0-SLRX)
        SGPX = MIN( MAX( (SGX-SGRX)/(1.D+0-SGRX),0.D+0 ),1.D+0 )
        RKGX = SQRT(SGPX)*((1.D+0-((1.D+0-SGPX)**(1.D+0/RPGC(2,IZN)))
     &      **RPGC(2,IZN))**2)
!
!---  Tabular gas relative permeability versus gas saturation
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( IRPG(IZN).EQ.10 ) THEN
        SGZ = MAX( SGX-SGTX,0.D+0 )
        ITBX = 0
        RKGX = FNTBLY( SGZ,IRGTBL(1,IZN),IRGTBL(2,IZN),ITBX )
!
!---  Tabular gas relative permeability versus gas saturation
!     w/ cubic spline interpolation  ---
!
      ELSEIF( IRPG(IZN).EQ.11 ) THEN
        SGZ = MAX( SGX-SGTX,0.D+0 )
        RKGX = FSPLNY( SGZ,IRGTBL(1,IZN),IRGTBL(2,IZN) )
!
!---  Tabular gas relative permeability versus ln(capillary head)
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( IRPG(IZN).EQ.12 ) THEN
        ITBX = 0
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKGX = FNTBLY( HDGLX,IRGTBL(1,IZN),IRGTBL(2,IZN),ITBX )
!
!---  Tabular gas relative permeability versus ln(capillary head)
!     w/ cubic spline interpolation  ---
!
      ELSEIF( IRPG(IZN).EQ.13 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKGX = FSPLNY( HDGLX,IRGTBL(1,IZN),IRGTBL(2,IZN) )
!
!---  Tabular gas relative permeability versus capillary head
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( IRPG(IZN).EQ.14 ) THEN
        ITBX = 0
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKGX = FNTBLY( HDGL,IRGTBL(1,IZN),IRGTBL(2,IZN),ITBX )
!
!---  Tabular gas relative permeability versus capillary head
!     w/ cubic spline interpolation  ---
!
      ELSEIF( IRPG(IZN).EQ.15 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKGX = FSPLNY( HDGL,IRGTBL(1,IZN),IRGTBL(2,IZN) )
!
!---  Doughty gas relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.18 ) THEN
        RKGMX = RPGC(1,IZN)
        GAMMAX = RPGC(2,IZN)
        SLRX = RPGC(3,IZN)
        CMX = RPGC(4,IZN)
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ASLPX = MIN( MAX( (ESLX+ESGTX),0.D+0 ),1.D+0 )
        RKGX = RKGMX*((1.D+0-ASLPX)**GAMMAX)*
     &    (((1.D+0-ASLPX)**(1.D+0/CMX))**(2.D+0*CMX))
!
!---  Doughty drainage-imbibition gas relative permeability 
!     function  ---
!
      ELSEIF( IRPG(IZN).EQ.19 ) THEN
        RKGMX = RPGC(1,IZN)
        GAMMAX = RPGC(2,IZN)
        SLRDX = RPGC(3,IZN)
        CMX = RPGC(4,IZN)
        SLRIX = RPGC(5,IZN)
        SGTMX = RPGC(6,IZN)
        SLDX = MIN( SLX+SGTX,1.D+0 )
        SLIX = MIN( SLX+SGTMX,1.D+0 )
        ESLDX = MIN( MAX( (SLDX-SLRDX)/(1.D+0-SLRDX),0.D+0 ),1.D+0 )
        ESLIX = MIN( MAX( (SLIX-SLRIX)/(1.D+0-SLRIX),0.D+0 ),1.D+0 )
        RKGDX = RKGMX*((1.D+0-ESLDX)**GAMMAX)*
     &    (((1.D+0-ESLDX)**(1.D+0/CMX))**(2.D+0*CMX))
        RKGIX = RKGMX*((1.D+0-ESLIX)**GAMMAX)*
     &    (((1.D+0-ESLIX)**(1.D+0/CMX))**(2.D+0*CMX))
        RSLX = MIN( MAX( ASLMINX/ASLX,0.D+0 ),1.D+0 )
        RKGX = RSLX*RKGDX + (1.D+0-RSLX)*RKGIX
!
!---  Variable-Corey relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.20 ) THEN
        SLRX = RPGC(1,IZN)
        SGRX = RPGC(3,IZN)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        ASLPX = SLPX + ESGTX
        IF( SLRX.GT.EPSL ) THEN
          BX = (1.D+0-RPGC(1,IZN))/SLRX
          AX = 1.D+0 - BX
          DSLX = 5.D-2
          RKG1X = AX + BX*(1.D+0-SLX)
          RKG2X = RPGC(2,IZN)*((1.D+0-ASLPX)**2)*
     &      (1.D+0-ASLPX**RPGC(4,IZN))
          RKG1X = MIN( MAX( RKG1X,0.D+0 ),1.D+0 )
          RKG2X = MIN( MAX( RKG2X,0.D+0 ),1.D+0 )
          SIGMAX = 1.D+0/(1.D+0 + EXP((SLRX-SLX)/DSLX))
          RKGX = (1.D+0-SIGMAX)*RKG1X + SIGMAX*RKG2X
        ELSE
          RKGX = RPGC(2,IZN)*((1.D+0-ASLPX)**2)*
     &      (1.D+0-ASLPX**RPGC(4,IZN))
        ENDIF
        
!
!---  Extended power law relative permeability function  ---
!
      ELSEIF( IRPG(IZN).EQ.21 ) THEN
        SLRX = RPGC(3,IZN)
        SGRX = RPGC(4,IZN)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        ASLPX = SLPX + ESGTX
        SGPX = 1.D+0 - ASLPX
        IF( SLRX.GT.EPSL ) THEN
          BX = (1.D+0-RPGC(1,IZN))/SLRX
          AX = 1.D+0 - BX
          DSLX = 5.D-2
          RKG1X = AX + BX*(1.D+0-SLX)
          RKG2X = RPGC(1,IZN)*(SGPX**(RPGC(2,IZN)))
          RKG1X = MIN( MAX( RKG1X,0.D+0 ),1.D+0 )
          RKG2X = MIN( MAX( RKG2X,0.D+0 ),1.D+0 )
          SIGMAX = 1.D+0/(1.D+0 + EXP((SLRX-SLX)/DSLX))
          RKGX = (1.D+0-SIGMAX)*RKG1X + SIGMAX*RKG2X
        ELSE
          RKGX = RPGC(1,IZN)*(SGPX**(RPGC(2,IZN)))
        ENDIF
      ENDIF
!
!---  Klinkenberg effect  ---
!
      IF( IRPG(IZN).GE.100 .AND. IRPG(IZN).LT.200 ) THEN
        RKGX = RKGX*(1.D+0 + RPGC(5,IZN)*((PGX+PATM)**
     &    (RPGC(6,IZN)-1.D+0)))
      ENDIF
      RKGX = MIN( MAX( RKGX,0.D+0 ),1.D+0 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKGS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKLS_GT( ASL_F,ASL_M,ESLX,PGX,PLX,RKLX,SLX,IZN )
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
!     Aqueous relative permeability.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/RKLS_GT'
      ESLZ = MAX( ESLX,0.D+0 )
!
!---  Constant relative permeability function  ---
!
      IF( MOD( IRPL(IZN),100 ).EQ.0 ) THEN
        RKLX(1) = RPLC(2,IZN)
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---    Single-pressure dual-porosity saturation functions  ---
!
        IF( ISCHR(IZN).EQ.3 .OR. ISCHR(IZN).EQ.4 ) THEN
          RKLM = RPLC(2,IZN)
          RKLF = RPLC(1,IZN)
          RKLX(1) = ( PERM(4,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(7,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(4,IZN)*(1.D+0-POR(4,IZN)) + PERM(7,IZN)*POR(4,IZN)
     &      + SMALL )
          RKLX(2) = ( PERM(5,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(8,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(5,IZN)*(1.D+0-POR(4,IZN)) + PERM(8,IZN)*POR(4,IZN)
     &      + SMALL )
          RKLX(3) = ( PERM(6,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(9,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(6,IZN)*(1.D+0-POR(4,IZN)) + PERM(9,IZN)*POR(4,IZN)
     &      + SMALL )
          ENDIF
!
!---  Mualem-irreducible porosity distribution function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.21 ) THEN
        SLRX = RPLC(1,IZN)
        SLPX = MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 )
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 .OR.
     &    ISCHR(IZN).EQ.201 ) THEN
          RKLX(1) = SQRT(SLPX)*((1.D+0-(1.D+0-SLPX**(1.D+0/RPLC(2,IZN)))
     &      **RPLC(2,IZN))**2)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 ) THEN
          RKLX(1) = SLPX**(2.5D+0 + 2.0D+0/RPLC(2,IZN))
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
        ENDIF
!
!---  Mualem porosity distribution function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.1 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 .OR.
     &    ISCHR(IZN).EQ.201 ) THEN
          RKLX(1) = SQRT(ESLZ)*((1.D+0-(1.D+0-ESLZ**(1.D+0/RPLC(2,IZN)))
     &      **RPLC(2,IZN))**2)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 ) THEN
          RKLX(1) = ESLZ**(2.5D+0 + 2.0D+0/RPLC(2,IZN))
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
          ELSEIF( ISCHR(IZN).EQ.3  ) THEN
            RKLM = (1.D+0-(ASL_M**(1.D+0/RPLC(2,IZN))))
            RKLM = SQRT(ASL_M)*((1.D+0-RKLM**RPLC(2,IZN))**2)
            RKLF = (1.D+0-(ASL_F**(1.D+0/RPLC(1,IZN))))
            RKLF = SQRT(ASL_F)*((1.D+0-RKLF**RPLC(1,IZN))**2)
            RKLX(1) = ( PERM(4,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(7,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(4,IZN)*(1.D+0-POR(4,IZN)) + PERM(7,IZN)*POR(4,IZN)
     &      + SMALL )
            RKLX(2) = ( PERM(5,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(8,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(5,IZN)*(1.D+0-POR(4,IZN)) + PERM(8,IZN)*POR(4,IZN)
     &      + SMALL )
            RKLX(3) = ( PERM(6,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(9,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(6,IZN)*(1.D+0-POR(4,IZN)) + PERM(9,IZN)*POR(4,IZN)
     &      + SMALL )
!
!---    Single-pressure dual-porosity Brooks and Corey  ---
!
          ELSEIF( ISCHR(IZN).EQ.4  ) THEN
            RKLM = ASL_M**(2.5D+0 + 2.0D+0/RPLC(2,IZN))
            RKLF = ASL_F**(2.5D+0 + 2.0D+0/RPLC(1,IZN))
            RKLX(1) = ( PERM(4,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(7,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(4,IZN)*(1.D+0-POR(4,IZN)) + PERM(7,IZN)*POR(4,IZN)
     &      + SMALL )
            RKLX(2) = ( PERM(5,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(8,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(5,IZN)*(1.D+0-POR(4,IZN)) + PERM(8,IZN)*POR(4,IZN)
     &      + SMALL )
            RKLX(3) = ( PERM(6,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(9,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(6,IZN)*(1.D+0-POR(4,IZN)) + PERM(9,IZN)*POR(4,IZN)
     &      + SMALL )
        ENDIF
!
!---  Burdine porosity distribution function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
          RKLX(1) = (ESLZ**2)*(1.D+0-(1.D+0-ESLZ**(1.D+0/RPLC(2,IZN)))
     &      **RPLC(2,IZN))
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 ) THEN
          RKLX(1) = ESLZ**(3.0D+0 + 2.0D+0/RPLC(2,IZN))
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
        ELSEIF( ISCHR(IZN).EQ.3 ) THEN
          RKLM = (ASL_M**2)*(1.D+0-(1.D+0-ASL_M**(1.D+0/RPLC(2,IZN)))
     &      **RPLC(2,IZN))
          RKLF = (ASL_F**2)*(1.D+0-(1.D+0-ASL_F**(1.D+0/RPLC(1,IZN)))
     &      **RPLC(1,IZN))
          RKLX(1) = ( PERM(4,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(7,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(4,IZN)*(1.D+0-POR(4,IZN)) + PERM(7,IZN)*POR(4,IZN)
     &      + SMALL )
          RKLX(2) = ( PERM(5,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(8,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(5,IZN)*(1.D+0-POR(4,IZN)) + PERM(8,IZN)*POR(4,IZN)
     &      + SMALL )
          RKLX(3) = ( PERM(6,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(9,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(6,IZN)*(1.D+0-POR(4,IZN)) + PERM(9,IZN)*POR(4,IZN)
     &      + SMALL )
!
!---  Dual porosity Brooks and Corey  ---
!
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          RKLM = ASL_M**(3.0D+0 + 2.0D+0/RPLC(2,IZN))
          RKLF = ASL_F**(3.0D+0 + 2.0D+0/RPLC(1,IZN))
          RKLX(1) = ( PERM(4,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(7,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(4,IZN)*(1.D+0-POR(4,IZN)) + PERM(7,IZN)*POR(4,IZN)
     &      + SMALL )
          RKLX(2) = ( PERM(5,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(8,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(5,IZN)*(1.D+0-POR(4,IZN)) + PERM(8,IZN)*POR(4,IZN)
     &      + SMALL )
          RKLX(3) = ( PERM(6,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(9,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(6,IZN)*(1.D+0-POR(4,IZN)) + PERM(9,IZN)*POR(4,IZN)
     &      + SMALL )
        ENDIF
!
!---  Corey relative permeability function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.3 ) THEN
        RKLX(1) = ESLZ**4
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Fatt and Klikoff relative permeability function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.4 ) THEN
        RKLX(1) = ESLZ**3
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Haverkamp relative permeability function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.5 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 .OR. ISCHR(IZN).EQ.5 ) THEN
          HDEN = SCHR(1,IZN)
        ELSE
          HDEN = ZERO
        ENDIF
        IF( HDGL.LE.HDEN ) THEN
          RKLX(1) = 1.D+0
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
        ELSE
          RKLX(1) = RPLC(1,IZN)/
     &      (RPLC(1,IZN) + (((HDGL-HDEN)/SCHR(5,IZN))**RPLC(2,IZN)))
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
        ENDIF
!
!---  Touma and Vauclin relative permeability function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.6 ) THEN
        RKLX(1) = RPLC(1,IZN)*(ESLZ**RPLC(2,IZN))
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Free Corey relative permeability function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.7 ) THEN
        SLRX = RPLC(3,IZN)
        SGRX = RPLC(4,IZN)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        RKLX(1) = RPLC(1,IZN)*(SLPX**(RPLC(2,IZN)))
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Rijtema-Gardner modified exponential function  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.9 ) THEN
        RKLX(1) = EXP( RPLC(1,IZN)*HDGL + RPLC(2,IZN) )
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Tabular aqueous relative permeability versus gas saturation
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.10 ) THEN
        ITBX = 0
        RKLX(1) = FNTBLY( SLX,IRLTBL(1,IZN),IRLTBL(2,IZN),ITBX )
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Tabular aqueous relative permeability versus gas saturation
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.11 ) THEN
        RKLX(1) = FSPLNY( SLX,IRLTBL(1,IZN),IRLTBL(2,IZN) )
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Tabular aqueous relative permeability versus ln(capillary head)
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.12 ) THEN
        ITBX = 0
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKLX(1) = FNTBLY( HDGLX,IRLTBL(1,IZN),IRLTBL(2,IZN),ITBX )
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Tabular aqueous relative permeability versus ln(capillary head)
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.13 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKLX(1) = FSPLNY( HDGLX,IRLTBL(1,IZN),IRLTBL(2,IZN) )
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Tabular aqueous relative permeability versus capillary head
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.14 ) THEN
        ITBX = 0
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKLX(1) = FNTBLY( HDGL,IRLTBL(1,IZN),IRLTBL(2,IZN),ITBX )
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
!
!---  Tabular aqueous relative permeability versus capillary head
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPL(IZN),100 ).EQ.15 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKLX(1) = FSPLNY( HDGL,IRLTBL(1,IZN),IRLTBL(2,IZN) )
        RKLX(2) = RKLX(1)
        RKLX(3) = RKLX(1)
      ENDIF
!
!---  Polmann anisotropy permeability function  ---
!
      IF( IRPL(IZN).GE.100 .AND. IRPL(IZN).LT.200 ) THEN
        PSI = HDGL*1.D+2
        SKLX = RPLC(5,IZN) - RPLC(10,IZN)*PSI -
     &    RPLC(6,IZN)*RPLC(9,IZN)*( RPLC(7,IZN) -
     &    (RPLC(7,IZN)**2)*PSI - (RPLC(8,IZN)**2)*PSI)/
     &    (1.D+0 + RPLC(10,IZN)*RPLC(9,IZN))
        SIGMA = RPLC(6,IZN)*(((1.D+0 - RPLC(7,IZN)*PSI)**2) +
     &    (RPLC(8,IZN)**2)*(PSI**2))/
     &    (1.D+0 + RPLC(10,IZN)*RPLC(9,IZN))
        SKHX = EXP( SKLX + 5.D-1*SIGMA )
        SKVX = EXP( SKLX - 5.D-1*SIGMA )
        ANISOX = MIN( MAX( SKHX/SKVX,0.D+0 ),RPLC(11,IZN) )
        ANISOX = MAX( ANISOX,RPLC(12,IZN) )
        RKLX(1) = RKLX(3)*ANISOX
        RKLX(2) = RKLX(3)*ANISOX
!
!---  Pruess anisotropy permeability function  ---
!
      ELSEIF( IRPL(IZN).GE.200 .AND. IRPL(IZN).LT.300 ) THEN
        HDGL_CM = HDGL*1.D+2
        ANISOX = RPLC(5,IZN)*(RPLC(6,IZN)**(RPLC(7,IZN)**HDGL_CM))
        ANISOX = MAX( ANISOX,1.D+0 )
        ANISOX = MIN( ANISOX,RPLC(5,IZN) )
        RKLX(1) = RKLX(3)*ANISOX
        RKLX(2) = RKLX(3)*ANISOX
      ENDIF
      RKLX(1) = MIN( MAX( RKLX(1),0.D+0 ),1.D+0 )
      RKLX(2) = MIN( MAX( RKLX(2),0.D+0 ),1.D+0 )
      RKLX(3) = MIN( MAX( RKLX(3),0.D+0 ),1.D+0 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKLS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKLST_GT( ASL_F,ASL_M,ESLX,PGX,PLX,
     &  RKLX,SLX,IZN,ITX )
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
!     Aqueous relative permeability.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/RKLST_GT'
      ESLZ = MAX( ESLX,0.D+0 )
!
!---  Constant relative permeability function  ---
!
      IF( MOD( IRPLT(ITX,IZN),100 ).EQ.0 ) THEN
        RKLX = RPLT(ITX,2,IZN)
!
!---    Single-pressure dual-porosity saturation functions  ---
!
        IF( ISCHR(IZN).EQ.3 .OR. ISCHR(IZN).EQ.4 ) THEN
          RKLM = RPLT(ITX,2,IZN)
          RKLF = RPLT(ITX,1,IZN)
          IMX = ITX + 3
          IFX = ITX + 6
          RKLX = ( PERM(IMX,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(IFX,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(IMX,IZN)*(1.D+0-POR(4,IZN)) + 
     &      PERM(IFX,IZN)*POR(4,IZN) + SMALL )
          ENDIF
!
!---  Mualem-irreducible porosity distribution function  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.21 ) THEN
        SLRX = RPLT(ITX,1,IZN)
        SLPX = MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 )
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 .OR.
     &    ISCHR(IZN).EQ.201 ) THEN
          RKLX = SQRT(SLPX)*((1.D+0-(1.D+0-SLPX
     &      **(1.D+0/RPLT(ITX,2,IZN)))**RPLT(ITX,2,IZN))**2)
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 ) THEN
          RKLX = SLPX**(2.5D+0 + 2.0D+0/RPLT(ITX,2,IZN))
        ENDIF
!
!---  Mualem porosity distribution function  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.1 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 .OR.
     &    ISCHR(IZN).EQ.201 ) THEN
          RKLX = SQRT(ESLZ)*((1.D+0-(1.D+0-ESLZ
     &      **(1.D+0/RPLT(ITX,2,IZN)))**RPLT(ITX,2,IZN))**2)
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 ) THEN
          RKLX = ESLZ**(2.5D+0 + 2.0D+0/RPLT(ITX,2,IZN))
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
          ELSEIF( ISCHR(IZN).EQ.3  ) THEN
            RKLM = (1.D+0-(ASL_M**(1.D+0/RPLT(ITX,2,IZN))))
            RKLM = SQRT(ASL_M)*((1.D+0-RKLM**RPLT(ITX,2,IZN))**2)
            RKLF = (1.D+0-(ASL_F**(1.D+0/RPLT(ITX,1,IZN))))
            RKLF = SQRT(ASL_F)*((1.D+0-RKLF**RPLT(ITX,1,IZN))**2)
            IMX = ITX + 3
            IFX = ITX + 6
            RKLX = ( PERM(IMX,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(IFX,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(IMX,IZN)*(1.D+0-POR(4,IZN))
     &      + PERM(IFX,IZN)*POR(4,IZN) + SMALL )
!
!---    Single-pressure dual-porosity Brooks and Corey  ---
!
          ELSEIF( ISCHR(IZN).EQ.4  ) THEN
            RKLM = ASL_M**(2.5D+0 + 2.0D+0/RPLT(ITX,2,IZN))
            RKLF = ASL_F**(2.5D+0 + 2.0D+0/RPLT(ITX,1,IZN))
            IMX = ITX + 3
            IFX = ITX + 6
            RKLX = ( PERM(IMX,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(IFX,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(IMX,IZN)*(1.D+0-POR(4,IZN))
     &      + PERM(IFX,IZN)*POR(4,IZN) + SMALL )
        ENDIF
!
!---  Burdine porosity distribution function  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.2 ) THEN
!
!---    van Genuchten saturation function  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.101 ) THEN
          RKLX = (ESLZ**2)*(1.D+0-(1.D+0-ESLZ**(1.D+0/RPLT(ITX,2,IZN)))
     &      **RPLT(ITX,2,IZN))
!
!---    Brooks and Corey saturation function  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 ) THEN
          RKLX = ESLZ**(3.0D+0 + 2.0D+0/RPLT(ITX,2,IZN))
!
!---    Single-pressure dual-porosity van Genuchten  ---
!
        ELSEIF( ISCHR(IZN).EQ.3 ) THEN
          RKLM = (ASL_M**2)*(1.D+0-(1.D+0-ASL_M
     &      **(1.D+0/RPLT(ITX,2,IZN)))**RPLT(ITX,2,IZN))
          RKLF = (ASL_F**2)*(1.D+0-(1.D+0-ASL_F
     &      **(1.D+0/RPLT(ITX,1,IZN)))**RPLT(ITX,1,IZN))
          IMX = ITX + 3
          IFX = ITX + 6
          RKLX = ( PERM(IMX,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(IFX,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(IMX,IZN)*(1.D+0-POR(4,IZN))
     &      + PERM(IFX,IZN)*POR(4,IZN) + SMALL )
!
!---  Dual porosity Brooks and Corey  ---
!
        ELSEIF( ISCHR(IZN).EQ.4 ) THEN
          RKLM = ASL_M**(3.0D+0 + 2.0D+0/RPLT(ITX,2,IZN))
          RKLF = ASL_F**(3.0D+0 + 2.0D+0/RPLT(ITX,1,IZN))
          IMX = ITX + 3
          IFX = ITX + 6
          RKLX = ( PERM(IMX,IZN)*RKLM*(1.D+0-POR(4,IZN)) +
     &      PERM(IFX,IZN)*RKLF*POR(4,IZN) )/
     &      ( PERM(IMX,IZN)*(1.D+0-POR(4,IZN))
     &      + PERM(IFX,IZN)*POR(4,IZN) + SMALL )
        ENDIF
!
!---  Corey relative permeability function  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.3 ) THEN
        RKLX = ESLZ**4
!
!---  Fatt and Klikoff relative permeability function  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.4 ) THEN
        RKLX = ESLZ**3
!
!---  Haverkamp relative permeability function  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.5 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &    ISCHR(IZN).EQ.202 .OR. ISCHR(IZN).EQ.5 ) THEN
          HDEN = SCHR(1,IZN)
        ELSE
          HDEN = ZERO
        ENDIF
        IF( HDGL.LE.HDEN ) THEN
          RKLX = 1.D+0
        ELSE
          RKLX = RPLT(ITX,1,IZN)/(RPLT(ITX,1,IZN) + 
     &      (((HDGL-HDEN)/SCHR(5,IZN))**RPLT(ITX,2,IZN)))
        ENDIF
!
!---  Touma and Vauclin relative permeability function  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.6 ) THEN
        RKLX = RPLT(ITX,1,IZN)*(ESLZ**RPLT(ITX,2,IZN))
!
!---  Free Corey relative permeability function  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.7 ) THEN
        SLRX = RPLT(ITX,3,IZN)
        SGRX = RPLT(ITX,4,IZN)
        SLPX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX-SGRX),0.D+0 ),1.D+0 )
        RKLX = RPLT(ITX,1,IZN)*(SLPX**(RPLT(ITX,2,IZN)))
!
!---  Rijtema-Gardner modified exponential function  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.9 ) THEN
        RKLX = EXP( RPLT(ITX,1,IZN)*HDGL + RPLT(ITX,2,IZN) )
!
!---  Tabular aqueous relative permeability versus gas saturation
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.10 ) THEN
        ITBX = 0
        RKLX = FNTBLY( SLX,IRLTBLT(1,IZN,ITX),IRLTBLT(2,IZN,ITX),ITBX )
!
!---  Tabular aqueous relative permeability versus gas saturation
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.11 ) THEN
        RKLX = FSPLNY( SLX,IRLTBLT(1,IZN,ITX),IRLTBLT(2,IZN,ITX) )
!
!---  Tabular aqueous relative permeability versus ln(capillary head)
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.12 ) THEN
        ITBX = 0
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKLX = FNTBLY( HDGLX,IRLTBLT(1,IZN,ITX),
     &    IRLTBLT(2,IZN,ITX),ITBX )
!
!---  Tabular aqueous relative permeability versus ln(capillary head)
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.13 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        HDGLX = LOG(HDGL)
        RKLX = FSPLNY( HDGLX,IRLTBLT(1,IZN,ITX),IRLTBLT(2,IZN,ITX) )
!
!---  Tabular aqueous relative permeability versus capillary head
!     w/ linear interpolation and table truncation beyond limits  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.14 ) THEN
        ITBX = 0
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKLX = FNTBLY( HDGL,IRLTBLT(1,IZN,ITX),IRLTBLT(2,IZN,ITX),ITBX )
!
!---  Tabular aqueous relative permeability versus capillary head
!     w/ cubic spline interpolation  ---
!
      ELSEIF( MOD( IRPLT(ITX,IZN),100 ).EQ.15 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        RKLX = FSPLNY( HDGL,IRLTBLT(1,IZN,ITX),IRLTBLT(2,IZN,ITX) )
      ENDIF
      RKLX = MIN( MAX( RKLX,1.D-24 ),1.D+0 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RKLST_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_GT
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
!     Compute the maximum relative residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
      USE FDVT
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
      SUB_LOG(ISUB_LOG) = '/RSDL_GT'
!
!---  Zero maximum residuals  ---
!
      DO 100 M = 1,ISVC
        RSD(M) = 0.D+0
        NSD(M) = 0
  100 CONTINUE
!
!---  Update primary variables
!
      DO N = 1,NFLD
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
        IZN = IZ(N)
        N_DB = N
        NMD = IXP(N)
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) MPT = IM(IEQT,NMD)
        MPL = IM(IEQW,NMD)
        IF( ISLC(37).EQ.0 ) MPG = IM(IEQA,NMD)
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) MPS = IM(IEQS,NMD)
!
!---    Skip selected nodes in the residual calculation  ---
!
        IF( ISKP(IZN).EQ.1 ) CYCLE
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
!
!---      Energy equation  ---
!
          ACP = ((1.D+0-PORT(2,N))*RHOS(IZN)*CPS(IZN)*T(2,N) +
     &      (PORT(2,N)-PORD(2,N))*RHOL(2,N)*HL(2,N) +
     &      PORD(2,N)*(SL(2,N)*RHOL(2,N)*HL(2,N) +
     &      (YLS(2,N)-XLS(2,N))*SL(2,N)*RHOL(2,N)*HSP(2,N) +
     &      SG(2,N)*RHOG(2,N)*UEG(2,N)))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPT))/TABS,
     &      ABS(RSDL(IEQT,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQT) ) THEN
            RSD(IEQT) = RSDX
            NSD(IEQT) = N
          ENDIF
        ENDIF
!
!---    Saturated system w/o entrapped gas
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - aqueous-air mole fraction
!       NaCl mass - total NaCl brine mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
!
!---      Water mass equation  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPL))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
!
!---      Air mass equation, ignore residual for small aqueous-air  ---
!
          INDX = 0
          CALL REGION_4( T(2,N),PSWX,INDX )
          PWX = MAX( PSWX,PL(2,N)+PATM )
          CALL SOL_BRNS( T(2,N),PWX,XLSMX )
          XLSX = MIN( YLS(2,N),XLSMX )
          XLS(2,N) = XLSX
          CALL SP_B( T(2,N),XLS(2,N),PSBX )
          IF( ISLC(37).EQ.0 ) THEN
            PVAX = MAX( PWX-PSBX,0.D+0 )
            XMLAX = PVAX/HCAW
            XMLA(2,N) = PVA(2,N)/HCAW
            IF( XMLA(2,N).GT.(1.D-6*XMLAX) ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLA(2,N))*DTI*VOL(N)
              RSDX = MIN( ABS(BLU(MPG))/MAX( XMLAX,PATM/HCAW ),
     &          ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD(IEQA) ) THEN
                RSD(IEQA) = RSDX
                NSD(IEQA) = N
              ENDIF
            ENDIF
          ENDIF
!
!---      Salt mass equation, isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS(2,N)*DTI*VOL(N)
            RSDX = MIN( (ABS(BLU(MPS))/XLSMX),
     &        ABS(RSDL(IEQS,N)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            IF( RSDX.GT.RSD(IEQS) ) THEN
              RSD(IEQS) = RSDX
              NSD(IEQS) = N
            ENDIF
          ENDIF
!
!---    Unsaturated system w/ or w/o entrapped gas
!       
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
!
!---      Water mass equation  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPL))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
!
!---      Air mass equation  ---
!
          IF( ISLC(37).EQ.0 ) THEN
            IF( SG(2,N).GT.1.D-3 ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLA(2,N))*DTI*VOL(N)
              RSDX = MIN( ABS(BLU(MPG))/(ABS(PG(2,N))+PATM),
     &          ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD(IEQA) ) THEN
                RSD(IEQA) = RSDX
                NSD(IEQA) = N
              ENDIF
            ENDIF
          ENDIF
!
!---      Salt mass equation, isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS(2,N)*DTI*VOL(N)
            INDX = 0
            CALL REGION_4( T(2,N),PWX,INDX )
            CALL SOL_BRNS( T(2,N),PWX,XLSMX )
            RSDX = MIN( (ABS(BLU(MPS))/XLSMX),
     &        ABS(RSDL(IEQS,N)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            IF( RSDX.GT.RSD(IEQS) ) THEN
              RSD(IEQS) = RSDX
              NSD(IEQS) = N
            ENDIF
          ENDIF
!
!---    Saturated system w/ entrapped gas
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
!
!---      Water mass equation  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPL))/(ABS(PL(2,N))+PATM),
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
!
!---      Air mass equation  ---
!
          IF( ISLC(37).EQ.0 ) THEN
            IF( SG(2,N).GT.1.D-5 ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLA(2,N))*DTI*VOL(N)
              RSDX = MIN( ABS(BLU(MPG)/1.D+1),
     &          ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD(IEQA) ) THEN
                RSD(IEQA) = RSDX
                NSD(IEQA) = N
              ENDIF
            ENDIF
          ENDIF
!
!---      Salt mass equation, isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            ACP = TMS(2,N)*DTI*VOL(N)
            INDX = 0
            CALL REGION_4( T(2,N),PWX,INDX )
            CALL P_IAPWS( T(2,N),PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL SOL_BRNS( T(2,N),PWX,XLSMX )
            RSDX = MIN( (ABS(BLU(MPS))/XLSMX),
     &        ABS(RSDL(IEQS,N)/(ACP+SMALL)) )
            RSDX = RSDX*1.D-1
            IF( RSDX.GT.RSD(IEQS) ) THEN
              RSD(IEQS) = RSDX
              NSD(IEQS) = N
            ENDIF
          ENDIF
!
!---    Fully unsaturated conditions
!
!       Energy - temperature
!       Water mass - water vapor partial pressure
!       Air mass - gas pressure
!       NaCl mass - salt mass  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
!
!---      Water mass equation  ---
!
          INDX = 0
          CALL REGION_4( T(2,N),PWX,INDX )
          PWX = MIN( PWX,PCRW )
          ACP = PORD(2,N)*RHOG(2,N)*SG(2,N)*XGW(2,N)*DTI*VOL(N)
          RSDX = MIN( ABS(BLU(MPL))/PWX,
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
!
!---      Air mass equation  ---
!
          IF( ISLC(37).EQ.0 ) THEN
            PGX = PG(2,N) + PATM
            ACP = PORD(2,N)*RHOG(2,N)*SG(2,N)*XGA(2,N)*DTI*VOL(N)
            RSDX = MIN( ABS(BLU(MPG))/ABS(PGX),
     &        ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQA) ) THEN
              RSD(IEQA) = RSDX
              NSD(IEQA) = N
            ENDIF
          ENDIF
!
!---      Salt mass equation, isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            CALL DENS_S( T(2,N),PGX,RHOSPX )
            ACP = TMS(2,N)*DTI*VOL(N)
            RSDX = MIN( (ABS(BLU(MPS))/RHOSPX),
     &        ABS(RSDL(IEQS,N)/(ACP+SMALL)) )
            IF( RSDX.GT.RSD(IEQS) ) THEN
              RSD(IEQS) = RSDX
              NSD(IEQS) = N
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!
!---  Assign a convergence index  ---
!
      RSDX = 1.D-20
      DO M = 1,ISVC
        IF( RSD(M).GT.RSDMX ) ICNV = 2
        RSDX = MAX( RSD(M),RSDX )
      ENDDO
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) ICNV = 1
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.1 ) THEN
        IF( RSDX.GE.1.D+2 ) THEN
          WRITE(ISC,'(10X,A)') '---  Excessive Residual  ---'
          WRITE(IWR,'(10X,A)') '---  Excessive Residual  ---'
        ENDIF
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          NT = NSD(IEQT)
          IF( NT.GT.0 ) THEN
            NPT = NPHAZ(2,NT)
            NCHT = INDEX( PH_CND(NPT),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &        'Energy Equation Maximum Residual = ',RSD(IEQT),
     &        ': Node = ',NT,': I,J,K = ',ID(NT),',',JD(NT),',',KD(NT),
     &        ': Phase Condition = ',PH_CND(NPT)(1:NCHT)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &        'Energy Equation Maximum Residual = ',RSD(IEQT),
     &        ': Node = ',NT,': I,J,K = ',ID(NT),',',JD(NT),',',KD(NT),
     &        ': Phase Condition = ',PH_CND(NPT)(1:NCHT)
          ENDIF
        ENDIF
        NW = NSD(IEQW)
        IF( NW.GT.0 ) THEN
          NPW = NPHAZ(2,NW)
          NCHW = INDEX( PH_CND(NPW),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &      'Water Equation Maximum Residual = ',RSD(IEQW),
     &      ': Node = ',NW,': I,J,K = ',ID(NW),',',JD(NW),',',KD(NW),
     &      ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &      'Water Equation Maximum Residual = ',RSD(IEQW),
     &      ': Node = ',NW,': I,J,K = ',ID(NW),',',JD(NW),',',KD(NW),
     &      ': Phase Condition = ',PH_CND(NPW)(1:NCHW)
!
!---      Extended output on convergence failure  ---
!
          IF( ISLC(62).EQ.1 ) THEN
            WRITE(ISC,'(4X,2A)')
     &        'Current Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Current Aqueous Saturation = ',SL(2,NW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Current Trapped Gas Saturation = ',SGT(2,NW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Current Salt Saturation = ',SS(2,NW)
            WRITE(IWR,'(4X,2A)')
     &        'Current Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Current Aqueous Saturation = ',SL(2,NW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Current Trapped Gas Saturation = ',SGT(2,NW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Current Salt Saturation = ',SS(2,NW)
            NPW = NPHAZ(1,NW)
            NCHW = INDEX( PH_CND(NPW),'  ') - 1
            WRITE(ISC,'(4X,2A)')
     &        'Previous Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Previous Aqueous Saturation = ',SL(1,NW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Previous Trapped Gas Saturation = ',SGT(1,NW)
            WRITE(ISC,'(4X,A,1PE11.4)')
     &        'Previous Salt Saturation = ',SS(1,NW)
            WRITE(IWR,'(4X,2A)')
     &        'Previous Phase Condition = ',PH_CND(NPW)(1:NCHW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Previous Aqueous Saturation = ',SL(1,NW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Previous Trapped Gas Saturation = ',SGT(1,NW)
            WRITE(IWR,'(4X,A,1PE11.4)')
     &        'Previous Salt Saturation = ',SS(1,NW)
          ENDIF
        ENDIF
!
!---    Isoair option
!
        IF( ISLC(37).EQ.0 ) THEN
          NA = NSD(IEQA)
          IF( NA.GT.0 ) THEN
            NPA = NPHAZ(2,NA)
            NCHA = INDEX( PH_CND(NPA),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &        'Air Equation Maximum Residual = ',RSD(IEQA),
     &        ': Node = ',NA,': I,J,K = ',ID(NA),',',JD(NA),',',KD(NA),
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &        'Air Equation Maximum Residual = ',RSD(IEQA),
     &        ': Node = ',NA,': I,J,K = ',ID(NA),',',JD(NA),',',KD(NA),
     &        ': Phase Condition = ',PH_CND(NPA)(1:NCHA)
!
!---        Extended output on convergence failure  ---
!
            IF( ISLC(62).EQ.1 ) THEN
              WRITE(ISC,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL(2,NA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Trapped Gas Saturation = ',SGT(2,NA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS(2,NA)
              WRITE(IWR,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL(2,NA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Trapped Gas Saturation = ',SGT(2,NA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS(2,NA)
              NPA = NPHAZ(1,NA)
              NCHA = INDEX( PH_CND(NPA),'  ') - 1
              WRITE(ISC,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL(1,NA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Trapped Gas Saturation = ',SGT(1,NA)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS(1,NA)
              WRITE(IWR,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPA)(1:NCHA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL(1,NA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Trapped Gas Saturation = ',SGT(1,NA)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS(1,NA)
            ENDIF
          ENDIF
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          NS = NSD(IEQS)
          IF( NS.GT.0 ) THEN
            NPS = NPHAZ(2,NS)
            NCHS = INDEX( PH_CND(NPS),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD(IEQS),
     &        ': Node = ',NS,': I,J,K = ',ID(NS),',',JD(NS),',',KD(NS),
     &        ': Phase Condition = ',PH_CND(NPS)(1:NCHS)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6,A,I6,A,I6,2A)')
     &        'Salt Equation Maximum Residual = ',RSD(IEQS),
     &        ': Node = ',NS,': I,J,K = ',ID(NS),',',JD(NS),',',KD(NS),
     &        ': Phase Condition = ',PH_CND(NPS)(1:NCHS)
!
!---        Extended output on convergence failure  ---
!
            IF( ISLC(62).EQ.1 ) THEN
              WRITE(ISC,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL(2,NS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Trapped Gas Saturation = ',SGT(2,NS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS(2,NS)
              WRITE(IWR,'(4X,2A)')
     &          'Current Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Aqueous Saturation = ',SL(2,NS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Trapped Gas Saturation = ',SGT(2,NS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Current Salt Saturation = ',SS(2,NS)
              NPS = NPHAZ(1,NS)
              NCHS = INDEX( PH_CND(NPS),'  ') - 1
              WRITE(ISC,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL(1,NS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Trapped Gas Saturation = ',SGT(1,NS)
              WRITE(ISC,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS(1,NS)
              WRITE(IWR,'(4X,2A)')
     &          'Previous Phase Condition = ',PH_CND(NPS)(1:NCHS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Aqueous Saturation = ',SL(1,NS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Trapped Gas Saturation = ',SGT(1,NS)
              WRITE(IWR,'(4X,A,1PE11.4)')
     &          'Previous Salt Saturation = ',SS(1,NS)
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
!---    Number of time step reductions failure: stop simulation  ---
!
        ELSE
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit ' // 
     &      'Exceeded  ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit ' // 
     &      'Exceeded  ---'
          ICNV = 4
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SBND_GT( NSL )
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
!     Modify the Jacobian matrix for the solute transport equation
!     to incorporate boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, October 1995.
!






!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PORMED
      USE JACOB
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

      REAL*8 BCX(LSPBC+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SBND_GT'
      NBCT = MIN( NSL+LUK,NSOLU+LUK+1 )
!
!---  Loop over number of specified boundary conditions  ---
!
      DO NB = 1,NBC
!
!---    Zero flux boundary condition  ---
!
        IF( IBCT(NSL+LUK,NB).EQ.3 ) CYCLE
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) CYCLE
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) CYCLE
        IF( IBCM(NB).EQ.1 ) THEN
!
!---      Solute transport  ---
!
          IF( NSL.LE.NSOLU ) THEN
            BCX(1) = BC(NSL+LBCU,1,MB)
            IF( IBCT(NBCT,NB).EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---      Reactive species transport  ---
!
          ELSE
            BCX(1) = 0.D+0
            DO NSPX = 1,IBCSP(1,NB)
              MX = NSOLU+LBCU+NSPX
              BCX(NSPX+1) = BC(MX,1,MB)
              IF( IBCT(NBCT,NB).EQ.12 ) BCX(NSPX+1) = SP_CBO(NB,NSP)
            ENDDO
          ENDIF
        ELSE
          DO M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
              TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
              DTBC = MIN( BC(1,M,MB)-TMZ,DT )
              TFBC = (TMZ-5.D-1*DTBC-BC(1,M-1,MB))/TDBC
!
!---          Solute transport  ---
!
              IF( NSL.LE.NSOLU ) THEN
                BCX(1) = BC(NSL+LBCU,M-1,MB) +
     &            TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
                IF( IBCT(NBCT,NB).EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---          Reactive species transport  ---
!
              ELSE
                BCX(1) = 0.D+0
                DO NSPX = 1,IBCSP(1,NB)
                  NSP = IBCSP(NSPX+1,NB)
                  MX = NSOLU+LBCU+NSPX
                  BCX(NSPX+1) = BC(MX,M-1,MB) +
     &              TFBC*(BC(MX,M,MB)-BC(MX,M-1,MB))
                  IF( IBCT(NBCT,NB).EQ.12 ) BCX(NSPX+1) = SP_CBO(NB,NSP)
                ENDDO
              ENDIF
              GOTO 110
            ENDIF
          ENDDO
          CYCLE
        ENDIF
  110   CONTINUE
        N = IBCN(NB)
        IZN = IZ(N)
        MF = 1
        MP = IXP(N)
        I = ID(N)
        J = JD(N)
        K = KD(N)
        IF( ILES.EQ.1 ) THEN
          MCOL = MP
          MROW = MDT
        ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
          MA = 1
          MCOL = KLUC(MP,MA)
          MA = MA + 1




        ENDIF
!
!---    Diffusion coefficients at node adjacent to boundary  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SMDLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
!
!---    Kemper and van Schaik empirical model  ---
!
        IF( IEDL(NSL).EQ.2 ) THEN
          DLP = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &      EXP(SL(2,N)*PORD(2,N)*SDCL(3,IZN,NSL))
!
!---    Temperature dependent molecular diffusion coefficient  ---
!
        ELSEIF( IEDL(NSL).EQ.3 ) THEN
          DLP = TORL(2,N)*SL(2,N)*PORD(2,N)*SMDL(NSL)
!
!---    Constant molecular diffusion coefficient  ---
!
        ELSE
          DLP = TORL(2,N)*SL(2,N)*PORD(2,N)*SMDLP
        ENDIF
        PCOR = (PG(2,N)+PATM)/PATM
        SMDGP = SMDG(NSL)*(TCOR**1.75)/PCOR
        DGP = TORG(2,N)*SG(2,N)*PORD(2,N)*SMDGP
!
!---    Phase fraction factors at node adjacent to boundary  ---
!
        FCLP = 0.D+0
        IF( SL(2,N).GT.SMALL ) FCLP = YL(N,NSL)/(SL(2,N)*PORD(2,N))
        FCGP = 0.D+0
        IF( SG(2,N).GT.SMALL ) FCGP = YG(N,NSL)/(SG(2,N)*PORD(2,N))
!
!---    Phase fraction factors at boundary  ---
!
        IF( IPCL(NSL).EQ.2 ) THEN
          XVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))*SLB(2,NB)
        ELSE
          XVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
        ENDIF
        XVLB = SLB(2,NB)*PORDB(2,NB)
        XVGB = SGB(2,NB)*PORDB(2,NB)
!
!---    Constant gas-aqueous partition coefficient  ---
!
        IF( IPCGL(NSL).EQ.0 ) THEN
          PCGLX = PCGL(1,NSL)
!
!---    Temperature dependent gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.1 ) THEN
          TK = TB(2,NB)+TABS
          PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &      + PCGL(3,NSL)*LOG(TK)
     &      + PCGL(4,NSL)*TK + PCGL(5,NSL)*TK**2 )
!
!---    Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.2 ) THEN
          PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
        ENDIF
        PCGLX = MAX( PCGLX,1.D-20 )
        PCGLX = MIN( PCGLX,1.D+20 )
!
!---    Solute transport only, skip calculations for reactive
!       species transport  ---
!
        IF( NSL.LE.NSOLU ) THEN
!
!---      Phase-volumetric concentration ratios  ---
!
          FCL = 1.D+0/(XVSB + XVLB + XVGB*PCGLX)
          FCG = 1.D+0/((XVSB + XVLB)/PCGLX + XVGB)
!
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = XVLB*FCL
          YGB(NB,NSL) = XVGB*FCG
!
!---      Convert boundary phase concentrations to
!         volumetric concentrations  ---
!
          IF( IBCT(NSL+LUK,NB).EQ.8 .OR. IBCT(NSL+LUK,NB).EQ.14 ) THEN
            BCX(1) = BCX(1)/FCL
          ELSEIF( IBCT(NSL+LUK,NB).EQ.9 .OR. 
     &      IBCT(NSL+LUK,NB).EQ.15 ) THEN
            BCX(1) = BCX(1)/FCG
          ENDIF
          CB(NB,NSL) = BCX(1)

        ELSE
!
!---      Convert species concentrations to total-component
!         concentrations  ---
!
          IF( NSL.LE.NSOLU+NEQC ) THEN
            NEQ = NSL-NSOLU
            DO NSP = 1,IEQ_C(1,NEQ)
              DO NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_C(NSP+1,NEQ) ) THEN
                  BCX(1) = BCX(1) + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
              ENDDO
            ENDDO
!
!---      Convert species concentrations to total-kinetic
!         concentrations  ---
!
          ELSEIF( NSL.LE.NSOLU+NEQC+NEQK ) THEN
            NEQ = NSL-NSOLU-NEQC
            DO NSP = 1,IEQ_K(1,NEQ)
              DO NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_K(NSP+1,NEQ) ) THEN
                  BCX(1) = BCX(1) + EQ_K(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
!
!---      Phase-volumetric concentration ratios  ---
!
          FCL = 0.D+0
          IF( XVLB/EPSL.GT.EPSL ) FCL = 1.D+0/XVLB
!
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = 1.D+0
!
!---      Convert boundary phase concentrations to
!         volumetric concentrations  ---
!
          IF( IBCT(NBCT,NB).EQ.8 .OR. IBCT(NBCT,NB).EQ.14 ) THEN
            BCX(1) = BCX(1)/FCL
          ELSEIF( IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.15 ) THEN
            BCX(1) = BCX(1)/FCG
          ENDIF
!
!---      Load boundary concentration  ---
!
          CB(NB,NSL) = BCX(1)

        ENDIF
!
!---  Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
          NPZ = NSZ(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULBX,VLBX,WLBX,N,MF )
            CALL SHDP( WLBX,ULBX,VLBX,DISPL(IZN),DISPT(IZN),DPLB )
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGBX,VGBX,WGBX,N,MF )
            CALL SHDP( WGBX,UGBX,VGBX,DISPL(IZN),DISPT(IZN),DPGB )
          ELSE
            DPLB = 0.D+0
            DPGB = 0.D+0
          ENDIF
          FLB = AFZ(NPZ)*WL(1,NPZ)
          FGB = AFZ(NPZ)*WG(1,NPZ)
          CRLB = ABS( WL(1,NPZ) )*DT/(DZGF(N)*XVLB+SMALL)
          CRGB = ABS( WG(1,NPZ) )*DT/(DZGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8
     &       .OR. IBCT(NSL+LUK,NB).EQ.9  .OR. IBCT(NSL+LUK,NB).EQ.10
     &       .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            IF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SLB(2,NB)*PORDB(2,NB)*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDL(NSL)
            ELSE
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            ENDIF
            INDX = 16
            DLZ = DIFMN(DLB,DLP,DZGF(N),DZGF(N),WL(1,NPZ),INDX)
            DLZ = AFZ(NPZ)*(DLZ+DPLB)/(5.D-1*DZGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGZ = DIFMN(DGB,DGP,DZGF(N),DZGF(N),WG(1,NPZ),INDX)
            DGZ = AFZ(NPZ)*(DGZ+DPGB)/(5.D-1*DZGF(N))
!
!---        Patankar transport for the boundary surface  ---
!
            ALB = MAX( FLB,0.D+0 ) +
     &        DLZ*MAX((ONE-(TENTH*ABS(FLB)/(DLZ+SMALL)))**5,0.D+0)
            AGB = MAX( FGB,0.D+0 ) +
     &        DGZ*MAX((ONE-(TENTH*ABS(FGB)/(DGZ+SMALL)))**5,0.D+0)
            AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP
            AB = ALB*FCL + AGB*FCG
            BLU(MP) = BLU(MP) + AB*BCX(1)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLB = MIN( FLB,0.D+0 )
            FGB = MIN( FGB,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALB = MAX( FLB,0.D+0 )
            AGB = MAX( FGB,0.D+0 )
            AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP
            AB = ALB*FCL + AGB*FCG
            BLU(MP) = BLU(MP) + AB*BCX(1)
!
!---      Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 ) THEN
            FLB = MAX( FLB,0.D+0 )
            FGB = MAX( FGB,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALB = MAX( FLB,0.D+0 )
            AGB = MAX( FGB,0.D+0 )
            AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP
            AB = ALB*FCL + AGB*FCG
            BLU(MP) = BLU(MP) + AB*BCX(1)
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLB = MAX( FLB,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALB = MAX( FLB,0.D+0 )
            AP = (ALB-FLB)*FCLP
            AB = ALB*FCL
            BLU(MP) = BLU(MP) + AB*BCX(1)
!
!---      Inflow gas ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.15 ) THEN
            FGB = MAX( FGB,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            AGB = MAX( FGB,0.D+0 )
            AP = (AGB-FGB)*FCGP
            AB = AGB*FCG
            BLU(MP) = BLU(MP) + AB*BCX(1)
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          NPY = NSY(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULSX,VLSX,WLSX,N,MF )
            CALL SHDP( VLSX,WLSX,ULSX,DISPL(IZN),DISPT(IZN),DPLS )
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGSX,VGSX,WGSX,N,MF )
            CALL SHDP( VGSX,WGSX,UGSX,DISPL(IZN),DISPT(IZN),DPGS )
          ELSE
            DPLS = 0.D+0
            DPGS = 0.D+0
          ENDIF
          FLS = AFY(NPY)*VL(1,NPY)
          FGS = AFY(NPY)*VG(1,NPY)
          CRLS = ABS( VL(1,NPY) )*DT/(RP(I)*DYGF(N)*XVLB+SMALL)
          CRGS = ABS( VG(1,NPY) )*DT/(RP(I)*DYGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8
     &       .OR. IBCT(NSL+LUK,NB).EQ.9  .OR. IBCT(NSL+LUK,NB).EQ.10
     &       .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            IF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SLB(2,NB)*PORDB(2,NB)*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDL(NSL)
            ELSE
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            ENDIF
            INDX = 16
            DLY = DIFMN(DLB,DLP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
            DLY = AFY(NPY)*(DLY+DPLS)/RP(I)/(5.D-1*DYGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGY = DIFMN(DGB,DGP,DYGF(N),DYGF(N),VG(1,NPY),INDX)
            DGY = AFY(NPY)*(DGY+DPGS)/RP(I)/(5.D-1*DYGF(N))
!
!---        Patankar transport for the boundary surface  ---
!
            ALS = MAX( FLS,0.D+0 ) +
     &        DLY*MAX((ONE-(TENTH*ABS(FLS)/(DLY+SMALL)))**5,0.D+0)
            AGS = MAX( FGS,0.D+0 ) +
     &        DGY*MAX((ONE-(TENTH*ABS(FGS)/(DGY+SMALL)))**5,0.D+0)
            AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP
            AS = ALS*FCL + AGS*FCG
            BLU(MP) = BLU(MP) + AS*BCX(1)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLS = MIN( FLS,0.D+0 )
            FGS = MIN( FGS,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALS = MAX( FLS,0.D+0 )
            AGS = MAX( FGS,0.D+0 )
            AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP
            AS = ALS*FCL + AGS*FCG
            BLU(MP) = BLU(MP) + AS*BCX(1)
!
!---      Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 ) THEN
            FLS = MAX( FLS,0.D+0 )
            FGS = MAX( FGS,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALS = MAX( FLS,0.D+0 )
            AGS = MAX( FGS,0.D+0 )
            AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP
            AS = ALS*FCL + AGS*FCG
            BLU(MP) = BLU(MP) + AS*BCX(1)
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLS = MAX( FLS,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALS = MAX( FLS,0.D+0 )
            AP = (ALS-FLS)*FCLP
            AS = ALS*FCL
            BLU(MP) = BLU(MP) + AS*BCX(1)
!
!---      Inflow gas ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.15 ) THEN
            FGS = MAX( FGS,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            AGS = MAX( FGS,0.D+0 )
            AP = (AGS-FGS)*FCGP
            AS = AGS*FCG
            BLU(MP) = BLU(MP) + AS*BCX(1)
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          NPX = NSX(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULWX,VLWX,WLWX,N,MF )
            CALL SHDP( ULWX,VLWX,WLWX,DISPL(IZN),DISPT(IZN),DPLW )
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGWX,VGWX,WGWX,N,MF )
            CALL SHDP( UGWX,VGWX,WGWX,DISPL(IZN),DISPT(IZN),DPGW )
          ELSE
            DPLW = 0.D+0
            DPGW = 0.D+0
          ENDIF
          FLW = AFX(NPX)*UL(1,NPX)
          FGW = AFX(NPX)*UG(1,NPX)
          CRLW = ABS( UL(1,NPX) )*DT/(DXGF(N)*XVLB+SMALL)
          CRGW = ABS( UG(1,NPX) )*DT/(DXGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8
     &       .OR. IBCT(NSL+LUK,NB).EQ.9  .OR. IBCT(NSL+LUK,NB).EQ.10
     &       .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            IF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SLB(2,NB)*PORDB(2,NB)*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDL(NSL)
            ELSE
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            ENDIF
            INDX = 16
            DLX = DIFMN(DLB,DLP,DXGF(N),DXGF(N),UL(1,NPX),INDX)
            DLX = AFX(NPX)*(DLX+DPLW)/(5.D-1*DXGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGX = DIFMN(DGB,DGP,DXGF(N),DXGF(N),UG(1,NPX),INDX)
            DGX = AFX(NPX)*(DGX+DPGW)/(5.D-1*DXGF(N))
!
!---        Patankar transport for the boundary surface  ---
!
            ALW = MAX(FLW,0.D+0)
     &        + DLX*MAX((ONE-(TENTH*ABS(FLW)/(DLX+SMALL)))**5,0.D+0)
            AGW = MAX(FGW,0.D+0)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGW)/(DGX+SMALL)))**5,0.D+0)
            AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
            AW = ALW*FCL + AGW*FCG
            BLU(MP) = BLU(MP) + AW*BCX(1)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLW = MIN( FLW,0.D+0 )
            FGW = MIN( FGW,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALW = MAX(FLW,0.D+0)
            AGW = MAX(FGW,0.D+0)
            AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
            AW = ALW*FCL + AGW*FCG
            BLU(MP) = BLU(MP) + AW*BCX(1)
!
!---      Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 ) THEN
            FLW = MAX( FLW,0.D+0 )
            FGW = MAX( FGW,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALW = MAX(FLW,0.D+0)
            AGW = MAX(FGW,0.D+0)
            AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
            AW = ALW*FCL + AGW*FCG
            BLU(MP) = BLU(MP) + AW*BCX(1)
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLW = MAX( FLW,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALW = MAX(FLW,0.D+0)
            AP = (ALW-FLW)*FCLP
            AW = ALW*FCL
            BLU(MP) = BLU(MP) + AW*BCX(1)
!
!---      Inflow gas ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.15 ) THEN
            FGW = MAX( FGW,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            AGW = MAX(FGW,0.D+0)
            AP = (AGW-FGW)*FCGP
            AW = AGW*FCG
            BLU(MP) = BLU(MP) + AW*BCX(1)
          ENDIF
!
!---    East boundary
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          NQX = NSX(N) + 1
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULEX,VLEX,WLEX,N,MF )
            CALL SHDP( ULEX,VLEX,WLEX,DISPL(IZN),DISPT(IZN),DPLE )
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGEX,VGEX,WGEX,N,MF )
            CALL SHDP( UGEX,VGEX,WGEX,DISPL(IZN),DISPT(IZN),DPGE )
          ELSE
            DPLE = 0.D+0
            DPGE = 0.D+0
          ENDIF
          FLE = AFX(NQX)*UL(1,NQX)
          FGE = AFX(NQX)*UG(1,NQX)
          CRLE = ABS( UL(1,NQX) )*DT/(DXGF(N)*XVLB+SMALL)
          CRGE = ABS( UG(1,NQX) )*DT/(DXGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.4 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            IF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SLB(2,NB)*PORDB(2,NB)*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDL(NSL)
            ELSE
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            ENDIF
            INDX = 16
            DLX = DIFMN(DLP,DLB,DXGF(N),DXGF(N),UL(1,NQX),INDX)
            DLX = AFX(NQX)*(DLX+DPLE)/(5.D-1*DXGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGX = DIFMN(DGP,DGB,DXGF(N),DXGF(N),UG(1,NQX),INDX)
            DGX = AFX(NQX)*(DGX+DPGE)/(5.D-1*DXGF(N))
!
!---        Patankar transport for the boundary surface  ---
!
            ALE = MAX( -FLE,0.D+0 ) +
     &        DLX*MAX((ONE-(TENTH*ABS(FLE)/(DLX+SMALL)))**5,0.D+0)
            AGE = MAX( -FGE,0.D+0 ) +
     &        DGX*MAX((ONE-(TENTH*ABS(FGE)/(DGX+SMALL)))**5,0.D+0)
            AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP
            AE = ALE*FCL + AGE*FCG
            BLU(MP) = BLU(MP) + AE*BCX(1)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLE = MAX( FLE,0.D+0 )
            FGE = MAX( FGE,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALE = MAX( -FLE,0.D+0 )
            AGE = MAX( -FGE,0.D+0 )
            AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP
            AE = ALE*FCL + AGE*FCG
            BLU(MP) = BLU(MP) + AE*BCX(1)
!
!---      Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 ) THEN
            FLE = MIN( FLE,0.D+0 )
            FGE = MIN( FGE,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALE = MAX( -FLE,0.D+0 )
            AGE = MAX( -FGE,0.D+0 )
            AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP
            AE = ALE*FCL + AGE*FCG
            BLU(MP) = BLU(MP) + AE*BCX(1)
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLE = MIN( FLE,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALE = MAX( -FLE,0.D+0 )
            AP = (ALE+FLE)*FCLP
            AE = ALE*FCL
            BLU(MP) = BLU(MP) + AE*BCX(1)
!
!---      Inflow gas ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.15 ) THEN
            FGE = MIN( FGE,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            AGE = MAX( -FGE,0.D+0 )
            AP = (AGE+FGE)*FCGP
            AE = AGE*FCG
            BLU(MP) = BLU(MP) + AE*BCX(1)
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          NQY = NSY(N) + IFLD
!
!---      Hydraulic dispersion  ---
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULNX,VLNX,WLNX,N,MF )
            CALL SHDP( VLNX,WLNX,ULNX,DISPL(IZN),DISPT(IZN),DPLN )
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGNX,VGNX,WGNX,N,MF )
            CALL SHDP( VGNX,WGNX,UGNX,DISPL(IZN),DISPT(IZN),DPGN )
          ELSE
            DPLN = 0.D+0
            DPGN = 0.D+0
          ENDIF
          FLN = AFY(NQY)*VL(1,NQY)
          FGN = AFY(NQY)*VG(1,NQY)
          CRLN = ABS( VL(1,NQY) )*DT/(RP(I)*DYGF(N)*XVLB+SMALL)
          CRGN = ABS( VG(1,NQY) )*DT/(RP(I)*DYGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8
     &       .OR. IBCT(NSL+LUK,NB).EQ.9  .OR. IBCT(NSL+LUK,NB).EQ.10
     &       .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            IF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SLB(2,NB)*PORDB(2,NB)*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDL(NSL)
            ELSE
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            ENDIF
            INDX = 16
            DLY = DIFMN(DLP,DLB,DYGF(N),DYGF(N),VL(1,NQY),INDX)
            DLY = AFY(NQY)*(DLY+DPLN)/RP(I)/(5.D-1*DYGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGY = DIFMN(DGP,DGB,DYGF(N),DYGF(N),VG(1,NQY),INDX)
            DGY = AFY(NQY)*(DGY+DPGN)/RP(I)/(5.D-1*DYGF(N))
!
!---        Patankar transport for the boundary surface  ---
!
            ALN = MAX( -FLN,0.D+0 ) +
     &        DLY*MAX((ONE-(TENTH*ABS(FLN)/(DLY+SMALL)))**5,0.D+0)
            AGN = MAX( -FGN,0.D+0 ) +
     &        DGY*MAX((ONE-(TENTH*ABS(FGN)/(DGY+SMALL)))**5,0.D+0)
            AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP
            AN = ALN*FCL + AGN*FCG
            BLU(MP) = BLU(MP) + AN*BCX(1)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLN = MAX( FLN,0.D+0 )
            FGN = MAX( FGN,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALN = MAX( -FLN,0.D+0 )
            AGN = MAX( -FGN,0.D+0 )
            AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP
            AN = ALN*FCL + AGN*FCG
            BLU(MP) = BLU(MP) + AN*BCX(1)
!
!---      Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 ) THEN
            FLN = MIN( FLN,0.D+0 )
            FGN = MIN( FGN,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALN = MAX( -FLN,0.D+0 )
            AGN = MAX( -FGN,0.D+0 )
            AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP
            AN = ALN*FCL + AGN*FCG
            BLU(MP) = BLU(MP) + AN*BCX(1)
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLN = MIN( FLN,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALN = MAX( -FLN,0.D+0 )
            AP = (ALN+FLN)*FCLP
            AN = ALN*FCL
            BLU(MP) = BLU(MP) + AN*BCX(1)
!
!---      Inflow gas ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.15 ) THEN
            FGN = MIN( FGN,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            AGN = MAX( -FGN,0.D+0 )
            AP = (AGN+FGN)*FCGP
            AN = AGN*FCG
            BLU(MP) = BLU(MP) + AN*BCX(1)
          ENDIF
!
!---    Top boundary
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          NQZ = NSZ(N) + IJFLD
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULTX,VLTX,WLTX,N,MF )
            CALL SHDP( WLTX,ULTX,VLTX,DISPL(IZN),DISPT(IZN),DPLT )
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGTX,VGTX,WGTX,N,MF )
            CALL SHDP( WGTX,UGTX,VGTX,DISPL(IZN),DISPT(IZN),DPGT )
          ELSE
            DPLT = 0.D+0
            DPGT = 0.D+0
          ENDIF
          FLT = AFZ(NQZ)*WL(1,NQZ)
          FGT = AFZ(NQZ)*WG(1,NQZ)
          CRLT = ABS( WL(1,NQZ) )*DT/(DZGF(N)*XVLB+SMALL)
          CRGT = ABS( WG(1,NQZ) )*DT/(DZGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8
     &       .OR. IBCT(NSL+LUK,NB).EQ.9  .OR. IBCT(NSL+LUK,NB).EQ.10
     &       .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            IF( IEDL(NSL).EQ.2 ) THEN
              DLB = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &          EXP(SLB(2,NB)*PORDB(2,NB)*SDCL(3,IZN,NSL))
            ELSEIF( IEDL(NSL).EQ.3 ) THEN
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDL(NSL)
            ELSE
              DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            ENDIF
            INDX = 16
            DLZ = DIFMN(DLP,DLB,DZGF(N),DZGF(N),WL(1,NQZ),INDX)
            DLZ = AFZ(NQZ)*(DLZ+DPLT)/(5.D-1*DZGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGZ = DIFMN(DGP,DGB,DZGF(N),DZGF(N),WG(1,NQZ),INDX)
            DGZ = AFZ(NQZ)*(DGZ+DPGT)/(5.D-1*DZGF(N))
!
!---        Patankar transport for the boundary surface  ---
!
            ALT = MAX( -FLT,0.D+0 ) +
     &        DLZ*MAX((ONE-(TENTH*ABS(FLT)/(DLZ+SMALL)))**5,0.D+0)
            AGT = MAX( -FGT,0.D+0 ) +
     &        DGZ*MAX((ONE-(TENTH*ABS(FGT)/(DGZ+SMALL)))**5,0.D+0)
            AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP
            AT = ALT*FCL + AGT*FCG
            BLU(MP) = BLU(MP) + AT*BCX(1)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLT = MAX( FLT,0.D+0 )
            FGT = MAX( FGT,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALT = MAX( -FLT,0.D+0 )
            AGT = MAX( -FGT,0.D+0 )
            AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP
            AT = ALT*FCL + AGT*FCG
            BLU(MP) = BLU(MP) + AT*BCX(1)
!
!---      Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 ) THEN
            FLT = MIN( FLT,0.D+0 )
            FGT = MIN( FGT,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALT = MAX( -FLT,0.D+0 )
            AGT = MAX( -FGT,0.D+0 )
            AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP
            AT = ALT*FCL + AGT*FCG
            BLU(MP) = BLU(MP) + AT*BCX(1)
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLT = MIN( FLT,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            ALT = MAX( -FLT,0.D+0 )
            AP = (ALT+FLT)*FCLP
            AT = ALT*FCL
            BLU(MP) = BLU(MP) + AT*BCX(1)
!
!---      Inflow gas ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.15 ) THEN
            FGT = MIN( FGT,0.D+0 )
!
!---        Patankar transport for the boundary surface  ---
!
            AGT = MAX( -FGT,0.D+0 )
            AP = (AGT+FGT)*FCGP
            AT = AGT*FCG
            BLU(MP) = BLU(MP) + AT*BCX(1)
          ENDIF
        ENDIF
        IF( ILES.EQ.1 ) THEN
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + AP
        ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
          DLU(MCOL) = DLU(MCOL) + AP




        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SBND_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_GT( ASL_F,ASL_M,ASLX,ASLMINX,ESLX,ESGTX,ESGTMX,
     &  PGX,PLX,BTGLX,SGX,SGTX,SLX,SLFX,SLMX,SLRX,INDX,IZN )
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
!     Aqueous and gas saturations.
!
!     INDX = 0 : Trapped-gas saturation computed.
!     INDX = 1 : Trapped-gas saturation given.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 June 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/SP_GT'
!
!---  Zero scaling factor  ---
!
      IF( BTGLX.LT.EPSL ) THEN
        SGX = 1.D+0
        SGTX = 0.D+0
        SLX = 0.D+0
        ESLX = 0.D+0
        ESGTX = 0.D+0
        ASLX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  van Genuchten saturation function w/o gas entrapment,
!     w/ or w/o Webb extension  ---
!
      IF( ISCHR(IZN).EQ.1 ) THEN
!
!---    van Genuchten saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(9,IZN)
          SLRX = SCHR(4,IZN)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(8,IZN)
            HDGL = MIN( HDGL,SCHR(12,IZN) )
            DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
            SLX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPX
            ASLX = SLX
!
!---      Capillary head at or below the matching point head,
!         use van Genuchten function
!
          ELSE
            CN = MAX( SCHR(3,IZN),SMALL )
            IF( SCHR(14,IZN).LE.ZERO ) THEN
              IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                CM = 1.D+0 - 2.D+0/CN
              ELSE
                CM = 1.D+0 - 1.D+0/CN
              ENDIF
            ELSE
              CM = SCHR(14,IZN)
            ENDIF
            ASLX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGL)**CN))**CM
            SLX = ASLX*(1.D+0-SLRX) + SLRX
            ASLX = SLX
          ENDIF
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          ESGTX = 0.D+0
          SGTX = 0.D+0
!
!---    van Genuchten saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL  ------
!
        ELSE
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          CN = MAX( SCHR(3,IZN),SMALL )
          IF( SCHR(14,IZN).LE.ZERO ) THEN
            IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(14,IZN)
          ENDIF
          ASLX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGL)**CN))**CM
          SLRX = SCHR(4,IZN)
          SLX = ASLX*(1.D+0-SLRX) + SLRX
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          ESGTX = 0.D+0
          SGTX = 0.D+0
        ENDIF
!
!---  Brooks and Corey saturation function w/o gas entrapment,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 ) THEN
!
!---    Brooks and Corey saturation function w/o gas entrapment, 
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(9,IZN)
          SLRX = SCHR(4,IZN)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(8,IZN)
            HDGL = MIN( HDGL,SCHR(12,IZN) )
            DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
            SLX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPX
            ASLX = SLX
!
!---      Capillary head at or below the matching point head,
!         use Brooks-Corey function
!
          ELSE
            CL = MAX( SCHR(3,IZN),SMALL )
            IF( HDGL.LE.(SCHR(1,IZN)+1.D-9) ) THEN
              ASLX = 1.D+0
            ELSE
              ASLX = (SCHR(1,IZN)/HDGL)**CL
            ENDIF
            SLRX = SCHR(4,IZN)
            SLX = ASLX*(1.D+0-SLRX) + SLRX
            ASLX = SLX
          ENDIF
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          ESGTX = 0.D+0
          SGTX = 0.D+0
!
!---    Brooks and Corey saturation function w/o gas entrapment,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          CL = MAX( SCHR(3,IZN),SMALL )
          IF( HDGL.LE.(SCHR(1,IZN)+1.D-9) ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (SCHR(1,IZN)/HDGL)**CL
          ENDIF
          SLRX = SCHR(4,IZN)
          SGRX = SCHR(15,IZN)
          SLX = ASLX*(1.D+0-SLRX-SGRX) + SLRX
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          ESGTX = 0.D+0
          SGTX = 0.D+0
        ENDIF
!
!---  Dual porosity van Genuchten saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.3 ) THEN
!
!---    Matrix van Genuchten saturation function  
!       w/o gas entrapment, w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(9,IZN)
          SLR_M = SCHR(4,IZN)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(8,IZN)
            HDGL = MIN( HDGL,SCHR(12,IZN) )
            DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
            SLMX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPX
            ASL_M = SLMX
!
!---      Capillary head at or below the matching point head,
!         use van Genuchten function
!
          ELSE
            CN = MAX( SCHR(3,IZN),SMALL )
            IF( SCHR(14,IZN).LE.ZERO ) THEN
              IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                CM = 1.D+0 - 2.D+0/CN
              ELSE
                CM = 1.D+0 - 1.D+0/CN
              ENDIF
            ELSE
              CM = SCHR(14,IZN)
            ENDIF
            ASL_M = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGL)**CN))**CM
            SLMX = ASL_M*(1.D+0-SLR_M) + SLR_M
            ASL_M = SLMX
          ENDIF
!
!---    Matrix van Genuchten saturation function 
!       w/o gas entrapment, w/o Webb extension ASL = ESL  ------
!
        ELSE
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          CN = MAX( SCHR(3,IZN),SMALL )
          IF( SCHR(14,IZN).LE.ZERO ) THEN
            IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(14,IZN)
          ENDIF
          ASL_M = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGL)**CN))**CM
          SLR_M = SCHR(4,IZN)
          SLMX = ASL_M*(1.D+0-SLR_M) + SLR_M
        ENDIF
!
!---    Fracture van Genuchten saturation function  
!       w/o gas entrapment, w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(11,IZN)
          SLR_F = SCHR(7,IZN)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(10,IZN)
            HDGL = MIN( HDGL,SCHR(12,IZN) )
            DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
            SLFX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPX
            ASL_F = SLFX
!
!---      Capillary head at or below the matching point head,
!         use van Genuchten function
!
          ELSE
            CN = MAX( SCHR(6,IZN),SMALL )
            IF( SCHR(15,IZN).LE.ZERO ) THEN
              IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                CM = 1.D+0 - 2.D+0/CN
              ELSE
                CM = 1.D+0 - 1.D+0/CN
              ENDIF
            ELSE
              CM = SCHR(15,IZN)
            ENDIF
            ASL_F = (1.D+0/(1.D+0 + (SCHR(5,IZN)*HDGL)**CN))**CM
            SLFX = ASL_F*(1.D+0-SLR_F) + SLR_F
            ASL_F = SLFX
          ENDIF
!
!---    Fracture van Genuchten saturation function 
!       w/o gas entrapment, w/o Webb extension ASL = ESL  ------
!
        ELSE
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          CN = MAX( SCHR(6,IZN),SMALL )
          IF( SCHR(15,IZN).LE.ZERO ) THEN
            IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(15,IZN)
          ENDIF
          ASL_F = (1.D+0/(1.D+0 + (SCHR(5,IZN)*HDGL)**CN))**CM
          SLR_F = SCHR(7,IZN)
          SLFX = ASL_F*(1.D+0-SLR_F) + SLR_F
        ENDIF
!
!---    Combination van Genuchten saturation function 
!       w/o gas entrapment, w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          PORD_MX = (1.D+0-POR(4,IZN))*POR(2,IZN)/
     &      ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
          PORD_FX = POR(4,IZN)/
     &      ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
          SLX = SLFX*PORD_FX + SLMX*PORD_MX
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
          ASLX = SLX
          ESLX = MIN( MAX( (SLFX-SLR_F)/(1.D+0-SLR_F),0.D+0 ),1.D+0 )
     &      *PORD_FX + MIN( MAX( (SLMX-SLR_M)/(1.D+0-SLR_M),0.D+0 ),
     &      1.D+0 )*PORD_MX
          ESGTX = 0.D+0
          SGTX = 0.D+0
!
!---    Combination van Genuchten saturation function
!       w/o gas entrapment, w/o Webb extension ASL = ESL  ------
!
        ELSE
          PORD_MX = (1.D+0-POR(4,IZN))*POR(2,IZN)/
     &      ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
          PORD_FX = POR(4,IZN)/
     &      ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
          SLX = SLFX*PORD_FX + SLMX*PORD_MX
          SLRX = SLR_F*PORD_FX + SLR_M*PORD_MX
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
          ASLX = ASL_F*PORD_FX + ASL_M*PORD_MX
          ESLX = MIN( MAX( (SLFX-SLR_F)/(1.D+0-SLR_F),0.D+0 ),1.D+0 )
     &      *PORD_FX + MIN( MAX( (SLMX-SLR_M)/(1.D+0-SLR_M),0.D+0 ),
     &      1.D+0 )*PORD_MX
          ESGTX = 0.D+0
          SGTX = 0.D+0
        ENDIF
!
!---  Dual porosity Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.4 ) THEN
!
!---    Matrix Brooks and Corey saturation function  
!       w/o gas entrapment, w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(9,IZN)
          SLR_M = SCHR(4,IZN)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(8,IZN)
            HDGL = MIN( HDGL,SCHR(12,IZN) )
            DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
            SLMX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPX
            ASL_M = SLMX
!
!---      Capillary head at or below the matching point head,
!         use Brooks and Corey function
!
          ELSE
            CL = MAX( SCHR(3,IZN),SMALL )
            IF( HDGL.LE.(SCHR(1,IZN)+1.D-9) ) THEN
              ASL_M = 1.D+0
            ELSE
              ASL_M = (SCHR(1,IZN)/HDGL)**CL
            ENDIF
            SLR_M = SCHR(4,IZN)
            SLMX = ASL_M*(1.D+0-SLR_M) + SLR_M
            ASL_M = SLMX
          ENDIF
!
!---    Matrix Brooks and Corey saturation function 
!       w/o gas entrapment, w/o Webb extension ASL = ESL  ------
!
        ELSE
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          CL = MAX( SCHR(3,IZN),SMALL )
          IF( HDGL.LE.(SCHR(1,IZN)+1.D-9) ) THEN
            ASL_M = 1.D+0
          ELSE
            ASL_M = (SCHR(1,IZN)/HDGL)**CL
          ENDIF
          SLR_M = SCHR(4,IZN)
          SGR_M = SCHR(14,IZN)
          SLMX = ASL_M*(1.D+0-SLR_M-SGR_M) + SLR_M
        ENDIF
!
!---    Fracture Brooks and Corey saturation function  
!       w/o gas entrapment, w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(11,IZN)
          SLR_F = SCHR(7,IZN)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(10,IZN)
            HDGL = MIN( HDGL,SCHR(12,IZN) )
            DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
            SLFX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPX
            ASL_F = SLFX
!
!---      Capillary head at or below the matching point head,
!         use Brooks and Corey function
!
          ELSE
            CL = MAX( SCHR(6,IZN),SMALL )
            IF( HDGL.LE.SCHR(5,IZN) ) THEN
              ASL_F = 1.D+0
            ELSE
              ASL_F = (SCHR(5,IZN)/HDGL)**CL
            ENDIF
            SLR_F = SCHR(7,IZN)
            SLFX = ASL_F*(1.D+0-SLR_F) + SLR_F
            ASL_F = SLFX
          ENDIF
!
!---    Fracture Brooks and Corey saturation function 
!       w/o gas entrapment, w/o Webb extension ASL = ESL  ------
!
        ELSE
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          CL = MAX( SCHR(6,IZN),SMALL )
          IF( HDGL.LE.SCHR(5,IZN) ) THEN
            ASL_F = 1.D+0
          ELSE
            ASL_F = (SCHR(5,IZN)/HDGL)**CL
          ENDIF
          SLR_F = SCHR(7,IZN)
          SGR_F = SCHR(15,IZN)
          SLFX = ASL_F*(1.D+0-SLR_F-SGR_F) + SLR_F
        ENDIF
!
!---    Combination Brooks and Corey saturation function 
!       w/o gas entrapment, w/ Webb extension ASL = SL  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          PORD_MX = (1.D+0-POR(4,IZN))*POR(2,IZN)/
     &      ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
          PORD_FX = POR(4,IZN)/
     &      ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
          SLX = SLFX*PORD_FX + SLMX*PORD_MX
          SLRX = SLR_F*PORD_FX + SLR_M*PORD_MX
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
          ASLX = SLX
          ESLX = MIN( MAX( (SLFX-SLR_F)/(1.D+0-SLR_F),0.D+0 ),1.D+0 )
     &      *PORD_FX + MIN( MAX( (SLMX-SLR_M)/(1.D+0-SLR_M),0.D+0 ),
     &      1.D+0 )*PORD_MX
          ESGTX = 0.D+0
          SGTX = 0.D+0
!
!---    Combination Brooks and Corey saturation function
!       w/o gas entrapment, w/o Webb extension ASL = ESL  ------
!
        ELSE
          PORD_MX = (1.D+0-POR(4,IZN))*POR(2,IZN)/
     &      ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
          PORD_FX = POR(4,IZN)/
     &      ( POR(4,IZN) + (1.D+0-POR(4,IZN))*POR(2,IZN) + SMALL )
          SLX = SLFX*PORD_FX + SLMX*PORD_MX
          SLRX = SLR_F*PORD_FX + SLR_M*PORD_MX
          SGRX = SGR_F*PORD_FX + SGR_M*PORD_MX
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
          ASLX = ASL_F*PORD_FX + ASL_M*PORD_MX
          ESLX = MIN( MAX( (SLFX-SLR_F)/(1.D+0-SLR_F),0.D+0 ),1.D+0 )
     &      *PORD_FX + MIN( MAX( (SLMX-SLR_M)/(1.D+0-SLR_M),0.D+0 ),
     &      1.D+0 )*PORD_MX
          ESGTX = 0.D+0
          SGTX = 0.D+0
        ENDIF
!
!---  Haverkamp saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.5 .OR. ISCHR(IZN).EQ.12 ) THEN
        HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
        IF( HDGL.LE.(SCHR(1,IZN)+1.D-9) ) THEN
          ASLX = 1.D+0
        ELSE
          ALPHAX = SCHR(2,IZN)/SCHR(5,IZN)
          IF( ISCHR(IZN).EQ.12 ) THEN
            HDGLX = LOG(HDGL/SCHR(5,IZN))
          ELSE
            HDGLX = HDGL/SCHR(5,IZN)
          ENDIF
          ASLX = ALPHAX/(ALPHAX + (HDGLX**SCHR(3,IZN)))
        ENDIF
        SLRX = SCHR(4,IZN)
        SLX = ASLX*(1.D+0-SLRX) + SLRX
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ESGTX = 0.D+0
        SGTX = 0.D+0
!
!---  Rossi-Nimmo-Sum saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.40 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        PSIOX = SCHR(1,IZN)
        PSIIX = SCHR(2,IZN)
        PSIDX = SCHR(4,IZN)
        CCX = SCHR(5,IZN)
        CLX = SCHR(6,IZN)
        CAX = SCHR(7,IZN)
        IF( HDGL.LE.PSIIX ) THEN
          ASLX = 1.D+0 - CCX*((HDGL/PSIOX)**2)
        ELSEIF( HDGL.LE.PSIDX ) THEN
          ASLX = (((PSIOX/HDGL)**CLX)-((PSIOX/PSIDX)**CLX))
     &      + CAX*LOG(PSIDX/HDGL)
        ELSEIF( HDGL.GT.PSIDX ) THEN
          ASLX = 0.D+0
        ELSE
          ASLX = 1.D+0
        ENDIF
        SLRX = 0.D+0
        SLX = ASLX
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ESGTX = 0.D+0
        SGTX = 0.D+0
!
!---  Rossi-Nimmo-Junction saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.41 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
        PSIOX = SCHR(1,IZN)
        PSIIX = SCHR(2,IZN)
        PSIJX = SCHR(3,IZN)
        PSIDX = SCHR(4,IZN)
        CCX = SCHR(5,IZN)
        CLX = SCHR(6,IZN)
        CAX = SCHR(7,IZN)
        IF( HDGL.LE.PSIIX ) THEN
          ASLX = 1.D+0 - CCX*((HDGL/PSIOX)**2)
        ELSEIF( HDGL.LE.PSIJX ) THEN
          ASLX = (PSIOX/HDGL)**CLX
        ELSEIF( HDGL.LE.PSIDX ) THEN
          ASLX = CAX*LOG(PSIDX/HDGL)
        ELSEIF( HDGL.GT.PSIDX ) THEN
          ASLX = 0.D+0
        ELSE
          ASLX = 1.D+0
        ENDIF
        SLRX = 0.D+0
        SLX = ASLX
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ESGTX = 0.D+0
        SGTX = 0.D+0
!
!---  Linear or linear-log interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.10 .OR. ISCHR(IZN).EQ.11 ) THEN
        HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
        IF( ISCHR(IZN).EQ.11 ) HDGL = LOG(HDGL)
        ITBX = 0
        ASLX = FNTBLY( HDGL,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        SLX = ASLX
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ESGTX = 0.D+0
        SGTX = 0.D+0
!
!---  Drainage-imbibition linear or linear-log 
!     interpolation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.12 .OR. ISCHR(IZN).EQ.13 ) THEN
        HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
        IF( ISCHR(IZN).EQ.13 ) HDGL = LOG(HDGL)
        ITBX = 0
        SLDX = FNTBLY( HDGL,ISLTBL(1,IZN),ISLTBL(2,IZN),ITBX )
        ITBX = 0
        SLIX = FNTBLY( HDGL,ISLTBL(3,IZN),ISLTBL(4,IZN),ITBX )
        ASLM = MAX( MIN( SLDX,ASLMINX ),0.D+0 )
        IF( SLDX.GE.SLIX ) THEN
          ASLX = 5.D-1*(SLIX + SQRT((SLIX**2) + 4.D+0*ASLM*SLDX -
     &      4.D+0*SLIX*ASLM))
        ELSE
          ASLX = 5.D-1*(SLIX - SQRT((SLIX**2) + 4.D+0*ASLM*SLDX -
     &      4.D+0*SLIX*ASLM))
        ENDIF
!
!---  Polynomial function  ---
!
      ELSEIF( ISCHR(IZN).EQ.19 ) THEN
        HDGL = MAX( (PGX-PLX)/RHORL/GRAV,1.D-14 )
!
!---    Convert head units for polynomial function basis  ---
!
        HDGLU = HDGL/SCHR(1,IZN)
        IF( HDGLU.LT.CPLY_SL(1,1,IZN) ) THEN
          ASLX = 1.D+0
        ELSEIF( HDGLU.GE.CPLY_SL(2,NPLY_SL(IZN),IZN) ) THEN
          ASLX = 0.D+0
        ELSE
          DO 192 NP = 1,NPLY_SL(IZN)
            IF( HDGLU.GE.CPLY_SL(1,NP,IZN) .AND.
     &        HDGLU.LT.CPLY_SL(2,NP,IZN) ) THEN
              ASLX = 0.D+0
              NPOLYC = LPOLYC
              DO 190 NC = 5,NPOLYC
                ASLX = ASLX + CPLY_SL(NC,NP,IZN)*(LOG10(HDGLU)**(NC-5))
  190         CONTINUE
              GOTO 194
            ENDIF
  192     CONTINUE
  194     CONTINUE
        ENDIF
        SLX = ASLX
        SGX = 1.D+0-SLX
        IF( SGX.LT.EPSL ) SGX = 0.D+0
        ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        ESGTX = 0.D+0
        SGTX = 0.D+0
!
!---  van Genuchten saturation function w/ gas entrapment,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.101 ) THEN
!
!---    van Genuchten saturation function w/ gas entrapment, 
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(9,IZN)
          SLRX = SCHR(4,IZN)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(8,IZN)
            HDGL = MIN( HDGL,SCHR(12,IZN) )
            DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
            SLX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPX
            ASLX = SLX
!
!---      Capillary head at or below the matching point head,
!         use van Genuchten function
!
          ELSE
            CN = MAX( SCHR(3,IZN),SMALL )
            IF( SCHR(14,IZN).LE.ZERO ) THEN
              IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                CM = 1.D+0 - 2.D+0/CN
              ELSE
                CM = 1.D+0 - 1.D+0/CN
              ENDIF
            ELSE
              CM = SCHR(14,IZN)
            ENDIF
            ASLX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGL)**CN))**CM
            SLX = ASLX*(1.D+0-SLRX) + SLRX
            ASLX = SLX
          ENDIF
!
!---      Compute trapped gas saturation, using the minimum
!         apparent aqueous saturation  ---
!
          IF( INDX.EQ.0 ) THEN
            ASLM = MAX( MIN( ASLX,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
              SGTMX = ESGTMX
              R = 1.D+0/SGTMX - 1.D+0
              ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
              IF( ESGTX.LT.EPSL ) ESGTX = 0.D+0
            ELSE
              ESGTX = 0.D+0
            ENDIF
!
!---      Trapped-gas saturation given, compute the minimum
!         apparent aqueous saturation  ---
!
          ELSEIF( INDX.EQ.1 ) THEN
            IF( ESGTMX.GT.EPSL ) THEN
              ESGTX = SGX
              SGTMX = ESGTMX
              ASLM = (SGTMX - SGTX)/(SGTMX + SGTMX*SGTX - SGTX)
              ASLMINX = MAX( MIN( ASLM,ASLX ),0.D+0 )
            ELSE
              ESGTX = 0.D+0
              ASLMINX = MAX( ASLX,0.D+0 )
            ENDIF
!
!---      Unrecognized option  ---
!
          ELSE
            CHMSG = 'Unrecognized Trapped-Gas Option'
            INDX = 12
            IMSG = INDX
            CALL WRMSGS( INDX )
          ENDIF
          SLX = ASLX - ESGTX
          SGTX = ESGTX
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
!
!---    van Genuchten saturation function w/ gas entrapment, 
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          CN = MAX( SCHR(3,IZN),SMALL )
          IF( SCHR(14,IZN).LE.ZERO ) THEN
            IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
              CM = 1.D+0 - 2.D+0/CN
            ELSE
              CM = 1.D+0 - 1.D+0/CN
            ENDIF
          ELSE
            CM = SCHR(14,IZN)
          ENDIF
          ASLX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGL)**CN))**CM
          SLRX = SCHR(4,IZN)
!
!---      Compute trapped gas saturation, using the minimum
!         apparent aqueous saturation  ---
!
          IF( INDX.EQ.0 ) THEN
            ASLM = MAX( MIN( ASLX,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
              R = 1.D+0/ESGTMX - 1.D+0
              ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
              IF( ESGTX.LT.EPSL ) ESGTX = 0.D+0
            ELSE
              ESGTX = 0.D+0
            ENDIF
!
!---      Trapped-gas saturation given, compute the minimum
!         apparent aqueous saturation  ---
!
          ELSEIF( INDX.EQ.1 ) THEN
            IF( ESGTMX.GT.EPSL ) THEN
              ESGTX = SGX/(1.D+0-SLRX)
              ASLM = (ESGTMX - ESGTX)/(ESGTMX + ESGTMX*ESGTX - ESGTX)
              ASLMINX = MAX( MIN( ASLM,ASLX ),0.D+0 )
            ELSE
              ESGTX = 0.D+0
              ASLMINX = MAX( ASLX,0.D+0 )
            ENDIF
!
!---      Unrecognized option  ---
!
          ELSE
            CHMSG = 'Unrecognized Trapped-Gas Option'
            INDX = 12
            IMSG = INDX
            CALL WRMSGS( INDX )
          ENDIF
          ESLX = ASLX - ESGTX
          SLX = ESLX*(1.D+0-SLRX) + SLRX
          SGTX = ESGTX*(1.D+0-SLRX)
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
        ENDIF
!
!---  Brooks and Corey saturation function w/ gas entrapment,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.102 ) THEN
!
!---    Brooks and Corey saturation function w/ gas entrapment, 
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          HMPX = SCHR(9,IZN)
          SLRX = SCHR(4,IZN)
!
!---      Capillary head above the matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPX ) THEN
            SMPX = SCHR(8,IZN)
            HDGL = MIN( HDGL,SCHR(12,IZN) )
            DMPX = SMPX/(LOG10(SCHR(12,IZN))-LOG10(HMPX))
            SLX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPX
            ASLX = SLX
!
!---      Capillary head at or below the matching point head,
!         use van Genuchten function
!
          ELSE
            CL = MAX( SCHR(3,IZN),SMALL )
            IF( HDGL.LE.(SCHR(1,IZN)+1.D-9) ) THEN
              ASLX = 1.D+0
            ELSE
              ASLX = (SCHR(1,IZN)/HDGL)**CL
            ENDIF
            SLRX = SCHR(4,IZN)
            SLX = ASLX*(1.D+0-SLRX) + SLRX
            ASLX = SLX
          ENDIF
!
!---      Compute trapped gas saturation, using the minimum
!         apparent aqueous saturation  ---
!
          IF( INDX.EQ.0 ) THEN
            ASLM = MAX( MIN( ASLX,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
              SGTMX = ESGTMX
              R = 1.D+0/SGTMX - 1.D+0
              ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
              IF( ESGTX.LT.EPSL ) ESGTX = 0.D+0
            ELSE
              ESGTX = 0.D+0
            ENDIF
!
!---      Trapped-gas saturation given, compute the minimum
!         apparent aqueous saturation  ---
!
          ELSEIF( INDX.EQ.1 ) THEN
            IF( ESGTMX.GT.EPSL ) THEN
              ESGTX = SGX
              SGTMX = ESGTMX
              ASLM = (SGTMX - SGTX)/(SGTMX + SGTMX*SGTX - SGTX)
              ASLMINX = MAX( MIN( ASLM,ASLX ),0.D+0 )
            ELSE
              ESGTX = 0.D+0
              ASLMINX = MAX( ASLX,0.D+0 )
            ENDIF
!
!---      Unrecognized option  ---
!
          ELSE
            CHMSG = 'Unrecognized Trapped-Gas Option'
            INDX = 12
            IMSG = INDX
            CALL WRMSGS( INDX )
          ENDIF
          SLX = ASLX - ESGTX
          SGTX = ESGTX
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
!
!---    Brooks and Corey saturation function w/ gas entrapment, 
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          CL = MAX( SCHR(3,IZN),SMALL )
          IF( HDGL.LE.(SCHR(1,IZN)+1.D-9) ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (SCHR(1,IZN)/HDGL)**CL
          ENDIF
          SLRX = SCHR(4,IZN)
!
!---      Compute trapped gas saturation, using the minimum
!         apparent aqueous saturation  ---
!
          IF( INDX.EQ.0 ) THEN
            ASLM = MAX( MIN( ASLX,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
              R = 1.D+0/ESGTMX - 1.D+0
              ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
              IF( ESGTX.LT.EPSL ) ESGTX = 0.D+0
            ELSE
              ESGTX = 0.D+0
            ENDIF
!
!---      Trapped-gas saturation given, compute the minimum
!         apparent aqueous saturation  ---
!
          ELSEIF( INDX.EQ.1 ) THEN
            IF( ESGTMX.GT.EPSL ) THEN
              ESGTX = SGX/(1.D+0-SLRX)
              ASLM = (ESGTMX - ESGTX)/(ESGTMX + ESGTMX*ESGTX - ESGTX)
              ASLMINX = MAX( MIN( ASLM,ASLX ),0.D+0 )
            ELSE
              ESGTX = 0.D+0
              ASLMINX = MAX( ASLX,0.D+0 )
            ENDIF
!
!---      Unrecognized option  ---
!
          ELSE
            CHMSG = 'Unrecognized Trapped-Gas Option'
            INDX = 12
            IMSG = INDX
            CALL WRMSGS( INDX )
          ENDIF
          ESLX = ASLX - ESGTX
          SLX = ESLX*(1.D+0-SLRX) + SLRX
          SGTX = ESGTX*(1.D+0-SLRX)
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
        ENDIF
!
!---  van Genuchten drainage-imbibition saturation function,
!     w/ or w/o extensions---
!
      ELSEIF( ISCHR(IZN).EQ.201 ) THEN
!
!---    van Genuchten saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          HMPDX = SCHR(9,IZN)
          HMPIX = SCHR(11,IZN)
          SLRX = SCHR(4,IZN)
          HDGL = MIN( HDGL,SCHR(12,IZN) )
!
!---      Capillary head above the drainage matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPDX ) THEN
              SMPDX = SCHR(8,IZN)
              DMPDX = SMPDX/(LOG10(SCHR(12,IZN))-LOG10(HMPDX))
              SLDX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPDX
              ASLDX = SLDX
              ASLMX = MAX( MIN( ASLDX,ASLMINX ),0.D+0 )
!
!---        Capillary head above the imbibition matching point head,
!           use Webb extension  ---
!
            IF( HDGL.GT.HMPIX ) THEN
              SMPIX = SCHR(10,IZN)
              DMPIX = SMPIX/(LOG10(SCHR(12,IZN))-LOG10(HMPIX))
              SLIX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPIX
              ASLIX = SLIX
              IF( ASLDX.GE.ASLIX ) THEN
                ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ELSE
                ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ENDIF
!
!---        Capillary head at or below the imbibition matching 
!           point head, use van Genuchten function
!
            ELSE
              CNI = MAX( SCHR(5,IZN),SMALL )
              IF( SCHR(13,IZN).LE.ZERO ) THEN
                IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                  CMI = 1.D+0 - 2.D+0/CNI
                ELSE
                  CMI = 1.D+0 - 1.D+0/CNI
                ENDIF
              ELSE
                CMI = SCHR(13,IZN)
              ENDIF
              ASLIX = (1.D+0/(1.D+0 + (SCHR(2,IZN)*HDGL)**CNI))**CMI
              ASLIX = ASLIX*(1.D+0-SLRX) + SLRX
              IF( ASLDX.GE.ASLIX ) THEN
                ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ELSE
                ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ENDIF
            ENDIF
!
!---      Capillary head at or below the drainage matching point head,
!         use van Genuchten function
!
          ELSE
            CND = MAX( SCHR(3,IZN),SMALL )
            IF( SCHR(14,IZN).LE.ZERO ) THEN
              IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                CMD = 1.D+0 - 2.D+0/CND
              ELSE
                CMD = 1.D+0 - 1.D+0/CND
              ENDIF
            ELSE
              CMD = SCHR(14,IZN)
            ENDIF
            ASLDX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGL)**CND))**CMD
            ASLDX = ASLDX*(1.D+0-SLRX) + SLRX
            ASLMX = MAX( MIN( ASLDX,ASLMINX ),0.D+0 )
!
!---        Capillary head above the imbibition matching point head,
!           use Webb extension  ---
!
            IF( HDGL.GT.HMPIX ) THEN
              SMPIX = SCHR(10,IZN)
              DMPIX = SMPIX/(LOG10(SCHR(12,IZN))-LOG10(HMPIX))
              SLIX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPIX
              ASLIX = SLIX
              IF( ASLDX.GE.ASLIX ) THEN
                ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ELSE
                ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ENDIF
            ELSE
              CNI = MAX( SCHR(5,IZN),SMALL )
              IF( SCHR(13,IZN).LE.ZERO ) THEN
                IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
                  CMI = 1.D+0 - 2.D+0/CNI
                ELSE
                  CMI = 1.D+0 - 1.D+0/CNI
                ENDIF
              ELSE
                CMI = SCHR(13,IZN)
              ENDIF
              ASLIX = (1.D+0/(1.D+0 + (SCHR(2,IZN)*HDGL)**CNI))**CMI
              ASLIX = ASLIX*(1.D+0-SLRX) + SLRX
              IF( ASLDX.GE.ASLIX ) THEN
                ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ELSE
                ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ENDIF
            ENDIF
          ENDIF
!
!---      Compute trapped gas saturation, using the minimum
!         apparent aqueous saturation  ---
!
          IF( INDX.EQ.0 ) THEN
            ASLM = MAX( MIN( ASLX,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
              SGTMX = ESGTMX
              R = 1.D+0/SGTMX - 1.D+0
              ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
              IF( ESGTX.LT.EPSL ) ESGTX = 0.D+0
            ELSE
              ESGTX = 0.D+0
            ENDIF
!
!---      Trapped-gas saturation given, compute the minimum
!         apparent aqueous saturation  ---
!
          ELSEIF( INDX.EQ.1 ) THEN
            IF( ESGTMX.GT.EPSL ) THEN
              ESGTX = SGX
              SGTMX = ESGTMX
              ASLM = (SGTMX - SGTX)/(SGTMX + SGTMX*SGTX - SGTX)
              ASLMINX = MAX( MIN( ASLM,ASLX ),0.D+0 )
            ELSE
              ESGTX = 0.D+0
              ASLMINX = MAX( ASLX,0.D+0 )
            ENDIF
!
!---      Unrecognized option  ---
!
          ELSE
            CHMSG = 'Unrecognized Trapped-Gas Option'
            INDX = 12
            IMSG = INDX
            CALL WRMSGS( INDX )
          ENDIF
          SLX = ASLX - ESGTX
          SGTX = ESGTX
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
!
!---    van Genuchten saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          CND = MAX( SCHR(3,IZN),SMALL )
          IF( SCHR(14,IZN).LE.ZERO ) THEN
            IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
              CMD = 1.D+0 - 2.D+0/CND
            ELSE
              CMD = 1.D+0 - 1.D+0/CND
            ENDIF
          ELSE
            CMD = SCHR(14,IZN)
          ENDIF
          ASLDX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HDGL)**CND))**CMD
          CNI = MAX( SCHR(5,IZN),SMALL )
          IF( SCHR(13,IZN).LE.ZERO ) THEN
            IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
              CMI = 1.D+0 - 2.D+0/CNI
            ELSE
              CMI = 1.D+0 - 1.D+0/CNI
            ENDIF
          ELSE
            CMI = SCHR(13,IZN)
          ENDIF
          ASLIX = (1.D+0/(1.D+0 + (SCHR(2,IZN)*HDGL)**CNI))**CMI
          ASLMX = MAX( MIN( ASLDX,ASLMINX ),0.D+0 )
          IF( ASLDX.GE.ASLIX ) THEN
            ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &         + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
          ELSE
            ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &         + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
          ENDIF
!
!---      Compute trapped gas saturation, using the minimum
!         apparent aqueous saturation  ---
!
          IF( INDX.EQ.0 ) THEN
            ASLM = MAX( MIN( ASLX,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
              R = 1.D+0/ESGTMX - 1.D+0
              ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
              IF( ESGTX.LT.EPSL ) ESGTX = 0.D+0
            ELSE
              ESGTX = 0.D+0
            ENDIF
!
!---      Trapped-gas saturation given, compute the minimum
!         apparent aqueous saturation  ---
!
          ELSEIF( INDX.EQ.1 ) THEN
            IF( ESGTMX.GT.EPSL ) THEN
              ESGTX = SGX/(1.D+0-SLRX)
              ASLM = (ESGTMX - ESGTX)/(ESGTMX + ESGTMX*ESGTX - ESGTX)
              ASLMINX = MAX( MIN( ASLM,ASLX ),0.D+0 )
            ELSE
              ESGTX = 0.D+0
              ASLMINX = MAX( ASLX,0.D+0 )
            ENDIF
!
!---      Unrecognized option  ---
!
          ELSE
            CHMSG = 'Unrecognized Trapped-Gas Option'
            INDX = 12
            IMSG = INDX
            CALL WRMSGS( INDX )
          ENDIF
          ESLX = ASLX - ESGTX
          SLX = ESLX*(1.D+0-SLRX) + SLRX
          SGTX = ESGTX*(1.D+0-SLRX)
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
        ENDIF
!
!---  Brooks and Corey drainage-imbibition saturation function,
!     w/ or w/o Webb extension  ---
!
      ELSEIF( ISCHR(IZN).EQ.202 ) THEN
!
!---    Brooks and Corey drainage-imbibition saturation function,
!       w/ Webb extension ASL = SL + SGT  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HDGL = MAX( (PGX-PLX)/BTGLX/RHORL/GRAV,1.D-14 )
          HMPDX = SCHR(9,IZN)
          HMPIX = SCHR(11,IZN)
          SLRX = SCHR(4,IZN)
          HDGL = MIN( HDGL,SCHR(12,IZN) )
!
!---      Capillary head above the drainage matching point head,
!         use Webb extension  ---
!
          IF( HDGL.GT.HMPDX ) THEN
            SMPDX = SCHR(8,IZN)
            DMPDX = SMPDX/(LOG10(SCHR(12,IZN))-LOG10(HMPDX))
            SLDX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPDX
            ASLDX = SLDX
            ASLMX = MAX( MIN( ASLDX,ASLMINX ),0.D+0 )
!
!---        Capillary head above the imbibition matching point head,
!           use Webb extension  ---
!
            IF( HDGL.GT.HMPIX ) THEN
              SMPIX = SCHR(10,IZN)
              DMPIX = SMPIX/(LOG10(SCHR(12,IZN))-LOG10(HMPIX))
              SLIX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPIX
              ASLIX = SLIX
              IF( ASLDX.GE.ASLIX ) THEN
                ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ELSE
                ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ENDIF
!
!---        Capillary head at or below the imbibition matching 
!           point head, use Brooks and Corey function
!
            ELSE
              CLI = MAX( SCHR(6,IZN),SMALL )
              IF( HDGL.LE.SCHR(5,IZN) ) THEN
                ASLIX = 1.D+0
              ELSE
                ASLIX = (SCHR(5,IZN)/HDGL)**CLI
              ENDIF
              ASLIX = ASLIX*(1.D+0-SLRX) + SLRX
              IF( ASLDX.GE.ASLIX ) THEN
                ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ELSE
                ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ENDIF
            ENDIF
!
!---      Capillary head at or below the drainage matching point head,
!         use Brooks and Corey function
!
          ELSE
            CLD = MAX( SCHR(3,IZN),SMALL )
            IF( HDGL.LE.(SCHR(1,IZN)+1.D-9) ) THEN
              ASLDX = 1.D+0
            ELSE
              ASLDX = (SCHR(1,IZN)/HDGL)**CLD
            ENDIF
            ASLDX = ASLDX*(1.D+0-SLRX) + SLRX
            ASLMX = MAX( MIN( ASLDX,ASLMINX ),0.D+0 )
!
!---        Capillary head above the imbibition matching point head,
!           use Webb extension  ---
!
            IF( HDGL.GT.HMPIX ) THEN
              SMPIX = SCHR(10,IZN)
              DMPIX = SMPIX/(LOG10(SCHR(12,IZN))-LOG10(HMPIX))
              SLIX = -(LOG10(HDGL)-LOG10(SCHR(12,IZN)))*DMPIX
              ASLIX = SLIX
              IF( ASLDX.GE.ASLIX ) THEN
                ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ELSE
                ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ENDIF
!
!---        Capillary head at or below the imbibition matching 
!           point head, use Brooks and Corey function
!
            ELSE
              CLI = MAX( SCHR(6,IZN),SMALL )
              IF( HDGL.LE.SCHR(5,IZN) ) THEN
                ASLIX = 1.D+0
              ELSE
                ASLIX = (SCHR(5,IZN)/HDGL)**CLI
              ENDIF
              ASLIX = ASLIX*(1.D+0-SLRX) + SLRX
              IF( ASLDX.GE.ASLIX ) THEN
                ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ELSE
                ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
              ENDIF
            ENDIF
          ENDIF
!
!---      Compute trapped gas saturation, using the minimum
!         apparent aqueous saturation  ---
!
          IF( INDX.EQ.0 ) THEN
            ASLM = MAX( MIN( ASLX,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
              SGTMX = ESGTMX
              R = 1.D+0/SGTMX - 1.D+0
              ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
              IF( ESGTX.LT.EPSL ) ESGTX = 0.D+0
            ELSE
              ESGTX = 0.D+0
            ENDIF
!
!---      Trapped-gas saturation given, compute the minimum
!         apparent aqueous saturation  ---
!
          ELSEIF( INDX.EQ.1 ) THEN
            IF( ESGTMX.GT.EPSL ) THEN
              ESGTX = SGX
              SGTMX = ESGTMX
              ASLM = (SGTMX - SGTX)/(SGTMX + SGTMX*SGTX - SGTX)
              ASLMINX = MAX( MIN( ASLM,ASLX ),0.D+0 )
            ELSE
              ESGTX = 0.D+0
              ASLMINX = MAX( ASLX,0.D+0 )
            ENDIF
!
!---      Unrecognized option  ---
!
          ELSE
            CHMSG = 'Unrecognized Trapped-Gas Option'
            INDX = 12
            IMSG = INDX
            CALL WRMSGS( INDX )
          ENDIF
          SLX = ASLX - ESGTX
          SGTX = ESGTX
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
!
!---    Brooks and Corey drainage-imbibition saturation function,
!       w/o Webb extension ASL = ESL + ESGT  ---
!
        ELSE
          CLD = MAX( SCHR(3,IZN),SMALL )
          IF( HDGL.LE.(SCHR(1,IZN)+1.D-9) ) THEN
            ASLDX = 1.D+0
          ELSE
            ASLDX = (SCHR(1,IZN)/HDGL)**CLD
          ENDIF
          CLI = MAX( SCHR(6,IZN),SMALL )
          IF( HDGL.LE.SCHR(5,IZN) ) THEN
            ASLIX = 1.D+0
          ELSE
            ASLIX = (SCHR(5,IZN)/HDGL)**CLI
          ENDIF
          ASLMX = MAX( MIN( ASLDX,ASLMINX ),0.D+0 )
          IF( ASLDX.GE.ASLIX ) THEN
            ASLX = 5.D-1*(ASLIX + SQRT((ASLIX**2)
     &         + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
          ELSE
            ASLX = 5.D-1*(ASLIX - SQRT((ASLIX**2)
     &             + 4.D+0*ASLMX*ASLDX - 4.D+0*ASLIX*ASLMX))
          ENDIF
!
!---      Compute trapped gas saturation, using the minimum
!         apparent aqueous saturation  ---
!
          IF( INDX.EQ.0 ) THEN
            ASLM = MAX( MIN( ASLX,ASLMINX ),0.D+0 )
            IF( ESGTMX.GT.EPSL .AND. ASLX.GT.ASLM ) THEN
              R = 1.D+0/ESGTMX - 1.D+0
              ESGTX = (1.D+0-ASLM)/(1.D+0 + R*(1.D+0-ASLM)) -
     &          (1.D+0-ASLX)/(1.D+0 + R*(1.D+0-ASLX))
              IF( ESGTX.LT.EPSL ) ESGTX = 0.D+0
            ELSE
              ESGTX = 0.D+0
            ENDIF
!
!---      Trapped-gas saturation given, compute the minimum
!         apparent aqueous saturation  ---
!
          ELSEIF( INDX.EQ.1 ) THEN
            IF( ESGTMX.GT.EPSL ) THEN
              ESGTX = SGX/(1.D+0-SLRX)
              ASLM = (ESGTMX - ESGTX)/(ESGTMX + ESGTMX*ESGTX - ESGTX)
              ASLMINX = MAX( MIN( ASLM,ASLX ),0.D+0 )
            ELSE
              ESGTX = 0.D+0
              ASLMINX = MAX( ASLX,0.D+0 )
            ENDIF
!
!---      Unrecognized option  ---
!
          ELSE
            CHMSG = 'Unrecognized Trapped-Gas Option'
            INDX = 12
            IMSG = INDX
            CALL WRMSGS( INDX )
          ENDIF
          ESLX = ASLX - ESGTX
          SLX = ESLX*(1.D+0-SLRX) + SLRX
          SGTX = ESGTX*(1.D+0-SLRX)
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          SGX = 1.D+0-SLX
          IF( SGX.LT.EPSL ) SGX = 0.D+0
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORC_GT
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
!     Compute source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE PORMED
      USE OUTPU
      USE JACOB
      USE GRID
      USE FDVT
      USE FDVS
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
      REAL*8 SRX(8+LSOLU)
      REAL*8 HCFG(LFZ)
      REAL*8 MAXX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC_GT'
!
!---  Zero source terms  ---
!
      DO N = 1,NFLD
        IF( IXP(N).EQ.0 ) CYCLE
        DO M = 2,ISVC+2
          SRCT(M,N) = 0.D+0
          SRCS(M,N) = 0.D+0
          SRCA(M,N) = 0.D+0
          SRCW(M,N) = 0.D+0
        ENDDO
      ENDDO
!
!---  Loop over sources  ---
!
      DO NS = 1,NSR
        IF( TM.LE.SRC(1,1,NS) ) CYCLE
        SRX(1) = TM
        IF( ISRM(NS).EQ.1 .OR. ISRT(NS).GT.100 ) THEN
          DO N = 2,8+NSOLU
            SRX(N) = SRC(N,1,NS)
          ENDDO
        ELSE
          DO M = 2,ISRM(NS)
            IF( TM.LE.SRC(1,M,NS) ) THEN
              DTSR = MIN( SRC(1,M,NS)-TM,DT )
              TFSR = (TM-0.5D+0*DTSR-SRC(1,M-1,NS))/
     &          (SRC(1,M,NS)-SRC(1,M-1,NS))
              DO N = 2,8+NSOLU
                SRX(N) = SRC(N,M-1,NS)+TFSR*(SRC(N,M,NS)-SRC(N,M-1,NS))
              ENDDO
              GOTO 110
            ENDIF
          ENDDO
          CYCLE
        ENDIF
  110   CONTINUE
!---  Skip over domain loop for well sources  ---
!
      IF( ISRT(NS).GE.23 .AND. ISRT(NS).LE.25 ) GOTO 510
!
!---  Loop source domain  ---
!
        DO I = ISRDM(1,NS),ISRDM(2,NS)
          DO J = ISRDM(3,NS),ISRDM(4,NS)
            DO K = ISRDM(5,NS),ISRDM(6,NS)
              N = ND(I,J,K)
              IZN = IZ(N)
              IF( IXP(N).EQ.0 ) CYCLE
              DO 400 M = 2,ISVC+2
                PGX = PG(M,N) + PATM
                PLX = PL(M,N) + PATM
!
!---            Specified Power  ---
!
                IF( ISRT(NS).EQ.1 ) THEN
                  SRCT(M,N) = SRCT(M,N) + SRX(4)
!
!---            Power Density  ---
!
                ELSEIF( ISRT(NS).EQ.2 ) THEN
                  SRCT(M,N) = SRCT(M,N) + SRX(4)*VOL(N)
!
!---            Setting of Portland Cement  ---
!
                ELSEIF( ISRT(NS).EQ.11 ) THEN
!
!---              Temperature dependent Sigmoid function parameters
!                 for Portland cement  ---
!
!                  XHALFX = 3.9653D+3 + 9.7147D+4*EXP(-4.7852D-2*T(M,N))
!                  RATEX = 8.2143D+2 + 3.0896D+5*EXP(-1.227D-1*T(M,N))
!
!---              Temperature dependent Sigmoid function parameters
!                 for Portland cement, with 2x rate and xhalf
!                 shifted later by 2000 s  ---
!
!                  XHALFX = 5.9653D+3 + 9.7147D+4*EXP(-4.7852D-2*T(M,N))
!                  RATEX = 16.429D+2 + 6.1792D+5*EXP(-1.227D-1*T(M,N))
                  XHALFX = SRX(2) + SRX(3)*EXP(-SRX(5)*T(M,N))
                  RATEX = SRX(6) + SRX(7)*EXP(-SRX(8)*T(M,N))
                  BASEX = 0.D+0
                  MAXX = SRX(4)
                  VARX = EXP((TM-XHALFX)/RATEX)
!
!---              Heat evolution rate, W/kg  ---
!
                  HERX = MAXX*VARX/(((1.D+0 + VARX)**2)*RATEX)
                  SRCT(M,N) = SRCT(M,N) + HERX*VOL(N)*RHOS(IZN)
!
!---            Aqueous Volumetric Rate w/ 
!               Dissolved-Air Relative Saturation ---
!
                ELSEIF( ISRT(NS).EQ.3 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    TX = SRX(2)
                    PX = SRX(3)
                    WLAX = SRX(5)
                    XLSX = SRX(6)
                    IF( ISLC(32).EQ.0 ) XLSX = XLS(2,N)
                    CALL SOL_BRNS( TX,PX,XLSMX )
                    XLSX = XLSX*XLSMX
                    INDX = 0
                    CALL REGION_4( TX,PSWX,INDX )
                    PX = MAX( PGX,PLX,PX,PSWX )
                    PWX = PX
                    CALL P_IAPWS( TX,PWX,RHOX,RHOLWX,HX,HLWX,UX,ULWX )
                    PVAX = WLAX*MAX( PX-PSWX,0.D+0 )
                    XMLAX = PVAX/HCAW
                    XMLWX = MAX( 1.D+0-XMLAX,0.D+0 )
                    WTMX = XMLWX*WTMW + XMLAX*WTMA
                    XLAX = (XMLAX*WTMA)/WTMX
                    XLWX = (XMLWX*WTMW)/WTMX
                    CALL DENS_B( XLSX,RHOLX,RHOLWX,TX )
                    CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                    HLX = MAX(1.D+0-XLAX,0.D+0)*HLWX + XLAX*HGAX
                    SRCS(M,N) = SRCS(M,N) + SRX(4)*RHOLX*XLSX
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOLX*XLAX
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOLX*XLWX
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*RHOLX*HLX
                  ELSE
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*RHOL(M,N)*HL(M,N)
                    SRCS(M,N) = SRCS(M,N) + SRX(4)*RHOL(M,N)*XLS(M,N)
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOL(M,N)*XLA(M,N)
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOL(M,N)*XLW(M,N)
                  ENDIF
!
!---            Gas Volumetric Rate w/ Component Relative Humidity ---
!
                ELSEIF( ISRT(NS).EQ.4 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    PX = MAX( PGX,PLX,SRX(3) )
                    INDX = 0
                    CALL REGION_4( SRX(2),PSWX,INDX )
                    PVWX = PSWX*SRX(5)
                    PVAX = MAX( PX-PVWX,0.D+0 )
                    INDX = 0
                    CALL WATGSD( SRX(2),PVWX,RHOW,INDX )
                    CALL AIRGSD( SRX(2),PVAX,RHOA )
                    CALL WATGSH( SRX(2),PSWX,RHOW,HGWX,UGWX )
                    CALL AIRGSH( SRX(2),PVAX,HGAX,UGAX )
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*(RHOW*HGWX+RHOA*HGAX)
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOA
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOW
                  ELSE
                    HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*RHOG(M,N)*HGX
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOG(M,N)*XGA(M,N)
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOG(M,N)*XGW(M,N)
                  ENDIF
!
!---            Gas Volumetric Rate w/ Component Mass Fractions ---
!
                ELSEIF( ISRT(NS).EQ.5 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    PX = MAX( PGX,PLX,SRX(3) )
                    XGWX = SRX(5)
                    XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                    XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
                    INDX = 0
                    CALL REGION_4( SRX(2),PSWX,INDX )
                    PVWX = MIN( PX*XMGWX,PSWX )
                    PVAX = MAX( PX-PVWX,0.D+0 )
                    INDX = 0
                    CALL WATGSD( SRX(2),PVWX,RHOW,INDX )
                    CALL AIRGSD( SRX(2),PVAX,RHOA )
                    CALL WATGSH( SRX(2),PSWX,RHOW,HGWX,UGWX )
                    CALL AIRGSH( SRX(2),PVAX,HGAX,UGAX )
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*(RHOW*HGWX+RHOA*HGAX)
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOA
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOW
                  ELSE
                    HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*RHOG(M,N)*HGX
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOG(M,N)*XGA(M,N)
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOG(M,N)*XGW(M,N)
                  ENDIF
!
!---            Undocumented Condensate Source ---
!
                ELSEIF( ISRT(NS).EQ.6 ) THEN
                  NSFG = 0
                  DO 160 NSU = 1,NSF
                    IF( ISFT(NSU).EQ.-6 ) NSFG = NSFG + 1
                    IF( NSFG.EQ.INT(SRX(4)) ) THEN
                      CMSX = MAX(SF(1,NSU),0.D+0)*VOL(N)/SRX(6)
                      GOTO 161
                    ENDIF
  160             CONTINUE
  161             CONTINUE
                  CMSX = CMSX*MAX( 1.D+0-SRX(7),0.D+0 )
                  PX = MAX( PGX,PLX,SRX(3) )
                  CALL WATLQH( SRX(2),PX,HLX )
                  INDX = 0
                  CALL REGION_4( SRX(2),PSWX,INDX )
                  PVAX = MAX( PX-PSWX,0.D+0 )*SRX(5)
                  XMLAX = PVAX/HCAW
                  XMLWX = MAX( 1.D+0-XMLAX,0.D+0 )
                  XLAX = (XMLAX*WTMA)/(XMLAX*WTMA+XMLWX*WTMW)
                  SRCT(M,N) = SRCT(M,N) + CMSX*HLX
                  SRCA(M,N) = SRCA(M,N) + CMSX*XLAX
                  SRCW(M,N) = SRCW(M,N) + CMSX*(1.D+0-XLAX)
!
!---            Aqueous Mass Rate w/ Dissolved-Air Relative Saturation  ---
!
                ELSEIF( ISRT(NS).EQ.7 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    TX = SRX(2)
                    PX = SRX(3)
                    WLAX = SRX(5)
                    INDX = 0
                    CALL REGION_4( TX,PSWX,INDX )
                    PX = MAX( PGX,PLX,PX,PSWX )
                    PWX = PX
                    CALL P_IAPWS( TX,PWX,RHOX,RHOLWX,HX,HLWX,UX,ULWX )
                    PVAX = WLAX*MAX( PX-PSWX,0.D+0 )
                    XMLAX = PVAX/HCAW
                    XMLWX = MAX( 1.D+0-XMLAX,0.D+0 )
                    WTMX = XMLWX*WTMW + XMLAX*WTMA
                    XLAX = (XMLAX*WTMA)/WTMX
                    XLWX = (XMLWX*WTMW)/WTMX
                    CALL AIRGSH( TX,PVAX,HGAX,UEGAX )
                    HLX = MAX(1.D+0-XLAX,0.D+0)*HLWX + XLAX*HGAX
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*XLAX
                    SRCS(M,N) = SRCS(M,N) + SRX(4)*XLSX
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*XLWX
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*HLX
                  ELSE
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*HL(M,N)
                    SRCS(M,N) = SRCS(M,N) + SRX(4)*XLS(M,N)
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*XLA(M,N)
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*XLW(M,N)
                  ENDIF
!
!---            Gas Mass Rate w/ Component Relative Humidity ---
!
                ELSEIF( ISRT(NS).EQ.8 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    PX = MAX( PGX,PLX,SRX(3) )
                    INDX = 0
                    CALL REGION_4( SRX(2),PSWX,INDX )
                    PVWX = PSWX*SRX(5)
                    PVAX = MAX( PGX-PVWX,0.D+0 )
                    INDX = 0
                    CALL WATGSD( SRX(2),PVWX,RHOW,INDX )
                    CALL AIRGSD( SRX(2),PVAX,RHOA )
                    CALL WATGSH( SRX(2),PSWX,RHOW,HGWX,UGWX )
                    CALL AIRGSH( SRX(2),PVAX,HGAX,UGAX )
                    XGWX = RHOW/(RHOW+RHOA)
                    XGAX = RHOA/(RHOW+RHOA)
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*(XGWX*HGWX+XGAX*HGAX)
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*XGAX
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*XGWX
                  ELSE
                    HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*HGX
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*XGA(M,N)
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*XGW(M,N)
                  ENDIF
!
!---            Gas Mass Rate w/ Component Mass Fractions ---
!
                ELSEIF( ISRT(NS).EQ.9 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    PX = MAX( PLX,PGX,SRX(3) )
                    XGWX = SRX(5)
                    XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                    XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
                    INDX = 0
                    CALL REGION_4( SRX(2),PSWX,INDX )
                    PVWX = MIN( PX*XMGWX,PSWX )
                    PVAX = MAX( PX-PVWX,0.D+0 )
                    INDX = 0
                    CALL WATGSD( SRX(2),PVWX,RHOW,INDX )
                    CALL AIRGSD( SRX(2),PVAX,RHOA )
                    CALL WATGSH( SRX(2),PSWX,RHOW,HGWX,UGWX )
                    CALL AIRGSH( SRX(2),PVAX,HGAX,UGAX )
                    XGWX = RHOW/(RHOW+RHOA)
                    XGAX = RHOA/(RHOW+RHOA)
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*(XGWX*HGWX+XGAX*HGAX)
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*XGAX
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*XGWX
                  ELSE
                    HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*HGX
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*XGA(M,N)
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*XGW(M,N)
                  ENDIF
!
!---            Mass Sink ---
!
                ELSEIF( ISRT(NS).EQ.10 ) THEN
                  RKLX = (RKL(1,M,N)*RKL(2,M,N)*RKL(3,M,N))**
     &              (1.D+0/3.D+0)
                  RMLX = (RKLX*RHOL(M,N)/VISL(M,N))
                  RMGX = (RKG(M,N)*RHOG(M,N)/VISG(M,N))
                  RMLX = RMLX/(RMLX+RMGX)
                  RMGX = 1.D+0-RMLX
                  SRCT(M,N) = SRCT(M,N) + SRX(4)*(RMLX*HL(M,N) + 
     &              RMGX*HG(M,N))
                  SRCA(M,N) = SRCA(M,N) + SRX(4)*(RMLX*XLA(M,N) + 
     &              RMGX*XGA(M,N))
                  SRCW(M,N) = SRCW(M,N) + SRX(4)*(RMLX*XLW(M,N) + 
     &              RMGX*XGW(M,N))
!
!---            Z-direction vertical injection well  ---
!
                ELSEIF( ISRT(NS).GE.13 .AND. ISRT(NS).LE.15  ) THEN
!
!---              Geometric factors  ---
!
                  RDW = SRX(3)
                  RDE = SQRT( AFZ(NSZ(N))/GPI/SRX(4) )
                  DRD2 = (RDE**2-RDW**2)
!
!---              Well pressure  ---
!
                  IF( M.EQ.2 ) THEN
                    IF( K.EQ.ISRDM(5,NS) ) THEN
                      PGWX = SRX(2)
                      PX = PGWX+PATM
                      IF( ISRT(NS).EQ.14 ) THEN
                        XGWX = SRX(5)
                        XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                        XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
                      ENDIF
                      INDX = 0
                      CALL REGION_4( SRX(7),PSWX,INDX )
                      IF( ISRT(NS).EQ.14 ) THEN
                        PVWX = MIN( PX*XMGWX,PSWX )
                      ELSE
                        PVWX = PSWX*SRX(5)
                      ENDIF
                      PVAX = MAX( PX-PVWX,0.D+0 )
                      INDX = 0
                      CALL WATGSD( SRX(7),PVWX,RHOGWX,INDX )
                      CALL AIRGSD( SRX(7),PVAX,RHOGAX )
                      CALL WATGSH( SRX(7),PSWX,RHOGWX,HGWX,UGWX )
                      CALL AIRGSH( SRX(7),PVAX,HGAX,UGAX )
                    ELSE
                      PX = PGWX+PATM
                      IF( ISRT(NS).EQ.14 ) THEN
                        XGWX = SRX(5)
                        XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                        XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
                      ENDIF
                      INDX = 0
                      CALL REGION_4( SRX(7),PSWX,INDX )
                      IF( ISRT(NS).EQ.14 ) THEN
                        PVWX = MIN( PX*XMGWX,PSWX )
                      ELSE
                        PVWX = PSWX*SRX(5)
                      ENDIF
                      PVAX = MAX( PX-PVWX,0.D+0 )
                      INDX = 0
                      CALL WATGSD( SRX(7),PVWX,RHOGWX,INDX )
                      CALL AIRGSD( SRX(7),PVAX,RHOGAX )
                      CALL WATGSH( SRX(7),PSWX,RHOGWX,HGWX,UGWX )
                      CALL AIRGSH( SRX(7),PVAX,HGAX,UGAX )
                      NX = ND(I,J,K-1)
                      GB = (ZP(N)-ZP(NX))*GRAVZ
                      RHOGX = RHOGWX+RHOGAX
                      PGWX = PGWX - RHOGX*GB
                    ENDIF
                  ENDIF
                  IF( (PERM(1,IZN)/EPSL).GT.EPSL ) THEN
                    IF( (PERM(2,IZN)/EPSL).GT.EPSL ) THEN
                      PERMX = SQRT( PERM(1,IZN)*PERM(2,IZN) )
                    ELSE
                      PERMX = PERM(1,IZN)
                    ENDIF
                  ELSE
                    PERMX = PERM(2,IZN)
                  ENDIF
!
!---              If permeability is too low, skip
!
                  IF ( PERMX.GT.SRX(6) ) THEN
!
!---                Injection  ---
!
                    IF( PGWX-PG(M,N).GT.EPSL ) THEN
                      XGAX = RHOGAX/(RHOGAX+RHOGWX)
                      XGWX = RHOGWX/(RHOGAX+RHOGWX)
                      WTMGX = 1.D+0/(XGAX/WTMA + XGWX/WTMW)
                      XMGAX = XGAX*WTMGX/WTMA
                      XMGWX = XGWX*WTMGX/WTMW
                      CALL WATGSV( SRX(7),VISGWX )
                      CALL AIRGSV( SRX(7),VISGAX )
                      CALL GASVIS( XMGWX,ZERO,XMGAX,VISGWX,SMALL,VISGAX,
     &                  VISGX )
                      HCGX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
     &                  (VISGX*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
                      QGX = SRX(4)*(PGWX-PG(M,N))*HCGX
                      SRCT(M,N) = SRCT(M,N) + QGX*(RHOGWX*HGWX+RHOGAX
     &                  *HGAX)
                      SRCW(M,N) = SRCW(M,N) + QGX*RHOGWX
                      SRCA(M,N) = SRCA(M,N) + QGX*RHOGAX
!
!---                Withdrawl  ---
!
                    ELSE
                      HCGX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
     &                  (VISG(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
                      QGX = SRX(4)*(PGWX-PG(M,N))*RKG(M,N)*HCGX
                      HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                      SRCT(M,N) = SRCT(M,N) + QGX*RHOG(M,N)*HGX
                      SRCW(M,N) = SRCW(M,N) + QGX*RHOG(M,N)*XGW(M,N)
                      SRCA(M,N) = SRCA(M,N) + QGX*RHOG(M,N)*XGA(M,N)
                    ENDIF
                  ENDIF
!
!---            X-direction horizontal injection well  ---
!
                ELSEIF( ISRT(NS).GE.16 .AND. ISRT(NS).LE.18  ) THEN
!
!---              Geometric factors  ---
!
                  RDW = SRX(3)
                  RDE = SQRT( AFX(NSX(N))/GPI/SRX(4) )
                  DRD2 = (RDE**2-RDW**2)
!
!---              Well pressure  ---
!
                  IF( M.EQ.2 ) THEN
                    IF( I.EQ.ISRDM(1,NS) ) THEN
                      PGWX = SRX(2)
                      PX = PGWX+PATM
                      IF( ISRT(NS).EQ.14 ) THEN
                        XGWX = SRX(5)
                        XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                        XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
                      ENDIF
                      INDX = 0
                      CALL REGION_4( SRX(7),PSWX,INDX )
                      IF( ISRT(NS).EQ.14 ) THEN
                        PVWX = MIN( PX*XMGWX,PSWX )
                      ELSE
                        PVWX = PSWX*SRX(5)
                      ENDIF
                      PVAX = MAX( PX-PVWX,0.D+0 )
                      INDX = 0
                      CALL WATGSD( SRX(7),PVWX,RHOGWX,INDX )
                      CALL AIRGSD( SRX(7),PVAX,RHOGAX )
                      CALL WATGSH( SRX(7),PSWX,RHOGWX,HGWX,UGWX )
                      CALL AIRGSH( SRX(7),PVAX,HGAX,UGAX )
                    ELSE
                      PX = PGWX+PATM
                      IF( ISRT(NS).EQ.14 ) THEN
                        XGWX = SRX(5)
                        XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                        XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
                      ENDIF
                      INDX = 0
                      CALL REGION_4( SRX(7),PSWX,INDX )
                      IF( ISRT(NS).EQ.14 ) THEN
                        PVWX = MIN( PX*XMGWX,PSWX )
                      ELSE
                        PVWX = PSWX*SRX(5)
                      ENDIF
                      PVAX = MAX( PX-PVWX,0.D+0 )
                      INDX = 0
                      CALL WATGSD( SRX(7),PVWX,RHOGWX,INDX )
                      CALL AIRGSD( SRX(7),PVAX,RHOGAX )
                      CALL WATGSH( SRX(7),PSWX,RHOGWX,HGWX,UGWX )
                      CALL AIRGSH( SRX(7),PVAX,HGAX,UGAX )
                      NX = ND(I-1,J,K)
                      GB = (XP(N)-XP(NX))*GRAVX
                      RHOGX = RHOGWX+RHOGAX
                      PGWX = PGWX - RHOGX*GB
                    ENDIF
                  ENDIF
                  IF( (PERM(2,IZN)/EPSL).GT.EPSL ) THEN
                    IF( (PERM(3,IZN)/EPSL).GT.EPSL ) THEN
                      PERMX = SQRT( PERM(2,IZN)*PERM(3,IZN) )
                    ELSE
                      PERMX = PERM(2,IZN)
                    ENDIF
                  ELSE
                    PERMX = PERM(3,IZN)
                  ENDIF
!
!---              If permeability is too low, skip
!
                  IF ( PERMX.GT.SRX(6) ) THEN
!
!---                Injection  ---
!
                    IF( PGWX-PG(M,N).GT.EPSL ) THEN
                      XGAX = RHOGAX/(RHOGAX+RHOGWX)
                      XGWX = RHOGWX/(RHOGAX+RHOGWX)
                      WTMGX = 1.D+0/(XGAX/WTMA + XGWX/WTMW)
                      XMGAX = XGAX*WTMGX/WTMA
                      XMGWX = XGWX*WTMGX/WTMW
                      CALL WATGSV( SRX(7),VISGWX )
                      CALL AIRGSV( SRX(7),VISGAX )
                      CALL GASVIS( XMGWX,ZERO,XMGAX,VISGWX,SMALL,VISGAX,
     &                  VISGX )
                      HCGX = 2.D+0*GPI*PERMX*DRD2*DXGF(N)/
     &                  (VISGX*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
                      QGX = SRX(4)*(PGWX-PG(M,N))*HCGX
                      SRCT(M,N) = SRCT(M,N) + QGX*(RHOGWX*HGWX+RHOGAX
     &                  *HGAX)
                      SRCW(M,N) = SRCW(M,N) + QGX*RHOGWX
                      SRCA(M,N) = SRCA(M,N) + QGX*RHOGAX
!
!---                Withdrawl  ---
!
                    ELSE
                      HCGX = 2.D+0*GPI*PERMX*DRD2*DXGF(N)/
     &                  (VISG(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
                      QGX = SRX(4)*(PGWX-PG(M,N))*RKG(M,N)*HCGX
                      HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                      SRCT(M,N) = SRCT(M,N) + QGX*RHOG(M,N)*HGX
                      SRCW(M,N) = SRCW(M,N) + QGX*RHOG(M,N)*XGW(M,N)
                      SRCA(M,N) = SRCA(M,N) + QGX*RHOG(M,N)*XGA(M,N)
                    ENDIF
                  ENDIF
!
!---            Y-direction horizontal injection well  ---
!
                ELSEIF( ISRT(NS).GE.19 .AND. ISRT(NS).LE.21  ) THEN
!
!---              Geometric factors  ---
!
                  RDW = SRX(3)
                  RDE = SQRT( AFY(NSY(N))/GPI/SRX(4) )
                  DRD2 = (RDE**2-RDW**2)
!
!---              Well pressure  ---
!
                  IF( M.EQ.2 ) THEN
                    IF( J.EQ.ISRDM(3,NS) ) THEN
                      PGWX = SRX(2)
                      PX = PGWX+PATM
                      IF( ISRT(NS).EQ.14 ) THEN
                        XGWX = SRX(5)
                        XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                        XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
                      ENDIF
                      INDX = 0
                      CALL REGION_4( SRX(7),PSWX,INDX )
                      IF( ISRT(NS).EQ.14 ) THEN
                        PVWX = MIN( PX*XMGWX,PSWX )
                      ELSE
                        PVWX = PSWX*SRX(5)
                      ENDIF
                      PVAX = MAX( PX-PVWX,0.D+0 )
                      INDX = 0
                      CALL WATGSD( SRX(7),PVWX,RHOGWX,INDX )
                      CALL AIRGSD( SRX(7),PVAX,RHOGAX )
                      CALL WATGSH( SRX(7),PSWX,RHOGWX,HGWX,UGWX )
                      CALL AIRGSH( SRX(7),PVAX,HGAX,UGAX )
                    ELSE
                      PX = PGWX+PATM
                      IF( ISRT(NS).EQ.14 ) THEN
                        XGWX = SRX(5)
                        XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                        XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
                      ENDIF
                      INDX = 0
                      CALL REGION_4( SRX(7),PSWX,INDX )
                      IF( ISRT(NS).EQ.14 ) THEN
                        PVWX = MIN( PX*XMGWX,PSWX )
                      ELSE
                        PVWX = PSWX*SRX(5)
                      ENDIF
                      PVAX = MAX( PX-PVWX,0.D+0 )
                      INDX = 0
                      CALL WATGSD( SRX(7),PVWX,RHOGWX,INDX )
                      CALL AIRGSD( SRX(7),PVAX,RHOGAX )
                      CALL WATGSH( SRX(7),PSWX,RHOGWX,HGWX,UGWX )
                      CALL AIRGSH( SRX(7),PVAX,HGAX,UGAX )
                      NX = ND(I,J-1,K)
                      GB = (YP(N)-YP(NX))*RP(I)*GRAVY
                      RHOGX = RHOGWX+RHOGAX
                      PGWX = PGWX - RHOGX*GB
                    ENDIF
                  ENDIF
                  IF( (PERM(1,IZN)/EPSL).GT.EPSL ) THEN
                    IF( (PERM(3,IZN)/EPSL).GT.EPSL ) THEN
                      PERMX = SQRT( PERM(1,IZN)*PERM(3,IZN) )
                    ELSE
                      PERMX = PERM(1,IZN)
                    ENDIF
                  ELSE
                    PERMX = PERM(3,IZN)
                  ENDIF
!
!---              If permeability is too low, skip
!
                  IF ( PERMX.GT.SRX(6) ) THEN
!
!---                Injection  ---
!
                    IF( PGWX-PG(M,N).GT.EPSL ) THEN
                      XGAX = RHOGAX/(RHOGAX+RHOGWX)
                      XGWX = RHOGWX/(RHOGAX+RHOGWX)
                      WTMGX = 1.D+0/(XGAX/WTMA + XGWX/WTMW)
                      XMGAX = XGAX*WTMGX/WTMA
                      XMGWX = XGWX*WTMGX/WTMW
                      CALL WATGSV( SRX(7),VISGWX )
                      CALL AIRGSV( SRX(7),VISGAX )
                      CALL GASVIS( XMGWX,ZERO,XMGAX,VISGWX,SMALL,VISGAX,
     &                  VISGX )
                      HCGX = 2.D+0*GPI*PERMX*DRD2*DYGF(N)/
     &                  (VISGX*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
                      QGX = SRX(4)*(PGWX-PG(M,N))*HCGX
                      SRCT(M,N) = SRCT(M,N) + QGX*(RHOGWX*HGWX+RHOGAX
     &                  *HGAX)
                      SRCW(M,N) = SRCW(M,N) + QGX*RHOGWX
                      SRCA(M,N) = SRCA(M,N) + QGX*RHOGAX
!
!---                Withdrawl  ---
!
                    ELSE
                      HCGX = 2.D+0*GPI*PERMX*DRD2*DYGF(N)/
     &                  (VISG(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
                      QGX = SRX(4)*(PGWX-PG(M,N))*RKG(M,N)*HCGX
                      HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                      SRCT(M,N) = SRCT(M,N) + QGX*RHOG(M,N)*HGX
                      SRCW(M,N) = SRCW(M,N) + QGX*RHOG(M,N)*XGW(M,N)
                      SRCA(M,N) = SRCA(M,N) + QGX*RHOG(M,N)*XGA(M,N)
                    ENDIF
                  ENDIF
                ENDIF
  400         CONTINUE
            ENDDO
          ENDDO
        ENDDO
!
!---  Volumetric injection wells  ---
!
  510   CONTINUE
!
!---    Z-direction Gas Volumetric injection well  ---
!
        IF( ISRT(NS).EQ.23 .OR. ISRT(NS).EQ.24 ) THEN
          I = ISRDM(1,NS)
          J = ISRDM(3,NS)
          K1X = ISRDM(5,NS)
          K2X = ISRDM(6,NS)
!
!---      Partition the applied flux according to the air
!         permeabilities of each grid block, weighted by the
!         length of the well bore contained in each grid block.
!
          DO 518 M = 2,ISVC+2
            SUMKX = 0.D+0
            DO 512 K = K1X,K2X
              N = ND(I,J,K)
              IZN = IZ(N)
              IF( (PERM(1,IZN)/EPSL).GT.EPSL ) THEN
                IF( (PERM(2,IZN)/EPSL).GT.EPSL ) THEN
                  PERMX = SQRT( PERM(1,IZN)*PERM(2,IZN) )
                ELSE
                  PERMX = PERM(1,IZN)
                ENDIF
              ELSE
                PERMX = PERM(2,IZN)
              ENDIF
              HCFG(K) = PERMX*RKG(M,N)*DZGF(N)
              SUMKX = SUMKX + HCFG(K)
  512       CONTINUE
            DO 516 K = K1X,K2X
              N = ND(I,J,K)
              PGX = PG(M,N) + PATM
              PLX = PL(M,N) + PATM
              TX = T(M,N)
              WTX = HCFG(K)/SUMKX
!
!---          Geometric factors  ---
!
              RDW = SRX(6)
              RDE = SQRT( AFZ(NSZ(N))/GPI/SRX(7) )
              DRD2 = (RDE**2-RDW**2)
              HCGX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
     &          (VISG(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
!
!---  Gas Volumetric Rate w/ Component Relative Humidity ---
!
              IF( ISRT(NS).EQ.23 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    PGWX = SRX(4)*WTX/SRX(7)/HCGX+PG(M,N)
                    PX = MAX( PGX,PLX,PGWX+PATM )
                    INDX = 0
                    CALL REGION_4( SRX(2),PSWX,INDX )
                    PVWX = PSWX*SRX(5)
                    PVAX = MAX( PX-PVWX,0.D+0 )
                    INDX = 0
                    CALL WATGSD( SRX(2),PVWX,RHOW,INDX )
                    CALL AIRGSD( SRX(2),PVAX,RHOA )
                    CALL WATGSH( SRX(2),PSWX,RHOW,HGWX,UGWX )
                    CALL AIRGSH( SRX(2),PVAX,HGAX,UGAX )
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*(RHOW*HGWX+RHOA*HGAX)
     &                *WTX
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOA*WTX
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOW*WTX
                  ELSE
                    DPGX = MIN( SRX(3)-PG(M,N),0.D+0 )
                    QGX = MAX( DPGX*RKG(M,N)*HCGX,WTX*SRX(4) )
                    HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                    SRCT(M,N) = SRCT(M,N) + QGX*RHOG(M,N)*HGX
                    SRCA(M,N) = SRCA(M,N) + QGX*RHOG(M,N)*XGA(M,N)
                    SRCW(M,N) = SRCW(M,N) + QGX*RHOG(M,N)*XGW(M,N)
                  ENDIF
                ENDIF
!
!---  Gas Volumetric Rate w/ Component Mass Fractions ---
!
                IF( ISRT(NS).EQ.24 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    PGWX = SRX(4)*WTX/SRX(7)/HCGX+PG(M,N)
                    PX = MAX( PGX,PLX,PGWX+PATM )
                    XGWX = SRX(5)
                    XGAX = MAX( 1.D+0-XGWX,0.D+0 )
                    XMGWX = (XGWX/WTMW)/((XGAX/WTMA)+(XGWX/WTMW))
                    INDX = 0
                    CALL REGION_4( SRX(2),PSWX,INDX )
                    PVWX = MIN( PX*XMGWX,PSWX )
                    PVAX = MAX( PX-PVWX,0.D+0 )
                    INDX = 0
                    CALL WATGSD( SRX(2),PVWX,RHOW,INDX )
                    CALL AIRGSD( SRX(2),PVAX,RHOA )
                    CALL WATGSH( SRX(2),PSWX,RHOW,HGWX,UGWX )
                    CALL AIRGSH( SRX(2),PVAX,HGAX,UGAX )
                    SRCT(M,N) = SRCT(M,N) + SRX(4)*(RHOW*HGWX+RHOA*HGAX)
     &                *WTX
                    SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOA*WTX
                    SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOW*WTX
                  ELSE
                    DPGX = MIN( SRX(3)-PG(M,N),0.D+0 )
                    QGX = MAX( DPGX*RKG(M,N)*HCGX,WTX*SRX(4) )
                    HGX = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
                    SRCT(M,N) = SRCT(M,N) + QGX*RHOG(M,N)*HGX
                    SRCA(M,N) = SRCA(M,N) + QGX*RHOG(M,N)*XGA(M,N)
                    SRCW(M,N) = SRCW(M,N) + QGX*RHOG(M,N)*XGW(M,N)
                  ENDIF
                ENDIF
  516       CONTINUE
  518     CONTINUE
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORC_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORIT_GT( NSL )
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
!     Compute solute transport source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE GRID
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
      REAL*8 SRX(8+LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORIT_GT'
!
!---  Loop over sources  ---
!
      DO NS = 1,NSR
        IF( TM.LE.SRC(1,1,NS) ) CYCLE
        SRX(1) = TM
        IF( ISRM(NS).EQ.1 ) THEN
          DO N = 2,8+NSOLU
            SRX(N) = SRC(N,1,NS)
          ENDDO
        ELSE
          DO M = 2,ISRM(NS)
            IF( TM.LE.SRC(1,M,NS) ) THEN
              DTSR = MIN( SRC(1,M,NS)-TM,DT )
              TFSR = (TM-0.5D+0*DTSR-SRC(1,M-1,NS))/
     &          (SRC(1,M,NS)-SRC(1,M-1,NS))
              DO N = 2,8+NSOLU
                SRX(N) = SRC(N,M-1,NS)+TFSR*(SRC(N,M,NS)-SRC(N,M-1,NS))
              ENDDO
             GOTO 110
            ENDIF
          ENDDO
          CYCLE
        ENDIF
  110   CONTINUE
!
!---    Loop over source domain  ---
!
        DO I = ISRDM(1,NS),ISRDM(2,NS)
        DO J = ISRDM(3,NS),ISRDM(4,NS)
        DO K = ISRDM(5,NS),ISRDM(6,NS)
          N = ND(I,J,K)
          IF( IXP(N).EQ.0 ) CYCLE
!
!---      Aqueous Volumetric Sink  ---
!
          IF( ISRT(NS).EQ.3 .AND. SRX(4).LT.0.D+0 ) THEN
            SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &      YL(N,NSL)*DT/(PORD(2,N)*SL(2,N))
!
!---      Gas Volumetric Sink  ---
!
          ELSEIF( (ISRT(NS).EQ.4 .OR. ISRT(NS).EQ.5)
     &      .AND. SRX(4).LT.0.D+0 ) THEN
            SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &      YG(N,NSL)*DT/(PORD(2,N)*SG(2,N))
!
!---      Aqueous Mass Sink  ---
!
          ELSEIF( ISRT(NS).EQ.7 .AND. SRX(4).LT.0.D+0 ) THEN
            SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &        YL(N,NSL)*DT/(RHOL(2,N)*PORD(2,N)*SL(2,N))
!
!---      Gas Mass Sink  ---
!
          ELSEIF( (ISRT(NS).EQ.8 .OR. ISRT(NS).EQ.9)
     &      .AND. SRX(4).LT.0.D+0 ) THEN
            SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &        YG(N,NSL)*DT/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---      Solute source  ---
!
          ELSEIF( ISRT(NS).LT.0 .AND. ISRT(NS).GE.-NSOLU ) THEN
            SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(4)*DT
!
!---      Solute source  ---
!
          ELSEIF( ISRT(NS).LT.-NSOLU .AND.
     &      ISRT(NS).GE.-2*NSOLU ) THEN
            SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(4)*DT*VOL(N)
!
!---      Diffusion-dominated solute release model  ---
!
          ELSEIF( ISRT(NS).EQ.-(NSL+5*NSOLU) .OR.
     &      ISRT(NS).EQ.-(NSL+8*NSOLU) ) THEN
            SRCIC(N,NSL) = SRCIC(N,NSL) + CNL(N,NSL)*DT
!
!---      Solubility-controlled solute release model or
!         solubility-controlled salt-cake release model  ---
!
          ELSEIF( ISRT(NS).EQ.-(NSL+6*NSOLU) .OR.
     &      ISRT(NS).EQ.-(NSL+7*NSOLU) ) THEN
!
!---        SRX(2): nodal solute inventory  ---
!
            IF ( CNL(N,NSL).GT.ZERO ) THEN
!
!---          Integrated solute mass has exceeded solute inventory
!             and is being truncated to solute inventory  ---
!
              SRCIC(N,NSL) = SRCIC(N,NSL) + CNL(N,NSL)*DT
            ELSEIF ( SRCIC(N,NSL).LT.SRX(2) ) THEN
!
!---          Integrated solute mass still below solute inventory  ---
!
              NPX = NSX(N)
              NPY = NSY(N)
              NPZ = NSZ(N)
              NQX = NPX + 1
              NQY = NPY + IFLD
              NQZ = NPZ + IJFLD
              DCDT =  (C(N,NSL)-CO(N,NSL))*VOL(N)*DTI
              CLX = ZERO
              IF ( I.NE.1 ) THEN
                CLX =  CLX + UC(NPX,NSL)*AFX(NPX)
              ENDIF
              IF ( I.NE.IFLD ) THEN
                CLX =  CLX - UC(NQX,NSL)*AFX(NQX)
              ENDIF
              IF ( J.NE.1 ) THEN
                CLX =  CLX + VC(NPY,NSL)*AFY(NPY)
              ENDIF
              IF ( J.NE.JFLD ) THEN
                CLX =  CLX - VC(NQY,NSL)*AFY(NQY)
              ENDIF
              IF ( K.NE.1 ) THEN
                CLX =  CLX + WC(NPZ,NSL)*AFZ(NPZ)
              ENDIF
              IF ( K.NE.KFLD ) THEN
                CLX =  CLX - WC(NQZ,NSL)*AFZ(NQZ)
              ENDIF
              SRCIC(N,NSL) = SRCIC(N,NSL) + (DCDT-CLX)*DT
            ENDIF
          ENDIF
        ENDDO
        ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORIT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_GT( NSL )
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
!     Compute solute transport source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!






!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE PORMED
      USE JACOB
      USE GRID
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





      REAL*8 SRX(8+LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORT_GT'
!
!---  Loop over sources  ---
!
      DO NS = 1,NSR
        IF( TM.LE.SRC(1,1,NS) ) CYCLE
        SRX(1) = TM
        IF( ISRM(NS).EQ.1 ) THEN
          DO N = 2,8+NSOLU
            SRX(N) = SRC(N,1,NS)
          ENDDO
          DTMAX = TM-SRC(1,1,NS)
        ELSE
          DO M = 2,ISRM(NS)
            IF( TM.LE.SRC(1,M,NS) ) THEN
              DTSR = MIN( SRC(1,M,NS)-TM,DT )
              TFSR = (TM-0.5D+0*DTSR-SRC(1,M-1,NS))/
     &          (SRC(1,M,NS)-SRC(1,M-1,NS))
              DO N = 2,8+NSOLU
                SRX(N) = SRC(N,M-1,NS)+TFSR*(SRC(N,M,NS)-SRC(N,M-1,NS))
              ENDDO
              DTMAX = TM-SRC(1,M-1,NS)
              GOTO 110
            ENDIF
          ENDDO
          CYCLE
        ENDIF
  110   CONTINUE








!
!---  Loop over source domain  ---
!
        DO I = ISRDM(1,NS),ISRDM(2,NS)
        DO J = ISRDM(3,NS),ISRDM(4,NS)
        DO K = ISRDM(5,NS),ISRDM(6,NS)
          N = ND(I,J,K)
          IF( IXP(N).EQ.0 ) CYCLE
          MP = IXP(N)
          IF( ILES.EQ.1 ) THEN
            MCOL = MP
            MROW = MDT
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1




          ENDIF
          IADDVAL = 1
          SORTX = 0.D+0
!
!---      Aqueous Volumetric Sink  ---
!
          IF( ISRT(NS).EQ.3 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YL(N,NSL)/(PORD(2,N)*SL(2,N))
!
!---      Gas Volumetric Sink  ---
!
          ELSEIF( ISRT(NS).EQ.4 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YG(N,NSL)/(PORD(2,N)*SG(2,N))
!
!---      Gas Volumetric Sink  ---
!
          ELSEIF( ISRT(NS).EQ.5 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YG(N,NSL)/(PORD(2,N)*SG(2,N))
!
!---      Aqueous Mass Sink  ---
!
          ELSEIF( ISRT(NS).EQ.7 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YL(N,NSL)/(RHOL(2,N)*PORD(2,N)*SL(2,N))
!
!---      Gas Mass Sink  ---
!
          ELSEIF( ISRT(NS).EQ.8 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YG(N,NSL)/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---      Gas Mass Sink  ---
!
          ELSEIF( ISRT(NS).EQ.9 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YG(N,NSL)/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---      Solute source  ---
!
          ELSEIF( ISRT(NS).EQ.-NSL ) THEN
            BLU(MP) = BLU(MP) + SRX(4)
!
!---      Solute density source  ---
!
          ELSEIF( ISRT(NS).EQ.-(NSL+NSOLU) ) THEN
            BLU(MP) = BLU(MP) + SRX(4)*VOL(N)
!
!---      Diffusion-dominated solute release model  ---
!
          ELSEIF( ISRT(NS).EQ.-(NSL+5*NSOLU) ) THEN
            QSRIX = 2.D+0*SRX(3)*SQRT(SRX(5)*DTMAX/GPI)/SRX(4)
            QSRIX = MIN( SRX(3),QSRIX )
            QSRRX = MAX( (QSRIX-SRCIC(N,NSL))*DTI,0.D+0 )
!
!---        Hold solute release rate in variable CNL  ---
!
            CNL(N,NSL) = QSRRX
            BLU(MP) = BLU(MP) + QSRRX
!
!---      Solubility-controlled solute release models or
!         solubility-controlled salt cake release models  ---
!
          ELSEIF( ISRT(NS).EQ.-(NSL+6*NSOLU) .OR.
     &      ISRT(NS).EQ.-(NSL+7*NSOLU) ) THEN
!
!---        SRX(2): nodal solute inventory  ---
!
            IF ( SRCIC(N,NSL).LT.SRX(2) ) THEN
!
!---          Forecast whether the cummulative mass dissolved will exceed
!             the nodal inventory within the current time step.  This is
!             an approximation since the fluxes were computed at the 
!             previous time step.  However, the calculation does provide
!             an estimate for minimizing the error due to overshooting the
!             mass dissolved.  ---
!
              NPX = NSX(N)
              NPY = NSY(N)
              NPZ = NSZ(N)
              NQX = NPX + 1
              NQY = NPY + IFLD
              NQZ = NPZ + IJFLD
              DCDT =  (C(N,NSL)-CO(N,NSL))*VOL(N)*DTI
              CLX = ZERO
              IF ( I.NE.1 ) THEN
                CLX =  CLX + UC(NPX,NSL)*AFX(NPX)
              ENDIF
              IF ( I.NE.IFLD ) THEN
                CLX =  CLX - UC(NQX,NSL)*AFX(NQX)
              ENDIF
              IF ( J.NE.1 ) THEN
                CLX =  CLX + VC(NPY,NSL)*AFY(NPY)
              ENDIF
              IF ( J.NE.JFLD ) THEN
                CLX =  CLX - VC(NQY,NSL)*AFY(NQY)
              ENDIF
              IF ( K.NE.1 ) THEN
                CLX =  CLX + WC(NPZ,NSL)*AFZ(NPZ)
              ENDIF
              IF ( K.NE.KFLD ) THEN
                CLX =  CLX - WC(NQZ,NSL)*AFZ(NQZ)
              ENDIF
              CNL(N,NSL) = ZERO
              IF ( SRCIC(N,NSL)+(DCDT-CLX)*DT.GT.SRX(2) ) THEN
!
!---            If the forecast predicts overshoot, set the source mass
!               as the difference between the current cummulative mass
!               and the nodal inventory.  ---
!
                QSRRX = (SRX(2)-SRCIC(N,NSL))*DTI
                CNL(N,NSL) = QSRRX
                BLU(MP) = BLU(MP) + QSRRX
              ELSE
!
!---            Otherwise, set diagonal matrix entry to 1 and all 
!               off-diagonal to 0, to force the concentration to be the 
!               solubility limit.  ---
!
                IF ( K.NE.1 ) THEN
                  IF( ILES.EQ.1 ) THEN
                    MCB = IXP(N-IJFLD)
                    MCOL = MCB
                    MROW = MP-MCB+MDT
                    ALU(MROW,MCOL) = 0.D+0
                  ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
                    MROW = KLUC(MP,MA)
                    DLU(MROW) = 0.D+0
                    MA = MA + 1







                  ENDIF
                ENDIF
                IF ( J.NE.1 ) THEN
                  IF( ILES.EQ.1 ) THEN
                    MCS = IXP(N-IFLD)
                    MCOL = MCS
                    MROW = MP-MCS+MDT
                    ALU(MROW,MCOL) = 0.D+0
                  ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
                    MROW = KLUC(MP,MA)
                    DLU(MROW) = 0.D+0
                    MA = MA + 1







                  ENDIF
                ENDIF
                IF( ILES.EQ.1 ) THEN
                  ALU(MROW,MCOL) = 1.D+0
                ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
                  DLU(MCOL) = 1.D+0






                ENDIF
                IF ( I.NE.1 ) THEN
                  IF( ILES.EQ.1 ) THEN
                    MCW = IXP(N-1)
                    MCOL = MCW
                    MROW = MP-MCW+MDT
                    ALU(MROW,MCOL) = 0.D+0
                  ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
                    MROW = KLUC(MP,MA)
                    DLU(MROW) = 0.D+0
                    MA = MA + 1







                  ENDIF
                ENDIF
                IF ( I.NE.IFLD ) THEN
                  IF( ILES.EQ.1 ) THEN
                    MCE = IXP(N+1)
                    MCOL = MCE
                    MROW = MP-MCE+MDT
                    ALU(MROW,MCOL) = ZERO
                  ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
                    MROW = KLUC(MP,MA)
                    DLU(MROW) = 0.D+0
                    MA = MA + 1







                  ENDIF
                ENDIF
                IF ( J.NE.JFLD ) THEN
                  IF( ILES.EQ.1 ) THEN
                    MCN = IXP(N+IFLD)
                    MCOL = MCN
                    MROW = MP-MCN+MDT
                    ALU(MROW,MCOL) = 0.D+0
                  ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
                    MROW = KLUC(MP,MA)
                    DLU(MROW) = 0.D+0
                    MA = MA + 1







                  ENDIF
                ENDIF
                IF ( K.NE.KFLD ) THEN
                  IF( ILES.EQ.1 ) THEN
                    MCT = IXP(N+IJFLD)
                    MCOL = MCT
                    MROW = MP-MCT+MDT
                    ALU(MROW,MCOL) = 0.D+0
                  ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
                    MROW = KLUC(MP,MA)
                    DLU(MROW) = 0.D+0
                    MA = MA + 1







                  ENDIF
                ENDIF
!
!---            Solubility-controlled solute release model  ---
!
                IF( ISRT(NS).EQ.-(NSL+6*NSOLU) ) THEN
!
!---              SRX(2): nodal solute inventory
!                 SRX(3): aqueous solubility  ---
!
                  BLU(MP) = SRX(3)*SL(2,N)*PORD(2,N)/YL(N,NSL)
!
!---              Solubility-controlled salt cake release model  ---
!
                ELSEIF ( ISRT(NS).EQ.-(NSL+7*NSOLU) ) THEN
!
!---              SRX(2): nodal solute inventory
!                 SRX(3): nodal salt cake inventory
!                 SRX(4): salt cake solubility  ---
!
                  BLU(MP) = SRX(4)*SRX(2)/SRX(3)*SL(2,N)*PORD(2,N)/
     &                    YL(N,NSL)
                ENDIF
              ENDIF
            ENDIF
            IADDVAL = 0
!
!---      Diffusion-dominated solute release model (w/variable diffusion)  ---
!
          ELSEIF( ISRT(NS).EQ.-(NSL+8*NSOLU) ) THEN
!
!---        SRX(3): nodal solute inventory
!           SRX(4): vertical depth of residual waste
!           SRX(5): diffusion coeficient
!           SRX(6): constrictivity  ---
!
            TORLX = TORL(2,N)
            DIFFIX = SRX(5)*PORD(2,N)*SRX(6)/(TORLX*TORLX)
            IZN = IZ(N)
            ALPHAX = PORD(2,N) + RHOS(IZN)*PCSL(1,IZN,NSL)*
     &        (1.D+0-PORT(2,N))
            DIFFAX = DIFFIX/ALPHAX
            QSRIX = 2.D+0*SRX(3)*SQRT(DIFFAX*DTMAX/GPI)/SRX(4)
            QSRIX = MIN( SRX(3),QSRIX )
            QSRRX = MAX( (QSRIX-SRCIC(N,NSL))*DTI,0.D+0 )
!
!---        Hold solute release rate in variable CNL  ---
!
            CNL(N,NSL) = QSRRX
            BLU(MP) = BLU(MP) + QSRRX
          ENDIF
!
!---      Load Jacobian  ---
!
          IF( IADDVAL.EQ.1 ) THEN
            IF( ILES.EQ.1 ) THEN
              ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
            ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
              DLU(MCOL) = DLU(MCOL) + SORTX





            ENDIF
          ENDIF
        ENDDO
        ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORT_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP_GT( NSL )
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
!     Calculates the aqueous- and gas-phase solute
!     mole fractions from user-specified partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE GRID
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
      SUB_LOG(ISUB_LOG) = '/SPRP_GT'
!
!---  Loop over all nodes  ---
!
      DO 900 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 900
        IZN = IZ(N)
        IF( IPCL(NSL).EQ.2 ) THEN
          XVS = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))*SL(2,N)
        ELSE
          XVS = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
        ENDIF
        XVL = SL(2,N)*PORD(2,N)
        XVG = SG(2,N)*PORD(2,N)
!
!---    Constant gas-aqueous partition coefficient  ---
!
        IF( IPCGL(NSL).EQ.0 ) THEN
          PCGLX = PCGL(1,NSL)
!
!---    Temperature dependent gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.1 ) THEN
          TK = T(2,N)+TABS
          PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &      + PCGL(3,NSL)*LOG(TK)
     &      + PCGL(4,NSL)*TK + PCGL(5,NSL)*TK**2 )
!
!---    Water-vapor gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.2 ) THEN
          PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
        ENDIF
        PCGLX = MAX( PCGLX,1.D-20 )
        PCGLX = MIN( PCGLX,1.D+20 )
!
!---    Phase-volumetric concentration ratios  ---
!
        YVL = 1.D+0/(XVS + XVL + XVG*PCGLX)
        YVG = 1.D+0/((XVS + XVL)/PCGLX + XVG)
!
!---    Phase mole fractions  ---
!
        YL(N,NSL) = XVL*YVL
        YG(N,NSL) = XVG*YVG
  900 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPRP_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TPORT_GT( NSL )
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
!     Solute/Reactive Species Transport Shell.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 September 2005.
!

!
!----------------------PETSc Modules-----------------------------------!
!
      USE STOMP_LIS_MODULE







!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/TPORT_GT'
!
!---  Zero Jacobian matrix  ---
!



      INDX = 6
      CALL JCBZ( ISVT,MUT,MLT,MKT,INDX )
!
!---  Compute solute sources ---
!
      CALL SORT_GT( NSL )
!
!---  Zero solute transport fluxes  ---
!
      CALL SFXZ( NSL )
!
!---  Load Jacobian matrix (aqueous-phase transport)  ---
!
      CALL SJCBL( NSL )
!
!---  Load Jacobian matrix (gas-phase transport)  ---
!
      CALL SJCBG( NSL )
!
!---  Modify Jacobian matrix for boundary conditions ---
!
      CALL SBND_GT( NSL )
!
!---  Fracture solute sources ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL SORT_FRC_GT( NSL )
!
!---  Borehole solute sources ---
!
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL SORT_BH_GT( NSL )
!
!---  Zero solute transport fluxes for fractures ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL SFXZ_FRC_GT( NSL )
!
!---  Zero solute transport fluxes for boreholes ---
!
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL SFXZ_BH_GT( NSL )
!
!---  Load Jacobian matrix (aqueous transport) for fractures  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL SJCBL_FRC_GT( NSL )
!
!---  Load Jacobian matrix (aqueous transport) for boreholes  ---
!
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL SJCBL_BH_GT( NSL )
!
!---  Load Jacobian matrix (gas transport) for fractures  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL SJCBG_FRC_GT( NSL )
!
!---  Load Jacobian matrix (gas transport) for boreholes  ---
!
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL SJCBG_BH_GT( NSL )
!
!---  Transfer functions between fracture and matrix, fracture
!     and matrix equations ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL TRNSC_FM_GT( NSL )
!
!---  Transfer functions between borehole and matrix, borehole
!     and matrix equations ---
!
      IF( ISLC(74).EQ.3 .OR. ISLC(74).EQ.3 ) CALL TRNSC_BM_GT( NSL )
!
!---  Transfer functions between borehole and fracture, borehole
!     and fracture equations ---
!
      IF( ISLC(74).EQ.3 ) CALL TRNSC_BF_GT( NSL )
!
!---  Linear equation solver  ---
!
      IF( ILES.EQ.1 ) THEN
        INDX = 1
        CALL BAND( 0,MUT,MLT,INDX )
      ELSEIF( ILES.EQ.3 ) THEN
        INDX = 1
        CALL PSPLIB( 0,INDX )

      ELSEIF( ILES.EQ.4 ) THEN
        INDX = 1
        CALL STOMP_LIS_SOLVE(-1,T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,INDX)







      ENDIF
!
!---  Update solute concentrations ---
!
      CALL UPDTC( NSL )
!
!---  Update solute concentrations for fractures  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL UPDTC_FRC( NSL )
!
!---  Update solute concentrations for boreholes  ---
!
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL UPDTC_BH( NSL )
!
!---  Compute solute aqueous-phase fluxes (interior nodes)  ---
!
      CALL SFXL( NSL )
!
!---  Compute solute aqueous-phase fluxes (boundary surfaces)  ---
!
      CALL SFXLB( NSL )
!
!---  Compute solute gas-phase fluxes (interior nodes)  ---
!
      CALL SFXG( NSL )
!
!---  Compute solute gas-phase fluxes (boundary surfaces)  ---
!
      CALL SFXGB( NSL )
!
!---  Integrate solute sources  ---
!
      CALL SORIT_GT( NSL )
!
!---  Solute flux for fracture triangle to fracture triangle flow  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL SFX_FRC_GT( NSL )
!
!---  Solute flux for borehole node to borehole node flow  ---
!
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL SFX_BH_GT( NSL )
!
!---  Solute flux for fracture triangle to matrix grid cell flow  ---
!
      IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 ) CALL SFX_FM_GT( NSL )
!
!---  Solute flux for borehole node to matrix grid cell flow  ---
!
      IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 ) CALL SFX_BM_GT( NSL )
!
!---  Solute flux for borehole node to fracture triangle flow  ---
!
      IF( ISLC(74).EQ.3 ) CALL SFX_BF_GT( NSL )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TPORT_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_GT
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
!     Update the primary variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_GT'
      IF( ICNV.EQ.1 ) GOTO 300
      IERR = 0
!
!---  Update primary variables
!
      DO N = 1,NFLD
        N_DB = N
        IZN = IZ(N)
        IF( IERR.NE.0 ) CYCLE
!
!---    Skip inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
        NMD = IXP(N)
        MPL = IM(IEQW,NMD)
        DPL = BLU(MPL)
!
!---    Isoair option  ---
!
        IF( ISLC(37).EQ.0 ) THEN
          MPG = IM(IEQA,NMD)
          DPG = BLU(MPG)
        ELSE
          DPG = 0.D+0
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          MPS = IM(IEQS,NMD)
          DPS = BLU(MPS)
        ELSE
          DPS = 0.D+0
        ENDIF
!
!---    Nonisothermal simulations  ---
!
        IF( ISLC(30).EQ.0 ) THEN
          MPT = IM(IEQT,NMD)
          DPT = BLU(MPT)
          DPX = 2.5D+0
          DPT = SIGN( MIN( DPX,ABS(DPT) ),DPT )
          T(2,N) = T(2,N) + DPT
        ENDIF
!
!---    Saturated system w/o entrapped gas
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - gas air partial pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          INDX = 0
          CALL REGION_4( T(2,N),PSWX,INDX )
!          DPX = 2.5D-1*(PL(2,N)+PATM)
          DPX = 1.D+6
          DPL = SIGN( MIN( ABS(DPL),DPX ),DPL )
          PL(2,N) = PL(2,N) + DPL
!
!---      Zero negative corrections for zero aqueous air  ---
!
          IF( ISLC(37).EQ.0 ) THEN
!            IF( XMLA(2,N)/EPSL.LT.EPSL.AND.BLU(MPG)/EPSL.LT.EPSL ) THEN
!              BLU(MPG) = 0.D+0
!              DPG = 0.D+0
!            ENDIF
!            XMLA(2,N) = XMLA(2,N) + DPG
!            IF( XMLA(2,N).LT.1.D-12 ) XMLA(2,N) = 0.D+0
            PVA(2,N) = MAX( PVA(2,N)+DPG,0.D+0 )
            IF( PVA(2,N).LT.EPSL ) PVA(2,N) = 0.D+0
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Limit salt mass fraction changes to 0.25 of the
!           maximum value if salt mass fraction is less than
!           the maximum   ---
!
            PWX = MAX( PSWX,PL(2,N)+PATM )
            CALL SOL_BRNS( T(2,N),PWX,XLSMX )
            IF( YLS(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero aqueous salt  ---
!
            IF( YLS(2,N)/EPSL.LT.EPSL .AND. BLU(MPS)/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS(2,N) = YLS(2,N)+DPS
            IF( YLS(2,N).LT.EPSL ) YLS(2,N) = 0.D+0
            XLS(2,N) = MIN( YLS(2,N),XLSMX )
          ENDIF
!
!---    Unsaturated system w/ or w/o entrapped gas
!     
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - gas pressure
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          INDX = 0
          CALL REGION_4( T(2,N),PSWX,INDX )
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Limit salt mass fraction changes to 0.25 of the
!           maximum value if salt mass fraction is less than
!           the maximum   ---
!
            CALL SOL_BRNS( T(2,N),PSWX,XLSMX )
            IF( YLS(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero aqueous salt  ---
!
            IF( YLS(2,N)/EPSL.LT.EPSL .AND. BLU(MPS)/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS(2,N) = YLS(2,N)+DPS
            IF( YLS(2,N).LT.EPSL ) YLS(2,N) = 0.D+0
            XLS(2,N) = MIN( YLS(2,N),XLSMX )
          ENDIF
!
!---      Assign gas-entry pressure for non Brooks-Corey;
!         Brooks-Corey; Brooks-Corey, Dual Porosity; and
!         Brooks-Corey, Entrapment  ---
!
          IF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.102 .OR.
     &        ISCHR(IZN).EQ.202 ) THEN
            ENPR = SCHRV(1,N)*RHORL*GRAV
          ELSEIF( ISCHR(IZN).EQ.4 ) THEN
            ENPR = MIN( SCHRV(1,N),SCHR(5,IZN) )*RHORL*GRAV
          ELSE
            ENPR = 0.D+0
          ENDIF
!!
!!---      Limit changes in pressure  ---
!!
!          DPX = MAX( 1.D+6,2.5D-1*(PG(2,N)-PL(2,N)) )
!          DPCX = MAX( 1.D+4,2.5D-1*(PG(2,N)-PL(2,N)) )
!          IF( (DPG-DPL).GT.DPCX ) THEN
!            DPG = DPG - 5.D-1*DPCX
!            DPL = DPL + 5.D-1*DPCX
!          ELSE
!            DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
!            DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
!          ENDIF
!
!---      Limit changes in pressure, and trap large changes in 
!         gas pressure  ---
!
          DPX = MAX( 1.D+6,2.5D-1*(PG(2,N)-PL(2,N)) )
          DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
          IF( ISLC(37).EQ.0 ) THEN
            DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
!
!---       Relax pressure updates when transitioning to unsaturated
!           conditions  ---
!
            IF( (PG(2,N)+DPG)-(PL(2,N)+DPL).LT.ENPR ) THEN
              DPG = 6.D-1*DPG
              DPL = 6.D-1*DPL
            ENDIF
            PG(2,N) = PG(2,N) + DPG
            PG(2,N) = MIN( PG(2,N),5.D+8 )
          ENDIF
          PL(2,N) = PL(2,N) + DPL
          PL(2,N) = MAX( PL(2,N),(PG(2,N)-SCHR(12,IZN)*RHORL*GRAV) )
!
!---      Maintain the gas pressure above or at the water vapor
!         pressure  ---
!
          CALL P_IAPWS( T(2,N),PSWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL DENS_B( XLS(2,N),RHOBX,RHOLWX,T(2,N) )
          CALL SP_B( T(2,N),XLS(2,N),PSBX )
          PGX = PG(2,N) + PATM
          PLX = PL(2,N) + PATM
          PCX = PGX - PLX
!          PCX = MIN( PCX,PCMX )
          CALL VPL_B( T(2,N),PSBX,PCX,RHOBX,PVBX )
          PG(2,N) = MAX( PG(2,N),(PVBX-PATM) )
!!
!!---      Maintain the aqueous pressure within a maximum
!!         capillary pressure, based on a 0.05 change in aqueous
!!         saturation  ---
!!
!          PL(2,N) = MAX( PL(2,N),PG(2,N)-PCMX )
!
!---    Saturated system w/ entrapped gas
!
!       Energy - temperature
!       Water mass - aqueous pressure
!       Air mass - trapped gas saturation
!       NaCl mass - total NaCl brine mass fraction  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
!
!---      Limit changes in pressure  ---
!
          DPX = 2.5D-1*MAX(PL(2,N)+PATM,1.D+6)
          DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
          PL(2,N) = PL(2,N) + DPL
          PL(2,N) = MIN( PL(2,N),5.D+8 )
!
!---      Limit changes in trapped gas  ---
!
          IF( ISLC(37).EQ.0 ) THEN
            DPX = 1.D-1*SCHR(15,IZ(N))
            DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
            SG(2,N) = SG(2,N) + DPG
            IF( SG(2,N).LT.EPSL ) SG(2,N) = 0.D+0
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Limit salt mass fraction changes to 0.25 of the
!           maximum value if salt mass fraction is less than
!           the maximum   ---
!
            INDX = 0
            CALL REGION_4( T(2,N),PWX,INDX )
            CALL P_IAPWS( T(2,N),PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
            CALL SOL_BRNS( T(2,N),PWX,XLSMX )
            IF( YLS(2,N).LT.XLSMX ) THEN
              DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
            ENDIF
!
!---        Zero negative corrections for zero aqueous salt  ---
!
            IF( YLS(2,N)/EPSL.LT.EPSL .AND. BLU(MPS)/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            YLS(2,N) = YLS(2,N)+DPS
            IF( YLS(2,N).LT.EPSL ) YLS(2,N) = 0.D+0
            XLS(2,N) = MIN( YLS(2,N),XLSMX )
          ENDIF
!
!---    Fully unsaturated conditions
!
!       Energy - temperature
!       Water mass - water vapor partial pressure
!       Air mass - gas pressure
!       NaCl mass - salt mass  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
!
!---      Limit changes in water vapor partial pressure  ---
!
          INDX = 0
          CALL REGION_4( T(2,N),PWX,INDX )
          PWX = MIN( PWX,PCRW )
          DPX = 2.5D-1*PWX
          DPL = SIGN( MIN(ABS(DPX),ABS(DPL)),DPL )
          PVW(2,N) = MAX( PVW(2,N)+DPL,1.D-6 )
!
!---      Limit changes in gas pressure  ---
!
          IF( ISLC(37).EQ.0 ) THEN
            DPX = 2.5D-1*MAX(PVA(2,N),1.D+6)
            DPG = SIGN( MIN(ABS(DPX),ABS(DPG)),DPG )
            PVA(2,N) = MAX( PVA(2,N)+DPG,0.D+0 )
            PG(2,N) = PVW(2,N) + PVA(2,N) - PATM
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Zero negative corrections for salt volumetric conc.  ---
!
            IF( TMS(2,N)/EPSL.LT.EPSL .AND. BLU(MPS)/EPSL.LT.EPSL ) THEN
              BLU(MPS) = 0.D+0
              DPS = 0.D+0
            ENDIF
            TMS(2,N) = TMS(2,N) + DPS
            IF( TMS(2,N).LT.1.D-9 ) TMS(2,N) = 0.D+0
          ENDIF
        ENDIF 
!
!---    Check for excessive pressure or temperature   ---
!
        PGX = PG(2,N)+PATM
        PLX = PL(2,N)+PATM
        TKX = T(2,N)+TABS
        IF( IERR.EQ.0 ) THEN
          IF( PGX.GT.100.D+6 .OR. PGX.LT.1.001D-6 ) IERR = N
          IF( PLX.GT.100.D+6 ) IERR = N
          IF( TKX.GT.1473.15D+0 .OR. TKX.LT.TABS ) IERR = N
        ENDIF
        ENDDO
!
!---  Reduce time step for excessive changes in primary variables   ---
!
      IF( IERR.EQ.1 ) THEN
        N = IERR
        ICNV = 1
        WRITE(ISC,'(10X,A)') '---  Excessive Primary Variable Change
     &---'
        WRITE(IWR,'(10X,A)') '---  Excessive Primary Variable Change
     &---'
        WRITE(ISC,'(4X,A,I6)') 'Node = ',N
        WRITE(IWR,'(4X,A,I6)') 'Node = ',N
        WRITE(ISC,'(4X,2A)') 'Phase Condition = ',PH_CND(NPHAZ(2,N))
        WRITE(IWR,'(4X,2A)') 'Phase Condition = ',PH_CND(NPHAZ(2,N))
        WRITE(ISC,'(4X,A,1PE12.5)') 'Temperature, C = ',T(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Temperature, C = ',T(2,N)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Aqu. Pressure, Pa = ',PL(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Aqu. Pressure, Pa = ',PL(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Pressure, Pa = ',PG(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Pressure, Pa = ',PG(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Gas Saturation = ',SG(2,N)
        WRITE(ISC,'(4X,A,1PE12.5)')
     &    'Aqueous-Air Mole Fraction = ',XMLA(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)')
     &    'Aqueous-Air Mole Fraction = ',XMLA(2,N)
        WRITE(ISC,'(4X,A,1PE12.5)')
     &    'Water Vapor Pressure = ',PVW(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)')
     &    'Water Vapor Pressure = ',PVW(2,N)
        WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &    'Total-Salt Aqu. Mass Fraction = ',YLS(2,N)
        WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &    'Total-Salt Aqu. Mass Fraction = ',YLS(2,N)
        WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &    'Salt Volumetric Concentration = ',TMS(2,N)
        WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &    'Salt Volumetric Concentration = ',TMS(2,N)
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
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
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
!---      Coupled-well pressure  ---
!
          IF( L_CW.EQ.1 ) THEN
            DO NCW = 1,N_CW
              P_CW(2,NCW) = P_CW(1,NCW)
            ENDDO
          ENDIF
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
!---  End of UPDT_GT group
!
      RETURN
      END


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
!----------------------------------------------------------------------!
!
!     STOMP: Subsurface Transport Over Multiple Phases
!
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     This engineering program numerically simulates the transport
!     of H2O, NaCl, CH4 and O2 through multifluid subsurface 
!     environments under non-isothermal conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
#ifdef lis
!
!----------------------Lis Modules-----------------------------------!
!
      USE STOMP_LIS_MODULE
#endif
#ifdef petsc
!
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE OUTPU
      USE JACOB
      USE GEO_MECH
      USE FILES
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
#ifdef lis
#include "include/lisf.h"
#endif
!
!----------------------Type Declarations-------------------------------!
!
      LOGICAL HALT,PLOT,RESTART
#ifdef petsc
      PetscErrorCode :: IERR
#endif
#ifdef lis
      LIS_INTEGER :: IERR
#endif
!
!----------------------Executable Lines--------------------------------!
!
#if defined(intel) || defined(gnu)
!
!---  Initialize system clock ---
!
      CALL SYSTEM_CLOCK(ICOUNT1,ICRATE,ICMAX)
#endif
#ifdef petsc
!
!---  Initialize PETSc ---
!
      CALL PetscInitialize(PETSC_NULL_CHARACTER,IERR)
#endif
#ifdef lis
!
!---  Initialize Lis ---
!
      CALL lis_initialize(IERR)
#endif
!
!---  Read input file to determine memory requirements  ---
!
      CALL STEP
!
!---  Allocate memory  ---
!
      CALL ALLOC
#ifdef demo
!
!---  Number of active nodes exceeds demonstration limit  ---
!
      IF( LAN.GT.1000 ) THEN
        INDX = 3
        CHMSG = 'Number of Active Nodes > Demonstration Limit of 1000'
        CALL WRMSGS( INDX )
      ENDIF
#endif
!
!---  Number of grid cells exceeds binary file size limit  ---
!
      IF( ISLC(67).EQ.1 .AND. LFLD.GT.12000000 ) THEN
        INDX = 3
        CHMSG = 'Number of Grid Cells Yields Binary Files > 2.4 GB'
        CALL WRMSGS( INDX )
      ENDIF
      ISUB_LOG = 1
      SUB_LOG(1) = 'STOMP-HYDT-KE'
      ICODE = 39
!
!---  Intialize variables in common blocks and open files  ---
!
      CALL INTLZ
!
!---  Print banner on screen and output file  ---
!
      CALL BANNER
!
!---  Read user input and restart files  ---
!
      CALL RDINPT_HYD_O2_CH4_N2
!
!---  Create a node connection map  ---
!
      CALL CONNMAP
!
!---  Check for internal boundary surfaces and write connectivity
!     list file  --
!
      CALL CONNLST
!
!---  For geomechanics simulations create a finite-element node map  --
!
      IF( ISLC(50).NE.0 ) CALL CONNFEN
!
!---  For geomechanics simulations check and preprocess boundary
!     conditions, and set the reference volumetric stress from
!     the initial displacements stored in the restart file  ---
!
      IF( ISLC(50).NE.0 ) CALL CHK_GM
!
!---  Check thermodynamic and hydrologic initial states  ---
!
      CALL CHK_HYD_O2_CH4_N2
!
!---  Write bin input files for parallel processing  ---
!
      IF( ISLC(67).EQ.1 ) THEN
        CALL INCRM_HYD_O2_CH4_N2
        CALL PROP_HYD_O2_CH4_N2
!
!---    For geomechanics simulations compute Jacobian 
!       matrix pointers  --
!
        IF( ISLC(50).NE.0 ) CALL JCBP_GM_PPC
#ifdef ecke
!
!---    Sequence reaction equations  ---
!
        IF( ISLC(40).EQ.1 ) CALL SEQEQ
#endif
        CALL WRITE_BIN_HYD_O2_CH4_N2
        WRITE(IWR,'(A)') 'NOTE: Preprocessing for STOMPX-HYD-O2-CH4-N2'
        WRITE(ISC,'(A)') 'NOTE: Preprocessing for STOMPX-HYD-O2-CH4-N2'
!
!---    End parallel pre-processing unless output of the flow or
!       transport linear system is requested  ---
!
        IF( ISLC(34).EQ.0 ) THEN
          WRITE(IWR,'(/,A,/)') '---  End of STOMPX-HYD-O2-CH4-N2' //
     &      ' Preprocessing  ---'
          WRITE(ISC,'(/,A,/)') '---  End of STOMPX-HYD-O2-CH4-N2' //
     &      ' Preprocessing  ---'
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
      CALL INCRM_HYD_O2_CH4_N2
!
!---  Initial hydrologic and thermodynamic properties on interior nodes
!     and boundary surfaces  ---
!
      CALL PROP_HYD_O2_CH4_N2
      CALL BCP_HYD_O2_CH4_N2
!
!---  Compute initial solute concentrations  ---
!
      CALL CISC_HYD_O2_CH4_N2
!
!---  Initialize SPLIB routines  ---
!
      IF( ILES.EQ.3 ) THEN
        INDX = -1
        CALL PSPLIB( 0,INDX )
      ENDIF
#ifdef lis
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
#ifdef ecke
        NSL = NEQ + NSOLU
#else
        NSL = NSOLU
#endif
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
#endif
#ifdef petsc
!
!---  Initialize PETSc routines  ---
!
      IF( ILES.EQ.5 ) THEN
!
!---    Solver for coupled flow  ---
!
        INDX = 0
        CALL PETSC_SOLVER_INIT(ISVC,F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,
     &    NLU,INDX)
!
!---    Solver for solute/species transport  ---
!
#ifdef ecke
        NSL = NEQ + NSOLU
#else
        NSL = NSOLU
#endif
        IF( NSL.GT.0 ) THEN
          INDX = 1
          CALL PETSC_SOLVER_INIT(0,T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,
     &      NLUC,INDX)
        ENDIF
!
!---    Solver for geomechanics  ---
!
        IF( ISLC(50).NE.0 ) THEN
          INDX = 2
          CALL PETSC_SOLVER_INIT(0,G_KSP,G_MAT,G_RHS_VEC,G_SOL_VEC,
     &      NLU_GM,INDX)
        ENDIF
      ENDIF
#endif
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
!---  Compute initial fluxes on non-boundary and boundary surfaces  ---
!
      ISVF = 1
!
!---  Aqueous volumetric flux (non-boundary surfaces)  ---
!
      CALL DRCVL
!
!---  Gas volumetric flux (non-boundary surfaces)  ---
!
      CALL DRCVG
!
!---  Mobile O2 diffusive flux through gas (non-boundary surfaces)
!     inactive for iso-O2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFGA
!
!---  Mobile CH4 diffusive flux through gas (non-boundary surfaces)
!     inactive for iso-CH4 option  ---
!
      IF( ISLC(64).EQ.0 ) CALL DFFGO
!
!---  Water diffusive flux through gas (non-boundary surfaces)  ---
!
      CALL DFFGW
!
!---  Mobile N2 diffusive flux through gas (non-boundary surfaces)
!     inactive for iso-N2 option (last in gas diffusion order as
!     molar N2 diffusion in gas phase equal to minus the sum of the 
!     molar diffusion of H2O, CH4, and O2)  ---
!
      IF( ISLC(65).EQ.0 ) CALL DFFGN
!
!---  Mobile O2 diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-O2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFLA
!
!---  Mobile CH4 diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-Ch4 option  ---
!
      IF( ISLC(64).EQ.0 ) CALL DFFLO
!
!---  Mobile N2 diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-N2 option  ---
!
      IF( ISLC(65).EQ.0 ) CALL DFFLN
!
!---  Salt diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-brine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS
!
!---  Water diffusive flux through aqueous (non-boundary surfaces)
!     (last in aqueous diffusion order as molar H2O diffusion 
!     in aqueous phase equal to minus the sum of the 
!     molar diffusion of salt, N2, CH4, and O2)  ---
!
      CALL DFFLW
!
!---  Thermal diffusive flux (non-boundary surfaces)  ---
!
      CALL THD_HYDT_KE
!
!---  Thermal advective flux for aqueous (non-boundary surfaces)  ---
!
      CALL THL_HYDT_KE
!
!---  Thermal advective flux for gas (non-boundary surfaces)  ---
!
      CALL THG_HYDT_KE
!
!---  Volumetric fluxes, component diffusive fluxes, and thermal
!     fluxes (boundary surfaces)  ---
!
      CALL BCF_HYD_O2_CH4_N2
      ISVF = 2*ISVC+1
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
      CALL LDO_HYD_O2_CH4_N2
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
#if defined(intel) || defined(gnu)
!
!---  Stop simulation if system clock time is exceeded
!
      IF( CLKMX .GT. 1.D-6 ) THEN
        CALL SYSTEM_CLOCK(ICOUNT2)
        CLOCK_TIME = FLOAT(ICOUNT2-ICOUNT1)/FLOAT(ICRATE)
        IF( CLOCK_TIME.GT.CLKMX ) THEN
          INDX = 1
          CHMSG = 'Simulation Stopped:  Clock Time Limit'
          CALL WRMSGS( INDX )
          GOTO 900
        ENDIF
      ENDIF
#endif
!
!---  Generate plot file if file "plot_stomp" exists  ---
!
      INQUIRE( FILE="plot_stomp", EXIST=PLOT )
      IF( PLOT ) THEN
        OPEN( UNIT=19, FILE="plot_stomp" )
        CLOSE( UNIT=19, STATUS='DELETE' )
        CALL WRPLOT
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
        IF( ISLC(18).LT.1 ) CALL WRRST
      ENDIF
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
        CALL BCP_HYD_O2_CH4_N2
        GOTO 600
      ENDIF
      NTSR = 0
!
!---  Top of sequential flow and transport and geomechanics  ---
!
      K_GM(1) = 0
      K_GM(2) = 0
  190 CONTINUE
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
!---  Compute boundary saturation, relative permeability, and
!     thermodynamic properties  ---
!
      CALL BCP_HYD_O2_CH4_N2
!
!---  Compute source contributions  ---
!
      CALL SORC_HYD_O2_CH4_N2
!
!---  Aqueous volumetric flux (non-boundary surfaces)  ---
!
      CALL DRCVL
!
!---  Gas volumetric flux (non-boundary surfaces)  ---
!
      CALL DRCVG
!
!---  Mobile O2 diffusive flux through gas (non-boundary surfaces)
!     inactive for iso-O2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFGA
!
!---  Mobile CH4 diffusive flux through gas (non-boundary surfaces)
!     inactive for iso-CH4 option  ---
!
      IF( ISLC(64).EQ.0 ) CALL DFFGO
!
!---  Water diffusive flux through gas (non-boundary surfaces)  ---
!
      CALL DFFGW
!
!---  Mobile N2 diffusive flux through gas (non-boundary surfaces)
!     inactive for iso-N2 option (last in gas diffusion order as
!     molar N2 diffusion in gas phase equal to minus the sum of the 
!     molar diffusion of H2O, CH4, and O2)  ---
!
      IF( ISLC(65).EQ.0 ) CALL DFFGN
!
!---  Mobile O2 diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-O2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFLA
!
!---  Mobile CH4 diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-Ch4 option  ---
!
      IF( ISLC(64).EQ.0 ) CALL DFFLO
!
!---  Mobile N2 diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-N2 option  ---
!
      IF( ISLC(65).EQ.0 ) CALL DFFLN
!
!---  Salt diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-brine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS
!
!---  Water diffusive flux through aqueous (non-boundary surfaces)
!     (last in aqueous diffusion order as molar H2O diffusion 
!     in aqueous phase equal to minus the sum of the 
!     molar diffusion of salt, N2, CH4, and O2)  ---
!
      CALL DFFLW
!
!---  Thermal diffusive flux (non-boundary surfaces)  ---
!
      CALL THD_HYDT_KE
!
!---  Thermal advective flux for aqueous (non-boundary surfaces)  ---
!
      CALL THL_HYDT_KE
!
!---  Thermal advective flux for gas (non-boundary surfaces)  ---
!
      CALL THG_HYDT_KE
!
!---  Volumetric fluxes, component diffusive fluxes, and thermal
!     fluxes (boundary surfaces)  ---
!
      CALL BCF_HYD_O2_CH4_N2
!
!---  Zero Jacobian matrix  ---
!
#ifdef petsc
      IF( ILES.EQ.5 ) CALL MatZeroEntries(F_MAT,IERR)
#endif
      INDX = 0
      CALL JCBZ( ISVC,MUC,MLC,MKC,INDX )
!
!---  Load Jacobian matrix for the energy equation
!     (zero flux boundary)  ---
!
      CALL JCBT_HYDT_KE
!
!---  Load Jacobian matrix for the water equation
!     (zero flux boundary)  ---
!
      CALL JCBW_HYDT_KE
!
!---  Load Jacobian matrix for the mobile-O2 equation
!     (zero flux boundary)  ---
!
      IF( ISLC(45).EQ.0 ) CALL JCBA_HYDT_KE
!
!---  Load Jacobian matrix for the mobile-CH4 equation
!     (zero flux boundary)  ---
!
      IF( ISLC(64).EQ.0 ) CALL JCBO_HYDT_KE
!
!---  Load Jacobian matrix for the mobile-N2 equation
!     (zero flux boundary)  ---
!
      IF( ISLC(65).EQ.0 ) CALL JCBN_HYDT_KE
!
!---  Load Jacobian matrix for the hydrate-O2 equation
!     (zero flux boundary)  ---
!
      IF( ISLC(45).EQ.0 ) CALL JCBHA_HYDT_KE
!
!---  Load Jacobian matrix for the hydrate-CH4 equation
!     (zero flux boundary)  ---
!
      IF( ISLC(64).EQ.0 ) CALL JCBHO_HYDT_KE
!
!---  Load Jacobian matrix for the hydrate-N2 equation
!     (zero flux boundary)  ---
!
      IF( ISLC(65).EQ.0 ) CALL JCBHN_HYDT_KE
!
!---  Load Jacobian matrix for the salt equation
!     (zero flux boundary), isobrine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL JCBS_HYDT_KE
!
!---  Modify the Jacobian matrix for boundary conditions  ---
!
      CALL BCJ_HYD_O2_CH4_N2
!
!---  Linear equation solver  ---
!
      IF( ILES.EQ.1 ) THEN
        INDX = 0
        CALL BAND( ISVC,MUC,MLC,INDX )
      ELSEIF( ILES.EQ.3 ) THEN
        INDX = 0
        CALL PSPLIB( ISVC,INDX )
#ifdef lis
      ELSEIF( ILES.EQ.4 ) THEN
        INDX = 0
        CALL STOMP_LIS_SOLVE(ISVC,F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC,INDX)
#endif
#ifdef petsc
      ELSEIF( ILES.EQ.5 ) THEN
        CALL PETSC_PUT_RHS(BLU,F_RHS_VEC)
        CALL PETSC_SOLVE(F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC)
        CALL PETSC_GET_SOLUTION(BLU,F_SOL_VEC)
#endif
      ENDIF
!
!---  Update primary variables  ---
!
      CALL UPDT_HYD_O2_CH4_N2
!
!---  Compute convergence from maximum relative residuals  ---
!
      CALL RSDL_HYD_O2_CH4_N2
!
!---  Compute primary variable increments, saturation,
!     relative permeability, porosity, tortuosity,
!     thermodynamic properties for interior nodes,
!     except immediately after a new time step  ---
!
      CALL INCRM_HYD_O2_CH4_N2
      CALL PROP_HYD_O2_CH4_N2
!
!---  For geomechanics simulations alter permeability with
!     porosity  --
!
      IF( ISLC(50).NE.0 ) CALL PERMRF_GM
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
        CALL VOLSS_GM( INDX )
!
!---    Static porothermoelastic geomechanics  ---
!
        CALL STATIC_GM
!
!---    Convergence check for sequential coupled flow and transport
!       and geomechanics  ---
!
        CALL RSDL_GM
        IF( RSD_GM.GT.RSDM_GM(IEPD) .AND. K_GM(1).LT.8 ) THEN
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
          GOTO 190
        ENDIF
!
!---    Update porosity and permeability for geomechical stress  ---
!
        CALL PORSTY_GM
        CALL PERMRF_GM
      ENDIF
!
!---  Update fluxes for transport solutions or flux integrations  ---
!
      ISVF = 1
!
!---  Aqueous volumetric flux (non-boundary surfaces)  ---
!
      CALL DRCVL
!
!---  Gas volumetric flux (non-boundary surfaces)  ---
!
      CALL DRCVG
!
!---  Mobile O2 diffusive flux through gas (non-boundary surfaces)
!     inactive for iso-O2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFGA
!
!---  Mobile CH4 diffusive flux through gas (non-boundary surfaces)
!     inactive for iso-CH4 option  ---
!
      IF( ISLC(64).EQ.0 ) CALL DFFGO
!
!---  Water diffusive flux through gas (non-boundary surfaces)  ---
!
      CALL DFFGW
!
!---  Mobile N2 diffusive flux through gas (non-boundary surfaces)
!     inactive for iso-N2 option (last in gas diffusion order as
!     molar N2 diffusion in gas phase equal to minus the sum of the 
!     molar diffusion of H2O, CH4, and O2)  ---
!
      IF( ISLC(65).EQ.0 ) CALL DFFGN
!
!---  Mobile O2 diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-O2 option  ---
!
      IF( ISLC(45).EQ.0 ) CALL DFFLA
!
!---  Mobile CH4 diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-Ch4 option  ---
!
      IF( ISLC(64).EQ.0 ) CALL DFFLO
!
!---  Mobile N2 diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-N2 option  ---
!
      IF( ISLC(65).EQ.0 ) CALL DFFLN
!
!---  Salt diffusive flux through aqueous (non-boundary surfaces)
!     inactive for iso-brine option  ---
!
      IF( ISLC(32).EQ.0 ) CALL DFFLS
!
!---  Water diffusive flux through aqueous (non-boundary surfaces)
!     (last in aqueous diffusion order as molar H2O diffusion 
!     in aqueous phase equal to minus the sum of the 
!     molar diffusion of salt, N2, CH4, and O2)  ---
!
      CALL DFFLW
!
!---  Thermal diffusive flux (non-boundary surfaces)  ---
!
      CALL THD_HYDT_KE
!
!---  Thermal advective flux for aqueous (non-boundary surfaces)  ---
!
      CALL THL_HYDT_KE
!
!---  Thermal advective flux for gas (non-boundary surfaces)  ---
!
      CALL THG_HYDT_KE
!
!---  Volumetric fluxes, component diffusive fluxes, and thermal
!     fluxes (boundary surfaces)  ---
!
      CALL BCF_HYD_O2_CH4_N2
!
!---  Compute Local Courant Numbers  ---
!
      IF( ICRNT.EQ.1 ) CALL CRNTNB
      ISVF = 2*ISVC+1
!
!---  Beginning of transport equation solution  ---
!
      IF( IEQC.EQ.0 ) GOTO 800
!
!---  Loop over number of solutes  ---
!
      DO 700 NSL = 1,NSOLU
!
!---  Courant number limiting  ---
!
        N_CRN(NSL) = 1
        IF( ISLC(17).NE.0 ) CALL CRN_LIM( NSL )
        DO 690 NC = 1,N_CRN(NSL)
          IF( ISLC(17).NE.0 ) TM = MIN( TM+DT,TM_CRN )
!
!---      Compute solute mole fractions ---
!
          CALL SPRP_HYD_O2_CH4_N2( NSL )
!
!---      Zero Jacobian matrix  ---
!
#ifdef petsc
          IF( ILES.EQ.5 ) CALL MatZeroEntries(T_MAT,IERR)
#endif
          INDX = 1
          CALL JCBZ( ISVT,MUT,MLT,MKT,INDX )
!
!---      Compute solute sources ---
!
          CALL SORT_HYD_O2_CH4_N2( NSL )
!
!---      Zero solute transport fluxes  ---
!
          CALL SFXZ( NSL )
!
!---      Load Jacobian matrix (aqueous-phase transport)  ---
!
          CALL SJCBL( NSL )
!
!---      Load Jacobian matrix (gas-phase transport)  ---
!
          CALL SJCBG( NSL )
!
!---      Modify Jacobian matrix for boundary conditions ---
!
          CALL SBND_HYD_O2_CH4_N2( NSL )
!
!---      Linear equation solver  ---
!
          IF( ILES.EQ.1 ) THEN
            INDX = 1
            CALL BAND( 0,MUT,MLT,INDX )
          ELSEIF( ILES.EQ.3 ) THEN
            INDX = 1
            CALL PSPLIB( 0,INDX )
#ifdef lis
          ELSEIF( ILES.EQ.4 ) THEN
            INDX = 1
            CALL STOMP_LIS_SOLVE(-1,T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC,INDX)
#endif
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            CALL PETSC_PUT_RHS(BLU,T_RHS_VEC)
            CALL PETSC_SOLVE(T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC)
            CALL PETSC_GET_SOLUTION(BLU,T_SOL_VEC)
#endif
          ENDIF
!
!---      Update solute concentrations ---
!
          CALL UPDTC( NSL )
!
!---      Compute solute aqueous-phase fluxes (interior nodes)  ---
!
          CALL SFXL( NSL )
!
!---      Compute solute aqueous-phase fluxes (boundary surfaces)  ---
!
          CALL SFXLB( NSL )
!
!---      Compute solute gas-phase fluxes (interior nodes)  ---
!
          CALL SFXG( NSL )
!
!---      Compute solute gas-phase fluxes (boundary surfaces)  ---
!
          CALL SFXGB( NSL )
!
!---      Integrate solute sources  ---
!
          CALL SORIT_HYD_O2_CH4_N2( NSL )
!
!---  Load old sub-time-step concentrations  ---
!
          IF( ISLC(17).NE.0 ) CALL UPDTCO( NSL)
!
!---  Bottom of sub-time step loop  ---
!
  690   CONTINUE
!
!---  Courant number limiting, reset time stepping  ---
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
#ifdef radioactive
!
!---  Decay matrix, fracture, and borehole solutes via Bateman
!     chain decay solution  ---
!
      CALL CHAIN_DECAY
#endif
#ifdef reactive
!
!---  Solute reactive chemistry  ---
!
        CALL SRCHM
#endif
  800 CONTINUE
!
!---  Surface flux integrator  ---
!
      CALL SFIN
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
      IF( ISLC(18).LT.2 ) CALL WRRST
      WRITE(IWR,'(/,A)') '---  End of STOMP Simulation ---'
      WRITE(ISC,'(/,A)') '---  End of STOMP Simulation ---'
#ifdef lis
!
!---  Finalize Lis execution  ---
!
      CALL lis_finalize(IERR)
#endif
#ifdef petsc
!
!---  Finalize PETSc execution  ---
!
      CALL PETSC_SOLVER_DESTROY(F_KSP,F_MAT,F_RHS_VEC,F_SOL_VEC)
      IF( NSOLU.GT.0 ) THEN
        CALL PETSC_SOLVER_DESTROY(T_KSP,T_MAT,T_RHS_VEC,T_SOL_VEC)
      ENDIF
      CALL PetscFinalize(IERR)
#endif
      STOP
!
!---  End of STOMP program  ---
!
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCF_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Compute boundary surface fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
      USE FLUXD
      USE CONST
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
      REAL*8 BCX(LBCV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCF_HYD_O2_CH4_N2'
!
!---  Zero boundary fluxes  ---
!
      DO 70 NB = 1,NBC
        N = IBCN(NB)
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N)+1
        NQY = NSY(N)+IFLD
        NQZ = NSZ(N)+IJFLD
        IF( IBCD(NB).EQ.-3 ) THEN
          DO 10 M = 1,ISVF
            WL(M,NPZ) = 0.D+0
            WG(M,NPZ) = 0.D+0
            WLA(M,NPZ) = 0.D+0
            WLN(M,NPZ) = 0.D+0
            WLO(M,NPZ) = 0.D+0
            WLW(M,NPZ) = 0.D+0
            WGA(M,NPZ) = 0.D+0
            WGN(M,NPZ) = 0.D+0
            WGO(M,NPZ) = 0.D+0
            WGW(M,NPZ) = 0.D+0
            WDLA(M,NPZ) = 0.D+0
            WDLN(M,NPZ) = 0.D+0
            WDLO(M,NPZ) = 0.D+0
            WDLW(M,NPZ) = 0.D+0
            WDGA(M,NPZ) = 0.D+0
            WDGN(M,NPZ) = 0.D+0
            WDGO(M,NPZ) = 0.D+0
            WDGW(M,NPZ) = 0.D+0
            WQ(M,NPZ) = 0.D+0
            WS(M,NPZ) = 0.D+0
            WDS(M,NPZ) = 0.D+0
   10     CONTINUE
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DO 20 M = 1,ISVF
            VL(M,NPY) = 0.D+0
            VG(M,NPY) = 0.D+0
            VLA(M,NPY) = 0.D+0
            VLN(M,NPY) = 0.D+0
            VLO(M,NPY) = 0.D+0
            VLW(M,NPY) = 0.D+0
            VGA(M,NPY) = 0.D+0
            VGN(M,NPY) = 0.D+0
            VGO(M,NPY) = 0.D+0
            VGW(M,NPY) = 0.D+0
            VDLA(M,NPY) = 0.D+0
            VDLN(M,NPY) = 0.D+0
            VDLO(M,NPY) = 0.D+0
            VDLW(M,NPY) = 0.D+0
            VDGA(M,NPY) = 0.D+0
            VDGN(M,NPY) = 0.D+0
            VDGO(M,NPY) = 0.D+0
            VDGW(M,NPY) = 0.D+0
            VQ(M,NPY) = 0.D+0
            VS(M,NPY) = 0.D+0
            VDS(M,NPY) = 0.D+0
   20     CONTINUE
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DO 30 M = 1,ISVF
            UL(M,NPX) = 0.D+0
            UG(M,NPX) = 0.D+0
            ULA(M,NPX) = 0.D+0
            ULN(M,NPX) = 0.D+0
            ULO(M,NPX) = 0.D+0
            ULW(M,NPX) = 0.D+0
            UGA(M,NPX) = 0.D+0
            UGN(M,NPX) = 0.D+0
            UGO(M,NPX) = 0.D+0
            UGW(M,NPX) = 0.D+0
            UDLA(M,NPX) = 0.D+0
            UDLN(M,NPX) = 0.D+0
            UDLO(M,NPX) = 0.D+0
            UDLW(M,NPX) = 0.D+0
            UDGA(M,NPX) = 0.D+0
            UDGN(M,NPX) = 0.D+0
            UDGO(M,NPX) = 0.D+0
            UDGW(M,NPX) = 0.D+0
            UQ(M,NPX) = 0.D+0
            US(M,NPX) = 0.D+0
            UDS(M,NPX) = 0.D+0
   30     CONTINUE
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DO 40 M = 1,ISVF
            UL(M,NQX) = 0.D+0
            UG(M,NQX) = 0.D+0
            ULA(M,NQX) = 0.D+0
            ULN(M,NQX) = 0.D+0
            ULO(M,NQX) = 0.D+0
            ULW(M,NQX) = 0.D+0
            UGA(M,NQX) = 0.D+0
            UGN(M,NQX) = 0.D+0
            UGO(M,NQX) = 0.D+0
            UGW(M,NQX) = 0.D+0
            UDLA(M,NQX) = 0.D+0
            UDLN(M,NQX) = 0.D+0
            UDLO(M,NQX) = 0.D+0
            UDLW(M,NQX) = 0.D+0
            UDGA(M,NQX) = 0.D+0
            UDGN(M,NQX) = 0.D+0
            UDGO(M,NQX) = 0.D+0
            UDGW(M,NQX) = 0.D+0
            UQ(M,NQX) = 0.D+0
            US(M,NQX) = 0.D+0
            UDS(M,NQX) = 0.D+0
   40     CONTINUE
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DO 50 M = 1,ISVF
            VL(M,NQY) = 0.D+0
            VG(M,NQY) = 0.D+0
            VLA(M,NQY) = 0.D+0
            VLN(M,NQY) = 0.D+0
            VLO(M,NQY) = 0.D+0
            VLW(M,NQY) = 0.D+0
            VGA(M,NQY) = 0.D+0
            VGN(M,NQY) = 0.D+0
            VGO(M,NQY) = 0.D+0
            VGW(M,NQY) = 0.D+0
            VDLA(M,NQY) = 0.D+0
            VDLN(M,NQY) = 0.D+0
            VDLO(M,NQY) = 0.D+0
            VDLW(M,NQY) = 0.D+0
            VDGA(M,NQY) = 0.D+0
            VDGN(M,NQY) = 0.D+0
            VDGO(M,NQY) = 0.D+0
            VDGW(M,NQY) = 0.D+0
            VQ(M,NQY) = 0.D+0
            VS(M,NQY) = 0.D+0
            VDS(M,NQY) = 0.D+0
   50     CONTINUE
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DO 60 M = 1,ISVF
            WL(M,NQZ) = 0.D+0
            WG(M,NQZ) = 0.D+0
            WLA(M,NQZ) = 0.D+0
            WLN(M,NQZ) = 0.D+0
            WLO(M,NQZ) = 0.D+0
            WLW(M,NQZ) = 0.D+0
            WGA(M,NQZ) = 0.D+0
            WGN(M,NQZ) = 0.D+0
            WGO(M,NQZ) = 0.D+0
            WGW(M,NQZ) = 0.D+0
            WDLA(M,NQZ) = 0.D+0
            WDLN(M,NQZ) = 0.D+0
            WDLO(M,NQZ) = 0.D+0
            WDLW(M,NQZ) = 0.D+0
            WDGA(M,NQZ) = 0.D+0
            WDGN(M,NQZ) = 0.D+0
            WDGO(M,NQZ) = 0.D+0
            WDGW(M,NQZ) = 0.D+0
            WQ(M,NQZ) = 0.D+0
            WS(M,NQZ) = 0.D+0
            WDS(M,NQZ) = 0.D+0
   60     CONTINUE
        ENDIF
   70 CONTINUE
!
!---  Loop over boundary conditions  ---
!
      DO 200 NB = 1,NBC
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 200
        IF( IBCM(NB).EQ.1 ) THEN
          DO 80 N = 1,LBCV
            BCX(N) = BC(N,1,MB)
   80     CONTINUE
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
             TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
             DTBC = MIN( BC(1,M,MB)-TMZ,DT )
             TFBC = (TMZ-5.D-1*DTBC-BC(1,M-1,MB))/TDBC
             DO 90 N = 1,LBCV
               BCX(N) = BC(N,M-1,MB) + TFBC*(BC(N,M,MB)-BC(N,M-1,MB))
   90        CONTINUE
             GOTO 105
            ENDIF
  100     CONTINUE
          GOTO 200
        ENDIF
  105   CONTINUE
        N = IBCN(NB)
        I = ID(N)
        J = JD(N)
        K = KD(N)
        NPZ = NSZ(N)
        NPY = NSY(N)
        NPX = NSX(N)
        NQX = NPX+1
        NQY = NPY+IFLD
        NQZ = NPZ+IJFLD
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---      Neumann flow boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 110 M = 1,ISVF
              WG(M,NPZ) = BCX(3)*SGB(2,NB)
              WL(M,NPZ) = BCX(3)*SLB(2,NB)
  110       CONTINUE
            IF( ISLC(45).EQ.0 ) CALL DFFLAB( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOB( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNB( N,NB )
            CALL DFFLWB( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAB( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOB( N,NB )
            CALL DFFGWB( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNB( N,NB )
!
!---      Aqueous-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.21 .OR. IBCT(2,NB).EQ.27 ) THEN
            CALL DRCVLB( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAB( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOB( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNB( N,NB )
            CALL DFFLWB( N,NB )
!
!---      Gas-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.31 .OR. IBCT(2,NB).EQ.37 ) THEN
            CALL DRCVGB( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAB( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOB( N,NB )
            CALL DFFGWB( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNB( N,NB )
!
!---      Non zero-flux flow boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLB( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAB( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOB( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNB( N,NB )
            CALL DFFLWB( N,NB )
            CALL DRCVGB( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAB( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOB( N,NB )
            CALL DFFGWB( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNB( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(3,NB).NE.3 ) THEN
              CALL DFFLSB( N,NB )
            ENDIF
          ENDIF
!
!---      Energy Neumann  ---
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 114 M = 1,ISVF
              WQ(M,NPZ) = BCX(2)
  114       CONTINUE
!
!---      Energy outflow, inflow, or inflow-outflow  ---
!
          ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13
     &      .OR. IBCT(1,NB).EQ.14  ) THEN
            CALL THLB_HYDT_KE( N,NB )
            CALL THGB_HYDT_KE( N,NB )
!
!---      Energy zero flux  ---
!
          ELSEIF( IBCT(1,NB).NE.3 ) THEN
            CALL THDB_HYDT_KE( N,NB )
            CALL THLB_HYDT_KE( N,NB )
            CALL THGB_HYDT_KE( N,NB )
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Neumann flow boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 120 M = 1,ISVF
              VG(M,NPY) = BCX(3)*SGB(2,NB)
              VL(M,NPY) = BCX(3)*SLB(2,NB)
  120       CONTINUE
            IF( ISLC(45).EQ.0 ) CALL DFFLAS( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOS( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNS( N,NB )
            CALL DFFLWS( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAS( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOS( N,NB )
            CALL DFFGWS( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNS( N,NB )
!
!---      Aqueous-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.21 .OR. IBCT(2,NB).EQ.27 ) THEN
            CALL DRCVLS( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAS( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOS( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNS( N,NB )
            CALL DFFLWS( N,NB )
!
!---      Gas-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.31 .OR. IBCT(2,NB).EQ.37 ) THEN
            CALL DRCVGS( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAS( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOS( N,NB )
            CALL DFFGWS( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNS( N,NB )
!
!---      Non zero-flux flow boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLS( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAS( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOS( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNS( N,NB )
            CALL DFFLWS( N,NB )
            CALL DRCVGS( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAS( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOS( N,NB )
            CALL DFFGWS( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNS( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(3,NB).NE.3 ) THEN
              CALL DFFLSS( N,NB )
            ENDIF
          ENDIF
!
!---      Energy Neumann  ---
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 124 M = 1,ISVF
              VQ(M,NPY) = BCX(2)
  124       CONTINUE
!
!---      Energy outflow, inflow, or inflow-outflow  ---
!
          ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13
     &      .OR. IBCT(1,NB).EQ.14  ) THEN
            CALL THLS_HYDT_KE( N,NB )
            CALL THGS_HYDT_KE( N,NB )
!
!---      Energy not zero flux  ---
!
          ELSEIF( IBCT(1,NB).NE.3 ) THEN
            CALL THDS_HYDT_KE( N,NB )
            CALL THLS_HYDT_KE( N,NB )
            CALL THGS_HYDT_KE( N,NB )
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---     Neumann flow boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 130 M = 1,ISVF
              UG(M,NPX) = BCX(3)*SGB(2,NB)
              UL(M,NPX) = BCX(3)*SLB(2,NB)
  130       CONTINUE
            IF( ISLC(45).EQ.0 ) CALL DFFLAW( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOW( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNW( N,NB )
            CALL DFFLWW( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAW( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOW( N,NB )
            CALL DFFGWW( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNW( N,NB )
!
!---      Aqueous-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.21 .OR. IBCT(2,NB).EQ.27 ) THEN
            CALL DRCVLW( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAW( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOW( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNW( N,NB )
            CALL DFFLWW( N,NB )
!
!---      Gas-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.31 .OR. IBCT(2,NB).EQ.37 ) THEN
            CALL DRCVGW( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAW( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOW( N,NB )
            CALL DFFGWW( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNW( N,NB )
!
!---      Non zero-flux flow boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLW( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAW( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOW( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNW( N,NB )
            CALL DFFLWW( N,NB )
            CALL DRCVGW( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAW( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOW( N,NB )
            CALL DFFGWW( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNW( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(3,NB).NE.3 ) THEN
              CALL DFFLSW( N,NB )
            ENDIF
          ENDIF
!
!---      Energy Neumann  ---
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 134 M = 1,ISVF
              UQ(M,NPX) = BCX(2)
  134       CONTINUE
!
!---      Energy outflow, inflow, or inflow-outflow  ---
!
          ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13
     &      .OR. IBCT(1,NB).EQ.14  ) THEN
            CALL THLW_HYDT_KE( N,NB )
            CALL THGW_HYDT_KE( N,NB )
!
!---      Energy not zero flux  ---
!
          ELSEIF( IBCT(1,NB).NE.3 ) THEN
            CALL THDW_HYDT_KE( N,NB )
            CALL THLW_HYDT_KE( N,NB )
            CALL THGW_HYDT_KE( N,NB )
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Neumann flow boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 140 M = 1,ISVF
              UG(M,NQX) = BCX(3)*SGB(2,NB)
              UL(M,NQX) = BCX(3)*SLB(2,NB)
  140       CONTINUE
            IF( ISLC(45).EQ.0 ) CALL DFFLAE( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOE( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNE( N,NB )
            CALL DFFLWE( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAE( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOE( N,NB )
            CALL DFFGWE( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNE( N,NB )
!
!---      Aqueous-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.21 .OR. IBCT(2,NB).EQ.27 ) THEN
            CALL DRCVLE( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAE( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOE( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNE( N,NB )
            CALL DFFLWE( N,NB )
!
!---      Gas-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.31 .OR. IBCT(2,NB).EQ.37 ) THEN
            CALL DRCVGE( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAE( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOE( N,NB )
            CALL DFFGWE( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNE( N,NB )
!
!---      Non zero-flux flow boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLE( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAE( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOE( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNE( N,NB )
            CALL DFFLWE( N,NB )
            CALL DRCVGE( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAE( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOE( N,NB )
            CALL DFFGWE( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNE( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(3,NB).NE.3 ) THEN
              CALL DFFLSE( N,NB )
            ENDIF
          ENDIF
!
!---      Energy Neumann  ---
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 144 M = 1,ISVF
              UQ(M,NQX) = BCX(2)
  144       CONTINUE
!
!---      Energy outflow, inflow, or inflow-outflow  ---
!
          ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13
     &       .OR. IBCT(1,NB).EQ.14  ) THEN
            CALL THLE_HYDT_KE( N,NB )
            CALL THGE_HYDT_KE( N,NB )
!
!---      Energy not zero flux  ---
!
          ELSEIF( IBCT(1,NB).NE.3 ) THEN
            CALL THDE_HYDT_KE( N,NB )
            CALL THLE_HYDT_KE( N,NB )
            CALL THGE_HYDT_KE( N,NB )
          ENDIF
!
!---    North boundary
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Neumann flow boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 150 M = 1,ISVF
              VG(M,NQY) = BCX(3)*SGB(2,NB)
              VL(M,NQY) = BCX(3)*SLB(2,NB)
  150       CONTINUE
            IF( ISLC(45).EQ.0 ) CALL DFFLAN( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLON( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNN( N,NB )
            CALL DFFLWN( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAN( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGON( N,NB )
            CALL DFFGWN( N,NB )
!
!---      Aqueous-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.21 .OR. IBCT(2,NB).EQ.27 ) THEN
            CALL DRCVLN( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAN( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLON( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNN( N,NB )
            CALL DFFLWN( N,NB )
!
!---      Gas-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.31 .OR. IBCT(2,NB).EQ.37 ) THEN
            CALL DRCVGN( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAN( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGON( N,NB )
            CALL DFFGWN( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNN( N,NB )
!
!---      Non zero-flux flow boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLN( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAN( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLON( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNN( N,NB )
            CALL DFFLWN( N,NB )
            CALL DRCVGN( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAN( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGON( N,NB )
            CALL DFFGWN( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNN( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(3,NB).NE.3 ) THEN
              CALL DFFLSN( N,NB )
            ENDIF
          ENDIF
!
!---      Energy Neumann  ---
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 154 M = 1,ISVF
              VQ(M,NQY) = BCX(2)
  154       CONTINUE
!
!---      Energy outflow, inflow, or inflow-outflow  ---
!
          ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13
     &      .OR. IBCT(1,NB).EQ.14  ) THEN
            CALL THLN_HYDT_KE( N,NB )
            CALL THGN_HYDT_KE( N,NB )
!
!---      Energy not zero flux  ---
!
          ELSEIF( IBCT(1,NB).NE.3 ) THEN
            CALL THDN_HYDT_KE( N,NB )
            CALL THLN_HYDT_KE( N,NB )
            CALL THGN_HYDT_KE( N,NB )
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Neumann flow boundary  ---
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 160 M = 1,ISVF
              WG(M,NQZ) = BCX(3)*SGB(2,NB)
              WL(M,NQZ) = BCX(3)*SLB(2,NB)
  160       CONTINUE
            IF( ISLC(45).EQ.0 ) CALL DFFLAT( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOT( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNT( N,NB )
            CALL DFFLWT( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAT( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOT( N,NB )
            CALL DFFGWT( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNT( N,NB )
!
!---      Aqueous-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.21 .OR. IBCT(2,NB).EQ.27 ) THEN
            CALL DRCVLT( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAT( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOT( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNT( N,NB )
            CALL DFFLWT( N,NB )
!
!---      Gas-only Dirichlet boundary  ---
!
          ELSEIF( IBCT(2,NB).EQ.31 .OR. IBCT(2,NB).EQ.37 ) THEN
            CALL DRCVGT( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAT( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOT( N,NB )
            CALL DFFGWT( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNT( N,NB )
!
!---      Non zero-flux flow boundary  ---
!
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVLT( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFLAT( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFLOT( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFLNT( N,NB )
            CALL DFFLWT( N,NB )
            CALL DRCVGT( N,NB )
            IF( ISLC(45).EQ.0 ) CALL DFFGAT( N,NB )
            IF( ISLC(64).EQ.0 ) CALL DFFGOT( N,NB )
            CALL DFFGWT( N,NB )
            IF( ISLC(65).EQ.0 ) CALL DFFGNT( N,NB )
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Salt not zero flux boundary  ---
!
            IF( IBCT(3,NB).NE.3 ) THEN
              CALL DFFLST( N,NB )
            ENDIF
          ENDIF
!
!---      Energy Neumann  ---
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 164 M = 1,ISVF
              WQ(M,NQZ) = BCX(2)
  164       CONTINUE
!
!---      Energy outflow, inflow, or inflow-outflow  ---
!
          ELSEIF( IBCT(1,NB).EQ.7 .OR. IBCT(1,NB).EQ.13
     &      .OR. IBCT(1,NB).EQ.14  ) THEN
            CALL THLT_HYDT_KE( N,NB )
            CALL THGT_HYDT_KE( N,NB )
!
!---      Energy not zero flux  ---
!
          ELSEIF( IBCT(1,NB).NE.3 ) THEN
            CALL THDT_HYDT_KE( N,NB )
            CALL THLT_HYDT_KE( N,NB )
            CALL THGT_HYDT_KE( N,NB )
          ENDIF
        ENDIF
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCF_HYD_O2_CH4_N2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCJ_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Modify the Jacobian matrix for boundary conditions
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
      SUB_LOG(ISUB_LOG) = '/BCJ_HYD_O2_CH4_N2'
!
!---  Loop over boundary conditions  ---
!
      DO 100 NB = 1,NBC
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 100
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) GOTO 100
        N = IBCN(NB)
        IF( IXP(N).EQ.0 ) GOTO 100
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX + 1
        NQY = NPY + IFLD
        NQZ = NPZ + IJFLD
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---      Flow  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLWB( N,NB,NPZ )
            IF( ISLC(45).EQ.0 ) CALL JCBLAB( N,NB,NPZ )
            IF( ISLC(64).EQ.0 ) CALL JCBLOB( N,NB,NPZ )
            IF( ISLC(65).EQ.0 ) CALL JCBLNB( N,NB,NPZ )
            CALL JCBGWB( N,NB,NPZ )
            IF( ISLC(45).EQ.0 ) CALL JCBGAB( N,NB,NPZ )
            IF( ISLC(64).EQ.0 ) CALL JCBGOB( N,NB,NPZ )
            IF( ISLC(65).EQ.0 ) CALL JCBGNB( N,NB,NPZ )
          ENDIF
!
!---      Energy  ---
!
          IF( IBCT(1,NB).NE.3 ) THEN
            CALL JCBTB( N,NB,NPZ )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(3,NB).NE.3 ) CALL JCBSB( N,NB,NPZ )
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Flow  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLWS( N,NB,NPY )
            IF( ISLC(45).EQ.0 ) CALL JCBLAS( N,NB,NPY )
            IF( ISLC(64).EQ.0 ) CALL JCBLOS( N,NB,NPY )
            IF( ISLC(65).EQ.0 ) CALL JCBLNS( N,NB,NPY )
            CALL JCBGWS( N,NB,NPY )
            IF( ISLC(45).EQ.0 ) CALL JCBGAS( N,NB,NPY )
            IF( ISLC(64).EQ.0 ) CALL JCBGOS( N,NB,NPY )
            IF( ISLC(65).EQ.0 ) CALL JCBGNS( N,NB,NPY )
          ENDIF
!
!---      Energy  ---
!
          IF( IBCT(1,NB).NE.3 ) THEN
            CALL JCBTS( N,NB,NPY )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(3,NB).NE.3 ) CALL JCBSS( N,NB,NPY )
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Flow  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLWW( N,NB,NPX )
            IF( ISLC(45).EQ.0 ) CALL JCBLAW( N,NB,NPX )
            IF( ISLC(64).EQ.0 ) CALL JCBLOW( N,NB,NPX )
            IF( ISLC(65).EQ.0 ) CALL JCBLNW( N,NB,NPX )
            CALL JCBGWW( N,NB,NPX )
            IF( ISLC(45).EQ.0 ) CALL JCBGAW( N,NB,NPX )
            IF( ISLC(64).EQ.0 ) CALL JCBGOW( N,NB,NPX )
            IF( ISLC(65).EQ.0 ) CALL JCBGNW( N,NB,NPX )
          ENDIF
!
!---      Energy  ---
!
          IF( IBCT(1,NB).NE.3 ) THEN
             CALL JCBTW( N,NB,NPX )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(3,NB).NE.3 ) CALL JCBSW( N,NB,NPX )
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Flow  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLWE( N,NB,NQX )
            IF( ISLC(45).EQ.0 ) CALL JCBLAE( N,NB,NQX )
            IF( ISLC(64).EQ.0 ) CALL JCBLOE( N,NB,NQX )
            IF( ISLC(65).EQ.0 ) CALL JCBLNE( N,NB,NQX )
            CALL JCBGWE( N,NB,NQX )
            IF( ISLC(45).EQ.0 ) CALL JCBGAE( N,NB,NQX )
            IF( ISLC(64).EQ.0 ) CALL JCBGOE( N,NB,NQX )
            IF( ISLC(65).EQ.0 ) CALL JCBGNE( N,NB,NQX )
          ENDIF
!
!---      Energy  ---
!
          IF( IBCT(1,NB).NE.3 ) THEN
            CALL JCBTE( N,NB,NQX )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(3,NB).NE.3 ) CALL JCBSE( N,NB,NQX )
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Flow  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLWN( N,NB,NQY )
            IF( ISLC(45).EQ.0 ) CALL JCBLAN( N,NB,NQY )
            IF( ISLC(64).EQ.0 ) CALL JCBLON( N,NB,NQY )
            IF( ISLC(65).EQ.0 ) CALL JCBLNN( N,NB,NQY )
            CALL JCBGWN( N,NB,NQY )
            IF( ISLC(45).EQ.0 ) CALL JCBGAN( N,NB,NQY )
            IF( ISLC(64).EQ.0 ) CALL JCBGON( N,NB,NQY )
            IF( ISLC(65).EQ.0 ) CALL JCBGNN( N,NB,NQY )
          ENDIF
!
!---      Energy  ---
!
          IF( IBCT(1,NB).NE.3 ) THEN
            CALL JCBTN( N,NB,NQY )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(3,NB).NE.3 ) CALL JCBSN( N,NB,NQY )
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Flow  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBLWT( N,NB,NQZ )
            IF( ISLC(45).EQ.0 ) CALL JCBLAT( N,NB,NQZ )
            IF( ISLC(64).EQ.0 ) CALL JCBLOT( N,NB,NQZ )
            IF( ISLC(65).EQ.0 ) CALL JCBLNT( N,NB,NQZ )
            CALL JCBGWT( N,NB,NQZ )
            IF( ISLC(45).EQ.0 ) CALL JCBGAT( N,NB,NQZ )
            IF( ISLC(64).EQ.0 ) CALL JCBGOT( N,NB,NQZ )
            IF( ISLC(65).EQ.0 ) CALL JCBGNT( N,NB,NQZ )
          ENDIF
!
!---      Energy  ---
!
          IF( IBCT(1,NB).NE.3 ) THEN
            CALL JCBTT( N,NB,NQZ )
          ENDIF
!
!---      Salt  ---
!
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(3,NB).NE.3 ) CALL JCBST( N,NB,NQZ )
          ENDIF
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJ_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCP_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Compute saturation, relative permeability and thermodynamic
!     properties for boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
      USE JACOB
      USE GRID
      USE FDVS
      USE FDVP
      USE FDVH
      USE CONST
      USE BCVT
      USE BCVS
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVG
      USE BCVA
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
      REAL*8 PHIX(LHF_HT+1,LHF_HT+1)
      REAL*8 ZMCX(LHF_HT)
      REAL*8 HIDX(LHF_HT),HRESX(LHF_HT),XMGX(LHF_HT),XMNX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCP_HYD_O2_CH4_N2'
      EPSLX = 1.D-4
!
!---  Assign values for initial condition type boundary conditions  ---
!
      IF( NSTEP-NRST.LE.1 .AND. NITER.LE.1 ) THEN
        DO 50 NB = 1,NBC
          N = IBCN(NB)
          IF( IBCT(2,NB).EQ.12 .OR. IBCT(2,NB).EQ.14 .OR.
     &      IBCT(2,NB).EQ.15 ) THEN
            IF( IBCD(NB).EQ.-3 ) THEN
              DB = 0.5D+0*DZGF(N)
              NPZ = NSZ(N)
              DB = DZGP(NPZ)
              GB = GRVZ(NPZ)*DB
            ELSEIF( IBCD(NB).EQ.-2 ) THEN
              DB = 0.5D+0*DYGF(N)*RP(ID(N))
              NPY = NSY(N)
              DB = DYGP(NPY)*RP(ID(N))
              GB = GRVY(NPY)*DB
            ELSEIF( IBCD(NB).EQ.-1 ) THEN
              DB = 0.5D+0*DXGF(N)
              NPX = NSX(N)
              DB = DXGP(NPX)
              GB = GRVX(NPX)*DB
            ELSEIF( IBCD(NB).EQ.1 ) THEN
              DB = -0.5D+0*DXGF(N)
              NQX = NSX(N)+1
              IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
              DB = -DXGP(NQX)
              GB = GRVX(NQX)*DB
            ELSEIF( IBCD(NB).EQ.2 ) THEN
              DB = -0.5D+0*DYGF(N)*RP(ID(N))
              NQY = NSY(N)+IFLD
              IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
              DB = -DYGP(NQY)*RP(ID(N))
              GB = GRVY(NQY)*DB
            ELSEIF( IBCD(NB).EQ.3 ) THEN
              DB = -0.5D+0*DZGF(N)
              NQZ = NSZ(N)+IJFLD
              IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
              DB = -DZGP(NQZ)
              GB = GRVZ(NQZ)*DB
            ENDIF
          ENDIF
          IF( IBCT(1,NB).EQ.12 ) TB(1,NB) = T(2,N)
          IF( IBCT(2,NB).EQ.12 .OR. IBCT(2,NB).EQ.14 .OR.
     &      IBCT(2,NB).EQ.15 ) THEN
            IF( SG(2,N).GT.EPSL ) THEN
              PGB(1,NB) = PG(2,N) + RHOG(2,N)*GB
              IF( SN(2,N).GT.EPSL ) THEN
                PNB(1,NB) = PN(2,N) + RHON(2,N)*GB
                IF( SL(2,N).GT.EPSL ) THEN
                  PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
                ELSE
                  PLB(1,NB) = PNB(1,NB) + PL(2,N) - PN(2,N)
                ENDIF
              ELSE
                IF( SL(2,N).GT.EPSL ) THEN
                  PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
                ELSE
                  PLB(1,NB) = PGB(1,NB) + PL(2,N) - PG(2,N)
                  PNB(1,NB) = PGB(1,NB) + PN(2,N) - PG(2,N)
                ENDIF
              ENDIF
            ELSE
              IF( SN(2,N).GT.EPSL ) THEN
                PNB(1,NB) = PN(2,N) + RHON(2,N)*GB
                PGB(1,NB) = PNB(1,NB) + PG(2,N) - PN(2,N)
                IF( SL(2,N).GT.EPSL ) THEN
                  PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
                ELSE
                  PLB(1,NB) = PNB(1,NB) + PL(2,N) - PN(2,N)
                ENDIF
              ELSE
                PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
                PNB(1,NB) = PLB(1,NB) + PN(2,N) - PL(2,N)
                PGB(1,NB) = PLB(1,NB) + PG(2,N) - PL(2,N)
              ENDIF
            ENDIF
            PSOB(1,NB) = MAX( PGB(1,NB),PLB(1,NB),PNB(1,NB) )
            SLB(1,NB) = SL(2,N)
            SHB(1,NB) = SH(2,N)
            SGB(1,NB) = SG(2,N)
            SNB(1,NB) = SN(2,N)
            XMLAB(1,NB) = ZMCA(2,N)
            XMLOB(1,NB) = ZMCO(2,N)
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 )  THEN
            IF( IBCT(3,NB).EQ.12 ) YLSB(1,NB) = YLS(2,N)
          ENDIF
   50   CONTINUE
      ENDIF
!
!---  Loop over boundary conditions  ---
!
      DO 400 NB = 1,NBC
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 400
!
!---  Assign local boundary condition variables  ---
!
        IF( IBCM(NB).EQ.1 ) THEN
          DO 80 N = 1,LBCV
            BCX(N) = BC(N,1,MB)
   80     CONTINUE
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
             TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
             DTBC = MIN( BC(1,M,MB)-TMZ,DT )
             TFBC = (TMZ-BC(1,M-1,MB))/TDBC
             DO 90 N = 1,LBCV
               BCX(N) = BC(N,M-1,MB) + TFBC*(BC(N,M,MB)-BC(N,M-1,MB))
   90        CONTINUE
!
!---         Energy boundary  ---
!
             IF( IBCT(1,NB).EQ.2 ) THEN
               BCX(2) = BCX(2)-5.D-1*DTBC*(BC(2,M,MB)-BC(2,M-1,MB))/TDBC
             ENDIF
!
!---         Flow boundary  ---
!
             IF( IBCT(2,NB).EQ.2 ) THEN
               BCX(3) = BCX(3)-5.D-1*DTBC*(BC(3,M,MB)-BC(3,M-1,MB))/TDBC
             ENDIF
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 400
        ENDIF
  110   CONTINUE
!
!---    Initial condition boundary condition  ---
!
        IF( IBCT(1,NB).EQ.12 ) BCX(2) = TB(1,NB)
        IF( IBCT(2,NB).EQ.12 .OR. IBCT(2,NB).EQ.14 .OR.
     &      IBCT(2,NB).EQ.15 ) THEN
          BCX(3) = PSOB(1,NB)
          BCX(4) = SLB(1,NB)
          BCX(5) = SHB(1,NB)
          BCX(6) = SIB(1,NB)
          BCX(10) = XMLAB(1,NB)
          BCX(11) = XMLOB(1,NB)
          BCX(12) = 1.D+0
          BCX(8) = 1.D+0
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 )  THEN
          IF( IBCT(3,NB).EQ.12 ) BCX(6) = YLSB(1,NB)
        ENDIF
        N = IBCN(NB)
        N_DB = NB
        IBD = ABS(IBCD(NB))
        IZN = IZ(N)
!
!---    Boundary Direction  ---
!
        I = ID(N)
        J = JD(N)
        K = KD(N)
        I1 = I+1
        J1 = J+1
        K1 = K+1
        IF( IBCD(NB).EQ.-3 ) THEN
          DB = 0.5D+0*DZGF(N)
          GB = GRVZ(NSZ(N))*DB
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DB = 0.5D+0*DYGF(N)*RP(ID(N))
          GB = GRVY(NSY(N))*DB
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DB = 0.5D+0*DXGF(N)
          GB = GRVX(NSX(N))*DB
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DB = -0.5D+0*DXGF(N)
          GB = GRVX(NSX(N)+1)*DB
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DB = -0.5D+0*DYGF(N)*RP(ID(N))
          GB = GRVY(NSY(N)+IFLD)*DB
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DB = -0.5D+0*DZGF(N)
          GB = GRVZ(NSZ(N)+IJFLD)*DB
        ENDIF
!
!---    Loop over secondary variable indices  ---
!
        DO 300 M = 2,ISVC+2
!
!---      Energy: Dirichlet, initial conditions, inflow, outflow,
!         or inflow-outflow  ---
!
          IF( IBCT(1,NB).EQ.1 .OR. IBCT(1,NB).EQ.12 .OR. 
     &      IBCT(1,NB).EQ.13 .OR. IBCT(1,NB).EQ.14 ) THEN
            TB(M,NB) = BCX(2)
!
!---      Energy: Neumann  ---
!
          ELSEIF( IBCT(1,NB).EQ.2 ) THEN
            TB(M,NB) = T(M,N)
!
!---      Energy: zero flux  ---
!
          ELSEIF( IBCT(1,NB).EQ.3 ) THEN
            TB(M,NB) = T(M,N)
!
!---      Energy: outflow  ---
!
          ELSEIF( IBCT(1,NB).EQ.7 ) THEN
            TB(M,NB) = T(M,N)
          ENDIF
!
!---      Flow: Dirichlet or Dirichlet inflow  ---
!
          IF( MOD(IBCT(2,NB),10).EQ.1.OR.MOD(IBCT(2,NB),10).EQ.7 ) THEN
            PSOB(M,NB) = BCX(3)
            SLB(M,NB) = BCX(4)
            SHB(M,NB) = 0.D+0
            SIB(M,NB) = 0.D+0
            ZMCX(1) = BCX(10)
            ZMCX(2) = BCX(11)
            ZMCX(3) = MAX( 1.D+0-ZMCX(1)-ZMCX(2),0.D+0 )
            WSGWX = BCX(12)
            WSLCX = BCX(8)
!
!---      Flow: initial conditions or initial conditions inflow  ---
!
          ELSEIF( IBCT(2,NB).EQ.12 .OR. IBCT(2,NB).EQ.14 ) THEN
            PSOB(M,NB) = BCX(3)
            SLB(M,NB) = BCX(4)
            SHB(M,NB) = BCX(5)
            SIB(M,NB) = BCX(6)
            ZMCX(1) = BCX(10)
            ZMCX(2) = BCX(11)
            ZMCX(3) = MAX( 1.D+0-ZMCX(1)-ZMCX(2),0.D+0 )
            WSGWX = BCX(12)
            WSLCX = BCX(8)
!
!---      Flow: hydrostatic, hydrostatic inflow  ---
!
          ELSEIF( IBCT(2,NB).EQ.28 .OR. IBCT(2,NB).EQ.29 ) THEN
            PSOB(M,NB) = BCX(3) + BCX(5)*(ZPBC(NB)-BCX(6))
            SLB(M,NB) = BCX(4)
            SHB(M,NB) = 0.D+0
            SIB(M,NB) = 0.D+0
            ZMCX(1) = BCX(10)
            ZMCX(2) = BCX(11)
            ZMCX(3) = MAX( 1.D+0-ZMCX(1)-ZMCX(2),0.D+0 )
            WSGWX = BCX(12)
            WSLCX = BCX(8)
!
!---      Flow: initial condition outflow, Dirichlet outflow  ---
!
          ELSEIF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 ) THEN
            PSOB(M,NB) = BCX(3)
            SLB(M,NB) = SL(M,N)
            SGB(M,NB) = SG(M,N)
            SNB(M,NB) = SN(M,N)
            SHB(M,NB) = SH(M,N)
            SIB(M,NB) = SI(M,N)
            ZMCX(1) = ZMCA(M,N)
            ZMCX(2) = ZMCO(M,N)
            ZMCX(3) = ZMCN(M,N)
            WSGWX = 1.D+0
            WSLCX = 1.D+0
!
!---      Flow: hydrostatic outflow  ---
!
          ELSEIF( IBCT(2,NB).EQ.30 ) THEN
            PSOB(M,NB) = BCX(3) + BCX(5)*(ZPBC(NB)-BCX(6))
            SLB(M,NB) = SL(M,N)
            SGB(M,NB) = SG(M,N)
            SNB(M,NB) = SN(M,N)
            SHB(M,NB) = SH(M,N)
            SIB(M,NB) = SI(M,N)
            ZMCX(1) = ZMCA(M,N)
            ZMCX(2) = ZMCO(M,N)
            ZMCX(3) = ZMCN(M,N)
            WSGWX = 1.D+0
            WSLCX = 1.D+0
!
!---      Flow: Neumann  ---
!
          ELSEIF( IBCT(2,NB).EQ.2 ) THEN
            PSOB(M,NB) = PSO(M,N)
            SLB(M,NB) = BCX(4)
            SHB(M,NB) = 0.D+0
            SIB(M,NB) = 0.D+0
            ZMCX(1) = BCX(10)
            ZMCX(2) = BCX(11)
            ZMCX(3) = MAX( 1.D+0-ZMCX(1)-ZMCX(2),0.D+0 )
            WSGWX = BCX(12)
            WSLCX = BCX(8)
!
!---      Flow: zero flux ---
!
          ELSEIF( IBCT(2,NB).EQ.3 ) THEN
            PGBX = PG(M,N)
            PLBX = PL(M,N)
            PNBX = PN(M,N)
            IF( SG(M,N).GT.EPSL ) THEN
              PGBX = PG(M,N) + RHOG(M,N)*GB
              IF( SN(M,N).GT.EPSL ) THEN
                PNBX = PN(M,N) + RHON(M,N)*GB
                IF( SL(M,N).GT.EPSL ) THEN
                  PLBX = PL(M,N) + RHOL(M,N)*GB
                ELSE
                  PLBX = PNBX + PL(M,N) - PN(M,N)
                ENDIF
              ELSE
                IF( SL(M,N).GT.EPSL ) THEN
                  PLBX = PL(M,N) + RHOL(M,N)*GB
                ELSE
                  PLBX = PGBX + PL(M,N) - PG(M,N)
                  PNBX = PGBX + PN(M,N) - PG(M,N)
                ENDIF
              ENDIF
            ELSE
              IF( SN(M,N).GT.EPSL ) THEN
                PNBX = PN(M,N) + RHON(M,N)*GB
                PGBX = PNBX + PG(M,N) - PN(M,N)
                IF( SL(M,N).GT.EPSL ) THEN
                  PLBX = PL(M,N) + RHOL(M,N)*GB
                ELSE
                  PLBX = PNBX + PL(M,N) - PN(M,N)
                ENDIF
              ELSE
                PLBX = PL(M,N) + RHOL(M,N)*GB
                PNBX = PLBX + PN(M,N) - PL(M,N)
                PGBX = PLBX + PG(M,N) - PL(M,N)
              ENDIF
            ENDIF
            PSOB(M,NB) = MAX( PGBX,PLBX,PNBX )
            SLB(M,NB) = SL(M,N)
            SGB(M,NB) = SG(M,N)
            SNB(M,NB) = SN(M,N)
            SHB(M,NB) = SH(M,N)
            SIB(M,NB) = SI(M,N)
            ZMCX(1) = ZMCA(M,N)
            ZMCX(2) = ZMCO(M,N)
            ZMCX(3) = MAX( 1.D+0-ZMCX(1)-ZMCX(2),0.D+0 )
            WSGWX = 1.D+0
            WSLCX = 1.D+0
          ENDIF
!
!---      Salt concentration  ---
!
          CALL SOL_LS( TB(M,NB),XLSMX )
          IF( ISLC(32).EQ.0 ) THEN
            IF( IBCT(3,NB).EQ.34 .OR. IBCT(3,NB).EQ.35 ) THEN
              YLSB(M,NB) = BCX(6)*XLSMX
              XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
            ELSEIF( IBCT(3,NB).EQ.36 .OR. IBCT(3,NB).EQ.37 ) THEN
              YLSB(M,NB) = BCX(6)
              XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
            ENDIF
          ELSE
            YLSB(M,NB) = 0.D+0
            XLSB(M,NB) = MIN( YLSB(M,NB),XLSMX )
          ENDIF
!
!---      Phase equilibria for boundary surface  ---
!
          CALL FLH_BC( PGB(M,NB),PLB(M,NB),PNB(M,NB),PPELX,PPEUX,PSBX,
     &      PSOB(M,NB),PVAB(M,NB),PVNB(M,NB),PVOB(M,NB),PVWB(M,NB),
     &      RHOGB(M,NB),RHONB(M,NB),RKGB(M,NB),RKLB(1,M,NB),
     &      RKNB(M,NB),SGB(M,NB),SHB(M,NB),SIB(M,NB),SLB(M,NB),
     &      SNB(M,NB),TB(M,NB),TCRX,TCTX,WSGWX,WSLCX,XGAB(M,NB),
     &      XGNB(M,NB),XGOB(M,NB),XGWB(M,NB),XNAB(M,NB),XNNB(M,NB),
     &      XNOB(M,NB),XNWB(M,NB),XMLAB(M,NB),XMLNB(M,NB),XMLOB(M,NB),
     &      XLSB(M,NB),YLSB(M,NB),YMGAX,YMGNX,YMGOX,YMNAX,YMNNX,YMNOX,
     &      ZMCX(1),ZMCX(2),ZMCX(3),IZN,N,NP_HT )
!
!---      Absolute pressures  ---
!
          PX = PSOB(M,NB) + PATM
!
!---      Aqueous component fractions and density  ---
!
          WTMLX = (XMLAB(M,NB)*WTMA + XMLOB(M,NB)*WTMO + 
     &      XMLNB(M,NB)*WTMN + 
     &      WTMW*(1.D+0-XMLAB(M,NB)-XMLOB(M,NB)-XMLNB(M,NB)))/
     &      (1.D+0-XLSB(M,NB)+WTMW*XLSB(M,NB)/WTMS)
          XLAB(M,NB) = XMLAB(M,NB)*WTMA/WTMLX
          XLOB(M,NB) = XMLOB(M,NB)*WTMO/WTMLX
          XLNB(M,NB) = XMLNB(M,NB)*WTMN/WTMLX
          CALL DENS_B( TB(M,NB),PX,XLSB(M,NB),RHOBX )
          CALL DENS_L( TB(M,NB),RHOBX,XLAB(M,NB),XLOB(M,NB),XLNB(M,NB),
     &      RHOLB(M,NB) )
          XLWB(M,NB) = 1.D+0-XLAB(M,NB)-XLOB(M,NB)-XLNB(M,NB)-XLSB(M,NB)
          IF( XLWB(M,NB).LT.EPSL ) XLWB(M,NB) = 0.D+0
          XMLSB(M,NB) = WTMLX*XLSB(M,NB)/WTMS
          XMLWB(M,NB) = WTMLX*XLWB(M,NB)/WTMW
!
!---      Aqueous molar density  ---
!
          RHOMLB(M,NB) = RHOLB(M,NB)/WTMLX
!
!---      Porous-media porosity  ---
!
          PPX = MAX( PGB(M,NB),PLB(M,NB),PNB(M,NB) ) + PATM
          CALL PORSTY_HYDT_KE( N,PPX,PCMP(N),PORDB(M,NB),
     &      PORTB(M,NB),SHB(M,NB) )
          PORDB(M,NB) = MAX( PORDB(M,NB),EPSL )
          PORTB(M,NB) = MAX( PORTB(M,NB),PORDB(M,NB) )
!
!---      Gas molar fractions and density  ---
!
          RHOGAX = XGAB(M,NB)*RHOGB(M,NB)
          RHOGOX = XGOB(M,NB)*RHOGB(M,NB)
          RHOGNX = XGNB(M,NB)*RHOGB(M,NB)
          RHOGWX = XGWB(M,NB)*RHOGB(M,NB)
          WTMGX = 1.D+0/(XGAB(M,NB)/WTMA + XGOB(M,NB)/WTMO
     &      + XGNB(M,NB)/WTMN + XGWB(M,NB)/WTMW)
          XMGAB(M,NB) = XGAB(M,NB)*WTMGX/WTMA
          XMGOB(M,NB) = XGOB(M,NB)*WTMGX/WTMO
          XMGNB(M,NB) = XGNB(M,NB)*WTMGX/WTMN
          XMGWB(M,NB) = XGWB(M,NB)*WTMGX/WTMW
          RHOMGB(M,NB) = RHOGB(M,NB)/WTMGX
!
!---      Gas viscosity, Pa s  ---
!
          RHOGAX = RHOMGB(M,NB)*WTMA
          CALL VISC_O2( TB(M,NB),RHOGAX,VISGAX )
          RHOGOX = RHOMGB(M,NB)*WTMO
          CALL VISC_CH4( TB(M,NB),RHOGOX,VISGOX )
          RHOGNX = RHOMGB(M,NB)*WTMN
          CALL VISC_N2( TB(M,NB),RHOGNX,VISGNX )
          RHOGWX = RHOMGB(M,NB)*WTMW
          CALL VISC_W( TB(M,NB),PX,RHOGWX,VISGWX )
          CALL VISC_G( PHIX,VISGAX,VISGOX,VISGNX,VISGWX,XMGAB(M,NB),
     &      XMGOB(M,NB),XMGNB(M,NB),XMGWB(M,NB),VISGB(M,NB) )
!
!---      Gas diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFGAB(M,NB) = DFGAC
            DFGOB(M,NB) = DFGOC
            DFGNB(M,NB) = DFGNC
            DFGWB(M,NB) = DFGWC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL DIFC_G( TB(M,NB),PX,DFGAB(M,NB),DFGOB(M,NB),
     &        DFGNB(M,NB),DFGWB(M,NB),XMGAB(M,NB),XMGOB(M,NB),
     &        XMGNB(M,NB),XMGWB(M,NB) )
          ENDIF
!
!---      Aqueous viscosity  ---
!
          ISRX = 1
          CALL DENS_W( TB(M,NB),PX,RHOLWX,RHOX,ISRX )
          CALL VISC_W( TB(M,NB),PX,RHOLWX,VISLWX )
          CALL VISC_B( TB(M,NB),XLSB(M,NB),VISLWX,VISBX )
          CALL VISC_L( XMLAB(M,NB),XMLOB(M,NB),XMLNB(M,NB),
     &      VISBX,VISGAX,VISGOX,VISGNX,VISLB(M,NB) )
!
!---      Aqueous diffusion coefficients  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLAB(M,NB) = DFLAC
            DFLOB(M,NB) = DFLOC
            DFLNB(M,NB) = DFLNC
            DFLSB(M,NB) = DFLSC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            CALL DIFC_LA( TB(M,NB),VISLB(M,NB),DFLAB(M,NB) )
            CALL DIFC_LO( TB(M,NB),VISLB(M,NB),DFLOB(M,NB) )
            CALL DIFC_LN( TB(M,NB),VISLB(M,NB),DFLNB(M,NB) )
            CALL DIFC_LS( TB(M,NB),XLSB(M,NB),VISLB(M,NB),DFLSB(M,NB) )
          ENDIF
!
!---      Nonaqueous-liquid molar fractions and density  ---
!
          RHONAX = XNAB(M,NB)*RHONB(M,NB)
          RHONOX = XNOB(M,NB)*RHONB(M,NB)
          RHONNX = XNNB(M,NB)*RHONB(M,NB)
          RHONWX = XNWB(M,NB)*RHONB(M,NB)
          WTMNX = 1.D+0/(XNAB(M,NB)/WTMA + XNOB(M,NB)/WTMO
     &       + XNNB(M,NB)/WTMN + XNWB(M,NB)/WTMW)
          XMNAB(M,NB) = XNAB(M,NB)*WTMNX/WTMA
          XMNOB(M,NB) = XNOB(M,NB)*WTMNX/WTMO
          XMNNB(M,NB) = XNNB(M,NB)*WTMNX/WTMN
          XMNWB(M,NB) = XNWB(M,NB)*WTMNX/WTMW
          RHOMNB(M,NB) = RHONB(M,NB)/WTMNX
!
!---      Nonaqueous-liquid viscosity, Pa s  ---
!
          RHONAX = RHOMNB(M,NB)*WTMA
          CALL VISC_O2( TB(M,NB),RHONAX,VISNAX )
          RHONOX = RHOMNB(M,NB)*WTMO
          CALL VISC_CH4( TB(M,NB),RHONOX,VISNOX )
          RHONNX = RHOMNB(M,NB)*WTMN
          CALL VISC_N2( TB(M,NB),RHONNX,VISNNX )
          RHONWX = RHOMNB(M,NB)*WTMW
          CALL VISC_W( TB(M,NB),PX,RHONWX,VISNWX )
          CALL VISC_G( PHIX,VISNAX,VISNOX,VISNNX,VISNWX,XMNAB(M,NB),
     &      XMNOB(M,NB),XMNNB(M,NB),XMNWB(M,NB),VISNB(M,NB) )
!
!---      Nonaqueous-liquid diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFNAB(M,NB) = DFNAC
            DFNOB(M,NB) = DFNOC
            DFNNB(M,NB) = DFNNC
            DFNWB(M,NB) = DFNWC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL DIFC_G( TB(M,NB),PX,DFNAB(M,NB),DFNOB(M,NB),
     &        DFNNB(M,NB),DFNWB(M,NB),XMNAB(M,NB),XMNOB(M,NB),
     &        XMNNB(M,NB),XMNWB(M,NB) )
          ENDIF
!
!---      Gas thermal conductivity, W/m K  ---
!
          RHOGAX = RHOMGB(M,NB)*WTMA
          CALL THK_O2( TB(M,NB),RHOGAX,THKGAX )
          RHOGOX = RHOMGB(M,NB)*WTMO
          CALL THK_CH4( TB(M,NB),RHOGOX,THKGOX )
          RHOGNX = RHOMGB(M,NB)*WTMN
          CALL THK_N2( TB(M,NB),RHOGNX,THKGNX )
          RHOGWX = RHOMGB(M,NB)*WTMW
          CALL THK_W( TB(M,NB),PX,RHOGWX,THKGWX )
          CALL THK_G( PHIX,THKGAX,THKGOX,THKGNX,THKGWX,XMGAB(M,NB),
     &      XMGOB(M,NB),XMGNB(M,NB),XMGWB(M,NB),THKGB(M,NB) )
!  
!---      Aqueous thermal conductivity  ---
!
          CALL THK_W( TB(M,NB),PX,RHOLWX,THKLWX )
          CALL THK_B( TB(M,NB),XLSB(M,NB),THKLWX,THKLB(M,NB) )
!
!---      Nonaqueous thermal conductivity, W/m K  ---
!
          RHONAX = RHOMNB(M,NB)*WTMA
          CALL THK_O2( TB(M,NB),RHONAX,THKNAX )
          RHONOX = RHOMNB(M,NB)*WTMO
          CALL THK_CH4( TB(M,NB),RHONOX,THKNOX )
          RHONNX = RHOMNB(M,NB)*WTMN
          CALL THK_N2( TB(M,NB),RHONNX,THKNNX )
          RHONWX = RHOMNB(M,NB)*WTMW
          CALL THK_W( TB(M,NB),PX,RHONWX,THKNWX )
          CALL THK_G( PHIX,THKNAX,THKNOX,THKNNX,THKNWX,XMNAB(M,NB),
     &      XMNOB(M,NB),XMNNB(M,NB),XMNWB(M,NB),THKNB(M,NB) )
!
!---      Pure water vapor enthalpy and internal energy  ---
!
          CALL SP_W( TB(M,NB),PSWX )
          PEX = MIN( PSWX,PX )
          ISRX = 2
          CALL ENTH_W( TB(M,NB),PEX,HLWX,HGWB(M,NB),ISRX )
!
!---      Gas enthalpy and internal energy  ---
!
          XMGX(1) = YMGAX
          XMGX(2) = YMGOX
          XMGX(3) = YMGNX
          CALL ENTH_ID( TB(M,NB),HIDX )
          INDX = 1
          CALL ENTH_RES( PX,TB(M,NB),HRESX,XMGX,INDX )
          HGAB(M,NB) = (HIDX(1) + HRESX(1))/WTMA
          HGOB(M,NB) = (HIDX(2) + HRESX(2))/WTMO
          HGNX = (HIDX(3) + HRESX(3))/WTMN
          HGB(M,NB) = XGAB(M,NB)*HGAB(M,NB) + XGOB(M,NB)*HGOB(M,NB) + 
     &      XGNB(M,NB)*HGNX + XGWB(M,NB)*HGWB(M,NB)
!
!---      Nonaqueous-liquid enthalpy and internal energy  ---
!
          XMNX(1) = YMNAX
          XMNX(2) = YMNOX
          XMNX(3) = YMNNX
          CALL ENTH_ID( TB(M,NB),HIDX )
          INDX = 2
          CALL ENTH_RES( PX,TB(M,NB),HRESX,XMNX,INDX )
          HNAX = (HIDX(1) + HRESX(1))/WTMA
          HNOX = (HIDX(2) + HRESX(2))/WTMO
          HNNX = (HIDX(3) + HRESX(3))/WTMN
          HNB(M,NB) = XNAB(M,NB)*HNAX + XNOB(M,NB)*HNOX + 
     &      XNNB(M,NB)*HNNX + XNWB(M,NB)*HGWB(M,NB)
!
!---      Pure water liquid enthalpy  ---
!
          CALL SP_W( TB(M,NB),PSWX )
          PEX = MAX( PSWX,PX )
          ISRX = 1
          CALL ENTH_W( TB(M,NB),PEX,HLWB(M,NB),HGWX,ISRX )
!
!---      Aqueous enthalpy  ---
!
          CALL ENTH_B( TB(M,NB),XLSB(M,NB),HLWB(M,NB),HBX )
          FGAX = SGB(M,NB)*RHOGB(M,NB)*XGAB(M,NB)
          FNAX = SNB(M,NB)*RHONB(M,NB)*XNAB(M,NB)
          IF( (FGAX+FNAX)/EPSL.LT.EPSL ) THEN
            HGNAX = 0.D+0
          ELSE
            HGNAX = (HGAB(M,NB)*FGAX + HNAX*FNAX)/(FGAX + FNAX)
          ENDIF
          FGOX = SGB(M,NB)*RHOGB(M,NB)*XGOB(M,NB)
          FNOX = SNB(M,NB)*RHONB(M,NB)*XNOB(M,NB)
          IF( (FGOX+FNOX)/EPSL.LT.EPSL ) THEN
            HGNOX = 0.D+0
          ELSE
            HGNOX = (HGOB(M,NB)*FGOX + HNOX*FNOX)/(FGOX + FNOX)
          ENDIF
          FGNX = SGB(M,NB)*RHOGB(M,NB)*XGNB(M,NB)
          FNNX = SNB(M,NB)*RHONB(M,NB)*XNNB(M,NB)
          IF( (FGNX+FNNX)/EPSL.LT.EPSL ) THEN
            HGNNX = 0.D+0
          ELSE
            HGNNX = (HGNX*FGNX + HNNX*FNNX)/(FGNX + FNNX)
          ENDIF
          CALL ENTH_L( TB(M,NB),XLSB(M,NB),XLAB(M,NB),XLOB(M,NB),
     &      XLNB(M,NB),HBX,HGNAX,HGNOX,HGNNX,HLB(M,NB) )
!
!---      Aqueous, gas, and nonaqueous-liquid tortuosity  ---
!
          PORDX = PORDB(M,NB)*(1.D+0-SHB(M,NB)-SIB(M,NB))
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SLB(M,NB),SGB(M,NB),
     &      SNB(M,NB),PORDX,TORLB(M,NB),TORGB(M,NB),TORNB(M,NB) )
!
!---      NaCl volumetric concentration  ---
!
          TMSB(M,NB) = YLSB(M,NB)*RHOLB(M,NB)*SLB(M,NB)*PORDB(M,NB)
  300   CONTINUE
  400 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCP_HYD_O2_CH4_N2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP_HYD_O2_CH4_N2( CPGLX,CPGNX,CPNLX,SHX,SLIX,SNX,
     &  IZN,N )
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
!     H2O-O2-CH4-NaCl-E-HYD1 Mode
!
!     Compute the gas-aqueous and hydrate aqueous capillary pressures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CAP_HYD_O2_CH4_N2'
      N_DB = N
!
!---  van Genuchten saturation function
!
      IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 .OR. 
     &  ISCHR(IZN).EQ.13 ) THEN
        CNX = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( MOD( IRPL(IZN),100 ).EQ.2 ) THEN
            CMX = 1.D+0 - 2.D+0/CNX
          ELSE
            CMX = 1.D+0 - 1.D+0/CNX
          ENDIF
        ELSE
          CMX = SCHR(14,IZN)
        ENDIF
        SRX = SCHR(4,IZN)
!
!---    Effective saturations  ---
!
        ESLX = SLIX/(1.D+0-SHX)
        ESNX = SNX/(1.D+0-SHX)
        ESTX = ESLX+ESNX
!
!---    Webb extension  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HMPX = SCHR(17,IZN)
          SMPX = SCHR(16,IZN)
!
!---      Total-liquid saturation below the matching point,
!         use Webb extension  ---
!
          IF( ESTX.LT.SMPX ) THEN
            DMPX = -(LOG10(HDOD)-LOG10(HMPX))/SMPX
            HDGN = 1.D+1**(DMPX*(ESTX-SMPX) + LOG10(HMPX))
!
!---      Total-liquid saturation at or above the matching point,
!         use van Genuchten function
!
          ELSE
            ASTX = (ESTX-SRX)/(1.D+0-SRX)
            VARX = (1.D+0/ASTX)**(1.D+0/CMX) - 1.D+0
            IF( VARX.LT.EPSL ) VARX = 0.D+0
            HDGN = (VARX**(1.D+0/CNX))/SCHR(1,IZN)
          ENDIF
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( ESLX.LT.SMPX ) THEN
            DMPX = -(LOG10(HDOD)-LOG10(HMPX))/SMPX
            HDNL = 1.D+1**(DMPX*(ESLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use van Genuchten function
!
          ELSE
            ASLX = (ESLX-SRX)/(1.D+0-SRX)
            VARX = (1.D+0/ASLX)**(1.D+0/CMX) - 1.D+0
            IF( VARX.LT.EPSL ) VARX = 0.D+0
            HDNL = (VARX**(1.D+0/CNX))/SCHR(1,IZN)
          ENDIF
!
!---    No extension  ---
!
        ELSE
          ASLX = (ESLX-SRX)/(1.D+0-SRX)
          ASNX = ESNX/(1.D+0-SRX)
          ASTX = ASLX+ASNX
          IF( ASTX.GT.EPSL ) THEN
            VARX = (1.D+0/ASTX)**(1.D+0/CMX) - 1.D+0
            IF( VARX.LT.EPSL ) VARX = 0.D+0
            HDGN = (VARX**(1.D+0/CNX))/SCHR(1,IZN)
          ELSE
            HDGN = 0.D+0
          ENDIF
          IF( ASLX.GT.EPSL ) THEN
            VARX = (1.D+0/ASLX)**(1.D+0/CMX) - 1.D+0
            IF( VARX.LT.EPSL ) VARX = 0.D+0
            HDNL = (VARX**(1.D+0/CNX))/SCHR(1,IZN)
          ELSE
            HDNL = 0.D+0
          ENDIF
        ENDIF
        HCPBFX = 1.D+0
!
!---    Entry pressure with hydrate barrier adjustment  ---
!
        IF( ISCHR(IZN).EQ.11 ) THEN
          HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &      SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
!
!---    Entry pressure with hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.13 ) THEN
          HCPBFX = (1.D+0-SHX)**(-(SCHR(10,IZN)-1.D+0)/SCHR(10,IZN))
        ENDIF
        CPGNX = HCPBFX*HDGN*RHORL*GRAV/BGN
        CPNLX = HCPBFX*HDNL*RHORL*GRAV/BNL
        CPGLX = CPGNX+CPNLX
!
!---  Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 .OR. 
     &  ISCHR(IZN).EQ.14 ) THEN
        CLX = MAX( SCHR(3,IZN),SMALL )
        SRX = SCHR(4,IZN)
!
!---    Effective saturations  ---
!
        ESLX = SLIX/(1.D+0-SHX)
        ESNX = SNX/(1.D+0-SHX)
        ESTX = ESLX+ESNX
!
!---    Webb extension  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HMPX = SCHR(17,IZN)
          SMPX = SCHR(16,IZN)
!
!---      Total-liquid saturation below the matching point,
!         use Webb extension  ---
!
          IF( ESTX.LT.SMPX ) THEN
            DMPX = -(LOG10(HDOD)-LOG10(HMPX))/SMPX
            HDGN = 1.D+1**(DMPX*(ESTX-SMPX) + LOG10(HMPX))
!
!---      Total-liquid saturation at or above the matching point,
!         use Brooks and Corey function
!
          ELSE
            ASTX = (ESTX-SRX)/(1.D+0-SRX)
            IF( ASTX.GT.EPSL ) THEN
              HDGN = SCHR(2,IZN)*(1.D+0/ASTX)**(1.D+0/CLX)
            ELSE
              HDGN = 0.D+0
            ENDIF
          ENDIF
!
!---      Aqueous saturation below the matching point,
!         use Webb extension  ---
!
          IF( ESLX.LT.SMPX ) THEN
            DMPX = -(LOG10(HDOD)-LOG10(HMPX))/SMPX
            HDNL = 1.D+1**(DMPX*(ESLX-SMPX) + LOG10(HMPX))
!
!---      Aqueous saturation at or above the matching point,
!         use Brooks and Corey function
!
          ELSE
            ASLX = (ESLX-SRX)/(1.D+0-SRX)
            IF( ASLX.GT.EPSL ) THEN
              HDNL = SCHR(2,IZN)*(1.D+0/ASLX)**(1.D+0/CLX)
            ELSE
              HDNL = 0.D+0
            ENDIF
          ENDIF
!
!---    No extension  ---
!
        ELSE
          ASLX = (ESLX-SRX)/(1.D+0-SRX)
          ASNX = ESNX/(1.D+0-SRX)
          ASTX = ASLX+ASNX
          IF( ASTX.GT.EPSL ) THEN
            HDGN = SCHR(2,IZN)*(1.D+0/ASTX)**(1.D+0/CLX)
          ELSE
            HDGN = 0.D+0
          ENDIF
          IF( ASLX.GT.EPSL ) THEN
            HDNL = SCHR(2,IZN)*(1.D+0/ASLX)**(1.D+0/CLX)
          ELSE
            HDNL = 0.D+0
          ENDIF
        ENDIF
        HCPBFX = 1.D+0
!
!---    Entry pressure with hydrate barrier adjustment  ---
!
        IF( ISCHR(IZN).EQ.12 ) THEN
          HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &      SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
!
!---    Entry pressure with hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.14 ) THEN
          HCPBFX = (1.D+0-SHX)**(-(SCHR(10,IZN)-1.D+0)/SCHR(10,IZN))
        ENDIF
        CPGNX = HCPBFX*HDGN*RHORL*GRAV/BGN
        CPNLX = HCPBFX*HDNL*RHORL*GRAV/BNL
        CPGLX = CPGNX+CPNLX
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Check the thermodynamic and hydrologic states declared through
!     user inputs.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE GEO_MECH
      USE FDVS
      USE FDVP
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/CHK_HYD_O2_CH4_N2'
      EPSLX = 1.D-4
!
!---  Component name  ---
!
      HCNM(1) = 'o2'
      HCNM(2) = 'ch4'
      HCNM(3) = 'n2'
      HCNM(4) = 'h2o'
!
!---  Component molecular weight, kg/kmol  ---
!
      HCPP(1,1) = 31.999D+0
      HCPP(1,2) = 16.0428D+0
      HCPP(1,3) = 28.01348D+0
      HCPP(1,4) = 18.015D+0
      WTMA = HCPP(1,1)
      WTMO = HCPP(1,2)
      WTMN = HCPP(1,3)
      WTMW = HCPP(1,4)
!
!---  Component normal boiling point, K  ---
!
      HCPP(2,1) = 90.188D+0
      HCPP(2,2) = 111.667D+0
      HCPP(2,3) = 77.355D+0
      HCPP(2,4) = 373.2D+0
!
!---  Component critical temperature, K  ---
!
      HCPP(3,1) = 154.581D+0
      HCPP(3,2) = 190.564D+0
      HCPP(3,3) = 126.192D+0
      HCPP(3,4) = 647.3D+0
!
!---  Component critical pressure, Pa  ---
!
      HCPP(4,1) = 5042.786D+3
      HCPP(4,2) = 4599.2D+3
      HCPP(4,3) = 3395.8D+3
      HCPP(4,4) = 221.2D+5
!
!---  Component triple-point temperature, K  ---
!
      HCPP(16,1) = 54.361D+0
      HCPP(16,2) = 90.694D+0
      HCPP(16,3) = 63.1526D+0
      HCPP(16,4) = 273.16D+0
!
!---  Component triple-point pressure, Pa  ---
!
      HCPP(17,1) = 0.1462776D+3
      HCPP(17,2) = 11.696D+3
      HCPP(17,3) = 12.523D+3
      HCPP(17,4) = 0.6117D+3
!
!---  Component critical molar volume, m^3/kmol ---
!
      HCPP(5,1) = 1.D+0/13.63D+0
      HCPP(5,2) = 1.D+0/10.139128D+0
      HCPP(5,3) = 1.D+0/11.1839D+0
      HCPP(5,4) = 57.1D-3
!
!---  Component critical compressibility  ---
!
      HCPP(6,1) = 0.288D+0
      HCPP(6,2) = 0.288D+0
      HCPP(6,3) = 0.292D+0
      HCPP(6,4) = 0.235D+0
!
!---  Component Pitzer acentric factor  ---
!
      HCPP(7,1) = 0.025D+0
      HCPP(7,2) = 0.01142D+0
      HCPP(7,3) = 0.0372D+0
      HCPP(7,4) = 0.344D+0
!
!---  Component dipole moment, debyes  ---
!
      HCPP(8,1) = 0.D+0
      HCPP(8,2) = 0.D+0
      HCPP(8,3) = 0.D+0
      HCPP(8,4) = 1.8D+0
!
!---  Liquid molar volume, m^3/kmol  ---
!
      HCPP(9,1) = 28.045D-3
      HCPP(9,2) = 24.429D-3
      HCPP(9,3) = 17.335D-3
      HCPP(9,4) = 18.016D-3
!
!---  Ideal gas heat capacity coefficient A  ---
!
      HCPP(18,1) = 2.811D+4
      HCPP(18,2) = 1.925D+4
      HCPP(18,3) = 3.115D+4
      HCPP(18,4) = 3.224D+4
!
!---  Ideal gas heat capacity coefficient B  ---
!
      HCPP(19,1) = -3.680D-3
      HCPP(19,2) = 5.213D+1
      HCPP(19,3) = -1.357D+1
      HCPP(19,4) = 1.924D+1
!
!---  Ideal gas heat capacity coefficient C  ---
!
      HCPP(20,1) = 1.746D-2
      HCPP(20,2) = 1.197D-2
      HCPP(20,3) = 2.680D-2
      HCPP(20,4) = 1.055D-2
!
!---  Ideal gas heat capacity coefficient D  ---
!
      HCPP(21,1) = -1.065D-5
      HCPP(21,2) = -1.132D-5
      HCPP(21,3) = -1.168D-5
      HCPP(21,4) = -3.596D-6
!
!---  Peng-Robinson pure component parameters (a, b, fw)  ---
!
      DO 80 I = 1,4
        HCPP(13,I) = 0.45724D+0*(RCU**2)*(HCPP(3,I)**2)/HCPP(4,I)
        HCPP(14,I) = 0.07780D+0*RCU*HCPP(3,I)/HCPP(4,I)
        HCPP(15,I) = 0.37464D+0 + 1.54226D+0*HCPP(7,I) -
     &    0.26992D+0*(HCPP(7,I)**2)
   80 CONTINUE
!
!---  Read O2-CH4-N2 equilibria file  ---
!
      CALL RDPF_HYD_O2_CH4_N2
!
!---  Check initial temperature, aqueous pressure, gas pressure,
!     and aqueous saturation  ---
!
      INDX = 0
      DO 100 N = 1,NFLD
        N_DB = N
        IF( IXP(N).EQ.0 ) GOTO 100
        IF( T(2,N).GT.374.14D+0 .OR. T(2,N).LT.-50.D+0 ) THEN
          INDX = 16
          IMSG = N
          RLMSG = T(2,N)
          CHMSG = 'Out of Range Initial Temperature(C) @ Node'
          CALL WRMSGS( INDX )
        ENDIF
        IF( PSO(2,N).GT.8.D+8-PATM ) THEN
          INDX = 16
          IMSG = N
          CHMSG = 'Out of Range Initial Pressure(Pa) @ Node'
          RLMSG = PSO(2,N)+PATM
          CALL WRMSGS( INDX )
        ENDIF
  100 CONTINUE
      IF( INDX.GT.0 ) STOP
!
!---  Webb saturation and capillary pressure matching points  ---
!
      DO 120 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 120
        IZN = IZ(N)
        IF( ISM(IZN).EQ.2 ) THEN
!
!---      van Genuchten moisture retension function  ---
!
          IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11
     &      .OR. ISCHR(IZN).EQ.13 ) THEN
            IF( SCHR(14,IZN).LE.0.D+0 ) THEN
              CNX = MAX( SCHR(3,IZN),SMALL )
              IF( IRPN(IZN).EQ.2 ) THEN
                SCHR(14,IZN) = 1.D+0 - 2.D+0/CNX
              ELSE
                SCHR(14,IZN) = 1.D+0 - 1.D+0/CNX
              ENDIF
            ENDIF
            CALL WEBB_VG( N )
!
!---      Brooks and Corey moisture retension function  ---
!
          ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12
     &      .OR. ISCHR(IZN).EQ.14 ) THEN
            CALL WEBB_BC( N )
          ENDIF
        ENDIF
  120 CONTINUE
!
!---  Skip flash calculations for restart simulations  ---
!
      IF( IEO.EQ.2 ) GOTO 400
!
!---  Convert initial conditions through flash calculations  ---
!
      DO 300 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 300
        IZN = IZ(N)
        N_DB = N
!
!---    Restart simulation with overwrite  ---
!
        IF( IEO.EQ.6 ) THEN
!
!---      Phase conditions 21 through 23 (20 series)  ---
!
!         SH = 0, SG + SN > 0, SL + SI < 1
!         ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #21
!         ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #22
!         ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #23
!
!         energy - temperature, C - T
!         water mass - aqueous + ice saturation - POSM
!         maximum mobile mass component - pressure, Pa (gauge) - PSO
!         mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!         mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!         mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!         hydrate O2 mass - hydrate O2 mass, kg - TMHA
!         hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!         hydrate N2 mass - hydrate N2 mass, kg - TMHN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          IF( NPHAZ(2,N).GE.21 .AND. NPHAZ(2,N).LE.23 ) THEN
!
!---        Nonaqueous O2 or CH4 or N2 mole fraction of formers 
!           specified  ---
!
            IF( IC_OPT(8,N).EQ.2 .OR. IC_OPT(9,N).EQ.2 .OR.
     &        IC_OPT(10,N).EQ.2 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: Restart w/ Overwrite: ' // 
     &          '20 Series Phase Condition: Nonaqueous ' //
     &          'O2 or CH4 or N2 Mole Fraction of Formers Specified.'
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Mobile O2 + CH4 + N2 mole fraction of formers specified  ---
!
            IF( IC_OPT(8,N).EQ.1 .AND. IC_OPT(9,N).EQ.1 .AND.
     &        IC_OPT(10,N).EQ.1 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: Restart w/ Overwrite: ' // 
     &          '20 Series Phase Condition: Mobile ' //
     &          'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Mobile O2 or CH4 mole fraction of formers unspecified  ---
!
            IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: Restart w/ Overwrite: ' // 
     &          '20 Series Phase Condition: Mobile ' //
     &          'O2 or CH4 Mole Fraction of Formers Unspecified.'
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Mobile O2 or N2 mole fraction of formers unspecified  ---
!
            IF( (IC_OPT(8,N)+IC_OPT(10,N)).EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: Restart w/ Overwrite: ' // 
     &          '20 Series Phase Condition: Mobile ' //
     &          'O2 or N2 Mole Fraction of Formers Unspecified.'
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Mobile CH4 or N2 mole fraction of formers unspecified  ---
!
            IF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.0 ) THEN
              INDX = 4
              CHMSG = 'Specification Error: Restart w/ Overwrite: ' // 
     &          '20 Series Phase Condition: Mobile ' //
     &          'CH4 or N2 Mole Fraction of Formers Unspecified.'
              CALL WRMSGS( INDX )
            ENDIF
            IF( IC_OPT(6,N).EQ.1 ) POSM(2,N) = SL(2,N)
            IF( IC_OPT(8,N).EQ.1 ) ZMCA(2,N) = XMGA(2,N)
            IF( IC_OPT(9,N).EQ.1 ) ZMCO(2,N) = XMGO(2,N)
            IF( IC_OPT(10,N).EQ.1 ) ZMCN(2,N) = XMGN(2,N)
            IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.2 ) THEN
              ZMCN(2,N) = MAX( 1.D+0-ZMCA(2,N)-ZMCO(2,N),0.D+0 )
              IF( ZMCN(2,N).LT.EPSL ) ZMCN(2,N) = 0.D+0
            ELSEIF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.2 ) THEN
              ZMCA(2,N) = MAX( 1.D+0-ZMCO(2,N)-ZMCN(2,N),0.D+0 )
              IF( ZMCA(2,N).LT.EPSL ) ZMCA(2,N) = 0.D+0
            ELSEIF( (IC_OPT(10,N)+IC_OPT(8,N)).EQ.2 ) THEN
              ZMCO(2,N) = MAX( 1.D+0-ZMCN(2,N)-ZMCA(2,N),0.D+0 )
              IF( ZMCO(2,N).LT.EPSL ) ZMCO(2,N) = 0.D+0
            ENDIF
          ENDIF
          CYCLE
        ENDIF
!
!---    Initial condition options
!
!       IC_OPT(1,N) initial condition state option
!       IC_OPT(2,N) 0 hydrate saturation not specified
!       IC_OPT(2,N) 1 hydrate saturation specified
!       IC_OPT(3,N) 0 hydrate O2 mole frac. of formers not specified
!       IC_OPT(3,N) 1 hydrate O2 mole frac. of formers specified
!       IC_OPT(4,N) 0 hydrate CH4 mole frac. of formers not specified
!       IC_OPT(4,N) 1 hydrate CH4 mole frac. of formers specified
!       IC_OPT(5,N) 0 hydrate N2 mole frac. of formers not specified
!       IC_OPT(5,N) 1 hydrate N2 mole frac. of formers specified
!       IC_OPT(6,N) 0 aqueous saturation not specified
!       IC_OPT(6,N) 1 aqueous saturation specified
!       IC_OPT(8,N) 0 mobile O2 mole frac. of formers not specified
!       IC_OPT(8,N) 1 mobile O2 mole frac. of formers specified
!       IC_OPT(9,N) 0 mobile CH4 mole frac. of formers not specified
!       IC_OPT(9,N) 1 mobile CH4 mole frac. of formers specified
!       IC_OPT(10,N) 0 mobile N2 mole frac. of formers not specified
!       IC_OPT(10,N) 1 mobile N2 mole frac. of formers specified
!       IC_OPT(8,N) 0 total nonaqueous O2 mole fraction not specified
!       IC_OPT(8,N) 2 total nonaqueous O2 mole fraction specified
!       IC_OPT(9,N) 0 total nonaqueous CH4 mole fraction not specified
!       IC_OPT(9,N) 2 total nonaqueous CH4 mole fraction specified
!       IC_OPT(10,N) 0 total nonaqueous N2 mole fraction not specified
!       IC_OPT(10,N) 2 total nonaqueous N2 mole fraction specified
!       IC_OPT(11,N) 0 aqu. relative saturation of formers not specified
!       IC_OPT(11,N) 1 aqu. relative saturation of formers specified
!       IC_OPT(12,N) salt concentration option
!       IC_OPT(13,N) 0 temperature not specified
!       IC_OPT(13,N) 1 temperature specified
!       IC_OPT(14,N) 0 pressure not specified
!       IC_OPT(14,N) 1 pressure specified
!       IC_OPT(15,N) 0 total nonaqueous saturation not specified
!       IC_OPT(15,N) 1 total nonaqueous saturation specified
!
!---    Initial condition state not declared  ---
!
        IF( IC_OPT(1,N).LT.1 .OR. IC_OPT(1,N).GT.7 ) THEN
          INDX = 7
          IMSG = N_DB
          CHMSG = 'Specification Error: ' // 
     &      'Unspecified Initial Condition State: Node'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Salt concentration undeclared  ---
!
        IF( IC_OPT(12,N).EQ.0 ) IC_OPT(12,N) = 3
!
!---    Aqueous salt relative saturation or
!       aqueous salt mass fraction  ---
!
        IF( IC_OPT(12,N).EQ.2 .OR. IC_OPT(12,N).EQ.3 ) THEN
!
!---      Aqueous salt relative saturation  ---
!
          IF( IC_OPT(12,N).EQ.2 ) THEN
            PHILSX = TMS(2,N)
            CALL SOL_LS( T(2,N),XLSMX )
            YLS(2,N) = PHILSX*XLSMX
            XLS(2,N) = MIN( YLS(2,N),XLSMX )
!
!---      Aqueous salt mass fraction  ---
!
          ELSEIF( IC_OPT(12,N).EQ.3 ) THEN
            YLS(2,N) = TMS(2,N)
            CALL SOL_LS( T(2,N),XLSMX )
            XLS(2,N) = MIN( YLS(2,N),XLSMX )
          ENDIF
        ELSE
          INDX = 4
          CHMSG = 'Uncoded Initial Condition Option:  ' //
     &      'Dissolved-Salt Volumetric Concentration.'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Hydrate saturation declared for a no-hydrate rock  ---
!
        IF( INHYD(IZN).EQ.1 .AND. SH(2,N).GT.0.D+0 ) THEN
          INDX = 21
          RLMSG = SH(2,N)
          CHMSG = 'Hydrate Saturation in No-Hydrate Rock: ' //
     &      ROCK(IZN) // 'SH = '
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Initial condition state #1
!
!       SHX = 0.0
!       SGX + SNX = 0.0
!       SLX = 1.0
!       Hydrate-mobile phase non-equilibrium
!
!       Declared variables:
!
!       T - temperature, C (T/T)
!       PSO - pressure, Pa (PSO/PSO)
!       2/3 YMMA, YMMO, YMMN
!       YMMA - mobile O2 mole fraction of formers (XMGA/PVA)
!       YMMO - mobile CH4 mole fraction of formers (XMGO/PVO)
!       YMMN - mobile N2 mole fraction of formers (XMGN/PVN)
!       or
!       ZMCA - total nonaqueous O2 mole fraction (ZMCA/PVA)
!       ZMCO - total nonaqueous CH4 mole fraction (ZMCO/PVO)
!       ZMCN - total nonaqueous N2 mole fraction (ZMCN/PVN)
!       WLF - aqueous relative saturation of formers (XMLO/PVW)
!
        IF( IC_OPT(1,N).EQ.1 ) THEN
!
!---      Hydrate saturation specified  ---
!
          IF( IC_OPT(2,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: ' //
     &        'Hydrate Saturation Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 mole fraction of formers specified  ---
!
          IF( IC_OPT(3,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: ' //
     &        'Hydrate O2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate CH4 mole fraction of formers specified  ---
!
          IF( IC_OPT(4,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: ' //
     &        'Hydrate CH4 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate N2 mole fraction of formers specified  ---
!
          IF( IC_OPT(5,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: ' //
     &        'Hydrate N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous saturation specified  ---
!
          IF( IC_OPT(6,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: ' //
     &        'Aqueous Saturation Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 + CH4 + N2 mole fraction of formers specified  ---
!
          IF( IC_OPT(8,N).EQ.1 .AND. IC_OPT(9,N).EQ.1 .AND.
     &      IC_OPT(10,N).EQ.1 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: Mobile ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous O2 + CH4 + N2 mole fraction of formers 
!         specified  ---
!
          IF( IC_OPT(8,N).EQ.2 .AND. IC_OPT(9,N).EQ.2 .AND.
     &      IC_OPT(10,N).EQ.2 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: Nonaqueous ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile or nonaqueous O2 or CH4 mole fraction of formers 
!         unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: Mobile ' //
     &        'O2 or CH4 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile or nonaqueous O2 or N2 mole fraction of formers 
!         unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: Mobile ' //
     &        'O2 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile or nonaqueous CH4 or N2 mole fraction of formers 
!         unspecified  ---
!
          IF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: Mobile ' //
     &        'CH4 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Relative saturation of dissolved formers not specified  ---
!
          IF( IC_OPT(11,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: ' //
     &        'Relative Saturation of Dissolved Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Temperature not specified  ---
!
          IF( IC_OPT(13,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: ' //
     &        'Temperature Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Pressure not specified  ---
!
          IF( IC_OPT(14,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 1: ' //
     &        'Pressure Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Flash calculation for initial condition #1  ---
!
          PVA(2,N) = XMGA(2,N)
          PVO(2,N) = XMGO(2,N)
          PVN(2,N) = XMGN(2,N)
          PVW(2,N) = XMLO(2,N)
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.2 ) THEN
            PVN(2,N) = MAX( 1.D+0-PVA(2,N)-PVO(2,N),0.D+0 )
            IF( PVN(2,N).LT.EPSL ) PVN(2,N) = 0.D+0
          ELSEIF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.2 ) THEN
            PVA(2,N) = MAX( 1.D+0-PVO(2,N)-PVN(2,N),0.D+0 )
            IF( PVA(2,N).LT.EPSL ) PVA(2,N) = 0.D+0
          ELSEIF( (IC_OPT(10,N)+IC_OPT(8,N)).EQ.2 ) THEN
            PVO(2,N) = MAX( 1.D+0-PVN(2,N)-PVA(2,N),0.D+0 )
            IF( PVO(2,N).LT.EPSL ) PVO(2,N) = 0.D+0
          ENDIF
          SL(2,N) = 1.D+0
          CALL FLH_IC1( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),PVN(2,N),
     &      PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),RHOH(2,N),
     &      SH(2,N),SI(2,N),SG(2,N),SL(2,N),SN(2,N),
     &      T(2,N),TCR(N),TCT(N),TEQHX,
     &      XHA(2,N),XHN(2,N),XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),
     &      YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),YMHGN(2,N),
     &      YMHGO(2,N),ZMCA(2,N),ZMCN(2,N),ZMCO(2,N),
     &      IZN,N,NP_HT )
          POSM(2,N) = SL(2,N) + SI(2,N)
          IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
            NPHAZ(2,N) = 10 + (100*NP_HT)
          ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
            NPHAZ(2,N) = 21 + (100*NP_HT)
          ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
            NPHAZ(2,N) = 22 + (100*NP_HT)
          ELSE
            NPHAZ(2,N) = 23 + (100*NP_HT)
          ENDIF
!
!---    Initial condition state #2
!
!       SHX = 0.0
!       SGX + SNX > 0.0
!       SLX + SIX < 1.0
!       Hydrate-mobile phase non-equilibrium
!
!       Declared variables:
!
!       T - temperature, C (T/T)
!       PSO - pressure, Pa (PSO/PSO)
!       SL - aqueous saturation (SL/SL)
!       2/3 YMMA, YMMO, YMMN
!       YMMA - mobile O2 mole fraction of formers (XMGA/PVA)
!       YMMO - mobile CH4 mole fraction of formers (XMGO/PVO)
!       YMMN - mobile N2 mole fraction of formers (XMGN/PVN)
!
        ELSEIF( IC_OPT(1,N).EQ.2 ) THEN
!
!---      Hydrate saturation specified  ---
!
          IF( IC_OPT(2,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: ' //
     &        'Hydrate Saturation Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 mole fraction of formers specified  ---
!
          IF( IC_OPT(3,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: ' //
     &        'Hydrate O2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate CH4 mole fraction of formers specified  ---
!
          IF( IC_OPT(4,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: ' //
     &        'Hydrate CH4 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate N2 mole fraction of formers specified  ---
!
          IF( IC_OPT(5,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: ' //
     &        'Hydrate N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous saturation not specified  ---
!
          IF( IC_OPT(6,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: ' //
     &        'Aqueous Saturation Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous O2 or CH4 or N2 mole fraction of formers 
!         specified  ---
!
          IF( IC_OPT(8,N).EQ.2 .OR. IC_OPT(9,N).EQ.2 .OR.
     &      IC_OPT(10,N).EQ.2 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: Nonaqueous ' //
     &        'O2 or CH4 or N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 + CH4 + N2 mole fraction of formers specified  ---
!
          IF( IC_OPT(8,N).EQ.1 .AND. IC_OPT(9,N).EQ.1 .AND.
     &      IC_OPT(10,N).EQ.1 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: Mobile ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 or CH4 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: Mobile ' //
     &        'O2 or CH4 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: Mobile ' //
     &        'O2 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile CH4 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: Mobile ' //
     &        'CH4 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Relative saturation of dissolved formers specified  ---
!
          IF( IC_OPT(11,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: ' //
     &        'Relative Saturation of Dissolved Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Temperature not specified  ---
!
          IF( IC_OPT(13,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: ' //
     &        'Temperature Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Pressure not specified  ---
!
          IF( IC_OPT(14,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 2: ' //
     &        'Pressure Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Flash calculation for initial condition #2 ---
!
          PVA(2,N) = XMGA(2,N)
          PVO(2,N) = XMGO(2,N)
          PVN(2,N) = XMGN(2,N)
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.2 ) THEN
            PVN(2,N) = MAX( 1.D+0-PVA(2,N)-PVO(2,N),0.D+0 )
            IF( PVN(2,N).LT.EPSL ) PVN(2,N) = 0.D+0
          ELSEIF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.2 ) THEN
            PVA(2,N) = MAX( 1.D+0-PVO(2,N)-PVN(2,N),0.D+0 )
            IF( PVA(2,N).LT.EPSL ) PVA(2,N) = 0.D+0
          ELSEIF( (IC_OPT(10,N)+IC_OPT(8,N)).EQ.2 ) THEN
            PVO(2,N) = MAX( 1.D+0-PVN(2,N)-PVA(2,N),0.D+0 )
            IF( PVO(2,N).LT.EPSL ) PVO(2,N) = 0.D+0
          ENDIF
          CALL FLH_IC2( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),PVN(2,N),
     &      PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),RHOH(2,N),
     &      SH(2,N),SI(2,N),SG(2,N),SL(2,N),SN(2,N),
     &      T(2,N),TCR(N),TCT(N),TEQHX,
     &      XHA(2,N),XHN(2,N),XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),
     &      YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),YMHGN(2,N),
     &      YMHGO(2,N),ZMCA(2,N),ZMCN(2,N),ZMCO(2,N),
     &      IZN,N,NP_HT )
          POSM(2,N) = SL(2,N) + SI(2,N)
          IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
            NPHAZ(2,N) = 10 + (100*NP_HT)
          ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
            NPHAZ(2,N) = 21 + (100*NP_HT)
          ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
            NPHAZ(2,N) = 22 + (100*NP_HT)
          ELSE
            NPHAZ(2,N) = 23 + (100*NP_HT)
          ENDIF
!
!---    Initial condition state #3
!
!       SHX > 0.0
!       SGX + SNX = 0.0
!       SLX < 1.0
!       Hydrate-mobile phase non-equilibrium
!
!       Declared variables:
!
!       T - temperature, C (T/T)
!       PSO - pressure, Pa (PSO/PSO)
!       2/3 YMMA, YMMO, YMMN
!       YMMA - mobile O2 mole fraction of formers (XMGA/PVA)
!       YMMO - mobile CH4 mole fraction of formers (XMGO/PVO)
!       YMMN - mobile N2 mole fraction of formers (XMGN/PVN)
!       or
!       ZMCA - total nonaqueous O2 mole fraction (ZMCA/PVA)
!       ZMCO - total nonaqueous CH4 mole fraction (ZMCO/PVO)
!       ZMCN - total nonaqueous N2 mole fraction (ZMCN/PVN)
!       WLF - aqueous relative saturation of formers (XMLO/PVA)
!       SH - hydrate saturation (SH/SH)
!       2/3 YMHA, YMHO, YMHN
!       YMHA - hydrate O2 mole fraction of formers (XMGA/PVA)
!       YMHO - hydrate CH4 mole fraction of formers (XMGO/PVO)
!       YMHN - hydrate N2 mole fraction of formers (XMGN/PVN)
!
        ELSEIF( IC_OPT(1,N).EQ.3 ) THEN
!
!---      Hydrate saturation not specified  ---
!
          IF( IC_OPT(2,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: ' //
     &        'Hydrate Saturation Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 + CH4 + N2 mole fraction of formers specified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(4,N)+IC_OPT(5,N)).EQ.3 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: Hydrate ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 or CH4 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(4,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: Hydrate ' //
     &        'O2 or CH4 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(5,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: Hydrate ' //
     &        'O2 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate CH4 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(4,N)+IC_OPT(5,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: Hydrate ' //
     &        'CH4 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous saturation specified  ---
!
          IF( IC_OPT(6,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: ' //
     &        'Aqueous Saturation Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 + CH4 + N2 mole fraction of formers specified  ---
!
          IF( IC_OPT(8,N).EQ.1 .AND. IC_OPT(9,N).EQ.1 .AND.
     &      IC_OPT(10,N).EQ.1 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: Mobile ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous O2 + CH4 + N2 mole fraction of formers 
!         specified  ---
!
          IF( IC_OPT(8,N).EQ.2 .AND. IC_OPT(9,N).EQ.2 .AND.
     &      IC_OPT(10,N).EQ.2 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: Nonaqueous ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile or nonaqueous O2 or CH4 mole fraction of formers 
!         unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: Mobile ' //
     &        'O2 or CH4 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile or nonaqueous O2 or N2 mole fraction of formers 
!         unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: Mobile ' //
     &        'O2 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile or nonaqueous CH4 or N2 mole fraction of formers 
!         unspecified  ---
!
          IF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: Mobile ' //
     &        'CH4 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Relative saturation of dissolved formers not specified  ---
!
          IF( IC_OPT(11,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: ' //
     &        'Relative Saturation of Dissolved Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Temperature not specified  ---
!
          IF( IC_OPT(13,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: ' //
     &        'Temperature Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Pressure not specified  ---
!
          IF( IC_OPT(14,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 3: ' //
     &        'Pressure Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Flash calculation for initial condition #3 ---
!
          PVA(2,N) = XMGA(2,N)
          PVO(2,N) = XMGO(2,N)
          PVN(2,N) = XMGN(2,N)
          PVW(2,N) = XMLO(2,N)
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.2 ) THEN
            PVN(2,N) = MAX( 1.D+0-PVA(2,N)-PVO(2,N),0.D+0 )
            IF( PVN(2,N).LT.EPSL ) PVN(2,N) = 0.D+0
          ELSEIF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.2 ) THEN
            PVA(2,N) = MAX( 1.D+0-PVO(2,N)-PVN(2,N),0.D+0 )
            IF( PVA(2,N).LT.EPSL ) PVA(2,N) = 0.D+0
          ELSEIF( (IC_OPT(10,N)+IC_OPT(8,N)).EQ.2 ) THEN
            PVO(2,N) = MAX( 1.D+0-PVN(2,N)-PVA(2,N),0.D+0 )
            IF( PVO(2,N).LT.EPSL ) PVO(2,N) = 0.D+0
          ENDIF
          PVHA(2,N) = XHA(2,N)
          PVHO(2,N) = XHO(2,N)
          PVHN(2,N) = XHN(2,N)
          IF( (IC_OPT(3,N)+IC_OPT(4,N)).EQ.2 ) THEN
            PVHN(2,N) = MAX( 1.D+0-PVHA(2,N)-PVHO(2,N),0.D+0 )
            IF( PVHN(2,N).LT.EPSL ) PVHN(2,N) = 0.D+0
          ELSEIF( (IC_OPT(4,N)+IC_OPT(5,N)).EQ.2 ) THEN
            PVHA(2,N) = MAX( 1.D+0-PVHO(2,N)-PVHN(2,N),0.D+0 )
            IF( PVHA(2,N).LT.EPSL ) PVHA(2,N) = 0.D+0
          ELSEIF( (IC_OPT(5,N)+IC_OPT(3,N)).EQ.2 ) THEN
            PVHO(2,N) = MAX( 1.D+0-PVHN(2,N)-PVHA(2,N),0.D+0 )
            IF( PVHO(2,N).LT.EPSL ) PVHO(2,N) = 0.D+0
          ENDIF
          CALL FLH_IC3( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),PVN(2,N),
     &      PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),RHOH(2,N),
     &      SH(2,N),SI(2,N),SG(2,N),SL(2,N),SN(2,N),
     &      T(2,N),TCR(N),TCT(N),TEQHX,TMHA(2,N),TMHN(2,N),TMHO(2,N),
     &      XHA(2,N),XHN(2,N),XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),
     &      YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),YMHGN(2,N),
     &      YMHGO(2,N),ZMCA(2,N),ZMCN(2,N),ZMCO(2,N),
     &      IZN,N,NP_HT )
          POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      Zero hydrate saturation  ---
!
          IF( SH(2,N).LT.EPSL ) THEN
            IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
              NPHAZ(2,N) = 10 + (100*NP_HT)
            ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              NPHAZ(2,N) = 21 + (100*NP_HT)
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              NPHAZ(2,N) = 22 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 23 + (100*NP_HT)
            ENDIF
!
!---      Non-zero hydrate saturation  ---
!
          ELSE
            IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 41 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 42 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 43 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 31 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 32 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 33 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 34 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 35 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 36 + (100*NP_HT)
              ENDIF
            ELSE
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 37 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 38 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 39 + (100*NP_HT)
              ENDIF
            ENDIF
          ENDIF
!
!---    Initial condition state #4
!
!       SHX > 0.0
!       SGX + SNX > 0.0
!       SLX + SIX < 1.0
!       Hydrate-mobile phase non-equilibrium
!
!       Declared variables:
!
!       T - temperature, C (T/T)
!       PSO - pressure, Pa (PSO/PSO)
!       SL - aqueous saturation (SL/SL)
!       2/3 YMMA, YMMO, YMMN
!       YMMA - mobile O2 mole fraction of formers (XMGA/PVA)
!       YMMO - mobile CH4 mole fraction of formers (XMGO/PVO)
!       YMMN - mobile N2 mole fraction of formers (XMGN/PVN)
!       SH - hydrate saturation (SH/SH)
!       2/3 YMHA, YMHO, YMHN
!       YMHA - hydrate O2 mole fraction of formers (XMGA/PVA)
!       YMHO - hydrate CH4 mole fraction of formers (XMGO/PVO)
!       YMHN - hydrate N2 mole fraction of formers (XMGN/PVN)
!
        ELSEIF( IC_OPT(1,N).EQ.4 ) THEN
!
!---      Hydrate saturation not specified  ---
!
          IF( IC_OPT(2,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: ' //
     &        'Hydrate Saturation Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 + CH4 + N2 mole fraction of formers specified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(4,N)+IC_OPT(5,N)).EQ.3 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: Hydrate ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 or CH4 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(4,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: Hydrate ' //
     &        'O2 or CH4 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(5,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: Hydrate ' //
     &        'O2 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate CH4 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(4,N)+IC_OPT(5,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: Hydrate ' //
     &        'CH4 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous saturation unspecified  ---
!
          IF( IC_OPT(6,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: ' //
     &        'Aqueous Saturation Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous O2 or CH4 or N2 mole fraction of formers 
!         specified  ---
!
          IF( IC_OPT(8,N).EQ.2 .OR. IC_OPT(9,N).EQ.2 .OR.
     &      IC_OPT(10,N).EQ.2 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: Nonaqueous ' //
     &        'O2 or CH4 or N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 + CH4 + N2 mole fraction of formers specified  ---
!
          IF( IC_OPT(8,N).EQ.1 .AND. IC_OPT(9,N).EQ.1 .AND.
     &      IC_OPT(10,N).EQ.1 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: Mobile ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 or CH4 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: Mobile ' //
     &        'O2 or CH4 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: Mobile ' //
     &        'O2 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile CH4 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: Mobile ' //
     &        'CH4 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Relative saturation of dissolved formers specified  ---
!
          IF( IC_OPT(11,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: ' //
     &        'Relative Saturation of Dissolved Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Temperature not specified  ---
!
          IF( IC_OPT(13,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: ' //
     &        'Temperature Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Pressure not specified  ---
!
          IF( IC_OPT(14,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 4: ' //
     &        'Pressure Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Flash calculation for initial condition #4 ---
!
          PVA(2,N) = XMGA(2,N)
          PVO(2,N) = XMGO(2,N)
          PVN(2,N) = XMGN(2,N)
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.2 ) THEN
            PVN(2,N) = MAX( 1.D+0-PVA(2,N)-PVO(2,N),0.D+0 )
            IF( PVN(2,N).LT.EPSL ) PVN(2,N) = 0.D+0
          ELSEIF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.2 ) THEN
            PVA(2,N) = MAX( 1.D+0-PVO(2,N)-PVN(2,N),0.D+0 )
            IF( PVA(2,N).LT.EPSL ) PVA(2,N) = 0.D+0
          ELSEIF( (IC_OPT(10,N)+IC_OPT(8,N)).EQ.2 ) THEN
            PVO(2,N) = MAX( 1.D+0-PVN(2,N)-PVA(2,N),0.D+0 )
            IF( PVO(2,N).LT.EPSL ) PVO(2,N) = 0.D+0
          ENDIF
          PVHA(2,N) = XHA(2,N)
          PVHO(2,N) = XHO(2,N)
          PVHN(2,N) = XHN(2,N)
          IF( (IC_OPT(3,N)+IC_OPT(4,N)).EQ.2 ) THEN
            PVHN(2,N) = MAX( 1.D+0-PVHA(2,N)-PVHO(2,N),0.D+0 )
            IF( PVHN(2,N).LT.EPSL ) PVHN(2,N) = 0.D+0
          ELSEIF( (IC_OPT(4,N)+IC_OPT(5,N)).EQ.2 ) THEN
            PVHA(2,N) = MAX( 1.D+0-PVHO(2,N)-PVHN(2,N),0.D+0 )
            IF( PVHA(2,N).LT.EPSL ) PVHA(2,N) = 0.D+0
          ELSEIF( (IC_OPT(5,N)+IC_OPT(3,N)).EQ.2 ) THEN
            PVHO(2,N) = MAX( 1.D+0-PVHN(2,N)-PVHA(2,N),0.D+0 )
            IF( PVHO(2,N).LT.EPSL ) PVHO(2,N) = 0.D+0
          ENDIF
          CALL FLH_IC4( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),PVN(2,N),
     &      PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),RHOH(2,N),
     &      SH(2,N),SI(2,N),SG(2,N),SL(2,N),SN(2,N),
     &      T(2,N),TCR(N),TCT(N),TEQHX,
     &      XHA(2,N),XHN(2,N),XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),
     &      YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),YMHGN(2,N),
     &      YMHGO(2,N),YMNAX,YMNNX,YMNOX,ZMCA(2,N),ZMCN(2,N),ZMCO(2,N),
     &      IZN,N,NP_HT )
          POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      Zero hydrate saturation  ---
!
          IF( SH(2,N).LT.EPSL ) THEN
            IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
              NPHAZ(2,N) = 10 + (100*NP_HT)
            ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              NPHAZ(2,N) = 21 + (100*NP_HT)
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              NPHAZ(2,N) = 22 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 23 + (100*NP_HT)
            ENDIF
!
!---      Non-zero hydrate saturation  ---
!
          ELSE
            IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 41 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 42 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 43 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 31 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 32 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 33 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 34 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 35 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 36 + (100*NP_HT)
              ENDIF
            ELSE
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 37 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 38 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 39 + (100*NP_HT)
              ENDIF
            ENDIF
          ENDIF
!
!---    Initial condition state #5
!
!       SHX > 0.0
!       SGX + SNX > 0.0
!       SLX + SIX < 1.0
!       Hydrate-mobile phase equilibrium
!
!       Declared variables:
!
!       T - temperature, C (T/T)
!       PSO - pressure, Pa (PSO/PSO) (optional)
!       SL - aqueous saturation (SL/SL)
!       SH - hydrate saturation (SH/SH)
!       2/3 YMHA, YMHO, YMHN
!       YMHA - hydrate O2 mole fraction of formers (XMGA/PVA)
!       YMHO - hydrate CH4 mole fraction of formers (XMGO/PVO)
!       YMHN - hydrate N2 mole fraction of formers (XMGN/PVN)
!
        ELSEIF( IC_OPT(1,N).EQ.5 ) THEN
!
!---      Hydrate saturation not specified  ---
!
          IF( IC_OPT(2,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: ' //
     &        'Hydrate Saturation Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 + CH4 + N2 mole fraction of formers specified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(4,N)+IC_OPT(5,N)).EQ.3 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: Hydrate ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 or CH4 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(4,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: Hydrate ' //
     &        'O2 or CH4 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(5,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: Hydrate ' //
     &        'O2 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate CH4 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(4,N)+IC_OPT(5,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: Hydrate ' //
     &        'CH4 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous saturation unspecified  ---
!
          IF( IC_OPT(6,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: ' //
     &        'Aqueous Saturation Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 mole fraction of formers specified  ---
!
          IF( IC_OPT(8,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: ' //
     &        'Mobile O2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile CH4 mole fraction of formers specified  ---
!
          IF( IC_OPT(9,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: ' //
     &        'Mobile CH4 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile N2 mole fraction of formers specified  ---
!
          IF( IC_OPT(10,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: ' //
     &        'Mobile N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Relative saturation of dissolved formers specified  ---
!
          IF( IC_OPT(11,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: ' //
     &        'Relative Saturation of Dissolved Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Temperature not specified  ---
!
          IF( IC_OPT(13,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 5: ' //
     &        'Temperature Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!!
!!---      Pressure specified  ---
!!
!          IF( IC_OPT(14,N).EQ.0 ) THEN
!            INDX = 1
!            CHMSG = 'Specification Note: ' // 
!     &        'Initial Condition State 5: ' //
!     &        'Pressure Unspecified.'
!            CALL WRMSGS( INDX )
!          ENDIF
!
!---      Flash calculation for initial condition #5 ---
!
          PVHA(2,N) = XHA(2,N)
          PVHO(2,N) = XHO(2,N)
          PVHN(2,N) = XHN(2,N)
          IF( (IC_OPT(3,N)+IC_OPT(4,N)).EQ.2 ) THEN
            PVHN(2,N) = MAX( 1.D+0-PVHA(2,N)-PVHO(2,N),0.D+0 )
            IF( PVHN(2,N).LT.EPSL ) PVHN(2,N) = 0.D+0
          ELSEIF( (IC_OPT(4,N)+IC_OPT(5,N)).EQ.2 ) THEN
            PVHA(2,N) = MAX( 1.D+0-PVHO(2,N)-PVHN(2,N),0.D+0 )
            IF( PVHA(2,N).LT.EPSL ) PVHA(2,N) = 0.D+0
          ELSEIF( (IC_OPT(5,N)+IC_OPT(3,N)).EQ.2 ) THEN
            PVHO(2,N) = MAX( 1.D+0-PVHN(2,N)-PVHA(2,N),0.D+0 )
            IF( PVHO(2,N).LT.EPSL ) PVHO(2,N) = 0.D+0
          ENDIF
          CALL FLH_IC5( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),PVN(2,N),
     &      PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),RHOH(2,N),
     &      SH(2,N),SI(2,N),SG(2,N),SL(2,N),SN(2,N),
     &      T(2,N),TCR(N),TCT(N),TEQHX,TMHA(2,N),TMHN(2,N),TMHO(2,N),
     &      XHA(2,N),XHN(2,N),XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),
     &      YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),YMHGN(2,N),
     &      YMHGO(2,N),ZMCA(2,N),ZMCN(2,N),ZMCO(2,N),
     &      IZN,N,NP_HT )
          POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      Zero hydrate saturation  ---
!
          IF( SH(2,N).LT.EPSL ) THEN
            IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
              NPHAZ(2,N) = 10 + (100*NP_HT)
            ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              NPHAZ(2,N) = 21 + (100*NP_HT)
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              NPHAZ(2,N) = 22 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 23 + (100*NP_HT)
            ENDIF
!
!---      Non-zero hydrate saturation  ---
!
          ELSE
            IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 41 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 42 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 43 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 31 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 32 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 33 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 34 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 35 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 36 + (100*NP_HT)
              ENDIF
            ELSE
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 37 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 38 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 39 + (100*NP_HT)
              ENDIF
            ENDIF
          ENDIF
!
!---    Initial condition state #6
!
!       SHX = 0.0
!       SGX + SNX > 0.0
!       SLX + SIX < 1.0
!       Hydrate-mobile phase non-equilibrium
!
!       Declared variables:
!
!       T - temperature, C (T/T)
!       PSO - pressure, Pa (PSO/PSO)
!       SL - aqueous saturation (SL/SL)
!       2/3 ZMCA, ZMCO, ZMCN
!       ZMCA - total nonaqueous O2 mole fraction (XMGA/PVA)
!       ZMCO - total nonaqueous CH4 mole fraction (XMGO/PVO)
!       ZMCN - total nonaqueous N2 mole fraction (XMGN/PVN)
!
        ELSEIF( IC_OPT(1,N).EQ.6 ) THEN
!
!---      Hydrate saturation specified  ---
!
          IF( IC_OPT(2,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: ' //
     &        'Hydrate Saturation Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 mole fraction of formers specified  ---
!
          IF( IC_OPT(3,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: ' //
     &        'Hydrate O2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate CH4 mole fraction of formers specified  ---
!
          IF( IC_OPT(4,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: ' //
     &        'Hydrate CH4 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate N2 mole fraction of formers specified  ---
!
          IF( IC_OPT(5,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: ' //
     &        'Hydrate N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous saturation not specified  ---
!
          IF( IC_OPT(6,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: ' //
     &        'Aqueous Saturation Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 or CH4 or N2 mole fraction of formers 
!         specified  ---
!
          IF( IC_OPT(8,N).EQ.1 .OR. IC_OPT(9,N).EQ.1 .OR.
     &      IC_OPT(10,N).EQ.1 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: Mobile ' //
     &        'O2 or CH4 or N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous O2 + CH4 + N2 mole fraction of formers 
!         specified  ---
!
          IF( IC_OPT(8,N).EQ.2 .AND. IC_OPT(9,N).EQ.2 .AND.
     &      IC_OPT(10,N).EQ.2 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: Nonaqueous ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous O2 or CH4 mole fraction of formers 
!         unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: Nonaqueous ' //
     &        'O2 or CH4 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous O2 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: Nonaqueous ' //
     &        'O2 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous CH4 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: Nonaqueous ' //
     &        'CH4 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Relative saturation of dissolved formers specified  ---
!
          IF( IC_OPT(11,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: ' //
     &        'Relative Saturation of Dissolved Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Temperature not specified  ---
!
          IF( IC_OPT(13,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: ' //
     &        'Temperature Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Pressure not specified  ---
!
          IF( IC_OPT(14,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 6: ' //
     &        'Pressure Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Flash calculation for initial condition #6 ---
!
          PVA(2,N) = XMGA(2,N)
          PVO(2,N) = XMGO(2,N)
          PVN(2,N) = XMGN(2,N)
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.2 ) THEN
            PVN(2,N) = MAX( 1.D+0-PVA(2,N)-PVO(2,N),0.D+0 )
            IF( PVN(2,N).LT.EPSL ) PVN(2,N) = 0.D+0
          ELSEIF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.2 ) THEN
            PVA(2,N) = MAX( 1.D+0-PVO(2,N)-PVN(2,N),0.D+0 )
            IF( PVA(2,N).LT.EPSL ) PVA(2,N) = 0.D+0
          ELSEIF( (IC_OPT(10,N)+IC_OPT(8,N)).EQ.2 ) THEN
            PVO(2,N) = MAX( 1.D+0-PVN(2,N)-PVA(2,N),0.D+0 )
            IF( PVO(2,N).LT.EPSL ) PVO(2,N) = 0.D+0
          ENDIF
          CALL FLH_IC6( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),PVN(2,N),
     &      PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),RHOH(2,N),
     &      SH(2,N),SI(2,N),SG(2,N),SL(2,N),SN(2,N),
     &      T(2,N),TCR(N),TCT(N),TEQHX,
     &      XHA(2,N),XHN(2,N),XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),
     &      YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),YMHGN(2,N),
     &      YMHGO(2,N),ZMCA(2,N),ZMCN(2,N),ZMCO(2,N),
     &      IZN,N,NP_HT )
          POSM(2,N) = SL(2,N) + SI(2,N)
          IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
            NPHAZ(2,N) = 10 + (100*NP_HT)
          ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
            NPHAZ(2,N) = 21 + (100*NP_HT)
          ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
            NPHAZ(2,N) = 22 + (100*NP_HT)
          ELSE
            NPHAZ(2,N) = 23 + (100*NP_HT)
          ENDIF
!
!---    Initial condition state #7
!
!       SHX > 0.0
!       SGX + SNX > 0.0
!       SLX + SIX < 1.0
!       Hydrate-mobile phase non-equilibrium
!
!       Declared variables:
!
!       T - temperature, C (T/T)
!       PSO - pressure, Pa (PSO/PSO)
!       SL - aqueous saturation (SL/SL)
!       2/3 ZMCA, ZMCO, ZMCN
!       ZMCA - total nonaqueous O2 mole fraction (XMGA/PVA)
!       ZMCO - total nonaqueous CH4 mole fraction (XMGO/PVO)
!       ZMCN - total nonaqueous N2 mole fraction (XMGN/PVN)
!       SH - hydrate saturation (SH/SH)
!       2/3 YMHA, YMHO, YMHN
!       YMHA - hydrate O2 mole fraction of formers (XMGA/PVA)
!       YMHO - hydrate CH4 mole fraction of formers (XMGO/PVO)
!       YMHN - hydrate N2 mole fraction of formers (XMGN/PVN)
!
        ELSEIF( IC_OPT(1,N).EQ.7 ) THEN
!
!---      Hydrate saturation not specified  ---
!
          IF( IC_OPT(2,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: ' //
     &        'Hydrate Saturation Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 + CH4 + N2 mole fraction of formers specified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(4,N)+IC_OPT(5,N)).EQ.3 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: Hydrate ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 or CH4 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(4,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: Hydrate ' //
     &        'O2 or CH4 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate O2 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(3,N)+IC_OPT(5,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: Hydrate ' //
     &        'O2 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Hydrate CH4 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(4,N)+IC_OPT(5,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: Hydrate ' //
     &        'CH4 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Aqueous saturation unspecified  ---
!
          IF( IC_OPT(6,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: ' //
     &        'Aqueous Saturation Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Mobile O2 or CH4 or N2 mole fraction of formers 
!         specified  ---
!
          IF( IC_OPT(8,N).EQ.1 .OR. IC_OPT(9,N).EQ.1 .OR.
     &      IC_OPT(10,N).EQ.1 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: Mobile ' //
     &        'O2 or CH4 or N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous O2 + CH4 + N2 mole fraction of formers 
!         specified  ---
!
          IF( IC_OPT(8,N).EQ.2 .AND. IC_OPT(9,N).EQ.2 .AND.
     &      IC_OPT(10,N).EQ.2 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: Nonaqueous ' //
     &        'O2 + CH4 + N2 Mole Fraction of Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous O2 or CH4 mole fraction of formers 
!         unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: Nonaqueous ' //
     &        'O2 or CH4 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous O2 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(8,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: Nonaqueous ' //
     &        'O2 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Nonaqueous CH4 or N2 mole fraction of formers unspecified  ---
!
          IF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: Mobile ' //
     &        'CH4 or N2 Mole Fraction of Formers Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Relative saturation of dissolved formers specified  ---
!
          IF( IC_OPT(11,N).NE.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: ' //
     &        'Relative Saturation of Dissolved Formers Specified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Temperature not specified  ---
!
          IF( IC_OPT(13,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: ' //
     &        'Temperature Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Pressure not specified  ---
!
          IF( IC_OPT(14,N).EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Specification Error: ' // 
     &        'Initial Condition State 7: ' //
     &        'Pressure Unspecified.'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Flash calculation for initial condition #7 ---
!
          PVA(2,N) = XMGA(2,N)
          PVO(2,N) = XMGO(2,N)
          PVN(2,N) = XMGN(2,N)
          IF( (IC_OPT(8,N)+IC_OPT(9,N)).EQ.2 ) THEN
            PVN(2,N) = MAX( 1.D+0-PVA(2,N)-PVO(2,N),0.D+0 )
            IF( PVN(2,N).LT.EPSL ) PVN(2,N) = 0.D+0
          ELSEIF( (IC_OPT(9,N)+IC_OPT(10,N)).EQ.2 ) THEN
            PVA(2,N) = MAX( 1.D+0-PVO(2,N)-PVN(2,N),0.D+0 )
            IF( PVA(2,N).LT.EPSL ) PVA(2,N) = 0.D+0
          ELSEIF( (IC_OPT(10,N)+IC_OPT(8,N)).EQ.2 ) THEN
            PVO(2,N) = MAX( 1.D+0-PVN(2,N)-PVA(2,N),0.D+0 )
            IF( PVO(2,N).LT.EPSL ) PVO(2,N) = 0.D+0
          ENDIF
          PVHA(2,N) = XHA(2,N)
          PVHO(2,N) = XHO(2,N)
          PVHN(2,N) = XHN(2,N)
          IF( (IC_OPT(3,N)+IC_OPT(4,N)).EQ.2 ) THEN
            PVHN(2,N) = MAX( 1.D+0-PVHA(2,N)-PVHO(2,N),0.D+0 )
            IF( PVHN(2,N).LT.EPSL ) PVHN(2,N) = 0.D+0
          ELSEIF( (IC_OPT(4,N)+IC_OPT(5,N)).EQ.2 ) THEN
            PVHA(2,N) = MAX( 1.D+0-PVHO(2,N)-PVHN(2,N),0.D+0 )
            IF( PVHA(2,N).LT.EPSL ) PVHA(2,N) = 0.D+0
          ELSEIF( (IC_OPT(5,N)+IC_OPT(3,N)).EQ.2 ) THEN
            PVHO(2,N) = MAX( 1.D+0-PVHN(2,N)-PVHA(2,N),0.D+0 )
            IF( PVHO(2,N).LT.EPSL ) PVHO(2,N) = 0.D+0
          ENDIF
          CALL FLH_IC7( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),PVN(2,N),
     &      PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),RHOH(2,N),
     &      SH(2,N),SI(2,N),SG(2,N),SL(2,N),SN(2,N),
     &      T(2,N),TCR(N),TCT(N),TEQHX,TMHA(2,N),TMHN(2,N),TMHO(2,N),
     &      XHA(2,N),XHN(2,N),XHO(2,N),XHW(2,N),XLS(2,N),YLS(2,N),
     &      YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),YMHGN(2,N),
     &      YMHGO(2,N),ZMCA(2,N),ZMCN(2,N),ZMCO(2,N),
     &      IZN,N,NP_HT )
          POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      Zero hydrate saturation  ---
!
          IF( SH(2,N).LT.EPSL ) THEN
            IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
              NPHAZ(2,N) = 10 + (100*NP_HT)
            ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              NPHAZ(2,N) = 21 + (100*NP_HT)
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              NPHAZ(2,N) = 22 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 23 + (100*NP_HT)
            ENDIF
!
!---      Non-zero hydrate saturation  ---
!
          ELSE
            IF( SG(2,N)+SN(2,N).LT.EPSL ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 41 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 42 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 43 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 31 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 32 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 33 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 34 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 35 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 36 + (100*NP_HT)
              ENDIF
            ELSE
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 37 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 38 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 39 + (100*NP_HT)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
  300 CONTINUE
!
!---    Establish reference pressure for soil compressibility  ---
!
      DO 320 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 320
        TCMP(N) = T(2,N)
        IZN = IZ(N)
        IF( CMP(3,IZN).GT.PATM ) THEN
          PCMP(N) = CMP(3,IZN)
        ELSEIF( ISLC(61).EQ.0 ) THEN
          PCMP(N) = PSO(2,N) + PATM
        ENDIF
  320 CONTINUE
  400 CONTINUE
!
!---  Initializations for normal and restart simulations  ---
!
      DO 500 N = 1,NFLD
!
!---    Aqueous and ice saturation  ---
!
        POSM(2,N) = SL(2,N) + SI(2,N)
!
!---    Porous-media porosity  ---
!
        PPX = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
        IF( ISLC(50).NE.0 ) P_GM(2,N) = PPX
        CALL PORSTY_HYDT_KE( N,PPX,PCMP(N),PORD(2,N),PORT(2,N),
     &    SH(2,N) )
        PORD(2,N) = MAX( PORD(2,N),EPSL )
        PORT(2,N) = MAX( PORT(2,N),PORD(2,N) )
!
!---    Initialize old-time-step values  ---
!
        T(1,N) = T(2,N)
        PL(1,N) = PL(2,N)
        PG(1,N) = PG(2,N)
        PN(1,N) = PN(2,N)
        PSO(1,N) = PSO(2,N)
        ZMCA(1,N) = ZMCA(2,N)
        ZMCO(1,N) = ZMCO(2,N)
        ZMCN(1,N) = ZMCN(2,N)
        SG(1,N) = SG(2,N)
        SL(1,N) = SL(2,N)
        SH(1,N) = SH(2,N)
        SI(1,N) = SI(2,N)
        SN(1,N) = SN(2,N)
        POSM(1,N) = POSM(2,N)
        YLS(1,N) = YLS(2,N)
        YMGA(1,N) = YMGA(2,N)
        YMGO(1,N) = YMGO(2,N)
        YMGN(1,N) = YMGN(2,N)
        YMHGA(1,N) = YMHGA(2,N)
        YMHGO(1,N) = YMHGO(2,N)
        YMHGN(1,N) = YMHGN(2,N)
        PORD(1,N) = PORD(2,N)
        TMHA(1,N) = TMHA(2,N)
        TMHO(1,N) = TMHO(2,N)
        TMHN(1,N) = TMHN(2,N)
        NPHAZ(1,N) = NPHAZ(2,N)
!
!---    Assign dissolved-salt mass fractions,
!       isobrine option  ---
!
        IF( ISLC(32).EQ.1 ) THEN
          DO 410 M = 3,ISVC+2
            XLS(M,N) = XLS(2,N)
            YLS(M,N) = YLS(2,N)
  410     CONTINUE
        ENDIF
  500 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Compute initial solute concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
      SUB_LOG(ISUB_LOG) = '/CISC_HYD_O2_CH4_N2'
      DO 140 NSL = 1,NSOLU
        DO 110 N = 1,NFLD
          N_DB = N
          IF( IXP(N).EQ.0 ) GOTO 110
          IZN = IZ(N)
          IF( IPCL(NSL).EQ.2 ) THEN
            XVS = SL(2,N)*RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
          ELSE
            XVS = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
          ENDIF
          XVL = SL(2,N)*PORD(2,N)
          XVG = SG(2,N)*PORD(2,N)
!
!---      Constant gas-aqueous partition coefficient  ---
!
          IF( IPCGL(NSL).EQ.0 ) THEN
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
!---  Phase-volumetric concentration ratios  ---
!
          IF( ICT(N,NSL).EQ.2 ) THEN
            C(N,NSL) = C(N,NSL)*(XVS + XVL + XVG*PCGLX)
          ELSEIF( ICT(N,NSL).EQ.3 ) THEN
            C(N,NSL) = C(N,NSL)*((XVS + XVL)/PCGLX + XVG)
          ENDIF
  110   CONTINUE
!
!---  Assign boundary solute concentrations for initial condition
!     type boundary conditions  ---
!
        DO 130 NB = 1,NBC
          IF( IBCT(NSL+LUK+LPH,NB).EQ.12 ) THEN
            N = IBCN(NB)
            CBO(NB,NSL) = C(N,NSL)
          ENDIF
  130   CONTINUE
  140 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CISC_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFINA_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Reset temporarily inactive nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      SUB_LOG(ISUB_LOG) = '/DFINA_HYD_O2_CH4_N2'
      IF( ICNV.EQ.3 ) THEN
        DO 100 N = 1,NFLD
          IXP(N) = ABS(IXP(N))
  100   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFINA_HYD_O2_CH4_N2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Compute primary variable increments.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE FDVS
      USE FDVP
      USE FDVN
      USE FDVH
      USE FDVG
      USE FDVA
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
      SUB_LOG(ISUB_LOG) = '/INCRM_HYD_O2_CH4_N2'
      TOLX = 1.D-4
!
!---  Phase options, compute phase condition   ---
!
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 100
        IZN = IZ(N)
        N_DB = N
        IF( NITER.LE.3 ) THEN
          ZMIH_HT(1) = YMHGA(2,N)
          ZMIH_HT(2) = YMHGO(2,N)
          ZMIH_HT(3) = YMHGN(2,N)
          SNR(1,N) = YMHGA(2,N)
          SNR(2,N) = YMHGO(2,N)
          SNR(3,N) = YMHGN(2,N)
        ELSE
          ZMIH_HT(1) = SNR(1,N)
          ZMIH_HT(2) = SNR(2,N)
          ZMIH_HT(3) = SNR(3,N)
        ENDIF
!
!---    Absolute temperature and pressure  ---
!
        TKX = T(2,N) + TABS
        PX = PSO(2,N) + PATM
!
!---    Core phase condition  ---
!
        NPHZX = MOD(NPHAZ(2,N),100)
!
!---    Phase condition #10 (10 series)  ---
!
!       SH = 0, SG + SN = 0, SL+SI = 1
!
!       energy - temperature, C - T
!       water mass - pressure, Pa (gauge) - PSO
!       mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!       mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!       mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!       hydrate O2 mass - hydrate O2 mass, kg - TMHA
!       hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!       hydrate N2 mass - hydrate N2 mass, kg - TMHN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        IF( NPHZX.EQ.10 ) THEN
          CALL FLH_PC10( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),PVN(2,N),
     &      PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),RHOG(2,N),RHOH(2,N),
     &      RHON(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),
     &      SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,SN(2,N),T(2,N),TCR(N),
     &      TCT(N),TEQHX,XGA(2,N),
     &      XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),XHN(2,N),XHO(2,N),
     &      XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),XNW(2,N),XMLA(2,N),
     &      XMLN(2,N),XMLO(2,N),XLS(2,N),YLS(2,N),YMGA(2,N),YMGN(2,N),
     &      YMGO(2,N),YMHGA(2,N),YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,
     &      YMNOX,ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
          POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      Absolute total vapor pressure, mobile component 
!         vapor pressure, hydrate component vapor pressure,
!         gas pressure, and nonaqueous-liquid pressure  ---
!
          PVX = PVA(2,N) + PVO(2,N) + PVN(2,N) + PVW(2,N)
          PVMX = PVA(2,N) + PVO(2,N) + PVN(2,N)
          PVHX = PVHA(2,N) + PVHO(2,N) + PVHN(2,N)
          PGX = PG(2,N) + PATM
          PNX = PN(2,N) + PATM
!
!---      Transition to hydrate conditions with nonaqueous phases
!         (10 series -> 30 series)  ---
!
          IF( PVHX.GE.PH(2,N) .AND. PVX.GT.(PX*(1.D+0+TOLX)) .AND. 
     &      INHYD(IZN).EQ.0 ) THEN
            IF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 31 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 32 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 33 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 34 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 35 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 36 + (100*NP_HT)
              ENDIF
            ELSE
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 37 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 38 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 39 + (100*NP_HT)
              ENDIF
            ENDIF
            SH(2,N) = 1.D-7
            POSM(2,N) = POSM(2,N) - 1.D-4*POSM(2,N)
            SLIX = POSM(2,N)
            CALL FLH_PC30( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RHON(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),
     &        SLIX,SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,TMHA(2,N),
     &        TMHN(2,N),TMHO(2,N),XGA(2,N),XGN(2,N),XGO(2,N),
     &        XGW(2,N),XHA(2,N),XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),
     &        XNN(2,N),XNO(2,N),XNW(2,N),XMLA(2,N),XMLN(2,N),
     &        XMLO(2,N),XLS(2,N),YLS(2,N),YMGA(2,N),YMGN(2,N),
     &        YMGO(2,N),YMHGA(2,N),YMHGN(2,N),YMHGO(2,N),YMNAX,
     &        YMNNX,YMNOX,ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),
     &        IZN,N,NP_HT )
!
!---      Transition to hydrate conditions with no nonaqueous phases
!         (10 series -> 40 series)  ---
!
          ELSEIF( PVHX.GE.PH(2,N) .AND. INHYD(IZN).EQ.0 ) THEN
            IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
              NPHAZ(2,N) = 41 + (100*NP_HT)
            ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
              NPHAZ(2,N) = 42 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 43 + (100*NP_HT)
            ENDIF
            SH(2,N) = 1.D-7
            CALL FLH_PC40( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RHON(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &        SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,TMHA(2,N),TMHN(2,N),
     &        TMHO(2,N),XGA(2,N),XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),
     &        XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),
     &        XNW(2,N),XMLA(2,N),XMLN(2,N),XMLO(2,N),XLS(2,N),
     &        YLS(2,N),YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),
     &        YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,YMNOX,
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
            POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      Transition to no hydrate conditions with nonaqueous phases
!         (10 series -> 20 series)  ---
!
          ELSEIF( PVX.GT.(PX*(1.D+0+TOLX)) ) THEN
            IF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              NPHAZ(2,N) = 21 + (100*NP_HT)
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              NPHAZ(2,N) = 22 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 23 + (100*NP_HT)
            ENDIF
            POSM(2,N) = POSM(2,N) - 1.D-4*POSM(2,N)
            SLIX = POSM(2,N)
            CALL FLH_PC20( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &        SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,
     &        XGA(2,N),XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),
     &        XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),
     &        XNW(2,N),XMLA(2,N),XMLN(2,N),XMLO(2,N),XLS(2,N),
     &        YLS(2,N),YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),
     &        YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,YMNOX,
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
          ENDIF
!
!---    Phase conditions 21 through 23 (20 series)  ---
!
!       SH = 0, SG + SN > 0, SL + SI < 1
!       ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #21
!       ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #22
!       ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #23
!
!       energy - temperature, C - T
!       water mass - aqueous + ice saturation - POSM
!       maximum mobile mass component - pressure, Pa (gauge) - PSO
!       mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!       mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!       mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!       hydrate O2 mass - hydrate O2 mass, kg - TMHA
!       hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!       hydrate N2 mass - hydrate N2 mass, kg - TMHN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.21 .AND. NPHZX.LE.23  ) THEN
          IF( NPHZX.EQ.21 ) THEN
            ZMCA(2,N) = MIN( 1.D+0-ZMCO(2,N)-ZMCN(2,N),1.D+0 )
            IF( ZMCA(2,N).LT.EPSL ) ZMCA(2,N) = 0.D+0
          ELSEIF( NPHZX.EQ.22 ) THEN
            ZMCO(2,N) = MIN( 1.D+0-ZMCA(2,N)-ZMCN(2,N),1.D+0 )
            IF( ZMCO(2,N).LT.EPSL ) ZMCO(2,N) = 0.D+0
          ELSE
            ZMCN(2,N) = MIN( 1.D+0-ZMCA(2,N)-ZMCO(2,N),1.D+0 )
            IF( ZMCN(2,N).LT.EPSL ) ZMCN(2,N) = 0.D+0
          ENDIF
          SLIX = POSM(2,N)
          CALL FLH_PC20( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),PVN(2,N),
     &      PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),RHOG(2,N),RHOH(2,N),
     &      RKG(2,N),RKL(1,2,N),RKN(2,N),SH(2,N),
     &      SI(2,N),SG(2,N),SL(2,N),SLIX,SN(2,N),T(2,N),TCR(N),
     &      TCT(N),TEQHX,XGA(2,N),
     &      XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),XHN(2,N),XHO(2,N),
     &      XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),XNW(2,N),XMLA(2,N),
     &      XMLN(2,N),XMLO(2,N),XLS(2,N),YLS(2,N),YMGA(2,N),YMGN(2,N),
     &      YMGO(2,N),YMHGA(2,N),YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,
     &      YMNOX,ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT ) 
!
!---      Absolute total vapor pressure, mobile component 
!         vapor pressure, hydrate component vapor pressure,
!         gas pressure, and nonaqueous-liquid pressure  ---
!
          PVX = PVA(2,N) + PVO(2,N) + PVN(2,N) + PVW(2,N)
          PVMX = PVA(2,N) + PVO(2,N) + PVN(2,N)
          PVHX = PVHA(2,N) + PVHO(2,N) + PVHN(2,N)
          PGX = PG(2,N) + PATM
          PNX = PN(2,N) + PATM
!
!---      Transition to hydrate conditions with nonaqueous phases
!         (20 series -> 30 series)  ---
!
          IF( PVHX.GE.PH(2,N) .AND. (SG(2,N)+SN(2,N)).GT.EPSL .AND. 
     &      INHYD(IZN).EQ.0 ) THEN
            IF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 31 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 32 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 33 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 34 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 35 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 36 + (100*NP_HT)
              ENDIF
            ELSE
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 37 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 38 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 39 + (100*NP_HT)
              ENDIF
            ENDIF
            SH(2,N) = 1.D-7
            SLIX = POSM(2,N)
            CALL FLH_PC30( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RHON(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),
     &        SLIX,SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,TMHA(2,N),
     &        TMHN(2,N),TMHO(2,N),XGA(2,N),XGN(2,N),XGO(2,N),
     &        XGW(2,N),XHA(2,N),XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),
     &        XNN(2,N),XNO(2,N),XNW(2,N),XMLA(2,N),XMLN(2,N),
     &        XMLO(2,N),XLS(2,N),YLS(2,N),YMGA(2,N),YMGN(2,N),
     &        YMGO(2,N),YMHGA(2,N),YMHGN(2,N),YMHGO(2,N),YMNAX,
     &        YMNNX,YMNOX,ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),
     &        IZN,N,NP_HT )
!
!---      Transition to hydrate conditions with no nonaqueous phases
!         (20 series -> 40 series)  ---
!
          ELSEIF( PVHX.GE.PH(2,N) .AND. (SG(2,N)+SN(2,N)).LE.EPSL .AND. 
     &      INHYD(IZN).EQ.0 ) THEN
            IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
              NPHAZ(2,N) = 41 + (100*NP_HT)
            ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
              NPHAZ(2,N) = 42 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 43 + (100*NP_HT)
            ENDIF
            SH(2,N) = 1.D-7
            CALL FLH_PC40( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RHON(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &        SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,TMHA(2,N),TMHN(2,N),
     &        TMHO(2,N),XGA(2,N),XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),
     &        XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),
     &        XNW(2,N),XMLA(2,N),XMLN(2,N),XMLO(2,N),XLS(2,N),
     &        YLS(2,N),YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),
     &        YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,YMNOX,
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
            POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      Transition to no hydrate conditions with no nonaqueous phases
!         (20 series -> 10 series)  ---
!
          ELSEIF( (SG(2,N)+SN(2,N)).LE.EPSL ) THEN
            NPHAZ(2,N) = 10 + (100*NP_HT)
            CALL FLH_PC10( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RHON(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &        SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,
     &        XGA(2,N),XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),
     &        XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),
     &        XNW(2,N),XMLA(2,N),XMLN(2,N),XMLO(2,N),XLS(2,N),
     &        YLS(2,N),YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),
     &        YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,YMNOX,
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
            POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      No transition from no hydrate with nonaqueous phases  ---
!
          ELSE
            IF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              NPHAZ(2,N) = 21 + (100*NP_HT)
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              NPHAZ(2,N) = 22 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 23 + (100*NP_HT)
            ENDIF
          ENDIF
!
!---    Phase conditions #31 through #39 (30 series)  ---
!
!       SH > 0, SG + SN > 0, SL + SI < 1
!       ZMCA = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #31
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #32
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #33
!       ZMCO = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #34
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #35
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #36
!       ZMCN = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #37
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #38
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #39
!
!       energy - temperature, C - T
!       water mass - aqueous + ice saturation, POSM
!       maximum mobile mass component - pressure, Pa (gauge) - PSO
!       mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!       mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!       mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!       maximum hydrate mass component - hydrate saturation - SH
!       hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!       hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!       hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.31 .AND. NPHZX.LE.39 ) THEN
          IF( NPHZX.GE.31 .AND. NPHZX.LE.33 ) THEN
            ZMCA(2,N) = MIN( 1.D+0-ZMCO(2,N)-ZMCN(2,N),1.D+0 )
            IF( ZMCA(2,N).LT.EPSL ) ZMCA(2,N) = 0.D+0
          ELSEIF( NPHZX.GE.34 .AND. NPHZX.LE.36 ) THEN
            ZMCO(2,N) = MIN( 1.D+0-ZMCA(2,N)-ZMCN(2,N),1.D+0 )
            IF( ZMCO(2,N).LT.EPSL ) ZMCO(2,N) = 0.D+0
          ELSE
            ZMCN(2,N) = MIN( 1.D+0-ZMCA(2,N)-ZMCO(2,N),1.D+0 )
            IF( ZMCN(2,N).LT.EPSL ) ZMCN(2,N) = 0.D+0
          ENDIF
          IF( NPHZX.EQ.31 .OR. NPHZX.EQ.34 .OR. NPHZX.EQ.37 ) THEN
            YMHGA(2,N) = MIN( 1.D+0-YMHGO(2,N)-YMHGN(2,N),1.D+0 )
            IF( YMHGA(2,N).LT.EPSL ) YMHGA(2,N) = 0.D+0
          ELSEIF( NPHZX.EQ.32 .OR. NPHZX.EQ.35 .OR. NPHZX.EQ.38 )THEN
            YMHGO(2,N) = MIN( 1.D+0-YMHGA(2,N)-YMHGN(2,N),1.D+0 )
            IF( YMHGO(2,N).LT.EPSL ) YMHGO(2,N) = 0.D+0
          ELSE
            YMHGN(2,N) = MIN( 1.D+0-YMHGA(2,N)-YMHGO(2,N),1.D+0 )
            IF( YMHGN(2,N).LT.EPSL ) YMHGN(2,N) = 0.D+0
          ENDIF
          SLIX = POSM(2,N)
          CALL FLH_PC30( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),PVN(2,N),
     &      PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),RHOG(2,N),RHOH(2,N),
     &      RHON(2,N),RKG(2,N),RKL(1,2,N),RKN(2,N),
     &      SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,SN(2,N),T(2,N),TCR(N),
     &      TCT(N),TEQHX,TMHA(2,N),TMHN(2,N),TMHO(2,N),XGA(2,N),
     &      XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),XHN(2,N),XHO(2,N),
     &      XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),XNW(2,N),XMLA(2,N),
     &      XMLN(2,N),XMLO(2,N),XLS(2,N),YLS(2,N),YMGA(2,N),YMGN(2,N),
     &      YMGO(2,N),YMHGA(2,N),YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,
     &      YMNOX,ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT ) 
!
!---      Absolute total vapor pressure, mobile component 
!         vapor pressure, hydrate component vapor pressure,
!         gas pressure, and nonaqueous-liquid pressure  ---
!
          PVX = PVA(2,N) + PVO(2,N) + PVN(2,N) + PVW(2,N)
          PVMX = PVA(2,N) + PVO(2,N) + PVN(2,N)
          PVHX = PVHA(2,N) + PVHO(2,N) + PVHN(2,N)
          PGX = PG(2,N) + PATM
          PNX = PN(2,N) + PATM
!
!---      Transition to no-hydrate conditions + nonaqueous phases 
!         (30 series -> 20 series)  ---
!
          IF( SH(2,N).LE.EPSL .AND. PVMX.LT.PVHX  ) THEN
            IF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              NPHAZ(2,N) = 21 + (100*NP_HT)
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              NPHAZ(2,N) = 22 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 23 + (100*NP_HT)
            ENDIF
            SLIX = POSM(2,N)
            CALL FLH_PC20( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &        SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,
     &        XGA(2,N),XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),
     &        XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),
     &        XNW(2,N),XMLA(2,N),XMLN(2,N),XMLO(2,N),XLS(2,N),
     &        YLS(2,N),YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),
     &        YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,YMNOX,
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
!
!---      Transition to no-hydrate conditions + no nonaqueous phases 
!         (30 series -> 10 series)  ---
!
          ELSEIF( SH(2,N).LE.EPSL .AND. PVMX.LT.PVHX .AND. 
     &      (SG(2,N)+SN(2,N)).LT.EPSL  ) THEN
            NPHAZ(2,N) = 10 + (100*NP_HT)
            CALL FLH_PC10( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RHON(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &        SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,
     &        XGA(2,N),XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),
     &        XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),
     &        XNW(2,N),XMLA(2,N),XMLN(2,N),XMLO(2,N),XLS(2,N),
     &        YLS(2,N),YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),
     &        YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,YMNOX,
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
            POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      Transition to hydrate conditions + no nonaqueous phases 
!         (30 series -> 40 series)  ---
!
          ELSEIF( (SG(2,N)+SN(2,N)).LT.EPSL  ) THEN
            IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
              NPHAZ(2,N) = 41 + (100*NP_HT)
            ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
              NPHAZ(2,N) = 42 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 43 + (100*NP_HT)
            ENDIF
            CALL FLH_PC40( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RHON(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &        SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,TMHA(2,N),TMHN(2,N),
     &        TMHO(2,N),XGA(2,N),XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),
     &        XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),
     &        XNW(2,N),XMLA(2,N),XMLN(2,N),XMLO(2,N),XLS(2,N),
     &        YLS(2,N),YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),
     &        YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,YMNOX,
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
            POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      No transition from hydrate conditions + nonaqueous phases  ---
!
          ELSE
            IF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 31 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 32 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 33 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 34 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 35 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 36 + (100*NP_HT)
              ENDIF
            ELSE
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 37 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 38 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 39 + (100*NP_HT)
              ENDIF
            ENDIF
          ENDIF
!
!---    Phase conditions #41 through #43 (40 series)  ---
!
!       SH > 0, SG + SN = 0, SL+SI < 1
!       YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #41
!       YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #42
!       YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #43
!
!       energy - temperature, C - T
!       water mass - pressure, Pa (gauge) - PSO
!       mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!       mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!       mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!       maximum hydrate mass component - hydrate saturation - SH
!       hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!       hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!       hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.41 .AND. NPHZX.LE.43 ) THEN
          IF( NPHZX.EQ.41 ) THEN
            YMHGA(2,N) = MIN( 1.D+0-YMHGO(2,N)-YMHGN(2,N),1.D+0 )
            IF( YMHGA(2,N).LT.EPSL ) YMHGA(2,N) = 0.D+0
          ELSEIF( NPHZX.EQ.42 )THEN
            YMHGO(2,N) = MIN( 1.D+0-YMHGA(2,N)-YMHGN(2,N),1.D+0 )
            IF( YMHGO(2,N).LT.EPSL ) YMHGO(2,N) = 0.D+0
          ELSE
            YMHGN(2,N) = MIN( 1.D+0-YMHGA(2,N)-YMHGO(2,N),1.D+0 )
            IF( YMHGN(2,N).LT.EPSL ) YMHGN(2,N) = 0.D+0
          ENDIF
          CALL FLH_PC40( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &      PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &      PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &      RHOG(2,N),RHOH(2,N),RHON(2,N),RKG(2,N),
     &      RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &      SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,TMHA(2,N),TMHN(2,N),
     &      TMHO(2,N),XGA(2,N),XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),
     &      XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),
     &      XNW(2,N),XMLA(2,N),XMLN(2,N),XMLO(2,N),XLS(2,N),
     &      YLS(2,N),YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),
     &      YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,YMNOX,
     &      ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
          POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      Absolute total vapor pressure, mobile component 
!         vapor pressure, hydrate component vapor pressure,
!         gas pressure, and nonaqueous-liquid pressure  ---
!
          PVX = PVA(2,N) + PVO(2,N) + PVN(2,N) + PVW(2,N)
          PVMX = PVA(2,N) + PVO(2,N) + PVN(2,N)
          PVHX = PVHA(2,N) + PVHO(2,N) + PVHN(2,N)
          PGX = PG(2,N) + PATM
          PNX = PN(2,N) + PATM
!
!---      Transition to hydrate conditions + nonaqueous phases 
!         (40 series -> 30 series)  ---
!
          IF( PVHX.GE.PH(2,N) .AND. PVX.GT.(PX*(1.D+0+TOLX)) ) THEN
            IF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 31 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 32 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 33 + (100*NP_HT)
              ENDIF
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 34 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 35 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 36 + (100*NP_HT)
              ENDIF
            ELSE
              IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
                NPHAZ(2,N) = 37 + (100*NP_HT)
              ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
                NPHAZ(2,N) = 38 + (100*NP_HT)
              ELSE
                NPHAZ(2,N) = 39 + (100*NP_HT)
              ENDIF
            ENDIF
            POSM(2,N) = POSM(2,N) - 1.D-4*POSM(2,N)
            SLIX = POSM(2,N)
            CALL FLH_PC30( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RHON(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &        SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,TMHA(2,N),
     &        TMHN(2,N),TMHO(2,N),XGA(2,N),XGN(2,N),XGO(2,N),
     &        XGW(2,N),XHA(2,N),XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),
     &        XNN(2,N),XNO(2,N),XNW(2,N),XMLA(2,N),XMLN(2,N),
     &        XMLO(2,N),XLS(2,N),YLS(2,N),YMGA(2,N),YMGN(2,N),
     &        YMGO(2,N),YMHGA(2,N),YMHGN(2,N),YMHGO(2,N),YMNAX,
     &        YMNNX,YMNOX,ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),
     &        IZN,N,NP_HT )
!
!---      Transition to no hydrate conditions + nonaqueous phases 
!         (40 series -> 20 series)  ---
!
          ELSEIF( SH(2,N).LE.EPSL .AND. PVMX.LT.PVHX .AND.
     &       PVX.GT.(PX*(1.D+0+TOLX)) ) THEN
            IF( ZMCA(2,N).GE.MAX( ZMCO(2,N),ZMCN(2,N) ) ) THEN
              NPHAZ(2,N) = 21 + (100*NP_HT)
            ELSEIF( ZMCO(2,N).GE.MAX( ZMCN(2,N),ZMCA(2,N) ) ) THEN
              NPHAZ(2,N) = 22 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 23 + (100*NP_HT)
            ENDIF
            POSM(2,N) = POSM(2,N) - 1.D-4*POSM(2,N)
            SLIX = POSM(2,N)
            CALL FLH_PC20( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &        SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,
     &        XGA(2,N),XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),
     &        XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),
     &        XNW(2,N),XMLA(2,N),XMLN(2,N),XMLO(2,N),XLS(2,N),
     &        YLS(2,N),YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),
     &        YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,YMNOX,
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
!
!---      Transition to no hydrate conditions + no nonaqueous phases 
!         (40 series -> 10 series)  ---
!
          ELSEIF( SH(2,N).LE.EPSL .AND. PVMX.LT.PVHX  ) THEN
            NPHAZ(2,N) = 10 + (100*NP_HT)
            CALL FLH_PC10( PH(2,N),PG(2,N),PI(2,N),PL(2,N),PN(2,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(2,N),PVA(2,N),PVHA(2,N),
     &        PVN(2,N),PVHN(2,N),PVO(2,N),PVHO(2,N),PVW(2,N),
     &        RHOG(2,N),RHOH(2,N),RHON(2,N),RKG(2,N),
     &        RKL(1,2,N),RKN(2,N),SH(2,N),SI(2,N),SG(2,N),SL(2,N),SLIX,
     &        SN(2,N),T(2,N),TCR(N),TCT(N),TEQHX,
     &        XGA(2,N),XGN(2,N),XGO(2,N),XGW(2,N),XHA(2,N),
     &        XHN(2,N),XHO(2,N),XHW(2,N),XNA(2,N),XNN(2,N),XNO(2,N),
     &        XNW(2,N),XMLA(2,N),XMLN(2,N),XMLO(2,N),XLS(2,N),
     &        YLS(2,N),YMGA(2,N),YMGN(2,N),YMGO(2,N),YMHGA(2,N),
     &        YMHGN(2,N),YMHGO(2,N),YMNAX,YMNNX,YMNOX,
     &        ZMCA(2,N),ZMCO(2,N),ZMCN(2,N),IZN,N,NP_HT )
            POSM(2,N) = SL(2,N) + SI(2,N)
!
!---      No transition from hydrate conditions + no nonaqueous phases 
!         (40 series -> 40 series)  ---
!
          ELSE
            IF( YMHGA(2,N).GE.MAX( YMHGO(2,N),YMHGN(2,N) ) ) THEN
              NPHAZ(2,N) = 41 + (100*NP_HT)
            ELSEIF( YMHGO(2,N).GE.MAX( YMHGN(2,N),YMHGA(2,N) ) ) THEN
              NPHAZ(2,N) = 42 + (100*NP_HT)
            ELSE
              NPHAZ(2,N) = 43 + (100*NP_HT)
            ENDIF
          ENDIF
        ENDIF
  100 CONTINUE
!
!--- Compute increments  ---
!
      DO 300 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 300
        IZN = IZ(N)
        N_DB = N
!
!---    Core phase condition  ---
!
        NPHZX = MOD(NPHAZ(2,N),100)
!
!---    Energy (temperature)  ---
!
        DNR(IEQT,N) = -1.D-7
!
!---    Assign dissolved-salt mass fraction increments,
!       isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          CALL SOL_LS( T(2,N),XLSMX )
          XLSX = MIN( YLS(2,N),XLSMX )
          DNR(IEQS,N) = 1.D-6*XLSMX
        ENDIF
!
!---    Phase condition #10 (10 series)  ---
!
!       SH = 0, SG + SN = 0, SL+SI = 1
!
!       energy - temperature, C - T
!       water mass - pressure, Pa (gauge) - PSO
!       mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!       mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!       mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!       hydrate O2 mass - hydrate O2 mass, kg - TMHA
!       hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!       hydrate N2 mass - hydrate N2 mass, kg - TMHN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        IF( NPHZX.EQ.10 ) THEN
!
!---      Water mass (pressure), increment  ---
!
          DNR(IEQW,N) = MAX( 1.D-2,1.D-6*PSO(2,N) )
!
!---      Mobile O2 mass (O2 vapor pressure), increment  ---
!
          IF( ISLC(45).EQ.0 ) DNR(IEQA,N) = MAX( 1.D-2,1.D-6*PVA(2,N) )
!
!---      Mobile CH4 mass (CH4 vapor pressure), increment  ---
!
          IF( ISLC(64).EQ.0 ) DNR(IEQO,N) = MAX( 1.D-2,1.D-6*PVO(2,N) )
!
!---      Mobile N2 mass (N2 vapor pressure), increment  ---
!
          IF( ISLC(65).EQ.0 ) DNR(IEQN,N) = MAX( 1.D-2,1.D-6*PVN(2,N) )
!
!---      Hydrate O2 mass (hydrate O2 mass)  ---
!
          IF( ISLC(45).EQ.0 ) DNR(IEQHA,N) = 1.D-6
!
!---      Hydrate CH4 mass (hydrate CH4 mass)  ---
!
          IF( ISLC(64).EQ.0 ) DNR(IEQHO,N) = 1.D-6
!
!---      Hydrate N2 mass (hydrate N2 mass)  ---
!
          IF( ISLC(65).EQ.0 ) DNR(IEQHN,N) = 1.D-6
!
!---    Phase conditions 21 through 23 (20 series)  ---
!
!       SH = 0, SG + SN > 0, SL + SI < 1
!       ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #21
!       ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #22
!       ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #23
!
!       energy - temperature, C - T
!       water mass - aqueous + ice saturation - POSM
!       maximum mobile mass component - pressure, Pa (gauge) - PSO
!       mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!       mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!       mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!       hydrate O2 mass - hydrate O2 mass, kg - TMHA
!       hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!       hydrate N2 mass - hydrate N2 mass, kg - TMHN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.21 .AND. NPHZX.LE.23  ) THEN
!
!---      Water mass (aqueous + ice saturation), decrement for
!         high saturation, increment for low saturation  ---
!
          DNR(IEQW,N) = SIGN( 1.D-6,(5.D-1-POSM(2,N)) )
!
!---      Mobile O2 mass (pressure), increment  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( NPHZX.EQ.21 ) THEN
              DNRX = MAX( 1.D-2,1.D-6*ABS(PG(2,N)-PL(2,N)),
     &          1.D-6*ABS(PG(2,N)-PN(2,N)) )
              DNR(IEQA,N) = DNRX
!
!---        Mobile O2 mass (nonaqueous O2 mole fraction), increment  ---
!
            ELSE
              DNR(IEQA,N) = SIGN(1.D-7,5.D-1-ZMCA(2,N))
            ENDIF
          ENDIF
!
!---      Mobile CH4 mass (pressure), increment  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            IF( NPHZX.EQ.22 ) THEN
              DNRX = MAX( 1.D-2,1.D-6*ABS(PG(2,N)-PL(2,N)),
     &          1.D-6*ABS(PG(2,N)-PN(2,N)) )
              DNR(IEQO,N) = DNRX
!
!---        Mobile CH4 mass (nonaqueous CH4 mole fraction), increment  ---
!
            ELSE
              DNR(IEQO,N) = SIGN(1.D-7,5.D-1-ZMCO(2,N))
            ENDIF
          ENDIF
!
!---      Mobile N2 mass (pressure), increment  ---
!
          IF( ISLC(65).EQ.0 ) THEN
            IF( NPHZX.EQ.23  ) THEN
              DNRX = MAX( 1.D-2,1.D-6*ABS(PG(2,N)-PL(2,N)),
     &          1.D-6*ABS(PG(2,N)-PN(2,N)) )
              DNR(IEQN,N) = DNRX
!
!---        Mobile N2 mass (nonaqueous N2 mole fraction), increment  ---
!
            ELSE
              DNR(IEQN,N) = SIGN(1.D-7,5.D-1-ZMCN(2,N))
            ENDIF
          ENDIF
!
!---      Hydrate O2 mass (hydrate O2 mass)  ---
!
          IF( ISLC(45).EQ.0 ) DNR(IEQHA,N) = 1.D-6
!
!---      Hydrate CH4 mass (hydrate CH4 mass)  ---
!
          IF( ISLC(64).EQ.0 ) DNR(IEQHO,N) = 1.D-6
!
!---      Hydrate N2 mass (hydrate N2 mass)  ---
!
          IF( ISLC(65).EQ.0 ) DNR(IEQHN,N) = 1.D-6
!
!---    Phase conditions #31 through #39 (30 series)  ---
!
!       SH > 0, SG + SN > 0, SL + SI < 1
!       ZMCA = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #31
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #32
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #33
!       ZMCO = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #34
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #35
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #36
!       ZMCN = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #37
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #38
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #39
!
!       energy - temperature, C - T
!       water mass - aqueous + ice saturation, POSM
!       maximum mobile mass component - pressure, Pa (gauge) - PSO
!       mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!       mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!       mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!       maximum hydrate mass component - hydrate saturation - SH
!       hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!       hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!       hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.31 .AND. NPHZX.LE.39 ) THEN
!
!---      Water mass (aqueous saturation), decrement for
!         high saturation, increment for low saturation  ---
!
          DNR(IEQW,N) = SIGN( 1.D-6,(5.D-1*(1.D+0-SH(2,N))-POSM(2,N)) )
!
!---      Mobile O2 mass (pressure), increment  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( NPHZX.GE.31 .AND. NPHZX.LE.33  ) THEN
              DNRX = MAX( 1.D-2,1.D-7*ABS(PG(2,N)-PL(2,N)),
     &          1.D-6*ABS(PG(2,N)-PN(2,N)) )
              DNR(IEQA,N) = DNRX
!
!---        Mobile O2 mass (nonaqueous O2 mole fraction), increment  ---
!
            ELSE
              DZMCX = MAX( MIN( 1.D-7,1.D-4*ZMCA(2,N) ),1.D-14 )
              DNR(IEQA,N) = SIGN(DZMCX,5.D-1-ZMCA(2,N))
            ENDIF
          ENDIF
!
!---      Mobile CH4 mass (pressure), increment  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            IF( NPHZX.GE.34 .AND. NPHZX.LE.36  ) THEN
              DNRX = MAX( 1.D-2,1.D-7*ABS(PG(2,N)-PL(2,N)),
     &          1.D-6*ABS(PG(2,N)-PN(2,N)) )
              DNR(IEQO,N) = DNRX
!
!---        Mobile CH4 mass (nonaqueous CH4 mole fraction), increment  ---
!
            ELSE
              DZMCX = MAX( MIN( 1.D-7,1.D-4*ZMCO(2,N) ),1.D-14 )
              DNR(IEQO,N) = SIGN(DZMCX,5.D-1-ZMCO(2,N))
            ENDIF
          ENDIF
!
!---      Mobile N2 mass (pressure), increment  ---
!
          IF( ISLC(65).EQ.0 ) THEN
            IF( NPHZX.GE.37 .AND. NPHZX.LE.39  ) THEN
              DNRX = MAX( 1.D-2,1.D-7*ABS(PG(2,N)-PL(2,N)),
     &          1.D-6*ABS(PG(2,N)-PN(2,N)) )
              DNR(IEQN,N) = DNRX
!
!---        Mobile N2 mass (nonaqueous N2 mole fraction), increment  ---
!
            ELSE
              DZMCX = MAX( MIN( 1.D-7,1.D-4*ZMCN(2,N) ),1.D-14 )
              DNR(IEQN,N) = SIGN(DZMCX,5.D-1-ZMCN(2,N))
            ENDIF
          ENDIF
!
!---      Hydrate O2 mass (hydrate saturation)  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( NPHZX.EQ.31 .OR. NPHZX.EQ.34 .OR. NPHZX.EQ.37 ) THEN
              DSHX = MAX( MIN( 1.D-7,1.D-4*SH(2,N) ),1.D-14 )
              DNR(IEQHA,N) = SIGN(DSHX,5.D-1-SH(2,N))
!
!---        Hydrate O2 mass (hydrate O2 mole fraction of formers)  ---
!
            ELSE
              DYMHX = MAX( MIN( 1.D-7,1.D-4*YMHGA(2,N) ),1.D-14 )
              DNR(IEQHA,N) = SIGN(DYMHX,5.D-1-YMHGA(2,N))
            ENDIF
          ENDIF
!
!---      Hydrate CH4 mass (hydrate saturation)  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            IF( NPHZX.EQ.32 .OR. NPHZX.EQ.35 .OR. NPHZX.EQ.38 ) THEN
              DSHX = MAX( MIN( 1.D-7,1.D-4*SH(2,N) ),1.D-14 )
              DNR(IEQHO,N) = SIGN(DSHX,5.D-1-SH(2,N))
!
!---        Hydrate CH4 mass (hydrate O2 mole fraction of formers)  ---
!
            ELSE
              DYMHX = MAX( MIN( 1.D-7,1.D-4*YMHGO(2,N) ),1.D-14 )
              DNR(IEQHO,N) = SIGN(DYMHX,5.D-1-YMHGO(2,N))
            ENDIF
          ENDIF
!
!---      Hydrate N2 mass (hydrate saturation)  ---
!
          IF( ISLC(65).EQ.0 ) THEN
            IF( NPHZX.EQ.33 .OR. NPHZX.EQ.36 .OR. NPHZX.EQ.39 ) THEN
              DSHX = MAX( MIN( 1.D-7,1.D-4*SH(2,N) ),1.D-14 )
              DNR(IEQHN,N) = SIGN(DSHX,5.D-1-SH(2,N))
!
!---        Hydrate N2 mass (hydrate O2 mole fraction of formers)  ---
!
            ELSE
              DYMHX = MAX( MIN( 1.D-7,1.D-4*YMHGN(2,N) ),1.D-14 )
              DNR(IEQHN,N) = SIGN(DYMHX,5.D-1-YMHGN(2,N))
            ENDIF
          ENDIF
!
!---    Phase conditions #41 through #43 (40 series)  ---
!
!       SH > 0, SG + SN = 0, SL + SI < 1
!       YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #41
!       YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #42
!       YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #43
!
!       energy - temperature, C - T
!       water mass - pressure, Pa (gauge) - PSO
!       mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!       mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!       mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!       maximum hydrate mass component - hydrate saturation - SH
!       hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!       hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!       hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.41 .AND. NPHZX.LE.43 ) THEN
!
!---      Water mass (pressure), increment  ---
!
          DNR(IEQW,N) = MAX( 1.D-2,1.D-8*PSO(2,N) )
!
!---      Mobile O2 mass (O2 vapor pressure), increment  ---
!
          IF( ISLC(45).EQ.0 ) DNR(IEQA,N) = MAX( 1.D-2,1.D-8*PVA(2,N) )
!
!---      Mobile CH4 mass (CH4 vapor pressure), increment  ---
!
          IF( ISLC(64).EQ.0 ) DNR(IEQO,N) = MAX( 1.D-2,1.D-8*PVO(2,N) )
!
!---      Mobile N2 mass (N2 vapor pressure), increment  ---
!
          IF( ISLC(65).EQ.0 ) DNR(IEQN,N) = MAX( 1.D-2,1.D-8*PVN(2,N) )
!
!---      Hydrate O2 mass (hydrate saturation)  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( NPHZX.EQ.41 ) THEN
              DNR(IEQHA,N) = SIGN(1.D-7,5.D-1-SH(2,N))
!
!---        Hydrate O2 mass (hydrate O2 mole fraction of formers)  ---
!
            ELSE
              DYMHX = MAX( MIN( 1.D-7,1.D-4*YMHGA(2,N) ),1.D-14 )
              DNR(IEQHA,N) = SIGN(DYMHX,5.D-1-YMHGA(2,N))
            ENDIF
          ENDIF
!
!---      Hydrate CH4 mass (hydrate saturation)  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            IF( NPHZX.EQ.42 ) THEN
              DNR(IEQHO,N) = SIGN(1.D-7,5.D-1-SH(2,N))
!
!---        Hydrate CH4 mass (hydrate O2 mole fraction of formers)  ---
!
            ELSE
              DYMHX = MAX( MIN( 1.D-7,1.D-4*YMHGO(2,N) ),1.D-14 )
              DNR(IEQHO,N) = SIGN(DYMHX,5.D-1-YMHGO(2,N))
            ENDIF
          ENDIF
!
!---      Hydrate N2 mass (hydrate saturation)  ---
!
          IF( ISLC(65).EQ.0 ) THEN
            IF( NPHZX.EQ.43 ) THEN
              DNR(IEQHN,N) = SIGN(1.D-7,5.D-1-SH(2,N))
!
!---        Hydrate N2 mass (hydrate O2 mole fraction of formers)  ---
!
            ELSE
              DYMHX = MAX( MIN( 1.D-7,1.D-4*YMHGN(2,N) ),1.D-14 )
              DNR(IEQHN,N) = SIGN(DYMHX,5.D-1-YMHGN(2,N))
            ENDIF
          ENDIF
        ENDIF
!
!---    Increment the primary variables  ---
!
        DO 200 M = 3,ISVC+2
          T(M,N) = T(2,N)
          PL(M,N) = PL(2,N)
          PSO(M,N) = PSO(2,N)
          PVO(M,N) = PVO(2,N)
          PVA(M,N) = PVA(2,N)
          PVN(M,N) = PVN(2,N)
          SH(M,N) = SH(2,N)
          SL(M,N) = SL(2,N)
          POSM(M,N) = POSM(2,N)
          YLS(M,N) = YLS(2,N)
          YMHGA(M,N) = YMHGA(2,N)
          YMHGO(M,N) = YMHGO(2,N)
          YMHGN(M,N) = YMHGN(2,N)
          TMHA(M,N) = TMHA(2,N)
          TMHO(M,N) = TMHO(2,N)
          TMHN(M,N) = TMHN(2,N)
          ZMCA(M,N) = ZMCA(2,N)
          ZMCO(M,N) = ZMCO(2,N)
          ZMCN(M,N) = ZMCN(2,N)
!
!---      Phase condition #10 (10 series)  ---
!
!         SH = 0, SG + SN = 0, SL + SI = 1
!
!         energy - temperature, C - T
!         water mass - pressure, Pa (gauge) - PSO
!         mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!         mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!         mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!         hydrate O2 mass - hydrate O2 mass, kg - TMHA
!         hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!         hydrate N2 mass - hydrate N2 mass, kg - TMHN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          IF( NPHZX.EQ.10 ) THEN
            IF( M.EQ.IEQT+2 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PSO(M,N) = PSO(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQA+2 ) THEN
              PVA(M,N) = PVA(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQO+2 ) THEN
              PVO(M,N) = PVO(M,N) + DNR(IEQO,N)
            ELSEIF( M.EQ.IEQN+2  ) THEN
              PVN(M,N) = PVN(M,N) + DNR(IEQN,N)
            ELSEIF( M.EQ.IEQHA+2 ) THEN
              TMHA(M,N) = TMHA(M,N) + DNR(IEQHA,N)
            ELSEIF( M.EQ.IEQHO+2 ) THEN
              TMHO(M,N) = TMHO(M,N) + DNR(IEQHO,N)
            ELSEIF( M.EQ.IEQHN+2 ) THEN
              TMHN(M,N) = TMHN(M,N) + DNR(IEQHN,N)
            ELSEIF( M.EQ.IEQS+2 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Phase conditions 21 through 23 (20 series)  ---
!
!         SH = 0, SG + SN > 0, SL + SI < 1
!         ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #21
!         ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #22
!         ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #23
!
!         energy - temperature, C - T
!         water mass - aqueous + ice saturation - POSM
!         maximum mobile mass component - pressure, Pa (gauge) - PSO
!         mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!         mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!         mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!         hydrate O2 mass - hydrate O2 mass, kg - TMHA
!         hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!         hydrate N2 mass - hydrate N2 mass, kg - TMHN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          ELSEIF( NPHZX.GE.21 .AND. NPHZX.LE.23  ) THEN
            IF( M.EQ.IEQT+2 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              POSM(M,N) = POSM(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQA+2 ) THEN
              IF( NPHZX.EQ.21 ) THEN
                PSO(M,N) = PSO(M,N) + DNR(IEQA,N)
              ELSE
                ZMCA(M,N) = ZMCA(M,N) + DNR(IEQA,N)
              ENDIF
            ELSEIF( M.EQ.IEQO+2 ) THEN
              IF( NPHZX.EQ.22 ) THEN
                PSO(M,N) = PSO(M,N) + DNR(IEQO,N)
              ELSE
                ZMCO(M,N) = ZMCO(M,N) + DNR(IEQO,N)
              ENDIF
            ELSEIF( M.EQ.IEQN+2 ) THEN
              IF( NPHZX.EQ.23 ) THEN
                PSO(M,N) = PSO(M,N) + DNR(IEQN,N)
              ELSE
                ZMCN(M,N) = ZMCN(M,N) + DNR(IEQN,N)
              ENDIF
            ELSEIF( M.EQ.IEQHA+2 ) THEN
              TMHA(M,N) = TMHA(M,N) + DNR(IEQHA,N)
            ELSEIF( M.EQ.IEQHO+2 ) THEN
              TMHO(M,N) = TMHO(M,N) + DNR(IEQHO,N)
            ELSEIF( M.EQ.IEQHN+2 ) THEN
              TMHN(M,N) = TMHN(M,N) + DNR(IEQHN,N)
            ELSEIF( M.EQ.IEQS+2 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Phase conditions #31 through #39 (30 series)  ---
!
!         SH > 0, SG + SN > 0, SL + SI < 1
!         ZMCA = MAX( ZMCA,ZMCO,ZMCN )
!           YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #31
!           YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #32
!           YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #33
!         ZMCO = MAX( ZMCA,ZMCO,ZMCN )
!           YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #34
!           YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #35
!           YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #36
!         ZMCN = MAX( ZMCA,ZMCO,ZMCN )
!           YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #37
!           YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #38
!           YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #39
!
!         energy - temperature, C - T
!         water mass - aqueous + ice saturation, POSM
!         maximum mobile mass component - pressure, Pa (gauge) - PSO
!         mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!         mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!         mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!         maximum hydrate mass component - hydrate saturation - SH
!         hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!         hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!         hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          ELSEIF( NPHZX.GE.31 .AND. NPHZX.LE.39 ) THEN
            IF( M.EQ.IEQT+2 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              POSM(M,N) = POSM(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQA+2 ) THEN
              IF( NPHZX.GE.31 .AND. NPHZX.LE.33 ) THEN
                PSO(M,N) = PSO(M,N) + DNR(IEQA,N)
              ELSE
                ZMCA(M,N) = ZMCA(M,N) + DNR(IEQA,N)
              ENDIF
            ELSEIF( M.EQ.IEQO+2 ) THEN
              IF( NPHZX.GE.34 .AND. NPHZX.LE.36 ) THEN
                PSO(M,N) = PSO(M,N) + DNR(IEQO,N)
              ELSE
                ZMCO(M,N) = ZMCO(M,N) + DNR(IEQO,N)
              ENDIF
            ELSEIF( M.EQ.IEQN+2 ) THEN
              IF( NPHZX.GE.37 .AND. NPHZX.LE.39 ) THEN
                PSO(M,N) = PSO(M,N) + DNR(IEQN,N)
              ELSE
                ZMCN(M,N) = ZMCN(M,N) + DNR(IEQN,N)
              ENDIF
            ELSEIF( M.EQ.IEQHA+2 ) THEN
              IF( NPHZX.EQ.31 .OR. NPHZX.EQ.34 .OR. NPHZX.EQ.37 ) THEN
                SH(M,N) = SH(M,N) + DNR(IEQHA,N)
              ELSE
                YMHGA(M,N) = YMHGA(M,N) + DNR(IEQHA,N)
              ENDIF
            ELSEIF( M.EQ.IEQHO+2 ) THEN
              IF( NPHZX.EQ.32 .OR. NPHZX.EQ.35 .OR. NPHZX.EQ.38 ) THEN
                SH(M,N) = SH(M,N) + DNR(IEQHO,N)
              ELSE
                YMHGO(M,N) = YMHGO(M,N) + DNR(IEQHO,N)
              ENDIF
            ELSEIF( M.EQ.IEQHN+2 ) THEN
              IF( NPHZX.EQ.33 .OR. NPHZX.EQ.36 .OR. NPHZX.EQ.39 ) THEN
                SH(M,N) = SH(M,N) + DNR(IEQHN,N)
              ELSE
                YMHGN(M,N) = YMHGN(M,N) + DNR(IEQHN,N)
              ENDIF           
            ELSEIF( M.EQ.IEQS+2 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
!
!---      Phase conditions #41 through #43 (40 series)  ---
!
!         SH > 0, SG + SN = 0, SL + SI < 1
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #41
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #42
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #43
!
!         energy - temperature, C - T
!         water mass - pressure, Pa (gauge) - PSO
!         mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!         mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!         mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!         maximum hydrate mass component - hydrate saturation - SH
!         hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!         hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!         hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          ELSEIF( NPHZX.GE.41 .AND. NPHZX.LE.43 ) THEN
            IF( M.EQ.IEQT+2 ) THEN
              T(M,N) = T(M,N) + DNR(IEQT,N)
            ELSEIF( M.EQ.IEQW+2 ) THEN
              PSO(M,N) = PSO(M,N) + DNR(IEQW,N)
            ELSEIF( M.EQ.IEQA+2 ) THEN
              PVA(M,N) = PVA(M,N) + DNR(IEQA,N)
            ELSEIF( M.EQ.IEQO+2 ) THEN
              PVO(M,N) = PVO(M,N) + DNR(IEQO,N)
            ELSEIF( M.EQ.IEQN+2  ) THEN
              PVN(M,N) = PVN(M,N) + DNR(IEQN,N)
            ELSEIF( M.EQ.IEQHA+2 ) THEN
              IF( NPHZX.EQ.41 ) THEN
                SH(M,N) = SH(M,N) + DNR(IEQHA,N)
              ELSE
                YMHGA(M,N) = YMHGA(M,N) + DNR(IEQHA,N)
              ENDIF
            ELSEIF( M.EQ.IEQHO+2 ) THEN
              IF( NPHZX.EQ.42 ) THEN
                SH(M,N) = SH(M,N) + DNR(IEQHO,N)
              ELSE
                YMHGO(M,N) = YMHGO(M,N) + DNR(IEQHO,N)
              ENDIF
            ELSEIF( M.EQ.IEQHN+2 ) THEN
              IF( NPHZX.EQ.43 ) THEN
                SH(M,N) = SH(M,N) + DNR(IEQHN,N)
              ELSE
                YMHGN(M,N) = YMHGN(M,N) + DNR(IEQHN,N)
              ENDIF           
            ELSEIF( M.EQ.IEQS+2 ) THEN
              YLS(M,N) = YLS(M,N) + DNR(IEQS,N)
            ENDIF
          ENDIF
  200   CONTINUE
  300 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM_HYD_O2_CH4_N2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Load the current time step values into the old time step
!     variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE HYST
      USE GRID
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVN
      USE FDVI
      USE FDVH
      USE FDVG
      USE FDVA
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LDO_HYD_O2_CH4_N2'
!
!---  Assign old time step values  ---
!
      DO 100 N = 1,NFLD
        BTGL(1,N) = BTGL(2,N)
        DFGA(1,N) = DFGA(2,N)
        DFGO(1,N) = DFGO(2,N)
        DFGW(1,N) = DFGW(2,N)
        DFLA(1,N) = DFLA(2,N)
        DFLO(1,N) = DFLO(2,N)
        DFLN(1,N) = DFLN(2,N)
        DFNO(1,N) = DFNO(2,N)
        DFNN(1,N) = DFNN(2,N)
        DFNW(1,N) = DFNW(2,N)
        NPHAZ(1,N) = NPHAZ(2,N)
        PSO(1,N) = PSO(2,N)
        PG(1,N) = PG(2,N)
        PH(1,N) = PH(2,N)
        PI(1,N) = PI(2,N)
        PL(1,N) = PL(2,N)
        PN(1,N) = PN(2,N)
        PORD(1,N) = PORD(2,N)
        PORT(1,N) = PORT(2,N)
        PVA(1,N) = PVA(2,N)
        PVO(1,N) = PVO(2,N)
        PVN(1,N) = PVN(2,N)
        PVW(1,N) = PVW(2,N)
        PVHA(1,N) = PVHA(2,N)
        PVHO(1,N) = PVHO(2,N)
        PVHN(1,N) = PVHN(2,N)
        TMHA(1,N) = TMHA(2,N)
        TMHO(1,N) = TMHO(2,N)
        TMHN(1,N) = TMHN(2,N)
        TMBP_A(1,N) = TMBP_A(2,N)
        TMBP_O(1,N) = TMBP_O(2,N)
        TMBP_N(1,N) = TMBP_N(2,N)
        RHOG(1,N) = RHOG(2,N)
        RHOH(1,N) = RHOH(2,N)
        RHOI(1,N) = RHOI(2,N)
        RHOL(1,N) = RHOL(2,N)
        RHON(1,N) = RHON(2,N)
        RHOSP(1,N) = RHOSP(2,N)
        RKG(1,N) = RKG(2,N)
        RKL(1,1,N) = RKL(1,2,N)
        RKL(2,1,N) = RKL(2,2,N)
        RKL(3,1,N) = RKL(3,2,N)
        RKN(1,N) = RKN(2,N)
        TMS(1,N) = TMS(2,N)
        SH(1,N) = SH(2,N)
        SI(1,N) = SI(2,N)
        SG(1,N) = SG(2,N)
        SL(1,N) = SL(2,N)
        POSM(1,N) = POSM(2,N)
        SS(1,N) = SS(2,N)
        SN(1,N) = SN(2,N)
        T(1,N) = T(2,N)
        TORG(1,N) = TORG(2,N)
        TORL(1,N) = TORL(2,N)
        TORN(1,N) = TORN(2,N)
        VISG(1,N) = VISG(2,N)
        VISL(1,N) = VISL(2,N)
        VISN(1,N) = VISN(2,N)
        XGA(1,N) = XGA(2,N)
        XGO(1,N) = XGO(2,N)
        XGN(1,N) = XGN(2,N)
        XGW(1,N) = XGW(2,N)
        XHA(1,N) = XHA(2,N)
        XHO(1,N) = XHO(2,N)
        XHN(1,N) = XHN(2,N)
        XHW(1,N) = XHW(2,N)
        XLA(1,N) = XLA(2,N)
        XLO(1,N) = XLO(2,N)
        XLN(1,N) = XLN(2,N)
        XLS(1,N) = XLS(2,N)
        XLW(1,N) = XLW(2,N)
        XNA(1,N) = XNA(2,N)
        XNN(1,N) = XNN(2,N)
        XNO(1,N) = XNO(2,N)
        XNW(1,N) = XNW(2,N)
        XMGA(1,N) = XMGA(2,N)
        XMGO(1,N) = XMGO(2,N)
        XMGN(1,N) = XMGN(2,N)
        XMGW(1,N) = XMGW(2,N)
        XMLA(1,N) = XMLA(2,N)
        XMLO(1,N) = XMLO(2,N)
        XMLN(1,N) = XMLN(2,N)
        XMLW(1,N) = XMLW(2,N)
        XMNA(1,N) = XMNA(2,N)
        XMNN(1,N) = XMNN(2,N)
        XMNO(1,N) = XMNO(2,N)
        XMNW(1,N) = XMNW(2,N)
        YLS(1,N) = YLS(2,N)
        YMGA(1,N) = YMGA(2,N)
        YMGO(1,N) = YMGO(2,N)
        YMGN(1,N) = YMGN(2,N)
        YMHGA(1,N) = YMHGA(2,N)
        YMHGO(1,N) = YMHGO(2,N)
        YMHGN(1,N) = YMHGN(2,N)
        HGA(1,N) = HGA(2,N)
        HGO(1,N) = HGO(2,N)
        HGW(1,N) = HGW(2,N)
        HG(1,N) = HG(2,N)
        HLW(1,N) = HLW(2,N)
        HL(1,N) = HL(2,N)
        HH(1,N) = HH(2,N)
        HI(1,N) = HI(2,N)
        HN(1,N) = HN(2,N)
        HSP(1,N) = HSP(2,N)
        THKL(1,N) = THKL(2,N)
        THKG(1,N) = THKG(2,N)
        THKN(1,N) = THKN(2,N)
        THKH(1,N) = THKH(2,N)
        THKI(1,N) = THKI(2,N)
        UEG(1,N) = UEG(2,N)
        UEL(1,N) = UEL(2,N)
        UEN(1,N) = UEN(2,N)
        ZMCA(1,N) = ZMCA(2,N)
        ZMCO(1,N) = ZMCO(2,N)
        ZMCN(1,N) = ZMCN(2,N)
        DO 90 NSL = 1,NSOLU
          CO(N,NSL) = C(N,NSL)
   90   CONTINUE
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO_HYD_O2_CH4_N2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Compute hydrologic, thermodynamic and physical properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NAPL
      USE JACOB
      USE HYST
      USE HYDT
      USE GRID
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVN
      USE FDVI
      USE FDVH
      USE FDVG
      USE FDVA
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 PHIX(LHF_HT+1,LHF_HT+1)
      REAL*8 HIDX(LHF_HT),HRESX(LHF_HT),XMGX(LHF_HT),XMNX(LHF_HT)
      REAL*8 SP_AREAX,RDMX
!
!----------------------Data Statements---------------------------------!
!
      DATA HFLH / 5.4D+7 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROP_HYD_O2_CH4_N2'
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        IZN = IZ(N)
        N_DB = N
        IF( NITER.LE.3 ) THEN
          ZMIH_HT(1) = YMHGA(2,N)
          ZMIH_HT(2) = YMHGO(2,N)
          ZMIH_HT(3) = YMHGN(2,N)
          SNR(1,N) = YMHGA(2,N)
          SNR(2,N) = YMHGO(2,N)
          SNR(3,N) = YMHGN(2,N)
        ELSE
          ZMIH_HT(1) = SNR(1,N)
          ZMIH_HT(2) = SNR(2,N)
          ZMIH_HT(3) = SNR(3,N)
        ENDIF
!
!---    Core phase condition  ---
!
        NPHZX = MOD(NPHAZ(2,N),100)
!
!---    Loop over increment indices  ---
!
        DO 100 M = 2,ISVC+2
          TKX = T(M,N) + TABS
!
!---      Phase condition #10 (10 series)  ---
!
!         SH = 0, SG + SN = 0, SL + SI = 1
!
!         energy - temperature, C - T
!         water mass - pressure, Pa (gauge) - PSO
!         mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!         mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!         mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!         hydrate O2 mass - hydrate O2 mass, kg - TMHA
!         hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!         hydrate N2 mass - hydrate N2 mass, kg - TMHN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          IF( NPHZX.EQ.10 ) THEN
            CALL FLH_PC10( PH(M,N),PG(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(M,N),PVA(M,N),PVHA(M,N),PVN(M,N),
     &        PVHN(M,N),PVO(M,N),PVHO(M,N),PVW(M,N),RHOG(M,N),RHOH(M,N),
     &        RHON(M,N),RKG(M,N),RKL(1,M,N),RKN(M,N),SH(M,N),
     &        SI(M,N),SG(M,N),SL(M,N),SLIX,SN(M,N),T(M,N),TCR(N),
     &        TCT(N),TEQHX,XGA(M,N),
     &        XGN(M,N),XGO(M,N),XGW(M,N),XHA(M,N),XHN(M,N),XHO(M,N),
     &        XHW(M,N),XNA(M,N),XNN(M,N),XNO(M,N),XNW(M,N),XMLA(M,N),
     &        XMLN(M,N),XMLO(M,N),XLS(M,N),YLS(M,N),YMGA(M,N),
     &        YMGN(M,N),YMGO(M,N),YMHGA(M,N),YMHGN(M,N),YMHGO(M,N),
     &        YMNAX,YMNNX,YMNOX,ZMCA(M,N),ZMCO(M,N),ZMCN(M,N),
     &        IZN,N,NP_HT )
            POSM(M,N) = SL(M,N) + SI(M,N)
!
!---      Phase conditions 21 through 23 (20 series)  ---
!
!         SH = 0, SG + SN > 0, SL + SI < 1
!         ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #21
!         ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #22
!         ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #23
!
!         energy - temperature, C - T
!         water mass - aqueous + ice saturation - POSM
!         maximum mobile mass component - pressure, Pa (gauge) - PSO
!         mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!         mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!         mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!         hydrate O2 mass - hydrate O2 mass, kg - TMHA
!         hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!         hydrate N2 mass - hydrate N2 mass, kg - TMHN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          ELSEIF( NPHZX.GE.21 .AND. NPHZX.LE.23  ) THEN
            IF( NPHZX.EQ.21 ) THEN
              ZMCA(M,N) = MIN( 1.D+0-ZMCO(M,N)-ZMCN(M,N),1.D+0 )
              IF( ZMCA(M,N).LT.EPSL ) ZMCA(M,N) = 0.D+0
            ELSEIF( NPHZX.EQ.22 ) THEN
              ZMCO(M,N) = MIN( 1.D+0-ZMCA(M,N)-ZMCN(M,N),1.D+0 )
              IF( ZMCO(M,N).LT.EPSL ) ZMCO(M,N) = 0.D+0
            ELSE
              ZMCN(M,N) = MIN( 1.D+0-ZMCA(M,N)-ZMCO(M,N),1.D+0 )
              IF( ZMCN(M,N).LT.EPSL ) ZMCN(M,N) = 0.D+0
            ENDIF
            SLIX = POSM(M,N)
            CALL FLH_PC20( PH(M,N),PG(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(M,N),PVA(M,N),PVHA(M,N),PVN(M,N),
     &        PVHN(M,N),PVO(M,N),PVHO(M,N),PVW(M,N),RHOG(M,N),RHOH(M,N),
     &        RKG(M,N),RKL(1,M,N),RKN(M,N),SH(M,N),
     &        SI(M,N),SG(M,N),SL(M,N),SLIX,SN(M,N),T(M,N),TCR(N),
     &        TCT(N),TEQHX,XGA(M,N),
     &        XGN(M,N),XGO(M,N),XGW(M,N),XHA(M,N),XHN(M,N),XHO(M,N),
     &        XHW(M,N),XNA(M,N),XNN(M,N),XNO(M,N),XNW(M,N),XMLA(M,N),
     &        XMLN(M,N),XMLO(M,N),XLS(M,N),YLS(M,N),YMGA(M,N),
     &        YMGN(M,N),YMGO(M,N),YMHGA(M,N),YMHGN(M,N),YMHGO(M,N),
     &        YMNAX,YMNNX,YMNOX,ZMCA(M,N),ZMCO(M,N),ZMCN(M,N),
     &        IZN,N,NP_HT )
!
!---      Phase conditions #31 through #39 (30 series)  ---
!
!         SH > 0, SG + SN > 0, SL + SI < 1
!         ZMCA = MAX( ZMCA,ZMCO,ZMCN )
!           YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #31
!           YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #32
!           YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #33
!         ZMCO = MAX( ZMCA,ZMCO,ZMCN )
!           YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #34
!           YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #35
!           YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #36
!         ZMCN = MAX( ZMCA,ZMCO,ZMCN )
!           YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #37
!           YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #38
!           YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #39
!
!         energy - temperature, C - T
!         water mass - aqueous + ice saturation, POSM
!         maximum mobile mass component - pressure, Pa (gauge) - PSO
!         mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!         mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!         mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!         maximum hydrate mass component - hydrate saturation - SH
!         hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!         hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!         hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          ELSEIF( NPHZX.GE.31 .AND. NPHZX.LE.39 ) THEN
            IF( NPHZX.GE.31 .AND. NPHZX.LE.33 ) THEN
              ZMCA(M,N) = MIN( 1.D+0-ZMCO(M,N)-ZMCN(M,N),1.D+0 )
              IF( ZMCA(M,N).LT.EPSL ) ZMCA(M,N) = 0.D+0
            ELSEIF( NPHZX.GE.34 .AND. NPHZX.LE.36 ) THEN
              ZMCO(M,N) = MIN( 1.D+0-ZMCA(M,N)-ZMCN(M,N),1.D+0 )
              IF( ZMCO(M,N).LT.EPSL ) ZMCO(M,N) = 0.D+0
            ELSE
              ZMCN(M,N) = MIN( 1.D+0-ZMCA(M,N)-ZMCO(M,N),1.D+0 )
              IF( ZMCN(M,N).LT.EPSL ) ZMCN(M,N) = 0.D+0
            ENDIF
            IF( NPHZX.EQ.31 .OR. NPHZX.EQ.34 .OR. NPHZX.EQ.37 ) THEN
              YMHGA(M,N) = MIN( 1.D+0-YMHGO(M,N)-YMHGN(M,N),1.D+0 )
              IF( YMHGA(M,N).LT.EPSL ) YMHGA(M,N) = 0.D+0
            ELSEIF( NPHZX.EQ.32 .OR. NPHZX.EQ.35 .OR. NPHZX.EQ.38 )THEN
              YMHGO(M,N) = MIN( 1.D+0-YMHGA(M,N)-YMHGN(M,N),1.D+0 )
              IF( YMHGO(M,N).LT.EPSL ) YMHGO(M,N) = 0.D+0
            ELSE
              YMHGN(M,N) = MIN( 1.D+0-YMHGA(M,N)-YMHGO(M,N),1.D+0 )
              IF( YMHGN(M,N).LT.EPSL ) YMHGN(M,N) = 0.D+0
            ENDIF
            SLIX = POSM(M,N)
            CALL FLH_PC30( PH(M,N),PG(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(M,N),PVA(M,N),PVHA(M,N),PVN(M,N),
     &        PVHN(M,N),PVO(M,N),PVHO(M,N),PVW(M,N),RHOG(M,N),RHOH(M,N),
     &        RHON(M,N),RKG(M,N),RKL(1,M,N),RKN(M,N),SH(M,N),
     &        SI(M,N),SG(M,N),SL(M,N),SLIX,SN(M,N),T(M,N),TCR(N),
     &        TCT(N),TEQHX,TMHA(M,N),TMHN(M,N),TMHO(M,N),XGA(M,N),
     &        XGN(M,N),XGO(M,N),XGW(M,N),XHA(M,N),XHN(M,N),XHO(M,N),
     &        XHW(M,N),XNA(M,N),XNN(M,N),XNO(M,N),XNW(M,N),XMLA(M,N),
     &        XMLN(M,N),XMLO(M,N),XLS(M,N),YLS(M,N),YMGA(M,N),
     &        YMGN(M,N),YMGO(M,N),YMHGA(M,N),YMHGN(M,N),YMHGO(M,N),
     &        YMNAX,YMNNX,YMNOX,ZMCA(M,N),ZMCO(M,N),ZMCN(M,N),
     &        IZN,N,NP_HT )
!
!---      Phase conditions #41 through #43 (40 series)  ---
!
!         SH > 0, SG + SN = 0, SL + SI < 1
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #41
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #42
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #43
!
!         energy - temperature, C - T
!         water mass - pressure, Pa (gauge) - PSO
!         mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!         mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!         mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!         maximum hydrate mass component - hydrate saturation - SH
!         hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!         hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!         hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!         NaCl mass - total NaCl aqueous mass fraction - YLS
!
          ELSEIF( NPHZX.GE.41 .AND. NPHZX.LE.43 ) THEN
            IF( NPHZX.EQ.41 ) THEN
              YMHGA(M,N) = MIN( 1.D+0-YMHGO(M,N)-YMHGN(M,N),1.D+0 )
              IF( YMHGA(M,N).LT.EPSL ) YMHGA(M,N) = 0.D+0
            ELSEIF( NPHZX.EQ.42 )THEN
              YMHGO(M,N) = MIN( 1.D+0-YMHGA(M,N)-YMHGN(M,N),1.D+0 )
              IF( YMHGO(M,N).LT.EPSL ) YMHGO(M,N) = 0.D+0
            ELSE
              YMHGN(M,N) = MIN( 1.D+0-YMHGA(M,N)-YMHGO(M,N),1.D+0 )
              IF( YMHGN(M,N).LT.EPSL ) YMHGN(M,N) = 0.D+0
            ENDIF
            CALL FLH_PC40( PH(M,N),PG(M,N),PI(M,N),PL(M,N),PN(M,N),
     &        PPEL(N),PPEU(N),PSBX,PSO(M,N),PVA(M,N),PVHA(M,N),PVN(M,N),
     &        PVHN(M,N),PVO(M,N),PVHO(M,N),PVW(M,N),RHOG(M,N),RHOH(M,N),
     &        RHON(M,N),RKG(M,N),RKL(1,M,N),RKN(M,N),SH(M,N),
     &        SI(M,N),SG(M,N),SL(M,N),SLIX,SN(M,N),T(M,N),TCR(N),
     &        TCT(N),TEQHX,TMHA(M,N),TMHN(M,N),TMHO(M,N),XGA(M,N),
     &        XGN(M,N),XGO(M,N),XGW(M,N),XHA(M,N),XHN(M,N),XHO(M,N),
     &        XHW(M,N),XNA(M,N),XNN(M,N),XNO(M,N),XNW(M,N),XMLA(M,N),
     &        XMLN(M,N),XMLO(M,N),XLS(M,N),YLS(M,N),YMGA(M,N),
     &        YMGN(M,N),YMGO(M,N),YMHGA(M,N),YMHGN(M,N),YMHGO(M,N),
     &        YMNAX,YMNNX,YMNOX,ZMCA(M,N),ZMCO(M,N),ZMCN(M,N),
     &        IZN,N,NP_HT )
            POSM(M,N) = SL(M,N) + SI(M,N)
          ENDIF
!
!---      Absolute system and gas pressure  ---
!
          PX = PSO(M,N) + PATM
          PVX = PVA(M,N)+PVO(M,N)+PVN(M,N)+PVW(M,N)
          PVMX = PVA(M,N)+PVO(M,N)+PVN(M,N)
          PVHX = PVHA(M,N)+PVHO(M,N)+PVHN(M,N)
!
!---      Aqueous component fractions and density  ---
!
          WTMLX = (XMLA(M,N)*WTMA + XMLO(M,N)*WTMO + XMLN(M,N)*WTMN +
     &      WTMW*(1.D+0-XMLA(M,N)-XMLO(M,N)-XMLN(M,N)))/
     &      (1.D+0-XLS(M,N)+WTMW*XLS(M,N)/WTMS)
          XLA(M,N) = XMLA(M,N)*WTMA/WTMLX
          XLO(M,N) = XMLO(M,N)*WTMO/WTMLX
          XLN(M,N) = XMLN(M,N)*WTMN/WTMLX
          CALL DENS_B( T(M,N),PX,XLS(M,N),RHOBX )
          CALL DENS_L( T(M,N),RHOBX,XLA(M,N),XLO(M,N),XLN(M,N),
     &      RHOL(M,N) )
          XLW(M,N) = 1.D+0-XLA(M,N)-XLO(M,N)-XLN(M,N)-XLS(M,N)
          IF( XLW(M,N).LT.EPSL ) XLW(M,N) = 0.D+0
          XMLS(M,N) = WTMLX*XLS(M,N)/WTMS
          XMLW(M,N) = WTMLX*XLW(M,N)/WTMW
!
!---      Aqueous molar density  ---
!
          RHOML(M,N) = RHOL(M,N)/WTMLX
!
!---      Ice density kg/m^3  ---
!
          CALL DENS_I( T(M,N),RHOI(M,N) )
!
!---      Porous-media porosity  ---
!
          PPX = MAX( PG(M,N),PL(M,N),PN(M,N) ) + PATM
          CALL PORSTY_HYDT_KE( N,PPX,PCMP(N),PORD(M,N),
     &      PORT(M,N),SH(M,N) )
          PORD(M,N) = MAX( PORD(M,N),EPSL )
          PORT(M,N) = MAX( PORT(M,N),PORD(M,N) )
!
!---      Gas molar fractions and density  ---
!
          RHOGAX = XGA(M,N)*RHOG(M,N)
          RHOGOX = XGO(M,N)*RHOG(M,N)
          RHOGNX = XGN(M,N)*RHOG(M,N)
          RHOGWX = XGW(M,N)*RHOG(M,N)
          WTMGX = 1.D+0/(XGA(M,N)/WTMA + XGO(M,N)/WTMO + XGN(M,N)/WTMN
     &      + XGW(M,N)/WTMW)
          XMGA(M,N) = XGA(M,N)*WTMGX/WTMA
          XMGO(M,N) = XGO(M,N)*WTMGX/WTMO
          XMGN(M,N) = XGN(M,N)*WTMGX/WTMN
          XMGW(M,N) = XGW(M,N)*WTMGX/WTMW
          RHOMG(M,N) = RHOG(M,N)/WTMGX
!
!---      Gas viscosity, Pa s  ---
!
          RHOGAX = RHOMG(M,N)*WTMA
          CALL VISC_O2( T(M,N),RHOGAX,VISGAX )
          RHOGOX = RHOMG(M,N)*WTMO
          CALL VISC_CH4( T(M,N),RHOGOX,VISGOX )
          RHOGNX = RHOMG(M,N)*WTMN
          CALL VISC_N2( T(M,N),RHOGNX,VISGNX )
          RHOGWX = RHOMG(M,N)*WTMW
          CALL VISC_W( T(M,N),PX,RHOGWX,VISGWX )
          CALL VISC_G( PHIX,VISGAX,VISGOX,VISGNX,VISGWX,XMGA(M,N),
     &      XMGO(M,N),XMGN(M,N),XMGW(M,N),VISG(M,N) )
!
!---      Gas diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFGA(M,N) = DFGAC
            DFGO(M,N) = DFGOC
            DFGN(M,N) = DFGNC
            DFGW(M,N) = DFGWC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL DIFC_G( T(M,N),PX,DFGA(M,N),DFGO(M,N),DFGN(M,N),
     &        DFGW(M,N),XMGA(M,N),XMGO(M,N),XMGN(M,N),XMGW(M,N) )
          ENDIF
!
!---      Aqueous viscosity  ---
!
          ISRX = 1
          CALL DENS_W( T(M,N),PX,RHOLWX,RHOX,ISRX )
          CALL VISC_W( T(M,N),PX,RHOLWX,VISLWX )
          CALL VISC_B( T(M,N),XLS(M,N),VISLWX,VISBX )
          CALL VISC_L( XMLA(M,N),XMLO(M,N),XMLN(M,N),
     &      VISBX,VISGAX,VISGOX,VISGNX,VISL(M,N) )
!
!---      Aqueous diffusion coefficients  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLA(M,N) = DFLAC
            DFLO(M,N) = DFLOC
            DFLN(M,N) = DFLNC
            DFLS(M,N) = DFLSC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            CALL DIFC_LA( T(M,N),VISL(M,N),DFLA(M,N) )
            CALL DIFC_LO( T(M,N),VISL(M,N),DFLO(M,N) )
            CALL DIFC_LN( T(M,N),VISL(M,N),DFLN(M,N) )
            CALL DIFC_LS( T(M,N),XLS(M,N),VISL(M,N),DFLS(M,N) )
          ENDIF
!
!---      Nonaqueous-liquid molar fractions and density  ---
!
          RHONAX = XNA(M,N)*RHON(M,N)
          RHONOX = XNO(M,N)*RHON(M,N)
          RHONNX = XNN(M,N)*RHON(M,N)
          RHONWX = XNW(M,N)*RHON(M,N)
          WTMNX = 1.D+0/(XNA(M,N)/WTMA + XNO(M,N)/WTMO + XNN(M,N)/WTMN
     &      + XNW(M,N)/WTMW)
          XMNA(M,N) = XNA(M,N)*WTMNX/WTMA
          XMNO(M,N) = XNO(M,N)*WTMNX/WTMO
          XMNN(M,N) = XNN(M,N)*WTMNX/WTMN
          XMNW(M,N) = XNW(M,N)*WTMNX/WTMW
          RHOMN(M,N) = RHON(M,N)/WTMNX
!
!---      Nonaqueous-liquid viscosity, Pa s  ---
!
          RHONAX = RHOMN(M,N)*WTMA
          CALL VISC_O2( T(M,N),RHONAX,VISNAX )
          RHONOX = RHOMN(M,N)*WTMO
          CALL VISC_CH4( T(M,N),RHONOX,VISNOX )
          RHONNX = RHOMN(M,N)*WTMN
          CALL VISC_N2( T(M,N),RHONNX,VISNNX )
          RHONWX = RHOMN(M,N)*WTMW
          CALL VISC_W( T(M,N),PX,RHONWX,VISNWX )
          CALL VISC_G( PHIX,VISNAX,VISNOX,VISNNX,VISNWX,XMNA(M,N),
     &      XMNO(M,N),XMNN(M,N),XMNW(M,N),VISN(M,N) )
!
!---      Nonaqueous-liquid diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFNA(M,N) = DFNAC
            DFNO(M,N) = DFNOC
            DFNN(M,N) = DFNNC
            DFNW(M,N) = DFNWC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL DIFC_G( T(M,N),PX,DFNA(M,N),DFNO(M,N),DFNN(M,N),
     &        DFNW(M,N),XMNA(M,N),XMNO(M,N),XMNN(M,N),XMNW(M,N) )
          ENDIF
!
!---      Gas thermal conductivity, W/m K  ---
!
          RHOGAX = RHOMG(M,N)*WTMA
          CALL THK_O2( T(M,N),RHOGAX,THKGAX )
          RHOGOX = RHOMG(M,N)*WTMO
          CALL THK_CH4( T(M,N),RHOGOX,THKGOX )
          RHOGNX = RHOMG(M,N)*WTMN
          CALL THK_N2( T(M,N),RHOGNX,THKGNX )
          RHOGWX = RHOMG(M,N)*WTMW
          CALL THK_W( T(M,N),PX,RHOGWX,THKGWX )
          CALL THK_G( PHIX,THKGAX,THKGOX,THKGNX,THKGWX,
     &      XMGA(M,N),XMGO(M,N),XMGN(M,N),XMGW(M,N),THKG(M,N) )
!  
!---      Aqueous thermal conductivity  ---
!
          CALL THK_W( T(M,N),PX,RHOLWX,THKLWX )
          CALL THK_B( T(M,N),XLS(M,N),THKLWX,THKL(M,N) )
!
!---      Nonaqueous thermal conductivity, W/m K  ---
!
          RHONAX = RHOMN(M,N)*WTMA
          CALL THK_O2( T(M,N),RHONAX,THKNAX )
          RHONOX = RHOMN(M,N)*WTMO
          CALL THK_CH4( T(M,N),RHONOX,THKNOX )
          RHONNX = RHOMN(M,N)*WTMN
          CALL THK_N2( T(M,N),RHONNX,THKNNX )
          RHONWX = RHOMN(M,N)*WTMW
          CALL THK_W( T(M,N),PX,RHONWX,THKNWX )
          CALL THK_G( PHIX,THKNAX,THKNOX,THKNNX,THKNWX,
     &      XMNA(M,N),XMNO(M,N),XMNN(M,N),XMNW(M,N),THKN(M,N) )
!
!---      Ice thermal conductivity  ---
!
          CALL THK_I( T(M,N),THKI(M,N) )
!
!---      Hydrate thermal conductivity  ---
!
          CALL THK_H( T(M,N),THKH(M,N) )
!
!---      Pure water vapor enthalpy and internal energy  ---
!
          CALL SP_W( T(M,N),PSWX )
          PEX = MIN( PSWX,PX )
          ISRX = 2
          CALL ENTH_W( T(M,N),PEX,HLWX,HGW(M,N),ISRX )
          UGWX = HGW(M,N) - PX/(RHOGWX+SMALL)
!
!---      Gas enthalpy and internal energy  ---
!
          XMGX(1) = YMGA(M,N)
          XMGX(2) = YMGO(M,N)
          XMGX(3) = YMGN(M,N)
          CALL ENTH_ID( T(M,N),HIDX )
          INDX = 1
          CALL ENTH_RES( PX,T(M,N),HRESX,XMGX,INDX )
          HGA(M,N) = (HIDX(1) + HRESX(1))/WTMA
          HGO(M,N) = (HIDX(2) + HRESX(2))/WTMO
          HGNX = (HIDX(3) + HRESX(3))/WTMN
          HG(M,N) = XGA(M,N)*HGA(M,N) + XGO(M,N)*HGO(M,N) + 
     &      XGN(M,N)*HGNX + XGW(M,N)*HGW(M,N)
          UEG(M,N) = HG(M,N) - PX/(RHOG(M,N)+SMALL)
!
!---      Nonaqueous-liquid enthalpy and internal energy  ---
!
          XMNX(1) = YMNAX
          XMNX(2) = YMNOX
          XMNX(3) = YMNNX
          CALL ENTH_ID( T(M,N),HIDX )
          INDX = 2
          CALL ENTH_RES( PX,T(M,N),HRESX,XMNX,INDX )
          HNAX = (HIDX(1) + HRESX(1))/WTMA
          HNOX = (HIDX(2) + HRESX(2))/WTMO
          HNNX = (HIDX(3) + HRESX(3))/WTMN
          HN(M,N) = XNA(M,N)*HNAX + XNO(M,N)*HNOX + 
     &      XNN(M,N)*HNNX + XNW(M,N)*HGW(M,N)
          UEN(M,N) = HN(M,N) - PX/(RHON(M,N)+SMALL)
!
!---      Pure water liquid enthalpy  ---
!
          CALL SP_W( T(M,N),PSWX )
          PEX = MAX( PSWX,PX )
          ISRX = 1
          CALL ENTH_W( T(M,N),PEX,HLW(M,N),HGWX,ISRX )
!
!---      Aqueous enthalpy  ---
!
          CALL ENTH_B( T(M,N),XLS(M,N),HLW(M,N),HBX )
          FGAX = SG(M,N)*RHOG(M,N)*XGA(M,N)
          FNAX = SN(M,N)*RHON(M,N)*XNA(M,N)
          IF( (FGAX+FNAX)/EPSL.LT.EPSL ) THEN
            HGNAX = 0.D+0
          ELSE
            HGNAX = (HGA(M,N)*FGAX + HNAX*FNAX)/(FGAX + FNAX)
          ENDIF
          FGOX = SG(M,N)*RHOG(M,N)*XGO(M,N)
          FNOX = SN(M,N)*RHON(M,N)*XNO(M,N)
          IF( (FGOX+FNOX)/EPSL.LT.EPSL ) THEN
            HGNOX = 0.D+0
          ELSE
            HGNOX = (HGO(M,N)*FGOX + HNOX*FNOX)/(FGOX + FNOX)
          ENDIF
          FGNX = SG(M,N)*RHOG(M,N)*XGN(M,N)
          FNNX = SN(M,N)*RHON(M,N)*XNN(M,N)
          IF( (FGNX+FNNX)/EPSL.LT.EPSL ) THEN
            HGNNX = 0.D+0
          ELSE
            HGNNX = (HGNX*FGNX + HNNX*FNNX)/(FGNX + FNNX)
          ENDIF
          CALL ENTH_L( T(M,N),XLS(M,N),XLA(M,N),XLO(M,N),XLN(M,N),HBX,
     &      HGNAX,HGNOX,HGNNX,HL(M,N) )
          UEL(M,N) = HL(M,N) - PX/RHOL(M,N)
!
!---      Ice enthalpy  ---
!
          CALL ENTH_I( T(M,N),HI(M,N) )
!
!---      Hydrate enthalpy  ---
!
          XMGX(1) = YMHGA(M,N)
          XMGX(2) = YMHGO(M,N)
          XMGX(3) = YMHGN(M,N)
          INDX = 1
          CALL ENTH_RES( PX,T(M,N),HRESX,XMGX,INDX )
          HGAX = (HIDX(1) + HRESX(1))/WTMA
          HGOX = (HIDX(2) + HRESX(2))/WTMO
          HGNX = (HIDX(3) + HRESX(3))/WTMN
          CALL ENTH_H( PVHX,T(M,N),XHA(M,N),XHO(M,N),XHN(M,N),XHW(M,N),
     &      XLS(M,N),HGAX,HGOX,HGNX,HH(M,N),HLW(M,N) )
!
!---      Precipitated NaCl density, saturation, and enthalpy  ---
!
          CALL DENS_S( T(M,N),PX,RHOSP(M,N) )
          SS(M,N) = MAX(YLS(M,N)-XLS(M,N),0.D+0)*RHOL(M,N)*SL(M,N)/
     &      RHOSP(M,N)
!
!---      Precipitated NaCl enthalpy  ---
!
          CALL ENTH_S( T(M,N),HSP(M,N) )
!
!---      Civan model for permeability reduction
!         Civan, F.C., Scale effect on porosity and permeability: 
!         kinetics, model and correlation, AIChE Journal, 
!         2001, 47(2), 271-287.  ---
!
          IF( IPRF(IZN).EQ.3 ) THEN
!            PORDX = PORD(M,N)*(1.D+0-SH(M,N)-SI(M,N)-SS(M,N))
            PORDX = PORD(M,N)*(1.D+0-SH(M,N))
            PERMRF(M,N) = (PORDX/PORD(M,N))*(((PORDX*(1.D+0-PORD(M,N)))/
     &        (PORD(M,N)*(1.D+0-PORDX)))**(2.D+0*PERM(6,IZN)))
!
!---      Hydrate-porosity permeability reduction model from 
!         IGHCCS2 Benchmark Problem 2  ---
!
          ELSEIF( IPRF(IZN).EQ.4 ) THEN
            VARX = (5.D+0*PERM(4,IZN) + 4.D+0)/(2.D+0*PERM(4,IZN))
            PERMRF1X = (1.D+0-SH(M,N))**VARX
            VARX = (1.D+0-PORD(M,N))/(1.D+0-PERM(6,IZN))
            PERMRF2X = (1.D+0/(VARX**PERM(5,IZN)))**2
            PERMRF(M,N) = PERMRF1X*PERMRF2X
!
!---      Hydrate-porosity permeability reduction model from 
!         IGHCCS2 Benchmark Problem 5  ---
!
          ELSEIF( IPRF(IZN).EQ.5 ) THEN
            VARX = PERM(6,IZN)*((PORD(M,N)/PERM(5,IZN))-1.D+0)
            PERMRF(M,N) = EXP(VARX)*(1.D+0-SH(M,N))**PERM(4,IZN)
          ENDIF
!
!---      Change in porosity as a function of precipitated salt  ---
!
          PORD(M,N) = PORD(M,N)*(1.D+0-SS(M,N))
!
!---      Aqueous, gas, and nonaqueous-liquid tortuosity  ---
!
          PORDX = PORD(M,N)*(1.D+0-SH(M,N)-SI(M,N))
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL(M,N),SG(M,N),SN(M,N),
     &      PORDX,TORL(M,N),TORG(M,N),TORN(M,N) )
!
!---      NaCl volumetric concentration  ---
!
          TMS(M,N) = YLS(M,N)*RHOL(M,N)*SL(M,N)*PORD(M,N)
!
!---      Zero net mass transfer of formers into hydrate  ---
!
          TMBP_A(M,N) = 0.D+0
          TMBP_O(M,N) = 0.D+0
          TMBP_N(M,N) = 0.D+0
!
!---      Effective porosity and specific surface area of hydrate
!
!         Sun, X., Nanchary, N., Mohanty, K.K., 2005. 1-D modeling 
!         of hydrate depressurization in porous media. 
!         Transp. Porous Media 58, 315e338. 
!         http://dx.doi.org/10.1007/s11242-004-1410-x ---
!
          PORDX = PORD(M,N)*(1.D+0-SH(M,N))
          PERMX = MAX( PERM(1,IZN),PERM(2,IZN),PERM(3,IZN) )
          SP_AREAX = SQRT((PORDX**3.D+0)/(2.D+0*PERMX))
          ASL(N) = SP_AREAX*VOL(N)
!
!---      Activation energy for hydrate dissociation, J/kmol  ---
!
          EAX = 81.D+6
!
!---      Hydrate dissociation rate coefficient, kmol/m^2 Pa s ---
!
          CHKND = CHKN(3)*EXP(-(EAX/RCU)/(T(M,N)+273.15D+0))
!
!---      Bound-water saturation ---
!
          SLBX = MIN( SCHR(15,IZN),0.99999D+0 )  
!!
!!---      Independent hydrate formation rates  ---
!!
!          RFMAX = CHKN(2)*PORD(M,N)*VOL(N)*SP_AREAX*
!     &      (PVA(M,N)-PVHA(M,N))
!          RFMNX = CHKN(2)*PORD(M,N)*VOL(N)*SP_AREAX*
!     &      (PVN(M,N)-PVHN(M,N))
!          RFMOX = CHKN(2)*PORD(M,N)*VOL(N)*SP_AREAX*
!     &      (PVO(M,N)-PVHO(M,N))
!!
!!---      Hydrate formation rate, limited by bound water, kmol/s  ---
!!
!          IF( SL(M,N).GT.SLBX ) THEN
!            RFMAX = RFMAX*(SL(M,N)-SLBX)/(1.D+0-SLBX)
!            RFMNX = RFMNX*(SL(M,N)-SLBX)/(1.D+0-SLBX)
!            RFMOX = RFMOX*(SL(M,N)-SLBX)/(1.D+0-SLBX)
!          ELSE
!            RFMAX = 0.D+0
!            RFMNX = 0.D+0
!            RFMOX = 0.D+0
!          ENDIF
!!
!!---      Independent hydrate dissociation rate, kmol/s  ---
!!
!          IF( SH(2,N).GT.EPSL ) THEN
!            RDMAX = CHKND*PORD(M,N)*VOL(N)*SP_AREAX*
!     &        (PVA(M,N)-PVHA(M,N))
!            RDMNX = CHKND*PORD(M,N)*VOL(N)*SP_AREAX*
!     &        (PVN(M,N)-PVHN(M,N))
!            RDMOX = CHKND*PORD(M,N)*VOL(N)*SP_AREAX*
!     &        (PVO(M,N)-PVHO(M,N))
!          ELSE
!            RDMAX = 0.D+0
!            RDMNX = 0.D+0
!            RDMOX = 0.D+0
!          ENDIF
!!
!!---      O2 mobile pressure greater than O2 hydrate 
!!         equilibrium pressure toward hydrate formation  ---
!!
!          IF( (PVA(M,N)-PVHA(M,N)).GE.0.D+0 ) THEN
!            TMBP_A(M,N) = RFMAX*YMGA(M,N)
!!
!!---      O2 mobile pressure less than O2 hydrate 
!!         equilibrium pressure toward hydrate dissociation  ---
!!
!          ELSE 
!            TMBP_A(M,N) = RDMAX*YMHGA(M,N)
!          ENDIF              
!!
!!---      N2 mobile pressure greater than N2 hydrate 
!!         equilibrium pressure toward hydrate formation  ---
!!
!          IF( (PVN(M,N)-PVHN(M,N)).GE.0.D+0 ) THEN
!            TMBP_N(M,N) = RFMNX*YMGN(M,N)
!!
!!---      N2 mobile pressure less than N2 hydrate 
!!         equilibrium pressure toward hydrate dissociation  ---
!!
!          ELSE 
!            TMBP_N(M,N) = RDMNX*YMHGN(M,N)
!          ENDIF              
!!
!!---      CH4 mobile pressure greater than CH4 hydrate 
!!         equilibrium pressure toward hydrate formation  ---
!!
!          IF( (PVO(M,N)-PVHO(M,N)).GE.0.D+0 ) THEN
!            TMBP_O(M,N) = RFMOX*YMGO(M,N)
!!
!!---      CH4 mobile pressure less than CH4 hydrate 
!!         equilibrium pressure toward hydrate dissociation  ---
!!
!          ELSE 
!            TMBP_O(M,N) = RDMOX*YMHGO(M,N)
!          ENDIF              
!          TMBP_A(M,N) = CHKN(6)*TMBP_A(M,N)*WTMA
!          TMBP_N(M,N) = CHKN(7)*TMBP_N(M,N)*WTMN
!          TMBP_O(M,N) = CHKN(8)*TMBP_O(M,N)*WTMO
!
!---      Total mobile and hydrate equilibriumn pressures  ---
!
          PVMX = PVA(M,N) + PVO(M,N) + PVN(M,N)
          PVHX = PVHA(M,N) + PVHO(M,N) + PVHN(M,N)
!
!---      Hydrate formation rate, limited by bound water, kmol/s  ---
!
          RFMX = CHKN(2)*PORD(M,N)*VOL(N)*SP_AREAX*(PVMX-PVHX)
          IF( SL(M,N).GT.SLBX ) THEN
            RFMX = RFMX*(SL(M,N)-SLBX)/(1.D+0-SLBX)
          ELSE
            RFMX = 0.D+0
          ENDIF
!
!---      Hydrate dissociation rate, kmol/s  ---
!
          IF( SH(2,N).GT.EPSL ) THEN
            RDMX = CHKND*PORD(M,N)*VOL(N)*SP_AREAX*(PVMX-PVHX)
          ELSE
            RDMX = 0.D+0
          ENDIF
!
!---      Total mobile pressure greater than total hydrate 
!         equilibrium pressure toward hydrate formation  ---
!
          IF(PVMX-PVHX.GE.0.D+0) THEN
            TMBP_A(M,N) = RFMX*YMGA(M,N)
            TMBP_O(M,N) = RFMX*YMGO(M,N)
            TMBP_N(M,N) = RFMX*YMGN(M,N)
            AST(N) = RFMX/(VOL(N)*SP_AREAX)
!
!---      Total mobile pressure less than total hydrate 
!         equilibrium pressure toward hydrate dissociation  ---
!
          ELSE 
            TMBP_A(M,N) = RDMX*YMHGA(M,N)
            TMBP_O(M,N) = RDMX*YMHGO(M,N)
            TMBP_N(M,N) = RDMX*YMHGN(M,N)  
            AST(N) = RDMX/(VOL(N)*SP_AREAX)   
          ENDIF              
!
!---      Exchange of hydrate formers driven by difference in 
!         mobile and hydrate equilibrium mole fractions  ---
!
          IF( ISLC(73).EQ.0 ) THEN
            IF( SH(2,N).GT.EPSL ) THEN
              PVMHX = 5.D-1*(PVMX+PVHX)
              REMX = CHKN(4)*PVMHX*PORD(M,N)*VOL(N)*SP_AREAX
            ELSE
              REMX = 0.D+0
            ENDIF
!
!---        Total molar rate from mobile to hydrate  ---
!
            TMH_AX = REMX*MAX( (YMGA(M,N)-YMHGA(M,N)),0.D+0 )
            TMH_OX = REMX*MAX( (YMGO(M,N)-YMHGO(M,N)),0.D+0 )
            TMH_NX = REMX*MAX( (YMGN(M,N)-YMHGN(M,N)),0.D+0 )
            TMH_X = TMH_AX + TMH_OX + TMH_NX
!
!---        Factored total molar rate from mobile to hydrate  ---
!
            FMH_AX = CHKN(6)*TMH_AX
            FMH_OX = TMH_OX
            FMH_NX = CHKN(7)*TMH_NX
            FMH_X = FMH_AX + FMH_OX + FMH_NX
!
!---        Total molar rate from hydrate to mobile  ---
!
            THM_AX = REMX*MAX( (YMHGA(M,N)-YMGA(M,N)),0.D+0 )
            THM_OX = REMX*MAX( (YMHGO(M,N)-YMGO(M,N)),0.D+0 )
            THM_NX = REMX*MAX( (YMHGN(M,N)-YMGN(M,N)),0.D+0 )
            THM_X = THM_AX + THM_OX + THM_NX
!
!---        Factored total molar rate from mobile to hydrate  ---
!
            FHM_AX = CHKN(6)*THM_AX
            FHM_OX = THM_OX
            FHM_NX = CHKN(7)*THM_NX
            FHM_X = FHM_AX + FHM_OX + FHM_NX
!
!---        Molar exchange rates  ---
!
            TMBP_AX = 0.D+0
            TMBP_OX = 0.D+0
            TMBP_NX = 0.D+0
            IF( FMH_X/EPSL.GT.EPSL ) THEN
              TMBP_AX = TMBP_AX + FMH_AX*TMH_X/FMH_X
              TMBP_OX = TMBP_OX + FMH_OX*TMH_X/FMH_X
              TMBP_NX = TMBP_NX + FMH_NX*TMH_X/FMH_X
            ENDIF
            IF( FHM_X/EPSL.GT.EPSL ) THEN
              TMBP_AX = TMBP_AX - FHM_AX*THM_X/FHM_X
              TMBP_OX = TMBP_OX - FHM_OX*THM_X/FHM_X
              TMBP_NX = TMBP_NX - FHM_NX*THM_X/FHM_X
            ENDIF
!            TMBP_AX = REMX*(YMGA(M,N)-YMHGA(M,N))
!            TMBP_OX = REMX*(YMGO(M,N)-YMHGO(M,N))
!!
!!---        Balance molar exchange rates  ---
!!
!            TMBP_NX = -(TMBP_AX+TMBP_OX)
!
!---      Exchange of hydrate formers driven by difference in 
!         mobile and hydrate equilibrium partial pressures  ---
!
          ELSEIF( ISLC(73).EQ.1 ) THEN
            IF( SH(2,N).GT.EPSL ) THEN
              TMBP_AX = CHKN(4)*(PVA(M,N)-PVHA(M,N))*
     &          PORD(M,N)*VOL(N)*SP_AREAX
              TMBP_NX = CHKN(4)*(PVN(M,N)-PVHN(M,N))*
     &          PORD(M,N)*VOL(N)*SP_AREAX
              TMBP_OX = CHKN(4)*(PVO(M,N)-PVHO(M,N))*
     &          PORD(M,N)*VOL(N)*SP_AREAX
            ELSE
              TMBP_AX = 0.D+0
              TMBP_NX = 0.D+0
              TMBP_OX = 0.D+0
            ENDIF
!
!---        Balance molar exchange rates  ---
!
            AX = 0.D+0
            BX = 0.D+0
            IF( TMBP_AX.GE.0.D+0 ) THEN
              AX = AX + TMBP_AX
            ELSE
              BX = BX - TMBP_AX
            ENDIF
            IF( TMBP_NX.GE.0.D+0 ) THEN
              AX = AX + TMBP_NX
            ELSE
              BX = BX - TMBP_NX
            ENDIF
            IF( TMBP_OX.GE.0.D+0 ) THEN
              AX = AX + TMBP_OX
            ELSE
              BX = BX - TMBP_OX
            ENDIF
            FX = (1.D+0/AX)*SQRT(AX*BX)
            IF( TMBP_AX.GE.0.D+0 ) THEN
              TMBP_AX = TMBP_AX*FX
            ELSE
              TMBP_AX = TMBP_AX/FX
            ENDIF
            IF( TMBP_NX.GE.0.D+0 ) THEN
              TMBP_NX = TMBP_NX*FX
            ELSE
              TMBP_NX = TMBP_NX/FX
            ENDIF
            IF( TMBP_OX.GE.0.D+0 ) THEN
              TMBP_OX = TMBP_OX*FX
            ELSE
              TMBP_OX = TMBP_OX/FX
            ENDIF
          ENDIF
!
!---      Combined formation/dissocation and exchange  ---
!
          TMBP_A(M,N) = (TMBP_A(M,N)+TMBP_AX)*WTMA
          TMBP_O(M,N) = (TMBP_O(M,N)+TMBP_OX)*WTMO
          TMBP_N(M,N) = (TMBP_N(M,N)+TMBP_NX)*WTMN
!
!---      Rate limits for O2  ---
!
          IF(TMBP_A(M,N).GT.0.D+0) THEN
!
!---        Initial mass of mobile O2
!
            TMMAX0 = PORD(M,N)*VOL(N)*(SG(M,N)*RHOG(M,N)*XGA(M,N)
     &        + SN(M,N)*RHON(M,N)*XNA(M,N) 
     &        + SL(M,N)*RHOL(M,N)*XLA(M,N))
!
!---        Estimate the mass of mobile O2 at the end of the time
!           step, assuming a closed grid cell, with the formulation
!           of Kim et al. 1987. "Kinetics of methane hydrate
!           decomposition," Chemical Engineering Science,
!           42(7):1645-1653  ---
!
            IF( TMMAX0.GT.0.D+0 ) THEN
              TMMAX = TMMAX0*((1.D+0-(TMBP_A(M,N)*DT)/(3.D+0*
     &          TMMAX0))**3.D+0)
!
!---          Limit mass transfer to zero
!
              IF (TMMAX.LT.0.D+0) TMMAX = 0.D+0
!
!---          Convert to linear rate
!
              TMBP_A(M,N) = (TMMAX0-TMMAX)*DTI
            ENDIF
          ELSEIF( TMBP_A(M,N).LT.0.D+0 ) THEN
!
!---        Initial kmol of hydrate O2
!
            TMHAX0 = TMHA(M,N)*PORD(M,N)*VOL(N)
!
!---        Estimate the mass of hydrate O2 at the end of the time
!           step, assuming a closed grid cell, with the formulation
!           of Kim et al. 1987. "Kinetics of methane hydrate
!           decomposition," Chemical Engineering Science,
!           42(7):1645-1653  ---
!
            IF( TMHAX0.GT.0.D+0 ) THEN
              TMHAX = TMHAX0*((1.D+0+(TMBP_A(M,N)*DT)/(3.D+0*
     &          TMHAX0))**3.D+0)
!
!---          Limit mole transfer to zero
!
              IF( TMHAX.LT.0.D+0 ) TMHAX= 0.D+0
!
!---          Convert to linear rate
!
              TMBP_A(M,N)=(TMHAX-TMHAX0)*DTI
            ENDIF
          ENDIF
!
!---      Rate limits for CH4  ---
!
          IF(TMBP_O(M,N).GT.0.D+0) THEN
!
!---        Initial mass of mobile CH4
!
            TMMOX0 = PORD(M,N)*VOL(N)*(SG(M,N)*RHOG(M,N)*XGO(M,N)
     &        + SN(M,N)*RHON(M,N)*XNO(M,N) 
     &        + SL(M,N)*RHOL(M,N)*XLO(M,N))
!
!---        Estimate the mass of mobile CH4 at the end of the time
!           step, assuming a closed grid cell, with the formulation
!           of Kim et al. 1987. "Kinetics of methane hydrate
!           decomposition," Chemical Engineering Science,
!           42(7):1645-1653  ---
!
            IF( TMMOX0.GT.0.D+0 ) THEN
              TMMOX = TMMOX0*((1.D+0-(TMBP_O(M,N)*DT)/(3.D+0*
     &          TMMOX0))**3.D+0)
!
!---          Limit mass transfer to zero
!
              IF (TMMOX.LT.0.D+0) TMMOX = 0.D+0
!
!---          Convert to linear rate
!
              TMBP_O(M,N) = (TMMOX0-TMMOX)*DTI
            ENDIF
          ELSEIF( TMBP_O(M,N).LT.0.D+0 ) THEN
!
!---        Initial kmol of hydrate CH4
!
            TMHOX0 = TMHO(M,N)*PORD(M,N)*VOL(N)
!
!---        Estimate the mass of hydrate CH4 at the end of the time
!           step, assuming a closed grid cell, with the formulation
!           of Kim et al. 1987. "Kinetics of methane hydrate
!           decomposition," Chemical Engineering Science,
!           42(7):1645-1653  ---
!
            IF( TMHOX0.GT.0.D+0 ) THEN
              TMHOX = TMHOX0*((1.D+0+(TMBP_O(M,N)*DT)/(3.D+0*
     &          TMHOX0))**3.D+0)
!
!---          Limit mole transfer to zero
!
              IF( TMHOX.LT.0.D+0 ) TMHOX= 0.D+0
!
!---          Convert to linear rate
!
              TMBP_O(M,N)=(TMHOX-TMHOX0)*DTI
            ENDIF
          ENDIF
!
!---      Rate limits for N2  ---
!
          IF(TMBP_N(M,N).GT.0.D+0) THEN
!
!---        Initial mass of mobile N2
!
            TMMNX0 = PORD(M,N)*VOL(N)*(SG(M,N)*RHOG(M,N)*XGN(M,N)
     &        + SN(M,N)*RHON(M,N)*XNN(M,N) 
     &        + SL(M,N)*RHOL(M,N)*XLN(M,N))
!
!---        Estimate the mass of mobile N2 at the end of the time
!           step, assuming a closed grid cell, with the formulation
!           of Kim et al. 1987. "Kinetics of methane hydrate
!           decomposition," Chemical Engineering Science,
!           42(7):1645-1653  ---
!
            IF( TMMNX0.GT.0.D+0 ) THEN
              TMMNX = TMMNX0*((1.D+0-(TMBP_N(M,N)*DT)/(3.D+0*
     &          TMMNX0))**3.D+0)
!
!---          Limit mass transfer to zero
!
              IF (TMMNX.LT.0.D+0) TMMNX = 0.D+0
!
!---          Convert to linear rate
!
              TMBP_N(M,N) = (TMMNX0-TMMNX)*DTI
            ENDIF
          ELSEIF( TMBP_N(M,N).LT.0.D+0 ) THEN
!
!---        Initial kmol of hydrate N2
!
            TMHNX0 = TMHN(M,N)*PORD(M,N)*VOL(N)
!
!---        Estimate the mass of hydrate N2 at the end of the time
!           step, assuming a closed grid cell, with the formulation
!           of Kim et al. 1987. "Kinetics of methane hydrate
!           decomposition," Chemical Engineering Science,
!           42(7):1645-1653  ---
!
            IF( TMHNX0.GT.0.D+0 ) THEN
              TMHNX = TMHNX0*((1.D+0+(TMBP_N(M,N)*DT)/(3.D+0*
     &          TMHNX0))**3.D+0)
!
!---          Limit mole transfer to zero
!
              IF( TMHNX.LT.0.D+0 ) TMHNX= 0.D+0
!
!---          Convert to linear rate
!
              TMBP_N(M,N)=(TMHNX-TMHNX0)*DTI
            ENDIF
          ENDIF
  100   CONTINUE
!
!---    Aqueous and gas compressibility  ---
!
        IF( NPHZX.EQ.10 ) THEN
          TORL(1,N) = -RHOL(2,N)*
     &     (1.D+0/RHOL(4,N)-1.D+0/RHOL(2,N))/(PL(4,N)-PL(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &     (1.D+0/RHOG(5,N)-1.D+0/RHOG(2,N))/(PG(5,N)-PG(2,N)+SMALL)
        ELSEIF( NPHZX.EQ.21  ) THEN
          TORL(1,N) = -RHOL(2,N)*
     &     (1.D+0/RHOL(5,N)-1.D+0/RHOL(2,N))/(PL(5,N)-PL(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &     (1.D+0/RHOG(5,N)-1.D+0/RHOG(2,N))/(PG(5,N)-PG(2,N)+SMALL)
        ELSEIF( NPHZX.EQ.22  ) THEN
          TORL(1,N) = -RHOL(2,N)*
     &     (1.D+0/RHOL(6,N)-1.D+0/RHOL(2,N))/(PL(6,N)-PL(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &     (1.D+0/RHOG(6,N)-1.D+0/RHOG(2,N))/(PG(6,N)-PG(2,N)+SMALL)
        ELSEIF( NPHZX.EQ.23  ) THEN
          TORL(1,N) = -RHOL(2,N)*
     &     (1.D+0/RHOL(7,N)-1.D+0/RHOL(2,N))/(PL(7,N)-PL(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &     (1.D+0/RHOG(7,N)-1.D+0/RHOG(2,N))/(PG(7,N)-PG(2,N)+SMALL)
        ELSEIF( NPHZX.GE.31 .AND. NPHZX.LE.33 ) THEN
          TORL(1,N) = -RHOL(2,N)*
     &     (1.D+0/RHOL(5,N)-1.D+0/RHOL(2,N))/(PL(5,N)-PL(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &     (1.D+0/RHOG(5,N)-1.D+0/RHOG(2,N))/(PG(5,N)-PG(2,N)+SMALL)
        ELSEIF( NPHZX.GE.34 .AND. NPHZX.LE.36 ) THEN
          TORL(1,N) = -RHOL(2,N)*
     &     (1.D+0/RHOL(6,N)-1.D+0/RHOL(2,N))/(PL(6,N)-PL(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &     (1.D+0/RHOG(6,N)-1.D+0/RHOG(2,N))/(PG(6,N)-PG(2,N)+SMALL)
        ELSEIF( NPHZX.GE.37 .AND. NPHZX.LE.39 ) THEN
          TORL(1,N) = -RHOL(2,N)*
     &     (1.D+0/RHOL(7,N)-1.D+0/RHOL(2,N))/(PL(7,N)-PL(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &     (1.D+0/RHOG(7,N)-1.D+0/RHOG(2,N))/(PG(7,N)-PG(2,N)+SMALL)
        ELSEIF( NPHZX.GE.41 .AND. NPHZX.LE.43 ) THEN
          TORL(1,N) = -RHOL(2,N)*
     &     (1.D+0/RHOL(4,N)-1.D+0/RHOL(2,N))/(PL(4,N)-PL(2,N)+SMALL)
          TORG(1,N) = -RHOG(2,N)*
     &     (1.D+0/RHOG(5,N)-1.D+0/RHOG(2,N))/(PG(5,N)-PG(2,N)+SMALL)
        ENDIF
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDBC_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Read input file for boundary condition information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
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
      CHARACTER*64 ADUM,BDUM(LUK+LSOLU)
      CHARACTER*64 UNTS
      CHARACTER*512 CHDUM
      REAL*8 VAR(LBTM,LBCV)
      REAL*8 XVX(4),YVX(4),ZVX(4)
      INTEGER ITYP(LUK+LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDBC_HYD_O2_CH4_N2'
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
        ENDIF
!
!---  Read boundary types  ---
!
        VARB = 'Boundary Condition Type'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(2))
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 )
     &    CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(3))
!
!---    Solute boundary condition types  ---
!
        IF( IEQC.GT.0 ) THEN
          DO 5 NSL = 1,NSOLU
            BDUM(NSL+LUK) = 'zero flux'
            IDFLT = 1
            VARB = 'Solute Boundary Type'
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM(NSL+LUK))
    5     CONTINUE
        ENDIF
!
!---    Write boundary condition type(s) ---
!
        WRITE(IWR,'(A)') 'Boundary Condition Type: '
!
!---    Read energy boundary condition type ---
!
        IF( INDEX(BDUM(1)(1:),'dirichlet').NE.0 ) THEN
          ITYP(1) = 1
          WRITE(IWR,'(2X,2A)') 'Energy: Dirichlet'
        ELSEIF( INDEX(BDUM(1)(1:),'inflow').NE.0 .AND.
     &    INDEX(BDUM(1)(1:),'outflow').NE.0 ) THEN
          ITYP(1) = 14
          WRITE(IWR,'(2X,2A)') 'Energy: Inflow-Outflow'
        ELSEIF( INDEX(BDUM(1)(1:),'neumann').NE.0 ) THEN
          ITYP(1) = 2
          ITYP(2) = 3
          ITYP(3) = 3
          ITYP(4) = 3
          WRITE(IWR,'(2X,2A)') 'Energy: Neumann'
        ELSEIF( INDEX(BDUM(1)(1:),'zero flux').NE.0 ) THEN
          ITYP(1) = 3
          WRITE(IWR,'(2X,2A)') 'Energy: Zero Flux'
        ELSEIF( INDEX(BDUM(1)(1:),'outflow').NE.0 ) THEN
          ITYP(1) = 7
          WRITE(IWR,'(2X,2A)') 'Energy: Outflow'
        ELSEIF( INDEX(BDUM(1)(1:),'initial cond').NE.0 ) THEN
          ITYP(1) = 12
          WRITE(IWR,'(2X,2A)') 'Energy: Initial Conditions'
        ELSEIF( INDEX(BDUM(1)(1:),'inflow').NE.0 ) THEN
          ITYP(1) = 13
          WRITE(IWR,'(2X,2A)') 'Energy: Inflow'
        ELSE
          INDX = 4
          CHMSG = 'Unknown Energy Boundary Condition: '//BDUM(1)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read flow boundary condition type  ---
!
        IF( ITYP(1).NE.2 ) THEN
          IF( INDEX(BDUM(2)(1:),'hydrostatic').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'outflow').NE.0 ) THEN
            ITYP(2) = 30
            WRITE(IWR,'(2X,2A)') 'Flow: Hydrostatic Outflow'
          ELSEIF( INDEX(BDUM(2)(1:),'hydrostatic').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'inflow').NE.0 ) THEN
            ITYP(2) = 29
            WRITE(IWR,'(2X,2A)') 'Flow: Hydrostatic Inflow'
          ELSEIF( INDEX(BDUM(2)(1:),'hydrostatic').NE.0 ) THEN
            ITYP(2) = 28
          ELSEIF( INDEX(BDUM(2)(1:),'dirichlet').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'outflow').NE.0 ) THEN
            ITYP(2) = 26
            WRITE(IWR,'(2X,2A)') 'Flow: Dirichlet Outflow'
          ELSEIF( INDEX(BDUM(2)(1:),'dirichlet').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'inflow').NE.0 ) THEN
            IF( INDEX(BDUM(2)(1:),'aqueous-only').NE.0 ) THEN
              ITYP(2) = 27
              WRITE(IWR,'(2X,2A)') 'Flow: Dirichlet Aqueous-Only'
            ELSEIF( INDEX(BDUM(2)(1:),'gas-only').NE.0 ) THEN
              ITYP(2) = 37
              WRITE(IWR,'(2X,2A)') 'Flow: Dirichlet Gas-Only'
            ELSEIF( INDEX(BDUM(2)(1:),'nonaqueous-only').NE.0 ) THEN
              ITYP(2) = 47
              WRITE(IWR,'(2X,2A)') 'Flow: Dirichlet Nonaqueous-Only'
            ELSE
              ITYP(2) = 7
              WRITE(IWR,'(2X,2A)') 'Flow: Dirichlet Inflow'
            ENDIF
          ELSEIF( INDEX(BDUM(2)(1:),'dirichlet').NE.0 ) THEN
            IF( INDEX(BDUM(2)(1:),'aqueous-only').NE.0 ) THEN
              ITYP(2) = 21
              WRITE(IWR,'(2X,2A)') 'Flow: Dirichlet Aqueous-Only'
            ELSEIF( INDEX(BDUM(2)(1:),'gas-only').NE.0 ) THEN
              ITYP(2) = 31
              WRITE(IWR,'(2X,2A)') 'Flow: Dirichlet Gas-Only'
            ELSEIF( INDEX(BDUM(2)(1:),'nonaqueous-only').NE.0 ) THEN
              ITYP(2) = 41
              WRITE(IWR,'(2X,2A)') 'Flow: Dirichlet Nonaqueous-Only'
            ELSE
              ITYP(2) = 1
              WRITE(IWR,'(2X,2A)') 'Flow: Dirichlet'
            ENDIF
          ELSEIF( INDEX(BDUM(2)(1:),'neumann').NE.0 ) THEN
            ITYP(2) = 2
            WRITE(IWR,'(2X,2A)') 'Flow: Neumann'
          ELSEIF( INDEX(BDUM(2)(1:),'zero flux').NE.0 ) THEN
            ITYP(2) = 3
            WRITE(IWR,'(2X,2A)') 'Flow: Zero Flux'
          ELSEIF( INDEX(BDUM(2)(1:),'initial').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'condition').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'outflow').NE.0 ) THEN
            ITYP(2) = 15
            WRITE(IWR,'(2X,2A)') 'Flow: Initial Condition Outflow'
          ELSEIF( INDEX(BDUM(2)(1:),'initial').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'condition').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'inflow').NE.0 ) THEN
            ITYP(2) = 14
            WRITE(IWR,'(2X,2A)') 'Flow: Initial Condition Inflow'
          ELSEIF( INDEX(BDUM(2)(1:),'hydrostatic').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'inflow').NE.0 ) THEN
            ITYP(2) = 29
            WRITE(IWR,'(2X,2A)') 'Flow: Hydrostatic Inflow'
          ELSEIF( INDEX(BDUM(2)(1:),'initial').NE.0 .AND.
     &      INDEX(BDUM(2)(1:),'condition').NE.0 ) THEN
            ITYP(2) = 12
            WRITE(IWR,'(2X,2A)') 'Flow: Initial Condition'
          ELSE
            INDX = 4
            CHMSG = 'Unknown Flow Boundary Condition: '
     &        // BDUM(2)
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Read salt boundary condition type ---
!
          IF( ISLC(32).EQ.0 .AND. ITYP(2).NE.3 ) THEN
            IF( INDEX(BDUM(3)(1:),'aqu').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'rel').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'sat').NE.0 ) THEN
              ITYP(3) = 34
              WRITE(IWR,'(2X,2A)') 'Salt: Aqueous Relative Saturation'
            ELSEIF( INDEX(BDUM(3)(1:),'inflow').NE.0 .AND.
     &        (INDEX(BDUM(3)(1:),'rel').NE.0 .OR.
     &        INDEX(BDUM(3)(1:),'sat').NE.0) ) THEN
              ITYP(3) = 35
              WRITE(IWR,'(2X,2A)') 'Salt: Inflow Relative Saturation'
            ELSEIF( INDEX(BDUM(3)(1:),'aqu').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'mass').NE.0 .AND.
     &        INDEX(BDUM(3)(1:),'frac').NE.0 ) THEN
               ITYP(3) = 36
              WRITE(IWR,'(2X,2A)') 'Salt: Aqueous Mass Fraction'
            ELSEIF( INDEX(BDUM(3)(1:),'inflow').NE.0 .AND.
     &        (INDEX(BDUM(3)(1:),'mass').NE.0 .OR.
     &        INDEX(BDUM(3)(1:),'frac').NE.0) ) THEN
              ITYP(3) = 37
              WRITE(IWR,'(2X,2A)') 'Salt: Inflow Mass Fraction'
            ELSEIF( INDEX(BDUM(3)(1:),'zero flux').NE.0 ) THEN
               ITYP(3) = 3
              WRITE(IWR,'(2X,2A)') 'Salt: Zero Flux'
            ELSEIF( INDEX(BDUM(3)(1:),'outflow').NE.0 ) THEN
               ITYP(3) = 7
              WRITE(IWR,'(2X,2A)') 'Salt: Outflow'
            ELSEIF( INDEX(BDUM(3)(1:),'initial cond').NE.0 ) THEN
               ITYP(3) = 12
              WRITE(IWR,'(2X,2A)') 'Salt: Initial Condition'
            ELSE
              INDX = 4
              CHMSG = 'Unknown Salt Boundary Condition: '//BDUM(3)
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
!
!---    Read solute boundary condition type(s) ---
!
          IF( IEQC.GT.0 ) THEN
            DO NSL = 1,NSOLU
              IDB = INDEX( SOLUT(NSL)(1:),'  ') - 1
              IF(INDEX(BDUM(NSL+LUK)(1:),'volumetric conc').NE.0 ) THEN
                 ITYP(NSL+LUK) = 1
                 WRITE(IWR,'(2X,2A)') SOLUT(NSL)(1:IDB),
     &             ': Volumetric Concentration'
              ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'aqueous conc').NE.0 )THEN
                 ITYP(NSL+LUK) = 8
                 WRITE(IWR,'(2X,2A)') SOLUT(NSL)(1:IDB),
     &             ': Aqueous Concentration'
              ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'gas conc').NE.0 ) THEN
                 ITYP(NSL+LUK) = 9
                 WRITE(IWR,'(2X,2A)') SOLUT(NSL)(1:IDB),
     &             ': Gas Concentration'
              ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'id-o2 conc').NE.0 ) THEN
                 ITYP(NSL+LUK) = 10
                 WRITE(IWR,'(2X,2A)') SOLUT(NSL)(1:IDB),
     &             ': ID-O2 Concentration'
              ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'outflow').NE.0 ) THEN
                 ITYP(NSL+LUK) = 7
                 WRITE(IWR,'(2X,2A)') SOLUT(NSL)(1:IDB),
     &             ': Outflow'
              ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'zero flux').NE.0 ) THEN
                 ITYP(NSL+LUK) = 3
                 WRITE(IWR,'(2X,2A)') SOLUT(NSL)(1:IDB),
     &             ': Zero Flux'
              ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'initial con').NE.0 ) THEN
                 ITYP(NSL+LUK) = 12
                 WRITE(IWR,'(2X,2A)') SOLUT(NSL)(1:IDB),
     &             ': Initial Concentration'
              ELSE
                INDX = 4
                CHMSG = 'Unrecognized Solute Boundary Condition: ' //
     &            BDUM(NSL+LUK)
                CALL WRMSGS( INDX )
              ENDIF
            ENDDO
          ENDIF
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
!---  Check boundary domain  ---
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
!---  Read number of boundary times  ---
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
          DO 40 M = 1,LBCV
            VAR(NTM,M) = 0.D+0
   40     CONTINUE
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
     &,': ',VAR(NTM,1)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR(NTM,1),INDX)
!
!---      Dirichlet or inflow energy ---
!
          IF( ITYP(1).EQ.1 .OR. ITYP(1).EQ.13 .OR. ITYP(1).EQ.14 ) THEN
            VARB = 'Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
!
!---      Neumann energy ---
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
          ENDIF
!
!---      Dirichlet flow, Dirichlet-inflow flow, or
!         Dirichlet-outflow flow  ---
!
          IF( MOD(ITYP(2),10).EQ.1 .OR. ITYP(2).EQ.26
     &      .OR. MOD(ITYP(2),10).EQ.7 ) THEN
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
!---      Hydrostatic flow, hydrostatic-inflow flow, or
!         hydrostatic-outflow flow  ---
!
          ELSEIF( ITYP(2).EQ.28 .OR. ITYP(2).EQ.29
     &      .OR. ITYP(2).EQ.30 ) THEN
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
            VARB = 'Gradient'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,5)
            INDX = 0
            IUNM = -2
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,5),INDX)
            VARB = 'Reference Elevation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,6)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
!
!---      Neumann flow  ---
!
          ELSEIF( ITYP(2).EQ.2 ) THEN
            VARB = 'Volumetric Flux'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
            INDX = 0
            IUNM = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
          ENDIF
!
!---      Flow compositions  ---
!
          IF( MOD(ITYP(2),10).EQ.1 .OR. ITYP(2).EQ.2 .OR. ITYP(2).EQ.28
     &      .OR. ITYP(2).EQ.29 .OR. MOD(ITYP(2),10).EQ.7 ) THEN
!
!---        Aqueous saturation  ---
!
            VARB = 'Aqueous saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,4)
            IF( VAR(NTM,4).GT.1.D+0 .OR. VAR(NTM,4).LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Out of Range Value: ' // 
     &          'Aqueous Saturation'
              RLMSG = VAR(NTM,4)
              CALL WRMSGS( INDX )
            ENDIF
!
!---        O2 Mole Fraction of Nonaqueous Components ---
!
            VARB = 'O2 Mole Fraction of Nonaqueous Components'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,10))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,10)
            IF( VAR(NTM,10).GT.1.D+0 .OR. VAR(NTM,10).LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Out of Range Value: Nonaqueous Boundary:' // 
     &          'O2 Mole Fraction of Nonaqueous Components'
              RLMSG = VAR(NTM,10)
              CALL WRMSGS( INDX )
            ENDIF
!
!---        CH4 Mole Fraction of Nonaqueous Components ---
!
            VARB = 'CH4 Mole Fraction of Nonaqueous Components'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,11))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,11)
            IF( VAR(NTM,11).GT.1.D+0 .OR. VAR(NTM,11).LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Out of Range Value: Nonaqueous Boundary: ' // 
     &          'CH4 Mole Fraction of Nonaqueous Components'
              RLMSG = VAR(NTM,11)
              CALL WRMSGS( INDX )
            ENDIF
            VARX = VAR(NTM,10)+VAR(NTM,11)
            IF( VARX.GT.1.D+0 .OR. VARX.LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Out of Range Value: Nonaqueous Boundary: ' // 
     &          'O2+CH4 Mole Fraction of Nonaqueous Components'
              RLMSG = VARX
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Nonaqueous: Water-vapor relative saturation ---
!
            VARB = 'Nonaqueous: Water-Vapor Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,12))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,12)
            IF( VAR(NTM,12).GT.1.D+0 .OR. VAR(NTM,12).LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Out of Range Value: Nonaqueous: ' // 
     &          'Water-Vapor Relative Saturation'
              RLMSG = VAR(NTM,12)
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Aqueous: Relative saturation of nonaqueous components ---
!
            VARB = 'Aqueous: Relative Saturation of ' // 
     &        'Nonaqueous Components'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,8)
            IF( VAR(NTM,8).GT.1.D+0 .OR. VAR(NTM,8).LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Out of Range Value: Aqueous ' // 
     &          'Relative Saturation of Nonaqueous Components'
              RLMSG = VAR(NTM,8)
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
!
!---      Isobrine option  ---
!
          IF( ISLC(32).EQ.0 ) THEN
!
!---        Read dissolved-NaCl aqueous relative saturation ---
!
            IF( ITYP(3).EQ.34 .OR. ITYP(3).EQ.35 ) THEN
              VARB = 'Dissolved-NaCl Aqueous Relative Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,9))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,9)
!
!---        Read dissolved-NaCl aqueous mass fraction ---
!
            ELSEIF( ITYP(3).EQ.36 .OR. ITYP(3).EQ.37 ) THEN
              VARB = 'Dissolved-NaCl Aqueous Mass Fraction'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,9))
              WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &          ': ',VAR(NTM,9)
!
!---        Initial condition, zero flux or outflow  ---
!
            ELSE
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            ENDIF
          ENDIF
!
!---      Read solute transport boundary condition variables ---
!
          IF( IEQC.GT.0 ) THEN
            DO 50 NSL = 1,NSOLU
              IF( ITYP(NSL+LUK).EQ.1 ) THEN
                VARB = 'Volumetric Concentration'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+LBCU))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL),VARB(1:IVR),
     &            ', ',UNTS(1:NCH),': ',VAR(NTM,NSL+LBCU)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VAR(NTM,NSL+LBCU),INDX)
              ELSEIF( ITYP(NSL+LUK).EQ.8 ) THEN
                VARB = 'Aqueous-Phase Volumetric Concentration'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+LBCU))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL),VARB(1:IVR),
     &            ', ',UNTS(1:NCH),': ',VAR(NTM,NSL+LBCU)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VAR(NTM,NSL+LBCU),INDX)
              ELSEIF( ITYP(NSL+LUK).EQ.9 ) THEN
                VARB = 'Gas-Phase Volumetric Concentration'
                CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+LBCU))
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL),VARB(1:IVR),
     &            ', ',UNTS(1:NCH),': ',VAR(NTM,NSL+LBCU)
                INDX = 0
                IUNM = -3
                CALL RDUNIT(UNTS,VAR(NTM,NSL+LBCU),INDX)
              ELSE
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              ENDIF
   50       CONTINUE
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
        DO 106 NTM = 1,IBCMX
          DO 102 M = 1,LBCU
            BC(M,NTM,NB) = VAR(NTM,M)
  102     CONTINUE
          DO 104 NSL = 1,NSOLU
            BC(NSL+LBCU,NTM,NB) = VAR(NTM,NSL+LBCU)
  104     CONTINUE
  106   CONTINUE
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
              IF( ISLC(32).EQ.0 ) IBCT(3,NBC) = ITYP(3)
              IF( IEQC.GT.0 ) THEN
                DO 110 NSL = 1,NSOLU
                  IBCT(NSL+LUK,NBC) = ITYP(NSL+LUK)
  110           CONTINUE
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
  300       CONTINUE
  310     CONTINUE
  320   CONTINUE
        IF( INDEX(ADUM(1:),'file').NE.0 ) CLOSE(UNIT=26)
  400 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDBC_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIC_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Read input file for initial conditions information.
!
!       IC_OPT(1,N) initial condition state option
!       IC_OPT(2,N) 0 hydrate saturation not specified
!       IC_OPT(2,N) 1 hydrate saturation specified
!       IC_OPT(3,N) 0 hydrate CH4 mole frac. of formers not specified
!       IC_OPT(3,N) 1 hydrate CH4 mole frac. of formers specified
!       IC_OPT(4,N) 0 hydrate CH4 mole frac. of formers not specified
!       IC_OPT(4,N) 1 hydrate CH4 mole frac. of formers specified
!       IC_OPT(5,N) 0 hydrate N2 mole frac. of formers not specified
!       IC_OPT(5,N) 1 hydrate N2 mole frac. of formers specified
!       IC_OPT(6,N) 0 aqueous saturation not specified
!       IC_OPT(6,N) 1 aqueous saturation specified
!       IC_OPT(8,N) 0 mobile O2 mole frac. of formers not specified
!       IC_OPT(8,N) 1 mobile O2 mole frac. of formers specified
!       IC_OPT(9,N) 0 mobile CH4 mole frac. of formers not specified
!       IC_OPT(9,N) 1 mobile CH4 mole frac. of formers specified
!       IC_OPT(10,N) 0 mobile N2 mole frac. of formers not specified
!       IC_OPT(10,N) 1 mobile N2 mole frac. of formers specified
!       IC_OPT(8,N) 0 total nonaqueous O2 mole fraction not specified
!       IC_OPT(8,N) 2 total nonaqueous O2 mole fraction specified
!       IC_OPT(9,N) 0 total nonaqueous CH4 mole fraction not specified
!       IC_OPT(9,N) 2 total nonaqueous CH4 mole fraction specified
!       IC_OPT(10,N) 0 total nonaqueous N2 mole fraction not specified
!       IC_OPT(10,N) 2 total nonaqueous N2 mole fraction specified
!       IC_OPT(11,N) 0 aqu. relative saturation of formers not specified
!       IC_OPT(11,N) 1 aqu. relative saturation of formers specified
!       IC_OPT(12,N) salt concentration option
!       IC_OPT(13,N) 0 temperature not specified
!       IC_OPT(13,N) 1 temperature specified
!       IC_OPT(14,N) 0 pressure not specified
!       IC_OPT(14,N) 1 pressure specified
!       IC_OPT(15,N) 0 total nonaqueous saturation not specified
!       IC_OPT(15,N) 1 total nonaqueous saturation specified
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE FILES
      USE FDVP
      USE FDVH
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
      INTEGER IDOM(6)
      REAL*8 VAR(5)
      LOGICAL FCHK
!
!----------------------Data Statements---------------------------------!
!
      DATA CHLB /'X-Direction Gradient, ','Y-Direction Gradient, ',
     &           'Z-Direction Gradient, '/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDIC_HYD_O2_CH4_N2'
      IERR = 0
!
!---  Write card information to ouput file  ---
!
      CARD = 'Initial Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Restart file will be read for initial conditions  ---
!
      IF( IEO.EQ.2 ) THEN
        INDX = 2
        CALL RDRST(INDX)
      ENDIF
!
!---  Read initial conditions  ---
!
      WRITE(IWR,'(/,A)') 'Initial Condition Variable(s) and Domain(s)'
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
        ICX = 0
        IF( INDEX( ADUM(1:),'overwrite').EQ.0 .AND.
     &    ( IEO.EQ.2 .OR. IEO.EQ.6 ) ) GOTO 1000
        IF( INDEX( ADUM(1:),'overwrite').NE.0 .AND.
     &    ( IEO.EQ.2 ) ) IEO = 6 
        IF( INDEX( ADUM(1:),'temperature' ).NE.0 ) THEN
          VARB = 'Initial Temperature'
          IUNK = 1
        ELSEIF( INDEX( ADUM(1:),'pressure' ).NE.0 ) THEN
          VARB = 'Initial Pressure'
          IUNM = -1
          IUNKG = 1
          IUNS = -2
        ELSEIF( INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Hydrate Saturation'
        ELSEIF( INDEX( ADUM(1:),'aqueous' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous Saturation'
        ELSEIF( INDEX( ADUM(1:),'o2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'total' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'nonaqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Total Nonaqueous O2 Mole Fraction'
        ELSEIF( INDEX( ADUM(1:),'o2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mobile' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Mobile O2 Mole Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'o2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mobile' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA) // ': Only Mole Fractions Allowed'
          CALL WRMSGS( INDX )
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'total' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'nonaqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Total Nonaqueous CH4 Mole Fraction'
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mobile' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Mobile CH4 Mole Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mobile' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA) // ': Only Mole Fractions Allowed'
          CALL WRMSGS( INDX )
        ELSEIF( INDEX( ADUM(1:),'n2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'total' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'nonaqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Total Nonaqueous N2 Mole Fraction'
        ELSEIF( INDEX( ADUM(1:),'n2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mobile' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Mobile N2 Mole Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'n2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mobile' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA) // ': Only Mole Fractions Allowed'
          CALL WRMSGS( INDX )
        ELSEIF( INDEX( ADUM(1:),'formers' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqueous' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Aqueous Hydrate Formers ' // 
     &      'Relative Saturation'
        ELSEIF( INDEX( ADUM(1:),'o2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Hydrate O2 Mole Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'o2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA) // ': Only Mole Fractions Allowed'
          CALL WRMSGS( INDX )
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Hydrate CH4 Mole Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA) // ': Only Mole Fractions Allowed'
          CALL WRMSGS( INDX )
        ELSEIF( INDEX( ADUM(1:),'n2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          VARB = 'Initial Hydrate N2 Mole Fraction of Formers'
        ELSEIF( INDEX( ADUM(1:),'n2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mass' ).NE.0 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA) // ': Only Mole Fractions Allowed'
          CALL WRMSGS( INDX )
        ELSEIF( INDEX( ADUM(1:),'ic1' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic-1' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic:1' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic#1' ).NE.0 ) THEN
          VARB = 'Initial Condition State #1'
          ICX = 1
        ELSEIF( INDEX( ADUM(1:),'ic2' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic-2' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic:2' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic#2' ).NE.0 ) THEN
          VARB = 'Initial Condition State #2'
          ICX = 2
        ELSEIF( INDEX( ADUM(1:),'ic3' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic-3' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic:3' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic#3' ).NE.0 ) THEN
          VARB = 'Initial Condition State #3'
          ICX = 3
        ELSEIF( INDEX( ADUM(1:),'ic4' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic-4' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic:4' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic#4' ).NE.0 ) THEN
          VARB = 'Initial Condition State #4'
          ICX = 4
        ELSEIF( INDEX( ADUM(1:),'ic5' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic-5' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic:5' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic#5' ).NE.0 ) THEN
          VARB = 'Initial Condition State #5'
          ICX = 5
        ELSEIF( INDEX( ADUM(1:),'ic6' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic-6' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic:6' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic#6' ).NE.0 ) THEN
          VARB = 'Initial Condition State #6'
          ICX = 6
        ELSEIF( INDEX( ADUM(1:),'ic7' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic-7' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic:7' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'ic#7' ).NE.0 ) THEN
          VARB = 'Initial Condition State #7'
          ICX = 7
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Salt Relative Saturation'
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    ( INDEX( ADUM(1:),'mass' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'frac' ).NE.0 ) ) THEN
          VARB = 'Initial Dissolved-Salt Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'conc' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Salt Concentration'
          IUNKG = 1
          IUNM = -3
        ELSEIF( INDEX( ADUM(1:),'salt' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'vol' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'conc' ).NE.0 ) THEN
          VARB = 'Initial Volumetric Concentration'
          IUNKG = 1
          IUNM = -3
        ELSEIF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          VARB = 'Initial Solute Concentration'
          IUNM = -3
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
        IF( ICX.EQ.0 ) THEN
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(1))
          CALL RDCHR(ISTART,ICOMMA,NCHU,CHDUM,UNTS)
        ENDIF
!
!---    Read initial conditions input from an external file  ---
!
        IF( INDEX( ADUM(1:),'file' ).NE.0 ) THEN
          IF( INDEX( ADUM(1:),'binary' ).NE.0 ) THEN
            IF( ICX.EQ.0 ) THEN
              WRITE(IWR,'(2X,3A)') ADUM(1:NCHA),',',UNTS(1:NCHU)
            ENDIF
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
            IF( ICX.EQ.0 ) THEN
              WRITE(IWR,'(2X,3A)') ADUM(1:NCHA),',',UNTS(1:NCHU)
            ENDIF
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
            IF( ICX.EQ.0 ) THEN
              WRITE(IWR,'(2X,4A,1PE11.4)') ADUM(1:NCHA),
     &          ' (Default Value), ',UNTS(1:NCHU),': ',VAR(1)
            ENDIF
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
            IF( ICX.EQ.0 ) THEN
              INDX = 0
              CALL RDUNIT( UNTS,VAR(1),INDX )
            ENDIF
          ENDIF
!
!---    Read initial conditions according to rock/soil zonations  ---
!
        ELSEIF( INDEX( ADUM(1:),'rock' ).NE.0 .OR.
     &    INDEX( ADUM(1:),'zonation' ).NE.0 ) THEN
          VARB = 'Rock/Soil Name'
          CALL RDCHR(ISTART,ICOMMA,NCHF,CHDUM,FDUM)
!
!---      Search known rock types for a matching type ---
!
          DO 20 M = 1, NROCK
            IF( FDUM .EQ. ROCK(M)) THEN
            IROCK = M
            GOTO 30
          ENDIF
   20     CONTINUE
          INDX = 2
          CHMSG = 'Unrecognized Rock/Soil Type: '//FDUM
          CALL WRMSGS( INDX )
          GOTO 1000
   30     CONTINUE
          IF( ICX.EQ.0 ) THEN
            WRITE(IWR,'(2X,3A,1PE11.4,2A)') ADUM(1:NCHA),UNTS(1:NCHU),
     &        ': ',VAR(1),' Rock/Soil Type: ',FDUM(1:NCHF)
            INDX = 0
            CALL RDUNIT( UNTS,VAR(1),INDX )
          ELSE
            WRITE(IWR,'(2X,3A)') ADUM(1:NCHA),
     &        ': Rock/Soil Type: ',FDUM(1:NCHF)
          ENDIF
!
!---    Read initial condtions input from the input file  ---
!
        ELSE
          IF( ICX.EQ.0 ) THEN
            WRITE(IWR,'(2X,4A,1PE11.4)') ADUM(1:NCHA),', ',
     &        UNTS(1:NCHU),': ',VAR(1)
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
     &          ': ',VAR(I)
              INDX = 0
              IUNM = -1
              CALL RDUNIT( UNTS,VAR(I),INDX )
              WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(I),', 1/m)'
  100       CONTINUE
          ELSE
            WRITE(IWR,'(2X,2A)') 'Initial Condition Option: ',
     &        ADUM(1:NCHA)
          ENDIF
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
!---    Read temperature input  ---
!
        IF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNK = 1
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( T,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( T,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( T,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
            DO 202 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 202
              IC_OPT(13,N) = 1
  202       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( T,VAR(1),ADDER,IROCK,INDX )
            DO 204 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 204
              IF( IZ(N).EQ.IROCK ) IC_OPT(13,N) = 1
  204       CONTINUE
          ELSE
            CALL RDINIS( T,VAR,ADDER,IDOM,INDX )
            DO 206 K = IDOM(5),IDOM(6)
            DO 206 J = IDOM(3),IDOM(4)
            DO 206 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 206
              IC_OPT(13,N) = 1
  206       CONTINUE
          ENDIF
!
!---    Read pressure input  ---
!
        ELSEIF( INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          ADDER = -PATM
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( PSO,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( PSO,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( PSO,VAR,ADDER,UNTS,INDX )
            ENDIF
            CLOSE(UNIT=26)
            DO 212 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 212
              IF( IZ(N).EQ.IROCK ) IC_OPT(14,N) = 1
  212       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PSO,VAR(1),ADDER,IROCK,INDX )
            DO 214 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 214
              IC_OPT(14,N) = 1
  214       CONTINUE
          ELSE
            CALL RDINIS( PSO,VAR,ADDER,IDOM,INDX )
            DO 216 K = IDOM(5),IDOM(6)
            DO 216 J = IDOM(3),IDOM(4)
            DO 216 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 216
              IC_OPT(14,N) = 1
  216       CONTINUE
          ENDIF
!
!---    Read hydrate saturation  ---
!
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'sat').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( SH,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( SH,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( SH,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 222 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 222
              IC_OPT(2,N) = 1
  222       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( SH,VAR(1),ADDER,IROCK,INDX )
            DO 224 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 224
              IF( IZ(N).EQ.IROCK ) IC_OPT(2,N) = 1
  224       CONTINUE
          ELSE
            CALL RDINIS( SH,VAR,ADDER,IDOM,INDX )
            DO 226 K = IDOM(5),IDOM(6)
            DO 226 J = IDOM(3),IDOM(4)
            DO 226 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 226
              IC_OPT(2,N) = 1
  226       CONTINUE
          ENDIF
!
!---    Read aqueous saturation  ---
!
        ELSEIF( INDEX( ADUM(1:),'aqueous' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).EQ.0 ) THEN
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
            DO 242 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 242
              IC_OPT(6,N) = 1
  242       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( SL,VAR(1),ADDER,IROCK,INDX )
            DO 244 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 244
              IF( IZ(N).EQ.IROCK ) IC_OPT(6,N) = 1
  244       CONTINUE
          ELSE
            CALL RDINIS( SL,VAR,ADDER,IDOM,INDX )
            DO 246 K = IDOM(5),IDOM(6)
            DO 246 J = IDOM(3),IDOM(4)
            DO 246 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 246
              IC_OPT(6,N) = 1
  246       CONTINUE
          ENDIF
!
!---    Read total nonaqueous O2 mole fraction  ---
!
        ELSEIF( INDEX( ADUM(1:),'o2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'total' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'nonaqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( XMGA,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( XMGA,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( XMGA,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 252 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 252
              IC_OPT(8,N) = 2
  252       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XMGA,VAR(1),ADDER,IROCK,INDX )
            DO 254 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 254
              IF( IZ(N).EQ.IROCK ) IC_OPT(8,N) = 2
  254       CONTINUE
          ELSE
            CALL RDINIS( XMGA,VAR,ADDER,IDOM,INDX )
            DO 256 K = IDOM(5),IDOM(6)
            DO 256 J = IDOM(3),IDOM(4)
            DO 256 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 256
              IC_OPT(8,N) = 2
  256       CONTINUE
          ENDIF
!
!---    Read mobile O2 mole-fraction of formers  ---
!
        ELSEIF( INDEX( ADUM(1:),'o2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mobile' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( XMGA,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( XMGA,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( XMGA,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 262 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 262
              IC_OPT(8,N) = 1
  262       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XMGA,VAR(1),ADDER,IROCK,INDX )
            DO 264 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 264
              IF( IZ(N).EQ.IROCK ) IC_OPT(8,N) = 1
  264       CONTINUE
          ELSE
            CALL RDINIS( XMGA,VAR,ADDER,IDOM,INDX )
            DO 266 K = IDOM(5),IDOM(6)
            DO 266 J = IDOM(3),IDOM(4)
            DO 266 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 266
              IC_OPT(8,N) = 1
  266       CONTINUE
          ENDIF
!
!---    Read total nonaqueous CH4 mole fraction  ---
!
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'total' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'nonaqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( XMGO,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( XMGO,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( XMGO,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 272 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 272
              IC_OPT(9,N) = 2
  272       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XMGO,VAR(1),ADDER,IROCK,INDX )
            DO 274 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 274
              IF( IZ(N).EQ.IROCK ) IC_OPT(9,N) = 2
  274       CONTINUE
          ELSE
            CALL RDINIS( XMGO,VAR,ADDER,IDOM,INDX )
            DO 276 K = IDOM(5),IDOM(6)
            DO 276 J = IDOM(3),IDOM(4)
            DO 276 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 276
              IC_OPT(9,N) = 2
  276       CONTINUE
          ENDIF
!
!---    Read mobile CH4 mole-fraction of formers  ---
!
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mobile' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( XMGO,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( XMGO,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( XMGO,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 282 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 282
              IC_OPT(9,N) = 1
  282       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XMGO,VAR(1),ADDER,IROCK,INDX )
            DO 284 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 284
              IF( IZ(N).EQ.IROCK ) IC_OPT(9,N) = 1
  284       CONTINUE
          ELSE
            CALL RDINIS( XMGO,VAR,ADDER,IDOM,INDX )
            DO 286 K = IDOM(5),IDOM(6)
            DO 286 J = IDOM(3),IDOM(4)
            DO 286 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 286
              IC_OPT(9,N) = 1
  286       CONTINUE
          ENDIF
!
!---    Read total nonaqueous N2 mole fraction  ---
!
        ELSEIF( INDEX( ADUM(1:),'n2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'total' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'nonaqu' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( XMGN,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( XMGN,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( XMGN,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 292 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 292
              IC_OPT(10,N) = 2
  292       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XMGN,VAR(1),ADDER,IROCK,INDX )
            DO 294 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 294
              IF( IZ(N).EQ.IROCK ) IC_OPT(10,N) = 2
  294       CONTINUE
          ELSE
            CALL RDINIS( XMGN,VAR,ADDER,IDOM,INDX )
            DO 296 K = IDOM(5),IDOM(6)
            DO 296 J = IDOM(3),IDOM(4)
            DO 296 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 296
              IC_OPT(10,N) = 2
  296       CONTINUE
          ENDIF
!
!---    Read mobile N2 mole-fraction of formers  ---
!
        ELSEIF( INDEX( ADUM(1:),'n2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mobile' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( XMGN,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( XMGN,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( XMGN,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 302 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 302
              IC_OPT(10,N) = 1
  302       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XMGN,VAR(1),ADDER,IROCK,INDX )
            DO 304 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 304
              IF( IZ(N).EQ.IROCK ) IC_OPT(10,N) = 1
  304       CONTINUE
          ELSE
            CALL RDINIS( XMGN,VAR,ADDER,IDOM,INDX )
            DO 306 K = IDOM(5),IDOM(6)
            DO 306 J = IDOM(3),IDOM(4)
            DO 306 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 306
              IC_OPT(10,N) = 1
  306       CONTINUE
          ENDIF
!
!---    Read O2 hydrate mole-fraction of formers  ---
!
        ELSEIF( INDEX( ADUM(1:),'o2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( XHA,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( XHA,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( XHA,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 312 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 312
              IC_OPT(3,N) = 1
  312       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XHA,VAR(1),ADDER,IROCK,INDX )
            DO 314 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 314
              IF( IZ(N).EQ.IROCK ) IC_OPT(3,N) = 1
  314       CONTINUE
          ELSE
            CALL RDINIS( XHA,VAR,ADDER,IDOM,INDX )
            DO 316 K = IDOM(5),IDOM(6)
            DO 316 J = IDOM(3),IDOM(4)
            DO 316 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 316
              IC_OPT(3,N) = 1
  316       CONTINUE
          ENDIF
!
!---    Read CH4 hydrate mole-fraction of formers  ---
!
        ELSEIF( INDEX( ADUM(1:),'ch4' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( XHO,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( XHO,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( XHO,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 322 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 322
              IC_OPT(4,N) = 1
  322       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XHO,VAR(1),ADDER,IROCK,INDX )
            DO 324 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 324
              IF( IZ(N).EQ.IROCK ) IC_OPT(4,N) = 1
  324       CONTINUE
          ELSE
            CALL RDINIS( XHO,VAR,ADDER,IDOM,INDX )
            DO 326 K = IDOM(5),IDOM(6)
            DO 326 J = IDOM(3),IDOM(4)
            DO 326 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 326
              IC_OPT(4,N) = 1
  326       CONTINUE
          ENDIF
!
!---    Read N2 hydrate mole-fraction of formers  ---
!
        ELSEIF( INDEX( ADUM(1:),'n2' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'hydrate' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'form' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'mole' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( XHN,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( XHN,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( XHN,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 332 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 332
              IC_OPT(5,N) = 1
  332       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XHN,VAR(1),ADDER,IROCK,INDX )
            DO 334 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 334
              IF( IZ(N).EQ.IROCK ) IC_OPT(5,N) = 1
  334       CONTINUE
          ELSE
            CALL RDINIS( XHN,VAR,ADDER,IDOM,INDX )
            DO 336 K = IDOM(5),IDOM(6)
            DO 336 J = IDOM(3),IDOM(4)
            DO 336 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 336
              IC_OPT(5,N) = 1
  336       CONTINUE
          ENDIF
!
!---    Read aqueous dissolved formers relative saturation  ---
!
        ELSEIF( INDEX( ADUM(1:),'formers' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'aqueous' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'rel' ).NE.0 .AND.
     &    INDEX( ADUM(1:),'sat' ).NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( XMLO,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( XMLO,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( XMLO,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 342 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 342
              IC_OPT(11,N) = 1
  342       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XMLO,VAR(1),ADDER,IROCK,INDX )
            DO 344 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 344
              IF( IZ(N).EQ.IROCK ) IC_OPT(11,N) = 1
  344       CONTINUE
          ELSE
            CALL RDINIS( XMLO,VAR,ADDER,IDOM,INDX )
            DO 346 K = IDOM(5),IDOM(6)
            DO 346 J = IDOM(3),IDOM(4)
            DO 346 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 346
              IC_OPT(11,N) = 1
  346       CONTINUE
          ENDIF
!
!---    Read initial condition state, using the ice pressure
!       variable as a dummy  ---
!
        ELSEIF( INDEX(ADUM(1:),'ic').NE.0 .AND.
     &    ICX.GE.1 .AND. ICX.LE.10 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'binary').NE.0 ) THEN
              CALL RDINBS( PI,ADDER,UNTS,INDX )
            ELSEIF( INDEX(ADUM(1:),'ascii').NE.0 ) THEN
              CALL RDINAS( PI,ADDER,UNTS,INDX )
            ELSE
              CALL RDINFS( PI,VAR,ADDER,UNTS,INDX )
            ENDIF
            DO 352 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 352
              IC_OPT(1,N) = ICX
  352       CONTINUE
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PI,VAR(1),ADDER,IROCK,INDX )
            DO 354 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 354
              IF( IZ(N).EQ.IROCK ) IC_OPT(1,N) = ICX
  354       CONTINUE
          ELSE
            CALL RDINIS( PI,VAR,ADDER,IDOM,INDX )
            DO 356 K = IDOM(5),IDOM(6)
            DO 356 J = IDOM(3),IDOM(4)
            DO 356 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 356
              IC_OPT(1,N) = ICX
  356       CONTINUE
          ENDIF
!
!---    Read salt/inhibitor concentration  ---
!
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'vol').NE.0 .AND.
     &      INDEX(ADUM(1:),'conc').NE.0 ) THEN
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
          ELSE
            IVAR = 3
          ENDIF
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IF( IVAR.EQ.0 .OR. IVAR.EQ.1 ) THEN
              IUNM = -3
              IUKG = 1
            ENDIF
            CALL RDINFS( TMS,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 362 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 362
              IC_OPT(12,N) = IVAR
  362       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( TMS,VAR(1),ADDER,IROCK,INDX )
            DO 364 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 364
              IF( IZ(N).EQ.IROCK ) IC_OPT(12,N) = IVAR
  364       CONTINUE
          ELSE
            CALL RDINIS( TMS,VAR,ADDER,IDOM,INDX )
            DO 366 K = IDOM(5),IDOM(6)
            DO 366 J = IDOM(3),IDOM(4)
            DO 366 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 366
              IC_OPT(12,N) = IVAR
  366       CONTINUE
          ENDIF
!
!---    Read solute concentration  ---
!
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
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: ' //
     &      ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIC_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINPT_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Read input file cards.
!     Direct control to card reader subroutines.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      SUB_LOG(ISUB_LOG) = '/RDINPT_HYD_O2_CH4_N2'
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
      INDX = 4
      CHMSG = 'Missing Simulation Title Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
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
          INDX = 4
          CHMSG = 'Incompatible Operational Mode'
          CALL WRMSGS( INDX )
        ENDIF
        REWIND(IRD)
        GOTO 300
      ELSE
        GOTO 200
      ENDIF
  210 CONTINUE
      INDX = 4
      CHMSG = 'Missing Solution Control Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for grid card  ---
!
  300 CONTINUE
  309 READ(IRD,'(A)', END=310) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 309
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'grid').NE.0 ) THEN
        CALL RDGRID
        REWIND(IRD)
        GOTO 400
      ELSE
        GOTO 300
      ENDIF
  310 CONTINUE
      INDX = 4
      CHMSG = 'Missing Grid Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for rock/soil zonation card  ---
!
  400 CONTINUE
  409 READ(IRD,'(A)', END=410) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 409
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'rock/soil').NE.0 ) THEN
        CALL RDROCK
        REWIND(IRD)
        GOTO 500
      ELSE
        GOTO 400
      ENDIF
  410 CONTINUE
      INDX = 4
      CHMSG = 'Missing Rock/Soil Zonation Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for inactive nodes card  ---
!
  500 CONTINUE
  509 READ(IRD,'(A)', END=510) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 509
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'inactive').NE.0 ) THEN
        CALL RDINAC
        REWIND(IRD)
        GOTO 600
      ELSE
        GOTO 500
      ENDIF
  510 CONTINUE
      INDX = 1
      CHMSG = 'Missing Inactive Nodes Card'
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
      INDX = 4
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
      INDX = 4
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
        CALL RDSP_HYD_O2_CH4_N2
        REWIND(IRD)
        GOTO 900
      ELSE
        GOTO 800
      ENDIF
  810 CONTINUE
      INDX = 4
      CHMSG = 'Missing Saturation Function Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for aqueous relative permeability card  ---
!
  900 CONTINUE
  909 READ(IRD,'(A)', END=910) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 909
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'aqueous rel').NE.0 ) THEN
        CALL RDLRP
        REWIND(IRD)
        GOTO 1000
      ELSE
        GOTO 900
      ENDIF
  910 CONTINUE
      INDX = 4
      CHMSG = 'Missing Aqueous Relative Permeability Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for nonaqueous liquid relative 
!     permeability card  ---
!
 1000 CONTINUE
 1009 READ(IRD,'(A)', END=1010) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 1009
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  (INDEX(CHDUM(2:),'nonaqueous').NE.0 .OR.
     &  INDEX(CHDUM(2:),'non-aqueous').NE.0) .AND.
     &  INDEX(CHDUM(2:),'liq').NE.0 .AND.
     &  INDEX(CHDUM(2:),'rel').NE.0 ) THEN
        CALL RDNRP
        REWIND(IRD)
        GOTO 1100
      ELSE
        GOTO 1000
      ENDIF
 1010 CONTINUE
      INDX = 4
      CHMSG = 'Missing Nonaqueous Liquid Relative Permeability Card'
      CALL WRMSGS( INDX )
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
      INDX = 4
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
        GOTO 2000
      ELSE
        GOTO 1200
      ENDIF
 1250 CONTINUE
      INDX = 4
      CHMSG = 'Missing Thermal Properties Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for solute/fluid interaction card --
!
 2000 CONTINUE
 2009 READ(IRD,'(A)', END=2010) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2009
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'solute/fluid').NE.0 ) THEN
        CALL RDTF_HYD_O2_CH4_N2
        REWIND(IRD)
        GOTO 2100
      ELSE
        GOTO 2000
      ENDIF
 2010 CONTINUE
      IF( IEQC.EQ.0 ) THEN
        REWIND(IRD)
      ELSE
        INDX = 4
        CHMSG = 'Missing Solute/Fluid Interaction Card'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Search input file for solute/rock interaction card --
!
 2100 CONTINUE
 2109 READ(IRD,'(A)', END=2110) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2109
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'solute/porous').NE.0 ) THEN
        CALL RDTP_HYD_O2_CH4_N2
        REWIND(IRD)
        GOTO 2200
      ELSE
        GOTO 2100
      ENDIF
 2110 CONTINUE
      IF( IEQC.EQ.0 ) THEN
       REWIND(IRD)
      ELSE
        INDX = 4
        CHMSG = 'Missing Solute/Porous Media Interaction Card'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Search input file for salt transport card --
!
 2200 CONTINUE
 2209 READ(IRD,'(A)', END=2210) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2209
      CALL LCASE( CHDUM )
      IF( CHDUM(1:15).EQ.'~salt transport' ) THEN
        CALL RDST_HYD_O2_CH4_N2
        REWIND(IRD)
        GOTO 4000
      ELSE
        GOTO 2200
      ENDIF
 2210 CONTINUE
      INDX = 4
      CHMSG = 'Missing Salt Transport Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for initial conditions card --
!
 4000 CONTINUE
 4009 READ(IRD,'(A)', END=4010) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4009
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'initial').NE.0 ) THEN
        CALL RDIC_HYD_O2_CH4_N2
        REWIND(IRD)
        GOTO 4100
      ELSE
        GOTO 4000
      ENDIF
 4010 CONTINUE
      IF( IEO.EQ.2 ) THEN
        INDX = 1
        CHMSG = 'Missing Initial Conditions Card'
        CALL WRMSGS( INDX )
        INDX = 2
        CALL RDRST(INDX)
        REWIND(IRD)
      ELSE
        INDX = 4
        CHMSG = 'Missing Initial Conditions Card'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Search input file for boundary conditions card --
!
 4100 CONTINUE
 4109 READ(IRD,'(A)', END=4110) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4109
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &  INDEX(CHDUM(2:),'boundary').NE.0 .AND.
     &  INDEX(CHDUM(2:),'geomech').EQ.0 ) THEN
        CALL RDBC_HYD_O2_CH4_N2
        REWIND(IRD)
        GOTO 4200
      ELSE
        GOTO 4100
      ENDIF
 4110 CONTINUE
      INDX = 1
      CHMSG = 'Missing Boundary Conditions Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for source card --
!
 4200 CONTINUE
 4209 READ(IRD,'(A)', END= 4210) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4209
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'source').NE.0 ) THEN
        CALL RDSR_HYD_O2_CH4_N2
        REWIND(IRD)
        GOTO 4300
      ELSE
        GOTO 4200
      ENDIF
 4210 CONTINUE
      INDX = 1
      CHMSG = 'Missing Source Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for output control card --
!
 4300 CONTINUE
 4309 READ(IRD,'(A)', END=4310) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4309
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'output').NE.0 ) THEN
        CALL RDOU_HYD_O2_CH4_N2
        REWIND(IRD)
        GOTO 4400
      ELSE
        GOTO 4300
      ENDIF
 4310 CONTINUE
      INDX = 1
      CHMSG = 'Missing Output Control Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
!
!---  Search input file for surface flux card --
!
 4400 CONTINUE
 4409 READ(IRD,'(A)', END=4410) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4409
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'surface').NE.0 ) THEN
        CALL RDSF_HYD_O2_CH4_N2
        REWIND(IRD)
        GOTO 4500
      ELSE
        GOTO 4400
      ENDIF
 4410 CONTINUE
      INDX = 1
      CHMSG = 'Missing Surface Flux Card'
      CALL WRMSGS( INDX )
      REWIND(IRD)
 4500 CONTINUE
!
!---  Geomechanics  ---
!
      IF( ISLC(50).NE.0 ) THEN
 7100   CONTINUE
!
!---  Search input file for inactive nodes card  ---
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
!---    Search input file for geomechanics property card --
!
 7209   READ(IRD,'(A)', END=7210) CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 7209
        CALL LCASE( CHDUM )
          IF( CHDUM(1:1).EQ.'~' .AND.
     &      INDEX(CHDUM(2:),'geomech').NE.0 .AND.
     &      INDEX(CHDUM(2:),'prop').NE.0 ) THEN
          CALL RDGMP
          REWIND(IRD)
          GOTO 7300
        ELSE
          GOTO 7200
        ENDIF
 7210   CONTINUE
!
!---    Geomechanical simulations  ---
!
        IF( ISLC(50).NE.0 ) THEN
          INDX = 4
          CHMSG = 'Missing Geomechanical Properties Card'
          CALL WRMSGS( INDX )
        ELSE
          REWIND(IRD)
        ENDIF
 7300   CONTINUE
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
          GOTO 7300
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
!---  End of input record --
!
      CARD = 'End of Input Record'
      ICD = INDEX( CARD,'  ' )-1
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINPT_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDOU_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Read input file for output information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE OUTPU
      USE GRID
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE FDVA
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
      CHARACTER*64 ADUM,UNTS,SOLNM
      CHARACTER*512 CHDUM
      CHARACTER*6 FORM
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM / '(I6,$)' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDOU_HYD_O2_CH4_N2'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Output Control Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read reference node information  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Reference Nodes'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NREF)
      IF( NREF.GT.LREF ) THEN
        INDX = 5
        CHMSG = 'Number of Reference Nodes > Parameter LREF'
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(IWR,'(/,A,I6)') 'Reference Node No. and Indices: ',NREF
      DO 100 N = 1,NREF
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Reference Node Index'
        CALL RDINT(ISTART,ICOMMA,CHDUM,IRF)
        CALL RDINT(ISTART,ICOMMA,CHDUM,JRF)
        CALL RDINT(ISTART,ICOMMA,CHDUM,KRF)
        IF( IRF.LT.1 .OR. IRF.GT.IFLD ) THEN
          INDX = 7
          CHMSG = 'Unrecognized Reference Node I Index'
          IMSG = IRF
          CALL WRMSGS( INDX )
        ENDIF
        IF( JRF.LT.1 .OR. JRF.GT.JFLD ) THEN
          INDX = 7
          CHMSG = 'Unrecognized Reference Node J Index'
          IMSG = JRF
          CALL WRMSGS( INDX )
        ENDIF
        IF( KRF.LT.1 .OR. KRF.GT.KFLD) THEN
          INDX = 7
          CHMSG = 'Unrecognized Reference Node K Index'
          IMSG = KRF
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
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
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
      IF( ICS.EQ.2 .OR. ICS.EQ.6 .OR. ICS.EQ.8 ) THEN
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
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Reference Node Variables'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NVREF)
      WRITE( IWR,'(/,A,I6)') 'Reference Node Variables: ',NVREF
      NVC = 0
      DO 200 NV = 1,NVREF
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Reference Node Variable'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Reference Node Variable: Solute Name'
          CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SOLNM)
          DO 110 NSL = 1,NSOLU
            IF( SOLNM.EQ.SOLUT(NSL) ) GOTO 120
  110     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute Name: '//SOLNM
          CALL WRMSGS( INDX )
          NVC = NVC -1
          GOTO 200
  120     CONTINUE
        ENDIF
        IF( (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'liq').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IREF(NV) = 3
        ELSEIF( INDEX(ADUM(1:),'aqueous compress').NE.0 ) THEN
          IREF(NV) = 103
        ELSEIF( INDEX(ADUM(1:),'gas compress').NE.0 ) THEN
          IREF(NV) = 104
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IREF(NV) = 1
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IREF(NV) = 2
        ELSEIF( INDEX(ADUM(1:),'critical').NE.0 .AND.
     &    INDEX(ADUM(1:),'temp').NE.0 ) THEN
          IREF(NV) = 206
          CHREF(206) = 'TCR'
        ELSEIF( INDEX(ADUM(1:),'cricondentherm').NE.0 .AND.
     &    INDEX(ADUM(1:),'temp').NE.0 ) THEN
          IREF(NV) = 213
          CHREF(213) = 'TCT'
        ELSEIF( INDEX(ADUM(1:),'temp').NE.0 ) THEN
          IREF(NV) = 4
        ELSEIF( INDEX(ADUM(1:),'phase').NE.0 .AND.
     &    INDEX(ADUM(1:),'cond').NE.0 ) THEN
          IREF(NV) = 5
        ELSEIF( (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'liq').NE.0 .AND.
     &    INDEX(ADUM(1:),'gauge').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IREF(NV) = 8
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'gauge').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IREF(NV) = 6
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'gauge').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IREF(NV) = 7
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'area').NE.0 ) THEN
          IREF(NV) = 9
          CHREF(9) = 'HYDA'
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'reaction').NE.0 ) THEN
          IREF(NV) = 10
          CHREF(10) = 'HYDR'
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'equil').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IREF(NV) = 386
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'formers').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IREF(NV) = 387
        ELSEIF( INDEX(ADUM(1:),'apparent').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IREF(NV) = 9
        ELSEIF( (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'liq').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IREF(NV) = 13
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IREF(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IREF(NV) = 12
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'moisture').NE.0 .AND.
     &    INDEX(ADUM(1:),'content').NE.0 ) THEN
          IREF(NV) = 15
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 ) THEN
          IREF(NV) = 165
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IREF(NV) = 294
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 295
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 ) THEN
          IREF(NV) = 221
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IREF(NV) = 379
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IREF(NV) = 380
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 ) THEN
          IREF(NV) = 162
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IREF(NV) = 381
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IREF(NV) = 382
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IREF(NV) = 383
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IREF(NV) = 384
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IREF(NV) = 385
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          CHREF(204) = 'XMLA'
          IREF(NV) = 204
        ELSEIF( (INDEX(ADUM(1:),'salt').NE.0 .OR.
     &    INDEX(ADUM(1:),'nacl').NE.0) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          CHREF(205) = 'XMLS'
          IREF(NV) = 205
        ELSEIF( INDEX(ADUM(1:),'effective').NE.0 .AND.
     &    INDEX(ADUM(1:),'trapped').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IREF(NV) = 19
        ELSEIF( INDEX(ADUM(1:),'diffusive').NE.0 .AND.
     &    INDEX(ADUM(1:),'porosity').NE.0 ) THEN
          IREF(NV) = 20
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 21
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 22
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 23
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 113
          CHREF(113) = 'XGN'
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 163
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 223
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 166
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0 ) .AND.
     &    (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 222
          CHREF(222) = 'XNN'
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 24
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 25
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 26
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 152
          CHREF(152) = 'XLN'
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydraulic').NE.0 .AND.
     &    INDEX(ADUM(1:),'head').NE.0 ) THEN
          IREF(NV) = 27
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydraulic').NE.0 .AND.
     &    INDEX(ADUM(1:),'head').NE.0 ) THEN
          IREF(NV) = 28
        ELSEIF ( INDEX(ADUM(1:),'rock').NE.0 .AND.
     &    INDEX(ADUM(1:),'soil').NE.0 .AND.
     &    INDEX(ADUM(1:),'type').NE.0 ) THEN
          IREF(NV) = 30
        ELSEIF( INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 .AND.
     &    INDEX(ADUM(1:),'relative').NE.0 .AND.
     &    INDEX(ADUM(1:),'permeability').NE.0 ) THEN
          IREF(NV) = 33
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'relative').NE.0 .AND.
     &    INDEX(ADUM(1:),'permeability').NE.0 ) THEN
          IREF(NV) = 31
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'relative').NE.0 .AND.
     &    INDEX(ADUM(1:),'permeability').NE.0 ) THEN
          IREF(NV) = 32
        ELSEIF( INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 .AND.
     &    INDEX(ADUM(1:),'density').NE.0 ) THEN
          IREF(NV) = 36
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'density').NE.0 ) THEN
          IREF(NV) = 34
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'density').NE.0 ) THEN
          IREF(NV) = 35
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 38
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 39
        ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          IREF(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          IREF(NV) = 41
        ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          IREF(NV) = 42
        ELSEIF( INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &    (INDEX(ADUM(1:),'salt').NE.0 .OR.
     &    INDEX(ADUM(1:),'nacl').NE.0) .AND.
     &    INDEX(ADUM(1:),'concentration').NE.0 ) THEN
          IREF(NV) = 47
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    (INDEX(ADUM(1:),'salt').NE.0 .OR.
     &    INDEX(ADUM(1:),'nacl').NE.0) .AND.
     &    INDEX(ADUM(1:),'concentration').NE.0 ) THEN
          IREF(NV) = 48
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'courant').NE.0 ) THEN
          ICRNT = 1
          IREF(NV) = 49
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    (INDEX(ADUM(1:),'salt').NE.0 .OR.
     &    INDEX(ADUM(1:),'nacl').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 50
        ELSEIF( INDEX(ADUM(1:),'xnc').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 87
        ELSEIF( INDEX(ADUM(1:),'ync').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 88
        ELSEIF( INDEX(ADUM(1:),'znc').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 89
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 51
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 52
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 53
        ELSEIF( INDEX(ADUM(1:),'xnc').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 90
        ELSEIF( INDEX(ADUM(1:),'ync').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 91
        ELSEIF( INDEX(ADUM(1:),'znc').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 92
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 54
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 55
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IREF(NV) = 56
        ELSEIF( INDEX(ADUM(1:),'xnc').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 96
        ELSEIF( INDEX(ADUM(1:),'ync').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 97
        ELSEIF( INDEX(ADUM(1:),'znc').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 98
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 60
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 61
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 62
        ELSEIF( INDEX(ADUM(1:),'xnc').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 67
        ELSEIF( INDEX(ADUM(1:),'ync').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 68
        ELSEIF( INDEX(ADUM(1:),'znc').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 69
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 64
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 65
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 66
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 70
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 71
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 72
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IREF(NV) = 73
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IREF(NV) = 74
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IREF(NV) = 75
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IREF(NV) = 76
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IREF(NV) = 77
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IREF(NV) = 78
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'courant').NE.0 ) THEN
          ICRNT = 1
          IREF(NV) = 79
        ELSEIF( INDEX(ADUM(1:),'ice').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IREF(NV) = 80
        ELSEIF( INDEX(ADUM(1:),'ice').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IREF(NV) = 81
        ELSEIF( INDEX(ADUM(1:),'ice').NE.0 .AND.
     &    INDEX(ADUM(1:),'density').NE.0 ) THEN
          IREF(NV) = 82
        ELSEIF( (INDEX(ADUM(1:),'nacl').NE.0 .OR.
     &    INDEX(ADUM(1:),'salt').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IREF(NV) = 110
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'vapor').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IREF(NV) = 128
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'vapor').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IREF(NV) = 129
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'vapor').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IREF(NV) = 130
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'scaling').NE.0 ) THEN
          IREF(NV) = 131
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IREF(NV) = 140
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IREF(NV) = 141
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IREF(NV) = 142
        ELSEIF( (INDEX(ADUM(1:),'nacl').NE.0 .OR.
     &    INDEX(ADUM(1:),'salt').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IREF(NV) = 147
        ELSEIF( (INDEX(ADUM(1:),'nacl').NE.0 .OR.
     &    INDEX(ADUM(1:),'salt').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          IREF(NV) = 148
        ELSEIF( INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 .AND.
     &    INDEX(ADUM(1:),'viscosity').NE.0 ) THEN
          IREF(NV) = 290
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'viscosity').NE.0 ) THEN
          IREF(NV) = 176
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'viscosity').NE.0 ) THEN
          IREF(NV) = 289
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 277
          CHREF(277) = 'DMHW'
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 278
          CHREF(278) = 'DMHA'
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 279
          CHREF(279) = 'DMHO'
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 280
          CHREF(280) = 'DMHN'
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'mobile').NE.0 ) THEN
          IREF(NV) = 281
          CHREF(281) = 'DMMW'
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'mobile').NE.0 ) THEN
          IREF(NV) = 282
          CHREF(282) = 'DMMA'
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mobile').NE.0 ) THEN
          IREF(NV) = 283
          CHREF(283) = 'DMMO'
        ELSEIF( INDEX(ADUM(1:),'differential').NE.0 .AND.
     &    INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'mobile').NE.0 ) THEN
          IREF(NV) = 284
          CHREF(284) = 'DMMN'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 256
          CHREF(256) = 'IMHW'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 257
          CHREF(257) = 'IMHA'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 258
          CHREF(258) = 'IMHO'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 171
          CHREF(171) = 'IMHN'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
          IREF(NV) = 195
          CHREF(195) = 'IMLA'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
          IREF(NV) = 259
          CHREF(259) = 'IMLO'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
          IREF(NV) = 239
          CHREF(239) = 'IMLN'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IREF(NV) = 197
          CHREF(197) = 'IMGW'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IREF(NV) = 198
          CHREF(198) = 'IMGA'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IREF(NV) = 260
          CHREF(260) = 'IMGO'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IREF(NV) = 240
          CHREF(240) = 'IMGN'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 190
          CHREF(190) = 'IMN'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 191
          CHREF(191) = 'IMW'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 192
          CHREF(192) = 'IMA'
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IREF(NV) = 193
          CHREF(193) = 'IMO'
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 ) THEN
          IREF(NV) = 250
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 ) THEN
          IREF(NV) = 251
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 ) THEN
          IREF(NV) = 252
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 ) THEN
          IREF(NV) = 241
          CHREF(241) = 'XHN'
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'density').NE.0 ) THEN
          IREF(NV) = 253
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IREF(NV) = 254
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IREF(NV) = 255
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IREF(NV) = 264
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 265
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 266
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IREF(NV) = 267
        ELSEIF( INDEX(ADUM(1:),'integrated').NE.0 .AND.
     &    INDEX(ADUM(1:),'energy').NE.0 ) THEN
          IREF(NV) = 378
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
     &    INDEX(ADUM(1:),'volum').NE.0 .AND.
     &    INDEX(ADUM(1:),'refer').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 372
          CHREF(372) = 'EPSV-R'
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'volum').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 200
          ENDIF
          IREF(NV) = 373
          CHREF(373) = 'EPSV'
!        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
!     &    INDEX(ADUM(1:),'yz').NE.0 ) THEN
!          IF( ISLC(50).EQ.0 ) THEN
!            NVC = NVC - 1
!            GOTO 200
!          ENDIF
!          IREF(NV) = 372
!          CHREF(372) = 'EPS-YZ'
!        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
!     &    INDEX(ADUM(1:),'xz').NE.0 ) THEN
!          IF( ISLC(50).EQ.0 ) THEN
!            NVC = NVC - 1
!            GOTO 200
!          ENDIF
!          IREF(NV) = 373
!          CHREF(373) = 'EPS-XZ'
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
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 1
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 2
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 3
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mol').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 5
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 6
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 8
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 9
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 10
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 11
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Reference Node Variable: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for duplicate reference node variables  ---
!
        IF( IREF(NV).NE.375 .AND. IREF(NV).NE.376 .AND.
     &    IREF(NV).NE.377 ) THEN
          DO NX = 1,NV-1
            IF( IREF(NV).EQ.IREF(NX) ) THEN
              INDX = 4
              CHMSG = 'Duplicate Reference Node Variable: '//ADUM
              CALL WRMSGS( INDX )
            ENDIF
          ENDDO
        ENDIF
!
!---    Reference node variable units  ---
!
        IDFLT = 1
        VARB = 'Reference Node Variable Unit'
        CALL RDCHR(ISTART,ICOMMA,NCU,CHDUM,UNREF(IREF(NV)))
        IF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          WRITE( IWR,'(2X,3A,2X,2A,I2,A)' ) ADUM(1:NCH),', ',
     &      UNREF(IREF(NV))(1:NCU),SOLNM(1:NCS),' Solute(',NSL,')'
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
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
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
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Plot File Variables'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NVPLOT)
      NVC = 0
      DO 400 NV = 1,NVPLOT
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Plot File Variable'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Plot File Variable: Solute Name'
          CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SOLNM)
          DO 310 NSL = 1,NSOLU
            IF( SOLNM.EQ.SOLUT(NSL) ) GOTO 320
  310     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute Name: '//SOLNM
          CALL WRMSGS( INDX )
          NVC = NVC -1
          GOTO 400
  320     CONTINUE
        ENDIF
        IF( INDEX(ADUM(1:),'final').NE.0 .AND.
     &    INDEX(ADUM(1:),'restart').NE.0 ) THEN
          ISLC(18) = 1
          IPLOT(NV) = 200
        ELSEIF( INDEX(ADUM(1:),'no').NE.0 .AND.
     &    INDEX(ADUM(1:),'restart').NE.0 ) THEN
          ISLC(18) = 2
          IPLOT(NV) = 200
        ELSEIF( (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'liq').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IPLOT(NV) = 3
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IPLOT(NV) = 1
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IPLOT(NV) = 2
        ELSEIF( INDEX(ADUM(1:),'critical').NE.0 .AND.
     &    INDEX(ADUM(1:),'temp').NE.0 ) THEN
          IPLOT(NV) = 206
        ELSEIF( INDEX(ADUM(1:),'cricondentherm').NE.0 .AND.
     &    INDEX(ADUM(1:),'temp').NE.0 ) THEN
          IPLOT(NV) = 213
        ELSEIF( INDEX(ADUM(1:),'temp').NE.0 ) THEN
          IPLOT(NV) = 4
        ELSEIF( INDEX(ADUM(1:),'phase').NE.0 .AND.
     &    INDEX(ADUM(1:),'cond').NE.0 ) THEN
          IPLOT(NV) = 5
        ELSEIF( (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'liq').NE.0 .AND.
     &    INDEX(ADUM(1:),'gauge').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 8
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'gauge').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IPLOT(NV) = 6
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'gauge').NE.0 .AND.
     &    INDEX(ADUM(1:),'press').NE.0 ) THEN
          IPLOT(NV) = 7
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'equil').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 386
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'formers').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 387
        ELSEIF( INDEX(ADUM(1:),'apparent').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IPLOT(NV) = 9
        ELSEIF( (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'liq').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IPLOT(NV) = 13
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IPLOT(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IPLOT(NV) = 12
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'moisture').NE.0 .AND.
     &    INDEX(ADUM(1:),'content').NE.0 ) THEN
          IPLOT(NV) = 15
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 ) THEN
          IPLOT(NV) = 165
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IPLOT(NV) = 294
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IPLOT(NV) = 295
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 ) THEN
          IPLOT(NV) = 221
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IPLOT(NV) = 379
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IPLOT(NV) = 380
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 ) THEN
          IPLOT(NV) = 162
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IPLOT(NV) = 381
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IPLOT(NV) = 382
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IPLOT(NV) = 383
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IPLOT(NV) = 384
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'total').NE.0 .AND.
     &    INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'former').NE.0 ) THEN
          IPLOT(NV) = 385
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 204
        ELSEIF( (INDEX(ADUM(1:),'salt').NE.0 .OR.
     &    INDEX(ADUM(1:),'nacl').NE.0) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 205
        ELSEIF( INDEX(ADUM(1:),'effective').NE.0 .AND.
     &    INDEX(ADUM(1:),'trapped').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IPLOT(NV) = 19
        ELSEIF( INDEX(ADUM(1:),'diffusive').NE.0 .AND.
     &    INDEX(ADUM(1:),'porosity').NE.0 ) THEN
          IPLOT(NV) = 20
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 21
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 22
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 113
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 23
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 163
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 223
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 166
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0 ) .AND.
     &    (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 222
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 24
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 25
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 26
        ELSEIF( (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 152
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydraulic').NE.0 .AND.
     &    INDEX(ADUM(1:),'head').NE.0 ) THEN
          IPLOT(NV) = 27
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydraulic').NE.0 .AND.
     &    INDEX(ADUM(1:),'head').NE.0 ) THEN
          IPLOT(NV) = 28
        ELSEIF ( INDEX(ADUM(1:),'rock').NE.0 .AND.
     &    INDEX(ADUM(1:),'soil').NE.0 .AND.
     &    INDEX(ADUM(1:),'type').NE.0 ) THEN
          IPLOT(NV) = 30
        ELSEIF( (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 .AND.
     &    INDEX(ADUM(1:),'relative').NE.0 .AND.
     &    INDEX(ADUM(1:),'permeability').NE.0 ) THEN
          IPLOT(NV) = 33
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'relative').NE.0 .AND.
     &    INDEX(ADUM(1:),'permeability').NE.0 ) THEN
          IPLOT(NV) = 31
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'relative').NE.0 .AND.
     &    INDEX(ADUM(1:),'permeability').NE.0 ) THEN
          IPLOT(NV) = 32
        ELSEIF( (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 .AND.
     &    INDEX(ADUM(1:),'density').NE.0 ) THEN
          IPLOT(NV) = 36
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'density').NE.0 ) THEN
          IPLOT(NV) = 34
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'density').NE.0 ) THEN
          IPLOT(NV) = 35
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IPLOT(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IPLOT(NV) = 38
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IPLOT(NV) = 39
        ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          IPLOT(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          IPLOT(NV) = 41
        ELSEIF( INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          IPLOT(NV) = 42
        ELSEIF( INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &    (INDEX(ADUM(1:),'salt').NE.0 .OR.
     &    INDEX(ADUM(1:),'nacl').NE.0) .AND.
     &    INDEX(ADUM(1:),'concentration').NE.0 ) THEN
          IPLOT(NV) = 47
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    (INDEX(ADUM(1:),'salt').NE.0 .OR.
     &    INDEX(ADUM(1:),'nacl').NE.0) .AND.
     &    INDEX(ADUM(1:),'concentration').NE.0 ) THEN
          IPLOT(NV) = 48
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'courant').NE.0 ) THEN
          ICRNT = 1
          IPLOT(NV) = 49
        ELSEIF( INDEX(ADUM(1:),'total').NE.0 .AND.
     &    (INDEX(ADUM(1:),'salt').NE.0 .OR.
     &    INDEX(ADUM(1:),'nacl').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          IPLOT(NV) = 50
        ELSEIF( INDEX(ADUM(1:),'xnc').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 87
        ELSEIF( INDEX(ADUM(1:),'ync').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 88
        ELSEIF( INDEX(ADUM(1:),'znc').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 89
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 51
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 52
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 53
        ELSEIF( INDEX(ADUM(1:),'xnc').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 90
        ELSEIF( INDEX(ADUM(1:),'ync').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 91
        ELSEIF( INDEX(ADUM(1:),'znc').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 92
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 54
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 55
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 ) THEN
          IPLOT(NV) = 56
        ELSEIF( INDEX(ADUM(1:),'xnc').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 96
        ELSEIF( INDEX(ADUM(1:),'ync').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 97
        ELSEIF( INDEX(ADUM(1:),'znc').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 98
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 60
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 61
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'heat').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 62
        ELSEIF( INDEX(ADUM(1:),'xnc').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 67
        ELSEIF( INDEX(ADUM(1:),'ync').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 68
        ELSEIF( INDEX(ADUM(1:),'znc').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 69
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 64
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 65
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 66
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 70
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 71
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 72
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 73
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 74
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 75
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 76
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 77
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 78
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'courant').NE.0 ) THEN
          ICRNT = 1
          IPLOT(NV) = 79
        ELSEIF( INDEX(ADUM(1:),'ice').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 80
        ELSEIF( INDEX(ADUM(1:),'ice').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IPLOT(NV) = 81
        ELSEIF( INDEX(ADUM(1:),'ice').NE.0 .AND.
     &    INDEX(ADUM(1:),'density').NE.0 ) THEN
          IPLOT(NV) = 82
        ELSEIF( (INDEX(ADUM(1:),'nacl').NE.0 .OR.
     &    INDEX(ADUM(1:),'salt').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 ) THEN
          IPLOT(NV) = 110
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'vapor').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 128
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'vapor').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 129
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'vapor').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 130
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'scaling').NE.0 ) THEN
          IPLOT(NV) = 131
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IPLOT(NV) = 140
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IPLOT(NV) = 141
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IPLOT(NV) = 142
        ELSEIF( (INDEX(ADUM(1:),'nacl').NE.0 .OR.
     &    INDEX(ADUM(1:),'salt').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'rate').NE.0 ) THEN
          IPLOT(NV) = 147
        ELSEIF( (INDEX(ADUM(1:),'nacl').NE.0 .OR.
     &    INDEX(ADUM(1:),'salt').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 .AND.
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          IPLOT(NV) = 148
        ELSEIF( (INDEX(ADUM(1:),'nonaqueous').NE.0 .OR.
     &    INDEX(ADUM(1:),'non-aqueous').NE.0) .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 .AND.
     &    INDEX(ADUM(1:),'viscosity').NE.0 ) THEN
          IPLOT(NV) = 290
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'viscosity').NE.0 ) THEN
          IPLOT(NV) = 176
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'viscosity').NE.0 ) THEN
          IPLOT(NV) = 289
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 ) THEN
          IPLOT(NV) = 250
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 ) THEN
          IPLOT(NV) = 251
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0 ) .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 ) THEN
          IPLOT(NV) = 252
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND. 
     &    (INDEX(ADUM(1:),'n2').NE.0 .OR.
     &    INDEX(ADUM(1:),'nitrogen').NE.0) .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    INDEX(ADUM(1:),'fraction').NE.0 ) THEN
          IPLOT(NV) = 241
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'density').NE.0 ) THEN
          IPLOT(NV) = 253
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.00 .AND.
     &    INDEX(ADUM(1:),'equil').NE.0 ) THEN
          IPLOT(NV) = 256
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IPLOT(NV) = 254
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'pressure').NE.0 ) THEN
          IPLOT(NV) = 255
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0 .AND.
     &    INDEX(ADUM(1:),'saturation').NE.0 ) THEN
          IPLOT(NV) = 264
        ELSEIF( (INDEX(ADUM(1:),'h2o').NE.0 .OR.
     &    INDEX(ADUM(1:),'water').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IPLOT(NV) = 265
        ELSEIF( (INDEX(ADUM(1:),'o2').NE.0 .OR.
     &    INDEX(ADUM(1:),'carbon').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IPLOT(NV) = 266
        ELSEIF( (INDEX(ADUM(1:),'ch4').NE.0 .OR.
     &    INDEX(ADUM(1:),'methane').NE.0) .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 .AND.
     &    INDEX(ADUM(1:),'frac').NE.0 .AND.
     &    INDEX(ADUM(1:),'hydrate').NE.0 ) THEN
          IPLOT(NV) = 267
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
     &    INDEX(ADUM(1:),'volum').NE.0 .AND.
     &    INDEX(ADUM(1:),'refer').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 372
        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
     &    INDEX(ADUM(1:),'volum').NE.0 ) THEN
          IF( ISLC(50).EQ.0 ) THEN
            NVC = NVC - 1
            GOTO 400
          ENDIF
          IPLOT(NV) = 373
!        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
!     &    INDEX(ADUM(1:),'yz').NE.0 ) THEN
!          IF( ISLC(50).EQ.0 ) THEN
!            NVC = NVC - 1
!            GOTO 400
!          ENDIF
!          IPLOT(NV) = 372
!        ELSEIF( INDEX(ADUM(1:),'strain').NE.0 .AND.
!     &    INDEX(ADUM(1:),'xz').NE.0 ) THEN
!          IF( ISLC(50).EQ.0 ) THEN
!            NVC = NVC - 1
!            GOTO 400
!          ENDIF
!          IPLOT(NV) = 373
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
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'volumetric').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 1
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 2
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 3
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'mol').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 5
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'gas').NE.0 .AND.
     &    INDEX(ADUM(1:),'mole').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 6
        ELSEIF( INDEX(ADUM(1:),'x').NE.0 .AND.
     &    INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 8
        ELSEIF( INDEX(ADUM(1:),'y').NE.0 .AND.
     &    INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 9
        ELSEIF( INDEX(ADUM(1:),'z').NE.0 .AND.
     &    INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 10
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 .AND.
     &    INDEX(ADUM(1:),'source').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 11
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Plot File Variable: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for duplicate plot file variables  ---
!
        IF( IPLOT(NV).NE.375 .AND. IPLOT(NV).NE.376 .AND.
     &    IPLOT(NV).NE.377 ) THEN
          DO NX = 1,NV-1
            IF( IPLOT(NV).EQ.IPLOT(NX) ) THEN
              INDX = 4
              CHMSG = 'Duplicate Plot File Variable: '//ADUM
              CALL WRMSGS( INDX )
            ENDIF
          ENDDO
        ENDIF
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDOU_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSF_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Read input file surface flux information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE OUTPU
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
      CHARACTER*64 ADUM,BDUM
      CHARACTER*512 CHDUM,CHDUMX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSF_HYD_O2_CH4_N2'
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
      NC = 0
      DO 100 NS = 1, NSF
        IF( NS.NE.1 ) WRITE(IWR, '(/)')
        CALL RDINPL( CHDUM )
        CHDUMX = CHDUM
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---  Check for specified surface flux filename  ---
!
        NSFF = 0
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
          IF( NSFGP.GT.LSF ) THEN
            INDX = 4
            CHMSG = 'Number of Surface Flux Files > LSF'
            CALL WRMSGS( INDX )
          ENDIF
          NC = 0
          FNSF(NSFGP) = ADUM
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
        ENDIF
        NC = NC + 1
        IF( (NSFF-NC).GE.0 )THEN
          ISFF(NS) = NSFGP
        ELSE
          ISFF(NS) = 1
          ISFGP(1) = ISFGP(1) + 1
        ENDIF
!
!---  Read surface flux type  ---
!
        VARB = 'Surface Flux Type'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'heat flux').NE.0) THEN
          ISFT(NS) = 1
          WRITE(IWR,'(A)') 'Heat Flux Surface'
          UNSF(1,NS) = 'w'
          UNSF(2,NS) = 'j'
        ELSEIF( INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'liquid').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'volum').NE.0) THEN
            ISFT(NS) = 4
            WRITE(IWR,'(A)') 'Nonaqueous-Liquid Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ELSEIF( INDEX(ADUM(1:),'carbon').NE.0 .OR.
     &      INDEX(ADUM(1:),'o2').NE.0 ) THEN
            ISFT(NS) = 55
            WRITE(IWR,'(A)') 'Nonaqueous-Liquid O2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'methane').NE.0 .OR.
     &      INDEX(ADUM(1:),'ch4').NE.0 ) THEN
            ISFT(NS) = 56
            WRITE(IWR,'(A)') 'Nonaqueous-Liquid CH4 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'nitrogen').NE.0 .OR.
     &      INDEX(ADUM(1:),'n2').NE.0 ) THEN
            ISFT(NS) = 61
            WRITE(IWR,'(A)') 'Nonaqueous-Liquid N2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'water').NE.0 .OR.
     &      INDEX(ADUM(1:),'h2o').NE.0 ) THEN
            ISFT(NS) = 57
            WRITE(IWR,'(A)') 'Nonaqueous-Liquid Water Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 7
            WRITE(IWR,'(A)') 'Nonaqueous-Liquid Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0) THEN
          IF( INDEX(ADUM(1:),'volum').NE.0) THEN
            ISFT(NS) = 2
            WRITE(IWR,'(A)') 'Aqueous Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ELSEIF( INDEX(ADUM(1:),'carbon').NE.0 .OR.
     &      INDEX(ADUM(1:),'o2').NE.0 ) THEN
            ISFT(NS) = 29
            WRITE(IWR,'(A)') 'Aqueous O2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'methane').NE.0 .OR.
     &      INDEX(ADUM(1:),'ch4').NE.0 ) THEN
            ISFT(NS) = 12
            WRITE(IWR,'(A)') 'Aqueous CH4 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'nitrogen').NE.0 .OR.
     &      INDEX(ADUM(1:),'n2').NE.0 ) THEN
            ISFT(NS) = 59
            WRITE(IWR,'(A)') 'Aqueous N2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'water').NE.0 .OR.
     &      INDEX(ADUM(1:),'h2o').NE.0 ) THEN
            ISFT(NS) = 43
            WRITE(IWR,'(A)') 'Aqueous Water Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 5
            WRITE(IWR,'(A)') 'Aqueous Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'volum').NE.0) THEN
            ISFT(NS) = 3
            WRITE(IWR,'(A)') 'Gas Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ELSEIF( INDEX(ADUM(1:),'carbon').NE.0 .OR.
     &      INDEX(ADUM(1:),'o2').NE.0 ) THEN
            ISFT(NS) = 28
            WRITE(IWR,'(A)') 'Gas O2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'methane').NE.0 .OR.
     &      INDEX(ADUM(1:),'ch4').NE.0 ) THEN
            ISFT(NS) = 11
            WRITE(IWR,'(A)') 'Gas CH4 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'nitrogen').NE.0 .OR.
     &      INDEX(ADUM(1:),'n2').NE.0 ) THEN
            ISFT(NS) = 60
            WRITE(IWR,'(A)') 'Gas N2 Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSEIF( INDEX(ADUM(1:),'water').NE.0 .OR.
     &      INDEX(ADUM(1:),'h2o').NE.0 ) THEN
            ISFT(NS) = 44
            WRITE(IWR,'(A)') 'Gas Water Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 6
            WRITE(IWR,'(A)') 'Gas Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'salt').NE.0) THEN
          ISFT(NS) = 8
          WRITE(IWR,'(A)') 'Salt-Mass Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'carbon').NE.0 .OR.
     &    INDEX(ADUM(1:),'o2').NE.0 ) THEN
          ISFT(NS) = 30
          WRITE(IWR,'(A)') 'O2-Mass Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'methane').NE.0 .OR.
     &    INDEX(ADUM(1:),'ch4').NE.0 ) THEN
          ISFT(NS) = 13
          WRITE(IWR,'(A)') 'CH4-Mass Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'nitrogen').NE.0 .OR.
     &    INDEX(ADUM(1:),'n2').NE.0 ) THEN
          ISFT(NS) = 58
          WRITE(IWR,'(A)') 'N2-Mass Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'water').NE.0 .OR.
     &    INDEX(ADUM(1:),'h2o').NE.0 ) THEN
          ISFT(NS) = 45
          WRITE(IWR,'(A)') 'Water-Mass Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
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
              GOTO 20
            ENDIF
   10     CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Solute Name: '//BDUM
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
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(2A,$)') VARB(1:IVR),': '
        ISFSN(NS) = 0
        IF( INDEX(ADUM(1:),'surface normal').NE.0 )  ISFSN(NS) = 1
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
        ENDIF
!
!---  Read surface domain  ---
!
        VARB = 'Surface Flux Domain'
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(1,NS))
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(2,NS))
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(3,NS))
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(4,NS))
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(5,NS))
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(6,NS))
        IF( ISFC(1,NS).LT.1 .OR. ISFC(1,NS).GT.IFLD .OR.
     &    ISFC(2,NS).LT.1 .OR. ISFC(2,NS).GT.IFLD .OR.
     &    ISFC(1,NS).GT.ISFC(2,NS) ) THEN
          INDX = 4
          CHMSG = 'Illegal Surface Flux Domain: I Indices'
          CALL WRMSGS( INDX )
        ENDIF
        IF( ISFC(3,NS).LT.1 .OR. ISFC(3,NS).GT.JFLD .OR.
     &    ISFC(4,NS).LT.1 .OR. ISFC(4,NS).GT.JFLD .OR.
     &    ISFC(3,NS).GT.ISFC(4,NS) ) THEN
          INDX = 4
          CHMSG = 'Illegal Surface Flux Domain: J Indices'
          CALL WRMSGS( INDX )
        ENDIF
        IF( ISFC(5,NS).LT.1 .OR. ISFC(5,NS).GT.KFLD .OR.
     &    ISFC(6,NS).LT.1 .OR. ISFC(6,NS).GT.KFLD .OR.
     &    ISFC(5,NS).GT.ISFC(6,NS) ) THEN
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
        WRITE(IWR,'(2X,2(A,I6))') 'I = ',ISFC(1,NS),' to ',ISFC(2,NS)
        WRITE(IWR,'(2X,2(A,I6))') 'J = ',ISFC(3,NS),' to ',ISFC(4,NS)
        WRITE(IWR,'(2X,2(A,I6))') 'K = ',ISFC(5,NS),' to ',ISFC(6,NS)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSF_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSP_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Read input file for rock/soil saturation function information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NAPL
      USE HYST
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSP_HYD_O2_CH4_N2'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Rock/Soil Saturation Function Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read new line  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
!
!---  Read fluid pair interfacial tensions  ---
!
      VARB = 'Gas-Aqueous Interfacial Tension'
      SCHR(5,1) = -1.D+0
      IDFLT = 1
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(5,1))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &  ': ',SCHR(5,1)
      INDX = 0
      IUNKG = 1
      IUNS = -2
      CALL RDUNIT(UNTS,SCHR(5,1),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(5,1),', N/m)'
      VARB ='Nonaqueous Liquid-Aqueous Interfacial Tension'
      SCHR(6,1) = -1.D+0
      IDFLT = 1
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(6,1))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &  ': ',SCHR(6,1)
      INDX = 0
      IUNKG = 1
      IUNS = -2
      CALL RDUNIT(UNTS,SCHR(6,1),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(6,1),', N/m)'
      SFCSF(3) = SCHR(6,1)
      VARB ='Gas-Nonaqueous Liquid Interfacial Tension'
      SCHR(7,1) = -1.D+0
      IDFLT = 1
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(7,1))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &  ': ',SCHR(7,1)
      INDX = 0
      IUNKG = 1
      IUNS = -2
      CALL RDUNIT(UNTS,SCHR(7,1),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(7,1),', N/m)'
      VARB ='Hydrate-Aqueous Interfacial Tension'
      SCHR(9,1) = 0.0267D+0
      IDFLT = 1
      CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(9,1))
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &  ': ',SCHR(9,1)
      INDX = 0
      IUNKG = 1
      IUNS = -2
      CALL RDUNIT(UNTS,SCHR(9,1),INDX)
      WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(9,1),', N/m)'
!
!---  Check for under specified interfacial tension systems  ---
!
      IF( SCHR(5,1).LT.ZERO .AND. SCHR(6,1).LT.ZERO ) THEN
        INDX = 4
        CHMSG = 'Under Specified Interfacial Tensions (Two Minimum)'
        CALL WRMSGS( INDX )
      ELSEIF( SCHR(6,1).LT.ZERO .AND. SCHR(7,1).LT.ZERO ) THEN
        INDX = 4
        CHMSG = 'Under Specified Interfacial Tensions (Two Minimum)'
        CALL WRMSGS( INDX )
      ELSEIF( SCHR(7,1).LT.ZERO .AND. SCHR(5,1).LT.ZERO ) THEN
        INDX = 4
        CHMSG = 'Under Specified Interfacial Tensions (Two Minimum)'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Assign defaulted interfacial tensions  ---
!
      IF( SCHR(5,1).LT.ZERO ) THEN
        SCHR(5,1) = SCHR(6,1) + SCHR(7,1)
      ELSEIF( SCHR(6,1).LT.ZERO ) THEN
        SCHR(6,1) = SCHR(5,1) - SCHR(7,1)
      ELSEIF( SCHR(7,1).LT.ZERO ) THEN
        SCHR(7,1) = SCHR(5,1) - SCHR(6,1)
      ENDIF
!
!---  Check for indistinguishable interfacial tensions  ---
!
      IF( ABS(SCHR(5,1)-SCHR(6,1)).LT.EPSL ) THEN
        INDX = 4
        CHMSG = 'Gas-Aqueous Surface Tension = ' //
     &    'Nonaqueous Liquid-Aqueous Surface Tension'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Check for zero gas-NAPL interfacial tension  ---
!
      IF( SCHR(7,1).LT.EPSL ) THEN
        INDX = 4
        CHMSG = 'Zero Gas-Nonaqueous Liquid Interfacial Tension'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write spreading coefficient  ---
!
      C_SPRD = MIN( (SCHR(5,1)-SCHR(6,1)-SCHR(7,1)),0.D+0 )
      WRITE(IWR,'(2X,A,1PE11.4)') 
     &  'Nonaqueous Liquid Spreading Coefficient: ',C_SPRD
      IF( C_SPRD.LT.-EPSL ) THEN
        CA_NL = ACOS( ((SCHR(6,1)**2) - (SCHR(7,1)**2) + (SCHR(5,1)**2))
     &    /(2.D+0*SCHR(5,1)*SCHR(6,1)) )
        CA_GN = ACOS( ((SCHR(7,1)**2) - (SCHR(6,1)**2) + (SCHR(5,1)**2))
     &    /(2.D+0*SCHR(5,1)*SCHR(7,1)) )
!
!---    Write gas-NAPL and NAPL-aqueous contact angles  ---
!

        WRITE(IWR,'(2X,A,1PE11.4)') 'Gas-Nonaqueous Liquid Contact ' //
     &    'Angle, deg: ',(CA_GN*1.8D+2/GPI)
        WRITE(IWR,'(2X,A,1PE11.4)') 'Nonaqueous Liquid-Aqueous ' //
     &    'Contact Angle, deg: ',(CA_NL*1.8D+2/GPI)
        CA_GN = COS(CA_GN)
        CA_NL = COS(CA_NL)
      ENDIF
!
!---  Assign interfacial tensions and hydrate formation-dissociation
!     rate to rock/soil characteristic parameters  ---
!
      DO 10 N = 2,NROCK
        SCHR(5,N) = SCHR(5,1)
        SCHR(6,N) = SCHR(6,1)
        SCHR(7,N) = SCHR(7,1)
        SCHR(9,N) = SCHR(9,1)
  10  CONTINUE
!
!---  Loop over the rock/soil saturation information lines  ---
!
      N = 0
      IJK = 0
      IROCK = 1
   20 CONTINUE
      IF( N.GE.NROCK .OR. IJK.GT.0 ) GOTO 500
      ISTART = 1
      VARB = 'Saturation Function: Rock Name: '
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
        GOTO 20
      ENDIF
      INDX = 2
      CHMSG = 'Unrecognized Rock/Soil Type: ' // RDUM(1:NCH)
      CALL WRMSGS( INDX )
      GOTO 20
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
      N = N + 1
  220 CONTINUE
!
!---  Read saturation/capillary pressure function  ---
!
      VARB = 'Saturation Function Type: '
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
      ISMX = 0
      IF( INDEX(ADUM(1:),'extended').NE.0 ) ISMX = 1
      IF( INDEX(ADUM(1:),'webb').NE.0 ) ISMX = 2
      IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
        IF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'barrier').NE.0 ) THEN
          ISCHRX = 11
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'porosity').NE.0 ) THEN
          ISCHRX = 13
        ELSE
          ISCHRX = 1
        ENDIF
      ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &    INDEX(ADUM(1:),'corey').NE.0 ) THEN
        IF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'barrier').NE.0 ) THEN
          ISCHRX = 12
        ELSEIF( INDEX(ADUM(1:),'hydrate').NE.0 .AND.
     &    INDEX(ADUM(1:),'porosity').NE.0 ) THEN
          ISCHRX = 14
        ELSE
          ISCHRX = 2
        ENDIF
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Saturation Function: '//ADUM
        CALL WRMSGS( INDX )
      ENDIF
      IF( IJK.GT.0 ) THEN
        DO 230 N = 1,NFLD
          ISCHR(IZ(N)) = ISCHRX
          ISM(IZ(N)) = ISMX
  230   CONTINUE
      ELSE
        ISCHR(IROCK) = ISCHRX
        ISM(IROCK) = ISMX
      ENDIF
!
!---  van Genuchten Function  ---
!
      IF( ISCHRX.EQ.1 .OR. ISCHRX.EQ.11 .OR. ISCHRX.EQ.13 ) THEN
        IF( ISCHRX.EQ.1 ) THEN
          WRITE(IWR,'(A)') 'van Genuchten Function'
        ELSEIF( ISCHRX.EQ.11 ) THEN
         WRITE(IWR,'(A)')'van Genuchten  Function w/ Hydrate Barrier'
        ELSEIF( ISCHRX.EQ.13 ) THEN
         WRITE(IWR,'(A)')'van Genuchten  Function w/ Hydrate-Porosity'
        ENDIF
        VARB = 'van Genuchten (alpha)'
        IF( IJK.GT.0 ) THEN
          INDX = 1
          LNDX = LSCHR
          UNTS = '1/m'
          IUNM = -1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(1,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHR(1,IROCK)
          INDX = 0
          IUNM = -1
          CALL RDUNIT(UNTS,SCHR(1,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(1,IROCK),', 1/m)'
        ENDIF
        VARB = 'van Genuchten (n): '
        IF( IJK.GT.0 ) THEN
          INDX = 3
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(3,IROCK))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(3,IROCK)
        ENDIF
        VARB = 'van Genuchten (residual saturation): '
        IF( IJK.GT.0 ) THEN
          INDX = 4
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(4,IROCK))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(4,IROCK)
        ENDIF
        VARB = 'van Genuchten (m)'
        IF( IJK.GT.0 ) THEN
          INDX = 14
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(14,IROCK))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(14,IROCK)
        ENDIF
        VARB = 'van Genuchten (Reference Interfacial Tension)'
        IDFLT = 1
        IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
          TX = 20.D+0
          XLSX = 0.D+0
          CALL SFT_L( TX,XLSX,SCHR(8,IROCK) )
        ELSE
          SCHR(8,IROCK) = SCHR(5,IROCK)
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 8
          LNDX = LSCHR
          UNTS = 'n/m'
          IUNKG = 1
          IUNS = -2
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(8,IROCK))
          UNTS = 'n/m'
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',SCHR(8,IROCK)
          INDX = 0
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,SCHR(8,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(8,IROCK),', N/m)'
        ENDIF
!
!---    Hydrate barrier model  ---
!
        IF( ISCHRX.EQ.11 ) THEN
          VARB = 'Hydrate Barrier Model W0 Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 10
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(10,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(10,IROCK)
          ENDIF
          VARB = 'Hydrate Barrier Model W1 Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 11
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(11,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(11,IROCK)
          ENDIF
          VARB = 'Hydrate Barrier Model W2 Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 12
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(12,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(12,IROCK)
          ENDIF
          VARB = 'Hydrate Barrier Model W3 Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 13
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(13,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(13,IROCK)
          ENDIF
        ENDIF
!
!---    Hydrate-porosity model  ---
!
        IF( ISCHRX.EQ.13 ) THEN
          VARB = 'Hydrate-Porosity Model: Hydrate m Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 10
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(10,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(10,IROCK)
          ENDIF
          VARB = 'Hydrate-Porosity Model: Porosity a Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 11
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(11,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(11,IROCK)
          ENDIF
          VARB = 'Hydrate-Porosity Model: Initial Porosity Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 12
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(12,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(12,IROCK)
          ENDIF
        ENDIF
!
!---  Brooks and Corey Function  ---
!
      ELSEIF( ISCHRX.EQ.2 .OR. ISCHRX.EQ.12 .OR. ISCHRX.EQ.14 ) THEN
        IF( ISCHRX.EQ.2 ) THEN
          WRITE(IWR,'(A)') 'Brooks and Corey Function'
        ELSEIF( ISCHRX.EQ.12 ) THEN
         WRITE(IWR,'(A)')'Brooks and Corey Function ' // 
     &     'w/ Hydrate Barrier'
        ELSEIF( ISCHRX.EQ.14 ) THEN
         WRITE(IWR,'(A)')'Brooks and Corey Function ' // 
     &     'w/ Hydrate-Porosity'
        ENDIF
        VARB = 'Brooks and Corey (psi): '
        IF( IJK.GT.0 ) THEN
          INDX = 1
          LNDX = LSCHR
          UNTS = 'm'
          IUNM = 1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(2,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHR(2,IROCK)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SCHR(2,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(2,IROCK),', m)'
        ENDIF
        VARB = 'Brooks and Corey (lambda): '
        IF( IJK.GT.0 ) THEN
          INDX = 3
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(3,IROCK))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(3,IROCK)
        ENDIF
        VARB = 'Brooks and Corey (residual saturation): '
        IF( IJK.GT.0 ) THEN
          INDX = 4
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(4,IROCK))
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),SCHR(4,IROCK)
        ENDIF
        VARB = 'Brooks and Corey (Reference Interfacial Tension)'
        IDFLT = 1
        IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 ) THEN
          TX = 20.D+0
          XLSX = 0.D+0
          CALL SFT_L( TX,XLSX,SCHR(8,IROCK) )
        ELSE
          SCHR(8,IROCK) = SCHR(5,IROCK)
        ENDIF
        IF( IJK.GT.0 ) THEN
          INDX = 8
          LNDX = LSCHR
          UNTS = 'n/m'
          IUNKG = 1
          IUNS = -2
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(8,IROCK))
          UNTS = 'n/m'
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',SCHR(8,IROCK)
          INDX = 0
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,SCHR(8,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(8,IROCK),', N/m)'
        ENDIF
!
!---    Hydrate barrier model  ---
!
        IF( ISCHRX.EQ.12 ) THEN
          VARB = 'Hydrate Barrier Model W0 Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 10
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(10,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(10,IROCK)
          ENDIF
          VARB = 'Hydrate Barrier Model W1 Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 11
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(11,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(11,IROCK)
          ENDIF
          VARB = 'Hydrate Barrier Model W2 Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 12
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(12,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(12,IROCK)
          ENDIF
          VARB = 'Hydrate Barrier Model W3 Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 13
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(13,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(13,IROCK)
          ENDIF
        ENDIF
!
!---    Hydrate-porosity model  ---
!
        IF( ISCHRX.EQ.14 ) THEN
          VARB = 'Hydrate-Porosity Model: Hydrate m Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 10
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(10,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(10,IROCK)
          ENDIF
          VARB = 'Hydrate-Porosity Model: Porosity a Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 11
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(11,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(11,IROCK)
          ENDIF
          VARB = 'Hydrate-Porosity Model: Initial Porosity Coefficient'
          IF( IJK.GT.0 ) THEN
            INDX = 12
            LNDX = LSCHR
            UNTS = 'null'
            CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
          ELSE
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(12,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &        SCHR(12,IROCK)
          ENDIF
        ENDIF
      ENDIF
!
!---  Bound-water aqueous saturation  ---
!
      VARB = 'Bound-Water Aqueous Saturation'
      IF( IJK.GT.0 ) THEN
        INDX = 15
        LNDX = LSCHR
        UNTS = 'null'
        CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
      ELSE
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(15,IROCK))
        WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',
     &    SCHR(15,IROCK)
      ENDIF
!
!---  Loop over remaining rock/soils within scaling group  ---
!
      IF( ISLC(19).EQ.1 .AND. IROCK.LT.NROCK ) THEN
        DO 490 M = IROCK+1,NROCK
          IF( ISCALE(M).EQ.ISGRP ) THEN
            N = N+1
            ISM(M) = ISM(IROCK)
            ISCHR(M) = ISCHR(IROCK)
            DO 480 L = 1,LSCHR
              SCHR(L,M) = SCHR(L,IROCK)
  480       CONTINUE
            DO 482 L = 1,2
              ISLTBL(L,M) = ISLTBL(L,IROCK)
  482       CONTINUE
          ENDIF
  490   CONTINUE
      ENDIF
!
!---  Read next rock/soil type or scaling group  ---
!
      IF( N.LT.NROCK ) WRITE(IWR,'(/)')
      GOTO 20
 500  CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSP_HYD_O2_CH4_N2 group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSR_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Read input file for source information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
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
      REAL*8 VAR(LSTM,8+LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSR_HYD_O2_CH4_N2'
      I1X = 0
      I2X = 0
      J1X = 0
      J2X = 0
      K1X = 0
      K2X = 0
!
!---  Write card information to ouput file  ---
!
      CARD = 'Source Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      NSR = 0
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Sources'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO 140 NS = 1, NLIN
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
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
        ELSEIF( INDEX(ADUM(1:),'nonaqueous').NE.0 .AND.
     &    INDEX(ADUM(1:),'fluid').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Nonaqueous-Fluid Mass Source'
          ISRTX = 4
        ELSEIF( INDEX(ADUM(1:),'fluid').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'ic').NE.0 .OR.
     &    INDEX(ADUM(1:),'initial').NE.0 ) ) THEN
          WRITE(IWR,'(2X,A)') 'Fluid-Mass Source'
          ISRTX = 5
        ELSEIF( INDEX(ADUM(1:),'fluid').NE.0 .AND.
     &    INDEX(ADUM(1:),'mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Fluid-Mass Source'
          ISRTX = 3
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
              GOTO 32
            ENDIF
   30     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute Name: '//ADUM
          CALL WRMSGS( INDX )
   32     CONTINUE
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Source Type: '//ADUM
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read source domain indices  ---
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
!---    Check for ill-defined source domains  ---
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
!---    Check for sources applied to inactive nodes  ---
!
        DO 50 K = K1X,K2X
        DO 50 J = J1X,J2X
        DO 50 I = I1X,I2X
          IF( IXP(ND(I,J,K)).EQ.0 ) THEN
            INDX = 7
            IMSG = ND(I,J,K)
            CHMSG = 'Source Applied to an Inactive Node: '
            CALL WRMSGS( INDX )
          ENDIF
   50   CONTINUE
!
!---    Read number of source times  ---
!
        VARB = 'Number of Source Times'
        CALL RDINT(ISTART,ICOMMA,CHDUM,ISRM(NS))
        IF( ISRM(NS).GT.LSTM ) THEN
          INDX = 5
          CHMSG = 'Number of Source Times > Parameter LSTM'
          CALL WRMSGS( INDX )
        ENDIF
        SRTMO = -SMALL
!
!---    Loop over number of source times  ---
!
        DO 100 NTM = 1,ISRM(NS)
          DO 70 M = 1,8+NSOLU
            VAR(NTM,M) = 0.D+0
   70     CONTINUE
!
!---      Read and write source values and units  ---
!
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
          VARB = 'Source Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',VAR(NTM,1)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR(NTM,1),INDX)
!
!---      Power Source  ---
!
          IF( ISRTX.EQ.1 ) THEN
            VARB = 'Source: Power'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 2
            IUNKG = 1
            IUNS = -3
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
!
!---      Power-Density Source  ---
!
          ELSEIF( ISRTX.EQ.2 ) THEN
            VARB = 'Source: Power Density'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -3
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
!
!---      Fluid Mass Sink  ---
!
          ELSEIF( ISRTX.EQ.3 ) THEN
            VARB = 'Source: Fluid Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            IF( VAR(NTM,4).GT.0.D+0 ) THEN
              INDX = 4
              CHMSG = 'Fluid Mass Rate > 0.0'
              CALL WRMSGS( INDX )
            ENDIF
!
!---      Nonaqueous Fluid Mass Source  ---
!
          ELSEIF( ISRTX.EQ.4 ) THEN
            VARB = 'Source: Nonaqueous Fluid Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            WRITE(IWR,'(2X,A,1PE11.4,A)') ' (',VAR(NTM,4),', kg/s)'
            VARB = 'Source: Temperature'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,2)
            INDX = 0
            IUNK = 1
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,2),', C)'
!
!---        O2 Mole Fraction of Nonaqueous Components ---
!
            VARB = 'Source: O2 Mole Fraction of Nonaqueous Components'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,5)
            IF( VAR(NTM,5).GT.1.D+0 .OR. VAR(NTM,5).LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Out of Range Value: ' // 
     &          'Source: O2 Mole Fraction of Nonaqueous Components'
              RLMSG = VAR(NTM,5)
              CALL WRMSGS( INDX )
            ENDIF
!
!---        CH4 Mole Fraction of Nonaqueous Components ---
!
            VARB = 'Source: CH4 Mole Fraction of Nonaqueous Components'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,6)
            IF( VAR(NTM,6).GT.1.D+0 .OR. VAR(NTM,6).LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Out of Range Value: ' // 
     &          'Source: CH4 Mole Fraction of Nonaqueous Components'
              RLMSG = VAR(NTM,6)
              CALL WRMSGS( INDX )
            ENDIF
            VARX = VAR(NTM,5)+VAR(NTM,6)
            IF( VARX.GT.1.D+0 .OR. VARX.LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Out of Range Value: ' // 
     &          'Gas: O2+CH4 Mole Fraction of Nonaqueous Components'
              RLMSG = VARX
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Water relative saturation  ---
!
            VARB = 'Source: Water-Vapor Relative Saturation'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',VAR(NTM,3)
            IF( VAR(NTM,3).GT.1.D+0 .OR. VAR(NTM,3).LT.0.D+0 ) THEN
              INDX = 9
              CHMSG = 'Out of Range Value: ' // 
     &          'Source: Water-Vapor Relative Saturation'
              RLMSG = VAR(NTM,3)
              CALL WRMSGS( INDX )
            ENDIF
!
!---      Initial Condition Fluid Mass Sink  ---
!
          ELSEIF( ISRTX.EQ.5 ) THEN
            VARB = 'Source: Fluid Mass Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            IF( VAR(NTM,4).GT.0.D+0 ) THEN
              INDX = 4
              CHMSG = 'Fluid Mass Rate > 0.0'
              CALL WRMSGS( INDX )
            ENDIF
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
!
!---  Check for nonascending source times  ---
!
          IF( VAR(NTM,1).LT.SRTMO ) THEN
            INDX = 4
            CHMSG = 'Source Time Sequencing'
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
          DO 120 M = 1,8+NSOLU
            SRC(M,NTM,NSR) = VAR(NTM,M)
  120     CONTINUE
  130   CONTINUE
  140 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSR_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDST_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Reads the salt transport card dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
      SUB_LOG(ISUB_LOG) = '/RDST_HYD_O2_CH4_N2'
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
      N = N + 1
  220 CONTINUE
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
     &    DPLGS(IROCK)
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
     &    DPTRS(IROCK)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DPTRS(IROCK),INDX)
        IF( DPTRS(IROCK).GE.SMALL ) IDSPS = 1
      ENDIF
!
!---  Read next rock/soil type or scaling group  ---
!
      IF( N.LT.NROCK ) WRITE(IWR,'(/)')
      GOTO 10
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDST_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTF_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Reads solute/fluid interaction card for diffusion and partition
!     coefficients, and internodal diffusion term averaging scheme for
!     single phase (aqueous) solute transport equation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICDSX
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDTF_HYD_O2_CH4_N2'
!
!---  Write card information to output file  ---
!
      CARD = 'Solute/Fluid Interaction Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of different solutes  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Solutes'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      NSOLU = 0
      DO 200 NL = 1, NLIN
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
!---  Aqueous-phase molecular diffusion coefficient  ---
!
        VARB = 'Aqueous-Phase Molecular Diffusion Coefficient'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDL(NSL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',SMDL(NSL)
        INDX = 0
        IUNM = 2
        IUNS = -1
        CALL RDUNIT(UNTS,SMDL(NSL),INDX)
!
!---  Gas-phase molecular diffusion coefficient  ---
!
        VARB = 'Gas-Phase Molecular Diffusion Coefficient'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDG(NSL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',
     &    UNTS(1:NCH),': ',SMDG(NSL)
        INDX = 0
        IUNM = 2
        IUNS = -1
        CALL RDUNIT(UNTS,SMDG(NSL),INDX)
!
!---  Gas-aqueous partition coefficient option  ---
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
          WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),': '
     &      ,PCGL(1,NSL)
          INDX = 0
          CALL RDUNIT(UNTS,PCGL(1,NSL),INDX)
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 )  THEN
          IPCGL(NSL) = 1
          WRITE( IWR,'(A)' ) ': Temperature Dependent'
          WRITE( IWR,'(A)' ) 'ln( Kgl ) = a + b/T + c ln(T) + dT + eT^2'
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
     &    INDEX(ADUM(1:),'vapor').NE.0 )  THEN
          IPCGL(NSL) = 2
          WRITE( IWR,'(A)' ) ': Water Vapor Equilibrium'
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Gas-Aqueous Partition Option: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Solid-Aqueous Partition option  ---
!
        VARB = 'Solid-Aqueous Partition Option'
        WRITE( IWR,'(/,2A,$)' ) VARB(1:IVR),': '
        WRITE( IWR,'(A)' ) 'Continuous Solid Wetting'
#ifdef radioactive
!
!---  Half-life  ---
!
          IDFLT = 1
          VARB = 'Radioactive Half-Life'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,HLF(NSL))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),': '
     &,HLF(NSL)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,HLF(NSL),INDX)
          HLF(NSL) = MAX( HLF(NSL),SMALL )
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
          CALL WRMSGP( INDX )
        ENDIF
        ALLOCATE( ICLX(1:NSOLU),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: ICLX'
          CALL WRMSGP( INDX )
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
#endif
#ifdef reactive
!
!---  First-order reaction rates  ---
!
        VARB = 'Number of reactions'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NCHEM(NSL))
        IF( NCHEM(NSL).GT.LCHEM ) THEN
          INDX = 5
          CHMSG = 'Number of Reactions > Parameter LCHEM'
          CALL WRMSGS( INDX )
        ENDIF
        DO 120 NCHM = 1,NCHEM(NSL)
          IDFLT = 1
          VARB = 'First-Order Reaction Rate Constant'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,RHLF(NSL,NCHM))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',RHLF(NSL,NCHM)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,RHLF(NSL,NCHM),INDX)
          RHLF(NSL,NCHM) = MAX( RHLF(NSL,NCHM),SMALL )
          RCHDF(NSL,NSL,NCHM) = -1.D+0
  120   CONTINUE
  200 CONTINUE
!
!---  Read number of lines of reaction information  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Reaction Lines'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      IF( NLIN.GT.0 ) THEN
        WRITE(IWR,'(/,A)') 'Reaction Data:'
      ENDIF
      DO 400 NL = 1, NLIN
!
!---  Read new line  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        ADUM(1:) = ' '
        VARB = 'Parent Solute Name'
        NPSL = 0
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        DO 300 NSL = 1,NSOLU
          IF( SOLUT(NSL).EQ.ADUM ) NPSL = NSL
  300   CONTINUE
        IF( NCHEM(NPSL).LE.0 ) THEN
          INDX = 4
          CHMSG = 'Unrecognized Solute Name: '//SOLUT(NSL)
          CALL WRMSGS( INDX )
        ENDIF
        BDUM(1:) = ' '
        VARB = 'Daughter Solute Name'
        NDSL = 0
        CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
        DO 310 NSL = 1,NSOLU
          IF( SOLUT(NSL).EQ.BDUM ) NDSL = NSL
  310   CONTINUE
        IF( NPSL.EQ.0 .OR. NDSL.EQ.0 ) THEN
          INDX = 4
          CHMSG = 'Reaction Solute: '//
     &      ADUM(1:NCHA)//': '//BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ELSEIF( NPSL.EQ.NDSL ) THEN
          INDX = 4
          CHMSG = 'Parent Solute = Progeny Solute: '//
     &      ADUM(1:NCHA)//': '//BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ENDIF
        VARB = 'Reaction Number'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NCHM)
        VARB = 'Reaction Stoichiometry (Fraction)'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,RCHDF(NPSL,NDSL,NCHM))
        WRITE(IWR,'(2X,4A,I2,A,1PE11.4)') ADUM(1:NCHA),' to ',
        WRITE(IWR,'(2X,4A,I2,A,1PE11.4)') ADUM(1:NCHA),' to ',
     &    BDUM(1:NCHB),'; Reaction = #',NCHM,
     &    '; Stoichiometry = ',RCHDF(NPSL,NDSL,NCHM)
  400 CONTINUE
#endif
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTF_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTP_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Reads the solute/porous media interaction card for the
!     dispersivities, half-lives, and partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
      SUB_LOG(ISUB_LOG) = '/RDTP_HYD_O2_CH4_N2'
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
        VARB = 'Longitudinal Dispersivity'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DISPL(IROCK))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),': '
     &,DISPL(IROCK)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DISPL(IROCK),INDX)
        IF( DISPL(IROCK).GE.SMALL ) IDISP = 1
!
!---  Transverse dispersivity  ---
!
        VARB = 'Transverse Dispersivity'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DISPT(IROCK))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),': '
     &,DISPT(IROCK)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DISPT(IROCK),INDX)
        IF( DISPT(IROCK).GE.SMALL ) IDISP = 1
!
!---  Loop over number of solutes or radionuclides  ---
!
        DO 500 NS = 1,NSOLU
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
          VARB = 'Solute Name'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---  Search known solutes for matching name  ---
!
          DO 300 NSL = 1,NSOLU
            IF( ADUM.EQ.SOLUT(NSL)) GOTO 400
  300     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute Name: '//ADUM
          CALL WRMSGS( INDX )
  400     CONTINUE
          WRITE(IWR,'(/,2A)') 'Solute Name:',SOLUT(NSL)
!
!---  Solid-aqueous partition coefficient  ---
!
          IDFLT = 1
          VARB = 'Solid-Aqueous Partition Coefficient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',PCSL(1,IROCK,NSL)
          INDX = 0
          IUNM = 3
          IUNKG = -1
          CALL RDUNIT(UNTS,PCSL(1,IROCK,NSL),INDX)
          PCSL(1,IROCK,NSL) = MAX( PCSL(1,IROCK,NSL),1.D-12 )
  500   CONTINUE
!
!---  Read next rock/soil type or scaling group  ---
!
      IF( N.LT.NROCK ) WRITE(IWR,'(/)')
      GOTO 10
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTP_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKG_HYD_O2_CH4_N2( RKGX,SGX,SHX,SLX,SNX,IZN )
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Gas relative permeability
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RKG_HYD_O2_CH4_N2'
      ISUB_LOG = ISUB_LOG-1
!
!---  Effective saturations for relative permeability
!     calculations  ---
!
      IF( SHX.LT.1.D+0 ) THEN
        SLZ = SLX/(1.D+0-SHX)
        SNZ = SNX/(1.D+0-SHX)
        SGZ = SGX/(1.D+0-SHX)
      ELSE
        SLZ = 0.D+0
        SNZ = 0.D+0
        SGZ = 0.D+0
      ENDIF
      SLZ = SLX
      SNZ = SNX
      SGZ = SGX
      SLRX = SCHR(4,IZN)
      ASLX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
      ASNX = MIN( MAX( SNZ/(1.D+0-SLRX),0.D+0 ),1.D+0 )
      ASGX = MIN( MAX( SGZ/(1.D+0-SLRX),0.D+0 ),1.D+0 )
      ASTX = MIN( MAX( ASLX+ASNX,0.D+0 ),1.D+0 )
!
!---  Constant aqueous relative permeability  ---
!
      IF( IRPG(IZN).EQ.0 ) THEN
        RKGX = RPGC(3,IZN)
!
!---  Mualem aqueous relative permeability  ---
!
      ELSEIF( IRPG(IZN).EQ.1 ) THEN
!
!---    van Genuchten saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 ) THEN
          CMX = RPGC(3,IZN)
          RKGX = SQRT(ASGX)*((1.D+0-ASTX**(1.D+0/CMX))**CMX)**2
!
!---    Brooks and Corey saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 ) THEN
          CLX = RPGC(3,IZN)
          RKGX = SQRT(ASGX)*(1.D+0-ASTX**(1.D+0+1.D+0/CLX))**2
!
!---    van Genuchten saturation function w/
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.13 ) THEN
          CMX = RPGC(3,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SGPX = MIN( MAX( SGZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SNPX = MIN( MAX( SNZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKGX = SQRT(SGPX)*((1.D+0-(SLPX+SNPX)**(1.D+0/CMX))**CMX)**2
!
!---    Brooks and Corey saturation function w/
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.14 ) THEN
          CLX = RPGC(3,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SGPX = MIN( MAX( SGZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SNPX = MIN( MAX( SNZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKGX = SQRT(SGPX)*(1.D+0-(SLPX+SNPX)**(1.D+0+1.D+0/CLX))**2
        ENDIF
!
!---  Burdine aqueous relative permeability  ---
!
      ELSEIF( IRPG(IZN).EQ.2 ) THEN
!
!---    van Genuchten saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 ) THEN
          CMX = RPGC(3,IZN)
          RKGX = (ASGX**2)*((1.D+0-ASTX**(1.D+0/CMX))**CMX)
!
!---    Brooks and Corey saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 ) THEN
          CLX = RPGC(3,IZN)
          RKGX = (ASGX**2)*(1.D+0-ASTX**(1.D+0 + 2.D+0/CLX))
!
!---    van Genuchten saturation function w/
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.13 ) THEN
          CMX = RPGC(3,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SGPX = MIN( MAX( SGZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SNPX = MIN( MAX( SNZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKGX = (SGPX**2)*((1.D+0-(SLPX+SNPX)**(1.D+0/CMX))**CMX)
!
!---    Brooks and Corey saturation function w/
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.14 ) THEN
          CLX = RPGC(3,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SGPX = MIN( MAX( SGZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SNPX = MIN( MAX( SNZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKGX = (SGPX**2)*(1.D+0-(SLPX+SNPX)**(1.D+0 + 2.D+0/CLX))
        ENDIF
!
!---  Free Corey relative permeability function,
!     modified from STOMP-O2 for the Moridis
!     code comparison problem  ---
!
      ELSEIF( IRPG(IZN).EQ.7 ) THEN
        SLRX = RPGC(3,IZN)
        SGRX = RPGC(4,IZN)
        SGPX = MIN( MAX( (SGZ-SGRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        RKGX = RPGC(1,IZN)*(SGPX**(RPGC(2,IZN)))
!
!---  Unrecognized aqueous relative permeability  ---
!
      ELSE
        CHMSG = 'Unrecognized Gas Relative Permeability' //
     &      ' Function: IRPG(IZN) = '
        INDX = 12
        IMSG = IRPG(IZN)
        CALL WRMSGS( INDX )
      ENDIF
      RKGX = MAX( RKGX,0.D+0 )
!
!---  End of RKG_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKL_HYD_O2_CH4_N2( RKLX,SGX,SHX,SLX,SNX,IZN )
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Aqueous relative permeability
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/RKL_HYD_O2_CH4_N2'
      ISUB_LOG = ISUB_LOG-1
!
!---  Effective saturations for relative permeability
!     calculations  ---
!
      IF( SHX.LT.1.D+0 ) THEN
        SLZ = SLX/(1.D+0-SHX)
        SNZ = SNX/(1.D+0-SHX)
        SGZ = SGX/(1.D+0-SHX)
      ELSE
        SLZ = 0.D+0
        SNZ = 0.D+0
        SGZ = 0.D+0
      ENDIF
      SLZ = SLX
      SNZ = SNX
      SGZ = SGX
      SLRX = SCHR(4,IZN)
      ASLX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
      ASNX = MIN( MAX( SNZ/(1.D+0-SLRX),0.D+0 ),1.D+0 )
      ASGX = MIN( MAX( SGZ/(1.D+0-SLRX),0.D+0 ),1.D+0 )
      ASTX = MIN( MAX( ASLX+ASNX,0.D+0 ),1.D+0 )
!
!---  Constant aqueous relative permeability  ---
!
      IF( IRPL(IZN).EQ.0 ) THEN
        RKLX(1) = RPLC(2,IZN)
!
!---  Mualem aqueous relative permeability  ---
!
      ELSEIF( IRPL(IZN).EQ.1 ) THEN
!
!---    van Genuchten saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 ) THEN
          CMX = RPLC(2,IZN)
          RKLX(1) = SQRT(ASLX)*(1.D+0 -
     &    (1.D+0-ASLX**(1.D+0/CMX))**CMX)**2
!
!---    Brooks and Corey saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 ) THEN
          CLX = RPLC(2,IZN)
          RKLX(1) = ASLX**(2.5D+0 + 2.D+0/CLX)
!
!---    van Genuchten saturation function w/ the
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.13 ) THEN
          CMX = RPLC(2,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKLX(1) = SQRT(SLPX)*(1.D+0 -
     &    (1.D+0-SLPX**(1.D+0/CMX))**CMX)**2
!
!---    Brooks and Corey saturation function w/ the
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.14 ) THEN
          CLX = RPLC(2,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKLX(1) = SLPX**(2.5D+0 + 2.D+0/CLX)
        ENDIF
!
!---  Burdine aqueous relative permeability  ---
!
      ELSEIF( IRPL(IZN).EQ.2 ) THEN
!
!---    van Genuchten saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 ) THEN
          CMX = RPLC(2,IZN)
          RKLX(1) = (ASLX**2)*(1.D+0 -
     &      (1.D+0-ASLX**(1.D+0/CMX))**CMX)
!
!---    Brooks and Corey saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 ) THEN
          CLX = RPLC(2,IZN)
          RKLX(1) = ASLX**(3.0D+0 + 2.D+0/CLX)
!
!---    van Genuchten saturation function w/ the
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.13 ) THEN
          CMX = RPLC(2,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKLX(1) = (SLPX**2)*(1.D+0 -
     &      (1.D+0-SLPX**(1.D+0/CMX))**CMX)
!
!---    Brooks and Corey saturation function w/ the
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.14 ) THEN
          CLX = RPLC(2,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKLX(1) = SLPX**(3.0D+0 + 2.D+0/CLX)
        ENDIF
!
!---  Free Corey relative permeability function,
!     modified from STOMP-O2 for the Moridis
!     code comparison problem  ---
!
      ELSEIF( IRPL(IZN).EQ.7 ) THEN
        SLRX = RPLC(3,IZN)
        SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        RKLX(1) = RPLC(1,IZN)*(SLPX**(RPLC(2,IZN)))
!
!---  Modified Mualem aqueous relative permeability  ---
!
      ELSEIF( IRPL(IZN).EQ.22 ) THEN
!
!---    van Genuchten saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 ) THEN
          RKLX(1) = (ASLX**RPLC(1,IZN))*((1.D+0-(1.D+0-ASLX**
     &      (1.D+0/RPLC(2,IZN)))**RPLC(2,IZN))**2)
!
!---    Brooks and Corey saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 ) THEN
          RKLX(1) = (ASLX**RPLC(1,IZN))*
     &      ASLX**(2.D+0 + 2.0D+0/RPLC(2,IZN))
!
!---    van Genuchten saturation function w/ the
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.13 ) THEN
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKLX(1) = (SLPX**RPLC(1,IZN))*((1.D+0-(1.D+0-SLPX**
     &      (1.D+0/RPLC(2,IZN)))**RPLC(2,IZN))**2)
!
!---    Brooks and Corey saturation function w/ the
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.14 ) THEN
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKLX(1) = (SLPX**RPLC(1,IZN))*
     &      SLPX**(2.D+0 + 2.0D+0/RPLC(2,IZN))
        ENDIF
!
!---  Unrecognized aqueous relative permeability  ---
!
      ELSE
        CHMSG = 'Unrecognized Aqueous Relative Permeability' //
     &    ' Function: IRPL(IZN) = '
        INDX = 12
        IMSG = IRPL(IZN)
        CALL WRMSGS( INDX )
      ENDIF
      RKLX(1) = MAX( RKLX(1),0.D+0 )
      RKLX(2) = MAX( RKLX(1),0.D+0 )
      RKLX(3) = MAX( RKLX(1),0.D+0 )
!
!---  End of RKL_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RKN_HYD_O2_CH4_N2( RKNX,SGX,SHX,SLX,SNX,IZN )
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Nonaqueous-liquid relative permeability
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RKN_HYD_O2_CH4_N2'
      ISUB_LOG = ISUB_LOG-1
!
!---  Effective saturations for relative permeability
!     calculations  ---
!
      IF( SHX.LT.1.D+0 ) THEN
        SLZ = SLX/(1.D+0-SHX)
        SNZ = SNX/(1.D+0-SHX)
        SGZ = SGX/(1.D+0-SHX)
      ELSE
        SLZ = 0.D+0
        SNZ = 0.D+0
        SGZ = 0.D+0
      ENDIF
      SLZ = SLX
      SNZ = SNX
      SGZ = SGX
      SLRX = SCHR(4,IZN)
      ASLX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
      ASNX = MIN( MAX( SNZ/(1.D+0-SLRX),0.D+0 ),1.D+0 )
      ASGX = MIN( MAX( SGZ/(1.D+0-SLRX),0.D+0 ),1.D+0 )
      ASTX = MIN( MAX( ASLX+ASNX,0.D+0 ),1.D+0 )
!
!---  Constant aqueous relative permeability  ---
!
      IF( IRPN(IZN).EQ.0 ) THEN
        RKNX = RPNC(1,IZN)
!
!---  Mualem aqueous relative permeability  ---
!
      ELSEIF( IRPN(IZN).EQ.1 ) THEN
!
!---    van Genuchten saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 ) THEN
          CMX = RPNC(1,IZN)
          RKNX = SQRT(ASNX)*((1.D+0 - (ASTX-ASNX)**(1.D+0/CMX))**CMX -
     &      (1.D+0 - (ASTX)**(1.D+0/CMX))**CMX)**2
!
!---    Brooks and Corey saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 ) THEN
          CLX = RPNC(1,IZN)
          RKNX = SQRT(ASNX)*( ASTX**(2.0D+0 + 2.D+0/CLX) -
     &        (ASTX-ASNX)**(2.0D+0 + 2.D+0/CLX) )
!
!---    van Genuchten saturation function w/
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.13 ) THEN
          CMX = RPNC(1,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SNPX = MIN( MAX( SNZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKNX = SQRT(SNPX)*((1.D+0 - (SLPX)**(1.D+0/CMX))**CMX -
     &      (1.D+0 - (SLPX+SNPX)**(1.D+0/CMX))**CMX)**2
!
!---    Brooks and Corey saturation function w/
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.14 ) THEN
          CLX = RPNC(1,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SNPX = MIN( MAX( SNZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKNX = SQRT(SNPX)*( (SLPX+SNPX)**(2.0D+0 + 2.D+0/CLX) -
     &      SLPX**(2.0D+0 + 2.D+0/CLX) )
        ENDIF
!
!---  Burdine aqueous relative permeability  ---
!
      ELSEIF( IRPN(IZN).EQ.2 ) THEN
!
!---    van Genuchten saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 ) THEN
          CMX = RPNC(1,IZN)
          RKNX = (ASNX**2)*((1.D+0 - (ASTX-ASNX)**(1.D+0/CMX))**CMX -
     &      (1.D+0 - (ASTX)**(1.D+0/CMX))**CMX)
!
!---    Brooks and Corey saturation function w/o or w/ the
!       hydrate capillary pressure barrier  ---
!
        ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 ) THEN
          CLX = RPNC(1,IZN)
          RKNX = (ASNX**2)*( ASTX**(1.0D+0 + 2.D+0/CLX) -
     &      (ASTX-ASNX)**(1.0D+0 + 2.D+0/CLX) )
!
!---    van Genuchten saturation function w/
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.13 ) THEN
          CMX = RPNC(1,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SNPX = MIN( MAX( SNZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKNX = (SNPX**2)*((1.D+0 - (SLPX)**(1.D+0/CMX))**CMX -
     &      (1.D+0 - (SLPX+SNPX)**(1.D+0/CMX))**CMX)
!
!---    Brooks and Corey saturation function w/
!       hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.14 ) THEN
          CLX = RPNC(1,IZN)
          SLRX = SCHR(4,IZN)
          SLPX = MIN( MAX( (SLZ-SLRX)/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          SNPX = MIN( MAX( SNZ/(1.D+0-SHX-SLRX),0.D+0 ),1.D+0 )
          RKNX = (SNPX**2)*( (SLPX+SNPX)**(1.0D+0 + 2.D+0/CLX) -
     &      SLPX**(1.0D+0 + 2.D+0/CLX) )
        ENDIF
!
!---  Free Corey relative permeability function,
!     modified from STOMP-O2 for the Moridis
!     code comparison problem  ---
!
      ELSEIF( IRPN(IZN).EQ.7 ) THEN
        SLRX = RPNC(3,IZN)
        SNRX = RPNC(4,IZN)
        SNPX = MIN( MAX( (SNZ-SNRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
        RKNX = RPNC(1,IZN)*(SNPX**(RPNC(2,IZN)))
!
!---  Unrecognized aqueous relative permeability  ---
!
      ELSE
        CHMSG = 'Unrecognized NAPL Relative Permeability' //
     &      ' Function: IRPN(IZN) = '
        INDX = 12
        IMSG = IRPN(IZN)
        CALL WRMSGS( INDX )
      ENDIF
      RKNX = MAX( RKNX,0.D+0 )
!
!---  End of RKN_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Compute the maximum relative residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE OUTPU
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVI
      USE FDVH
      USE FDVA
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*128 PH_CND(43)
      INTEGER IRSDX(LFD)
!
!----------------------Data Statements---------------------------------!
!
      DATA PH_CND /9*'','SH = 0, SG = 0, SN = 0',10*'',
     &  'SH = 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN )',
     &  'SH = 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN )',
     &  'SH = 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO )',7*'',
     &  'SH > 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN ), YMHGA >= MAX(YMH
     &GO,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN ), YMHGA >= MAX(YMH
     &GO,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO ), YMHGA >= MAX(YMH
     &GO,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN ), YMHGO >= MAX(YMH
     &GA,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN ), YMHGO >= MAX(YMH
     &GA,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO ), YMHGO >= MAX(YMH
     &GA,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN ), YMHGN >= MAX(YMH
     &GA,YMHGO )',
     &  'SH > 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN ), YMHGN >= MAX(YMH
     &GA,YMHGO )',
     &  'SH > 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO ), YMHGN >= MAX(YMH
     &GA,YMHGO )','',
     &  'SH > 0, SG + SN = 0, YMHGA >= MAX(YMHGO,YMHGN )',
     &  'SH > 0, SG + SN = 0, YMHGO >= MAX(YMHGA,YMHGN )',
     &  'SH > 0, SG + SN = 0, YMHGN >= MAX(YMHGA,YMHGO )'/
!
!----------------------Executable Lines--------------------------------!
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL_HYD_O2_CH4_N2'
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
      DO 200 N = 1,NFLD
        N_DB = N
        IF( IXP(N).EQ.0 ) GOTO 200
        NMD = IXP(N)
        IZN = IZ(N)
        MPT = IM(IEQT,NMD)
        DPT = BLU(MPT)
        MPW = IM(IEQW,NMD)
        DPW = BLU(MPW)
        IF( ISLC(45).EQ.0 ) THEN
          MPA = IM(IEQA,NMD)
          DPA = BLU(MPA)
        ELSE
          DPA = 0.D+0
        ENDIF
        IF( ISLC(64).EQ.0 ) THEN
          MPO = IM(IEQO,NMD)
          DPO = BLU(MPO)
        ELSE
          DPO = 0.D+0
        ENDIF
        IF( ISLC(65).EQ.0 ) THEN
          MPN = IM(IEQN,NMD)
          DPN = BLU(MPN)
        ELSE
          DPN = 0.D+0
        ENDIF
        IF( ISLC(45).EQ.0 ) THEN
          MPHA = IM(IEQHA,NMD)
          DPHA = BLU(MPHA)
        ELSE
          DPHA = 0.D+0
        ENDIF
        IF( ISLC(64).EQ.0 ) THEN
          MPHO = IM(IEQHO,NMD)
          DPHO = BLU(MPHO)
        ELSE
          DPHO = 0.D+0
        ENDIF
        IF( ISLC(65).EQ.0 ) THEN
          MPHN = IM(IEQHN,NMD)
          DPHN = BLU(MPHN)
        ELSE
          DPHN = 0.D+0
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
!---    Skip selected nodes in the residual calculation  ---
!
        IF( ISKP(IZN).EQ.1 ) GOTO 200
!
!---    Salt mass equation, isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          ACP = TMS(2,N)*DTI*VOL(N)
          CALL SOL_LS( T(2,N),XLSMX )
          RSDX = MIN( (ABS(DPS)/XLSMX),
     &      ABS(RSDL(IEQS,N)/(ACP+SMALL)) )
          RSDX = RSDX*1.D-1
          IF( RSDX.GT.RSD(IEQS) ) THEN
            RSD(IEQS) = RSDX
            NSD(IEQS) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
        ENDIF
!
!---    Core phase condition  ---
!
        NPHZX = MOD(NPHAZ(2,N),100)
!
!---    Phase condition #10 (10 series)  ---
!
!       SH = 0, SG + SN = 0, SL + SI = 1
!
!       energy - temperature, C - T
!       water mass - pressure, Pa (gauge) - PSO
!       mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!       mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!       mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!       hydrate O2 mass - hydrate O2 mass, kg - TMHA
!       hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!       hydrate N2 mass - hydrate N2 mass, kg - TMHN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        IF( NPHZX.EQ.10 ) THEN
          PX = PSO(2,N) + PATM
          PX = MAX( PX,PG(2,N)-PL(2,N) )
!
!---      Energy equation  ---
!
          RSDX = ABS(DPT)/TABS
          IF( RSDX.GT.RSD(IEQT) ) THEN
            RSD(IEQT) = RSDX
            NSD(IEQT) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      Water residual  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) + 
     &      RHOI(2,N)*SI(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N) + 
     &      RHON(2,N)*SN(2,N)*XNW(2,N) + 
     &      RHOH(2,N)*SH(2,N)*XHW(2,N) )*DTI*VOL(N)
          RSDX = MIN( ABS(DPW)/PX,
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( POSM(2,N).LT.1.D-6 ) RSDX = 0.D+0
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      Mobile O2 residual  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( PVA(2,N).GT.(1.D-6*PX) ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLA(2,N) + 
     &          RHON(2,N)*SN(2,N)*XNA(2,N))*DTI*VOL(N)
              RSDX = MIN( ABS(DPA)/PX,
     &          ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD(IEQA) ) THEN
                RSD(IEQA) = RSDX
                NSD(IEQA) = N
              ENDIF
              IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
            ENDIF
          ENDIF
!
!---      Mobile CH4 residual  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            IF( PVO(2,N).GT.(1.D-6*PX) ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGO(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLO(2,N) + 
     &          RHON(2,N)*SN(2,N)*XNO(2,N))*DTI*VOL(N)
              RSDX = MIN( ABS(DPO)/PX,
     &          ABS(RSDL(IEQO,N)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD(IEQO) ) THEN
                RSD(IEQO) = RSDX
                NSD(IEQO) = N
              ENDIF
              IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
            ENDIF
          ENDIF
!
!---      Mobile N2 residual  ---
!
          IF( ISLC(65).EQ.0 ) THEN
            IF( PVN(2,N).GT.(1.D-6*PX) ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGN(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLN(2,N) + 
     &          RHON(2,N)*SN(2,N)*XNN(2,N))*DTI*VOL(N)
              RSDX = MIN( ABS(DPN)/PX,
     &          ABS(RSDL(IEQN,N)/(ACP+SMALL)) )
              IF( RSDX.GT.RSD(IEQN) ) THEN
                RSD(IEQN) = RSDX
                NSD(IEQN) = N
              ENDIF
              IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
            ENDIF
          ENDIF
!
!---    Phase conditions 21 through 23 (20 series)  ---
!
!       SH = 0, SG + SN > 0, SL + SI < 1
!       ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #21
!       ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #22
!       ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #23
!
!       energy - temperature, C - T
!       water mass - aqueous + ice saturation - POSM
!       maximum mobile mass component - pressure, Pa (gauge) - PSO
!       mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!       mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!       mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!       hydrate O2 mass - hydrate O2 mass, kg - TMHA
!       hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!       hydrate N2 mass - hydrate N2 mass, kg - TMHN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.21 .AND. NPHZX.LE.23  ) THEN
          PX = PSO(2,N) + PATM
          PX = MAX( PX,PG(2,N)-PL(2,N) )
!
!---      Energy equation  ---
!
          RSDX = ABS(DPT)/TABS
          IF( RSDX.GT.RSD(IEQT) ) THEN
            RSD(IEQT) = RSDX
            NSD(IEQT) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      Water residual  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) + 
     &      RHOI(2,N)*SI(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N) + 
     &      RHON(2,N)*SN(2,N)*XNW(2,N) + 
     &      RHOH(2,N)*SH(2,N)*XHW(2,N) )*DTI*VOL(N)
          RSDX = MIN( ABS(DPW)/PX,
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( POSM(2,N).LT.1.D-6 ) RSDX = 0.D+0
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      Mobile O2 residual  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N) + 
     &        RHON(2,N)*SN(2,N)*XNA(2,N))*DTI*VOL(N)
            IF( NPHZX.EQ.21 ) THEN
              IF( PSO(2,N).GT.(1.D-6*PX) ) THEN
                RSDX = MIN( ABS(DPA)/PX,
     &            ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              ENDIF
            ELSE
              RSDX = MIN( ABS(DPA),
     &          ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              IF( ZMCA(2,N).LT.EPSL .AND. DPA.LT.0.D+0 ) RSDX = 0.D+0
            ENDIF
            IF( RSDX.GT.RSD(IEQA) ) THEN
              RSD(IEQA) = RSDX
              NSD(IEQA) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Mobile CH4 residual  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGO(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLO(2,N) + 
     &        RHON(2,N)*SN(2,N)*XNO(2,N))*DTI*VOL(N)
            IF( NPHZX.EQ.22 ) THEN
              IF( PSO(2,N).GT.(1.D-6*PX) ) THEN
                RSDX = MIN( ABS(DPO)/PX,
     &            ABS(RSDL(IEQO,N)/(ACP+SMALL)) )
              ENDIF
            ELSE
              RSDX = MIN( ABS(DPO),
     &          ABS(RSDL(IEQO,N)/(ACP+SMALL)) )
              IF( ZMCO(2,N).LT.EPSL .AND. DPO.LT.0.D+0 ) RSDX = 0.D+0
            ENDIF
            IF( RSDX.GT.RSD(IEQO) ) THEN
              RSD(IEQO) = RSDX
              NSD(IEQO) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Mobile N2 residual  ---
!
          IF( ISLC(65).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGN(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLN(2,N) + 
     &        RHON(2,N)*SN(2,N)*XNN(2,N))*DTI*VOL(N)
            IF( NPHZX.EQ.22 ) THEN
              IF( PSO(2,N).GT.(1.D-6*PX) ) THEN
                RSDX = MIN( ABS(DPN)/PX,
     &            ABS(RSDL(IEQN,N)/(ACP+SMALL)) )
              ENDIF
            ELSE
              RSDX = MIN( ABS(DPN),
     &          ABS(RSDL(IEQN,N)/(ACP+SMALL)) )
              IF( ZMCN(2,N).LT.EPSL .AND. DPN.LT.0.D+0 ) RSDX = 0.D+0
            ENDIF
            IF( RSDX.GT.RSD(IEQN) ) THEN
              RSD(IEQN) = RSDX
              NSD(IEQN) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---    Phase conditions #31 through #39 (30 series)  ---
!
!       SH > 0, SG + SN > 0, SL + SI < 1
!       ZMCA = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #31
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #32
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #33
!       ZMCO = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #34
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #35
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #36
!       ZMCN = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #37
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #38
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #39
!
!       energy - temperature, C - T
!       water mass - aqueous + ice saturation, POSM
!       maximum mobile mass component - pressure, Pa (gauge) - PSO
!       mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!       mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!       mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!       maximum hydrate mass component - hydrate saturation - SH
!       hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!       hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!       hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.31 .AND. NPHZX.LE.39 ) THEN
          PX = PSO(2,N) + PATM
          PX = MAX( PX,PG(2,N)-PL(2,N) )
!
!---      Energy equation  ---
!
          RSDX = ABS(DPT)/TABS
          IF( RSDX.GT.RSD(IEQT) ) THEN
            RSD(IEQT) = RSDX
            NSD(IEQT) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      Water residual  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) + 
     &      RHOI(2,N)*SI(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N) + 
     &      RHON(2,N)*SN(2,N)*XNW(2,N) + 
     &      RHOH(2,N)*SH(2,N)*XHW(2,N) )*DTI*VOL(N)
          RSDX = MIN( ABS(DPW)/PX,
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( POSM(2,N).LT.1.D-6 ) RSDX = 0.D+0
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      Mobile O2 residual  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLA(2,N) + 
     &        RHON(2,N)*SN(2,N)*XNA(2,N))*DTI*VOL(N)
            IF( NPHZX.GE.31 .AND. NPHZX.LE.33 ) THEN
              IF( PSO(2,N).GT.(1.D-6*PX) ) THEN
                RSDX = MIN( ABS(DPA)/PX,
     &            ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              ENDIF
            ELSE
              RSDX = MIN( ABS(DPA),
     &          ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              IF( ZMCA(2,N).LT.EPSL .AND. DPA.LT.0.D+0 ) RSDX = 0.D+0
            ENDIF
            RSDX = 1.D-1*RSDX
            IF( RSDX.GT.RSD(IEQA) ) THEN
              RSD(IEQA) = RSDX
              NSD(IEQA) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Mobile CH4 residual  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGO(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLO(2,N) + 
     &        RHON(2,N)*SN(2,N)*XNO(2,N))*DTI*VOL(N)
            IF( NPHZX.GE.34 .AND. NPHZX.LE.36 ) THEN
              IF( PSO(2,N).GT.(1.D-6*PX) ) THEN
                RSDX = MIN( ABS(DPO)/PX,
     &            ABS(RSDL(IEQO,N)/(ACP+SMALL)) )
              ENDIF
            ELSE
              RSDX = MIN( ABS(DPO),
     &          ABS(RSDL(IEQO,N)/(ACP+SMALL)) )
              IF( ZMCO(2,N).LT.EPSL .AND. DPO.LT.0.D+0 ) RSDX = 0.D+0
            ENDIF
            RSDX = 1.D-1*RSDX
            IF( RSDX.GT.RSD(IEQO) ) THEN
              RSD(IEQO) = RSDX
              NSD(IEQO) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Mobile N2 residual  ---
!
          IF( ISLC(65).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGN(2,N) +
     &        RHOL(2,N)*SL(2,N)*XLN(2,N) + 
     &        RHON(2,N)*SN(2,N)*XNN(2,N))*DTI*VOL(N)
            IF( NPHZX.GE.37 .AND. NPHZX.LE.39 ) THEN
              IF( PSO(2,N).GT.(1.D-6*PX) ) THEN
                RSDX = MIN( ABS(DPN)/PX,
     &            ABS(RSDL(IEQN,N)/(ACP+SMALL)) )
              ENDIF
            ELSE
              RSDX = MIN( ABS(DPN),
     &          ABS(RSDL(IEQN,N)/(ACP+SMALL)) )
              IF( ZMCN(2,N).LT.EPSL .AND. DPN.LT.0.D+0 ) RSDX = 0.D+0
            ENDIF
            RSDX = 1.D-1*RSDX
            IF( RSDX.GT.RSD(IEQN) ) THEN
              RSD(IEQN) = RSDX
              NSD(IEQN) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Hydrate-O2 residual  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOH(2,N)*SH(2,N)*XHA(2,N))*DTI*VOL(N)
            RSDX = 1.D-3*MIN( ABS(DPHA),
     &        ABS(RSDL(IEQHA,N)/(ACP+SMALL)) )
            IF( YMHGA(2,N).LT.1.D-9 ) RSDX = 0.D+0
            IF( RSDX.GT.RSD(IEQHA) ) THEN
              RSD(IEQHA) = RSDX
              NSD(IEQHA) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Hydrate-CH4 residual  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOH(2,N)*SH(2,N)*XHO(2,N))*DTI*VOL(N)
            RSDX = 1.D-3*MIN( ABS(DPHO),
     &        ABS(RSDL(IEQHO,N)/(ACP+SMALL)) )
            IF( YMHGO(2,N).LT.1.D-9 ) RSDX = 0.D+0
            IF( RSDX.GT.RSD(IEQHO) ) THEN
              RSD(IEQHO) = RSDX
              NSD(IEQHO) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Hydrate-N2 residual  ---
!
          IF( ISLC(65).EQ.0 ) THEN 
            ACP = PORD(2,N)*(RHOH(2,N)*SH(2,N)*XHN(2,N))*DTI*VOL(N)
            RSDX = 1.D-3*MIN( ABS(DPHN),
     &        ABS(RSDL(IEQHN,N)/(ACP+SMALL)) )
            IF( YMHGN(2,N).LT.1.D-9 ) RSDX = 0.D+0
            IF( RSDX.GT.RSD(IEQHN) ) THEN
              RSD(IEQHN) = RSDX
              NSD(IEQHN) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---    Phase conditions #41 through #43 (40 series)  ---
!
!       SH > 0, SG + SN = 0, SL + SI < 1
!       YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #41
!       YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #42
!       YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #43
!
!       energy - temperature, C - T
!       water mass - pressure, Pa (gauge) - PSO
!       mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!       mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!       mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!       maximum hydrate mass component - hydrate saturation - SH
!       hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!       hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!       hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.41 .AND. NPHZX.LE.43 ) THEN
          PX = PSO(2,N) + PATM
          PX = MAX( PX,PG(2,N)-PL(2,N) )
!
!---      Energy equation  ---
!
          RSDX = ABS(DPT)/TABS
          IF( RSDX.GT.RSD(IEQT) ) THEN
            RSD(IEQT) = RSDX
            NSD(IEQT) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      Water residual  ---
!
          ACP = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) + 
     &      RHOI(2,N)*SI(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N) + 
     &      RHON(2,N)*SN(2,N)*XNW(2,N) + 
     &      RHOH(2,N)*SH(2,N)*XHW(2,N) )*DTI*VOL(N)
          RSDX = MIN( ABS(DPW)/PX,
     &      ABS(RSDL(IEQW,N)/(ACP+SMALL)) )
          IF( POSM(2,N).LT.1.D-6 ) RSDX = 0.D+0
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
!
!---      Mobile O2 residual  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( PVA(2,N).GT.(1.D-6*PX) ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGA(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLA(2,N) + 
     &          RHON(2,N)*SN(2,N)*XNA(2,N))*DTI*VOL(N)
              RSDX = MIN( ABS(DPA)/PX,
     &          ABS(RSDL(IEQA,N)/(ACP+SMALL)) )
              RSDX = 1.D-1*RSDX
              IF( RSDX.GT.RSD(IEQA) ) THEN
                RSD(IEQA) = RSDX
                NSD(IEQA) = N
              ENDIF
              IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
            ENDIF
          ENDIF
!
!---      Mobile CH4 residual  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            IF( PVO(2,N).GT.(1.D-6*PX) ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGO(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLO(2,N) + 
     &          RHON(2,N)*SN(2,N)*XNO(2,N))*DTI*VOL(N)
              RSDX = MIN( ABS(DPO)/PX,
     &          ABS(RSDL(IEQO,N)/(ACP+SMALL)) )
              RSDX = 1.D-1*RSDX
              IF( RSDX.GT.RSD(IEQO) ) THEN
                RSD(IEQO) = RSDX
                NSD(IEQO) = N
              ENDIF
              IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
            ENDIF
          ENDIF
!
!---      Mobile N2 residual  ---
!
          IF( ISLC(65).EQ.0 ) THEN
            IF( PVN(2,N).GT.(1.D-6*PX) ) THEN
              ACP = PORD(2,N)*(RHOG(2,N)*SG(2,N)*XGN(2,N) +
     &          RHOL(2,N)*SL(2,N)*XLN(2,N) + 
     &          RHON(2,N)*SN(2,N)*XNN(2,N))*DTI*VOL(N)
              RSDX = MIN( ABS(DPN)/PX,
     &          ABS(RSDL(IEQN,N)/(ACP+SMALL)) )
              RSDX = 1.D-1*RSDX
              IF( RSDX.GT.RSD(IEQN) ) THEN
                RSD(IEQN) = RSDX
                NSD(IEQN) = N
              ENDIF
              IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
            ENDIF
          ENDIF
!
!---      Hydrate-O2 residual  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOH(2,N)*SH(2,N)*XHA(2,N))*DTI*VOL(N)
            RSDX = 1.D-3*MIN( ABS(DPHA),
     &        ABS(RSDL(IEQHA,N)/(ACP+SMALL)) )
            IF( YMHGA(2,N).LT.1.D-9 ) RSDX = 0.D+0
            IF( RSDX.GT.RSD(IEQHA) ) THEN
              RSD(IEQHA) = RSDX
              NSD(IEQHA) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Hydrate-CH4 residual  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            ACP = PORD(2,N)*(RHOH(2,N)*SH(2,N)*XHO(2,N))*DTI*VOL(N)
            RSDX = 1.D-3*MIN( ABS(DPHO),
     &        ABS(RSDL(IEQHO,N)/(ACP+SMALL)) )
            IF( YMHGO(2,N).LT.1.D-9 ) RSDX = 0.D+0
            IF( RSDX.GT.RSD(IEQHO) ) THEN
              RSD(IEQHO) = RSDX
              NSD(IEQHO) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
!
!---      Hydrate-N2 residual  ---
!
          IF( ISLC(65).EQ.0 ) THEN 
            ACP = PORD(2,N)*(RHOH(2,N)*SH(2,N)*XHN(2,N))*DTI*VOL(N)
            RSDX = 1.D-3*MIN( ABS(DPHN),
     &        ABS(RSDL(IEQHN,N)/(ACP+SMALL)) )
            IF( YMHGN(2,N).LT.1.D-9 ) RSDX = 0.D+0
            IF( RSDX.GT.RSD(IEQHN) ) THEN
              RSD(IEQHN) = RSDX
              NSD(IEQHN) = N
            ENDIF
            IF( RSDX.GT.RSDMX ) IRSDX(N) = 1
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  Assign a convergence index  ---
!
      DO 300 M = 1,ISVC
        IF( RSD(M).GT.RSDMX ) ICNV = 2
  300 CONTINUE
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) ICNV = 1
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.1 ) THEN
        WRITE(ISC,'(10X,A)') '---  Convergence Failure  ---'
        WRITE(IWR,'(10X,A)') '---  Convergence Failure  ---'
!
!---    Energy equation  ---
!
        NX = NSD(IEQT)
        IF( NX.GT.0 ) THEN
          NPX = MOD(NPHAZ(2,NX),100)
          NCHX = INDEX( PH_CND(NPX),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &      'Energy Equation Maximum Residual = ',
     &      RSD(IEQT),': Node = ',NX,
     &      ': Temperature = ',T(2,NX),
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &      'Energy Equation Maximum Residual = ',
     &      RSD(IEQT),': Node = ',NX,
     &      ': Temperature = ',T(2,NX),
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
        ENDIF
!
!---    Water mass equation  ---
!
        NX = NSD(IEQW)
        IF( NX.GT.0 ) THEN
          NPX = MOD(NPHAZ(2,NX),100)
          NCHX = INDEX( PH_CND(NPX),'  ') - 1
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &      'Water Mass Equation Maximum Residual = ',
     &      RSD(IEQW),': Node = ',NX,
     &      ': Aqueous Saturation = ',SL(2,NX),
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &      'Water Mass Equation Maximum Residual = ',
     &      RSD(IEQW),': Node = ',NX,
     &      ': Aqueous Saturation = ',SL(2,NX),
     &      ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
        ENDIF
!
!---    Mobile-O2 mass equation  ---
!
        IF( ISLC(45).EQ.0 ) THEN
          NX = NSD(IEQA)
          IF( NX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NX),100)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Mobile-O2 Mass Equation Maximum Residual = ',
     &        RSD(IEQA),': Node = ',NX,
     &        ': Gas Saturation = ',SG(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Mobile-O2 Mass Equation Maximum Residual = ',
     &        RSD(IEQA),': Node = ',NX,
     &        ': Gas Saturation = ',SG(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---    Mobile-CH4 mass equation  ---
!
        IF( ISLC(64).EQ.0 ) THEN
          NX = NSD(IEQO)
          IF( NX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NX),100)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Mobile-CH4 Mass Equation Maximum Residual = ',
     &        RSD(IEQO),': Node = ',NX,
     &        ': Gas Saturation = ',SG(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Mobile-CH4 Mass Equation Maximum Residual = ',
     &        RSD(IEQO),': Node = ',NX,
     &        ': Gas Saturation = ',SG(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---    Mobile-N2 mass equation  ---
!
        IF( ISLC(65).EQ.0 ) THEN
          NX = NSD(IEQN)
          IF( NX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NX),100)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Mobile-N2 Mass Equation Maximum Residual = ',
     &        RSD(IEQN),': Node = ',NX,
     &        ': Gas Saturation = ',SL(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Mobile-N2 Mass Equation Maximum Residual = ',
     &        RSD(IEQN),': Node = ',NX,
     &        ': Gas Saturation = ',SL(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---    Hydrate-O2 mass equation  ---
!
        IF( ISLC(45).EQ.0 ) THEN
          NX = NSD(IEQHA)
          IF( NX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NX),100)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Hydrate-O2 Mass Equation Maximum Residual = ',
     &        RSD(IEQHA),': Node = ',NX,
     &        ': Hydrate Saturation = ',SH(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Hydrate-O2 Mass Equation Maximum Residual = ',
     &        RSD(IEQHA),': Node = ',NX,
     &        ': Hydrate Saturation = ',SH(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---    Hydrate-CH4 mass equation  ---
!
        IF( ISLC(64).EQ.0 ) THEN
          NX = NSD(IEQHO)
          IF( NX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NX),100)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Hydrate-CH4 Mass Equation Maximum Residual = ',
     &        RSD(IEQHO),': Node = ',NX,
     &        ': Hydrate Saturation = ',SH(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Hydrate-CH4 Mass Equation Maximum Residual = ',
     &        RSD(IEQHO),': Node = ',NX,
     &        ': Hydrate Saturation = ',SH(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---    Hydrate-N2 mass equation  ---
!
        IF( ISLC(65).EQ.0 ) THEN 
          NX = NSD(IEQHN)
          IF( NX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NX),100)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        ': Hydrate-N2 Mass Equation Maximum Residual = ',
     &        RSD(IEQHN),': Node = ',NX,
     &        ': Hydrate Saturation = ',SH(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Hydrate-N2 Mass Equation Maximum Residual = ',
     &        RSD(IEQHN),': Node = ',NX,
     &        'Hydrate Saturation = ',SH(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
          ENDIF
        ENDIF
!
!---    Isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          NX = NSD(IEQS)
          IF( NX.GT.0 ) THEN
            NPX = MOD(NPHAZ(2,NX),100)
            NCHX = INDEX( PH_CND(NPX),'  ') - 1
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Salt Equation Maximum Residual = ',
     &        RSD(IEQS),': Node = ',NX,
     &        ': Salt Mass Fraction = ',YLS(2,NX),
     &        ': Phase Condition = ',PH_CND(NPX)(1:NCHX)
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,1PE11.4,2A)')
     &        'Salt Equation Maximum Residual = ',
     &        RSD(IEQS),': Node = ',NX,
     &        ': Salt Mass Fraction = ',YLS(2,NX),
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
          DO 400 N = 1,NFLD
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            PSO(2,N) = PSO(1,N)
            ZMCA(2,N) = ZMCA(1,N)
            ZMCO(2,N) = ZMCO(1,N)
            ZMCN(2,N) = ZMCN(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            POSM(2,N) = POSM(1,N)
            SH(2,N) = SH(1,N)
            SI(2,N) = SI(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            YMGA(2,N) = YMGA(1,N)
            YMGO(2,N) = YMGO(1,N)
            YMGN(2,N) = YMGN(1,N)
            YMHGA(2,N) = YMHGA(1,N)
            YMHGO(2,N) = YMHGO(1,N)
            YMHGN(2,N) = YMHGN(1,N)
            TMHA(2,N) = TMHA(1,N)
            TMHO(2,N) = TMHO(1,N)
            TMHN(2,N) = TMHN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  400     CONTINUE
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
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL_HYD_O2_CH4_N2 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SBND_HYD_O2_CH4_N2( NSL )
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Modify the Jacobian matrix for the solute transport equation
!     to incorporate boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
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
#ifdef petsc
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#endif
#ifdef petsc
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: IROW
      PetscErrorCode :: IERR
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SBND_HYD_O2_CH4_N2'
!
!---  Loop over number of specified boundary conditions  ---
!
      DO 200 NB = 1,NBC
!
!---    Zero flux boundary condition  ---
!
        IF( IBCT(NSL+LUK,NB).EQ.3 ) GOTO 200
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 200
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) GOTO 200
        IF( IBCM(NB).EQ.1 ) THEN
          BCX = BC(NSL+LBCU,1,MB)
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
             TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
             DTBC = MIN( BC(1,M,MB)-TMZ,DT )
             TFBC = (TMZ-5.D-1*DTBC-BC(1,M-1,MB))/TDBC
             BCX = BC(NSL+LBCU,M-1,MB) +
     &         TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
             IF( IBCT(NSL+LUK,NB).EQ.12 ) BCX = CBO(NB,NSL)
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 200
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
#ifdef petsc
        ELSEIF( ILES.EQ.5 ) THEN
          IROW = petsc_offset(N)
#endif
        ENDIF
!
!---  Diffusion coefficients at node adjacent to boundary  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SMDLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
        DLP = TORL(2,N)*SL(2,N)*PORD(2,N)*SMDLP
        PCOR = (PG(2,N)+PATM)/PATM
        SMDGP = SMDG(NSL)*(TCOR**1.75)/PCOR
        DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SMDGP
!
!---  Phase fraction factors at node adjacent to boundary  ---
!
        FCLP = 0.D+0
        IF( SL(2,N).GT.SMALL ) FCLP = YL(N,NSL)/(SL(2,N)*PORD(2,N))
        FCGP = 0.D+0
        IF( SG(2,N).GT.SMALL ) FCGP = YG(N,NSL)/(SG(2,N)*PORD(2,N))
!
!---  Phase fraction factors at boundary  ---
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
!---  Phase-volumetric concentration ratios  ---
!
        FCL = 1.D+0/(XVSB + XVLB + XVGB*PCGLX)
        FCG = 1.D+0/((XVSB + XVLB)/PCGLX + XVGB)
!
!---  Phase mole fractions  ---
!
        YLB(NB,NSL) = XVLB*FCL
        YGB(NB,NSL) = XVGB*FCG
!
!---  Convert boundary phase concentrations to
!     volumetric concentrations  ---
!
        IF( IBCT(NSL+LUK,NB).EQ.8 ) THEN
          BCX = BCX/FCL
        ELSEIF( IBCT(NSL+LUK,NB).EQ.9 ) THEN
          BCX = BCX/FCG
        ENDIF
        CB(NB,NSL) = BCX
!
!---  Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
          NPZ = NSZ(N)
!
!---      Hydraulic dispersion  ---
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
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR.
     &      IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
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
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              IF( FLB.GE.ZERO ) THEN
                WCLZ = BCX*FCL*FLB
              ELSEIF( FLB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCLT = YL(NBT,NSL)/(SL(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBT,NSL)*FCLT)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                WCLZ = BCX*FLB*THETA*FCL
     &            + C(N,NSL)*FLB*(1.D+0-THETA)*FCLP
              ELSEIF( FLB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCLZ = C(N,NSL)*FLB*FCLP
              ENDIF
              IF( FGB.GE.ZERO ) THEN
                WCGZ = BCX*FCG*FGB
              ELSEIF( FGB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCGT = YG(NBT,NSL)/(SG(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBT,NSL)*FCGT)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                WCGZ = BCX*FGB*THETA*FCG
     &            + C(N,NSL)*FGB*(1.D+0-THETA)*FCGP
              ELSEIF( FGB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCGZ = C(N,NSL)*FGB*FCGP
              ENDIF
              AB = DLZ*FCL + DGZ*FCG
              AP = DGZ*FCGP + DLZ*FCLP
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCGZ)/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCGZ
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVLX = SL(2,NBT)*PORD(2,NBT)
                FCLT = YL(NBT,NSL)/(XVLX+SMALL)
                CRLT = ABS( WL(1,NQZ) )*DT/(DZGP(NQZ)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBT,NSL)*FCLT-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCLZ = C(N,NSL)*FLT*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBT,NSL)*FLT*THETA*DZF*FCLT
                WCLZF = CO(N,NSL)*FLT*FCLP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ-WCLZF)/AFZ(NQZ)
                BLU(MP) = BLU(MP) - WCLZ + WCLZF
                BLU(IXP(NBT)) = BLU(IXP(NBT)) + WCLZ - WCLZF
              ENDIF
              FGT = AFZ(NQZ)*WG(1,NQZ)
              IF( FGT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVGX = SG(2,NBT)*PORD(2,NBT)
                FCGT = YG(NBT,NSL)/(XVGX+SMALL)
                CRGT = ABS( WG(1,NQZ) )*DT/(DZGP(NQZ)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBT,NSL)*FCGT-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCGZ = C(N,NSL)*FGT*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBT,NSL)*FGT*THETA*DZF*FCGT
                WCGZF = CO(N,NSL)*FGT*FCGP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCGZ-WCGZF)/AFZ(NQZ)
                BLU(MP) = BLU(MP) - WCGZ + WCGZF
                BLU(IXP(NBT)) = BLU(IXP(NBT)) + WCGZ - WCGZF
              ENDIF
            ELSE
              ALB = MAX( FLB,ZERO ) +
     &          DLZ*MAX((ONE-(TENTH*ABS(FLB)/(DLZ+SMALL)))**5,ZERO)
              AGB = MAX( FGB,ZERO ) +
     &          DGZ*MAX((ONE-(TENTH*ABS(FGB)/(DGZ+SMALL)))**5,ZERO)
              AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP
              AB = ALB*FCL + AGB*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AB*BCX
!
!---  Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLB = MIN( FLB,0.D+0 )
            FGB = MIN( FGB,0.D+0 )
            FNB = MIN( FNB,0.D+0 )
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              WCLZ = 0.D+0
              IF( FLB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCLT = YL(NBT,NSL)/(SL(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBT,NSL)*FCLT)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                WCLZ = BCX*FLB*THETA*FCL
     &            + C(N,NSL)*FLB*(1.D+0-THETA)*FCLP
              ELSEIF( FLB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCLZ = C(N,NSL)*FLB*FCLP
              ENDIF
              WCGZ = 0.D+0
              IF( FGB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCGT = YG(NBT,NSL)/(SG(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBT,NSL)*FCGT)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                WCGZ = BCX*FGB*THETA*FCG
     &            + C(N,NSL)*FGB*(1.D+0-THETA)*FCGP
              ELSEIF( FGB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCGZ = C(N,NSL)*FGB*FCGP
              ENDIF
              AB = 0.D+0
              AP = 0.D+0
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCGZ)/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCGZ
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVLX = SL(2,NBT)*PORD(2,NBT)
                FCLT = YL(NBT,NSL)/(XVLX+SMALL)
                CRLT = ABS( WL(1,NQZ) )*DT/(DZGP(NQZ)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBT,NSL)*FCLT-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCLZ = C(N,NSL)*FLT*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBT,NSL)*FLT*THETA*DZF*FCLT
                WCLZF = CO(N,NSL)*FLT*FCLP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ-WCLZF)/AFZ(NQZ)
                BLU(MP) = BLU(MP) - WCLZ + WCLZF
                BLU(IXP(NBT)) = BLU(IXP(NBT)) + WCLZ - WCLZF
              ENDIF
              FGT = AFZ(NQZ)*WG(1,NQZ)
              IF( FGT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVGX = SG(2,NBT)*PORD(2,NBT)
                FCGT = YG(NBT,NSL)/(XVGX+SMALL)
                CRGT = ABS( WG(1,NQZ) )*DT/(DZGP(NQZ)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBT,NSL)*FCGT-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCGZ = C(N,NSL)*FGT*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBT,NSL)*FGT*THETA*DZF*FCGT
                WCGZF = CO(N,NSL)*FGT*FCGP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCGZ-WCGZF)/AFZ(NQZ)
                BLU(MP) = BLU(MP) - WCGZ + WCGZF
                BLU(IXP(NBT)) = BLU(IXP(NBT)) + WCGZ - WCGZF
              ENDIF
            ELSE
              ALB = MAX( FLB,ZERO )
              AGB = MAX( FGB,ZERO )
              AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP
              AB = ALB*FCL + AGB*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AB*BCX
!
!--- Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 .OR.
     &      IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLB = MAX( FLB,0.D+0 )
            FGB = MAX( FGB,0.D+0 )
            FNB = MAX( FNB,0.D+0 )
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              WCLZ = 0.D+0
              IF( FLB.GE.ZERO ) WCLZ = BCX*FCL*FLB
              WCGZ = 0.D+0
              IF( FGB.GE.ZERO ) WCGZ = BCX*FCG*FGB
              AB = 0.D+0
              AP = 0.D+0
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCGZ)/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCGZ
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVLX = SL(2,NBT)*PORD(2,NBT)
                FCLT = YL(NBT,NSL)/(XVLX+SMALL)
                CRLT = ABS( WL(1,NQZ) )*DT/(DZGP(NQZ)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBT,NSL)*FCLT-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCLZ = C(N,NSL)*FLT*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBT,NSL)*FLT*THETA*DZF*FCLT
                WCLZF = CO(N,NSL)*FLT*FCLP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ-WCLZF)/AFZ(NQZ)
                BLU(MP) = BLU(MP) - WCLZ + WCLZF
                BLU(IXP(NBT)) = BLU(IXP(NBT)) + WCLZ - WCLZF
              ENDIF
              FGT = AFZ(NQZ)*WG(1,NQZ)
              IF( FGT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVGX = SG(2,NBT)*PORD(2,NBT)
                FCGT = YG(NBT,NSL)/(XVGX+SMALL)
                CRGT = ABS( WG(1,NQZ) )*DT/(DZGP(NQZ)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBT,NSL)*FCGT-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCGZ = C(N,NSL)*FGT*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBT,NSL)*FGT*THETA*DZF*FCGT
                WCGZF = CO(N,NSL)*FGT*FCGP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCGZ-WCGZF)/AFZ(NQZ)
                BLU(MP) = BLU(MP) - WCGZ + WCGZF
                BLU(IXP(NBT)) = BLU(IXP(NBT)) + WCGZ - WCGZF
              ENDIF
            ELSE
              ALB = MAX( FLB,ZERO )
              AGB = MAX( FGB,ZERO )
              AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP
              AB = ALB*FCL + AGB*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AB*BCX
          ENDIF
!
!---  South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          NPY = NSY(N)
!
!---      Hydraulic dispersion  ---
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
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR.
     &      IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
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
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              IF( FLS.GE.ZERO ) THEN
                VCLY = BCX*FCL*FLS
              ELSEIF( FLS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCLN = YL(NBN,NSL)/(SL(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBN,NSL)*FCLN)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                VCLY = BCX*FLS*THETA*FCL
     &            + C(N,NSL)*FLS*(1.D+0-THETA)*FCLP
              ELSEIF( FLS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                 VCLY = C(N,NSL)*FLS*FCLP
             ENDIF
              IF( FGS.GE.ZERO ) THEN
                VCGY = BCX*FCG*FGS
              ELSEIF( FGS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCGN = YG(NBN,NSL)/(SG(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBN,NSL)*FCGN)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                VCGY = BCX*FGS*THETA*FCG
     &            + C(N,NSL)*FGS*(1.D+0-THETA)*FCGP
              ELSEIF( FGS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCGY = C(N,NSL)*FGS*FCGP
              ENDIF
              AS = DLY*FCL + DGY*FCG
              AP = DGY*FCGP + DLY*FCLP
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCGY)/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCGY
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQY = NPY+IFLD
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVLX = SL(2,NBN)*PORD(2,NBN)
                FCLN = YL(NBN,NSL)/(XVLX+SMALL)
                CRLN = ABS( VL(1,NQY) )*DT/(DYGP(NQY)*RP(I)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBN,NSL)*FCLN-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCLY = C(N,NSL)*FLN*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBN,NSL)*FLN*THETA*DYF*FCLN
                VCLYF = CO(N,NSL)*FLN*FCLP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCLY-VCLYF)/AFY(NQY)
                BLU(MP) = BLU(MP) - VCLY + VCLYF
                BLU(IXP(NBN)) = BLU(IXP(NBN)) + VCLY - VCLYF
              ENDIF
              FGN = AFY(NQY)*VG(1,NQY)
              IF( FGN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVGX = SG(2,NBN)*PORD(2,NBN)
                FCGN = YG(NBN,NSL)/(XVGX+SMALL)
                CRGN = ABS( VG(1,NQY) )*DT/(DYGP(NQY)*RP(I)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBN,NSL)*FCGN-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCGY = C(N,NSL)*FGN*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBN,NSL)*FGN*THETA*DYF*FCGN
                VCGYF = CO(N,NSL)*FGN*FCGP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCGY-VCGYF)/AFY(NQY)
                BLU(MP) = BLU(MP) - VCGY + VCGYF
                BLU(IXP(NBN)) = BLU(IXP(NBN)) + VCGY - VCGYF
              ENDIF
            ELSE
              ALS = MAX( FLS,ZERO ) +
     &          DLY*MAX((ONE-(TENTH*ABS(FLS)/(DLY+SMALL)))**5,ZERO)
              AGS = MAX( FGS,ZERO ) +
     &          DGY*MAX((ONE-(TENTH*ABS(FGS)/(DGY+SMALL)))**5,ZERO)
              AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP
              AS = ALS*FCL + AGS*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AS*BCX
!
!---  Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLS = MIN( FLS,0.D+0 )
            FGS = MIN( FGS,0.D+0 )
            FNS = MIN( FNS,0.D+0 )
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              VCLY = 0.D+0
              IF( FLS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCLN = YL(NBN,NSL)/(SL(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBN,NSL)*FCLN)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                VCLY = BCX*FLS*THETA*FCL
     &            + C(N,NSL)*FLS*(1.D+0-THETA)*FCLP
              ELSEIF( FLS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                 VCLY = C(N,NSL)*FLS*FCLP
              ENDIF
              VCGY = 0.D+0
              IF( FGS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCGN = YG(NBN,NSL)/(SG(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBN,NSL)*FCGN)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                VCGY = BCX*FGS*THETA*FCG
     &            + C(N,NSL)*FGS*(1.D+0-THETA)*FCGP
              ELSEIF( FGS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCGY = C(N,NSL)*FGS*FCGP
              ENDIF
              AS = 0.D+0
              AP = 0.D+0
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCGY)/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCGY
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQY = NPY+IFLD
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVLX = SL(2,NBN)*PORD(2,NBN)
                FCLN = YL(NBN,NSL)/(XVLX+SMALL)
                CRLN = ABS( VL(1,NQY) )*DT/(DYGP(NQY)*RP(I)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBN,NSL)*FCLN-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCLY = C(N,NSL)*FLN*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBN,NSL)*FLN*THETA*DYF*FCLN
                VCLYF = CO(N,NSL)*FLN*FCLP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCLY-VCLYF)/AFY(NQY)
                BLU(MP) = BLU(MP) - VCLY + VCLYF
                BLU(IXP(NBN)) = BLU(IXP(NBN)) + VCLY - VCLYF
              ENDIF
              FGN = AFY(NQY)*VG(1,NQY)
              IF( FGN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVGX = SG(2,NBN)*PORD(2,NBN)
                FCGN = YG(NBN,NSL)/(XVGX+SMALL)
                CRGN = ABS( VG(1,NQY) )*DT/(DYGP(NQY)*RP(I)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBN,NSL)*FCGN-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCGY = C(N,NSL)*FGN*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBN,NSL)*FGN*THETA*DYF*FCGN
                VCGYF = CO(N,NSL)*FGN*FCGP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCGY-VCGYF)/AFY(NQY)
                BLU(MP) = BLU(MP) - VCGY + VCGYF
                BLU(IXP(NBN)) = BLU(IXP(NBN)) + VCGY - VCGYF
              ENDIF
            ELSE
              ALS = MAX( FLS,ZERO )
              AGS = MAX( FGS,ZERO )
              AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP
              AS = ALS*FCL + AGS*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AS*BCX
!
!--- Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 .OR.
     &      IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLS = MAX( FLS,0.D+0 )
            FGS = MAX( FGS,0.D+0 )
            FNS = MAX( FNS,0.D+0 )
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              VCLY = 0.D+0
              IF( FLS.GE.ZERO ) VCLY = BCX*FCL*FLS
              VCGY = 0.D+0
              IF( FGS.GE.ZERO ) VCGY = BCX*FCG*FGS
              AS = 0.D+0
              AP = 0.D+0
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCGY)/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCGY
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQY = NPY+IFLD
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVLX = SL(2,NBN)*PORD(2,NBN)
                FCLN = YL(NBN,NSL)/(XVLX+SMALL)
                CRLN = ABS( VL(1,NQY) )*DT/(DYGP(NQY)*RP(I)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBN,NSL)*FCLN-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCLY = C(N,NSL)*FLN*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBN,NSL)*FLN*THETA*DYF*FCLN
                VCLYF = CO(N,NSL)*FLN*FCLP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCLY-VCLYF)/AFY(NQY)
                BLU(MP) = BLU(MP) - VCLY + VCLYF
                BLU(IXP(NBN)) = BLU(IXP(NBN)) + VCLY - VCLYF
              ENDIF
              FGN = AFY(NQY)*VG(1,NQY)
              IF( FGN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVGX = SG(2,NBN)*PORD(2,NBN)
                FCGN = YG(NBN,NSL)/(XVGX+SMALL)
                CRGN = ABS( VG(1,NQY) )*DT/(DYGP(NQY)*RP(I)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBN,NSL)*FCGN-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCGY = C(N,NSL)*FGN*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBN,NSL)*FGN*THETA*DYF*FCGN
                VCGYF = CO(N,NSL)*FGN*FCGP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCGY-VCGYF)/AFY(NQY)
                BLU(MP) = BLU(MP) - VCGY + VCGYF
                BLU(IXP(NBN)) = BLU(IXP(NBN)) + VCGY - VCGYF
              ENDIF
            ELSE
              ALS = MAX( FLS,ZERO )
              AGS = MAX( FGS,ZERO )
              AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP
              AS = ALS*FCL + AGS*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AS*BCX
          ENDIF
!
!---  West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          NPX = NSX(N)
!
!---      Hydraulic dispersion  ---
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
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR.
     &      IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
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
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              IF( FLW.GE.ZERO ) THEN
                UCLX = BCX*FCL*FLW
              ELSEIF( FLW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCLE = YL(NBE,NSL)/(SL(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBE,NSL)*FCLE)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA)*FCLP
     &            + BCX*FLW*THETA*FCL
              ELSEIF( FLW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCLX = C(N,NSL)*FLW*FCLP
              ENDIF
              IF( FGW.GE.ZERO ) THEN
                UCGX = BCX*FCG*FGW
              ELSEIF( FGW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCGE = YG(NBE,NSL)/(SG(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBE,NSL)*FCGE)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA)*FCGP
     &            + BCX*FGW*THETA*FCG
              ELSEIF( FGW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCGX = C(N,NSL)*FGW*FCGP
              ENDIF
              AW = DLX*FCL + DGX*FCG
              AP = DLX*FCLP + DGX*FCGP
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCGX)/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCGX
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQX = NPX+1
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                XVLX = SL(2,NBE)*PORD(2,NBE)
                FCLE = YL(NBE,NSL)/(XVLX+SMALL)
                CRLE = ABS( UL(1,NQX) )*DT/(DXGP(NQX)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBE,NSL)*FCLE-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCLX = C(N,NSL)*FLE*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBE,NSL)*FLE*THETA*DXF*FCLE
                UCLXF = CO(N,NSL)*FLE*FCLP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCLX-UCLXF)/AFX(NQX)
                BLU(MP) = BLU(MP) - UCLX + UCLXF
                BLU(IXP(NBE)) = BLU(IXP(NBE)) + UCLX - UCLXF
              ENDIF
              FGE = AFX(NQX)*UG(1,NQX)
              IF( FGE.GE.ZERO ) THEN
                NBE = N+1
                XVGX = SG(2,NBE)*PORD(2,NBE)
                FCGE = YG(NBE,NSL)/(XVGX+SMALL)
                CRGE = ABS( UG(1,NQX) )*DT/(DXGP(NQX)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBE,NSL)*FCGE-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCGX = C(N,NSL)*FGE*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBE,NSL)*FGE*THETA*DXF*FCGE
                UCGXF = CO(N,NSL)*FGE*FCGP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCGX-UCGXF)/AFX(NQX)
                BLU(MP) = BLU(MP) - UCGX + UCGXF
                BLU(IXP(NBE)) = BLU(IXP(NBE)) + UCGX - UCGXF
              ENDIF
            ELSE
              ALW = MAX(FLW,ZERO)
     &          + DLX*MAX((ONE-(TENTH*ABS(FLW)/(DLX+SMALL)))**5,ZERO)
              AGW = MAX(FGW,ZERO)
     &          + DGX*MAX((ONE-(TENTH*ABS(FGW)/(DGX+SMALL)))**5,ZERO)
              AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
              AW = ALW*FCL + AGW*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AW*BCX
!
!---  Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLW = MIN( FLW,0.D+0 )
            FGW = MIN( FGW,0.D+0 )
            FNW = MIN( FNW,0.D+0 )
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              UCLX = 0.D+0
              IF( FLW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCLE = YL(NBE,NSL)/(SL(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBE,NSL)*FCLE)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA)*FCLP
     &            + BCX*FLW*THETA*FCL
              ELSEIF( FLW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCLX = C(N,NSL)*FLW*FCLP
              ENDIF
              UCGX = 0.D+0
              IF( FGW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCGE = YG(NBE,NSL)/(SG(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBE,NSL)*FCGE)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA)*FCGP
     &            + BCX*FGW*THETA*FCG
              ELSEIF( FGW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCGX = C(N,NSL)*FGW*FCGP
              ENDIF
              AW = 0.D+0
              AP = 0.D+0
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCGX)/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCGX
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQX = NPX+1
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                XVLX = SL(2,NBE)*PORD(2,NBE)
                FCLE = YL(NBE,NSL)/(XVLX+SMALL)
                CRLE = ABS( UL(1,NQX) )*DT/(DXGP(NQX)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBE,NSL)*FCLE-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCLX = C(N,NSL)*FLE*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBE,NSL)*FLE*THETA*DXF*FCLE
                UCLXF = CO(N,NSL)*FLE*FCLP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCLX-UCLXF)/AFX(NQX)
                BLU(MP) = BLU(MP) - UCLX + UCLXF
                BLU(IXP(NBE)) = BLU(IXP(NBE)) + UCLX - UCLXF
              ENDIF
              FGE = AFX(NQX)*UG(1,NQX)
              IF( FGE.GE.ZERO ) THEN
                NBE = N+1
                XVGX = SG(2,NBE)*PORD(2,NBE)
                FCGE = YG(NBE,NSL)/(XVGX+SMALL)
                CRGE = ABS( UG(1,NQX) )*DT/(DXGP(NQX)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBE,NSL)*FCGE-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCGX = C(N,NSL)*FGE*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBE,NSL)*FGE*THETA*DXF*FCGE
                UCGXF = CO(N,NSL)*FGE*FCGP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCGX-UCGXF)/AFX(NQX)
                BLU(MP) = BLU(MP) - UCGX + UCGXF
                BLU(IXP(NBE)) = BLU(IXP(NBE)) + UCGX - UCGXF
              ENDIF
            ELSE
              ALW = MAX(FLW,ZERO)
              AGW = MAX(FGW,ZERO)
              AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
              AW = ALW*FCL + AGW*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AW*BCX
!
!--- Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 .OR.
     &      IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLW = MAX( FLW,0.D+0 )
            FGW = MAX( FGW,0.D+0 )
            FNW = MAX( FNW,0.D+0 )
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              UCLX = 0.D+0
              IF( FLW.GE.ZERO ) UCLX = BCX*FCL*FLW
              UCGX = 0.D+0
              IF( FGW.GE.ZERO ) UCGX = BCX*FCG*FGW
              AW = 0.D+0
              AP = 0.D+0
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCGX)/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCGX
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQX = NPX+1
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                XVLX = SL(2,NBE)*PORD(2,NBE)
                FCLE = YL(NBE,NSL)/(XVLX+SMALL)
                CRLE = ABS( UL(1,NQX) )*DT/(DXGP(NQX)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBE,NSL)*FCLE-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCLX = C(N,NSL)*FLE*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBE,NSL)*FLE*THETA*DXF*FCLE
                UCLXF = CO(N,NSL)*FLE*FCLP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCLX-UCLXF)/AFX(NQX)
                BLU(MP) = BLU(MP) - UCLX + UCLXF
                BLU(IXP(NBE)) = BLU(IXP(NBE)) + UCLX - UCLXF
              ENDIF
              FGE = AFX(NQX)*UG(1,NQX)
              IF( FGE.GE.ZERO ) THEN
                NBE = N+1
                XVGX = SG(2,NBE)*PORD(2,NBE)
                FCGE = YG(NBE,NSL)/(XVGX+SMALL)
                CRGE = ABS( UG(1,NQX) )*DT/(DXGP(NQX)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBE,NSL)*FCGE-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCGX = C(N,NSL)*FGE*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBE,NSL)*FGE*THETA*DXF*FCGE
                UCGXF = CO(N,NSL)*FGE*FCGP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCGX-UCGXF)/AFX(NQX)
                BLU(MP) = BLU(MP) - UCGX + UCGXF
                BLU(IXP(NBE)) = BLU(IXP(NBE)) + UCGX - UCGXF
              ENDIF
            ELSE
              ALW = MAX(FLW,ZERO)
              AGW = MAX(FGW,ZERO)
              AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
              AW = ALW*FCL + AGW*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AW*BCX
          ENDIF
!
!---  East boundary
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          NQX = NSX(N) + 1
!
!---      Hydraulic dispersion  ---
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
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR.
     &      IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
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
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              IF( FLE.LT.ZERO ) THEN
                UCLX = BCX*FCL*FLE
              ELSEIF( FLE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCLW = YL(NBW,NSL)/(SL(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBW,NSL)*FCLW)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                UCLX =  C(N,NSL)*FLE*(1.D+0-THETA)*FCLP
     &          +  BCX*FLE*THETA*FCL
              ELSEIF( FLE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCLX =  C(N,NSL)*FLE*FCLP
              ENDIF
              IF( FGE.LT.ZERO ) THEN
                UCGX = BCX*FCG*FGE
              ELSEIF( FGE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCGW = YG(NBW,NSL)/(SG(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBW,NSL)*FCGW)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                UCGX =  C(N,NSL)*FGE*(1.D+0-THETA)*FCGP
     &          +  BCX*FGE*THETA*FCG
              ELSEIF( FGE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCGX =  C(N,NSL)*FGE*FCGP
              ENDIF
              AE = DLX*FCL + DGX*FCG
              AP = DLX*FCLP + DGX*FCGP
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCGX)/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCGX
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPX = NSX(N)
              FLW = AFX(NPX)*UL(1,NPX)
              IF( FLW.LT.ZERO ) THEN
                NBW = N-1
                XVLX = SL(2,NBW)*PORD(2,NBW)
                CRLW = ABS( UL(1,NPX) )*DT/(DXGP(NPX)*XVLX+SMALL)
                FCLW = YL(NBW,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBW,NSL)*FCLW-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBW,NSL)*FLW*THETA*DXF*FCLW
                UCLXF = CO(N,NSL)*FLW*FCLP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCLX-UCLXF)/AFX(NPX)
                BLU(MP) = BLU(MP) + UCLX - UCLXF
                BLU(IXP(NBW)) = BLU(IXP(NBW)) - UCLX + UCLXF
              ENDIF
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/(DXGP(NPX)*XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBW,NSL)*FCGW-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBW,NSL)*FGW*THETA*DXF*FCGW
                UCGXF = CO(N,NSL)*FGW*FCGP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCGX-UCGXF)/AFX(NPX)
                BLU(MP) = BLU(MP) + UCGX - UCGXF
                BLU(IXP(NBW)) = BLU(IXP(NBW)) - UCGX + UCGXF
              ENDIF
            ELSE
              ALE = MAX( -FLE,ZERO ) +
     &          DLX*MAX((ONE-(TENTH*ABS(FLE)/(DLX+SMALL)))**5,ZERO)
              AGE = MAX( -FGE,ZERO ) +
     &          DGX*MAX((ONE-(TENTH*ABS(FGE)/(DGX+SMALL)))**5,ZERO)
              AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP
              AE = ALE*FCL + AGE*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AE*BCX
!
!---  Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLE = MAX( FLE,0.D+0 )
            FGE = MAX( FGE,0.D+0 )
            FNE = MAX( FNE,0.D+0 )
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              UCLX = 0.D+0
              IF( FLE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCLW = YL(NBW,NSL)/(SL(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBW,NSL)*FCLW)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                UCLX =  C(N,NSL)*FLE*(1.D+0-THETA)*FCLP
     &          +  BCX*FLE*THETA*FCL
              ELSEIF( FLE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCLX =  C(N,NSL)*FLE*FCLP
              ENDIF
              UCGX = 0.D+0
              IF( FGE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCGW = YG(NBW,NSL)/(SG(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBW,NSL)*FCGW)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                UCGX =  C(N,NSL)*FGE*(1.D+0-THETA)*FCGP
     &          +  BCX*FGE*THETA*FCG
              ELSEIF( FGE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCGX =  C(N,NSL)*FGE*FCGP
              ENDIF
              AE = 0.D+0
              AP = 0.D+0
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCGX)/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCGX
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPX = NSX(N)
              FLW = AFX(NPX)*UL(1,NPX)
              IF( FLW.LT.ZERO ) THEN
                NBW = N-1
                XVLX = SL(2,NBW)*PORD(2,NBW)
                CRLW = ABS( UL(1,NPX) )*DT/(DXGP(NPX)*XVLX+SMALL)
                FCLW = YL(NBW,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBW,NSL)*FCLW-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBW,NSL)*FLW*THETA*DXF*FCLW
                UCLXF = CO(N,NSL)*FLW*FCLP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCLX-UCLXF)/AFX(NPX)
                BLU(MP) = BLU(MP) + UCLX - UCLXF
                BLU(IXP(NBW)) = BLU(IXP(NBW)) - UCLX + UCLXF
              ENDIF
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/(DXGP(NPX)*XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBW,NSL)*FCGW-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBW,NSL)*FGW*THETA*DXF*FCGW
                UCGXF = CO(N,NSL)*FGW*FCGP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCGX-UCGXF)/AFX(NPX)
                BLU(MP) = BLU(MP) + UCGX - UCGXF
                BLU(IXP(NBW)) = BLU(IXP(NBW)) - UCGX + UCGXF
              ENDIF
            ELSE
              ALE = MAX( -FLE,ZERO )
              AGE = MAX( -FGE,ZERO )
              AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP
              AE = ALE*FCL + AGE*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AE*BCX
!
!--- Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 .OR.
     &      IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLE = MIN( FLE,0.D+0 )
            FGE = MIN( FGE,0.D+0 )
            FNE = MIN( FNE,0.D+0 )
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              UCLX = 0.D+0
              IF( FLE.LT.ZERO ) UCLX = BCX*FCL*FLE
              UCGX = 0.D+0
              IF( FGE.LT.ZERO ) UCGX = BCX*FCG*FGE
              AE = 0.D+0
              AP = 0.D+0
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCGX)/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCGX
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPX = NSX(N)
              FLW = AFX(NPX)*UL(1,NPX)
              IF( FLW.LT.ZERO ) THEN
                NBW = N-1
                XVLX = SL(2,NBW)*PORD(2,NBW)
                CRLW = ABS( UL(1,NPX) )*DT/(DXGP(NPX)*XVLX+SMALL)
                FCLW = YL(NBW,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBW,NSL)*FCLW-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBW,NSL)*FLW*THETA*DXF*FCLW
                UCLXF = CO(N,NSL)*FLW*FCLP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCLX-UCLXF)/AFX(NPX)
                BLU(MP) = BLU(MP) + UCLX - UCLXF
                BLU(IXP(NBW)) = BLU(IXP(NBW)) - UCLX + UCLXF
              ENDIF
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/(DXGP(NPX)*XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBW,NSL)*FCGW-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBW,NSL)*FGW*THETA*DXF*FCGW
                UCGXF = CO(N,NSL)*FGW*FCGP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCGX-UCGXF)/AFX(NPX)
                BLU(MP) = BLU(MP) + UCGX - UCGXF
                BLU(IXP(NBW)) = BLU(IXP(NBW)) - UCGX + UCGXF
              ENDIF
            ELSE
              ALE = MAX( -FLE,ZERO )
              AGE = MAX( -FGE,ZERO )
              AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP
              AE = ALE*FCL + AGE*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AE*BCX
          ENDIF
!
!---  North boundary  ---
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
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR.
     &      IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
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
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              IF( FLN.LT.ZERO ) THEN
                VCLY = BCX*FCL*FLN
              ELSEIF( FLN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCLS = YL(NBS,NSL)/(SL(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBS,NSL)*FCLS)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                VCLY =  BCX*FLN*THETA*FCL
     &            + C(N,NSL)*FLN*(1.D+0-THETA)*FCLP
              ELSEIF( FLN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCLY =  C(N,NSL)*FLN*FCLP
              ENDIF
              IF( FGN.LT.ZERO ) THEN
                VCGY = BCX*FCG*FGN
              ELSEIF( FGN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCGS = YG(NBS,NSL)/(SG(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBS,NSL)*FCGS)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                VCGY =  BCX*FGN*THETA*FCG
     &            + C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ELSEIF( FGN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCGY =  C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ENDIF
              AN = DLY*FCL + DGY*FCG
              AP = DLY*FCLP + DGY*FCGP
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCGY)/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCGY
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPY = NSY(N)
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/(DYGP(NPY)*XVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBS,NSL)*FCLS-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCLY = C(N,NSL)*FLS*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBS,NSL)*FLS*THETA*DYF*FCLS
                VCLYF = CO(N,NSL)*FLS*FCLP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCLY-VCLYF)/AFY(NPY)
                BLU(MP) = BLU(MP) + VCLY - VCLYF
                BLU(IXP(NBS)) = BLU(IXP(NBS)) - VCLY + VCLYF
              ENDIF
              FGS = AFY(NPY)*VG(1,NPY)
              IF( FGS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/(DYGP(NPY)*XVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBS,NSL)*FCGS-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCGY = C(N,NSL)*FGS*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBS,NSL)*FGS*THETA*DYF*FCGS
                VCGYF = CO(N,NSL)*FGS*FCGP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCGY-VCGYF)/AFY(NPY)
                BLU(MP) = BLU(MP) + VCGY - VCGYF
                BLU(IXP(NBS)) = BLU(IXP(NBS)) - VCGY + VCGYF
              ENDIF
            ELSE
              ALN = MAX( -FLN,ZERO ) +
     &          DLY*MAX((ONE-(TENTH*ABS(FLN)/(DLY+SMALL)))**5,ZERO)
              AGN = MAX( -FGN,ZERO ) +
     &          DGY*MAX((ONE-(TENTH*ABS(FGN)/(DGY+SMALL)))**5,ZERO)
              AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP
              AN = ALN*FCL + AGN*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AN*BCX
!
!---  Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLN = MAX( FLN,0.D+0 )
            FGN = MAX( FGN,0.D+0 )
            FNN = MAX( FNN,0.D+0 )
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              VCLY = 0.D+0
              IF( FLN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCLS = YL(NBS,NSL)/(SL(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBS,NSL)*FCLS)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                VCLY =  BCX*FLN*THETA*FCL
     &            + C(N,NSL)*FLN*(1.D+0-THETA)*FCLP
              ELSEIF( FLN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCLY =  C(N,NSL)*FLN*FCLP
              ENDIF
              VCGY = 0.D+0
              IF( FGN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCGS = YG(NBS,NSL)/(SG(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBS,NSL)*FCGS)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                VCGY =  BCX*FGN*THETA*FCG
     &            + C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ELSEIF( FGN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCGY =  C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ENDIF
              AN = 0.D+0
              AP = 0.D+0
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCGY)/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCGY
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPY = NSY(N)
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/(DYGP(NPY)*XVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBS,NSL)*FCLS-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCLY = C(N,NSL)*FLS*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBS,NSL)*FLS*THETA*DYF*FCLS
                VCLYF = CO(N,NSL)*FLS*FCLP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCLY-VCLYF)/AFY(NPY)
                BLU(MP) = BLU(MP) + VCLY - VCLYF
                BLU(IXP(NBS)) = BLU(IXP(NBS)) - VCLY + VCLYF
              ENDIF
              FGS = AFY(NPY)*VG(1,NPY)
              IF( FGS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/(DYGP(NPY)*XVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBS,NSL)*FCGS-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCGY = C(N,NSL)*FGS*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBS,NSL)*FGS*THETA*DYF*FCGS
                VCGYF = CO(N,NSL)*FGS*FCGP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCGY-VCGYF)/AFY(NPY)
                BLU(MP) = BLU(MP) + VCGY - VCGYF
                BLU(IXP(NBS)) = BLU(IXP(NBS)) - VCGY + VCGYF
              ENDIF
            ELSE
              ALN = MAX( -FLN,ZERO )
              AGN = MAX( -FGN,ZERO )
              AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP
              AN = ALN*FCL + AGN*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AN*BCX
!
!--- Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 .OR.
     &      IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLN = MIN( FLN,0.D+0 )
            FGN = MIN( FGN,0.D+0 )
            FNN = MIN( FNN,0.D+0 )
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              VCLY = 0.D+0
              IF( FLN.LT.ZERO ) VCLY = BCX*FCL*FLN
              VCGY = 0.D+0
              IF( FGN.LT.ZERO ) VCGY = BCX*FCG*FGN
              AN = 0.D+0
              AP = 0.D+0
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCGY)/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCGY
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPY = NSY(N)
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/(DYGP(NPY)*XVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBS,NSL)*FCLS-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCLY = C(N,NSL)*FLS*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBS,NSL)*FLS*THETA*DYF*FCLS
                VCLYF = CO(N,NSL)*FLS*FCLP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCLY-VCLYF)/AFY(NPY)
                BLU(MP) = BLU(MP) + VCLY - VCLYF
                BLU(IXP(NBS)) = BLU(IXP(NBS)) - VCLY + VCLYF
              ENDIF
              FGS = AFY(NPY)*VG(1,NPY)
              IF( FGS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/(DYGP(NPY)*XVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBS,NSL)*FCGS-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCGY = C(N,NSL)*FGS*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBS,NSL)*FGS*THETA*DYF*FCGS
                VCGYF = CO(N,NSL)*FGS*FCGP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCGY-VCGYF)/AFY(NPY)
                BLU(MP) = BLU(MP) + VCGY - VCGYF
                BLU(IXP(NBS)) = BLU(IXP(NBS)) - VCGY + VCGYF
              ENDIF
            ELSE
              ALN = MAX( -FLN,ZERO )
              AGN = MAX( -FGN,ZERO )
              AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP
              AN = ALN*FCL + AGN*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AN*BCX
          ENDIF
!
!---  Top boundary
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          NQZ = NSZ(N) + IJFLD
!
!---      Hydraulic dispersion  ---
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
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR.
     &      IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
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
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              IF( FLT.LT.ZERO ) THEN
                WCLZ = BCX*FCL*FLT
              ELSEIF( FLT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCLB = YL(NBB,NSL)/(SL(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBB,NSL)*FCLB)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                WCLZ =  C(N,NSL)*FLT*(1.D+0-THETA)*FCLP
     &            + BCX*FLT*THETA*FCL
              ELSEIF( FLT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCLZ =  C(N,NSL)*FLT*FCLP
              ENDIF
              IF( FGT.LT.ZERO ) THEN
                WCGZ = BCX*FCG*FGT
              ELSEIF( FGT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCGB = YG(NBB,NSL)/(SG(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBB,NSL)*FCGB)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                WCGZ =  C(N,NSL)*FGT*(1.D+0-THETA)*FCGP
     &            + BCX*FGT*THETA*FCG
              ELSEIF( FGT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCGZ =  C(N,NSL)*FGT*FCGP
              ENDIF
              AT = DLZ*FCL + DGZ*FCG
              AP = DLZ*FCLP + DGZ*FCGP
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCGZ)/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCGZ
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPZ = NSZ(N)
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/(DZGP(NPZ)*XVLX+SMALL)
                FCLB = YL(NBB,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBB,NSL)*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCLZ = C(N,NSL)*FLB*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBB,NSL)*FLB*THETA*DZF*FCLB
                WCLZF = CO(N,NSL)*FLB*FCLP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ-WCLZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCLZ - WCLZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCLZ + WCLZF
              ENDIF
              FGB = AFZ(NPZ)*WG(1,NPZ)
              IF( FGB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/(DZGP(NPZ)*XVGX+SMALL)
                FCGB = YG(NBB,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBB,NSL)*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCGZ = C(N,NSL)*FGB*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBB,NSL)*FGB*THETA*DZF*FCGB
                WCGZF = CO(N,NSL)*FGB*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCGZ-WCGZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCGZ - WCGZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCGZ + WCGZF
              ENDIF
            ELSE
              ALT = MAX( -FLT,ZERO ) +
     &          DLZ*MAX((ONE-(TENTH*ABS(FLT)/(DLZ+SMALL)))**5,ZERO)
              AGT = MAX( -FGT,ZERO ) +
     &          DGZ*MAX((ONE-(TENTH*ABS(FGT)/(DGZ+SMALL)))**5,ZERO)
              AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP
              AT = ALT*FCL + AGT*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AT*BCX
!
!---  Outflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.7 ) THEN
            FLT = MAX( FLT,0.D+0 )
            FGT = MAX( FGT,0.D+0 )
            FNT = MAX( FNT,0.D+0 )
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              WCLZ = 0.D+0
              IF( FLT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCLB = YL(NBB,NSL)/(SL(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBB,NSL)*FCLB)
     &            /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                WCLZ =  C(N,NSL)*FLT*(1.D+0-THETA)*FCLP
     &            + BCX*FLT*THETA*FCL
              ELSEIF( FLT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCLZ =  C(N,NSL)*FLT*FCLP
              ENDIF
              WCGZ = 0.D+0
              IF( FGT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCGB = YG(NBB,NSL)/(SG(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBB,NSL)*FCGB)
     &            /(BCX*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                WCGZ =  C(N,NSL)*FGT*(1.D+0-THETA)*FCGP
     &            + BCX*FGT*THETA*FCG
              ELSEIF( FGT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCGZ =  C(N,NSL)*FGT*FCGP
              ENDIF
              AT = 0.D+0
              AP = 0.D+0
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCGZ)/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCGZ
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPZ = NSZ(N)
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/(DZGP(NPZ)*XVLX+SMALL)
                FCLB = YL(NBB,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBB,NSL)*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCLZ = C(N,NSL)*FLB*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBB,NSL)*FLB*THETA*DZF*FCLB
                WCLZF = CO(N,NSL)*FLB*FCLP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ-WCLZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCLZ - WCLZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCLZ + WCLZF
              ENDIF
              FGB = AFZ(NPZ)*WG(1,NPZ)
              IF( FGB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/(DZGP(NPZ)*XVGX+SMALL)
                FCGB = YG(NBB,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBB,NSL)*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCGZ = C(N,NSL)*FGB*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBB,NSL)*FGB*THETA*DZF*FCGB
                WCGZF = CO(N,NSL)*FGB*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCGZ-WCGZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCGZ - WCGZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCGZ + WCGZF
              ENDIF
            ELSE
              ALT = MAX( -FLT,ZERO )
              AGT = MAX( -FGT,ZERO )
              AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP
              AT = ALT*FCL + AGT*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AT*BCX
!
!--- Inflow ---
!
          ELSEIF( IBCT(NSL+LUK,NB).EQ.13 .OR.
     &      IBCT(NSL+LUK,NB).EQ.14 ) THEN
            FLT = MIN( FLT,0.D+0 )
            FGT = MIN( FGT,0.D+0 )
            FNT = MIN( FNT,0.D+0 )
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              WCLZ = 0.D+0
              IF( FLT.LT.ZERO ) WCLZ = BCX*FCL*FLT
              WCGZ = 0.D+0
              IF( FGT.LT.ZERO ) WCGZ = BCX*FCG*FGT
              AT = 0.D+0
              AP = 0.D+0
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCGZ)/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCGZ
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPZ = NSZ(N)
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/(DZGP(NPZ)*XVLX+SMALL)
                FCLB = YL(NBB,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCL)
     &            /(C(NBB,NSL)*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCLZ = C(N,NSL)*FLB*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBB,NSL)*FLB*THETA*DZF*FCLB
                WCLZF = CO(N,NSL)*FLB*FCLP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ-WCLZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCLZ - WCLZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCLZ + WCLZF
              ENDIF
              FGB = AFZ(NPZ)*WG(1,NPZ)
              IF( FGB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/(DZGP(NPZ)*XVGX+SMALL)
                FCGB = YG(NBB,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCG)
     &            /(C(NBB,NSL)*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCGZ = C(N,NSL)*FGB*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBB,NSL)*FGB*THETA*DZF*FCGB
                WCGZF = CO(N,NSL)*FGB*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCGZ-WCGZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCGZ - WCGZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCGZ + WCGZF
              ENDIF
            ELSE
              ALT = MAX( -FLT,ZERO )
              AGT = MAX( -FGT,ZERO )
              AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP
              AT = ALT*FCL + AGT*FCG
            ENDIF
            BLU(MP) = BLU(MP) + AT*BCX
          ENDIF
        ENDIF
        IF( ILES.EQ.1 ) THEN
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + AP
        ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
          DLU(MCOL) = DLU(MCOL) + AP
#ifdef petsc
        ELSEIF( ILES.EQ.5 ) THEN
          CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,ADD_VALUES,IERR)
#endif
        ENDIF
  200 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SBND_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORC_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Compute source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE NAPL
      USE JACOB
      USE GRID
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVH
      USE FDVA
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 SRX(8)
      REAL*8 HIDX(LHF_HT),HRESX(LHF_HT)
      REAL*8 YMGX(LHF_HT),YMNX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC_HYD_O2_CH4_N2'
!
!---  Zero source terms  ---
!
      DO 60 NS = 1,NSR
        DO 50 I = ISRDM(1,NS),ISRDM(2,NS)
          DO 50 J = ISRDM(3,NS),ISRDM(4,NS)
            DO 50 K = ISRDM(5,NS),ISRDM(6,NS)
              N = ND(I,J,K)
              DO 40 M = 2,ISVC+2
                SRCA(M,N) = 0.D+0
                SRCO(M,N) = 0.D+0
                SRCN(M,N) = 0.D+0
                SRCW(M,N) = 0.D+0
                SRCS(M,N) = 0.D+0
                SRCT(M,N) = 0.D+0
   40         CONTINUE
   50   CONTINUE
   60 CONTINUE
!
!---  Loop over sources  ---
!
      DO 600 NS = 1,NSR
        IF( TM.LE.SRC(1,1,NS) ) GOTO 600
        SRX(1) = TM
        IF( ISRM(NS).EQ.1 ) THEN
          DO 70 N = 1,8
            SRX(N) = SRC(N,1,NS)
   70     CONTINUE
        ELSE
          DO 100 M = 2,ISRM(NS)
            IF( TM.LE.SRC(1,M,NS) ) THEN
             DTSR = MIN( SRC(1,M,NS)-TM,DT )
             TFSR = (TM-0.5D+0*DTSR-SRC(1,M-1,NS))/
     &         (SRC(1,M,NS)-SRC(1,M-1,NS))
             DO 80 N = 1,8
               SRX(N) = SRC(N,M-1,NS) + TFSR*(SRC(N,M,NS)-SRC(N,M-1,NS))
   80        CONTINUE
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 600
        ENDIF
  110   CONTINUE
!
!---  Loop source domain  ---
!
        DO 500 I = ISRDM(1,NS),ISRDM(2,NS)
          DO 500 J = ISRDM(3,NS),ISRDM(4,NS)
            DO 500 K = ISRDM(5,NS),ISRDM(6,NS)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 500
              DO 400 M = 2,ISVC+2
                PGX = PG(M,N) + PATM
                PNX = PN(M,N) + PATM
                PLX = PL(M,N) + PATM
                TX = T(M,N)
!
!---            Power  ---
!
                IF( ISRT(NS).EQ.1 ) THEN
                  SRCT(M,N) = SRCT(M,N) + SRX(4)
!
!---            Power Density  ---
!
                ELSEIF( ISRT(NS).EQ.2 ) THEN
                  SRCT(M,N) = SRCT(M,N) + SRX(4)*VOL(N)
!
!---            Fluid mass rate  ---
!
                ELSEIF( ISRT(NS).EQ.3 ) THEN
                  XFAGX = XGA(M,N)*RHOG(M,N)*SG(M,N)
                  XFALX = XLA(M,N)*RHOL(M,N)*SL(M,N)
                  XFANX = XNA(M,N)*RHON(M,N)*SN(M,N)
                  XFAX = XFAGX + XFALX + XFANX
                  XFOGX = XGO(M,N)*RHOG(M,N)*SG(M,N)
                  XFOLX = XLO(M,N)*RHOL(M,N)*SL(M,N)
                  XFONX = XNO(M,N)*RHON(M,N)*SN(M,N)
                  XFOX = XFOGX + XFOLX + XFONX
                  XFNGX = XGN(M,N)*RHOG(M,N)*SG(M,N)
                  XFNLX = XLN(M,N)*RHOL(M,N)*SL(M,N)
                  XFNNX = XNN(M,N)*RHON(M,N)*SN(M,N)
                  XFNX = XFNGX + XFNLX + XFNNX
                  XFSLX = XLS(M,N)*RHOL(M,N)*SL(M,N)
                  XFSX = XFSLX
                  XFWGX = XGW(M,N)*RHOG(M,N)*SG(M,N)
                  XFWLX = XLW(M,N)*RHOL(M,N)*SL(M,N)
                  XFWNX = XNW(M,N)*RHON(M,N)*SN(M,N)
                  XFWX = XFWGX + XFWLX + XFWNX
                  XFTX = XFAX+XFOX+XFSX+XFWX+XFNX
                  SRCA(M,N) = SRCA(M,N) + SRX(4)*XFAX/XFTX
                  SRCO(M,N) = SRCO(M,N) + SRX(4)*XFOX/XFTX
                  SRCN(M,N) = SRCN(M,N) + SRX(4)*XFNX/XFTX
                  SRCW(M,N) = SRCW(M,N) + SRX(4)*XFWX/XFTX
                  SRCS(M,N) = SRCS(M,N) + SRX(4)*XFSX/XFTX
                  XFGX = (XFAGX+XFOGX+XFWGX+XFNGX)/XFTX
                  XFLX = (XFALX+XFOLX+XFSLX+XFWLX+XFNLX)/XFTX
                  XFNX = (XFANX+XFONX+XFNNX)/XFTX
                  SRCT(M,N) = SRCT(M,N) + SRX(4)*
     &              (XFGX*HG(M,N) + XFLX*HL(M,N) + XFNX*HN(M,N))
!
!---            Nonaqueous fluid mass rate  ---
!
                ELSEIF( ISRT(NS).EQ.4 ) THEN
                  PX = PSO(M,N) + PATM
                  PSOX = PSO(M,N)
                  TX = SRX(2)
                  WGWX = SRX(3)
                  CALL SP_W( TX,PSWX )
                  PVWX = WGWX*PSWX          
!
!---              Total nonaqueous component mole fractions  ---
!
                  ZMCAX = SRX(5)
                  ZMCOX = SRX(6)
                  ZMCNX = MIN( MAX( 1.D+0-ZMCAX-ZMCOX,0.D+0 ),1.D+0 )
!
!---              Nonaqueous fluid state and compositions  ---
!
                  CALL FLH_SR1( PPELX,PPEUX,PSOX,PVAX,PVNX,PVOX,PVWX,
     &              RHOGX,RHONX,SGX,SNX,TX,TCRX,TCTX,XGAX,XGNX,
     &              XGOX,XGWX,XNAX,XNNX,XNOX,XNWX,YMGAX,YMGNX,YMGOX,
     &              YMNAX,YMNNX,YMNOX,ZMCAX,ZMCOX,ZMCNX )
!
!---              Pure water vapor enthalpy  ---
!
                  CALL SP_W( TX,PSWX )
                  PEX = MIN( PSWX,PSOX )
                  ISRX = 2
                  CALL ENTH_W( TX,PEX,HLWX,HGWX,ISRX )
!
!---              Gas enthalpy and internal energy  ---
!
                  YMGX(1) = YMGAX
                  YMGX(2) = YMGOX
                  YMGX(3) = YMGNX
                  CALL ENTH_ID( TX,HIDX )
                  INDX = 1
                  CALL ENTH_RES( PSOX,TX,HRESX,YMGX,INDX )
                  HGAX = (HIDX(1) + HRESX(1))/WTMA
                  HGOX = (HIDX(2) + HRESX(2))/WTMO
                  HGNX = (HIDX(3) + HRESX(3))/WTMN
                  HGX = XGAX*HGAX + XGOX*HGOX + XGNX*HGNX + XGWX*HGWX
!
!---              Nonaqueous-liquid enthalpy and internal energy  ---
!
                  YMNX(1) = YMNAX
                  YMNX(2) = YMNOX
                  YMNX(3) = YMNNX
                  CALL ENTH_ID( TX,HIDX )
                  INDX = 2
                  CALL ENTH_RES( PSOX,TX,HRESX,YMNX,INDX )
                  HNAX = (HIDX(1) + HRESX(1))/WTMA
                  HNOX = (HIDX(2) + HRESX(2))/WTMO
                  HNNX = (HIDX(3) + HRESX(3))/WTMN
                  HNX = XNAX*HNAX + XNOX*HNOX + XNNX*HNNX + XNWX*HGWX
!
!---              Total density  ---
!
                  RHOX = SGX*RHOGX + SNX*RHONX
!
!---              O2 mass fractions  ---
!
                  XFGAX = XGAX*SGX*RHOGX/(RHOX+SMALL)
                  XFNAX = XNAX*SNX*RHONX/(RHOX+SMALL)
                  XFAX = XFGAX + XFNAX
!
!---              CH4 mass fractions  ---
!
                  XFGOX = XGOX*SGX*RHOGX/(RHOX+SMALL)
                  XFNOX = XNOX*SNX*RHONX/(RHOX+SMALL)
                  XFOX = XFGOX + XFNOX
!
!---              N2 mass fractions  ---
!
                  XFGNX = XGNX*SGX*RHOGX/(RHOX+SMALL)
                  XFNNX = XNNX*SNX*RHONX/(RHOX+SMALL)
                  XFNX = XFGNX + XFNNX
!
!---              Water mass fractions  ---
!
                  XFGWX = XGWX*SGX*RHOGX/(RHOX+SMALL)
                  XFNWX = XNWX*SNX*RHONX/(RHOX+SMALL)
                  XFWX = XFGWX + XFNWX
!
!---              Energy fractions  ---
!
                  XFGEX = HGX*SGX*RHOGX/(RHOX+SMALL)
                  XFNEX = HNX*SNX*RHONX/(RHOX+SMALL)
                  XFEX = XFGEX + XFNEX
!
!---              Mass and energy sources  ---
!
                  SRCA(M,N) = SRCA(M,N) + SRX(4)*XFAX
                  SRCO(M,N) = SRCO(M,N) + SRX(4)*XFOX
                  SRCN(M,N) = SRCN(M,N) + SRX(4)*XFNX
                  SRCW(M,N) = SRCW(M,N) + SRX(4)*XFWX
                  SRCT(M,N) = SRCT(M,N) + SRX(4)*XFEX
!
!---            Initial-condition fluid mass rate  ---
!
                ELSEIF( ISRT(NS).EQ.5 ) THEN
                  IF( (NSTEP-NRST).EQ.1 .AND. NITER.EQ.1
     &              .AND. M.EQ.2 ) THEN
                    DO 200 MX = 1,ISRM(NS)
                      SRC(3,MX,NS) = HL(2,N)
                      SRC(5,MX,NS) = XLA(2,N)
                      SRC(6,MX,NS) = XLO(2,N)
                      SRC(7,MX,NS) = XLS(2,N)
                      SRC(8,MX,NS) = XLW(2,N)
                      SRX(3) = HL(2,N)
                      SRX(5) = XLA(2,N)
                      SRX(6) = XLO(2,N)
                      SRX(7) = XLS(2,N)
                      SRX(8) = XLW(2,N)
  200               CONTINUE
                  ENDIF
                  SRCT(M,N) = SRCT(M,N) + SRX(4)*SRX(3)
                  SRCA(M,N) = SRCA(M,N) + SRX(4)*SRX(5)
                  SRCO(M,N) = SRCO(M,N) + SRX(4)*SRX(6)
                  SRCS(M,N) = SRCS(M,N) + SRX(4)*SRX(7)
                  SRCW(M,N) = SRCW(M,N) + SRX(4)*SRX(8)
                ENDIF
  400         CONTINUE
  500   CONTINUE
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORC_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORIT_HYD_O2_CH4_N2( NSL )
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Compute solute transport source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE GRID
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 SRX(8)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORIT_HYD_O2_CH4_N2'
!
!---  Loop over sources  ---
!
      DO 600 NS = 1,NSR
        IF( TM.LE.SRC(1,1,NS) ) GOTO 600
        SRX(1) = TM
        IF( ISRM(NS).EQ.1 ) THEN
          SRX(4) = SRC(4,1,NS)
        ELSE
          DO 100 M = 2,ISRM(NS)
            IF( TM.LE.SRC(1,M,NS) ) THEN
             DTSR = MIN( SRC(1,M,NS)-TM,DT )
             TFSR = (TM-0.5D+0*DTSR-SRC(1,M-1,NS))/
     &         (SRC(1,M,NS)-SRC(1,M-1,NS))
             SRX(4) = SRC(4,M-1,NS) + TFSR*(SRC(4,M,NS)-SRC(4,M-1,NS))
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 600
        ENDIF
  110   CONTINUE
!
!---  Loop over source domain  ---
!
        DO 500 I = ISRDM(1,NS),ISRDM(2,NS)
          DO 500 J = ISRDM(3,NS),ISRDM(4,NS)
            DO 500 K = ISRDM(5,NS),ISRDM(6,NS)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 500
!
!---  Aqueous Volumetric Sink  ---
!
              IF( ISRT(NS).EQ.3 .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YL(N,NSL)*DT/(PORD(2,N)*SL(2,N))
!
!---  Gas Volumetric Sink  ---
!
              ELSEIF( (ISRT(NS).EQ.4 .OR. ISRT(NS).EQ.5)
     &          .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YG(N,NSL)*DT/(PORD(2,N)*SG(2,N))
!
!---  Aqueous Mass Sink  ---
!
              ELSEIF( ISRT(NS).EQ.7 .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YL(N,NSL)*DT/(RHOL(2,N)*PORD(2,N)*SL(2,N))
!
!---  Gas Mass Sink  ---
!
              ELSEIF( (ISRT(NS).EQ.8 .OR. ISRT(NS).EQ.9)
     &          .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YG(N,NSL)*DT/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---  Solute source  ---
!
              ELSEIF( ISRT(NS).LT.0 .AND. ISRT(NS).GE.-NSOLU ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(4)*DT
!
!---  Solute source  ---
!
              ELSEIF( ISRT(NS).LT.-NSOLU .AND.
     &          ISRT(NS).GE.-2*NSOLU ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(4)*DT*VOL(N)
              ENDIF
  500   CONTINUE
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORIT_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT_HYD_O2_CH4_N2( NSL )
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Compute solute transport source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE JACOB
      USE GRID
      USE FDVP
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
#ifdef petsc
      INTEGER :: IROW
      REAL*8 :: BUFFER
      PetscErrorCode :: IERR
#endif
      REAL*8 SRX(8)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORT_HYD_O2_CH4_N2'
!
!---  Loop over sources  ---
!
      DO 600 NS = 1,NSR
        IF( TM.LE.SRC(1,1,NS) ) GOTO 600
        SRX(1) = TM
        IF( ISRM(NS).EQ.1 ) THEN
          SRX(4) = SRC(4,1,NS)
        ELSE
          DO 100 M = 2,ISRM(NS)
            IF( TM.LE.SRC(1,M,NS) ) THEN
             DTSR = MIN( SRC(1,M,NS)-TM,DT )
             TFSR = (TM-0.5D+0*DTSR-SRC(1,M-1,NS))/
     &         (SRC(1,M,NS)-SRC(1,M-1,NS))
             SRX(4) = SRC(4,M-1,NS) + TFSR*(SRC(4,M,NS)-SRC(4,M-1,NS))
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 600
        ENDIF
  110   CONTINUE
!
!---    Loop over source domain  ---
!
        DO 500 I = ISRDM(1,NS),ISRDM(2,NS)
        DO 500 J = ISRDM(3,NS),ISRDM(4,NS)
        DO 500 K = ISRDM(5,NS),ISRDM(6,NS)
          N = ND(I,J,K)
          IF( IXP(N).EQ.0 ) GOTO 500
          MP = IXP(N)
          IF( ILES.EQ.1 ) THEN
            MCOL = MP
            MROW = MDT
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            MA = 1
            MCOL = KLUC(MP,MA)
            MA = MA + 1
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            IROW = petsc_offset(N)
#endif
          ENDIF
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
          ENDIF
!
!---      Load Jacobian  ---
!
          IF( ILES.EQ.1 ) THEN
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + SORTX
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCOL) = DLU(MCOL) + SORTX
#ifdef petsc
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,SORTX,
     &                        ADD_VALUES,IERR)
#endif
          ENDIF
  500   CONTINUE
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORT_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_HYD_O2_CH4_N2( PIX,PLX,PNX,SHX,SIX,SLX,IZN,N )
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
!----------------------Data Statements---------------------------------!
!
      DATA TOL /1.D-12/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SP_HYD_O2_CH4_N2'
      N_DB = N
!
!---  Residual aqueous saturation for relative permeability  ---
!
      SLRX = SCHR(4,IZN)
!
!---  van Genuchten saturation function w/o or w/ the
!     hydrate capillary pressure barrier or w/ the
!     hydrate-porosity model  ---
!
      IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 .OR.
     &  ISCHR(IZN).EQ.13 ) THEN
        IF( SHX.LT.TOL ) SHX = 0.D+0
        HCPBFX = 1.D+0
!
!---    Entry pressure with hydrate barrier adjustment  ---
!
        IF( ISCHR(IZN).EQ.11 ) THEN
          HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &      SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
!
!---    Entry pressure with hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.13 ) THEN
          HCPBFX = (1.D+0-SHX)**(-(SCHR(10,IZN)-1.D+0)/SCHR(10,IZN))
        ENDIF
        HNL = MAX( 0.D+0,BNL*(PNX-PLX)/RHORL/GRAV )/HCPBFX
        HIL = MAX( 0.D+0,BIL*(PIX-PLX)/RHORL/GRAV )/HCPBFX
        CNX = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( IRPL(IZN).EQ.2 ) THEN
            CMX = 1.D+0 - 2.D+0/CNX
          ELSE
            CMX = 1.D+0 - 1.D+0/CNX
         ENDIF
        ELSE
          CMX = SCHR(14,IZN)
        ENDIF
        SRX = SCHR(4,IZN)
!
!---    Webb extension  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HMPX = SCHR(17,IZN)
!
!---      Gas-aqueous capillary head above matching-point
!         head, use Webb extension  ---
!
          IF( HNL.GT.HMPX ) THEN
            SMPX = SCHR(16,IZN)
            HNL = MIN( HNL,HDOD )
            DMPX = SMPX/(LOG10(HDOD)-LOG10(HMPX))
            ASLX = -(LOG10(HNL)-LOG10(HDOD))*DMPX
            ESLX = ASLX
!
!---      Gas-aqueous capillary head below 
!         matching-point head  ---
!
          ELSE
            ASLX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HNL)**CNX))**CMX
            ESLX = (ASLX*(1.D+0-SRX) + SRX)
          ENDIF
!
!---      Ice-aqueous capillary head above matching-point
!         head, use Webb extension  ---
!
          IF( HIL.GT.HMPX ) THEN
            SMPX = SCHR(16,IZN)
            HIL = MIN( HIL,HDOD )
            DMPX = SMPX/(LOG10(HDOD)-LOG10(HMPX))
            ASLIX = -(LOG10(HIL)-LOG10(HDOD))*DMPX
            ESLIX = ASLIX
            ESLIX = MIN( ESLX,ESLIX )
!
!---      Ice-aqueous capillary head below 
!         matching-point head  ---
!
          ELSE
            ASLIX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HIL)**CNX))**CMX
            ESLIX = (ASLIX*(1.D+0-SRX) + SRX)
            ESLIX = MIN( ESLX,ESLIX )
          ENDIF
!
!---    No extension  ---
!
        ELSE
          ASLX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HNL)**CNX))**CMX
          ASLIX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HIL)**CNX))**CMX
          ESLX = (ASLX*(1.D+0-SRX) + SRX)
          ESLIX = (ASLIX*(1.D+0-SRX) + SRX)
          ESLIX = MIN( ESLX,ESLIX )
        ENDIF
!
!---    Ice saturation  ---
!
        ESIX = ESLX-ESLIX
        IF( ESIX.LT.TOL ) ESIX = 0.D+0
        ESLX = ESLX-ESIX
        IF( ESLX.LT.TOL ) ESLX = 0.D+0
        ESGNX = 1.D+0-ESLX-ESIX
        IF( ESGNX.LT.TOL ) ESGNX = 0.D+0
!
!---    Actual saturations  ---
!
!        SLX = ESLX*(1.D+0-SHX)
        SIX = ESIX*(1.D+0-SHX)
        SLX = MAX( SLX-SIX,0.D+0 )
!
!---    Effective saturations for relative permeability
!       calculations  ---
!
        ASLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
!
!---  Brooks and Corey saturation function w/o or w/ the
!     hydrate capillary pressure barrier or w/ the
!     hydrate-porosity model  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 .OR.
     &  ISCHR(IZN).EQ.14 ) THEN
        IF( SHX.LT.TOL ) SHX = 0.D+0
        HCPBFX = 1.D+0
!
!---    Entry pressure with hydrate barrier adjustment  ---
!
        IF( ISCHR(IZN).EQ.12 ) THEN
          HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &      SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
!
!---    Entry pressure with hydrate-porosity model  ---
!
        ELSEIF( ISCHR(IZN).EQ.14 ) THEN
          HCPBFX = (1.D+0-SHX)**(-(SCHR(10,IZN)-1.D+0)/SCHR(10,IZN))
        ENDIF
        HNL = MAX( 0.D+0,BNL*(PNX-PLX)/RHORL/GRAV )/HCPBFX
        HIL = MAX( ZERO,BIL*(PIX-PLX)/RHORL/GRAV )/HCPBFX
        CLX = MAX( SCHR(3,IZN),SMALL )
        SRX = SCHR(4,IZN)
!
!---    Webb extension  ---
!
        IF( ISM(IZN).EQ.2 ) THEN
          HMPX = SCHR(17,IZN)
!
!---      Gas-aqueous capillary head above matching-point
!         head, use Webb extension  ---
!
          IF( HNL.GT.HMPX ) THEN
            SMPX = SCHR(16,IZN)
            HNL = MIN( HNL,HDOD )
            DMPX = SMPX/(LOG10(HDOD)-LOG10(HMPX))
            ASLX = -(LOG10(HNL)-LOG10(HDOD))*DMPX
            ESLX = ASLX
!
!---      Gas-aqueous capillary head below 
!          matching-point head  ---
!
          ELSE
            IF( HNL-SCHR(2,IZN).LE.1.D-12 ) THEN
              ASLX = 1.D+0
            ELSE
              ASLX = (SCHR(2,IZN)/HNL)**CLX
            ENDIF
            ESLX = (ASLX*(1.D+0-SRX) + SRX)
          ENDIF
!
!---      Ice-aqueous capillary head above matching-point
!         head, use Webb extension  ---
!
          IF( HIL.GT.HMPX ) THEN
            SMPX = SCHR(16,IZN)
            HIL = MIN( HIL,HDOD )
            DMPX = SMPX/(LOG10(HDOD)-LOG10(HMPX))
            ASLIX = -(LOG10(HIL)-LOG10(HDOD))*DMPX
            ESLIX = ASLIX
            ESLIX = MIN( ESLX,ESLIX )
!
!---      Ice-aqueous capillary head below 
!         matching-point head  ---
!
          ELSE
            IF( HIL-SCHR(2,IZN).LE.1.D-12 ) THEN
              ASLIX = 1.D+0
            ELSE
              ASLIX = (SCHR(2,IZN)/HNL)**CLX
            ENDIF
            ESLIX = (ASLIX*(1.D+0-SRX) + SRX)
            ESLIX = MIN( ESLX,ESLIX )
          ENDIF
!
!---    No extension  ---
!
        ELSE
          IF( HNL-SCHR(2,IZN).LE.1.D-12 ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (SCHR(2,IZN)/HNL)**CLX
          ENDIF
          IF( HIL-SCHR(2,IZN).LE.1.D-12 ) THEN
            ASLIX = 1.D+0
          ELSE
            ASLIX = (SCHR(2,IZN)/HNL)**CLX
          ENDIF
          ESLX = (ASLX*(1.D+0-SRX) + SRX)
          ESLIX = (ASLIX*(1.D+0-SRX) + SRX)
          ESLIX = MIN( ESLX,ESLIX )
        ENDIF
!
!---    Ice saturation  ---
!
        ESIX = ESLX-ESLIX
        IF( ESIX.LT.TOL ) ESIX = 0.D+0
        ESLX = ESLX-ESIX
        IF( ESLX.LT.TOL ) ESLX = 0.D+0
        ESGNX = 1.D+0-ESLX-ESIX
        IF( ESGNX.LT.TOL ) ESGNX = 0.D+0
!
!---    Actual saturations  ---
!
!        SLX = ESLX*(1.D+0-SHX)
        SIX = ESIX*(1.D+0-SHX)
        SLX = MAX( SLX-SIX,0.D+0 )
!
!---    Effective saturations for relative permeability
!       calculations  ---
!
        ASLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP_HYD_O2_CH4_N2( NSL )
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Calculates the aqueous- and gas-phase solute
!     mole fractions from user-specified partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
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
      SUB_LOG(ISUB_LOG) = '/SPRP_HYD_O2_CH4_N2'
!
!---  Loop over all nodes  ---
!
      DO 900 N = 1,NFLD
        N_DB = N
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
     &      + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK + PCGL(5,NSL)*TK**2 )
!
!---    Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.2 ) THEN
          PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
        ENDIF
        PCGLX = MAX( PCGLX,1.D-20 )
        PCGLX = MIN( PCGLX,1.D+20 )
!
!---  Phase-volumetric concentration ratios  ---
!
        YVL = 1.D+0/(XVS + XVL + XVG*PCGLX)
        YVG = 1.D+0/((XVS + XVL)/PCGLX + XVG)
!
!---  Phase mole fractions  ---
!
        YL(N,NSL) = XVL*YVL
        YG(N,NSL) = XVG*YVG
!
  900 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPRP_HYD_O2_CH4_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT_HYD_O2_CH4_N2
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
!     HYDT-KE (Ternary Gas Hydrate w/ Kinetic Exchange) Mode
!
!     Update the primary variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 July 2019
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE OUTPU
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
      USE FDVS
      USE FDVP
      USE FDVH
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*128 PH_CND(43)
!
!----------------------Data Statements---------------------------------!
!
      DATA PH_CND /9*'','SH = 0, SG = 0, SN = 0',10*'',
     &  'SH = 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN )',
     &  'SH = 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN )',
     &  'SH = 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO )',7*'',
     &  'SH > 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN ), YMHGA >= MAX(YMH
     &GO,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN ), YMHGA >= MAX(YMH
     &GO,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO ), YMHGA >= MAX(YMH
     &GO,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN ), YMHGO >= MAX(YMH
     &GA,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN ), YMHGO >= MAX(YMH
     &GA,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO ), YMHGO >= MAX(YMH
     &GA,YMHGN )',
     &  'SH > 0, SG + SN > 0, ZMCA >= MAX( ZMCO,ZMCN ), YMHGN >= MAX(YMH
     &GA,YMHGO )',
     &  'SH > 0, SG + SN > 0, ZMCO >= MAX( ZMCA,ZMCN ), YMHGN >= MAX(YMH
     &GA,YMHGO )',
     &  'SH > 0, SG + SN > 0, ZMCN >= MAX( ZMCA,ZMCO ), YMHGN >= MAX(YMH
     &GA,YMHGO )','',
     &  'SH > 0, SG + SN = 0, YMHGA >= MAX(YMHGO,YMHGN )',
     &  'SH > 0, SG + SN = 0, YMHGO >= MAX(YMHGA,YMHGN )',
     &  'SH > 0, SG + SN = 0, YMHGN >= MAX(YMHGA,YMHGO )'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDT_HYD_O2_CH4_N2'
      IF( ICNV.EQ.1 ) GOTO 300
      IERR = 0
!
!---  Update primary variables
!
      DO 200 N = 1,NFLD
        N_DB = N
        IF( IXP(N).EQ.0 ) GOTO 200
        NMD = IXP(N)
        IZN = IZ(N)
        MPT = IM(IEQT,NMD)
        DPT = BLU(MPT)
        MPW = IM(IEQW,NMD)
        DPW = BLU(MPW)
        IF( ISLC(45).EQ.0 ) THEN
          MPA = IM(IEQA,NMD)
          DPA = BLU(MPA)
        ENDIF
        IF( ISLC(64).EQ.0 ) THEN
          MPO = IM(IEQO,NMD)
          DPO = BLU(MPO)
        ENDIF
        IF( ISLC(65).EQ.0 ) THEN
          MPN = IM(IEQN,NMD)
          DPN = BLU(MPN)
        ELSE
          DPN = 0.D+0
        ENDIF
        IF( ISLC(45).EQ.0 ) THEN
          MPHA = IM(IEQHA,NMD)
          DPHA = BLU(MPHA)
        ENDIF
        IF( ISLC(64).EQ.0 ) THEN
          MPHO = IM(IEQHO,NMD)
          DPHO = BLU(MPHO)
        ENDIF
        IF( ISLC(65).EQ.0 ) THEN
          MPHN = IM(IEQHN,NMD)
          DPHN = BLU(MPHN)
        ELSE
          DPHN = 0.D+0
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
!---    Excess changes in primary variables, skip updates  ---
!
        IF( ABS(DPT).GT.1.D+2 ) IERR = N
        IF( IERR.NE.0 ) CYCLE
!
!---    Limit salt mass fraction changes to 0.25 of the
!       maximum value if salt mass fraction is less than
!       the maximum, isobrine option   ---
!
        IF( ISLC(32).EQ.0 ) THEN
          CALL SOL_LS( T(2,N),XLSMX )
          IF( YLS(2,N).LT.XLSMX ) THEN
            DPS = SIGN( MIN( ABS(2.5D-1*XLSMX),ABS(DPS) ),DPS )
          ENDIF
!
!---      Zero negative corrections for zero dissolved salt  ---
!
          IF( YLS(2,N)/EPSL.LT.EPSL .AND. DPS/EPSL.LT.EPSL ) THEN
            BLU(MPS) = 0.D+0
            DPS = 0.D+0
          ENDIF
          YLS(2,N) = YLS(2,N)+DPS
          IF( YLS(2,N).LT.1.D-12 ) YLS(2,N) = 0.D+0
          XLS(2,N) = MIN( YLS(2,N),XLSMX )
        ENDIF
!
!---    Freezing point  ---
!
        CALL ICE( T(2,N),PL(2,N),XLS(2,N),TFPX,PIX )
!
!---    Limit changes in temperature to 0.25 C without hydrate, 
!       and to 0.1 C with hydrate   ---
!
        IF( SH(2,N).GT.EPSL ) THEN
          DPT = SIGN( MIN( 1.D-1,ABS(DPT) ),DPT )
        ELSE
          DPT = SIGN( MIN( 2.5D-1,ABS(DPT) ),DPT )
        ENDIF
!
!---    Relax temperature change in transition from no-ice to ice,
!       relax temperature change in transition from ice to no-ice  ---
!
        DTFPX = T(2,N)-TFPX
        IF( DTFPX.GE.0.D+0 .AND. (DTFPX+DPT).LT.0.D+0 ) THEN
          T(2,N) = TFPX + 6.D-1*DPT
        ELSEIF( DTFPX.LT.0.D+0 .AND. (DTFPX+DPT).GE.0.D+0 ) THEN
          T(2,N) = TFPX + 6.D-1*DPT
        ELSE
          T(2,N) = T(2,N) + DPT
        ENDIF
!
!---    Core phase condition  ---
!
        NPHZX = MOD(NPHAZ(2,N),100)
!
!---    Phase condition #10 (10 series)  ---
!
!       SH = 0, SG + SN = 0, SL + SI = 1
!
!       energy - temperature, C - T
!       water mass - pressure, Pa (gauge) - PSO
!       mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!       mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!       mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!       hydrate O2 mass - hydrate O2 mass, kg - TMHA
!       hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!       hydrate N2 mass - hydrate N2 mass, kg - TMHN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        IF( NPHZX.EQ.10 ) THEN
!
!---      Limit changes in pressure  ---
!
          DPX = 2.D-1*(PSO(2,N)+PATM)
          DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
          PSO(2,N) = PSO(2,N) + DPW
!
!---      Zero negative corrections for zero O2 partial pressure  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( PVA(2,N)/EPSL.LT.EPSL .AND. DPA/EPSL.LT.EPSL ) THEN
              BLU(MPA) = 0.D+0
              DPA = 0.D+0
            ENDIF
            DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
            PVA(2,N) = MAX( (PVA(2,N)+DPA),0.D+0 )
            IF( PVA(2,N).LT.1.D-6 ) PVA(2,N) = 0.D+0
          ENDIF
!
!---      Zero negative corrections for zero CH4 partial pressure  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            IF( PVO(2,N)/EPSL.LT.EPSL .AND. DPO/EPSL.LT.EPSL ) THEN
              BLU(MPO) = 0.D+0
              DPO = 0.D+0
            ENDIF
            DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
            PVO(2,N) = MAX( (PVO(2,N)+DPO),0.D+0 )
            IF( PVO(2,N).LT.1.D-6 ) PVO(2,N) = 0.D+0
          ENDIF
!
!---      Zero negative corrections for zero N2 partial pressure  ---
!
          IF( ISLC(65).EQ.0 ) THEN
            IF( PVN(2,N)/EPSL.LT.EPSL .AND. DPN/EPSL.LT.EPSL ) THEN
              BLU(MPN) = 0.D+0
              DPN = 0.D+0
            ENDIF
            DPN = SIGN( MIN( ABS(DPX),ABS(DPN) ),DPN )
            PVN(2,N) = MAX( (PVN(2,N)+DPN),0.D+0 )
            IF( PVN(2,N).LT.1.D-6 ) PVN(2,N) = 0.D+0
          ENDIF
!
!---      Limit changes in O2 hydrate mass  ---
!
          TMHA(2,N) = 0.D+0
!
!---      Limit changes in CH4 hydrate mass  ---
!
          TMHO(2,N) = 0.D+0
!
!---      Limit changes in N2 hydrate mass  ---
!
          TMHN(2,N) = 0.D+0
!
!---    Phase conditions 21 through 23 (20 series)  ---
!
!       SH = 0, SG + SN > 0, SL + SI < 1
!       ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #21
!       ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #22
!       ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #23
!
!       energy - temperature, C - T
!       water mass - aqueous + ice saturation - POSM
!       maximum mobile mass component - pressure, Pa (gauge) - PSO
!       mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!       mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!       mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!       hydrate O2 mass - hydrate O2 mass, kg - TMHA
!       hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!       hydrate N2 mass - hydrate N2 mass, kg - TMHN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.21 .AND. NPHZX.LE.23  ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          POSM(2,N) = POSM(2,N) + DPW
          IF( POSM(2,N).LT.EPSL ) POSM(2,N) = 0.D+0
!
!---      ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #21
!
          IF( NPHZX.EQ.21 ) THEN
!
!---        Gas conditions  ---
!
            IF( SG(2,N).GT.EPSL ) THEN
!
!---          Limit changes in pressure  ---
!
              DPX = MAX( 1.D+6,1.D-1*(PG(2,N)-PL(2,N)) )
              DPA = SIGN( MIN(ABS(DPX),ABS(DPA)),DPA )
              PSO(2,N) = MAX( (PSO(2,N)+DPA),0.D+0 )
              IF( PSO(2,N).LT.1.D-6 ) PSO(2,N) = 0.D+0
            ELSE
!
!---          Limit changes in pressure  ---
!
              DPX = 2.D-1*(PSO(2,N)+PATM)
              DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
              PSO(2,N) = PSO(2,N) + DPA
            ENDIF
          ELSE
!
!---        Zero negative corrections total nonaqueous 
!           O2 mole fraction  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              DPX = 2.5D-1
              DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
              ZMCA(2,N) = MAX( (ZMCA(2,N)+DPA),0.D+0 )
              IF( ZMCA(2,N).LT.EPSL ) ZMCA(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---      ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #22
!
          IF( NPHZX.EQ.22 ) THEN
!
!---        Gas conditions  ---
!
            IF( SG(2,N).GT.EPSL ) THEN
!
!---          Limit changes in pressure  ---
!
              DPX = MAX( 1.D+6,1.D-1*(PG(2,N)-PL(2,N)) )
              DPO = SIGN( MIN(ABS(DPX),ABS(DPO)),DPO )
              PSO(2,N) = MAX( (PSO(2,N)+DPO),0.D+0 )
              IF( PSO(2,N).LT.1.D-6 ) PSO(2,N) = 0.D+0
            ELSE
!
!---          Limit changes in pressure  ---
!
              DPX = 2.D-1*(PSO(2,N)+PATM)
              DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
              PSO(2,N) = PSO(2,N) + DPO
            ENDIF
          ELSE
!
!---        Zero negative corrections total nonaqueous 
!           CH4 mole fraction  ---
!
            IF( ISLC(64).EQ.0 ) THEN
              DPX = 2.5D-1
              DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
              ZMCO(2,N) = MAX( (ZMCO(2,N)+DPO),0.D+0 )
              IF( ZMCO(2,N).LT.EPSL ) ZMCO(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---      ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #23
!
          IF( NPHZX.EQ.23 ) THEN
!
!---        Gas conditions  ---
!
            IF( SG(2,N).GT.EPSL ) THEN
!
!---          Limit changes in pressure  ---
!
              DPX = MAX( 1.D+6,1.D-1*(PG(2,N)-PL(2,N)) )
              DPA = SIGN( MIN(ABS(DPX),ABS(DPN)),DPN )
              PSO(2,N) = MAX( (PSO(2,N)+DPN),0.D+0 )
              IF( PSO(2,N).LT.1.D-6 ) PSO(2,N) = 0.D+0
            ELSE
!
!---          Limit changes in pressure  ---
!
              DPX = 2.D-1*(PSO(2,N)+PATM)
              DPN = SIGN( MIN( ABS(DPX),ABS(DPN) ),DPN )
              PSO(2,N) = PSO(2,N) + DPN
            ENDIF
          ELSE
!
!---        Zero negative corrections total nonaqueous 
!           N2 mole fraction  ---
!
            IF( ISLC(65).EQ.0 ) THEN
              DPX = 2.5D-1
              DPN = SIGN( MIN( ABS(DPX),ABS(DPN) ),DPN )
              ZMCN(2,N) = MAX( (ZMCN(2,N)+DPN),0.D+0 )
              IF( ZMCN(2,N).LT.EPSL ) ZMCN(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---      Limit changes in O2 hydrate mass  ---
!
          TMHA(2,N) = 0.D+0
!
!---      Limit changes in CH4 hydrate mass  ---
!
          TMHO(2,N) = 0.D+0
!
!---      Limit changes in N2 hydrate mass  ---
!
          TMHN(2,N) = 0.D+0
!
!---    Phase conditions #31 through #39 (30 series)  ---
!
!       SH > 0, SG + SN > 0, SL + SI < 1
!       ZMCA = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #31
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #32
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #33
!       ZMCO = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #34
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #35
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #36
!       ZMCN = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #37
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #38
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #39
!
!       energy - temperature, C - T
!       water mass - aqueous + ice saturation, POSM
!       maximum mobile mass component - pressure, Pa (gauge) - PSO
!       mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!       mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!       mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!       maximum hydrate mass component - hydrate saturation - SH
!       hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!       hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!       hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.31 .AND. NPHZX.LE.39 ) THEN
!
!---      Limit changes in aqueous saturation  ---
!
          DPX = 2.5D-2
          DPW = SIGN( MIN(ABS(DPX),ABS(DPW)),DPW )
          POSM(2,N) = POSM(2,N) + DPW
          IF( POSM(2,N).LT.EPSL ) POSM(2,N) = 0.D+0
!
!---      ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #31, #32, or #33
!
          IF( NPHZX.GE.31 .AND. NPHZX.LE.33 ) THEN
!
!---        Gas conditions  ---
!
            IF( SG(2,N).GT.EPSL ) THEN
!
!---          Limit changes in pressure  ---
!
              DPX = MAX( 1.D+6,1.D-1*(PG(2,N)-PL(2,N)) )
              DPA = SIGN( MIN(ABS(DPX),ABS(DPA)),DPA )
              PSO(2,N) = MAX( (PSO(2,N)+DPA),0.D+0 )
              IF( PSO(2,N).LT.1.D-6 ) PSO(2,N) = 0.D+0
            ELSE
!
!---          Limit changes in pressure  ---
!
              DPX = 2.D-1*(PSO(2,N)+PATM)
              DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
              PSO(2,N) = PSO(2,N) + DPA
            ENDIF
          ELSE
!
!---        Zero negative corrections total nonaqueous 
!           O2 mole fraction  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              DPX = 2.5D-1
              DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
              ZMCA(2,N) = MAX( (ZMCA(2,N)+DPA),0.D+0 )
              IF( ZMCA(2,N).LT.EPSL ) ZMCA(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---      ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #34, #35, or #36
!
          IF( NPHZX.GE.34 .AND. NPHZX.LE.36 ) THEN
!
!---        Gas conditions  ---
!
            IF( SG(2,N).GT.EPSL ) THEN
!
!---          Limit changes in pressure  ---
!
              DPX = MAX( 1.D+6,1.D-1*(PG(2,N)-PL(2,N)) )
              DPO = SIGN( MIN(ABS(DPX),ABS(DPO)),DPO )
              PSO(2,N) = MAX( (PSO(2,N)+DPO),0.D+0 )
              IF( PSO(2,N).LT.1.D-6 ) PSO(2,N) = 0.D+0
            ELSE
!
!---          Limit changes in pressure  ---
!
              DPX = 2.D-1*(PSO(2,N)+PATM)
              DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
              PSO(2,N) = PSO(2,N) + DPO
            ENDIF
          ELSE
!
!---        Zero negative corrections total nonaqueous 
!           CH4 mole fraction  ---
!
            IF( ISLC(64).EQ.0 ) THEN
              DPX = 2.5D-1
              DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
              ZMCO(2,N) = MAX( (ZMCO(2,N)+DPO),0.D+0 )
              IF( ZMCO(2,N).LT.EPSL ) ZMCO(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---      ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #37, #38, or #39
!
          IF( NPHZX.GE.37 .AND. NPHZX.LE.39 ) THEN
!
!---        Gas conditions  ---
!
            IF( SG(2,N).GT.EPSL ) THEN
!
!---          Limit changes in pressure  ---
!
              DPX = MAX( 1.D+6,1.D-1*(PG(2,N)-PL(2,N)) )
              DPA = SIGN( MIN(ABS(DPX),ABS(DPN)),DPN )
              PSO(2,N) = MAX( (PSO(2,N)+DPN),0.D+0 )
              IF( PSO(2,N).LT.1.D-6 ) PSO(2,N) = 0.D+0
            ELSE
!
!---          Limit changes in pressure  ---
!
              DPX = 2.D-1*(PSO(2,N)+PATM)
              DPN = SIGN( MIN( ABS(DPX),ABS(DPN) ),DPN )
              PSO(2,N) = PSO(2,N) + DPN
            ENDIF
          ELSE
!
!---        Zero negative corrections total nonaqueous 
!           N2 mole fraction  ---
!
            IF( ISLC(65).EQ.0 ) THEN
              DPX = 2.5D-1
              DPN = SIGN( MIN( ABS(DPX),ABS(DPN) ),DPN )
              ZMCN(2,N) = MAX( (ZMCN(2,N)+DPN),0.D+0 )
              IF( ZMCN(2,N).LT.EPSL ) ZMCN(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---      YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #31, #34, or #37
!
          IF( NPHZX.EQ.31 .OR. NPHZX.EQ.34 .OR. NPHZX.EQ.37 ) THEN
!
!---        Limit changes in hydrate saturation  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHA = SIGN( MIN( ABS(DPHX),ABS(DPHA) ),DPHA )
              SH(2,N) = MIN( MAX( SH(2,N)+DPHA,0.D+0 ),0.999D+0 )
            ENDIF
          ELSE
!
!---        Limit changes in hydrate-gas O2 mole fraction of formers  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHA = SIGN( MIN( ABS(DPHX),ABS(DPHA) ),DPHA )
              YMHGA(2,N) = MIN( MAX( YMHGA(2,N)+DPHA,0.D+0 ),1.D+0 )
              IF( YMHGA(2,N).LT.EPSL ) YMHGA(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---      YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #32, #35, or #38
!
          IF( NPHZX.EQ.32 .OR. NPHZX.EQ.35 .OR. NPHZX.EQ.38 ) THEN
!
!---        Limit changes in hydrate saturation  ---
!
            IF( ISLC(64).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHO = SIGN( MIN( ABS(DPHX),ABS(DPHO) ),DPHO )
              SH(2,N) = MIN( MAX( SH(2,N)+DPHO,0.D+0 ),0.999D+0 )
            ENDIF
          ELSE
!
!---        Limit changes in hydrate-gas CH4 mole fraction of formers  ---
!
            IF( ISLC(64).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHO = SIGN( MIN( ABS(DPHX),ABS(DPHO) ),DPHO )
              YMHGO(2,N) = MIN( MAX( YMHGO(2,N)+DPHO,0.D+0 ),1.D+0 )
              IF( YMHGO(2,N).LT.EPSL ) YMHGO(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---      YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #33, #36, or #39
!
          IF( NPHZX.EQ.33 .OR. NPHZX.EQ.36 .OR. NPHZX.EQ.39 ) THEN
!
!---        Limit changes in hydrate saturation  ---
!
            IF( ISLC(65).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHN = SIGN( MIN( ABS(DPHX),ABS(DPHN) ),DPHN )
              SH(2,N) = MIN( MAX( SH(2,N)+DPHN,0.D+0 ),0.999D+0 )
            ENDIF
          ELSE
!
!---        Limit changes in hydrate-gas N2 mole fraction of formers  ---
!
            IF( ISLC(65).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHN = SIGN( MIN( ABS(DPHX),ABS(DPHN) ),DPHN )
              YMHGN(2,N) = MIN( MAX( YMHGN(2,N)+DPHN,0.D+0 ),1.D+0 )
              IF( YMHGN(2,N).LT.EPSL ) YMHGN(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---    Phase conditions #41 through #43 (40 series)  ---
!
!       SH > 0, SG + SN = 0, SL + SI < 1
!       YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #41
!       YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #42
!       YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #43
!
!       energy - temperature, C - T
!       water mass - pressure, Pa (gauge) - PSO
!       mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!       mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!       mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!       maximum hydrate mass component - hydrate saturation - SH
!       hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!       hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!       hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.41 .AND. NPHZX.LE.43 ) THEN
!
!---      Limit changes in pressure  ---
!
          DPX = 2.D-1*(PSO(2,N)+PATM)
          DPW = SIGN( MIN( ABS(DPX),ABS(DPW) ),DPW )
          PSO(2,N) = PSO(2,N) + DPW
!
!---      Zero negative corrections for zero O2 partial pressure  ---
!
          IF( ISLC(45).EQ.0 ) THEN
            IF( PVA(2,N)/EPSL.LT.EPSL .AND. DPA/EPSL.LT.EPSL ) THEN
              BLU(MPA) = 0.D+0
              DPA = 0.D+0
            ENDIF
            DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
            PVA(2,N) = MAX( (PVA(2,N)+DPA),0.D+0 )
            IF( PVA(2,N).LT.1.D-6 ) PVA(2,N) = 0.D+0
          ENDIF
!
!---      Zero negative corrections for zero CH4 partial pressure  ---
!
          IF( ISLC(64).EQ.0 ) THEN
            IF( PVO(2,N)/EPSL.LT.EPSL .AND. DPO/EPSL.LT.EPSL ) THEN
              BLU(MPO) = 0.D+0
              DPO = 0.D+0
            ENDIF
            DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
            PVO(2,N) = MAX( (PVO(2,N)+DPO),0.D+0 )
            IF( PVO(2,N).LT.1.D-6 ) PVO(2,N) = 0.D+0
          ENDIF
!
!---      Zero negative corrections for zero N2 partial pressure  ---
!
          IF( ISLC(65).EQ.0 ) THEN
            IF( PVN(2,N)/EPSL.LT.EPSL .AND. DPN/EPSL.LT.EPSL ) THEN
              BLU(MPN) = 0.D+0
              DPN = 0.D+0
            ENDIF
            DPN = SIGN( MIN( ABS(DPX),ABS(DPN) ),DPN )
            PVN(2,N) = MAX( (PVN(2,N)+DPN),0.D+0 )
            IF( PVN(2,N).LT.1.D-6 ) PVN(2,N) = 0.D+0
          ENDIF
!
!---      YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #41
!
          IF( NPHZX.EQ.41 ) THEN
!
!---        Limit changes in hydrate saturation  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHA = SIGN( MIN( ABS(DPHX),ABS(DPHA) ),DPHA )
              SH(2,N) = MIN( MAX( SH(2,N)+DPHA,0.D+0 ),0.999D+0 )
            ENDIF
          ELSE
!
!---        Limit changes in hydrate-gas O2 mole fraction of formers  ---
!
            IF( ISLC(45).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHA = SIGN( MIN( ABS(DPHX),ABS(DPHA) ),DPHA )
              YMHGA(2,N) = MIN( MAX( YMHGA(2,N)+DPHA,0.D+0 ),1.D+0 )
              IF( YMHGA(2,N).LT.EPSL ) YMHGA(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---      YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #42
!
          IF( NPHZX.EQ.42 ) THEN
!
!---        Limit changes in hydrate saturation  ---
!
            IF( ISLC(64).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHO = SIGN( MIN( ABS(DPHX),ABS(DPHO) ),DPHO )
              SH(2,N) = MIN( MAX( SH(2,N)+DPHO,0.D+0 ),0.999D+0 )
            ENDIF
          ELSE
!
!---        Limit changes in hydrate-gas CH4 mole fraction of formers  ---
!
            IF( ISLC(64).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHO = SIGN( MIN( ABS(DPHX),ABS(DPHO) ),DPHO )
              YMHGO(2,N) = MIN( MAX( YMHGO(2,N)+DPHO,0.D+0 ),1.D+0 )
              IF( YMHGO(2,N).LT.EPSL ) YMHGO(2,N) = 0.D+0
            ENDIF
          ENDIF
!
!---      YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #43
!
          IF( NPHZX.EQ.43 ) THEN
!
!---        Limit changes in hydrate saturation  ---
!
            IF( ISLC(65).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHN = SIGN( MIN( ABS(DPHX),ABS(DPHN) ),DPHN )
              SH(2,N) = MIN( MAX( SH(2,N)+DPHN,0.D+0 ),0.999D+0 )
            ENDIF
          ELSE
!
!---        Limit changes in hydrate-gas N2 mole fraction of formers  ---
!
            IF( ISLC(65).EQ.0 ) THEN
              DPHX = 1.D-1
              DPHN = SIGN( MIN( ABS(DPHX),ABS(DPHN) ),DPHN )
              YMHGN(2,N) = MIN( MAX( YMHGN(2,N)+DPHN,0.D+0 ),1.D+0 )
              IF( YMHGN(2,N).LT.EPSL ) YMHGN(2,N) = 0.D+0
            ENDIF
          ENDIF
        ENDIF
!
!---    Reduce time step for excessive changes in primary variables  ---
!
  200 CONTINUE
      IF( IERR.NE.0 ) THEN
        N = IERR
        ICNV = 1
        WRITE(ISC,'(10X,A)') '---  Excessive Primary Variable Change
     &---'
        WRITE(IWR,'(10X,A)') '---  Excessive Primary Variable Change
     &---'
        WRITE(ISC,'(4X,A,I6)') 'Node = ',N
        WRITE(IWR,'(4X,A,I6)') 'Node = ',N
        NPHZX = MOD(NPHAZ(2,N),100)
        WRITE(ISC,'(4X,2A)') 'Phase Condition = ',PH_CND(NPHZX)
        WRITE(IWR,'(4X,2A)') 'Phase Condition = ',PH_CND(NPHZX)
        WRITE(ISC,'(4X,A,1PE12.5)') 'DPT = ',DPT
        WRITE(IWR,'(4X,A,1PE12.5)') 'DPT = ',DPT
        WRITE(ISC,'(4X,A,1PE12.5)') 'DPW = ',DPW
        WRITE(IWR,'(4X,A,1PE12.5)') 'DPW = ',DPW
        IF( ISLC(45).EQ.0 ) THEN
          WRITE(ISC,'(4X,A,1PE12.5)') 'DPA = ',DPA
          WRITE(IWR,'(4X,A,1PE12.5)') 'DPA = ',DPA
        ENDIF
        IF( ISLC(64).EQ.0 ) THEN
          WRITE(ISC,'(4X,A,1PE12.5)') 'DPO = ',DPO
          WRITE(IWR,'(4X,A,1PE12.5)') 'DPO = ',DPO
        ENDIF
        IF( ISLC(65).EQ.0 ) THEN
          WRITE(ISC,'(4X,A,1PE12.5)') 'DPN = ',DPN
          WRITE(IWR,'(4X,A,1PE12.5)') 'DPN = ',DPN
        ENDIF
        IF( ISLC(45).EQ.0 ) THEN
          WRITE(ISC,'(4X,A,1PE12.5)') 'DPHA = ',DPHA
          WRITE(IWR,'(4X,A,1PE12.5)') 'DPHA = ',DPHA
        ENDIF
        IF( ISLC(64).EQ.0 ) THEN
          WRITE(ISC,'(4X,A,1PE12.5)') 'DPHO = ',DPHO
          WRITE(IWR,'(4X,A,1PE12.5)') 'DPHO = ',DPHO
        ENDIF
        IF( ISLC(65).EQ.0 ) THEN
          WRITE(ISC,'(4X,A,1PE12.5)') 'DPHN = ',DPHN
          WRITE(IWR,'(4X,A,1PE12.5)') 'DPHN = ',DPHN
        ENDIF
        WRITE(ISC,'(4X,A,1PE12.5)') 'Temperature, C = ',T(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Temperature, C = ',T(2,N)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Pressure, Pa = ',PSO(2,N)+PATM
        WRITE(IWR,'(4X,A,1PE12.5)') 'Pressure, Pa = ',PSO(2,N)+PATM
        WRITE(ISC,'(4X,A,1PE12.5)') 'Aqueous Saturation = ',SL(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Aqueous Saturation = ',SL(2,N)
        WRITE(ISC,'(4X,A,1PE12.5)') 'Ice Saturation = ',SI(2,N)
        WRITE(IWR,'(4X,A,1PE12.5)') 'Ice Saturation = ',SI(2,N)
!
!---    Core phase condition  ---
!
        NPHZX = MOD(NPHAZ(2,N),100)
!
!---    Phase condition #10 (10 series)  ---
!
!       SH = 0, SG + SN = 0, SL + SI = 1
!
!       energy - temperature, C - T
!       water mass - pressure, Pa (gauge) - PSO
!       mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!       mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!       mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!       hydrate O2 mass - hydrate O2 mass, kg - TMHA
!       hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!       hydrate N2 mass - hydrate N2 mass, kg - TMHN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        IF( NPHZX.EQ.10 ) THEN
          IF( ISLC(45).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'O2 vapor pressure, Pa absolute = ',PVA(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'O2 vapor pressure, Pa absolute = ',PVA(2,N)
          ENDIF
          IF( ISLC(64).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'CH4 vapor pressure, Pa absolute = ',PVO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'CH4 vapor pressure, Pa absolute = ',PVO(2,N)
          ENDIF
          IF( ISLC(65).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'N2 vapor pressure, Pa absolute = ',PVN(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'N2 vapor pressure, Pa absolute = ',PVN(2,N)
          ENDIF
          IF( ISLC(45).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'O2 Hydrate Mass = ',TMHA(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'O2 Hydrate Mass = ',TMHA(2,N)
          ENDIF
          IF( ISLC(64).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'CH4 Hydrate Mass = ',TMHO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'CH4 Hydrate Mass = ',TMHO(2,N)
          ENDIF
          IF( ISLC(65).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'N2 Hydrate Mass = ',TMHN(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'N2 Hydrate Mass = ',TMHN(2,N)
          ENDIF
!
!---    Phase conditions 21 through 23 (20 series)  ---
!
!       SH = 0, SG + SN > 0, SL + SI < 1
!       ZMCA = MAX( ZMCA,ZMCO,ZMCN ) #21
!       ZMCO = MAX( ZMCA,ZMCO,ZMCN ) #22
!       ZMCN = MAX( ZMCA,ZMCO,ZMCN ) #23
!
!       energy - temperature, C - T
!       water mass - aqueous + ice saturation - POSM
!       maximum mobile mass component - pressure, Pa (gauge) - PSO
!       mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!       mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!       mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!       hydrate O2 mass - hydrate O2 mass, kg - TMHA
!       hydrate CH4 mass - hydrate CH4 mass, kg - TMHO
!       hydrate N2 mass - hydrate N2 mass, kg - TMHN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.21 .AND. NPHZX.LE.23  ) THEN
          IF( ISLC(45).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous O2 Mole Fraction = ',ZMCA(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous O2 Mole Fraction = ',ZMCA(2,N)
          ENDIF
          IF( ISLC(64).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous CH4 Mole Fraction = ',ZMCO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous CH4 Mole Fraction = ',ZMCO(2,N)
          ENDIF
          IF( ISLC(65).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous N2 Mole Fraction = ',ZMCN(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous N2 Mole Fraction = ',ZMCN(2,N)
          ENDIF
          IF( ISLC(45).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'O2 Hydrate Mass = ',TMHA(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'O2 Hydrate Mass = ',TMHA(2,N)
          ENDIF
          IF( ISLC(64).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'CH4 Hydrate Mass = ',TMHO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'CH4 Hydrate Mass = ',TMHO(2,N)
          ENDIF
          IF( ISLC(65).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'N2 Hydrate Mass = ',TMHN(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'N2 Hydrate Mass = ',TMHN(2,N)
          ENDIF
!
!---    Phase conditions #31 through #39 (30 series)  ---
!
!       SH > 0, SG + SN > 0, SL + SI < 1
!       ZMCA = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #31
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #32
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #33
!       ZMCO = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #34
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #35
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #36
!       ZMCN = MAX( ZMCA,ZMCO,ZMCN )
!         YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #37
!         YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #38
!         YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #39
!
!       energy - temperature, C - T
!       water mass - aqueous + ice saturation, POSM
!       maximum mobile mass component - pressure, Pa (gauge) - PSO
!       mobile O2 mass - total nonaqueous CH4 mole fraction - ZMCA
!       mobile CH4 mass - total nonaqueous CH4 mole fraction - ZMCO
!       mobile N2 mass - total nonaqueous N2 mole fraction - ZMCN
!       maximum hydrate mass component - hydrate saturation - SH
!       hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!       hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!       hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
      ELSEIF( NPHZX.GE.31 .AND. NPHZX.LE.39 ) THEN
          IF( ISLC(45).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous O2 Mole Fraction = ',ZMCA(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous O2 Mole Fraction = ',ZMCA(2,N)
          ENDIF
          IF( ISLC(64).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous CH4 Mole Fraction = ',ZMCO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous CH4 Mole Fraction = ',ZMCO(2,N)
          ENDIF
          IF( ISLC(65).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous N2 Mole Fraction = ',ZMCN(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Total Nonaqueous N2 Mole Fraction = ',ZMCN(2,N)
          ENDIF
          WRITE(ISC,'(4X,A,1PE12.5)')
     &      'Hydrate Saturation = ',SH(2,N)
          WRITE(IWR,'(4X,A,1PE12.5)')
     &      'Hydrate Saturation = ',SH(2,N)
          IF( ISLC(45).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Hydrate Gas O2 Mole Fraction of Formers = ',YMHGA(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Hydrate Gas O2 Mole Fraction of Formers = ',YMHGA(2,N)
          ENDIF
          IF( ISLC(64).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Hydrate Gas CH4 Mole Fraction of Formers = ',YMHGO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Hydrate Gas CH4 Mole Fraction of Formers = ',YMHGO(2,N)
          ENDIF
          IF( ISLC(65).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Hydrate Gas N2 Mole Fraction of Formers = ',YMHGN(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Hydrate Gas N2 Mole Fraction of Formers = ',YMHGN(2,N)
          ENDIF
!
!---    Phase conditions #41 through #43 (40 series)  ---
!
!       SH > 0, SG + SN = 0, SL + SI < 1
!       YMHGA = MAX( YMHGA,YMHGO,YMHGN ) #41
!       YMHGO = MAX( YMHGA,YMHGO,YMHGN ) #42
!       YMHGN = MAX( YMHGA,YMHGO,YMHGN ) #43
!
!       energy - temperature, C - T
!       water mass - pressure, Pa (gauge) - PSO
!       mobile O2 mass - O2 vapor pressure, Pa absolute - PVA
!       mobile CH4 mass - CH4 vapor pressure, Pa absolute - PVO
!       mobile N2 mass - N2 vapor pressure, Pa absolute - PVN
!       maximum hydrate mass component - hydrate saturation - SH
!       hydrate O2 mass - hyd-gas CH4 mole frac. of formers - YMHGA
!       hydrate CH4 mass - hyd-gas CH4 mole frac. of formers - YMHGO
!       hydrate N2 mass - hyd-gas N2 mole frac. of formers - YMHGN
!       NaCl mass - total NaCl aqueous mass fraction - YLS
!
        ELSEIF( NPHZX.GE.41 .AND. NPHZX.LE.43 ) THEN
          IF( ISLC(45).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'O2 vapor pressure, Pa absolute = ',PVA(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'O2 vapor pressure, Pa absolute = ',PVA(2,N)
          ENDIF
          IF( ISLC(64).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'CH4 vapor pressure, Pa absolute = ',PVO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'CH4 vapor pressure, Pa absolute = ',PVO(2,N)
          ENDIF
          IF( ISLC(65).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'N2 vapor pressure, Pa absolute = ',PVN(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'N2 vapor pressure, Pa absolute = ',PVN(2,N)
          ENDIF
          WRITE(ISC,'(4X,A,1PE12.5)')
     &      'Hydrate Saturation = ',SH(2,N)
          WRITE(IWR,'(4X,A,1PE12.5)')
     &      'Hydrate Saturation = ',SH(2,N)
          IF( ISLC(45).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Hydrate Gas O2 Mole Fraction of Formers = ',YMHGA(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Hydrate Gas O2 Mole Fraction of Formers = ',YMHGA(2,N)
          ENDIF
          IF( ISLC(64).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Hydrate Gas CH4 Mole Fraction of Formers = ',YMHGO(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Hydrate Gas CH4 Mole Fraction of Formers = ',YMHGO(2,N)
          ENDIF
          IF( ISLC(65).EQ.0 ) THEN
            WRITE(ISC,'(4X,A,1PE12.5)')
     &        'Hydrate Gas N2 Mole Fraction of Formers = ',YMHGN(2,N)
            WRITE(IWR,'(4X,A,1PE12.5)')
     &        'Hydrate Gas N2 Mole Fraction of Formers = ',YMHGN(2,N)
          ENDIF
        ENDIF
        WRITE(ISC,'(4X,A,1PE12.5,A,I6)')
     &    'Total-Salt Aqu. Mass Fraction = ',YLS(2,N)
        WRITE(IWR,'(4X,A,1PE12.5,A,I6)')
     &    'Total-Salt Aqu. Mass Fraction = ',YLS(2,N)
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
          DO 400 N = 1,NFLD
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            PSO(2,N) = PSO(1,N)
            ZMCA(2,N) = ZMCA(1,N)
            ZMCO(2,N) = ZMCO(1,N)
            ZMCN(2,N) = ZMCN(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            POSM(2,N) = POSM(1,N)
            SH(2,N) = SH(1,N)
            SI(2,N) = SI(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            YMGA(2,N) = YMGA(1,N)
            YMGO(2,N) = YMGO(1,N)
            YMGN(2,N) = YMGN(1,N)
            YMHGA(2,N) = YMHGA(1,N)
            YMHGO(2,N) = YMHGO(1,N)
            YMHGN(2,N) = YMHGN(1,N)
            TMHA(2,N) = TMHA(1,N)
            TMHO(2,N) = TMHO(1,N)
            TMHN(2,N) = TMHN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  400     CONTINUE
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
          DO 410 N = 1,NFLD
            T(2,N) = T(1,N)
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            PSO(2,N) = PSO(1,N)
            ZMCA(2,N) = ZMCA(1,N)
            ZMCO(2,N) = ZMCO(1,N)
            ZMCN(2,N) = ZMCN(1,N)
            SG(2,N) = SG(1,N)
            SL(2,N) = SL(1,N)
            POSM(2,N) = POSM(1,N)
            SH(2,N) = SH(1,N)
            SI(2,N) = SI(1,N)
            SN(2,N) = SN(1,N)
            YLS(2,N) = YLS(1,N)
            YMGA(2,N) = YMGA(1,N)
            YMGO(2,N) = YMGO(1,N)
            YMGN(2,N) = YMGN(1,N)
            YMHGA(2,N) = YMHGA(1,N)
            YMHGO(2,N) = YMHGO(1,N)
            YMHGN(2,N) = YMHGN(1,N)
            TMHA(2,N) = TMHA(1,N)
            TMHO(2,N) = TMHO(1,N)
            TMHN(2,N) = TMHN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  410     CONTINUE
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          ICNV = 4
        ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT_HYD_O2_CH4_N2 group
!
      RETURN
      END

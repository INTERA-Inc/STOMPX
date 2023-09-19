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
!     Water-Oil-Air Mode
!
!     This engineering program numerically simulates thermal
!     and hydrologic transport phenomena in variably saturated
!     subsurface environments, contaminated with a water immiscible
!     volatile organic compound.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
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
      USE CONST
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
      ISUB_LOG = 1
      SUB_LOG(1) = 'STOMP5'
      ICODE = 5
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
      CALL RDINPT5
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
      CALL CHK5
!
!---  For geomechanics set k iterate value of pore pressure  ---
!
      IF( ISLC(50).NE.0 ) THEN
        INDX = 2
        CALL PRESS_GM( INDX )
      ENDIF
!
!---  Compute primary variable increments  ---
!
      CALL INCRM5
!
!---  Saturation, relative permeability, porosity, and tortuosity  ---
!
      CALL PROP5
!
!---  Thermodynamic properties and equations of state  ---
!
      CALL TMPR5
!
!---  Compute boundary saturation, relative permeability, and
!     thermodynamic properties  ---
!
      CALL BCP5
!
!---  Compute Jacobian matrix pointers  ---
!
      CALL JCBP
!
!---  For geomechanics simulations compute Jacobian matrix pointers  --
!
      IF( ISLC(50).NE.0 .AND. ISLC(67).EQ.0 ) CALL JCBP_GM
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
!---  Compute initial fluxes  ---
!
      ISVF = 1
!
!---  Compute aqueous-phase volumetric flux (interior surfaces)  ---
!
      CALL DRCVL
!
!---  Compute gas-phase volumetric flux (interior surfaces)  ---
!
      CALL DRCVG
!
!---  Compute NAPL volumetric flux (interior surfaces)  ---
!
      CALL DRCVN
!
!---  Compute air diffusion flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLA
!
!---  Compute oil diffusion flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLO
!
!---  Compute oil diffusion flux through the napl phase
!     (interior surfaces)  ---
!
      CALL DFFNO
!
!---  Compute water mass flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLW
!
!---  Compute water vapor diffusion flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGW
!
!---  Compute oil vapor diffusion flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGO
!
!---  Compute air mass flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGA
!
!---  Compute aqueous-phase, and NAPL volumetric flux
!     (boundary surfaces)  ---
!
      CALL BCF5
!
!---  Average gas and aqueous pore velocities  ---
!
      IF( ISLC(46).GE.1 ) CALL AV5
!
!---  Reset flux increment range  ---
!
      ISVF = 2*ISVC+1
!
!---  Surface flux integrator for zero time step  ---
!
      DTX = DT
      DT = 0.D+0
      CALL SFIN
      DT = DTX
!
!---  Compute initial solute concentrations  ---
!
      CALL CISC5
!
!---  New Time Step ---
!
  100 CONTINUE
!
!---  Update porosity and permeability in response to geomechanical
!     stress  ---
!
      IF( ISLC(50).NE.0 ) THEN
        CALL PORSTY_GM
        CALL PERMRF_GM
      ENDIF
!
!---  Load old time step arrays  ---
!
      CALL LDO5
!
!---  Load old time step arrays for the volumetric stress
!     and pore pressure  ---
!
      IF( ISLC(50).NE.0 ) THEN
        INDX = 1
        CALL LD_GM( INDX )
      ENDIF
!
!---  Compute trapping number  ---
!
      CALL TRPNL5
!
!---  Reference node(s) output  ---
!
      IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 .OR.
     &  MOD( (NSTEP-NRST),IFQO ).EQ.0) CALL REFNOD
!
!---  End of initial conditions simulations  ---
!
      IF( IEO.EQ.3 ) GOTO 900
!
!---  Stop simulation if simulation time exceeds limit  ---
!
      IF( ABS(TMMX-TM).LE.1.D-6 ) THEN
        INDX = 1
        CHMSG = 'Simulation Stopped: Simulation Time Limit'
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
      IF( ABS(TMPR-TM)/(ABS(TM)+EPSL).LE.EPSL ) THEN
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
        CALL BCP5
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
!---  Courant number limit check  ---
!
      IF( CRNTMXC.GT.EPSL ) THEN
        ISVF = 1
        CALL DRCVL
        CALL DRCVG
        CALL DRCVN
        CALL BCF5
        ISVF = 2*ISVC+1
        ICNV = 3
        CALL CRN_G
        CALL CRN_L
        CALL CRN_N
!
!---    Cut time step  ---
!
        IF( ICNV.EQ.1 ) GOTO 200
      ENDIF
!
!---  Compute boundary saturation, relative permeability, and
!     thermodynamic properties  ---
!
      CALL BCP5
!
!---  Compute source contributions  ---
!
      CALL SORC5
!
!---  Compute aqueous volumetric flux (interior surfaces)  ---
!
      CALL DRCVL
!
!---  Compute gas volumetric flux (interior surfaces)  ---
!
      CALL DRCVG
!
!---  Compute NAPL volumetric flux (interior surfaces)  ---
!
      CALL DRCVN
!
!---  Compute air mass flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLA
!
!---  Compute oil mass flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLO
!
!---  Compute oil diffusion flux through the napl phase
!     (interior surfaces)  ---
!

      CALL DFFNO
!
!---  Compute water mass flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLW
!
!---  Compute water mass flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGW
!
!---  Compute oil mass flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGO
!
!---  Compute air mass flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGA
!
!---  Compute aqueous, gas, and NAPL phase volumetric fluxes,
!     and water, air, and oil component mass fluxes.
!     (boundary surfaces)  ---
!
      CALL BCF5
!
!---  Kinetic volatilzation rates  ---
!
      IF( ISLC(46).GE.1 ) CALL KV5
!
!---  Zero Jacobian matrix  ---
!
#ifdef petsc
      IF( ILES.EQ.5 ) CALL MatZeroEntries(F_MAT,IERR)
#endif
      INDX = 0
      CALL JCBZ( ISVC,MUC,MLC,MKC,INDX )
!
!---  Load Jacobian matrix for the water equation
!     (zero flux boundary)  ---
!
      CALL JCBW5
!
!---  Load Jacobian matrix for the air equation
!     (zero flux boundary)  ---
!
      CALL JCBA5
!
!---  Kinetic volatilization option  ---
!
      IF( ISLC(46).GE.1 ) THEN
!
!---    Load Jacobian matrix for the NAPL oil equation
!       (zero flux boundary)  ---
!
        CALL JCBNO5
!
!---    Load Jacobian matrix for the gas-aqueous oil equation
!       (zero flux boundary)  ---
!
        CALL JCBDO5
      ELSE
!
!---    Load Jacobian matrix for the oil equation
!       (zero flux boundary)  ---
!
        CALL JCBO5
      ENDIF
!
!---  Modify the Jacobian matrix for the boundary conditions  ---
!
      CALL BCJ5
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
      CALL UPDT5
!
!---  Compute convergence from maximum relative residuals  ---
!
      CALL RSDL5
!
!---  Compute primary variable increments, saturation,
!     relative permeability, porosity, tortuosity,
!     thermodynamic properties for interior nodes,
!     except immediately after a new time step  ---
!
      CALL INCRM5
      CALL PROP5
      CALL TMPR5
!
!---  For geomechanics simulations alter permeability with
!     porosity  --
!
      IF( ISLC(50).NE.0 ) CALL PERMRF_GM
      GOTO( 200,300,600,900 ) ICNV
!
!---  Beginning of transport equation solution  ---
!
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
          GOTO 190
        ENDIF
      ENDIF
!
!---  Temporarily set flux index  ---
!
      ISVF = 1
!
!---  Compute aqueous volumetric flux (interior surfaces)  ---
!
      CALL DRCVL
!
!---  Compute gas volumetric flux (interior surfaces)  ---
!
      CALL DRCVG
!
!---  Compute NAPL volumetric flux (interior surfaces)  ---
!
      CALL DRCVN
!
!---  Compute air mass flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLA
!
!---  Compute oil mass flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLO
!
!---  Compute oil diffusion flux through the napl phase
!     (interior surfaces)  ---
!
      CALL DFFNO
!
!---  Compute water mass flux through the aqueous phase
!     (interior surfaces)  ---
!
      CALL DFFLW
!
!---  Compute water mass flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGW
!
!---  Compute oil mass flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGO
!
!---  Compute air mass flux through the gas phase
!     (interior surfaces)  ---
!
      CALL DFFGA
!
!---  Compute aqueous, gas, and NAPL phase volumetric fluxes,
!     and water, air, and oil component mass fluxes.
!     (boundary surfaces)  ---
!
      CALL BCF5
!
!---  Average gas and aqueous pore velocities  ---
!
      IF( ISLC(46).GE.1 ) CALL AV5
!
!---  Compute Local Courant Numbers  ---
!
      IF( ICRNT.EQ.1 ) CALL CRNTNB
!
!---  Reset flux index  ---
!
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
          CALL SPRP5( NSL )
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
          CALL SORT5( NSL )
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
!---      Load Jacobian matrix (NAPL transport)  ---
!
          CALL SJCBN( NSL )
!
!---      Modify Jacobian matrix for boundary conditions ---
!
          CALL SBND5( NSL )
!
!---  Linear equation solver  ---
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
!---      Compute solute NAPL fluxes (interior nodes)  ---
!
          CALL SFXN( NSL )
!
!---      Compute solute NAPL fluxes (boundary surfaces)  ---
!
          CALL SFXNB( NSL )
!
!---      Integrate solute sources ---
!
          CALL SORIT5( NSL )
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
!---  Correct aqueous liquid density and viscosity for electrolyte
!     solute concentration  ---
!
      IF( ISLC(16).EQ.1 ) CALL ELC5
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
      WRITE(IWR,'(/,A)') '---  End of STOMP Simulation  ---'
      WRITE(ISC,'(/,A)') '---  End of STOMP Simulation  ---'
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
      SUBROUTINE AV5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute average gas and aqueous velocities for the kinetic
!     volatilization option.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 May 2008.
!     Last Modified by MD White, PNNL, 23 May 2008.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVD
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
      SUB_LOG(ISUB_LOG) = '/AV5'
!
!---  Loop over nodes  ---
!
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 100
        AVPVG(N) = 0.D+0
        AVPVL(N) = 0.D+0
!
!---    Skip calculation for nodes without NAPL  ---
!
        IF( SN(2,N).LT.EPSL ) GOTO 100
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX+1
        NQY = NPY+IFLD
        NQZ = NPZ+IJFLD
!
!---    Average pore gas velocity, m/s,
!       skipping calculation for nodes without
!       gas  ---
!
        IF( SG(2,N).LT.EPSL ) GOTO 50
        NC = 0
        IF( ABS(UG(1,NPX)).GT.EPSL ) THEN
          AVPVG(N) = AVPVG(N) + ABS(UG(1,NPX))
          NC = NC + 1
        ENDIF
        IF( ABS(UG(1,NQX)).GT.EPSL ) THEN
          AVPVG(N) = AVPVG(N) + ABS(UG(1,NQX))
          NC = NC + 1
        ENDIF
        IF( ABS(VG(1,NPY)).GT.EPSL ) THEN
          AVPVG(N) = AVPVG(N) + ABS(VG(1,NPY))
          NC = NC + 1
        ENDIF
        IF( ABS(VG(1,NQY)).GT.EPSL ) THEN
          AVPVG(N) = AVPVG(N) + ABS(VG(1,NQY))
          NC = NC + 1
        ENDIF
        IF( ABS(WG(1,NPZ)).GT.EPSL ) THEN
          AVPVG(N) = AVPVG(N) + ABS(WG(1,NPZ))
          NC = NC + 1
        ENDIF
        IF( ABS(WG(1,NQZ)).GT.EPSL ) THEN
          AVPVG(N) = AVPVG(N) + ABS(WG(1,NQZ))
          NC = NC + 1
        ENDIF
        IF( NC.GT.0 ) THEN
          AVPVG(N) = AVPVG(N)/REAL(NC)/PORD(2,N)/SG(2,N)
        ENDIF
   50 CONTINUE
!
!---    Average pore aqueous velocity, m/s,
!       skipping calculation for nodes without
!       aqueous  ---
!
        IF( SL(2,N).LT.EPSL ) GOTO 100
        AVPVL(N) = 0.D+0
        NC = 0
        IF( ABS(UL(1,NPX)).GT.EPSL ) THEN
          AVPVL(N) = AVPVL(N) + ABS(UL(1,NPX))
          NC = NC + 1
        ENDIF
        IF( ABS(UL(1,NQX)).GT.EPSL ) THEN
          AVPVL(N) = AVPVL(N) + ABS(UL(1,NQX))
          NC = NC + 1
        ENDIF
        IF( ABS(VL(1,NPY)).GT.EPSL ) THEN
          AVPVL(N) = AVPVL(N) + ABS(VL(1,NPY))
          NC = NC + 1
        ENDIF
        IF( ABS(VL(1,NQY)).GT.EPSL ) THEN
          AVPVL(N) = AVPVL(N) + ABS(VL(1,NQY))
          NC = NC + 1
        ENDIF
        IF( ABS(WL(1,NPZ)).GT.EPSL ) THEN
          AVPVL(N) = AVPVL(N) + ABS(WL(1,NPZ))
          NC = NC + 1
        ENDIF
        IF( ABS(WL(1,NQZ)).GT.EPSL ) THEN
          AVPVL(N) = AVPVL(N) + ABS(WL(1,NQZ))
          NC = NC + 1
        ENDIF
        IF( NC.GT.0 ) THEN
          AVPVL(N) = AVPVL(N)/REAL(NC)/PORD(2,N)/SL(2,N)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of AV5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCF5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute boundary surface fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     Last Modified by ML Rockhold, PNNL, 29 January 2003.
!     Last Modified by M Oostrom, PNNL, 20 March 2003.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXD
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
      SUB_LOG(ISUB_LOG) = '/BCF5'
!
!---  Zero boundary fluxes  ---
!
      DO 70 NB = 1,NBC
        N = IBCN(NB)
        IF( IBCD(NB).EQ.-3 ) THEN
          DO 10 M = 1,ISVF
            WL(M,NSZ(N)) = 0.D+0
            WG(M,NSZ(N)) = 0.D+0
            WN(M,NSZ(N)) = 0.D+0
            WDGO(M,NSZ(N)) = 0.D+0
            WDGW(M,NSZ(N)) = 0.D+0
            WDLO(M,NSZ(N)) = 0.D+0
            WDLA(M,NSZ(N)) = 0.D+0
            WLA(M,NSZ(N)) = 0.D+0
            WGA(M,NSZ(N)) = 0.D+0
            WLO(M,NSZ(N)) = 0.D+0
            WGO(M,NSZ(N)) = 0.D+0
            WNA(M,NSZ(N)) = 0.D+0
            WLW(M,NSZ(N)) = 0.D+0
            WGW(M,NSZ(N)) = 0.D+0
   10     CONTINUE
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          DO 20 M = 1,ISVF
            VL(M,NSY(N)) = 0.D+0
            VG(M,NSY(N)) = 0.D+0
            VN(M,NSY(N)) = 0.D+0
            VDGO(M,NSY(N)) = 0.D+0
            VDGW(M,NSY(N)) = 0.D+0
            VDLO(M,NSY(N)) = 0.D+0
            VDLA(M,NSY(N)) = 0.D+0
            VLA(M,NSY(N)) = 0.D+0
            VGA(M,NSY(N)) = 0.D+0
            VLO(M,NSY(N)) = 0.D+0
            VGO(M,NSY(N)) = 0.D+0
            VNA(M,NSY(N)) = 0.D+0
            VLW(M,NSY(N)) = 0.D+0
            VGW(M,NSY(N)) = 0.D+0
   20     CONTINUE
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          DO 30 M = 1,ISVF
            UL(M,NSX(N)) = 0.D+0
            UG(M,NSX(N)) = 0.D+0
            UN(M,NSX(N)) = 0.D+0
            UDGO(M,NSX(N)) = 0.D+0
            UDGW(M,NSX(N)) = 0.D+0
            UDLO(M,NSX(N)) = 0.D+0
            UDLA(M,NSX(N)) = 0.D+0
            ULA(M,NSX(N)) = 0.D+0
            UGA(M,NSX(N)) = 0.D+0
            ULO(M,NSX(N)) = 0.D+0
            UGO(M,NSX(N)) = 0.D+0
            UNA(M,NSX(N)) = 0.D+0
            ULW(M,NSX(N)) = 0.D+0
            UGW(M,NSX(N)) = 0.D+0
   30     CONTINUE
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          DO 40 M = 1,ISVF
            UL(M,NSX(N)+1) = 0.D+0
            UG(M,NSX(N)+1) = 0.D+0
            UN(M,NSX(N)+1) = 0.D+0
            UDGO(M,NSX(N)+1) = 0.D+0
            UDGW(M,NSX(N)+1) = 0.D+0
            UDLO(M,NSX(N)+1) = 0.D+0
            UDLA(M,NSX(N)+1) = 0.D+0
            ULA(M,NSX(N)+1) = 0.D+0
            UGA(M,NSX(N)+1) = 0.D+0
            ULO(M,NSX(N)+1) = 0.D+0
            UGO(M,NSX(N)+1) = 0.D+0
            UNA(M,NSX(N)+1) = 0.D+0
            ULW(M,NSX(N)+1) = 0.D+0
            UGW(M,NSX(N)+1) = 0.D+0
   40     CONTINUE
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          DO 50 M = 1,ISVF
            VL(M,NSY(N)+IFLD) = 0.D+0
            VG(M,NSY(N)+IFLD) = 0.D+0
            VN(M,NSY(N)+IFLD) = 0.D+0
            VDGO(M,NSY(N)+IFLD) = 0.D+0
            VDGW(M,NSY(N)+IFLD) = 0.D+0
            VDLO(M,NSY(N)+IFLD) = 0.D+0
            VDLA(M,NSY(N)+IFLD) = 0.D+0
            VLA(M,NSY(N)+IFLD) = 0.D+0
            VGA(M,NSY(N)+IFLD) = 0.D+0
            VLO(M,NSY(N)+IFLD) = 0.D+0
            VGO(M,NSY(N)+IFLD) = 0.D+0
            VNA(M,NSY(N)+IFLD) = 0.D+0
            VLW(M,NSY(N)+IFLD) = 0.D+0
            VGW(M,NSY(N)+IFLD) = 0.D+0
   50     CONTINUE
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          DO 60 M = 1,ISVF
            WL(M,NSZ(N)+IJFLD) = 0.D+0
            WG(M,NSZ(N)+IJFLD) = 0.D+0
            WN(M,NSZ(N)+IJFLD) = 0.D+0
            WDGO(M,NSZ(N)+IJFLD) = 0.D+0
            WDGW(M,NSZ(N)+IJFLD) = 0.D+0
            WDLO(M,NSZ(N)+IJFLD) = 0.D+0
            WDLA(M,NSZ(N)+IJFLD) = 0.D+0
            WLA(M,NSZ(N)+IJFLD) = 0.D+0
            WGA(M,NSZ(N)+IJFLD) = 0.D+0
            WLO(M,NSZ(N)+IJFLD) = 0.D+0
            WGO(M,NSZ(N)+IJFLD) = 0.D+0
            WNA(M,NSZ(N)+IJFLD) = 0.D+0
            WLW(M,NSZ(N)+IJFLD) = 0.D+0
            WGW(M,NSZ(N)+IJFLD) = 0.D+0
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
!---  Bottom boundary
!
        IF( IBCD(NB).EQ.-3 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 110 M = 1,ISVF
              WL(M,NSZ(N)) = BCX(2)
  110       CONTINUE
            CALL DFFLOB( N,NB )
            CALL DFFLAB( N,NB )
            CALL DFFLWB( N,NB )
          ELSEIF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL DRCVLB( N,NB )
            CALL DFFLOB( N,NB )
            CALL DFFLAB( N,NB )
            CALL DFFLWB( N,NB )
          ENDIF
!
!---  Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 111 M = 1,ISVF
              WG(M,NSZ(N)) = BCX(3)
  111       CONTINUE
            CALL DFFGOB( N,NB )
            CALL DFFGWB( N,NB )
            CALL DFFGAB( N,NB )
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVGB( N,NB )
            CALL DFFGOB( N,NB )
            CALL DFFGWB( N,NB )
            CALL DFFGAB( N,NB )
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 112 M = 1,ISVF
              WN(M,NSZ(N)) = BCX(4)
  112       CONTINUE
            CALL DFFNOB( N,NB )
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVNB( N,NB )
            CALL DFFNOB( N,NB )
          ENDIF
!
!---  South boundary
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 120 M = 1,ISVF
              VL(M,NSY(N)) = BCX(2)
  120       CONTINUE
            CALL DFFLOS( N,NB )
            CALL DFFLAS( N,NB )
            CALL DFFLWS( N,NB )
          ELSEIF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL DRCVLS( N,NB )
            CALL DFFLOS( N,NB )
            CALL DFFLAS( N,NB )
            CALL DFFLWS( N,NB )
          ENDIF
!
!---  Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 121 M = 1,ISVF
              VG(M,NSY(N)) = BCX(3)
  121       CONTINUE
            CALL DFFGOS( N,NB )
            CALL DFFGWS( N,NB )
            CALL DFFGAS( N,NB )
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVGS( N,NB )
            CALL DFFGOS( N,NB )
            CALL DFFGWS( N,NB )
            CALL DFFGAS( N,NB )
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 122 M = 1,ISVF
              VN(M,NSY(N)) = BCX(4)
  122       CONTINUE
            CALL DFFNOS( N,NB )
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVNS( N,NB )
            CALL DFFNOS( N,NB )
          ENDIF
!
!---  West boundary
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 130 M = 1,ISVF
              UL(M,NSX(N)) = BCX(2)
  130       CONTINUE
            CALL DFFLOW( N,NB )
            CALL DFFLAW( N,NB )
            CALL DFFLWW( N,NB )
          ELSEIF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL DRCVLW( N,NB )
            CALL DFFLOW( N,NB )
            CALL DFFLAW( N,NB )
            CALL DFFLWW( N,NB )
          ENDIF
!
!---  Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 131 M = 1,ISVF
              UG(M,NSX(N)) = BCX(3)
  131       CONTINUE
            CALL DFFGOW( N,NB )
            CALL DFFGWW( N,NB )
            CALL DFFGAW( N,NB )
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVGW( N,NB )
            CALL DFFGOW( N,NB )
            CALL DFFGWW( N,NB )
            CALL DFFGAW( N,NB )
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 132 M = 1,ISVF
              UN(M,NSX(N)) = BCX(4)
  132       CONTINUE
            CALL DFFNOW( N,NB )
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVNW( N,NB )
            CALL DFFNOW( N,NB )
          ENDIF
!
!---  East boundary
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 140 M = 1,ISVF
              UL(M,NSX(N)+1) = BCX(2)
  140       CONTINUE
            CALL DFFLOE( N,NB )
            CALL DFFLAE( N,NB )
            CALL DFFLWE( N,NB )
          ELSEIF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL DRCVLE( N,NB )
            CALL DFFLOE( N,NB )
            CALL DFFLAE( N,NB )
            CALL DFFLWE( N,NB )
          ENDIF
!
!---  Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 141 M = 1,ISVF
              UG(M,NSX(N)+1) = BCX(3)
  141       CONTINUE
            CALL DFFGOE( N,NB )
            CALL DFFGWE( N,NB )
            CALL DFFGAE( N,NB )
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVGE( N,NB )
            CALL DFFGOE( N,NB )
            CALL DFFGWE( N,NB )
            CALL DFFGAE( N,NB )
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 142 M = 1,ISVF
              UN(M,NSX(N)+1) = BCX(4)
  142       CONTINUE
            CALL DFFNOE( N,NB )
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVNE( N,NB )
            CALL DFFNOE( N,NB )
          ENDIF
!
!---  North boundary
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 150 M = 1,ISVF
              VL(M,NSY(N)+IFLD) = BCX(2)
  150       CONTINUE
            CALL DFFLON( N,NB )
            CALL DFFLAN( N,NB )
            CALL DFFLWN( N,NB )
          ELSEIF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL DRCVLN( N,NB )
            CALL DFFLON( N,NB )
            CALL DFFLAN( N,NB )
            CALL DFFLWN( N,NB )
          ENDIF
!
!---  Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 151 M = 1,ISVF
              VG(M,NSY(N)+IFLD) = BCX(3)
  151       CONTINUE
            CALL DFFGON( N,NB )
            CALL DFFGWN( N,NB )
            CALL DFFGAN( N,NB )
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVGN( N,NB )
            CALL DFFGON( N,NB )
            CALL DFFGWN( N,NB )
            CALL DFFGAN( N,NB )
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 152 M = 1,ISVF
              VN(M,NSY(N)+IFLD) = BCX(4)
  152       CONTINUE
            CALL DFFNON( N,NB )
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVNN( N,NB )
            CALL DFFNON( N,NB )
          ENDIF
!
!---  Top boundary
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---  Aqueous Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(1,NB).EQ.2 ) THEN
            DO 160 M = 1,ISVF
              WL(M,NSZ(N)+IJFLD) = BCX(2)
  160       CONTINUE
            CALL DFFLOT( N,NB )
            CALL DFFLAT( N,NB )
            CALL DFFLWT( N,NB )
          ELSEIF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL DRCVLT( N,NB )
            CALL DFFLOT( N,NB )
            CALL DFFLAT( N,NB )
            CALL DFFLWT( N,NB )
          ENDIF
!
!---  Gas Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(2,NB).EQ.2 ) THEN
            DO 161 M = 1,ISVF
              WG(M,NSZ(N)+IJFLD) = BCX(3)
  161       CONTINUE
            CALL DFFGOT( N,NB )
            CALL DFFGWT( N,NB )
            CALL DFFGAT( N,NB )
          ELSEIF( IBCT(2,NB).NE.3 ) THEN
            CALL DRCVGT( N,NB )
            CALL DFFGOT( N,NB )
            CALL DFFGWT( N,NB )
            CALL DFFGAT( N,NB )
          ENDIF
!
!---  NAPL Neumann else Dirichlet, Saturated, Unit Gradient
!
          IF( IBCT(3,NB).EQ.2 ) THEN
            DO 162 M = 1,ISVF
              WN(M,NSZ(N)+IJFLD) = BCX(4)
  162       CONTINUE
            CALL DFFNOT( N,NB )
          ELSEIF( IBCT(3,NB).NE.3 ) THEN
            CALL DRCVNT( N,NB )
            CALL DFFNOT( N,NB )
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCF5 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCJ5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Modify the Jacobian matrix for boundary conditions
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/BCJ5'
!
!---  Loop over boundary conditions  ---
!
      DO 100 NB = 1,NBC
!
!---  Special decreasing pond boundary condition for the MOFHYS
!     verification problem  ---
!
!        IF( NB.EQ.2 .AND. BC(5,1,NB).LE.SMALL ) GOTO 100
        TMZ = TM
        IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 100
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) GOTO 100
        N = IBCN(NB)
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
!---      Aqueous  ---
!
          IF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL JCBLWB( N,NB,NPZ )
            CALL JCBLOB( N,NB,NPZ )
            CALL JCBLAB( N,NB,NPZ )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBGWB( N,NB,NPZ )
            CALL JCBGOB( N,NB,NPZ )
            CALL JCBGAB( N,NB,NPZ )
          ENDIF
!
!---      NAPL  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBNOB( N,NB,NPZ )
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL JCBLWS( N,NB,NPY )
            CALL JCBLOS( N,NB,NPY )
            CALL JCBLAS( N,NB,NPY )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBGWS( N,NB,NPY )
            CALL JCBGOS( N,NB,NPY )
            CALL JCBGAS( N,NB,NPY )
          ENDIF
!
!---      NAPL
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBNOS( N,NB,NPY )
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL JCBLWW( N,NB,NPX )
            CALL JCBLOW( N,NB,NPX )
            CALL JCBLAW( N,NB,NPX )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBGWW( N,NB,NPX )
            CALL JCBGOW( N,NB,NPX )
            CALL JCBGAW( N,NB,NPX )
          ENDIF
!
!---      NAPL  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBNOW( N,NB,NPX )
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL JCBLWE( N,NB,NQX )
            CALL JCBLOE( N,NB,NQX )
            CALL JCBLAE( N,NB,NQX )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBGWE( N,NB,NQX )
            CALL JCBGOE( N,NB,NQX )
            CALL JCBGAE( N,NB,NQX )
          ENDIF
!
!---      NAPL  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBNOE( N,NB,NQX )
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL JCBLWN( N,NB,NQY )
            CALL JCBLON( N,NB,NQY )
            CALL JCBLAN( N,NB,NQY )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBGWN( N,NB,NQY )
            CALL JCBGON( N,NB,NQY )
            CALL JCBGAN( N,NB,NQY )
          ENDIF
!
!---      NAPL  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBNON( N,NB,NQY )
          ENDIF
!
!---    Top boundary  ---
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
!
!---      Aqueous  ---
!
          IF( IBCT(1,NB).NE.3 .AND.
     &      ABS(IBCT(1,NB)).NE.318 ) THEN
            CALL JCBLWT( N,NB,NQZ )
            CALL JCBLOT( N,NB,NQZ )
            CALL JCBLAT( N,NB,NQZ )
          ENDIF
!
!---      Gas  ---
!
          IF( IBCT(2,NB).NE.3 ) THEN
            CALL JCBGWT( N,NB,NQZ )
            CALL JCBGOT( N,NB,NQZ )
            CALL JCBGAT( N,NB,NQZ )
          ENDIF
!
!---      NAPL  ---
!
          IF( IBCT(3,NB).NE.3 ) THEN
            CALL JCBNOT( N,NB,NQZ )
          ENDIF
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJ5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCP5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute saturation, relative permeability and thermodynamic
!     properties for boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     Last Modified by ML Rockhold, PNNL, 29 January 2003.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FDVP
      USE FDVN
      USE FDVG
      USE FDVD
      USE CONST
      USE BCVP
      USE BCVN
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCP5'
      INDX = 1
!
!---  Assign values for initial condition type boundary conditions  ---
!
      IF( NSTEP-NRST.LE.1 .AND. NITER.LE.1 ) THEN
       DO 50 NB = 1,NBC
        IF( IBCT(1,NB).EQ.12 .OR. IBCT(2,NB).EQ.12
     &    .OR. IBCT(3,NB).EQ.12 ) THEN
          N = IBCN(NB)
          IF( IBCD(NB).EQ.-3 ) THEN
             DB = 0.5D+0*DZGF(N)
             NPZ = NSZ(N)
             GB = GRVZ(NPZ)*DB
          ELSEIF( IBCD(NB).EQ.-2 ) THEN
             DB = 0.5D+0*DYGF(N)*RP(ID(N))
             NPY = NSY(N)
             GB = GRVY(NPY)*DB
          ELSEIF( IBCD(NB).EQ.-1 ) THEN
             DB = 0.5D+0*DXGF(N)
             NPX = NSX(N)
             GB = GRVX(NPX)*DB
          ELSEIF( IBCD(NB).EQ.1 ) THEN
             DB = -0.5D+0*DXGF(N)
             NQX = NSX(N)+1
             GB = GRVX(NQX)*DB
          ELSEIF( IBCD(NB).EQ.2 ) THEN
             DB = -0.5D+0*DYGF(N)*RP(ID(N))
             NQY = NSY(N)+IFLD
             GB = GRVY(NQY)*DB
          ELSEIF( IBCD(NB).EQ.3 ) THEN
             DB = -0.5D+0*DZGF(N)
             NQZ = NSZ(N)+IJFLD
             GB = GRVZ(NQZ)*DB
          ENDIF
        ENDIF
        IF( IBCT(1,NB).EQ.12 ) PLB(1,NB) = PL(2,N) + RHOL(2,N)*GB
        IF( IBCT(2,NB).EQ.12 ) PGB(1,NB) = PG(2,N) + RHOG(2,N)*GB
        IF( IBCT(3,NB).EQ.12 ) PNB(1,NB) = PN(2,N) + RHON(2,N)*GB
   50  CONTINUE
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
             IF( IBCT(3,NB).EQ.2 ) THEN
               BCX(4) = BCX(4)-5.D-1*DTBC*(BC(4,M,MB)-BC(4,M-1,MB))/TDBC
             ENDIF
             IF( IBCT(1,NB).EQ.2 ) THEN
               BCX(2) = BCX(2)-5.D-1*DTBC*(BC(2,M,MB)-BC(2,M-1,MB))/TDBC
             ENDIF
             IF( IBCT(2,NB).EQ.2 ) THEN
               BCX(3) = BCX(3)-5.D-1*DTBC*(BC(3,M,MB)-BC(3,M-1,MB))/TDBC
             ENDIF
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 400
        ENDIF
  110   CONTINUE
        N = IBCN(NB)
        IZN = IZ(N)
        IBD = ABS(IBCD(NB))
!
!---    Initial condition boundary condition  ---
!
        IF( IBCT(3,NB).EQ.12 ) BCX(4) = PNB(1,NB)
        IF( IBCT(1,NB).EQ.12 ) BCX(2) = PLB(1,NB)
        IF( IBCT(2,NB).EQ.12 ) BCX(3) = PGB(1,NB)
!
!---  Boundary Direction  ---
!
        I = ID(N)
        J = JD(N)
        K = KD(N)
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX+1
        NQY = NPY+IFLD
        NQZ = NPZ+IJFLD
        IF( IBCD(NB).EQ.-3 ) THEN
           DB = 0.5D+0*DZGF(N)
           GB = GRVZ(NPZ)*DB
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
           DB = 0.5D+0*DYGF(N)*RP(ID(N))
           GB = GRVY(NPY)*DB
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
           DB = 0.5D+0*DXGF(N)
           GB = GRVX(NPX)*DB
        ELSEIF( IBCD(NB).EQ.1 ) THEN
           DB = -0.5D+0*DXGF(N)
           GB = GRVX(NQX)*DB
        ELSEIF( IBCD(NB).EQ.2 ) THEN
           DB = -0.5D+0*DYGF(N)*RP(ID(N))
           GB = GRVY(NQY)*DB
        ELSEIF( IBCD(NB).EQ.3 ) THEN
           DB = -0.5D+0*DZGF(N)
           GB = GRVZ(NQZ)*DB
        ENDIF
!
!---    Fluid-pair scaling factors and entry pressure  ---
!
        BGL = SCHR(8,IZN)/SCHR(5,IZN)
        BNL = SCHR(8,IZN)/(SCHR(6,IZN)*CA_NL)
        BGN = SCHR(8,IZN)/(SCHR(7,IZN)*CA_GN)
        ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---    Loop over secondary variable indices  ---
!
        DO 300 M = 2,ISVC+2
          TX = T(2,N)
          PLX = PL(M,N)
          PGX = PG(M,N)
          PNX = PN(M,N)
!
!---      Aqueous Dirichlet  ---
!
          IF( IBCT(1,NB).EQ.1 .OR. IBCT(1,NB).EQ.12 ) THEN
            PLX = BCX(2)
!
!---      Aqueous Neumann  ---
!
          ELSEIF( IBCT(1,NB).EQ.2 ) THEN
            PLX = PLX + BCX(2)*DB*VISL(M,N)/PERM(IBD,IZN)
     &          + RHOL(M,N)*GB
!
!---      Aqueous Zero Flux  ---
!
          ELSEIF( IBCT(1,NB).EQ.3 ) THEN
            IF( ABS(BCX(2)+PATM).GT.EPSL ) THEN
              PLX = BCX(2)
            ELSE
              PLX = PLX + RHOL(M,N)*GB
            ENDIF
!
!---      Aqueous Unit Gradient  ---
!
          ELSEIF( IBCT(1,NB).EQ.5 ) THEN
            PLX = PLX
!
!---      Aqueous Hydraulic Gradient Base Node  ---
!
          ELSEIF( IBCT(1,NB).EQ.-11 ) THEN
            PLX = BCX(2)
!
!---      Aqueous Hydraulic Gradient  ---
!
          ELSEIF( IBCT(1,NB).EQ.11 ) THEN
            IF( M.EQ.2 ) THEN
              PLX = FNHGBL( NBHG(1,NB),NB,M )
            ELSE
              PLX = PLB(2,NB)
            ENDIF
!
!---      Aqueous Fluctuating Water Table Base Node  ---
!
          ELSEIF( IBCT(1,NB).EQ.-18 ) THEN
            PLX = BCX(2)
!
!---      Aqueous Fluctuating Water Table  ---
!
          ELSEIF( IBCT(1,NB).EQ.18 ) THEN
            IF( M.EQ.2 ) THEN
              PLX = FNHGBL( NBHG(1,NB),NB,M )
            ELSE
              PLX = PLB(2,NB)
            ENDIF
          ENDIF
!
!---      Gas Dirichlet or Dirichlet-Outflow  ---
!
          IF( IBCT(2,NB).EQ.1 .OR. IBCT(2,NB).EQ.19
     &      .OR. IBCT(2,NB).EQ.12 ) THEN
            PGX = BCX(3)
!
!---      Gas Neumann  ---
!
          ELSEIF( IBCT(2,NB).EQ.2 ) THEN
            PGX = PGX + BCX(3)*DB*VISG(M,N)/PERM(IBD,IZN)
     &          + RHOG(M,N)*GB
!
!---      Gas Zero Flux  ---
!
          ELSEIF( IBCT(2,NB).EQ.3 ) THEN
            IF( ABS(BCX(3)+PATM).GT.EPSL ) THEN
              PGX = BCX(3)
            ELSE
              PGX = PGX + RHOG(M,N)*GB
            ENDIF
!
!---      Gas Unit Gradient  ---
!
          ELSEIF( IBCT(2,NB).EQ.5 ) THEN
            PGX = PGX
!
!---      Gas Hydraulic Gradient Base Node  ---
!
          ELSEIF( IBCT(2,NB).EQ.-11 ) THEN
            PGX = BCX(3)
!
!---      Gas Hydraulic Gradient  ---
!
          ELSEIF( IBCT(2,NB).EQ.11 ) THEN
            IF( M.EQ.2 ) THEN
              PGX = FNHGBG( NBHG(2,NB),NB,M )
            ELSE
              PGX = PGB(2,NB)
            ENDIF
          ENDIF
!
!---      NAPL Dirichlet  ---
!
          IF( IBCT(3,NB).EQ.1 .OR. IBCT(3,NB).EQ.12 ) THEN
            PNX = BCX(4)
!
!---      NAPL Neumann; Logic prevents NAPL flow 
!         from two-phase nodes  ---
!
          ELSEIF( IBCT(3,NB).EQ.2 ) THEN
            IF( NPHAZ(2,N).LE.2 .AND. BCX(4)*DB.LT.ZERO ) THEN
              INDX = 12
              CHMSG = 'NAPL Boundary Flux: Boundary Node: '
              IMSG = N
              CALL WRMSGS( INDX )
            ELSE
              PNX = PNX + BCX(4)*DB*VISN(M,N)/PERM(IBD,IZN)
     &          + RHON(M,N)*GB
            ENDIF
!
!---      NAPL Zero Flux  ---
!
          ELSEIF( IBCT(3,NB).EQ.3 ) THEN
            IF( ABS(BCX(4)+PATM).GT.EPSL ) THEN
              PNX = BCX(4)
            ELSE
              PNX = PNX + RHON(M,N)*GB
            ENDIF
!
!---      NAPL Unit Gradient  ---
!
          ELSEIF( IBCT(3,NB).EQ.5 ) THEN
            PNX = PNX
!
!---      NAPL Hydraulic Gradient Base Node  ---
!
          ELSEIF( IBCT(3,NB).EQ.-11 ) THEN
            PNX = BCX(4)
!
!---      NAPL Hydraulic Gradient  ---
!
          ELSEIF( IBCT(3,NB).EQ.11 ) THEN
            IF( M.EQ.2 ) THEN
              PNX = FNHGBN( NBHG(3,NB),NB,M )
            ELSE
              PNX = PNB(2,NB)
            ENDIF
          ENDIF
!
!---      NAPL critical pressure  ---
!
          PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),PLX+ENPR/BNL )
!
!----     Check aqueous pressure and switch between
!         hydraulic gradient (for saturated nodes) and zero flux
!         (for unsaturated nodes)  ---
!
          IF( MOD( ABS(IBCT(1,NB)),100 ).EQ.18 ) THEN
            HGN = MAX(ZERO,BGN*(PGX-MAX(PNC,PNX)))
            IF( HGN.GT.ENPR ) THEN
              PLX = PL(M,N) + RHOL(M,N)*GB
              PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),PLX+ENPR/BNL )
              IF( IBCT(1,NB).LT.0 ) THEN
                IBCT(1,NB) = -318
              ELSE
                IBCT(1,NB) = 318
              ENDIF
            ELSE
              IF( IBCT(1,NB).LT.0 ) THEN
                IBCT(1,NB) = -18
              ELSE
                IBCT(1,NB) = 18
              ENDIF
            ENDIF
          ENDIF
          PNX = MAX( PNC,PNX )
!
!---      Compute phase saturation and relative permeability
!         assuming no entrapment on boundary surfaces  ---
!
          INDX = 2
          TRPNLX = SCHR(9,IZN)/(SCHR(9,IZN)+TRPNL(2,N))
          ASNTMX = SCHR(11,IZN)*TRPNLX
          ASNRMX = SCHR(10,IZN)
          CALL KSP5( PGX,PNX,PLX,SLB(M,NB),SNB(M,NB),SGB(M,NB),
     &      RKLB(1,M,NB),RKNB(M,NB),RKGB(M,NB),ASTX,ASLX,
     &      ASLMINX,ASTMINX,ASTMAXX,ASNTX,ASNRX,ASGTX,ASGTLX,
     &      ASGTNX,ASNTMX,ASNRMX,ASGTMX,INDX,IZN,NB )
!
!---  Compute porosity and tortuosity  ---
!
          PX = MAX( PGX,PLX,PNX )+PATM
          CALL PORSTY( N,PX,PCMP(N),PORDB(M,NB),PORTB(M,NB) )
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SLB(M,NB),SGB(M,NB),
     &      SNB(M,NB),PORDB(M,NB),TORLB(M,NB),TORGB(M,NB),TORNB(M,NB) )
!
!---  Convert pressure to absolute prior to computing physical
!     properties  ---
!
          PLX = PLX + PATM
          PGX = PGX + PATM
          PNX = PNX + PATM
!
!---  Compute the liquid water density and viscosity  ---
!
          PX = MAX( PLX,PGX,PNX )
          CALL WATLQD( TX,PX,RHOLB(M,NB) )
          CALL WATLQV( TX,PX,PSWB(2,NB),VISLB(M,NB) )
!
!---  Compute the saturated water pressure, water vapor partial
!     pressure, saturated oil pressure, oil vapor partial pressure,
!     air partial pressure  ---
!
          CALL WATSP( TX,PSWB(2,NB) )
          PVWB(M,NB) = PSWB(2,NB)*BCX(7)
          CALL VOCSP( TX,PSOB(2,NB) )
          PVOB(M,NB) = PSOB(2,NB)*BCX(8)
          PVAB(M,NB) = MAX( PGX-PVWB(2,NB)-PVOB(2,NB),ZERO )
!
!---  Compute the gas-phase component densities, density,
!     component mass fractions, and component mole fractions  ---
!
          INDX = 0
          CALL WATGSD( TX,PVWB(M,NB),RHOW,INDX )
          CALL VOCGSD( TX,PVOB(M,NB),RHOO,INDX )
          CALL AIRGSD( TX,PVAB(M,NB),RHOA )
          RHOGB(M,NB) = RHOW+RHOO+RHOA
          XGWB(M,NB) = RHOW/RHOGB(M,NB)
          XGOB(M,NB) = RHOO/RHOGB(M,NB)
          XGAB(M,NB) = RHOA/RHOGB(M,NB)
          XMGWB(M,NB) = PVWB(M,NB)/PGX
          XMGAB(M,NB) = PVAB(M,NB)/PGX
          XMGOB(M,NB) = PVOB(M,NB)/PGX
          WTMX = XMGWB(M,NB)*WTMW + XMGAB(M,NB)*WTMA + XMGOB(M,NB)*WTMO
          RHOMGB(M,NB) = RHOGB(M,NB)/WTMX
!
!---  Compute the aqueous-phase mole and mass fractions  ---
!
          XMLAB(M,NB) = BCX(5)*PVAB(M,NB)/HCAW
          XMLOB(M,NB) = BCX(6)*PVOB(M,NB)/HCOW
          XMLWB(M,NB) = MAX( 1.D+0-XMLAB(M,NB)-XMLOB(M,NB),0.D+0 )
          WTMX = XMLWB(M,NB)*WTMW + XMLAB(M,NB)*WTMA + XMLOB(M,NB)*WTMO
          XLWB(M,NB) = XMLWB(M,NB)*WTMW/WTMX
          XLAB(M,NB) = XMLAB(M,NB)*WTMA/WTMX
          XLOB(M,NB) = XMLOB(M,NB)*WTMO/WTMX
          RHOMLB(M,NB) = RHOLB(M,NB)/WTMX
!
!---  Compute the NAPL mass fractions  ---
!
          XNOB(M,NB) = 1.D+0
!
!---  Compute the NAPL density and viscosity  ---
!
          CALL VOCLQD( TX,PX,PVOB(M,NB),RHONB(M,NB) )
          CALL VOCLQV( TX,VISNB(M,NB),VISCO,IVISO )
!
!---  Compute the gas component viscosity and phase viscosity  ---
!
          CALL WATGSV( TX,VISW )
          CALL AIRGSV( TX,VISA )
          CALL VOCGSV( TX,VISO )
          CALL GASVIS( XMGWB(M,NB),XMGOB(M,NB),XMGAB(M,NB),
     &      VISW,VISO,VISA,VISGB(M,NB) )
!
!---  Compute vapor diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFGOB(M,NB) = DFGOC
            DFGWB(M,NB) = DFGWC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL BNDFAW( TX,PGX,DFGAW )
            CALL BNDFAO( TX,PGX,DFGAO )
            CALL BNDFOW( TX,PGX,DFGOW )
            CALL VOCDFG( XMGWB(M,NB),XMGOB(M,NB),XMGAB(M,NB),
     &        DFGAO,DFGOW,DFGOB(M,NB) )
            CALL WATDFG( XMGWB(M,NB),XMGOB(M,NB),XMGAB(M,NB),
     &        DFGAW,DFGOW,DFGWB(M,NB) )
          ENDIF
!
!---  Compute liquid diffusion coefficients  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLOB(M,NB) = DFLOC
            DFLAB(M,NB) = DFLAC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            CALL AIRDFL( TX,VISLB(M,NB),DFLAB(M,NB) )
            CALL VOCDFL( TX,VISLB(M,NB),DFLOB(M,NB) )
          ENDIF
!
!---  Convert boundary pressure to gauge pressure  ---
!
          TB(M,NB) = TX
          PLB(M,NB) = PLX - PATM
          PNB(M,NB) = PNX - PATM
          PGB(M,NB) = PGX - PATM
  300   CONTINUE
  400 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCP5 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CAP5( ASLMINX,PGX,PLX,PNX,ASNTMX,ASNX,INDX,IZN,NX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!
!     Compute the NAPL pressure from the NAPL saturation and
!     aqueous pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1999.
!     Last Modified by MD White on January 11, 1999.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE FILES
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
      SUB_LOG(ISUB_LOG) = '/CAP5'
!
!---  Entrapment van Genuchten saturation function  ---
!
      IF( ISCHR(IZN).EQ.32 ) THEN
        CN = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( IRPL(IZN).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(14,IZN)
        ENDIF
        ENPR = SCHR(2,IZN)*RHORL*GRAV
        PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
        IF( ASNX.LT.1.D-6 ) THEN
          PNX = PNC
          GOTO 152
        ENDIF
        PLND = 1.D+0/ASNTMX - 1.D+0
        TNFZ = MAX( 1.D-3,MIN( (SNT(2,NX)/(SN(2,NX)+SMALL)),0.999D+0 ) )
        ASNTZ = ASNX*TNFZ
        ASNFZ = ASNX*(1.D+0-TNFZ)
        SNRZ = (1.D+0-ASLMINX)/(1.D+0+PLND*(1.D+0-ASLMINX))
        ASLZ = (1.D+0+(1.D+0+PLND)*(ASNTZ-SNRZ))/
     &    (1.D+0+PLND*(ASNTZ-SNRZ))
        ASLZ = MIN( ASLZ,1.D+0-ASNFZ )
        IF( (1.D+0-ASLZ).GT.1.D-9 ) THEN
          VAR = (1.D+0/ASLZ)
          IF( VAR.LE.0.D+0 ) THEN
            PRINT *,'ASLZ: ',ASLZ
            PRINT *,'ASNX: ',ASNX
            PRINT *,'TNFZ: ',TNFZ
            PRINT *,'ASNTZ: ',ASNTZ
            PRINT *,'ASNFZ: ',ASNFZ
            PRINT *,'ASLMINX: ',ASLMINX
            PRINT *,'PLND: ',PLND
            PRINT *,'SNRZ: ',SNRZ
            PRINT *,'SNT: ',SNT(2,NX)
            PRINT *,'SN2: ',SN(2,NX)
            PRINT *,'SN1: ',SN(1,NX)
            PRINT *,'NPHAZ2: ',NPHAZ(2,NX)
            PRINT *,'NPHAZ1: ',NPHAZ(1,NX)
            STOP
          ENDIF
          HNLZ = ((-1.D+0 + VAR**(1.D+0/CM))**(1.D+0/CN))/SCHR(1,IZN)
          PNZ = PLX + HNLZ*GRAV*RHORL/BNL
        ELSE
          PNZ = PLX
        ENDIF
        PNZI = PNZ
        NC = 0
  150   CONTINUE
        HNLZ = MAX( ZERO,BNL*(PNZ-PLX)/RHORL/GRAV )
        HGNZ = MAX( ZERO,BGN*(PGX-PNZ)/RHORL/GRAV )
        ASTZ = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HGNZ)**CN))**CM
        ASLZ = MIN(ASTZ,(1.D+0/(1.D+0 + (SCHR(1,IZN)*HNLZ)**CN))**CM)
        ASLMINZ = MIN( ASLZ,ASLMINX )
        SNRZ = (1.D+0-ASLMINZ)/(1.D+0+PLND*(1.D+0-ASLMINZ))
     &    - (1.D+0-ASLZ)/(1.D+0+PLND*(1.D+0-ASLZ))
        ASNTZ = MIN( ASNX,SNRZ )
        F = ASNX - ASTZ + ASLZ - ASNTZ
        DASTZ = CM*CN*SCHR(1,IZN)*BGN*((SCHR(1,IZN)*HGNZ)**(CN-1.D+0))/
     &    ((1.D+0 + (SCHR(1,IZN)*HGNZ)**CN)**(CM+1.D+0))/RHORL/GRAV
        DASLZ = -CM*CN*SCHR(1,IZN)*BNL*((SCHR(1,IZN)*HNLZ)**(CN-1.D+0))/
     &    ((1.D+0 + (SCHR(1,IZN)*HNLZ)**CN)**(CM+1.D+0))/RHORL/GRAV
        IF( ASLZ.GT.ASLMINZ ) THEN
          DASNTZ = DASLZ*((1.D+0/(1.D+0+(1.D+0-ASLZ)*PLND)) -
     &      ((1.D+0-ASLZ)*PLND)/((1.D+0+(1.D+0-ASLZ)*PLND)**2))
        ELSE
          DASNTZ = 0.D+0
        ENDIF
        DF = DASLZ - DASTZ - DASNTZ
        IF( ABS(DF).LT.EPSL ) GOTO 151
        DPNZ = -F/DF
        IF( PNZ+DPNZ.LT.PNC .AND. NC.LT.3 ) THEN
          PNZ = PNC + 1.D+2/(1.D+1**NC)
        ELSE
          PNZ = MAX( PNZ+DPNZ,PNC )
        ENDIF
        NC = NC + 1
        IF( NC.GT.32 ) THEN
          ICNV = 1
          WRITE(ISC,'(A)') 'Execution Note: CAP5 Convergence Failure'
          WRITE(IWR,'(A)') 'Execution Note: CAP5 Convergence Failure'
          WRITE(ISC,'(A,I6,$)') '  Node: ',NX
          WRITE(IWR,'(A,I6,$)') '  Node: ',NX
          WRITE(ISC,'(A,I1,$)') '  Phase: ',NPHAZ(2,NX)
          WRITE(IWR,'(A,I1,$)') '  Phase: ',NPHAZ(2,NX)
          WRITE(ISC,'(A,I2,$)') '  Iterations: ',NC
          WRITE(IWR,'(A,I2,$)') '  Iterations: ',NC
          WRITE(ISC,'(A,1PE12.5,$)') '  ASNX: ',ASNX
          WRITE(IWR,'(A,1PE12.5,$)') '  ASNX: ',ASNX
          WRITE(ISC,'(A,1PE12.5,$)') '  PLX: ',PLX
          WRITE(IWR,'(A,1PE12.5,$)') '  PLX: ',PLX
          WRITE(ISC,'(A,1PE12.5)') '  PNZI: ',PNZI
          WRITE(IWR,'(A,1PE12.5)') '  PNZI: ',PNZI
          WRITE(ISC,'(A,1PE12.5,$)') '  PNZ: ',PNZ
          WRITE(ISC,'(A,1PE12.5,$)') '  BGN: ',BGN
          WRITE(ISC,'(A,1PE12.5,$)') '  BNL: ',BNL
          WRITE(ISC,'(A,1PE12.5,$)') '  BGL: ',BGL
          WRITE(ISC,'(A,1PE12.5,$)') '  CN: ',CN
          WRITE(ISC,'(A,1PE12.5)') '  CM: ',CM
          WRITE(ISC,'(A,1PE12.5,$)') '  ALPHA: ',SCHR(1,IZN)
          WRITE(ISC,'(A,1PE12.5,$)') '  ASLMINX: ',ASLMINX
          WRITE(ISC,'(A,1PE12.5)') '  PLND: ',PLND
          CALL UPDT5
          INDX = -INDX
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
        IF( ABS(DPNZ).GT.1.D-6 ) GOTO 150
  151   CONTINUE
        IF( NC.LT.1 .AND. ASNX.GT.1.D-6 ) THEN
          ICNV = 1
          WRITE(ISC,'(A)') 'Execution Note: CAP5 Convergence Failure'
          WRITE(IWR,'(A)') 'Execution Note: CAP5 Convergence Failure'
          WRITE(ISC,'(A,I6,$)') '  Node: ',NX
          WRITE(IWR,'(A,I6,$)') '  Node: ',NX
          WRITE(ISC,'(A,I1,$)') '  Phase: ',NPHAZ(2,NX)
          WRITE(IWR,'(A,I1,$)') '  Phase: ',NPHAZ(2,NX)
          WRITE(ISC,'(A,I2,$)') '  Iterations: ',NC
          WRITE(IWR,'(A,I2,$)') '  Iterations: ',NC
          WRITE(ISC,'(A,1PE12.5,$)') '  ASNX: ',ASNX
          WRITE(IWR,'(A,1PE12.5,$)') '  ASNX: ',ASNX
          WRITE(ISC,'(A,1PE12.5,$)') '  PLX: ',PLX
          WRITE(IWR,'(A,1PE12.5,$)') '  PLX: ',PLX
          WRITE(ISC,'(A,1PE12.5)') '  PNZI: ',PNZI
          WRITE(IWR,'(A,1PE12.5)') '  PNZI: ',PNZI
          CALL UPDT5
          INDX = -INDX
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
        PNX = MAX( PNC,PNZ )
  152   CONTINUE
!
!---  Entrapment Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.36 ) THEN
        CL = MAX( SCHR(3,IZN),SMALL )
!        SMX = SCHR(4,IZN)
!        ASNX = SNX/(1.D+0-SMX)
        ENPR = SCHR(2,IZN)*RHORL*GRAV
        PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
        IF( ASNX.LT.1.D-6 ) THEN
          PNX = PNC
          GOTO 162
        ENDIF
        PLND = 1.D+0/ASNTMX - 1.D+0
        TNFZ = MAX( 1.D-3,MIN( (SNT(2,NX)/(SN(2,NX)+SMALL)),0.999D+0 ) )
        ASNTZ = ASNX*TNFZ
        ASNFZ = ASNX*(1.D+0-TNFZ)
        SNRZ = (1.D+0-ASLMINX)/(1.D+0+PLND*(1.D+0-ASLMINX))
        ASLZ = (1.D+0+(1.D+0+PLND)*(ASNTZ-SNRZ))/
     &    (1.D+0+PLND*(ASNTZ-SNRZ))
        ASLZ = MIN( ASLZ,1.D+0-ASNFZ )
        IF( (1.D+0-ASLZ).GT.1.D-9 ) THEN
          VAR = (1.D+0/ASLZ)
          IF( VAR.LE.0.D+0 ) THEN
            PRINT *,'ASLZ: ',ASLZ
            PRINT *,'ASNX: ',ASNX
            PRINT *,'TNFZ: ',TNFZ
            PRINT *,'ASNTZ: ',ASNTZ
            PRINT *,'ASNFZ: ',ASNFZ
            PRINT *,'ASLMINX: ',ASLMINX
            PRINT *,'PLND: ',PLND
            PRINT *,'SNRZ: ',SNRZ
            PRINT *,'SNT: ',SNT(2,NX)
            PRINT *,'SN2: ',SN(2,NX)
            PRINT *,'SN1: ',SN(1,NX)
            PRINT *,'NPHAZ2: ',NPHAZ(2,NX)
            PRINT *,'NPHAZ1: ',NPHAZ(1,NX)
            STOP
          ENDIF
          HNLZ = SCHR(2,IZN)*VAR**(1./CL)
          PNZ = PLX + HNLZ*GRAV*RHORL/BNL
        ELSE
          PNZ = PLX
        ENDIF
        PNZI = PNZ
        NC = 0
  160   CONTINUE
        HNLZ = MAX( ZERO,BNL*(PNZ-PLX)/RHORL/GRAV )
        HGNZ = MAX( ZERO,BGN*(PGX-PNZ)/RHORL/GRAV )
        IF( HGNZ-SCHR(2,IZN).LE.EPSL ) THEN
          ASTZ = 1.D+0
        ELSE
          ASTZ = (SCHR(2,IZN)/HGNZ)**CL
        ENDIF
        IF( HNLZ-SCHR(2,IZN).LE.EPSL ) THEN
          ASLZ = 1.D+0
        ELSE
          ASLZ = (SCHR(2,IZN)/HNLZ)**CL
        ENDIF
        ASLZ = MIN( ASTZ,ASLZ )
        ASLMINZ = MIN( ASLZ,ASLMINX )
        SNRZ = (1.D+0-ASLMINZ)/(1.D+0+PLND*(1.D+0-ASLMINZ))
     &    - (1.D+0-ASLZ)/(1.D+0+PLND*(1.D+0-ASLZ))
        ASNTZ = MIN( ASNX,SNRZ )
        F = ASNX - ASTZ + ASLZ - ASNTZ
        DASTZ = -CL*ASTZ*BGN/HGNZ/RHORL/GRAV
        DASLZ = -CL*ASLZ*BNL/HNLZ/RHORL/GRAV
        IF( ASLZ.GT.ASLMINZ ) THEN
          DASNTZ = DASLZ*((1.D+0/(1.D+0+(1.D+0-ASLZ)*PLND)) -
     &      ((1.D+0-ASLZ)*PLND)/((1.D+0+(1.D+0-ASLZ)*PLND)**2))
        ELSE
          DASNTZ = 0.D+0
        ENDIF
        DF = DASLZ - DASTZ - DASNTZ
        IF( ABS(DF).LT.EPSL ) GOTO 161
        DPNZ = -F/DF
        IF( PNZ+DPNZ.LT.PNC .AND. NC.LT.3 ) THEN
          PNZ = PNC + 1.D+2/(1.D+1**NC)
        ELSE
          PNZ = MAX( PNZ+DPNZ,PNC )
        ENDIF
        NC = NC + 1
        IF( NC.GT.32 ) THEN
          ICNV = 1
          WRITE(ISC,'(A)') 'Execution Note: CAP5 Convergence Failure'
          WRITE(IWR,'(A)') 'Execution Note: CAP5 Convergence Failure'
          WRITE(ISC,'(A,I6,$)') '  Node: ',NX
          WRITE(IWR,'(A,I6,$)') '  Node: ',NX
          WRITE(ISC,'(A,I1,$)') '  Phase: ',NPHAZ(2,NX)
          WRITE(IWR,'(A,I1,$)') '  Phase: ',NPHAZ(2,NX)
          WRITE(ISC,'(A,I2,$)') '  Iterations: ',NC
          WRITE(IWR,'(A,I2,$)') '  Iterations: ',NC
          WRITE(ISC,'(A,1PE12.5,$)') '  ASNX: ',ASNX
          WRITE(IWR,'(A,1PE12.5,$)') '  ASNX: ',ASNX
          WRITE(ISC,'(A,1PE12.5,$)') '  PLX: ',PLX
          WRITE(IWR,'(A,1PE12.5,$)') '  PLX: ',PLX
          WRITE(ISC,'(A,1PE12.5)') '  PNZI: ',PNZI
          WRITE(IWR,'(A,1PE12.5)') '  PNZI: ',PNZI
          WRITE(ISC,'(A,1PE12.5,$)') '  PNZ: ',PNZ
          WRITE(ISC,'(A,1PE12.5,$)') '  BGN: ',BGN
          WRITE(ISC,'(A,1PE12.5,$)') '  BNL: ',BNL
          WRITE(ISC,'(A,1PE12.5,$)') '  BGL: ',BGL
          WRITE(ISC,'(A,1PE12.5,$)') '  CN: ',CN
          WRITE(ISC,'(A,1PE12.5)') '  CM: ',CM
          WRITE(ISC,'(A,1PE12.5,$)') '  ALPHA: ',SCHR(1,IZN)
          WRITE(ISC,'(A,1PE12.5,$)') '  ASLMINX: ',ASLMINX
          WRITE(ISC,'(A,1PE12.5)') '  PLND: ',PLND
          CALL UPDT5
          INDX = -INDX
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
        IF( ABS(DPNZ).GT.1.D-6 ) GOTO 160
  161   CONTINUE
        IF( NC.LT.1 .AND. ASNX.GT.1.D-6 ) THEN
          ICNV = 1
          WRITE(ISC,'(A)') 'Execution Note: CAP5 Convergence Failure'
          WRITE(IWR,'(A)') 'Execution Note: CAP5 Convergence Failure'
          WRITE(ISC,'(A,I6,$)') '  Node: ',NX
          WRITE(IWR,'(A,I6,$)') '  Node: ',NX
          WRITE(ISC,'(A,I1,$)') '  Phase: ',NPHAZ(2,NX)
          WRITE(IWR,'(A,I1,$)') '  Phase: ',NPHAZ(2,NX)
          WRITE(ISC,'(A,I2,$)') '  Iterations: ',NC
          WRITE(IWR,'(A,I2,$)') '  Iterations: ',NC
          WRITE(ISC,'(A,1PE12.5,$)') '  ASNX: ',ASNX
          WRITE(IWR,'(A,1PE12.5,$)') '  ASNX: ',ASNX
          WRITE(ISC,'(A,1PE12.5,$)') '  PLX: ',PLX
          WRITE(IWR,'(A,1PE12.5,$)') '  PLX: ',PLX
          WRITE(ISC,'(A,1PE12.5)') '  PNZI: ',PNZI
          WRITE(IWR,'(A,1PE12.5)') '  PNZI: ',PNZI
          CALL UPDT5
          INDX = -INDX
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
        PNX = MAX( PNC,PNZ )
  162   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CAP5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHK5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Check the thermodynamic and hydrologic states declared through
!     user inputs.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, October 15, 1997.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE FDVP
      USE FDVD
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
      SUB_LOG(ISUB_LOG) = '/CHK5'
!
!---  Pressure tolerance (0.1 Pa)  ---
!
      PTOLX = 0.1
!
!---  Pure liquid water reference density  ---
!
      TX = 20.D+0
      PX = 101325.0
      CALL WATLQD( TX,PX,RHORL )
!
!---  Relative permeability function compatibility  ---
!
      DO 10 N = 1,NROCK
        IF( ISCHR(N).EQ.13 .OR. ISCHR(N).EQ.15 .OR.
     &    ISCHR(N).EQ.17 ) THEN
          IF( SCHR(14,N).LE.ZERO ) THEN
            IF( IRPL(N).EQ.2 ) THEN
              SCHR(14,N) = 1.D+0 - 2.D+0/SCHR(3,N)
            ELSE
              SCHR(14,N) = 1.D+0 - 1.D+0/SCHR(3,N)
            ENDIF
          ENDIF
          IF( SCHR(14,N).LE.ZERO ) THEN
            NCH = INDEX( ROCK(N),'  ' )-1
            CHMSG = 'Negative or Zero van Genuchten m Parameter: ' //
     &        ROCK(N)(1:NCH)
            INDX = 3
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
   10 CONTINUE
!
!---  Compute saturated vapor pressures and initialize apparent
!     minimum saturations for restart simulations  ---
!
      DO 50 N = 1,NFLD
        CALL WATSP( T(2,N),PSW(2,N) )
        CALL VOCSP( T(2,N),PSO(2,N) )
        ASLMIN(1,N) = ASLMIN(2,N)
        ASTMAX(1,N) = ASTMAX(2,N)
        ASTMIN(1,N) = ASTMIN(2,N)
   50 CONTINUE
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 100
        DO 90 M = 1,ISVC+2
          T(M,N) = T(2,N)
          PG(M,N) = PG(2,N)
          PL(M,N) = PL(2,N)
          PN(M,N) = PN(2,N)
          SG(M,N) = SG(2,N)
          SL(M,N) = SL(2,N)
          SN(M,N) = SN(2,N)
   90   CONTINUE
!
!---  Check initial pressures and temperature  ---
!
        IF( T(2,N).GT.TMX .OR. T(2,N).LT.TMN ) THEN
          CHMSG = 'Temperature, C: Node: '
          INDX = 16
          IMSG = N
          RLMSG = T(2,N)
          CALL WRMSGS( INDX )
        ENDIF
        IF( PL(2,N).GT.PMX-PATM ) THEN
          CHMSG = 'Aqueous Pressure, Pa: Node: '
          INDX = 16
          IMSG = N
          RLMSG = PL(2,N)+PATM
          CALL WRMSGS( INDX )
        ENDIF
        IF( PN(2,N).GT.PMX-PATM ) THEN
          CHMSG = 'NAPL Pressure, Pa: Node: '
          INDX = 16
          IMSG = N
          RLMSG = PN(2,N)+PATM
          CALL WRMSGS( INDX )
       ENDIF
        IF( PG(2,N).GT.PMX-PATM .OR. PG(2,N)+PATM.LT.0.D+0 ) THEN
          CHMSG = 'Gas Pressure, Pa: Node: '
          INDX = 16
          IMSG = N
          RLMSG = PG(2,N)+PATM
          CALL WRMSGS( INDX )
        ENDIF
  100 CONTINUE
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 110 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 110
        PX = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
        PX = MAX( PX,PSW(2,N)+PSO(2,N) )
        CALL WATLQD( T(2,N),PX,RHOLWX )
!
!---    Restart  ---
!
        IF( IEO.EQ.2 ) THEN
          XLA(2,N) = -XMLA(2,N)
          XLO(2,N) = -XMLO(2,N)
        ENDIF
!
!---    Dissolved-air initial condition provided as concentration  ---
!
        IF( XLA(2,N).LT.-EPSL .AND. XMLA(2,N).LT.-EPSL ) THEN
          RHOLAX = -XLA(2,N)
!
!---      Dissolved-oil initial condition provided as concentration  ---
!
          IF( XLO(2,N).LT.-EPSL .AND. XMLO(2,N).LT.-EPSL ) THEN
            RHOLOX = -XLO(2,N)
            RHOL(2,N) = RHOLWX + RHOLOX + RHOLAX
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Dissolved-oil initial condition provided as mass fraction  ---
!
          ELSEIF( XLO(2,N).GT.EPSL .AND. XMLO(2,N).LT.-EPSL ) THEN
            RHOL(2,N) = (RHOLAX+RHOLWX)/(1.D+0-XLO(2,N))
            RHOLOX = RHOL(2,N) - RHOLAX - RHOLWX
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Dissolved-oil initial condition provided as mole fraction  ---
!
          ELSEIF( XLO(2,N).LT.-EPSL .AND. XMLO(2,N).GT.EPSL ) THEN
            RHOLOX = XMLO(2,N)*(RHOLAX*WTMO*WTMW + RHOLWX*WTMO*WTMA)/
     &        (WTMA*WTMW*(1.D+0-XMLO(2,N)))
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Zero dissolved-oil initial condition  ---
!
          ELSE
            RHOLOX = 0.D+0
            RHOL(2,N) = RHOLWX + RHOLOX + RHOLAX
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
          ENDIF
!
!---    Dissolved-air initial condition provided as mass fraction  ---
!
        ELSEIF( XLA(2,N).GT.EPSL .AND. XMLA(2,N).LT.-EPSL ) THEN
!
!---      Dissolved-oil initial condition provided as concentration  ---
!
          IF( XLO(2,N).LT.-EPSL .AND. XMLO(2,N).LT.-EPSL ) THEN
            RHOLOX = -XLO(2,N)
            RHOLAX = XLA(2,N)*(RHOLOX+RHOLAX)/(1.D+0-XLA(2,N))
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Dissolved-oil initial condition provided as mass fraction  ---
!
          ELSEIF( XLO(2,N).GT.EPSL .AND. XMLO(2,N).LT.-EPSL ) THEN
            RHOL(2,N) = RHOLWX/(1.D+0-XLO(2,N)-XLA(2,N))
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Dissolved-oil initial condition provided as mole fraction  ---
!
          ELSEIF( XLO(2,N).LT.-EPSL .AND. XMLO(2,N).GT.EPSL ) THEN
            RHOL(2,N) = 
     &        (RHOLWX*WTMA*(XMLO(2,N)*WTMW-XMLO(2,N)*WTMO-WTMW))/
     &        (WTMW*(XMLO(2,N)*WTMA-XMLO(2,N)*WTMA*XLA(2,N)+
     &        XMLO(2,N)*XLA(2,N)*WTMO-WTMA+WTMA*XLA(2,N)))
            RHOLAX = XLA(2,N)*RHOL(2,N)
            RHOLOX = RHOL(2,N)-RHOLAX-RHOLWX
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Zero dissolved-oil initial condition  ---
!
          ELSE
            RHOLOX = 0.D+0
            RHOL(2,N) = RHOLWX/(1.D+0-XLA(2,N))
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
          ENDIF
!
!---    Dissolved-air initial condition provided as mole fraction  ---
!
        ELSEIF( XLA(2,N).LT.-EPSL .AND. XMLA(2,N).GT.EPSL ) THEN
!
!---      Dissolved-oil initial condition provided as concentration  ---
!
          IF( XLO(2,N).LT.-EPSL .AND. XMLO(2,N).LT.-EPSL ) THEN
            RHOLOX = -XLO(2,N)
            RHOLAX = (XMLA(2,N)*RHOLOX*WTMA*WTMW +
     &        XMLA(2,N)*RHOLWX*WTMO*WTMA)/((1.D+0-XMLA(2,N))*WTMO*WTMW)
            RHOL(2,N) = RHOLAX+RHOLOX+RHOLWX
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Dissolved-oil initial condition provided as mass fraction  ---
!
          ELSEIF( XLO(2,N).GT.EPSL .AND. XMLO(2,N).LT.-EPSL ) THEN
            RHOL(2,N) = 
     &        (RHOLWX*WTMO*(XMLA(2,N)*WTMW-XMLA(2,N)*WTMA-WTMW))/
     &        (WTMW*(XMLA(2,N)*WTMO-XMLA(2,N)*WTMO*XLO(2,N)+
     &        XMLA(2,N)*XLO(2,N)*WTMA-WTMO+WTMO*XLO(2,N)))
            RHOLOX = XLO(2,N)*RHOL(2,N)
            RHOLAX = RHOL(2,N)-RHOLOX-RHOLWX
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Dissolved-oil initial condition provided as mole fraction  ---
!
          ELSEIF( XLO(2,N).LT.-EPSL .AND. XMLO(2,N).GT.EPSL ) THEN
            VARX = RHOLWX*(XMLO(2,N)+XMLA(2,N))/
     &        (WTMW*(1.D+0-XMLO(2,N)-XMLA(2,N)))
            RHOLOX = XMLO(2,N)*WTMO*(VARX+(RHOLWX/WTMW))
            RHOLAX = XMLA(2,N)*WTMA*(VARX+(RHOLWX/WTMW))
            RHOL(2,N) = RHOLAX+RHOLOX+RHOLWX
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Zero dissolved-oil initial condition  ---
!
          ELSE
            RHOLOX = 0.D+0
            RHOLAX = (XMLA(2,N)*RHOLWX*WTMA)/((1.D+0-XMLA(2,N))*WTMW)
            RHOL(2,N) = RHOLAX+RHOLOX+RHOLWX
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
          ENDIF
        ELSE
!
!---      Dissolved-oil initial condition provided as concentration  ---
!
          IF( XLO(2,N).LT.-EPSL .AND. XMLO(2,N).LT.-EPSL ) THEN
            RHOLOX = -XLO(2,N)
            RHOLAX = 0.D+0
            RHOL(2,N) = RHOLWX + RHOLOX + RHOLAX
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Dissolved-oil initial condition provided as mass fraction  ---
!
          ELSEIF( XLO(2,N).GT.EPSL .AND. XMLO(2,N).LT.-EPSL ) THEN
            RHOLAX = 0.D+0
            RHOL(2,N) = RHOLWX/(1.D+0-XLO(2,N))
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
!
!---      Dissolved-oil initial condition provided as mole fraction  ---
!
          ELSEIF( XLO(2,N).LT.-EPSL .AND. XMLO(2,N).GT.EPSL ) THEN
            RHOLAX = 0.D+0
            RHOLOX = (XMLO(2,N)*RHOLWX*WTMO)/((1.D+0-XMLO(2,N))*WTMW)
            RHOL(2,N) = RHOLAX+RHOLOX+RHOLWX
            XLA(2,N) = RHOLAX/RHOL(2,N)
            XLO(2,N) = RHOLOX/RHOL(2,N)
            XLW(2,N) = RHOLWX/RHOL(2,N)
            WTMLX = XLW(2,N)/WTMW + XLO(2,N)/WTMO + XLA(2,N)/WTMA
            XMLA(2,N) = XLA(2,N)/WTMA/WTMLX
            XMLO(2,N) = XLO(2,N)/WTMO/WTMLX
            XMLW(2,N) = XLW(2,N)/WTMW/WTMLX
          ELSE
            XLA(2,N) = 0.D+0
            XLO(2,N) = 0.D+0
            XLW(2,N) = 0.D+0
            XMLA(2,N) = 0.D+0
            XMLO(2,N) = 0.D+0
            XMLW(2,N) = 0.D+0
          ENDIF
        ENDIF
        IF( XMLO(2,N).LT.EPSL ) XMLO(2,N) = 0.D+0
        IF( XMLA(2,N).LT.EPSL ) XMLA(2,N) = 0.D+0
        IF( XMLW(2,N).LT.EPSL ) XMLW(2,N) = 0.D+0
        IF( XLO(2,N).LT.EPSL ) XLO(2,N) = 0.D+0
        IF( XLA(2,N).LT.EPSL ) XLA(2,N) = 0.D+0
        IF( XLW(2,N).LT.EPSL ) XLW(2,N) = 0.D+0
  110 CONTINUE
!
!---  Return for restart simulations  ---
!
      IF( IEO.EQ.2 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---    Establish reference pressure for soil compressibility  ---
!
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        TCMP(N) = T(2,N)
        IZN = IZ(N)
        IF( CMP(3,IZN).GT.PATM ) THEN
          PCMP(N) = CMP(3,IZN)
        ELSEIF( ISLC(61).EQ.0 ) THEN
          PCMP(N) = MAX( PL(2,N),PG(2,N),PN(2,N) )+PATM
        ENDIF
  200 CONTINUE
!
!---  Determine initial phase conditions  ---
!
      INDX = 0
      DO 210 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 210
        IZN = IZ(N)
        SIG_GL = SCHR(5,IZN)
        SIG_NL = SCHR(6,IZN)
        SIG_GN = SCHR(7,IZN)
        RHOLOX = XLO(2,N)*RHOL(2,N)
        CALL SFT5( RHOLOX )
        BGL = SCHR(8,IZN)/SIG_GL
        BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
        BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
        ENPR = SCHR(2,IZN)*RHORL*GRAV
        IF( PN(2,N).LE.-1.D+20 ) THEN
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        PL(2,N)+ENPR/BNL )
          PN(2,N) = PN(2,N) + 1.D+4
        ENDIF
        PX = MAX( PL(2,N),PG(2,N),PN(2,N) ) + PATM
        HGL = MAX( ZERO,BGL*(PG(2,N)-PL(2,N)) )
        HNL = MAX( ZERO,BNL*(PN(2,N)-PL(2,N)) )
        HGN = MAX( ZERO,BGN*(PG(2,N)-PN(2,N)) )
        XMLOX = PSO(2,N)/HCOW
!
!---    No free NAPL conditions  ---
!
        IF( HGN-HNL.GE.-1.D-9 .OR. HNL.LE.ENPR ) THEN
!
!---      Aqueous-gas w/ no free NAPL conditions  ---
!
          IF( HGL.GT.ENPR ) THEN
!
!---        Aqueous-gas w/trapped NAPL conditions  ---
!
            IF( SCHR(11,IZN).GT.EPSL .AND. SNT(2,N).GT.EPSL ) THEN
              NPHAZ(2,N) = 6
!
!---          Kinetic volatilization option  ---
!
              IF( ISLC(46).GE.1 ) THEN
                PVO(2,N) = XMLO(2,N)*HCOW
              ELSE
                PVO(2,N) = PSO(2,N)
                XMLO(2,N) = PVO(2,N)/HCOW
              ENDIF
              PVW(2,N) = PSW(2,N)
              PVA(2,N) = MAX( PG(2,N)+PATM-PVO(2,N)-PVW(2,N),0.D+0 )
              IF( PG(2,N)+PATM+PTOLX.LT.PVO(2,N)+PVW(2,N) ) THEN
                CHMSG = 'Phase Condition #6: Gas Pressure < ' //
     &            'Oil+Water Vapor Pressure: Node: '
                INDX = 12
                IMSG = N
                CALL WRMSGS( INDX )
              ENDIF
              XMLA(2,N) = PVA(2,N)/HCAW
!
!---        Aqueous-gas w/dissolved-oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
              NPHAZ(2,N) = 5
              PVO(2,N) = XMLO(2,N)*HCOW
              PVW(2,N) = PSW(2,N)
              PVA(2,N) = MAX( PG(2,N)+PATM-PVO(2,N)-PVW(2,N),0.D+0 )
              IF( XMLO(2,N).GT.(VSLC(1)*XMLOX+PTOLX/HCOW) ) THEN
                CHMSG = 'Phase Condition #5: Initial Dissolved-Oil ' //
     &            'Concentration > Solubility Limit: Node: '
                INDX = 12
                IMSG = N
                CALL WRMSGS( INDX )
              ENDIF
              IF( PG(2,N)+PATM+PTOLX.LT.PVO(2,N)+PVW(2,N) ) THEN
                CHMSG = 'Phase Condition #5: Gas Pressure < ' //
     &            'Oil+Water Vapor Pressure: Node: '
                INDX = 12
                IMSG = N
                CALL WRMSGS( INDX )
              ENDIF
              XMLA(2,N) = PVA(2,N)/HCAW
!
!---        Aqueous-gas w/no oil dissolution  ---
!
            ELSE
              NPHAZ(2,N) = 9
              PVO(2,N) = XMLO(2,N)*HCOW
              PVW(2,N) = PSW(2,N)
              PVA(2,N) = MAX( PG(2,N)+PATM-PVO(2,N)-PVW(2,N),0.D+0 )
              IF( XMLO(2,N).GT.(VSLC(1)*XMLOX+PTOLX/HCOW) ) THEN
                CHMSG = 'Phase Condition #9: Initial Dissolved-Oil ' //
     &            'Concentration > Solubility Limit: Node: '
                INDX = 12
                IMSG = N
                CALL WRMSGS( INDX )
              ENDIF
              IF( PG(2,N)+PATM+PTOLX.LT.PVO(2,N)+PVW(2,N) ) THEN
                CHMSG = 'Phase Condition #9: Gas Pressure < ' //
     &            'Oil+Water Vapor Pressure: Node: '
                INDX = 12
                IMSG = N
                CALL WRMSGS( INDX )
              ENDIF
              XMLA(2,N) = PVA(2,N)/HCAW
            ENDIF
            PN(2,N) = (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN)
!
!---      Aqueous w/ no free NAPL conditions  ---
!
          ELSE
!
!---        Aqueous w/ trapped NAPL conditions  ---
!
            IF( SCHR(11,IZN).GT.EPSL .AND. SNT(2,N).GT.EPSL ) THEN
!
!---          Aqueous w/ trapped gas and trapped NAPL conditions  ---
!
              IF( SCHR(17,IZN).GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
                NPHAZ(2,N) = 4
!
!---            Kinetic volatilization option  ---
!
                IF( ISLC(46).GE.1 ) THEN
                  PVO(2,N) = XMLO(2,N)*HCOW
                ELSE
                  PVO(2,N) = PSO(2,N)
                  XMLO(2,N) = PVO(2,N)/HCOW
                ENDIF
                PVW(2,N) = PSW(2,N)
                PG(2,N) = PL(2,N) + ENPR/BGL
                PVA(2,N) = MAX( PG(2,N)+PATM-PVO(2,N)-PVW(2,N),0.D+0 )
                PGX = PVO(2,N)+PVW(2,N)+PVA(2,N)-PATM
                IF( PG(2,N)+PTOLX.LT.PGX ) THEN
                  CHMSG = 'Phase Condition #4: Gas Pressure < ' //
     &              'Oil+Water Vapor Pressure: Node: '
                  INDX = 12
                  IMSG = N
                  CALL WRMSGS( INDX )
                ENDIF
                XMLA(2,N) = PVA(2,N)/HCAW
!
!---          Aqueous w/ dissolved-air and trapped-NAPL conditions  ---
!
              ELSE
                NPHAZ(2,N) = 3
!
!---            Kinetic volatilization option  ---
!
                IF( ISLC(46).GE.1 ) THEN
                  PVO(2,N) = XMLO(2,N)*HCOW
                ELSE
                  PVO(2,N) = PSO(2,N)
                  XMLO(2,N) = PVO(2,N)/HCOW
                ENDIF
                PVW(2,N) = PSW(2,N)
                PVA(2,N) = XMLA(2,N)*HCAW
                PG(2,N) = PL(2,N) + ENPR/BGL
                PGX = PVO(2,N)+PVW(2,N)+PVA(2,N)-PATM
                IF( PG(2,N)+PTOLX.LT.PGX ) THEN
                  CHMSG = 'Phase Condition #3: Gas Pressure < ' //
     &              'Oil+Water Vapor Pressure: Node: '
                  INDX = 12
                  IMSG = N
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
!
!---        Aqueous w/ dissolved oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
!
!---          Aqueous w/ trapped gas and dissolved oil conditions  ---
!
              IF( SCHR(17,IZN).GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
                NPHAZ(2,N) = 2
                PVO(2,N) = XMLO(2,N)*HCOW
                PVW(2,N) = PSW(2,N)
                PG(2,N) = PL(2,N) + ENPR/BGL
                PVA(2,N) = MAX( PG(2,N)+PATM-PVO(2,N)-PVW(2,N),0.D+0 )
                PGX = PVO(2,N)+PVW(2,N)+PVA(2,N)
                IF( XMLO(2,N).GT.(VSLC(1)*XMLOX+PTOLX/HCOW) ) THEN
                  CHMSG = 'Phase Condition #2: Initial Dissolved-Oil '//
     &              'Concentration > Solubility Limit: Node: '
                  INDX = 12
                  IMSG = N
                  CALL WRMSGS( INDX )
                ENDIF
                IF( PG(2,N)+PTOLX.LT.PGX ) THEN
                  CHMSG = 'Phase Condition #2: Gas Pressure < ' //
     &              'Oil+Water Vapor Pressure: Node: '
                  INDX = 12
                  IMSG = N
                  CALL WRMSGS( INDX )
                ENDIF
                XMLA(2,N) = PVA(2,N)/HCAW
!
!---          Aqueous w/ dissolved air and dissolved oil conditions  ---
!
              ELSE
                NPHAZ(2,N) = 1
                PVO(2,N) = XMLO(2,N)*HCOW
                PVW(2,N) = PSW(2,N)
                PVA(2,N) = XMLA(2,N)*HCAW
                PG(2,N) = PL(2,N) + ENPR/BGL
                PGX = PVO(2,N)+PVW(2,N)+PVA(2,N)-PATM
                IF( XMLO(2,N).GT.(VSLC(1)*XMLOX+PTOLX/HCOW) ) THEN
                  CHMSG = 'Phase Condition #1: Initial Dissolved-Oil '//
     &              'Concentration > Solubility Limit: Node: '
                  INDX = 12
                  IMSG = N
                  CALL WRMSGS( INDX )
                ENDIF
                IF( PG(2,N)+PTOLX.LT.PGX ) THEN
                  CHMSG = 'Phase Condition #1: Gas Pressure < ' //
     &              'Oil+Water Vapor Pressure: Node: '
                  INDX = 12
                  IMSG = N
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
!
!---        Aqueous w/ no oil dissolution  ---
!
            ELSE
!
!---          Aqueous w/ trapped gas and no oil dissolution
!             conditions  ---
!
              IF( SCHR(17,IZN).GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
                NPHAZ(2,N) = 8
                PVO(2,N) = XMLO(2,N)*HCOW
                PVW(2,N) = PSW(2,N)
                PG(2,N) = PN(2,N) + ENPR/BGN
                PGX = PVO(2,N)+PVW(2,N)-PATM
                PVA(2,N) = MAX( PG(2,N)+PATM-PVO(2,N)-PVW(2,N),0.D+0 )
                IF( XMLO(2,N).GT.(VSLC(1)*XMLOX+PTOLX/HCOW) ) THEN
                  CHMSG = 'Phase Condition #8: Initial Dissolved-Oil '//
     &              'Concentration > Solubility Limit: Node: '
                  INDX = 12
                  IMSG = N
                  CALL WRMSGS( INDX )
                ENDIF
                IF( PG(2,N)+PTOLX.LT.PGX ) THEN
                  CHMSG = 'Phase Condition #8: Gas Pressure < ' //
     &              'Oil+Water Vapor Pressure: Node: '
                  INDX = 12
                  IMSG = N
                  CALL WRMSGS( INDX )
                ENDIF
                XMLA(2,N) = PVA(2,N)/HCAW
!
!---          Aqueous w/ dissolved air and no oil dissolution
!             conditions  ---
!
              ELSE
                NPHAZ(2,N) = 7
                PVO(2,N) = XMLO(2,N)*HCOW
                PVW(2,N) = PSW(2,N)
                PVA(2,N) = XMLA(2,N)*HCAW
                PG(2,N) = PN(2,N) + ENPR/BGN
                PGX = PVO(2,N)+PVW(2,N)+PVA(2,N)-PATM
                IF( XMLO(2,N).GT.(VSLC(1)*XMLOX+PTOLX/HCOW) ) THEN
                  CHMSG = 'Phase Condition #7: Initial Dissolved-Oil '//
     &              'Concentration > Solubility Limit: Node: '
                  INDX = 12
                  IMSG = N
                  CALL WRMSGS( INDX )
                ENDIF
                IF( PG(2,N)+PTOLX.LT.PGX ) THEN
                  CHMSG = 'Phase Condition #7: Gas Pressure < ' //
     &              'Oil+Water Vapor Pressure: Node: '
                  INDX = 12
                  IMSG = N
                  CALL WRMSGS( INDX )
                ENDIF
              ENDIF
            ENDIF
            PN(2,N) = PL(2,N)+ENPR/BNL
          ENDIF
!
!---    Aqueous w/ free NAPL conditions  ---
!
        ELSE
!
!---      Aqueous w/ free gas and free NAPL conditions  ---
!
          IF( HGN.GT.ENPR ) THEN
            NPHAZ(2,N) = 9
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              PVO(2,N) = XMLO(2,N)*HCOW
            ELSE
              PVO(2,N) = PSO(2,N)
              XMLO(2,N) = PVO(2,N)/HCOW
            ENDIF
            PVW(2,N) = PSW(2,N)
            PVA(2,N) = MAX( PG(2,N)+PATM-PVO(2,N)-PVW(2,N),0.D+0 )
            IF( PG(2,N)+PATM+PTOLX.LT.PVO(2,N)+PVW(2,N) ) THEN
              CHMSG = 'Phase Condition #9: Gas Pressure < ' //
     &          'Oil+Water Vapor Pressure: Node: '
              INDX = 12
              IMSG = N
              CALL WRMSGS( INDX )
            ENDIF
            XMLA(2,N) = PVA(2,N)/HCAW
!
!---      Aqueous w/ no free gas and free NAPL conditions  ---
!
          ELSE
!
!---        Aqueous w/ trapped gas and free NAPL conditions  ---
!
            IF( SCHR(17,IZN).GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
              NPHAZ(2,N) = 8
!
!---          Kinetic volatilization option  ---
!
              IF( ISLC(46).GE.1 ) THEN
                PVO(2,N) = XMLO(2,N)*HCOW
              ELSE
                PVO(2,N) = PSO(2,N)
                XMLO(2,N) = PVO(2,N)/HCOW
              ENDIF
              PVW(2,N) = PSW(2,N)
              PG(2,N) = PN(2,N) + ENPR/BGN
              PGX = PVO(2,N)+PVW(2,N)-PATM
              PVA(2,N) = MAX( PG(2,N)+PATM-PVO(2,N)-PVW(2,N),0.D+0 )
              IF( PG(2,N)+PTOLX.LT.PGX ) THEN
                CHMSG = 'Phase Condition #8: Gas Pressure < ' //
     &            'Oil+Water Vapor Pressure: Node: '
                INDX = 12
                IMSG = N
                CALL WRMSGS( INDX )
              ENDIF
              XMLA(2,N) = PVA(2,N)/HCAW
!
!---        Aqueous w/ dissolved air and free NAPL conditions  ---
!
            ELSE
              NPHAZ(2,N) = 7
!
!---          Kinetic volatilization option  ---
!
              IF( ISLC(46).GE.1 ) THEN
                PVO(2,N) = XMLO(2,N)*HCOW
              ELSE
                PVO(2,N) = PSO(2,N)
                XMLO(2,N) = PVO(2,N)/HCOW
              ENDIF
              PVW(2,N) = PSW(2,N)
              PVA(2,N) = XMLA(2,N)*HCAW
              PG(2,N) = PN(2,N) + ENPR/BGN
              PGX = PVO(2,N)+PVW(2,N)+PVA(2,N)-PATM
              IF( PG(2,N)+PTOLX.LT.PGX ) THEN
                CHMSG = 'Phase Condition #7: Gas Pressure < ' //
     &            'Oil+Water Vapor Pressure: Node: '
                INDX = 12
                IMSG = N
                CALL WRMSGS( INDX )
              ENDIF
            ENDIF
          ENDIF
        ENDIF
!
!---    Call the saturation function to initialize the
!       apparent total-liquid and apparent aqueous saturations  ---
!
        INDX = 2
        ASNT(N) = SNT(2,N)/(1.D+0-SCHR(4,IZN))
        ASNR(N) = SNR(2,N)/(1.D+0-SCHR(4,IZN))
        ASGT(N) = SGT(2,N)/(1.D+0-SCHR(4,IZN))
        TRPNLX = SCHR(9,IZN)/(SCHR(9,IZN)+TRPNL(2,N))
        ASNTMX = SCHR(11,IZN)*TRPNLX
        ASNRMX = SCHR(10,IZN)
        CALL KSP5( PG(2,N),PN(2,N),PL(2,N),SLX,SNX,SGX,
     &    RKLX,RKNX,RKGX,AST(N),ASL(N),
     &    ASLMINX,ASTMINX,ASTMAXX,ASNT(N),ASNR(N),ASGT(N),
     &    ASGTL(N),ASGTN(N),ASNTMX,ASNRMX,ASGTMX,
     &    INDX,IZN,N )
!
!---    Residual NAPL check  ---
!
        IF( ASNRMX.GE.EPSL ) THEN
          ASTMAXX = 1.D+0
          ASNRX = ESNR( ASL(N),AST(N),ASTMAXX,ASNRMX )
          ASNRX = MIN( AST(N)-ASL(N),ASNRX )
        ELSE
          ASNRX = 0.D+0
        ENDIF
!
!---    Residual-NAPL saturation greater than mobile-residual-NAPL
!       saturation stop simulation with error message  ---
!
        IF( ASNR(N).GT.MAX((AST(N)-ASL(N)),EPSL) ) THEN
          CHMSG = 'Residual-NAPL Saturation > ' //
     &      'Free-NAPL Saturation: Node'
          INDX = 16
          IMSG = N
          RLMSG = ASNR(N)
          CALL WRMSGS( INDX )
!
!---    Residual-NAPL saturation greater than maximum,
!       stop simulation with error message  ---
!
        ELSEIF( ASNR(N).GT.MAX(ASNRX,EPSL) ) THEN
          CHMSG = 'Residual-NAPL Saturation > ' //
     &      'Maximum Possible: Node'
          INDX = 16
          IMSG = N
          RLMSG = ASNR(N)
          CALL WRMSGS( INDX )
!
!---    No residual NAPL, set maximum apparent-total-liquid
!       saturation  ---
!
        ELSEIF( ASNR(N).LE.EPSL ) THEN
          ASTMAX(2,N) = ASL(N)
!
!---    Residual NAPL within limits, set maximum apparent-total-liquid
!       saturation  ---
!
        ELSE
          ASTMAX(2,N) = ASTMX( ASL(N),AST(N),ASNR(N),ASNRMX )
!
!---      Apparent maximum total-liquid saturation less than
!         apparent total-liquid saturation, calculate
!         residual-NAPL saturation using apparent
!         total liquid saturation; default, where the residual-NAPL
!         saturation is uninitialized  ---
!
          IF( AST(N).GT.ASTMAX(2,N) ) THEN
            ASTMAX(2,N) = AST(N)
            ASNR(N) = ESNR( ASL(N),AST(N),ASTMAX(2,N),ASNRMX )
            ASNRX = MIN( AST(N)-ASL(N),ASNRX )
          ENDIF
        ENDIF
  210 CONTINUE
!
!---  Compute the total trapping number  ---
!
!      CALL TRPNL5
!
!---  Compute the initial effective aqueous and total-liquid
!     minimum saturations, and check for excessive trapped-NAPL and
!     trapped-gas saturations  ---
!
      DO 220 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 220
        IZN = IZ(N)
        SIG_GL = SCHR(5,IZN)
        SIG_NL = SCHR(6,IZN)
        SIG_GN = SCHR(7,IZN)
        RHOLOX = XLO(2,N)*RHOL(2,N)
        CALL SFT5( RHOLOX )
        BGL = SCHR(8,IZN)/SIG_GL
        BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
        BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
!
!---  Trapped NAPL check  ---
!
        ASNTMX = SCHR(11,IZN)
        IF( ASNTMX.GE.EPSL ) THEN
          PLND = 1.D+0/ASNTMX - 1.D+0
          ASNTX = ASNTMX - (1.D+0-ASL(N))/(1.D+0 + PLND*(1.D+0-ASL(N)))
          ASNTX = MAX( ASNTX,ZERO )
        ELSE
          ASNTX = 0.D+0
        ENDIF
        ASNT(N) = SNT(2,N)/(1.D+0-SCHR(4,IZN))
        IF( ASNT(N).GT.ASNTX ) THEN
          CHMSG = 'Trapped NAPL Saturation: Node: '
          INDX = 16
          IMSG = N
          RLMSG = ASNT(N)
          CALL WRMSGS( INDX )
        ELSEIF( ASNT(N).LE.EPSL ) THEN
          ASLMIN(2,N) = 1.D+0
        ELSE
          ASNRX = (1.D+0-ASL(N))/(1.D+0+PLND*(1.D+0-ASL(N)))
          ASLMIN(2,N) = (1.D+0-ASNRX-PLND*ASNRX-ASNT(N)-PLND*ASNT(N))/
     &      (1.D+0-PLND*ASNRX-PLND*ASNT(N))
        ENDIF
        IF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.4 .OR. NPHAZ(2,N).EQ.6 )
     &    SN(2,N) = SNT(2,N)
!
!---  Trapped gas check  ---
!
        ASGTMX = SCHR(17,IZN)
        IF( ASGTMX.GE.EPSL ) THEN
          PLND = 1.D+0/ASGTMX - 1.D+0
          ASGTX = ASGTMX - (1.D+0-AST(N))/(1.D+0+PLND*(1.D+0-AST(N)))
          ASGTX = MAX( ASGTX,ZERO )
        ELSE
          ASGTX = 0.D+0
        ENDIF
        ASGTL(N) = SGTL(N)/(1.D+0-SCHR(4,IZN))
        ASGTN(N) = SGTN(N)/(1.D+0-SCHR(4,IZN))
        IF( ASGTL(N)+ASGTN(N).GT.ASGTX ) THEN
          CHMSG = 'Trapped Gas Saturation: Node: '
          INDX = 16
          IMSG = N
          RLMSG = ASGTL(N)+ASGTN(N)
          CALL WRMSGS( INDX )
        ELSEIF( ASGTL(N)+ASGTN(N).LE.EPSL ) THEN
          ASTMIN(2,N) = AST(N)
        ELSE
          ASGTX = ASGTL(N) + ASGTN(N)
          ASGR = (1.D+0-AST(N))/(1.D+0+PLND*(1.D+0-AST(N)))
          ASTMIN(2,N) = (1.D+0-ASGR-PLND*ASGR-ASGTX-PLND*ASGTX)/
     &      (1.D+0-PLND*ASGR-PLND*ASGTX)
        ENDIF
        IF( NPHAZ(2,N).EQ.2 .OR. NPHAZ(2,N).EQ.4 .OR. NPHAZ(2,N).EQ.8 )
     &    SG(2,N) = SGT(2,N)
        ASLMIN(1,N) = ASLMIN(2,N)
        ASTMIN(1,N) = ASTMIN(2,N)
  220 CONTINUE
!
!---    Establish reference pressure for soil compressibility  ---
!
      DO 240 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 240
        IZN = IZ(N)
        IF( CMP(3,IZN).GT.PATM ) THEN
          PCMP(N) = CMP(3,IZN)
        ELSEIF( ISLC(61).EQ.0 ) THEN
          PCMP(N) = MAX( PL(2,N),PG(2,N),PN(2,N) )+PATM
        ENDIF
  240 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHK5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CISC5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute initial solute concentrations.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, May 1994.
!     Last Modified by MD White, Battelle, June 30, 1994.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/CISC5'
!
!---  Loop over all solutes  ---
!
      DO 140 NSL = 1,NSOLU
        DO 110 N = 1,NFLD
          IF( IXP(N).EQ.0 ) GOTO 110
          IZN = IZ(N)
          XVS = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
          XVL = SL(2,N)*PORD(2,N)
          XVN = SN(2,N)*PORD(2,N)
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
            C(N,NSL) = C(N,NSL)*(XVS + XVL + XVN/PCLN(1,NSL) +XVG*PCGLX)
          ELSEIF( ICT(N,NSL).EQ.3 ) THEN
            C(N,NSL) = C(N,NSL)*((XVS + XVL + XVN)/PCGLX + XVG)
          ELSEIF( ICT(N,NSL).EQ.4 ) THEN
            C(N,NSL) = C(N,NSL)*((XVS +XVL +XVG*PCGLX)*PCLN(1,NSL) +XVN)
          ENDIF
  110   CONTINUE
!
!---  Assign boundary solute concentrations for initial condition
!     type boundary conditions  ---
!
        DO 130 NB = 1,NBC
          IF( IBCT(NSL+LUK,NB).EQ.12 ) THEN
            N = IBCN(NB)
            CBO(NB,NSL) = C(N,NSL)
          ENDIF
  130   CONTINUE
  140 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CISC5 group  ---
!
      RETURN
      END
!
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ELC5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
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
!     Correct aqueous liquid density and viscosity for electrolyte
!     solute concentration
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, October 2000.
!     Last Modified by MD White, PNNL, October 11, 2000.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ELC5'
      DO 400 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 400
        PGX = PG(2,N) + PATM
        DO 300 M = 2,ISVC+2
          PLX = PL(M,N) + PATM
          CALL WATSP( T(2,N),PVW(M,N) )
          PX = MAX( PLX,PGX,PVW(M,N) )
          CALL WATLQD( T(2,N),PX,RHOL(M,N) )
          CALL WATLQV( T(2,N),PX,PSW(2,N),VISL(M,N) )
          CLX = C(N,NSL_ELC)*YL(N,NSL_ELC)/(SL(M,N)*PORD(M,N)+SMALL)
          XLW(M,N) = RHOL(M,N)
          CALL ELC_DEN( RHOL(M,N),CLX,ELC_DCF )
          XLW(M,N) = XLW(M,N)/RHOL(M,N)
          CALL ELC_VIS( VISL(M,N),CLX,ELC_VCF )
  300   CONTINUE
  400 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ELC5 group  ---
!
      RETURN
      END
!
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRM5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute primary variable increments
!
!     Phase Conditions
!     Number              Description                   Primary Unknowns
!       1    Aqueous w/ Dissolved Oil w/ Dissolved Air    PL,XMLO,XMLA
!            PG = MAX( PL+ENPR/BGL, PVO+PVW+PVA-PATM )
!            PG > PL+ENPR/BGL --> Free Gas
!       2    Aqueous w/ Dissolved Oil w/ Trapped Gas      PL,XMLO,SG
!            PG = PL+ENPR/BGL
!            PVA = MAX( PG+PATM-PVO-PVW, 0. )
!       3    Aqueous w/ Trapped NAPL w/ Dissolved Air     PL,SN,XMLA
!            PG = MAX( PL+ENPR/BGL, PSO+PVW+PVA-PATM )
!            PG > PL+ENPR/BGL --> Free Gas
!       4    Aqueous w/ Trapped NAPL w/ Trapped Gas       PL,SN,SG
!            PG = PL+ENPR/BGL
!            PVA = MAX( PG+PATM-PSO-PVW, 0. )
!       5    Aqueous-Gas w/ Dissolved Oil                 PL,XMLO,PG
!            PVA = MAX( PG+PATM-PVO-PVW, 0. )
!       6    Aqueous-Gas w/ Trapped NAPL                  PL,SN,PG
!            PVA = MAX( PG+PATM-PSO-PVW, 0. )
!       7    Aqueous-NAPL w/ Dissolved Air                PL,PN,XMLA
!            PG = MAX( PN+ENPR/BGN, PSO+PVW+PVA-PATM )
!            PG > PN+ENPR/BGN --> Free Gas
!       8    Aqueous-NAPL w/ Trapped Gas                  PL,PN,SG
!            PG = MAX( PN+ENPR/BGN, PSO+PVW+PVA-PATM )
!            PVA = MAX( PG+PATM-PSO-PVW, 0. )
!       9    Aqueous-NAPL-Gas                             PL,PN,PG
!            PVA = MAX( PG+PATM-PSO-PVW, 0. )
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, October 20, 1997.
!     Last Modified by MD White, PNNL, 27 June 2002.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FDVP
      USE FDVD
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
!----------------------Data Statements---------------------------------!
!
      SAVE TOL1,TOL2
      DATA TOL1,TOL2 /0.D+0,0.D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INCRM5'
!
!---  Phase options, compute phase condition   ---
!
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 100
        IZN = IZ(N)
        XMLOX = MIN( XMLO(2,N),PSO(2,N)/HCOW )
        XMLWX = MAX( 1.D+0-XMLOX,0.D+0 )
        XLOX = XMLOX*WTMO/(XMLOX*WTMO + XMLWX*WTMW)
        PX = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
        CALL WATLQD( T(2,N),PX,RHOLX )
        SIG_GL = SCHR(5,IZN)
        SIG_NL = SCHR(6,IZN)
        SIG_GN = SCHR(7,IZN)
        RHOLOX = XLOX*RHOLX
        CALL SFT5( RHOLOX )
        BGL = SCHR(8,IZN)/SIG_GL
        BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
        BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
        ENPR = SCHR(2,IZN)*RHORL*GRAV
        HGL = MAX( ZERO,BGL*(PG(2,N)-PL(2,N)) )
        HNL = MAX( ZERO,BNL*(PN(2,N)-PL(2,N)) )
        HGN = MAX( ZERO,BGN*(PG(2,N)-PN(2,N)) )
        XMLOX = PSO(2,N)/HCOW
!
!---    Aqueous w/ dissolved-oil, dissolved-air conditions
!       Unknowns: PL,XMLA,XMLO
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,PN,XMLO  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          PG(2,N) = PL(2,N)+ENPR/BGL
          PVOX = MIN( XMLO(2,N)*HCOW,PSO(2,N) )
          PVAX = XMLA(2,N)*HCAW
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
          PGX = PVOX + PVAX + PVW(2,N)
!
!---      Gas- and vapor-partial pressures exceed gas pressure,
!         transition from aqueous-saturated to free-gas conditions  ---
!
          IF( PGX.GT.VSLC(1)*(PG(2,N)+PATM)  ) THEN
!
!---        Dissolved-oil mole fraction exceeds solubility limit,
!           transition from dissolved-oil to free-NAPL conditions,
!           bump gas and NAPL pressure by 0.1 Pa to create gas
!           and NAPL saturations  --
!
            IF( XMLO(2,N).GT.XMLOX ) THEN
              PG(2,N) = PG(2,N) + 1.D-1
              PN(2,N) = PN(2,N) + 1.D-1
              NPHAZ(2,N) = 9
!
!---        Dissolved-oil mole fraction less than solubility limit,
!           no transition from dissolved-oil conditions,
!           bump gas pressure by 0.1 Pa to create gas saturation  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
              PG(2,N) = PG(2,N) + 1.D-1
              NPHAZ(2,N) = 5
!
!---        Transition from dissolved-oil to zero-solubility conditions,
!           bump gas pressure by 0.1 Pa to create gas saturation
!           (transition should not be possible)  ---
!
            ELSE
              PG(2,N) = PG(2,N) + 1.D-1
              NPHAZ(2,N) = 9
            ENDIF
!
!---        Gas- and vapor-partial pressures less than gas pressure,
!           no transition from aqueous-saturated conditions  ---
!
          ELSE
!
!---        Dissolved-oil mole fraction exceeds solubility limit,
!           transition from dissolved-oil to free-NAPL conditions,
!           bump NAPL pressure by 0.1 Pa to create NAPL saturation  --
!
            IF( XMLO(2,N).GT.XMLOX ) THEN
              NPHAZ(2,N) = 7
              PN(2,N) = PN(2,N) + 1.D-1
!
!---        Dissolved-oil mole fraction less than solubility limit,
!           no transition from dissolved-oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
              NPHAZ(2,N) = 1
!
!---        Transition from dissolved-oil to zero-solubility 
!           conditions (transition should not be possible)  ---
!
            ELSE
              NPHAZ(2,N) = 7
            ENDIF
          ENDIF
!
!---    Aqueous w/ dissolved-oil, trapped-gas conditions
!       Unknowns: PL,SG,XMLO
!
!       Kinetic volatilization option
!       Unknowns: PL,SG,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          PG(2,N) = PL(2,N) + ENPR/BGL
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
          ASGTMX = SCHR(17,IZN)
          IF( ASGTMX.GE.EPSL ) THEN
            PLND = 1.D+0/(ASGTMX+SMALL) - 1.D+0
            ASGTX = (1.D+0-ASTMIN(1,N))/(1.D+0+PLND*(1.D+0-ASTMIN(1,N)))
          ELSE
            ASGTX = 0.D+0
          ENDIF
          ASGX = SG(2,N)/(1.D+0-SCHR(4,IZN))
!
!---      Transition from trapped-gas to free-gas conditions  ---
!
          IF( (ASGX-ASGTX)/(ASGTMX+EPSL).GT.1.D-1 ) THEN
!
!---        Transition from dissolved-oil to free-NAPL conditions  --
!
            IF( XMLO(2,N).GT.XMLOX ) THEN
              PG(2,N) = PG(2,N) + 1.D-1
              PN(2,N) = PN(2,N) + 1.D-1
              NPHAZ(2,N) = 9
!
!---        No transition from dissolved-oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
              PG(2,N) = PG(2,N) + 1.D-1
              NPHAZ(2,N) = 5
!
!---        Transition from dissolved-oil to zero-solubility
!           conditions  ---
!
            ELSE
              PG(2,N) = PG(2,N) + 1.D-1
              NPHAZ(2,N) = 9
            ENDIF
!
!---      Transition from trapped-gas to dissolved-air conditions  ---
!
          ELSEIF( ASGX.LT.EPSL ) THEN
!
!---        Transition from dissolved-oil to free-NAPL conditions  --
!
            IF( XMLO(2,N).GT.XMLOX ) THEN
              PN(2,N) = PN(2,N) + 1.D-1
              NPHAZ(2,N) = 7
!
!---        No transition from dissolved-oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
              NPHAZ(2,N) = 1
!
!---        Transition from dissolved-oil to zero-solubility
!           conditions  ---
!
            ELSE
              NPHAZ(2,N) = 7
            ENDIF
!
!---      Trapped-gas to trapped-gas transition
!
          ELSE
!
!---        Transition from dissolved-oil to free-NAPL conditions  --
!
            IF( XMLO(2,N).GT.XMLOX ) THEN
              PN(2,N) = PN(2,N) + 1.D-1
              NPHAZ(2,N) = 8
!
!---        No transition from dissolved-oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
              NPHAZ(2,N) = 2
!
!---        Transition from dissolved-oil to zero-solubility
!           conditions  ---
!
            ELSE
              NPHAZ(2,N) = 8
            ENDIF
          ENDIF
!
!---    Aqueous w/ trapped-NAPL, dissolved-air conditions
!       Unknowns: PL,XMLA,SN
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,SN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
!
!---      Set gas pressure to the entry pressure,
!         set NAPL pressure to the critical pressure  ---
!
          PG(2,N) = PL(2,N)+ENPR/BGL
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
!
!---      Sum gas partial pressures  ---
!
          PVAX = XMLA(2,N)*HCAW
          PGX = PSO(2,N) + PVAX + PVW(2,N)
!
!---      Effective trapped-NAPL saturation  ---
!
          ASNX = SN(2,N)/(1.D+0-SCHR(4,IZN))
!
!---      Maximum effective trapped-NAPL saturation and
!         minimum effective aqueous saturation  ---
!
          TRPNLX = (SCHR(9,IZN)+TRPNL(2,N)*SCHR(8,IZN)/
     &      (BNL*SIG_NL))/(SCHR(9,IZN)+TRPNL(2,N))
          ASNTMX = SCHR(11,IZN)*TRPNLX
          IF( ASNTMX.GE.EPSL ) THEN
            PLND = 1.D+0/ASNTMX - 1.D+0
            ASNTX = (1.D+0-ASLMIN(1,N))/(1.D+0+PLND*(1.D+0-ASLMIN(1,N)))
          ELSE
            ASNTX = 0.D+0
          ENDIF
!
!---      Transition from aqueous-saturated to free-gas conditions  ---
!
          IF( PGX.GT.(PG(2,N)+PATM)  ) THEN
!
!---        Trapped-NAPL saturation exceeds maximum trapped-NAPL
!           saturation, transition to mobile-residual NAPL conditions,
!           bump gas and NAPL pressure by 0.1 Pa to create gas
!           and NAPL saturations  ---
!
            IF( (ASNX-ASNTX)/(ASNTMX+EPSL).GT.1.D-1 ) THEN
              PG(2,N) = PG(2,N) + 1.D-1
              PN(2,N) = PN(2,N) + 1.D-1
              NPHAZ(2,N) = 9
!
!---        Transition from trapped-NAPL to dissolved-oil
!           conditions  ---
!
            ELSEIF( ASNX.LT.EPSL ) THEN
!
!---          Transition from trapped-NAPL to dissolved-oil conditions,
!             bump gas pressure by 0.1 Pa to create gas saturation  ---
!
              IF( XMLOX.GT.1.D-12 ) THEN
                PG(2,N) = PG(2,N) + 1.D-1
                NPHAZ(2,N) = 5
!
!---          Transition from trapped-NAPL to zero-solubility condition,
!             bump gas pressure by 0.1 Pa to create gas saturation  ---
!
              ELSE
                PG(2,N) = PG(2,N) + 1.D-1
                NPHAZ(2,N) = 9
              ENDIF
!
!---        No transition from trapped-NAPL conditions,
!           bump gas pressure by 0.1 Pa to create gas saturation  ---
!
            ELSE
              PG(2,N) = PG(2,N) + 1.D-1
              NPHAZ(2,N) = 6
            ENDIF
!
!---      Gas- and vapor-partial pressures less than gas pressure,
!         no transition from aqueous-saturated conditions  ---
!
          ELSE
!
!---        Trapped-NAPL saturation exceeds maximum trapped-NAPL
!           saturation, transition to mobile-residual NAPL conditions,
!           bump NAPL pressure by 0.1 Pa to create NAPL saturation  ---
!
            IF( (ASNX-ASNTX)/(ASNTMX+EPSL).GT.1.D-1 ) THEN
              PN(2,N) = PN(2,N) + 1.D-1
              NPHAZ(2,N) = 7
!
!---        Transition from trapped-NAPL to dissolved-oil
!           conditions  ---
!
            ELSEIF( ASNX.LT.EPSL ) THEN
!
!---          Transition from trapped-NAPL to dissolved-oil
!             conditions  ----
!
              IF( XMLOX.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 1
!
!---          Transition from trapped-NAPL to zero-solubility
!             conditions  ----
!
              ELSE
                NPHAZ(2,N) = 7
              ENDIF
!
!---        No transition from trapped-NAPL conditions  ---
!
            ELSE
              NPHAZ(2,N) = 3
            ENDIF
          ENDIF
!
!---    Aqueous w/ trapped-NAPL, trapped-gas conditions
!       Unknowns: PL,SG,SN
!
!       Kinetic volatilization option
!       Unknowns: PL,SG,SN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
          PG(2,N) = PL(2,N) + ENPR/BGL
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
          ASGTMX = SCHR(17,IZN)
          IF( ASGTMX.GE.EPSL ) THEN
            PLND = 1.D+0/ASGTMX - 1.D+0
            ASGTX = (1.D+0-ASTMIN(1,N))/(1.D+0+PLND*(1.D+0-ASTMIN(1,N)))
          ELSE
            ASGTX = 0.D+0
          ENDIF
          ASGX = SG(2,N)/(1.D+0-SCHR(4,IZN))
          ASNTMX = SCHR(11,IZN)*(SCHR(9,IZN)+TRPNL(2,N)*SCHR(8,IZN)/
     &      (BNL*SIG_NL))/(SCHR(9,IZN)+TRPNL(2,N))
          IF( ASNTMX.GE.EPSL ) THEN
            PLND = 1.D+0/ASNTMX - 1.D+0
            ASNTX = (1.D+0-ASLMIN(1,N))/(1.D+0+PLND*(1.D+0-ASLMIN(1,N)))
          ELSE
            ASNTX = 0.D+0
          ENDIF
          ASNX = SN(2,N)/(1.D+0-SCHR(4,IZN))
!
!---      Transition from trapped-gas to free-gas conditions  ---
!
          IF( (ASGX-ASGTX)/(ASGTMX+EPSL).GT.1.D-1 ) THEN
!
!---        Transition from trapped-NAPL to free-NAPL conditions  ---
!
            IF( (ASNTX-ASNX).LT.EPSL ) THEN
              PG(2,N) = PG(2,N) + 1.D-1
              PN(2,N) = PN(2,N) + 1.D-1
              NPHAZ(2,N) = 9
!
!---        Transition from trapped-NAPL to dissolved-oil
!           conditions  ---
!
            ELSEIF( ASNX.LT.EPSL ) THEN
!
!---          Transition from trapped-NAPL to dissolved-oil
!             conditions  ---
!
              IF( XMLOX.GT.1.D-12 ) THEN
                PG(2,N) = PG(2,N) + 1.D-1
                NPHAZ(2,N) = 5
!
!---          Transition from trapped-NAPL to zero-solubility
!             conditions  ---
!
              ELSE
                PG(2,N) = PG(2,N) + 1.D-1
                NPHAZ(2,N) = 9
              ENDIF
!
!---        No transition from trapped-NAPL conditions  ---
!
            ELSE
              PG(2,N) = PG(2,N) + 1.D-1
              NPHAZ(2,N) = 6
            ENDIF
!
!---      Transition from trapped-gas to dissolved-air conditions  ---
!
          ELSEIF( ASGX.LT.EPSL ) THEN
!
!---        Transition from trapped-NAPL to free-NAPL conditions  ---
!
            IF( (ASNTX-ASNX).LT.EPSL ) THEN
              PN(2,N) = PN(2,N) + 1.D-1
              NPHAZ(2,N) = 7
!
!---        Transition from trapped-NAPL to dissolved-oil
!           conditions  ---
!
            ELSEIF( ASNX.LT.EPSL ) THEN
!
!---          Transition from trapped-NAPL to dissolved-oil
!             conditions  ---
!
              IF( XMLOX.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 1
!
!---          Transition from trapped-NAPL to
!             zero-solubility conditions  ---
!
              ELSE
                NPHAZ(2,N) = 7
              ENDIF
!
!---        No transition from trapped-NAPL conditions  ---
!
            ELSE
              NPHAZ(2,N) = 3
            ENDIF
!
!---      No transition from trapped-gas conditions  ---
!
          ELSE
!
!---        Transition from trapped-NAPL to free-NAPL conditions  ---
!
            IF( (ASNTX-ASNX).LT.EPSL ) THEN
              PN(2,N) = PN(2,N) + 1.D-1
              NPHAZ(2,N) = 8
!
!---        Transition from trapped-NAPL to dissolved-oil
!           conditions  ---
!
            ELSEIF( ASNX.LT.EPSL ) THEN
!
!---          Transition from trapped-NAPL to dissolved-oil
!             conditions  ---
!
              IF( XMLOX.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 2
!
!---          Transition from trapped-NAPL to zero-solubility
!             conditions  ---
!
              ELSE
                NPHAZ(2,N) = 8
              ENDIF
!
!---        No transition from trapped-NAPL conditions  ---
!
            ELSE
              NPHAZ(2,N) = 4
            ENDIF
          ENDIF
!
!---    Aqueous-gas w/ dissolved-oil conditions
!       Unknowns:  PL,PG,XMLO
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
          ASGTMX = SCHR(17,IZN)
!
!---      No transition from free-gas conditions  ---
!
          IF( HGL.GT.ENPR ) THEN
!
!---        Dissolved-oil mole fraction exceeds solubility limit,
!           transition from dissolved-oil to free-NAPL conditions,
!           bump NAPL pressure by 0.1 Pa to create NAPL saturation  ---
!
            IF( XMLO(2,N).GT.XMLOX ) THEN
              NPHAZ(2,N) = 9
              PN(2,N) = PN(2,N) + 1.D-1
!
!---          Dissolved-oil mole fraction less than solubility limit,
!             no transition from dissolved-oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
              NPHAZ(2,N) = 5
!
!---        Transition from dissolved-oil to zero-solubility
!           conditions  ---
!
            ELSE
              NPHAZ(2,N) = 9
            ENDIF
!
!---      Transition from free-gas to trapped-gas conditions  ---
!
          ELSEIF( ASGTMX.GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
            PG(2,N) = PL(2,N) + ENPR/BGL
            PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        PL(2,N)+ENPR/BNL )
!
!---        Transition from dissolved-oil to free-NAPL conditions  ---
!
            IF( XMLO(2,N).GT.XMLOX ) THEN
              NPHAZ(2,N) = 8
              PN(2,N) = PN(2,N) + 1.D-1
!
!---        No transition from dissolved-oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
              NPHAZ(2,N) = 2
!
!---        Transition from dissolved-oil to zero-solubility conditions
!
            ELSE
              NPHAZ(2,N) = 8
            ENDIF
!
!---        Transition from free-gas to dissolved-air conditions,
!           set gas pressure to entry pressure and set NAPL pressure
!           to the critical pressure  ---
!
          ELSE
            PG(2,N) = PL(2,N) + ENPR/BGL
            PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        PL(2,N)+ENPR/BNL )
!
!---        Transition from free-gas to dissolved-air conditions  ---
!
            IF( PATM/HCAW.GT.EPSL ) THEN
!
!---          Transition from dissolved-oil to free-NAPL conditions  ---
!
              IF( XMLO(2,N).GT.XMLOX ) THEN
                PN(2,N) = PN(2,N) + 1.D-1
                NPHAZ(2,N) = 7
!
!---          No transition from dissolved-oil conditions  ---
!
              ELSEIF( XMLOX.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 1
!
!---          Transition from dissolved-oil to zero-solubility
!             conditions  ---
!
              ELSE
                NPHAZ(2,N) = 7
              ENDIF
            ENDIF
          ENDIF
!
!---    Aqueous-gas w/ trapped-NAPL conditions
!       Unknowns: PL,PG,SN
!
!       Kinetic volatilization option
!       Unknowns: PL,PG,SN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
!
!---      Set NAPL pressure to the critical pressure  ---
!
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
          ASNTMX = SCHR(11,IZN)*(SCHR(9,IZN)+TRPNL(2,N)*SCHR(8,IZN)/
     &      (BNL*SIG_NL))/(SCHR(9,IZN)+TRPNL(2,N))
          ASGTMX = SCHR(17,IZN)
!
!---      Scaled gas-aqueous capillary pressure greater than
!         entry pressure, transition from free-gas conditions  ---
!
          IF( HGL.GT.ENPR ) THEN
            IF( ASNTMX.GT.EPSL ) THEN
              INDX = 4
              ASLMINX = ASLMIN(1,N)
              ASTMINX = ASTMIN(1,N)
              ASTMAXX = ASTMAX(1,N)
              CALL KSP5( PG(2,N),PN(2,N),PL(2,N),SLX,SNX,SGX,
     &          RKLX,RKNX,RKGX,ASTX,ASLX,
     &          ASLMINX,ASTMINX,ASTMAXX,ASNTX,ASNRX,ASGTX,
     &          ASGTLX,ASGTNX,ASNTMX,ASNRMX,ASGTMX,
     &          INDX,IZN,N )
              PLND = 1.D+0/ASNTMX - 1.D+0
              ASNTX = (1.D+0-ASLMIN(1,N))/(1.D+0+PLND*
     &          (1.D+0-ASLMIN(1,N)))- (1.D+0-ASLX)/
     &          (1.D+0+PLND*(1.D+0-ASLX))
            ELSE
              ASNTX = 0.D+0
            ENDIF
            ASNX = SN(2,N)/(1.D+0-SCHR(4,IZN))
!
!---        Trapped-NAPL saturation exceeds maximum trapped-NAPL
!           saturation, transition to mobile-residual NAPL conditions,
!           bump NAPL pressure by 0.1 Pa to create NAPL saturation  ---
!
            IF( (ASNX-ASNTX)/(ASNTMX+EPSL).GT.1.D-1 ) THEN
              NPHAZ(2,N) = 9
              PN(2,N) = PN(2,N) + 1.D-1
!
!---        Transition from trapped-NAPL to dissolved-oil
!           conditions  ---
!
            ELSEIF( ASNX.LT.EPSL ) THEN
!
!---          Transition from trapped-NAPL to dissolved-oil
!             conditions  ---
!
              IF( XMLOX.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 5
!
!---          Transition from trapped-NAPL to zero-solubility
!             conditions  ---
!
              ELSE
                NPHAZ(2,N) = 9
              ENDIF
!
!---        No transition from trapped-NAPL conditions  ---
!
            ELSE
              NPHAZ(2,N) = 6
            ENDIF
!
!---      Transition from free-gas to trapped-gas conditions  ---
!
          ELSEIF( ASGTMX.GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
            PG(2,N) = PL(2,N) + ENPR/BGL
            PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        PL(2,N)+ENPR/BNL )
            IF( ASNTMX.GE.EPSL ) THEN
              PLND = 1.D+0/ASNTMX - 1.D+0
              ASNTX = (1.D+0-ASLMIN(1,N))/(1.D+0+PLND*
     &          (1.D+0-ASLMIN(1,N)))
            ELSE
              ASNTX = 0.D+0
            ENDIF
            ASNX = SN(2,N)/(1.D+0-SCHR(4,IZN))
!
!---        Transition from trapped-NAPL to free-NAPL conditions  ---
!
            IF( (ASNX-ASNTX)/(ASNTMX+EPSL).GT.1.D-1 ) THEN
              NPHAZ(2,N) = 8
              PN(2,N) = PN(2,N) + 1.D-1
!
!---        Transition from trapped-NAPL to dissolved-oil
!           conditions  ---
!
            ELSEIF( ASNX.LT.EPSL ) THEN
!
!---          Transition from trapped-NAPL to dissolved-oil
!             conditions  ---
!
              IF( XMLOX.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 2
!
!---          Transition from trapped-NAPL to zero-solubility
!             conditions  ---
!
              ELSE
                NPHAZ(2,N) = 8
              ENDIF
!
!---        No transition from trapped-NAPL conditions  ---
!
            ELSE
              NPHAZ(2,N) = 4
            ENDIF
!
!---      Scaled gas-aqueous capillary pressure less than entry
!         pressure, transition from free-gas to dissolved-air 
!         conditions  ---
!
          ELSE
!
!---        Set gas pressure to entry pressure and set NAPL pressure
!           to the critical pressure  ---
!
            PG(2,N) = PL(2,N) + ENPR/BGL
            PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        PL(2,N)+ENPR/BNL )
            IF( ASNTMX.GE.EPSL ) THEN
              PLND = 1.D+0/ASNTMX - 1.D+0
              ASNTX = (1.D+0-ASLMIN(1,N))/(1.D+0+PLND*
     &          (1.D+0-ASLMIN(1,N)))
            ELSE
              ASNTX = 0.D+0
            ENDIF
            ASNX = SN(2,N)/(1.D+0-SCHR(4,IZN))
!
!---        Trapped-NAPL saturation exceeds maximum trapped-NAPL
!           saturation, transition to mobile-residual NAPL conditions,
!           bump NAPL pressure by 0.1 Pa to create NAPL saturation  ---
!
            IF( (ASNX-ASNTX)/(ASNTMX+EPSL).GT.1.D-1 ) THEN
              NPHAZ(2,N) = 7
              PN(2,N) = PN(2,N) + 1.D-1
!
!---        Transition from trapped-NAPL to dissolved-oil
!           conditions  ---
!
            ELSEIF( ASNX.LT.EPSL ) THEN
!
!---          Transition from trapped-NAPL to dissolved-oil
!             conditions  ---
!
              IF( XMLOX.GT.1.D-12 ) THEN
                NPHAZ(2,N) = 1
!
!---          Transition from trapped-NAPL to zero-solubility
!             conditions  ---
!
              ELSE
                NPHAZ(2,N) = 7
              ENDIF
!
!---        No transition from trapped-NAPL conditions
!
            ELSE
              NPHAZ(2,N) = 3
            ENDIF
          ENDIF
!
!---    Aqueous-NAPL w/ dissolved-air conditions
!       Unknowns: PL,XMLA,PN
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.7 ) THEN
!
!---      Set gas pressure to entry pressure  ---
!
          PG(2,N) = PN(2,N)+ENPR/BGN
!
!---      Maximum effective trapped-NAPL saturation and
!         minimum effective aqueous saturation  ---
!
          ASNTMX = SCHR(11,IZN)*(SCHR(9,IZN)+TRPNL(2,N)*SCHR(8,IZN)/
     &      (BNL*SIG_NL))/(SCHR(9,IZN)+TRPNL(2,N))
!
!---      Dissolved-oil vapor pressure, dissolved-air partial 
!         pressure, water vapor pressure  ---
!
          PGX = PVW(2,N) + XMLA(2,N)*HCAW + PSO(2,N)
!
!---      Scaled NAPL-aqueous capillary pressure less than the 
!         entry pressure, transition from free-NAPL to trapped-NAPL or 
!         dissolved-oil conditions  ---
!
          IF( HNL.LE.ENPR ) THEN
!
!---        Set NAPL pressure to the critical point  ---
!
            PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        PL(2,N)+ENPR/BNL )
!
!---        Transition from free-NAPL to trapped-NAPL conditions  ---
!
            IF( ASNTMX.GT.EPSL .AND. SNT(2,N).GT.EPSL ) THEN
!
!---          Transition from dissolved-air to free-gas conditions,
!             bump gas and NAPL pressures by 0.1 Pa to create
!             gas and NAPL saturations  ---
!
              IF( PGX.GT.(PG(2,N)+PATM) ) THEN
                PG(2,N) = PG(2,N) + 1.D-1
                PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &            PL(2,N)+ENPR/BNL )
                NPHAZ(2,N) = 6
!
!---          No transition from dissolved-air conditions  ---
!
              ELSE
                NPHAZ(2,N) = 3
              ENDIF
!
!---        Transition from free-NAPL to dissolved-oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
!
!---          Transition from dissolved-air to free-gas conditions,
!             bump gas pressure by 0.1 Pa to create gas saturation  ---
!
              IF( PGX.GT.(PG(2,N)+PATM) ) THEN
                PG(2,N) = PG(2,N) + 1.D-1
                NPHAZ(2,N) = 5
!
!---          No transition from dissolved-air conditions  ---
!
              ELSE
                NPHAZ(2,N) = 1
              ENDIF
!
!---        Transition from free-NAPL to zero-solubility conditions  ---
!
            ELSE
!
!---          Transition from dissolved-air to free-gas conditions,
!             bump gas pressure by 0.1 Pa to create gas saturation  ---
!
              IF( PGX.GT.(PG(2,N)+PATM) ) THEN
                PG(2,N) = PG(2,N) + 1.D-1
                NPHAZ(2,N) = 9
!
!---          No transition from dissolved-air conditions  ---
!
              ELSE
                NPHAZ(2,N) = 7
              ENDIF
            ENDIF
!
!---      No transition from free-NAPL conditions  ---
!
          ELSE
!
!---        Transition from dissolved-air to free-gas conditions,
!           bump gas pressure by 0.1 Pa to create gas saturation  ---
!
            IF( PGX.GT.(PG(2,N)+PATM) ) THEN
              PG(2,N) = PG(2,N) + 1.D-1
              NPHAZ(2,N) = 9
!
!---        No transition from dissolved-air conditions  ---
!
            ELSE
              NPHAZ(2,N) = 7
            ENDIF
          ENDIF
!
!---    Aqueous-NAPL w/ trapped-gas conditions
!       Unknowns: PL,SG,PN
!
!       Kinetic volatilization option
!       Unknowns: PL,SG,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 ) THEN
          PG(2,N) = PN(2,N)+ENPR/BGN
          ASNTMX = SCHR(11,IZN)*(SCHR(9,IZN)+TRPNL(2,N)*SCHR(8,IZN)/
     &      (BNL*SIG_NL))/(SCHR(9,IZN)+TRPNL(2,N))
!
!---      Transition from free-NAPL to trapped-NAPL or
!         dissolved-oil conditions  ---
!
          IF( HNL.LE.ENPR ) THEN
            PG(2,N) = PL(2,N)+ENPR/BGL
            PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &        PL(2,N)+ENPR/BNL )
            ASGTMX = SCHR(17,IZN)
            IF( ASGTMX.GE.EPSL ) THEN
            PLND = 1.D+0/ASGTMX - 1.D+0
            ASGTX = (1.D+0-ASTMIN(1,N))/(1.D+0+PLND*(1.D+0-ASTMIN(1,N)))
            ELSE
              ASGTX = 0.D+0
            ENDIF
            ASGX = SG(2,N)/(1.D+0-SCHR(4,IZN))
!
!---        Transition from free-NAPL to trapped-NAPL conditions  ---
!
            IF( ASNTMX.GT.EPSL .AND. SNT(2,N).GT.EPSL ) THEN
!
!---          Transition from trapped-gas to free-gas conditions  ---
!
              IF( (ASGX-ASGTX)/(ASGTMX+EPSL).GT.1.D-1 ) THEN
                PG(2,N) = PG(2,N) + 1.D-1
                NPHAZ(2,N) = 6
!
!---          Transition from trapped-gas to dissolved-air
!             conditions  ---
!
              ELSEIF( ASGX.LT.EPSL ) THEN
                NPHAZ(2,N) = 3
!
!---          No transition from trapped-gas conditions  ---
!
              ELSE
                NPHAZ(2,N) = 4
              ENDIF
!
!---        Transition from free-NAPL to dissolved-oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
!
!---          Transition from trapped-gas to free-gas conditions  ---
!
              IF( (ASGX-ASGTX)/(ASGTMX+EPSL).GT.1.D-1 ) THEN
                PG(2,N) = PG(2,N) + 1.D-1
                NPHAZ(2,N) = 5
!
!---          Transition from trapped-gas to dissolved-air
!             conditions  ---
!
              ELSEIF( ASGX.LT.EPSL ) THEN
                NPHAZ(2,N) = 1
!
!---          No transition from trapped-gas conditions  ---
!
              ELSE
                NPHAZ(2,N) = 2
              ENDIF
!
!---        Transition from free-NAPL to zero-solubility conditions  ---
!
            ELSE
!
!---          Transition from trapped-gas to free-gas conditions  ---
!
              IF( (ASGX-ASGTX)/(ASGTMX+EPSL).GT.1.D-1 ) THEN
                PG(2,N) = PG(2,N) + 1.D-1
                NPHAZ(2,N) = 9
!
!---          Transition from trapped-gas to dissolved-air
!             conditions  ---
!
              ELSEIF( ASGX.LT.EPSL ) THEN
                NPHAZ(2,N) = 7
!
!---          No transition from trapped-gas conditions  ---
!
              ELSE
                NPHAZ(2,N) = 8
              ENDIF
            ENDIF
!
!---      No transition from free-NAPL conditions  ---
!
          ELSE
            IF( ASGTMX.GT.EPSL ) THEN
            PLND = 1.D+0/ASGTMX - 1.D+0
            ASGTX = (1.D+0-ASTMIN(1,N))/(1.D+0+PLND*(1.D+0-ASTMIN(1,N)))
            ELSE
              ASGTX = 0.D+0
            ENDIF
            ASGX = SG(2,N)/(1.D+0-SCHR(4,IZN))
!
!---        Transition from trapped-gas to free-gas conditions  ---
!
            IF( (ASGX-ASGTX)/(ASGTMX+EPSL).GT.1.D-1 ) THEN
              PG(2,N) = PG(2,N) + 1.D-1
              NPHAZ(2,N) = 9
!
!---        Transition from trapped-gas to dissolved-air conditions  ---
!
            ELSEIF( ASGX.LT.EPSL ) THEN
                NPHAZ(2,N) = 7
!
!---        No transition from trapped-gas conditions  ---
!
            ELSE
              NPHAZ(2,N) = 8
            ENDIF
          ENDIF
!
!---    Aqueous-NAPL-gas conditions
!       Unknowns: PL,PG,PN
!
!       Kinetic volatilization option
!       Unknowns: PL,PG,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.9 ) THEN
          ASNTMX = SCHR(11,IZN)*(SCHR(9,IZN)+TRPNL(2,N)*SCHR(8,IZN)/
     &      (BNL*SIG_NL))/(SCHR(9,IZN)+TRPNL(2,N))
          ASGTMX = SCHR(17,IZN)
!
!---      Transition from free-NAPL to trapped-NAPL or
!         dissolved-oil conditions  ---
!
          IF( HGN.GE.HNL .OR. HNL.LE.ENPR ) THEN
!
!---        Transition from free-NAPL to trapped-NAPL conditions  ---
!
            IF( ASNTMX.GT.EPSL .AND. SNT(2,N).GT.EPSL ) THEN
!
!---          No transition from free-gas conditions  ---
!
              IF( HGL.GT.ENPR ) THEN
                PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &            PL(2,N)+ENPR/BNL )
                NPHAZ(2,N) = 6
!
!---          Transition from free-gas to trapped-gas conditions  ---
!
              ELSEIF( ASGTMX.GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
                PG(2,N) = PL(2,N) + ENPR/BGL
                PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &            PL(2,N)+ENPR/BNL )
                NPHAZ(2,N) = 4
!
!---          Transition from free-gas to dissolved-air conditions  ---
!
              ELSE
                PG(2,N) = PL(2,N) + ENPR/BGL
                PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &            PL(2,N)+ENPR/BNL )
                NPHAZ(2,N) = 3
              ENDIF
!
!---        Transition from free-NAPL to dissolved-oil conditions  ---
!
            ELSEIF( XMLOX.GT.1.D-12 ) THEN
!
!---          No transition from free-gas conditions  ---
!
              IF( HGL.GT.ENPR ) THEN
                PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &            PL(2,N)+ENPR/BNL )
                NPHAZ(2,N) = 5
!
!---          Transition from free-gas to trapped-gas conditions  ---
!
              ELSEIF( ASGTMX.GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
                PG(2,N) = PL(2,N) + ENPR/BGL
                PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &            PL(2,N)+ENPR/BNL )
                NPHAZ(2,N) = 2
!
!---          Transition from free-gas to dissolved-air conditions  ---
!
              ELSE
                PG(2,N) = PL(2,N)+ENPR/BGL
                PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &            PL(2,N)+ENPR/BNL )
                NPHAZ(2,N) = 1
              ENDIF
!
!---        Transition from free-NAPL to no dissolved-oil
!           conditions  ---
!
            ELSE
!
!---          No transition from free-gas conditions  ---
!
              IF( HGL.GT.ENPR ) THEN
                NPHAZ(2,N) = 9
!
!---          Transition from free-gas to trapped-gas conditions  ---
!
              ELSEIF( ASGTMX.GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
                PG(2,N) = PN(2,N)+ENPR/BGN
                NPHAZ(2,N) = 8
!
!---          Transition from free-gas to dissolved-air conditions  ---
!
              ELSE
                PG(2,N) = PN(2,N)+ENPR/BGN
                NPHAZ(2,N) = 7
              ENDIF
            ENDIF
!
!---      No transition from free-NAPL conditions  ---
!
          ELSE
!
!---        No transition from free-gas conditions  ---
!
            IF( HGL.GT.ENPR ) THEN
              NPHAZ(2,N) = 9
!
!---        Transition from free-gas to trapped-gas conditions  ---
!
            ELSEIF( ASGTMX.GT.EPSL .AND. SGT(2,N).GT.EPSL ) THEN
              PG(2,N) = PN(2,N)+ENPR/BGN
              NPHAZ(2,N) = 8
!
!---        Transition from free-gas to dissolved-air conditions  ---
!
            ELSE
              PG(2,N) = PN(2,N)+ENPR/BGN
              NPHAZ(2,N) = 7
            ENDIF
          ENDIF
        ENDIF
  100 CONTINUE
!
!--- Load previous time step condition for initial time step  ---
!
      IF( NSTEP-NRST.EQ.0 ) THEN
        DO 110 N = 1,NFLD
          IF( IXP(N).EQ.0 ) GOTO 110
          NPHAZ(1,N) = NPHAZ(2,N)
  110   CONTINUE
      ENDIF
!
!--- Compute increments  ---
!
      DO 300 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 300
        IZN = IZ(N)
        XMLOX = MIN( XMLO(2,N),PSO(2,N)/HCOW )
        XMLWX = MAX( 1.D+0-XMLOX,0.D+0 )
        XLOX = XMLOX*WTMO/(XMLOX*WTMO + XMLWX*WTMW)
        PX = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
        CALL WATLQD( T(2,N),PX,RHOLX )
        SIG_GL = SCHR(5,IZN)
        SIG_NL = SCHR(6,IZN)
        SIG_GN = SCHR(7,IZN)
        RHOLOX = XLOX*RHOLX
        CALL SFT5( RHOLOX )
        BGL = SCHR(8,IZN)/SIG_GL
        BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
        BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
        ENPR = SCHR(2,IZN)*RHORL*GRAV
        XMLOX = PSO(2,N)/HCOW
!
!---    Aqueous w/ dissolved-oil, dissolved-air conditions
!       Unknowns: PL,XMLA,XMLO
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,PN,XMLO  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          PG(2,N) = PL(2,N) + ENPR/BGL
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
          PVAX = MAX( PG(2,N)+PATM-PSO(2,N)-PSW(2,N),0.D+0 )
          XMLAX = MAX(PVAX/HCAW,PATM/HCAW)
          DNR(IEQW,N) = MAX( 1.D-1,1.D-6*(PL(2,N)+PATM) )
!          DNR(IEQA,N) = SIGN( 1.D-5*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
!          IF( XMLO(2,N).LT.1.D-4*XMLOX ) THEN
!            DNR(IEQO,N) = 1.D-7*XMLOX
!          ELSE
!            DNR(IEQO,N) = SIGN( 1.D-6*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
!          ENDIF
          IF( XMLA(2,N).GT.(1.D-2*XMLAX) ) THEN
            DNR(IEQA,N) = SIGN( 1.D-5*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
          ELSE
            DNR(IEQA,N) = SIGN( 1.D-5*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
            DNR(IEQO,N) = 1.D-2
          ELSE
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
          ENDIF
!
!---    Aqueous w/ trapped-NAPL, dissolved-air conditions
!       Unknowns: PL,XMLA,SN
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,SN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          PG(2,N) = PL(2,N) + ENPR/BGL
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
          DNR(IEQW,N) = MAX( 1.D-1,1.D-6*(PL(2,N)+PATM) )
          DNR(IEQA,N) = SIGN( 1.D-5,5.D-1-SG(2,N) )
!          IF( XMLO(2,N).LT.1.D-4*XMLOX ) THEN
!            DNR(IEQO,N) = 1.D-7*XMLOX
!          ELSE
!            DNR(IEQO,N) = SIGN( 1.D-6*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
!          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
            DNR(IEQO,N) = 1.D-2
          ELSE
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
          ENDIF
!
!---    Aqueous w/ trapped-NAPL, dissolved-air conditions
!       Unknowns: PL,XMLA,SN
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,SN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
          PG(2,N) = PL(2,N) + ENPR/BGL
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
          DNR(IEQW,N) = -1.D-3
          PVAX = MAX( PG(2,N)+PATM-PSO(2,N)-PSW(2,N),ZERO )
          XMLAX = PVAX/HCAW
!          DNR(IEQA,N) = SIGN( 1.D-5*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
          IF( XMLA(2,N).GT.(1.D-2*XMLAX) ) THEN
            DNR(IEQA,N) = SIGN( 1.D-5*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
          ELSE
            DNR(IEQA,N) = SIGN( 1.D-5*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
            DNRO = MAX(1.D-9,1.D-4*SN(2,N))
            DNR(IEQO,N) = SIGN( DNRO,5.D-1-SN(2,N) )
          ELSE
            DNRO = MAX(1.D-9,1.D-4*SN(2,N))
            DNR(IEQO,N) = SIGN( DNRO,5.D-1-SN(2,N) )
          ENDIF
!
!---    Aqueous w/ trapped-NAPL, trapped-gas conditions
!       Unknowns: PL,SG,SN
!
!       Kinetic volatilization option
!       Unknowns: PL,SG,SN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
          PG(2,N) = PL(2,N) + ENPR/BGL
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
          DNR(IEQW,N) = -1.D-3
          DNR(IEQA,N) = SIGN( 1.D-5,5.D-1-SG(2,N) )
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
            DNR(IEQO,N) = SIGN( 1.D-5,5.D-1-SN(2,N) )
          ELSE
            DNR(IEQO,N) = SIGN( 1.D-5,5.D-1-SN(2,N) )
          ENDIF
!
!---    Aqueous-gas w/ dissolved-oil conditions
!       Unknowns:  PL,PG,XMLO
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
!          DNR(IEQW,N) = -MAX( 1.D-1,1.D-6*ABS(PG(2,N)-PL(2,N)) )
!          DNR(IEQA,N) = MAX( 1.D-1,1.D-6*ABS(PG(2,N)-PL(2,N)) )
!          IF( XMLO(2,N).LT.1.D-4*XMLOX ) THEN
!            DNR(IEQO,N) = 1.D-7*XMLOX
!          ELSE
!            DNR(IEQO,N) = SIGN( 1.D-6*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
!          ENDIF
          DNRX = MAX( 1.D-2,1.D-6*ABS(PG(2,N)-PL(2,N)) )
          DNR(IEQW,N) = SIGN( DNRX,(5.D-1-SL(2,N)) )
          DNR(IEQA,N) = -DNR(IEQW,N)
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
            DNR(IEQO,N) = 1.D-2
          ELSE
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
          ENDIF
!
!---    Aqueous-gas w/ trapped-NAPL conditions
!       Unknowns: PL,PG,SN
!
!       Kinetic volatilization option
!       Unknowns: PL,PG,SN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
          PN(2,N) = MAX( (BNL*PL(2,N)+BGN*PG(2,N))/(BNL+BGN),
     &      PL(2,N)+ENPR/BNL )
          DNRX = MAX( 1.D-2,1.D-6*ABS(PG(2,N)-PL(2,N)) )
          DNR(IEQW,N) = SIGN( DNRX,(5.D-1-SL(2,N)) )
          DNR(IEQA,N) = -DNR(IEQW,N)
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
            DNRO = MAX(1.D-9,1.D-4*SN(2,N))
            DNR(IEQO,N) = SIGN( DNRO,5.D-1-SN(2,N) )
          ELSE
            DNRO = MAX(1.D-9,1.D-4*SN(2,N))
            DNR(IEQO,N) = SIGN( DNRO,5.D-1-SN(2,N) )
          ENDIF
!
!---    Aqueous-NAPL w/ dissolved-air conditions
!       Unknowns: PL,XMLA,PN
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.7 ) THEN
          PG(2,N) = PN(2,N) + ENPR/BGN
          DNR(IEQW,N) = -MAX( 1.D-2,1.D-6*ABS(PN(2,N)-PL(2,N))  )
          IF( SN(2,N)/(1.D+0-SCHR(4,IZN)).GT.5.D-1 ) THEN
            DNR(IEQO,N) = -DNR(IEQO,N)
          ENDIF
          PVAX = MAX( PG(2,N)+PATM-PSO(2,N)-PSW(2,N),ZERO )
          XMLAX = PVAX/HCAW
!          DNR(IEQA,N) = SIGN( 1.D-5*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
          IF( XMLA(2,N).GT.(1.D-2*XMLAX) ) THEN
            DNR(IEQA,N) = SIGN( 1.D-5*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
          ELSE
            DNR(IEQA,N) = SIGN( 1.D-5*XMLAX,5.D-1*XMLAX-XMLA(2,N) )
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
            DNR(IEQO,N) = MAX( 1.D-2,1.D-5*ABS(PN(2,N)-PL(2,N))  )
            IF( SN(2,N)/(1.D+0-SCHR(4,IZN)).GT.5.D-1 ) THEN
              DNR(IEQO,N) = -DNR(IEQO,N)
            ENDIF
          ELSE
            DNR(IEQO,N) = MAX( 1.D-2,1.D-5*ABS(PN(2,N)-PL(2,N))  )
            IF( SN(2,N)/(1.D+0-SCHR(4,IZN)).GT.5.D-1 ) THEN
              DNR(IEQO,N) = -DNR(IEQO,N)
            ENDIF
          ENDIF
!
!---    Aqueous-NAPL w/ trapped-gas conditions
!       Unknowns: PL,SG,PN
!
!       Kinetic volatilization option
!       Unknowns: PL,SG,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 ) THEN
          PG(2,N) = PN(2,N) + ENPR/BGN
          DNR(IEQW,N) = -MAX( 1.D-2,1.D-6*ABS(PN(2,N)-PL(2,N))  )
          DNR(IEQO,N) = MAX( 1.D-2,1.D-6*ABS(PN(2,N)-PL(2,N))  )
          IF( SN(2,N)/(1.D+0-SCHR(4,IZN)).GT.5.D-1 ) THEN
            DNR(IEQW,N) = -DNR(IEQW,N)
            DNR(IEQO,N) = -DNR(IEQO,N)
          ENDIF
          DNR(IEQA,N) = SIGN( 1.D-5,5.D-1-SG(2,N) )
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
            DNR(IEQO,N) = MAX( 1.D-2,1.D-5*ABS(PN(2,N)-PL(2,N))  )
            IF( SN(2,N)/(1.D+0-SCHR(4,IZN)).GT.5.D-1 ) THEN
              DNR(IEQO,N) = -DNR(IEQO,N)
            ENDIF
          ELSE
            DNR(IEQO,N) = MAX( 1.D-2,1.D-5*ABS(PN(2,N)-PL(2,N))  )
            IF( SN(2,N)/(1.D+0-SCHR(4,IZN)).GT.5.D-1 ) THEN
              DNR(IEQO,N) = -DNR(IEQO,N)
            ENDIF
          ENDIF
!
!---    Aqueous-NAPL-gas conditions
!       Unknowns: PL,PG,PN
!
!       Kinetic volatilization option
!       Unknowns: PL,PG,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.9 ) THEN
          DNR(IEQW,N) = -MAX( 1.D-2,1.D-6*ABS(PN(2,N)-PL(2,N)) )
          DNR(IEQA,N) = MAX( 1.D-2,1.D-6*ABS(PG(2,N)-PN(2,N)) )
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IF( XMLO(2,N).GT.(1.D-2*XMLOX) ) THEN
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ELSE
              DNR(IEQDO,N) = SIGN( 1.D-5*XMLOX,5.D-1*XMLOX-XMLO(2,N) )
            ENDIF
            DNR(IEQO,N) = MAX( 1.D-2,1.D-5*ABS(PG(2,N)-PN(2,N)) )
          ELSE
            DNR(IEQO,N) = MAX( 1.D-2,1.D-5*ABS(PG(2,N)-PN(2,N)) )
          ENDIF
        ENDIF
!
!--- Increment the primary variables  ---
!
        DO 200 M = 2,ISVC+2
          PL(M,N) = PL(2,N)
          PN(M,N) = PN(2,N)
          PG(M,N) = PG(2,N)
          XMLO(M,N) = XMLO(2,N)
          XMLA(M,N) = XMLA(2,N)
          SN(M,N) = SN(2,N)
          SG(M,N) = SG(2,N)
!
!---      Aqueous w/ dissolved-oil, dissolved-air conditions
!         Unknowns: PL,XMLA,XMLO
!
!         Kinetic volatilization option
!         Unknowns: PL,XMLA,PN,XMLO  ---
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
            IF( M.EQ.IEQW+2 ) PL(M,N) = PL(M,N) + DNR(IEQW,N)
            IF( M.EQ.IEQA+2 ) XMLA(M,N) = XMLA(M,N) + DNR(IEQA,N)
            PG(M,N) = PL(M,N)+ENPR/BGL
            PN(M,N) = MAX( (BNL*PL(M,N)+BGN*PG(M,N))/(BNL+BGN),
     &        PL(M,N)+ENPR/BNL )
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              IF( M.EQ.IEQDO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQDO,N)
              IF( M.EQ.IEQO+2 ) PN(M,N) = PN(M,N) + DNR(IEQO,N)
            ELSE
              IF( M.EQ.IEQO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQO,N)
            ENDIF
!
!---      Aqueous w/ dissolved-oil, trapped-gas conditions
!         Unknowns: PL,SG,XMLO
!
!         Kinetic volatilization option
!         Unknowns: PL,SG,PN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
            IF( M.EQ.IEQW+2 ) PL(M,N) = PL(M,N) + DNR(IEQW,N)
            IF( M.EQ.IEQA+2 ) SG(M,N) = SG(M,N) + DNR(IEQA,N)
            PG(M,N) = PL(M,N)+ENPR/BGL
            PN(M,N) = MAX( (BNL*PL(M,N)+BGN*PG(M,N))/(BNL+BGN),
     &        PL(M,N)+ENPR/BNL )
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              IF( M.EQ.IEQDO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQDO,N)
              IF( M.EQ.IEQO+2 ) PN(M,N) = PN(M,N) + DNR(IEQO,N)
            ELSE
              IF( M.EQ.IEQO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQO,N)
            ENDIF
!
!---      Aqueous w/ trapped-NAPL, dissolved-air conditions
!         Unknowns: PL,XMLA,SN
!
!         Kinetic volatilization option
!         Unknowns: PL,XMLA,SN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
            IF( M.EQ.IEQW+2 ) PL(M,N) = PL(M,N) + DNR(IEQW,N)
            IF( M.EQ.IEQA+2 ) XMLA(M,N) = XMLA(M,N) + DNR(IEQA,N)
            PG(M,N) = PL(M,N)+ENPR/BGL
            PN(M,N) = MAX( (BNL*PL(M,N)+BGN*PG(M,N))/(BNL+BGN),
     &        PL(M,N)+ENPR/BNL )
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              IF( M.EQ.IEQDO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQDO,N)
              IF( M.EQ.IEQO+2 ) SN(M,N) = SN(M,N) + DNR(IEQO,N)
            ELSE
              IF( M.EQ.IEQO+2 ) SN(M,N) = SN(M,N) + DNR(IEQO,N)
            ENDIF
!
!---      Aqueous w/ trapped-NAPL, trapped-gas conditions
!         Unknowns: PL,SG,SN
!
!         Kinetic volatilization option
!         Unknowns: PL,SG,SN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
            IF( M.EQ.IEQW+2 ) PL(M,N) = PL(M,N) + DNR(IEQW,N)
            IF( M.EQ.IEQA+2 ) SG(M,N) = SG(M,N) + DNR(IEQA,N)
            IF( M.EQ.IEQO+2 ) SN(M,N) = SN(M,N) + DNR(IEQO,N)
            PG(M,N) = PL(M,N)+ENPR/BGL
            PN(M,N) = MAX( (BNL*PL(M,N)+BGN*PG(M,N))/(BNL+BGN),
     &        PL(M,N)+ENPR/BNL )
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              IF( M.EQ.IEQDO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQDO,N)
              IF( M.EQ.IEQO+2 ) SN(M,N) = SN(M,N) + DNR(IEQO,N)
            ELSE
              IF( M.EQ.IEQO+2 ) SN(M,N) = SN(M,N) + DNR(IEQO,N)
            ENDIF
!
!---      Aqueous-gas w/ dissolved-oil conditions
!         Unknowns:  PL,PG,XMLO
!
!         Kinetic volatilization option
!         Unknowns: PL,XMLA,PN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
            IF( M.EQ.IEQW+2 ) PL(M,N) = PL(M,N) + DNR(IEQW,N)
            IF( M.EQ.IEQA+2 ) PG(M,N) = PG(M,N) + DNR(IEQA,N)
            PN(M,N) = MAX( (BNL*PL(M,N)+BGN*PG(M,N))/(BNL+BGN),
     &        PL(M,N)+ENPR/BNL )
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              IF( M.EQ.IEQDO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQDO,N)
              IF( M.EQ.IEQO+2 ) PN(M,N) = PN(M,N) + DNR(IEQO,N)
            ELSE
              IF( M.EQ.IEQO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQO,N)
            ENDIF
!
!---      Aqueous-gas w/ trapped-NAPL conditions
!         Unknowns: PL,PG,SN
!
!         Kinetic volatilization option
!         Unknowns: PL,PG,SN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
            IF( M.EQ.IEQW+2 ) PL(M,N) = PL(M,N) + DNR(IEQW,N)
            IF( M.EQ.IEQA+2 ) PG(M,N) = PG(M,N) + DNR(IEQA,N)
            PN(M,N) = MAX( (BNL*PL(M,N)+BGN*PG(M,N))/(BNL+BGN),
     &        PL(M,N)+ENPR/BNL )
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              IF( M.EQ.IEQDO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQDO,N)
              IF( M.EQ.IEQO+2 ) SN(M,N) = SN(M,N) + DNR(IEQO,N)
            ELSE
              IF( M.EQ.IEQO+2 ) SN(M,N) = SN(M,N) + DNR(IEQO,N)
            ENDIF
!
!---      Aqueous-NAPL w/ dissolved-air conditions
!         Unknowns: PL,XMLA,PN
!
!         Kinetic volatilization option
!         Unknowns: PL,XMLA,PN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.7 ) THEN
            IF( M.EQ.IEQW+2 ) PL(M,N) = PL(M,N) + DNR(IEQW,N)
            IF( M.EQ.IEQA+2 ) XMLA(M,N) = XMLA(M,N) + DNR(IEQA,N)
            PG(M,N) = PN(M,N)+ENPR/BGN
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              IF( M.EQ.IEQDO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQDO,N)
              IF( M.EQ.IEQO+2 ) PN(M,N) = PN(M,N) + DNR(IEQO,N)
            ELSE
              IF( M.EQ.IEQO+2 ) PN(M,N) = PN(M,N) + DNR(IEQO,N)
            ENDIF
!
!---     Aqueous-NAPL w/ trapped-gas conditions
!         Unknowns: PL,SG,PN
!
!         Kinetic volatilization option
!         Unknowns: PL,SG,PN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.8 ) THEN
            IF( M.EQ.IEQW+2 ) PL(M,N) = PL(M,N) + DNR(IEQW,N)
            IF( M.EQ.IEQA+2 ) SG(M,N) = SG(M,N) + DNR(IEQA,N)
            IF( M.EQ.IEQO+2 ) PN(M,N) = PN(M,N) + DNR(IEQO,N)
            PG(M,N) = PN(M,N)+ENPR/BGN
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              IF( M.EQ.IEQDO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQDO,N)
              IF( M.EQ.IEQO+2 ) PN(M,N) = PN(M,N) + DNR(IEQO,N)
            ELSE
              IF( M.EQ.IEQO+2 ) PN(M,N) = PN(M,N) + DNR(IEQO,N)
            ENDIF
!
!---    Aqueous-NAPL-gas conditions
!       Unknowns: PL,PG,PN
!
!       Kinetic volatilization option
!       Unknowns: PL,PG,PN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.9 ) THEN
            IF( M.EQ.IEQW+2 ) PL(M,N) = PL(M,N) + DNR(IEQW,N)
            IF( M.EQ.IEQA+2 ) PG(M,N) = PG(M,N) + DNR(IEQA,N)
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              IF( M.EQ.IEQDO+2 ) XMLO(M,N) = XMLO(M,N) + DNR(IEQDO,N)
              IF( M.EQ.IEQO+2 ) PN(M,N) = PN(M,N) + DNR(IEQO,N)
            ELSE
              IF( M.EQ.IEQO+2 ) PN(M,N) = PN(M,N) + DNR(IEQO,N)
            ENDIF
          ENDIF
  200   CONTINUE
  300 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRM5 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KSP5( PGX,PNX,PLX,SLX,SNX,SGX,RKLX,RKNX,RKGX,
     &                 ASTX,ASLX,ASLMINX,ASTMINX,ASTMAXX,ASNTX,
     &                 ASNRX,ASGTX,ASGTLX,ASGTNX,ASNTMX,ASNRMX,ASGTMX,
     &                 INDX,IZN,NX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!
!     Compute the aqueous saturation from the oil/aqueous capillary
!     pressure, and compute the aqueous relative permeability from the
!     aqueous saturation.
!     Compute the total saturation from the air/oil capillary
!     pressure, and compute the napl relative permeability from the
!     napl saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1997.
!     Last Modified by MD White on October 16, 1997.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
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
!----------------------Data Statements---------------------------------!
!
      SAVE TOL
      DATA TOL /1.D-12/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/KSP5'
!
!---  Residual aqueous saturation  ---
!
      SLRX = SCHR(4,IZN)
!
!---  No entrapment van Genuchten saturation and
!     relative permeability functions  ---
!
      IF( ISCHR(IZN).EQ.31 .OR. ISCHR(IZN).EQ.33 ) THEN
        HNL = MAX( ZERO,BNL*(PNX-PLX)/RHORL/GRAV )
        HGN = MAX( ZERO,BGN*(PGX-PNX)/RHORL/GRAV )
        CN = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( IRPL(IZN).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(14,IZN)
        ENDIF
        ASTX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HGN)**CN))**CM
        ASLX = MIN(ASTX,(1.D+0/(1.D+0 + (SCHR(1,IZN)*HNL)**CN))**CM)
        ASNX = ASTX-ASLX
        IF( ASNX.LT.TOL ) ASNX = 0.D+0
        ASGX = 1.D+0-ASTX
        IF( ASGX.LT.TOL ) ASGX = 0.D+0
        ASNTX = 0.D+0
        ASGTX = 0.D+0
        ASGTLX = 0.D+0
        ASGTNX = 0.D+0
        SLX = ASLX*(1.D+0-SLRX) + SLRX
        SNX = ASNX*(1.D+0-SLRX)
        SGX = ASGX*(1.D+0-SLRX)
        IF( SLX.LT.TOL ) SLX = 0.D+0
        IF( SNX.LT.TOL ) SNX = 0.D+0
        IF( SGX.LT.TOL ) SGX = 0.D+0
!
!---    Assign effective minimum aqueous saturation and apparent
!       maximum total-liquid saturation for conditions with
!       mobile-residual NAPL  ---
!
        IF( (ASTX-ASLX)/EPSL.GT.EPSL ) THEN
          ASLMINX = MIN( ASLMINX,ASLX )
          ASTMAXX = MAX( ASTMAXX,ASTX )
        ENDIF
!
!---    Residual NAPL function  ---
!
        IF( ISCHR(IZN).EQ.33 ) THEN
!
!---      Constant residual NAPL option  ---
!
          IF( ISNR.EQ.7 ) THEN
            SNRMX = ASNRMX*(1.D+0-SLRX)
            SNRX = MIN( SNRMX,SNX )
            ASNRX = SNRX/(1.D+0-SLRX)
!
!---      Default residual NAPL option  ---
!
          ELSE
            IF( (ASTMAXX-ASLX).LT.TOL .OR. ASNRMX.LT.TOL ) THEN
              ASNRX = 0.D+0
            ELSE
              ASNRX = ESNR( ASLX,ASTX,ASTMAXX,ASNRMX )
              IF( ASNRX.GT.(ASTX-ASLX) ) THEN
                ASNRX = ASTX-ASLX
                ASTMAXX = ASTMX( ASLX,ASTX,ASNRX,ASNRMX )
              ENDIF
            ENDIF
          ENDIF
          ASNRX = MAX( ASNRX,0.D+0 )
!
!---    No residual NAPL Function  ---
!
        ELSE
          ASNRX = 0.D+0
        ENDIF
!
!---    Mobile NAPL  ---
!
        ASNMX = ASTX-ASLX-ASNRX
        IF( ASNMX.LT.TOL ) ASNMX = 0.D+0
!
!---    Constant aqueous relative permeability  ---
!
        IF( IRPL(IZN).EQ.0 ) THEN
          RKLX(1) = RPLC(2,IZN)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Mualem aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.1 ) THEN
          CM = RPLC(2,IZN)
          RKLX(1) = SQRT(ASLX)*(1.D+0 -
     &      (1.D+0-ASLX**(1.D+0/CM))**CM)**2
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Burdine aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.2 ) THEN
          CM = RPLC(2,IZN)
          RKLX(1) = (ASLX**2)*(1.D+0 -
     &      (1.D+0-ASLX**(1.D+0/CM))**CM)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Unrecognized aqueous relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized Aqueous Relative Permeability' //
     &      ' Function: IRPL(IZN) = '
          INDX = 12
          IMSG = IRPL(IZN)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Constant gas relative permeability  ---
!
        IF( MOD(IRPG(IZN),100).EQ.0 ) THEN
          RKGX = RPGC(3,IZN)
!
!---    Mualem gas relative permeability  ---
!
        ELSEIF( MOD(IRPG(IZN),100).EQ.1 ) THEN
          CM = RPGC(3,IZN)
          RKGX = SQRT(ASGX)*((1.D+0-ASTX**(1.D+0/CM))**CM)**2
!
!---    Burdine gas relative permeability  ---
!
        ELSEIF( MOD(IRPG(IZN),100).EQ.2 ) THEN
          CM = RPGC(3,IZN)
          RKGX = (ASGX**2)*((1.D+0-ASTX**(1.D+0/CM))**CM)
!
!---    Unrecognized gas relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized Gas Relative Permeability' //
     &      ' Function: IRPG(IZN) = '
          INDX = 12
          IMSG = IRPG(IZN)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Constant NAPL relative permeability  ---
!
        IF( IRPN(IZN).EQ.0 ) THEN
          RKNX = RPNC(1,IZN)
!
!---    Mualem NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.1 ) THEN
          CM = RPNC(1,IZN)
          ASX = MAX( MIN( (ASLX+ASNRX),1.D+0 ),0.D+0 )
          RKNX = SQRT(ASNMX)*((1.D+0 - ASX**(1.D+0/CM))**CM -
     &      (1.D+0 - (ASTX)**(1.D+0/CM))**CM)**2
!
!---    Burdine NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.2 ) THEN
          CM = RPNC(1,IZN)
          ASX = MAX( MIN( (ASLX+ASNRX),1.D+0 ),0.D+0 )
          RKNX = (ASNMX**2)*((1.D+0 - ASX**(1.D+0/CM))**CM -
     &      (1.D+0 - (ASTX)**(1.D+0/CM))**CM)
!
!---    Unrecognized NAPL relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized NAPL Relative Permeability' //
     &      ' Function: IRPN(IZN) = '
          INDX = 12
          IMSG = IRPN(IZN)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  No entrapment Brooks and Corey saturation and
!     relative permeability functions  ---
!
      ELSEIF( ISCHR(IZN).EQ.35 .OR. ISCHR(IZN).EQ.37 ) THEN
        HNL = MAX( ZERO,BNL*(PNX-PLX)/RHORL/GRAV )
        HGN = MAX( ZERO,BGN*(PGX-PNX)/RHORL/GRAV )
        CL = MAX( SCHR(3,IZN),SMALL )
        IF( HGN-SCHR(2,IZN).LE.EPSL ) THEN
          ASTX = 1.D+0
        ELSE
          ASTX = (SCHR(2,IZN)/HGN)**CL
        ENDIF
        IF( HNL-SCHR(2,IZN).LE.EPSL ) THEN
          ASLX = 1.D+0
        ELSE
          ASLX = (SCHR(2,IZN)/HNL)**CL
        ENDIF
        ASNX = ASTX-ASLX
        IF( ASTX-ASLX.LT.TOL ) ASNX = 0.D+0
        ASGX = 1.D+0-ASTX
        IF( 1.D+0-ASTX.LT.TOL ) ASGX = 0.D+0
        ASNTX = 0.D+0
        SLX = ASLX*(1.D+0-SLRX) + SLRX
        SNX = ASNX*(1.D+0-SLRX)
        SGX = ASGX*(1.D+0-SLRX)
        IF( SLX.LT.TOL ) SLX = 0.D+0
        IF( SNX.LT.TOL ) SNX = 0.D+0
        IF( SGX.LT.TOL ) SGX = 0.D+0
!
!---    Assign effective minimum aqueous saturation and apparent
!       maximum total-liquid saturation for conditions with
!       mobile-residual NAPL  ---
!
        IF( (ASTX-ASLX)/EPSL.GT.EPSL ) THEN
          ASLMINX = MIN( ASLMINX,ASLX )
          ASTMAXX = MAX( ASTMAXX,ASTX )
        ENDIF
!
!---    Residual NAPL function ---
!
        IF( ISCHR(IZN).EQ.37 ) THEN
!
!---      Constant residual NAPL option  ---
!
          IF( ISNR.EQ.7 ) THEN
            SNRMX = ASNRMX*(1.D+0-SLRX)
            SNRX = MIN( SNRMX,SNX )
            ASNRX = SNRX/(1.D+0-SLRX)
            ASNRKX = MAX( (SNX-SNRX)/(1.D+0-SLRX-SNRX),0.D+0 )
!
!---      Default residual NAPL option  ---
!
          ELSE
            IF( (ASTMAXX-ASLX).LT.TOL .OR. ASNRMX.LT.TOL ) THEN
              ASNRX = 0.D+0
            ELSE
              ASNRX = ESNR( ASLX,ASTX,ASTMAXX,ASNRMX )
              IF( ASNRX.GT.(ASTX-ASLX) ) THEN
                ASNRX = ASTX-ASLX
                ASTMAXX = ASTMX( ASLX,ASTX,ASNRX,ASNRMX )
              ENDIF
            ENDIF
          ENDIF
          ASNRX = MAX( ASNRX,0.D+0 )
!
!---    No residual NAPL Function  ---
!
        ELSE
          ASNRX = 0.D+0
        ENDIF
!
!---    Mobile NAPL  ---
!
        ASNMX = ASTX-ASLX-ASNRX
        IF( ASNMX.LT.TOL ) ASNMX = 0.D+0
!
!---    Constant aqueous relative permeability  ---
!
        IF( IRPL(IZN).EQ.0 ) THEN
          RKLX(1) = RPLC(2,IZN)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Mualem aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.1 ) THEN
          CL = RPLC(2,IZN)
          RKLX(1) = ASLX**(2.5D+0 + 2.D+0/CL)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Burdine aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.2 ) THEN
          CL = RPLC(2,IZN)
          RKLX(1) = ASLX**(3.0D+0 + 2.D+0/CL)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Unrecognized aqueous relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized Aqueous Relative Permeability' //
     &      ' Function: IRPL(IZN) = '
          INDX = 12
          IMSG = IRPL(IZN)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Constant gas relative permeability  ---
!
        IF( MOD(IRPG(IZN),100).EQ.0 ) THEN
          RKGX = RPGC(3,IZN)
!
!---    Mualem gas relative permeability  ---
!
        ELSEIF( MOD(IRPG(IZN),100).EQ.1 ) THEN
          CL = RPGC(3,IZN)
          RKGX = SQRT(ASGX)*(1.D+0-ASTX**(1.D+0+1.D+0/CL))**2
!
!---    Burdine gas relative permeability  ---
!
        ELSEIF( MOD(IRPG(IZN),100).EQ.2 ) THEN
          CL = RPGC(3,IZN)
          RKGX = (ASGX**2)*(1.D+0-ASTX**(1.D+0 + 2.D+0/CL))
!
!---    Unrecognized gas relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized Gas Relative Permeability' //
     &      ' Function: IRPG(IZN) = '
          INDX = 12
          IMSG = IRPG(IZN)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Constant NAPL relative permeability  ---
!
        IF( IRPN(IZN).EQ.0 ) THEN
          RKNX = RPNC(1,IZN)
!
!---    Mualem NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.1 ) THEN
          CL = RPNC(1,IZN)
          ASX = MAX( MIN( (ASLX+ASNRX),1.D+0 ),0.D+0 )
          RKNX = SQRT(ASNMX)*( ASTX**(2.0D+0 + 2.D+0/CL) -
     &        ASX**(2.0D+0 + 2.D+0/CL) )
!
!---    Burdine NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.2 ) THEN
          CL = RPNC(1,IZN)
          ASX = MAX( MIN( (ASLX+ASNRX),1.D+0 ),0.D+0 )
          RKNX = (ASNMX**2)*( ASTX**(1.0D+0 + 2.D+0/CL) -
     &        ASX**(1.0D+0 + 2.D+0/CL) )
!
!---    Unrecognized NAPL relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized NAPL Relative Permeability' //
     &      ' Function: IRPN(IZN) = '
          INDX = 12
          IMSG = IRPN(IZN)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Entrapment van Genuchten saturation and
!     relative permeability functions  ---
!
      ELSEIF( ISCHR(IZN).EQ.32 .OR. ISCHR(IZN).EQ.34 ) THEN
        HNL = MAX( ZERO,BNL*(PNX-PLX)/RHORL/GRAV )
        HGN = MAX( ZERO,BGN*(PGX-PNX)/RHORL/GRAV )
        CN = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( IRPL(IZN).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(14,IZN)
        ENDIF
        ASTX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HGN)**CN))**CM
        ASLX = MIN(ASTX,(1.D+0/(1.D+0 + (SCHR(1,IZN)*HNL)**CN))**CM)
        ASTX = MAX(ASTX,0.D+0)
        ASLX = MAX(ASLX,0.D+0)
!
!---    No entrapment, used for initial condition checks and
!       boundary conditions  ---
!
        IF( INDX.EQ.2 ) THEN
          ASNTX = 0.D+0
          ASNFX = ASTX-ASLX
          IF( ASTX-ASLX.LT.TOL ) ASNFX = 0.D+0
          ASGTX = 0.D+0
          ASGTLX = 0.D+0
          ASGTNX = 0.D+0
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
!
!---    Dissolved oil condition, assume no trapped NAPL,
!       no residual NAPL, no free NAPL,
!       zero NAPL relative permeability,
!       entrapment of gas by aqueous, and
!       no entrapment of gas by NAPL  ---
!
        ELSEIF( INDX.EQ.3 ) THEN
          ASNTX = 0.D+0
          ASNFX = 0.D+0
          ASTMNX = MIN( ASTX,ASTMINX )
          IF( (ASTX-ASTMNX).LT.EPSL .OR. ASGTMX.LT.EPSL ) THEN
            ASGTLX = 0.D+0
          ELSE
            PLND = 1.D+0/ASGTMX - 1.D+0
            SGRLX = (1.D+0-ASTMNX)/(1.D+0 + PLND*(1.D+0-ASTMNX))
            ASGTLX = SGRLX-(1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ENDIF
          ASGTNX = 0.D+0
          ASGTX = ASGTLX+ASGTNX
          IF( ASGTX.LT.EPSL ) ASGTX = 0.D+0
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.EPSL ) ASGFX = 0.D+0
          ASTMINX = MIN( ASTMINX,ASTX )
!
!---    Entrapped NAPL phase condition, compute an estimate
!       of the NAPL pressure, assume no residual NAPL,
!       zero NAPL relative permeability, 
!       entrapment of gas by aqueous, and 
!       no entrapment of gas by NAPL  ---
!
        ELSEIF( INDX.EQ.4 ) THEN
          ASNX = SNX/(1.D+0-SLRX)
          CALL CAP5( ASLMINX,PGX,PLX,PNX,ASNTMX,ASNX,INDX,IZN,NX )
          IF( INDX.LT.0 ) THEN
            ISUB_LOG = ISUB_LOG-1
            RETURN
          ENDIF
          HNLZ = MAX( ZERO,BNL*(PNX-PLX)/RHORL/GRAV )
          HGNZ = MAX( ZERO,BGN*(PGX-PNX)/RHORL/GRAV )
          ASTX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HGNZ)**CN))**CM
          ASLX = MIN(ASTX,(1.D+0/(1.D+0 + (SCHR(1,IZN)*HNLZ)**CN))**CM)
          ASLMINX = MIN( ASLX,ASLMINX )
          PLND = 1.D+0/ASNTMX - 1.D+0
          ASNTX = (1.D+0-ASLMINX)/(1.D+0+PLND*(1.D+0-ASLMINX))
     &    - (1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ASNTX = MIN( ASNX,ASNTX )
          ASNFX = ASNX - ASNTX
          ASTMNX = MIN( ASTX,ASTMINX )
          IF( (ASTX-ASTMNX).LT.EPSL .OR. ASGTMX.LT.EPSL ) THEN
            ASGTLX = 0.D+0
          ELSE
            PLND = 1.D+0/ASGTMX - 1.D+0
            SGRLX = (1.D+0-ASTMNX)/(1.D+0 + PLND*(1.D+0-ASTMNX))
            ASGTLX = SGRLX-(1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ENDIF
          ASGTNX = 0.D+0
          ASGTX = ASGTLX+ASGTNX
          IF( ASGTX.LT.TOL ) ASGTX = 0.D+0
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
          ASTMINX = MIN( ASTX,ASTMINX )
!
!---    Full entrapment, NAPL by aqueous, gas by aqueous and
!       gas by NAPL  ---
!
        ELSE
          IF( (ASLX-ASLMINX).LT.TOL .OR. ASNTMX.LT.TOL ) THEN
            ASNTX = 0.D+0
          ELSE
            PLND = 1.D+0/ASNTMX - 1.D+0
            ASNTX = (1.D+0-ASLMINX)/(1.D+0 + PLND*(1.D+0-ASLMINX))
            ASNTX = ASNTX-(1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ENDIF
          ASTMNX = MIN( ASTX,ASTMINX )
          IF( (ASLX-ASTMNX).LT.TOL .OR. ASGTMX.LT.TOL ) THEN
            ASGTLX = 0.D+0
          ELSE
            PLND = 1.D+0/ASGTMX - 1.D+0
            SGRLX = (1.D+0-ASTMNX)/(1.D+0 + PLND*(1.D+0-ASTMNX))
            ASGTLX = SGRLX-(1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ENDIF
          ASLMNX = MAX( ASTMINX,ASLX )
          IF( (ASTX-ASLMNX).LT.TOL .OR. ASGTMX.LT.TOL ) THEN
            ASGTNX = 0.D+0
          ELSE
            PLND = 1.D+0/ASGTMX - 1.D+0
            SGRNX = (1.D+0-ASLMNX)/(1.D+0 + PLND*(1.D+0-ASLMNX))
            ASGTNX = SGRNX-(1.D+0-ASTX)/(1.D+0+PLND*(1.D+0-ASTX))
          ENDIF
          ASGTX = ASGTLX+ASGTNX
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
          ASNFX = ASTX-ASLX-ASGTNX
          IF( ASNFX.LT.TOL ) ASNFX = 0.D+0
        ENDIF
        ESLX = ASLX-ASNTX-ASGTLX
        SMX = SCHR(4,IZN)
        SLX = (ESLX)*(1.D+0-SMX) + SMX
        SNX = (ASNFX+ASNTX)*(1.D+0-SMX)
        SGX = (ASGFX+ASGTX)*(1.D+0-SMX)
        IF( SLX.LT.TOL ) SLX = 0.D+0
        IF( SNX.LT.TOL ) SNX = 0.D+0
        IF( SGX.LT.TOL ) SGX = 0.D+0
!
!---    Assign effective minimum aqueous saturation and apparent
!       maximum total-liquid saturation for conditions with
!       mobile-residual NAPL  ---
!
        IF( (ASTX-ASLX)/EPSL.GT.EPSL ) THEN
          ASLMINX = MIN( ASLMINX,ASLX )
          ASTMAXX = MAX( ASTMAXX,ASTX )
        ENDIF
!
!---    Residual NAPL function  ---
!
        IF( ISCHR(IZN).EQ.34 ) THEN
!
!---      Constant residual NAPL option  ---
!
          IF( ISNR.EQ.7 ) THEN
            SNRMX = ASNRMX*(1.D+0-SLRX)
            SNRX = MIN( SNRMX,SNX )
            ASNRX = SNRX/(1.D+0-SLRX)
!
!---      Default residual NAPL option  ---
!
          ELSE
            IF( (ASTMAXX-ASLX).LT.TOL .OR. ASNRMX.LT.TOL ) THEN
              ASNRX = 0.D+0
            ELSE
              ASNRX = ESNR( ASLX,ASTX,ASTMAXX,ASNRMX )
              IF( ASNRX.GT.(ASTX-ASLX) ) THEN
                ASNRX = ASTX-ASLX
                ASTMAXX = ASTMX( ASLX,ASTX,ASNRX,ASNRMX )
              ENDIF
            ENDIF
          ENDIF
          ASNRX = MAX( ASNRX,0.D+0 )
!
!---    No residual NAPL Function  ---
!
        ELSE
          ASNRX = 0.D+0
        ENDIF
!
!---    Mobile NAPL  ---
!
        ASNMX = ASTX-ASLX-ASNRX
        IF( ASNMX.LT.TOL ) ASNMX = 0.D+0
!
!---    Constant aqueous relative permeability  ---
!
        IF( IRPL(IZN).EQ.0 ) THEN
          RKLX(1) = RPLC(2,IZN)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Mualem aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.1 ) THEN
          CM = RPLC(2,IZN)
          RKLX(1) = SQRT(ESLX)*(1.D+0 -
     &      (1.D+0-ESLX**(1.D+0/CM))**CM)**2
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Burdine aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.2 ) THEN
          CM = RPLC(2,IZN)
          RKLX(1) = (ESLX**2)*(1.D+0 -
     &      (1.D+0-ESLX**(1.D+0/CM))**CM)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Unrecognized aqueous relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized Aqueous Relative Permeability' //
     &      ' Function: IRPL(IZN) = '
          INDX = 12
          IMSG = IRPL(IZN)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Constant gas relative permeability  ---
!
        IF( MOD(IRPG(IZN),100).EQ.0 ) THEN
          RKGX = RPGC(3,IZN)
!
!---    Mualem gas relative permeability  ---
!
        ELSEIF( MOD(IRPG(IZN),100).EQ.1 ) THEN
          CM = RPGC(3,IZN)
          RKGX = SQRT(ASGFX)*((1.D+0-ASTX**(1.D+0/CM))**CM)**2
!
!---    Burdine gas relative permeability  ---
!
        ELSEIF( MOD(IRPG(IZN),100).EQ.2 ) THEN
          CM = RPGC(3,IZN)
          RKGX = (ASGFX**2)*((1.D+0-ASTX**(1.D+0/CM))**CM)
!
!---    Unrecognized gas relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized Gas Relative Permeability' //
     &      ' Function: IRPG(IZN) = '
          INDX = 12
          IMSG = IRPG(IZN)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Constant NAPL relative permeability  ---
!
        IF( IRPN(IZN).EQ.0 ) THEN
          RKNX = RPNC(1,IZN)
!
!---    Mualem NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.1 ) THEN
          CM = RPNC(1,IZN)
          ASX = MAX( MIN( (ASLX+ASNRX),1.D+0 ),0.D+0 )
          RKNX = SQRT(ASNMX)*((1.D+0 - ASX**(1.D+0/CM))**CM -
     &      (1.D+0 - (ASTX)**(1.D+0/CM))**CM)**2
!
!---    Burdine NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.2 ) THEN
          CM = RPNC(1,IZN)
          ASX = MAX( MIN( (ASLX+ASNRX),1.D+0 ),0.D+0 )
          RKNX = (ASNMX**2)*((1.D+0 - ASX**(1.D+0/CM))**CM -
     &      (1.D+0 - (ASTX)**(1.D+0/CM))**CM)
!
!---    Unrecognized NAPL relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized NAPL Relative Permeability' //
     &      ' Function: IRPL(IZN) = '
          INDX = 12
          IMSG = IRPN(IZN)
          CALL WRMSGS( INDX )
        ENDIF
        IF( INDX.EQ.4 ) RKNX = 0.D+0
!
!---  Entrapment Brooks and Corey saturation and
!     relative permeability functions  ---
!
      ELSEIF( ISCHR(IZN).EQ.36 .OR. ISCHR(IZN).EQ.38 ) THEN
        HNL = MAX( ZERO,BNL*(PNX-PLX)/RHORL/GRAV )
        HGN = MAX( ZERO,BGN*(PGX-PNX)/RHORL/GRAV )
        CL = MAX( SCHR(3,IZN),SMALL )
        IF( HGN-SCHR(2,IZN).LT.EPSL ) THEN
          ASTX = 1.D+0
        ELSE
          ASTX = (SCHR(2,IZN)/HGN)**CL
        ENDIF
        IF( HNL-SCHR(2,IZN).LT.EPSL ) THEN
          ASLX = 1.D+0
        ELSE
          ASLX = (SCHR(2,IZN)/HNL)**CL
        ENDIF
!
!---    No entrapment, used for initial condition checks and
!       boundary conditions  ---
!
        IF( INDX.EQ.2 ) THEN
          ASNTX = 0.D+0
          ASNFX = ASTX-ASLX
          IF( ASTX-ASLX.LT.TOL ) ASNFX = 0.D+0
          ASGTX = 0.D+0
          ASGTLX = 0.D+0
          ASGTNX = 0.D+0
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
!
!---    Dissolved oil condition, assume no trapped NAPL,
!       no residual NAPL, no free NAPL,
!       zero NAPL relative permeability,
!       entrapment of gas by aqueous, and
!       no entrapment of gas by NAPL  ---
!
        ELSEIF( INDX.EQ.3 ) THEN
          ASNTX = 0.D+0
          ASNFX = 0.D+0
          ASTMNX = MIN( ASTX,ASTMINX )
          IF( (ASTX-ASTMNX).LT.TOL .OR. ASGTMX.LT.TOL ) THEN
            ASGTLX = 0.D+0
          ELSE
            PLND = 1.D+0/ASGTMX - 1.D+0
            SGRLX = (1.D+0-ASTMNX)/(1.D+0 + PLND*(1.D+0-ASTMNX))
            ASGTLX = SGRLX-(1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ENDIF
          ASGTNX = 0.D+0
          ASGTX = ASGTLX+ASGTNX
          IF( ASGTX.LT.TOL ) ASGTX = 0.D+0
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
          ASLMINX = ASLX
          ASTMINX = ASTX
!
!---    Entrapped NAPL phase condition, compute an estimate
!       of the NAPL pressure, assume no residual NAPL,
!       zero NAPL relative permeability, 
!       entrapment of gas by aqueous, and 
!       no entrapment of gas by NAPL  ---
!
        ELSEIF( INDX.EQ.4 ) THEN
          SMX = SCHR(4,IZN)
          ASNX = SNX/(1.D+0-SMX)
          CALL CAP5( ASLMINX,PGX,PLX,PNX,ASNTMX,ASNX,INDX,IZN,NX )
          IF( INDX.LT.0 ) THEN
            ISUB_LOG = ISUB_LOG-1
            RETURN
          ENDIF
          HNLZ = MAX( ZERO,BNL*(PNX-PLX)/RHORL/GRAV )
          HGNZ = MAX( ZERO,BGN*(PGX-PNX)/RHORL/GRAV )
          IF( HGNZ-SCHR(2,IZN).LE.EPSL ) THEN
            ASTX = 1.D+0
          ELSE
            ASTX = (SCHR(2,IZN)/HGNZ)**CL
          ENDIF
          IF( HNLZ-SCHR(2,IZN).LE.EPSL ) THEN
            ASLX = 1.D+0
          ELSE
            ASLX = (SCHR(2,IZN)/HNLZ)**CL
          ENDIF
          ASLX = MIN( ASTX,ASLX )
          ASLMINX = MIN( ASLX,ASLMINX )
          PLND = 1.D+0/ASNTMX - 1.D+0
          ASNTX = (1.D+0-ASLMINX)/(1.D+0+PLND*(1.D+0-ASLMINX))
     &    - (1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ASNTX = MIN( ASNX,ASNTX )
          ASNFX = ASNX - ASNTX
          ASTMNX = MIN( ASTX,ASTMINX )
          IF( (ASTX-ASTMNX).LT.TOL .OR. ASGTMX.LT.TOL ) THEN
            ASGTLX = 0.D+0
          ELSE
            PLND = 1.D+0/ASGTMX - 1.D+0
            SGRLX = (1.D+0-ASTMNX)/(1.D+0 + PLND*(1.D+0-ASTMNX))
            ASGTLX = SGRLX-(1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ENDIF
          ASGTNX = 0.D+0
          ASGTX = ASGTLX+ASGTNX
          IF( ASGTX.LT.TOL ) ASGTX = 0.D+0
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
          ASTMINX = ASTX
!
!---  Full entrapment, NAPL by aqueous, gas by aqueous and
!     gas by NAPL  ---
!
        ELSE
          IF( (ASLX-ASLMINX).LT.TOL .OR. ASNTMX.LT.TOL ) THEN
            ASNTX = 0.D+0
          ELSE
            PLND = 1.D+0/ASNTMX - 1.D+0
            SNRX = (1.D+0-ASLMINX)/(1.D+0 + PLND*(1.D+0-ASLMINX))
            ASNTX = SNRX-(1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ENDIF
          ASTMNX = MIN( ASTX,ASTMINX )
          IF( (ASLX-ASTMNX).LT.TOL .OR. ASGTMX.LT.TOL ) THEN
            ASGTLX = 0.D+0
          ELSE
            PLND = 1.D+0/ASGTMX - 1.D+0
            SGRLX = (1.D+0-ASTMNX)/(1.D+0 + PLND*(1.D+0-ASTMNX))
            ASGTLX = SGRLX-(1.D+0-ASLX)/(1.D+0+PLND*(1.D+0-ASLX))
          ENDIF
          ASLMNX = MAX( ASTMINX,ASLX )
          IF( (ASTX-ASLMNX).LT.TOL .OR. ASGTMX.LT.TOL ) THEN
            ASGTNX = 0.D+0
          ELSE
            PLND = 1.D+0/ASGTMX - 1.D+0
            SGRNX = (1.D+0-ASLMNX)/(1.D+0 + PLND*(1.D+0-ASLMNX))
            ASGTNX = SGRNX-(1.D+0-ASTX)/(1.D+0+PLND*(1.D+0-ASTX))
          ENDIF
          ASGTX = ASGTLX+ASGTNX
          ASGFX = 1.D+0-ASTX
          IF( ASGFX.LT.TOL ) ASGFX = 0.D+0
          ASNFX = ASTX-ASLX-ASGTNX
          IF( ASNFX.LT.TOL ) ASNFX = 0.D+0
        ENDIF
        ESLX = ASLX-ASNTX-ASGTLX
        SMX = SCHR(4,IZN)
        SLX = (ESLX)*(1.D+0-SMX) + SMX
        SNX = (ASNFX+ASNTX)*(1.D+0-SMX)
        SGX = (ASGFX+ASGTX)*(1.D+0-SMX)
        IF( SLX.LT.TOL ) SLX = 0.D+0
        IF( SNX.LT.TOL ) SNX = 0.D+0
        IF( SGX.LT.TOL ) SGX = 0.D+0
!
!---    Assign effective minimum aqueous saturation and apparent
!       maximum total-liquid saturation for conditions with
!       mobile-residual NAPL  ---
!
        IF( (ASTX-ASLX)/EPSL.GT.EPSL ) THEN
          ASLMINX = MIN( ASLMINX,ASLX )
          ASTMAXX = MAX( ASTMAXX,ASTX )
        ENDIF
!
!---    Residual NAPL function  ---
!
        IF( ISCHR(IZN).EQ.38 ) THEN
!
!---      Constant residual NAPL option  ---
!
          IF( ISNR.EQ.7 ) THEN
            SNRMX = ASNRMX*(1.D+0-SLRX)
            SNRX = MIN( SNRMX,SNX )
            ASNRX = SNRX/(1.D+0-SLRX)
!
!---      Default residual NAPL option  ---
!
          ELSE
            IF( (ASTMAXX-ASLX).LT.TOL .OR. ASNRMX.LT.TOL ) THEN
              ASNRX = 0.D+0
            ELSE
              ASNRX = ESNR( ASLX,ASTX,ASTMAXX,ASNRMX )
              IF( ASNRX.GT.(ASTX-ASLX) ) THEN
                ASNRX = ASTX-ASLX
                ASTMAXX = ASTMX( ASLX,ASTX,ASNRX,ASNRMX )
              ENDIF
            ENDIF
          ENDIF
          ASNRX = MAX( ASNRX,0.D+0 )
!
!---    No residual NAPL Function  ---
!
        ELSE
          ASNRX = 0.D+0
        ENDIF
!
!---    Mobile NAPL  ---
!
        ASNMX = ASTX-ASLX-ASNRX
        IF( ASNMX.LT.TOL ) ASNMX = 0.D+0
!
!---    Constant aqueous relative permeability  ---
!
        CL = RPLC(2,IZN)
        IF( IRPL(IZN).EQ.0 ) THEN
          RKLX(1) = CL
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Mualem aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.1 ) THEN
          RKLX(1) = ESLX**(2.5D+0 + 2.D+0/CL)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Burdine aqueous relative permeability  ---
!
        ELSEIF( IRPL(IZN).EQ.2 ) THEN
          RKLX(1) = ESLX**(3.0D+0 + 2.D+0/CL)
          RKLX(2) = RKLX(1)
          RKLX(3) = RKLX(1)
!
!---    Unrecognized aqueous relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized Aqueous Relative Permeability' //
     &      ' Function: IRPL(IZN) = '
          INDX = 12
          IMSG = IRPL(IZN)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Constant gas relative permeability  ---
!
        CL = RPGC(3,IZN)
        IF( MOD(IRPG(IZN),100).EQ.0 ) THEN
          RKGX = CL
!
!---    Mualem gas relative permeability  ---
!
        ELSEIF( MOD(IRPG(IZN),100).EQ.1 ) THEN
          RKGX = SQRT(ASGFX)*( 1.D+0 - ASTX**(2.D+0 + 2.D+0/CL) )
!
!---    Burdine gas relative permeability  ---
!
        ELSEIF( MOD(IRPG(IZN),100).EQ.2 ) THEN
          RKGX = (ASGFX**2)*(1.D+0 - ASTX**(1.D+0 + 2.D+0/CL))
!
!---    Unrecognized gas relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized Gas Relative Permeability' //
     &      ' Function: IRPG(IZN) = '
          INDX = 12
          IMSG = IRPG(IZN)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Constant NAPL relative permeability  ---
!
        CL = RPNC(1,IZN)
        IF( IRPN(IZN).EQ.0 ) THEN
          RKNX = CL
!
!---    Mualem NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.1 ) THEN
          RKNX = SQRT(ASNMX)*( ASTX**(2.0D+0 + 2.D+0/CL) -
     &        ASLX**(2.0D+0 + 2.D+0/CL) )
!
!---    Burdine NAPL relative permeability  ---
!
        ELSEIF( IRPN(IZN).EQ.2 ) THEN
          RKNX = (ASNMX**2)*( ASTX**(1.0D+0 + 2.D+0/CL) -
     &        ASLX**(1.0D+0 + 2.D+0/CL) )
!
!---    Unrecognized NAPL relative permeability  ---
!
        ELSE
          CHMSG = 'Unrecognized NAPL Relative Permeability' //
     &      ' Function: IRPL(IZN) = '
          INDX = 12
          IMSG = IRPN(IZN)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Klinkenberg effect  ---
!
      IF( IRPG(IZN).GE.100 .AND. IRPG(IZN).LT.200 ) THEN
        RKGX = RKGX*(1.D+0 + RPGC(5,IZN)*((PGX+PATM)**
     &    (RPGC(6,IZN)-1.D+0)))
      ENDIF
      RKLX(1) = MAX( RKLX(1),1.D-6 )
      RKLX(2) = MAX( RKLX(2),1.D-6 )
      RKLX(3) = MAX( RKLX(3),1.D-6 )
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KSP5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KV5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute kinetic volatilization rate (kg/s), where mass transfer
!     of oil between phases is defined as positive from NAPL-oil
!     to dissolved-oil.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 21 May 2008.
!     Last Modified by MD White, PNNL, 21 May 2008.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE NAPL
      USE JACOB
      USE GRID
      USE FDVP
      USE FDVG
      USE FDVD
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
      SUB_LOG(ISUB_LOG) = '/KV5'
!
!---  Loop over nodes and secondary variable indices  ---
!
      DO 300 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 300
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX+1
        NQY = NPY+IFLD
        NQZ = NPZ+IJFLD
        IZN = IZ(N)
        D50X = D50(IZN)
        DO 200 M = 2,ISVC+2
!
!---      No NAPL conditions  ---
!
          IF( SN(2,N).LT.EPSL ) THEN
            TMBP_O(M,N) = 0.D+0
            GOTO 200
          ENDIF
!
!---      Average gas velocity, cm/s  ---
!
          ZGX = 1.D+2*AVPVG(N)
!
!---      Limit gas velocity to 1D-8 cm/s  ---
!
          ZGX = MAX( ZGX,1.0D-8 )
!
!---      Oil gas diffusion coefficient, cm^2/s  ---
!
          DFGOX = 1.D+4*DFGO(M,N)
!
!---      Limit gas Peclet number to 0.02  ---
!
!          ZGX = MAX( ZGX,2.D-2*DFGOX/(1.D+2*D50X) )
!
!---      Kinetic volatilization coeff. w/ Wilkins et al. model  ---
!
          IF( ISLC(46).EQ.1 ) THEN
            VOLKX = 3.80189D-1*(ZGX**6.2D-1)*
     &        (DFGOX**3.8D-1)*((1.D+2*D50X)**(4.4D-1))
!
!---      Kinetic volatilization coeff. w/ Van der Ham and 
!         Brouwers model  ---
!
          ELSEIF( ISLC(46).EQ.2 ) THEN
            VOLKX = 2.17771D-1*(ZGX**8.8D-1)*
     &        (DFGOX**1.2D-1)*((1.D+2*D50X)**(7.0D-1))
!
!---      Kinetic volatilization coeff. w/ Yoon et al. model  ---
!
          ELSEIF( ISLC(46).EQ.3 ) THEN
            VOLKX = 2.63027D-1*(ZGX**6.8D-1)*
     &        (DFGOX**3.2D-1)*((1.D+2*D50X)**3.6D-1)
!
!---      Kinetic volatilization coeff. w/ constant model  ---
!
          ELSEIF( ISLC(46).EQ.4) THEN
!
!---        Back calculate D50 using volatilization coeff  ---
!
            VOLKX = D50(IZN)
            PEX = 2.D-2
            D50X = (3.783D-1*(PEX**6.2D-1)*DFGOX/VOLKX)**5.556D+0
            D50X = D50X/1.0D+2
          ENDIF
!
!---      Kinetic volatilization mass transfer of oil between
!         phases, positive being defined from the NAPL-oil to
!         dissolved-oil, kg/s  ---
!
          CALL VOCGSD( T(2,N),PSO(2,N),RHOGOX,INDX )
          TMBP_O(M,N) = VOLKX*PORD(M,N)*SG(M,N)*
     &      (RHOGOX-XGO(M,N)*RHOG(M,N))*VOL(N)
!
!---      Aqueous pore-water Reynolds number  ---
!
!          RELX = MAX( AVPVL(N)*RHOL(M,N)*D50(IZN)/VISL(M,N),1.D-4 )
          RELX = MAX( AVPVL(N),1.0D-10 )*RHOL(M,N)*D50X/VISL(M,N)
!
!---      Sherwood number  ---
!
          SHWX = 3.72D+1*(RELX**6.1D-1)*(SN(M,N)**1.24D+0)
!
!---      Kinetic dissolution coefficient, 1/s  ---
!
          DISKX = SHWX*DFLO(M,N)/(D50X**2)
!
!---      Kinetic dissolution mass transfer of oil between
!         phases, positive being defined from the NAPL-oil to
!         dissolved-oil, kg/s  ---
!
          XMLOX = PSO(2,N)/HCOW
          XMLWX = MAX( 1.D+0-XMLA(M,N)-XMLOX,0.D+0 )
          WTMX = XMLWX*WTMW + XMLA(M,N)*WTMA + XMLOX*WTMO
          XLOX = XMLOX*WTMO/WTMX
          TMBP_O(M,N) = TMBP_O(M,N) + DISKX*PORD(M,N)*SL(M,N)*
     &      (XLOX*RHOL(M,N)-XLO(M,N)*RHOL(M,N))*VOL(N)
  200   CONTINUE
  300 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KV5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDO5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Load the current time step values into the old time step
!     variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE HYST
      USE GRID
      USE FDVP
      USE FDVN
      USE FDVG
      USE FDVD
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
      SUB_LOG(ISUB_LOG) = '/LDO5'
!
!---  Assign old time step values  ---
!
      DO 100 N = 1,NFLD
        PG(1,N) = PG(2,N)
        T(1,N) = T(2,N)
        PORD(1,N) = PORD(2,N)
        PORT(1,N) = PORT(2,N)
        SL(1,N) = SL(2,N)
        SN(1,N) = SN(2,N)
        SG(1,N) = SG(2,N)
        PL(1,N) = PL(2,N)
        PN(1,N) = PN(2,N)
        PVA(1,N) = PVA(2,N)
        PVW(1,N) = PVW(2,N)
        PVO(1,N) = PVO(2,N)
        XGA(1,N) = XGA(2,N)
        XGW(1,N) = XGW(2,N)
        XGO(1,N) = XGO(2,N)
        XMGA(1,N) = XMGA(2,N)
        XMGW(1,N) = XMGW(2,N)
        XMGO(1,N) = XMGO(2,N)
        XLA(1,N) = XLA(2,N)
        XLW(1,N) = XLW(2,N)
        XLO(1,N) = XLO(2,N)
        XMLA(1,N) = XMLA(2,N)
        XMLW(1,N) = XMLW(2,N)
        XMLO(1,N) = XMLO(2,N)
        XSO(1,N) = XSO(2,N)
        RHOL(1,N) = RHOL(2,N)
        RHOML(1,N) = RHOML(2,N)
        RHOMG(1,N) = RHOMG(2,N)
        RHON(1,N) = RHON(2,N)
        RHOG(1,N) = RHOG(2,N)
        VISL(1,N) = VISL(2,N)
        VISG(1,N) = VISG(2,N)
        VISN(1,N) = VISN(2,N)
        RKL(1,1,N) = RKL(1,2,N)
        RKL(2,1,N) = RKL(2,2,N)
        RKL(3,1,N) = RKL(3,2,N)
        RKG(1,N) = RKG(2,N)
        RKN(1,N) = RKN(2,N)
        TORL(1,N) = TORL(2,N)
        TORG(1,N) = TORG(2,N)
        TORN(1,N) = TORN(2,N)
        NPHAZ(1,N) = NPHAZ(2,N)
        TRPNL(1,N) = TRPNL(2,N)
        IF( SN(2,N).GT.EPSL ) ASLMIN(1,N) = MIN( ASL(N),ASLMIN(2,N) )
        ASTMIN(1,N) = MIN( AST(N),ASTMIN(2,N) )
        ASLMIN(2,N) = ASLMIN(1,N)
        ASTMIN(2,N) = ASTMIN(1,N)
        DO 90 NSL = 1,NSOLU
          CO(N,NSL) = C(N,NSL)
   90   CONTINUE
  100 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDO5 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROP5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Control saturation, relative permeability, porosity, and
!     tortuosity calculations
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FDVP
      USE FDVN
      USE FDVG
      USE FDVD
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
      SUB_LOG(ISUB_LOG) = '/PROP5'
!
!---  Loop over all nodes  ---
!
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        IZN = IZ(N)
        DO 100 M = ISVC+2,2,-1
          INDX = 1
!
!---      Interphase scaling factors   ---
!
          XMLOX = MIN( XMLO(M,N),PSO(2,N)/HCOW )
          XMLWX = MAX( 1.D+0-XMLOX,0.D+0 )
          XLOX = XMLOX*WTMO/(XMLOX*WTMO + XMLWX*WTMW)
          PX = MAX( PG(M,N),PL(M,N),PN(M,N) ) + PATM
          CALL WATLQD( T(2,N),PX,RHOLX )
          SIG_GL = SCHR(5,IZN)
          SIG_NL = SCHR(6,IZN)
          SIG_GN = SCHR(7,IZN)
          RHOLOX = XLOX*RHOLX
          CALL SFT5( RHOLOX )
          BGL = SCHR(8,IZN)/SIG_GL
          BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
          BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
!
!---      Maximum effective trapped-NAPL saturation  ---
!
          TRPNLX = SCHR(9,IZN)/(SCHR(9,IZN)+TRPNL(2,N))
          ASNTMX = SCHR(11,IZN)*TRPNLX
!
!---      Maximum effective trapped-gas saturation  ---
!
          ASGTMX = SCHR(17,IZN)
!
!---      Maximum effective residual-NAPL saturation  ---
!
          ASNRMX = SCHR(10,IZN)
!
!---  For phase conditions with dissolved-oil/air compute saturations
!     assuming no NAPL/gas entrapment  ---
!
          INDX = 1
          IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.2 .OR.
     &      NPHAZ(2,N).EQ.5 ) INDX = 3
          IF( NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.4 .OR.
     &      NPHAZ(2,N).EQ.6 ) INDX = 4
          ASLMINX = ASLMIN(1,N)
          ASTMINX = ASTMIN(1,N)
          ASTMAXX = ASTMAX(1,N)
          CALL KSP5( PG(M,N),PN(M,N),PL(M,N),SL(M,N),SN(M,N),SG(M,N),
     &      RKL(1,M,N),RKN(M,N),RKG(M,N),ASTX,ASLX,
     &      ASLMINX,ASTMINX,ASTMAXX,ASNTX,ASNRX,ASGTX,
     &      ASGTLX,ASGTNX,ASNTMX,ASNRMX,ASGTMX,INDX,IZN,N )
          IF( M.EQ.2 ) THEN
            AST(N) = ASTX
            ASL(N) = ASLX
            ASNT(N) = ASNTX
            ASGT(N) = ASGTX
            ASNR(N) = ASNRX
            ASLMIN(2,N) = ASLMINX
            ASTMIN(2,N) = ASTMINX
            ASTMAX(2,N) = ASTMAXX
          ENDIF
          SNT(M,N) = ASNTX*(1.D+0-SCHR(4,IZN))
          SGT(M,N) = ASGTX*(1.D+0-SCHR(4,IZN))
!
!---  Compute porosity and tortuosity  ---
!
          PX = MAX( PG(M,N),PL(M,N),PN(M,N) )+PATM
          CALL PORSTY( N,PX,PCMP(N),PORD(M,N),PORT(M,N) )
          IF( ISLC(3).EQ.1 ) CALL TORTU( IZN,SL(M,N),SG(M,N),
     &      SN(M,N),PORD(M,N),TORL(M,N),TORG(M,N),TORN(M,N) )
  100   CONTINUE
  200 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROP5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE QSAT5( PLX,PGX,PNX,DPGLX,N )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!
!     Calculate pressure differential required to yield a differential
!     change in saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 2001.
!     Last Modified by MD White, Battelle, February 9, 2001.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/QSAT5'
!
!---  Interphase scaling factors   ---
!
      IZN = IZ(N)
      BGL = SCHR(8,IZN)/SCHR(5,IZN)
      BNL = SCHR(8,IZN)/(SCHR(6,IZN)*CA_NL)
      BGN = SCHR(8,IZN)/(SCHR(7,IZN)*CA_GN)
      ENPR = SCHR(2,IZN)*RHORL*GRAV
      HNL = MAX( 0.D+0,BNL*(PNX-PLX)/RHORL/GRAV )
      HGN = MAX( 0.D+0,BGN*(PGX-PNX)/RHORL/GRAV )
!
!---  van Genuchten saturation functions  ---
!
      IF( ISCHR(IZN).EQ.13 .OR. ISCHR(IZN).EQ.15 .OR.
     &  ISCHR(IZN).EQ.17 ) THEN
        CN = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( IRPL(IZN).EQ.2 ) THEN
            CM = 1.D+0 - 2.D+0/CN
          ELSE
            CM = 1.D+0 - 1.D+0/CN
          ENDIF
        ELSE
          CM = SCHR(14,IZN)
        ENDIF
        ASTX = (1.D+0/(1.D+0 + (SCHR(1,IZN)*HGN)**CN))**CM
        ASLX = MIN(ASTX,(1.D+0/(1.D+0 + (SCHR(1,IZN)*HNL)**CN))**CM)
        ASNX = ASTX-ASLX
        IF( ASNX.LE.EPSL ) ASNX = 0.D+0
        ASGX = 1.D+0-ASTX
        IF( ASGX.LE.EPSL ) ASGX = 0.D+0
!
!---  No free NAPL  ---
!
        IF( ASNX.LE.EPSL ) THEN
          IF( ASLX.LT.5.D-1 ) THEN
            ASLY = ASLX + 1.D-6
          ELSE
            ASLY = ASLX - 1.D-6
          ENDIF
          HGL = ((-1.D+0 + (1.D+0/ASLY)**(1.D+0/CM))**(1.D+0/CN))/
     &      SCHR(1,IZN)
          DPGLX = HGL*RHORL*GRAV/BGL + PLX - PGX
        ENDIF
!
!---  Brooks and Corey saturation functions  ---
!
      ELSEIF( ISCHR(IZN).EQ.14 .OR. ISCHR(IZN).EQ.16 .OR.
     &  ISCHR(IZN).EQ.18 ) THEN
        CL = MAX( SCHR(3,IZN),SMALL )
        IF( HGN-SCHR(2,IZN).LE.EPSL ) THEN
          ASTX = 1.D+0
        ELSE
          ASTX = (SCHR(2,IZN)/HGN)**CL
        ENDIF
        IF( HNL-SCHR(2,IZN).LE.EPSL ) THEN
          ASLX = 1.D+0
        ELSE
          ASLX = (SCHR(2,IZN)/HNL)**CL
        ENDIF
        ASNX = ASTX-ASLX
        IF( ASNX.LE.EPSL ) ASNX = 0.D+0
        ASGX = 1.D+0-ASTX
        IF( ASGX.LE.EPSL ) ASGX = 0.D+0
!
!---  No free NAPL  ---
!
        IF( ASNX.LE.EPSL ) THEN
          DPGLY = 1.D-2
          IF( ASLX.LT.5.D-1 ) DPGLY = -1.D-2
          DO 200 I = 1,10
            DPGLY = DPGLY*1.D+1
            HGL = MAX( 0.D+0,BGL*(PGX-PLX+DPGLY)/RHORL/GRAV )
            IF( HGL-SCHR(2,IZN).LE.EPSL ) THEN
              ASLY = 1.D+0
            ELSE
              ASLY = (SCHR(2,IZN)/HGL)**CL
            ENDIF
            IF( ABS(ASLY-ASLX).GT.1.D-7 ) GOTO 210
  200     CONTINUE
  210     CONTINUE
          DO 220 I = 1,10
            DPGLZ = DPGLY/(1.D+1**(0.5**(10-I)))
            HGL = MAX( 0.D+0,BGL*(PGX-PLX+DPGLZ)/RHORL/GRAV )
            IF( HGL-SCHR(2,IZN).LE.EPSL ) THEN
              ASLY = 1.D+0
            ELSE
              ASLY = (SCHR(2,IZN)/HGL)**CL
            ENDIF
            IF( ABS(ASLY-ASLX).GT.1.D-7 .AND.
     &        ABS(ASLY-ASLX).LE.1.D-5 ) THEN
              DPGLX = DPGLZ
              GOTO 900
            ENDIF
  220     CONTINUE
          INDX = 20
          IMSG = N
          CHMSG = 'Saturation Increment Failure, Node: '
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
  900 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of QSAT5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDBC5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Read input file for boundary condition information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     Last Modified by MD White, PNNL, 24 June 2002.
!     Last Modified by ML Rockhold, PNNL, 29 January 2003.
!     Last Modified by M Oostrom, PNNL, 20 March 2003.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      CHARACTER*64 ADUM,BDUM(LUK+LSOLU),FDUM
      CHARACTER*64 UNTS
      CHARACTER*32 CHTYP(19)
      CHARACTER*512 CHDUM
      REAL*8 VAR(LBTM,LBCV)
      INTEGER ITYP(LUK+LSOLU),NCHB(LUK+LSOLU)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CHTYP
      DATA CHTYP /'Dirichlet ','Neumann ','Zero Flux ','Saturated ',
     &            'Unit Gradient ','Free Gradient ','Outflow ',
     &            'Aqueous Concentration ','Gas Concentration ',
     &            'NAPL Concentration ','Hydraulic Gradient ',
     &            'Initial Condition ','Inflow ',
     &            'Inflow Aqueous-Phase ','Inflow Gas-Phase ',
     &            'Inflow NAPL ','Seepage Face ',
     &            'Fluctuating Water Table ','Dirichlet-Outflow'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDBC5'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Boundary Conditions Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      NBC = 0
   91 READ(IRD,'(A)') CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 91
      ISTART = 1
      VARB = 'Number of Boundary Condition Cards'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO 400 NB = 1, NLIN
   92   READ(IRD,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 92
        CALL LCASE( CHDUM )
!
!---  Read boundary orientation  ---
!
        ISTART = 1
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
        VARB = 'Boundary Type'
        CALL RDCHR(ISTART,ICOMMA,NCHB(IEQW),CHDUM,BDUM(IEQW))
        CALL RDCHR(ISTART,ICOMMA,NCHB(IEQA),CHDUM,BDUM(IEQA))
        CALL RDCHR(ISTART,ICOMMA,NCHB(IEQO),CHDUM,BDUM(IEQO))
        IF( IEQC.GT.0 ) THEN
          DO 5 NSL = 1,NSOLU
            BDUM(NSL+LUK) = 'zero flux'
            IDFLT = 1
            VARB = 'Solute Boundary Type'
            CALL RDCHR(ISTART,ICOMMA,NCHB(NSL+LUK),CHDUM,BDUM(NSL+LUK))
    5     CONTINUE
        ENDIF
!
!---  Read aqueous boundary condition type ---
!
        IF( INDEX(BDUM(IEQW)(1:),'dirichlet').NE.0 ) THEN
           ITYP(IEQW) = 1
        ELSEIF( INDEX(BDUM(IEQW)(1:),'neumann').NE.0 ) THEN
           ITYP(IEQW) = 2
        ELSEIF( INDEX(BDUM(IEQW)(1:),'zero flux').NE.0 ) THEN
           ITYP(IEQW) = 3
        ELSEIF( INDEX(BDUM(IEQW)(1:),'unit gradient').NE.0 ) THEN
           ITYP(IEQW) = 5
        ELSEIF( INDEX(BDUM(IEQW)(1:),'outflow').NE.0 ) THEN
           ITYP(IEQW) = 7
        ELSEIF( INDEX(BDUM(IEQW)(1:),'hydraulic gradient').NE.0 ) THEN
           ITYP(IEQW) = 11
        ELSEIF( INDEX(BDUM(IEQW)(1:),'initial cond').NE.0 ) THEN
           ITYP(IEQW) = 12
        ELSEIF( INDEX(BDUM(IEQW)(1:),'fluctuating water').NE.0 ) THEN
           ITYP(IEQW) = 18
        ELSE
          INDX = 4
          NCH = NCHB(IEQW)
          CHMSG = 'Unrecognized Boundary Condition: '//BDUM(IEQW)(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Read gas boundary condition type ---
!
        IF( INDEX(BDUM(IEQA)(1:),'dirichlet-outflow').NE.0 ) THEN
           ITYP(IEQA) = 19
        ELSEIF( INDEX(BDUM(IEQA)(1:),'dirichlet').NE.0 ) THEN
           ITYP(IEQA) = 1
        ELSEIF( INDEX(BDUM(IEQA)(1:),'neumann').NE.0 ) THEN
           ITYP(IEQA) = 2
        ELSEIF( INDEX(BDUM(IEQA)(1:),'zero flux').NE.0 ) THEN
           ITYP(IEQA) = 3
        ELSEIF( INDEX(BDUM(IEQA)(1:),'unit gradient').NE.0 ) THEN
           ITYP(IEQA) = 5
        ELSEIF( INDEX(BDUM(IEQA)(1:),'outflow').NE.0 ) THEN
           ITYP(IEQA) = 7
        ELSEIF( INDEX(BDUM(IEQA)(1:),'hydraulic gradient').NE.0 ) THEN
           ITYP(IEQA) = 11
        ELSEIF( INDEX(BDUM(IEQA)(1:),'initial cond').NE.0 ) THEN
           ITYP(IEQA) = 12
        ELSE
          INDX = 4
          NCH = NCHB(IEQA)
          CHMSG = 'Unrecognized Boundary Condition: '//BDUM(IEQA)(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Read NAPL boundary condition type ---
!
        IF( INDEX(BDUM(IEQO)(1:),'dirichlet').NE.0 ) THEN
           ITYP(IEQO) = 1
        ELSEIF( INDEX(BDUM(IEQO)(1:),'neumann').NE.0 ) THEN
           ITYP(IEQO) = 2
        ELSEIF( INDEX(BDUM(IEQO)(1:),'zero flux').NE.0 ) THEN
           ITYP(IEQO) = 3
        ELSEIF( INDEX(BDUM(IEQO)(1:),'unit gradient').NE.0 ) THEN
           ITYP(IEQO) = 5
        ELSEIF( INDEX(BDUM(IEQO)(1:),'outflow').NE.0 ) THEN
           ITYP(IEQO) = 7
        ELSEIF( INDEX(BDUM(IEQO)(1:),'hydraulic gradient').NE.0 ) THEN
           ITYP(IEQO) = 11
        ELSEIF( INDEX(BDUM(IEQO)(1:),'initial cond').NE.0 ) THEN
           ITYP(IEQO) = 12
        ELSE
          INDX = 4
          NCH = NCHB(IEQO)
          CHMSG = 'Unrecognized Boundary Condition: '//BDUM(IEQO)(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Read solute boundary condition type(s) ---
!
        IF( IEQC.GT.0 ) THEN
          DO 25 NSL = 1,NSOLU
            IF( INDEX(BDUM(NSL+LUK)(1:),'outflow').NE.0 ) THEN
               ITYP(NSL+LUK) = 7
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'initial cond').NE.0 ) THEN
               ITYP(NSL+LUK) = 12
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'inflow vol').NE.0 ) THEN
               ITYP(NSL+LUK) = 13
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'inflow aqu').NE.0 ) THEN
               ITYP(NSL+LUK) = 14
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'inflow gas').NE.0 ) THEN
               ITYP(NSL+LUK) = 15
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'inflow napl').NE.0 ) THEN
               ITYP(NSL+LUK) = 16
            ELSEIF(INDEX(BDUM(NSL+LUK)(1:),
     &        'volumetric conc').NE.0 ) THEN
               ITYP(NSL+LUK) = 1
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'aqueous conc').NE.0 )THEN
               ITYP(NSL+LUK) = 8
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'gas conc').NE.0 ) THEN
               ITYP(NSL+LUK) = 9
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'napl conc').NE.0 ) THEN
               ITYP(NSL+LUK) = 10
            ELSEIF( INDEX(BDUM(NSL+LUK)(1:),'zero flux').NE.0 ) THEN
               ITYP(NSL+LUK) = 3
            ELSE
              INDX = 4
              NCH = NCHB(NSL+LUK)
              CHMSG = 'Unrecognized Solute Boundary Condition: '//
     &          BDUM(NSL+LUK)(1:NCH)
              CALL WRMSGS( INDX )
            ENDIF
   25     CONTINUE
        ENDIF
!
!---  Write boundary condition type(s) ---
!
        WRITE(IWR,'(A)') 'Boundary Condition Type: '
        WRITE(IWR,'(2X,2A)') 'Aqueous: ',CHTYP(ITYP(IEQW))
        WRITE(IWR,'(2X,2A)') 'Gas :    ',CHTYP(ITYP(IEQA))
        WRITE(IWR,'(2X,2A)') 'NAPL:    ',CHTYP(ITYP(IEQO))
        IF( IEQC.GT.0 ) THEN
          DO 30 NSL = 1,NSOLU
            IDB = INDEX( SOLUT(NSL)(1:),'  ') - 1
            WRITE(IWR,'(2X,2A)') CHTYP(ITYP(NSL+LUK)),SOLUT(NSL)(1:IDB)
   30     CONTINUE
        ENDIF
!
!---  Read and write boundary domain indices  ---
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
          INDX = 7
          CHMSG = 'Zero Boundary Condition Times: Boundary Number: '
          IMSG = NB
          CALL WRMSGS( INDX )
        ELSE
          INDX = 7
          CHMSG = 'Number of Cyclic Boundary Condition Times < 3 : '
     &      //'Boundary Number: '
          IMSG = NB
          CALL WRMSGS( INDX )
        ENDIF
        IF( IBCMX.GT.LBTM ) THEN
          INDX = 6
          CHMSG = 'Number of Boundary Condition Times > Parameter LBTM'
          CALL WRMSGS( INDX )
        ENDIF
        BCTMO = -SMALL
        WRITE(IWR,'(A)') 'Boundary Condition Times and Variables:'
        DO 100 NTM = 1,IBCMX
          DO 40 M = 1,LBCV
            VAR(NTM,M) = 0.D+0
   40     CONTINUE
!
!---  Read, write, and convert boundary condition time, variables,
!     and units  ---
!
   94     READ(IRD,'(A)') CHDUM
          IF( CHDUM(1:1).EQ.'#' ) GOTO 94
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
!---      Aqueous Dirichlet boundary condition  ---
!
          IF( ITYP(IEQW).EQ.1 ) THEN
            VARB = 'Aqueous Pressure'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
              VAR(NTM,2) = VAR(NTM,2) - PATM
!
!---      Aqueous Neumann boundary condition  ---
!
          ELSEIF( ITYP(IEQW).EQ.2 ) THEN
            VARB = 'Volumetric Aqueous Flux'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNM = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
!
!---      Aqueous zero-flux boundary condition  ---
!
          ELSEIF( ITYP(IEQW).EQ.3 ) THEN
            VARB = 'Aqueous Pressure'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
              VAR(NTM,2) = VAR(NTM,2) - PATM
!
!---      Aqueous hydraulic-gradient boundary condition  ---
!
          ELSEIF( ITYP(IEQW).EQ.11 ) THEN
            VARB = 'Base Aqueous Pressure'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
              VAR(NTM,2) = VAR(NTM,2) - PATM
!
!---      Aqueous fluctuating water table boundary condition  ---
!         (Hydraulic gradient for liquid saturated nodes with
!          zero-flux for overlying unsaturated nodes)
!
          ELSEIF( ITYP(IEQW).EQ.18 ) THEN
            VARB = 'Base Aqueous Pressure'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,2)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
              VAR(NTM,2) = VAR(NTM,2) - PATM
!
!---      Read dummy variables  ---
!
          ELSE
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          ENDIF
!
!---      Aqueous dissolved-air relative saturation  ---
!
          VARB = 'Aqueous Dissolved-Air Relative Saturation'
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),': '
          ISX = ISTART
          ICX = ICOMMA
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(1PE11.4)') VAR(NTM,5)
!
!---      Aqueous dissolved-oil relative saturation  ---
!
          VARB = 'Aqueous Dissolved-Oil Relative Saturation'
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),': '
          ISX = ISTART
          ICX = ICOMMA
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(1PE11.4)') VAR(NTM,6)
!
!---      Gas Dirichlet or Dirichlet-Outflow boundary condition  ---
!
          IF( ITYP(IEQA).EQ.1 .OR. ITYP(IEQA).EQ.19 ) THEN
            VARB = 'Gas Pressure'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
              VAR(NTM,3) = VAR(NTM,3) - PATM
!
!---      Gas Neumann boundary condition  ---
!
          ELSEIF( ITYP(IEQA).EQ.2 ) THEN
            VARB = 'Volumetric Gas Flux'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
              INDX = 0
              IUNM = 1
              IUNS = -1
              CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
!
!---      Gas zero-flux boundary condition  ---
!
          ELSEIF( ITYP(IEQA).EQ.3 ) THEN
            VARB = 'Gas Pressure'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
              VAR(NTM,3) = VAR(NTM,3) - PATM
!
!---      Gas hydraulic-gradient boundary condition  ---
!
          ELSEIF( ITYP(IEQA).EQ.11 ) THEN
            VARB = 'Base Gas Pressure'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,3)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
              VAR(NTM,3) = VAR(NTM,3) - PATM
          ELSE
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          ENDIF
!
!---      Gas water-vapor relative humidity  ---
!
          VARB = 'Gas Water-Vapor Relative Humidity'
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),': '
          ISX = ISTART
          ICX = ICOMMA
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,7))
            WRITE(IWR,'(1PE11.4)') VAR(NTM,7)
!
!---      Gas oil-vapor relative humidity  ---
!
          VARB = 'Gas Oil-Vapor Relative Humidity'
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),': '
          ISX = ISTART
          ICX = ICOMMA
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
            ISTART = ISX
            ICOMMA = ICX
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,8))
            WRITE(IWR,'(1PE11.4)') VAR(NTM,8)
!
!---      NAPL Dirichlet boundary condition  ---
!
          IF( ITYP(IEQO).EQ.1 ) THEN
            VARB = 'NAPL Pressure'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,4)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
              VAR(NTM,4) = VAR(NTM,4) - PATM
          ELSEIF( ITYP(IEQO).EQ.2 ) THEN
            VARB = 'Volumetric NAPL Flux, '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,A,$)') VARB(1:IVR)
            WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,4)
            INDX = 0
            IUNM = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
!
!---      NAPL zero-flux boundary condition  ---
!
          ELSEIF( ITYP(IEQO).EQ.3 ) THEN
            VARB = 'NAPL Pressure'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,4)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
              VAR(NTM,4) = VAR(NTM,4) - PATM
!
!---      NAPL hydraulic-gradient boundary condition  ---
!
          ELSEIF( ITYP(IEQO).EQ.11 ) THEN
            VARB = 'Base NAPL Pressure'
            WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),', '
            ISX = ISTART
            ICX = ICOMMA
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
              ISTART = ISX
              ICOMMA = ICX
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
              CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
              WRITE(IWR,'(2A,1PE11.4)') UNTS(1:NCH),': ',VAR(NTM,4)
              INDX = 0
              IUNM = -1
              IUNKG = 1
              IUNS = -2
              CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
              VAR(NTM,4) = VAR(NTM,4) - PATM
          ELSE
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          ENDIF
!
!---      Solute boundary conditions  ---
!
          IF( IEQC.GT.0 ) THEN
            DO 50 NSL = 1,NSOLU
              IDB = INDEX( SOLUT(NSL)(1:),'  ') - 1
!
!---          Volumetric concentration or inflow boundary condition  ---
!
              IF( ITYP(NSL+LUK).EQ.1 .OR. ITYP(NSL+LUK).EQ.13 ) THEN
                VARB = 'Volumetric Concentration'
                ISX = ISTART
                ICX = ICOMMA
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
                  ISTART = ISX
                  ICOMMA = ICX
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+LBCU))
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL)(1:IDB),
     &              VARB(1:IVR),', ',UNTS(1:NCH),': ',VAR(NTM,NSL+LBCU)
                  INDX = 0
                  IUNM = -3
                  CALL RDUNIT(UNTS,VAR(NTM,NSL+LBCU),INDX)
!
!---          Aqueous concentration or inflow boundary condition  ---
!
              ELSEIF( ITYP(NSL+LUK).EQ.8 .OR. ITYP(NSL+LUK).EQ.14 ) THEN
                VARB = 'Aqueous-Phase Volumetric Concentration'
                ISX = ISTART
                ICX = ICOMMA
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
                  ISTART = ISX
                  ICOMMA = ICX
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+LBCU))
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL)(1:IDB),
     &              VARB(1:IVR),', ',UNTS(1:NCH),': ',VAR(NTM,NSL+LBCU)
                  INDX = 0
                  IUNM = -3
                  CALL RDUNIT(UNTS,VAR(NTM,NSL+LBCU),INDX)
!
!---          Gas concentration or inflow boundary condition  ---
!
              ELSEIF( ITYP(NSL+LUK).EQ.9 .OR. ITYP(NSL+LUK).EQ.15 ) THEN
                VARB = 'Gas-Phase Volumetric Concentration'
                ISX = ISTART
                ICX = ICOMMA
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
                  ISTART = ISX
                  ICOMMA = ICX
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+LBCU))
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL)(1:IDB),
     &              VARB(1:IVR),', ',UNTS(1:NCH),': ',VAR(NTM,NSL+LBCU)
                  INDX = 0
                  IUNM = -3
                  CALL RDUNIT(UNTS,VAR(NTM,NSL+LBCU),INDX)
!
!---          NAPL concentration or inflow boundary condition  ---
!
              ELSEIF( ITYP(NSL+LUK).EQ.10 .OR. ITYP(NSL+LUK).EQ.16 )THEN
                VARB = 'NAPL Volumetric Concentration'
                ISX = ISTART
                ICX = ICOMMA
                CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,FDUM)
                  ISTART = ISX
                  ICOMMA = ICX
                  CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,NSL+LBCU))
                  CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
                  WRITE(IWR,'(2X,5A,1PE11.4)') SOLUT(NSL)(1:IDB),
     &              VARB(1:IVR),', ',UNTS(1:NCH),': ',VAR(NTM,NSL+LBCU)
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
            INDX = 7
            CHMSG = 'Boundary Condition Time Sequence: '//
     &        'Boundary Number: '
            IMSG = NB
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
!---    Assign values to boundary variables  ---
!
        NBCL = 0
        DO 320 K = K1X, K2X
          DO 310 J = J1X, J2X
            DO 300 I = I1X, I2X
              N = ND(I,J,K)
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
              IF( IBCDX.EQ.-3 .AND. K.NE.1) THEN
                IF( IXP(N-IJFLD).NE.0 .AND. INBS(1,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'Bottom Boundary'
                  WRITE(IWR,'(A)') 'Bottom Boundary'
                ENDIF
              ELSEIF( IBCDX.EQ.-2 .AND. J.NE.1) THEN
                IF( IXP(N-IFLD).NE.0 .AND. INBS(2,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'South Boundary'
                  WRITE(IWR,'(A)') 'South Boundary'
                ENDIF
              ELSEIF( IBCDX.EQ.-1 .AND. I.NE.1) THEN
                IF( IXP(N-1).NE.0 .AND. INBS(3,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'West Boundary'
                  WRITE(IWR,'(A)') 'West Boundary'
                ENDIF
              ELSEIF( IBCDX.EQ.1 .AND. I.NE.IFLD) THEN
                IF( IXP(N+1).NE.0 .AND. INBS(4,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'East Boundary'
                  WRITE(IWR,'(A)') 'East Boundary'
                ENDIF
              ELSEIF( IBCDX.EQ.2 .AND. J.NE.JFLD) THEN
                IF( IXP(N+IFLD).NE.0 .AND. INBS(5,N).EQ.0 ) THEN
                  IERR = 1
                  WRITE(ISC,'(A)') 'North Boundary'
                  WRITE(IWR,'(A)') 'North Boundary'
                ENDIF
              ELSEIF( IBCDX.EQ.3 .AND. K.NE.KFLD) THEN
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
              IBCT(1,NBC) = ITYP(IEQW)
              IBCT(2,NBC) = ITYP(IEQA)
              IBCT(3,NBC) = ITYP(IEQO)
              IF( IEQC.GT.0 ) THEN
                DO 110 NSL = 1,NSOLU
                  IBCT(NSL+LUK,NBC) = ITYP(NSL+LUK)
  110           CONTINUE
              ENDIF
              IBCM(NBC) = IBCMX
              IBCIN(NBC) = NB
!
!---          Assign hydraulic gradient type boundary condition
!             for the base surface  ---
!
              IF( ITYP(IEQW).EQ.11 ) THEN
                IF( NBCL.EQ.1 ) THEN
                  IBCT(IEQW,NBC) = -11
                  NBHG(1,NBC) = NBC
                ELSE
                  NBHG(1,NBC) = NBHG(1,NBC-1)
                ENDIF
              ENDIF
              IF( ITYP(IEQA).EQ.11 ) THEN
                IF( NBCL.EQ.1 ) THEN
                  IBCT(IEQA,NBC) = -11
                  NBHG(2,NBC) = NBC
                ELSE
                  NBHG(2,NBC) = NBHG(2,NBC-1)
                ENDIF
              ENDIF
              IF( ITYP(IEQO).EQ.11 ) THEN
                IF( NBCL.EQ.1 ) THEN
                  IBCT(IEQO,NBC) = -11
                  NBHG(3,NBC) = NBC
                ELSE
                  NBHG(3,NBC) = NBHG(3,NBC-1)
                ENDIF
              ENDIF
              IF( ITYP(IEQW).EQ.18 ) THEN
                IF( NBCL.EQ.1 ) THEN
                  IBCT(IEQW,NBC) = -18
                  NBHG(1,NBC) = NBC
                ELSE
                  NBHG(1,NBC) = NBHG(1,NBC-1)
                ENDIF
              ENDIF
!
!---  Check for double boundary conditions  ---
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
  300       CONTINUE
  310     CONTINUE
  320   CONTINUE
  400 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDBC5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDOT5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!
!     Reads the oil transport card for solid-aqueous 
!     partition coefficients and hydraulic dispersivities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1996.
!     Last Modified by MD White, PNNL, August 1, 1996.
!     Last Modified by MD White, PNNL, 8 January 2004.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/RDOT5'
!
!---  Write card information to output file  ---
!
      CARD = 'Oil Transport Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      ISLC(12) = 1
      IDSPD = 0
!
!---  Kinetic volatilization option  ---
!
      IF( ISLC(46).GE.1 ) THEN
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---    Read volatilization model  ---
!
        VARB = 'Kinetic Volatilization Model'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX(ADUM(1:),'wilkins').NE.0 ) THEN
          WRITE(IWR,'(/,2A)') VARB(1:IVR),': Wilkins et al.'
          ISLC(46) = 1
        ELSEIF( INDEX(ADUM(1:),'van der ham').NE.0 .OR.
     &    INDEX(ADUM(1:),'brouwers').NE.0 ) THEN
          WRITE(IWR,'(/,2A)') VARB(1:IVR),': Van der Ham and Brouwers'
          ISLC(46) = 2
        ELSEIF( INDEX(ADUM(1:),'yoon').NE.0 ) THEN
          WRITE(IWR,'(/,2A)') VARB(1:IVR),': Yoon et al.'
          ISLC(46) = 3
        ELSEIF( INDEX(ADUM(1:),'constant').NE.0 ) THEN
          WRITE(IWR,'(/,2A)') VARB(1:IVR),': Constant Coefficient'
          ISLC(46) = 4
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Kinetic Volatilization Model: ' 
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Check for advection-dispersion-diffusion
!     transport scheme, default to upwind  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      IADTSX = 0
      VARB = 'Oil Advection-Dispersion-Diffusion Transport Scheme'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
      IF( INDEX(ADUM(1:),'patankar').NE.0 ) THEN
        ISLC(8) = 0
        WRITE(IWR,'(A)') 'Patankar Advection-Dispersion-' // 
     &    'Diffusion Transport Scheme'
        IADTSX = 1
      ELSEIF( INDEX(ADUM(1:),'tvd').NE.0 .OR.
     &  INDEX(ADUM(1:),'leonard').NE.0 ) THEN
        ISLC(8) = 1
        WRITE(IWR,'(A)') 'Leonard-TVD Advection-Dispersion-' //
     &    'Diffusion Transport Scheme'
        IADTSX = 1
      ELSEIF( INDEX(ADUM(1:),'roe').NE.0 .OR.
     &  INDEX(ADUM(1:),'superbee').NE.0 ) THEN
        ISLC(8) = 2
        WRITE(IWR,'(A)') 'Roe''s Superbee Advection-Dispersion-' //
     &    'Diffusion Transport Scheme'
        IADTSX = 1
      ELSEIF( INDEX(ADUM(1:),'upwind').NE.0 ) THEN
        ISLC(8) = 3
        WRITE(IWR,'(A)') 'Upwind Advection-Dispersion-' //
     &    'Diffusion Transport Scheme'
        IADTSX = 1
      ELSE
        ISLC(8) = 3
        WRITE(IWR,'(A)') 'Upwind Advection-Dispersion-' // 
     &    'Diffusion Transport Scheme'
      ENDIF
!
!---  Maximum Courant number  ---
!
      IF( IADTSX.EQ.1 ) THEN
        CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
        IF( INDX.EQ.1 ) THEN
          VARB = 'Maximum Courant Number'
          IDFLT = 1
          CRNTMXC = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,CRNTMXC)
          WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),CRNTMXC
        ENDIF
      ENDIF
!
!---  Loop over the rock/soil saturation information lines  ---
!
      N = 0
      IJK = 0
   10 CONTINUE
        IF( N.GE.NROCK .OR. IJK.GT.0 ) GOTO 600
        IF( IADTSX.EQ.1 ) THEN
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
        ENDIF
        IADTSX = 1
        ISTART = 1
        VARB = 'Rock Name: '
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,RDUM)
!
!---    IJK, KIJ, or JKI indexing  ---
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
!---    Search known rock types for a matching type ---
!
        DO 100 M = 1, NROCK
          IF( RDUM.EQ.ROCK(M)) THEN
            IROCK = M
            GOTO 200
          ENDIF
  100   CONTINUE
!
!---    Search known scaling groups for a matching type ---
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
!---    Loop over rock/soils within scaling group  ---
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
!---    Longitudinal dispersivity  ---
!
        VARB = 'Longitudinal Dispersivity: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DPLD(IROCK))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),
     &    UNTS(1:NCH),': ',DPLD(IROCK)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DPLD(IROCK),INDX)
        IF( DPLD(IROCK).GE.SMALL ) IDSPD = 1
!
!---    Transverse dispersivity  ---
!
        VARB = 'Transverse Dispersivity: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DPTD(IROCK))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),
     &    UNTS(1:NCH),': ',DPTD(IROCK)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DPTD(IROCK),INDX)
        IF( DPTD(IROCK).GE.SMALL ) IDSPD = 1
!
!---    Read adsorption function  ---
!
        VARB = 'Solid-Aqueous Oil Adsorption Function'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'linear').NE.0 .AND.
     &    INDEX(ADUM(1:),'kd').NE.0 ) THEN
          IPCSLD(IROCK) = 1
          WRITE(IWR,'(A)') 'Linear Isotherm w/ Kd'
          WRITE(IWR,'(A)') 'Kd = ' //
     &      '(m^3 aqu/kg soil)(kg oil sorbed/kg oil aqu)'
          VARB = 'Kd parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSLD(1,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',PCSLD(1,IROCK)
          INDX = 0
          IUNKG = -1
          IUNM = 3
          CALL RDUNIT(UNTS,PCSLD(1,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PCSLD(1,IROCK),', m^3/kg)'
        ELSEIF( INDEX(ADUM(1:),'linear').NE.0 ) THEN
          IPCSLD(IROCK) = 2
          WRITE(IWR,'(A)') 'Linear Isotherm'
          WRITE(IWR,'(A)') 'k = ' //
     &      '(m^3 aqu/m^3 node)(kg oil sorbed/kg oil aqu)'
          VARB = 'k parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSLD(1,IROCK))
          WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),': ',PCSLD(1,IROCK)
        ELSEIF( INDEX(ADUM(1:),'freundlich').NE.0 ) THEN
          IPCSLD(IROCK) = 3
          WRITE(IWR,'(A)') 'Freundlich Isotherm'
          WRITE(IWR,'(A)') 'k = ((m^3 aqu/kg oil aqu)^(1/n))' //
     &      '(kg oil sorbed/m^3 node)'
          VARB = 'k parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSLD(1,IROCK))
          WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),': ',PCSLD(1,IROCK)
          VARB = 'n Parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSLD(2,IROCK))
          WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PCSLD(2,IROCK)
        ELSEIF( INDEX(ADUM(1:),'langmuir').NE.0 ) THEN
          IPCSLD(IROCK) = 4
          WRITE(IWR,'(A)') 'Langmuir Isotherm'
          WRITE(IWR,'(A)') 'a = ' //
     &      '(m^3 aqu/m^3 node)(kg oil sorbed/kg oil aqu)'
          WRITE(IWR,'(A)') 'b = (m^3 aqu/kg oil aqu)'
          VARB = 'a parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSLD(1,IROCK))
          WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),': ',PCSLD(1,IROCK)
          VARB = 'b parameter'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSLD(2,IROCK))
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',PCSLD(2,IROCK)
          INDX = 0
          IUNKG = -1
          IUNM = 3
          CALL RDUNIT(UNTS,PCSLD(2,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PCSLD(2,IROCK),', m^3/kg)'
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Oil Adsorption Function: ' 
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Kinetic volatilization coeff.:
!       1 - w/ Wilkins et al. model
!       2 - w/ Van der Ham and Brouwers model
!       3 - w/ Yoon et al. model
!
        IF( ISLC(46).GE.1 .AND. ISLC(46).LE.3 ) THEN
!
!---      (d50) Mean particle diameter  ---
!
          VARB = 'Mean Particle Diameter (d50)'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,D50(IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',D50(IROCK)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,D50(IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',D50(IROCK),', m)'
!
!---    Kinetic volatilization coeff.: constant  ---
!
        ELSEIF( ISLC(46).EQ.4 ) THEN
!
!---      Volatilization Coefficient (stored in d50 variable)  ---
!
          VARB = 'Volatilization Coefficient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,D50(IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',D50(IROCK)
          INDX = 0
          IUNS = -1
          CALL RDUNIT(UNTS,D50(IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',D50(IROCK),', 1/s)'
        ENDIF
!
!---  Read next rock/soil type or scaling group  ---
!
        IF( N.LT.NROCK ) WRITE(IWR,'(/)')
        GOTO 10
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDOT5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDIC5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Read input file for initial conditions information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!     Last Modified by MD White, PNNL, October 15, 1997.
!     Last Modified by MD White, PNNL, 24 June 2002.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE HYST
      USE GRID
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
      CHARACTER*64 ADUM,BDUM,FDUM,FMDUM,UNTS
      CHARACTER*24 CHLB(3)
      CHARACTER*512 CHDUM
      INTEGER IDOM(6)
      REAL*8 VAR(5)
      LOGICAL FCHK
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
      SUB_LOG(ISUB_LOG) = '/RDIC5'
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
        ISIC = 3
      ENDIF
!
!---  Read initial conditions  ---
!
      WRITE(IWR,'(/,A)') 'Initial Condition Variable(s) and Domain(s)'
   92 READ(IRD,'(A)') CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 92
      ISTART = 1
      VARB = 'Number of Initial Condition Cards: '
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO 1000 NL = 1, NLIN
   93   READ(IRD,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 93
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Initial Condition Variable: '
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'overwrite').EQ.0 .AND.
     &    ( IEO.EQ.2 ) ) GOTO 1000
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
        ELSEIF( INDEX( ADUM(1:),'napl pres' ).NE.0 ) THEN
          VARB = 'Initial NAPL Pressure'
          IUNM = -1
          IUNKG = 1
          IUNS = -2
        ELSEIF( INDEX( ADUM(1:),'temperature' ).NE.0 ) THEN
          VARB = 'Initial Temperature'
          IUNK = 1
        ELSEIF( INDEX( ADUM(1:),'trapped napl' ).NE.0 ) THEN
          VARB = 'Initial Trapped NAPL Saturation'
        ELSEIF( INDEX( ADUM(1:),'trapped gas' ).NE.0 ) THEN
          VARB = 'Initial Trapped Gas Saturation'
        ELSEIF( INDEX( ADUM(1:),'oil conc' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Oil Concentratin'
          IUNM = -3
          IUNKG = 1
        ELSEIF( INDEX( ADUM(1:),'oil mole' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Oil Mole Fraction'
        ELSEIF( INDEX( ADUM(1:),'oil mass' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Oil Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'air conc' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Air Concentration'
          IUNM = -3
          IUNKG = 1
        ELSEIF( INDEX( ADUM(1:),'air mole' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Air Mole Fraction'
        ELSEIF( INDEX( ADUM(1:),'air mass' ).NE.0 ) THEN
          VARB = 'Initial Dissolved-Air Mass Fraction'
        ELSEIF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
          VARB = 'Initial Solute Concentration'
          IUNM = -3
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: '//
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
          DO 20 M = 1, NROCK
            IF( FDUM .EQ. ROCK(M)) THEN
            IROCK = M
            GOTO 30
          ENDIF
   20     CONTINUE
          INDX = 2
          CHMSG = 'Unrecognized Rock/Soil Type: '// FDUM(1:NCHF)
          CALL WRMSGS( INDX )
          GOTO 1000
   30     CONTINUE
          WRITE(IWR,'(2X,3A,1PE11.4,2A)') ADUM(1:NCHA),UNTS(1:NCHU),
     &      ': ',VAR(1),' Rock/Soil Type: ',FDUM(1:NCHF)
          INDX = 0
          CALL RDUNIT( UNTS,VAR(1),INDX )
        ELSE
!
!---  Read initial condtions input from the input file  ---
!
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
!---    Read aqueous pressure  ---
!
        IF( INDEX(ADUM(1:),'aqueous pres').NE.0 ) THEN
          ADDER = -PATM
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDINFS( PL,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PL,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( PL,VAR,ADDER,IDOM,INDX )
          ENDIF
!
!---    Read gas pressure  ---
!
        ELSEIF( INDEX(ADUM(1:),'gas pres').NE.0 ) THEN
          ADDER = -PATM
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDINFS( PG,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PG,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( PG,VAR,ADDER,IDOM,INDX )
          ENDIF
!
!---    Read NAPL pressure  ---
!
        ELSEIF( INDEX(ADUM(1:),'napl pres').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'auto prime').NE.0 ) THEN
            VAR(1) = -1.D+22
            ADDER = 0.D+0
          ELSE
            ADDER = -PATM
          ENDIF
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDINFS( PN,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( PN,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( PN,VAR,ADDER,IDOM,INDX )
          ENDIF
!
!---    Read trapped-NAPL saturation  ---
!
        ELSEIF( INDEX(ADUM(1:),'trapped napl').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            CALL RDINFS( SNT,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( SNT,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( SNT,VAR,ADDER,IDOM,INDX )
          ENDIF
!
!---    Read temperature  ---
!
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNK = 1
            CALL RDINFS( T,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( T,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( T,VAR,ADDER,IDOM,INDX )
          ENDIF
!
!---    Read trapped-gas saturation  ---
!
        ELSEIF( INDEX(ADUM(1:),'trapped gas').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            CALL RDINFS( SGT,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( SGT,VAR(1),ADDER,IROCK,INDX )
          ELSE
            CALL RDINIS( SGT,VAR,ADDER,IDOM,INDX )
          ENDIF
!
!---    Read dissolved-oil concentration  ---
!
        ELSEIF( INDEX(ADUM(1:),'oil conc').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
         IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -3
            IUNKG = 1
            VAR(1) = -VAR(1)
            CALL RDINFS( XLO,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 210 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 210
              XMLO(2,N) = XLO(2,N)
  210       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            VAR(1) = -VAR(1)
            CALL RDINZS( XLO,VAR(1),ADDER,IROCK,INDX )
            DO 211 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 211
              IF( IZ(N).EQ.IROCK ) XMLO(2,N) = XLO(2,N)
  211       CONTINUE
          ELSE
            VAR(1) = -VAR(1)
            CALL RDINIS( XLO,VAR,ADDER,IDOM,INDX )
            DO 212 K = IDOM(5),IDOM(6)
            DO 212 J = IDOM(3),IDOM(4)
            DO 212 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 212
              XMLO(2,N) = XLO(2,N)
  212       CONTINUE
          ENDIF
!
!---    Read dissolved-oil mass fraction  ---
!
        ELSEIF( INDEX(ADUM(1:),'oil mass').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            CALL RDINFS( XLO,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 220 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 220
              XMLO(2,N) = -XLO(2,N)
  220       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XLO,VAR(1),ADDER,IROCK,INDX )
            DO 221 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 221
              IF( IZ(N).EQ.IROCK ) XMLO(2,N) = -XLO(2,N)
  221       CONTINUE
          ELSE
            CALL RDINIS( XLO,VAR,ADDER,IDOM,INDX )
            DO 222 K = IDOM(5),IDOM(6)
            DO 222 J = IDOM(3),IDOM(4)
            DO 222 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 222
              XMLO(2,N) = -XLO(2,N)
  222       CONTINUE
          ENDIF
!
!---    Read dissolved-oil mole fraction  ---
!
        ELSEIF( INDEX(ADUM(1:),'oil mole').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            CALL RDINFS( XMLO,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 230 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 230
              XLO(2,N) = -XMLO(2,N)
  230       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XMLO,VAR(1),ADDER,IROCK,INDX )
            DO 231 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 231
              IF( IZ(N).EQ.IROCK ) XLO(2,N) = -XMLO(2,N)
  231       CONTINUE
          ELSE
            CALL RDINIS( XMLO,VAR,ADDER,IDOM,INDX )
            DO 232 K = IDOM(5),IDOM(6)
            DO 232 J = IDOM(3),IDOM(4)
            DO 232 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 232
              XLO(2,N) = -XMLO(2,N)
  232       CONTINUE
          ENDIF
!
!---    Read dissolved-air concentration  ---
!
        ELSEIF( INDEX(ADUM(1:),'air conc').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
         IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            IUNM = -3
            IUNKG = 1
            VAR(1) = -VAR(1)
            CALL RDINFS( XLA,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 240 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 240
              XMLA(2,N) = XLA(2,N)
  240       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            VAR(1) = -VAR(1)
            CALL RDINZS( XLA,VAR(1),ADDER,IROCK,INDX )
            DO 241 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 241
              IF( IZ(N).EQ.IROCK ) XMLA(2,N) = XLA(2,N)
  241       CONTINUE
          ELSE
            VAR(1) = -VAR(1)
            CALL RDINIS( XLA,VAR,ADDER,IDOM,INDX )
            DO 242 K = IDOM(5),IDOM(6)
            DO 242 J = IDOM(3),IDOM(4)
            DO 242 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 242
              XMLA(2,N) = XLA(2,N)
  242       CONTINUE
          ENDIF
!
!---    Read dissolved-air mass fraction  ---
!
        ELSEIF( INDEX(ADUM(1:),'air mass').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            CALL RDINFS( XLA,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 250 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 250
              XMLA(2,N) = -XLA(2,N)
  250       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XLA,VAR(1),ADDER,IROCK,INDX )
            DO 251 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 251
              IF( IZ(N).EQ.IROCK ) XMLA(2,N) = -XLA(2,N)
  251       CONTINUE
          ELSE
            CALL RDINIS( XLA,VAR,ADDER,IDOM,INDX )
            DO 252 K = IDOM(5),IDOM(6)
            DO 252 J = IDOM(3),IDOM(4)
            DO 252 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 252
              XMLA(2,N) = -XLA(2,N)
  252       CONTINUE
          ENDIF
!
!---    Read dissolved-air mole fraction  ---
!
        ELSEIF( INDEX(ADUM(1:),'air mole').NE.0 ) THEN
          ADDER = 0.D+0
          INDX = 2
          IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
            CALL RDINFS( XMLA,VAR,ADDER,UNTS,INDX )
            CLOSE(UNIT=26)
            DO 260 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 260
              XLA(2,N) = -XMLA(2,N)
  260       CONTINUE
          ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &      INDEX(ADUM(1:),'zonation').NE.0 )  THEN
            CALL RDINZS( XMLA,VAR(1),ADDER,IROCK,INDX )
            DO 261 N = 1,NFLD
              IF( IXP(N).EQ.0 ) GOTO 261
              IF( IZ(N).EQ.IROCK ) XLA(2,N) = -XMLA(2,N)
  261       CONTINUE
          ELSE
            CALL RDINIS( XMLA,VAR,ADDER,IDOM,INDX )
            DO 262 K = IDOM(5),IDOM(6)
            DO 262 J = IDOM(3),IDOM(4)
            DO 262 I = IDOM(1),IDOM(2)
              N = ND(I,J,K)
              IF( IXP(N).EQ.0 ) GOTO 262
              XLA(2,N) = -XMLA(2,N)
  262       CONTINUE
          ENDIF
!
!---    Read solute initial conditions  ---
!
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'gas').NE.0 ) THEN
            IVAR = 4
          ELSEIF( INDEX(ADUM(1:),'napl').NE.0 ) THEN
            IVAR = 3
          ELSEIF( INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
            IVAR = 2
          ELSE
            IVAR = 1
          ENDIF
          IF( INDEX( UNTS(1:),'bd' ).NE.0 ) IVAR = -IVAR
          DO 300 NSL = 1,NSOLU
            IDB = INDEX(SOLUT(NSL)(1:),'  ') - 1
            IF( BDUM(1:NCHB).EQ.SOLUT(NSL)(1:IDB) ) THEN
              ADDER = 0.D+0
              IF( INDEX(ADUM(1:),'file').NE.0 ) THEN
                IUNM = -3
                CALL RDINFP( C(1,NSL),VAR,ADDER,ICT(1,NSL),IVAR,UNTS )
                CLOSE(UNIT=26)
              ELSEIF( INDEX(ADUM(1:),'rock').NE.0 .OR.
     &          INDEX(ADUM(1:),'zonation').NE.0 )  THEN
                CALL RDINZP( C(1,NSL),VAR(1),ADDER,ICT(1,NSL),
     &            IVAR,IROCK )
              ELSE
                CALL RDINIP( C(1,NSL),VAR,ADDER,ICT(1,NSL),IVAR,IDOM )
              ENDIF
              GOTO 310
            ENDIF
  300     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute: '//BDUM(1:NCHB)
          CALL WRMSGS( INDX )
  310     CONTINUE
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Initial Condition Variable: '//
     &      ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
 1000 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDIC5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINPT5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Read input file cards.
!     Direct control to card reader subroutines.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     Last Modified by MD White, PNNL, 24 June 2002.
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
      SUB_LOG(ISUB_LOG) = '/RDINTP5'
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
          INDX = 12
          CHMSG = 'Incompatible Operational Mode'
          IMSG = IOM
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
        CALL RDSP5
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
!---  Search input file for napl relative permeability card  ---
!
 1000 CONTINUE
 1009 READ(IRD,'(A)', END=1010) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 1009
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'napl rel').NE.0 ) THEN
        CALL RDNRP
        REWIND(IRD)
        GOTO 1100
      ELSE
        GOTO 1000
      ENDIF
 1010 CONTINUE
      INDX = 4
      CHMSG = 'Missing NAPL Relative Permeability Card'
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
        GOTO 1300
      ELSE
        GOTO 1100
      ENDIF
 1110 CONTINUE
      INDX = 4
      CHMSG = 'Missing Gas Relative Permeability Card'
      CALL WRMSGS( INDX )
!
!---  Search input file for volatile organic compound properties card --
!
 1300 CONTINUE
 1309 READ(IRD,'(A)', END=1310) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 1309
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    ( INDEX(CHDUM(2:),'volatile').NE.0 .OR.
     &      INDEX(CHDUM(2:),'organic').NE.0 .OR.
     &      INDEX(CHDUM(2:),'oil').NE.0 .OR.
     &      INDEX(CHDUM(2:),'voc').NE.0 ) .AND.
     &      INDEX(CHDUM(2:),'prop').NE.0 ) THEN
        CALL RDVOC
        REWIND(IRD)
        GOTO 1400
      ELSE
        GOTO 1300
      ENDIF
 1310 CONTINUE
      INDX = 4
      CHMSG = 'Missing Volatile Organic Compound Properties Card'
      CALL WRMSGS( INDX )
 1400 CONTINUE
!
!---  Search input file for solute/fluid interaction card --
!
 2000 CONTINUE
 2009 READ(IRD,'(A)', END=2010) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2009
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'solute/fluid').NE.0 ) THEN
        CALL RDTF5
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
        CALL RDTP5
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
        CHMSG = 'Missing Solute/Rock Interaction Card'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Search input file for oil transport card --
!
 2200 CONTINUE
 2209 READ(IRD,'(A)', END=2210) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 2209
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    ( INDEX(CHDUM(2:),'oil transport').NE.0 .OR.
     &    INDEX(CHDUM(2:),'voc transport').NE.0 )) THEN
        CALL RDOT5
        REWIND(IRD)
        GOTO 4000
      ELSE
        GOTO 2200
      ENDIF
 2210 CONTINUE
!
!---  Kinetic volatilization option --
!
      IF( ISLC(46).GE.1 ) THEN
        INDX = 4
        CHMSG = 'Missing Oil Transport Card'
        CALL WRMSGS( INDX )
      ELSE
        INDX = 1
        CHMSG = 'Missing Oil Transport Card'
        CALL WRMSGS( INDX )
        REWIND(IRD)
      ENDIF
!
!---  Search input file for initial conditions card --
!
 4000 CONTINUE
 4009 READ(IRD,'(A)', END=4010) CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 4009
      CALL LCASE( CHDUM )
      IF( CHDUM(1:1).EQ.'~' .AND.
     &    INDEX(CHDUM(2:),'initial').NE.0 ) THEN
        CALL RDIC5
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
        ISIC = 3
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
        CALL RDBC5
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
        CALL RDSR5
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
        CALL RDOU5
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
        CALL RDSF5
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
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINPT5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDOU5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Read input file for output information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE OUTPU
      USE GRID
      USE FILES
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
      SAVE FORM
      DATA FORM / '(I6,$)' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDOU5'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Output Control Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read reference node information  ---
!
   90 READ(IRD,'(A)') CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 90
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Reference Nodes'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NREF)
      WRITE(IWR,'(/,A,I6)') 'Reference Node No. and Indices: ',NREF
      DO 100 N = 1,NREF
   91   READ(IRD,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 91
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Reference Node Index'
        CALL RDINT(ISTART,ICOMMA,CHDUM,IRF)
        CALL RDINT(ISTART,ICOMMA,CHDUM,JRF)
        CALL RDINT(ISTART,ICOMMA,CHDUM,KRF)
        IF( IRF.LT.1 .OR. IRF.GT.IFLD .OR. JRF.LT.1 .OR.
     &    JRF.GT.JFLD. OR. KRF.LT.1 .OR. KRF.GT.KFLD) THEN
          INDX = 4
          CHMSG = 'Unrecognized Reference Node Index: '//CHDUM(1:ISTART)
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
   92 READ(IRD,'(A)') CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 92
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
   93 READ(IRD,'(A)') CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 93
      ISTART = 1
      VARB = 'Number of Reference Node Variables'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NVREF)
      WRITE( IWR,'(/,A,I6)') 'Reference Node Variables: ',NVREF
      NVC = 0
      DO 200 NV = 1,NVREF
   94   READ(IRD,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 94
        CALL LCASE( CHDUM )
        ISTART=1
        VARB = 'Reference Node Variable'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Reference Node Variable: Solute Name'
          CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SOLNM)
          DO 110 NSL = 1,NSOLU
            IF( SOLNM.EQ.SOLUT(NSL) ) GOTO 120
  110     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute: '//SOLNM(1:NCS)
          CALL WRMSGS( INDX )
  120     CONTINUE
        ENDIF
        IF( INDEX(ADUM(1:),'aqueous pressure').NE.0 ) THEN
          IREF(NV) = 1
        ELSEIF( INDEX(ADUM(1:),'gas pressure').NE.0 ) THEN
          IREF(NV) = 2
        ELSEIF( INDEX(ADUM(1:),'napl pressure').NE.0 ) THEN
          IREF(NV) = 3
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          IREF(NV) = 4
        ELSEIF( INDEX(ADUM(1:),'phase condition').NE.0 ) THEN
          IREF(NV) = 5
        ELSEIF( INDEX(ADUM(1:),'aqueous gauge pressure').NE.0 ) THEN
          IREF(NV) = 6
        ELSEIF( INDEX(ADUM(1:),'gas gauge pressure').NE.0 ) THEN
          IREF(NV) = 7
        ELSEIF( INDEX(ADUM(1:),'napl gauge pressure').NE.0 ) THEN
          IREF(NV) = 8
        ELSEIF( INDEX(ADUM(1:),'apparent aqueous sat').NE.0 ) THEN
          IREF(NV) = 9
        ELSEIF( INDEX(ADUM(1:),'apparent total sat').NE.0 ) THEN
          IREF(NV) = 10
        ELSEIF( INDEX(ADUM(1:),'effective trapped napl').NE.0 ) THEN
          IREF(NV) = 18
        ELSEIF( INDEX(ADUM(1:),'effective trapped gas').NE.0 ) THEN
          IREF(NV) = 19
        ELSEIF( INDEX(ADUM(1:),'aqueous trapped gas').NE.0 ) THEN
          IREF(NV) = 107
        ELSEIF( INDEX(ADUM(1:),'napl trapped gas').NE.0 ) THEN
          IREF(NV) = 108
        ELSEIF( INDEX(ADUM(1:),'trapped gas sat').NE.0 ) THEN
          IREF(NV) = 105
        ELSEIF( INDEX(ADUM(1:),'trapped napl sat').NE.0 ) THEN
          IREF(NV) = 106
        ELSEIF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
          IREF(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'gas saturation').NE.0 ) THEN
          IREF(NV) = 12
        ELSEIF( INDEX(ADUM(1:),'napl saturation').NE.0 ) THEN
          IREF(NV) = 13
        ELSEIF( INDEX(ADUM(1:),'total saturation').NE.0 ) THEN
          IREF(NV) = 14
        ELSEIF( INDEX(ADUM(1:),'aqueous moisture cont').NE.0 ) THEN
          IREF(NV) = 15
        ELSEIF( INDEX(ADUM(1:),'napl moisture cont').NE.0 ) THEN
          IREF(NV) = 16
        ELSEIF( INDEX(ADUM(1:),'total moisture cont').NE.0 ) THEN
          IREF(NV) = 17
        ELSEIF( INDEX(ADUM(1:),'diffusive porosity').NE.0 ) THEN
          IREF(NV) = 20
        ELSEIF( INDEX(ADUM(1:),'water gas mass frac').NE.0 ) THEN
          IREF(NV) = 21
        ELSEIF( INDEX(ADUM(1:),'air gas mass frac').NE.0 ) THEN
          IREF(NV) = 22
        ELSEIF( INDEX(ADUM(1:),'oil gas mass frac').NE.0 ) THEN
          IREF(NV) = 23
        ELSEIF( INDEX(ADUM(1:),'water aqueous mass frac').NE.0 ) THEN
          IREF(NV) = 24
        ELSEIF( INDEX(ADUM(1:),'air aqueous mass frac').NE.0 ) THEN
          IREF(NV) = 25
        ELSEIF( INDEX(ADUM(1:),'oil aqueous mass frac').NE.0 ) THEN
          IREF(NV) = 26
        ELSEIF( INDEX(ADUM(1:),'aqueous hydraulic head').NE.0 ) THEN
          IREF(NV) = 27
        ELSEIF( INDEX(ADUM(1:),'gas hydraulic head').NE.0 ) THEN
          IREF(NV) = 28
        ELSEIF( INDEX(ADUM(1:),'napl hydraulic head').NE.0 ) THEN
          IREF(NV) = 29
        ELSEIF( INDEX(ADUM(1:),'rock/soil type').NE.0 ) THEN
          IREF(NV) = 30
        ELSEIF( INDEX(ADUM(1:),'aqueous relative perm').NE.0 ) THEN
          IREF(NV) = 31
        ELSEIF( INDEX(ADUM(1:),'gas relative perm').NE.0 ) THEN
          IREF(NV) = 32
        ELSEIF( INDEX(ADUM(1:),'napl relative perm').NE.0 ) THEN
          IREF(NV) = 33
        ELSEIF( INDEX(ADUM(1:),'aqueous density').NE.0 ) THEN
          IREF(NV) = 34
        ELSEIF( INDEX(ADUM(1:),'gas density').NE.0 ) THEN
          IREF(NV) = 35
        ELSEIF( INDEX(ADUM(1:),'napl density').NE.0 ) THEN
          IREF(NV) = 36
        ELSEIF( INDEX(ADUM(1:),'total water mass').NE.0 ) THEN
          IREF(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total air mass').NE.0 ) THEN
          IREF(NV) = 38
        ELSEIF( INDEX(ADUM(1:),'total oil mass').NE.0 ) THEN
          IREF(NV) = 39
        ELSEIF( INDEX(ADUM(1:),'water mass source int').NE.0 ) THEN
          IREF(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'air mass source int').NE.0 ) THEN
          IREF(NV) = 41
        ELSEIF( INDEX(ADUM(1:),'oil mass source int').NE.0 ) THEN
          IREF(NV) = 42
        ELSEIF( INDEX(ADUM(1:),'aqueous courant').NE.0 ) THEN
          ICRNT = 1
          IREF(NV) = 49
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
        ELSEIF( INDEX(ADUM(1:),'x napl vol').NE.0 ) THEN
          IREF(NV) = 57
        ELSEIF( INDEX(ADUM(1:),'y napl vol').NE.0 ) THEN
          IREF(NV) = 58
        ELSEIF( INDEX(ADUM(1:),'z napl vol').NE.0 ) THEN
          IREF(NV) = 59
        ELSEIF( INDEX(ADUM(1:),'water gas mole').NE.0 ) THEN
          IREF(NV) = 70
        ELSEIF( INDEX(ADUM(1:),'air gas mole').NE.0 ) THEN
          IREF(NV) = 71
        ELSEIF( INDEX(ADUM(1:),'oil gas mole').NE.0 ) THEN
          IREF(NV) = 72
        ELSEIF( INDEX(ADUM(1:),'water gas conc').NE.0 ) THEN
          IREF(NV) = 73
        ELSEIF( INDEX(ADUM(1:),'air gas conc').NE.0 ) THEN
          IREF(NV) = 74
        ELSEIF( INDEX(ADUM(1:),'oil gas conc').NE.0 ) THEN
          IREF(NV) = 75
        ELSEIF( INDEX(ADUM(1:),'water aqueous conc').NE.0 ) THEN
          IREF(NV) = 76
        ELSEIF( INDEX(ADUM(1:),'air aqueous conc').NE.0 ) THEN
          IREF(NV) = 77
        ELSEIF( INDEX(ADUM(1:),'oil aqueous conc').NE.0 ) THEN
          IREF(NV) = 78
        ELSEIF( INDEX(ADUM(1:),'gas courant').NE.0 ) THEN
          ICRNT = 1
          IREF(NV) = 79
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
        ELSEIF( INDEX(ADUM(1:),'xnc napl vol').NE.0 ) THEN
          IREF(NV) = 93
        ELSEIF( INDEX(ADUM(1:),'ync napl vol').NE.0 ) THEN
          IREF(NV) = 94
        ELSEIF( INDEX(ADUM(1:),'znc napl vol').NE.0 ) THEN
          IREF(NV) = 95
        ELSEIF( INDEX(ADUM(1:),'napl courant').NE.0 ) THEN
          ICRNT = 1
          IREF(NV) = 99
        ELSEIF( INDEX(ADUM(1:),'water mass source rate').NE.0 ) THEN
          IREF(NV) = 140
        ELSEIF( INDEX(ADUM(1:),'air mass source rate').NE.0 ) THEN
          IREF(NV) = 141
        ELSEIF( INDEX(ADUM(1:),'oil mass source rate').NE.0 ) THEN
          IREF(NV) = 142
        ELSEIF( INDEX(ADUM(1:),'integrated water mass').NE.0 ) THEN
          IREF(NV) = 191
        ELSEIF( INDEX(ADUM(1:),'integrated air mass').NE.0 ) THEN
          IREF(NV) = 192
        ELSEIF( INDEX(ADUM(1:),'integrated oil mass').NE.0 ) THEN
          IREF(NV) = 193
        ELSEIF( INDEX(ADUM(1:),'integrated aqueous water').NE.0 ) THEN
          IREF(NV) = 194
        ELSEIF( INDEX(ADUM(1:),'integrated aqueous air').NE.0 ) THEN
          IREF(NV) = 195
        ELSEIF( INDEX(ADUM(1:),'integrated aqueous oil').NE.0 ) THEN
          IREF(NV) = 196
        ELSEIF( INDEX(ADUM(1:),'integrated gas water').NE.0 ) THEN
          IREF(NV) = 197
        ELSEIF( INDEX(ADUM(1:),'integrated gas air').NE.0 ) THEN
          IREF(NV) = 198
        ELSEIF( INDEX(ADUM(1:),'integrated gas oil').NE.0 ) THEN
          IREF(NV) = 199
        ELSEIF( INDEX(ADUM(1:),'sorbed oil mass frac').NE.0 ) THEN
          IREF(NV) = 241
        ELSEIF( INDEX(ADUM(1:),'sorbed oil mass').NE.0 ) THEN
          IREF(NV) = 240
        ELSEIF( INDEX(ADUM(1:),'sorbed oil volumetric conc').NE.0 ) THEN
          IREF(NV) = 242
        ELSEIF( INDEX(ADUM(1:),'solute volumetric conc').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 1
        ELSEIF( INDEX(ADUM(1:),'solute aqueous conc').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 2
        ELSEIF( INDEX(ADUM(1:),'solute gas conc').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 3
        ELSEIF( INDEX(ADUM(1:),'solute napl conc').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 4
        ELSEIF( INDEX(ADUM(1:),'solute aqueous mol').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 5
        ELSEIF( INDEX(ADUM(1:),'solute gas mol').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 6
        ELSEIF( INDEX(ADUM(1:),'solute napl mol').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 7
        ELSEIF( INDEX(ADUM(1:),'x solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 8
        ELSEIF( INDEX(ADUM(1:),'y solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 9
        ELSEIF( INDEX(ADUM(1:),'z solute flux').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 10
        ELSEIF( INDEX(ADUM(1:),'solute source').NE.0 ) THEN
          IREF(NV) = 400 + (NSL-1)*33 + 11
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Reference Node Variable: '//ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Check for duplicate reference node variables  ---
!
        DO 190 NX = 1,NV-1
          IF( IREF(NV).EQ.IREF(NX) ) THEN
            INDX = 4
            CHMSG = 'Duplicate Reference Node Variable: '//ADUM
            CALL WRMSGS( INDX )
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
   97 READ(IRD,'(A)') CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 97
      ISTART = 1
      VARB = 'Number of Plot File Variables'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NVPLOT)
      NVC = 0
      DO 400 NV = 1,NVPLOT
   98   READ(IRD,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 98
        CALL LCASE( CHDUM )
        ISTART=1
        VARB = 'Plot File Variable'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
          VARB = 'Plot File Variable: Solute Name'
          CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SOLNM)
          DO 310 NSL = 1,NSOLU
            IF( SOLNM.EQ.SOLUT(NSL) ) GOTO 320
  310     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Solute: '//SOLNM(1:NCS)
          CALL WRMSGS( INDX )
  320     CONTINUE
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
        ELSEIF( INDEX(ADUM(1:),'napl pressure').NE.0 ) THEN
          IPLOT(NV) = 3
        ELSEIF( INDEX(ADUM(1:),'temperature').NE.0 ) THEN
          IPLOT(NV) = 4
        ELSEIF( INDEX(ADUM(1:),'phase condition').NE.0 ) THEN
          IPLOT(NV) = 5
        ELSEIF( INDEX(ADUM(1:),'aqueous gauge pressure').NE.0 ) THEN
          IPLOT(NV) = 6
        ELSEIF( INDEX(ADUM(1:),'gas gauge pressure').NE.0 ) THEN
          IPLOT(NV) = 7
        ELSEIF( INDEX(ADUM(1:),'napl gauge pressure').NE.0 ) THEN
          IPLOT(NV) = 8
        ELSEIF( INDEX(ADUM(1:),'apparent aqueous sat').NE.0 ) THEN
          IPLOT(NV) = 9
        ELSEIF( INDEX(ADUM(1:),'apparent total sat').NE.0 ) THEN
          IPLOT(NV) = 10
        ELSEIF( INDEX(ADUM(1:),'effective trapped napl').NE.0 ) THEN
          IPLOT(NV) = 18
        ELSEIF( INDEX(ADUM(1:),'effective trapped gas').NE.0 ) THEN
          IPLOT(NV) = 19
        ELSEIF( INDEX(ADUM(1:),'aqueous trapped gas').NE.0 ) THEN
          IPLOT(NV) = 107
        ELSEIF( INDEX(ADUM(1:),'napl trapped gas').NE.0 ) THEN
          IPLOT(NV) = 108
        ELSEIF( INDEX(ADUM(1:),'trapped gas sat').NE.0 ) THEN
          IPLOT(NV) = 105
        ELSEIF( INDEX(ADUM(1:),'trapped napl sat').NE.0 ) THEN
          IPLOT(NV) = 106
        ELSEIF( INDEX(ADUM(1:),'aqueous saturation').NE.0 ) THEN
          IPLOT(NV) = 11
        ELSEIF( INDEX(ADUM(1:),'gas saturation').NE.0 ) THEN
          IPLOT(NV) = 12
        ELSEIF( INDEX(ADUM(1:),'napl saturation').NE.0 ) THEN
          IPLOT(NV) = 13
        ELSEIF( INDEX(ADUM(1:),'total saturation').NE.0 ) THEN
          IPLOT(NV) = 14
        ELSEIF( INDEX(ADUM(1:),'aqueous moisture cont').NE.0 ) THEN
          IPLOT(NV) = 15
        ELSEIF( INDEX(ADUM(1:),'napl moisture cont').NE.0 ) THEN
          IPLOT(NV) = 16
        ELSEIF( INDEX(ADUM(1:),'total moisture cont').NE.0 ) THEN
          IPLOT(NV) = 17
        ELSEIF( INDEX(ADUM(1:),'diffusive porosity').NE.0 ) THEN
          IPLOT(NV) = 20
        ELSEIF( INDEX(ADUM(1:),'water gas mass frac').NE.0 ) THEN
          IPLOT(NV) = 21
        ELSEIF( INDEX(ADUM(1:),'air gas mass frac').NE.0 ) THEN
          IPLOT(NV) = 22
        ELSEIF( INDEX(ADUM(1:),'oil gas mass frac').NE.0 ) THEN
          IPLOT(NV) = 23
        ELSEIF( INDEX(ADUM(1:),'water aqueous mass frac').NE.0 ) THEN
          IPLOT(NV) = 24
        ELSEIF( INDEX(ADUM(1:),'air aqueous mass frac').NE.0 ) THEN
          IPLOT(NV) = 25
        ELSEIF( INDEX(ADUM(1:),'oil aqueous mass frac').NE.0 ) THEN
          IPLOT(NV) = 26
        ELSEIF( INDEX(ADUM(1:),'aqueous hydraulic head').NE.0 ) THEN
          IPLOT(NV) = 27
        ELSEIF( INDEX(ADUM(1:),'gas hydraulic head').NE.0 ) THEN
          IPLOT(NV) = 28
        ELSEIF( INDEX(ADUM(1:),'napl hydraulic head').NE.0 ) THEN
          IPLOT(NV) = 29
        ELSEIF( INDEX(ADUM(1:),'rock/soil type').NE.0 ) THEN
          IPLOT(NV) = 30
        ELSEIF( INDEX(ADUM(1:),'aqueous relative perm').NE.0 ) THEN
          IPLOT(NV) = 31
        ELSEIF( INDEX(ADUM(1:),'gas relative perm').NE.0 ) THEN
          IPLOT(NV) = 32
        ELSEIF( INDEX(ADUM(1:),'napl relative perm').NE.0 ) THEN
          IPLOT(NV) = 33
        ELSEIF( INDEX(ADUM(1:),'aqueous density').NE.0 ) THEN
          IPLOT(NV) = 34
        ELSEIF( INDEX(ADUM(1:),'gas density').NE.0 ) THEN
          IPLOT(NV) = 35
        ELSEIF( INDEX(ADUM(1:),'napl density').NE.0 ) THEN
          IPLOT(NV) = 36
        ELSEIF( INDEX(ADUM(1:),'total water mass').NE.0 ) THEN
          IPLOT(NV) = 37
        ELSEIF( INDEX(ADUM(1:),'total air mass').NE.0 ) THEN
          IPLOT(NV) = 38
        ELSEIF( INDEX(ADUM(1:),'total oil mass').NE.0 ) THEN
          IPLOT(NV) = 39
        ELSEIF( INDEX(ADUM(1:),'water mass source int').NE.0 ) THEN
          IPLOT(NV) = 40
        ELSEIF( INDEX(ADUM(1:),'air mass source int').NE.0 ) THEN
          IPLOT(NV) = 41
        ELSEIF( INDEX(ADUM(1:),'oil mass source int').NE.0 ) THEN
          IPLOT(NV) = 42
        ELSEIF( INDEX(ADUM(1:),'aqueous courant').NE.0 ) THEN
          ICRNT = 1
          IPLOT(NV) = 49
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
        ELSEIF( INDEX(ADUM(1:),'x napl vol').NE.0 ) THEN
          IPLOT(NV) = 57
        ELSEIF( INDEX(ADUM(1:),'y napl vol').NE.0 ) THEN
          IPLOT(NV) = 58
        ELSEIF( INDEX(ADUM(1:),'z napl vol').NE.0 ) THEN
          IPLOT(NV) = 59
        ELSEIF( INDEX(ADUM(1:),'water gas mole').NE.0 ) THEN
          IPLOT(NV) = 70
        ELSEIF( INDEX(ADUM(1:),'air gas mole').NE.0 ) THEN
          IPLOT(NV) = 71
        ELSEIF( INDEX(ADUM(1:),'oil gas mole').NE.0 ) THEN
          IPLOT(NV) = 72
        ELSEIF( INDEX(ADUM(1:),'water gas conc').NE.0 ) THEN
          IPLOT(NV) = 73
        ELSEIF( INDEX(ADUM(1:),'air gas conc').NE.0 ) THEN
          IPLOT(NV) = 74
        ELSEIF( INDEX(ADUM(1:),'oil gas conc').NE.0 ) THEN
          IPLOT(NV) = 75
        ELSEIF( INDEX(ADUM(1:),'water aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 76
        ELSEIF( INDEX(ADUM(1:),'air aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 77
        ELSEIF( INDEX(ADUM(1:),'oil aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 78
        ELSEIF( INDEX(ADUM(1:),'gas courant').NE.0 ) THEN
          ICRNT = 1
          IPLOT(NV) = 79
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
        ELSEIF( INDEX(ADUM(1:),'xnc napl vol').NE.0 ) THEN
          IPLOT(NV) = 93
        ELSEIF( INDEX(ADUM(1:),'ync napl vol').NE.0 ) THEN
          IPLOT(NV) = 94
        ELSEIF( INDEX(ADUM(1:),'znc napl vol').NE.0 ) THEN
          IPLOT(NV) = 95
        ELSEIF( INDEX(ADUM(1:),'napl courant').NE.0 ) THEN
          ICRNT = 1
          IPLOT(NV) = 99
        ELSEIF( INDEX(ADUM(1:),'water mass source rate').NE.0 ) THEN
          IPLOT(NV) = 140
        ELSEIF( INDEX(ADUM(1:),'air mass source rate').NE.0 ) THEN
          IPLOT(NV) = 141
        ELSEIF( INDEX(ADUM(1:),'oil mass source rate').NE.0 ) THEN
          IPLOT(NV) = 142
        ELSEIF( INDEX(ADUM(1:),'sorbed oil mass frac').NE.0 ) THEN
          IPLOT(NV) = 241
        ELSEIF( INDEX(ADUM(1:),'sorbed oil mass').NE.0 ) THEN
          IPLOT(NV) = 240
        ELSEIF( INDEX(ADUM(1:),'sorbed oil volumetric conc').NE.0 ) THEN
          IPLOT(NV) = 242
        ELSEIF( INDEX(ADUM(1:),'solute volumetric conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 1
        ELSEIF( INDEX(ADUM(1:),'solute aqueous conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 2
        ELSEIF( INDEX(ADUM(1:),'solute gas conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 3
        ELSEIF( INDEX(ADUM(1:),'solute napl conc').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 4
        ELSEIF( INDEX(ADUM(1:),'solute aqueous mol').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 5
        ELSEIF( INDEX(ADUM(1:),'solute gas mol').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 6
        ELSEIF( INDEX(ADUM(1:),'solute napl mol').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 7
        ELSEIF( INDEX(ADUM(1:),'x solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 8
        ELSEIF( INDEX(ADUM(1:),'y solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 9
        ELSEIF( INDEX(ADUM(1:),'z solute flux').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 10
        ELSEIF( INDEX(ADUM(1:),'solute source').NE.0 ) THEN
          IPLOT(NV) = 400 + (NSL-1)*33 + 11
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Plot File Variable: '//ADUM(1:NCH)
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
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDOU5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSF5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Read input file surface flux information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     Last Modified by MD White, PNNL, 24 June 2002.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/RDSF5'
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
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'aqueous').NE.0) THEN
          IF( INDEX(ADUM(1:),'mass').NE.0) THEN
            ISFT(NS) = 5
            WRITE(IWR,'(A)') 'Aqueous-Phase Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 2
            WRITE(IWR,'(A)') 'Aqueous-Phase Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'gas').NE.0) THEN
          IF( INDEX(ADUM(1:),'oil').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'advective').NE.0 ) THEN
              ISFT(NS) = 31
              WRITE(IWR,'(A)') 'Gas-Advective Oil Mass Flux Surface'
              UNSF(1,NS) = 'kg/s'
              UNSF(2,NS) = 'kg'
            ELSEIF( INDEX(ADUM(1:),'diffusive').NE.0 ) THEN
              ISFT(NS) = 32
              WRITE(IWR,'(A)') 'Gas-Diffusive Oil Mass Flux Surface'
              UNSF(1,NS) = 'kg/s'
              UNSF(2,NS) = 'kg'
            ELSEIF( INDEX(ADUM(1:),'total').NE.0 ) THEN
              ISFT(NS) = 33
              WRITE(IWR,'(A)') 'Gas-Total Oil Mass Flux Surface'
              UNSF(1,NS) = 'kg/s'
              UNSF(2,NS) = 'kg'
            ELSE
              INDX = 4
              CHMSG = 'Unrecognized Surface Flux Type: '//ADUM(1:NCHA)
              CALL WRMSGS( INDX )
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'mass').NE.0) THEN
            ISFT(NS) = 6
            WRITE(IWR,'(A)') 'Gas-Phase Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 3
            WRITE(IWR,'(A)') 'Gas-Phase Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'napl').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'mass').NE.0) THEN
            ISFT(NS) = 7
            WRITE(IWR,'(A)') 'NAPL Mass Flux Surface'
            UNSF(1,NS) = 'kg/s'
            UNSF(2,NS) = 'kg'
          ELSE
            ISFT(NS) = 4
            WRITE(IWR,'(A)') 'NAPL Volumetric Flux Surface'
            UNSF(1,NS) = 'm^3/s'
            UNSF(2,NS) = 'm^3'
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'dissolved oil').NE.0 .OR.
     &    (INDEX(ADUM(1:),'aqueous').NE.0 .AND.
     &    (INDEX(ADUM(1:),'oil').NE.0)) ) THEN
          ISFT(NS) = 9
          WRITE(IWR,'(A)') 'Aqueous Oil Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'water mass').NE.0) THEN
          ISFT(NS) = 10
          WRITE(IWR,'(A)') 'Water Mass Flux Surface'
          UNSF(1,NS) = 'kg/s'
          UNSF(2,NS) = 'kg'
        ELSEIF( INDEX(ADUM(1:),'solute').NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
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
            CHMSG = 'Unrecognized Solute: '//BDUM(1:NCHB)
            CALL WRMSGS( INDX )
   20     CONTINUE
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Surface Flux Type: '//ADUM(1:NCHA)
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
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
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
        WRITE(IWR,'(2A)') VARB(1:IVR),': '
        WRITE (IWR,'(2X,2(A,I6))') 'I = ',ISFC(1,NS),' to ',ISFC(2,NS)
        WRITE (IWR,'(2X,2(A,I6))') 'J = ',ISFC(3,NS),' to ',ISFC(4,NS)
        WRITE (IWR,'(2X,2(A,I6))') 'K = ',ISFC(5,NS),' to ',ISFC(6,NS)
  100 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSF5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSP5

!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Read input file for rock/soil saturation function information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, October 20, 1997.
!     Last Modified by MD White, PNNL, 5 October 2001.
!     Last Modified by MD White, PNNL, 24 June 2002.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!      REAL*8 DRCOSX(2),DRSINX(2),RCOSX(3),RSINX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSP5'
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
      VARB ='NAPL-Aqueous Interfacial Tension'
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
      VARB ='Gas-NAPL Interfacial Tension'
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
        CHMSG = 'Gas-Aqu. Surface Tension = NAPL-Aqu. Surface Tension'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Check for zero gas-NAPL interfacial tension  ---
!
      IF( SCHR(7,1).LT.EPSL ) THEN
        INDX = 4
        CHMSG = 'Zero Gas-NAPL Interfacial Tension'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Write spreading coefficient  ---
!
      C_SPRD = MIN( (SCHR(5,1)-SCHR(6,1)-SCHR(7,1)),0.D+0 )
      WRITE(IWR,'(2X,A,1PE11.4)') 'NAPL Spreading Coefficient: ',C_SPRD
      IF( C_SPRD.LT.-EPSL ) THEN
!
!---    Interfacial contact angles  ---
!
!        NC = 0
!        CA_NL = GPI/3.D+0
!        CA_GN = GPI/3.D+0
!   20   CONTINUE
!        DO 30 M = 1,3
!          CA_OWX = CA_NL
!          CA_AOX = CA_GN
!          IF( M.EQ.2 ) CA_OWX = CA_NL + 1.D-6
!          IF( M.EQ.3 ) CA_AOX = CA_GN + 1.D-6
!          RCOSX(M) = SCHR(5,1) - SCHR(6,1)*COS(CA_OWX) -
!     &      SCHR(7,1)*COS(CA_AOX)
!          RSINX(M) = SCHR(6,1)*SIN(CA_OWX) - SCHR(7,1)*SIN(CA_AOX)
!   30   CONTINUE
!        DRCOSX(1) = 1.D+6*(RCOSX(2)-RCOSX(1))
!        DRCOSX(2) = 1.D+6*(RCOSX(3)-RCOSX(1))
!        DRSINX(1) = 1.D+6*(RSINX(2)-RSINX(1))
!        DRSINX(2) = 1.D+6*(RSINX(3)-RSINX(1))
!        DCA_AOX = (DRSINX(1)*RCOSX(1) - DRCOSX(1)*RSINX(1))/
!     &    (DRSINX(2)*DRCOSX(1)-DRSINX(1)*DRCOSX(2))
!        DCA_OWX = (DRCOSX(2)*RSINX(1) - DRSINX(2)*RCOSX(1))/
!     &    (DRSINX(2)*DRCOSX(1)-DRSINX(1)*DRCOSX(2))
!        CA_GN = CA_GN + DCA_AOX
!        CA_NL = CA_NL + DCA_OWX
!        NC = NC + 1
!        IF( NC.GT.32 ) THEN
!          INDX = 3
!          CHMSG = 'Unconverged Interfacial Contact Angle:'
!          CALL WRMSGS( INDX )
!        ENDIF
!        IF( ABS(DCA_AOX).GT.1.D-9 .OR. ABS(DCA_OWX).GT.1.D-9 ) GOTO 20
        CA_NL = ACOS( ((SCHR(6,1)**2) - (SCHR(7,1)**2) + (SCHR(5,1)**2))
     &    /(2.D+0*SCHR(5,1)*SCHR(6,1)) )
        CA_GN = ACOS( ((SCHR(7,1)**2) - (SCHR(6,1)**2) + (SCHR(5,1)**2))
     &    /(2.D+0*SCHR(5,1)*SCHR(7,1)) )
!
!---    Write gas-NAPL and NAPL-aqueous contact angles  ---
!

        WRITE(IWR,'(2X,A,1PE11.4)') 'Gas-NAPL Contact Angle, deg: ',
     &    (CA_GN*1.8D+2/GPI)
        WRITE(IWR,'(2X,A,1PE11.4)') 'NAPL-Aqueous Contact Angle, deg: ',
     &    (CA_NL*1.8D+2/GPI)
        CA_GN = COS(CA_GN)
        CA_NL = COS(CA_NL)
      ENDIF
!
!---  Assign interfacial tensions to rock/soil characteristic
!     parameters  ---
!
      DO 100 N = 2,NROCK
        SCHR(5,N) = SCHR(5,1)
        SCHR(6,N) = SCHR(6,1)
        SCHR(7,N) = SCHR(7,1)
 100  CONTINUE
!
!---  Loop over the rock/soil saturation information lines  ---
!
      N = 0
      IJK = 0
  110 CONTINUE
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
        GOTO 320
      ENDIF
!
!---  Search known rock types for a matching type ---
!
      DO 200 M = 1, NROCK
        IF( RDUM.EQ.ROCK(M)) THEN
          IROCK = M
          GOTO 300
        ENDIF
  200 CONTINUE
!
!---  Search known scaling groups for a matching type ---
!
      IF( ISLC(19).EQ.1 ) THEN
        DO 210 M = 1,NSCALE
           IF( RDUM.EQ.SCALNM(M) ) THEN
              ISGRP = M
              IROCK = 1
              GOTO 300
           ENDIF
  210   CONTINUE
        INDX = 2
        CHMSG = 'Unrecognized Rock/Soil Type or Scaling Group: '
     &    // RDUM(1:NCH)
        CALL WRMSGS( INDX )
        GOTO 110
      ENDIF
      INDX = 2
      CHMSG = 'Unrecognized Rock/Soil Type: ' // RDUM(1:NCH)
      CALL WRMSGS( INDX )
      GOTO 110
  300 CONTINUE
!
!---  Loop over rock/soils within scaling group  ---
!
      IF( ISLC(19).EQ.1 .AND. ISGRP.NE.0 ) THEN
        DO 302 M = IROCK,NROCK
          IF( ISCALE(M).EQ.ISGRP ) THEN
            IROCK = M
            GOTO 304
          ENDIF
  302   CONTINUE
      ENDIF
  304 CONTINUE
!
!---  Write rock/soil name  ---
!
      WRITE (IWR,'(/,2A)') 'Rock/Soil Name: ',ROCK(IROCK)
      N = N + 1
  320 CONTINUE
!
!---  Read saturation/capillary pressure function  ---
!
        VARB = 'Saturation Function Type'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
        IF( INDEX(ADUM(1:),'entrap').NE.0 .AND.
     &    INDEX(ADUM(1:),'resid').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            ISCHRX = 34
          ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &      INDEX(ADUM(1:),'corey').NE.0 ) THEN
            ISCHRX = 38
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Saturation Function: '//ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
          IF( INDEX(ADUM(1:),'geel').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'linear').NE.0 ) THEN
              IF( INDEX(ADUM(1:),'mod').NE.0 ) THEN
                ISNR = 2
              ELSE
                ISNR = 1
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'land').NE.0 ) THEN
              IF( INDEX(ADUM(1:),'mod').NE.0 ) THEN
                ISNR = 4
              ELSE
                ISNR = 3
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'lenhard').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'mod').NE.0 ) THEN
              ISNR = 6
            ELSE
              ISNR = 5
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'constant').NE.0 ) THEN
            ISNR = 7
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'entrap').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            ISCHRX = 32
          ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &      INDEX(ADUM(1:),'corey').NE.0 ) THEN
            ISCHRX = 36
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Saturation Function: '//ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'resid').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            ISCHRX = 33
          ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &      INDEX(ADUM(1:),'corey').NE.0 ) THEN
            ISCHRX = 37
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Saturation Function: '//ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
          ISNR = 0
          IF( INDEX(ADUM(1:),'geel').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'linear').NE.0 ) THEN
              IF( INDEX(ADUM(1:),'mod').NE.0 ) THEN
                ISNR = 2
              ELSE
                ISNR = 1
              ENDIF
            ELSEIF( INDEX(ADUM(1:),'land').NE.0 ) THEN
              IF( INDEX(ADUM(1:),'mod').NE.0 ) THEN
                ISNR = 4
              ELSE
                ISNR = 3
              ENDIF
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'lenhard').NE.0 ) THEN
            IF( INDEX(ADUM(1:),'mod').NE.0 ) THEN
              ISNR = 6
            ELSE
              ISNR = 5
            ENDIF
          ELSEIF( INDEX(ADUM(1:),'constant').NE.0 ) THEN
            ISNR = 7
          ENDIF
        ELSE
          IF( INDEX(ADUM(1:),'van genuchten').NE.0 ) THEN
            ISCHRX = 31
          ELSEIF( INDEX(ADUM(1:),'brooks').NE.0 .AND.
     &      INDEX(ADUM(1:),'corey').NE.0 ) THEN
            ISCHRX = 35
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Saturation Function: '//ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
        IF( IJK.GT.0 ) THEN
          DO 330 N = 1,NFLD
            ISCHR(IZ(N)) = ISCHRX
            IF( ISCHRX.EQ.7 ) IRPL(IZ(N)) = 7
  330     CONTINUE
        ELSE
          ISCHR(IROCK) = ISCHRX
          IF( ISCHRX.EQ.7 ) IRPL(IROCK) = 7
        ENDIF
!
!---    van Genuchten Function  ---
!
        IF( ISCHRX.GE.31 .AND. ISCHRX.LE.34 ) THEN
          IF( ISCHRX.EQ.31 ) THEN
            WRITE(IWR,'(A)') 'van Genuchten s-P Function'
          ELSEIF( ISCHRX.EQ.32 ) THEN
            WRITE(IWR,'(2A)') 'van Genuchten w/ NAPL Entrapment',
     &      ' s-P Function'
          ELSEIF( ISCHRX.EQ.33 ) THEN
            WRITE(IWR,'(2A)') 'van Genuchten w/ NAPL Residual',
     &      ' s-P Function'
          ELSEIF( ISCHRX.EQ.34 ) THEN
            WRITE(IWR,'(2A)') 'van Genuchten w/ Combined',
     &      ' NAPL Entrapment and NAPL Residual s-P Function'
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
     &        ': ',SCHR(1,IROCK)
            INDX = 0
            IUNM = -1
            CALL RDUNIT(UNTS,SCHR(1,IROCK),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(1,IROCK),', 1/m)'
          ENDIF
          SCHR(2,IROCK) = 0.D+0
          VARB = 'van Genuchten (n)'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(3,IROCK))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(3,IROCK)
          VARB = 'van Genuchten (minimum saturation)'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(4,IROCK))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(4,IROCK)
          VARB = 'van Genuchten (Reference Interfacial Tension)'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(8,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHR(8,IROCK)
          INDX = 0
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,SCHR(8,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(8,IROCK),', N/m)'
          SCHR(9,IROCK) = 1.D+20
          VARB = 'van Genuchten (m)'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(14,IROCK))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(14,IROCK)
          IF( ISCHRX.EQ.32 ) THEN
            VARB = 'van Genuchten (Effective Maximum Trapped NAPL ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(11,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(11,IROCK)
            VARB = 'van Genuchten (Critical Trapping Number)'
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(9,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(9,IROCK)
            IF( SCHR(9,IROCK)/EPSL.LT.EPSL ) THEN
              INDX = 4
              CHMSG = 'Negative or Zero Critical Trapping Number'
              CALL WRMSGS( INDX )
            ENDIF
            VARB = 'van Genuchten (Effective Maximum Trapped Gas ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(17,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(17,IROCK)
          ELSEIF( ISCHRX.EQ.33 ) THEN
            VARB = 'van Genuchten (Effective Maximum Residual NAPL ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(10,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(10,IROCK)
          ELSEIF( ISCHRX.EQ.34 ) THEN
            VARB = 'van Genuchten (Effective Maximum Trapped NAPL ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(11,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(11,IROCK)
            VARB = 'van Genuchten (Effective Maximum Residual NAPL ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(10,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(10,IROCK)
            VARB = 'van Genuchten (Critical Trapping Number)'
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(9,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(9,IROCK)
            IF( SCHR(9,IROCK)/EPSL.LT.EPSL ) THEN
              INDX = 4
              CHMSG = 'Negative or Zero Critical Trapping Number'
              CALL WRMSGS( INDX )
            ENDIF
            VARB = 'van Genuchten (Effective Maximum Trapped Gas ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(17,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(17,IROCK)
          ENDIF
!
!---    Brooks and Corey Function  ---
!
        ELSEIF( ISCHRX.GE.35 .AND. ISCHRX.LE.38 ) THEN
          IF( ISCHRX.EQ.35 ) THEN
            WRITE(IWR,'(A)') 'Brooks and Corey s-P Function'
          ELSEIF( ISCHRX.EQ.36 ) THEN
            WRITE(IWR,'(2A)') 'Brooks and Corey w/ NAPL Entrapment',
     &        ' s-P Function'
          ELSEIF( ISCHRX.EQ.37 ) THEN
            WRITE(IWR,'(2A)') 'Brooks and Corey w/ NAPL Residual',
     &        ' s-P Function'
          ELSEIF( ISCHRX.EQ.38 ) THEN
            WRITE(IWR,'(2A)') 'Brooks and Corey w/ Combined',
     &      ' NAPL Entrapment and NAPL Residual s-P Function'
          ENDIF
          VARB = 'Brooks-Corey (entry head)'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(2,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHR(2,IROCK)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SCHR(2,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(2,IROCK),', m)'
          VARB = 'Brooks-Corey (lambda parameter)'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(3,IROCK))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(3,IROCK)
          VARB = 'Brooks-Corey (residual saturation)'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(4,IROCK))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(4,IROCK)
          VARB = 'Brooks-Corey (Reference Interfacial Tension)'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(8,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHR(8,IROCK)
          INDX = 0
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,SCHR(8,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(8,IROCK),', N/m)'
          SCHR(9,IROCK) = 1.D+20
          IF( ISCHRX.EQ.36 ) THEN
            VARB = 'Brooks-Corey (Effective Maximum Trapped NAPL ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(11,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(11,IROCK)
            VARB = 'Brooks-Corey (Critical Trapping Number)'
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(9,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(9,IROCK)
            IF( SCHR(9,IROCK)/EPSL.LT.EPSL ) THEN
              INDX = 4
              CHMSG = 'Negative or Zero Critical Trapping Number'
              CALL WRMSGS( INDX )
            ENDIF
            VARB = 'Brooks-Corey (Effective Maximum Trapped Gas ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(17,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(17,IROCK)
          ELSEIF( ISCHRX.EQ.37 ) THEN
            VARB = 'Brooks-Corey (Effective Maximum Resdiual NAPL ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(10,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(10,IROCK)
          ELSEIF( ISCHRX.EQ.38 ) THEN
            VARB = 'Brooks-Corey (Effective Maximum Trapped NAPL ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(11,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(11,IROCK)
            VARB = 'Brooks-Corey (Effective Maximum Residual NAPL ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(10,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(10,IROCK)
            VARB = 'Brooks-Corey (Critical Trapping Number)'
            IDFLT = 1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(9,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(9,IROCK)
            IF( SCHR(9,IROCK)/EPSL.LT.EPSL ) THEN
              INDX = 4
              CHMSG = 'Negative or Zero Critical Trapping Number'
              CALL WRMSGS( INDX )
            ENDIF
            VARB = 'Brooks-Corey (Effective Maximum Trapped Gas ' //
     &        'Saturation)'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(17,IROCK))
            WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(17,IROCK)
          ENDIF
        ENDIF
!
!---  Loop over remaining rock/soils within scaling group  ---
!
      IF( ISLC(19).EQ.1 .AND. IROCK.LT.NROCK ) THEN
        DO 490 M = IROCK+1,NROCK
          IF( ISCALE(M).EQ.ISGRP ) THEN
            N = N+1
            ISCHR(M) = ISCHR(IROCK)
            IF( ISCHR(M).EQ.7 ) IRPL(M) = 7
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
      GOTO 110
 500  CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSP5 group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSR5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Read input file for source information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     Last Modified by MD White, Battelle, 2 October 2003.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 VAR(LSTM,6),SRCPX(8)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSR5'
      EPSLX = 1.D-12
!
!---  Write card information to ouput file  ---
!
      CARD = 'Source Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
      NSR = 0
   90 READ(IRD,'(A)') CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 90
      ISTART = 1
      VARB = 'Number of Sources'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      DO 140 NS = 1, NLIN
   91   READ(IRD,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 91
        CALL LCASE( CHDUM )
        ISTART = 1
!
!---  Read source type  ---
!
        VARB = 'Source Type'
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'aqueous volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Volumetric Source w/ Mass Frac.'
          ISRTX = 3
        ELSEIF( INDEX(ADUM(1:),'gas volu').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'mass frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Volumetric Source w/ Mass Frac.'
            ISRTX = 4
          ELSEIF( INDEX(ADUM(1:),'relative').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Volumetric Source w/ Rel. Humidity'
            ISRTX = 5
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'napl volu').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'NAPL Volumetric Source'
          ISRTX = 6
        ELSEIF( INDEX(ADUM(1:),'aqueous mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Mass Source w/ Mass Frac.'
          ISRTX = 7
        ELSEIF( INDEX(ADUM(1:),'gas mass').NE.0 ) THEN
          IF( INDEX(ADUM(1:),'mass frac').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Mass Source w/ Mass Frac.'
            ISRTX = 8
          ELSEIF( INDEX(ADUM(1:),'relative').NE.0 ) THEN
            WRITE(IWR,'(2X,A)') 'Gas Mass Source w/ Rel. Humidity'
            ISRTX = 9
          ENDIF
        ELSEIF( INDEX(ADUM(1:),'napl mass').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'NAPL Mass Source'
          ISRTX = 10
        ELSEIF( INDEX(ADUM(1:),'sve').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Soil Vapor Extraction Well'
          ISRTX = 22
        ELSEIF( INDEX(ADUM(1:),'bioslurping').NE.0 .AND.
     &    INDEX(ADUM(1:),'well').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Bioslurping Well'
          ISRTX = 30
        ELSEIF( IEQC.NE.0 .AND. INDEX(ADUM(1:),'solute').NE.0 ) THEN
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHB,CHDUM,BDUM)
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
            CHMSG = 'Unrecognized Solute: '//BDUM(1:NCHB)
            CALL WRMSGS( INDX )
   40     CONTINUE
        ELSE
          INDX = 4
          CHMSG = 'Unrecognized Source Type: '//ADUM(1:NCHA)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read source domain indices, unless source type is bioslurping
!       well  ---
!
        IF( ISRTX.NE.30 ) THEN
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
!---      Check for ill-defined source domains  ---
!
          IF( I1X.LT.1 .OR. I1X.GT.IFLD .OR. I2X.LT.1 .OR.
     &      I2X.GT.IFLD .OR. I2X.LT.I1X ) THEN
            INDX = 4
            CHMSG = 'Invalid Source Domain: ' // CHDUM(ISX:ICX)
            CALL WRMSGS( INDX )
          ENDIF
          IF( J1X.LT.1 .OR. J1X.GT.JFLD .OR. J2X.LT.1 .OR.
     &      J2X.GT.JFLD .OR. J2X.LT.J1X ) THEN
            INDX = 4
            CHMSG = 'Invalid Source Domain: ' // CHDUM(ISX:ICX)
            CALL WRMSGS( INDX )
          ENDIF
          IF( K1X.LT.1 .OR. K1X.GT.KFLD .OR. K2X.LT.1 .OR.
     &      K2X.GT.KFLD .OR. K2X.LT.K1X ) THEN
            INDX = 4
            CHMSG = 'Invalid Source Domain: ' // CHDUM(ISX:ICX)
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Check for z-direction well sources  ---
!
          IF( ISRTX.GE.20 .AND. ISRTX.LE.29 ) THEN
            IF( ((I2X-I1X).GT.1) .OR. ((J2X-J1X).GT.1) ) THEN
              INDX = 4
              CHMSG = 'Invalid Well-Domain Indices'
              CALL WRMSGS( INDX )
            ENDIF
            DO 45 M = 1,NS-1
              IF( ISRT(M).GE.20 .AND. ISRT(M).LE.29 ) THEN
                IF( I1X.EQ.ISRDM(1,M) .AND. J1X.EQ.ISRDM(3,M) ) THEN
                  IF( K1X.GE.ISRDM(5,M) .AND. K1X.LE.ISRDM(6,M) ) THEN
                    INDX = 4
                    CHMSG = 'Overlapping Well-Source Domains'
                    CALL WRMSGS( INDX )
                  ELSEIF( K2X.GE.ISRDM(5,M) .AND. 
     &              K2X.LE.ISRDM(6,M) ) THEN
                    INDX = 4
                    CHMSG = 'Overlapping Well-Source Domains'
                    CALL WRMSGS( INDX )
                  ENDIF
                ENDIF
              END IF
   45       CONTINUE
          ENDIF
!
!---    Check for sources applied to inactive nodes  ---
!
          DO 50 K = K1X,K2X
            DO 50 J = J1X,J2X
              DO 50 I = I1X,I2X
                IF( IXP(ND(I,J,K)).EQ.0 ) THEN
                  INDX = 4
                  CHMSG = 'Source Applied to Inactive Node: '//
     &              CHDUM(ISX:ISTART)
                  CALL WRMSGS( INDX )
                ENDIF
   50     CONTINUE
!
!---    Read x,y,z-coordinate location of the bioslurping tube,
!       x,y,z-component lengths of the well screened interval,
!       well radius, and well skin factor  ---
!
        ELSE
          VARB = 'Bioslurping Tube X-Direction Coordinate Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SRCPX(1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SRCPX(1)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SRCPX(1),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SRCPX(1),', m)'
          VARB = 'Bioslurping Tube Y-Direction Coordinate Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SRCPX(2))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SRCPX(2)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SRCPX(2),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SRCPX(2),', m)'
!
!---      Cylindrical coordinates with azimuthal symmetry  ---
!
          IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
            IF( ABS(SRCPX(1))/EPSL.GT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Non-Zero Bioslurping Tube X-Direction ' // 
     &          'Coordinate Location for Radially Symmetric Domain'
              RLMSG = SRCPX(1)
              CALL WRMSGS( INDX )
            ENDIF
            IF( ABS(SRCPX(2))/EPSL.GT.EPSL ) THEN
              INDX = 9
              CHMSG = 'Non-Zero Bioslurping Tube X-Direction ' // 
     &          'Coordinate Location for Radially Symmetric Domain'
              RLMSG = SRCPX(2)
              CALL WRMSGS( INDX )
            ENDIF
          ENDIF
          VARB = 'Bioslurping Tube Z-Direction Coordinate Location'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SRCPX(3))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SRCPX(3)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SRCPX(3),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SRCPX(3),', m)'
          VARB = 'Bioslurping Well Screened Interval X-Component Length'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SRCPX(4))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SRCPX(4)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SRCPX(4),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SRCPX(4),', m)'
          VARB = 'Bioslurping Well Screened Interval Y-Component Length'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SRCPX(5))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SRCPX(5)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SRCPX(5),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SRCPX(5),', m)'
          VARB = 'Bioslurping Well Screened Interval Z-Component Length'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SRCPX(6))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SRCPX(6)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SRCPX(6),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SRCPX(6),', m)'
          VARB = 'Bioslurping Well Radius'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SRCPX(7))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SRCPX(7)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SRCPX(7),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SRCPX(7),', m)'
          VARB = 'Bioslurping Well Skin Factor'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SRCPX(8))
          WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),': ',SRCPX(8)
!
!---      Loop over active nodes to find node number for the 
!         bioslurping tube location ---
!
          DO 70 N = 1,NFLD
            IF( IXP(N).EQ.0 ) GOTO 70
            I = ID(N)
            J = JD(N)
            K = KD(N)
!
!---        Cylindrical coordinates with azimuthal symmetry,
!           centrally located wells  ---
!
            IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. LFY.EQ.1
     &        .AND. I.EQ.1 ) THEN
!
!---          Node height greater than EPSLX  ---
!
              IF( ABS(ZE(1,N)-ZE(5,N)).GT.EPSLX ) THEN
                DZPX1 = SRCPX(3)-ZE(1,N)
                DZPX2 = ZE(5,N)-SRCPX(3)
                IF( ABS(DZPX1).LT.EPSLX ) DZPX1 = 0.D+0
                IF( ABS(DZPX2).LT.EPSLX ) DZPX2 = 0.D+0
!
!---            Transition point within vertical limits of node  ---
!
                IF( DZPX1.GE.0.D+0 .AND. DZPX2.GE.0.D+0 ) GOTO 80
              ENDIF
!
!---        Non-cylindrical coordinates with azimuthal symmetry,
!           centrally located wells  ---
!
            ELSE
!
!---          Check for point with hexahedron  ---
!
              CALL WITHIN( SRCPX(1),SRCPX(2),SRCPX(3),ICWX,N )
!
!---          Opposing rotations found, point outside hexahedron  ---
!
              IF( ICWX.EQ.0 ) GOTO 70
              GOTO 80
            ENDIF
   70     CONTINUE
          INDX = 4
          NCH = INDEX( CHDUM,'  ' ) - 1
          CHMSG = 'Bioslurping Point Outside Active Domain: ' // 
     &      CHDUM(1:NCH)
          CALL WRMSGS( INDX )
   80     CONTINUE
!
!---      No opposing rotations found, point inside hexahedron
!         store and record source node  ---
!
          I1X = ID(N)
          I2X = ID(N)
          J1X = JD(N)
          J2X = JD(N)
          K1X = KD(N)
          K2X = KD(N)
          WRITE(IWR,'(/,2X,A)') 'Source Domain:'
          WRITE(IWR,'(4X,A,I6,A,I6)') 'I = ',I1X,' to ',I2X
          WRITE(IWR,'(4X,A,I6,A,I6)') 'J = ',J1X,' to ',J2X
          WRITE(IWR,'(4X,A,I6,A,I6)') 'K = ',K1X,' to ',K2X
        ENDIF
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
!---    Minimum well gas pressures and well diameters for 
!       soil-vapor-extraction wells  ---
!
        IF( ISRTX.EQ.22 ) THEN
          CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
          IF( INDX.EQ.1 ) THEN
            VARB = 'SVE Well Minimum Gas Pressure: '
            IDFLT = 1
            SRCPX(1) = 0.D+0
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SRCPX(1))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SRCPX(1)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,SRCPX(1),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SRCPX(1),', Pa)'
            VARB = 'SVE Well Diameter: '
            IDFLT = 1
            SRCPX(2) = 1.7D-1
            CALL RDDPR(ISTART,ICOMMA,CHDUM,SRCPX(2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',SRCPX(2)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,SRCPX(2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',SRCPX(2),', m)'
!
!---        Convert minimum well pressure to guage and well diameter
!           to well radius  ---
!
            SRCPX(1) = SRCPX(1)-PATM
            SRCPX(2) = 5.D-1*SRCPX(2)
          ELSE
            SRCPX(1) = -PATM
            SRCPX(2) = 8.5D-2
          ENDIF
        ENDIF
!
!---    Initialize input variables  ---
!
        DO 100 NTM = 1,ISRM(NS)
          DO 60 M = 1,6
            VAR(NTM,M) = 0.D+0
   60     CONTINUE
!
!---  Read and write source values and units  ---
!
   92     READ(IRD,'(A)') CHDUM
          IF( CHDUM(1:1).EQ.'#' ) GOTO 92
          CALL LCASE( CHDUM )
          ISTART = 1
          VARB = 'Source Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,1))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4)')VARB(1:IVR),', ',UNTS(1:NCH),': ',
     &      VAR(NTM,1)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,VAR(NTM,1),INDX)
          IF( ISRTX.EQ.3 ) THEN
            VARB = 'Source Aqueous Volumetric Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &        VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Mass Fraction Dissolved Air: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,5)
            VARB = 'Source Mass Fraction Dissolved Oil: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,6)
          ELSEIF( ISRTX.EQ.4 ) THEN
            VARB = 'Source Gas Volumetric Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &        VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Water Vapor Mass Fraction: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,5)
            VARB = 'Source Oil Vapor Mass Fraction: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,6)
          ELSEIF( ISRTX.EQ.5 ) THEN
            VARB = 'Source Gas Volumetric Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &        VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Water Vapor Relative Humidity: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,5)
            VARB = 'Source Oil Vapor Relative Humidity: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,6)
          ELSEIF( ISRTX.EQ.6 ) THEN
            VARB = 'Source NAPL Volumetric Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ELSEIF( ISRTX.EQ.7 ) THEN
            VARB = 'Source Aqueous Mass Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &        VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Mass Fraction Dissolved Air: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,5)
            VARB = 'Source Mass Fraction Dissolved Oil: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,6)
          ELSEIF( ISRTX.EQ.8 ) THEN
            VARB = 'Source Gas Mass Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &        VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Water Vapor Mass Fraction: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,5)
            VARB = 'Source Oil Vapor Mass Fraction: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,6)
          ELSEIF( ISRTX.EQ.9 ) THEN
            VARB = 'Source Gas Mass Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &        VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            VARB = 'Source Water Vapor Relative Humidity: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,5)
            VARB = 'Source Oil Vapor Relative Humidity: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,6)
          ELSEIF( ISRTX.EQ.10 ) THEN
            VARB = 'Source NAPL Mass Rate: '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(/,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &        VAR(NTM,4)
            INDX = 0
            IUNKG = 1
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
!
!---      Soil-vapor-extraction well  ---
!
          ELSEIF( ISRTX.EQ.22 ) THEN
            VARB = 'Source Volumetric Gas Flow Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,4),', m^3/s)'
            VARB = 'Top of Screened Interval for SVE Well'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,5))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,5)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,5),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,5),', m)'
            VARB = 'Bottom of Screened Interval for SVE Well'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,6))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,6)
            IDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,VAR(NTM,6),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,6),', m)'
!
!---        Well used for injection, read water-vapor and oil-vapor
!           mass fractions  ---
!
            IF( VAR(NTM,4).GT.EPSL ) THEN
              VARB = 'Source Water Vapor Mass Fraction: '
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
              WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,2)
              VARB = 'Source Oil Vapor Mass Fraction: '
              CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
              WRITE(IWR,'(/,A,1PE11.4)') VARB(1:IVR),VAR(NTM,3)
            ENDIF
!
!---        Check screened interval against source domain  ---
!
            IF( (VAR(NTM,5)-VAR(NTM,6)).LT.EPSL ) THEN
              INDX = 4
              CHMSG = 'SVE Screened Interval'
              CALL WRMSGS( INDX )
            ENDIF
            IF( VAR(NTM,5).LT.ZE(1,ND(I1X,J1X,K1X)) .OR.
     &        VAR(NTM,5).GT.ZE(5,ND(I1X,J1X,K2X)) ) THEN
              INDX = 4
              CHMSG = 'SVE Top of Screened Interval'
              CALL WRMSGS( INDX )
            ENDIF
            IF( VAR(NTM,6).LT.ZE(1,ND(I1X,J1X,K1X)) .OR.
     &        VAR(NTM,6).GT.ZE(5,ND(I1X,J1X,K2X)) ) THEN
              INDX = 4
              CHMSG = 'SVE Bottom of Screened Interval'
              CALL WRMSGS( INDX )
            ENDIF
!
!---      Bioslurping well  ---
!
          ELSEIF( ISRTX.EQ.30 ) THEN
            VARB = 'Bioslurping Well: Gas Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,2))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,2)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,VAR(NTM,2),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,2),', Pa)'
            VARB = 'Bioslurping Well: Maximum Liquid Volumetric Flow'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,3))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,3)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,3),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',VAR(NTM,3),', m^3/2)'
         ELSEIF( ISRTX.LT.0 .AND. ISRTX.GE.-NSOLU ) THEN
            VARB = 'Source Solute Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
            INDX = 0
            IUNS = -1
            CALL RDUNIT(UNTS,VAR(NTM,4),INDX)
          ELSEIF( ISRTX.LT.-NSOLU .AND. ISRTX.GE.-2*NSOLU ) THEN
            VARB = 'Source Solute Density Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR(NTM,4))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',VAR(NTM,4)
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
        DO 110 M = 1,8
          SRCP(M,NSR) = SRCPX(M)
  110   CONTINUE
        DO 130 NTM = 1,ISRM(NS)
          DO 120 M = 1,6
            SRC(M,NTM,NSR) = VAR(NTM,M)
  120     CONTINUE
  130   CONTINUE
  140 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSR5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTF5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Reads solute/fluid interaction card for diffusion and partition
!     coefficients, and internodal diffusion term averaging scheme for
!     two phase (aqueous-NAPL) convective-dispersive transport equation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/RDTF5'
!
!---  Write card information to output file  ---
!
      CARD = 'Solute/Fluid Interaction Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of different solutes  ---
!
   92 READ(IRD,'(A)') CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 92
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Solutes: '
      CALL RDINT(ISTART,ICOMMA,CHDUM,NLIN)
      NSOLU = 0
      DO 200 NL = 1, NLIN
   93   READ(IRD,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 93
        CALL LCASE(CHDUM)
        ISTART = 1
        ADUM(1:) = ' '
        VARB = 'Solute Name: '
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        DO 100 NS = 1,NSOLU
          NCHS = INDEX( SOLUT(NS),'  ' ) -1
          IF( SOLUT(NS)(1:NCHS).EQ.ADUM(1:NCHA) ) THEN
            NSL = NS
            GOTO 110
          ENDIF
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
        WRITE(IWR,'(/,2A)') VARB(1:IVR),ADUM(1:NCHA)
!
!---  Aqueous molecular diffusion coefficient  ---
!
        VARB = 'Aqueous Molecular Diffusion Coefficient: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDL(NSL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': '
     &,SMDL(NSL)
        INDX = 0
        IUNM = 2
        IUNS = -1
        CALL RDUNIT(UNTS,SMDL(NSL),INDX)
!
!---  Gas molecular diffusion coefficient  ---
!
        VARB = 'Gas Molecular Diffusion Coefficient: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDG(NSL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': '
     &,SMDG(NSL)
        INDX = 0
        IUNM = 2
        IUNS = -1
        CALL RDUNIT(UNTS,SMDG(NSL),INDX)
!
!---  NAPL-phase molecular diffusion coefficient  ---
!
        VARB = 'NAPL Molecular Diffusion Coefficient: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SMDN(NSL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': '
     &,SMDN(NSL)
        INDX = 0
        IUNM = 2
        IUNS = -1
        CALL RDUNIT(UNTS,SMDN(NSL),INDX)
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
!---  Aqueous-NAPL partition coefficient  ---
!
        IDFLT = 1
        VARB = 'Aqueous-NAPL Partition Coefficient: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PCLN(1,NSL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': '
     &,PCLN(1,NSL)
        INDX = 0
        CALL RDUNIT(UNTS,PCLN(1,NSL),INDX)
        PCLN(1,NSL) = MAX( PCLN(1,NSL),1.D-20 )
        PCLN(1,NSL) = MIN( PCLN(1,NSL),1.D+20 )
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
!---  Electrolyte density option  ---
!
      IF( ISLC(16).EQ.1 ) THEN
        NSL_ELC = 0
        DO 210 MSL = 1,NSOLU
          IF( SOLUT(MSL).EQ.ELC_SOL ) NSL_ELC = MSL
  210   CONTINUE
        IF( NSL_ELC.EQ.0 ) THEN
          NCH = INDEX( ELC_SOL(1:),'  ' )-1
          INDX = 4
          CHMSG = 'Electrolyte Solute Not Listed' // ELC_SOL(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
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
          INDX = 6
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
   94 READ(IRD,'(A)') CHDUM
      IF( CHDUM(1:1).EQ.'#' ) GOTO 94
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
   95   READ(IRD,'(A)') CHDUM
        IF( CHDUM(1:1).EQ.'#' ) GOTO 95
        CALL LCASE(CHDUM)
        ISTART = 1
        ADUM(1:) = ' '
        VARB = 'Parent Solute Name'
        NPSL = 0
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        DO 300 NSL = 1,NSOLU
          NCHS = INDEX( SOLUT(NSL),'  ' ) - 1
          IF( SOLUT(NSL)(1:NCHS).EQ.ADUM(1:NCHA) ) NPSL = NSL
  300   CONTINUE
        IF( NCHEM(NPSL).LE.0 ) THEN
          INDX = 4
          CHMSG = 'Non-parent Reactive Solute: '//ADUM(1:NCHA)
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
          CHMSG = 'Chain Decay Solute: '//ADUM(1:NCHA)//': '
     &      //BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ELSEIF( NPSL.EQ.NDSL ) THEN
          INDX = 4
          CHMSG = 'Parent Solute = Progeny Solute: '//ADUM(1:NCHA)//': '
     &      //BDUM(1:NCHB)
          CALL WRMSGS( INDX )
        ENDIF
        VARB = 'Reaction Number'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NCHM)
        VARB = 'Reaction Stoichiometry (Fraction)'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,RCHDF(NPSL,NDSL,NCHM))
        WRITE(IWR,'(2X,4A,I2,A,1PE11.4)') ADUM(1:NCHA),' to ',
     &    BDUM(1:NCHB),'; Reaction = #',NCHM,
     &    '; Stoichiometry = ',RCHDF(NPSL,NDSL,NCHM)
  400 CONTINUE
#endif
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTF5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTP5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Reads the solute/porous media interaction card for the
!     dispersivities, species half-lives, and partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, September 2, 1993.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/RDTP5'
!
!---  Write card information to output file  ---
!
      CARD = 'Solute/Porous Media Interaction Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
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
        VARB = 'Longitudinal Dispersivity: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DISPL(IROCK))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),', ',
     &DISPL(IROCK)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DISPL(IROCK),INDX)
        IF( DISPL(IROCK).GE.SMALL ) IDISP = 1
!
!---  Transverse dispersivity  ---
!
        VARB = 'Transverse Dispersivity: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,DISPT(IROCK))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),', ',
     &DISPT(IROCK)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,DISPT(IROCK),INDX)
        IF( DISPT(IROCK).GE.SMALL ) IDISP = 1
!
!---  Loop over number of solutes or radionuclides  ---
!
        DO 500 NS = 1,NSOLU
   93     READ(IRD,'(A)') CHDUM
          IF( CHDUM(1:1).EQ.'#' ) GOTO 93
          CALL LCASE( CHDUM )
          ISTART = 1
          VARB = 'Solute Name: '
          CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
!
!---  Search known solutes for matching name  ---
!
          DO 300 NSS = 1,NSOLU
            NCHS = INDEX( SOLUT(NSS),'  ' ) - 1
            IF( ADUM(1:NCHA).EQ.SOLUT(NSS)(1:NCHS)) THEN
              NSL = NSS
              GOTO 400
            ENDIF
  300     CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Solute: '//ADUM(1:NCHA)
            CALL WRMSGS( INDX )
  400     CONTINUE
          WRITE (IWR,'(/,2A)') 'Solute Name:',SOLUT(NSL)
!
!---  Solid-aqueous partition coefficient  ---
!
          IDFLT = 1
          VARB = 'Solid-Aqueous Partition Coefficient'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PCSL(1,IROCK,NSL))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',
     &      PCSL(1,IROCK,NSL)
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
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTP5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDL5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute the maximum relative residuals
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     Last Modified by MD White, PNNL, 27 June 2002.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE OUTPU
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
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
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDL5'
!
!---  Initialize maximum residuals  ---
!
      DO 100 M = 1,ISVC
        RSD(M) = 0.D+0
        NSD(M) = 0
  100 CONTINUE
!
!---  Loop over all nodes
!
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        IZN = IZ(N)
        NMD = IXP(N)
        MWEX = IM(IEQW,NMD)
        MAEX = IM(IEQA,NMD)
        MOEX = IM(IEQO,NMD)
        IF( ISKP(IZN).EQ.1 ) GOTO 200
        XMLOX = PSO(2,N)/HCOW
!
!---    Aqueous w/ dissolved-oil, dissolved-air system
!       prior to iteration  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          RSDX = ABS(BLU(MWEX))/(ABS(PL(2,N))+PATM)
          ACCWX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQW,N))/ACCWX,RSDX )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          PTX = MAX( PL(2,N),PG(2,N),PN(2,N),0.D+0 )+PATM
          XMLAX = PTX/HCAW
          RSDX = ABS(BLU(MAEX))/XMLAX
          ACCAX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGA(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQA,N))/ACCAX,RSDX )
          IF( RSDX.GT.RSD(IEQA) ) THEN
            RSD(IEQA) = RSDX
            NSD(IEQA) = N
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            IEQNOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
            MNOEX = IM(IEQNOX,NMD)
            RSD(IEQNOX) = 0.D+0
            NSD(IEQNOX) = 0
          ELSE
            IEQDOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
          ENDIF
          RSDX = ABS(BLU(MDOEX))/MAX( XMLOX,1.D-9 )
          ACCOX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGO(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQDOX,N))/ACCOX,RSDX )
          IF( RSDX.GT.RSD(IEQDOX) ) THEN
            RSD(IEQDOX) = RSDX
            NSD(IEQDOX) = N
          ENDIF
!
!---    Aqueous w/ dissolved-oil, trapped-gas system
!       prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          RSDX = ABS(BLU(MWEX))/(ABS(PL(2,N))+PATM)
          ACCWX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQW,N))/ACCWX,RSDX )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          RSDX = ABS(BLU(MAEX))
          IF( RSDX.GT.RSD(IEQA) ) THEN
            RSD(IEQA) = RSDX
            NSD(IEQA) = N
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            IEQNOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
            MNOEX = IM(IEQNOX,NMD)
            RSD(IEQNOX) = 0.D+0
            NSD(IEQNOX) = 0
          ELSE
            IEQDOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
          ENDIF
          RSDX = ABS(BLU(MDOEX))/MAX( XMLOX,1.D-9 )
          ACCOX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGO(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(3,N))/ACCOX,RSDX )
          IF( RSDX.GT.RSD(IEQDOX) ) THEN
            RSD(IEQDOX) = RSDX
            NSD(IEQDOX) = N
          ENDIF
!
!---  Aqueous w/ trapped-NAPL, dissolved-air system
!     prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
          RSDX = ABS(BLU(MWEX))/(ABS(PL(2,N))+PATM)
          ACCWX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQW,N))/ACCWX,RSDX )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          PTX = MAX( PL(2,N),PG(2,N),PN(2,N),0.D+0 )+PATM
          XMLAX = PTX/HCAW
          RSDX = 1.D-2*ABS(BLU(MAEX))/XMLAX
          IF( RSDX.GT.RSD(IEQA) ) THEN
            RSD(IEQA) = RSDX
            NSD(IEQA) = N
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            IEQNOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
            MNOEX = IM(IEQNOX,NMD)
            RSDX = ABS(BLU(MDOEX))/MAX( XMLOX,1.D-9 )
            ACCOX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &        RHOG(2,N)*SG(2,N)*XGO(2,N))*DTI*VOL(N) + EPSL
            RSDX = MIN( ABS(RSDL(3,N))/ACCOX,RSDX )
            IF( RSDX.GT.RSD(IEQDOX) ) THEN
              RSD(IEQDOX) = RSDX
              NSD(IEQDOX) = N
            ENDIF
          ELSE
            IEQNOX = IEQO
            MNOEX = IM(IEQNOX,NMD)
          ENDIF
          RSDX = ABS(BLU(MNOEX))
          IF( RSDX.GT.RSD(IEQNOX) ) THEN
            RSD(IEQNOX) = RSDX
            NSD(IEQNOX) = N
          ENDIF
!
!---  Aqueous w/trapped-NAPL, trapped-gas system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
          RSDX = ABS(BLU(MWEX))/(ABS(PL(2,N))+PATM)
          ACCWX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQW,N))/ACCWX,RSDX )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          RSDX = ABS(BLU(MAEX))
          IF( RSDX.GT.RSD(IEQA) ) THEN
            RSD(IEQA) = RSDX
            NSD(IEQA) = N
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            IEQNOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
            MNOEX = IM(IEQNOX,NMD)
            RSDX = ABS(BLU(MDOEX))/MAX( XMLOX,1.D-9 )
            ACCOX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &        RHOG(2,N)*SG(2,N)*XGO(2,N))*DTI*VOL(N) + EPSL
            RSDX = MIN( ABS(RSDL(3,N))/ACCOX,RSDX )
            IF( RSDX.GT.RSD(IEQDOX) ) THEN
              RSD(IEQDOX) = RSDX
              NSD(IEQDOX) = N
            ENDIF
          ELSE
            IEQNOX = IEQO
            MNOEX = IM(IEQNOX,NMD)
          ENDIF
          RSDX = ABS(BLU(MNOEX))
          IF( RSDX.GT.RSD(IEQNOX) ) THEN
            RSD(IEQNOX) = RSDX
            NSD(IEQNOX) = N
          ENDIF
!
!---  Aqueous-gas w/ dissolved-oil system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
          RSDX = ABS(BLU(MWEX))/(ABS(PL(2,N))+PATM)
          ACCWX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQW,N))/ACCWX,RSDX )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          RSDX = ABS(BLU(MAEX))/(ABS(PG(2,N))+PATM)
          ACCAX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGA(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQA,N))/ACCAX,RSDX )
          IF( RSDX.GT.RSD(IEQA) ) THEN
            RSD(IEQA) = RSDX
            NSD(IEQA) = N
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            IEQNOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
            MNOEX = IM(IEQNOX,NMD)
            RSD(IEQNOX) = 0.D+0
            NSD(IEQNOX) = 0
          ELSE
            IEQDOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
          ENDIF
!
!---      Skip residual calculation for insignificant dissolved
!         oil concentrations  ---
!
          IF( XMLO(2,N).GT.1.D-5*XMLOX ) THEN
            RSDX = ABS(BLU(MDOEX))/MAX( XMLOX,1.D-9 )
            ACCOX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &        RHOG(2,N)*SG(2,N)*XGO(2,N))*DTI*VOL(N) + EPSL
            RSDX = MIN( ABS(RSDL(3,N))/ACCOX,RSDX )
            IF( RSDX.GT.RSD(IEQDOX) ) THEN
              RSD(IEQDOX) = RSDX
              NSD(IEQDOX) = N
            ENDIF
          ENDIF
!
!---  Aqueous-gas w/ trapped-NAPL system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
          RSDX = ABS(BLU(MWEX))/(ABS(PL(2,N))+PATM)
          ACCWX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLW(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGW(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQW,N))/ACCWX,RSDX )
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
          RSDX = ABS(BLU(MAEX))/(ABS(PG(2,N))+PATM)
          ACCAX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGA(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQA,N))/ACCAX,RSDX )
          IF( RSDX.GT.RSD(IEQA) ) THEN
            RSD(IEQA) = RSDX
            NSD(IEQA) = N
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            IEQNOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
            MNOEX = IM(IEQNOX,NMD)
            RSDX = ABS(BLU(MDOEX))/MAX( XMLOX,1.D-9 )
            ACCOX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &        RHOG(2,N)*SG(2,N)*XGO(2,N))*DTI*VOL(N) + EPSL
            RSDX = MIN( ABS(RSDL(3,N))/ACCOX,RSDX )
            IF( RSDX.GT.RSD(IEQDOX) ) THEN
              RSD(IEQDOX) = RSDX
              NSD(IEQDOX) = N
            ENDIF
          ELSE
            IEQNOX = IEQO
            MNOEX = IM(IEQNOX,NMD)
          ENDIF
          RSDX = ABS(BLU(MNOEX))
          ACCOX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGO(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(3,N))/ACCOX,RSDX )
          IF( RSDX.GT.RSD(IEQNOX) ) THEN
            RSD(IEQNOX) = RSDX
            NSD(IEQNOX) = N
          ENDIF
!
!---    Aqueous-NAPL w/ dissolved-air system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.7 ) THEN
          ACCWX = PORD(2,N)*RHOL(2,N)*SL(2,N)*DTI*VOL(N) + EPSL
          RSDX = ABS(RSDL(IEQW,N)/ACCWX)
          RSDX = MIN( RSDX, ABS(BLU(MWEX))/(ABS(PL(2,N))+PATM) )
          IF( 1.D+5-(PN(2,N)-PL(2,N)).LE.EPSL ) RSDX = 0.D+0
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            IEQNOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
            MNOEX = IM(IEQNOX,NMD)
            RSDX = ABS(BLU(MDOEX))/MAX( XMLOX,1.D-9 )
            ACCOX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &        RHOG(2,N)*SG(2,N)*XGO(2,N))*DTI*VOL(N) + EPSL
            RSDX = MIN( ABS(RSDL(3,N))/ACCOX,RSDX )
            IF( RSDX.GT.RSD(IEQDOX) ) THEN
              RSD(IEQDOX) = RSDX
              NSD(IEQDOX) = N
            ENDIF
          ELSE
            IEQNOX = IEQO
            MNOEX = IM(IEQNOX,NMD)
          ENDIF
          ACCOX = PORD(2,N)*RHON(2,N)*SN(2,N)*DTI*VOL(N) + EPSL
          RSDX = ABS(RSDL(IEQNOX,N)/ACCOX)
          RSDX = MIN( RSDX, ABS(BLU(MNOEX))/(ABS(PN(2,N))+PATM) )
          IF( RSDX.GT.RSD(IEQNOX) ) THEN
            RSD(IEQNOX) = RSDX
            NSD(IEQNOX) = N
          ENDIF
          PTX = MAX( PL(2,N),PG(2,N),PN(2,N),0.D+0 )+PATM
          XMLAX = PTX/HCAW
          RSDX = ABS(BLU(MAEX))/XMLAX
          ACCAX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGA(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQA,N))/ACCAX,RSDX )
          IF( RSDX.GT.RSD(IEQA) ) THEN
            RSD(IEQA) = RSDX
            NSD(IEQA) = N
          ENDIF
!
!---    Aqueous-NAPL w/ trapped-gas system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 ) THEN
          ACCWX = PORD(2,N)*RHOL(2,N)*SL(2,N)*DTI*VOL(N) + EPSL
          RSDX = ABS(RSDL(IEQW,N)/ACCWX)
          RSDX = MIN( RSDX, ABS(BLU(MWEX))/(ABS(PL(2,N))+PATM) )
          IF( 1.D+5-(PN(2,N)-PL(2,N)).LE.EPSL ) RSDX = 0.D+0
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            IEQNOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
            MNOEX = IM(IEQNOX,NMD)
            RSDX = ABS(BLU(MDOEX))/MAX( XMLOX,1.D-9 )
            ACCOX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &        RHOG(2,N)*SG(2,N)*XGO(2,N))*DTI*VOL(N) + EPSL
            RSDX = MIN( ABS(RSDL(3,N))/ACCOX,RSDX )
            IF( RSDX.GT.RSD(IEQDOX) ) THEN
              RSD(IEQDOX) = RSDX
              NSD(IEQDOX) = N
            ENDIF
          ELSE
            IEQNOX = IEQO
            MNOEX = IM(IEQNOX,NMD)
          ENDIF
          ACCOX = PORD(2,N)*RHON(2,N)*SN(2,N)*DTI*VOL(N) + EPSL
          RSDX = ABS(RSDL(IEQNOX,N)/ACCOX)
          RSDX = MIN( RSDX, ABS(BLU(MNOEX))/(ABS(PN(2,N))+PATM) )
          IF( RSDX.GT.RSD(IEQNOX) ) THEN
            RSD(IEQNOX) = RSDX
            NSD(IEQNOX) = N
          ENDIF
          RSDX = ABS(BLU(MAEX))
          ACCAX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLA(2,N) +
     &      RHOG(2,N)*SG(2,N)*XGA(2,N))*DTI*VOL(N) + EPSL
          RSDX = MIN( ABS(RSDL(IEQA,N))/ACCAX,RSDX )
          IF( RSDX.GT.RSD(IEQA) ) THEN
            RSD(IEQA) = RSDX
            NSD(IEQA) = N
          ENDIF
!
!---    Aqueous-NAPL-gas system prior to iteration  ---
!
        ELSEIF( NPHAZ(2,N).EQ.9 ) THEN
          ACCWX = PORD(2,N)*RHOL(2,N)*SL(2,N)*DTI*VOL(N) + EPSL
          RSDX = ABS(RSDL(IEQW,N)/ACCWX)
          RSDX = MIN( RSDX, ABS(BLU(MWEX))/(ABS(PL(2,N))+PATM) )
          IF( 1.D+5-(PN(2,N)-PL(2,N)).LE.EPSL ) RSDX = 0.D+0
          IF( RSDX.GT.RSD(IEQW) ) THEN
            RSD(IEQW) = RSDX
            NSD(IEQW) = N
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            IEQNOX = IEQO
            MDOEX = IM(IEQDOX,NMD)
            MNOEX = IM(IEQNOX,NMD)
            RSDX = ABS(BLU(MDOEX))/MAX( XMLOX,1.D-9 )
            ACCOX = PORD(2,N)*(RHOL(2,N)*SL(2,N)*XLO(2,N) +
     &        RHOG(2,N)*SG(2,N)*XGO(2,N))*DTI*VOL(N) + EPSL
            RSDX = MIN( ABS(RSDL(3,N))/ACCOX,RSDX )
            IF( RSDX.GT.RSD(IEQDOX) ) THEN
              RSD(IEQDOX) = RSDX
              NSD(IEQDOX) = N
            ENDIF
          ELSE
            IEQNOX = IEQO
            MNOEX = IM(IEQNOX,NMD)
          ENDIF
          ACCOX = PORD(2,N)*RHON(2,N)*SN(2,N)*DTI*VOL(N) + EPSL
          RSDX = ABS(RSDL(IEQNOX,N)/ACCOX)
          RSDX = MIN( RSDX, ABS(BLU(MNOEX))/(ABS(PN(2,N))+PATM) )
          IF( RSDX.GT.RSD(IEQNOX) ) THEN
            RSD(IEQNOX) = RSDX
            NSD(IEQNOX) = N
          ENDIF
          RSDX = ABS(BLU(MAEX))/(ABS(PG(2,N))+PATM)
          IF( RSDX.GT.RSD(IEQA) ) THEN
            RSD(IEQA) = RSDX
            NSD(IEQA) = N
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
        IF( NSD(IEQW).GT.0 ) THEN
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6)')
     &    'Water Equation Maximum Residual = ',RSD(IEQW),
     &    ' Node = ',NSD(IEQW),' Phase Condition = ',NPHAZ(2,NSD(IEQW))
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6)')
     &    'Water Equation Maximum Residual = ',RSD(IEQW),
     &    ' Node = ',NSD(IEQW),' Phase Condition = ',NPHAZ(2,NSD(IEQW))
        ENDIF
        IF( NSD(IEQA).GT.0 ) THEN
          WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6)')
     &    'Air Equation Maximum Residual = ',RSD(IEQA),
     &    ' Node = ',NSD(IEQA),' Phase Condition = ',NPHAZ(2,NSD(IEQA))
          WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6)')
     &    'Air Equation Maximum Residual = ',RSD(IEQA),
     &    ' Node = ',NSD(IEQA),' Phase Condition = ',NPHAZ(2,NSD(IEQA))
        ENDIF
!
!---    Kinetic volatilization option  ---
!
        IF( ISLC(46).GE.1 ) THEN
          IF( NSD(IEQDO).GT.0 ) THEN
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6)')
     &        'Dissolved Oil Equation Maximum Residual = ',RSD(IEQDO),
     &        ' Node = ',NSD(IEQDO),' Phase Condition = ',
     &        NPHAZ(2,NSD(IEQDO))
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6)')
     &        'Dissolved Oil Equation Maximum Residual = ',RSD(IEQDO),
     &        ' Node = ',NSD(IEQDO),' Phase Condition = ',
     &        NPHAZ(2,NSD(IEQDO))
          ENDIF
          IF( NSD(IEQO).GT.0 ) THEN
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6)')
     &        'NAPL Oil Equation Maximum Residual = ',RSD(IEQO),
     &        ' Node = ',NSD(IEQO),' Phase Condition = ',
     &        NPHAZ(2,NSD(IEQO))
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6)')
     &        'NAPL Oil Equation Maximum Residual = ',RSD(IEQO),
     &        ' Node = ',NSD(IEQO),' Phase Condition = ',
     &        NPHAZ(2,NSD(IEQO))
          ENDIF
        ELSE
          IF( NSD(IEQO).GT.0 ) THEN
            WRITE(ISC,'(4X,A,1PE11.4,A,I6,A,I6)')
     &        'Oil Equation Maximum Residual = ',RSD(IEQO),
     &        ' Node = ',NSD(IEQO),' Phase Condition = ',
     &        NPHAZ(2,NSD(IEQO))
            WRITE(IWR,'(4X,A,1PE11.4,A,I6,A,I6)')
     &        'Oil Equation Maximum Residual = ',RSD(IEQO),
     &        ' Node = ',NSD(IEQO),' Phase Condition = ',
     &        NPHAZ(2,NSD(IEQO))
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
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            XMLA(2,N) = XMLA(1,N)
            XMLO(2,N) = XMLO(1,N)
            SN(2,N) = SN(1,N)
            SG(2,N) = SG(1,N)
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
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDL5 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SBND5( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Modify the Jacobian matrix for the solute transport equation
!     to incorporate boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, October 16, 1995.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FDVP
      USE FDVN
      USE FDVG
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/SBND5'
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
!---  Diffusion coefficients at nodes adjacent to boundaries  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SMDLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
        SVLP = SL(2,N)*PORD(2,N)
        DLP = TORL(2,N)*SVLP*SMDLP
        PCOR = (PG(2,N)+PATM)/PATM
        SMDGP = SMDG(NSL)*(TCOR**1.75)/PCOR
        SVGP = SG(2,N)*PORD(2,N)
        DGP = TORG(2,N)*SVGP*SMDGP
        SMDNP = SMDN(NSL)*TCOR*(VISRN/VISN(2,N))
        SVNP = SN(2,N)*PORD(2,N)
        DNP = TORN(2,N)*SVNP*SMDNP
!
!---  Compute phase fractions  ---
!
        FCLP = 0.D+0
        IF( SL(2,N).GT.SMALL ) FCLP = YL(N,NSL)/(SVLP+SMALL)
        FCGP = 0.D+0
        IF( SG(2,N).GT.SMALL ) FCGP = YG(N,NSL)/(SVGP+SMALL)
        FCNP = 0.D+0
        IF( SN(2,N).GT.SMALL ) FCNP = YN(N,NSL)/(SVNP+SMALL)
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
        YLX = SVLP/(RHOS(IZN)*PCSL(1,IZN,NSL))
        YNX = SVNP/(RHOS(IZN)*PCSL(1,IZN,NSL)*PCLN(1,NSL))
        YGX = SVGP*PCGLX/(RHOS(IZN)*PCSL(1,IZN,NSL))
        YSX = 1.D+0-PORT(2,N)
        YLB(NB,NSL) = YLX/(YNX + YLX + YGX + YSX)
        YGB(NB,NSL) = YGX/(YNX + YLX + YGX + YSX)
        YNB(NB,NSL) = YNX/(YNX + YLX + YGX + YSX)
        FCLB = YLB(NB,NSL)/(SLB(2,NB)*PORDB(2,NB)+SMALL)
        FCGB = YGB(NB,NSL)/(SGB(2,NB)*PORDB(2,NB)+SMALL)
        FCNB = YNB(NB,NSL)/(SNB(2,NB)*PORDB(2,NB)+SMALL)
        IF( IBCT(NSL+LUK,NB).EQ.8 ) THEN
            BCX = BCX/( FCLB+SMALL )
        ELSEIF( IBCT(NSL+LUK,NB).EQ.9 ) THEN
            BCX = BCX/( FCGB+SMALL )
        ELSEIF( IBCT(NSL+LUK,NB).EQ.10 ) THEN
            BCX = BCX/( FCNB+SMALL )
        ENDIF
        CB(NB,NSL) = BCX
!
!-----  Phase fraction factors at boundary  ---
!
        SVLB = PORDB(2,NB)*SLB(2,NB)
        SVGB = PORDB(2,NB)*SGB(2,NB)
        SVNB = PORDB(2,NB)*SNB(2,NB)
!
!---  Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
          NPZ = NSZ(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULBX,VLBX,WLBX,N,MF )
            CALL SHDP( WLBX,ULBX,VLBX,DISPL(IZN),DISPT(IZN),DPLB )
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNBX,VNBX,WNBX,N,MF )
            CALL SHDP( WNBX,UNBX,VNBX,DISPL(IZN),DISPT(IZN),DPNB )
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGBX,VGBX,WGBX,N,MF )
            CALL SHDP( WGBX,UGBX,VGBX,DISPL(IZN),DISPT(IZN),DPGB )
          ELSE
            DPLB = 0.D+0
            DPNB = 0.D+0
            DPGB = 0.D+0
          ENDIF
          FLB = AFZ(NPZ)*WL(1,NPZ)
          FGB = AFZ(NPZ)*WG(1,NPZ)
          FNB = AFZ(NPZ)*WN(1,NPZ)
          CRLB = ABS( WL(1,NPZ) )*DT/(DZGF(N)*SVLB+SMALL)
          CRGB = ABS( WG(1,NPZ) )*DT/(DZGF(N)*SVGB+SMALL)
          CRNB = ABS( WN(1,NPZ) )*DT/(DZGF(N)*SVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.10 .OR.
     &      IBCT(NSL+LUK,NB).EQ.12 ) THEN
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
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SNB(2,NB)*PORDB(2,NB)*SMDNB
            INDX = 16
            DNZ = DIFMN(DNB,DNP,DZGF(N),DZGF(N),WN(1,NPZ),INDX)
            DNZ = AFZ(NPZ)*(DNZ+DPNB)/(5.D-1*DZGF(N))
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              IF( FLB.GE.ZERO ) THEN
                WCLZ = BCX*FCLB*FLB
              ELSEIF( FLB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCLT = YL(NBT,NSL)/(SL(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBT,NSL)*FCLT)
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                WCLZ = BCX*FLB*THETA*FCLB
     &            + C(N,NSL)*FLB*(1.D+0-THETA)*FCLP
              ELSEIF( FLB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCLZ = C(N,NSL)*FLB*FCLP
              ENDIF
              IF( FGB.GE.ZERO ) THEN
                WCGZ = BCX*FCGB*FGB
              ELSEIF( FGB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCGT = YG(NBT,NSL)/(SG(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBT,NSL)*FCGT)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                WCGZ = BCX*FGB*THETA*FCGB
     &            + C(N,NSL)*FGB*(1.D+0-THETA)*FCGP
              ELSEIF( FGB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCGZ = C(N,NSL)*FGB*FCGP
              ENDIF
              IF( FNB.GE.ZERO ) THEN
                WCNZ = BCX*FCNB*FNB
              ELSEIF( FNB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCNT = YN(NBT,NSL)/(SN(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBT,NSL)*FCNT)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRNB,ISLC(1) )
                WCNZ = BCX*FNB*THETA*FCNB
     &            + C(N,NSL)*FNB*(1.D+0-THETA)*FCNP
              ELSEIF( FNB.LT.ZERO .AND. K.LT.KFLD ) THEN
                WCNZ = C(N,NSL)*FNB*FCNP
              ENDIF
              AB = DLZ*FCLB + DGZ*FCGB + DNZ*FCNB
              AP = DGZ*FCGP + DLZ*FCLP + DNZ*FCNP
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCGZ+WCNZ)/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCGZ + WCNZ
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVLX = SL(2,NBT)*PORD(2,NBT)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(SVLX+SMALL)
                FCLT = YL(NBT,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBT)*PORD(2,NBT)
                CRGT = ABS( WG(1,NQZ) )*DT/DZGP(NQZ)/(SVGX+SMALL)
                FCGT = YG(NBT,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNT = AFZ(NQZ)*WN(1,NQZ)
              IF( FNT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVNX = SN(2,NBT)*PORD(2,NBT)
                CRNT = ABS( WN(1,NQZ) )*DT/DZGP(NQZ)/(SVNX+SMALL)
                FCNT = YN(NBT,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBT,NSL)*FCNT-C(N,NSL)*FCNP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRNT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCNZ = C(N,NSL)*FNT*(1.D+0-THETA*DZF)*FCNP
     &            + C(NBT,NSL)*FNT*THETA*DZF*FCNT
                WCNZF = CO(N,NSL)*FNT*FCNP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCNZ-WCNZF)/AFZ(NQZ)
                BLU(MP) = BLU(MP) - WCNZ + WCNZF
                BLU(IXP(NBT)) = BLU(IXP(NBT)) + WCNZ - WCNZF
              ENDIF
            ELSE
              ALB = MAX( FLB,ZERO ) +
     &          DLZ*MAX((ONE-(TENTH*ABS(FLB)/(DLZ+SMALL)))**5,ZERO)
              AGB = MAX( FGB,ZERO ) +
     &          DGZ*MAX((ONE-(TENTH*ABS(FGB)/(DGZ+SMALL)))**5,ZERO)
              ANB = MAX( FNB,ZERO ) +
     &          DNZ*MAX((ONE-(TENTH*ABS(FNB)/(DNZ+SMALL)))**5,ZERO)
              AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP + (ANB-FNB)*FCNP
              AB = ALB*FCLB + AGB*FCGB + ANB*FCNB
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
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                WCLZ = BCX*FLB*THETA*FCLB
     &            + C(N,NSL)*FLB*(1.D+0-THETA)*FCLP
              ELSEIF( FLB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCLZ = C(N,NSL)*FLB*FCLP
              ENDIF
              WCGZ = 0.D+0
              IF( FGB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCGT = YG(NBT,NSL)/(SG(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBT,NSL)*FCGT)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                WCGZ = BCX*FGB*THETA*FCGB
     &            + C(N,NSL)*FGB*(1.D+0-THETA)*FCGP
              ELSEIF( FGB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCGZ = C(N,NSL)*FGB*FCGP
              ENDIF
              WCNZ = 0.D+0
              IF( FNB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCNT = YN(NBT,NSL)/(SN(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBT,NSL)*FCNT)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRNB,ISLC(1) )
                WCNZ = BCX*FNB*THETA*FCNB
     &            + C(N,NSL)*FNB*(1.D+0-THETA)*FCNP
              ELSEIF( FNB.LT.ZERO .AND. K.LT.KFLD ) THEN
                WCNZ = C(N,NSL)*FNB*FCNP
              ENDIF
              AB = 0.D+0
              AP = 0.D+0
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCGZ+WCNZ)/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCGZ + WCNZ
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVLX = SL(2,NBT)*PORD(2,NBT)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(SVLX+SMALL)
                FCLT = YL(NBT,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBT)*PORD(2,NBT)
                CRGT = ABS( WG(1,NQZ) )*DT/DZGP(NQZ)/(SVGX+SMALL)
                FCGT = YG(NBT,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNT = AFZ(NQZ)*WN(1,NQZ)
              IF( FNT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVNX = SN(2,NBT)*PORD(2,NBT)
                CRNT = ABS( WN(1,NQZ) )*DT/DZGP(NQZ)/(SVNX+SMALL)
                FCNT = YN(NBT,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBT,NSL)*FCNT-C(N,NSL)*FCNP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRNT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCNZ = C(N,NSL)*FNT*(1.D+0-THETA*DZF)*FCNP
     &            + C(NBT,NSL)*FNT*THETA*DZF*FCNT
                WCNZF = CO(N,NSL)*FNT*FCNP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCNZ-WCNZF)/AFZ(NQZ)
                BLU(MP) = BLU(MP) - WCNZ + WCNZF
                BLU(IXP(NBT)) = BLU(IXP(NBT)) + WCNZ - WCNZF
              ENDIF
            ELSE
              ALB = MAX( FLB,ZERO )
              AGB = MAX( FGB,ZERO )
              ANB = MAX( FNB,ZERO )
              AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP + (ANB-FNB)*FCNP
              AB = ALB*FCLB + AGB*FCGB + ANB*FCNB
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
              IF( FLB.GE.ZERO ) WCLZ = BCX*FCLB*FLB
              WCGZ = 0.D+0
              IF( FGB.GE.ZERO ) WCGZ = BCX*FCGB*FGB
              WCNZ = 0.D+0
              IF( FNB.GE.ZERO ) WCNZ = BCX*FCNB*FNB
              AB = DLZ*FCLB + DGZ*FCGB + DNZ*FCNB
              AP = DGZ*FCGP + DLZ*FCLP + DNZ*FCNP
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCGZ+WCNZ)/AFZ(NPZ)
              BLU(MP) = BLU(MP) + WCLZ + WCGZ + WCNZ
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQZ = NPZ+IJFLD
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVLX = SL(2,NBT)*PORD(2,NBT)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(SVLX+SMALL)
                FCLT = YL(NBT,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBT)*PORD(2,NBT)
                CRGT = ABS( WG(1,NQZ) )*DT/DZGP(NQZ)/(SVGX+SMALL)
                FCGT = YG(NBT,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNT = AFZ(NQZ)*WN(1,NQZ)
              IF( FNT.GE.ZERO ) THEN
                NBT = N+IJFLD
                SVNX = SN(2,NBT)*PORD(2,NBT)
                CRNT = ABS( WN(1,NQZ) )*DT/DZGP(NQZ)/(SVNX+SMALL)
                FCNT = YN(NBT,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBT,NSL)*FCNT-C(N,NSL)*FCNP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRNT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCNZ = C(N,NSL)*FNT*(1.D+0-THETA*DZF)*FCNP
     &            + C(NBT,NSL)*FNT*THETA*DZF*FCNT
                WCNZF = CO(N,NSL)*FNT*FCNP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCNZ-WCNZF)/AFZ(NQZ)
                BLU(MP) = BLU(MP) - WCNZ + WCNZF
                BLU(IXP(NBT)) = BLU(IXP(NBT)) + WCNZ - WCNZF
              ENDIF
            ELSE
              ALB = MAX( FLB,ZERO )
              AGB = MAX( FGB,ZERO )
              ANB = MAX( FNB,ZERO )
              AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP + (ANB-FNB)*FCNP
              AB = ALB*FCLB + AGB*FCGB + ANB*FCNB
            ENDIF
            BLU(MP) = BLU(MP) + AB*BCX
          ENDIF
!
!---  South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          NPY = NSY(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULSX,VLSX,WLSX,N,MF )
            CALL SHDP( VLSX,WLSX,ULSX,DISPL(IZN),DISPT(IZN),DPLS )
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNSX,VNSX,WNSX,N,MF )
            CALL SHDP( VNSX,WNSX,UNSX,DISPL(IZN),DISPT(IZN),DPNS )
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGSX,VGSX,WGSX,N,MF )
            CALL SHDP( VGSX,WGSX,UGSX,DISPL(IZN),DISPT(IZN),DPGS )
          ELSE
            DPLS = 0.D+0
            DPNS = 0.D+0
            DPGS = 0.D+0
          ENDIF
          FLS = AFY(NPY)*VL(1,NPY)
          FGS = AFY(NPY)*VG(1,NPY)
          FNS = AFY(NPY)*VN(1,NPY)
          CRLS = ABS( VL(1,NPY) )*DT/(RP(I)*DYGF(N)*SVLB+SMALL)
          CRGS = ABS( VG(1,NPY) )*DT/(RP(I)*DYGF(N)*SVGB+SMALL)
          CRNS = ABS( VN(1,NPY) )*DT/(RP(I)*DYGF(N)*SVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.10 .OR.
     &      IBCT(NSL+LUK,NB).EQ.12 ) THEN
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
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SNB(2,NB)*PORDB(2,NB)*SMDNB
            INDX = 16
            DNY = DIFMN(DNB,DNP,DYGF(N),DYGF(N),VN(1,NPY),INDX)
            DNY = AFY(NPY)*(DNY+DPNS)/RP(I)/(5.D-1*DYGF(N))
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              IF( FLS.GE.ZERO ) THEN
                VCLY = BCX*FCLB*FLS
              ELSEIF( FLS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCLN = YL(NBN,NSL)/(SL(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBN,NSL)*FCLN)
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                VCLY = BCX*FLS*THETA*FCLB
     &            + C(N,NSL)*FLS*(1.D+0-THETA)*FCLP
              ELSEIF( FLS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCLY = C(N,NSL)*FLS*FCLP
              ENDIF
              IF( FGS.GE.ZERO ) THEN
                VCGY = BCX*FCGB*FGS
              ELSEIF( FGS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCGN = YG(NBN,NSL)/(SG(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBN,NSL)*FCGN)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                VCGY = BCX*FGS*THETA*FCGB
     &            + C(N,NSL)*FGS*(1.D+0-THETA)*FCGP
              ELSEIF( FGS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCGY = C(N,NSL)*FGS*FCGP
              ENDIF
              IF( FNS.GE.ZERO ) THEN
                VCNY = BCX*FCNB*FNS
              ELSEIF( FNS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCNN = YN(NBN,NSL)/(SN(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBN,NSL)*FCNN)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRNS,ISLC(1) )
                VCNY = BCX*FNS*THETA*FCNB
     &            + C(N,NSL)*FNS*(1.D+0-THETA)*FCNP
              ELSEIF( FNS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCNY = C(N,NSL)*FNS*FCNP
              ENDIF
              AS = DLY*FCLB + DGY*FCGB + DNY*FCNB
              AP = DGY*FCGP + DLY*FCLP + DNY*FCNP
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCGY+VCNY)/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCGY + VCNY
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQY = NPY+IFLD
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                SVLX = SL(2,NBN)*PORD(2,NBN)
                CRLN = ABS( VL(1,NQY) )*DT/DYGP(NQY)/(SVLX*RP(I)+SMALL)
                FCLN = YL(NBN,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBN)*PORD(2,NBN)
                CRGN = ABS( VG(1,NQY) )*DT/DYGP(NQY)/(SVGX*RP(I)+SMALL)
                FCGN = YG(NBN,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNN = AFY(NQY)*VN(1,NQY)
              IF( FNN.GE.ZERO ) THEN
                NBN = N+IFLD
                SVNX = SN(2,NBN)*PORD(2,NBN)
                CRNN = ABS( VN(1,NQY) )*DT/DYGP(NQY)/(SVNX*RP(I)+SMALL)
                FCNN = YN(NBN,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBN,NSL)*FCNN-C(N,NSL)*FCNP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRNN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCNY = C(N,NSL)*FNN*(1.D+0-THETA*DYF)*FCNP
     &            + C(NBN,NSL)*FNN*THETA*DYF*FCNN
                VCNYF = CO(N,NSL)*FNN*FCNP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCNY-VCNYF)/AFY(NQY)
                BLU(MP) = BLU(MP) - VCNY + VCNYF
                BLU(IXP(NBN)) = BLU(IXP(NBN)) + VCNY - VCNYF
              ENDIF
            ELSE
              ALS = MAX( FLS,ZERO ) +
     &          DLY*MAX((ONE-(TENTH*ABS(FLS)/(DLY+SMALL)))**5,ZERO)
              AGS = MAX( FGS,ZERO ) +
     &          DGY*MAX((ONE-(TENTH*ABS(FGS)/(DGY+SMALL)))**5,ZERO)
              ANS = MAX( FNS,ZERO ) +
     &          DNY*MAX((ONE-(TENTH*ABS(FNS)/(DNY+SMALL)))**5,ZERO)
              AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP + (ANS-FNS)*FCNP
              AS = ALS*FCLB + AGS*FCGB + ANS*FCNB
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
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                VCLY = BCX*FLS*THETA*FCLB
     &            + C(N,NSL)*FLS*(1.D+0-THETA)*FCLP
              ELSEIF( FLS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCLY = C(N,NSL)*FLS*FCLP
              ENDIF
              VCGY = 0.D+0
              IF( FGS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCGN = YG(NBN,NSL)/(SG(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBN,NSL)*FCGN)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                VCGY = BCX*FGS*THETA*FCGB
     &            + C(N,NSL)*FGS*(1.D+0-THETA)*FCGP
              ELSEIF( FGS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCGY = C(N,NSL)*FGS*FCGP
              ENDIF
              VCNY = 0.D+0
              IF( FNS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCNN = YN(NBN,NSL)/(SN(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBN,NSL)*FCNN)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRNS,ISLC(1) )
                VCNY = BCX*FNS*THETA*FCNB
     &            + C(N,NSL)*FNS*(1.D+0-THETA)*FCNP
              ELSEIF( FNS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCNY = C(N,NSL)*FNS*FCNP
              ENDIF
              AS = 0.D+0
              AP = 0.D+0
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCGY+VCNY)/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCGY + VCNY
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQY = NPY+IFLD
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                SVLX = SL(2,NBN)*PORD(2,NBN)
                CRLN = ABS( VL(1,NQY) )*DT/DYGP(NQY)/(SVLX*RP(I)+SMALL)
                FCLN = YL(NBN,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBN)*PORD(2,NBN)
                CRGN = ABS( VG(1,NQY) )*DT/DYGP(NQY)/(SVGX*RP(I)+SMALL)
                FCGN = YG(NBN,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNN = AFY(NQY)*VN(1,NQY)
              IF( FNN.GE.ZERO ) THEN
                NBN = N+IFLD
                SVNX = SN(2,NBN)*PORD(2,NBN)
                CRNN = ABS( VN(1,NQY) )*DT/DYGP(NQY)/(SVNX*RP(I)+SMALL)
                FCNN = YN(NBN,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBN,NSL)*FCNN-C(N,NSL)*FCNP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRNN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCNY = C(N,NSL)*FNN*(1.D+0-THETA*DYF)*FCNP
     &            + C(NBN,NSL)*FNN*THETA*DYF*FCNN
                VCNYF = CO(N,NSL)*FNN*FCNP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCNY-VCNYF)/AFY(NQY)
                BLU(MP) = BLU(MP) - VCNY + VCNYF
                BLU(IXP(NBN)) = BLU(IXP(NBN)) + VCNY - VCNYF
              ENDIF
            ELSE
              ALS = MAX( FLS,ZERO )
              AGS = MAX( FGS,ZERO )
              ANS = MAX( FNS,ZERO )
              AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP + (ANS-FNS)*FCNP
              AS = ALS*FCLB + AGS*FCGB + ANS*FCNB
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
              IF( FLS.GE.ZERO ) VCLY = BCX*FCLB*FLS
              VCGY = 0.D+0
              IF( FGS.GE.ZERO ) VCGY = BCX*FCGB*FGS
              VCNY = 0.D+0
              IF( FNS.GE.ZERO ) VCNY = BCX*FCNB*FNS
              AS = 0.D+0
              AP = 0.D+0
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCGY+VCNY)/AFY(NPY)
              BLU(MP) = BLU(MP) + VCLY + VCGY + VCNY
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQY = NPY+IFLD
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                SVLX = SL(2,NBN)*PORD(2,NBN)
                CRLN = ABS( VL(1,NQY) )*DT/DYGP(NQY)/(SVLX*RP(I)+SMALL)
                FCLN = YL(NBN,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBN)*PORD(2,NBN)
                CRGN = ABS( VG(1,NQY) )*DT/DYGP(NQY)/(SVGX*RP(I)+SMALL)
                FCGN = YG(NBN,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNN = AFY(NQY)*VN(1,NQY)
              IF( FNN.GE.ZERO ) THEN
                NBN = N+IFLD
                SVNX = SN(2,NBN)*PORD(2,NBN)
                CRNN = ABS( VN(1,NQY) )*DT/DYGP(NQY)/(SVNX*RP(I)+SMALL)
                FCNN = YN(NBN,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBN,NSL)*FCNN-C(N,NSL)*FCNP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRNN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCNY = C(N,NSL)*FNN*(1.D+0-THETA*DYF)*FCNP
     &            + C(NBN,NSL)*FNN*THETA*DYF*FCNN
                VCNYF = CO(N,NSL)*FNN*FCNP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCNY-VCNYF)/AFY(NQY)
                BLU(MP) = BLU(MP) - VCNY + VCNYF
                BLU(IXP(NBN)) = BLU(IXP(NBN)) + VCNY - VCNYF
              ENDIF
            ELSE
              ALS = MAX( FLS,ZERO )
              AGS = MAX( FGS,ZERO )
              ANS = MAX( FNS,ZERO )
              AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP + (ANS-FNS)*FCNP
              AS = ALS*FCLB + AGS*FCGB + ANS*FCNB
            ENDIF
            BLU(MP) = BLU(MP) + AS*BCX
          ENDIF
!
!---  West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          NPX = NSX(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULWX,VLWX,WLWX,N,MF )
            CALL SHDP( ULWX,VLWX,WLWX,DISPL(IZN),DISPT(IZN),DPLW )
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNWX,VNWX,WNWX,N,MF )
            CALL SHDP( UNWX,VNWX,WNWX,DISPL(IZN),DISPT(IZN),DPNW )
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGWX,VGWX,WGWX,N,MF )
            CALL SHDP( UGWX,VGWX,WGWX,DISPL(IZN),DISPT(IZN),DPGW )
          ELSE
            DPLW = 0.D+0
            DPNW = 0.D+0
            DPGW = 0.D+0
          ENDIF
          FLW = AFX(NPX)*UL(1,NPX)
          FGW = AFX(NPX)*UG(1,NPX)
          FNW = AFX(NPX)*UN(1,NPX)
          CRLW = ABS( UL(1,NPX) )*DT/(DXGF(N)*SVLB+SMALL)
          CRGW = ABS( UG(1,NPX) )*DT/(DXGF(N)*SVGB+SMALL)
          CRNW = ABS( UN(1,NPX) )*DT/(DXGF(N)*SVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.10 .OR.
     &      IBCT(NSL+LUK,NB).EQ.12 ) THEN
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
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SNB(2,NB)*PORDB(2,NB)*SMDNB
            INDX = 16
            DNX = DIFMN(DNB,DNP,DXGF(N),DXGF(N),UN(1,NPX),INDX)
            DNX = AFX(NPX)*(DNX+DPNW)/(5.D-1*DXGF(N))
!
!---  TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 )  THEN
              IF( FLW.GE.ZERO ) THEN
                UCLX = BCX*FCLB*FLW
              ELSEIF( FLW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCLE = YL(NBE,NSL)/(SL(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBE,NSL)*FCLE)
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA)*FCLP
     &            + BCX*FLW*THETA*FCLB
              ENDIF
              IF( FGW.GE.ZERO ) THEN
                UCGX = BCX*FCGB*FGW
              ELSEIF( FGW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCGE = YG(NBE,NSL)/(SG(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBE,NSL)*FCGE)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA)*FCGP
     &            + BCX*FGW*THETA*FCGB
              ELSEIF( FGW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCGX = C(N,NSL)*FGW*FCGP
              ENDIF
              IF( FNW.GE.ZERO ) THEN
                UCNX = BCX*FCNB*FNW
              ELSEIF( FNW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCNE = YN(NBE,NSL)/(SN(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBE,NSL)*FCNE)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRNW,ISLC(1) )
                UCNX = C(N,NSL)*FNW*(1.D+0-THETA)*FCNP
     &            + BCX*FNW*THETA*FCNB
              ELSEIF( FNW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCNX = C(N,NSL)*FNW*FCNP
              ENDIF
              AW = DLX*FCLB + DGX*FCGB + DNX*FCNB
              AP = DLX*FCLP + DGX*FCGP + DNX*FCNP
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCGX+UCNX)/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCGX + UCNX
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQX = NPX+1
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                SVLX = SL(2,NBE)*PORD(2,NBE)
                CRLE = ABS( UL(1,NQX) )*DT/DXGP(NQX)/(SVLX+SMALL)
                FCLE = YL(NBE,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBE)*PORD(2,NBE)
                CRGE = ABS( UG(1,NQX) )*DT/DXGP(NQX)/(SVGX+SMALL)
                FCGE = YG(NBE,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNE = AFX(NQX)*UN(1,NQX)
              IF( FNE.GE.ZERO ) THEN
                NBE = N+1
                SVNX = SN(2,NBE)*PORD(2,NBE)
                CRNE = ABS( UN(1,NQX) )*DT/DXGP(NQX)/(SVNX+SMALL)
                FCNE = YN(NBE,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBE,NSL)*FCNE-C(N,NSL)*FCNP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRNE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCNX = C(N,NSL)*FNE*(1.D+0-THETA*DXF)*FCNP
     &            + C(NBE,NSL)*FNE*THETA*DXF*FCNE
                UCNXF = CO(N,NSL)*FNE*FCNP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCNX-UCNXF)/AFX(NQX)
                BLU(MP) = BLU(MP) - UCNX + UCNXF
                BLU(IXP(NBE)) = BLU(IXP(NBE)) + UCNX - UCNXF
              ENDIF
            ELSE
              ALW = MAX(FLW,ZERO)
     &          + DLX*MAX((ONE-(TENTH*ABS(FLW)/(DLX+SMALL)))**5,ZERO)
              AGW = MAX(FGW,ZERO)
     &          + DGX*MAX((ONE-(TENTH*ABS(FGW)/(DGX+SMALL)))**5,ZERO)
              ANW = MAX(FNW,ZERO)
     &          + DNX*MAX((ONE-(TENTH*ABS(FNW)/(DNX+SMALL)))**5,ZERO)
              AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP + (ANW-FNW)*FCNP
              AW = ALW*FCLB + AGW*FCGB + ANW*FCNB
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
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA)*FCLP
     &            + BCX*FLW*THETA*FCLB
              ENDIF
              UCGX = 0.D+0
              IF( FGW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCGE = YG(NBE,NSL)/(SG(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBE,NSL)*FCGE)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA)*FCGP
     &            + BCX*FGW*THETA*FCGB
              ELSEIF( FGW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCGX = C(N,NSL)*FGW*FCGP
              ENDIF
              UCNX = 0.D+0
              IF( FNW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCNE = YN(NBE,NSL)/(SN(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBE,NSL)*FCNE)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRNW,ISLC(1) )
                UCNX = C(N,NSL)*FNW*(1.D+0-THETA)*FCNP
     &            + BCX*FNW*THETA*FCNB
              ELSEIF( FNW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCNX = C(N,NSL)*FNW*FCNP
              ENDIF
              AW = 0.D+0
              AP = 0.D+0
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCGX+UCNX)/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCGX + UCNX
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQX = NPX+1
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                SVLX = SL(2,NBE)*PORD(2,NBE)
                CRLE = ABS( UL(1,NQX) )*DT/DXGP(NQX)/(SVLX+SMALL)
                FCLE = YL(NBE,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBE)*PORD(2,NBE)
                CRGE = ABS( UG(1,NQX) )*DT/DXGP(NQX)/(SVGX+SMALL)
                FCGE = YG(NBE,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNE = AFX(NQX)*UN(1,NQX)
              IF( FNE.GE.ZERO ) THEN
                NBE = N+1
                SVNX = SN(2,NBE)*PORD(2,NBE)
                CRNE = ABS( UN(1,NQX) )*DT/DXGP(NQX)/(SVNX+SMALL)
                FCNE = YN(NBE,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBE,NSL)*FCNE-C(N,NSL)*FCNP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRNE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCNX = C(N,NSL)*FNE*(1.D+0-THETA*DXF)*FCNP
     &            + C(NBE,NSL)*FNE*THETA*DXF*FCNE
                UCNXF = CO(N,NSL)*FNE*FCNP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCNX-UCNXF)/AFX(NQX)
                BLU(MP) = BLU(MP) - UCNX + UCNXF
                BLU(IXP(NBE)) = BLU(IXP(NBE)) + UCNX - UCNXF
              ENDIF
            ELSE
              ALW = MAX(FLW,ZERO)
              AGW = MAX(FGW,ZERO)
              ANW = MAX(FNW,ZERO)
              AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP + (ANW-FNW)*FCNP
              AW = ALW*FCLB + AGW*FCGB + ANW*FCNB
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
              IF( FLW.GE.ZERO ) UCLX = BCX*FCLB*FLW
              UCGX = 0.D+0
              IF( FGW.GE.ZERO ) UCGX = BCX*FCGB*FGW
              UCNX = 0.D+0
              IF( FNW.GE.ZERO ) UCNX = BCX*FCNB*FNW
              AW = 0.D+0
              AP = 0.D+0
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCGX+UCNX)/AFX(NPX)
              BLU(MP) = BLU(MP) + UCLX + UCGX + UCNX
!
!---  TVD Transport for interior surface adjacent to boundary  ---
!
              NQX = NPX+1
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                SVLX = SL(2,NBE)*PORD(2,NBE)
                CRLE = ABS( UL(1,NQX) )*DT/DXGP(NQX)/(SVLX+SMALL)
                FCLE = YL(NBE,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBE)*PORD(2,NBE)
                CRGE = ABS( UG(1,NQX) )*DT/DXGP(NQX)/(SVGX+SMALL)
                FCGE = YG(NBE,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNE = AFX(NQX)*UN(1,NQX)
              IF( FNE.GE.ZERO ) THEN
                NBE = N+1
                SVNX = SN(2,NBE)*PORD(2,NBE)
                CRNE = ABS( UN(1,NQX) )*DT/DXGP(NQX)/(SVNX+SMALL)
                FCNE = YN(NBE,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBE,NSL)*FCNE-C(N,NSL)*FCNP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRNE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCNX = C(N,NSL)*FNE*(1.D+0-THETA*DXF)*FCNP
     &            + C(NBE,NSL)*FNE*THETA*DXF*FCNE
                UCNXF = CO(N,NSL)*FNE*FCNP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCNX-UCNXF)/AFX(NQX)
                BLU(MP) = BLU(MP) - UCNX + UCNXF
                BLU(IXP(NBE)) = BLU(IXP(NBE)) + UCNX - UCNXF
              ENDIF
            ELSE
              ALW = MAX(FLW,ZERO)
              AGW = MAX(FGW,ZERO)
              ANW = MAX(FNW,ZERO)
              AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP + (ANW-FNW)*FCNP
              AW = ALW*FCLB + AGW*FCGB + ANW*FCNB
            ENDIF
            BLU(MP) = BLU(MP) + AW*BCX
          ENDIF
!
!---  East boundary
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          NQX = NSX(N) + 1
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULEX,VLEX,WLEX,N,MF )
            CALL SHDP( ULEX,VLEX,WLEX,DISPL(IZN),DISPT(IZN),DPLE )
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNEX,VNEX,WNEX,N,MF )
            CALL SHDP( UNEX,VNEX,WNEX,DISPL(IZN),DISPT(IZN),DPNE )
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGEX,VGEX,WGEX,N,MF )
            CALL SHDP( UGEX,VGEX,WGEX,DISPL(IZN),DISPT(IZN),DPGE )
          ELSE
            DPLE = 0.D+0
            DPNE = 0.D+0
            DPGE = 0.D+0
          ENDIF
          FLE = AFX(NQX)*UL(1,NQX)
          FGE = AFX(NQX)*UG(1,NQX)
          FNE = AFX(NQX)*UN(1,NQX)
          CRLE = ABS( UL(1,NQX) )*DT/(DXGF(N)*SVLB+SMALL)
          CRGE = ABS( UG(1,NQX) )*DT/(DXGF(N)*SVGB+SMALL)
          CRNE = ABS( UN(1,NQX) )*DT/(DXGF(N)*SVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.10 .OR.
     &      IBCT(NSL+LUK,NB).EQ.12 ) THEN
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
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SNB(2,NB)*PORDB(2,NB)*SMDNB
            INDX = 16
            DNX = DIFMN(DNP,DNB,DXGF(N),DXGF(N),UN(1,NQX),INDX)
            DNX = AFX(NQX)*(DNX+DPNE)/(5.D-1*DXGF(N))
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              IF( FLE.LT.ZERO ) THEN
                UCLX = BCX*FCLB*FLE
              ELSEIF( FLE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCLW = YL(NBW,NSL)/(SL(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBW,NSL)*FCLW)
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                UCLX =  C(N,NSL)*FLE*(1.D+0-THETA)*FCLP
     &          +  BCX*FLE*THETA*FCLB
              ELSEIF( FLE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCLX =  C(N,NSL)*FLE*FCLP
              ENDIF
              IF( FGE.LT.ZERO ) THEN
                UCGX = BCX*FCGB*FGE
              ELSEIF( FGE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCGW = YG(NBW,NSL)/(SG(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBW,NSL)*FCGW)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                UCGX =  C(N,NSL)*FGE*(1.D+0-THETA)*FCGP
     &          +  BCX*FGE*THETA*FCGB
              ELSEIF( FGE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCGX =  C(N,NSL)*FGE*FCGP
              ENDIF
              IF( FNE.LT.ZERO ) THEN
                UCNX = BCX*FCNB*FNE
              ELSEIF( FNE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCNW = YN(NBW,NSL)/(SN(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBW,NSL)*FCNW)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRNE,ISLC(1) )
                UCNX =  C(N,NSL)*FNE*(1.D+0-THETA)*FCNP
     &          +  BCX*FNE*THETA*FCNB
              ELSEIF( FNE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCNX =  C(N,NSL)*FNE*FCNP
              ENDIF
              AE = DLX*FCLB + DGX*FCGB + DNX*FCNB
              AP = DLX*FCLP + DGX*FCGP + DNX*FCNP
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCGX+UCNX)/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCGX - UCNX
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPX = NSX(N)
              FLW = AFX(NPX)*UL(1,NPX)
              IF( FLW.LT.ZERO ) THEN
                NBW = N-1
                SVLX = SL(2,NBW)*PORD(2,NBW)
                CRLW = ABS( UL(1,NPX) )*DT/DXGP(NPX)/(SVLX+SMALL)
                FCLW = YL(NBW,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(SVGX+SMALL)
                FCGW = YG(NBW,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              IF( FNW.LT.ZERO ) THEN
                NBW = N-1
                SVNX = SN(2,NBW)*PORD(2,NBW)
                CRNW = ABS( UN(1,NPX) )*DT/DXGP(NPX)/(SVNX+SMALL)
                FCNW = YN(NBW,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBW,NSL)*FCNW-C(N,NSL)*FCNP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRNW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCNX = C(N,NSL)*FNW*(1.D+0-THETA*DXF)*FCNP
     &            + C(NBW,NSL)*FNW*THETA*DXF*FCNW
                UCNXF = CO(N,NSL)*FNW*FCNP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCNX-UCNXF)/AFX(NPX)
                BLU(MP) = BLU(MP) + UCNX - UCNXF
                BLU(IXP(NBW)) = BLU(IXP(NBW)) - UCNX + UCNXF
              ENDIF
            ELSE
              ALE = MAX( -FLE,ZERO ) +
     &          DLX*MAX((ONE-(TENTH*ABS(FLE)/(DLX+SMALL)))**5,ZERO)
              AGE = MAX( -FGE,ZERO ) +
     &          DGX*MAX((ONE-(TENTH*ABS(FGE)/(DGX+SMALL)))**5,ZERO)
              ANE = MAX( -FNE,ZERO ) +
     &          DNX*MAX((ONE-(TENTH*ABS(FNE)/(DNX+SMALL)))**5,ZERO)
              AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP + (ANE+FNE)*FCNP
              AE = ALE*FCLB + AGE*FCGB + ANE*FCNB
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
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                UCLX =  C(N,NSL)*FLE*(1.D+0-THETA)*FCLP
     &          +  BCX*FLE*THETA*FCLB
              ELSEIF( FLE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCLX =  C(N,NSL)*FLE*FCLP
              ENDIF
              UCGX = 0.D+0
              IF( FGE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCGW = YG(NBW,NSL)/(SG(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBW,NSL)*FCGW)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                UCGX =  C(N,NSL)*FGE*(1.D+0-THETA)*FCGP
     &          +  BCX*FGE*THETA*FCGB
              ELSEIF( FGE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCGX =  C(N,NSL)*FGE*FCGP
              ENDIF
              UCNX = 0.D+0
              IF( FNE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCNW = YN(NBW,NSL)/(SN(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBW,NSL)*FCNW)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRNE,ISLC(1) )
                UCNX =  C(N,NSL)*FNE*(1.D+0-THETA)*FCNP
     &          +  BCX*FNE*THETA*FCNB
              ELSEIF( FNE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCNX =  C(N,NSL)*FNE*FCNP
              ENDIF
              AE = 0.D+0
              AP = 0.D+0
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCGX+UCNX)/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCGX - UCNX
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPX = NSX(N)
              FLW = AFX(NPX)*UL(1,NPX)
              IF( FLW.LT.ZERO ) THEN
                NBW = N-1
                SVLX = SL(2,NBW)*PORD(2,NBW)
                CRLW = ABS( UL(1,NPX) )*DT/DXGP(NPX)/(SVLX+SMALL)
                FCLW = YL(NBW,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(SVGX+SMALL)
                FCGW = YG(NBW,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              IF( FNW.LT.ZERO ) THEN
                NBW = N-1
                SVNX = SN(2,NBW)*PORD(2,NBW)
                CRNW = ABS( UN(1,NPX) )*DT/DXGP(NPX)/(SVNX+SMALL)
                FCNW = YN(NBW,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBW,NSL)*FCNW-C(N,NSL)*FCNP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRNW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCNX = C(N,NSL)*FNW*(1.D+0-THETA*DXF)*FCNP
     &            + C(NBW,NSL)*FNW*THETA*DXF*FCNW
                UCNXF = CO(N,NSL)*FNW*FCNP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCNX-UCNXF)/AFX(NPX)
                BLU(MP) = BLU(MP) + UCNX - UCNXF
                BLU(IXP(NBW)) = BLU(IXP(NBW)) - UCNX + UCNXF
              ENDIF
            ELSE
              ALE = MAX( -FLE,ZERO )
              AGE = MAX( -FGE,ZERO )
              ANE = MAX( -FNE,ZERO )
              AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP + (ANE+FNE)*FCNP
              AE = ALE*FCLB + AGE*FCGB + ANE*FCNB
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
              IF( FLE.LT.ZERO ) UCLX = BCX*FCLB*FLE
              UCGX = 0.D+0
              IF( FGE.LT.ZERO ) UCGX = BCX*FCGB*FGE
              UCNX = 0.D+0
              IF( FNE.LT.ZERO )  UCNX = BCX*FCNB*FNE
              AE = 0.D+0
              AP = 0.D+0
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCGX+UCNX)/AFX(NQX)
              BLU(MP) = BLU(MP) - UCLX - UCGX - UCNX
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPX = NSX(N)
              FLW = AFX(NPX)*UL(1,NPX)
              IF( FLW.LT.ZERO ) THEN
                NBW = N-1
                SVLX = SL(2,NBW)*PORD(2,NBW)
                CRLW = ABS( UL(1,NPX) )*DT/DXGP(NPX)/(SVLX+SMALL)
                FCLW = YL(NBW,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(SVGX+SMALL)
                FCGW = YG(NBW,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              IF( FNW.LT.ZERO ) THEN
                NBW = N-1
                SVNX = SN(2,NBW)*PORD(2,NBW)
                CRNW = ABS( UN(1,NPX) )*DT/DXGP(NPX)/(SVNX+SMALL)
                FCNW = YN(NBW,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBW,NSL)*FCNW-C(N,NSL)*FCNP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRNW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCNX = C(N,NSL)*FNW*(1.D+0-THETA*DXF)*FCNP
     &            + C(NBW,NSL)*FNW*THETA*DXF*FCNW
                UCNXF = CO(N,NSL)*FNW*FCNP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCNX-UCNXF)/AFX(NPX)
                BLU(MP) = BLU(MP) + UCNX - UCNXF
                BLU(IXP(NBW)) = BLU(IXP(NBW)) - UCNX + UCNXF
              ENDIF
            ELSE
              ALE = MAX( -FLE,ZERO )
              AGE = MAX( -FGE,ZERO )
              ANE = MAX( -FNE,ZERO )
              AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP + (ANE+FNE)*FCNP
              AE = ALE*FCLB + AGE*FCGB + ANE*FCNB
            ENDIF
            BLU(MP) = BLU(MP) + AE*BCX
          ENDIF
!
!---  North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          NQY = NSY(N) + IFLD
!
!---  Hydraulic dispersion  ---
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULNX,VLNX,WLNX,N,MF )
            CALL SHDP( VLNX,WLNX,ULNX,DISPL(IZN),DISPT(IZN),DPLN )
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNNX,VNNX,WNNX,N,MF )
            CALL SHDP( VNNX,WNNX,UNNX,DISPL(IZN),DISPT(IZN),DPNN )
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGNX,VGNX,WGNX,N,MF )
            CALL SHDP( VGNX,WGNX,UGNX,DISPL(IZN),DISPT(IZN),DPGN )
          ELSE
            DPLN = 0.D+0
            DPNN = 0.D+0
            DPGN = 0.D+0
          ENDIF
          FLN = AFY(NQY)*VL(1,NQY)
          FGN = AFY(NQY)*VG(1,NQY)
          FNN = AFY(NQY)*VN(1,NQY)
          CRLN = ABS( VL(1,NQY) )*DT/(RP(I)*DYGF(N)*SVLB+SMALL)
          CRGN = ABS( VG(1,NQY) )*DT/(RP(I)*DYGF(N)*SVGB+SMALL)
          CRNN = ABS( VN(1,NQY) )*DT/(RP(I)*DYGF(N)*SVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.10 .OR.
     &      IBCT(NSL+LUK,NB).EQ.12 ) THEN
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
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SNB(2,NB)*PORDB(2,NB)*SMDNB
            INDX = 16
            DNY = DIFMN(DNP,DNB,DYGF(N),DYGF(N),VN(1,NQY),INDX)
            DNY = AFY(NQY)*(DNY+DPNN)/RP(I)/(5.D-1*DYGF(N))
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              IF( FLN.LT.ZERO ) THEN
                VCLY = BCX*FCLB*FLN
              ELSEIF( FLN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCLS = YL(NBS,NSL)/(SL(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBS,NSL)*FCLS)
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                VCLY =  BCX*FLN*THETA*FCLB
     &            + C(N,NSL)*FLN*(1.D+0-THETA)*FCLP
              ELSEIF( FLN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCLY = C(N,NSL)*FLN*FCLP
              ENDIF
              IF( FGN.LT.ZERO ) THEN
                VCGY = BCX*FCGB*FGN
              ELSEIF( FGN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCGS = YG(NBS,NSL)/(SG(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBS,NSL)*FCGS)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                VCGY =  BCX*FGN*THETA*FCGB
     &            + C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ELSEIF( FGN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCGY =  C(N,NSL)*FGN*FCGP
              ENDIF
              IF( FNN.LT.ZERO ) THEN
                VCNY = BCX*FCNB*FNN
              ELSEIF( FNN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCNS = YN(NBS,NSL)/(SN(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBS,NSL)*FCNS)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRNN,ISLC(1) )
                VCNY =  BCX*FNN*THETA*FCNB
     &            + C(N,NSL)*FNN*(1.D+0-THETA)*FCNP
              ELSEIF( FNN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCNY =  C(N,NSL)*FNN*FCNP
              ENDIF
              AN = DLY*FCLB + DGY*FCGB + DNY*FCNB
              AP = DLY*FCLP + DGY*FCGP + DNY*FCNP
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCGY+VCNY)/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCGY - VCNY
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPY = NSY(N)
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                SVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/DYGP(NPY)/(SVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/DYGP(NPY)/(SVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNS = AFY(NPY)*VN(1,NPY)
              IF( FNS.LT.ZERO ) THEN
                NBS = N-IFLD
                SVNX = SN(2,NBS)*PORD(2,NBS)
                CRNS = ABS( VN(1,NPY) )*DT/DYGP(NPY)/(SVNX*RP(I)+SMALL)
                FCNS = YN(NBS,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBS,NSL)*FCNS-C(N,NSL)*FCNP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRNS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCNY = C(N,NSL)*FNS*(1.D+0-THETA*DYF)*FCNP
     &            + C(NBS,NSL)*FNS*THETA*DYF*FCNS
                VCNYF = CO(N,NSL)*FNS*FCNP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCNY-VCNYF)/AFY(NPY)
                BLU(MP) = BLU(MP) + VCNY - VCNYF
                BLU(IXP(NBS)) = BLU(IXP(NBS)) - VCNY + VCNYF
              ENDIF
            ELSE
              ALN = MAX( -FLN,ZERO ) +
     &          DLY*MAX((ONE-(TENTH*ABS(FLN)/(DLY+SMALL)))**5,ZERO)
              AGN = MAX( -FGN,ZERO ) +
     &          DGY*MAX((ONE-(TENTH*ABS(FGN)/(DGY+SMALL)))**5,ZERO)
              ANN = MAX( -FNN,ZERO ) +
     &          DNY*MAX((ONE-(TENTH*ABS(FNN)/(DNY+SMALL)))**5,ZERO)
              AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP + (ANN+FNN)*FCNP
              AN = ALN*FCLB + AGN*FCGB + ANN*FCNB
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
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                VCLY =  BCX*FLN*THETA*FCLB
     &            + C(N,NSL)*FLN*(1.D+0-THETA)*FCLP
              ELSEIF( FLN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCLY = C(N,NSL)*FLN*FCLP
              ENDIF
              VCGY = 0.D+0
              IF( FGN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCGS = YG(NBS,NSL)/(SG(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBS,NSL)*FCGS)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                VCGY =  BCX*FGN*THETA*FCGB
     &            + C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ELSEIF( FGN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCGY =  C(N,NSL)*FGN*FCGP
              ENDIF
              VCNY = 0.D+0
              IF( FNN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCNS = YN(NBS,NSL)/(SN(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBS,NSL)*FCNS)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRNN,ISLC(1) )
                VCNY =  BCX*FNN*THETA*FCNB
     &            + C(N,NSL)*FNN*(1.D+0-THETA)*FCNP
              ELSEIF( FNN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCNY =  C(N,NSL)*FNN*FCNP
              ENDIF
              AN = 0.D+0
              AP = 0.D+0
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCGY+VCNY)/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCGY - VCNY
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPY = NSY(N)
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                SVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/DYGP(NPY)/(SVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/DYGP(NPY)/(SVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNS = AFY(NPY)*VN(1,NPY)
              IF( FNS.LT.ZERO ) THEN
                NBS = N-IFLD
                SVNX = SN(2,NBS)*PORD(2,NBS)
                CRNS = ABS( VN(1,NPY) )*DT/DYGP(NPY)/(SVNX*RP(I)+SMALL)
                FCNS = YN(NBS,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBS,NSL)*FCNS-C(N,NSL)*FCNP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRNS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCNY = C(N,NSL)*FNS*(1.D+0-THETA*DYF)*FCNP
     &            + C(NBS,NSL)*FNS*THETA*DYF*FCNS
                VCNYF = CO(N,NSL)*FNS*FCNP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCNY-VCNYF)/AFY(NPY)
                BLU(MP) = BLU(MP) + VCNY - VCNYF
                BLU(IXP(NBS)) = BLU(IXP(NBS)) - VCNY + VCNYF
              ENDIF
            ELSE
              ALN = MAX( -FLN,ZERO )
              AGN = MAX( -FGN,ZERO )
              ANN = MAX( -FNN,ZERO )
              AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP + (ANN+FNN)*FCNP
              AN = ALN*FCLB + AGN*FCGB + ANN*FCNB
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
              IF( FLN.LT.ZERO ) VCLY = BCX*FCLB*FLN
              VCGY = 0.D+0
              IF( FGN.LT.ZERO ) VCGY = BCX*FCGB*FGN
              VCNY = 0.D+0
              IF( FNN.LT.ZERO ) VCNY = BCX*FCNB*FNN
              AN = 0.D+0
              AP = 0.D+0
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCGY+VCNY)/AFY(NQY)
              BLU(MP) = BLU(MP) - VCLY - VCGY - VCNY
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPY = NSY(N)
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                SVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/DYGP(NPY)/(SVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
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
                SVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/DYGP(NPY)/(SVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
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
              FNS = AFY(NPY)*VN(1,NPY)
              IF( FNS.LT.ZERO ) THEN
                NBS = N-IFLD
                SVNX = SN(2,NBS)*PORD(2,NBS)
                CRNS = ABS( VN(1,NPY) )*DT/DYGP(NPY)/(SVNX*RP(I)+SMALL)
                FCNS = YN(NBS,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBS,NSL)*FCNS-C(N,NSL)*FCNP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRNS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCNY = C(N,NSL)*FNS*(1.D+0-THETA*DYF)*FCNP
     &            + C(NBS,NSL)*FNS*THETA*DYF*FCNS
                VCNYF = CO(N,NSL)*FNS*FCNP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCNY-VCNYF)/AFY(NPY)
                BLU(MP) = BLU(MP) + VCNY - VCNYF
                BLU(IXP(NBS)) = BLU(IXP(NBS)) - VCNY + VCNYF
              ENDIF
            ELSE
              ALN = MAX( -FLN,ZERO )
              AGN = MAX( -FGN,ZERO )
              ANN = MAX( -FNN,ZERO )
              AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP + (ANN+FNN)*FCNP
              AN = ALN*FCLB + AGN*FCGB + ANN*FCNB
            ENDIF
            BLU(MP) = BLU(MP) + AN*BCX
          ENDIF
!
!---  Top boundary
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          NQZ = NSZ(N) + IJFLD
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULTX,VLTX,WLTX,N,MF )
            CALL SHDP( WLTX,ULTX,VLTX,DISPL(IZN),DISPT(IZN),DPLT )
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNTX,VNTX,WNTX,N,MF )
            CALL SHDP( WNTX,UNTX,VNTX,DISPL(IZN),DISPT(IZN),DPNT )
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGTX,VGTX,WGTX,N,MF )
            CALL SHDP( WGTX,UGTX,VGTX,DISPL(IZN),DISPT(IZN),DPGT )
          ELSE
            DPLT = 0.D+0
            DPNT = 0.D+0
            DPGT = 0.D+0
          ENDIF
          FLT = AFZ(NQZ)*WL(1,NQZ)
          FGT = AFZ(NQZ)*WG(1,NQZ)
          FNT = AFZ(NQZ)*WN(1,NQZ)
          CRLT = ABS( WL(1,NQZ) )*DT/(DZGF(N)*SVLB+SMALL)
          CRGT = ABS( WG(1,NQZ) )*DT/(DZGF(N)*SVGB+SMALL)
          CRNT = ABS( WN(1,NQZ) )*DT/(DZGF(N)*SVNB+SMALL)
!
!---  Dirichlet ---
!
          IF( IBCT(NSL+LUK,NB).EQ.1 .OR. IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.9 .OR. IBCT(NSL+LUK,NB).EQ.10 .OR.
     &      IBCT(NSL+LUK,NB).EQ.12 ) THEN
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
            SMDNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SNB(2,NB)*PORDB(2,NB)*SMDNB
            INDX = 16
            DNZ = DIFMN(DNP,DNB,DZGF(N),DZGF(N),WN(1,NQZ),INDX)
            DNZ = AFZ(NQZ)*(DNZ+DPNT)/(5.D-1*DZGF(N))
!
!--- TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 ) THEN
              IF( FLT.LT.ZERO ) THEN
                WCLZ = BCX*FCLB*FLT
              ELSEIF( FLT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCLB = YL(NBB,NSL)/(SL(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBB,NSL)*FCLB)
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                WCLZ =  C(N,NSL)*FLT*(1.D+0-THETA)*FCLP
     &            + BCX*FLT*THETA*FCLB
              ELSEIF( FLT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCLZ =  C(N,NSL)*FLT*FCLP
              ENDIF
              IF( FGT.LT.ZERO ) THEN
                WCGZ = BCX*FCGB*FGT
              ELSEIF( FGT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCGB = YG(NBB,NSL)/(SG(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBB,NSL)*FCGB)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                WCGZ =  C(N,NSL)*FGT*(1.D+0-THETA)*FCGP
     &            + BCX*FGT*THETA*FCGB
              ELSEIF( FGT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCGZ =  C(N,NSL)*FGT*FCGP
              ENDIF
              IF( FNT.LT.ZERO ) THEN
                WCNZ = BCX*FCNB*FNT
              ELSEIF( FNT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCNB = YN(NBB,NSL)/(SN(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBB,NSL)*FCNB)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRNT,ISLC(1) )
                WCNZ =  C(N,NSL)*FNT*(1.D+0-THETA)*FCNP
     &            + BCX*FNT*THETA*FCNB
              ELSEIF( FNT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCNZ =  C(N,NSL)*FNT*FCNP
              ENDIF
              AT = DLZ*FCLB + DGZ*FCGB + DNZ*FCNB
              AP = DLZ*FCLP + DGZ*FCGP + DNZ*FCNP
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCGZ+WCNZ)/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCGZ - WCNZ
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPZ = NSZ(N)
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/DZGP(NPZ)/(SVLX+SMALL)
                FCLB = YL(NBB,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
     &            /(C(NBB,NSL)*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCLZ = C(N,NSL)*FLB*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBB,NSL)*FLB*THETA*DZF*FCLB
                WCLZF = CO(N,NSL)*FGT*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ-WCLZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCLZ - WCLZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCLZ + WCLZF
              ENDIF
              FGT = AFZ(NPZ)*WG(1,NPZ)
              IF( FGT.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/DZGP(NPZ)/(SVGX+SMALL)
                FCGB = YG(NBB,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
     &            /(C(NBB,NSL)*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCGZ = C(N,NSL)*FGT*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBB,NSL)*FGT*THETA*DZF*FCGB
                WCGZF = CO(N,NSL)*FGT*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCGZ-WCGZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCGZ - WCGZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCGZ + WCGZF
              ENDIF
              FNB = AFZ(NPZ)*WN(1,NPZ)
              IF( FNB.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVNX = SN(2,NBB)*PORD(2,NBB)
                CRNB = ABS( WN(1,NPZ) )*DT/DZGP(NPZ)/(SVNX+SMALL)
                FCNB = YN(NBB,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBB,NSL)*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRNB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCNZ = C(N,NSL)*FNB*(1.D+0-THETA*DZF)*FCNP
     &            + C(NBB,NSL)*FNB*THETA*DZF*FCNB
                WCNZF = CO(N,NSL)*FNB*FCNP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCNZ-WCNZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCNZ - WCNZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCNZ + WCNZF
              ENDIF
            ELSE
              ALT = MAX( -FLT,ZERO ) +
     &          DLZ*MAX((ONE-(TENTH*ABS(FLT)/(DLZ+SMALL)))**5,ZERO)
              AGT = MAX( -FGT,ZERO ) +
     &          DGZ*MAX((ONE-(TENTH*ABS(FGT)/(DGZ+SMALL)))**5,ZERO)
              ANT = MAX( -FNT,ZERO ) +
     &          DNZ*MAX((ONE-(TENTH*ABS(FNT)/(DNZ+SMALL)))**5,ZERO)
              AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP + (ANT+FNT)*FCNP
              AT = ALT*FCLB + AGT*FCGB + ANT*FCNB
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
     &            /(BCX*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                WCLZ =  C(N,NSL)*FLT*(1.D+0-THETA)*FCLP
     &            + BCX*FLT*THETA*FCLB
              ELSEIF( FLT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCLZ =  C(N,NSL)*FLT*FCLP
              ENDIF
              IF( FGT.LT.ZERO ) THEN
                WCGZ = BCX*FCGB*FGT
              ELSEIF( FGT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCGB = YG(NBB,NSL)/(SG(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBB,NSL)*FCGB)
     &            /(BCX*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                WCGZ =  C(N,NSL)*FGT*(1.D+0-THETA)*FCGP
     &            + BCX*FGT*THETA*FCGB
              ELSEIF( FGT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCGZ =  C(N,NSL)*FGT*FCGP
              ENDIF
              WCNZ = 0.D+0
              IF( FNT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCNB = YN(NBB,NSL)/(SN(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCNP-C(NBB,NSL)*FCNB)
     &            /(BCX*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRNT,ISLC(1) )
                WCNZ =  C(N,NSL)*FNT*(1.D+0-THETA)*FCNP
     &            + BCX*FNT*THETA*FCNB
              ELSEIF( FNT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCNZ =  C(N,NSL)*FNT*FCNP
              ENDIF
              AT = 0.D+0
              AP = 0.D+0
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCGZ+WCNZ)/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCGZ - WCNZ
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPZ = NSZ(N)
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/DZGP(NPZ)/(SVLX+SMALL)
                FCLB = YL(NBB,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
     &            /(C(NBB,NSL)*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCLZ = C(N,NSL)*FLB*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBB,NSL)*FLB*THETA*DZF*FCLB
                WCLZF = CO(N,NSL)*FGT*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ-WCLZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCLZ - WCLZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCLZ + WCLZF
              ENDIF
              FGT = AFZ(NPZ)*WG(1,NPZ)
              IF( FGT.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/DZGP(NPZ)/(SVGX+SMALL)
                FCGB = YG(NBB,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
     &            /(C(NBB,NSL)*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCGZ = C(N,NSL)*FGT*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBB,NSL)*FGT*THETA*DZF*FCGB
                WCGZF = CO(N,NSL)*FGT*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCGZ-WCGZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCGZ - WCGZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCGZ + WCGZF
              ENDIF
              FNB = AFZ(NPZ)*WN(1,NPZ)
              IF( FNB.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVNX = SN(2,NBB)*PORD(2,NBB)
                CRNB = ABS( WN(1,NPZ) )*DT/DZGP(NPZ)/(SVNX+SMALL)
                FCNB = YN(NBB,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBB,NSL)*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRNB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCNZ = C(N,NSL)*FNB*(1.D+0-THETA*DZF)*FCNP
     &            + C(NBB,NSL)*FNB*THETA*DZF*FCNB
                WCNZF = CO(N,NSL)*FNB*FCNP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCNZ-WCNZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCNZ - WCNZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCNZ + WCNZF
              ENDIF
            ELSE
              ALT = MAX( -FLT,ZERO )
              AGT = MAX( -FGT,ZERO )
              ANT = MAX( -FNT,ZERO )
              AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP + (ANT+FNT)*FCNP
              AT = ALT*FCLB + AGT*FCGB + ANT*FCNB
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
              IF( FLT.LT.ZERO ) WCLZ = BCX*FCLB*FLT
              WCGZ = 0.D+0
              IF( FGT.LT.ZERO )  WCGZ = BCX*FCGB*FGT
              WCNZ = 0.D+0
              IF( FNT.LT.ZERO ) WCNZ = BCX*FCNB*FNT
              AT = 0.D+0
              AP = 0.D+0
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCGZ+WCNZ)/AFZ(NQZ)
              BLU(MP) = BLU(MP) - WCLZ - WCGZ - WCNZ
!
!--- TVD Transport for interior surface adjacent to boundary  ---
!
              NPZ = NSZ(N)
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/DZGP(NPZ)/(SVLX+SMALL)
                FCLB = YL(NBB,NSL)/(SVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX*FCLB)
     &            /(C(NBB,NSL)*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCLZ = C(N,NSL)*FLB*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBB,NSL)*FLB*THETA*DZF*FCLB
                WCLZF = CO(N,NSL)*FGT*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ-WCLZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCLZ - WCLZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCLZ + WCLZF
              ENDIF
              FGT = AFZ(NPZ)*WG(1,NPZ)
              IF( FGT.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/DZGP(NPZ)/(SVGX+SMALL)
                FCGB = YG(NBB,NSL)/(SVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX*FCGB)
     &            /(C(NBB,NSL)*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCGZ = C(N,NSL)*FGT*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBB,NSL)*FGT*THETA*DZF*FCGB
                WCGZF = CO(N,NSL)*FGT*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCGZ-WCGZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCGZ - WCGZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCGZ + WCGZF
              ENDIF
              FNB = AFZ(NPZ)*WN(1,NPZ)
              IF( FNB.LT.ZERO ) THEN
                NBB = N-IJFLD
                SVNX = SN(2,NBB)*PORD(2,NBB)
                CRNB = ABS( WN(1,NPZ) )*DT/DZGP(NPZ)/(SVNX+SMALL)
                FCNB = YN(NBB,NSL)/(SVNX+SMALL)
                R = ((C(N,NSL)*FCNP-BCX*FCNB)
     &            /(C(NBB,NSL)*FCNB-C(N,NSL)*FCNP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRNB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCNZ = C(N,NSL)*FNB*(1.D+0-THETA*DZF)*FCNP
     &            + C(NBB,NSL)*FNB*THETA*DZF*FCNB
                WCNZF = CO(N,NSL)*FNB*FCNP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCNZ-WCNZF)/AFZ(NPZ)
                BLU(MP) = BLU(MP) + WCNZ - WCNZF
                BLU(IXP(NBB)) = BLU(IXP(NBB)) - WCNZ + WCNZF
              ENDIF
            ELSE
              ALT = MAX( -FLT,ZERO )
              AGT = MAX( -FGT,ZERO )
              ANT = MAX( -FNT,ZERO )
              AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP + (ANT+FNT)*FCNP
              AT = ALT*FCLB + AGT*FCGB + ANT*FCNB
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
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SBND5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFT5( RHOLOX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
#ifdef awtap
!
!     This software contains components that were developed solely by
!     and whose rights are retained exclusively by the Battelle
!     Memorial Institute.
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
#endif
!
!----------------------Description-------------------------------------!
!
!     Water-Air-Oil Mode
!
!     Gas-aqueous interfacial tension and contact angle as a function
!     of dissolved-oil concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 April 2005.
!     Last Modified by MD White, PNNL, 6 April 2005.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
      USE HYST
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
      SUB_LOG(ISUB_LOG) = '/SFT5'
!
!---  Gas-aqueous interfacial tension as a function of dissolved-oil
!     concentration  ---
!
      IF( ABS(SFCSF(1)).GT.EPSL ) THEN
        SIG_GL = SFCSF(1) + SFCSF(2)*EXP( -SFCSF(3)*RHOLOX )
        IF( ABS(SFCSF(4)).GT.EPSL ) THEN
          SIG_GL = SIG_GL + SFCSF(4)*EXP( -SFCSF(5)*RHOLOX )
        ENDIF
!
!---    NAPL-aqueous interfacial tension as a function of gas-NAPL
!       contact angle, gas-NAPL surface tension, and gas-aqueous
!       surface tension  ---
!
        SIG_NL = SQRT( -2.D+0*CA_GN*SIG_GL*SIG_GN + 
     &    (SIG_GL**2) + (SIG_GN**2) )
!
!---    NAPL-aqueous contact angle as a function of NAPL-aqueous,
!       gas-NAPL, and gas-aqueous interfacial tensions  ---
!
        CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &    /(2.D+0*SIG_GL*SIG_NL)
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFT5 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORC5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     Last Modified by MD White, Battelle, 2 October 2003.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 SRX(8),HCFG(LFZ)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORC5'
!
!---  Zero source terms  ---
!
      DO 60 NS = 1,NSR
        DO 50 I = ISRDM(1,NS),ISRDM(2,NS)
          DO 50 J = ISRDM(3,NS),ISRDM(4,NS)
            DO 50 K = ISRDM(5,NS),ISRDM(6,NS)
              N = ND(I,J,K)
              DO 40 M = 2,ISVC+2
                SRCW(M,N) = 0.D+0
                SRCA(M,N) = 0.D+0
                SRCO(M,N) = 0.D+0
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
!---  Skip over domain loop for well sources  ---
!
      IF( ISRT(NS).GT.20 .AND. ISRT(NS).LT.30 ) GOTO 510
!
!---  Loop source domain  ---
!
        DO 500 I = ISRDM(1,NS),ISRDM(2,NS)
        DO 500 J = ISRDM(3,NS),ISRDM(4,NS)
        DO 500 K = ISRDM(5,NS),ISRDM(6,NS)
          N = ND(I,J,K)
          IZN = IZ(N)
          IF( IXP(N).EQ.0 ) GOTO 500
          DO 400 M = 2,ISVC+2
            PGX = PG(M,N) + PATM
!
!---        Aqueous Volumetric Rate w/ Component Mass Fractions ---
!
            IF( ISRT(NS).EQ.3 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
               SRCW(M,N) = SRCW(M,N) +
     &           SRX(4)*RHOL(M,N)*(1.D+0-SRX(5)-SRX(6))
               SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOL(M,N)*SRX(5)
               SRCO(M,N) = SRCO(M,N) + SRX(4)*RHOL(M,N)*SRX(6)
              ELSE
               SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOL(M,N)*XLW(M,N)
               SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOL(M,N)*XLA(M,N)
               SRCO(M,N) = SRCO(M,N) + SRX(4)*RHOL(M,N)*XLO(M,N)
              ENDIF
!
!---         Gas Volumetric Rate w/ Component Mass Fractions ---
!
            ELSEIF( ISRT(NS).EQ.4 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                XGWX = SRX(5)
                XGOX = SRX(6)
                XGAX = MAX(1.D+0-XGWX-XGOX,ZERO)
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOG(M,N)*XGWX
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOG(M,N)*XGAX
                SRCO(M,N) = SRCO(M,N) + SRX(4)*RHOG(M,N)*XGOX
              ELSE
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOG(M,N)*XGW(M,N)
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOG(M,N)*XGA(M,N)
                SRCO(M,N) = SRCO(M,N) + SRX(4)*RHOG(M,N)*XGO(M,N)
              ENDIF
!
!---        Gas Volumetric Rate w/ Component Relative Humidity ---
!
            ELSEIF( ISRT(NS).EQ.5 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                PVWX = PSW(2,N)*SRX(5)
                PVOX = PSO(2,N)*SRX(6)
                PVAX = MAX( PGX-PVWX-PVOX,ZERO )
                INDX = 0
                CALL WATGSD( T(2,N),PVWX,RHOW,INDX )
                CALL VOCGSD( T(2,N),PVOX,RHOO,INDX )
                CALL AIRGSD( T(2,N),PVAX,RHOA )
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOW
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOA
                SRCO(M,N) = SRCO(M,N) + SRX(4)*RHOO
              ELSE
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOG(M,N)*XGW(M,N)
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOG(M,N)*XGA(M,N)
                SRCO(M,N) = SRCO(M,N) + SRX(4)*RHOG(M,N)*XGO(M,N)
              ENDIF
!
!---         NAPL Volumetric Rate  ---
!
            ELSEIF( ISRT(NS).EQ.6 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                SRCO(M,N) = SRCO(M,N) + SRX(4)*RHON(M,N)
              ELSE
                SRCO(M,N) = SRCO(M,N) + SRX(4)*RHON(M,N)
              ENDIF
!
!---        Aqueous Mass Rate w/ Component Mass Fractions ---
!
            ELSEIF( ISRT(NS).EQ.7 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                SRCW(M,N) = SRCW(M,N) + SRX(4)*(1.D+0-SRX(5)-SRX(6))
                SRCA(M,N) = SRCA(M,N) + SRX(4)*SRX(5)
                SRCO(M,N) = SRCO(M,N) + SRX(4)*SRX(6)
              ELSE
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XLW(M,N)
                SRCA(M,N) = SRCA(M,N) + SRX(4)*XLA(M,N)
                SRCO(M,N) = SRCO(M,N) + SRX(4)*XLO(M,N)
              ENDIF
!
!---        Gas Mass Rate w/ Component Mass Fractions ---
!
            ELSEIF( ISRT(NS).EQ.8 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                XGWX = SRX(5)
                XGOX = SRX(6)
                XGAX = MAX(1.D+0-XGWX-XGOX,ZERO)
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGW(M,N)
                SRCA(M,N) = SRCA(M,N) + SRX(4)*XGA(M,N)
                SRCO(M,N) = SRCO(M,N) + SRX(4)*XGO(M,N)
              ELSE
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGW(M,N)
                SRCA(M,N) = SRCA(M,N) + SRX(4)*XGA(M,N)
                SRCO(M,N) = SRCO(M,N) + SRX(4)*XGO(M,N)
              ENDIF
!
!---        Gas Mass Rate w/ Component Relative Humidity ---
!
            ELSEIF( ISRT(NS).EQ.9 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                PVWX = PSW(2,N)*SRX(5)
                PVOX = PSO(2,N)*SRX(6)
                PVAX = MAX( PGX-PVWX-PVOX,ZERO )
                INDX = 0
                CALL WATGSD( T(2,N),PVWX,RHOW,INDX )
                CALL VOCGSD( T(2,N),PVOX,RHOO,INDX )
                CALL AIRGSD( T(2,N),PVAX,RHOA )
                RHOGX = RHOW + RHOO + RHOA
                SRCW(M,N) = SRCW(M,N) + SRX(4)*RHOW/RHOGX
                SRCA(M,N) = SRCA(M,N) + SRX(4)*RHOA/RHOGX
                SRCO(M,N) = SRCO(M,N) + SRX(4)*RHOO/RHOGX
              ELSE
                SRCW(M,N) = SRCW(M,N) + SRX(4)*XGW(M,N)
                SRCA(M,N) = SRCA(M,N) + SRX(4)*XGA(M,N)
                SRCO(M,N) = SRCO(M,N) + SRX(4)*XGO(M,N)
              ENDIF
!
!---        NAPL Mass Rate ---
!
            ELSEIF( ISRT(NS).EQ.10 ) THEN
              IF( SRX(4).GE.0.D+0 ) THEN
                SRCO(M,N) = SRCO(M,N) + SRX(4)
              ELSE
                SRCO(M,N) = SRCO(M,N) + SRX(4)
              ENDIF
!
!---        Bioslurping well ---
!
            ELSEIF( ISRT(NS).EQ.30 ) THEN
!
!---          Bioslurping tube coordinates ---
!
              XPWX = SRCP(1,NS)
              YPWX = SRCP(2,NS)
              ZPWX = SRCP(3,NS)
!
!---          Bioslurping well screened interval component lengths ---
!
              XLX = SRCP(4,NS)
              YLX = SRCP(5,NS)
              ZLX = SRCP(6,NS)
!
!---          Well gas pressure ---
!
              PGWX = SRX(2) - PATM
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
                XPNX = 0.D+0
                YPNX = 0.D+0
                ZPNX = ZP(N)
!
!---          Cylindrical coordinates  ---
!
              ELSEIF( ICS.EQ.2 .OR. ICS.EQ.6 ) THEN
                XPNX = XP(N)*COS(YP(N))
                YPNX = XP(N)*SIN(YP(N))
                ZPNX = ZP(N)
!
!---          Cartesian or boundary-fitted orthogonal coordinates  ---
!
              ELSE
                XPNX = XP(N)
                YPNX = YP(N)
                ZPNX = ZP(N)
              ENDIF
!
!---          Formation gas pressure at the bioslurping point ---
!
              PGFX = PG(M,N) - (ZPWX-ZPNX)*GRAV)*RHOG(M,N)
!
!---          Formation aqueous pressure at the bioslurping point ---
!
              PLFX = PL(M,N) - (ZPWX-ZPNX)*GRAV)*RHOL(M,N)
!
!---          Formation nonaqueous pressure at the bioslurping point ---
!
              PNFX = PN(M,N) - (ZPWX-ZPNX)*GRAV)*RHON(M,N)
!
!---          Equivalent field node radius components  ---
!
              PERMX = MAX( PERM(1,IZN),1.D-20 )
              PERMY = MAX( PERM(2,IZN),1.D-20 )
              PERMZ = MAX( PERM(3,IZN),1.D-20 )
              RWX = MAX( SRCP(7,NS),1.D-20 )
!
!---          Cylindrical coordinates with azimuthal symmetry,
!             centrally located wells  ---
!
              IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. JFLD.EQ.1 ) THEN
                ROZ = RP(I)
                RWX = MIN( RWX,9.999D-1*ROZ )
                PERMX = PERMRF(M,N)*PERM(1,IZN)
                WI_CWX = 2.D+0*GPI*PERMX*ZLX/(LOG(ROZ/RWX)+SRCP(8,NS))
              ELSE
                PERMYZ = SQRT(PERMY/PERMZ)
                PERMZY = SQRT(PERMZ/PERMY)
                DXGFX = DXGF(N)
                DYGFX = DYGF(N)*RP(I)
                DZGFX = DZGF(N)
                ROX = 2.8D-1*SQRT(PERMYZ*(DZGFX**2) + 
     &            PERMZY*(DYGFX**2))/(SQRT(PERMYZ)+SQRT(PERMZY))
                PERMZX = SQRT(PERMZ/PERMX)
                PERMXZ = SQRT(PERMX/PERMZ)
                ROY = 2.8D-1*SQRT(PERMZX*(DXGFX**2) + 
     &            PERMXZ*(DZGFX**2))/(SQRT(PERMZX)+SQRT(PERMXZ))
                PERMYX = SQRT(PERMY/PERMX)
                PERMXY = SQRT(PERMX/PERMY)
                ROZ = 2.8D-1*SQRT(PERMYX*(DXGFX**2) + 
     &            PERMXY*(DYGFX**2))/(SQRT(PERMYX)+SQRT(PERMXY))
!
!---            Well index components  ---
!
                PERMX = PERMRF(M,N)*PERM(1,IZN)
                PERMY = PERMRF(M,N)*PERM(2,IZN)
                PERMZ = PERMRF(M,N)*PERM(3,IZN)
                WIX = 2.D+0*GPI*SQRT(PERMY*PERMZ)*XLX/
     &            (LOG(ROX/RWX)+SRCP(8,NS))
                WIY = 2.D+0*GPI*SQRT(PERMX*PERMZ)*YLX/
     &            (LOG(ROY/RWX)+SRCP(8,NS))
                WIZ = 2.D+0*GPI*SQRT(PERMX*PERMY)*ZLX/
     &            (LOG(ROZ/RWX)+SRCP(8,NS))
                WI_CWX = SQRT((WIX**2) + (WIY**2) + (WIZ**2))
              ENDIF
!
!---          Gas mass fluxes, positive into the well  ---
!
              DPGX = MAX( PGFX-PGWX,0.D+0 )
              DPGX = PGFX-PGWX
!
!---          Gas flux from reservoir to well  ---
!
              QGX = WI_CWX*RKG(M,N)*DPGX/VISG(M,N)
!
!---          Aqueous w/ dissolved-oil, dissolved-air conditions
!             Unknowns: PL,XMLA,XMLO
!
!---          Aqueous w/ dissolved-oil, trapped-gas conditions
!             Unknowns: PL,SG,XMLO
!
!---          Aqueous w/ trapped-NAPL, dissolved-air conditions
!             Unknowns: PL,XMLA,SN
!
!---          Aqueous w/ trapped-NAPL, trapped-gas conditions
!             Unknowns: PL,SG,SN
!
              IF( NPHAZ(2,N).EQ.1 .OR. NPHAZ(2,N).EQ.2 .OR.
     &          NPHAZ(2,N).EQ.3 .OR. NPHAZ(2,N).EQ.4 ) THEN
                DPLX = PLFX-PGWX
                RMLX = ((RKL(1,M,N)*RKL(2,M,N)*RKL(3,M,N))**
     &            (1.D+0/3.D+0))/VISL(M,N)
                QLX = WI_CWX*RMLX*DPLX
                QTX = QLX + QGX
                IF( ABS(QTX).GT.EPSL ) THEN
                  QLX = (QLX/QTX)*MIN( QTX,SRX(3) )
                  QGX = (QGX/QTX)*MIN( QGX,SRX(3) )
                ELSE
                  QLX = 0.D+0
                  QGX = 0.D+0
                ENDIF
!
!---            Source contributions from aqueous slurping
!               and gas flow ---
!
                SRCW(M,N) = SRCW(M,N) - QLX*RHOL(M,N)*XLW(M,N) -
     &            QGX*RHOG(M,N)*XGW(M,N)
                SRCA(M,N) = SRCA(M,N) - QLX*RHOL(M,N)*XLA(M,N) -
     &            QGX*RHOG(M,N)*XGA(M,N)
                SRCO(M,N) = SRCO(M,N) - QLX*RHOL(M,N)*XLO(M,N) -
     &            QGX*RHOG(M,N)*XGO(M,N)
!
!---          Aqueous-gas w/ dissolved-oil conditions
!             Unknowns:  PL,PG,XMLO
!
!---          Aqueous-gas w/ trapped-NAPL conditions
!             Unknowns: PL,PG,SN
!
              ELSEIF( NPHAZ(2,N).EQ.5 .OR. NPHAZ(2,N).EQ.6 ) THEN
                DPLX = MAX( PLFX-PGWX,0.D+0 )
                RMLX = ((RKL(1,M,N)*RKL(2,M,N)*RKL(3,M,N))**
     &            (1.D+0/3.D+0))/VISL(M,N)
                QLX = WI_CWX*RMLX*DPLX
                QTX = QLX + QGX
                IF( ABS(QTX).GT.EPSL ) THEN
                  QLX = (QLX/QTX)*MIN( QTX,SRX(3) )
                  QGX = (QGX/QTX)*MIN( QGX,SRX(3) )
                ELSE
                  QLX = 0.D+0
                  QGX = 0.D+0
                ENDIF
!
!---            Source contributions from aqueous slurping
!               and gas flow ---
!
                SRCW(M,N) = SRCW(M,N) - QLX*RHOL(M,N)*XLW(M,N) -
     &            QGX*RHOG(M,N)*XGW(M,N)
                SRCA(M,N) = SRCA(M,N) - QLX*RHOL(M,N)*XLA(M,N) -
     &            QGX*RHOG(M,N)*XGA(M,N)
                SRCO(M,N) = SRCO(M,N) - QLX*RHOL(M,N)*XLO(M,N) -
     &            QGX*RHOG(M,N)*XGO(M,N)
!
!---          Aqueous-NAPL w/ dissolved-air conditions
!             Unknowns: PL,XMLA,PN
!
!---          Aqueous-NAPL w/ trapped-gas conditions
!             Unknowns: PL,SG,PN
!
              ELSEIF( NPHAZ(2,N).EQ.7 .OR. NPHAZ(2,N).EQ.8 ) THEN
                RMLX = ((RKL(1,M,N)*RKL(2,M,N)*RKL(3,M,N))**
     &            (1.D+0/3.D+0))/VISL(M,N)
                RMNX = RKN(M,N)/VISN(M,N)
                DPLX = MAX( PLFX-PGWX,0.D+0 )
                DPLX = PLFX-PGWX
                QLX = WI_CWX*RMLX*DPLX
                DPNX = MAX( PNFX-PGWX,0.D+0 )
                DPNX = PNFX-PGWX
                QNX = WI_CWX*RMNX*DPNX
                QTX = QLX + QNX + QGX
                IF( ABS(QTX).GT.EPSL ) THEN
                  QLX = (QLX/QTX)*MIN( QTX,SRX(3) )
                  QNX = (QNX/QTX)*MIN( QTX,SRX(3) )
                  QGX = (QGX/QTX)*MIN( QTX,SRX(3) )
                ELSE
                  QLX = 0.D+0
                  QNX = 0.D+0
                  QGX = 0.D+0
                ENDIF
!
!---            Source contributions from aqueous slurping
!               and gas flow ---
!
                SRCW(M,N) = SRCW(M,N) - QLX*RHOL(M,N)*XLW(M,N) -
     &            QGX*RHOG(M,N)*XGW(M,N)
                SRCA(M,N) = SRCA(M,N) - QLX*RHOL(M,N)*XLA(M,N) -
     &            QGX*RHOG(M,N)*XGA(M,N)
                SRCO(M,N) = SRCO(M,N) - QLX*RHOL(M,N)*XLO(M,N) -
     &            QGX*RHOG(M,N)*XGO(M,N)
!
!---            Source contributions from NAPL slurping  ---
!
                SRCO(M,N) = SRCO(M,N) - QNX*RHON(M,N)
!
!---          Aqueous-NAPL-gas conditions
!             Unknowns: PL,PG,PN
!
              ELSEIF( NPHAZ(2,N).EQ.9 ) THEN
                RMLX = ((RKL(1,M,N)*RKL(2,M,N)*RKL(3,M,N))**
     &            (1.D+0/3.D+0))/VISL(M,N)
                RMNX = RKN(M,N)/VISN(M,N)
                DPLX = MAX( PLFX-PGWX,0.D+0 )
                QLX = WI_CWX*RMLX*DPLX
                DPNX = MAX( PNFX-PGWX,0.D+0 )
                QNX = WI_CWX*RMNX*DPNX
                QTX = QLX + QNX + QGX
                IF( ABS(QTX).GT.EPSL ) THEN
                  QLX = (QLX/QTX)*MIN( QTX,SRX(3) )
                  QNX = (QNX/QTX)*MIN( QTX,SRX(3) )
                  QGX = (QGX/QTX)*MIN( QTX,SRX(3) )
                ELSE
                  QLX = 0.D+0
                  QNX = 0.D+0
                  QGX = 0.D+0
                ENDIF
!
!---            Source contributions from aqueous slurping
!               and gas flow ---
!
                SRCW(M,N) = SRCW(M,N) - QLX*RHOL(M,N)*XLW(M,N) -
     &            QGX*RHOG(M,N)*XGW(M,N)
                SRCA(M,N) = SRCA(M,N) - QLX*RHOL(M,N)*XLA(M,N) -
     &            QGX*RHOG(M,N)*XGA(M,N)
                SRCO(M,N) = SRCO(M,N) - QLX*RHOL(M,N)*XLO(M,N) -
     &            QGX*RHOG(M,N)*XGO(M,N)
!
!---            Source contributions from NAPL slurping  ---
!
                SRCO(M,N) = SRCO(M,N) - QNX*RHON(M,N)
              ENDIF
            ENDIF
  400     CONTINUE
  500   CONTINUE
!
!---    Well sources  ---
!
  510   CONTINUE
!
!---    SVE Well Source  ---
!
        IF( ISRT(NS).EQ.22 ) THEN
          I = ISRDM(1,NS)
          J = ISRDM(3,NS)
          K1X = ISRDM(5,NS)
          K2X = ISRDM(6,NS)
          RDW = SRCP(2,NS)
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
              DZX = DZGF(N)
              IF( K1X.EQ.K2X ) THEN
                HCFG(K) = 1.D+0
                SUMKX = 1.D+0
                GOTO 514
              ENDIF
              IF( K.EQ.K1X ) THEN
                IF( SRX(6).GE.ZE(1,ND(I,J,K)) ) 
     &            DZX = ZE(5,ND(I,J,K)) - SRX(6)
              ELSEIF( K.EQ.K2X ) THEN
                IF( SRX(5).GT.ZE(1,ND(I,J,K)) ) 
     &            DZX = SRX(5) - ZE(1,ND(I,J,K))
              ENDIF
              HCFG(K) = PERMX*RKG(M,N)*DZX
              SUMKX = SUMKX + HCFG(K)
  512       CONTINUE
  514       CONTINUE
            DO 516 K = K1X,K2X
              N = ND(I,J,K)
              WTX = HCFG(K)/SUMKX
!
!---          Soil-vapor-extraction well used to inject gas  ---
!
              IF( SRX(4).GE.0.D+0 ) THEN
                XGWX = SRX(2)
                XGOX = SRX(3)
                XGAX = MAX(1.D+0-XGWX-XGOX,ZERO)
                SRCW(M,N) = SRCW(M,N) + WTX*SRX(4)*RHOG(M,N)*XGWX
                SRCA(M,N) = SRCA(M,N) + WTX*SRX(4)*RHOG(M,N)*XGAX
                SRCO(M,N) = SRCO(M,N) + WTX*SRX(4)*RHOG(M,N)*XGOX
!
!---          Soil-vapor-extraction well used to withdraw gas  ---
!
              ELSE
!
!---            Limit extraction rate by minimum gas pressure  ---
!
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
                IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. (I.EQ.1) ) THEN
                  RDE = XE(2,ND(I,J,K))
                ELSE
                  RDE = SQRT( AFZ(NSZ(N))/GPI )
                ENDIF
                ACWX = 2.D+0*GPI*RDW*DZGF(N)
                DRD2 = (RDE**2-RDW**2)
                HCGX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
     &            (VISG(M,N)*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
                DPGX = MIN( SRCP(1,NS)-PG(M,N),0.D+0 )
                QGX = MAX( DPGX*RKG(M,N)*HCGX,WTX*SRX(4) )
                SRCW(M,N) = SRCW(M,N) + QGX*RHOG(M,N)*XGW(M,N)
                SRCA(M,N) = SRCA(M,N) + QGX*RHOG(M,N)*XGA(M,N)
                SRCO(M,N) = SRCO(M,N) + QGX*RHOG(M,N)*XGO(M,N)
              ENDIF
  516       CONTINUE
  518     CONTINUE
        ENDIF
  600 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORC5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORIT5( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute solute transport source integrals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, June 1994.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/SORIT5'
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
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*YL(N,NSL)/
     &            (PORD(2,N)*SL(2,N))
!
!---  Gas Volumetric Sink  ---
!
              ELSEIF( ( ISRT(NS).EQ.4 .OR. ISRT(NS).EQ.5 ) .AND.
     &          SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YG(N,NSL)*DT/(PORD(2,N)*SG(2,N))
!
!---  NAPL Volumetric Sink  ---
!
              ELSEIF( ISRT(NS).EQ.6 .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YN(N,NSL)*DT/(PORD(2,N)*SN(2,N))
!
!---  Aqueous Mass Sink  ---
!
              ELSEIF( ISRT(NS).EQ.7 .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YL(N,NSL)*DT/(RHOL(2,N)*PORD(2,N)*SL(2,N))
!
!---  Gas Volumetric Sink  ---
!
              ELSEIF( ( ISRT(NS).EQ.8 .OR. ISRT(NS).EQ.9 ) .AND.
     &          SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YG(N,NSL)*DT/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---  NAPL Mass Sink  ---
!
              ELSEIF( ISRT(NS).EQ.10 .AND. SRX(4).LT.0.D+0 ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) - C(N,NSL)*SRX(4)*
     &            YN(N,NSL)*DT/(RHON(2,N)*PORD(2,N)*SN(2,N))
!
!---  Solute source  ---
!
              ELSEIF( ISRT(NS).LT.0 .AND. ISRT(NS).GE.-NSOLU ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(4)*DT
!
!---  Solute density source  ---
!
              ELSEIF( ISRT(NS).LT.-NSOLU .AND.
     &          ISRT(NS).GE.-2*NSOLU ) THEN
                SRCIC(N,NSL) = SRCIC(N,NSL) + SRX(4)*DT*VOL(N)
              ENDIF
  500   CONTINUE
  600 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORIT5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORT5( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute solute transport source terms.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/SORT5'
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
          ELSEIF( ( ISRT(NS).EQ.4 .OR. ISRT(NS).EQ.5 ) .AND.
     &      SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YG(N,NSL)/(PORD(2,N)*SG(2,N))
!
!---      NAPL Volumetric Sink  ---
!
          ELSEIF( ISRT(NS).EQ.6 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YN(N,NSL)/(PORD(2,N)*SN(2,N))
!
!---      Aqueous Mass Sink  ---
!
          ELSEIF( ISRT(NS).EQ.7 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YL(N,NSL)/(RHOL(2,N)*PORD(2,N)*SL(2,N))
!
!---      Gas Mass Sink  ---
!
          ELSEIF( ( ISRT(NS).EQ.8 .OR. ISRT(NS).EQ.9 ) .AND.
     &      SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YG(N,NSL)/(RHOG(2,N)*PORD(2,N)*SG(2,N))
!
!---      NAPL Mass Sink  ---
!
          ELSEIF( ISRT(NS).EQ.10 .AND. SRX(4).LT.0.D+0 ) THEN
            SORTX = -SRX(4)*YN(N,NSL)/(RHON(2,N)*PORD(2,N)*SN(2,N))
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
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORT5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPRP5( NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Calculates the aqueous-, NAPL-, and gas-phase solute
!     mole fractions from user-specified partition coefficients.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January 1995.
!     Last Modified by MD White, Battelle, PNL, January 19, 1995.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/SPRP5'
!
!---  Loop over all nodes  ---
!
      DO 900 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 900
        IZN = IZ(N)
        XVS = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
        XVL = SL(2,N)*PORD(2,N)
        XVG = SG(2,N)*PORD(2,N)
        XVN = SN(2,N)*PORD(2,N)
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
        YVL = 1.D+0/(XVS + XVL + XVN/PCLN(1,NSL) + XVG*PCGLX)
        YVG = 1.D+0/((XVS + XVL + XVN)/PCGLX + XVG)
        YVN = 1.D+0/((XVS + XVL + XVG*PCGLX)*PCLN(1,NSL) + XVN)
!
!---  Phase mole fractions  ---
!
        YL(N,NSL) = XVL*YVL
        YG(N,NSL) = XVG*YVG
        YN(N,NSL) = XVN*YVN
!
  900 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPRP5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TMPR5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Compute thermo-physical properties for aqueous, napl, and gas
!     phases.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FDVP
      USE FDVN
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
      SUB_LOG(ISUB_LOG) = '/TMPR5'
!
!---  Loop over nodes and secondary variable indices  ---
!
      DO 300 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 300
        IZN = IZ(N)
        DO 200 M = 2,ISVC+2
          PLX = PL(M,N) + PATM
          PNX = PN(M,N) + PATM
          PGX = PG(M,N) + PATM
          BGL = SCHR(8,IZN)/SCHR(5,IZN)
          BNL = SCHR(8,IZN)/(SCHR(6,IZN)*CA_NL)
          BGN = SCHR(8,IZN)/(SCHR(7,IZN)*CA_GN)
          ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---      Aqueous w/ dissolved-oil, dissolved-air conditions
!         Unknowns: PL,XMLA,XMLO
!
!         Kinetic volatilization option
!         Unknowns: PL,XMLA,PN,XMLO  ---
!
          IF( NPHAZ(2,N).EQ.1 ) THEN
            PCAP = MAX( PGX-PLX-ENPR/BGL,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),T(2,N),RHORL,RCW,PVW(M,N) )
            ENDIF
            PVO(M,N) = XMLO(M,N)*HCOW
            PVA(M,N) = XMLA(M,N)*HCAW
            PTX = MAX( PLX,PGX,PNX )
!
!---      Aqueous w/ dissolved-oil, trapped-gas conditions
!         Unknowns: PL,SG,XMLO
!
!         Kinetic volatilization option
!         Unknowns: PL,SG,PN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
            PCAP = MAX( PGX-PLX-ENPR/BGL,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),T(2,N),RHORL,RCW,PVW(M,N) )
            ENDIF
            PVO(M,N) = XMLO(M,N)*HCOW
            PVA(M,N) = MAX( PGX-PVO(M,N)-PVW(M,N),ZERO )
            XMLA(M,N) = PVA(M,N)/HCAW
            PTX = MAX( PLX,PGX,PNX )
!
!---      Aqueous w/ trapped-NAPL, dissolved-air conditions
!         Unknowns: PL,XMLA,SN
!
!         Kinetic volatilization option
!         Unknowns: PL,XMLA,SN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
            PCAP = MAX( PGX-PLX-ENPR/BGL,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),T(2,N),RHORL,RCW,PVW(M,N) )
            ENDIF
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              PVO(M,N) = XMLO(M,N)*HCOW
            ELSE
              PVO(M,N) = PSO(2,N)
              XMLO(M,N) = PVO(M,N)/HCOW
            ENDIF
            PVA(M,N) = XMLA(M,N)*HCAW
            PTX = MAX( PLX,PGX,PNX )
!
!---      Aqueous w/ trapped-NAPL, trapped-gas conditions
!         Unknowns: PL,SG,SN
!
!         Kinetic volatilization option
!         Unknowns: PL,SG,SN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
            PCAP = MAX( PGX-PLX-ENPR/BGL,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),T(2,N),RHORL,RCW,PVW(M,N) )
            ENDIF
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              PVO(M,N) = XMLO(M,N)*HCOW
            ELSE
              PVO(M,N) = PSO(2,N)
              XMLO(M,N) = PVO(M,N)/HCOW
            ENDIF
            PVA(M,N) = MAX( PGX-PVO(M,N)-PVW(M,N),ZERO )
            XMLA(M,N) = PVA(M,N)/HCAW
            PTX = MAX( PLX,PGX,PNX )
!
!---      Aqueous-gas w/ dissolved-oil conditions
!         Unknowns:  PL,PG,XMLO
!
!         Kinetic volatilization option
!         Unknowns: PL,XMLA,PN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
            PCAP = MAX( PGX-PLX-ENPR/BGL,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),T(2,N),RHORL,RCW,PVW(M,N) )
            ENDIF
            PVO(M,N) = XMLO(M,N)*HCOW
            PVA(M,N) = MAX( PGX-PVO(M,N)-PVW(M,N),ZERO )
            XMLA(M,N) = PVA(M,N)/HCAW
            PTX = MAX( PLX,PGX,PNX )
!
!---      Aqueous-gas w/ trapped-NAPL conditions
!         Unknowns: PL,PG,SN
!
!         Kinetic volatilization option
!         Unknowns: PL,PG,SN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
            PCAP = MAX( PGX-PLX-ENPR/BGL,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),T(2,N),RHORL,RCW,PVW(M,N) )
            ENDIF
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              PVO(M,N) = XMLO(M,N)*HCOW
            ELSE
              PVO(M,N) = PSO(2,N)
              XMLO(M,N) = PVO(M,N)/HCOW
            ENDIF
            PVA(M,N) = MAX( PGX-PVO(M,N)-PVW(M,N),ZERO )
            XMLA(M,N) = PVA(M,N)/HCAW
            PTX = MAX( PLX,PGX,PNX )
!
!---      Aqueous-NAPL w/ dissolved-air conditions
!         Unknowns: PL,XMLA,PN
!
!         Kinetic volatilization option
!         Unknowns: PL,XMLA,PN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.7 ) THEN
            PCAP = MAX( PNX-(PLX*BNL-ENPR)/BNL,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),T(2,N),RHORL,RCW,PVW(M,N) )
            ENDIF
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              PVO(M,N) = XMLO(M,N)*HCOW
            ELSE
              PVO(M,N) = PSO(2,N)
              XMLO(M,N) = PVO(M,N)/HCOW
            ENDIF
            PVA(M,N) = XMLA(M,N)*HCAW
            PTX = MAX( PLX,PGX,PNX )
!
!---      Aqueous-NAPL w/ trapped-gas conditions
!         Unknowns: PL,SG,PN
!
!         Kinetic volatilization option
!         Unknowns: PL,SG,PN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.8 ) THEN
            PCAP = MAX( PNX-(PLX*BNL-ENPR)/BNL,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),T(2,N),RHORL,RCW,PVW(M,N) )
            ENDIF
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              PVO(M,N) = XMLO(M,N)*HCOW
            ELSE
              PVO(M,N) = PSO(2,N)
              XMLO(M,N) = PVO(M,N)/HCOW
            ENDIF
            PVA(M,N) = MAX( PGX-PVO(M,N)-PVW(M,N),ZERO )
            XMLA(M,N) = PVA(M,N)/HCAW
            PTX = MAX( PLX,PGX,PNX )
!
!---      Aqueous-NAPL-gas conditions
!         Unknowns: PL,PG,PN
!
!         Kinetic volatilization option
!         Unknowns: PL,PG,PN,XMLO  ---
!
          ELSEIF( NPHAZ(2,N).EQ.9 ) THEN
            PCAP = MAX( PNX-(PLX*BNL-ENPR)/BNL,ZERO )
            IF( ISLC(44).EQ.1 ) THEN
              PVW(M,N) = PSW(2,N)
            ELSE
              CALL VPLWR( PCAP,PSW(2,N),T(2,N),RHORL,RCW,PVW(M,N) )
            ENDIF
            PCAP = MAX( PGX-PNX-ENPR/BGN,ZERO )
!            IF( ISLC(44).EQ.1 ) THEN
!              PVO(M,N) = PSO(2,N)
!            ELSE
!              CALL VPLWR( PCAP,PSO(2,N),T(2,N),RHORN,RCO,PVO(M,N) )
!            ENDIF
!
!---        Kinetic volatilization option  ---
!
            IF( ISLC(46).GE.1 ) THEN
              PVO(M,N) = XMLO(M,N)*HCOW
            ELSE
              PVO(M,N) = PSO(2,N)
              XMLO(M,N) = PVO(M,N)/HCOW
            ENDIF
            PVA(M,N) = MAX( PGX-PVO(M,N)-PVW(M,N),ZERO )
            XMLA(M,N) = PVA(M,N)/HCAW
            PTX = MAX( PLX,PGX,PNX )
          ENDIF
!
!---     Gas-phase component vapor concentration, mole and
!        mass fractions, density, molar density  ---
!
          PVX = PVW(M,N) + PVA(M,N) + PVO(M,N)
          INDX = 0
          CALL WATGSD( T(2,N),PVW(M,N),RHOW,INDX )
          CALL VOCGSD( T(2,N),PVO(M,N),RHOO,INDX )
          CALL AIRGSD( T(2,N),PVA(M,N),RHOA )
          RHOG(M,N) = RHOW+RHOO+RHOA
          XGW(M,N) = RHOW/RHOG(M,N)
          XGO(M,N) = RHOO/RHOG(M,N)
          XGA(M,N) = RHOA/RHOG(M,N)
          XMGW(M,N) = PVW(M,N)/PVX
          XMGA(M,N) = PVA(M,N)/PVX
          XMGO(M,N) = PVO(M,N)/PVX
          WTMX = XMGW(M,N)*WTMW + XMGA(M,N)*WTMA + XMGO(M,N)*WTMO
          RHOMG(M,N) = RHOG(M,N)/WTMX
!
!---      Aqueous-phase mass fraction  ---
!
          XMLW(M,N) = MAX( 1.D+0-XMLA(M,N)-XMLO(M,N),ZERO )
          WTMX = XMLW(M,N)*WTMW + XMLA(M,N)*WTMA + XMLO(M,N)*WTMO
          XLW(M,N) = XMLW(M,N)*WTMW/WTMX
          XLA(M,N) = XMLA(M,N)*WTMA/WTMX
          XLO(M,N) = XMLO(M,N)*WTMO/WTMX
!
!---      NAPL mass fraction  ---
!
          XNO(M,N) = 1.D+0
!
!---      Aqueous-phase density and viscosity  ---
!
          CALL WATLQD( T(2,N),PTX,RHOL(M,N) )
          RHOML(M,N) = RHOL(M,N)/WTMX
          CALL WATLQV( T(2,N),PTX,PSW(2,N),VISL(M,N) )
!
!---      NAPL density and viscosity  ---
!
          CALL VOCLQD( T(2,N),PTX,PVO(M,N),RHON(M,N) )
          CALL VOCLQV( T(2,N),VISN(M,N),VISCO,IVISO )
!
!---      Gas component viscosity and phase viscosity  ---
!
          CALL WATGSV( T(2,N),VISW )
          CALL AIRGSV( T(2,N),VISA )
          CALL VOCGSV( T(2,N),VISO )
          CALL GASVIS( XMGW(M,N),XMGO(M,N),XMGA(M,N),VISW,VISO,VISA,
     &      VISG(M,N) )
!
!---      Compute soil mass fraction of sorbed oil  ---
!
          CLOX = XLO(M,N)*RHOL(M,N)
          RHOBX = RHOS(IZN)*(1.D+0-PORT(M,N))
          CALL SORB( CLOX,RHOBX,XSO(M,N),IZN )
!
!---      Gas diffusion coefficients  ---
!
          IF( ISLC(2).EQ.1 ) THEN
            DFGO(M,N) = DFGOC
            DFGW(M,N) = DFGWC
          ELSEIF( ISLC(2).EQ.2 ) THEN
            CALL BNDFAW( T(2,N),PGX,DFGAW )
            CALL BNDFAO( T(2,N),PGX,DFGAO )
            CALL BNDFOW( T(2,N),PGX,DFGOW )
            CALL VOCDFG( XMGW(M,N),XMGO(M,N),XMGA(M,N),
     &        DFGAO,DFGOW,DFGO(M,N) )
            CALL WATDFG( XMGW(M,N),XMGO(M,N),XMGA(M,N),
     &        DFGAW,DFGOW,DFGW(M,N) )
          ENDIF
!
!---      Aqueous diffusion coefficients  ---
!
          IF( ISLC(4).EQ.1 ) THEN
            DFLO(M,N) = DFLOC
            DFLA(M,N) = DFLAC
          ELSEIF( ISLC(4).EQ.2 ) THEN
            CALL AIRDFL( T(2,N),VISL(M,N),DFLA(M,N) )
            CALL VOCDFL( T(2,N),VISL(M,N),DFLO(M,N) )
          ENDIF
  200   CONTINUE
  300 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TMPR5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRPNL5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Air-Oil Mode
!     Compute the total trapping number for NAPL entrapment in the
!     aqueous phase.
!
!     Pennell, K.D., G.A. Pope, L.M. Abriola.  1996.
!     "Influence of Viscous and Buoyancy Forces on the Mobilization
!     of Residual Tetrachloroethylene during Surfactant Flushing."
!     Environ. Sci. Technol.  30(4):1328-1335.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1997.
!     Last Modified by MD White on March 13, 1997.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVD
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
      SUB_LOG(ISUB_LOG) = '/TRPNL5'
!
!---  Loop over all nodes  ---
!
      DO 100 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 100
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N)+1
        NQY = NSY(N)+IFLD
        NQZ = NSZ(N)+IJFLD
        IZN = IZ(N)
        DFMLX = SQRT( ABS(UL(1,NSX(N)))*ABS(UL(1,NSX(N)+1)) +
     &    ABS(VL(1,NSY(N)))*ABS(VL(1,NSY(N)+IFLD)) +
     &    ABS(WL(1,NSZ(N)))*ABS(WL(1,NSZ(N)+IJFLD)))
        ULX = 0.5D+0*(UL(1,NSX(N))+UL(1,NSX(N)+1))
        VLX = 0.5D+0*(VL(1,NSY(N))+VL(1,NSY(N)+IFLD))
        WLX = 0.5D+0*(WL(1,NSZ(N))+WL(1,NSZ(N)+IJFLD))
        ULGX = 5.D-1*(UL(1,NPX)*GRVX(NPX)+UL(1,NQX)*GRVX(NQX))
        VLGX = 5.D-1*(VL(1,NPY)*GRVY(NPY)+VL(1,NQY)*GRVY(NQY))
        WLGX = 5.D-1*(WL(1,NPZ)*GRVZ(NPZ)+WL(1,NQZ)*GRVZ(NQZ))
        DFALX = ((ULGX + VLGX + WLGX)/GRAV)/
     &    ( SQRT( ULX**2 + VLX**2 + WLX**2 ) + SMALL )
        SKL = SQRT((PERM(1,IZN)*(5.D-1*(GRVX(NPX)+GRVX(NQX))))**2 + 
     &    (PERM(2,IZN)*5.D-1*(GRVY(NPY)+GRVY(NQY)))**2 +
     &    (PERM(3,IZN)*5.D-1*(GRVZ(NPZ)+GRVZ(NQZ)))**2)/GRAV
        BGL = SCHR(8,IZN)/SCHR(5,IZN)
        BNL = SCHR(8,IZN)/(SCHR(6,IZN)*CA_NL)
        BGN = SCHR(8,IZN)/(SCHR(7,IZN)*CA_GN)
        RKLMX = (RKL(1,2,N)*RKL(2,2,N)*RKL(3,2,N))**(THIRD)
        BNDX = (RHOL(2,N)-RHON(2,N))*GRAV*SKL*RKLMX*BNL/SCHR(8,IZN)
        CAPX = DFMLX*VISL(2,N)*SCHR(8,IZN)/BNL
        TRPNL(2,N) = SQRT(CAPX**2 + 2.D+0*CAPX*BNDX*DFALX + BNDX**2)
  100 CONTINUE
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TRPNL5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDT5
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Water-Oil-Air Mode
!     Update the primary variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by MD White, Battelle, June 21, 1994.
!     Last Modified by MD White, PNNL, 27 June 2002.
!     stomp_woa.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE OUTPU
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
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
      SUB_LOG(ISUB_LOG) = '/UPDT5'
!
!---  Skip for convergence error  ---
!
      IF( ICNV.EQ.1 ) GOTO 300
      IERR = 0
!
!---  Update primary variables
!
      DO 200 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 200
        IZN = IZ(N)
        NMD = IXP(N)
        MWEX = IM(IEQW,NMD)
        MAEX = IM(IEQA,NMD)
        BGL = SCHR(8,IZN)/SCHR(5,IZN)
        BNL = SCHR(8,IZN)/(SCHR(6,IZN)*CA_NL)
        BGN = SCHR(8,IZN)/(SCHR(7,IZN)*CA_GN)
        ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---    Aqueous dissolved-oil solubility limit
!
        XMLOX = PSO(2,N)/HCOW
!
!---    Aqueous w/ dissolved-oil, dissolved-air conditions
!       Unknowns: PL,XMLA,XMLO
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,PN,XMLO  ---
!
        IF( NPHAZ(2,N).EQ.1 ) THEN
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) IERR = 1
          IF( BGL*(PG(2,N)-PL(2,N)-BLU(MWEX)).GT.ENPR ) THEN
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MWEX)
          ELSE
            PL(2,N) = PL(2,N) + BLU(MWEX)
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            MDOEX = IM(IEQDO,NMD)
            MNOEX = IM(IEQO,NMD)
            PN(2,N) = PN(2,N) + BLU(MNOEX)
          ELSE
            IEQDOX = IEQO
            MNOEX = IM(IEQO,NMD)
            MDOEX = IM(IEQO,NMD)
          ENDIF
!
!---      Zero negative corrections for zero dissolved oil  ---
!
          IF( XMLO(2,N)/EPSL.LT.EPSL .AND. BLU(MDOEX)/EPSL.LT.EPSL )
     &      BLU(MDOEX) = 0.D+0
!
!---      Relax correction for sharp decrease in dissolved oil  ---
!
          IF( XMLO(2,N)+BLU(MDOEX).LT.1.D-1*XMLO(2,N) ) THEN
            XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MDOEX),0.D+0 )
          ELSE
            XMLO(2,N) = MAX( XMLO(2,N)+BLU(MDOEX),0.D+0 )
          ENDIF
          IF( XMLO(2,N).LT.1.D-12*XMLOX ) XMLO(2,N) = 0.D+0
!
!---      Zero negative corrections for zero dissolved air  ---
!
          IF( XMLA(2,N)/EPSL.LT.EPSL .AND. BLU(MAEX)/EPSL.LT.EPSL )
     &      BLU(MAEX) = 0.D+0
!
!---      Relax correction for sharp decrease in dissolved air  ---
!
          IF( XMLA(2,N)+BLU(MAEX).LT.1.D-1*XMLA(2,N) ) THEN
            XMLA(2,N) = MAX( XMLA(2,N)+6.D-1*BLU(MAEX),0.D+0 )
          ELSE
            XMLA(2,N) = MAX( XMLA(2,N)+BLU(MAEX),0.D+0 )
          ENDIF
!
!---    Aqueous w/ dissolved-oil, trapped-gas conditions
!       Unknowns: PL,SG,XMLO
!
!       Kinetic volatilization option
!       Unknowns: PL,SG,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.2 ) THEN
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) IERR = 1
          IF( BGL*(PG(2,N)-PL(2,N)-BLU(MWEX)).GT.ENPR ) THEN
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MWEX)
          ELSE
            PL(2,N) = PL(2,N) + BLU(MWEX)
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            MDOEX = IM(IEQDO,NMD)
            MNOEX = IM(IEQO,NMD)
            PN(2,N) = PN(2,N) + BLU(MNOEX)
          ELSE
            IEQDOX = IEQO
            MDOEX = IM(IEQO,NMD)
          ENDIF
!
!---      Zero negative corrections for zero dissolved oil  ---
!
          IF( XMLO(2,N)/EPSL.LT.EPSL .AND. BLU(MDOEX)/EPSL.LT.EPSL )
     &      BLU(MDOEX) = 0.D+0
          IF( ABS(RSDL(IEQDOX,N))/EPSL.GT.EPSL ) THEN
            IF( XMLO(2,N)+BLU(MDOEX).LT.1.D-1*XMLO(2,N) ) THEN
              XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MDOEX),0.D+0 )
            ELSE
              XMLO(2,N) = MAX( XMLO(2,N)+BLU(MDOEX),0.D+0 )
            ENDIF
          ENDIF
          IF( XMLO(2,N).LT.1.D-12*XMLOX ) XMLO(2,N) = 0.D+0
          IF( SG(2,N)+BLU(MAEX).LT.1.D-1*SG(2,N) ) THEN
            SG(2,N) = MAX( SG(2,N)+6.D-1*BLU(MAEX),0.D+0 )
          ELSE
            SG(2,N) = MAX( SG(2,N)+BLU(MAEX),0.D+0 )
          ENDIF
!
!---    Aqueous w/ trapped-NAPL, dissolved-air conditions
!       Unknowns: PL,XMLA,SN
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,SN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.3 ) THEN
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) IERR = 1
          IF( BGL*(PG(2,N)-PL(2,N)-BLU(MWEX)).GT.ENPR ) THEN
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MWEX)
          ELSE
            PL(2,N) = PL(2,N) + BLU(MWEX)
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            MNOEX = IM(IEQO,NMD)
            MDOEX = IM(IEQDO,NMD)
!
!---        Zero negative corrections for zero dissolved oil  ---
!
            IF( XMLO(2,N)/EPSL.LT.EPSL .AND. BLU(MDOEX)/EPSL.LT.EPSL )
     &        BLU(MDOEX) = 0.D+0
            IF( ABS(RSDL(IEQDOX,N))/EPSL.GT.EPSL ) THEN
              IF( XMLO(2,N)+BLU(MDOEX).LT.1.D-1*XMLO(2,N) ) THEN
                XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MDOEX),0.D+0 )
              ELSE
                XMLO(2,N) = MAX( XMLO(2,N)+BLU(MDOEX),0.D+0 )
              ENDIF
            ENDIF
            IF( XMLO(2,N).LT.1.D-12*XMLOX ) XMLO(2,N) = 0.D+0
          ELSE
            MNOEX = IM(IEQO,NMD)
          ENDIF
          IF( SN(2,N)+BLU(MNOEX).LT.1.D-1*SN(2,N) ) THEN
            SN(2,N) = MAX( SN(2,N)+6.D-1*BLU(MNOEX),0.D+0 )
          ELSEIF( SN(2,N)+BLU(MNOEX).GT.(1.D+0-SCHR(4,IZN)) ) THEN
            SN(2,N) = MIN(SN(2,N)+6.D-1*BLU(MNOEX),(1.D+0-SCHR(4,IZN)))
          ELSE
            SN(2,N) = MAX( SN(2,N)+BLU(MNOEX),0.D+0 )
            SN(2,N) = MIN( SN(2,N),(1.D+0-SCHR(4,IZN)) )
          ENDIF
          IF( SN(2,N).LT.EPSL ) SN(2,N) = 0.D+0
!
!---      Zero negative corrections for zero dissolved air  ---
!
          IF( XMLA(2,N)/EPSL.LT.EPSL .AND. BLU(MAEX)/EPSL.LT.EPSL )
     &      BLU(MAEX) = 0.D+0
          IF( XMLA(2,N)+BLU(MAEX).LT.1.D-1*XMLA(2,N) ) THEN
            XMLA(2,N) = MAX( XMLA(2,N)+6.D-1*BLU(MAEX),0.D+0 )
          ELSE
            XMLA(2,N) = MAX( XMLA(2,N)+BLU(MAEX),0.D+0 )
          ENDIF
!
!---    Aqueous w/ trapped-NAPL, trapped-gas conditions
!       Unknowns: PL,SG,SN
!
!       Kinetic volatilization option
!       Unknowns: PL,SG,SN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.4 ) THEN
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) IERR = 1
          IF( BGL*(PG(2,N)-PL(2,N)-BLU(MWEX)).GT.ENPR ) THEN
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MWEX)
          ELSE
            PL(2,N) = PL(2,N) + BLU(MWEX)
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            MNOEX = IM(IEQO,NMD)
            MDOEX = IM(IEQDO,NMD)
!
!---        Zero negative corrections for zero dissolved oil  ---
!
            IF( XMLO(2,N)/EPSL.LT.EPSL .AND. BLU(MDOEX)/EPSL.LT.EPSL )
     &        BLU(MDOEX) = 0.D+0
            IF( ABS(RSDL(IEQDOX,N))/EPSL.GT.EPSL ) THEN
              IF( XMLO(2,N)+BLU(MDOEX).LT.1.D-1*XMLO(2,N) ) THEN
                XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MDOEX),0.D+0 )
              ELSE
                XMLO(2,N) = MAX( XMLO(2,N)+BLU(MDOEX),0.D+0 )
              ENDIF
            ENDIF
            IF( XMLO(2,N).LT.1.D-12*XMLOX ) XMLO(2,N) = 0.D+0
          ELSE
            MNOEX = IM(IEQO,NMD)
          ENDIF
          IF( SN(2,N)+BLU(MNOEX).LT.1.D-1*SN(2,N) ) THEN
            SN(2,N) = MAX( SN(2,N)+6.D-1*BLU(MNOEX),0.D+0 )
          ELSEIF( SN(2,N)+BLU(MNOEX).GT.(1.D+0-SCHR(4,IZN)) ) THEN
            SN(2,N) = MIN(SN(2,N)+6.D-1*BLU(MNOEX),(1.D+0-SCHR(4,IZN)))
          ELSE
            SN(2,N) = MAX( SN(2,N)+BLU(MNOEX),0.D+0 )
            SN(2,N) = MIN( SN(2,N),(1.D+0-SCHR(4,IZN)) )
          ENDIF
          IF( SN(2,N).LT.EPSL ) SN(2,N) = 0.D+0
          IF( SG(2,N)+BLU(MAEX).LT.1.D-1*SG(2,N) ) THEN
            SG(2,N) = MAX( SG(2,N)+6.D-1*BLU(MAEX),0.D+0 )
          ELSE
            SG(2,N) = MAX( SG(2,N)+BLU(MAEX),0.D+0 )
          ENDIF
          IF( SG(2,N).LT.EPSL ) SG(2,N) = 0.D+0
!
!---    Aqueous-gas w/ dissolved-oil conditions
!       Unknowns:  PL,PG,XMLO
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.5 ) THEN
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) IERR = 1
          IF( ABS(BLU(MAEX)).GT.1.D+6 ) IERR = 1
          IF( BGL*(PG(2,N)+BLU(MAEX)-PL(2,N)-BLU(MWEX)).LT.ENPR ) THEN
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MWEX)
            PG(2,N) = PG(2,N) + 6.D-1*BLU(MAEX)
          ELSE
            PL(2,N) = PL(2,N) + BLU(MWEX)
            PG(2,N) = PG(2,N) + BLU(MAEX)
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            MDOEX = IM(IEQDO,NMD)
            MNOEX = IM(IEQO,NMD)
            PN(2,N) = PN(2,N) + BLU(MNOEX)
          ELSE
            IEQDOX = IEQO
            MDOEX = IM(IEQO,NMD)
          ENDIF
!
!---      Zero negative corrections for zero dissolved oil  ---
!
          IF( XMLO(2,N)/EPSL.LT.EPSL .AND. BLU(MDOEX)/EPSL.LT.EPSL )
     &      BLU(MDOEX) = 0.D+0
          IF( ABS(RSDL(IEQDOX,N))/EPSL.GT.EPSL ) THEN
            IF( XMLO(2,N)+BLU(MDOEX).LT.1.D-1*XMLO(2,N) ) THEN
              XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MDOEX),0.D+0 )
            ELSE
              XMLO(2,N) = XMLO(2,N)+BLU(MDOEX)
            ENDIF
          ENDIF
          IF( XMLO(2,N).LT.1.D-12*XMLOX ) XMLO(2,N) = 0.D+0
!
!---    Aqueous-gas w/ trapped-NAPL conditions
!       Unknowns: PL,PG,SN
!
!       Kinetic volatilization option
!       Unknowns: PL,PG,SN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.6 ) THEN
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) IERR = 1
          IF( ABS(BLU(MAEX)).GT.1.D+6 ) IERR = 1
          IF( BGL*(PG(2,N)+BLU(MAEX)-PL(2,N)-BLU(MWEX)).LT.ENPR ) THEN
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MWEX)
            PG(2,N) = PG(2,N) + 6.D-1*BLU(MAEX)
          ELSE
            PL(2,N) = PL(2,N) + BLU(MWEX)
            PG(2,N) = PG(2,N) + BLU(MAEX)
          ENDIF
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            MNOEX = IM(IEQO,NMD)
            MDOEX = IM(IEQDO,NMD)
!
!---        Zero negative corrections for zero dissolved oil  ---
!
            IF( XMLO(2,N)/EPSL.LT.EPSL .AND. BLU(MDOEX)/EPSL.LT.EPSL )
     &        BLU(MDOEX) = 0.D+0
            IF( ABS(RSDL(IEQDOX,N))/EPSL.GT.EPSL ) THEN
              IF( XMLO(2,N)+BLU(MDOEX).LT.1.D-1*XMLO(2,N) ) THEN
                XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MDOEX),0.D+0 )
              ELSE
                XMLO(2,N) = MAX( XMLO(2,N)+BLU(MDOEX),0.D+0 )
              ENDIF
            ENDIF
            IF( XMLO(2,N).LT.1.D-12*XMLOX ) XMLO(2,N) = 0.D+0
          ELSE
            MNOEX = IM(IEQO,NMD)
          ENDIF
          IF( SN(2,N)+BLU(MNOEX).LT.1.D-1*SN(2,N) ) THEN
            SN(2,N) = MAX( SN(2,N)+6.D-1*BLU(MNOEX),0.D+0 )
          ELSEIF( SN(2,N)+BLU(MNOEX).GT.(1.D+0-SCHR(4,IZN)) ) THEN
            SN(2,N) = MIN(SN(2,N)+6.D-1*BLU(MNOEX),(1.D+0-SCHR(4,IZN)))
          ELSE
            SN(2,N) = MAX( SN(2,N)+BLU(MNOEX),0.D+0 )
            SN(2,N) = MIN( SN(2,N),(1.D+0-SCHR(4,IZN)) )
          ENDIF
          IF( SN(2,N).LT.EPSL ) SN(2,N) = 0.D+0
!
!---    Aqueous-NAPL w/ dissolved-air conditions
!       Unknowns: PL,XMLA,PN
!
!       Kinetic volatilization option
!       Unknowns: PL,XMLA,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.7 ) THEN
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            MNOEX = IM(IEQO,NMD)
            MDOEX = IM(IEQDO,NMD)
!
!---        Zero negative corrections for zero dissolved oil  ---
!
            IF( XMLO(2,N)/EPSL.LT.EPSL .AND. BLU(MDOEX)/EPSL.LT.EPSL )
     &        BLU(MDOEX) = 0.D+0
            IF( ABS(RSDL(IEQDOX,N))/EPSL.GT.EPSL ) THEN
              IF( XMLO(2,N)+BLU(MDOEX).LT.1.D-1*XMLO(2,N) ) THEN
                XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MDOEX),0.D+0 )
              ELSE
                XMLO(2,N) = MAX( XMLO(2,N)+BLU(MDOEX),0.D+0 )
              ENDIF
            ENDIF
            IF( XMLO(2,N).LT.1.D-12*XMLOX ) XMLO(2,N) = 0.D+0
          ELSE
            MNOEX = IM(IEQO,NMD)
          ENDIF
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) BLU(MWEX) = 6.D-1*BLU(MWEX)
          IF( ABS(BLU(MNOEX)).GT.1.D+6 ) BLU(MNOEX) = 6.D-1*BLU(MNOEX)
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) IERR = 1
          IF( ABS(BLU(MNOEX)).GT.1.D+6 ) IERR = 1
          PLX = PL(2,N) + BLU(MWEX)
          PNX = PN(2,N) + BLU(MNOEX)
          PNC = MAX( (BNL*PLX+BGN*PG(2,N))/(BNL+BGN),
     &      (PLX*BNL+ENPR)/BNL )
          IF( PNX.LT.PNC ) THEN
            PN(2,N) = PN(2,N) + 6.D-1*BLU(MNOEX)
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MWEX)
          ELSE
            PN(2,N) = PNX
            PL(2,N) = PLX
          ENDIF
          PL(2,N) = MAX( PL(2,N),PN(2,N)-1.D+5 )
!
!---      Zero negative corrections for zero dissolved air  ---
!
          IF( XMLA(2,N)/EPSL.LT.EPSL .AND. BLU(MAEX)/EPSL.LT.EPSL )
     &      BLU(MAEX) = 0.D+0
          IF( XMLA(2,N)+BLU(MAEX).LT.1.D-1*XMLA(2,N) ) THEN
            XMLA(2,N) = MAX( XMLA(2,N)+6.D-1*BLU(MAEX),0.D+0 )
          ELSE
            XMLA(2,N) = MAX( XMLA(2,N)+BLU(MAEX),0.D+0 )
          ENDIF
!
!---    Aqueous-NAPL w/ trapped-gas conditions
!       Unknowns: PL,SG,PN
!
!       Kinetic volatilization option
!       Unknowns: PL,SG,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.8 ) THEN
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            MNOEX = IM(IEQO,NMD)
            MDOEX = IM(IEQDO,NMD)
!
!---        Zero negative corrections for zero dissolved oil  ---
!
            IF( XMLO(2,N)/EPSL.LT.EPSL .AND. BLU(MDOEX)/EPSL.LT.EPSL )
     &        BLU(MDOEX) = 0.D+0
            IF( ABS(RSDL(IEQDOX,N))/EPSL.GT.EPSL ) THEN
              IF( XMLO(2,N)+BLU(MDOEX).LT.1.D-1*XMLO(2,N) ) THEN
                XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MDOEX),0.D+0 )
              ELSE
                XMLO(2,N) = MAX( XMLO(2,N)+BLU(MDOEX),0.D+0 )
              ENDIF
            ENDIF
            IF( XMLO(2,N).LT.1.D-12*XMLOX ) XMLO(2,N) = 0.D+0
          ELSE
            MNOEX = IM(IEQO,NMD)
          ENDIF
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) BLU(MWEX) = 6.D-1*BLU(MWEX)
          IF( ABS(BLU(MNOEX)).GT.1.D+6 ) BLU(MNOEX) = 6.D-1*BLU(MNOEX)
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) IERR = 1
          IF( ABS(BLU(MNOEX)).GT.1.D+6 ) IERR = 1
          PLX = PL(2,N) + BLU(MWEX)
          PNX = PN(2,N) + BLU(MNOEX)
          PNC = MAX( (BNL*PLX+BGN*PG(2,N))/(BNL+BGN),
     &      (PLX*BNL+ENPR)/BNL )
          IF( PNX.LT.PNC ) THEN
            PN(2,N) = PN(2,N) + 6.D-1*BLU(MNOEX)
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MWEX)
          ELSE
            PN(2,N) = PNX
            PL(2,N) = PLX
          ENDIF
          PL(2,N) = MAX( PL(2,N),PN(2,N)-1.D+5 )
          IF( SG(2,N)+BLU(MAEX).LT.1.D-1*SG(2,N) ) THEN
            SG(2,N) = MAX( SG(2,N)+6.D-1*BLU(MAEX),0.D+0 )
          ELSE
            SG(2,N) = MAX( SG(2,N)+BLU(MAEX),0.D+0 )
          ENDIF
          IF( SG(2,N).LT.EPSL ) SG(2,N) = 0.D+0
!
!---    Aqueous-NAPL-gas conditions
!       Unknowns: PL,PG,PN
!
!       Kinetic volatilization option
!       Unknowns: PL,PG,PN,XMLO  ---
!
        ELSEIF( NPHAZ(2,N).EQ.9 ) THEN
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            IEQDOX = IEQDO
            MNOEX = IM(IEQO,NMD)
            MDOEX = IM(IEQDO,NMD)
!
!---        Zero negative corrections for zero dissolved oil  ---
!
            IF( XMLO(2,N)/EPSL.LT.EPSL .AND. BLU(MDOEX)/EPSL.LT.EPSL )
     &        BLU(MDOEX) = 0.D+0
            IF( ABS(RSDL(IEQDOX,N))/EPSL.GT.EPSL ) THEN
              IF( XMLO(2,N)+BLU(MDOEX).LT.1.D-1*XMLO(2,N) ) THEN
                XMLO(2,N) = MAX( XMLO(2,N)+6.D-1*BLU(MDOEX),0.D+0 )
              ELSE
                XMLO(2,N) = MAX( XMLO(2,N)+BLU(MDOEX),0.D+0 )
              ENDIF
            ENDIF
            IF( XMLO(2,N).LT.1.D-12*XMLOX ) XMLO(2,N) = 0.D+0
          ELSE
            MNOEX = IM(IEQO,NMD)
          ENDIF
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) BLU(MWEX) = 6.D-1*BLU(MWEX)
          IF( ABS(BLU(MNOEX)).GT.1.D+6 ) BLU(MNOEX) = 6.D-1*BLU(MNOEX)
          IF( ABS(BLU(MAEX)).GT.1.D+6 ) BLU(MAEX) = 6.D-1*BLU(MAEX)
          IF( ABS(BLU(MWEX)).GT.1.D+6 ) IERR = 1
          IF( ABS(BLU(MNOEX)).GT.1.D+6 ) IERR = 1
          IF( ABS(BLU(MAEX)).GT.1.D+6 ) IERR = 1
          PLX = PL(2,N) + BLU(MWEX)
          PNX = PN(2,N) + BLU(MNOEX)
          PGX = PG(2,N) + BLU(MAEX)
          PNC = MAX( (BNL*PLX+BGN*PG(2,N))/(BNL+BGN),
     &      (PLX*BNL+ENPR)/BNL )
          IF( PNX.LT.PNC ) THEN
            PN(2,N) = PN(2,N) + 6.D-1*BLU(MNOEX)
            PL(2,N) = PL(2,N) + 6.D-1*BLU(MWEX)
            PG(2,N) = PG(2,N) + 6.D-1*BLU(MAEX)
          ELSE
            PN(2,N) = PNX
            PL(2,N) = PLX
            PG(2,N) = PGX
          ENDIF
!          PL(2,N) = MAX( PL(2,N),PN(2,N)-1.D+5 )
        ENDIF
!
!---  Reduce time step for excessive changes in primary variables   ---
!
        IF( IERR.EQ.1 ) THEN
          ICNV = 1
          WRITE(ISC,'(10X,A)') '---  Excessive Primary Variable Change
     &---'
          WRITE(IWR,'(10X,A)') '---  Excessive Primary Variable Change
     &---'
          WRITE(ISC,'(2(A,I6))') '    Node = ',N,' Phase Condition = ',
     &      NPHAZ(2,N)
          WRITE(IWR,'(2(A,I6))') '    Node = ',N,' Phase Condition = ',
     &      NPHAZ(2,N)
          WRITE(ISC,'(2(A,1PE12.5))')
     &      '    Water Pressure = ',PL(2,N)+PATM,' Change = ',BLU(MWEX)
          WRITE(IWR,'(2(A,1PE12.5))')
     &      '    Water Pressure = ',PL(2,N)+PATM,' Change = ',BLU(MWEX)
          WRITE(ISC,'(2(A,1PE12.5))')
     &      '    Gas Pressure = ',PG(2,N)+PATM,' Change = ',BLU(MAEX)
          WRITE(IWR,'(2(A,1PE12.5))')
     &      '    Gas Pressure = ',PG(2,N)+PATM,' Change = ',BLU(MAEX)
          WRITE(ISC,'(2(A,1PE12.5))')
     &      '    NAPL Pressure = ',PN(2,N)+PATM,' Change = ',BLU(MNOEX)
          WRITE(IWR,'(2(A,1PE12.5))')
     &      '    NAPL Pressure = ',PN(2,N)+PATM,' Change = ',BLU(MNOEX)
!
!---      Kinetic volatilization option  ---
!
          IF( ISLC(46).GE.1 ) THEN
            WRITE(ISC,'(2(A,1PE12.5))')
     &        '    Aqueous Oil Mole Fraction = ',XMLO(2,N),
     &        ' Change = ',BLU(MDOEX)
            WRITE(IWR,'(2(A,1PE12.5))')
     &        '    Aqueous Oil Mole Fraction = ',XMLO(2,N),
     &        ' Change = ',BLU(MDOEX)
          ENDIF
          GOTO 300
        ENDIF
  200 CONTINUE
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
          DO 400 N = 1,NFLD
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            XMLA(2,N) = XMLA(1,N)
            XMLO(2,N) = XMLO(1,N)
            SN(2,N) = SN(1,N)
            SG(2,N) = SG(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  400     CONTINUE
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
          DO 410 N = 1,NFLD
            PL(2,N) = PL(1,N)
            PG(2,N) = PG(1,N)
            PN(2,N) = PN(1,N)
            XMLA(2,N) = XMLA(1,N)
            XMLO(2,N) = XMLO(1,N)
            SN(2,N) = SN(1,N)
            SG(2,N) = SG(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  410     CONTINUE
          NSTEP = NSTEP - 1
          TM = TM - DT
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit Exceeded
     & ---'
          ICNV = 4
        ENDIF
      ENDIF
!
!---  Reset subroutine name list  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDT5 group
!
      RETURN
      END
